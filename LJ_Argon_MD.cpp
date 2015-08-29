/*******************************************/
/*  Teapot を ArcBall で操作      前田 稔  */
/*  World, View, Projection, Color         */
/*******************************************/
#include "DXUT.h"
#include "SDKmisc.h"
#include "DXUTShapes.h"
#include "DXUTcamera.h"
#include "moleculardynamics/Ar_moleculardynamics.h"
#include "utility/utility.h"
#include <array>                // for std::array
#include <memory>               // for std::unique_ptr
#include <vector>               // std::vector

//! A global variable (constant).
/*!
    インデックスバッファの個数
*/
static auto const NUMINDEXBUFFER = 16U;

//! A global variable (constant).
/*!
    頂点バッファの個数
*/
static auto const NUMMESH = 2U;

//! A global variable (constant).
/*!
    頂点バッファの個数
*/
static auto const NUMVERTEXBUFFER = 8U;

//! A global variable.
/*!
    バッファー リソース
*/
D3D10_BUFFER_DESC bd;

//! A global variable.
/*!
    メッシュへのスマートポインタが格納された可変長配列
*/
std::vector<std::unique_ptr<ID3DX10Mesh, utility::Safe_Release<ID3DX10Mesh>>> pmeshvec(NUMMESH);

//! A global variable.
/*!
	エフェクト＝シェーダプログラムを読ませるところ
*/
std::unique_ptr<ID3D10Effect, utility::Safe_Release<ID3D10Effect>> pEffect;

//! A global variable.
/*!
	インデックスバッファ
*/
std::unique_ptr<ID3D10Buffer, utility::Safe_Release<ID3D10Buffer>> pIndexBuffer;

//! A global variable.
/*!
	入力レイアウト インターフェイス
*/
std::unique_ptr<ID3D10InputLayout, utility::Safe_Release<ID3D10InputLayout>> pInputLayout;

//! A global variable.
/*!
	頂点バッファ
*/
std::unique_ptr<ID3D10Buffer, utility::Safe_Release<ID3D10Buffer>> pVertexBuffer;

ID3D10EffectTechnique*      g_pRender = nullptr;
ID3D10EffectMatrixVariable* g_pWorldVariable = nullptr;
ID3D10EffectMatrixVariable* g_pViewVariable = nullptr;
ID3D10EffectMatrixVariable* g_pProjectionVariable = nullptr;
ID3D10EffectVectorVariable* g_pColorVariable = nullptr;
D3DXMATRIX                  g_View;
D3DXMATRIX                  g_Projection;
D3DXVECTOR4 g_Colors[2] = 
{
    D3DXVECTOR4( 1.0f, 1.0f, 1.0f, 1.0f ),
    D3DXVECTOR4( 1.0f, 0.3f, 0.3f, 1.0f ),
};
CModelViewerCamera          g_Camera;

moleculardynamics::Ar_moleculardynamics armd;

//--------------------------------------------------------------------------------------
// Structures
//--------------------------------------------------------------------------------------
struct SimpleVertex
{
    D3DXVECTOR3 Pos;
};

//--------------------------------------------------------------------------------------
// Render the scene using the D3D10 device
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D10FrameRender( ID3D10Device* pd3dDevice, double fTime, float fElapsedTime, void* pUserContext )
{
    armd.Calc_Forces();
    armd.Move_Atoms();

    // Clear render target and the depth stencil 
    std::array<float, 4> const ClearColor = { 0.176f, 0.196f, 0.667f, 0.0f };
    pd3dDevice->ClearRenderTargetView(DXUTGetD3D10RenderTargetView(), ClearColor.data());
    pd3dDevice->ClearDepthStencilView(DXUTGetD3D10DepthStencilView(), D3D10_CLEAR_DEPTH, 1.0, 0);
    
    // Update variables
    g_pWorldVariable->SetMatrix(reinterpret_cast<float *>(const_cast<D3DXMATRIX *>(&(*g_Camera.GetWorldMatrix()))));
    g_pViewVariable->SetMatrix(reinterpret_cast<float *>(const_cast<D3DXMATRIX *>(&(*g_Camera.GetViewMatrix()))));
    g_pProjectionVariable->SetMatrix(reinterpret_cast<float *>(const_cast<D3DXMATRIX *>(&(*g_Camera.GetProjMatrix()))));
    
    D3D10_TECHNIQUE_DESC techDesc;
    g_pRender->GetDesc(&techDesc);

    g_pColorVariable->SetFloatVector(reinterpret_cast<float *>(&g_Colors[0]));

    // Set vertex buffer
    auto const stride = sizeof(SimpleVertex);
    auto const offset = 0U;
    auto const pVertexBuffertmp = pVertexBuffer.get();
    pd3dDevice->IASetVertexBuffers(0, 1, &pVertexBuffertmp, &stride, &offset);

    // Set index buffer
    pd3dDevice->IASetIndexBuffer(pIndexBuffer.get(), DXGI_FORMAT_R32_UINT, 0);

    // Set primitive topology
    pd3dDevice->IASetPrimitiveTopology(D3D10_PRIMITIVE_TOPOLOGY_LINESTRIP);

    for (auto p = 0U; p < techDesc.Passes; p++)
    {
        g_pRender->GetPassByIndex(p)->Apply(0);
        pd3dDevice->DrawIndexed(NUMINDEXBUFFER, 0, 0);
    }

    auto i = 0;
    for (auto & pmesh : pmeshvec) {
        g_pColorVariable->SetFloatVector(reinterpret_cast<float *>(&g_Colors[1]));

        D3DXMATRIX  World;
        D3DXMatrixTranslation(&World, 1.0f * i, 0.0f, 0.0f);
        D3DXMatrixMultiply(&World, &(*g_Camera.GetWorldMatrix()), &World);

        // Update variables
        g_pWorldVariable->SetMatrix(World);

        UINT NumSubsets;
        pmesh->GetAttributeTable(nullptr, &NumSubsets);

        for (auto p = 0U; p < techDesc.Passes; p++)
        {
            g_pRender->GetPassByIndex(p)->Apply(0);
            for (auto s = 0U; s < NumSubsets; s++)
            {
                pmesh->DrawSubset(s);
            }
        }

        i++;
    }
}

//--------------------------------------------------------------------------------------
// Reject any D3D10 devices that aren't acceptable by returning false
//--------------------------------------------------------------------------------------
bool CALLBACK IsD3D10DeviceAcceptable( UINT Adapter, UINT Output, D3D10_DRIVER_TYPE DeviceType, DXGI_FORMAT BackBufferFormat, bool bWindowed, void* pUserContext )
{
    return true;
}

//--------------------------------------------------------------------------------------
// Called right before creating a D3D9 or D3D10 device, allowing the app to modify the device settings as needed
//--------------------------------------------------------------------------------------
bool CALLBACK ModifyDeviceSettings( DXUTDeviceSettings* pDeviceSettings, void* pUserContext )
{
    return true;
}

//--------------------------------------------------------------------------------------
// Create any D3D10 resources that aren't dependant on the back buffer
//--------------------------------------------------------------------------------------
HRESULT CALLBACK OnD3D10CreateDevice( ID3D10Device* pd3dDevice, const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc,
                                      void* pUserContext )
{
    HRESULT hr = S_OK;

    // Find the D3DX effect file
    std::array<WCHAR, MAX_PATH> str;
    V_RETURN( DXUTFindDXSDKMediaFileCch( str.data(), MAX_PATH, L"Normal.fx" ) );
    DWORD dwShaderFlags = D3D10_SHADER_ENABLE_STRICTNESS;
#if defined( DEBUG ) || defined( _DEBUG )
    dwShaderFlags |= D3D10_SHADER_DEBUG;
#endif

	ID3D10Effect * pEffecttmp;
	utility::v_return(
		D3DX10CreateEffectFromFile(
			str.data(),
			nullptr,
			nullptr,
			"fx_4_0",
			dwShaderFlags,
			0,
			pd3dDevice,
			nullptr,
			nullptr,
			&pEffecttmp,
			nullptr,
			nullptr ) );
	pEffect.reset(pEffecttmp);

    // Obtain the technique
    g_pRender = pEffect->GetTechniqueByName( "Render" );
    
	// Obtain the variables
    g_pWorldVariable = pEffect->GetVariableByName( "World" )->AsMatrix();
    g_pViewVariable = pEffect->GetVariableByName( "View" )->AsMatrix();
    g_pProjectionVariable = pEffect->GetVariableByName( "Projection" )->AsMatrix();
    g_pColorVariable = pEffect->GetVariableByName( "Color" )->AsVector();

    // Create an input layout
    D3D10_INPUT_ELEMENT_DESC const layout[] =
    {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0,  D3D10_INPUT_PER_VERTEX_DATA, 0 },
        { "NORMAL",   0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D10_INPUT_PER_VERTEX_DATA, 0 },
    };
    
	D3D10_PASS_DESC PassDesc;
    g_pRender->GetPassByIndex( 0 )->GetDesc( &PassDesc );

	ID3D10InputLayout * pInputLayouttmp;
    utility::v_return(
		pd3dDevice->CreateInputLayout(
			layout,
			sizeof(layout) / sizeof(layout[0]), 
            PassDesc.pIAInputSignature,
			PassDesc.IAInputSignatureSize, &pInputLayouttmp ) );
	pInputLayout.reset(pInputLayouttmp);

	pd3dDevice->IASetInputLayout(pInputLayout.get());

    for (auto & pmesh : pmeshvec) {
        ID3DX10Mesh * pmeshtmp;
        DXUTCreateSphere(pd3dDevice, 1.0f, 16, 16, &pmeshtmp);
        pmesh.reset(pmeshtmp);
    }

    // Create vertex buffer
    std::array<SimpleVertex, NUMVERTEXBUFFER> const vertices =
    {
        D3DXVECTOR3(-1.0f, 1.0f, -1.0f),
        D3DXVECTOR3(1.0f, 1.0f, -1.0f),
        D3DXVECTOR3(1.0f, 1.0f, 1.0f), 
        D3DXVECTOR3(-1.0f, 1.0f, 1.0f),

        D3DXVECTOR3(-1.0f, -1.0f, -1.0f), 
        D3DXVECTOR3(1.0f, -1.0f, -1.0f),
        D3DXVECTOR3(1.0f, -1.0f, 1.0f),
        D3DXVECTOR3(-1.0f, -1.0f, 1.0f),
    };

    bd.Usage = D3D10_USAGE_DEFAULT;
	bd.ByteWidth = sizeof(SimpleVertex) * NUMVERTEXBUFFER;
    bd.BindFlags = D3D10_BIND_VERTEX_BUFFER;
    bd.CPUAccessFlags = 0;
    bd.MiscFlags = 0;

    D3D10_SUBRESOURCE_DATA InitData;
    InitData.pSysMem = vertices.data();

	ID3D10Buffer * pVertexBuffertmp;
    utility::v_return(pd3dDevice->CreateBuffer(&bd, &InitData, &pVertexBuffertmp));
	pVertexBuffer.reset(pVertexBuffertmp);

    // Create index buffer
    // Create vertex buffer
    std::array<DWORD, NUMINDEXBUFFER> const indices =
    {
        0, 1, 2,
        3, 0, 4,

        5, 1, 2,
        6, 5, 4,

        7, 3, 7,
        6
    };

    bd.Usage = D3D10_USAGE_DEFAULT;
    bd.ByteWidth = sizeof(DWORD) * NUMINDEXBUFFER;
    bd.BindFlags = D3D10_BIND_INDEX_BUFFER;
    bd.CPUAccessFlags = 0;
    bd.MiscFlags = 0;
    InitData.pSysMem = indices.data();
	
	ID3D10Buffer * pIndexBuffertmp;
    utility::v_return(pd3dDevice->CreateBuffer(&bd, &InitData, &pIndexBuffertmp));
	pIndexBuffer.reset(pIndexBuffertmp);
	
    D3DXVECTOR3 vEye(0.0f, 4.0f, -4.0f);
    D3DXVECTOR3 vLook(0.0f, 0.0f, 0.0f);
    D3DXVECTOR3 const Up(0.0f, 1.0f, 0.0f);
    D3DXMatrixLookAtLH(&g_View, &vEye, &vLook, &Up);

    // Update Variables that never change
    g_pViewVariable->SetMatrix(reinterpret_cast<float *>(&g_View));

    g_Camera.SetViewParams(&vEye, &vLook);

    return S_OK;
}

//--------------------------------------------------------------------------------------
// Create any D3D10 resources that depend on the back buffer
//--------------------------------------------------------------------------------------
HRESULT CALLBACK OnD3D10ResizedSwapChain( ID3D10Device* pd3dDevice, IDXGISwapChain *pSwapChain, const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc, void* pUserContext )
{
    auto const fAspectRatio = static_cast<float>(pBackBufferSurfaceDesc->Width) /
        static_cast<float>(pBackBufferSurfaceDesc->Height);
    g_Camera.SetProjParams(D3DX_PI / 4, fAspectRatio, 0.1f, 1000.0f);
    g_Camera.SetWindow(pBackBufferSurfaceDesc->Width, pBackBufferSurfaceDesc->Height);

    return S_OK;
}

//--------------------------------------------------------------------------------------
// Handle updates to the scene.  This is called regardless of which D3D API is used
//--------------------------------------------------------------------------------------
void CALLBACK OnFrameMove( double fTime, float fElapsedTime, void* pUserContext )
{
    g_Camera.FrameMove(fElapsedTime);
}

//--------------------------------------------------------------------------------------
// Release D3D10 resources created in OnD3D10ResizedSwapChain 
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D10ReleasingSwapChain( void* pUserContext )
{
}

//--------------------------------------------------------------------------------------
// Release D3D10 resources created in OnD3D10CreateDevice 
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D10DestroyDevice( void* pUserContext )
{
    pEffect.reset();
	pInputLayout.reset();
    pIndexBuffer.reset();
    pVertexBuffer.reset();

    for (auto & pmesh : pmeshvec) {
        pmesh.reset();
    }
}

//--------------------------------------------------------------------------------------
// Handle messages to the application
//--------------------------------------------------------------------------------------
LRESULT CALLBACK MsgProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam, 
                          bool* pbNoFurtherProcessing, void* pUserContext )
{
    g_Camera.HandleMessages(hWnd, uMsg, wParam, lParam);

    return 0;
}

//--------------------------------------------------------------------------------------
// Handle key presses
//--------------------------------------------------------------------------------------
void CALLBACK OnKeyboard( UINT nChar, bool bKeyDown, bool bAltDown, void* pUserContext )
{
}

//--------------------------------------------------------------------------------------
// Handle mouse button presses
//--------------------------------------------------------------------------------------
void CALLBACK OnMouse( bool bLeftButtonDown, bool bRightButtonDown, bool bMiddleButtonDown, 
                       bool bSideButton1Down, bool bSideButton2Down, int nMouseWheelDelta, 
                       int xPos, int yPos, void* pUserContext )
{
}

//--------------------------------------------------------------------------------------
// Call if device was removed.  Return true to find a new device, false to quit
//--------------------------------------------------------------------------------------
bool CALLBACK OnDeviceRemoved( void* pUserContext )
{
    return true;
}

//--------------------------------------------------------------------------------------
// Initialize everything and go into a render loop
//--------------------------------------------------------------------------------------
int WINAPI wWinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance, LPWSTR lpCmdLine, int nCmdShow )
{
    // Enable run-time memory check for debug builds.
#if defined(DEBUG) | defined(_DEBUG)
    _CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif

    // Set general DXUT callbacks
    DXUTSetCallbackFrameMove( OnFrameMove );
    DXUTSetCallbackKeyboard( OnKeyboard );
    DXUTSetCallbackMouse( OnMouse );
    DXUTSetCallbackMsgProc( MsgProc );
    DXUTSetCallbackDeviceChanging( ModifyDeviceSettings );
    DXUTSetCallbackDeviceRemoved( OnDeviceRemoved );

    // Set the D3D10 DXUT callbacks. Remove these sets if the app doesn't need to support D3D10
    DXUTSetCallbackD3D10DeviceAcceptable( IsD3D10DeviceAcceptable );
    DXUTSetCallbackD3D10DeviceCreated( OnD3D10CreateDevice );
    DXUTSetCallbackD3D10SwapChainResized( OnD3D10ResizedSwapChain );
    DXUTSetCallbackD3D10FrameRender( OnD3D10FrameRender );
    DXUTSetCallbackD3D10SwapChainReleasing( OnD3D10ReleasingSwapChain );
    DXUTSetCallbackD3D10DeviceDestroyed( OnD3D10DestroyDevice );

    DXUTInit( true, true, nullptr ); // Parse the command line, show msgboxes on error, no extra command line params
    DXUTSetCursorSettings( true, true ); // Show the cursor and clip it when in full screen
    DXUTCreateWindow( L"Teapot ArcBall" );
    DXUTCreateDevice( true, 640, 480 );  
    DXUTMainLoop(); // Enter into the DXUT render loop

    return DXUTGetExitCode();
}

