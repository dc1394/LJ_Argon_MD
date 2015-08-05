/*******************************************/
/*  Teapot を ArcBall で操作      前田 稔  */
/*  World, View, Projection, Color         */
/*******************************************/
#include "DXUT.h"
#include "SDKmisc.h"
#include "DXUTShapes.h"
#include "DXUTcamera.h"

//! A global variable.
/*!
    バッファー リソース
*/
D3D10_BUFFER_DESC g_bd;

ID3DX10Mesh*                g_pMesh = NULL;
ID3D10Buffer*               g_pVertexBuffer = NULL;
ID3D10Buffer*               g_pIndexBuffer = NULL;
ID3D10InputLayout*          g_pLayout = NULL;
ID3D10Effect*               g_pEffect = NULL;
ID3D10EffectTechnique*      g_pRender = NULL;
ID3D10EffectMatrixVariable* g_pWorldVariable = NULL;
ID3D10EffectMatrixVariable* g_pViewVariable = NULL;
ID3D10EffectMatrixVariable* g_pProjectionVariable = NULL;
ID3D10EffectVectorVariable* g_pColorVariable = NULL;
D3DXMATRIX                  g_View;
D3DXMATRIX                  g_Projection;
D3DXVECTOR4 g_Colors[2] = 
{
    D3DXVECTOR4( 1.0f, 1.0f, 1.0f, 1.0f ),
    D3DXVECTOR4( 1.0f, 0.3f, 0.3f, 1.0f ),
};
CModelViewerCamera          g_Camera;

//--------------------------------------------------------------------------------------
// Structures
//--------------------------------------------------------------------------------------
struct SimpleVertex
{
    D3DXVECTOR3 Pos;
    D3DXVECTOR2 Tex;
};

//--------------------------------------------------------------------------------------
// Render the scene using the D3D10 device
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D10FrameRender( ID3D10Device* pd3dDevice, double fTime, float fElapsedTime, void* pUserContext )
{
    // Clear render target and the depth stencil 
    float ClearColor[4] = { 0.176f, 0.196f, 0.667f, 0.0f };
    pd3dDevice->ClearRenderTargetView(DXUTGetD3D10RenderTargetView(), ClearColor);
    pd3dDevice->ClearDepthStencilView(DXUTGetD3D10DepthStencilView(), D3D10_CLEAR_DEPTH, 1.0, 0);

    // Update variables
    g_pWorldVariable->SetMatrix((float*)&(*g_Camera.GetWorldMatrix()));
    g_pViewVariable->SetMatrix((float*)&(*g_Camera.GetViewMatrix()));
    g_pProjectionVariable->SetMatrix((float*)&(*g_Camera.GetProjMatrix()));

    g_pColorVariable->SetFloatVector((float*)&g_Colors[0]);

    D3D10_TECHNIQUE_DESC techDesc;
    g_pRender->GetDesc(&techDesc);

    UINT NumSubsets;
    g_pMesh->GetAttributeTable(NULL, &NumSubsets);

    pd3dDevice->IASetInputLayout(g_pLayout);

    for (UINT p = 0; p<techDesc.Passes; p++)
    {
        g_pRender->GetPassByIndex(p)->Apply(0);
        for (UINT s = 0; s<NumSubsets; s++)
        {
            g_pMesh->DrawSubset(s);
        }
    }

    g_pColorVariable->SetFloatVector((float*)&g_Colors[1]);
    
    // Set vertex buffer
    UINT stride = sizeof(SimpleVertex);
    UINT offset = 0;
    pd3dDevice->IASetVertexBuffers(0, 1, &g_pVertexBuffer, &stride, &offset);

    // Set index buffer
    pd3dDevice->IASetIndexBuffer(g_pIndexBuffer, DXGI_FORMAT_R32_UINT, 0);

    // Set primitive topology
    pd3dDevice->IASetPrimitiveTopology(D3D10_PRIMITIVE_TOPOLOGY_LINESTRIP);
        
    for (UINT p = 0; p<techDesc.Passes; p++)
    {
        g_pRender->GetPassByIndex(p)->Apply(0);
        pd3dDevice->DrawIndexed(9, 0, 0);
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
    WCHAR str[MAX_PATH];
    V_RETURN( DXUTFindDXSDKMediaFileCch( str, MAX_PATH, L"Normal.fx" ) );
    DWORD dwShaderFlags = D3D10_SHADER_ENABLE_STRICTNESS;
    #if defined( DEBUG ) || defined( _DEBUG )
    dwShaderFlags |= D3D10_SHADER_DEBUG;
    #endif
    V_RETURN( D3DX10CreateEffectFromFile( str, NULL, NULL, "fx_4_0", dwShaderFlags, 0, pd3dDevice, NULL, NULL, &g_pEffect, NULL, NULL ) );

    // Obtain the technique
    g_pRender = g_pEffect->GetTechniqueByName( "Render" );
    // Obtain the variables
    g_pWorldVariable = g_pEffect->GetVariableByName( "World" )->AsMatrix();
    g_pViewVariable = g_pEffect->GetVariableByName( "View" )->AsMatrix();
    g_pProjectionVariable = g_pEffect->GetVariableByName( "Projection" )->AsMatrix();
    g_pColorVariable = g_pEffect->GetVariableByName( "Color" )->AsVector();

    // Create an input layout
    const D3D10_INPUT_ELEMENT_DESC layout[] =
    {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0,  D3D10_INPUT_PER_VERTEX_DATA, 0 },
        { "NORMAL",   0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D10_INPUT_PER_VERTEX_DATA, 0 },
    };
    D3D10_PASS_DESC PassDesc;
    g_pRender->GetPassByIndex( 0 )->GetDesc( &PassDesc );
    V_RETURN( pd3dDevice->CreateInputLayout( layout, sizeof(layout)/sizeof(layout[0]), 
                          PassDesc.pIAInputSignature, PassDesc.IAInputSignatureSize, &g_pLayout ) );
    DXUTCreateSphere( pd3dDevice, 1.0f, 16, 16, &g_pMesh );

    // Create vertex buffer
    SimpleVertex vertices[] =
    {
        { D3DXVECTOR3(-1.0f, 1.0f, -1.0f), D3DXVECTOR2(0.0f, 0.0f) },
        { D3DXVECTOR3(1.0f, 1.0f, -1.0f), D3DXVECTOR2(1.0f, 0.0f) },
        { D3DXVECTOR3(1.0f, 1.0f, 1.0f), D3DXVECTOR2(1.0f, 1.0f) },
        { D3DXVECTOR3(-1.0f, 1.0f, 1.0f), D3DXVECTOR2(0.0f, 1.0f) },

        { D3DXVECTOR3(-1.0f, -1.0f, -1.0f), D3DXVECTOR2(0.0f, 0.0f) },
        { D3DXVECTOR3(1.0f, -1.0f, -1.0f), D3DXVECTOR2(1.0f, 0.0f) },
        { D3DXVECTOR3(1.0f, -1.0f, 1.0f), D3DXVECTOR2(1.0f, 1.0f) },
        { D3DXVECTOR3(-1.0f, -1.0f, 1.0f), D3DXVECTOR2(0.0f, 1.0f) },

        { D3DXVECTOR3(-1.0f, -1.0f, 1.0f), D3DXVECTOR2(0.0f, 0.0f) },
        { D3DXVECTOR3(-1.0f, -1.0f, -1.0f), D3DXVECTOR2(1.0f, 0.0f) },
        { D3DXVECTOR3(-1.0f, 1.0f, -1.0f), D3DXVECTOR2(1.0f, 1.0f) },
        { D3DXVECTOR3(-1.0f, 1.0f, 1.0f), D3DXVECTOR2(0.0f, 1.0f) },

        { D3DXVECTOR3(1.0f, -1.0f, 1.0f), D3DXVECTOR2(0.0f, 0.0f) },
        { D3DXVECTOR3(1.0f, -1.0f, -1.0f), D3DXVECTOR2(1.0f, 0.0f) },
        { D3DXVECTOR3(1.0f, 1.0f, -1.0f), D3DXVECTOR2(1.0f, 1.0f) },
        { D3DXVECTOR3(1.0f, 1.0f, 1.0f), D3DXVECTOR2(0.0f, 1.0f) },

        { D3DXVECTOR3(-1.0f, -1.0f, -1.0f), D3DXVECTOR2(0.0f, 0.0f) },
        { D3DXVECTOR3(1.0f, -1.0f, -1.0f), D3DXVECTOR2(1.0f, 0.0f) },
        { D3DXVECTOR3(1.0f, 1.0f, -1.0f), D3DXVECTOR2(1.0f, 1.0f) },
        { D3DXVECTOR3(-1.0f, 1.0f, -1.0f), D3DXVECTOR2(0.0f, 1.0f) },

        { D3DXVECTOR3(-1.0f, -1.0f, 1.0f), D3DXVECTOR2(0.0f, 0.0f) },
        { D3DXVECTOR3(1.0f, -1.0f, 1.0f), D3DXVECTOR2(1.0f, 0.0f) },
        { D3DXVECTOR3(1.0f, 1.0f, 1.0f), D3DXVECTOR2(1.0f, 1.0f) },
        { D3DXVECTOR3(-1.0f, 1.0f, 1.0f), D3DXVECTOR2(0.0f, 1.0f) },
    };

    g_bd.Usage = D3D10_USAGE_DEFAULT;
    g_bd.ByteWidth = sizeof(SimpleVertex) * 24;
    g_bd.BindFlags = D3D10_BIND_VERTEX_BUFFER;
    g_bd.CPUAccessFlags = 0;
    g_bd.MiscFlags = 0;

    D3D10_SUBRESOURCE_DATA InitData;
    InitData.pSysMem = vertices;
    /*V_RETURN(*/pd3dDevice->CreateBuffer(&g_bd, &InitData, &g_pVertexBuffer)/*)*/;

    // Create index buffer
    // Create vertex buffer
    DWORD indices[] =
    {
        0, 1, 2,
        3, 0, 4,

        5, 1, 2//6, 7,
        //7, 4, 6,

        //11, 9, 8,
        //10, 9, 11,

        //14, 12, 13,
        //15, 12, 14,

        //19, 17, 16,
        //18, 17, 19,

        //22, 20, 21,
        //23, 20, 22
    };

    g_bd.Usage = D3D10_USAGE_DEFAULT;
    g_bd.ByteWidth = sizeof(DWORD) * 9;
    g_bd.BindFlags = D3D10_BIND_INDEX_BUFFER;
    g_bd.CPUAccessFlags = 0;
    g_bd.MiscFlags = 0;
    InitData.pSysMem = indices;
    /*V_RETURN(*/pd3dDevice->CreateBuffer(&g_bd, &InitData, &g_pIndexBuffer)/*)*/;

    D3DXVECTOR3 vEye(0, 4, -4);
    D3DXVECTOR3 vLook(0, 0, 0);

    D3DXVECTOR3 Up(0.0f, 1.0f, 0.0f);
    D3DXMatrixLookAtLH(&g_View, &vEye, &vLook, &Up);

    // Update Variables that never change
    g_pViewVariable->SetMatrix((float*)&g_View);

    g_Camera.SetViewParams(&vEye, &vLook);

    return S_OK;
}

//--------------------------------------------------------------------------------------
// Create any D3D10 resources that depend on the back buffer
//--------------------------------------------------------------------------------------
HRESULT CALLBACK OnD3D10ResizedSwapChain( ID3D10Device* pd3dDevice, IDXGISwapChain *pSwapChain, const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc, void* pUserContext )
{
    float fAspectRatio = pBackBufferSurfaceDesc->Width / (FLOAT)pBackBufferSurfaceDesc->Height;
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
    SAFE_RELEASE(g_pEffect);
    SAFE_RELEASE(g_pLayout);
    SAFE_RELEASE(g_pIndexBuffer);
    SAFE_RELEASE(g_pVertexBuffer);
    SAFE_RELEASE(g_pMesh);
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

    DXUTInit( true, true, NULL ); // Parse the command line, show msgboxes on error, no extra command line params
    DXUTSetCursorSettings( true, true ); // Show the cursor and clip it when in full screen
    DXUTCreateWindow( L"Teapot ArcBall" );
    DXUTCreateDevice( true, 640, 480 );  
    DXUTMainLoop(); // Enter into the DXUT render loop

    return DXUTGetExitCode();
}

