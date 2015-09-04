/*! \file LJ_Argon_MD.cpp
    \brief 分子動力学シミュレーションを描画する

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "DXUT.h"
#include "SDKmisc.h"
#include "DXUTcamera.h"
#include "DXUTgui.h"
#include "DXUTsettingsDlg.h"
#include "DXUTShapes.h"
#include "moleculardynamics/Ar_moleculardynamics.h"
#include "utility/utility.h"
#include <array>									// for std::array
#include <memory>									// for std::unique_ptr
#include <vector>									// for std::vector
#include <boost/assert.hpp>                         // for BOOST_ASSERT
#include <boost/cast.hpp>							// for boost::numeric_cast
#include <boost/range/algorithm/max_element.hpp>	// for boost::max_element

//! A global variable (constant).
/*!
	箱の比率
*/
static auto const BOXRATIO = 1.1f;

//! A global variable (constant).
/*!
	色の比率
*/
static auto const COLORRATIO = 0.025f;

//! A global variable (constant).
/*!
    インデックスバッファの個数
*/
static auto const NUMINDEXBUFFER = 16U;

//! A global variable (constant).
/*!
    頂点バッファの個数
*/
static auto const NUMVERTEXBUFFER = 8U;

//! A global variable (constant).
/*!
    頂点バッファの個数
*/
static auto const TINY = 1.0E-30f;

//! A global variable (constant).
/*!
    画面サイズ（高さ）
*/
static auto const WINDOWHEIGHT = 960;

//! A global variable (constant).
/*!
    画面サイズ（幅）
*/
static auto const WINDOWWIDTH = 1280;

//! A global variable.
/*!
    分子動力学シミュレーションのオブジェクト
*/
moleculardynamics::Ar_moleculardynamics armd;

//! A global variable.
/*!
    バッファー リソース
*/
D3D10_BUFFER_DESC bd;

//! A global variable.
/*!
    Font for drawing text
*/
std::unique_ptr<ID3DX10Font, utility::Safe_Release<ID3DX10Font>> font;

//! A global variable.
/*!
    ブレンディング・ステート
*/
std::unique_ptr<ID3D10BlendState, utility::Safe_Release<ID3D10BlendState>> pBlendStateNoBlend;

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
    メッシュへのスマートポインタが格納された可変長配列
*/
std::vector<std::unique_ptr<ID3DX10Mesh, utility::Safe_Release<ID3DX10Mesh>>> pmeshvec;

//! A global variable.
/*!
	頂点バッファ
*/
std::unique_ptr<ID3D10Buffer, utility::Safe_Release<ID3D10Buffer>> pVertexBuffer;

//! A global variable.
/*!
    Sprite for batching text drawing
*/
std::unique_ptr<ID3DX10Sprite, utility::Safe_Release<ID3DX10Sprite>> sprite;

//! A global variable.
/*!
    テキスト表示用
*/
std::unique_ptr<CDXUTTextHelper, utility::Safe_Delete<CDXUTTextHelper>> txthelper;

//! A global variable.
/*!
    A model viewing camera
*/
CModelViewerCamera g_Camera;

//! A global variable.
/*!
    manager for shared resources of dialogs
*/
CDXUTDialogResourceManager g_DialogResourceManager;

//! A global variable.
/*!
    Device settings dialog
*/
CD3DSettingsDlg g_D3DSettingsDlg;

//! A global variable.
/*!
    manages the 3D UI
*/
CDXUTDialog g_HUD;

//! A global variable.
/*!
*/
ID3D10EffectVectorVariable* g_pColorVariable = nullptr;

//! A global variable.
/*!
*/
ID3D10EffectMatrixVariable* g_pProjectionVariable = nullptr;

//! A global variable.
/*!
*/
ID3D10EffectTechnique*      g_pRender = nullptr;

//! A global variable.
/*!
*/
ID3D10EffectMatrixVariable* g_pViewVariable = nullptr;

//! A global variable.
/*!
*/
ID3D10EffectMatrixVariable* g_pWorldVariable = nullptr;

//! A global variable.
/*!
    ビュー行列
*/
D3DXMATRIX                  g_View;

D3DXVECTOR4 g_Colors[2] = 
{
    D3DXVECTOR4( 1.0f, 1.0f, 1.0f, 1.0f ),
    D3DXVECTOR4( 1.0f, 0.0f, 1.0f, 1.0f ),
};

//--------------------------------------------------------------------------------------
// UI control IDs
//--------------------------------------------------------------------------------------
#define IDC_TOGGLEFULLSCREEN    1
#define IDC_CHANGEDEVICE        2
#define IDC_TOGGLEROTATION      3
#define IDC_RECALC              4
#define IDC_SLIDER				5

void CALLBACK OnGUIEvent(UINT nEvent, int nControlID, CDXUTControl* pControl, void* pUserContext);

//! A function.
/*!
    画面の左上に情報を表示する
*/
void RenderText(ID3D10Device* pd3dDevice, double fTime);

//! A function.
/*!
    UIを配置する
*/
void SetUI();

//--------------------------------------------------------------------------------------
// Initialize the app 
//--------------------------------------------------------------------------------------
void InitApp()
{
    g_D3DSettingsDlg.Init(&g_DialogResourceManager);
    g_HUD.Init(&g_DialogResourceManager);

    g_HUD.SetCallback(OnGUIEvent);

    SetUI();
}

//! A function.
/*!
    UIを配置する
*/
void SetUI();

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
    if (g_D3DSettingsDlg.IsActive())
    {
        std::array<float, 4> const ClearColor = { 0.176f, 0.196f, 0.667f, 1.0f };
        auto pRTV = DXUTGetD3D10RenderTargetView();
        pd3dDevice->ClearRenderTargetView(pRTV, ClearColor.data());

        g_D3DSettingsDlg.OnRender(fElapsedTime);
        return;
    }
    else
    {
        armd.Calc_Forces();
        armd.Move_Atoms();

        // Clear render target and the depth stencil 
        std::array<float, 4> const ClearColor = { 0.176f, 0.196f, 0.667f, 1.0f };
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

        auto const size = pmeshvec.size();
        static auto const origin = boost::numeric_cast<float>(*boost::max_element(armd.X())) / 2.0f;

        for (auto i = 0U; i < size; i++) {
            auto color = g_Colors[1];
			auto const rcolor = COLORRATIO * armd.getForce(i);
            color.x = rcolor > 1.0f ? 1.0f : rcolor;
            g_pColorVariable->SetFloatVector(reinterpret_cast<float *>(&color));

            D3DXMATRIX World;
            D3DXMatrixTranslation(
                &World,
                boost::numeric_cast<float>(armd.X()[i]) - origin,
                boost::numeric_cast<float>(armd.Y()[i]) - origin,
                boost::numeric_cast<float>(armd.Z()[i]) - origin);
            D3DXMatrixMultiply(&World, &(*g_Camera.GetWorldMatrix()), &World);

            // Update variables
            auto mWorld = World * (*g_Camera.GetWorldMatrix());
            g_pWorldVariable->SetMatrix(reinterpret_cast<float *>(&mWorld));

            UINT NumSubsets;
            pmeshvec[i]->GetAttributeTable(nullptr, &NumSubsets);

            for (auto p = 0U; p < techDesc.Passes; p++)
            {
                g_pRender->GetPassByIndex(p)->Apply(0);
                for (auto s = 0U; s < NumSubsets; s++)
                {
                    pmeshvec[i]->DrawSubset(s);
                }
            }
        }

        DXUT_BeginPerfEvent(DXUT_PERFEVENTCOLOR, L"HUD / Stats");
        g_HUD.OnRender(fElapsedTime);
        RenderText(pd3dDevice, fTime);
        DXUT_EndPerfEvent();
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

    V_RETURN(g_DialogResourceManager.OnD3D10CreateDevice(pd3dDevice));
    V_RETURN(g_D3DSettingsDlg.OnD3D10CreateDevice(pd3dDevice));

    ID3DX10Font * fonttemp;
    V_RETURN(D3DX10CreateFont(pd3dDevice, 15, 0, FW_BOLD, 1, FALSE, DEFAULT_CHARSET,
        OUT_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH | FF_DONTCARE,
        L"Arial", &fonttemp));
    font.reset(fonttemp);

    ID3DX10Sprite * spritetmp;
    V_RETURN(D3DX10CreateSprite(pd3dDevice, 512, &spritetmp));
    sprite.reset(spritetmp);
    txthelper.reset(new CDXUTTextHelper(nullptr, nullptr, font.get(), sprite.get(), 15));

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
    
    D3D10_BLEND_DESC BlendState;
    ZeroMemory(&BlendState, sizeof(D3D10_BLEND_DESC));
    BlendState.BlendEnable[0] = FALSE;
    BlendState.RenderTargetWriteMask[0] = D3D10_COLOR_WRITE_ENABLE_ALL;

    ID3D10BlendState * pBlendStateNoBlendtmp = nullptr;
    pd3dDevice->CreateBlendState(&BlendState, &pBlendStateNoBlendtmp);
    pBlendStateNoBlend.reset(pBlendStateNoBlendtmp);

    pmeshvec.resize(armd.X().size());
    for (auto & pmesh : pmeshvec) {
        ID3DX10Mesh * pmeshtmp = nullptr;
        DXUTCreateSphere(pd3dDevice, 0.2f, 16, 16, &pmeshtmp);
        pmesh.reset(pmeshtmp);
    }

    auto boundarylen = boost::numeric_cast<float>(*boost::max_element(armd.X()));
    BOOST_ASSERT(std::fabs(boost::numeric_cast<float>(*boost::max_element(armd.Y())) - boundarylen) < TINY);
    BOOST_ASSERT(std::fabs(boost::numeric_cast<float>(*boost::max_element(armd.Z())) - boundarylen) < TINY);

	auto const pos = boundarylen * BOXRATIO;

    // Create vertex buffer
    std::array<SimpleVertex, NUMVERTEXBUFFER> const vertices =
    {
        D3DXVECTOR3(-pos, pos, -pos),
        D3DXVECTOR3(pos, pos, -pos),
        D3DXVECTOR3(pos, pos, pos), 
        D3DXVECTOR3(-pos, pos, pos),

        D3DXVECTOR3(-pos, -pos, -pos), 
        D3DXVECTOR3(pos, -pos, -pos),
        D3DXVECTOR3(pos, -pos, pos),
        D3DXVECTOR3(-pos, -pos, pos),
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
	
    D3DXVECTOR3 vEye(0.0f, 10.0f, 10.0f);
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
    HRESULT hr;

    V_RETURN(g_DialogResourceManager.OnD3D10ResizedSwapChain(pd3dDevice, pBackBufferSurfaceDesc));
    V_RETURN(g_D3DSettingsDlg.OnD3D10ResizedSwapChain(pd3dDevice, pBackBufferSurfaceDesc));

    g_HUD.SetLocation(pBackBufferSurfaceDesc->Width - 170, 0);
    g_HUD.SetSize(170, 170);
    
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
// Handles the GUI events
//--------------------------------------------------------------------------------------
void CALLBACK OnGUIEvent(UINT nEvent, int nControlID, CDXUTControl* pControl, void* pUserContext)
{
    switch (nControlID)
    {
    case IDC_TOGGLEFULLSCREEN:
        DXUTToggleFullScreen();
        break;

    case IDC_CHANGEDEVICE:
        g_D3DSettingsDlg.SetActive(!g_D3DSettingsDlg.IsActive());
        break;

    case IDC_TOGGLEROTATION:
        //ROT_FLAG = !ROT_FLAG;
        break;

    case IDC_RECALC:
        armd.recalc();
        break;

   /* case IDC_READDATA:
        StopDraw();
        ReadData();
        SetUI();
        scene->Pgd = pgd;
        scene->Thread_end = false;
        scene->Redraw = true;
        first = true;
        ::SetWindowText(DXUTGetHWND(), CreateWindowTitle().c_str());
        break;

    case IDC_COMBOBOX:
    {
        auto const pItem = (static_cast<CDXUTComboBox *>(pControl))->GetSelectedItem();
        if (pItem)
        {
            drawdata = reinterpret_cast<std::uint32_t>(pItem->pData);
            RedrawFlagTrue();
        }
        break;
    }

    case IDC_RADIOA:
        reim = TDXScene::Re_Im_type::REAL;
        RedrawFlagTrue();
        break;

    case IDC_RADIOB:
        reim = TDXScene::Re_Im_type::IMAGINARY;
        RedrawFlagTrue();
        break;

    case IDC_SLIDER:
        scene->Vertexsize(static_cast<std::vector<TDXScene::SimpleVertex2>::size_type>((reinterpret_cast<CDXUTSlider*>(pControl))->GetValue()));
        RedrawFlagTrue();
        break;*/

    default:
        break;
    }
}

//--------------------------------------------------------------------------------------
// Release D3D10 resources created in OnD3D10ResizedSwapChain 
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D10ReleasingSwapChain( void* pUserContext )
{
    g_DialogResourceManager.OnD3D10ReleasingSwapChain();
}

//--------------------------------------------------------------------------------------
// Release D3D10 resources created in OnD3D10CreateDevice 
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D10DestroyDevice( void* pUserContext )
{
    font.reset();
    pBlendStateNoBlend.reset();
    pEffect.reset();
	pInputLayout.reset();
    pIndexBuffer.reset();
    pVertexBuffer.reset();
    sprite.reset();
    txthelper.reset();

    for (auto & pmesh : pmeshvec) {
        pmesh.reset();
    }

    g_DialogResourceManager.OnD3D10DestroyDevice();
    g_D3DSettingsDlg.OnD3D10DestroyDevice();
    DXUTGetGlobalResourceCache().OnDestroyDevice();
}

//--------------------------------------------------------------------------------------
// Handle messages to the application
//--------------------------------------------------------------------------------------
LRESULT CALLBACK MsgProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam, 
                          bool* pbNoFurtherProcessing, void* pUserContext )
{
    // Pass messages to dialog resource manager calls so GUI state is updated correctly
    *pbNoFurtherProcessing = g_DialogResourceManager.MsgProc(hWnd, uMsg, wParam, lParam);
    if (*pbNoFurtherProcessing)
        return 0;

    // Pass messages to settings dialog if its active
    if (g_D3DSettingsDlg.IsActive())
    {
        g_D3DSettingsDlg.MsgProc(hWnd, uMsg, wParam, lParam);
        return 0;
    }

    // Give the dialogs a chance to handle the message first
    *pbNoFurtherProcessing = g_HUD.MsgProc(hWnd, uMsg, wParam, lParam);
    if (*pbNoFurtherProcessing)
        return 0;

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

void SetUI()
{
    g_HUD.RemoveAllControls();

    std::int32_t iY = 10;

    g_HUD.AddButton(IDC_TOGGLEFULLSCREEN, L"Toggle full screen", 35, iY, 125, 22);
    g_HUD.AddButton(IDC_CHANGEDEVICE, L"Change device (F2)", 35, iY += 24, 125, 22, VK_F2);
    g_HUD.AddButton(IDC_TOGGLEROTATION, L"Toggle Rotaion Animation", 35, iY += 24, 125, 22);

    g_HUD.AddButton(IDC_RECALC, L"再計算", 35, iY += 34, 125, 22);

    // Combobox
    //CDXUTComboBox* pCombo;
    //g_HUD.AddComboBox(IDC_COMBOBOX, 35, iY += 34, 125, 22, L'O', false, &pCombo);
    //if (pCombo)
    //{
    //    pCombo->SetDropHeight(100);
    //    pCombo->RemoveAllItems();
    //    BOOST_ASSERT(pgd->N > static_cast<std::int32_t>(pgd->L));
    //    auto orbital(std::to_wstring(pgd->N));
    //    switch (pgd->L) {
    //    case 0:
    //    {
    //        orbital += L's';
    //        pCombo->AddItem(orbital.c_str(), reinterpret_cast<LPVOID>(0x11111111));
    //    }
    //    break;

    //    case 1:
    //    {
    //        pCombo->AddItem((orbital + L"px").c_str(), reinterpret_cast<LPVOID>(0x11111111));
    //        pCombo->AddItem((orbital + L"py").c_str(), reinterpret_cast<LPVOID>(0x12121212));
    //        pCombo->AddItem((orbital + L"pz").c_str(), reinterpret_cast<LPVOID>(0x13131313));
    //    }
    //    break;

    //    case 2:
    //    {
    //        pCombo->AddItem((orbital + L"dxy").c_str(), reinterpret_cast<LPVOID>(0x11111111));
    //        pCombo->AddItem((orbital + L"dyz").c_str(), reinterpret_cast<LPVOID>(0x12121212));
    //        pCombo->AddItem((orbital + L"dzx").c_str(), reinterpret_cast<LPVOID>(0x13131313));
    //        pCombo->AddItem((orbital + L"dx^2-y^2").c_str(), reinterpret_cast<LPVOID>(0x14141414));
    //        pCombo->AddItem((orbital + L"dz^2").c_str(), reinterpret_cast<LPVOID>(0x15151515));
    //    }
    //    break;

    //    case 3:
    //    {
    //        pCombo->AddItem((orbital + L"fxz^2").c_str(), reinterpret_cast<LPVOID>(0x11111111));
    //        pCombo->AddItem((orbital + L"fyz^2").c_str(), reinterpret_cast<LPVOID>(0x12121212));
    //        pCombo->AddItem((orbital + L"fz(x^2-y^2)").c_str(), reinterpret_cast<LPVOID>(0x13131313));
    //        pCombo->AddItem((orbital + L"fxyz").c_str(), reinterpret_cast<LPVOID>(0x14141414));
    //        pCombo->AddItem((orbital + L"fx(x^2-3y^2)").c_str(), reinterpret_cast<LPVOID>(0x15151515));
    //        pCombo->AddItem((orbital + L"fy(3x^2-y^2)").c_str(), reinterpret_cast<LPVOID>(0x16161616));
    //        pCombo->AddItem((orbital + L"fz^2").c_str(), reinterpret_cast<LPVOID>(0x17171717));
    //    }
    //    break;

    //    default:
    //        throw std::runtime_error("g以上の軌道には対応していません");
    //        break;
    //    }

    //    if (pgd->Rho_wf_type_ == getdata::GetData::Rho_Wf_type::WF) {
    //        // Radio buttons
    //        g_HUD.AddRadioButton(IDC_RADIOA, 1, L"実部", 35, iY += 34, 125, 22, true, L'1');
    //        g_HUD.AddRadioButton(IDC_RADIOB, 1, L"虚部", 35, iY += 28, 125, 22, false, L'2');
    //    }
    //}

    //// 角度の調整
    //g_HUD.AddStatic(IDC_OUTPUT, L"頂点数", 20, iY += 34, 125, 22);
    //g_HUD.GetStatic(IDC_OUTPUT)->SetTextColor(D3DCOLOR_ARGB(255, 255, 255, 255));
    //g_HUD.AddSlider(IDC_SLIDER, 35, iY += 24, 125, 22, 0, 1000000, TDXScene::VERTEXSIZE_FIRST);
}

//--------------------------------------------------------------------------------------
// Render the help and statistics text
//--------------------------------------------------------------------------------------
void RenderText(ID3D10Device* pd3dDevice, double fTime)
{
    txthelper->Begin();
    txthelper->SetInsertionPos(2, 0);
    txthelper->SetForegroundColor(D3DXCOLOR(1.000f, 0.945f, 0.059f, 1.000f));
    txthelper->DrawTextLine(DXUTGetFrameStats(DXUTIsVsyncEnabled()));
    txthelper->DrawTextLine(DXUTGetDeviceStats());
    //txthelper->DrawTextLine((boost::wformat(L"CPUスレッド数: %d") % cputhread).str().c_str());
    //txthelper->DrawTextLine((boost::wformat(L"頂点数 = %d") % scene->Vertexsize()).str().c_str());
    //txthelper->DrawTextLine(str.c_str());
    txthelper->End();
    pd3dDevice->IASetInputLayout(pInputLayout.get());

    auto const blendFactor = 0.0f;
    auto sampleMask = 0xffffffff;

    pd3dDevice->OMSetBlendState(pBlendStateNoBlend.get(), &blendFactor, sampleMask);
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
    
    InitApp();

    // ウィンドウを生成
    auto const dispx = ::GetSystemMetrics(SM_CXSCREEN);
    auto const dispy = ::GetSystemMetrics(SM_CYSCREEN);
    auto const xpos = (dispx - WINDOWWIDTH) >> 1;
    auto const ypos = (dispy - WINDOWHEIGHT) >> 1;
    DXUTCreateWindow(L"アルゴンの古典分子動力学シミュレーション", nullptr, nullptr, nullptr, xpos, ypos);
    DXUTCreateDevice(true, WINDOWWIDTH, WINDOWHEIGHT);

    // 垂直同期をオフにする
    auto ds = DXUTGetDeviceSettings();
    ds.d3d10.SyncInterval = 0;
    DXUTCreateDeviceFromSettings(&ds);

    DXUTMainLoop(); // Enter into the DXUT render loop

    return DXUTGetExitCode();
}

