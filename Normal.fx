//--------------------------------------------------------------------------------------
// File: Normal.fx    頂点座標＋法線ベクトル＋モデルの色 ＆ 光源    前田 稔
//--------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------
// Constant Buffer Variables
//--------------------------------------------------------------------------------------
matrix World;
matrix View;
matrix Projection;
float4 Color;

//--------------------------------------------------------------------------------------
// Global variables
//--------------------------------------------------------------------------------------
cbuffer cb0
{
//    float3 g_vLightDir = float3(-0.707,0.707,0);    // 光源の座標
    float3 g_vLightDir = float3(-0.707,0.707,-0.5); // 光源の座標
};

//--------------------------------------------------------------------------------------
// Vertex shader output structure
//--------------------------------------------------------------------------------------
struct VS_INPUT
{
    float4 Pos : POSITION;
    float3 Norm : NORMAL;
};

struct PS_INPUT
{
    float4 Pos : SV_POSITION;
    float3 Norm : NORMAL;
};

//--------------------------------------------------------------------------------------
// Vertex Shader
//--------------------------------------------------------------------------------------
PS_INPUT VS( VS_INPUT input )
{
    PS_INPUT output = (PS_INPUT)0;
    output.Pos = mul( input.Pos, World );
    output.Pos = mul( output.Pos, View );
    output.Pos = mul( output.Pos, Projection );
    output.Norm = mul( input.Norm, World );
    
    return output;
}

//--------------------------------------------------------------------------------------
float4 PS( PS_INPUT input ) : SV_Target
{ 
    //return Color * saturate( dot( normalize(input.Norm), g_vLightDir ) );
    return Color * ( saturate(dot(normalize(input.Norm),g_vLightDir)) * 0.7f + 0.3f);
}

//--------------------------------------------------------------------------------------
technique10 Render
{
    pass P0
    {
        SetVertexShader( CompileShader( vs_4_0, VS() ) );
        SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_4_0, PS() ) );
    }
}

