#define PIXEL_SHADER

#include "fp_consts.h"
#include "fp_structs.h"
#include "fp_funcs.h"




float4 ApplyTexProj(float3 V)
{
  return g_r2o_tex_proj_mat_0 * V.x +
         g_r2o_tex_proj_mat_1 * V.y +
         g_r2o_tex_proj_mat_2 * V.z +
         g_r2o_tex_proj_mat_3;
}


float GetWFromDepthBuffer(uniform texobj2D depth_buffer, float4 VP)
{
  float4 depth_texel = tex2Dproj(depth_buffer, VP.xyw);
  float  dot_prod = dot(depth_texel, g_r2o_zbuffer_constants);
  return 1.0/dot_prod;
}


half4 r2o_main(FpR2O input, float4 pos, float2 int_map_derivs,
               uniform texobj2D amb_map_coarse,
               uniform texobj2D amb_map_fine,
               uniform texobj2D depth_buffer,
               uniform texobj2D color_buffer,
               uniform texobj2D refl_texture,
               uniform texobj1D fresnel_texture)
{
  // get ambient map derivs
  float2 coarse_map_derivs = tex2D(amb_map_coarse, input.uv.xy).xy;
  float2 fine_map_derivs   = tex2D(amb_map_fine,   input.uv.zw).xy;

  // sum derivs
  float2 derivs = input.derivs.xy + g_r2o_normal_map_weights.x * coarse_map_derivs
                                  + g_r2o_normal_map_weights.y * fine_map_derivs
                                  + g_r2o_normal_map_weights.z * int_map_derivs;

  // V0 : point on water surface at pixel location
  float3 V0 = -input.viewvec.xyz;
  float  w0 = input.derivs.w;   //  = 1.0/pos.w

  //  normal
  float3 N = normalize(float3(-derivs.x, 1.0, -derivs.y));

  // transmitted vector (not even vaguely related to the physically correct one!)
  float3 T = float3(0,1,0) - N;

  // VT: transmission map lookup vector
  float  scale = w0 * g_r2o_misc_consts.x + g_r2o_misc_consts.y;
  float3 VT = V0 + scale * T;

  // derive homogeneous texcoords for this point
  float4 VTP = ApplyTexProj(VT);

  // look up in depth buffer
  float wT = GetWFromDepthBuffer(depth_buffer, VTP);

  // if looked-up depthbuffer value at perturbed position is closer than water surface
  if (wT < w0)
  {
    // just use current pixel location instead
    VTP = g_r2o_fb_scale * pos + g_r2o_fb_offset;
  }

  // look up transmitted colour in the reduced framebuffer texture
  float3 trans_color = tex2Dproj(color_buffer, VTP.xyw).xyz;
  wT = GetWFromDepthBuffer(depth_buffer, VTP);

  // view vector
  float3 V = normalize(input.viewvec.xyz);

  // reflected vector
  float  d = dot(N,V);
  float3 R = 2.0*d*N - V;

  // look up reflected colour
  float3 VR = float3(R.x, -R.y, R.z);
  float4 VRP = ApplyTexProj(VR);
  float3 refl_color = tex2Dproj(refl_texture, VRP.xyw).xyz;

  // prepare all the blend factors
  float thickness = abs(wT-w0);
  float transp    = exp(thickness * -g_r2o_water_color.w);
  float fresnel   = tex1D(fresnel_texture, 0.5*d+0.5).x;
  float fog       = input.viewvec.w;
  float bloom     = ComputeBloomContribution(refl_color) * (1.0 - fog);

  // a fudge-factor to soften the hard edges of reflections
  float fudgefactor = g_r2o_misc_consts.z * thickness * d;
  fudgefactor = saturate(abs(fudgefactor));
  fresnel *= fudgefactor;

  // dirty the transmitted colour using the computed blend factor
  float3 water_color = lerp(g_r2o_water_color.xyz, trans_color, transp);

  // apply Fresnel reflectance
  float3 final_color = lerp(water_color, refl_color, fresnel);

  // combine with bloom channel
  return half4(final_color, bloom);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


half4 r2o_opaque(FpR2O input, float4 pos, float2 int_map_derivs,
                 uniform texobj2D amb_map_coarse,
                 uniform texobj2D amb_map_fine,
                 uniform texobj2D refl_texture,
                 uniform texobj1D fresnel_texture)
{
  // get ambient map derivs
  float2 coarse_map_derivs = tex2D(amb_map_coarse, input.uv.xy).xy;
  float2 fine_map_derivs   = tex2D(amb_map_fine,   input.uv.zw).xy;

  // sum derivs
  float2 derivs = input.derivs.xy + g_r2o_normal_map_weights.x * coarse_map_derivs
                                  + g_r2o_normal_map_weights.y * fine_map_derivs
                                  + g_r2o_normal_map_weights.z * int_map_derivs;

  // V0 : point on water surface at pixel location
  float3 V0 = -input.viewvec.xyz;
  float  w0 = input.derivs.w;   //  = 1.0/pos.w

  //  normal
  float3 N = normalize(float3(-derivs.x, 1.0, -derivs.y));

  // view vector
  float3 V = normalize(input.viewvec.xyz);

  // reflected vector
  float  d = dot(N,V);
  float3 R = 2.0*d*N - V;

  // look up reflected colour
  float3 VR = float3(R.x, -R.y, R.z);
  float4 VRP = ApplyTexProj(VR);
  float3 refl_color = tex2Dproj(refl_texture, VRP.xyw).xyz;

  // prepare all the blend factors
  float fresnel = tex1D(fresnel_texture, 0.5*d+0.5).x;
  float fog = input.viewvec.w;
  float bloom = ComputeBloomContribution(refl_color) * (1.0 - fog);

  // apply Fresnel reflectance
  float3 final_color = lerp(g_r2o_water_color.xyz, refl_color, fresnel);

  // combine with bloom channel
  return half4(final_color, bloom);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float2 r2o_IntMapDerivs(float2  uv,
                        uniform texobj2D index_map,
                        uniform texobj3D norm_map)
{
  float tile_index   = tex2D(index_map, uv).x;
  float tile_indices = tile_index;
  tile_indices       = round(255.0*tile_indices) - 0.5;
  float2 uv_frac     = float2(16.0, 16.0) * frac(uv) + float2(4.0, 1.0);
  float3 coords      = float3(uv_frac, tile_indices.x);
  float2 derivs      = tex3D(norm_map, coords).xy;
  return derivs;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_0(FpR2O input, in float4 pos: WPOS,
                         uniform texobj2D amb_map_coarse  : TEXUNIT0,
                         uniform texobj2D amb_map_fine    : TEXUNIT1,
                         uniform texobj2D depth_buffer    : TEXUNIT2,
                         uniform texobj2D color_buffer    : TEXUNIT3,
                         uniform texobj2D refl_texture    : TEXUNIT4,
                         uniform texobj1D fresnel_texture : TEXUNIT5) : COLOR
{
  float2 int_map_derivs = float2(0.0, 0.0);

  return r2o_main(input, pos, int_map_derivs,
                  amb_map_coarse, amb_map_fine, depth_buffer, color_buffer, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_1(FpR2O input, in float4 pos: WPOS,
                         uniform texobj2D amb_map_coarse  : TEXUNIT0,
                         uniform texobj2D amb_map_fine    : TEXUNIT1,
                         uniform texobj2D depth_buffer    : TEXUNIT2,
                         uniform texobj2D color_buffer    : TEXUNIT3,
                         uniform texobj2D refl_texture    : TEXUNIT4,
                         uniform texobj1D fresnel_texture : TEXUNIT5,
                         uniform texobj2D int_index_map2  : TEXUNIT8,
                         uniform texobj3D int_norm_map2   : TEXUNIT12) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv2, int_index_map2, int_norm_map2);
  return r2o_main(input, pos, int_map_derivs,
                  amb_map_coarse, amb_map_fine, depth_buffer, color_buffer, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_2(FpR2O input, in float4 pos: WPOS,
                         uniform texobj2D amb_map_coarse  : TEXUNIT0,
                         uniform texobj2D amb_map_fine    : TEXUNIT1,
                         uniform texobj2D depth_buffer    : TEXUNIT2,
                         uniform texobj2D color_buffer    : TEXUNIT3,
                         uniform texobj2D refl_texture    : TEXUNIT4,
                         uniform texobj1D fresnel_texture : TEXUNIT5,
                         uniform texobj2D int_index_map1  : TEXUNIT7,
                         uniform texobj3D int_norm_map1   : TEXUNIT11) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.zw, int_index_map1, int_norm_map1);
  return r2o_main(input, pos, int_map_derivs,
                  amb_map_coarse, amb_map_fine, depth_buffer, color_buffer, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_3(FpR2O input, in float4 pos: WPOS,
                         uniform texobj2D amb_map_coarse  : TEXUNIT0,
                         uniform texobj2D amb_map_fine    : TEXUNIT1,
                         uniform texobj2D depth_buffer    : TEXUNIT2,
                         uniform texobj2D color_buffer    : TEXUNIT3,
                         uniform texobj2D refl_texture    : TEXUNIT4,
                         uniform texobj1D fresnel_texture : TEXUNIT5,
                         uniform texobj2D int_index_map1  : TEXUNIT7,
                         uniform texobj2D int_index_map2  : TEXUNIT8,
                         uniform texobj3D int_norm_map1   : TEXUNIT11,
                         uniform texobj3D int_norm_map2   : TEXUNIT12) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.zw, int_index_map1, int_norm_map1) +
                          r2o_IntMapDerivs(input.int_uv2,     int_index_map2, int_norm_map2);
  return r2o_main(input, pos, int_map_derivs,
                  amb_map_coarse, amb_map_fine, depth_buffer, color_buffer, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_4(FpR2O input, in float4 pos: WPOS,
                         uniform texobj2D amb_map_coarse  : TEXUNIT0,
                         uniform texobj2D amb_map_fine    : TEXUNIT1,
                         uniform texobj2D depth_buffer    : TEXUNIT2,
                         uniform texobj2D color_buffer    : TEXUNIT3,
                         uniform texobj2D refl_texture    : TEXUNIT4,
                         uniform texobj1D fresnel_texture : TEXUNIT5,
                         uniform texobj2D int_index_map0  : TEXUNIT6,
                         uniform texobj3D int_norm_map0   : TEXUNIT10) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.xy, int_index_map0, int_norm_map0);
  return r2o_main(input, pos, int_map_derivs,
                  amb_map_coarse, amb_map_fine, depth_buffer, color_buffer, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_5(FpR2O input, in float4 pos: WPOS,
                         uniform texobj2D amb_map_coarse  : TEXUNIT0,
                         uniform texobj2D amb_map_fine    : TEXUNIT1,
                         uniform texobj2D depth_buffer    : TEXUNIT2,
                         uniform texobj2D color_buffer    : TEXUNIT3,
                         uniform texobj2D refl_texture    : TEXUNIT4,
                         uniform texobj1D fresnel_texture : TEXUNIT5,
                         uniform texobj2D int_index_map0  : TEXUNIT6,
                         uniform texobj2D int_index_map2  : TEXUNIT8,
                         uniform texobj3D int_norm_map0   : TEXUNIT10,
                         uniform texobj3D int_norm_map2   : TEXUNIT12) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.xy, int_index_map0, int_norm_map0) +
                          r2o_IntMapDerivs(input.int_uv2,     int_index_map2, int_norm_map2);
  return r2o_main(input, pos, int_map_derivs,
                  amb_map_coarse, amb_map_fine, depth_buffer, color_buffer, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_6(FpR2O input, in float4 pos: WPOS,
                         uniform texobj2D amb_map_coarse  : TEXUNIT0,
                         uniform texobj2D amb_map_fine    : TEXUNIT1,
                         uniform texobj2D depth_buffer    : TEXUNIT2,
                         uniform texobj2D color_buffer    : TEXUNIT3,
                         uniform texobj2D refl_texture    : TEXUNIT4,
                         uniform texobj1D fresnel_texture : TEXUNIT5,
                         uniform texobj2D int_index_map0  : TEXUNIT6,
                         uniform texobj2D int_index_map1  : TEXUNIT7,
                         uniform texobj3D int_norm_map0   : TEXUNIT10,
                         uniform texobj3D int_norm_map1   : TEXUNIT11) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.xy, int_index_map0, int_norm_map0) +
                          r2o_IntMapDerivs(input.int_uv01.zw, int_index_map1, int_norm_map1);
  return r2o_main(input, pos, int_map_derivs,
                  amb_map_coarse, amb_map_fine, depth_buffer, color_buffer, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_7(FpR2O input, in float4 pos: WPOS,
                         uniform texobj2D amb_map_coarse  : TEXUNIT0,
                         uniform texobj2D amb_map_fine    : TEXUNIT1,
                         uniform texobj2D depth_buffer    : TEXUNIT2,
                         uniform texobj2D color_buffer    : TEXUNIT3,
                         uniform texobj2D refl_texture    : TEXUNIT4,
                         uniform texobj1D fresnel_texture : TEXUNIT5,
                         uniform texobj2D int_index_map0  : TEXUNIT6,
                         uniform texobj2D int_index_map1  : TEXUNIT7,
                         uniform texobj2D int_index_map2  : TEXUNIT8,
                         uniform texobj3D int_norm_map0   : TEXUNIT10,
                         uniform texobj3D int_norm_map1   : TEXUNIT11,
                         uniform texobj3D int_norm_map2   : TEXUNIT12) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.xy, int_index_map0, int_norm_map0) +
                          r2o_IntMapDerivs(input.int_uv01.zw, int_index_map1, int_norm_map1) +
                          r2o_IntMapDerivs(input.int_uv2,     int_index_map2, int_norm_map2);
  return r2o_main(input, pos, int_map_derivs,
                  amb_map_coarse, amb_map_fine, depth_buffer, color_buffer, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_opaque_0(FpR2O input, in float4 pos: WPOS,
                                uniform texobj2D amb_map_coarse  : TEXUNIT0,
                                uniform texobj2D amb_map_fine    : TEXUNIT1,
                                uniform texobj2D refl_texture    : TEXUNIT4,
                                uniform texobj1D fresnel_texture : TEXUNIT5) : COLOR
{
  float2 int_map_derivs = float2(0.0, 0.0);

  return r2o_opaque(input, pos, int_map_derivs,
                    amb_map_coarse, amb_map_fine, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_opaque_1(FpR2O input, in float4 pos: WPOS,
                                uniform texobj2D amb_map_coarse  : TEXUNIT0,
                                uniform texobj2D amb_map_fine    : TEXUNIT1,
                                uniform texobj2D refl_texture    : TEXUNIT4,
                                uniform texobj1D fresnel_texture : TEXUNIT5,
                                uniform texobj2D int_index_map2  : TEXUNIT8,
                                uniform texobj3D int_norm_map2   : TEXUNIT12) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv2, int_index_map2, int_norm_map2);
  return r2o_opaque(input, pos, int_map_derivs,
                    amb_map_coarse, amb_map_fine, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_opaque_2(FpR2O input, in float4 pos: WPOS,
                                uniform texobj2D amb_map_coarse  : TEXUNIT0,
                                uniform texobj2D amb_map_fine    : TEXUNIT1,
                                uniform texobj2D refl_texture    : TEXUNIT4,
                                uniform texobj1D fresnel_texture : TEXUNIT5,
                                uniform texobj2D int_index_map1  : TEXUNIT7,
                                uniform texobj3D int_norm_map1   : TEXUNIT11) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.zw, int_index_map1, int_norm_map1);
  return r2o_opaque(input, pos, int_map_derivs,
                    amb_map_coarse, amb_map_fine, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_opaque_3(FpR2O input, in float4 pos: WPOS,
                                uniform texobj2D amb_map_coarse  : TEXUNIT0,
                                uniform texobj2D amb_map_fine    : TEXUNIT1,
                                uniform texobj2D refl_texture    : TEXUNIT4,
                                uniform texobj1D fresnel_texture : TEXUNIT5,
                                uniform texobj2D int_index_map1  : TEXUNIT7,
                                uniform texobj2D int_index_map2  : TEXUNIT8,
                                uniform texobj3D int_norm_map1   : TEXUNIT11,
                                uniform texobj3D int_norm_map2   : TEXUNIT12) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.zw, int_index_map1, int_norm_map1) +
                          r2o_IntMapDerivs(input.int_uv2,     int_index_map2, int_norm_map2);
  return r2o_opaque(input, pos, int_map_derivs,
                    amb_map_coarse, amb_map_fine, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_opaque_4(FpR2O input, in float4 pos: WPOS,
                                uniform texobj2D amb_map_coarse  : TEXUNIT0,
                                uniform texobj2D amb_map_fine    : TEXUNIT1,
                                uniform texobj2D refl_texture    : TEXUNIT4,
                                uniform texobj1D fresnel_texture : TEXUNIT5,
                                uniform texobj2D int_index_map0  : TEXUNIT6,
                                uniform texobj3D int_norm_map0   : TEXUNIT10) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.xy, int_index_map0, int_norm_map0);
  return r2o_opaque(input, pos, int_map_derivs,
                    amb_map_coarse, amb_map_fine, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_opaque_5(FpR2O input, in float4 pos: WPOS,
                                uniform texobj2D amb_map_coarse  : TEXUNIT0,
                                uniform texobj2D amb_map_fine    : TEXUNIT1,
                                uniform texobj2D refl_texture    : TEXUNIT4,
                                uniform texobj1D fresnel_texture : TEXUNIT5,
                                uniform texobj2D int_index_map0  : TEXUNIT6,
                                uniform texobj2D int_index_map2  : TEXUNIT8,
                                uniform texobj3D int_norm_map0   : TEXUNIT10,
                                uniform texobj3D int_norm_map2   : TEXUNIT12) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.xy, int_index_map0, int_norm_map0) +
                          r2o_IntMapDerivs(input.int_uv2,     int_index_map2, int_norm_map2);
  return r2o_opaque(input, pos, int_map_derivs,
                    amb_map_coarse, amb_map_fine, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_opaque_6(FpR2O input, in float4 pos: WPOS,
                                uniform texobj2D amb_map_coarse  : TEXUNIT0,
                                uniform texobj2D amb_map_fine    : TEXUNIT1,
                                uniform texobj2D refl_texture    : TEXUNIT4,
                                uniform texobj1D fresnel_texture : TEXUNIT5,
                                uniform texobj2D int_index_map0  : TEXUNIT6,
                                uniform texobj2D int_index_map1  : TEXUNIT7,
                                uniform texobj3D int_norm_map0   : TEXUNIT10,
                                uniform texobj3D int_norm_map1   : TEXUNIT11) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.xy, int_index_map0, int_norm_map0) +
                          r2o_IntMapDerivs(input.int_uv01.zw, int_index_map1, int_norm_map1);
  return r2o_opaque(input, pos, int_map_derivs,
                    amb_map_coarse, amb_map_fine, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__ENTRY__ half4 FP_r2o_opaque_7(FpR2O input, in float4 pos: WPOS,
                                uniform texobj2D amb_map_coarse  : TEXUNIT0,
                                uniform texobj2D amb_map_fine    : TEXUNIT1,
                                uniform texobj2D refl_texture    : TEXUNIT4,
                                uniform texobj1D fresnel_texture : TEXUNIT5,
                                uniform texobj2D int_index_map0  : TEXUNIT6,
                                uniform texobj2D int_index_map1  : TEXUNIT7,
                                uniform texobj2D int_index_map2  : TEXUNIT8,
                                uniform texobj3D int_norm_map0   : TEXUNIT10,
                                uniform texobj3D int_norm_map1   : TEXUNIT11,
                                uniform texobj3D int_norm_map2   : TEXUNIT12) : COLOR
{
  float2 int_map_derivs = r2o_IntMapDerivs(input.int_uv01.xy, int_index_map0, int_norm_map0) +
                          r2o_IntMapDerivs(input.int_uv01.zw, int_index_map1, int_norm_map1) +
                          r2o_IntMapDerivs(input.int_uv2,     int_index_map2, int_norm_map2);
  return r2o_opaque(input, pos, int_map_derivs,
                    amb_map_coarse, amb_map_fine, refl_texture, fresnel_texture);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__ENTRY__ half4 FP_r2o_downsample_fb(FpR2OTex input,  in float4 pos: WPOS,
                                     uniform texobj2D frame_buffer : TEXUNIT0) : COLOR
{
  return 0.125 * (tex2D(frame_buffer, input.uv.xy + float2(-3,-1)) +
                  tex2D(frame_buffer, input.uv.xy + float2(-1,-1)) +
                  tex2D(frame_buffer, input.uv.xy + float2( 1,-1)) +
                  tex2D(frame_buffer, input.uv.xy + float2( 3,-1)) +
                  tex2D(frame_buffer, input.uv.xy + float2(-3, 1)) +
                  tex2D(frame_buffer, input.uv.xy + float2(-1, 1)) +
                  tex2D(frame_buffer, input.uv.xy + float2( 1, 1)) +
                  tex2D(frame_buffer, input.uv.xy + float2( 3, 1)));
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__ENTRY__ half4 FP_r2o_downsample_zb(FpR2OTex input,  in float4 pos: WPOS,
                                     uniform texobj2D depth_buffer : TEXUNIT0) : COLOR
{
  float4 rgba = tex2D(depth_buffer, input.uv.xy);
  return rgba;


  //float4 rgba0 = tex2D(depth_buffer, input.uv.xy + float(-1,-1));
  //float4 rgba1 = tex2D(depth_buffer, input.uv.xy + float(-1, 1));
  //float4 rgba2 = tex2D(depth_buffer, input.uv.xy + float( 1,-1));
  //float4 rgba3 = tex2D(depth_buffer, input.uv.xy + float( 1, 1));
  //
  //float w0 = 1.0 / dot(rgba0, g_r2o_zbuffer_constants);
  //float w1 = 1.0 / dot(rgba1, g_r2o_zbuffer_constants);
  //float w2 = 1.0 / dot(rgba2, g_r2o_zbuffer_constants);
  //float w3 = 1.0 / dot(rgba3, g_r2o_zbuffer_constants);
  //
  //float w  = 0.25 * (w0 + w1 + w2 + w3);

  // w = 1.0 / ((rgba.x * const.x) + (rgba.y * const.y) + (rgba.z * const.z) + (rgba.w * const.w))
  // 
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__ENTRY__ half4 FP_r2o_blur(FpR2OTex input,  in float4 pos: WPOS,
                            uniform texobj2D src_tex : TEXUNIT0) : COLOR
{
  return 0.0625 * (tex2D(src_tex, input.uv.xy + float2(-3,-3)) +
                   tex2D(src_tex, input.uv.xy + float2(-3,-1)) +
                   tex2D(src_tex, input.uv.xy + float2(-3, 1)) +
                   tex2D(src_tex, input.uv.xy + float2(-3, 3)) +
                   tex2D(src_tex, input.uv.xy + float2(-1,-3)) +
                   tex2D(src_tex, input.uv.xy + float2(-1,-1)) +
                   tex2D(src_tex, input.uv.xy + float2(-1, 1)) +
                   tex2D(src_tex, input.uv.xy + float2(-1, 3)) +
                   tex2D(src_tex, input.uv.xy + float2( 1,-3)) +
                   tex2D(src_tex, input.uv.xy + float2( 1,-1)) +
                   tex2D(src_tex, input.uv.xy + float2( 1, 1)) +
                   tex2D(src_tex, input.uv.xy + float2( 1, 3)) +
                   tex2D(src_tex, input.uv.xy + float2( 3,-3)) +
                   tex2D(src_tex, input.uv.xy + float2( 3,-1)) +
                   tex2D(src_tex, input.uv.xy + float2( 3, 1)) +
                   tex2D(src_tex, input.uv.xy + float2( 3, 3)));
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__ENTRY__ half4 FP_r2o_blend_cube_map(FpR2OBlend input,  in float4 pos: WPOS,
                                      uniform texobj2D src_tex : TEXUNIT0) : COLOR
{
  float4 src_rgba = tex2D(src_tex, input.uv.xy);
  float4 cube_map_rgba;
  cube_map_rgba.xyz = SampleSpecularCubeMap(input.viewvec, 0.0f);
  cube_map_rgba.w = 0;
  return lerp(src_rgba, cube_map_rgba, g_r2o_cube_map_blend_ratio);
}

