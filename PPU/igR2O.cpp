#include <stdlib.h>

#include "igg/igGraphics.h"
#include "igBucketer/PPU/igBucketer.h"
#include "igShadersRsx/igVertexPrograms.h"
#include "igShadersRsx/igFragmentPrograms.h"
#include "igGameplay//igGpAlloc.h"
#include "igg/igShaders.h"
#include "igg/igEffectDraw.h"
#include "igg/igProfile.h"
#include "igg/igScene.h"
#include "igRender/igRender.h"
#include "igR2O.h"



namespace IGG
{

R2OCon        g_R2OCon;
u32           g_NormalMapPointers[MAX_LODS] ALIGNED(16);
f32           g_R2OMaxHeight   = 0.5f;
f32           g_MaxWaterDepth  = 4.0f;

u32 g_x0, g_y0, g_w0, g_h0, g_w1_frame, g_h1_frame, g_w1_depth, g_h1_depth;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  one-time initialization of r2o controller
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2OCon::Init(IGG::BasicAllocator* alloc, u32 max_water_objects)
{
  // setup bitfield allocator
  {
    u32 bitfield_byte_size = ((max_water_objects + 63) >> 6) << 3;
    m_object_alloc_bot = (u64*)g_SystemMemoryAlloc.Alloc(bitfield_byte_size, 128);
    Memset(m_object_alloc_bot, 0x0, bitfield_byte_size);
    Memset(m_object_alloc_top, 0x0, sizeof(m_object_alloc_top));
  }

  // setup the object array
  u32 size = max_water_objects * sizeof(R2OWaterObject);
  m_ea_water_objects = (u32)g_SystemMemoryAlloc.Alloc(size, 128 );
  FastMemSet8((void *)m_ea_water_objects, 0, size);


  m_r2o_enable                = true;

  m_time                      = 0.0f;
  m_paused                    = 0;
  m_DD_coeff0                 = -0.075f;
  m_allocator                 = (u32)&g_GPUAllocator;
  m_step                      = 128.0f;
  m_near0                     = 25.0f * m_step;
  m_num_lods                  = 38;
  m_near_lod                  = 28;
  m_max_verts                 = 9216;
  m_suppress_rows             = 6;
  m_dd_extra                  = 9;
  m_height_scale              = R2O_HEIGHT_SCALE;
  m_ea_normal_map_pointers    = (u32)g_NormalMapPointers;
  m_debug_fp                  = 0;
  m_coarse_map_weight         = 0.6f;
  m_fine_map_weight           = 0.3f;
  m_debug_band                = -1;
  m_frustum_fudge_factor      = 10.0f;
  m_index_of_refraction       = 1.332986f;
  m_clip_plane_height         = 0.0f;
  m_tie_scale_adjust          = 1.0f;
  m_bsphere_min_angle_degrees = 10.0f;
  m_cube_map_blend_ratio      = 0.2f;
  m_misc_consts[0]            = 0.1f;
  m_misc_consts[1]            = 0.2f;
  m_misc_consts[2]            = 4.0f;
  m_misc_consts[3]            = 0.5f;
  m_caustics_lod              = 20;
  m_ea_caustics_data          = (u32)g_NormalBuffer;
  m_ea_caustics_envelope      = (u32)g_CausticsEnvelope;
  m_ea_query_groups           = (u32)g_QueryGroups;
  #if R2O_TASKS
  m_ea_task_list              = (u32)g_R2OTasks;
  #endif
  m_max_water_objects         = max_water_objects;
  m_num_water_objects         = 0;
  m_min_rad                   = 1.0f;
  m_object_last               = -1; 
  m_timer_line                = 0;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  reset r2o controller to default state
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2OCon::Reset()
{
  // reset water objects
  FastMemSet8((void *)m_ea_water_objects, 0, m_max_water_objects * sizeof(R2OWaterObject));

  // reset bitfield
  u32 bitfield_byte_size = ((m_max_water_objects + 63) >> 6) << 3;
  Memset(m_object_alloc_bot, 0x0, bitfield_byte_size);
  Memset(m_object_alloc_top, 0x0, sizeof(m_object_alloc_top));

  m_num_water_objects = 0;
  m_object_last       = -1;
}



////////////////////////////////////////////////////////////////////////////////////////////////
//
// water object management -- allocate a new object and return a pointer to it
//
////////////////////////////////////////////////////////////////////////////////////////////////

R2OWaterObject *R2OCon::WaterObjectCreate()
{
  printf("WaterObjectCreate()\n");
  i32 index = GpAlloc(m_object_alloc_top, m_object_alloc_bot, m_max_water_objects);

  // handle case where we run out of objects
  if (index < 0)
  {
    return NULL;
  }

  R2OWaterObject *p_obj = (R2OWaterObject *)m_ea_water_objects + index;
  p_obj->m_flags = R2O_WATER_OBJECT_FLAG_ACTIVE;
  m_num_water_objects++;
  m_object_last = MaxI32(m_object_last, index);    // update m_last if this is the highest allocated object

  return p_obj;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//
// water object management -- allocate a new object and initialize it based on a built object
//
////////////////////////////////////////////////////////////////////////////////////////////////

R2OWaterObject *R2OCon::WaterObjectCreate(R2OBuiltWaterObject* p_built_obj)
{
  printf("WaterObjectCreate(%08X)\n", (u32)p_built_obj);
  R2OWaterObject *p_obj = WaterObjectCreate();

  if (p_obj)
  {
    // copy data members
    p_obj->m_origin         = p_built_obj->m_origin;
    p_obj->m_dimensions     = p_built_obj->m_dimensions;
    p_obj->m_water_color    = p_built_obj->m_water_color;
    p_obj->m_amplitude      = p_built_obj->m_amplitude;
    p_obj->m_first_waveband = p_built_obj->m_first_waveband;
    p_obj->m_material       = p_built_obj->m_material;
    p_obj->m_cubemap_index  = p_built_obj->m_cubemap_index;

    // some data adjustments I need to make on the engine side
    p_obj->m_first_waveband += 6;   // lod 0 was changed from 16m to 128m steps, 6 semi-octaves difference
    if (p_obj->m_first_waveband > 17)
    {
      // prevent those naughty artists from creating silly water objects
      p_obj->m_first_waveband = 17;
    }
    p_obj->m_amplitude *= 2.0f;   // doubled because I removed the odd lod's
  }

  return p_obj;
}



////////////////////////////////////////////////////////////////////////////////////////////////
//
// water object management -- delete the specified object
//
////////////////////////////////////////////////////////////////////////////////////////////////

void R2OCon::WaterObjectDelete(R2OWaterObject *p_obj)
{
  R2O_DeleteInteractive(&p_obj->m_interactive_data);

  i64 index = p_obj - (R2OWaterObject *)m_ea_water_objects;
  if (!GpIsAlloced(m_object_alloc_bot, m_max_water_objects, index))
  {
    return;
  }

  // clear the slot - which therefore clears the 'active' flag
  FastMemSet8(p_obj, 0, sizeof(R2OWaterObject));

  // free the index in the allocation bitfield
  GpFree(m_object_alloc_top, m_object_alloc_bot, index);

  // decrement the highest object marker
  if (index >= m_object_last)
  {
    m_object_last = GpAllocLast(m_object_alloc_bot, m_object_last);
  }

  // put here any callbacks for notifying other systems that a water object got deleted
  if (g_R2OHeightQueryCallback)
  {
    g_R2OHeightQueryCallback(p_obj);
  }

  m_num_water_objects--;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_Init(u32 interactive_tile_count)
{
  R2O_InitJob();
  R2O_InitAmbient();
  R2O_InitInteractive(interactive_tile_count);
  R2O_GenerateFresnelTexture();
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_FrameInit()
{
  R2O_ResetQueryGroups();
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_AddWaterObjects(R2OBuiltWaterObject* built_water_objects, u32 num_water_objects, i64 zone_index)
{
  if (g_R2OCon.m_num_water_objects + num_water_objects > g_R2OCon.m_max_water_objects)
  {
    IG_ASSERT(0);
    return;
  }

  for (u32 i=0; i<num_water_objects; i++)
  {
    R2OWaterObject *p_obj = g_R2OCon.WaterObjectCreate(&built_water_objects[i]);
    IG_ASSERT(p_obj);

    p_obj->m_zone_index = zone_index;

    // disable the new water objects until they've had a chance to be updated once
    u32 idx = p_obj - (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
    g_R2OCon.m_water_object_disable_bits |= 1ULL << idx;
  }
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_DeleteUnusedWaterObjects(i64 zone_alloc_mask)
{
  R2OWaterObject *p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++, p_obj++)
  {
    if (p_obj->m_flags & R2O_WATER_OBJECT_FLAG_ACTIVE)
    {
      i64 zone_index = p_obj->m_zone_index;
      i64 zone_index_mask = 0x8000000000000000ULL >> zone_index;
      if (!(zone_index_mask & zone_alloc_mask))
      {
        g_R2OCon.WaterObjectDelete(p_obj);
      }
    }
  }
}

void R20_DisableZoneWaterObjects(i64 zone_index)
{
  R2OWaterObject *p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++, p_obj++)
  {
    if (p_obj->m_zone_index == zone_index)
    {
      p_obj->m_flags &= (~R2O_WATER_OBJECT_FLAG_ACTIVE);
    }
  }
}

void R20_EnableZoneWaterObjects(i64 zone_index)
{
  R2OWaterObject *p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++, p_obj++)
  {
    if (p_obj->m_zone_index == zone_index)
    {
      p_obj->m_flags |= R2O_WATER_OBJECT_FLAG_ACTIVE;
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_Update(f32 time_step, u32 num_viewports)
{
  //----------------------------------------------------------------------------

  g_R2OCon.m_time_step = g_R2OCon.m_paused ? 0.0f : time_step;

  // re-enable any water objects that were disabled until the next update
  g_R2OCon.m_water_object_disable_bits = 0;


  //----------------------------------------------------------------------------

  if (g_R2OCon.m_num_water_objects)
  {
    R2O_AssignImpulsesToWaterObjects();
  }


  //----------------------------------------------------------------------------

  g_R2OCon.m_num_viewports = num_viewports;

  // make some extra helper vars from the view data
  for (u32 v=0; v<num_viewports; v++)
  {
    R2OViewData *p_view_data = &g_R2OCon.m_view_data[v];

    vec4f &near_plane   = p_view_data->m_frustum_planes[FRUSTUM_PLANE_NEAR];
    vec4f &left_plane   = p_view_data->m_frustum_planes[FRUSTUM_PLANE_L];
    vec4f &right_plane  = p_view_data->m_frustum_planes[FRUSTUM_PLANE_R];
    vec4f &bottom_plane = p_view_data->m_frustum_planes[FRUSTUM_PLANE_B];
    vec4f &top_plane    = p_view_data->m_frustum_planes[FRUSTUM_PLANE_T];
    vec4f &far_plane    = p_view_data->m_frustum_planes[FRUSTUM_PLANE_FAR];

    // derive the near and far plane distances
    f32 n = -VecDot4(near_plane, p_view_data->m_camera_position);
    f32 f =  VecDot4(far_plane,  p_view_data->m_camera_position);

    // derive the projection matrix from the frustum
    mtx4f proj_mat;
    f32 l,r,b,t,dot;

    dot = VecDot3(near_plane, left_plane);
    l = n * dot / Sqrtf(1.0f - dot*dot);

    dot = VecDot3(near_plane, right_plane);
    r = -n * dot / Sqrtf(1.0f - dot*dot);

    dot = VecDot3(near_plane, bottom_plane);
    b = -n * dot / Sqrtf(1.0f - dot*dot);

    dot = VecDot3(near_plane, top_plane);
    t = n * dot / Sqrtf(1.0f - dot*dot);

    proj_mat.Set(
        (2*n)/(r-l),          0,                 0,                0,
             0,          (2*n)/(t-b),            0,                0,
        (r+l)/(l-r),     (t+b)/(b-t),          f /(f-n),           1,
             0,               0,              f*n/(n-f),           0);

    p_view_data->m_frag_texproj_matrix = proj_mat;

    // set frustum params
    p_view_data->m_frustum_near   = n;
    p_view_data->m_frustum_far    = f;
    p_view_data->m_frustum_hratio = (l-r)/(2.0f*n);
    p_view_data->m_frustum_vratio = (t-b)/(2.0f*n);

    // set the vp index
    p_view_data->m_vp_index = v;
  }

  //----------------------------------------------------------------------------
  #if R2O_TASKS
  // build the spu task list
  R2O_BuildTaskList();
  #endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define FRESNEL_TEXTURE_SIZE 256
IGTexture g_R2OFresnelTexture;
u16 *g_R2OFresnelVram;

void R2O_GenerateFresnelTexture()
{
  g_R2OFresnelVram = (u16 *)g_SystemMemoryAlloc.Alloc(FRESNEL_TEXTURE_SIZE * sizeof(u16));

  for (u32 i=0; i<FRESNEL_TEXTURE_SIZE; i++)
  {
    f32 d  = (f32)i * 2.0f/((f32)(FRESNEL_TEXTURE_SIZE-1)) - 1.0f;
    f32 eta;
    if (d<0)
    {
      eta = 1.0f/g_R2OCon.m_index_of_refraction;
      d = -d;
    }
    else
    {
      eta = g_R2OCon.m_index_of_refraction;
    }
    f32 e2 = eta * eta;
    f32 s2 = e2 - 1.0f + d*d;

    f32 f;
    if (s2 < 0)
    {
      // total internal reflection
      f = 1.0f;
    }
    else
    {
      f32 s  = Sqrtf(s2);
      f32 ax = e2*d - s;
      f32 bx = e2*d + s;
      f32 ay =    d - s;
      f32 by =    d + s;
      f32 cx = ax/bx;
      f32 cy = ay/by;

      f = 0.5f * (cx*cx + cy*cy);
    }

    g_R2OFresnelVram[i] = (u16)(f * 65535.0f);
  }

  // set IGTexture fields
  g_R2OFresnelTexture.m_baseOffset = MainMemoryAddressToOffset(g_R2OFresnelVram) & 0x7fffffff;
  g_R2OFresnelTexture.m_control2   = 0x80000000;
  g_R2OFresnelTexture.m_filter     = 0x02022000;
  g_R2OFresnelTexture.m_size1      = (FRESNEL_TEXTURE_SIZE<<16) | 1;
  g_R2OFresnelTexture.m_size2      = 0x00100000 | (FRESNEL_TEXTURE_SIZE << 1);
  g_R2OFresnelTexture.m_control1   = 0x00030303;
  g_R2OFresnelTexture.m_format     = 0x0001801A | (TEXTUREFORMAT_R16 << 8);
  g_R2OFresnelTexture.m_swizzle    = 0xAAEE;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define NORMAL_MAP_WIDTH   64
#define NORMAL_MAP_HEIGHT  64
#define NORMAL_MAP_FORMAT  TEXTUREFORMAT_G8B8
#define NORMAL_MAP_BYPP    2
#define NORMAL_MAP_SWIZZLE 0x1B

void R2O_SetNormalMapTexture(u32 tex_unit, u32 ea)
{
  if (ea)
  {
    u32 base = AddressToOffset((void *)ea) & 0x7FFFFFFF;
    u32 format   = 1<<16 | 0x7<<13 | NORMAL_MAP_FORMAT<<8 | 2<<4 | 0xA;
    u32 control1 = TEXTUREADDRESS_CLAMPTOEDGE<<16 | TEXTUREADDRESS_WRAP<<8 | TEXTUREADDRESS_WRAP;
    u32 control2 = 1<<31 | 0<<27 | 0<<15;
    u32 swizzle  = 0xAA<<8 | NORMAL_MAP_SWIZZLE;
    u32 filter   = 0xF<<28 | TEXTUREFILTER_MAG_LINEAR <<24 | TEXTUREFILTER_MIN_LINEAR <<16 | 1<<13 | 0x1C00;
    u32 size1    = NORMAL_MAP_WIDTH<<16 | NORMAL_MAP_HEIGHT;
    u32 size2    = 1<<20 | (NORMAL_MAP_BYPP * NORMAL_MAP_WIDTH);

    IGTexture tex = IGTexture(base, format, control1, control2, swizzle, filter, size1, size2);
    g_CurrentPB->IGPushBuffer::SetTexture(tex_unit, &tex);
  }
  else
  {
    g_CurrentPB->IGPushBuffer::SetTexture(tex_unit, &g_BlackTexture);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define DOWNSAMPLE_SHIFT_FRAME_X 2    // relative to resolved width
#define DOWNSAMPLE_SHIFT_FRAME_Y 2

#define DOWNSAMPLE_SHIFT_DEPTH_X 1    // relative to resolved width
#define DOWNSAMPLE_SHIFT_DEPTH_Y 1

#define DOWNSAMPLE_FRAME_TO_16BIT 0
#define BLUR_FRAME_TO_16BIT       0
#define READ_DEPTH_AS_G16R16      1

// allocate the memory for the reduced frame buffer and z buffer textures
bool R2O_AllocVram(u32 &fb_mem, u32 &zb_mem, u32 &blur_mem, u32 &blur_mem2)
{
  u32 down_bytes = DOWNSAMPLE_FRAME_TO_16BIT ? 2 : 4;
  u32 blur_bytes = BLUR_FRAME_TO_16BIT       ? 2 : 4;

  // alloc downsampled frame buffer
  void* fb_vram_addr = g_ScratchVram.Alloc(g_w1_frame * g_h1_frame * down_bytes, RSX_ALIGNMENT_RENDERTARGET);

  // alloc downsampled depth buffer
  void *zb_vram_addr=fb_vram_addr;    // initialize it to something sensible
  if (g_R2OCon.m_any_transparent_visible)
  {
    zb_vram_addr = g_ScratchVram.Alloc(g_w1_depth * g_h1_depth * 4, RSX_ALIGNMENT_RENDERTARGET);
  }

  // alloc mem for blurred versions of downsampled frame buffer
  void* blur_vram_addr  = g_ScratchVram.Alloc(128 * 64 * blur_bytes, RSX_ALIGNMENT_RENDERTARGET);
  void* blur_vram_addr2 = g_ScratchVram.Alloc( 64 * 32 * blur_bytes, RSX_ALIGNMENT_RENDERTARGET);

  if (g_ScratchVram.FreeMemory() == 0)
  {
    return false;
  }

  fb_mem    = AddressToOffset(fb_vram_addr);
  zb_mem    = AddressToOffset(zb_vram_addr);
  blur_mem  = AddressToOffset(blur_vram_addr);
  blur_mem2 = AddressToOffset(blur_vram_addr2);

  // allocs ok
  return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_DownsampleFrameBuffer(IGTexture *p_fb_tex, u32 fb_mem, IGTexture *p_src_tex)
{
  // A) Obtain the source surface as a texture
  IGTexture src_tex = *p_src_tex;
  src_tex.m_format |= 0x4000;
  src_tex.m_filter = 0x02022000;
  g_CurrentPB->IGPushBuffer::SetTexture(0, &src_tex);
  g_CurrentPB->InvalidateTextureCachedState(0);

  // B) Bind the destination memory as a render target
  IGRenderTarget rt;
  u32 stride = g_w1_frame << (DOWNSAMPLE_FRAME_TO_16BIT ? 1 : 2);
  IGColorBufferFormat format = DOWNSAMPLE_FRAME_TO_16BIT ? COLORBUFFER_RGB565 : COLORBUFFER_ARGB8888;
  rt.InitNoDepth(g_w1_frame, g_h1_frame, format, DEPTHBUFFER_D24S8, MULTISAMPLE_NONE, RENDERTARGET_LINEAR, fb_mem, stride);
  g_CurrentPB->SetRenderTarget(&rt);

  // C) Set the viewport, scissor rect and full screen scale factors.
  SetScreenSpaceScale(g_CurrentPB,  g_w1_frame, g_h1_frame);
  g_CurrentPB->SetViewport   (0, 0, g_w1_frame, g_h1_frame);
  g_CurrentPB->SetScissorRect(0, 0, g_w1_frame, g_h1_frame);

  // D) Set any render state you want to use [Z writes and Z test off is important if its just a color surface]
  g_CurrentPB->SetRenderState(RENDERSTATE_DEPTHTEST,  false);
  g_CurrentPB->SetRenderState(RENDERSTATE_DEPTHWRITE, false);
  g_CurrentPB->SetRenderState(RENDERSTATE_CULLFACE,   false);
  g_CurrentPB->SetRenderState(RENDERSTATE_BLEND,      false);

  // E) Bind the screen space vertex program
  g_CurrentPB->DisableVertexAttributes();
  g_CurrentPB->SetVertexProgram(g_VertexPrograms[VERTEX_PROG_SCREENSPACE]);

  // F) Bind the fragment program you want to use
  g_CurrentPB->SetFragmentProgram(g_FragmentPrograms[FRAGMENT_PROG_R2O_DOWNSAMPLE_FB]);

  // G) Draw a full screen quad.
  u32 ms = (g_RenderTarget[TARGET_RENDER_BUFFER].m_targetFormat >> 12) & 0xF;
  float ms_scl = ms ? 2.f : 1.f;    // must double the u-coord because the frame buffer is 2x sampled
  {
    g_CurrentPB->DrawBegin(PRIM_QUADLIST);

    // top left
    g_CurrentPB->SetVertexAttrib4f(1, ms_scl*g_x0, g_y0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0,  0, 0, 0, 1);

    // bottom left
    g_CurrentPB->SetVertexAttrib4f(1, ms_scl*g_x0, g_y0+g_h0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0,  0, g_h1_frame, 0, 1);

    // bottom right
    g_CurrentPB->SetVertexAttrib4f(1, ms_scl*(g_x0+g_w0), g_y0+g_h0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0, g_w1_frame, g_h1_frame, 0, 1);

    // top right
    g_CurrentPB->SetVertexAttrib4f(1, ms_scl*(g_x0+g_w0), g_y0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0, g_w1_frame, 0, 0, 1);

    g_CurrentPB->DrawEnd();
  }

  // set IGTexture fields
  p_fb_tex->m_baseOffset = fb_mem;
  p_fb_tex->m_control2   = 0x80000000;
  p_fb_tex->m_filter     = 0x02022000;
  p_fb_tex->m_size1      =  (g_w1_frame<<16)  | g_h1_frame;
  p_fb_tex->m_size2      = 0x00100000 | stride;
  p_fb_tex->m_control1   = 0x00720202;  // '7' to enable gamma
  #if DOWNSAMPLE_FRAME_TO_16BIT
  p_fb_tex->m_format     = 0x0001A029 | (TEXTUREFORMAT_R5G6B5 << 8);
  #else
  p_fb_tex->m_format     = 0x0001A029 | (TEXTUREFORMAT_A8R8G8B8 << 8);
  #endif
  p_fb_tex->m_swizzle    = p_src_tex->m_swizzle;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_DownsampleDepthBuffer(IGTexture *p_zb_tex, u32 zb_mem, IGTexture *p_src_tex)
{
  g_CurrentPB->SetRenderState(RENDERSTATE_GAMMACORRECTION, 0);

  u32 stride = g_w1_depth << 2;

  // A) Obtain the source surface as a texture
  IGTexture src_tex = *p_src_tex;
  src_tex.m_format  = 0x0001E529;
  src_tex.m_filter  = 0x01012000;
  src_tex.m_swizzle = 0x0000AAE4;
  g_CurrentPB->IGPushBuffer::SetTexture(0, &src_tex);
  g_CurrentPB->InvalidateTextureCachedState(0);

  // B) Bind the destination memory as a render target
  IGRenderTarget rt;
  rt.InitNoDepth(g_w1_depth, g_h1_depth, COLORBUFFER_ARGB8888, DEPTHBUFFER_D24S8, MULTISAMPLE_NONE, RENDERTARGET_LINEAR, zb_mem, stride);
  g_CurrentPB->SetRenderTarget(&rt);

  // C) Set the viewport, scissor rect and full screen scale factors.
  SetScreenSpaceScale(g_CurrentPB,  g_w1_depth, g_h1_depth);
  g_CurrentPB->SetViewport   (0, 0, g_w1_depth, g_h1_depth);
  g_CurrentPB->SetScissorRect(0, 0, g_w1_depth, g_h1_depth);

  // D) Set any render state you want to use [Z writes and Z test off is important if its just a color surface]
  g_CurrentPB->SetRenderState(RENDERSTATE_DEPTHTEST,  false);
  g_CurrentPB->SetRenderState(RENDERSTATE_DEPTHWRITE, false);
  g_CurrentPB->SetRenderState(RENDERSTATE_CULLFACE,   false);
  g_CurrentPB->SetRenderState(RENDERSTATE_BLEND,      false);

  // E) Bind the screen space vertex program
  g_CurrentPB->DisableVertexAttributes();
  g_CurrentPB->SetVertexProgram(g_VertexPrograms[VERTEX_PROG_SCREENSPACE]);

  // F) Bind the fragment program you want to use
  g_CurrentPB->SetFragmentProgram(g_FragmentPrograms[FRAGMENT_PROG_R2O_DOWNSAMPLE_ZB]);

  // G) Draw a full screen quad.
  u32 ms = (g_RenderTarget[TARGET_RENDER_BUFFER].m_targetFormat >> 12) & 0xF;
  float ms_scl = ms ? 2.f : 1.f;    // must double the u-coord because the frame buffer is 2x sampled
  {
    g_CurrentPB->DrawBegin(PRIM_QUADLIST);

    // top left
    g_CurrentPB->SetVertexAttrib4f(1, ms_scl*g_x0, g_y0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0,  0, 0, 0, 1);

    // bottom left
    g_CurrentPB->SetVertexAttrib4f(1, ms_scl*g_x0, g_y0+g_h0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0,  0, g_h1_depth, 0, 1);

    // bottom right
    g_CurrentPB->SetVertexAttrib4f(1, ms_scl*(g_x0+g_w0), g_y0+g_h0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0, g_w1_depth, g_h1_depth, 0, 1);

    // top right
    g_CurrentPB->SetVertexAttrib4f(1, ms_scl*(g_x0+g_w0), g_y0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0, g_w1_depth, 0, 0, 1);

    g_CurrentPB->DrawEnd();
  }

  // set IGTexture fields
  p_zb_tex->m_baseOffset = zb_mem;
  #if READ_DEPTH_AS_G16R16
  p_zb_tex->m_format     = TEXTURE_FORMAT(1, 0, 1, TEXTUREFORMAT_G16R16,   2, 1, 0, 0, 1);
  #else
  p_zb_tex->m_format     = TEXTURE_FORMAT(1, 0, 1, TEXTUREFORMAT_A8R8G8B8, 2, 1, 0, 0, 1);
  #endif
  p_zb_tex->m_control1   = TEXTURE_CONTROL1(0, 0, TEXTUREADDRESS_WRAP, 0, TEXTUREADDRESS_MIRROR, 0, TEXTUREADDRESS_MIRROR);
  p_zb_tex->m_control2   = TEXTURE_CONTROL2(1, 0, 0, 0, 0, 0);
  p_zb_tex->m_swizzle    = TEXTURE_SWIZZLE(0, 0x6A, 0xE4);
  p_zb_tex->m_filter     = TEXTURE_FILTER(0x0, TEXTUREFILTER_MAG_NEAREST, TEXTUREFILTER_MIN_NEAREST, 1, 0);
  p_zb_tex->m_size1      = TEXTURE_SIZE1(g_w1_depth, g_h1_depth);
  p_zb_tex->m_size2      = TEXTURE_SIZE2(1, stride);

  g_CurrentPB->SetRenderState(RENDERSTATE_GAMMACORRECTION, 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_RestoreViewport()
{
  g_CurrentPB->SetRenderState(RENDERSTATE_DEPTHTEST,  true);
  g_CurrentPB->SetRenderState(RENDERSTATE_DEPTHWRITE, true);
  g_CurrentPB->SetRenderState(RENDERSTATE_CULLFACE,   true);
  g_CurrentPB->SetRenderTarget(&g_RenderTarget[TARGET_RENDER_BUFFER]);
  g_CurrentPB->SetViewport(g_x0, g_y0, g_w0, g_h0,1.0f,false,0.0f,g_CurrentDrawBucketer->m_viewport->m_standard_depth_range);
  g_CurrentPB->SetScissorRect(g_x0, g_y0, g_w0, g_h0);
  SetScreenSpaceScale(g_CurrentPB, g_w0, g_h0);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_DownsampleBuffers(IGTexture *p_fb_tex, u32 fb_mem, IGTexture *p_zb_tex, u32 zb_mem)
{
  // downsample frame buffer
  R2O_DownsampleFrameBuffer(p_fb_tex, fb_mem, &g_RenderTargetTextures[TARGET_RENDER_BUFFER]);

  // z buffer
  if (g_R2OCon.m_any_transparent_visible)
  {
  R2O_DownsampleDepthBuffer(p_zb_tex, zb_mem, &g_DepthTexture);
  }

  // restore settings
  R2O_RestoreViewport();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_BlurTexture(IGTexture *p_blur_tex, u32 blur_mem, IGTexture *p_src_tex, u32 dest_width, u32 dest_height)
{
  u32 w0 = p_src_tex->m_size1 >> 16;
  u32 h0 = p_src_tex->m_size1 & 0x0000FFFF;
  u32 w1 = dest_width;
  u32 h1 = dest_height;

  // A) Obtain the source surface as a texture
  IGTexture src_tex = *p_src_tex;
  src_tex.m_format |= 0x4000;
  src_tex.m_filter = 0x02022000;
  src_tex.m_control1 = 0x00730303;  // '7' to enable gamma
  g_CurrentPB->IGPushBuffer::SetTexture(0, &src_tex);
  g_CurrentPB->InvalidateTextureCachedState(0);

  // B) Bind the destination memory as a render target
  IGRenderTarget rt;
  u32 stride = w1 << (BLUR_FRAME_TO_16BIT ? 1 : 2);
  IGColorBufferFormat color_format = BLUR_FRAME_TO_16BIT ? COLORBUFFER_RGB565 : COLORBUFFER_ARGB8888;
  rt.InitNoDepth(w1, h1, color_format, DEPTHBUFFER_D24S8, MULTISAMPLE_NONE, RENDERTARGET_LINEAR, blur_mem, stride);
  g_CurrentPB->SetRenderTarget(&rt);

  // C) Set the viewport, scissor rect and full screen scale factors.
  SetScreenSpaceScale(g_CurrentPB,  w1, h1);
  g_CurrentPB->SetViewport   (0, 0, w1, h1);
  g_CurrentPB->SetScissorRect(0, 0, w1, h1);

  // D) Set any render state you want to use [Z writes and Z test off is important if its just a color surface]
  g_CurrentPB->SetRenderState(RENDERSTATE_DEPTHTEST,  false);
  g_CurrentPB->SetRenderState(RENDERSTATE_DEPTHWRITE, false);
  g_CurrentPB->SetRenderState(RENDERSTATE_CULLFACE,   false);
  g_CurrentPB->SetRenderState(RENDERSTATE_BLEND,      false);

  // E) Bind the screen space vertex program
  g_CurrentPB->DisableVertexAttributes();
  g_CurrentPB->SetVertexProgram(g_VertexPrograms[VERTEX_PROG_SCREENSPACE]);

  // F) Bind the fragment program you want to use
  g_CurrentPB->SetFragmentProgram(g_FragmentPrograms[FRAGMENT_PROG_R2O_BLUR]);

  // G) Draw a full screen quad.
  {
    g_CurrentPB->DrawBegin(PRIM_QUADLIST);

    // top left
    g_CurrentPB->SetVertexAttrib4f(1,  0,  0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0,  0,  0, 0, 1);

    // bottom left
    g_CurrentPB->SetVertexAttrib4f(1,  0, h0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0,  0, h1, 0, 1);

    // bottom right
    g_CurrentPB->SetVertexAttrib4f(1, w0, h0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0, w1, h1, 0, 1);

    // top right
    g_CurrentPB->SetVertexAttrib4f(1, w0,  0, 0, 1);
    g_CurrentPB->SetVertexAttrib4f(0, w1,  0, 0, 1);

    g_CurrentPB->DrawEnd();
  }

  // restore settings
  R2O_RestoreViewport();

  // set IGTexture fields
  IGTextureFormat tex_format = BLUR_FRAME_TO_16BIT ? TEXTUREFORMAT_R5G6B5 : TEXTUREFORMAT_A8R8G8B8;
  p_blur_tex->m_baseOffset = blur_mem;
  p_blur_tex->m_control2   = 0x80000000;
  p_blur_tex->m_filter     = 0x02022000;
  p_blur_tex->m_size1      = (w1<<16)   | h1;
  p_blur_tex->m_size2      = 0x00100000 | stride;
  p_blur_tex->m_control1   = 0x00720202;  // '7' to enable gamma
  p_blur_tex->m_format     = 0x0001A029 | (tex_format << 8);
  p_blur_tex->m_swizzle    = p_src_tex->m_swizzle;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_BlendCubeMaps(IGTexture **pp_blended_texs, u8 *p_cubemap_idx, u32 num_blended_tex, IGTexture *p_src_tex,
                       R2OViewData *p_view_data)
                       
{
  u32 w0 = p_src_tex->m_size1 >> 16;
  u32 h0 = p_src_tex->m_size1 & 0x0000FFFF;
  u32 w1 = w0;
  u32 h1 = h0;

  // A) Obtain the source surface as a texture
  IGTexture src_tex = *p_src_tex;
  src_tex.m_format |= 0x4000;
  src_tex.m_filter = 0x02022000;
  src_tex.m_control1 = 0x00730303;  // '7' to enable gamma
  g_CurrentPB->IGPushBuffer::SetTexture(0, &src_tex);
  g_CurrentPB->InvalidateTextureCachedState(0);

  // C) Set the viewport, scissor rect and full screen scale factors.
  SetScreenSpaceScale(g_CurrentPB,  w1, h1);
  g_CurrentPB->SetViewport   (0, 0, w1, h1);
  g_CurrentPB->SetScissorRect(0, 0, w1, h1);

  // D) Set any render state you want to use [Z writes and Z test off is important if its just a color surface]
  g_CurrentPB->SetRenderState(RENDERSTATE_DEPTHTEST,  false);
  g_CurrentPB->SetRenderState(RENDERSTATE_DEPTHWRITE, false);
  g_CurrentPB->SetRenderState(RENDERSTATE_CULLFACE,   false);
  g_CurrentPB->SetRenderState(RENDERSTATE_BLEND,      false);

  // E) Bind the screen space vertex program
  g_CurrentPB->DisableVertexAttributes();
  g_CurrentPB->SetVertexProgram(g_VertexPrograms[VERTEX_PROG_SCREENSPACE]);

  // F) Bind the fragment program you want to use
  g_FragmentPrograms[FRAGMENT_PROG_R2O_BLEND_CUBE_MAP]->SetConstant(FP_CONST_R2O_CUBE_MAP_BLEND_RATIO,
                                                                    (float *)&g_R2OCon.m_cube_map_blend_ratio);
  g_CurrentPB->SetFragmentProgram(g_FragmentPrograms[FRAGMENT_PROG_R2O_BLEND_CUBE_MAP]);

  // compute corner direction vecs in camera space
  vec4f c0, c1, c2, c3;
  f32 h = p_view_data->m_frustum_hratio;
  f32 v = p_view_data->m_frustum_vratio;
  c0.SetXYZW( h, v, 1, 0);   // top left
  c1.SetXYZW( h,-v, 1, 0);   // bottom left
  c2.SetXYZW(-h,-v, 1, 0);   // bottom right
  c3.SetXYZW(-h, v, 1, 0);   // top right

  // rotate them into world space
  mtx4f camera_to_world_rot = p_view_data->m_world_to_camera_matrix;
  camera_to_world_rot.Transpose();
  MatMulVec(c0, camera_to_world_rot, c0);
  MatMulVec(c1, camera_to_world_rot, c1);
  MatMulVec(c2, camera_to_world_rot, c2);
  MatMulVec(c3, camera_to_world_rot, c3);

  for(u32 i = 0; i < num_blended_tex; ++i)
  {
    //Bind the destination memory as a render target
    IGRenderTarget  rt;
    u32             blended_mem;
    void*           blended_vram_addr;
    blended_vram_addr = g_ScratchVram.Alloc( 64 * 32 * 2, RSX_ALIGNMENT_RENDERTARGET);
    blended_mem       = AddressToOffset(blended_vram_addr);

    u32 stride = w1 << 1;
    rt.InitNoDepth(w1, h1, COLORBUFFER_RGB565, DEPTHBUFFER_D24S8, MULTISAMPLE_NONE, RENDERTARGET_LINEAR, blended_mem, stride);
    g_CurrentPB->SetRenderTarget(&rt);

    SetCubeMapTexture(g_CurrentPB, &g_ShaderDatabase.m_cube_maps[p_cubemap_idx[i]], CUBE_SPEC_VERY_TIGHT);

    // G) Draw a full screen quad.
    {
      g_CurrentPB->DrawBegin(PRIM_QUADLIST);
   
      // top left
      g_CurrentPB->SetVertexAttrib4f(1,  0,  0, 0, 1);
      g_CurrentPB->SetVertexAttrib4f(2, c0.x,-c0.y, c0.z, 1);
      g_CurrentPB->SetVertexAttrib4f(0,  0,  0, 0, 1);
   
      // bottom left
      g_CurrentPB->SetVertexAttrib4f(1,  0, h0, 0, 1);
      g_CurrentPB->SetVertexAttrib4f(2, c1.x,-c1.y, c1.z, 1);
      g_CurrentPB->SetVertexAttrib4f(0,  0, h1, 0, 1);
   
      // bottom right
      g_CurrentPB->SetVertexAttrib4f(1, w0, h0, 0, 1);
      g_CurrentPB->SetVertexAttrib4f(2, c2.x,-c2.y, c2.z, 1);
      g_CurrentPB->SetVertexAttrib4f(0, w1, h1, 0, 1);
   
      // top right
      g_CurrentPB->SetVertexAttrib4f(1, w0,  0, 0, 1);
      g_CurrentPB->SetVertexAttrib4f(2, c3.x,-c3.y, c3.z, 1);
      g_CurrentPB->SetVertexAttrib4f(0, w1,  0, 0, 1);
   
   
      g_CurrentPB->DrawEnd();
    }

    // set IGTexture fields
    pp_blended_texs[i]->m_baseOffset = blended_mem;
    pp_blended_texs[i]->m_control2   = 0x80000000;
    pp_blended_texs[i]->m_filter     = 0x02022000;
    pp_blended_texs[i]->m_size1      = (w1<<16)   | h1;
    pp_blended_texs[i]->m_size2      = 0x00100000 | stride;
    pp_blended_texs[i]->m_control1   = 0x00720202;  // '7' to enable gamma
    pp_blended_texs[i]->m_format     = 0x0001A029 | (TEXTUREFORMAT_R5G6B5 << 8);
    pp_blended_texs[i]->m_swizzle    = p_src_tex->m_swizzle;
  }

  // restore settings
  R2O_RestoreViewport();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SetBasis(vec4f &basis0, vec4f &basis1, vec4f &alt_basis0, vec4f &alt_basis1, i32 lod)
{
  const f32 r = 0.707106781187f;

  if (lod&1)
  {
    basis0.SetXYZ( r,0,r);
    basis1.SetXYZ(-r,0,r);

    alt_basis0.SetXYZ(1,0,0);
    alt_basis1.SetXYZ(0,0,1);
  }
  else
  {
    basis0.SetXYZ(1,0,0);
    basis1.SetXYZ(0,0,1);

    alt_basis0.SetXYZ( r,0,r);
    alt_basis1.SetXYZ(-r,0,r);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_VisibilityTests(Frustum *frustum_pointers[], u32 num_views)
{
  // loop over the frusta
  u32 any_visible = false;
  u32 any_transparent_visible = false;
  for (u32 vp_idx=0; vp_idx<num_views; vp_idx++)
  {
    // get the frustum pointer
    Frustum *p_frustum = frustum_pointers[vp_idx];

    // and a pointer to the R2O view data
    R2OViewData *p_view_data = &g_R2OCon.m_view_data[vp_idx];

    // copy any necessary info from the frustum to the R2O view data
    p_view_data->m_world_to_camera_matrix = p_frustum->m_view_matrix;
    p_view_data->m_camera_position        = p_frustum->m_position;
    for (u32 i=0; i<6; i++)
    {
      p_view_data->m_frustum_planes[i]    = p_frustum->m_planes[i];
    }


    // loop through the water objects and cull them against the frustum
    u32 any_visible_in_view = false;
    u32 any_transparent_visible_in_view = false;
    R2OWaterObject* p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
    for (i64 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++, p_obj++)
    {
      if (!(p_obj->m_flags & R2O_WATER_OBJECT_FLAG_ACTIVE))
      {
        continue;
      }

      vec4f aabb_min, aabb_max, aabb_dims;
      aabb_min = p_obj->m_origin;
      AddVector3(aabb_max, aabb_min, p_obj->m_dimensions);

      // fix me!!! (I am a nasty cheat)
      aabb_min.y -= 1.0f;
      aabb_max.y += 1.0f;

      u32 visible = p_frustum->FrustumTrueCullAABB(aabb_min, aabb_max);
      p_obj->m_render_data[vp_idx].m_visible = visible;
      any_visible_in_view |= visible;
      any_transparent_visible_in_view |= (visible && (p_obj->m_water_color.w < 100.0f));
    }

    g_R2OCon.m_view_data[vp_idx].m_any_visible = any_visible_in_view;
    any_visible |= any_visible_in_view;
    any_transparent_visible |= any_transparent_visible_in_view;
  }

  g_R2OCon.m_any_visible = any_visible;
  g_R2OCon.m_any_transparent_visible = any_transparent_visible;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ProcessLayerR2O(BucketEntry* entry, u32 entry_cnt)
{
  ProfilerGPUEvent(PROF_GPU_R2O);

  // get viewport context
  R2OViewData *p_view_data = (R2OViewData *)entry->m_environment;
  u32 vp_idx = p_view_data->m_vp_index;
  if (!p_view_data->m_any_visible)
  {
    return;
  }


  // determine source and dest dimensions for downsamples
  g_x0 = g_CurrentDrawBucketer->m_viewport->m_xpos;
  g_y0 = g_CurrentDrawBucketer->m_viewport->m_ypos;
  g_w0 = g_CurrentDrawBucketer->m_viewport->m_width;
  g_h0 = g_CurrentDrawBucketer->m_viewport->m_height;

  // frame buffer dest dimensions
  g_w1_frame = g_w0 >> DOWNSAMPLE_SHIFT_FRAME_X;
  g_h1_frame = g_h0 >> DOWNSAMPLE_SHIFT_FRAME_Y;

  // align dest width to multiple of 32
  g_w1_frame = (g_w1_frame + 31) & -32;

  // readjust source width accordingly
  g_w0 = g_w1_frame << DOWNSAMPLE_SHIFT_FRAME_X;

  // z buffer dest dimensions
  g_w1_depth = g_w0 >> DOWNSAMPLE_SHIFT_DEPTH_X;
  g_h1_depth = g_h0 >> DOWNSAMPLE_SHIFT_DEPTH_Y;



  // try to alloc vram buffers; quit if alloc fails
  u32 vram_state = g_ScratchVram.GetState();
  u32 fb_mem, zb_mem, blur_mem, blur_mem2;
  if (!R2O_AllocVram(fb_mem, zb_mem, blur_mem, blur_mem2))
  {
    g_ScratchVram.Reset(vram_state);
    return;
  }

  // downsample the frame & depth buffers
  IGTexture fb_texture, zb_texture, blur_texture, blur_texture2;
  R2O_DownsampleBuffers(&fb_texture, fb_mem, &zb_texture, zb_mem);

  // generate blurred copies of downsampled framebuffer
  R2O_BlurTexture(&blur_texture,  blur_mem,  &fb_texture,  128, 64);
  R2O_BlurTexture(&blur_texture2, blur_mem2, &blur_texture, 64, 32);

  u32  scratch_alloc_state = g_ScratchAlloc.GetState();
  u8 num_blendmaps = 0;
  u8* waterobj_to_blendmap = g_ScratchAlloc.Alloc(g_R2OCon.m_num_water_objects*sizeof(u8));
  u8* blendmap_to_cubemap = g_ScratchAlloc.Alloc(g_R2OCon.m_num_water_objects*sizeof(u8));
  IGTexture** blendmap_pointers = (IGTexture**)g_ScratchAlloc.Alloc(g_R2OCon.m_num_water_objects*sizeof(IGTexture*));

  // loop over water objects doing something or other
  R2OWaterObject* p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++, p_obj++)
  {
    if (!(p_obj->m_flags & R2O_WATER_OBJECT_FLAG_ACTIVE))
    {
      continue;
    }

    if (!p_obj->m_render_data[vp_idx].m_visible)
    {
      continue;
    }

    u8 blendmap_index = 0xFF;
    for(u32 j=0; j<num_blendmaps; j++)
    {
      if ( blendmap_to_cubemap[j] == p_obj->m_cubemap_index )
      {
        blendmap_index = j;
        break;
      }
    }

    if ( blendmap_index == 0xFF )
    {
      blendmap_index = num_blendmaps++;
      blendmap_to_cubemap[blendmap_index] = p_obj->m_cubemap_index;
      blendmap_pointers[blendmap_index] = (IGTexture*)g_ScratchAlloc.Alloc(sizeof(IGTexture), 64);
    }

    waterobj_to_blendmap[i_obj] = blendmap_index;
  }

  R2O_BlendCubeMaps(blendmap_pointers, blendmap_to_cubemap, num_blendmaps, &blur_texture2, p_view_data);

  // set the global textures

  // the texture unit usage is as follows:
  // 0 : coarse normal map
  // 1 : fine normal map
  // 2 : downsampled depth buffer
  // 3 : downsampled frame buffer
  // 4 : reflection texture
  // 5 : fresnel reflectance tex
  // 6 : interactive index map 0
  // 7 : interactive index map 1
  // 8 : interactive index map 2
  // 10: interactive normal maps 0
  // 11: interactive normal maps 1
  // 12: interactive normal maps 2

  if (g_R2OCon.m_any_transparent_visible)
  {
  g_CurrentPB->IGPushBuffer::SetTexture(2, &zb_texture);
  }
  g_CurrentPB->IGPushBuffer::SetTexture(3, &fb_texture);


  g_CurrentPB->IGPushBuffer::SetTexture(5, &g_R2OFresnelTexture);

  g_CurrentPB->SetVertexFormat(VERTEX_DESC_R2O);

  g_CurrentPB->SetRenderState(RENDERSTATE_BLEND, false);
  g_CurrentPB->SetRenderState(RENDERSTATE_DEPTHWRITE, true);
  g_CurrentPB->SetFrontFace(p_view_data->m_camera_underwater ? FRONTFACE_CCW : FRONTFACE_CW);

  g_CurrentPB->SetPrimitiveRestartIndex(0xFFFF);
  g_CurrentPB->SetRenderState(RENDERSTATE_PRIMITIVERESTART,true);



  float near = p_view_data->m_frustum_near;
  float far  = p_view_data->m_frustum_far;
  float k    = (near-far) / (near*far);

  // magic z-buffer constants
  f32 table[4];
  #if READ_DEPTH_AS_G16R16
  table[0] = 0.0f;
  table[1] = k * (65536.0f / 65537.0f);
  table[2] = k * (    1.0f / 65537.0f);
  #else // A8R8G8B8
  table[0] = k * (65536.0f/65793.0f);
  table[1] = k * (  256.0f/65793.0f);
  table[2] = k * (    1.0f/65793.0f);
  #endif
  table[3] = 1.0f / near;
  float zbuffer_consts[4] ALIGNED(16);
  u32 swizzle = zb_texture.m_swizzle;   // 0xE4  11100100
  u32 w_cmpnt = (swizzle >> 0) & 3;     // 0 
  u32 x_cmpnt = (swizzle >> 2) & 3;     // 1
  u32 y_cmpnt = (swizzle >> 4) & 3;     // 2
  u32 z_cmpnt = (swizzle >> 6) & 3;     // 3
  zbuffer_consts[0] = table[x_cmpnt];   // k * (65536.0f / 65537.0f);
  zbuffer_consts[1] = table[y_cmpnt];   // k * (    1.0f / 65537.0f);
  zbuffer_consts[2] = table[z_cmpnt];   // 1.0f / near;
  zbuffer_consts[3] = table[w_cmpnt];   // 0.0f

  // magic framebuffer constants
  float framebuffer_consts[8];
  framebuffer_consts[0] =  1.0f / g_w0;
  framebuffer_consts[1] =  1.0f / g_h0;
  framebuffer_consts[2] =  0.0f;
  framebuffer_consts[3] =  0.0f;
  framebuffer_consts[4] =  0.0f;
  framebuffer_consts[5] =  0.0f;
  framebuffer_consts[6] =  0.0f;
  framebuffer_consts[7] =  1.0f;

  // set constants for super-near clipping
  f32 g = 1.01f;  // special depth band is near to g*near
  f32 h = 0.005;  // super-near plane is h*near

  f32 near_consts[3];
  near_consts[0] = near * g;
  near_consts[1] = (g-1.0f) / (g-h);
  near_consts[2] = near * g * (1.0f-h) / (g-h);
  g_CurrentPB->SetVertexProgramConstant(VP_CONST_R2O_NEAR_CONSTANTS, (float *)near_consts);

  f32 w2c2[4];
  w2c2[0] = p_view_data->m_world_to_camera_matrix[0].z;
  w2c2[1] = p_view_data->m_world_to_camera_matrix[1].z;
  w2c2[2] = p_view_data->m_world_to_camera_matrix[2].z;
  w2c2[3] = p_view_data->m_world_to_camera_matrix[3].z;
  g_CurrentPB->SetVertexProgramConstant4(VP_CONST_R2O_W2C_COL_2, (float *)w2c2);


  // scale texproj matrix horizontally according to viewport
  mtx4f view_mat = p_view_data->m_world_to_camera_matrix;
  mtx4f proj_mat = p_view_data->m_frag_texproj_matrix;
  mtx4f texproj;
  view_mat[3].SetXYZW(0,0,0,1);
  float *p = (float *)&proj_mat;
  p[0]  = ( p[0]  + p[ 3]) * 0.5f;
  p[8]  = ( p[8]  + p[11]) * 0.5f;
  p[5]  = (-p[5]  + p[ 7]) * 0.5f;
  p[9]  = (-p[9]  + p[11]) * 0.5f;
  MatMulMat(texproj, view_mat, proj_mat);


  // settings common to all r2o fragment programs, for all lods of all water objects
  for (u32 fp=FRAGMENT_PROG_R2O_0; fp<=FRAGMENT_PROG_R2O_OPAQUE_7; fp++)
  {
    IGFragmentProgram* p_fragment_prog = g_FragmentPrograms[fp];

    // z buffer consts
    p_fragment_prog->SetConstant(FP_CONST_R2O_ZBUFFER_CONSTANTS, (float *)zbuffer_consts);

    // frame buffer consts
    p_fragment_prog->SetConstant(FP_CONST_R2O_FB_SCALE,          (float *)framebuffer_consts);
    p_fragment_prog->SetConstant(FP_CONST_R2O_FB_OFFSET,         (float *)framebuffer_consts + 4);

    // view projection matrix
    p_fragment_prog->SetConstant(FP_CONST_R2O_TEX_PROJ_MAT_0, (float *)&texproj + 0);
    p_fragment_prog->SetConstant(FP_CONST_R2O_TEX_PROJ_MAT_1, (float *)&texproj + 4);
    p_fragment_prog->SetConstant(FP_CONST_R2O_TEX_PROJ_MAT_2, (float *)&texproj + 8);
    p_fragment_prog->SetConstant(FP_CONST_R2O_TEX_PROJ_MAT_3, (float *)&texproj + 12);

    // index of refraction
    p_fragment_prog->SetConstant(FP_CONST_R2O_MISC_CONSTS, (float *)g_R2OCon.m_misc_consts);
  }




  IGTexture* tex_ptr = NULL;
  u32 total_verts   = 0;
  u32 total_indices = 0;

  vec4f zero;
  zero.SetXYZW(0,0,0,0);

  // loop over water objects
  p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++, p_obj++)
  {
    if (!(p_obj->m_flags & R2O_WATER_OBJECT_FLAG_ACTIVE))
    {
      continue;
    }

    if (AvailableDynamicRenderMem() < PUSH_BUFFER_PADDING)
      goto OutOfPushBufferSpace;

    if (!p_obj->m_render_data[vp_idx].m_visible)
    {
      continue;
    }

    // aux pointers
    R2ORenderData      *p_rdata = &p_obj->m_render_data[p_view_data->m_vp_index];
    R2OInteractiveData *p_idata = &p_obj->m_interactive_data;

    IGTexture *b_tex = blendmap_pointers[waterobj_to_blendmap[i_obj]];
    if (tex_ptr != b_tex)
    {
      tex_ptr = b_tex;
      g_CurrentPB->InvalidateTextureCachedState(4);
      g_CurrentPB->IGPushBuffer::SetTexture(4, b_tex);
      SetCubeMapTexture(g_CurrentPB, &g_ShaderDatabase.m_cube_maps[p_obj->m_cubemap_index], CUBE_SPEC_VERY_TIGHT);
    }

    // set for each water object in case the previous one set VERTEX_PROG_R2O_NEAR
    g_CurrentPB->SetVertexProgram(g_VertexPrograms[VERTEX_PROG_R2O]);

    // derive deriv scales based on amplitude
    f32 ratio = 0.2f / p_obj->m_amplitude;
    f32 vert_deriv_scale = 32768.0f * ratio;
    f32 height_scale     = g_R2OCon.m_height_scale * ratio;
    f32 map_deriv_scale  = 128.0f * ratio;
    f32 vert_deriv_invscale = 1.0f / vert_deriv_scale;
    f32 map_deriv_invscale  = 1.0f / map_deriv_scale;
    f32 height_invscale     = 1.0f / height_scale;
    f32 inv_scales[2];
    inv_scales[0] = vert_deriv_invscale;
    inv_scales[1] = height_invscale;
    g_CurrentPB->SetVertexProgramConstant(VP_CONST_R2O_INV_SCALES, (float *)inv_scales);

    // scale the normal map weights accordingly
    f32 map_weights[3];
    map_weights[0] = g_R2OCon.m_coarse_map_weight * 128.0f * map_deriv_invscale;
    map_weights[1] = g_R2OCon.m_fine_map_weight   * 128.0f * map_deriv_invscale;
    map_weights[2] = g_R2OCon.m_int_map_weight;


    // settings common to all r2o fragment programs, for all lods of the current water object
    for (u32 fp=FRAGMENT_PROG_R2O_0; fp<=FRAGMENT_PROG_R2O_OPAQUE_7; fp++)
    {
      IGFragmentProgram* p_fragment_prog = g_FragmentPrograms[fp];

      // normal map weights (dependent on amplitude)
      p_fragment_prog->SetConstant(FP_CONST_R2O_NORMAL_MAP_WEIGHTS, (float *)map_weights);

      // colour
      p_fragment_prog->SetConstant(FP_CONST_R2O_WATER_COLOR, p_view_data->m_camera_underwater ?
                                   (float *)&zero : (float *)&p_obj->m_water_color);
    }



    f32 step_even = g_R2OCon.m_step;

    // loop over lods submitting a set of verts for each (if nonempty)
    for (u32 lod=0; lod<g_R2OCon.m_num_lods; lod++)
    {
      // activate special vertex program to render stuff nearer than the near plane
      if (lod == g_R2OCon.m_near_lod)
      {
        g_CurrentPB->SetVertexProgram(g_VertexPrograms[VERTEX_PROG_R2O_NEAR]);
      }

      vec4f basis0, basis1, alt_basis0, alt_basis1;
      f32 step = (lod&1) ? step_even*0.707106781187f : step_even;
      SetBasis(basis0, basis1, alt_basis0, alt_basis1, lod);

      vec4f dvc, dvr;
      if (lod & 1)
      {
        f32 s = 0.5f * step_even;
        dvc.SetXYZ( s,0,s);
        dvr.SetXYZ(-s,0,s);
      }
      else
      {
        f32 s = step_even;
        dvc.SetXYZ(s,0,0);
        dvr.SetXYZ(0,0,s);
      }

      u32 num_verts   = p_rdata->m_num_verts[lod];
      u32 num_indices = p_rdata->m_num_indices[lod];

      if (num_verts && num_indices)
      {
        vec4f org      = p_rdata->m_origins[lod];

        // set normal map matrix
        vec4f norm_map_mat0, norm_map_mat1;

        if (!(lod & 1))
        {
          f32 coarse_scale =  4.0f / step;          // 4 semi-octaves above geometry lod
          f32 fine_scale   = 16.0f / step;          // 8 semi-octaves above geometry lod

          norm_map_mat0.x = basis0.x * coarse_scale;
          norm_map_mat0.y = basis1.x * coarse_scale;
          norm_map_mat0.z = basis0.x * fine_scale;
          norm_map_mat0.w = basis1.x * fine_scale;
          norm_map_mat1.x = basis0.z * coarse_scale;
          norm_map_mat1.y = basis1.z * coarse_scale;
          norm_map_mat1.z = basis0.z * fine_scale;
          norm_map_mat1.w = basis1.z * fine_scale;
        }
        else
        {
          f32 coarse_scale = Sqrtf(  8.0f) / step;   // 3 semi-octaves above geometry lod
          f32 fine_scale   = Sqrtf(128.0f) / step;   // 7 semi-octaves above geometry lod

          norm_map_mat0.x = alt_basis0.x * coarse_scale;
          norm_map_mat0.y = alt_basis1.x * coarse_scale;
          norm_map_mat0.z = alt_basis0.x * fine_scale;
          norm_map_mat0.w = alt_basis1.x * fine_scale;
          norm_map_mat1.x = alt_basis0.z * coarse_scale;
          norm_map_mat1.y = alt_basis1.z * coarse_scale;
          norm_map_mat1.z = alt_basis0.z * fine_scale;
          norm_map_mat1.w = alt_basis1.z * fine_scale;
        }

        // set vertex prog consts
        g_CurrentPB->SetVertexProgramConstant(VP_CONST_R2O_DVC_WORLD,     (float *)&dvc);
        g_CurrentPB->SetVertexProgramConstant(VP_CONST_R2O_DVR_WORLD,     (float *)&dvr);
        g_CurrentPB->SetVertexProgramConstant(VP_CONST_R2O_ORG_WORLD,     (float *)&org);
        g_CurrentPB->SetVertexProgramConstant(VP_CONST_R2O_NORM_MAP_MAT0, (float *)&norm_map_mat0);
        g_CurrentPB->SetVertexProgramConstant(VP_CONST_R2O_NORM_MAP_MAT1, (float *)&norm_map_mat1);

        //------------------------------------------------------------------------------------------
        // interactive water stuff

        u32 fp_mask = 0x0;

        f32 interactive_constants[4];
        f32 offset = 8.0f * step_even;
        f32 scale  = 0.5f / step_even;

        // loop over interactive layers for this lod
        for (u32 i=0,bit=4; i<3; i++,bit>>=1)
        {
          // set vertex prog consts
          vec4f *p_org = &p_rdata->m_origins[(lod-2*i) & -2];
          interactive_constants[0] = p_org->x + offset;
          interactive_constants[1] = 0.0f;
          interactive_constants[2] = p_org->z + offset;
          interactive_constants[3] = scale;
          g_CurrentPB->SetVertexProgramConstant(VP_CONST_R2O_INTERACTIVE_CONSTANTS0 + i, (float *)interactive_constants);

          // determine whether to set index & normal maps to anything non-zero
          u32 ea_idx=0;
          i32 l = (lod-2*i) >> 1;
          if (l>=0 && l<MAX_INT_LODS-3)
          {
            ea_idx  = p_rdata->m_ea_index_map[l];
          }

          if (ea_idx)
          {
            // set index map texture
            u16 cols_rows = p_rdata->m_cols_rows[(lod-2*i) & -2];
            u32 cols = ((u32)cols_rows >>  8 ) - 16;
            u32 rows = ((u32)cols_rows & 0xFF) - 16;
            cols >>= 1;
            rows >>= 1;
            R2O_SetInteractiveIndexMapTexture(6+i, ea_idx, cols, rows);

            // set normal map texture
            u32 ea_norm = p_idata->m_ea_normal_maps[l+3];
            R2O_SetInteractiveNormalMapTexture(10+i, ea_norm);

            fp_mask |= bit;
          }
          else
          {
            // both maps zero
            R2O_SetInteractiveIndexMapTexture(6+i, 0, 0, 0);
            R2O_SetInteractiveNormalMapTexture(10+i, 0);
          }

          // ready for next layer
          offset *= 2.0f;
          scale  *= 0.5f;
        }

        u32 fp_base = (p_obj->m_water_color.w >= 100.0f) ? FRAGMENT_PROG_R2O_OPAQUE_0 : FRAGMENT_PROG_R2O_0;
        u32 fp = fp_base + fp_mask;
        g_CurrentPB->SetFragmentProgram(g_FragmentPrograms[fp]);
        g_CurrentPB->RefreshFragmentProgram(g_FragmentPrograms[fp]);

        //------------------------------------------------------------------------------------------

        // set ambient normal map textures
        u32 coarse_map = (lod & 1) ? lod-1 : lod;
        u32 fine_map   = (lod & 1) ? lod+3 : lod+4;
        coarse_map = (coarse_map < 31) ? coarse_map : 31;
        fine_map   = (fine_map   < 31) ? fine_map   : 31;
        R2O_SetNormalMapTexture(0, g_NormalMapPointers[coarse_map]);
        R2O_SetNormalMapTexture(1, g_NormalMapPointers[fine_map  ]);

        // set the vertex streams
        g_CurrentPB->SetVertexAttribPointer(0, MainMemoryAddressToOffset((void *)p_rdata->m_ea_verts[lod]));
        g_CurrentPB->SetVertexAttribPointer(1, MainMemoryAddressToOffset((void *)p_rdata->m_ea_derivs[lod]));

        // set the index buffer
        g_CurrentPB->SetIndexBase(MainMemoryAddressToOffset((void *)p_rdata->m_ea_indices[lod]));

        // add the mesh
        IGPrimitive prim_type = g_R2OCon.m_wireframe ? PRIM_LINESTRIP : PRIM_TRIANGLEFAN;
        g_CurrentPB->DrawBegin(prim_type);
        g_CurrentPB->DrawIndexedPrim(0, num_indices);
        g_CurrentPB->DrawEnd();

        total_verts += num_verts;
        total_indices += num_indices;
      }

      if (lod & 1)
      {
        step_even *= 0.5f;
      }
    }
  }

OutOfPushBufferSpace:

  //Reset the PU scratch allocator
  g_ScratchAlloc.Reset(scratch_alloc_state);
  g_ScratchVram.Reset(vram_state);

  // reset any state we changed that is not reset by other render layers
  g_CurrentPB->SetRenderState(RENDERSTATE_PRIMITIVERESTART, false);
  g_CurrentPB->SetFrontFace(FRONTFACE_CCW);
  g_CurrentPB->DisableTexture(SAMPLER_CUBE_MAP);

  // invalidate the cached state of any texture units which were set using SetTexture
  g_CurrentPB->InvalidateTextureCachedState(0);
  g_CurrentPB->InvalidateTextureCachedState(1);
  g_CurrentPB->InvalidateTextureCachedState(2);
  g_CurrentPB->InvalidateTextureCachedState(3);
  g_CurrentPB->InvalidateTextureCachedState(4);
  g_CurrentPB->InvalidateTextureCachedState(5);
  g_CurrentPB->InvalidateTextureCachedState(6);
  g_CurrentPB->InvalidateTextureCachedState(7);
  g_CurrentPB->InvalidateTextureCachedState(8);
  g_CurrentPB->InvalidateTextureCachedState(10);
  g_CurrentPB->InvalidateTextureCachedState(11);
  g_CurrentPB->InvalidateTextureCachedState(12);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_RenderShell(BucketCon* bucketer, i32 vp_index)
{
  if (g_R2OCon.m_r2o_enable && g_R2OCon.m_num_water_objects && g_R2OCon.m_view_data[vp_index].m_any_visible)
  {
    bucketer->m_layers[ SORT_LAYER_WATER ].m_layer_func = &ProcessLayerR2O;
    IGG::BucketerAddEntry(bucketer, (u8*)&g_R2OCon, (u8 *)&g_R2OCon.m_view_data[vp_index], 0, SORT_LAYER_WATER, 0);
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

u32 R2O_TicksLastFrame()
{
  return g_R2OCon.m_spu_ticks;
}

} // namespace IGG

