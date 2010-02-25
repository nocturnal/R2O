#include "igg/igGraphics.h"
#include "igR2O.h"
#include "igR2ODebug.h"
#include "igImpulse/igImpulse.h"
#include "igEffectsConduit/api.h"


namespace IGG
{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// free a list of tiles back to the free list
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_FreeTileList(u16 head_idx, u16 cnt)
{
  if (!head_idx)
  {
    return;
  }

  Tile *tiles = (Tile *)g_R2OCon.m_int_ea_tiles;
  u16 tile_idx = head_idx;

  // find tail of list
  Tile *p_tile = 0;
  while (tile_idx)
  {
    p_tile = &tiles[tile_idx];
    tile_idx = p_tile->m_next_tile;
  }

  // append to front of free list
  p_tile->m_next_tile = g_R2OCon.m_int_free_tile;
  g_R2OCon.m_int_free_tile = head_idx;

  // reduce in-use count
  g_R2OCon.m_int_num_tiles -= cnt;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_InitInteractiveTiles(u32 max_tiles)
{
  u32 ideal_num_tiles = 3*max_tiles / 4;

  // allocate memory
  Tile *tiles = (Tile *)g_SystemMemoryAlloc.Alloc(max_tiles * sizeof(Tile), 128);
  i16  *maps  = (i16  *)g_SystemMemoryAlloc.Alloc(max_tiles * 256*2*sizeof(i16), 128);

  // init zero tile
  tiles[0].m_coords.m_u64 = 0x7F7FFFFF7F7FFFFFULL;
  tiles[0].m_next_tile = 0;
  tiles[0].m_cache_slot_upper = 0;
  tiles[0].m_cache_slot_lower = 0;
  tiles[0].m_flags = 0x01;  // ensure it persists

  // link all other tiles into the free list
  for (u32 i=1; i<max_tiles-1; i++)
  {
    tiles[i].m_next_tile = i+1;
  }
  tiles[max_tiles-1].m_next_tile = 0;

  // clear tile list head and tile count for each lod
  R2OWaterObject *p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (u32 i_obj=0; i_obj<g_R2OCon.m_max_water_objects; i_obj++, p_obj++)
  {
    for (u32 i=0; i<MAX_INT_LODS; i++)
    {
      p_obj->m_interactive_data.m_tile_list[i] = 0;
      p_obj->m_interactive_data.m_num_tiles[i] = 0;
    }
  }

  // clear all height maps
  for (u32 i=0; i<max_tiles*256*2; i++)
  {
    maps[i] = 0;
  }

  // set R2OCon vars
  g_R2OCon.m_int_ea_tiles        = (u32)tiles;
  g_R2OCon.m_int_ea_maps         = (u32)maps;
  g_R2OCon.m_int_max_tiles       = max_tiles;
  g_R2OCon.m_int_ideal_num_tiles = ideal_num_tiles;
  g_R2OCon.m_int_num_tiles       = 1;
  g_R2OCon.m_int_free_tile       = 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_InitInteractive( u32 interactive_tile_count )
{
  // set system constants
  u32 map_width  = 16;
  u32 ovlp_width = 16;
  u32 fft_width  = map_width + ovlp_width;
  u32 num_bands  = 32;
  u32 fft_area   = fft_width * fft_width;
  f32 fft_meters = (g_R2OCon.m_step != 0.0f) ? g_R2OCon.m_step * (f32)fft_width : 1.0f;

  u16 wavenumbers_squared[256];
  u32 pal_cnt  = R2O_GenerateWaveMap(wavenumbers_squared, 256, fft_width);
  u32 pal_cnt4 = ALIGN_4(pal_cnt);

  // allocs
  u8   *index_map         = (u8   *)g_SystemMemoryAlloc.Alloc(fft_area  * sizeof(u8), 128);
  f32  *frequency_palette = (f32  *)g_SystemMemoryAlloc.Alloc(pal_cnt4  * num_bands * sizeof(f32), 128);

  // generate index map & frequency palette
  R2O_GenerateIndexMap(index_map, wavenumbers_squared, fft_width, pal_cnt4);
  R2O_GenerateFrequencyPalette(frequency_palette, wavenumbers_squared, fft_meters, num_bands, pal_cnt4);

  // interactive tile budget
  R2O_InitInteractiveTiles( interactive_tile_count );

  // store constants to g_R2OCon
  g_R2OCon.m_int_map_width            = map_width;
  g_R2OCon.m_int_ovlp_width           = ovlp_width;
  g_R2OCon.m_int_fft_width            = fft_width;
  g_R2OCon.m_int_num_bands            = num_bands;
  g_R2OCon.m_int_palette_cnt          = pal_cnt4;
  g_R2OCon.m_int_half_life            = 5.0f;
  g_R2OCon.m_int_enable               = true;
  g_R2OCon.m_int_asm_line             = 10000;
  g_R2OCon.m_int_map_weight           = 0.2f;
  g_R2OCon.m_int_scale                = 5.0f;
  g_R2OCon.m_int_threshold_create     = 0.002;
  g_R2OCon.m_int_threshold_destroy    = 0.004;
  g_R2OCon.m_int_cull_distance        = 200.0f;
  g_R2OCon.m_int_cull_dist_ratio      = 0.5f;
  g_R2OCon.m_int_inner_radius         = 0.5f;
  g_R2OCon.m_int_inner_displacement   = -0.2f;
  g_R2OCon.m_int_outer_displacement   = -0.5f;
  g_R2OCon.m_int_attenuation_radius   = 1.0f;
  g_R2OCon.m_int_attenuation_min      = 0.25f;
  g_R2OCon.m_int_damping_coeff0       = 0.1f;
  g_R2OCon.m_int_damping_coeff1       = 0.5f;
  g_R2OCon.m_int_cull_dist_coeff0     = 1.1f;
  g_R2OCon.m_int_cull_dist_coeff1     = 0.99f;
  #if SANITY_CHECKS
  g_R2OCon.m_int_error                = 0;
  #endif


  // store ea's to g_R2OCon
  g_R2OCon.m_int_ea_frequency_palette = (u32)frequency_palette;
  g_R2OCon.m_int_ea_index_map         = (u32)index_map;
  g_R2OCon.m_int_ea_impulses          = (u32)IMPULSE::GetImpulseBuffer()->m_list;


  // get a consumer mask from igImpulse
#if (BUILD < PUBLISH)
  g_R2OCon.m_int_consumer_mask = IMPULSE::RequestConsumerMask("R2O");
#else
  g_R2OCon.m_int_consumer_mask = IMPULSE::RequestConsumerMask();
#endif
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// free up the interactive resources of a water object being streamed out
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_DeleteInteractive(R2OInteractiveData *p_int_data)
{
  for (u32 l=0; l<MAX_INT_LODS; l++)
  {
    u16 head_idx = p_int_data->m_tile_list[l];
    if (head_idx)
    {
      u16 cnt = p_int_data->m_num_tiles[l];
      R2O_FreeTileList(head_idx, cnt);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


IMPULSE::ImpulseData *g_R2OImpulseList;


static int CmpImpulse(const void *p0, const void *p1)
{
  u16 idx0 = *(u16 *)p0;
  u16 idx1 = *(u16 *)p1;

  f32 rad0 = g_R2OImpulseList[idx0].m_volume.w;
  f32 rad1 = g_R2OImpulseList[idx1].m_volume.w;
  
  return rad0<rad1 ? 1 : rad0>rad1 ? -1 : 0;
}


void R2O_AssignImpulsesToWaterObjects()
{
  // clear impulse counts
  R2OWaterObject *p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (u32 i_obj=0; i_obj<g_R2OCon.m_max_water_objects; i_obj++, p_obj++)
  {
    p_obj->m_interactive_data.m_num_impulses = 0;
    p_obj->m_interactive_data.m_impulse_lods = 0;
  }

  // get the impulse buffer and impulse list
  IMPULSE::ImpulseBuffer *impulse_buffer = IMPULSE::GetImpulseBuffer();
  g_R2OImpulseList = impulse_buffer->m_list;

  // loop over impulses
  for (u32 i_imp=0; i_imp<impulse_buffer->m_max_index; i_imp++)
  {
    IMPULSE::ImpulseData *p_data = &g_R2OImpulseList[i_imp];

    // skip this impulse if it's already been consumed by R2O
    if (p_data->m_consumer_mask & g_R2OCon.m_int_consumer_mask)
    {
      continue;
    }

    // set consumer mask
    p_data->m_consumer_mask |= g_R2OCon.m_int_consumer_mask;

    // get start & end spheres
    vec4f pos0, pos1, unit_vec;
    pos0 = p_data->m_volume;
    unit_vec = p_data->m_impulse;
    f32 length = unit_vec.w;
    VecMulAdd3(pos1, unit_vec, length, pos0);

    // generate bounding box of swept sphere
    vec4f imp_min, imp_max;
    f32 radius = pos0.w;

    VecMin3(imp_min, pos0, pos1);
    VecMax3(imp_max, pos0, pos1);
    ScalarSubVector3(imp_min, imp_min, radius);
    ScalarAddVector3(imp_max, imp_max, radius);

    // loop over water objects, performing bounding box rejection test
    R2OWaterObject *p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
    for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++, p_obj++)
    {
      if (!(p_obj->m_flags & R2O_WATER_OBJECT_FLAG_ACTIVE))
      {
        continue;
      }

      // get object extents
      vec4f obj_min, obj_max;
      obj_min = p_obj->m_origin;
      AddVector3(obj_max, obj_min, p_obj->m_dimensions);

      // test AABBs for overlap
      if (!(VecLessAny3(obj_max, imp_min) || VecLessAny3(imp_max, obj_min)))
      {
        //printf("registering impulse %d to water object %d\n", i_imp, i_obj);

        // overlap: register the given impulse to this water object (provided we haven't filled its buffer)
        u32 idx = p_obj->m_interactive_data.m_num_impulses;
        if (idx < MAX_IMPULSES)
        {
          p_obj->m_interactive_data.m_impulses[idx] = i_imp;
          p_obj->m_interactive_data.m_num_impulses  = idx+1;

          // record which lods it can potentially poke into
          f32 ratio = INT_STEP0 / radius;
          u32 min_lod = QUICK_LOG2(1.5f * ratio);
          p_obj->m_interactive_data.m_impulse_lods |= (0x7 << min_lod); // spans upto 2 lods (crossing upto 3)
        }
      }
    }
  }

  // for each water object, sort the impulse indices in order of decreasing radius (i.e. low to high lod)
  p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++, p_obj++)
  {
    if (!(p_obj->m_flags & R2O_WATER_OBJECT_FLAG_ACTIVE))
    {
      continue;
    }

    u16 *indices = p_obj->m_interactive_data.m_impulses;
    u32 cnt      = p_obj->m_interactive_data.m_num_impulses;
    IGG::Qsort(indices, cnt, sizeof(u16), CmpImpulse);
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// a rather kludgy way of passing back the spu-generated intersection information for triggering splash effects
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_TriggerEventsFromImpulses()
{
  // temphack - this might get called while the gameplay tier is torn down during a region transition so abort if so
  if (IGG::g_MobyInsts.m_array == NULL)
    return;

  // get the impulse buffer and impulse list
  IMPULSE::ImpulseBuffer *impulse_buffer = IMPULSE::GetImpulseBuffer();
  IMPULSE::ImpulseData   *impulse_list   = impulse_buffer->m_list;

  // loop over impulses
  for (u32 i_imp=0; i_imp<impulse_buffer->m_max_index; i_imp++)
  {
    IMPULSE::ImpulseData *p_data = &impulse_list[i_imp];

    // check whether to make a splash
    bool b_splash = (p_data->m_flags & IMPULSE::IMPULSE_FLAG_SPLASH);

    // splash (provided there's a callback function)
    if  (g_R2OInteractionCallback && b_splash)
    {
      g_R2OInteractionCallback(p_data);                   // notify the game of the interaction
      p_data->m_flags &= ~IMPULSE::IMPULSE_FLAG_SPLASH;   // clear the splash flag so the splash only happens once
    }

    // visualize impulse (with different colours for splash / no splash)
    R2O_VisualizeImpulse(p_data, b_splash);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// set the half-life for interactive waves
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_SetInteractiveHalfLife(f32 seconds)
{
  g_R2OCon.m_int_half_life = seconds;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Texture settings for interactive normal maps
//
// 3D, 2 bytes per texel, normalized texture coords
// dimensions 24 x 18 x (>=num_tiles)
// wrap mode should be irrelevant because no reads will occur outside the valid portion of each map
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_SetInteractiveNormalMapTexture(u32 tex_unit, u32 ea)
{
  if (!ea)
  {
    g_CurrentPB->IGPushBuffer::SetTexture(tex_unit, &g_BlackTexture);
    return;
  }

  const u32 width  = 24;
  const u32 height = 18;
  const u32 depth  = 512;
  const u32 stride = 2*24;
  //const u32 stride = (width+2) * (height+2) * 2;
  IGTexture tex;

  tex.m_baseOffset = AddressToOffset((void *)ea) & 0x7FFFFFFF;
  tex.m_format     = TEXTURE_FORMAT(1, 1, 1, TEXTUREFORMAT_G8B8, 3, 1, 0, 1, 0);
  tex.m_control1   = TEXTURE_CONTROL1(0, 0, TEXTUREADDRESS_BORDER, 0, TEXTUREADDRESS_BORDER, 0, TEXTUREADDRESS_BORDER);
  tex.m_control2   = TEXTURE_CONTROL2(1, 0, 0, 0, 0, 0);
  tex.m_swizzle    = TEXTURE_SWIZZLE(0, 0xAA, 0x1B);
  tex.m_filter     = TEXTURE_FILTER(0xF, TEXTUREFILTER_MAG_LINEAR, TEXTUREFILTER_MIN_LINEAR, 1, 0x1C00);
  tex.m_size1      = TEXTURE_SIZE1(width, height);
  tex.m_size2      = TEXTURE_SIZE2(depth, stride);

  g_CurrentPB->IGPushBuffer::SetTexture(tex_unit, &tex);

  g_CurrentPB->SetTextureBorderColor(tex_unit, 0x80808080);
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Texture settings for interactive index map
//
// 2D, 1 byte per texel, unnormalized texture coords
// dimensions 16 x 16 (temp)
// addresses the border colour for reads outside the index map
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_SetInteractiveIndexMapTexture(u32 tex_unit, u32 ea, u32 width, u32 height)
{
  if (!ea)
  {
    g_CurrentPB->IGPushBuffer::SetTexture(tex_unit, &g_BlackTexture);
    return;
  }

  const u32 stride = width;
  IGTexture tex;
  
  tex.m_baseOffset = AddressToOffset((void *)ea) & 0x7FFFFFFF;
  tex.m_format     = TEXTURE_FORMAT(1, 1, 1, TEXTUREFORMAT_B8, 2, 1, 0, 1, 0);
  tex.m_control1   = TEXTURE_CONTROL1(0, 0, TEXTUREADDRESS_BORDER, 0, TEXTUREADDRESS_BORDER, 0, TEXTUREADDRESS_BORDER);
  tex.m_control2   = TEXTURE_CONTROL2(1, 0, 0, 0, 0, 0);
  tex.m_swizzle    = TEXTURE_SWIZZLE(0, 0xAA, 0x1B);
  tex.m_filter     = TEXTURE_FILTER(0x0, TEXTUREFILTER_MAG_NEAREST, TEXTUREFILTER_MIN_NEAREST, 1, 0x1C00);
  tex.m_size1      = TEXTURE_SIZE1(width, height);
  tex.m_size2      = TEXTURE_SIZE2(1,     stride);
  
  g_CurrentPB->IGPushBuffer::SetTexture(tex_unit, &tex);

  g_CurrentPB->SetTextureBorderColor(tex_unit, 0x00000000);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Terry's callback
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

R2OInteractionCallback g_R2OInteractionCallback = NULL;

void R2O_SetInteractionCallback(R2OInteractionCallback function)
{
  g_R2OInteractionCallback = function;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// now just stubs!
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if 1
void R2O_SetInteractiveGrid(vec4f_arg origin, u32 lod)
{
}

void R2O_GetInteractiveGrid(vf32 &origin, u32 &lod)
{
}

void R2O_GetInteractiveGrid(vec4f &min, vec4f& max)
{
}
#endif

} // namespace IGG

