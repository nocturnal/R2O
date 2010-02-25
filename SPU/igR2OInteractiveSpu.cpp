#include "igR2OSpu.h"
#include "igImpulse/igImpulseStructs.h"

// R2O Interactive Wave System
// ---------------------------


// globals
Tile *g_Tiles;
i16  *g_MapCache;
u8    g_MapCacheMask[(NUM_MAP_CACHE_SLOTS+7) >> 3];
u32   g_NormIdx;
u32   g_NumTiles;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// alloc a tile from the free list (return NULL if none free)
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



Tile *R2O_AllocTile()
{
  // get the index of the first free tile
  u16 free_idx = g_R2OCon.m_int_free_tile;
  if (!free_idx)
  {
    IG_ASSERT(g_R2OCon.m_int_num_tiles == g_R2OCon.m_int_max_tiles);
    return NULL;
  }

  // get the tile pointer
  Tile *p_tile = &g_Tiles[free_idx];

  // point free list to next free tile
  g_R2OCon.m_int_free_tile = p_tile->m_next_tile;

  // clear the tile
  *(qword *)p_tile = (qword)(0);

  // track number tiles in use
  IG_ASSERT(g_R2OCon.m_int_num_tiles < g_R2OCon.m_int_max_tiles);
  g_R2OCon.m_int_num_tiles++;

  // return pointer to newly allocated tile
  return p_tile;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// free a tile back to the free list
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_FreeTile(Tile *p_tile)
{
  p_tile->m_next_tile = g_R2OCon.m_int_free_tile;
  g_R2OCon.m_int_free_tile = p_tile - g_Tiles;

  // track number tiles in use
  IG_ASSERT(g_R2OCon.m_int_num_tiles > 1);  // must keep 1 for the zero tile
  g_R2OCon.m_int_num_tiles--;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// add a tile with coords (x,z) immediately after *p_prev
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Tile *R2O_AddTile(f32 x, f32 z, Tile *p_prev)
{
  Tile *p_tile = R2O_AllocTile();

  if (p_tile)
  {
    // set coords
    p_tile->m_coords.m_f32[0] = x;
    p_tile->m_coords.m_f32[1] = z;

    // ensure it persists
    p_tile->m_flags = 0x01;

    // link the new tile in after the given prev tile
    p_tile->m_next_tile = p_prev->m_next_tile;
    p_prev->m_next_tile = p_tile - g_Tiles;
  }

  // return tile pointer, or NULL if alloc failed
  return p_tile;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Generate the NxN array of phase-shifts for the given lod and time-step
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_GeneratePhaseShifts(Complex phase_shifts[], u32 lod, f32 dt, f32 damping)
{
  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // get the angular frequency palette for the given lod
  u32 pal_cnt = g_R2OCon.m_int_palette_cnt;
  f32 *frequency_palette = (f32 *)R2O_AllocGetWait(g_R2OCon.m_int_ea_frequency_palette, lod, pal_cnt*sizeof(f32));

  // convert to a palette of complex phase-shifts
  // (scale down by 1/1024 to renormalize values after 2DFFT & inverse transform)
  Complex *phase_shift_palette = (Complex *)g_Scratchpad.Alloc(pal_cnt*sizeof(Complex));
  R2O_GenerateComplexExponentials(phase_shift_palette, frequency_palette, dt, pal_cnt, damping * (1.0f/1024.0f));

  // expand phase-shifts from paletted form to an explicit array
  u32 fft_area = g_R2OCon.m_int_fft_width * g_R2OCon.m_int_fft_width;
  R2O_DecompressComplexArraySeparate((f32 *)phase_shifts, (f32 *)phase_shifts+1024, (f32 *)phase_shift_palette, g_IndexMap, fft_area);

  // restore scratchpad state
  g_Scratchpad.Reset(state);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Zero out a (width x height) subrectangle of the fft buffer
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_ZeroRectangle(Complex *dest, u32 dest_col, u32 dest_row, u32 width, u32 height)
{
  u32 fft_width = g_R2OCon.m_int_fft_width;
  u32 dest_skip = fft_width - width;
  dest += (dest_row * fft_width) + dest_col;

  for (u32 r=0; r<height; r++)
  {
    for (u32 c=0; c<width; c++)
    {
      *dest++ = 0.0f;
    }

    dest += dest_skip;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Read into a (width x height) subrectangle of the fft buffer, from a region of the source map
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_ReadRectangle(Complex *dest, u32 dest_col, u32 dest_row, i16 *src, u32 src_col, u32 src_row, u32 width, u32 height)
{
  dest += (dest_row << 5) + dest_col;
  u32 cnt = width * height;

  src += (src_row << 4) + src_col;
  R2O_ReadRectangleAsm((f32 *)dest, src, cnt);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Write a (width x height) subrectangle of the fft buffer, starting at src, to a map, starting at dest
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_WriteRectangle(i16 *dest, u32 dest_col, u32 dest_row, Complex *src, u32 src_col, u32 src_row, u32 width, u32 height,
                        u32 b_half_map)
{
  u32 map_width  = 16;
  u32 fft_width  = 32;

  u32 src_skip  = fft_width - width;
  u32 dest_skip = map_width - width;

  src  += (src_row  * fft_width) + src_col;
  dest += (dest_row * map_width) + dest_col;

  // imaginary parts are offset halfway through the output map
  // (and the output map size is map_width * map_width * 2, or half of that for a half-map)
  const u32 re = 0;
  const u32 im = (map_width * map_width) >> b_half_map;

  for (u32 r=0; r<height; r++)
  {
    for (u32 c=0; c<width; c++)
    {
      dest[re] = si_to_short(si_cflts(si_from_float(src->re), 15));
      dest[im] = si_to_short(si_cflts(si_from_float(src->im), 15));
      src++;
      dest++;
    }

    src  += src_skip;
    dest += dest_skip;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Write the real parts of a (width x height) subrectangle of the fft buffer, starting at src, to a map, starting at dest
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_WriteRealRectangle(f32 *dest, u32 dest_stride, Complex *src, u32 src_stride, u32 width, u32 height)
{
  u32 src_skip  = src_stride  - width;
  u32 dest_skip = dest_stride - width;

  for (u32 r=0; r<height; r++)
  {
    for (u32 c=0; c<width; c++)
    {
      *dest++ = src++->re;
    }

    src  += src_skip;
    dest += dest_skip;
  }
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// threshold the first-order differences in the x and z directions
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

u32 R2O_ThresholdRectangle(f32 *src, u32 col, u32 row, u32 width, u32 height, f32 threshold)
{
  width  -= 1;  // only use differences fully contained in the rectangle
  height -= 1;

  u32 fft_width = g_R2OCon.m_int_fft_width;
  u32 src_skip  = fft_width - width;

  src += (row  * fft_width) + col;

  u32 b_exceeds = false;

  for (u32 r=0; r<height; r++)
  {
    for (u32 c=0; c<width; c++)
    {
      f32 x0 = src[0];
      f32 x1 = src[1];
      f32 x2 = src[fft_width];
      f32 dx1 = x1-x0;
      f32 dx2 = x2-x0;
      if ((fabsf(dx1) > threshold) || (fabsf(dx2) > threshold))
      {
        b_exceeds = true;
      }

      src++;
    }

    src  += src_skip;
  }

  return b_exceeds;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// return 0 if alloc fails (but note that InitMapCache makes the first alloc, which uses the return value 0 to mean 'slot 0')
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

u8 R2O_AllocMapCacheSlot()
{
  u32 byte_idx;
  u8  bit_mask=0;
  const u32 num_bytes = (NUM_MAP_CACHE_SLOTS+7) >> 3;

  // find the first byte with a clear bit
  for (byte_idx=0; byte_idx<num_bytes; byte_idx++)
  {
    u8 byte = g_MapCacheMask[byte_idx];
    if (byte != 0xFF)
    {
      bit_mask = ~byte & (byte+1);
      break;
    }
  }

  if (bit_mask)
  {
    // get the index of the bit within the byte
    u32 bit_idx=0;
    for (u8 b=0x01; b!=bit_mask; b<<=1)
    {
      bit_idx++;
    }

    // set the bit to indicate the cache slot is now in use
    g_MapCacheMask[byte_idx] |= bit_mask;

    // return the index of the slot
    return (u8)(byte_idx<<3 | bit_idx);
  }
  else
  {
    // return 'cache full'
    return (u8)0;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_FreeMapCacheSlot(u32 slot_idx)
{
  // we must never free slot zero
  if (slot_idx==0)
  {
    return;
  }

  u32 byte_idx = slot_idx >> 3;
  u32 bit_idx  = slot_idx &  7;
  u8  bit_mask = 0x01 << bit_idx;

  // clear the bit
  g_MapCacheMask[byte_idx] &= ~bit_mask;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_ClearMapCacheSlot(u32 slot_idx)
{
  i16 *map = &g_MapCache[slot_idx << 8];
  for (u32 i=0; i<256; i++)
  {
    map[i] = 0;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_InitMapCache()
{
  // alloc mem for cached half maps
  g_MapCache = (i16 *)g_Scratchpad.Alloc(NUM_MAP_CACHE_SLOTS * 256 * sizeof(i16));

  // initialize empty state
  const u32 num_bytes = (NUM_MAP_CACHE_SLOTS+7) >> 3;
  for (u32 i=0; i<num_bytes; i++)
  {
    g_MapCacheMask[i] = 0;
  }

  // fill the zero slot with zeros
  u32 slot = R2O_AllocMapCacheSlot();
  IG_ASSERT(slot==0);
  R2O_ClearMemAsm((qword *)&g_MapCache[slot << 8], 512);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline Tile *NextTile(Tile *p_tile)
{
  IG_ASSERT(p_tile);
  return &g_Tiles[p_tile->m_next_tile];
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Tile *AdvanceTile(Tile *p_tile, f32 x, f32 z)
{
  IG_ASSERT(p_tile);
  while (1)
  {
    Tile *p_next = &g_Tiles[p_tile->m_next_tile];
    IG_ASSERT(p_next);
    f32 tx = p_next->m_coords.m_f32[0];
    f32 tz = p_next->m_coords.m_f32[1];
    if (tz>z || (tz==z && tx>x))
    {
      return p_tile;
    }
    p_tile = p_next;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_CacheUpperHalf(Tile *p_tile)
{
  IG_ASSERT(p_tile);

  if (p_tile == g_Tiles)
  {
    return;
  }

  u32 slot = R2O_AllocMapCacheSlot();
  p_tile->m_cache_slot_upper = slot;

  if (slot)
  {
    u32 map_idx = p_tile-g_Tiles;
    u32 ls      = (u32)&g_MapCache[slot << 8];
    u32 ea      = g_R2OCon.m_int_ea_maps + (map_idx << 10);

    DmaGet((volatile void *)(ls+0x000), ea+0x000, 256, g_kDmaTag);
    DmaGet((volatile void *)(ls+0x100), ea+0x200, 256, g_kDmaTag);
    DmaWaitAll(1 << g_kDmaTag);
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_CacheLowerHalf(Tile *p_tile)
{
  IG_ASSERT(p_tile);

  if (p_tile == g_Tiles)
  {
    return;
  }

  u32 slot = R2O_AllocMapCacheSlot();
  p_tile->m_cache_slot_lower = slot;

  if (slot)
  {
    u32 map_idx = p_tile-g_Tiles;
    u32 ls      = (u32)&g_MapCache[slot << 8];
    u32 ea      = g_R2OCon.m_int_ea_maps + (map_idx << 10);

    DmaGet((volatile void *)(ls+0x000), ea+0x100, 256, g_kDmaTag);
    DmaGet((volatile void *)(ls+0x100), ea+0x300, 256, g_kDmaTag);
    DmaWaitAll(1 << g_kDmaTag);
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_UncacheUpperHalf(Tile *p_tile)
{

  IG_ASSERT(p_tile);

  if (p_tile == g_Tiles)
  {
    return;
  }

  u32 slot = p_tile->m_cache_slot_upper;
  R2O_FreeMapCacheSlot(slot);
  p_tile->m_cache_slot_upper = 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_UncacheLowerHalf(Tile *p_tile)
{
  IG_ASSERT(p_tile);

  if (p_tile == g_Tiles)
  {
    return;
  }

  u32 slot = p_tile->m_cache_slot_lower;
  R2O_FreeMapCacheSlot(slot);
  p_tile->m_cache_slot_lower = 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_WritebackUpperHalf(Tile *p_tile)
{
  IG_ASSERT(p_tile);

  if (p_tile == g_Tiles)
  {
    return;
  }

  u32 slot = p_tile->m_cache_slot_upper;
  u32 map_idx = p_tile-g_Tiles;

  u32 ls = (u32)&g_MapCache[slot << 8];
  u32 ea = g_R2OCon.m_int_ea_maps + (map_idx << 10);

  DmaPut((volatile void *)(ls+0x000), ea+0x000, 256, g_kDmaTag);
  DmaPut((volatile void *)(ls+0x100), ea+0x200, 256, g_kDmaTag);
  DmaWaitAll(1 << g_kDmaTag);

  if (slot)
  {
    R2O_FreeMapCacheSlot(slot);
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_WritebackLowerHalf(Tile *p_tile)
{
  IG_ASSERT(p_tile);

  if (p_tile == g_Tiles)
  {
    return;
  }

  u32 slot = p_tile->m_cache_slot_lower;
  u32 map_idx = p_tile-g_Tiles;

  u32 ls = (u32)&g_MapCache[slot << 8];
  u32 ea = g_R2OCon.m_int_ea_maps + (map_idx << 10);

  DmaPut((volatile void *)(ls+0x000), ea+0x100, 256, g_kDmaTag);
  DmaPut((volatile void *)(ls+0x100), ea+0x300, 256, g_kDmaTag);
  DmaWaitAll(1 << g_kDmaTag);

  if (slot)
  {
    R2O_FreeMapCacheSlot(slot);
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

u32 R2O_TilePrecedes(u16 i0, u16 i1)
{
  f32 z0 = g_Tiles[i0].m_coords.m_f32[1];
  f32 z1 = g_Tiles[i1].m_coords.m_f32[1];

  if (z0 < z1)
  {
    return true;
  }

  if (z0 > z1)
  {
    return false;
  }

  f32 x0 = g_Tiles[i0].m_coords.m_f32[0];
  f32 x1 = g_Tiles[i1].m_coords.m_f32[0];

  if (x0 < x1)
  {
    return true;
  }

  return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

u32 R2O_TilePrecedes(Tile *p0, Tile *p1)
{
  f32 z0 = p0->m_coords.m_f32[1];
  f32 z1 = p1->m_coords.m_f32[1];

  if (z0 < z1)
  {
    return true;
  }

  if (z0 > z1)
  {
    return false;
  }

  f32 x0 = p0->m_coords.m_f32[0];
  f32 x1 = p1->m_coords.m_f32[0];

  if (x0 < x1)
  {
    return true;
  }

  return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_AddImpulses(Complex *fft_buf, Tile *p_tile)
{
  // get tile origin
  f32 x0 = p_tile->m_coords.m_f32[0];
  f32 z0 = p_tile->m_coords.m_f32[1];

  // some attenuaton params
  f32 R2    = g_R2OCon.m_int_attenuation_radius * g_R2OCon.m_int_attenuation_radius;
  f32 Rf    = 1.0f / g_R2OCon.m_int_attenuation_radius;
  f32 min   = g_R2OCon.m_int_attenuation_min;

  f32 inv_scale = 1.0f / g_R2OCon.m_int_scale;

  vu32 mask_0100 = (vu32){0,0xFFFFFFFF,0,0};

  // loop over impulses registered to this tile
  u32 idx = p_tile->m_first_impulse;
  u32 cnt = p_tile->m_num_impulses;
  for (u32 i=0; i<cnt; i++)
  {
    u32 i_imp = g_InteractiveData.m_assignments[idx++];

    // get impulse disk
    IMPULSE::ImpulseData impulse;
    u32 mainmem_idx = g_InteractiveData.m_impulses[i_imp];
    R2O_GetWait(&impulse, g_R2OCon.m_int_ea_impulses, mainmem_idx, sizeof(IMPULSE::ImpulseData));
    vf32 disk = impulse.m_volume;
    
    vf32 impulse_center = spu_sel(disk, g_WaterObject.m_origin, mask_0100);
    f32 rad = spu_extract(disk, 3);

    // compute a scaling factor based on location relative to camera
    f32 scaling_factor = inv_scale;

    for (u32 v=0; v<g_R2OCon.m_num_viewports; v++)
    {
      vf32 cam_pos = g_R2OCon.m_view_data[v].m_camera_position;
      vf32 diff = impulse_center - cam_pos;
      f32 r2 = VecLenSquared3(diff);
      if (r2 < R2)
      {
        f32 r = sqrtf(r2);
        f32 t = r * Rf;
        f32 new_scaling_factor = scaling_factor * ((1.0f-t) * min + t);
        if (new_scaling_factor < scaling_factor)
        {
          scaling_factor = new_scaling_factor;
        }
      }
    }

    // take impulse magnitude into account (defaults to 1.0)
    scaling_factor *= impulse.m_impulse_mag;

    // multiply by inner and outer displacement scale factors
    f32 inner_scaling_factor = scaling_factor * g_R2OCon.m_int_inner_displacement;
    f32 outer_scaling_factor = scaling_factor * g_R2OCon.m_int_outer_displacement;

    // determine rectangular region to modify
    f32 xc = spu_extract(disk, 0);
    f32 zc = spu_extract(disk, 2);
    i32 c0 = (i32)((xc-rad-x0)*inv_step) + 8;
    i32 c1 = (i32)((xc+rad-x0)*inv_step) + 8;
    i32 r0 = (i32)((zc-rad-z0)*inv_step) + 8;
    i32 r1 = (i32)((zc+rad-z0)*inv_step) + 8;

    // clamp the impulse rectangle (though we shouldn't really have to)
    if (c0 <  0) c0= 0;
    if (c1 > 31) c1=31;
    if (r0 <  0) r0= 0;
    if (r1 > 31) r1=31;

    f32 rad_sqd = rad*rad;
    f32 inner   = rad_sqd * (1.0f - g_R2OCon.m_int_inner_radius * g_R2OCon.m_int_inner_radius);

    // loop over this region
    f32 z = z0 - zc + (f32)(r0-8)*step;       // z coord relative to disk centre
    for (i32 r=r0; r<=r1; r++, z+=step)
    {
      f32 r2_z2 = rad_sqd - z*z;
      f32 x = x0 - xc + (f32)(c0-8)*step;     // x coord relative to disk centre
      for (i32 c=c0; c<=c1; c++, x+=step)
      {
        f32 r2_z2_x2 = r2_z2 - x*x;
        if (r2_z2_x2 > 0)
        {
          // perturb height
          f32 perturbation = inv_step * sqrtf(r2_z2 - x*x);
          scaling_factor = (r2_z2_x2 > inner) ? inner_scaling_factor : outer_scaling_factor;
          ((f32 *)fft_buf)[r*32+c] = perturbation * scaling_factor;
        }
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Update a map using 'overlap-save'
//
// - rework ffts so the real & imaginary parts are in separate arrays
// - replace the real heights output with the actual byte derivatives
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_UpdateMap(i16 *p_out, f32 *p_real_out,
                   Tile *p_UL,  Tile *p_U0,  Tile *p_UR,
                   Tile *p_0L,  Tile *p_00,  Tile *p_0R,
                   Tile *p_DL,  Tile *p_D0,  Tile *p_DR,
                   Complex *fft_buf, Complex *phase_shifts, u16 neighbor_flags)
{
  //------------------------------------------------------------------------
  /* build the fft buffer from the input maps
  
    +-------+-------+-------+
    |       |       |       |
    |       |       |       |
    |   +---|-------|---+   |
    |   |UL |  U0   |UR |   |
    |   |   |       |   |   |
    +-------+-------+-------+
    |   |0L |  00   |0R |   |
    |   |   |       |   |   |
    |   +---|-------|---+   |
    |   |0L |  00   |0R |   |
    |   |   |       |   |   |
    +-------+-------+-------+
    |   |DL |  D0   |DR |   |
    |   |   |       |   |   |
    |   +---|-------|---+   |
    |       |       |       |
    |       |       |       |
    +-------+-------+------*/

  i16 *UL = &g_MapCache[( (neighbor_flags & 0x080) ? p_UL->m_cache_slot_lower : 0) << 8];
  i16 *U0 = &g_MapCache[( (neighbor_flags & 0x040) ? p_U0->m_cache_slot_lower : 0) << 8];
  i16 *UR = &g_MapCache[( (neighbor_flags & 0x020) ? p_UR->m_cache_slot_lower : 0) << 8];

  i16 *uL = &g_MapCache[( (neighbor_flags & 0x010) ? p_0L->m_cache_slot_upper : 0) << 8];
  i16 *u0 = &g_MapCache[(!(neighbor_flags & 0x100) ? p_00->m_cache_slot_upper : 0) << 8];
  i16 *uR = &g_MapCache[( (neighbor_flags & 0x008) ? p_0R->m_cache_slot_upper : 0) << 8];

  i16 *dL = &g_MapCache[( (neighbor_flags & 0x010) ? p_0L->m_cache_slot_lower : 0) << 8];
  i16 *d0 = &g_MapCache[(!(neighbor_flags & 0x100) ? p_00->m_cache_slot_lower : 0) << 8];
  i16 *dR = &g_MapCache[( (neighbor_flags & 0x008) ? p_0R->m_cache_slot_lower : 0) << 8];

  i16 *DL = &g_MapCache[( (neighbor_flags & 0x004) ? p_DL->m_cache_slot_upper : 0) << 8];
  i16 *D0 = &g_MapCache[( (neighbor_flags & 0x002) ? p_D0->m_cache_slot_upper : 0) << 8];
  i16 *DR = &g_MapCache[( (neighbor_flags & 0x001) ? p_DR->m_cache_slot_upper : 0) << 8];

  //------------------------------------------------------------------------
  // update the heightmap in the Fourier domain

  if (g_R2OCon.m_paused)
  {
    R2O_PrepareFftBufferSeparate((f32 *)fft_buf, (f32 *)fft_buf+1024, UL, U0, UR, uL, u0, uR, dL, d0, dR, DL, D0, DR);
  }
  else
  {
    // gather/transform each row
    R2O_GatherAndTransformRows((f32 *)fft_buf, (f32 *)fft_buf+1024, UL, U0, UR, uL, u0, uR, dL, d0, dR, DL, D0, DR);

    // transform/phase-shift/inv transform each column
    R2O_TransformColsAsm((f32 *)fft_buf, (f32 *)phase_shifts, 32);

    // inverse transform each row
    R2O_InvFFT32_DIT_rows_sep_asm((f32 *)fft_buf, (f32 *)fft_buf+1024, 32);
  }


  //------------------------------------------------------------------------
  // add impulses (must precede thresholding, in case it's a blank tile)
  R2O_AddImpulses(fft_buf, p_00);

  //------------------------------------------------------------------------
  // copy the central section of the fft buffer to the output map
  R2O_WriteRectangleCentralSeparate(p_out, (f32 *)fft_buf+8*32+8, (f32 *)fft_buf+8*32+8+1024, 16*16);
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

u16 R2O_Threshold(f32 *heights)
{
  // compute threshold flags for propogation and persistence
  f32 thresh_create  = g_R2OCon.m_int_threshold_create;
  f32 thresh_destroy = g_R2OCon.m_int_threshold_destroy;

  u16 flags;

  flags = R2O_ThresholdRectangleWidth16Asm(heights+( 8<<5)+ 8, 16, thresh_destroy) << 8
        | R2O_ThresholdRectangleWidth16Asm(heights+( 4<<5)+ 8,  4, thresh_create ) << 6
        | R2O_ThresholdRectangleWidth4Asm (heights+( 8<<5)+ 4, 16, thresh_create ) << 4
        | R2O_ThresholdRectangleWidth4Asm (heights+( 8<<5)+24, 16, thresh_create ) << 3
        | R2O_ThresholdRectangleWidth16Asm(heights+(24<<5)+ 8,  4, thresh_create ) << 1;


  return flags;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_PropagateToNeighbor(Tile * volatile &p_neighbor, Complex *fft_buf,
                             f32 dx, f32 dz, f32 x, f32 z,
                             u16 neighbor_flags, u16 new_neighbor_flags, u16 threshold_flags, u16 flag)
{
  // quit if there's already a neigbour there
  if (neighbor_flags & flag)
  {
    return;
  }

  // quit if we haven't reached the threshold
  if (!(threshold_flags & flag))
  {
    return;
  }

  // add a new neighbour if there isn't one already
  if (!(new_neighbor_flags & flag))
  {
    // don't exceed the per-lod max tile count
    if (g_NumTiles >= 255)
    {
      return;
    }

    // no tile, so create one
    Tile *p_prev = p_neighbor;
    Tile *p_new  = R2O_AddTile(x, z, p_prev);

    if (!p_new)
    {
      //PRINT("couldn't add neighbour\n");
      return;
    }
    p_neighbor = p_new;

    g_NumTiles++;
  }

  // ensure the appropriate cache slots have been allocated (clearing half-maps where new allocations are made)
  if ((dz != -1.0f) && !p_neighbor->m_cache_slot_upper)
  {
    // horizontally adjacent or lower neighbour; make sure the upper half is allocated
    u32 slot = R2O_AllocMapCacheSlot();
    if (!slot)
    {
      return;
    }

    R2O_ClearMemAsm((qword *)&g_MapCache[slot << 8], 512);
    p_neighbor->m_cache_slot_upper = slot;
  }

  if ((dz != 1.0f) && !p_neighbor->m_cache_slot_lower)
  {
    // horizontally adjacent or upper neighbour; make sure the lower half is allocated
    u32 slot = R2O_AllocMapCacheSlot();
    if (!slot)
    {
      return;
    }

    R2O_ClearMemAsm((qword *)&g_MapCache[slot << 8], 512);
    p_neighbor->m_cache_slot_lower = slot;
  }



  // copy the relevant part(s) of the fft buf
  i32 dst_col = dx==-1 ? 12 : 0;
  if (dz == 0.0f)
  {
    // alongside; need two 4x8 rectangles
    i32 src_col = dx==-1 ? 4 : 24;

    IG_ASSERT(p_neighbor->m_cache_slot_upper);
    IG_ASSERT(p_neighbor->m_cache_slot_lower);

    R2O_WriteRectangleAdjacentSeparate(&g_MapCache[p_neighbor->m_cache_slot_upper << 8] + dst_col,
                                       (f32 *)fft_buf+(8*32)+src_col, (f32 *)fft_buf+(8*32)+src_col+1024, 32);

    R2O_WriteRectangleAdjacentSeparate(&g_MapCache[p_neighbor->m_cache_slot_lower << 8] + dst_col,
                                       (f32 *)fft_buf+(16*32)+src_col, (f32 *)fft_buf+(16*32)+src_col+1024, 32);
  }
  else
  {
    // above or below, or a corner neighbour; only need a single rectangle
    i32 src_col = dx==-1 ?  4 : dx==1 ? 24 : 8;
    i32 src_row = dz==-1 ?  4 : 24;
    i32 dst_row = dz==-1 ?  4 : 0;
    u32 slot    = dz==-1 ? p_neighbor->m_cache_slot_lower : p_neighbor->m_cache_slot_upper;
    IG_ASSERT(slot);
    if (dx == 0.0f)
    {
      // above/below neighbour
      R2O_WriteRectangleCentralSeparate(&g_MapCache[slot<<8] + (dst_row << 4) + dst_col,
                                        (f32 *)fft_buf+(src_row<<5)+src_col, (f32 *)fft_buf+(src_row<<5)+src_col+1024, 16*4);
    }
    else
    {
      // corner neighbour
      R2O_WriteRectangleAdjacentSeparate(&g_MapCache[slot<<8] + (dst_row << 4) + dst_col,
                                         (f32 *)fft_buf+(src_row<<5)+src_col, (f32 *)fft_buf+(src_row<<5)+src_col+1024, 16);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_PropagateToNeighbors(                Tile *p_new_U0,
                              Tile *p_new_0L, Tile *p_tile_00, Tile *p_new_0R,
                                              Tile *p_new_D0,
                              Complex *fft_buf, u16 neighbor_flags, u16 new_neighbor_flags, u16 threshold_flags)
{
  f32 tile_step = step * 16.0f;
  f32 x0 = p_tile_00->m_coords.m_f32[0];
  f32 z0 = p_tile_00->m_coords.m_f32[1];
  f32 xL = x0 - tile_step;
  f32 xR = x0 + tile_step;
  f32 zU = z0 - tile_step;
  f32 zD = z0 + tile_step;

  R2O_PropagateToNeighbor(p_new_D0, fft_buf,  0,  1, x0, zD, neighbor_flags, new_neighbor_flags, threshold_flags, 0x02);
  R2O_PropagateToNeighbor(p_new_0R, fft_buf,  1,  0, xR, z0, neighbor_flags, new_neighbor_flags, threshold_flags, 0x08);
  R2O_PropagateToNeighbor(p_new_0L, fft_buf, -1,  0, xL, z0, neighbor_flags, new_neighbor_flags, threshold_flags, 0x10);
  R2O_PropagateToNeighbor(p_new_U0, fft_buf,  0, -1, x0, zU, neighbor_flags, new_neighbor_flags, threshold_flags, 0x40);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

u16 R2O_MergeTileLists(u16 idx0, u16 idx1)
{
  // treat degenerate cases
  if (!idx0)
  {
    return idx1;
  }
  if (!idx1)
  {
    return idx0;
  }

  // init working idx
  u16 idx;
  if (R2O_TilePrecedes(idx0, idx1))
  {
    idx  = idx0;
    idx0 = g_Tiles[idx0].m_next_tile;
  }
  else
  {
    idx  = idx1;
    idx1 = g_Tiles[idx1].m_next_tile;
  }

  // record merged head
  u16 h = idx;

  // merge lists until one of them becomes null
  while (idx0 && idx1)
  {
    if (R2O_TilePrecedes(idx0, idx1))
    {
      g_Tiles[idx].m_next_tile = idx0;
      idx  = idx0;
      idx0 = g_Tiles[idx0].m_next_tile;
    }
    else
    {
      g_Tiles[idx].m_next_tile = idx1;
      idx  = idx1;
      idx1 = g_Tiles[idx1].m_next_tile;
    }
  }

  // concatenate the rest of the non-null list
  g_Tiles[idx].m_next_tile = idx0 | idx1;

  // return the head of the merged result
  return h;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Tile *R2O_UpdateMaps(Tile *tile_list, u32 lod, u32 ea_norm, f32 dt)
{
  // quit on empty list
  if (!tile_list)
  {
    return 0;
  }

  u32 map_width   = 16;
  u32 map_height  = 16;
  u32 map_stride  = map_width * map_height * 4;

  u32 fft_width   = 32;
  u32 fft_height  = 32;
  u32 fft_area    = 1024;

  // get tile step size
  f32 tile_step = step * 16.0f;

  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // initialize map cache
  R2O_InitMapCache();

  // alloc space for 1 output map
  i16 *updated_map = (i16 *)g_Scratchpad.Alloc(map_stride);

  // alloc fft buffer & phase shift buffer
  Complex *fft_buf      = (Complex *)g_Scratchpad.Alloc(fft_area * sizeof(Complex));
  Complex *phase_shifts = (Complex *)g_Scratchpad.Alloc(fft_area * sizeof(Complex));

  // temp: alloc a buffer for the real components and a big derivs buffer
  f32 *heights_buf = (f32 *)g_Scratchpad.Alloc(fft_area * sizeof(f32));
  i8  *derivs  = (i8 *)g_Scratchpad.Alloc(INT_NORM_STRIDE);

  // generate phase shifts
  f32 damping = powf(0.5f, dt/(g_R2OCon.m_int_half_life * sqrtf(step)));
  R2O_GeneratePhaseShifts(phase_shifts, lod, dt, damping);

  // a tile to act as the head of the list of propagated tiles
  Tile head_tile_new;
  head_tile_new.m_coords.m_u64 = 0xFF7FFFFFFF7FFFFFULL; // -MAX_FLOAT for x and z
  head_tile_new.m_next_tile = 0;

  // set up tile pointer pipeline
  f32 x0 = tile_list->m_coords.m_f32[0];
  f32 z0 = tile_list->m_coords.m_f32[1];
  f32 xL = x0 - tile_step;
  f32 xR = x0 + tile_step;
  f32 zU = z0 - tile_step;
  f32 zD = z0 + tile_step;
  Tile * volatile p_tile_UL = tile_list;
  Tile * volatile p_tile_U0 = tile_list;
  Tile * volatile p_tile_UR = tile_list;
  Tile * volatile p_tile_0L = tile_list;
  Tile * volatile p_tile_00 = tile_list;
  Tile * volatile p_tile_0R = NextTile(tile_list);
  Tile * volatile p_tile_DL = AdvanceTile(p_tile_0R, xL, zD);
  Tile * volatile p_tile_D0 = AdvanceTile(p_tile_DL, x0, zD);
  Tile * volatile p_tile_DR = AdvanceTile(p_tile_D0, xR, zD);
  Tile * volatile p_tile_cache_upper = tile_list;
  Tile * volatile p_tile_uncache_lower = tile_list;


  Tile * volatile p_new_U0 = &head_tile_new;
  Tile * volatile p_new_0L = &head_tile_new;
  Tile * volatile p_new_0R = &head_tile_new;
  Tile * volatile p_new_D0 = &head_tile_new;
  Tile * volatile p_tile_writeback_upper = &head_tile_new;
  Tile * volatile p_tile_writeback_lower = &head_tile_new;

  // cache one initial upper and one lower half
  //("caching one upper and one lower half\n");
  R2O_CacheUpperHalf(p_tile_00);
  R2O_CacheLowerHalf(p_tile_00);

  g_NormIdx  = 0;

  // loop over tiles in scan-list order
  while (1)
  {
    // cache enough upper halves
    while (p_tile_cache_upper != p_tile_DR)
    {
      p_tile_cache_upper = NextTile(p_tile_cache_upper);
      R2O_CacheUpperHalf(p_tile_cache_upper);
    }

    // cache one more lower half
    IG_ASSERT(p_tile_0R);
    R2O_CacheLowerHalf(p_tile_0R);

    // compute neighbour flags
    u16 neighbor_flags     = ((p_tile_00->m_flags & 0x02) != 0)                                     << 8
                           | (p_tile_UL->m_coords.m_f32[0]==xL && p_tile_UL->m_coords.m_f32[1]==zU) << 7
                           | (p_tile_U0->m_coords.m_f32[0]==x0 && p_tile_U0->m_coords.m_f32[1]==zU) << 6
                           | (p_tile_UR->m_coords.m_f32[0]==xR && p_tile_UR->m_coords.m_f32[1]==zU) << 5
                           | (p_tile_0L->m_coords.m_f32[0]==xL && p_tile_0L->m_coords.m_f32[1]==z0) << 4
                           | (p_tile_0R->m_coords.m_f32[0]==xR && p_tile_0R->m_coords.m_f32[1]==z0) << 3
                           | (p_tile_DL->m_coords.m_f32[0]==xL && p_tile_DL->m_coords.m_f32[1]==zD) << 2
                           | (p_tile_D0->m_coords.m_f32[0]==x0 && p_tile_D0->m_coords.m_f32[1]==zD) << 1
                           | (p_tile_DR->m_coords.m_f32[0]==xR && p_tile_DR->m_coords.m_f32[1]==zD) << 0;

    u16 new_neighbor_flags = ( p_new_U0->m_coords.m_f32[0]==x0 &&  p_new_U0->m_coords.m_f32[1]==zU) << 6
                           | ( p_new_0L->m_coords.m_f32[0]==xL &&  p_new_0L->m_coords.m_f32[1]==z0) << 4
                           | ( p_new_0R->m_coords.m_f32[0]==xR &&  p_new_0R->m_coords.m_f32[1]==z0) << 3
                           | ( p_new_D0->m_coords.m_f32[0]==x0 &&  p_new_D0->m_coords.m_f32[1]==zD) << 1;


    // update the current map
    R2O_UpdateMap(updated_map, heights_buf,
                  p_tile_UL, p_tile_U0, p_tile_UR,
                  p_tile_0L, p_tile_00, p_tile_0R,
                  p_tile_DL, p_tile_D0, p_tile_DR,
                  fft_buf,   phase_shifts, neighbor_flags);

    // copy the updated map back to main memory
    R2O_PutWait(updated_map, g_R2OCon.m_int_ea_maps, p_tile_00-g_Tiles, map_stride);

    // --------------------------------------------------------------------------------------
    // create normal map and upload to main mem

    // generate derivs from the real parts
    heights_buf = (f32 *)fft_buf;
    f32 *heights = heights_buf + fft_width * (fft_height - INT_NORM_HEIGHT)/2 + (fft_width - INT_NORM_WIDTH)/2;
    f32 scale = inv_step * 128.0f;
    R2O_GenerateMapDerivs(derivs, heights, INT_NORM_WIDTH, INT_NORM_HEIGHT, fft_width, 2*INT_NORM_WIDTH, false, scale);



    // copy the updated normal map back to main memory
    R2O_PutWait(derivs, ea_norm, g_NormIdx, INT_NORM_STRIDE);
    IG_ASSERT(g_NormIdx < 255);
    p_tile_00->m_normal_map_idx = ++g_NormIdx;

    // --------------------------------------------------------------------------------------
    // threshold the current map
    u16 threshold_flags = R2O_Threshold(heights_buf);

    // copy the 'persist' flag into the tile
    // and clear the 'blank' flag because it's only for use when the tile has just been created
    p_tile_00->m_flags = (p_tile_00->m_flags & 0xFC) | (threshold_flags>>8 & 0x01);

    // --------------------------------------------------------------------------------------
    // propagate to new tiles
    R2O_PropagateToNeighbors(          p_new_U0,
                             p_new_0L, p_tile_00, p_new_0R,
                                       p_new_D0,
                             fft_buf, neighbor_flags, new_neighbor_flags, threshold_flags);

    // --------------------------------------------------------------------------------------

    // uncache one upper half (provided we've got off the ground, i.e. this isn't the first pass)
    if (p_tile_0L != p_tile_00)
    {
      R2O_UncacheUpperHalf(p_tile_0L);
    }

    // write back enough new upper halves
    while (p_tile_writeback_upper != p_new_0L)
    {
      p_tile_writeback_upper = NextTile(p_tile_writeback_upper);
      R2O_WritebackUpperHalf(p_tile_writeback_upper);
    }

    // write back enough new lower halves
    while (p_tile_writeback_lower != p_new_U0)
    {
      p_tile_writeback_lower = NextTile(p_tile_writeback_lower);
      R2O_WritebackLowerHalf(p_tile_writeback_lower);
    }

    // advance tile pointers
    p_tile_0L = p_tile_00;
    p_tile_00 = p_tile_0R;
    if (p_tile_00 == g_Tiles)
    {
      break;
    }
    IG_ASSERT(p_tile_0R);
    p_tile_0R = NextTile(p_tile_0R);

    // identify neighbour coords
    x0 = p_tile_00->m_coords.m_f32[0];
    z0 = p_tile_00->m_coords.m_f32[1];
    xL = x0 - tile_step;
    xR = x0 + tile_step;
    zU = z0 - tile_step;
    zD = z0 + tile_step;

    p_tile_UL = AdvanceTile(p_tile_UL, xL, zU);
    p_tile_U0 = AdvanceTile(p_tile_U0, x0, zU);
    p_tile_UR = AdvanceTile(p_tile_UR, xR, zU);

    p_tile_DL = AdvanceTile(p_tile_DL, xL, zD);
    p_tile_D0 = AdvanceTile(p_tile_D0, x0, zD);
    p_tile_DR = AdvanceTile(p_tile_DR, xR, zD);


    p_new_U0 = AdvanceTile(p_new_U0, x0, zU);
    p_new_0L = AdvanceTile(p_new_0L, xL, z0);
    p_new_0R = AdvanceTile(p_new_0R, xR, z0);
    p_new_D0 = AdvanceTile(p_new_D0, x0, zD);


    // uncache enough lower halves
    while (p_tile_uncache_lower != p_tile_UL)
    {
      R2O_UncacheLowerHalf(p_tile_uncache_lower);
      p_tile_uncache_lower = NextTile(p_tile_uncache_lower);
    }
  }
  

  // write back any new upper halves still in the cache
  while (p_tile_writeback_upper != g_Tiles)
  {
    p_tile_writeback_upper = NextTile(p_tile_writeback_upper);
    R2O_WritebackUpperHalf(p_tile_writeback_upper);
  }

  // write back any new lower halves still in the cache
  while (p_tile_writeback_lower != g_Tiles)
  {
    p_tile_writeback_lower = NextTile(p_tile_writeback_lower);
    R2O_WritebackLowerHalf(p_tile_writeback_lower);
  }

  // add newly created tiles
  u16 merged_head = R2O_MergeTileLists(head_tile_new.m_next_tile, tile_list-g_Tiles);

  // restore scratchpad state
  g_Scratchpad.Reset(state);

  // return the start of the merged list
  return g_Tiles + merged_head;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Tile *R2O_Prune(Tile *p_tile)
{
  // find the first persistent tile (the zero tile is the terminator - guaranteed to be persistent)
  while (!(p_tile->m_flags & 0x01))
  {
    Tile *p_next = &g_Tiles[p_tile->m_next_tile];
    R2O_FreeTile(p_tile);
    p_tile = p_next;
    g_NumTiles--;
  }

  // record pointer so we can return it
  Tile *p_head = p_tile;

  // loop over remainder of list, deleting non-persistent entries
  while (p_tile->m_next_tile)
  {
    Tile *p_next = &g_Tiles[p_tile->m_next_tile];
    if (!(p_next->m_flags & 0x01))
    {
      p_tile->m_next_tile = p_next->m_next_tile;
      R2O_FreeTile(p_next);
      g_NumTiles--;
    }
    else
    {
      p_tile = p_next;
    }
  }

  return p_head;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// format the impulses registered to the current water object (also rejects the ones that didn't actually hit the water object)
// then register impulses to tiles
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_PreprocessImpulses()
{
  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // alloc mem for local impulses
  R2OImpulse *local_impulses = (R2OImpulse *)g_Scratchpad.Alloc(MAX_IMPULSES * sizeof(R2OImpulse));

  f32 sealevel = spu_extract(g_WaterObject.m_origin, 1);

  // generate the tag corresponding to the current water object's origin
  vu32 coords_u32 = *(vu32 *)&g_WaterObject.m_origin;
  u32 tag = spu_extract(coords_u32, 0) ^ spu_extract(coords_u32, 2);

  // get the holes for this level
  u32 num_holes = g_R2OCon.m_num_holes;
  if (num_holes)
  {
    g_Holes = (R2OHole *)R2O_AllocGetWait(g_R2OCon.m_ea_holes, 0, num_holes * sizeof(R2OHole));
  }

  // loop over all impulses registered to the current water object
  const u32 cnt_in = g_InteractiveData.m_num_impulses;
  u32 cnt_out = 0;
  for (u32 i_imp=0; i_imp<cnt_in; i_imp++)
  {
    // get the impulse from main memory
    IMPULSE::ImpulseData impulse;
    u32 mainmem_idx = g_InteractiveData.m_impulses[i_imp];
    R2O_GetWait(&impulse, g_R2OCon.m_int_ea_impulses, mainmem_idx, sizeof(IMPULSE::ImpulseData));

    // extract swept sphere info
    vf32 p0 = impulse.m_volume;
    f32 sph_rad = spu_extract(p0, 3);
    vf32 dir = impulse.m_impulse;
    f32 length = spu_extract(dir, 3);
    vf32 p1 = VecMulAdd4(dir, length, p0);
    f32 y0 = spu_extract(p0, 1);
    f32 y1 = spu_extract(p1, 1);

    f32 dy0 = sealevel - y0;
    f32 dy1 = sealevel - y1;

    // convert swept sphere to a suitable disk
    f32 imp_rad;
    vf32 disk;
    if (dy0 * dy1 < 0)
    {
      // endpoints straddle; intersect sweep vector with plane
      f32 t = dy0/(y1-y0);
      disk = VecMulAdd4(p1-p0, t, p0);
      imp_rad = sph_rad;
    }
    else
    {
      // only one sphere in contact with plane
      f32 diff;
      if (fabsf(dy0) < fabsf(dy1))
      {
        // use p0
        disk = p0;
        diff = sph_rad*sph_rad - dy0*dy0;
      }
      else
      {
        // use p1
        disk = p1;
        diff = sph_rad*sph_rad - dy1*dy1;
      }
      imp_rad = diff>0.0f ? sqrtf(diff) : 0.0f;
    }
    disk = spu_insert(imp_rad, disk, 3);


    // 2D rejection test
    f32 x = spu_extract(disk, 0);
    f32 z = spu_extract(disk, 2);
    vf32 min = g_WaterObject.m_origin;
    vf32 max = g_WaterObject.m_dimensions + min;
    f32 x_min = spu_extract(min, 0);
    f32 z_min = spu_extract(min, 2);
    f32 x_max = spu_extract(max, 0);
    f32 z_max = spu_extract(max, 2);
    if (!(x>x_min && x<x_max && z>z_min && z<z_max))
    {
      continue;
    }

    // distance-based rejection test
    vf32 pos = spu_insert(sealevel, disk, 1);
    f32 dist_threshold = sph_rad * g_R2OCon.m_int_cull_distance;
    f32 min_dist = 1.0e30f;
    for (u32 v=0; v<g_R2OCon.m_num_viewports; v++)
    {
      vf32 rel = pos - g_R2OCon.m_view_data[v].m_camera_position;
      f32 dist = VecLen3(rel);
      if (dist < min_dist)
      {
        min_dist = dist;
      }
    }
    if (min_dist > dist_threshold)
    {
      continue;
    }

    // reject if the intersection point falls inside a cut hole
    u32 in_a_hole = false;
    R2OHole *p_hole = g_Holes;
    for (u32 i_hole=0; i_hole<num_holes; i_hole++,p_hole++)
    {
      // match the tag
      if (p_hole->m_tag != tag)
      {
        continue;
      }

      // compute the hole bounds
      f32 r = (f32)(0x10000 >> (p_hole->m_lod >> 1)) * 0.001953125f;
      f32 d = 2.0f * r;
      x_min = p_hole->m_xcoord - r;
      z_min = p_hole->m_zcoord - r;
      f32 cnt = (f32)p_hole->m_cnt;
      if (p_hole->m_dir)
      {
        x_max = x_min + d;
        z_max = z_min + d * cnt;
      }
      else
      {
        x_max = x_min + d * cnt;
        z_max = z_min + d;
      }

      // test (x,z) against hole
      if (x>x_min && x<x_max && z>z_min && z<z_max)
      {
        in_a_hole = true;
        break;
      }
    }
    if (in_a_hole)
    {
      continue;
    }

    // record this impulse
    R2OImpulse *p_impulse = &local_impulses[cnt_out++];
    p_impulse->m_disk  = disk;
    p_impulse->m_flags = 0;
    p_impulse->m_magnitude = impulse.m_impulse_mag;

    // compute lod range
    f32 ratio = INT_STEP0 / sph_rad;
    p_impulse->m_min_l = QUICK_LOG2(1.5f * ratio);
    if (min_dist > dist_threshold * g_R2OCon.m_int_cull_dist_ratio)
    {
      p_impulse->m_max_l = p_impulse->m_min_l + 1;  // too far; only span 1 lod: radius range is 1.5*step to 3.0*step
    }
    else
    {
      p_impulse->m_max_l = p_impulse->m_min_l + 2;  // span 2 lods: radius range is 1.5*step to 6.0*step
    }

    // write the intersection information to the main-mem impulse m_volume
    // (which is a very naughty thing to do, but the splash-effect callback needs this info)
    impulse.m_volume = disk;
    impulse.m_flags |= IMPULSE::IMPULSE_FLAG_SPLASH;
    R2O_PutWait(&impulse, g_R2OCon.m_int_ea_impulses, mainmem_idx, sizeof(IMPULSE::ImpulseData));

    // compact the list of impulse indices attached to this water object
    g_InteractiveData.m_impulses[cnt_out-1] = mainmem_idx;
  }

  // record the number of impulses we output
  g_InteractiveData.m_num_impulses = cnt_out;


  // ------------------------------------------------------------------------
  // register impulses to tiles


  // init step & inv_step
  f32 step=128.0f, inv_step=0.0078125f;

  u32 registration_idx = 0;

  u32 first_impulse_idx = 0;
  u32 last_impulse_idx  = 0;
  u32 num_impulses = g_InteractiveData.m_num_impulses;

  // loop over interactive lods for the current water object
  for (u32 l=0; l<MAX_INT_LODS; l++)
  {
    // tile step & inv
    vf32 tile_step     = spu_splats(16.0f   * step);
    vf32 inv_tile_step = spu_splats(0.0625f * inv_step);

    // set range of impulses to examine (the ppu kindly put them in sorted order for us)
    while ((first_impulse_idx < num_impulses) && (local_impulses[first_impulse_idx].m_max_l < l))
    {
      first_impulse_idx++;
    }
    while ((last_impulse_idx < num_impulses) && (local_impulses[last_impulse_idx].m_min_l <= l))
    {
      last_impulse_idx++;
    }

    // flag for this lod
    u32 flag = 1<<l;


    // loop over tiles for this lod
    u16 head_idx = g_InteractiveData.m_tile_list[l];
    u16 tile_idx = head_idx;
    while (tile_idx)
    {
      Tile *p_tile = g_Tiles + tile_idx;

      // set start index into array of registered impulses
      p_tile->m_first_impulse = registration_idx;

      // get bounds of tile
      f32 x0 = p_tile->m_coords.m_f32[0];
      f32 z0 = p_tile->m_coords.m_f32[1];

      f32 x1 = x0 + spu_extract(tile_step, 0);
      f32 z1 = z0 + spu_extract(tile_step, 0);

      // loop over range of impulses designated for this lod (the ppu kindly put them in sorted order for us)
      for (u32 i_imp=first_impulse_idx; i_imp<last_impulse_idx; i_imp++)
      {
        // get disk
        vf32 disk = local_impulses[i_imp].m_disk;
        f32 xc  = spu_extract(disk, 0);
        f32 zc  = spu_extract(disk, 2);
        f32 rad = spu_extract(disk, 3);

        // check for overlap with tile
        if (xc+rad>x0 && xc-rad<x1 && zc+rad>z0 && zc-rad<z1)
        {
          // overlap: register this disk to this tile
          IG_ASSERT(registration_idx < MAX_IMPULSES * 16);
          g_InteractiveData.m_assignments[registration_idx++] = i_imp;

          // mark the impulse as used by this lod
          local_impulses[i_imp].m_flags |= flag;
        }
      }

      // count how many impulses we registered to this tile
      p_tile->m_num_impulses = registration_idx - p_tile->m_first_impulse;

      // next tile
      tile_idx = p_tile->m_next_tile;
    }


    // clear new tiles array
    u16 new_tiles[MAX_IMPULSES];
    for (u32 i=0; i<MAX_IMPULSES; i++)
    {
      new_tiles[i] = 0;
    }
    u32 new_tile_cnt = 0;
    IG_ASSERT(g_InteractiveData.m_num_tiles[l] <= 255);
    u32 max_new_tiles = 255 - g_InteractiveData.m_num_tiles[l];

    // generate new tiles for impulses which didn't overlap any exisiting tiles on this lod
    // loop over range of impulses designated for this lod (the ppu kindly put them in sorted order for us)
    for (u32 i_imp=first_impulse_idx; i_imp<last_impulse_idx; i_imp++)
    {
      // make sure we can still allocate a tile if we need to
      if (new_tile_cnt >= max_new_tiles)
      {
        break;
      }

      // check whether already used in this lod
      if (local_impulses[i_imp].m_flags & flag)
      {
        continue;
      }

      // generate coords for new tile
      vf32 disk = local_impulses[i_imp].m_disk;
      vi32 magic = (vi32)(0x4B400000);
      vf32 vcoords = spu_madd(disk, inv_tile_step, (vf32)magic);
      vcoords = spu_sub(vcoords, (vf32)magic);
      vcoords = spu_mul(vcoords, tile_step);
      Coords coords;
      coords.m_f32[0] = spu_extract(vcoords, 0);
      coords.m_f32[1] = spu_extract(vcoords, 2);

      // check that another impulse hasn't already generated a new tile here
      u32 been_there_done_that = false;
      for (u32 i=0; i<new_tile_cnt; i++)
      {
        Tile *p_tile = g_Tiles + new_tiles[i];
        if (p_tile->m_coords.m_u64 == coords.m_u64)
        {
          been_there_done_that = true;
          break;
        }
      }
      if (been_there_done_that)
      {
        continue;
      }

      // this impulse hasn't been used by this lod, and is needed; alloc a new tile for it
      Tile *p_new = R2O_AllocTile();
      if (p_new)
      {
        p_new->m_coords = coords;
        p_new->m_next_tile = 0;
        p_new->m_flags = 0x03;    // blank + persist

        p_new->m_first_impulse = registration_idx;
        p_new->m_num_impulses  = 1;

        IG_ASSERT(registration_idx < MAX_IMPULSES * 16);
        g_InteractiveData.m_assignments[registration_idx++] = i_imp;

        new_tiles[new_tile_cnt++] = p_new - g_Tiles;
      }
      else
      {
        //PRINT("couldn't add new tile for impulse\n");
        break;
      }
    }

    // merge-sort the new tiles
    for (u32 stride=1; stride<new_tile_cnt; stride<<=1)
    {
      for (u32 i0=0; i0<new_tile_cnt; i0+=2*stride)
      {
        u32 i1 = i0+stride;
        new_tiles[i0] = R2O_MergeTileLists(new_tiles[i0], new_tiles[i1]);
      }
    }

    // merge-sort the new tiles with the existing ones
    u16 merged_head = R2O_MergeTileLists(new_tiles[0], head_idx);
    g_InteractiveData.m_tile_list[l] = merged_head;
    g_InteractiveData.m_num_tiles[l] += new_tile_cnt;

    // step step & inv_step
    step     *= 0.5f;
    inv_step *= 2.0f;
  }


  // restore scratchpad state
  g_Scratchpad.Reset(state);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_InitIndexMap(u8 *map, u32 width, u32 height)
{
  u32 cnt = width * height;
  for (u32 i=0; i<cnt; i++)
  {
    *map++ = 0;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

u32 R2O_WriteIndexMap(u8 *map, Tile *p_tile, f32 tile_step, f32 inv_tile_step,
                      u32 width, u32 height, f32 x_orig, f32 z_orig)
{
  u32 norm_idx = 0xFFFFFFFF;

  u32 idx = p_tile-g_Tiles;
  while (idx)
  {
    // get tile coords
    f32 x = p_tile->m_coords.m_f32[0];
    f32 z = p_tile->m_coords.m_f32[1];

    // convert to col/row
    i32 col = (i32)((x - x_orig) * inv_tile_step);
    i32 row = (i32)((z - z_orig) * inv_tile_step);

    // test for overlap with window
    if (col>=0 && row>=0 && col<(i32)width && row<(i32)height)
    {
      // compute pos within map
      u32 pos = row * width + col;

      // get normal map index
      norm_idx = p_tile->m_normal_map_idx;

      // write the index into the map entry
      map[pos] = norm_idx;
    }

    // next tile in list
    idx = p_tile->m_next_tile;
    p_tile = g_Tiles + idx;
  }

  // check whether we wrote anything to the map
  return (norm_idx != 0xFFFFFFFF);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_AddInteractiveMaps(f32 *heights, i32 width, i32 height, u32 l)
{
  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // get window extents
  f32 x_min = spu_extract(origin_world, 0);
  f32 z_min = spu_extract(origin_world, 2);
  f32 x_max = x_min + (f32)width  * step;
  f32 z_max = z_min + (f32)height * step;

  // tile width ( = height)
  f32 tile_step = 16.0f * step;

  // only need to get real parts
  const u32 real_size    = 16*16*sizeof(i16);
  const u32 complex_size = 2*real_size;

  // alloc mem for real parts of 1 interactive map
  i16 *real_map = (i16 *)g_Scratchpad.Alloc(real_size); 
  Tile tile;

  // set scale for this lod
  f32 scale = step * (1.0f/32768.0f) * g_R2OCon.m_int_scale;

  // loop over tiles for this lod
  u16 tile_idx = g_InteractiveData.m_tile_list[l];
  do
  {
    // get the tile
    R2O_GetWait(&tile, g_R2OCon.m_int_ea_tiles, tile_idx, sizeof(Tile));

    // get coords
    f32 x0 = tile.m_coords.m_f32[0];
    f32 z0 = tile.m_coords.m_f32[1];
    f32 x1 = x0 + tile_step;
    f32 z1 = z0 + tile_step;

    // test for overlap with the current window
    if (x1>x_min && x0<x_max && z1>z_min && z0<z_max)
    {
      // get the real parts of the interactive map
      R2O_GetWait(real_map, g_R2OCon.m_int_ea_maps, tile_idx, real_size, complex_size);

      // translate dest coords to col/row values
      i32 c0 = (i32)((x0-x_min) * inv_step);
      i32 r0 = (i32)((z0-z_min) * inv_step);
      i32 c1 = c0 + 16;
      i32 r1 = r0 + 16;

      // init pointers
      i16 *p_src = real_map;
      f32 *p_dst = heights + (r0 * width) + c0;

      // clip interactive map to current window
      if (c0 < 0)
      {
        p_src -= c0;
        p_dst -= c0;
        c0 = 0;
      }
      if (c1 > width)
      {
        c1 = width;
      }
      if (r0 < 0)
      {
        p_src -= r0*16;
        p_dst -= r0*width;
        r0 = 0;
      }
      if (r1 > height)
      {
        r1 = height;
      }

      // compute end-of-row skip values
      i32 src_skip = 16    - (c1-c0);
      i32 dst_skip = width - (c1-c0);

      // iterate over points
      for (i32 r=r0; r<r1; r++)
      {
        for (i32 c=c0; c<c1; c++)
        {
          *p_dst++ += (f32)*p_src++ * scale;
        }

        p_src += src_skip;
        p_dst += dst_skip;
      }
    }

    // next tile in list
    tile_idx = tile.m_next_tile;
  }
  while (tile_idx);

  // restore scratchpad state
  g_Scratchpad.Reset(state);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_UpdateInteractive()
{
  if (!g_R2OCon.m_int_enable)
  {
    return;
  }

  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // allocate memory for interactive normal maps
  //
  // If all the current maps are to be accessible within the same texture, and any one map is to be
  // uploaded in a single dma transfer, and we want more than 227 maps (=4096/18), then we have to
  // use a 3D texture. For bilinear continuity across tile boundaries each tile must have a 1-texel
  // border, implying 18x18 is the minimum size. The texture need not actually use border mode, but
  // calculating texture coords might be easier if it does. The start of each map must be qword
  // aligned for dma purposes, and the row stride must be 18 or more. This implies the smallest map
  // size we can use is 18x20, or 20x18 (at 2 bytes per texel). We'll choose 20x18 since the
  // storage alignment pattern cycles every 5 qwords, vs every 9. And, actually, we'll need to use
  // 24x18 until some extra code is written to handle 20x18.
  u32 size = g_R2OCon.m_int_max_tiles * INT_NORM_STRIDE;
  u32 ea_norm = IGG::DynamicRenderAlloc(g_R2OCon.m_allocator, size, 128, false);
  g_R2OCon.m_int_ea_normal_maps = ea_norm;

  // load the index map we use for decompressing paletted phasor arrays
  g_IndexMap = (u8 *)R2O_AllocGetWait(g_R2OCon.m_int_ea_index_map, 0, 1024);

  // get all the tiles
  size = g_R2OCon.m_int_max_tiles * sizeof(Tile);
  g_Tiles = (Tile *)R2O_AllocGetWait(g_R2OCon.m_int_ea_tiles, 0, size);

  // loop over water objects
  u32 ea  = g_R2OCon.m_ea_water_objects;
  u32 ofs = OFFSETOF(R2OWaterObject, m_interactive_data);
  for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++)
  {
    // get the water object and interactive data
    R2O_GetWait(&g_WaterObject, ea, i_obj, sizeof(R2OBasicWaterObject), sizeof(R2OWaterObject));
    if (!(g_WaterObject.m_flags & R2O_WATER_OBJECT_FLAG_ACTIVE))
    {
      continue;
    }
    R2O_GetWait(&g_InteractiveData, ea+ofs, i_obj, sizeof(R2OInteractiveData), sizeof(R2OWaterObject));

    // update interactive height maps & normal map


    // propagate
    step=128.0f, inv_step=0.0078125f;
    for (u32 l=0; l<MAX_INT_LODS; l++)
    {
      // prepare impulse data
      if (l==0)
      {
        R2O_PreprocessImpulses();
      }

      if (g_InteractiveData.m_num_tiles[l])
      {
        // set the normal map base pointer for this lod
        g_InteractiveData.m_ea_normal_maps[l] = ea_norm;

        // update the heightmaps (also generates normal maps)
        u32 lod = (l<<1);
        Tile *tile_list = g_Tiles + g_InteractiveData.m_tile_list[l];
        g_NumTiles = g_InteractiveData.m_num_tiles[l];
        g_InteractiveData.m_tile_list[l] = R2O_UpdateMaps(tile_list, lod, ea_norm, g_R2OCon.m_time_step) - g_Tiles;

        // step the normal map base pointer
        ea_norm += g_NumTiles * INT_NORM_STRIDE;

        // pruning
        g_InteractiveData.m_tile_list[l] = R2O_Prune(g_Tiles + g_InteractiveData.m_tile_list[l]) - g_Tiles;
        g_InteractiveData.m_num_tiles[l] = g_NumTiles;
      }

      step     *= 0.5f;
      inv_step *= 2.0f;
    }


    // put back the interactive data
    R2O_PutWait(&g_InteractiveData, ea+ofs, i_obj, sizeof(R2OInteractiveData), sizeof(R2OWaterObject));
  }

  // put all the tiles back
  size = g_R2OCon.m_int_max_tiles * sizeof(Tile);
  R2O_PutWait(g_Tiles, g_R2OCon.m_int_ea_tiles, 0, size);

  // balance the budget
  {
    i32 num   = g_R2OCon.m_int_num_tiles;
    i32 ideal = g_R2OCon.m_int_ideal_num_tiles;
    if (num > ideal)
    {
      // damp things a bit
      f32 surplus = (f32)(num - ideal) / (f32)ideal;
      f32 damping = 1.0f - g_R2OCon.m_int_damping_coeff0 * surplus;
      if (!g_R2OCon.m_int_free_tile)
      {
        damping = 1.0f - g_R2OCon.m_int_damping_coeff1 * surplus;
      }
      f32 reduced_half_life = g_R2OCon.m_int_half_life * damping;
      g_R2OCon.m_int_half_life = MIN(reduced_half_life, 3.0f);
      g_R2OCon.m_int_cull_dist_ratio *= g_R2OCon.m_int_cull_dist_coeff1;
    }
    else
    {
      // lengthen half-life, let cull distance creep up
      g_R2OCon.m_int_half_life = 10.0f;
      g_R2OCon.m_int_cull_dist_ratio = MIN(g_R2OCon.m_int_cull_dist_ratio * g_R2OCon.m_int_cull_dist_coeff0, 0.5f);
    }
  }

  // restore scratchpad state
  g_Scratchpad.Reset(state);
}


#if R2O_TASKS

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// task-based version: pre-update
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

u32 s_state;

void R2O_PreUpdateInteractive()
{
  // allocate memory for interactive normal maps
  u32 size = g_R2OCon.m_int_max_tiles * INT_NORM_STRIDE;
  g_R2OCon.m_ea_norm = IGG::DynamicRenderAlloc(g_R2OCon.m_allocator, size, 128, false);
  g_R2OCon.m_int_ea_normal_maps = g_R2OCon.m_ea_norm;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// task-based version: update pre-object
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_UpdateInteractivePreObject(u32 i_obj)
{
  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // get all the tiles
  u32 size = g_R2OCon.m_int_max_tiles * sizeof(Tile);
  g_Tiles = (Tile *)R2O_AllocGetWait(g_R2OCon.m_int_ea_tiles, 0, size);

  // get the water object and interactive data
  u32 ea  = g_R2OCon.m_ea_water_objects;
  u32 ofs = OFFSETOF(R2OWaterObject, m_interactive_data);
  R2O_GetWait(&g_WaterObject, ea, i_obj, sizeof(R2OBasicWaterObject), sizeof(R2OWaterObject));
  R2O_GetWait(&g_InteractiveData, ea+ofs, i_obj, sizeof(R2OInteractiveData), sizeof(R2OWaterObject));

  // update interactive height maps & normal map

  // prepare impulse data
  R2O_PreprocessImpulses();

  // put back the interactive data
  R2O_PutWait(&g_InteractiveData, ea+ofs, i_obj, sizeof(R2OInteractiveData), sizeof(R2OWaterObject));

  // put all the tiles back
  size = g_R2OCon.m_int_max_tiles * sizeof(Tile);
  R2O_PutWait(g_Tiles, g_R2OCon.m_int_ea_tiles, 0, size);

  // restore scratchpad state
  g_Scratchpad.Reset(state);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// task-based version: update
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_UpdateInteractive(u32 i_obj, u32 lod)
{
  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // load the index map we use for decompressing paletted phasor arrays
  g_IndexMap = (u8 *)R2O_AllocGetWait(g_R2OCon.m_int_ea_index_map, 0, 1024);

  // get all the tiles
  u32 size = g_R2OCon.m_int_max_tiles * sizeof(Tile);
  g_Tiles = (Tile *)R2O_AllocGetWait(g_R2OCon.m_int_ea_tiles, 0, size);

  // get the water object and interactive data
  u32 ea  = g_R2OCon.m_ea_water_objects;
  u32 ofs = OFFSETOF(R2OWaterObject, m_interactive_data);
  R2O_GetWait(&g_WaterObject, ea, i_obj, sizeof(R2OBasicWaterObject), sizeof(R2OWaterObject));
  R2O_GetWait(&g_InteractiveData, ea+ofs, i_obj, sizeof(R2OInteractiveData), sizeof(R2OWaterObject));

  // update interactive height maps & normal map
  u32 l = lod>>1;
  if (g_InteractiveData.m_num_tiles[l])
  {
    step     = 128.0f     * powf(0.5f, (f32)l);
    inv_step = 0.0078125f * powf(2.0f, (f32)l);

    // set the normal map base pointer for this lod
    g_InteractiveData.m_ea_normal_maps[l] = g_R2OCon.m_ea_norm;

    // update the heightmaps (also generates normal maps)
    Tile *tile_list = g_Tiles + g_InteractiveData.m_tile_list[l];
    g_NumTiles = g_InteractiveData.m_num_tiles[l];
    g_InteractiveData.m_tile_list[l] = R2O_UpdateMaps(tile_list, lod, g_R2OCon.m_ea_norm, g_R2OCon.m_time_step) - g_Tiles;

    // step the normal map base pointer
    g_R2OCon.m_ea_norm += g_NumTiles * INT_NORM_STRIDE;

    // pruning
    g_InteractiveData.m_tile_list[l] = R2O_Prune(g_Tiles + g_InteractiveData.m_tile_list[l]) - g_Tiles;
    g_InteractiveData.m_num_tiles[l] = g_NumTiles;
  }

  // put back the interactive data
  R2O_PutWait(&g_InteractiveData, ea+ofs, i_obj, sizeof(R2OInteractiveData), sizeof(R2OWaterObject));

  // put all the tiles back
  size = g_R2OCon.m_int_max_tiles * sizeof(Tile);
  R2O_PutWait(g_Tiles, g_R2OCon.m_int_ea_tiles, 0, size);

  // restore scratchpad state
  g_Scratchpad.Reset(state);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// task-based version: post-update
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_PostUpdateInteractive()
{
  // balance the budget
  {
    i32 num   = g_R2OCon.m_int_num_tiles;
    i32 ideal = g_R2OCon.m_int_ideal_num_tiles;
    if (num > ideal)
    {
      // damp things a bit
      f32 surplus = (f32)(num - ideal) / (f32)ideal;
      f32 damping = 1.0f - g_R2OCon.m_int_damping_coeff0 * surplus;
      if (!g_R2OCon.m_int_free_tile)
      {
        damping = 1.0f - g_R2OCon.m_int_damping_coeff1 * surplus;
      }
      f32 reduced_half_life = g_R2OCon.m_int_half_life * damping;
      g_R2OCon.m_int_half_life = MIN(reduced_half_life, 3.0f);
      g_R2OCon.m_int_cull_dist_ratio *= g_R2OCon.m_int_cull_dist_coeff1;
    }
    else
    {
      // lengthen half-life, let cull distance creep up
      g_R2OCon.m_int_half_life = 10.0f;
      g_R2OCon.m_int_cull_dist_ratio = MIN(g_R2OCon.m_int_cull_dist_ratio * g_R2OCon.m_int_cull_dist_coeff0, 0.5f);
    }
  }
}

#endif


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// generate index maps
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_GenerateIndexMaps()
{
  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // allocate spu mem for index map
  u8 *index_map = (u8 *)g_Scratchpad.Alloc(16 * 1024);

  // get all the tiles (note: we're not getting the heightmaps here, just the 16-byte tile structs)
  u32 size = g_R2OCon.m_int_max_tiles * sizeof(Tile);
  g_Tiles = (Tile *)R2O_AllocGetWait(g_R2OCon.m_int_ea_tiles, 0, size);

  // init step & inv
  step=128.0f, inv_step=0.0078125f;

  // default: no non-zero index maps
  for (u32 l=0; l<MAX_INT_LODS; l++)
  {
    g_RenderData.m_ea_index_map[l] = 0;
  }

  // loop over lods for the current water object
  for (u32 l=0; l<MAX_INT_LODS; l++)
  {
    // get dimensions
    u32 lod = (l<<1);
    u16 cols_rows = g_RenderData.m_cols_rows[lod];
    u32 cols = (u32)cols_rows >> 8;
    u32 rows = (u32)cols_rows & 0xFF;

    // quit the loop if no render output for this lod
    if (!(cols && rows))
    {
      break;
    }

    IG_ASSERT(cols>16 && rows>16);

    // get origin
    f32 x_orig = g_RenderData.m_origins[lod][0];
    f32 z_orig = g_RenderData.m_origins[lod][2];

    // strip off padding
    cols -= 16;
    rows -= 16;
    x_orig += 8.0f * step;
    z_orig += 8.0f * step;

    // halve the dimensions because each tile covers 2 qridsquares (at the lod where it first becomes visible)
    cols >>= 1;
    rows >>= 1;

    IG_ASSERT(cols*rows <= 16*1024);

    // clear the index map
    R2O_InitIndexMap(index_map, cols, rows);

    // write tiles (if any) into the index map
    u32 num_tiles = 0;
    if (l+3 < MAX_INT_LODS)
    {
      num_tiles = g_InteractiveData.m_num_tiles[l+3];
    }
    u32 b_used = false;
    if (num_tiles)
    {
      f32 tile_step     = 2.0f * step;
      f32 inv_tile_step = 0.5f * inv_step;
      Tile *tile_list   = g_Tiles + g_InteractiveData.m_tile_list[l+3];
      b_used = R2O_WriteIndexMap(index_map, tile_list, tile_step, inv_tile_step, cols, rows, x_orig, z_orig);
    }

    // only writeback non-trivial index maps
    if (b_used)
    {
      // allocate main memory for index map
      u32 size = cols * rows;
      size = (size + 15) & -16;   // round up to qword multiple
      u32 ea_idx = IGG::DynamicRenderAlloc(g_R2OCon.m_allocator, size, 128, false);
      g_RenderData.m_ea_index_map[l] = ea_idx;

      // pass back the index map
      DmaPut((volatile void *)index_map, ea_idx, size, g_kDmaTag);
      DmaWaitAll(1 << g_kDmaTag);
    }

    step     *= 0.5f;
    inv_step *= 2.0f;
  }

  // restore scratchpad state
  g_Scratchpad.Reset(state);
}

