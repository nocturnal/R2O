#pragma once


#define USE_INTERACTIVE         1
#define MAX_INT_LODS           20
#define INT_STEP0          128.0f
#define MAX_IMPULSES           32   // a power of 2
#define INT_NORM_WIDTH         24
#define INT_NORM_HEIGHT        18
#define NUM_MAP_CACHE_SLOTS   192
#define OPTIMIZE_INTERACTIVE    1
#define SANITY_CHECKS           0

#define INT_NORM_STRIDE (INT_NORM_WIDTH * INT_NORM_HEIGHT * 2)




union Coords
{
  u64 m_u64;
  f32 m_f32[2];
};


struct Tile
{
  Coords  m_coords;           // world (x,z)
  u16     m_next_tile;        // index of next tile in scan-list order
  u8      m_flags;
  u8      m_normal_map_idx;
  u8      m_cache_slot_upper;
  u8      m_cache_slot_lower;
  u8      m_first_impulse;
  u8      m_num_impulses;
}
ALIGNED(16);



// per-water object data
struct R2OInteractiveData
{
  u32   m_ea_normal_maps[MAX_INT_LODS];
  u16   m_tile_list     [MAX_INT_LODS];
  u16   m_num_tiles     [MAX_INT_LODS];
  u16   m_impulses      [MAX_IMPULSES];
  u8    m_assignments   [MAX_IMPULSES * 16];

  u32   m_impulse_lods;
  u16   m_num_impulses;
}
ALIGNED(16);

