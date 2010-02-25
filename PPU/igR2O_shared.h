#pragma once
#include "igR2OQuery_shared.h"
#include "igR2OInteractive_shared.h"

#ifdef R2O_PPU
#include "igMem/igMemory.h"
#endif

#define R2O_MAX_VIEWPORTS  2
#define MAX_LODS          40
#define R2O_TASKS          0


#define TEXTURE_FORMAT(LEVELS, UNNORM, LINEAR, FORMAT, DIM, NOBORD, CUBE, MAINMEM, VIDMEM)  \
    ( ((LEVELS)<<16) | (1<<15) | ((UNNORM)<<14) | ((LINEAR)<<13) | ((FORMAT)<<8) |          \
      ((DIM)<<4) | ((NOBORD)<<3) | ((CUBE)<<2) | ((MAINMEM)<<1) | ((VIDMEM)<<0) )

#define TEXTURE_CONTROL1(RCOMP, BGRA, RWRAP, NORMEXP, TWRAP, ANISO, SWRAP)                  \
    ( ((RCOMP)<<28) | ((BGRA)<<20) | ((RWRAP)<<16) | ((NORMEXP)<<12) | ((TWRAP)<<8) |       \
      ((ANISO)<<5) | ((SWRAP)<<0) )

#define TEXTURE_CONTROL2(ENABLE, MINLEV, MAXLEV, ANISO, ALPHATEST, COLORKEY)                \
    ( ((ENABLE)<<31) | ((MINLEV)<<27) | ((MAXLEV)<<15) | ((ANISO)<<4) | ((ALPHATEST)<<2) |  \
      ((COLORKEY)<<0) )

#define TEXTURE_SWIZZLE(Y16X16, CHANNELS, SWIZZLE)                                          \
    ( ((Y16X16)<<16) | ((CHANNELS)<<8) | ((SWIZZLE)<<0) )

#define TEXTURE_FILTER(SIGNS, MAG, MIN, KERN, LODBIAS)                                      \
    ( ((SIGNS)<<28) | ((MAG)<<24) | ((MIN)<<16) | ((KERN)<<13) | ((LODBIAS)<<0) )

#define TEXTURE_SIZE1(WIDTH, HEIGHT)                                                        \
    ( ((WIDTH)<<16) | ((HEIGHT)<<0) )

#define TEXTURE_SIZE2(DEPTH, ROWPITCH)                                                      \
    ( ((DEPTH)<<20) | ((ROWPITCH)<<0) )


#define PRINT Printf



#define R2O_WATER_OBJECT_FLAG_ACTIVE   0x01



enum AsmVer
{
  NO_ASM                                 = 0,

  REFINEMENT_COPY_MAJOR_EVEN_ASM         = 1,
  REFINEMENT_COPY_MINOR_EVEN_ASM         = 2,
  REFINEMENT_COPY_MAJOR_ODD_ASM          = 3,
  REFINEMENT_COPY_MINOR_ODD_ASM          = 4,
  INTERPOLATE_DD_MINUS_LINEAR_ASM        = 5,
  ADD_DD_ASM                             = 6,
  INTERP_LINEAR_ASM                      = 7,
  GEOMORPH_ASM                           = 8,
  REPLICATE_AMBIENT_ASM                  = 9,
  GENERATE_DERIVS_ASM                    = 10,
  CONVERT_HEIGHTS_TO_DERIVS_ASM          = 11,
  INTERP_TILE_BILINEAR_ASM               = 12,
  ROTATE_TILE_CLOCKWISE_ASM              = 13,
  ROTATE_TILE_ANTICLOCKWISE_ASM          = 14,
  COMPRESS_DERIVS_FOR_OUTPUT_ASM         = 15,
  COMPRESS_VERTS_FOR_OUTPUT_ASM          = 16,
  GENERATE_OUTCODES_ASM                  = 17,
  SET_NEAR_BITS_ASM                      = 18,
  GET_REFINEMENT_WINDOW_ASM              = 19,
  REFINE_OUTCODES_EVEN_ASM               = 20,
  INTERPOLATE_OUTCODES_HORIZONTAL_ASM    = 21,
  INTERPOLATE_OUTCODES_VERTICAL_EVEN_ASM = 22,
  INTERPOLATE_OUTCODES_VERTICAL_ODD_ASM  = 23,
  ASSIGN_INDICES_ASM                     = 24,
  FLAG_VERTS_RENDER_ASM                  = 25,
  FLAG_VERTS_KEEP_ASM                    = 26,
  FLAG_QUADS_RENDER_KEEP_ASM             = 27,
  ACCUMULATE_NORMAL_MAPS_ASM             = 28,

  ALL_ASM
};



// This struct will help determine the visibility of water objects
#define R2O_OCCL_QUERY_RESET   15
#define R2O_OCCL_QUERY_RECHECK 2
struct R2OOcclQuery
{
  u32 m_id;
  u16 m_visible;
  u16 m_time_stamp;
};



// per-viewport data
struct R2OViewData
{
  mtx4f       m_world_to_camera_matrix;
  mtx4f       m_frag_texproj_matrix;
  puspu_vec4  m_camera_position;
  puspu_vec4  m_frustum_planes[6];
  f32         m_frustum_near;
  f32         m_frustum_far;
  f32         m_frustum_hratio;
  f32         m_frustum_vratio;
  u8          m_camera_underwater;
  u8          m_vp_index;
  u8          m_any_visible;
}
ALIGNED(16);



// output render data for each water object
struct R2ORenderData
{
  puspu_vec4    m_origins     [MAX_LODS];
  u32           m_ea_verts    [MAX_LODS];
  u32           m_ea_derivs   [MAX_LODS];
  u32           m_ea_indices  [MAX_LODS];
  u32           m_ea_index_map[MAX_INT_LODS];
  u16           m_num_indices [MAX_LODS];
  u16           m_num_verts   [MAX_LODS];
  u16           m_cols_rows   [MAX_LODS];

  f32           m_camera_height_above_surface;
  u8            m_b_camera_over_water;
  u8            m_visible;
}
ALIGNED(16);


// built version of water object
// Be sure to keep R2OBuiltWaterObject in igCore/igHeaders/ps3strucst.h in sync!
struct R2OBuiltWaterObject
{
  puspu_vec4    m_origin;
  puspu_vec4    m_dimensions;
  puspu_vec4    m_water_color;      // linear

  f32           m_amplitude;
  u32           m_first_waveband;
  u32           m_material;
  u16           m_occl_id;
  u8            m_cubemap_index;
  u8            m_pad[1] ;
}
ALIGNED(16);


struct R2OBasicWaterObject
{
  puspu_vec4    m_origin;
  puspu_vec4    m_dimensions;
  puspu_vec4    m_water_color;    // linear

  u16           m_query_groups[MAX_QUERY_GROUPS_PER_OBJ];

  R2OOcclQuery  m_occl_query[R2O_MAX_VIEWPORTS];

  f32           m_amplitude;
  u32           m_first_waveband;
  u32           m_material;
  u16           m_occl_id;
  u16           m_num_registered_query_groups;
  u8            m_cubemap_index;
  u8            m_zone_index;
  u8            m_flags;
}
ALIGNED(16);


// runtime version of water object
struct R2OWaterObject : R2OBasicWaterObject
{
  // auxilliary structs
  R2OInteractiveData m_interactive_data;
  R2ORenderData      m_render_data[R2O_MAX_VIEWPORTS];
}
ALIGNED(16);




// for cutting holes!
struct R2OHole
{
  f32 m_xcoord;
  f32 m_zcoord;
  u32 m_tag;
  u8  m_lod;
  u8  m_dir;
  u8  m_cnt;
  u8  m_pad;
}
ALIGNED(16);



// for multiple spu support
enum R2OTaskType
{
  R2O_TASK_UPDATE_AMB,
  R2O_TASK_GENERATE_NORMS_AMB,
  R2O_TASK_PRE_UPDATE_INT,
  R2O_TASK_UPDATE_INT_PRE_OBJECT,
  R2O_TASK_UPDATE_INT,
  R2O_TASK_POST_UPDATE_INT,
  R2O_TASK_PRE_RENDER,
  R2O_TASK_RENDER,
};

struct R2OTask
{
  u32             m_type;
  u32             m_obj;
  u32             m_idx;    // lod for interactive update tasks, viewport index for render tasks
  u32             m_pad;
}
ALIGNED(16);




#ifdef R2O_PPU
namespace IGG
{
#endif

struct R2OCon
{
  #ifdef R2O_PPU
  void            Init(IGG::BasicAllocator* alloc, u32 num_water_objects);
  void            Reset();
  R2OWaterObject *WaterObjectCreate();
  R2OWaterObject *WaterObjectCreate(R2OBuiltWaterObject* p_built_obj);
  void            WaterObjectDelete(R2OWaterObject *p_obj);
  #endif

  R2OViewData m_view_data[R2O_MAX_VIEWPORTS];
  u32   m_num_viewports;

  u32   m_ea_water_objects;
  u64*  m_object_alloc_bot;             // bitfield allocator (one bit per object)
  u64   m_object_alloc_top[1];          // bitfield allocator (one bit per 64 objects)
  i32   m_object_last;                  // index of highest allocated object (-1 if none)

  // time measurements
  // please keep this member aligned, it is referenced from the spu
  u32   m_spu_ticks  ALIGNED(16);


  u32   m_ea_phase_shifts;
  u32   m_ea_normal_map_pointers;
  u32   m_ea_holes;
  u32   m_ea_caustics_data;
  u32   m_ea_caustics_envelope;
  u32   m_ea_query_groups;
  u32   m_allocator;
  f32   m_time;
  f32   m_time_step;
  f32   m_DD_coeff0;
  f32   m_near0;
  f32   m_step;
  f32   m_height_scale;
  i32   m_max_verts;
  u64   m_water_object_disable_bits;
  u8    m_any_visible;
  u8    m_any_transparent_visible;
  u8    m_dd_extra;
  u8    m_suppress_rows;
  u8    m_max_water_objects;
  u8    m_num_water_objects;
  u8    m_num_lods;
  u8    m_near_lod;
  u8    m_num_holes;
  u8    m_caustics_lod;

  // use this to switch off the whole bloody thing
  u8    m_r2o_enable;

  // ambient wave system
  u32   m_amb_ea_index_map;
  u32   m_amb_ea_spectrum;
  u32   m_amb_ea_heightmaps;
  u32   m_amb_ea_wave_scale_palette;
  u32   m_amb_ea_angle_palette;
  u32   m_amb_ea_frequency_palette;
  u32   m_amb_ea_phasor_palette;
  u8    m_amb_num_bands;
  u8    m_amb_fft_width;
  u8    m_amb_palette_cnt;

  // interactive wave system
  u32   m_int_ea_tiles;
  u32   m_int_ea_maps;
  u32   m_int_ea_normal_maps;
  u32   m_int_ea_impulses;
  u32   m_int_palette_cnt;
  u32   m_int_palette;
  u32   m_int_ea_frequency_palette;
  u32   m_int_ea_index_map;
  f32   m_int_half_life;
  f32   m_int_map_weight;
  f32   m_int_scale;
  f32   m_int_threshold_create;
  f32   m_int_threshold_destroy;
  f32   m_int_cull_distance;
  f32   m_int_cull_dist_ratio;
  f32   m_int_inner_radius;
  f32   m_int_inner_displacement;
  f32   m_int_outer_displacement;
  f32   m_int_attenuation_radius;
  f32   m_int_attenuation_min;
  f32   m_int_damping_coeff0;
  f32   m_int_damping_coeff1;
  f32   m_int_cull_dist_coeff0;
  f32   m_int_cull_dist_coeff1;
  u16   m_int_max_tiles;
  u16   m_int_ideal_num_tiles;
  u16   m_int_num_tiles;
  u16   m_int_free_tile;
  u16   m_int_consumer_mask;
  u16   m_int_asm_line;
  u8    m_int_num_bands;
  u8    m_int_fft_width;
  u8    m_int_map_width;
  u8    m_int_ovlp_width;
  i8    m_int_test_band;
  u8    m_int_enable;
  #if SANITY_CHECKS
  u8    m_int_error;
  #endif

  // debug variables
  f32   m_min_rad;
  f32   m_coarse_map_weight;
  f32   m_fine_map_weight;
  f32   m_frustum_fudge_factor;
  f32   m_index_of_refraction;
  f32   m_clip_plane_height;
  f32   m_tie_scale_adjust;
  f32   m_bsphere_min_angle_degrees;
  f32   m_cube_map_blend_ratio;
  f32   m_misc_consts[4];
  u32   m_timer_line;
  i16   m_debug_height;
  u8    m_paused;
  u8    m_wireframe;
  i8    m_debug_fp;
  i8    m_debug_band;

  // spu task-related
  u32   m_ea_task_list;
  u32   m_ea_norm;
  u32   m_num_tasks;
}
ALIGNED(16);

#ifdef R2O_PPU
} // namespace IGG
#endif



