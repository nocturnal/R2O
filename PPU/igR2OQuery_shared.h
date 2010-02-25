#pragma once


#define MAX_QUERY_GROUPS         256
#define MAX_QUERY_GROUPS_PER_OBJ 64
#define CAMERAS_USE_QUERY_SYSTEM 0

struct R2OWaterObject;


struct R2OQuery
{
  f32 m_x;
  f32 m_y;
  f32 m_z;

  // hands off!
  i8  m_lod;
  u8  m_flags;
  u8  m_pad[2];
}
ALIGNED(16);



struct R2OQueryGroup
{
  R2OWaterObject *m_water_object;
  void *m_queries;    // either R2OQuery * or R2OQuery **
  u16   m_num_queries;
  u16   m_indirect;
  u16   m_pad[2];
}
ALIGNED(16);

