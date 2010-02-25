#include "igR2O.h"


namespace IGG
{


R2OQueryGroup g_QueryGroups[MAX_QUERY_GROUPS];
u32 g_NumRegisteredQueryGroups;

#if CAMERAS_USE_QUERY_SYSTEM
// camera queries
R2OQuery       g_CameraQueries    [R2O_MAX_VIEWPORTS];
R2OQueryGroup *g_CameraQueryGroups[R2O_MAX_VIEWPORTS];
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Per-frame init
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_ResetQueryGroups()
{
  g_NumRegisteredQueryGroups = 0;

  R2OWaterObject *p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (u32 i=0; i<g_R2OCon.m_max_water_objects; i++,p_obj++)
  {
    p_obj->m_num_registered_query_groups = 0;
  }

  #if CAMERAS_USE_QUERY_SYSTEM
  // register a 1-entry query group for each camera
  for (u32 i_vp=0; i_vp<g_R2OCon.m_num_viewports; i_vp++)
  {
    vec4f cam_pos = g_R2OCon.m_view_data[i_vp].m_camera_position;
    *(vec4f *)&g_CameraQueries[i_vp] = cam_pos;
    R2OWaterObject *p_closest = R2O_ClosestWaterObjectAboveOrBelow(cam_pos);
    g_CameraQueryGroups[i_vp] = R2O_RegisterQueries(&g_CameraQueries[i_vp], 1, false, p_closest);
  }
  #endif
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Register a group of queries or query pointers
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

R2OQueryGroup *R2O_RegisterQueries(void *queries, u32 num_queries, u32 indirect, R2OWaterObject *p_obj)
{
  // return NULL if any input is NULL
  if (!queries || !num_queries || !p_obj)
  {
    return NULL;
  }

  // make sure we're in good standing
  IG_ASSERT(g_NumRegisteredQueryGroups <= MAX_QUERY_GROUPS);
  IG_ASSERT(p_obj->m_num_registered_query_groups <= MAX_QUERY_GROUPS_PER_OBJ);

  // set initial approximation (i.e. sealevel) for each result
  f32 sealevel = p_obj->m_origin.y;
  if (indirect)
  {
    for (u32 i=0; i<num_queries; i++)
    {
      ((R2OQuery **)queries)[i]->m_y   = sealevel;
      ((R2OQuery **)queries)[i]->m_lod = -1;
    }
  }
  else
  {
    for (u32 i=0; i<num_queries; i++)
    {
      ((R2OQuery *)queries)[i].m_y   = sealevel;
      ((R2OQuery *)queries)[i].m_lod = -1;
    }
  }

  // check there's room to register the group
  if (g_NumRegisteredQueryGroups == MAX_QUERY_GROUPS ||
      p_obj->m_num_registered_query_groups == MAX_QUERY_GROUPS_PER_OBJ)
  {
    // no room - the user will only get the sealevel approximation
    return NULL;
  }

  // use next available query group struct
  u16 global_iq = g_NumRegisteredQueryGroups++;
  u16 object_iq = p_obj->m_num_registered_query_groups++;
  R2OQueryGroup *p_group = &g_QueryGroups[global_iq];
  p_obj->m_query_groups[object_iq] = global_iq;

  p_group->m_water_object = p_obj;
  p_group->m_queries      = queries;
  p_group->m_num_queries  = num_queries;
  p_group->m_indirect     = indirect;

  return p_group;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// return pointer to the closest water object to the given point, or NULL if there aren't any water objects loaded
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

R2OWaterObject *R2O_ClosestWaterObject(vec4f &point)
{
  R2OWaterObject *p_closest = NULL;
  f32 min_sq_dist = 1e20f;

  // loop over water objects looking for closest water object
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

    // check whether it's closer to the point
    vec4f clamped, diff;
    VecClamp3(clamped, point, obj_min, obj_max);
    SubVector3(diff, point, clamped);
    f32 sq_dist = diff.LenSquared3();
    if (sq_dist < min_sq_dist)
    {
      // record new closest object
      p_closest = p_obj;
      min_sq_dist  = sq_dist;
    }
  }

  return p_closest;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// return pointer to the closest water object to the given point, or NULL if the point isn't over water (or under water)
// (also returns NULL if there aren't any water objects loaded)
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

R2OWaterObject *R2O_ClosestWaterObjectAboveOrBelow(vec4f &point)
{
  R2OWaterObject *p_closest = NULL;
  f32 min_abs_dy = 1e20f;

  // loop over water objects looking for closest water surface directly above or below point
  R2OWaterObject *p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++, p_obj++)
  {
    if (!(p_obj->m_flags & R2O_WATER_OBJECT_FLAG_ACTIVE))
    {
      continue;
    }

    // get object extents
    f32 x_min = p_obj->m_origin.x;
    f32 z_min = p_obj->m_origin.z;
    f32 x_max = x_min + p_obj->m_dimensions.x;
    f32 z_max = z_min + p_obj->m_dimensions.z;

    // skip this object if the camera isn't directly above or below it
    if (point.x<x_min || point.x>x_max || point.z<z_min || point.z>z_max)
    {
      continue;
    }

    // check whether it's the closest surface to the point
    f32 dy = point.y - p_obj->m_origin.y;
    f32 abs_dy = fabsf(dy);
    if (abs_dy < min_abs_dy)
    {
      // record new closest object
      p_closest = p_obj;
      min_abs_dy  = abs_dy;
    }
  }

  return p_closest;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_CameraUnderwaterTest()
{
  // test each camera
  for (u32 i_vp=0; i_vp<g_R2OCon.m_num_viewports; i_vp++)
  {
    R2OViewData *p_view_data  = &g_R2OCon.m_view_data[i_vp];
    vec4f cam_pos = p_view_data->m_camera_position;
    u32 b_underwater = false;

    #if CAMERAS_USE_QUERY_SYSTEM
    R2OQueryGroup *p_group = g_CameraQueryGroups[i_vp];
    if (p_group)
    {
      vec4f cam_pos = g_R2OCon.m_view_data[i_vp].m_camera_position;
      f32 dy = cam_pos.y - g_CameraQueries[i_vp].m_y;
      b_underwater = ((dy < 0.0f) && (dy > -g_MaxWaterDepth));
    }
    #else
    R2OWaterObject *p_obj = R2O_ClosestWaterObjectAboveOrBelow(cam_pos);
    if (p_obj && p_obj->m_render_data[i_vp].m_b_camera_over_water)
    {
      f32 dy = p_obj->m_render_data[i_vp].m_camera_height_above_surface;
      b_underwater = ((dy < 0.0f) && (dy > -g_MaxWaterDepth));
    }
    #endif
    
    // activate Mark's underwater settings
    UnderWaterSetState( b_underwater ?
                        (UnderWaterGetState() |  (ENABLE_UNDERWATER_FOG | ENABLE_UNDERWATER_CAUSTICS)) :
                        (UnderWaterGetState() & ~(ENABLE_UNDERWATER_FOG | ENABLE_UNDERWATER_CAUSTICS)) );

    p_view_data->m_camera_underwater = b_underwater;
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Eric's callback
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void (*g_R2OHeightQueryCallback)(R2OWaterObject *) = NULL;

void R2O_SetHeightQueryCallback(void (*function)(R2OWaterObject *))
{
  g_R2OHeightQueryCallback = function;
}


}; // namespace IGG

