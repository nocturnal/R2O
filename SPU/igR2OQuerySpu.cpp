#include "igR2OSpu.h"



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

f32 R2O_ApproxCameraHeightAboveSurface()
{
  return spu_extract(g_pViewData->m_camera_position - g_WaterObject.m_origin, 1);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_RefineCameraHeightAboveSurface()
{
  static u32 b_within_grid;

  if (lod==0)
  {
    b_within_grid = true;
  }
  else
  {
    if (!b_within_grid)
    {
      return;
    }
  }

  // Grid position (c,r) has world coords (origin_world + c*dvc_world + r*dvr_world).
  // When these coords equal the the camera position's x and z,
  // c and r must satisfy the following 2 equations:
  //
  // c*dvc_world.x + r*dvr_world.x = camera_pos.x - origin_world.x
  // c*dvc_world.z + r*dvr_world.z = camera_pos.z - origin_world.z
  //
  // A*c + B*r = E
  // C*c + D*r = F
  //
  // (A B)(c) = (E)
  // (C D)(r)   (F)
  //
  // (c) = ( D -B) (E) * 1/(AD-BC)
  // (r)   (-C  A) (F) 

  f32 A  = spu_extract(dvc_world, 0);
  f32 B  = spu_extract(dvr_world, 0);
  f32 C  = spu_extract(dvc_world, 2);
  f32 D  = spu_extract(dvr_world, 2);

  f32 recip_det = 1.0f / (A*D-B*C);

  vf32 d_pos = g_pViewData->m_camera_position - origin_world;
  f32 E = spu_extract(d_pos, 0);
  f32 F = spu_extract(d_pos, 2);

  f32 c_float = (D*E-B*F) * recip_det;
  f32 r_float = (A*F-C*E) * recip_det;

  u32 c_in_range = (c_float >= 4.0f) && (c_float < (f32)(nc-5));
  u32 r_in_range = (r_float >= 4.0f) && (r_float < (f32)(nr-5));

  if (!(c_in_range && r_in_range))
  {
    b_within_grid = false;
    return;
  }

  // get integer grid coords & fractional parts
  i32 c = (i32)c_float;
  i32 r = (i32)r_float;

  f32 c_frac = c_float - (f32)c;
  f32 r_frac = r_float - (f32)r;

  // get the heights of the 4 quad verts
  f32 y00 = g_Heightmap[ (c+0) + (r+0)*nc ];     // y00--y10
  f32 y10 = g_Heightmap[ (c+1) + (r+0)*nc ];     //  |    |
  f32 y01 = g_Heightmap[ (c+0) + (r+1)*nc ];     //  |    |
  f32 y11 = g_Heightmap[ (c+1) + (r+1)*nc ];     // y01--y11
  f32 y;

  // does the quad have a major or a minor diagonal?
  if (((c^r) & 1) == 0)
  {
    // major diagonal y00->y11
    if (c_float > r_float)
    {
      // use upper-right triangle
      y = y10 + (1.0f-c_frac) * (y00-y10) + r_frac * (y11-y10);
      //= (1-cf, cf-rf, 0, rf).(y00, y10, y01, y11)
      // (cf', cf-rf, 0, rf)
    }
    else
    {
      // use lower-left triangle
      y = y01 + c_frac * (y11-y01) + (1.0f-r_frac) * (y00-y01);
      //= (1-rf, 0, rf-cf, cf).(y00, y10, y01, y11)
      // (rf', 0, rf-cf, cf)
    }
  }
  else
  {
    // minor diagonal: (c+1,r), (c,r+1)
    if (c_float + r_float < 1.0f)
    {
      // use upper-left triangle
      y = y00 + c_frac * (y10-y00) + r_frac * (y01-y00);
      //= (1-cf-rf, cf, rf, 0).(y00, y10, y01, y11)
      // (cf'-rf, cf, rf, 0)
    }
    else
    {
      // use lower-right triangle
      y = y11 + (1.0-c_frac) * (y01-y11) + (1.0-r_frac) * (y10-y11);
      //= (0, 1-rf, 1-cf, cf+rf-1).(y00, y10, y01, y11)
      //= (0, rf', cf', rf-cf')
    }
  }

  f32 camy = spu_extract(g_pViewData->m_camera_position, 1);
  f32 orgy = spu_extract(g_WaterObject.m_origin, 1);
  g_RenderData.m_camera_height_above_surface = camy - (orgy + y);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_ProcessQueries()
{
  // set up some constants
  f32 A  = spu_extract(dvc_world, 0);
  f32 B  = spu_extract(dvr_world, 0);
  f32 C  = spu_extract(dvc_world, 2);
  f32 D  = spu_extract(dvr_world, 2);

  f32 recip_det = 1.0f / (A*D-B*C);

  f32 orgx = spu_extract(origin_world, 0);
  f32 orgy = spu_extract(origin_world, 1);
  f32 orgz = spu_extract(origin_world, 2);


  // loop over registered query groups
  for (u32 g=0; g<g_WaterObject.m_num_registered_query_groups; g++)
  {
    // get the query group index
    u32 idx = g_WaterObject.m_query_groups[g];

    // get the query group struct
    R2OQueryGroup group;
    u32 ea = g_R2OCon.m_ea_query_groups;
    R2O_GetWait(&group, ea, idx, sizeof(R2OQueryGroup));

    // loop over queries in the group
    for (u32 q=0; q<group.m_num_queries; q++)
    {
      // get the query
      R2OQuery query;
      if (!group.m_indirect)
      {
        ea = (u32)group.m_queries + q * sizeof(R2OQuery);
      }
      else
      {
        vu32 addr;
        ea = (u32)group.m_queries + q * sizeof(u32);
        u32 ls = (u32)&addr + (ea & 15);
        R2O_GetWait((void *)ls, ea, 0, sizeof(u32));
        ea = *(u32 *)ls;
      }
      R2O_GetWait(&query, ea, 0, sizeof(R2OQuery));

      // quit if the current level of refinement doesn't enclose the query point
      if (lod==0)
      {
        query.m_flags |= QUERY_FLAG_WITHIN_GRID;
      }
      else
      {
        if (!(query.m_flags & QUERY_FLAG_WITHIN_GRID))
        {
          continue;
        }
      }

      // don't bother recalculating this lod if another viewport has already done it
      if (lod <= query.m_lod)
      {
        continue;
      }

      f32 E = query.m_x - orgx;
      f32 F = query.m_z - orgz;

      f32 c_float = (D*E-B*F) * recip_det;
      f32 r_float = (A*F-C*E) * recip_det;

      u32 c_in_range = (c_float >= 4.0f) && (c_float < (f32)(nc-5));
      u32 r_in_range = (r_float >= 4.0f) && (r_float < (f32)(nr-5));

      if (!(c_in_range && r_in_range))
      {
        // we can't refine this query point any further
        query.m_flags &= ~QUERY_FLAG_WITHIN_GRID;
        continue;
      }

      // get integer grid coords & fractional parts
      i32 c = (i32)c_float;
      i32 r = (i32)r_float;

      f32 c_frac = c_float - (f32)c;
      f32 r_frac = r_float - (f32)r;

      // get the heights of the 4 quad verts
      f32 y00 = g_Heightmap[ (c+0) + (r+0)*nc ];     // y00--y10
      f32 y10 = g_Heightmap[ (c+1) + (r+0)*nc ];     //  |    |
      f32 y01 = g_Heightmap[ (c+0) + (r+1)*nc ];     //  |    |
      f32 y11 = g_Heightmap[ (c+1) + (r+1)*nc ];     // y01--y11
      f32 y;

      // does the quad have a major or a minor diagonal?
      if (((c^r) & 1) == 0)
      {
        // major diagonal y00->y11
        if (c_float > r_float)
        {
          // use upper-right triangle
          y = y10 + (1.0f-c_frac) * (y00-y10) + r_frac * (y11-y10);
          //= (1-cf, cf-rf, 0, rf).(y00, y10, y01, y11)
          // (cf', cf-rf, 0, rf)
        }
        else
        {
          // use lower-left triangle
          y = y01 + c_frac * (y11-y01) + (1.0f-r_frac) * (y00-y01);
          //= (1-rf, 0, rf-cf, cf).(y00, y10, y01, y11)
          // (rf', 0, rf-cf, cf)
        }
      }
      else
      {
        // minor diagonal: (c+1,r), (c,r+1)
        if (c_float + r_float < 1.0f)
        {
          // use upper-left triangle
          y = y00 + c_frac * (y10-y00) + r_frac * (y01-y00);
          //= (1-cf-rf, cf, rf, 0).(y00, y10, y01, y11)
          // (cf'-rf, cf, rf, 0)
        }
        else
        {
          // use lower-right triangle
          y = y11 + (1.0-c_frac) * (y01-y11) + (1.0-r_frac) * (y10-y11);
          //= (0, 1-rf, 1-cf, cf+rf-1).(y00, y10, y01, y11)
          //= (0, rf', cf', rf-cf')
        }
      }

      query.m_y = orgy + y;

      // put the query back for further refinement
      R2O_PutWait(&query, ea, 0, sizeof(R2OQuery));
    }
  }
}

