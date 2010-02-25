#include "igR2OSpu.h"


// globals

vf32 g_Planes[6];
vf32 origin_world, dvc_world, dvr_world;
vf32 origin_camera, dvc_camera, dvr_camera;
vf32 basis_col, basis_row;
vf32 clip_min, clip_max;
i32 nc, nr, true_nc, true_nr;
i32 lod;
i32 max_points;
f32 step, prev_step, inv_step;

f32 *g_Heightmap, *g_TempBuf;
u8 *g_Outcodes, *g_OutcodesAlt;
u16  *g_Indices;

vf32 even_basis_col, odd_basis_col;
vf32 even_basis_row, odd_basis_row;
f32  even_step, odd_step;
f32  even_inv_step, odd_inv_step;

R2OHole *g_Holes;

R2OViewData *g_pViewData;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_PreRender()
{
  // assume nothing is visibile
  g_R2OCon.m_any_visible = false;
  for (u32 i_vp=0; i_vp<g_R2OCon.m_num_viewports; i_vp++)
  {
    g_R2OCon.m_view_data[i_vp].m_any_visible = false;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// task-base render
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_RenderObject(u32 i_obj, u32 i_vp)
{
  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // get the water object and interactive data
  u32 ea  = g_R2OCon.m_ea_water_objects;
  u32 ofs = OFFSETOF(R2OWaterObject, m_interactive_data);
  R2O_GetWait(&g_WaterObject, ea, i_obj, sizeof(R2OBasicWaterObject), sizeof(R2OWaterObject));
  R2O_GetWait(&g_InteractiveData, ea+ofs, i_obj, sizeof(R2OInteractiveData), sizeof(R2OWaterObject));

  // render to viewport (if visible)
  R2OViewData *p_view_data = &g_R2OCon.m_view_data[0];
  ofs      = OFFSETOF(R2OWaterObject, m_render_data);
  u32 ofs1 = OFFSETOF(R2ORenderData, m_visible);
  R2O_GetWait(&g_RenderData.m_visible, ea+ofs+ofs1, i_obj, sizeof(g_RenderData.m_visible), sizeof(R2OWaterObject));
  if (g_RenderData.m_visible)
  {
    // render water object
    R2O_Render(p_view_data);

    // use the render data to generate the interactive index maps
    R2O_GenerateIndexMaps();

    // pass back the output for this object
    R2O_PutWait(&g_RenderData, ea+ofs, i_obj, sizeof(R2ORenderData), sizeof(R2OWaterObject));

    // accumulate visibility result
    g_R2OCon.m_view_data[i_vp].m_any_visible |= g_RenderData.m_visible;
    g_R2OCon.m_any_visible |= g_RenderData.m_visible;
  }

  // reset the scratch memory
  g_Scratchpad.Reset(state);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FlipBuffers()
{
  // flip outcode buffer base pointers
  u8 *temp_o = g_Outcodes;
  g_Outcodes = g_OutcodesAlt;
  g_OutcodesAlt = temp_o;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SetBasisEtc(i32 c_orig, i32 r_orig)
{
  if (lod & 1)
  {
    step      = odd_step;
    inv_step  = odd_inv_step;
    basis_col = odd_basis_col;
    basis_row = odd_basis_row;
  }
  else
  {
    step      = even_step;
    inv_step  = even_inv_step;
    basis_col = even_basis_col;
    basis_row = even_basis_row;
  }

  origin_world += spu_splats((f32)c_orig)*dvc_world + spu_splats((f32)r_orig)*dvr_world;

  if (lod & 1)
  {
    f32 s = 1.0f * even_step;
    dvc_world = (vf32){ s,0,s,0};
    dvr_world = (vf32){-s,0,s,0};
  }
  else
  {
    f32 s = even_step;
    dvc_world = (vf32){s,0,0,0};
    dvr_world = (vf32){0,0,s,0};
  }

  origin_camera = MatMulVec((const mtx4 &)g_pViewData->m_world_to_camera_matrix, origin_world);
  dvc_camera    = MatMulVec((const mtx4 &)g_pViewData->m_world_to_camera_matrix, dvc_world);
  dvr_camera    = MatMulVec((const mtx4 &)g_pViewData->m_world_to_camera_matrix, dvr_world);
}


void NextLod(i32 c_orig, i32 r_orig)
{
  if (lod & 1)
  {
    odd_step     *= 0.5f;
    odd_inv_step *= 2.0f;
  }
  else
  {
    even_step     *= 0.5f;
    even_inv_step *= 2.0f;
  }

  lod++;
  SetBasisEtc(c_orig, r_orig);
}


void InitBasisEtc()
{
  // Use a fixed initial step size for now; 128m for lod 0.
  // This yields an fft tile size of 32 x 128m = 4096m for lod 0
  // and maximum dimensions of 8192m x 8192m
  step = g_R2OCon.m_step;
  vf32 step_vec = (vf32){step, 0, step, 0};
  
  // get inverse-step using float magic (since taking the reciprocal of a power of 2 yields a 1-bit error)
  qword q_step  = si_from_float(step);
  qword q_magic = si_ilhu(0x7F00);
  inv_step = si_to_float(si_sf(q_step, q_magic));

  // set clip window
  clip_min = g_WaterObject.m_origin;
  clip_max = g_WaterObject.m_origin + g_WaterObject.m_dimensions;

  // set origin at gridpoint below clip min
  f32 magic_float = 1.5f * 8388608.0f * step;
  vf32 magic_vf32 = (vf32){magic_float, 0, magic_float, 0};
  origin_world = (clip_min + magic_vf32) - magic_vf32;

  // compute gridpoint above clip max
  vf32 max_corner = (clip_max + magic_vf32) - magic_vf32;
  max_corner += step_vec;

  // offset both corners by the necessary amount of padding
  origin_world -= step_vec * spu_splats(8.0f);
  max_corner   += step_vec * spu_splats(8.0f);

  // set num cols & num rows
  vf32 dims = max_corner - origin_world;
  nc = (i32)(spu_extract(dims,0) * inv_step) + 1;
  nr = (i32)(spu_extract(dims,2) * inv_step) + 1;

  // record true nc, nr
  true_nc = nc - 16;
  true_nr = nr - 16;

  // alignment requirements (ooh, that's a bit strict)
  nc = (nc + 7) & -8;
  nr = (nr + 7) & -8;

  // deal with large grids
  if (nc > 80)
  {
    nc = 80;
    true_nc = 64;
    dims = spu_insert((nc-1)*step, dims, 0);
  }
  if (nr > 80)
  {
    nr = 80;
    true_nr = 64;
    dims = spu_insert((nr-1)*step, dims, 2);
  }
  max_corner = origin_world + dims;


  even_step = step;
  even_inv_step = inv_step;
  even_basis_col = (vf32){1.0f, 0.0f, 0.0f, 0.0f};
  even_basis_row = (vf32){0.0f, 0.0f, 1.0f, 0.0f};

  const f32 r = 0.707106781187f;
  odd_step    = even_step * r;
  odd_inv_step= even_inv_step * r * 2.0f;
  odd_basis_col = (vf32){ r, 0.0f, r, 0.0f};
  odd_basis_row = (vf32){-r, 0.0f, r, 0.0f};

  basis_col = even_basis_col;
  basis_row = even_basis_row;
  dvc_world = spu_splats(step) * basis_col;
  dvr_world = spu_splats(step) * basis_row;

  // set base lod origin
  g_RenderData.m_origins[0]   = origin_world;
  g_RenderData.m_cols_rows[0] = nc<<8 | nr;
  c0_amb = 0;
  r0_amb = 0;

  SetBasisEtc(0,0);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GenerateFrustum()
{
  // g_R2OCon.m_frustum_planes has 6 planes in the following order:
  //
  // 0: near  (there is igg code that assumes this is first)
  // 1: left
  // 2: right
  // 3: bottom
  // 4: top
  // 5: far   (there is igg code that assumes this is last)

  // the plane normals point into the interior of the frustum

  // copy T,B,L,R planes directly
  for (u32 i=1; i<5; i++)
  {
    g_Planes[i] = g_pViewData->m_frustum_planes[i];

    // loosen frustum based on lod, amplitude, and a fudgefactor
    f32 nw = spu_extract(g_Planes[i], 3);
    f32 ny = spu_extract(g_Planes[i], 1);

    //f32 slop = g_R2OCon.m_frustum_fudge_factor * g_WaterObject.m_amplitude * powf(0.5f, 0.5f*(f32)lod)
    //            * (g_R2OCon.m_step * 0.0625f);    // readjust based on changed in lod 0 from 16m to 128m stepsize
    //f32 slop = g_R2OCon.m_frustum_fudge_factor * powf(0.5f, 0.5f*(f32)lod);
    f32 slop = g_R2OCon.m_frustum_fudge_factor * step;

    nw += slop * fabsf(ny);
    g_Planes[i] = spu_insert(nw, g_Planes[i], 3);
  }

  // compute near subfrustum plane based on lod
  //f32 d = g_R2OCon.m_near0 * powf(0.5f, 0.5f*(f32)lod);
  f32 d = g_R2OCon.m_near0 * step * (1.0f/128.0f);

  // test
  if (lod==g_R2OCon.m_num_lods-1)
  {
    d = step;
  }

  vf32 camera_direction = (vf32){0,0,1,0};
  camera_direction = spu_insert(g_pViewData->m_world_to_camera_matrix.m_v0.z, camera_direction, 0);
  camera_direction = spu_insert(g_pViewData->m_world_to_camera_matrix.m_v1.z, camera_direction, 1);
  camera_direction = spu_insert(g_pViewData->m_world_to_camera_matrix.m_v2.z, camera_direction, 2);

  f32 dot = spu_extract(spu_dot3(camera_direction, g_pViewData->m_camera_position), 0);
  f32 w = -(dot + d);

  g_Planes[0] = spu_insert(w, camera_direction, 3);

  // reverse the sense of the near plane
  g_Planes[0] = -g_Planes[0];
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ClipToRectangle, for even lods
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ClipToRectangle(vf32 clip_min, vf32 clip_max)
{
  // convert from world coords to integer rows & cols
  vf32 norm_min = spu_mul(clip_min - origin_world, spu_splats(inv_step));
  vf32 norm_max = spu_mul(clip_max - origin_world, spu_splats(inv_step));
  vi32 int_min  = VecFloor4(norm_min);
  vi32 int_max  = VecCeil4 (norm_max);

  // expand rectangle by 1 gridpoint because quads incident to the verts we're about to cull out will also be culled out
  i32 c_min = spu_extract(int_min, 0) - 1;
  i32 c_max = spu_extract(int_max, 0) + 1;
  i32 r_min = spu_extract(int_min, 2) - 1;
  i32 r_max = spu_extract(int_max, 2) + 1;

  // trim loop bounds to rectangle so we don't splat memory
  i32 c0 = c_min >= 0 ? c_min   : 0;
  i32 c1 = c_max < nc ? c_max+1 : nc;
  i32 r0 = r_min >= 0 ? r_min   : 0;
  i32 r1 = r_max < nr ? r_max+1 : nr;


  // cull left points
  if (c_min>=0 && c_min<nc)
  {
    u8 *p = &g_Outcodes[r0*nc+c_min];
    for (i32 r=r0; r<r1; r++,p+=nc)
    {
      *p |= 0x80;
    }
  }

  // cull right points
  if (c_max>=0 && c_max<nc)
  {
    u8 *p = &g_Outcodes[r0*nc+c_max];
    for (i32 r=r0; r<r1; r++,p+=nc)
    {
      *p |= 0x80;
    }
  }

  // cull upper points
  if (r_min>=0 && r_min<nr)
  {
    u8 *p = &g_Outcodes[r_min*nc+c0];
    for (i32 c=c0; c<c1; c++,p++)
    {
      *p |= 0x80;
    }
  }

  // cull lower points
  if (r_max>=0 && r_max<nr)
  {
    u8 *p = &g_Outcodes[r_max*nc+c0];
    for (i32 c=c0; c<c1; c++,p++)
    {
      *p |= 0x80;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_CutHoles(u8 *outcodes, i32 width, i32 height, u32 lod)
{
  // get window extents
  f32 x_min = spu_extract(origin_world, 0);
  f32 z_min = spu_extract(origin_world, 2);
  f32 x_max = x_min + (f32)width  * step;
  f32 z_max = z_min + (f32)height * step;

  // make sure we cut the right water object
  vu32 coords_u32 = *(vu32 *)&g_WaterObject.m_origin;
  u32 tag = spu_extract(coords_u32, 0) ^ spu_extract(coords_u32, 2);

  R2OHole *p_hole = g_Holes;
  for (u32 i=0; i<g_R2OCon.m_num_holes; i++, p_hole++)
  {
    if (tag==p_hole->m_tag && lod==p_hole->m_lod)
    {
      // set start of hole
      f32 x = p_hole->m_xcoord;
      f32 z = p_hole->m_zcoord;

      // set coords deltas
      f32 dx = p_hole->m_dir ? 0.0f : 2.0f*step;
      f32 dz = p_hole->m_dir ? 2.0f*step : 0.0f;

      // loop over length of hole
      for (u32 l=0; l<p_hole->m_cnt; l++)
      {
        // test for overlap with the current window
        if (x>=x_min && x<x_max && z>=z_min && z<z_max)
        {
          // translate coords to col/row values
          i32 c = (i32)((x-x_min) * inv_step);
          i32 r = (i32)((z-z_min) * inv_step);

          // punch out the hole
          outcodes[r*width + c] |= 0x80;

          // test camera against hole
          f32 x0 = x - step;
          f32 x1 = x + step;
          f32 z0 = z - step;
          f32 z1 = z + step;
          f32 x_cam = spu_extract(g_pViewData->m_camera_position, 0);
          f32 z_cam = spu_extract(g_pViewData->m_camera_position, 2);
          if (x_cam>=x0 && x_cam<x1 && z_cam>=z0 && z_cam<z1)
          {
            g_RenderData.m_b_camera_over_water = false;
          }
        }

        // step coords
        x += dx;
        z += dz;
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

u32 GetRefinementWindow(i32 &c_orig, i32 &r_orig, i32 &nc_dest, i32 &nr_dest, u32 &oflw,
                        f32 heights[], u8 outcodes[], i32 nc_src, i32 nr_src)
{
  oflw = false;

  // get window extents, assuming an initial origin coincident with the top-left corner of the source window
  i16 min_c_dest =  10000;
  i16 max_c_dest = -10000;
  i16 min_r_dest =  10000;
  i16 max_r_dest = -10000;

  GetRefinementWindowAsm(&min_c_dest, &max_c_dest, &min_r_dest, &max_r_dest, g_Outcodes, nc_src, nr_src, lod);

  // test for empty window
  if ((min_c_dest > max_c_dest) || (min_r_dest > max_r_dest))
  {
    return false;
  }

  // include bottom and right boundaries
  max_c_dest++;
  max_r_dest++;

  // extra rows & columns
  i16 extra = g_R2OCon.m_dd_extra;  //3 * DD_TAPS - 4;
  min_c_dest -= extra;
  min_r_dest -= extra;
  max_c_dest += extra;
  max_r_dest += extra;

  // parity constraints to ensure the dest window corners all lie on the degree-8 source verts
  min_c_dest &= -8; // should really only be a multiple of 4, but OutputFans can only handle multiples of 8 currently
  min_r_dest &= -4;
  max_c_dest = (max_c_dest+4) & -8; // ditto
  max_r_dest = (max_r_dest+3) & -4;

  // set origin and dimensions
  nc_dest = max_c_dest - min_c_dest;
  nr_dest = max_r_dest - min_r_dest;
  if (lod & 1)
  {
    c_orig = (min_r_dest + min_c_dest) >> 1;
    r_orig = (min_r_dest - min_c_dest) >> 1;
  }
  else
  {
    c_orig = (min_c_dest - min_r_dest) >> 1;
    r_orig = (min_c_dest + min_r_dest) >> 1;
  }

  // overflow check
  while (nc_dest * nr_dest > g_R2OCon.m_max_verts)
  {
    //Printf("refinement window overflow at t=%g!\n seconds", g_R2OCon.m_time);
    oflw = true;
    nr_dest -= 4;
  }

  // non-empty window
  return true;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SuppressRows(u8 outcodes[], u32 nc, u32 nr, u32 pad)
{
  u32 start = 0;
  u32 end   = nc * pad;
  for (u32 i=start; i<end; i++)
  {
    outcodes[i] = 0x80;
  }

  start = nc * (nr-pad);
  end   = nc * nr;
  for (u32 i=start; i<end; i++)
  {
    outcodes[i] = 0x80;
  }

  for (u32 r=pad; r<nr-pad; r++)
  {
    for (u32 c=0; c<pad; c++)
    {
      u32 i = r*nc + c;
      outcodes[i] = 0x80;
    }

    for (u32 c=nc-pad; c<nc; c++)
    {
      u32 i = r*nc + c;
      outcodes[i] = 0x80;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OutputMesh(f32 *heightmap, i16 *derivs)
{
  IG_ASSERT(lod < g_R2OCon.m_num_lods);
  
  // split g_TempBuf into 2 pieces for double-buffered dma
  u8 *buf0 = (u8 *)g_TempBuf;             // derivs currently reside in this half
  u8 *buf1 = (u8 *)g_TempBuf + (g_R2OCon.m_max_verts<<2);

  // mark quads with render/keep bits
  R2O_FlagQuadsRenderKeep(g_Outcodes, nc, nr);

  // flag verts with render bit if they're incident to a quad marked 'render'
  R2O_FlagVertsRender(g_Outcodes, nc, nr);

  // assign indices
  u32 num_verts = R2O_AssignIndices(g_Indices, g_Outcodes, nc*nr);

  // flag verts with cull bit if they're not incident to a quad marked 'keep'
  R2O_FlagVertsKeep(g_Outcodes, nc, nr);

  //-----------------------------------------------------------------
  // derivs
  CompressDerivsForOutputAsm(derivs, derivs, g_Indices, nc*nr);

  // size is used for both derivs and verts (each 4-bytes per, and in 1:1 correspondence)
  u32 size = num_verts * 2 * sizeof(i16);
  size = (size + 15) & -16;

  // allocate a buffer and dma into it
  u32 ea = IGG::DynamicRenderAlloc(g_R2OCon.m_allocator, size, 128, false);
  R2ODmaPutBig((volatile void *)derivs, ea, size, 0);

  // pass back the ea
  g_RenderData.m_ea_derivs[lod] = ea;


  //-----------------------------------------------------------------
  // verts
  u8 *verts_out = buf1;
  f32 height_scale = g_R2OCon.m_height_scale * (0.2f / g_WaterObject.m_amplitude);
  CompressVertsForOutputAsm(verts_out, heightmap, height_scale, g_Indices, nc, nr);

  // allocate a buffer and dma into it
  ea = IGG::DynamicRenderAlloc(g_R2OCon.m_allocator, size, 128, false);
  R2ODmaPutBig((volatile void *)verts_out, ea, size, 1);
  DmaWaitAll(1<<0 | 1<<1);

  // pass back the ea and count
  g_RenderData.m_ea_verts[lod]  = ea;
  g_RenderData.m_num_verts[lod] = num_verts;

  //-----------------------------------------------------------------
  // generate first half of indices
  u32 nr0 = ((nr>>1) + 3) & -4;
  DmaWaitAll(1<<0);
  u32 num_indices0 = R2O_GenerateFans((u16 *)buf0, g_Indices, nc, nr0, nc);
  
  // generate second half of indices
  u32 nr1 = nr - nr0;
  DmaWaitAll(1<<1);
  u32 num_indices1 = R2O_GenerateFans((u16 *)buf1, g_Indices+(nr0*nc), nc, nr1, nc);

  // allocate a buffer and dma into it
  u32 size0 = num_indices0 * sizeof(u16);
  u32 size1 = num_indices1 * sizeof(u16);
  ea = IGG::DynamicRenderAlloc(g_R2OCon.m_allocator, size0+size1, 128, false);
  R2ODmaPutBig((volatile void *)buf0, ea,       size0, g_kDmaTag);
  R2ODmaPutBig((volatile void *)buf1, ea+size0, size1, g_kDmaTag);

  // pass back the ea and count
  u32 num_indices = num_indices0 + num_indices1;
  g_RenderData.m_ea_indices[lod]  = ea;
  g_RenderData.m_num_indices[lod] = num_indices;

  // make sure we don't try to output indices if there weren't any verts to be indexed
  if (num_verts==0 && num_indices!=0)
  {
    Printf("no verts to index!\n");
  }

  if (num_verts && num_indices)
  {
    //g_RenderData.m_b_any_output = true;
    g_RenderData.m_visible = true;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_Render(R2OViewData *p_view_data)
{
  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // copy arg
  g_pViewData = p_view_data;

  // allocate buffers
  g_Heightmap   = (f32 *)g_Scratchpad.Alloc(g_R2OCon.m_max_verts << 2);
  g_TempBuf     = (f32 *)g_Scratchpad.Alloc(g_R2OCon.m_max_verts << 3);  // 64K!!! (the fans might need this much)
  g_Outcodes    = (u8  *)g_Scratchpad.Alloc(g_R2OCon.m_max_verts << 0);
  g_OutcodesAlt = (u8  *)g_Scratchpad.Alloc(g_R2OCon.m_max_verts << 0);
  g_Indices     = (u16 *)g_Scratchpad.Alloc(g_R2OCon.m_max_verts << 1);

  // get the holes for this level
  u32 num_holes = g_R2OCon.m_num_holes;
  if (num_holes)
  {
    g_Holes = (R2OHole *)R2O_AllocGetWait(g_R2OCon.m_ea_holes, 0, num_holes * sizeof(R2OHole));
  }

  lod = 0;
  i32 c_orig=0, r_orig=0, nc_dest, nr_dest;

  InitBasisEtc();

  // clear counts
  for (u32 i=0; i<MAX_LODS; i++)
  {
    g_RenderData.m_num_verts[i]   = 0;
    g_RenderData.m_num_indices[i] = 0;
  }
  //g_RenderData.m_b_any_output = false;
  g_RenderData.m_visible = false;

  // clear the grid
  for (i32 i=0; i<nc*nr; i++)
  {
    g_Heightmap[i] = 0.0f;
    g_Outcodes[i]  = 0x80;
  }
  for (i32 r=8; r<true_nr+8; r++)
  {
    for (i32 c=8; c<true_nc+8; c++)
    {
      i32 i = r*nc+c;
      g_Outcodes[i] = 0;
    }
  }

  // approx camera height
  g_RenderData.m_camera_height_above_surface = R2O_ApproxCameraHeightAboveSurface();
  g_RenderData.m_b_camera_over_water = true;


  // lod 0
  u32 window = true;
  u32 oflw   = false;
  {
    // do fft for current band
    if (g_WaterObject.m_first_waveband==0)
    {
      f32 *heightmap = (f32 *)g_OutcodesAlt; // grab 4K temp space
      u32 size = 32*32*sizeof(float);
      u32 ea   = g_R2OCon.m_amb_ea_heightmaps;
      DmaGet((volatile void *)heightmap, ea, size, g_kDmaTag);
      DmaWaitAll(1 << g_kDmaTag);
      R2O_AddAmbientWaves(g_Heightmap, heightmap, nc, nr, g_WaterObject.m_amplitude);
    }

    // generate the frustum for the current lod
    GenerateFrustum();

    // outcode verts relative to current lod's frustum
    GenerateOutcodesAsm(g_Outcodes, g_Heightmap, origin_world, dvc_world, dvr_world, g_Planes, nc, nr);
    SetNearBitsAsm(g_Outcodes, nc, nr);
    ClipToRectangle(clip_min, clip_max);

    // generate derivatives
    i16 *derivs = (i16 *)g_TempBuf;
    //R2O_GenerateVertexDerivs(derivs, g_Heightmap, step, basis_col, basis_row, nc, nr);

    f32 vert_deriv_scale = (0.2f / g_WaterObject.m_amplitude);
    vf32 vscale = spu_splats((0.5f / step) * vert_deriv_scale);    // why is there a 0.5 in there?
    basis_col *= vscale;
    basis_row *= vscale;

    vf32 basis_col_x = spu_splats(spu_extract(basis_col, 0));
    vf32 basis_col_z = spu_splats(spu_extract(basis_col, 2));
    vf32 basis_row_x = spu_splats(spu_extract(basis_row, 0));
    vf32 basis_row_z = spu_splats(spu_extract(basis_row, 2));

    R2O_GenerateVertexDerivs(derivs, g_Heightmap, basis_col_x, basis_col_z, basis_row_x, basis_row_z, nc, nr);


    OutputMesh(g_Heightmap, derivs);
    R2O_RefineCameraHeightAboveSurface();
    R2O_ProcessQueries();

    // determine refinement window; quit if empty
    window = GetRefinementWindow(c_orig, r_orig, nc_dest, nr_dest, oflw, g_Heightmap, g_Outcodes, nc, nr);
  }

  while (window)
  {
    // only g_Heightmap (32K) and g_Outcodes (8K) are in use here

    FlipBuffers();  // g_OutcodesAlt <- g_Outcodes
    NextLod(c_orig, r_orig);

    // send back origin for new lod
    g_RenderData.m_origins[lod]   = origin_world;

    // set new ambient tile origin
    R2O_ChangeAmbientTileOrigin(c_orig, r_orig);



    // -----------------------------------------------------------------------------------------------------
    // refinement copy operation (which rotates the refinement window by 45 degrees and interleaves the output points)
    f32 *refined = g_TempBuf + (g_R2OCon.m_max_verts>>1);  // max_verts/2 floats = 1/4 of the way through g_TempBuf
    DmaWaitAll(1<<0);  // sync to transfer going out of first half of g_TempBuf
    if (!(lod & 1))
    {
      R2O_RefinementCopyMajorEven(refined, nc_dest, nr_dest, g_Heightmap+(r_orig*nc+c_orig), nc);

      // DD-interpolation (minus the linear interpolant) on the original grid without changing the orientation
      R2O_InterpolateDDMinusLinear(g_Heightmap, nc, nr, g_R2OCon.m_DD_coeff0);

      // now do the refinement copy on the DD-interpd values
      // (need to counter-shift the origin to undo what the DD-interpolation did)
      f32 *dd_points = g_TempBuf;
      R2O_RefinementCopyMinorEven(dd_points, nc_dest, nr_dest, g_Heightmap+((r_orig-1)*nc+c_orig-1), nc);

      // the original heightmap is now free, so we can tile the ambient values in it
      // prepare the ambient waves for this lod
      // (these become the morph deltas)
      if (lod >= (i32)g_WaterObject.m_first_waveband)
      {
        f32 *heightmap = (f32 *)g_Outcodes; // grab 4K temp space
        u32 size = 32*32*sizeof(float);
        u32 ea   = g_R2OCon.m_amb_ea_heightmaps + (lod * size);
        DmaGet((volatile void *)heightmap, ea, size, g_kDmaTag);
        DmaWaitAll(1 << g_kDmaTag);
        R2O_ReplicateAmbient(g_Heightmap, nc_dest, nr_dest, heightmap, c0_amb, r0_amb, g_WaterObject.m_amplitude);

        // add DD-values into the morph deltas
        R2O_AddDD(g_Heightmap, dd_points, nc_dest, nr_dest);
      }
      else
      {
        // copy DD-values into the morph deltas
        R2O_CopyDD(g_Heightmap, dd_points, nc_dest, nr_dest);
      }


      // add interactive height maps
      i32 l = lod>>1;
      if (l>=0 && l<MAX_INT_LODS && g_InteractiveData.m_num_tiles[l])
      {
        R2O_AddInteractiveMaps(g_Heightmap, nc_dest, nr_dest, l);
      }
    }
    else
    {
      R2O_RefinementCopyMajorOdd(refined, nc_dest, nr_dest, g_Heightmap+(r_orig*nc+c_orig), nc);

      // DD-interpolation (minus the linear interpolant) on the original grid without changing the orientation
      R2O_InterpolateDDMinusLinear(g_Heightmap, nc, nr, g_R2OCon.m_DD_coeff0);

      // now do the refinement copy on the DD-interpd values
      // (need to counter-shift the origin to undo what the DD-interpolation did)
      f32 *dd_points = g_TempBuf;
      R2O_RefinementCopyMinorOdd(dd_points, nc_dest, nr_dest, g_Heightmap+((r_orig-1)*nc+c_orig-2), nc);

      // copy DD-values into the morph deltas
      R2O_CopyDD(g_Heightmap, dd_points, nc_dest, nr_dest);
    }





    // -----------------------------------------------------------------------------------------------------
    // linear interpolation from the packed points to the alt buffer
    f32 *linear = g_TempBuf + g_R2OCon.m_max_verts;  // max_verts floats = halfway through g_TempBuf
    DmaWaitAll(1<<1);  // sync to transfer going out of second half of g_TempBuf
    R2O_InterpLinear(linear, refined, nc_dest, nr_dest);


    // -----------------------------------------------------------------------------------------------------
    // geomorphing
    f32 near = g_R2OCon.m_near0 * powf(0.5f, 0.5f*(f32)lod);
    f32 far  = g_R2OCon.m_near0 * powf(0.5f, 0.5f*(f32)(lod-1));
    f32 z0   = near + 0.75f*(far-near);
    f32 z1   = near + 0.25f*(far-near);
    f32 *full_morph = linear;

    f32  scale    = 1.0f / (z1-z0);
    vf32 blend0   = spu_splats((spu_extract(origin_camera,2) - z0) * scale);
    vf32 dblend_c = spu_splats(spu_extract(dvc_camera,2) * scale);
    vf32 dblend_r = spu_splats(spu_extract(dvr_camera,2) * scale);
    vf32 dblend_y = spu_splats(spu_extract(*(vf32 *)&g_pViewData->m_world_to_camera_matrix.m_v1,2) * scale);
    
    blend0 = blend0 + (vf32){0,1,2,3} * dblend_c;
    dblend_c = dblend_c * spu_splats(4.0f);
    
    R2O_Geomorph(g_Heightmap, full_morph, linear, nc_dest, nr_dest, blend0, dblend_c, dblend_r, dblend_y);



    // -----------------------------------------------------------------------------------------------------
    // generate derivatives
    i16 *derivs = (i16 *)g_TempBuf;

    f32 vert_deriv_scale = (0.2f / g_WaterObject.m_amplitude);
    vf32 vscale = spu_splats((0.5f / step) * vert_deriv_scale);    // why is there a 0.5 in there?
    basis_col *= vscale;
    basis_row *= vscale;

    vf32 basis_col_x = spu_splats(spu_extract(basis_col, 0));
    vf32 basis_col_z = spu_splats(spu_extract(basis_col, 2));
    vf32 basis_row_x = spu_splats(spu_extract(basis_row, 0));
    vf32 basis_row_z = spu_splats(spu_extract(basis_row, 2));

    R2O_GenerateVertexDerivs(derivs, full_morph, basis_col_x, basis_col_z, basis_row_x, basis_row_z, nc_dest, nr_dest);


    // -----------------------------------------------------------------------------------------------------
    // outcode refinement
    u8 *temp_outcodes_0 = (u8 *)(g_TempBuf + g_R2OCon.m_max_verts);
    u8 *temp_outcodes_1 = (u8 *)(g_TempBuf + 3*g_R2OCon.m_max_verts/2);

    if (nc & 8)
    {
      R2O_PadOutcodesAsm(temp_outcodes_0, g_OutcodesAlt, nc, nr);
    }

    u8 *src        = (nc      & 8) ? temp_outcodes_0 : g_OutcodesAlt;
    i32 src_width  = (nc      & 8) ? nc+8            : nc;
    u8 *dest       = (nc_dest & 8) ? temp_outcodes_1 : g_Outcodes;
    i32 dest_width = (nc_dest & 8) ? nc_dest+8       : nc_dest;

    if (lod & 1)
    {
      RefineOutcodesOddAsm(dest, dest_width, nr_dest, src, src_width, nr, c_orig, r_orig);
    }
    else
    {
      RefineOutcodesEvenAsm(dest, dest_width, nr_dest, src, src_width, nr, c_orig, r_orig);
    }
    
    if (nc_dest & 8)
    {
      R2O_PackOutcodesAsm(g_Outcodes, nc_dest, nr_dest, temp_outcodes_1);
    }


    // -----------------------------------------------------------------------------------------------------
    // outcode interpolation

    // horizontally
    InterpolateOutcodesHorizontalAsm(g_Outcodes, nc_dest, nr_dest);

    // vertically
    if (nc_dest & 8)
    {
      InterpolateOutcodesVerticalOddAsm(g_Outcodes, nc_dest, nr_dest);
    }
    else
    {
      InterpolateOutcodesVerticalEvenAsm(g_Outcodes, nc_dest, nr_dest);
    }

    // right col
    for (i32 r=0; r<nr_dest; r++)
    {
      i32 i = r*nc_dest + nc_dest-1;
      g_Outcodes[i] = 0x80;
    }

    // bottom row
    for (i32 c=0; c<nc_dest; c++)
    {
      i32 i = (nr_dest-1)*nc_dest + c;
      g_Outcodes[i] = 0x80;
    }




    // -----------------------------------------------------------------------------------------------------
    // handle window overflows
    if (oflw)
    {
      SuppressRows(g_Outcodes, nc_dest, nr_dest, g_R2OCon.m_suppress_rows);
    }


    // -----------------------------------------------------------------------------------------------------
    // set new window
    nc = nc_dest;
    nr = nr_dest;
    g_RenderData.m_cols_rows[lod] = nc<<8 | nr;

    // -----------------------------------------------------------------------------------------------------
    // generate the frustum for the current lod
    GenerateFrustum();

    // outcode verts relative to current lod's frustum
    GenerateOutcodesAsm(g_Outcodes, g_Heightmap, origin_world, dvc_world, dvr_world, g_Planes, nc, nr);
    SetNearBitsAsm(g_Outcodes, nc, nr);

    // -----------------------------------------------------------------------------------------------------
    // clip to bounding rectangle and cut any holes for this lod
    if (!(lod&1))
    {
      ClipToRectangle(clip_min, clip_max);
      if (g_R2OCon.m_num_holes)
      {
        R2O_CutHoles(g_Outcodes, nc, nr, lod);
      }
    }

    OutputMesh(g_Heightmap, derivs);
    R2O_RefineCameraHeightAboveSurface();
    R2O_ProcessQueries();

    if (lod==g_R2OCon.m_num_lods-1)
    {
      break;
    }


    // determine refinement window; quit if empty
    window = GetRefinementWindow(c_orig, r_orig, nc_dest, nr_dest, oflw, g_Heightmap, g_Outcodes, nc, nr);
  }

  // wait for OutputMesh transfers to finish since we're about to free the buffers
  DmaWaitAll(1<<0 | 1<<1);

  // restore scratchpad state
  g_Scratchpad.Reset(state);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_Render()
{
  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // loop over water objects
  u32 ea  = g_R2OCon.m_ea_water_objects;
  for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++)
  {
    // get the water object and interactive data
    u32 ofs = OFFSETOF(R2OWaterObject, m_interactive_data);
    R2O_GetWait(&g_WaterObject, ea, i_obj, sizeof(R2OBasicWaterObject), sizeof(R2OWaterObject));
    if (!(g_WaterObject.m_flags & R2O_WATER_OBJECT_FLAG_ACTIVE))
    {
      continue;
    }
    R2O_GetWait(&g_InteractiveData, ea+ofs, i_obj, sizeof(R2OInteractiveData), sizeof(R2OWaterObject));

    // render to each viewport
    R2OViewData *p_view_data = &g_R2OCon.m_view_data[0];
    ofs      = OFFSETOF(R2OWaterObject, m_render_data);
    u32 ofs1 = OFFSETOF(R2ORenderData, m_visible);
    for (u32 i_vp=0; i_vp<g_R2OCon.m_num_viewports; i_vp++, ofs+=sizeof(R2ORenderData))
    {
      // test visibility
      R2O_GetWait(&g_RenderData.m_visible, ea+ofs+ofs1, i_obj, sizeof(g_RenderData.m_visible), sizeof(R2OWaterObject));
      if (g_RenderData.m_visible)
      {
        // render water object
        R2O_Render(p_view_data);

        // use the render data to generate the interactive index maps
        R2O_GenerateIndexMaps();

        // pass back the output for this object
        R2O_PutWait(&g_RenderData, ea+ofs, i_obj, sizeof(R2ORenderData), sizeof(R2OWaterObject));

        // accumulate visibility result
        g_R2OCon.m_view_data[i_vp].m_any_visible |= g_RenderData.m_visible;
        g_R2OCon.m_any_visible |= g_RenderData.m_visible;
      }
    }
  }

  // reset the scratch memory
  g_Scratchpad.Reset(state);
}
