#ifndef __IGR2ORENDERSPU_H
#define __IGR2ORENDERSPU_H


#include "..\PPU\igR2O_shared.h"
#include "igR2OAmbientSpu.h"



#define DD_TAPS  4

#define INTERP_O(O0,O1)  ( ( ((O0) | (O1)) & 0x80 ) | ( ((O0) & (O1)) & 0x1 ) << 7 )


void R2O_PreRender();
void R2O_RenderObject(u32 i_obj, u32 i_vp);
void R2O_Render(R2OViewData *p_view_data);
void R2O_Render();

extern i32 nc, nr, nc_dest, nr_dest;
extern i32 lod;
extern f32 *g_Heightmap;
extern f32 step, inv_step;
extern vf32 origin_camera, dvc_camera, dvr_camera;
extern vf32 origin_world,  dvc_world,  dvr_world;
extern R2OViewData *g_pViewData;
extern R2OHole *g_Holes;


#endif // __IGR2ORENDERSPU_H

