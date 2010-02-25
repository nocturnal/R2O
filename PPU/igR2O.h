#pragma once

#include "igg/igGraphics.h"

#define  R2O_PPU
#include "igR2O_shared.h"
#include "igR2OAmbient.h"
#include "igR2OInteractive.h"
#include "igR2OJobs.h"
#include "igR2OQuery.h"
#include "igR2OReflection.h"
#include "igUnderwater/PPU/igUnderWater.h"

namespace IGG
{

#define R2O_HEIGHT_SCALE 4096.0f


struct BucketCon;

void R2O_Init(u32 interactive_tile_count);
void R2O_FrameInit();
void R2O_AddWaterObjects(R2OBuiltWaterObject* built_water_objects, u32 num_water_objects, i64 zone_index);
void R2O_DeleteUnusedWaterObjects(i64 zone_alloc_mask);
void R2O_Update(BucketCon* bucketer, f32 time_step);
void R2O_Update(f32 time_step);
bool R2O_AllocDownsampleBuffers();
void R2O_VisibilityTests(Frustum *frustum_pointers[], u32 num_views);
void R2O_RenderShell(BucketCon* bucketer, i32 vp_idx);
void R2O_Update(f32 time_step, u32 num_viewports);
void R2O_GenerateFresnelTexture();
u32  R2O_TicksLastFrame();

void R20_DisableZoneWaterObjects(i64 zone_index);
void R20_EnableZoneWaterObjects(i64 zone_index);

inline void R2O_Start(f32 time_step, u32 num_viewports)
{
  R2O_Update(time_step, num_viewports);
  R2O_SendJob();
}

extern R2OCon g_R2OCon;
extern f32    g_MaxWaterDepth;


} // namespace IGG

