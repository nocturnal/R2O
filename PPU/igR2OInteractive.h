#pragma once

namespace IMPULSE
{
  struct ImpulseData;
}

namespace IGG
{

#include "igR2OQuery_shared.h" 
#include "igR2OInteractive_shared.h" 

inline u32 QUICK_LOG2(f32 x)
{
  IntFloat temp;
  temp.m_f32 = x;
  return (temp.m_i32 >> 23) - 127;
}


void R2O_InitInteractive(u32 interactive_tile_count);
void R2O_InitInteractiveTest();
void R2O_DeleteInteractive(R2OInteractiveData *p_int_data);
void R2O_AssignImpulsesToWaterObjects();
void R2O_TriggerEventsFromImpulses();
void R2O_SetInteractiveGrid(vec4f_arg origin, u32 lod);
void R2O_GetInteractiveGrid(vf32 &origin, u32 &lod);
void R2O_GetInteractiveGrid(vec4f &min, vec4f& max);
void R2O_SetInteractiveHalfLife(f32 seconds);
void R2O_SetInteractiveNormalMapTexture(u32 tex_unit, u32 ea);
void R2O_SetInteractiveIndexMapTexture(u32 tex_unit, u32 ea, u32 width, u32 height);
void R2O_HardCodeHeightMaps();
void R2O_HardCodeTileLists(R2OInteractiveData *p_int_data);

typedef void (*R2OInteractionCallback)(IMPULSE::ImpulseData* impulse_data);
extern R2OInteractionCallback g_R2OInteractionCallback;
void R2O_SetInteractionCallback(R2OInteractionCallback);

} // namespace IGG

