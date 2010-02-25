#pragma once


#include "..\igR2O\PPU\igR2OInteractive_shared.h"

#define QUICK_LOG2(X) ((si_to_int(si_from_float(X)) >> 23) - 127)


struct R2OImpulse
{
  vf32 m_disk;
  u32  m_flags;
  u32  m_min_l;
  u32  m_max_l;
  f32  m_magnitude;
}
ALIGNED(16);


void R2O_AddInteractiveMaps(f32 *heights, i32 width, i32 height, u32 lod);
void R2O_UpdateInteractive();
void R2O_PreUpdateInteractive();
void R2O_UpdateInteractivePreObject(u32 i_obj);
void R2O_UpdateInteractive(u32 i_obj, u32 lod);
void R2O_PostUpdateInteractive();
void R2O_GenerateIndexMaps();


extern u8 g_HelpPropagation, g_HelpImpulses;
extern Tile *g_Tiles;
extern R2OImpulse *g_Impulses;

