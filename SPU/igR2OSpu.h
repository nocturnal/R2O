#pragma once

#define DISABLE_WWSJOB

#if !(OPTIMIZED)
  #define DMA_VALIDATE
#endif//!(OPTIMIZED)

#include "igSpu\igDmaSpu.h"

#include "igSpu\igTypesSpu.h"
#include "igSpu\igMathSpu.h"
#include "igSpu\igRenderSpu.h"
#include "igspu\igMiscSpu.h"
#include "igSpurs\igSpuTimer.h"

#include "..\igR2O\PPU\igR2O_shared.h"
#include "igR2OAmbientSpu.h"
#include "igR2OAsmComplexSpu.h"
#include "igR2OAsmInteractiveSpu.h"
#include "igR2OAsmSpu.h"
#include "igR2OComplexSpu.h"
#include "igR2OInteractiveSpu.h"
#include "igR2ONormalMapSpu.h"
#include "igR2OQuerySpu.h"
#include "igR2ORenderSpu.h"



//-------------------------------------------------------------------------------
extern u32 g_TimerStart, g_TimerTotal;
#define TIMER_INIT  g_TimerStart=g_TimerTotal=0;
#define TIMER_CHECK if (__LINE__<=g_R2OCon.m_timer_line)                        \
                      g_TimerStart=TimerGetCount();                             \
                    if (__LINE__> g_R2OCon.m_timer_line && g_TimerStart)        \
                      g_TimerTotal+=TimerGetCount()-g_TimerStart, g_TimerStart=0;
#define TIMER_PRINT PRINT("%g ms\n", TimerTicksToMillis(g_TimerTotal));
//-------------------------------------------------------------------------------




struct R2O_Allocator
{
  u32 m_base;
  u32 m_size;
  u32 m_offset;


  ////////////////////////////////////////////////////////////////////////////////////////////////
  R2O_Allocator(): m_base( 0 ), m_size( 0 ), m_offset( 0 )
  {
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  inline void Init(void* base, u32 size)
  {
    m_base   = (u32)base;
    m_size   = size;
    m_offset = 0;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // we assume alignment is a power of 2, and beyond < alignment
  u8* Alloc(u32 size, u32 alignment=16, u32 beyond=0);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  inline void Reset(u32 state=0)
  {
    m_offset = state;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  inline u32 GetState()
  {
    return m_offset;
  }

} ALIGNED(16);



void R2ODmaGetBig(volatile void* ls, u32 ea, u32 size, u32 tag);
void R2ODmaPutBig(const volatile void* ls, u32 ea, u32 size, u32 tag);
void *R2O_AllocGetWait(    u32 ea, u32 idx, u32 size, u32 stride=0);
void R2O_GetWait(void *ls, u32 ea, u32 idx, u32 size, u32 stride=0);
void R2O_PutWait(void *ls, u32 ea, u32 idx, u32 size, u32 stride=0);

void R2O_Update();

extern R2OCon             g_R2OCon;
extern R2OBasicWaterObject g_WaterObject;
extern R2ORenderData      g_RenderData;
extern R2OInteractiveData g_InteractiveData;
extern R2O_Allocator      g_Scratchpad;

enum
{
  g_kDmaTag = 20,
};

