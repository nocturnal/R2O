
#include "igR2OSpu.h"



#define SCRATCHPAD_SIZE (160*1024)


// globals
R2OCon              g_R2OCon;
R2OBasicWaterObject g_WaterObject;
R2ORenderData       g_RenderData;
R2OInteractiveData  g_InteractiveData;

u8 g_ScratchpadMemory[SCRATCHPAD_SIZE]   ALIGNED(16) UNINITIALIZED;
R2O_Allocator g_Scratchpad;

u32 g_TimerStart, g_TimerTotal;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

u8 *R2O_Allocator::Alloc(u32 size, u32 alignment, u32 beyond)
{
  IG_ASSERT(alignment == (alignment & -alignment));
  IG_ASSERT(beyond < alignment);

  m_offset = ((m_offset + alignment-1 - beyond) & ~(alignment-1)) + beyond;
  u8* p_mem = (u8*)(m_base + m_offset);
  m_offset += size;

  // overflow check
  if (m_offset > m_size)
  {
    Printf("m_offset=%d, m_size=%d, size=%d, alignment=%d, beyond=%d\n",
            m_offset   , m_size   , size   , alignment   , beyond      );
  }
  IG_ASSERT(m_offset <= m_size);

  return p_mem;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Use these for dma's over 16K
void R2ODmaGetBig(volatile void* ls, u32 ea, u32 size, u32 tag)
{
  while (size > 16384)
  {
    DmaGet(ls, ea, 16384, tag);
    ls = (volatile void *)((u32)ls+16384);
    ea += 16384;
    size -= 16384;
  }
  DmaGet(ls, ea, size, tag);
}


void R2ODmaPutBig(const volatile void* ls, u32 ea, u32 size, u32 tag)
{
  while (size > 16384)
  {
    DmaPut(ls, ea, 16384, tag);
    ls = (volatile void *)((u32)ls+16384);
    ea += 16384;
    size -= 16384;
  }
  DmaPut(ls, ea, size, tag);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_GetWait(void *ls, u32 ea, u32 idx, u32 size, u32 stride)
{
  stride = stride ? stride : size;
  R2ODmaGetBig((volatile void *)ls, ea+idx*stride, size, g_kDmaTag);
  DmaWaitAll(1 << g_kDmaTag);
}


void R2O_PutWait(void *ls, u32 ea, u32 idx, u32 size, u32 stride)
{
  stride = stride ? stride : size;
  R2ODmaPutBig((volatile void *)ls, ea+idx*stride, size, g_kDmaTag);
  DmaWaitAll(1 << g_kDmaTag);
}


void *R2O_AllocGetWait(u32 ea, u32 idx, u32 size, u32 stride)
{
  stride = stride ? stride : size;
  u8 *ls = g_Scratchpad.Alloc(size);
  R2ODmaGetBig((volatile void *)ls, ea+idx*stride, size, g_kDmaTag);
  DmaWaitAll(1 << g_kDmaTag);
  return (void *)ls;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if R2O_TASKS

R2OTask s_Task;

void R2O_ProcessTasks(u32 ea_r2o_con)
{
  for (u32 i=0; i<g_R2OCon.m_num_tasks; i++)
  {
    R2O_GetWait(&s_Task, g_R2OCon.m_ea_task_list, i, sizeof(R2OTask));

    switch (s_Task.m_type)
    {
      case R2O_TASK_UPDATE_AMB:
      {
        R2O_UpdateAmbient(g_R2OCon.m_time_step);
        break;
      }

      case R2O_TASK_GENERATE_NORMS_AMB:
      {
        R2O_GenerateAmbientNormalMaps();
        break;
      }

      case R2O_TASK_PRE_UPDATE_INT:
      {
        R2O_PreUpdateInteractive();
        break;
      }

      case R2O_TASK_UPDATE_INT_PRE_OBJECT:
      {
        R2O_UpdateInteractivePreObject(s_Task.m_obj);
        break;
      }

      case R2O_TASK_UPDATE_INT:
      {
        R2O_UpdateInteractive(s_Task.m_obj, s_Task.m_idx);    // idx = lod
        break;
      }

      case R2O_TASK_POST_UPDATE_INT:
      {
        R2O_PostUpdateInteractive();
        break;
      }

      case R2O_TASK_PRE_RENDER:
      {
        R2O_PreRender();
        break;
      }

      case R2O_TASK_RENDER:
      {
        R2O_RenderObject(s_Task.m_obj, s_Task.m_idx);            // idx = viewport
        break;
      }

      // ...
      // (perhaps split off index map generation from the render task)
    }
  }
}

#endif



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_Update(u32 ea_r2o_con)
{
  #if R2O_TASKS
  // process the task list
  R2O_ProcessTasks(ea_r2o_con);
  #else
  // update ambient system
  R2O_UpdateAmbient(g_R2OCon.m_time_step);

  // create normal maps
  R2O_GenerateAmbientNormalMaps();

  // update interactive system
  R2O_UpdateInteractive();

  // init stuff for render
  R2O_PreRender();

  // render
  R2O_Render();
  #endif
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C" void JobMain(u32 ea_r2o_con)
{
  TIMER_INIT;

  // get the r2o data
  DmaGet(&g_R2OCon, ea_r2o_con, sizeof(R2OCon), g_kDmaTag);
  DmaWaitAll(1 << g_kDmaTag);

  // timing
  u32 r1 = spu_read_decrementer();

  // init scratchpad memory
  g_Scratchpad.Init(g_ScratchpadMemory, SCRATCHPAD_SIZE);

  // update
  R2O_Update(ea_r2o_con);

  // timing
  g_R2OCon.m_spu_ticks = r1 - spu_read_decrementer();

  DmaPut(&g_R2OCon, ea_r2o_con, sizeof(R2OCon), g_kDmaTag);
  DmaWaitAll(1 << g_kDmaTag);
}

