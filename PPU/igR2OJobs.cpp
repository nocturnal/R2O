#include <stdlib.h>

#include "igSys/igCommon.h"
#include "igMath/igMath.h"
#include "igSys/igDebug.h"
#include "igg/igColorEnum.h"
#include "igSpuJob/PPU/igSpuJobLoader.h"

#include "igR2O.h"


extern u8 _binary_igR2OSpu_elf_bin_start[];



namespace IGG
{

// Water module logic

#define s_R2OSpuId       4

volatile bool                  s_R2OBreak      = false;
static SpuJob::SpuModuleHandle s_bR2OModHandle (_binary_igR2OSpu_elf_bin_start, "R2O", COLOR_BLUE_LIGHT);
static SpuJob::SpuJobHandle    s_bR2OJobHandle = (SpuJob::SpuJobHandle)0;
static const char*             s_bR2OJobName   = "igR2O";



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_SendJob()
{
  if (!g_R2OCon.m_r2o_enable)
  {
    return;
  }

  SpuJob::SpuJobParams job_params;
  job_params.m_u32[0] = (u32)&g_R2OCon;

  s_bR2OJobHandle = SpuJob::RunSpuJob(s_bR2OModHandle, job_params, s_bR2OJobName, s_R2OSpuId, s_R2OBreak);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_InitJob()
{
  if (!g_R2OCon.m_r2o_enable)
  {
    return;
  }

  s_bR2OJobHandle = (SpuJob::SpuJobHandle)0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_Sync()
{
  if (!g_R2OCon.m_r2o_enable)
  {
    return;
  }

  SpuJob::SyncJob(s_bR2OJobHandle);

  R2O_TriggerEventsFromImpulses();

  R2O_CameraUnderwaterTest();
}






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#if R2O_TASKS

#define R2O_MAX_TASKS 256
R2OTask g_R2OTasks[R2O_MAX_TASKS];


void R2O_AddTask(R2OTaskType type, u32 obj, u32 lod)
{
  if (g_R2OCon.m_num_tasks < R2O_MAX_TASKS)
  {
    R2OTask *p_task = &g_R2OTasks[g_R2OCon.m_num_tasks];

    p_task->m_type = type;
    p_task->m_obj  = obj;
    p_task->m_idx  = lod;

    g_R2OCon.m_num_tasks++;
  }
}



void R2O_BuildTaskList()
{
  if (!g_R2OCon.m_r2o_enable)
  {
    return;
  }

  g_R2OCon.m_num_tasks = 0;


  // add a single task for ambient update
  R2O_AddTask(R2O_TASK_UPDATE_AMB, 0, 0);


  // add a single task for ambient normal map generation
  R2O_AddTask(R2O_TASK_GENERATE_NORMS_AMB, 0, 0);


  // add a single task for interactive pre-update
  R2O_AddTask(R2O_TASK_PRE_UPDATE_INT, 0, 0);


  // loop over active water objects
  R2OWaterObject* p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++, p_obj++)
  {
    if (!(p_obj->m_flags & R2O_WATER_OBJECT_FLAG_ACTIVE))
    {
      continue;
    }

    // create an interactive update 'pre-object' task for this object
    R2O_AddTask(R2O_TASK_UPDATE_INT_PRE_OBJECT, i_obj, 0);

    // create an interactive update task for each potentially active lod of each object
    for (u32 l=0,lod=0; l<MAX_INT_LODS; l++,lod+=2)
    {
      if (p_obj->m_interactive_data.m_num_tiles[l] ||
          (p_obj->m_interactive_data.m_impulse_lods & (0x1<<l)))
      {
        R2O_AddTask(R2O_TASK_UPDATE_INT, i_obj, lod);
      }
    }
  }
  

  // add a single task for interactive post-update
  R2O_AddTask(R2O_TASK_POST_UPDATE_INT, 0, 0);


  // add a single task for pre-render
  R2O_AddTask(R2O_TASK_PRE_RENDER, 0, 0);


  // loop over active water objects
  p_obj = (R2OWaterObject *)g_R2OCon.m_ea_water_objects;
  for (i32 i_obj=0; i_obj<=g_R2OCon.m_object_last; i_obj++, p_obj++)
  {
    if (!(p_obj->m_flags & R2O_WATER_OBJECT_FLAG_ACTIVE))
    {
      continue;
    }

    // create a render task for each visible object visible, in each viewport
    for (u32 i_vp=0; i_vp<g_R2OCon.m_num_viewports; i_vp++)
    {
      if (p_obj->m_render_data[i_vp].m_visible)
      {
        R2O_AddTask(R2O_TASK_RENDER, i_obj, i_vp);
      }
    }
  }
}

#endif


} // namespace IGG
