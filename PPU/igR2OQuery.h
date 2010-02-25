#pragma once

#include "igR2OQuery_shared.h"

namespace IGG
{

extern R2OQueryGroup g_QueryGroups[MAX_QUERY_GROUPS];
extern u32 g_NumRegisteredQueryGroups;

void R2O_ResetQueryGroups();
R2OQueryGroup *R2O_RegisterQueries(void *queries, u32 num_queries, u32 indirect, R2OWaterObject *p_obj);
R2OWaterObject *R2O_ClosestWaterObject(vec4f &point);
R2OWaterObject *R2O_ClosestWaterObjectAboveOrBelow(vec4f &point);
void R2O_CameraUnderwaterTest();
void R2O_SetHeightQueryCallback(void (*function)(R2OWaterObject *));

extern void (*g_R2OHeightQueryCallback)(R2OWaterObject *);

}; // namespace IGG
