#pragma once

#define  R2O_PPU
#include "igR2O_shared.h"

namespace IGG
{

void R2O_InitJob();
void R2O_SendJob();
void R2O_Sync();
void R2O_BuildTaskList();

extern R2OTask g_R2OTasks[];

} // namespace IGG
