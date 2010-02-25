#pragma once

extern "C"
{
  void R2O_MultiplyComplexArraysAsm(f32 *z, f32 *w, u32 cnt);
  void R2O_DecompressComplexArrayAsm(f32 *dst, f32 *palette, u8 *map, u32 cnt);
  void R2O_DecompressComplexArraySeparate(f32 *x, f32 *y, f32 *palette, u8 *map, u32 cnt);
  void R2O_ZeroRectangleAsm(f32 *dst, u32 cnt);
  void R2O_ReadRectangleAsm(f32 *dst, i16 *src, u32 cnt);
  void R2O_WriteRectangleCentralAsm(i16 *dest, f32 *src, u32 cnt);
  void R2O_WriteRectangleAlongsideAsm(i16 *dst, f32 *src, u32 cnt);
  void R2O_WriteRectangleCornerAsm(i16 *dst, f32 *src, u32 cnt);
  void R2O_WriteRectangleCentralSeparate(i16 *dst, f32 *x, f32 *y, u32 cnt);
  void R2O_WriteRectangleAdjacentSeparate(i16 *dst, f32 *x, f32 *y, u32 cnt);
  //void R2O_WriteRectangleCornerSeparate(i16 *dst, f32 *x, f32 *y, u32 cnt);
  void R2O_CopyRealsAsm(f32 *dst, f32 *src, u32 cnt);
  u32  R2O_ThresholdRectangleWidth4Asm (f32 *src, u32 height, f32 threshold);
  u32  R2O_ThresholdRectangleWidth16Asm(f32 *src, u32 height, f32 threshold);
  void R2O_ClearMemAsm(qword *mem, u32 cnt);
  void R2O_PostProcessNewHalfAsm(i16 *dst, u32 b_lower);
  void R2O_PrepareFftBuffer(f32 *z, i16 *UL, i16 *U0, i16 *UR, i16 *uL, i16 *u0, i16 *uR, i16 *dL, i16 *d0, i16 *dR, i16 *DL, i16 *D0, i16 *DR);
  void R2O_PrepareFftBufferSeparate(f32 *x, f32 *y, i16 *UL, i16 *U0, i16 *UR, i16 *uL, i16 *u0, i16 *uR, i16 *dL, i16 *d0, i16 *dR, i16 *DL, i16 *D0, i16 *DR);
  void R2O_GatherAndTransformRows(f32 *x, f32 *y, i16 *UL, i16 *U0, i16 *UR, i16 *uL, i16 *u0, i16 *uR, i16 *dL, i16 *d0, i16 *dR, i16 *DL, i16 *D0, i16 *DR);
}

