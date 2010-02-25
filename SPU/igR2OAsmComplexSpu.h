#pragma once

extern "C"
{
  void R2O_FFT32_DIF_rows_asm(f32 *z, u32 cnt);
  void R2O_InvFFT32_DIT_rows_asm(f32 *z, u32 cnt);
  void R2O_FFT32_DIF_cols_asm(f32 *z, u32 cnt);
  void R2O_InvFFT32_DIT_cols_asm(f32 *z, u32 cnt);
  void R2O_InvFFT32_DIT_cols_sep_asm(f32 *x, f32 *y, u32 cnt);
  void R2O_FFT32_DIF_cols_sep_asm(f32 *x, f32 *y, u32 cnt);
  void R2O_TransformColsAsm(f32 *z, f32 *phase_shifts, u32 cnt);
  void R2O_InvFFT32_DIT_rows_sep_asm(f32 *x, f32 *y, u32 cnt);
  void R2O_FFT32_DIF_rows_sep_asm(f32 *x, f32 *y, u32 cnt);
}
