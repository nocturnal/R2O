#include "igspu\igTypesSpu.h"

#include "igR2OSpu.h"




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 32x32 2D FFT
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_FFT2D(f32 z[], u32 b_inv)
{
  if (!b_inv)
  {
    // forward ffts

    // transform each row
    R2O_FFT32_DIF_rows_asm(z, 32);

    R2O_FFT32_DIF_cols_asm(z, 32);
  }
  else
  {
    // inverse ffts

    // transform each col
    R2O_InvFFT32_DIT_cols_asm(z, 32);

    // transform each row
    R2O_InvFFT32_DIT_rows_asm(z, 32);
  }
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 32x32 2D FFT using separate real and imaginary arrays
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_FFT2D_sep(f32 z[], u32 b_inv)
{
  if (!b_inv)
  {
    // forward ffts

    // transform each row
    R2O_FFT32_DIF_rows_sep_asm(z, z+32*32, 32);

    R2O_FFT32_DIF_cols_sep_asm(z, z+32*32, 32);
  }
  else
  {
    // inverse ffts

    // transform each col
    R2O_InvFFT32_DIT_cols_sep_asm(z, z+32*32, 32);

    // transform each row
    R2O_InvFFT32_DIT_rows_sep_asm(z, z+32*32, 32);
  }
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compute an array of complex exponentials z=exp(2*pi*i*w*dt) given the angular frequencies w and the time-step dt
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_GenerateComplexExponentials(Complex z[], f32 ang_freqs[], f32 dt, u32 cnt, f32 scale)
{
  f32 two_pi_dt = 2.0f * PI * dt;

  for (u32 i=0; i<cnt; i++)
  {
    z[i] = scale * Expi(two_pi_dt * ang_freqs[i]);
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Convert a paletted 8-bit 'image' of complex values int explicit form
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_DecompressComplexArray(Complex z[], Complex palette[], u8 map[], u32 cnt)
{
  for (u32 i=0; i<cnt; i++)
  {
    u8 idx = *map++;
    *z++ = palette[idx];
  }
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Multipy the elements of one complex array by those of another
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_MultiplyComplexArrays(Complex *z, Complex *w, u32 cnt)
{
  for (u32 i=0; i<cnt; i++)
  {
    *z++ *= *w++;
  }
}




