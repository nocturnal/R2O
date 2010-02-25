#pragma once

#include "igSpu\igMathSpu.h"

#define PI 3.14149265359f


struct Complex
{
  float re;
  float im;

  inline Complex(void);
  inline Complex(float x);
  inline Complex(float x, float y);

  inline void operator =(float x);
  inline void operator *=(float x);
  inline void operator *=(Complex z);
  inline void operator +=(Complex z);
};


inline Complex::Complex(void)
{
  re=im=0.0f;
}

inline Complex::Complex(float x)
{
  re=x;
  im=0.0f;
}

inline Complex::Complex(float x, float y)
{
  re=x;
  im=y;
}

inline Complex operator -(Complex z)
{
  z.re = -z.re;
  z.im = -z.im;
  return z;
}

inline Complex operator +(Complex z0, Complex z1)
{
  Complex z;
  z.re = z0.re + z1.re;
  z.im = z0.im + z1.im;
  return z;
}

inline Complex operator -(Complex z0, Complex z1)
{
  Complex z;
  z.re = z0.re - z1.re;
  z.im = z0.im - z1.im;
  return z;
}

inline Complex operator *(Complex z0, Complex z1)
{
  Complex z;
  z.re = z0.re * z1.re - z0.im * z1.im;
  z.im = z0.re * z1.im + z0.im * z1.re;
  return z;
}

inline Complex operator *(float x, Complex z)
{
  z.re *= x;
  z.im *= x;
  return z;
}

inline Complex operator *(Complex z, float x)
{
  z.re *= x;
  z.im *= x;
  return z;
}

inline Complex operator /(Complex z, float x)
{
  return z*(1.0f/x);
}

inline void Complex::operator =(float x)
{
  re = x;
  im = 0.0f;
}

inline void Complex::operator *=(float x)
{
  re *= x;
  im *= x;
}

inline void Complex::operator *=(Complex z)
{
  float r = re;
  re = re*z.re - im*z.im;
  im = r *z.im + im*z.re;
}

inline void Complex::operator +=(Complex z)
{
  re += z.re;
  im += z.im;
}

inline Complex Conjugate(Complex z0)
{
  Complex z;
  z.re =  z0.re;
  z.im = -z0.im;
  return z;
}

inline Complex Expi(float theta)
{
  Complex z;
  spu_sincos(z.im, z.re, theta);
  return z;
}




void R2O_FFT2D(f32 z[], u32 b_inv);
void R2O_FFT2D_sep(f32 z[], u32 b_inv);
void R2O_GenerateComplexExponentials(Complex z[], f32 ang_freqs[], f32 dt, u32 cnt, f32 scale=1.0f);
void R2O_DecompressComplexArray(Complex z[], Complex palette[], u8 map[], u32 cnt);
void R2O_MultiplyComplexArrays(Complex *z, Complex *w, u32 cnt);

