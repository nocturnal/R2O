#include "igg/igGraphics.h"
#include "igR2O.h"


namespace IGG
{


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Map out the wavenumbers to find how many distinct magnitudes there are (this number becomes the palette size)
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static int CmpU16(const void *pa, const void *pb)
{
  u16 a=*(u16 *)pa, b=*(u16 *)pb;
  return a<b ? -1 : a>b ? 1 : 0;
}

u32 R2O_GenerateWaveMap(u16 wavenumbers_squared[], u32 wavenum_cnt, u32 fft_width)
{
  // clear the table initially (although we don't actually need to)
  for (u32 i=0; i<wavenum_cnt; i++)
  {
    wavenumbers_squared[i] = 0;
  }

  // make a list of the distinct wavenumber magnitudes squared
  u32 cnt=0;
  u32 half_fft_width = fft_width >> 1;
  for (u32 r=0; r<fft_width; r++)
  {
    i32 kz = (r < half_fft_width) ? r : (i32)r-fft_width;
    for (u32 c=0; c<fft_width; c++)
    {
      i32 kx = (c < half_fft_width) ? c : (i32)c-fft_width;
      u32 sq = kx*kx + kz*kz;
      u32 i;
      for (i=0; i<cnt; i++)
      {
        if (wavenumbers_squared[i] == sq)
        {
          break;
        }
      }
      if (i == cnt)
      {
        wavenumbers_squared[cnt++] = sq;
      }
    }
  }

  // sort them (although we don't actually need to)
  IGG::Qsort(wavenumbers_squared, cnt, sizeof(u16), CmpU16);

  // return palette size
  return cnt;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// All the allocs required by the ambient wave system
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_InitAmbientMemory()
{
  u32 num_bands = g_R2OCon.m_amb_num_bands;
  u32 pal_cnt   = g_R2OCon.m_amb_palette_cnt;
  u32 map_cnt   = 32*32;

  g_R2OCon.m_amb_ea_angle_palette      = (u32)g_SystemMemoryAlloc.Alloc(pal_cnt * num_bands * sizeof(f32), 128);
  g_R2OCon.m_amb_ea_frequency_palette  = (u32)g_SystemMemoryAlloc.Alloc(pal_cnt * num_bands * sizeof(f32), 128);
  g_R2OCon.m_amb_ea_phasor_palette     = (u32)g_SystemMemoryAlloc.Alloc(pal_cnt * num_bands * 2*sizeof(f32), 128);
  g_R2OCon.m_amb_ea_wave_scale_palette = (u32)g_SystemMemoryAlloc.Alloc(pal_cnt * sizeof(f32), 128);
  g_R2OCon.m_amb_ea_index_map          = (u32)g_SystemMemoryAlloc.Alloc(map_cnt * sizeof(u8), 128);
  g_R2OCon.m_amb_ea_spectrum           = (u32)g_SystemMemoryAlloc.Alloc(num_bands * map_cnt * 2 * sizeof(i8), 128);
  g_R2OCon.m_amb_ea_heightmaps         = (u32)g_SystemMemoryAlloc.Alloc(num_bands * map_cnt * sizeof(f32), 128);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// bit-reverse an integer in the range 0 to n-1
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int R2O_BitReverse(int k, int n)
{
  int r=0;
  for (int bk=n>>1,br=1; bk; bk>>=1,br<<=1)
  {
    if (k & bk)
    {
      r |= br;
    }
  }
  return r;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// reorder the contents of byte array a in order of bit-reversed array subscripts
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_OrderBitReversed(u8 a[], u32 n, u32 s)
{
  for (u32 i=0; i<n; i++)
  {
    u32 r = R2O_BitReverse(i, n);
    if (i<r)
    {
      // swap element i with element r
      u8 temp = a[r*s];
      a[r*s] = a[i*s];
      a[i*s] = temp;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Generate the index map which will be used by the ambient wave palettes (same map for all lods)
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_GenerateIndexMap(u8 map[], u16 wavenumbers_squared[], u32 fft_width, u32 pal_cnt)
{
  // construct the index map for the palettes
  u32 half_fft_width = fft_width >> 1;
  u8 *p_idx = map;
  for (u32 r=0; r<fft_width; r++)
  {
    i32 kz = (r < half_fft_width) ? r : (i32)r-fft_width;
    for (u32 c=0; c<fft_width; c++)
    {
      i32 kx = (c < half_fft_width) ? c : (i32)c-fft_width;
      u32 sq = kx*kx + kz*kz;
      u32 i;
      for (i=0; i<pal_cnt; i++)
      {
        if (wavenumbers_squared[i] == sq)
        {
          *p_idx++ = (u8)i;
          break;
        }
      }
      IG_ASSERT(i<pal_cnt);
    }
  }

  // reorder each row in bit-reversed order
  for (u32 i=0; i<fft_width; i++)
  {
    R2O_OrderBitReversed(map+i*fft_width, fft_width, 1);
  }

  // reorder each col in bit-reversed order
  for (u32 i=0; i<fft_width; i++)
  {
    R2O_OrderBitReversed(map+i, fft_width, fft_width);
  }
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_GenerateFrequencyPalette(f32 frequency_palette[], u16 wavenumbers_squared[],
                                  f32 fft_meters, u32 num_bands, u32 pal_cnt)
{
  const f32 gravity = 9.80665f;
  const f32 density = 998.200f;
  const f32 tension = 0.07275f;

  f32 scale = 2.0f*PI/fft_meters;
  f32 next_scale = scale * SQRT_OF_2;

  f32 *freq = frequency_palette;
  for (u32 b=0; b<num_bands; b++)
  {
    for (u32 i=0; i<pal_cnt; i++)
    {
      f32 k = scale * Sqrtf((f32)wavenumbers_squared[i]);
      f32 ang_freq = Sqrtf(k*gravity + k*k*k*(tension/density));
      *freq++ = ang_freq * (0.5f/PI);
    }

    f32 temp = scale;
    scale = next_scale;
    next_scale = temp * 2.0f;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_GenerateScalePalette(f32 scale_palette[], u16 wavenumbers_squared[], u32 pal_cnt)
{
  // write palette of wave-number-determined scales
  f32 *scale = scale_palette;
  *scale++ = 0.0f;   // avoid divide-by-zero for 1st entry
  for (u32 i=1; i<pal_cnt; i++)
  {
    *scale++ = 1.0f / Sqrtf((f32)wavenumbers_squared[i]);
  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_InitAnglePalette(f32 angle_palette[], u32 num_bands, u32 pal_cnt)
{
  f32 *ang = angle_palette;
  for (u32 b=0; b<num_bands; b++)
  {
    for (u32 i=0; i<pal_cnt; i++)
    {
      *ang++ = 0.0f;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The ambient spectrum for a given band is just a square array of complex coefficients
// (these can be further scaled by a per-band and/or a global multiplier).
// If a component of the spectrum is set to 1, the corresponding sinusoid has amplitude proportional
// to its wavelength. Each complex number is compressed to 2 bytes.

void R2O_InitAmbientSpectrum(i8 *spectrum, u32 num_bands)
{
  for (u32 b=0; b<num_bands; b++)
  {
    for (i32 i=0; i<32; i++)
    {
      for (i32 j=0; j<32; j++)
      {
        i8 re, im;
        if ((i>=8 && i<24) || (j>=8 && j<24))
        {
          re = (i8)(rand() & 0xFF);
          im = (i8)(rand() & 0xFF);
        }
        else
        {
          re = im = 0;
        }
        *spectrum++ = re;
        *spectrum++ = im;
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_InitAmbient()
{
  u32 num_bands  = g_R2OCon.m_num_lods;
  u32 fft_width  = 32;
  u32 fft_area   = fft_width * fft_width;
  f32 fft_meters = (g_R2OCon.m_step != 0.0f) ? g_R2OCon.m_step * (f32)fft_width : 1.0f;

  u16 wavenumbers_squared[256];
  u32 pal_cnt  = R2O_GenerateWaveMap(wavenumbers_squared, 256, fft_width);
  u32 pal_cnt4 = (pal_cnt+3) & -4;    // rounded up to next multiple of 4

  // allocs
  u8  *map               = (u8  *)g_SystemMemoryAlloc.Alloc(fft_area * sizeof(u8), 128);
  f32 *frequency_palette = (f32 *)g_SystemMemoryAlloc.Alloc(pal_cnt4 * num_bands * sizeof(f32), 128);
  f32 *scale_palette     = (f32 *)g_SystemMemoryAlloc.Alloc(pal_cnt4 * sizeof(f32), 128);
  f32 *angle_palette     = (f32 *)g_SystemMemoryAlloc.Alloc(pal_cnt4 * num_bands * sizeof(f32), 128);
  f32 *phasor_palette    = (f32 *)g_SystemMemoryAlloc.Alloc(pal_cnt4 * num_bands * 2 * sizeof(f32), 128);
  i8  *spectrum          = (i8  *)g_SystemMemoryAlloc.Alloc(fft_area * num_bands * 2 * sizeof(i8), 128);
  f32 *heightmaps        = (f32 *)g_SystemMemoryAlloc.Alloc(fft_area * num_bands * sizeof(f32), 128);

  // generate map & palettes and initialize arrays
  R2O_GenerateIndexMap(map, wavenumbers_squared, fft_width, pal_cnt4);
  R2O_GenerateFrequencyPalette(frequency_palette, wavenumbers_squared, fft_meters, num_bands, pal_cnt4);
  R2O_GenerateScalePalette(scale_palette, wavenumbers_squared, pal_cnt4);
  R2O_InitAnglePalette(angle_palette, num_bands, pal_cnt4);
  R2O_InitAmbientSpectrum(spectrum, num_bands);

  // store constants to g_R2OCon
  g_R2OCon.m_amb_num_bands   = num_bands;
  g_R2OCon.m_amb_fft_width   = fft_width;
  g_R2OCon.m_amb_palette_cnt = pal_cnt4;

  // store ea's to g_R2OCon
  g_R2OCon.m_amb_ea_angle_palette      = (u32)angle_palette;
  g_R2OCon.m_amb_ea_frequency_palette  = (u32)frequency_palette;
  g_R2OCon.m_amb_ea_phasor_palette     = (u32)phasor_palette;
  g_R2OCon.m_amb_ea_wave_scale_palette = (u32)scale_palette;
  g_R2OCon.m_amb_ea_index_map          = (u32)map;
  g_R2OCon.m_amb_ea_spectrum           = (u32)spectrum;
  g_R2OCon.m_amb_ea_heightmaps         = (u32)heightmaps;
}



} // namespace IGG
