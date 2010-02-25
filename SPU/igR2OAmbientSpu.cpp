#include "igR2OSpu.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Here is the dma buffering model used for any processes requiring an input-process-output sequence.
// Traditionally this would either be done using 4 buffers (double-buffered input, double-buffered output)
// or 3 buffers via triple buffering.
// In the version here, only 2 buffers are needed due to the use of a dma fence.
// The fence is also preferred over having an explicit processor-controlled syncing of the 3 stages of
// standard triple-buffering. This is because the processor doesn't have to wait for the output buffer to
// be ready before kicking off the output dma: this sync is implicit by the use of a fence.
//
// The notation get(buf,n) will be used to denote a dma-get into buffer #buf (either 0 or 1) from the nth
// group of data in main memory. put(buf,n) is the corresponding dma-put back to main memory, and
// getf(buf,n) is the fenced version of the dma-get.
//
//
// get(0, 0);
// buf = 1;
// 
// for (i=0; i<n; i++)
// {
//   if (i+1 < n)
//     getf(buf, i+1); // get i+1 (with a fence)
//   buf ^= 1;
//   wait(buf);
//   process(buf);     // process i
//   put(buf, i);      // put i
// }
// 
// // and a final wait(buf) if necessary



// globals

u8 *g_IndexMap;
i32 c0_amb, r0_amb;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Rotate a set of phase angles according to the given time-step.
// Each phase angle corresponds to a sinusoid of known wavelength
// which is related to the angular frequency by the dispersion
// relation for water - the angular frequencies have been precomputed.
// A full cycle is represented by the 0.0->1.0 range, and each angle
// is reduced to this range after rotation. All the increments are
// positive, making the range-reduction very straightforward.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if 0
void R2O_RotatePhaseAngles(f32 angles[], f32 frequencies[], u32 cnt, f32 dt)
{
  for (u32 i=0; i<cnt; i++)
  {
    f32 a = angles[i];
    f32 da = frequencies[i] * dt;
    a += da;
    a -= (f32)(i32)a;
    angles[i] = a;
  }
}
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Convert the phase angle palette for a single band into the corresponding
// complex number representations (phasors).
// Also scale each phasor in the band by an lod-determined scale and a wave-number-determined scale.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if 0
void R2O_GeneratePhasorPalette(Complex phasors[], f32 angles[], f32 wave_scales[], f32 lod_scale, u32 cnt)
{
  for (u32 idx=0; idx<cnt; idx++)
  {
    phasors[idx] = lod_scale * wave_scales[idx] * Expi(2.0f*PI*angles[idx]);  // remember the angles use 0->1 for a full cycle
  }
}
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// R2O_UpdatePhases
//
// fenced double-buffered version
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_UpdatePhases(f32 dt)
{
  u32 state = g_Scratchpad.GetState();

  // determine buffer size
  u32 cnt  = g_R2OCon.m_amb_palette_cnt;
  u32 size = cnt * sizeof(f32);

  // single buffered array for palette of wave scales (uses same tag as first buffer load)
  f32 *wave_scales = (f32 *)g_Scratchpad.Alloc(size, 128);
  u32 ea_wave_scales = (u32)g_R2OCon.m_amb_ea_wave_scale_palette;
  DmaGet((volatile void *)wave_scales, ea_wave_scales, size, 0);

  // double buffered arrays
  f32 *angles[2], *frequencies[2];
  Complex *phasors[2];
  for (u32 buf=0; buf<2; buf++)
  {
    angles     [buf] = (f32     *)g_Scratchpad.Alloc(size,   128);
    frequencies[buf] = (f32     *)g_Scratchpad.Alloc(size,   128);
    phasors    [buf] = (Complex *)g_Scratchpad.Alloc(size*2, 128);
  }

  // set main memory addresses
  u32 ea_angles      = (u32)g_R2OCon.m_amb_ea_angle_palette;
  u32 ea_frequencies = (u32)g_R2OCon.m_amb_ea_frequency_palette;
  u32 ea_phasors     = (u32)g_R2OCon.m_amb_ea_phasor_palette;

  // get first buffer load of data
  DmaGet((volatile void *)angles[0],      ea_angles,      size, 0);
  DmaGet((volatile void *)frequencies[0], ea_frequencies, size, 0);

  u32 buf=1;
  u32 num_bands = g_R2OCon.m_amb_num_bands;
  for (u32 i=0; i<num_bands; i++)
  {
    // get i+1, with a fence
    if (i+1 < num_bands)
    {
      DmaGetFence((volatile void *)angles[buf],      ea_angles+size,      size, buf);
      DmaGetFence((volatile void *)frequencies[buf], ea_frequencies+size, size, buf);
    }

    // flip
    buf ^= 1;

    // wait before processing i
    DmaWaitAll(1<<buf);

    // process i
    {
      // rotate each angle according to frequency and timestep
      R2O_RotatePhaseAngles(angles[buf], frequencies[buf], cnt, dt);

      // convert phase angles into scaled phasors (the palette for this waveband)
      f32 lod_scale = powf(0.5f, 0.5*(f32)i) * g_R2OCon.m_step * 0.0625f; // was originally targetted to a 16m step-size
      R2O_GeneratePhasorPalette(phasors[buf], angles[buf], wave_scales, lod_scale, cnt);
    }

    // put i
    DmaPut(angles[buf],  ea_angles,  size,   buf);
    DmaPut(phasors[buf], ea_phasors, size*2, buf);

    // inc main mem addrs
    ea_angles      += size;
    ea_frequencies += size;
    ea_phasors     += size*2;
  }
 
  // restore scratchpad state
  g_Scratchpad.Reset(state);

  // final wait
  DmaWaitAll(1<<buf);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// reconstruct the phasors for a single waveband
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if 0
void R2O_ReconstructPhasors(f32 *x, f32 *y, Complex phasor_palette[], i8 spectrum[], u8 *index_map, u32 cnt)
{
  for (u32 i=0; i<cnt; i++)
  {
    u8 idx = *index_map++;
    Complex spec;
    spec.re = (f32)(*spectrum++) * (1.0f/128.0f);
    spec.im = (f32)(*spectrum++) * (1.0f/128.0f);
    Complex phasor =  spec * phasor_palette[idx];
    *x++ = phasor.re;
    *y++ = phasor.im;
  }
}
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if 0
void R2O_GenerateCausticsGradients(f32 *p_grad, f32 *p_height, f32 inv_step, u32 cnt)
{
  for (u32 i=0; i<cnt; i++)
  {
    p_grad[0] = (p_height[1 ] - p_height[0]) * inv_step;
    p_grad[1] = (p_height[32] - p_height[0]) * inv_step;

    p_grad += 2;
    p_height++;
  }
}
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_FixupCausticsBorder(f32 *p_grad, f32 *p_height, i32 ofs, u32 g_step, u32 h_step, f32 inv_step, u32 cnt)
{
  for (u32 i=0; i<cnt; i++)
  {
    *p_grad = (p_height[ofs] - p_height[0]) * inv_step;

    p_grad   += g_step;
    p_height += h_step;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// generate data for Mark's caustics system
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_GenerateCausticsData(f32 *caustics_buf, Complex *fft_buf, Complex *palette, i8 *spectrum, u32 lod)
{
  // upload caustics envelope
  Complex envelope[16];
  u32 ea = (u32)g_R2OCon.m_ea_caustics_envelope;
  DmaGet((volatile void *)envelope, ea, 16*sizeof(Complex), g_kDmaTag);
  DmaWaitAll(1 << g_kDmaTag);

  // use it to scale the main palette entries
  Complex e, *p=palette;
  for (u32 i=0; i<=16; i++)
  {
    u32 idx = i<=15 ? i : 15;
    e = envelope[idx];
    for (u32 j=0; j<8; j++)
    {
      *p++ *= e;
    }
  }

  // reconstruct the phasors from 'enveloped' palette
  R2O_ReconstructPhasors((f32 *)fft_buf, (f32 *)fft_buf+32*32, palette, spectrum, g_IndexMap, 32*32);

  // inverse fft to get the complex heights
  R2O_FFT2D_sep((f32 *)fft_buf, 1);

  // derive inv_step size from lod
  f32 inv_step = 0.0078125f * powf(2.0f, 0.5f*(f32)lod);

  // convert real parts to differences
  R2O_GenerateCausticsGradients(caustics_buf, (f32 *)fft_buf, inv_step, 32*32);

  // fix up right border x-gradients
  R2O_FixupCausticsBorder(caustics_buf+31*2, (f32 *)fft_buf+31, -31, 32*2, 32, inv_step, 32);

  // fix up bottom border z-gradients
  R2O_FixupCausticsBorder(caustics_buf+31*32*2+1, (f32 *)fft_buf+31*32, -31*32, 2, 1, inv_step, 32);

  // copy the caustics data back to main memory
  ea = (u32)g_R2OCon.m_ea_caustics_data;
  DmaPut(caustics_buf, ea, 1024*sizeof(Complex), g_kDmaTag);
  //DmaWaitAll(1 << g_kDmaTag);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// R2O_GenerateHeightmaps
//
// fenced double-buffered version
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

i64 g_Digits[16] =
{
  0x001C22222222221CULL,  // 0
  0x000818080808081CULL,  // 1
  0x001C22020C10203EULL,  // 2
  0x001C22020C02221CULL,  // 3
  0x00040C14243E0404ULL,  // 4
  0x003E203C0202221CULL,  // 5
  0x000C10203C22221CULL,  // 6
  0x003E020204081010ULL,  // 7
  0x001C22221C22221CULL,  // 8
  0x001C22221E020418ULL,  // 9
  0x00081422223E2222ULL,  // A
  0x003C22223C22223CULL,  // B
  0x001C22202020221CULL,  // C
  0x003C22222222223CULL,  // D
  0x003E20203C20203EULL,  // E
  0x003E20203C202020ULL   // F
};


void DrawDigit(f32 *buf, u32 digit, i16 height)
{
  i64 bits = g_Digits[digit];
  for (u32 i=0; i<8; i++,buf+=24)
  {
    for (u32 j=0; j<8; j++,buf++,bits<<=1)
    {
      *buf = (bits<0) ? height : 0;
    }
  }
}

void R2O_GenerateHeightmaps()
{
  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // buffer sizes
  u32 pal_size  = g_R2OCon.m_amb_palette_cnt * sizeof(Complex);
  u32 spec_size = 32 * 32 * 2 * sizeof(i8);
  u32 fft_size  = 32 * 32 * sizeof(Complex);
  u32 map_size  = 32 * 32 * sizeof(f32);

  // allocate double-buffered memory
  i8 *spectrum[2];
  Complex *palette[2], *fft_buf[2];
  for (u32 buf=0; buf<2; buf++)
  {
    palette [buf] = (Complex *)g_Scratchpad.Alloc(pal_size,  128);
    spectrum[buf] = (i8      *)g_Scratchpad.Alloc(spec_size, 128);
    fft_buf [buf] = (Complex *)g_Scratchpad.Alloc(fft_size,  128);
  }

  // allocate single-buffered memory
  f32 *caustics_buf = (f32 *)g_Scratchpad.Alloc(fft_size,  128);

  // main mem addresses
  u32 ea_pal  = g_R2OCon.m_amb_ea_phasor_palette;
  u32 ea_spec = g_R2OCon.m_amb_ea_spectrum;
  u32 ea_map  = g_R2OCon.m_amb_ea_heightmaps;

  // get first buffer load of data
  DmaGet((volatile void *)palette [0], ea_pal,  pal_size,  0);
  DmaGet((volatile void *)spectrum[0], ea_spec, spec_size, 0);

  // loop over wavebands
  u32 buf=1;
  u32 num_bands = g_R2OCon.m_amb_num_bands;
  for (u32 i=0; i<num_bands; i++)
  {
    // get i+1, with a fence
    if (i+1 < num_bands)
    {
      DmaGetFence((volatile void *)palette [buf], ea_pal +pal_size,  pal_size,  buf);
      DmaGetFence((volatile void *)spectrum[buf], ea_spec+spec_size, spec_size, buf);
    }

    // flip
    buf ^= 1;

    // wait before processing i
    DmaWaitAll(1<<buf);

    // process i
    #if 1
    if ((i&1)==0)
    {
      // reconstruct the full set of phasors for this band
      R2O_ReconstructPhasors((f32 *)fft_buf[buf], (f32 *)fft_buf[buf]+32*32, palette[buf], spectrum[buf], g_IndexMap, 32*32);

      // inverse fft to get the complex heights
      R2O_FFT2D_sep((f32 *)fft_buf[buf], 1);
    }
    else
    {
      // LOD TEMP
      qword *p = (qword *)fft_buf[buf];
      for (u32 j=0; j<512; j++)
      {
        *p++ = (qword){0,0,0,0};
      }
    }
    #else
    qword *p = (qword *)fft_buf[buf];
    for (u32 j=0; j<512; j++)
    {
      *p++ = (qword){0,0,0,0};
    }
    if ((i32)i==g_R2OCon.m_debug_band)
    {
      DrawDigit((f32 *)fft_buf[buf], i>>1, g_R2OCon.m_debug_height);
    }
    #endif

    // put i
    DmaPut((f32 *)fft_buf[buf], ea_map, map_size, buf);

    if (i == g_R2OCon.m_caustics_lod)
    {
      // generate data for Mark's caustics system
      R2O_GenerateCausticsData(caustics_buf, fft_buf[buf], palette[buf], spectrum[buf], i);
    }

    // inc main mem addrs
    ea_pal  += pal_size;
    ea_spec += spec_size;
    ea_map  += map_size;
  }

  // restore scratchpad state
  g_Scratchpad.Reset(state);

  // final wait
  DmaWaitAll(1<<buf | 1<<g_kDmaTag);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_UpdateAmbient(f32 dt)
{
  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  // init fft twiddle factors for ambient system
  u32 fft_width = g_R2OCon.m_amb_fft_width;

  // load the global index map
  u32 fft_area = fft_width * fft_width;
  g_IndexMap = (u8 *)R2O_AllocGetWait(g_R2OCon.m_amb_ea_index_map, 0, fft_area*sizeof(u8));

  // update the phases and generate the heightmaps which get stored out to main memory for the rendering phase
  R2O_UpdatePhases(dt);
  R2O_GenerateHeightmaps();

  // restore scratchpad state
  g_Scratchpad.Reset(state);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_ChangeAmbientTileOrigin(i32 c0, i32 r0)
{
  if (lod & 1)
  {
    i32 c = (c0_amb + c0) & 31;
    i32 r = (r0_amb + r0) & 31;
    c0_amb = (r + c) & 31;
    r0_amb = (r - c) & 31;
  }
  else
  {
    i32 c = (c0_amb + c0) & 31;
    i32 r = (r0_amb + r0) & 31;
    c0_amb = (c - r) & 31;
    r0_amb = (c + r) & 31;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_GetAmbientWaves(f32 heightmap[], u32 band)
{
  u32 size = 32*32*sizeof(float);
  u32 ea   = g_R2OCon.m_amb_ea_heightmaps + (band * size);
  DmaGet((volatile void *)heightmap, ea, size, g_kDmaTag);
  DmaWaitAll(1 << g_kDmaTag);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void R2O_AddAmbientWaves(f32 heightmap[], f32 ambient_waves[], u32 nc, u32 nr, f32 amplitude)
{
  // accumulate into to the heightmap heights
  for (u32 r=0,i=r0_amb; r<nr; r++,i=(i+1)&31)
  {
    for (u32 c=0,j=c0_amb; c<nc; c++,j=(j+1)&31)
    {
      u32 rc = r*nc+c;
      u32 ij = i*32+j;
      heightmap[rc] += amplitude * ambient_waves[ij];
    }
  }
}

