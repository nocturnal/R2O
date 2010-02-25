      
#include "igR2OSpu.h"


u32 g_NormalMapPointers[MAX_LODS] ALIGNED(16);


// ------------------------------------------------------------------------------------------------------------------------


void R2O_GenerateAmbientNormalMaps()
{
  // clear ea's
  for (u32 i=0; i<MAX_LODS; i++)
  {
    g_NormalMapPointers[i] = 0;
  }

  // save scratchpad state
  u32 state = g_Scratchpad.GetState();

  f32 amplitude = 0.2f; // hmm

  // allocate buffers
  f32 *heights    = (f32 *)g_Scratchpad.Alloc(33*512);  // 64x66 floats (2 rows duplicated)
  f32 *fft_buffer = (f32 *)g_Scratchpad.Alloc(8*1024);  // 32x32 floats x 2
  i8  *normal_map = (i8  *)g_Scratchpad.Alloc(8*1024);  // 64x64x2 i8's

  heights += 64;                                        // step past the duplicated top row
  vf32 *heights_vf32 = (vf32 *)heights;

  // set base pointers
  f32 *ffts[2];
  ffts[0] = fft_buffer + 0 * 32*32;
  ffts[1] = fft_buffer + 1 * 32*32;

  // get heightmaps 2 & 4
  u32 size = 32*32*sizeof(float);
  for (u32 i=0; i<2; i++)
  {
    u32 ea = g_R2OCon.m_amb_ea_heightmaps + size*2*(i+1);
    DmaGet((volatile void *)ffts[i], ea, size, 0);
  }
  
  // make 27 normal maps
  u32 norm_lod;
  for (norm_lod=0; norm_lod<MAX_LODS-4; norm_lod+=2)
  {
    f32 inv_step = powf(2.0f, 0.5f*(f32)(norm_lod+4)) * (1.0f/32.0f);  // map is 4 levels up from the geometry, and 32 is the fft size
    inv_step *= 16.0f / g_R2OCon.m_step;    // was originally targeted at 16m step-size

    // sum the 2 levels
    DmaWaitAll(1<<0);
    R2O_AccumulateNormalMaps(heights, ffts[0], ffts[1], amplitude, 64,  64);

    // get the heightmap for next level up
    u32 size = 32*32*sizeof(float);
    u32 ea   = g_R2OCon.m_amb_ea_heightmaps + size*(norm_lod+6);
    DmaGet((volatile void *)ffts[0], ea, size, 0);

    // duplicate the top row at the bottom of the heightmap, and the bottom row at the top
    // (which considerably simplifies the work done by GenerateMapDerivs)
    for (u32 i=0; i<16; i++)
    {
      heights_vf32[i-16]    = heights_vf32[i+63*16];
      heights_vf32[i+64*16] = heights_vf32[i];
    }

    // convert heights to derivatives
    f32 scale = inv_step * 128.0f;
    DmaWaitAll(1<<1);
    R2O_GenerateMapDerivs(normal_map, heights, 64, 64, 64, 128, true, scale);

    // allocate gpu memory for the normal map and pass back the ea
    size = 64*64*2;
    ea = IGG::DynamicRenderAlloc(g_R2OCon.m_allocator, size, 128, false);
    g_NormalMapPointers[norm_lod] = ea;

    // upload derivs to main mem
    DmaPut((volatile void *)normal_map, ea, size, 1);

    // cycle fft pointers
    f32 *tempf = ffts[0];
    ffts[0]    = ffts[1];
    ffts[1]    = tempf;
  }

  // wait for dma's, because we're about to free the buffers they're using
  DmaWaitAll(1<<0 | 1<<1);

  // restore scratchpad state
  g_Scratchpad.Reset(state);


  // pass back the normal map base pointers
  DmaPut((volatile void *)g_NormalMapPointers, g_R2OCon.m_ea_normal_map_pointers, sizeof(g_NormalMapPointers), g_kDmaTag);

  // wait for dma's
  DmaWaitAll(1 << g_kDmaTag);
}

// ------------------------------------------------------------------------------------------------------------------------


