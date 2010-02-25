#pragma once

       
namespace IGG
{

u32 R2O_GenerateWaveMap(u16 wavenumbers_squared[], u32 wavenum_cnt, u32 fft_width);
void R2O_GenerateIndexMap(u8 map[], u16 wavenumbers_squared[], u32 fft_width, u32 pal_cnt);
void R2O_GenerateFrequencyPalette(f32 frequency_palette[], u16 wavenumbers_squared[],
                                  f32 fft_meters, u32 num_bands, u32 pal_cnt);
void R2O_InitAmbient();

} // namespace IGG

