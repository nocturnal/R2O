#ifndef __IGR2OAMBIENTSPU_H
#define __IGR2OAMBIENTSPU_H


struct Complex;


void R2O_UpdateAmbient(f32 dt);
void R2O_ChangeAmbientTileOrigin(i32 c0, i32 r0);
void R2O_AddAmbientWaveforms(f32 heightmap[], u32 nc, u32 nr, u32 band);
void R2O_GetAmbientWaves(f32 heightmap[], u32 band);
void R2O_AddAmbientWaves(f32 heightmap[], f32 ambient_waves[], u32 nc, u32 nr, f32 amplitude);

// prototypes for asm functions
extern "C"
{
  void R2O_RotatePhaseAngles(f32 *angles, f32 *frequencies, u32 cnt, f32 dt);
  void R2O_GeneratePhasorPalette(Complex phasors[], f32 angles[], f32 wave_scales[], f32 lod_scale, u32 cnt);
  void R2O_ReconstructPhasors(f32 *x, f32 *y, Complex phasor_palette[], i8 spectrum[], u8 *index_map, u32 cnt);
  void R2O_GenerateCausticsGradients(f32 *p_grad, f32 *p_height, f32 inv_step, u32 cnt);
}

extern u8 *g_IndexMap;
extern i32 c0_amb, r0_amb;


#endif // #ifndef __IGR2OAMBIENTSPU_H

