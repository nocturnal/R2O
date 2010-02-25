#pragma once

extern "C"
{
  void R2O_RefinementCopyMajorEven(f32 *dest, u32 dc, u32 dr, f32 *src, u32 src_stride);
  void R2O_RefinementCopyMinorEven(f32 *dest, u32 dc, u32 dr, f32 *src, u32 src_stride);
  void R2O_RefinementCopyMajorOdd (f32 *dest, u32 dc, u32 dr, f32 *src, u32 src_stride);
  void R2O_RefinementCopyMinorOdd (f32 *dest, u32 dc, u32 dr, f32 *src, u32 src_stride);
  void R2O_InterpolateDDMinusLinear(f32 *data, u32 nc, u32 nr, f32 DD_coeff0);
  void R2O_AddDD(f32 *dest, f32 *src, u32 nc, u32 nr);
  void R2O_CopyDD(f32 *dest, f32 *src, u32 nc, u32 nr);
  void R2O_InterpLinear(f32 *dest, f32 *src, u32 nc, u32 nr);
  void R2O_Geomorph(f32 *dest, f32 *full_morph, f32 *src, u32 nc, u32 nr,
                       vf32 blend0, vf32 dblend_c, vf32 dblend_r, vf32 dblend_y);
  void R2O_ReplicateAmbient(f32 *dest, u32 dc, u32 dr, f32 *amb, i32 c_amb, i32 r_amb, f32 amplitude);
  u32  R2O_GenerateFans(u16 fans[], u16 indices[], u32 dc, u32 dr, u32 stride);
  u32  R2O_AssignIndices(u16 indices[], u8 outcodes[], u32 cnt);
  void R2O_FlagVertsRender(u8 outcodes[], u32 cnt, u32 stride);
  void R2O_FlagVertsKeep(u8 outcodes[], u32 cnt, u32 stride);
  void R2O_FlagQuadsRenderKeep(u8 outcodes[], i32 nc, i32 nr);
  void R2O_AccumulateNormalMaps(f32* dest, f32* lo, f32* hi, f32 weight, i32 nc_dest, i32 nr_dest);
  void R2O_GenerateMapDerivs(i8 *dest, f32 *src, i32 cols, i32 rows, u32 src_stride, u32 dest_stride, u32 wrap, f32 scale);
  void R2O_GenerateVertexDerivs(i16 *dest, f32 *src, vf32 basis_col_x, vf32 basis_col_z, vf32 basis_row_x, vf32 basis_row_z,
                                u32 nc, u32 nr);



  void CompressDerivsForOutputAsm(i16 *dest, i16 *src, u16 *indices, u32 cnt);
  void CompressVertsForOutputAsm(u8 *dest, f32 *src, f32 scale, u16 *indices, i32 cols, i32 rows);
  void InterpolateOutcodesHorizontalAsm(u8 outcodes[], u32 nc, u32 nr);
  void InterpolateOutcodesVerticalEvenAsm(u8 outcodes[], u32 nc, u32 nr);
  void InterpolateOutcodesVerticalOddAsm(u8 outcodes[], u32 nc, u32 nr);
  void GenerateOutcodesAsm(u8 *outcodes, f32 *heights, vf32 origin, vf32 dvc, vf32 dvr, vf32 *planes, u32 nc, u32 nr);
  void SetNearBitsAsm(u8 *outcodes, u32 nc, u32 nr);
  void GetRefinementWindowAsm(i16 *p_min_c_dest, i16 *p_max_c_dest, i16 *p_min_r_dest, i16 *p_max_r_dest,
                              u8 *outcodes, u32 nc_src, u32 nr_src, u32 lod);
  void R2O_PadOutcodesAsm (u8 *dst, u8 *src, u32 src_width, u32 src_height);
  void R2O_PackOutcodesAsm(u8 *dst, u32 dst_width, u32 dst_height, u8 *src);
  void RefineOutcodesEvenAsm(u8 *dst, u32 dst_width, u32 dst_height, u8 *src, u32 src_width, u32 src_height, i32 c_orig, i32 r_orig);
  void RefineOutcodesOddAsm (u8 *dst, u32 dst_width, u32 dst_height, u8 *src, u32 src_width, u32 src_height, i32 c_orig, i32 r_orig);
}
