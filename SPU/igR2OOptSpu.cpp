
#include "igR2OSpu.h"



// globals

DEF_SHUF(0000); DEF_SHUF(D000); DEF_SHUF(C000); DEF_SHUF(CD00);
DEF_SHUF(B000); DEF_SHUF(BD00); DEF_SHUF(BC00); DEF_SHUF(BCD0);
DEF_SHUF(A000); DEF_SHUF(AD00); DEF_SHUF(AC00); DEF_SHUF(ACD0);
DEF_SHUF(AB00); DEF_SHUF(ABD0); DEF_SHUF(ABC0); DEF_SHUF(ABCD);
DEF_SHUF(AAAA); DEF_SHUF(BBBB); DEF_SHUF(CCCC); DEF_SHUF(DDDD);
DEF_SHUF(BCDa); DEF_SHUF(CDab); DEF_SHUF(Dabc); DEF_SHUF(AABB);
DEF_SHUF(ABBC); DEF_SHUF(CDDa); DEF_SHUF(CCaa); DEF_SHUF(ACCa);
DEF_SHUF(CCDD); DEF_SHUF(BDDb); DEF_SHUF(aAcC); DEF_SHUF(ABab);
DEF_SHUF(CDcd); DEF_SHUF(DDbb); DEF_SHUF(Dbbd); DEF_SHUF(BBDD);
DEF_SHUF(ACac); DEF_SHUF(BDbd); DEF_SHUF(AaBb); DEF_SHUF(CcDd);
DEF_SHUF(AaaA); DEF_SHUF(BaDc);

DEF_SHUF(EFGHabcd); DEF_SHUF(b000aABC); DEF_SHUF(AcbaEFGH); DEF_SHUF(CcXXXXXX);
DEF_SHUF(BBBBBBBB); DEF_SHUF(d0cCDEe0); DEF_SHUF(AcCDEFGe); DEF_SHUF(dcXXXXXX);
DEF_SHUF(DDDDDDDD); DEF_SHUF(dCDEe000); DEF_SHUF(ABCDEedc); DEF_SHUF(cCXXXXXX);
DEF_SHUF(bCc000aA); DEF_SHUF(ABCcbaGH); DEF_SHUF(BCXXXXXX); DEF_SHUF(HHHHHHHH);
DEF_SHUF(BbDdFfHh); DEF_SHUF(BDFHbdfh); DEF_SHUF(AaCcEeGg); DEF_SHUF(ACEGaceg);         
DEF_SHUF(ABCDEFGH); DEF_SHUF(BCDEFGHa); DEF_SHUF(BFbf0000); DEF_SHUF(DHdh0000);
DEF_SHUF(0000BFbf); DEF_SHUF(0000DHdh); DEF_SHUF(BFbfBFbf); DEF_SHUF(DHdhDHdh);
DEF_SHUF(A0a0B0b0); DEF_SHUF(C0c0D0d0); DEF_SHUF(E0e0F0f0); DEF_SHUF(G0g0H0h0);
DEF_SHUF(AEae0000); DEF_SHUF(CGcg0000); DEF_SHUF(0000AEae); DEF_SHUF(0000CGcg);
DEF_SHUF(AEaeAEae); DEF_SHUF(CGcgCGcg); DEF_SHUF(A0B0C0D0); DEF_SHUF(E0F0G0H0);
DEF_SHUF(acegaceg);

DEF_SHUF(0A0B0C0D0E0F0G0H); DEF_SHUF(0I0J0K0L0M0N0O0P); DEF_SHUF(AABBCCDDEEFFGGHH); DEF_SHUF(IIJJKKLLMMNNOOPP);
DEF_SHUF(BCDEFGHIJKLMNOPa); DEF_SHUF(IJKLMNOPabcdefgh); DEF_SHUF(JKLMNOPabcdefghi); DEF_SHUF(BbDdFfHhJjLlNnPp);
DEF_SHUF(AaCcEeGgIiKkMmOo); DEF_SHUF(DCcdHGghLKklPOop); DEF_SHUF(AbCdEfGhIjKlMnOp); DEF_SHUF(AEIMAEIMAEIMAEIM);
DEF_SHUF(ABCD000000000000); DEF_SHUF(AaCcEeGgIJKLMNOP); DEF_SHUF(abcdefghIjKlMnOp); DEF_SHUF(Pabcdefghijklmno);
DEF_SHUF(HIJKLMNOPabcdefg); DEF_SHUF(abcdefghijklmnop); DEF_SHUF(ABCDEFGHIJKLMNOP); DEF_SHUF(Aa00Ee00Ii00Mm00);
DEF_SHUF(AaAeEeEeIiIiMmMm); DEF_SHUF(0Aa00Ee00Ii00Mm0); DEF_SHUF(a00Ae00Ai00Im00M); DEF_SHUF(0aA00eE00iI00mM0);
DEF_SHUF(aBCDeFGHiJKLmNOP); DEF_SHUF(ABCaEFGeIJKiMNOm); DEF_SHUF(A0a0E0e0I0i0M0m0); DEF_SHUF(ACacEGegIKikMOmo);

vu32 mask_F000 = (vu32){ 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000 };
vu32 mask_000F = (vu32){ 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF };
vu32 mask_F0F0 = (vu32){ 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0x00000000 };
vu32 mask_0F0F = (vu32){ 0x00000000, 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF };
vu32 mask_FFFF = (vu32){ 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF };
vu32 mask_00FF = (vu32){ 0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF };
vu32 mask_FF00 = (vu32){ 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000 };

vu8 vu8_0x40 = (vu8) { 0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40 };

vf32 zero    = (vf32){0.00f, 0.00f, 0.00f, 0.00f};
vf32 half    = (vf32){0.50f, 0.50f, 0.50f, 0.50f};
vf32 quarter = (vf32){0.25f, 0.25f, 0.25f, 0.25f};


//------------------------------------------------------------------------------------------------------------------------

/*            /\
             /..\
            /    \
    a0  a1 /a2  a3>
          /      /
         /..  ../               ->    (b0,b1,a2,a3)
        /      /
       <b0  b1/ b2  b3
        \    /
         \../
          \/
*/

#if 0
void R2O_RefinementCopyMajorEven(f32 *dest, u32 dc, u32 dr, f32 *src, u32 src_stride)
{
  // check alignments
  IG_ASSERT(((u32)dest  & 0xF) == 0);
  IG_ASSERT(((u32)src   & 0x3) == 0);
  IG_ASSERT((dc         & 0x3) == 0);
  IG_ASSERT((src_stride & 0x3) == 0);

  vf32 bL,bR,a03,b03=spu_splats(0.0f), out;
  u8 byte;
  vu8 shuf;

  u32 qdc = dc >> 2;

  // loop over destination columns
  for (u32 c=0; c<dc; c+=4)
  {
    vf32 *p_dest    = (vf32 *)&dest[c];
    f32  *p_src_f32 = &src[-(c>>1)*(src_stride-1)];

    // loop over destination rows
    for (u32 r=0; r<dr; r+=2)
    {
      // advance queue
      a03 = b03;

      // get alignment shuffle
      byte = (u32)p_src_f32 & 0xF;
      shuf = spu_splats(byte);
      shuf = (vu8)spu_add((vu32)shuf, (vu32)shuf_ABCD);

      // get next row of input
      bL = *(vf32 *)(p_src_f32+0);
      bR = *(vf32 *)(p_src_f32+4);
      b03 = spu_shuffle(bL, bR, shuf);

      // shuffle
      out = spu_sel(a03, b03, mask_FF00);

      // output
      *p_dest = out;

      // step pointers
      p_src_f32 += src_stride + 1;
      p_dest += qdc;
    }
  }
}
#endif

//------------------------------------------------------------------------------------------------------------------------

/*
              /\
  a0  a1  a2 /a3\
            /    \
           /..  ..>
          /      /
      b0 /b1  b2/ b3            ->    (c0,b1,a2,a3)
        /      /
       <..  ../
        \    /
         \c0/ c1  c2  c3
          \/
*/


#if 0
void R2O_RefinementCopyMinorEven(f32 *dest, u32 dc, u32 dr, f32 *src, u32 src_stride)
{
  // check alignments
  IG_ASSERT(((u32)dest  & 0xF) == 0);
  IG_ASSERT(((u32)src   & 0x3) == 0);
  IG_ASSERT((dc         & 0x3) == 0);
  IG_ASSERT((src_stride & 0x3) == 0);

  vf32 cL,cR,a03,b03=spu_splats(0.0f),c03,out;
  u8 byte;
  vu8 shuf;

  u32 qdc = dc >> 2;

  // loop over destination columns
  for (u32 c=0; c<dc; c+=4)
  {
    vf32 *p_dest    = (vf32 *)&dest[c];
    f32  *p_src_f32 = &src[-(c>>1)*(src_stride-1)] - (src_stride + 1);

    // get alignment shuffle
    byte = (u32)p_src_f32 & 0xF;
    shuf = spu_splats(byte);
    shuf = (vu8)spu_add((vu32)shuf, (vu32)shuf_ABCD);

    // initialize input queue
    cL = *(vf32 *)(p_src_f32+0);
    cR = *(vf32 *)(p_src_f32+4);
    c03 = spu_shuffle(cL, cR, shuf);

    // step source pointer
    p_src_f32 += src_stride + 1;

    // loop over destination rows
    for (u32 r=0; r<dr; r+=2)
    {
      // advance queue
      a03 = b03;
      b03 = c03;

      // get alignment shuffle
      byte = (u32)p_src_f32 & 0xF;
      shuf = spu_splats(byte);
      shuf = (vu8)spu_add((vu32)shuf, (vu32)shuf_ABCD);

      // get next row of input
      cL = *(vf32 *)(p_src_f32+0);
      cR = *(vf32 *)(p_src_f32+4);
      c03 = spu_shuffle(cL, cR, shuf);

      // shuffle
      out = spu_sel(b03, c03, mask_F000);
      out = spu_sel(out, a03, mask_000F);

      // output
      *p_dest = out;

      // step pointers
      p_src_f32 += src_stride + 1;
      p_dest += qdc;
    }
  }
}
#endif


//------------------------------------------------------------------------------------------------------------------------

/*        /\
         /a0\ a1  a2  a3
        /    \
       <..  ..\
        \      \ 
      b0 \b1  b2\ b3            ->    (a0,b1,b2,c3)
          \      \ 
           \..  ..>
            \    /
  c0  c1  c2 \c3/
              \/
              
*/


#if 0
void R2O_RefinementCopyMajorOdd(f32 *dest, u32 dc, u32 dr, f32 *src, u32 src_stride)
{
  // check alignments
  IG_ASSERT(((u32)dest  & 0xF) == 0);
  IG_ASSERT(((u32)src   & 0x3) == 0);
  IG_ASSERT((dc         & 0x3) == 0);
  IG_ASSERT((src_stride & 0x3) == 0);

  vf32 cL,cR,a03,b03=spu_splats(0.0f),c03,out;
  u8 byte;
  vu8 shuf;

  u32 qdc = dc >> 2;

  // loop over destination columns
  for (u32 c=0; c<dc; c+=4)
  {
    vf32 *p_dest    = (vf32 *)&dest[c];
    f32  *p_src_f32 = &src[(c>>1)*(src_stride+1)] + (src_stride - 1);

    // get alignment shuffle
    byte = (u32)p_src_f32 & 0xF;
    shuf = spu_splats(byte);
    shuf = (vu8)spu_add((vu32)shuf, (vu32)shuf_ABCD);

    // initialize input queue
    cL = *(vf32 *)(p_src_f32+0);
    cR = *(vf32 *)(p_src_f32+4);
    c03 = spu_shuffle(cL, cR, shuf);

    // step source pointer
    p_src_f32 += src_stride - 1;

    // loop over destination rows
    for (u32 r=0; r<dr; r+=2)
    {
      // advance queue
      a03 = b03;
      b03 = c03;

      // get alignment shuffle
      byte = (u32)p_src_f32 & 0xF;
      shuf = spu_splats(byte);
      shuf = (vu8)spu_add((vu32)shuf, (vu32)shuf_ABCD);

      // get next row of input
      cL = *(vf32 *)(p_src_f32+0);
      cR = *(vf32 *)(p_src_f32+4);
      c03 = spu_shuffle(cL, cR, shuf);

      // shuffle
      out = spu_sel(b03, a03, mask_F000);
      out = spu_sel(out, c03, mask_000F);

      // output
      *p_dest = out;

      // step pointers
      p_src_f32 += src_stride - 1;
      p_dest += qdc;
    }
  }
}
#endif


//------------------------------------------------------------------------------------------------------------------------

/*        /\
         /..\
        /    \
       <a0  a1\ a2  a3
        \      \ 
         \..  ..\               ->    (a0,a1,b2,b3)
          \      \ 
    b0  b1 \b2  b3>
            \    /
             \../
              \/
              
*/


#if 0
void R2O_RefinementCopyMinorOdd(f32 *dest, u32 dc, u32 dr, f32 *src, u32 src_stride)
{
  // check alignments
  IG_ASSERT(((u32)dest  & 0xF) == 0);
  IG_ASSERT(((u32)src   & 0x3) == 0);
  IG_ASSERT((dc         & 0x3) == 0);
  IG_ASSERT((src_stride & 0x3) == 0);

  vf32 bL,bR,a03,b03=spu_splats(0.0f),out;
  u8 byte;
  vu8 shuf;

  u32 qdc = dc >> 2;

  // loop over destination columns
  for (u32 c=0; c<dc; c+=4)
  {
    vf32 *p_dest    = (vf32 *)&dest[c];
    f32  *p_src_f32 = &src[(c>>1)*(src_stride+1)] + (src_stride - 1);

    // loop over destination rows
    for (u32 r=0; r<dr; r+=2)
    {
      // advance queue
      a03 = b03;

      // get alignment shuffle
      byte = (u32)p_src_f32 & 0xF;
      shuf = spu_splats(byte);
      shuf = (vu8)spu_add((vu32)shuf, (vu32)shuf_ABCD);

      // get next row of input
      bL = *(vf32 *)(p_src_f32+0);
      bR = *(vf32 *)(p_src_f32+4);
      b03 = spu_shuffle(bL, bR, shuf);

      // shuffle
      out = spu_sel(a03, b03, mask_00FF);

      // output
      *p_dest = out;

      // step pointers
      p_src_f32 += src_stride - 1;
      p_dest += qdc;
    }
  }
}
#endif


//------------------------------------------------------------------------------------------------------------------------

#if 0
void R2O_RefinementCopyMajor(f32 *dest, u32 dc, u32 dr, f32 *src, u32 src_stride)
{
  if (lod & 1)
  {
    R2O_RefinementCopyMajorOdd(dest, dc, dr, src, src_stride);
  }
  else
  {
    R2O_RefinementCopyMajorEven(dest, dc, dr, src, src_stride);
  }
}
#endif

//------------------------------------------------------------------------------------------------------------------------

#if 0
void R2O_RefinementCopyMinor(f32 *dest, u32 dc, u32 dr, f32 *src, u32 src_stride)
{
  if (lod & 1)
  {
    R2O_RefinementCopyMinorOdd(dest, dc, dr, src, src_stride);
  }
  else
  {
    R2O_RefinementCopyMinorEven(dest, dc, dr, src, src_stride);
  }
}
#endif

//--------------------------------------------------------------------------------------------------------------------------

#if 0
void R2O_InterpolateDDMinusLinear(f32 *data, u32 nc, u32 nr, f32 DD_coeff0)
{
  // check alignments
  IG_ASSERT(((u32)data  & 0xF) == 0);
  IG_ASSERT((nc         & 0x3) == 0);

  vf32 a03, a47, b03, b47, c03, c47, d03, d47, e03, e47, v03, v14, v25, v36, v47, h03, j03, k03, l03, out;

  // 4-tap DD-interpolation coeffs
  vf32 p = spu_splats(DD_coeff0);
  vf32 q = spu_splats(0.5f - DD_coeff0);

  u32  qstride = nc >> 2;
  u32  offset0 = 3*qstride;
  u32  offset1 = 3*qstride + 1;

  for (u32 c=0; c<nc; c+=4)
  {
    vf32 *p_data = (vf32 *)&data[c];

    // initialize input queue
    c03 = p_data[0];
    c47 = p_data[1];
    d03 = p_data[qstride];
    d47 = p_data[qstride+1];
    e03 = p_data[2*qstride];
    e47 = p_data[2*qstride+1];
    l03 = spu_shuffle(d03, d47, (vu8)shuf_BDDb);

    for (u32 r=0; r<nr; r+=2)
    {
      // advance input queue
      a03 = c03;
      a47 = c47;
      b03 = d03;
      b47 = d47;
      c03 = e03;
      c47 = e47;
      j03 = l03;

      // even row

      // get 2 qwords of input
      d03 = p_data[offset0];
      d47 = p_data[offset1];

      // apply vertical kernel
      v03 = p*a03 + q*b03 + q*c03 + p*d03;
      v47 = p*a47 + q*b47 + q*c47 + p*d47;
      v14 = spu_shuffle(v03, v47, (vu8)shuf_BCDa);
      v25 = spu_shuffle(v03, v47, (vu8)shuf_CDab);
      v36 = spu_shuffle(v03, v47, (vu8)shuf_Dabc);

      // apply horizontal kernel
      h03 = p*v03 + q*v14 + q*v25 + p*v36;

      // shuffles for linear interp
      //j03 = spu_shuffle(b03, b47, (vu8)shuf_BDDb);
      k03 = spu_shuffle(c03, c47, (vu8)shuf_CCaa);

      // subtract linear interpolant
      out = h03 - half * (j03+k03);

      // store 1 word of output
      p_data[0] = out;

      // step pointer
      p_data += qstride;

      // odd row

      // get 2 qwords of input
      e03 = p_data[offset0];
      e47 = p_data[offset1];

      // apply vertical kernel
      v03 = p*b03 + q*c03 + q*d03 + p*e03;
      v47 = p*b47 + q*c47 + q*d47 + p*e47;
      v14 = spu_shuffle(v03, v47, (vu8)shuf_BCDa);
      v25 = spu_shuffle(v03, v47, (vu8)shuf_CDab);
      v36 = spu_shuffle(v03, v47, (vu8)shuf_Dabc);

      // apply horizontal kernel
      h03 = p*v03 + q*v14 + q*v25 + p*v36;

      // shuffles for linear interp
      l03 = spu_shuffle(d03, d47, (vu8)shuf_BDDb);

      // subtract linear interpolant
      out = h03 - half * (k03+l03);

      // store 1 word of output
      p_data[0] = out;

      // step pointer
      p_data += qstride;
    }
  }
}
#endif


//--------------------------------------------------------------------------------------------------------------------------

#if 0
void R2O_AddDD(f32 *dest, f32 *src, u32 nc, u32 nr)
{
  // check alignments
  IG_ASSERT(((u32)dest & 0xF) == 0);
  IG_ASSERT(((u32)src  & 0xF) == 0);
  IG_ASSERT((nc        & 0x3) == 0);

  vf32 out0, out1, dd, dd0, dd1;

  u32  qstride = nc >> 2;

  for (u32 c=0; c<nc; c+=4)
  {
    vf32 *p_source = (vf32 *)&src[c];
    vf32 *p_dest   = (vf32 *)&dest[c];

    for (u32 r=0; r<nr; r+=2)
    {
      // read 1 qword from source array, 2 from dest array
      out0 = p_dest[0];
      out1 = p_dest[qstride];

      // configure dd interpolants
      dd   = p_source[0];
      dd0  = spu_sel(dd, zero, mask_F0F0);
      dd1  = spu_sel(dd, zero, mask_0F0F);

      // add
      out0 += dd0;
      out1 += dd1;

      // store back to dest array
      p_dest[0]       = out0;
      p_dest[qstride] = out1;

      // step pointers
      p_source += qstride;
      p_dest   += qstride * 2;
    }
  }
}
#endif


//--------------------------------------------------------------------------------------------------------------------------


#if 0
void R2O_InterpLinear(f32 *dest, f32 *src, u32 nc, u32 nr)
{
  // check alignments
  IG_ASSERT(((u32)dest & 0xF) == 0);
  IG_ASSERT(((u32)src  & 0xF) == 0);
  IG_ASSERT((nc        & 0x3) == 0);

  vf32 a03, a47, b03, a0a2a2a4, b0a0b2a2, avg, out0, out1;

  u32  qstride = nc >> 2;

  vf32 *p_source = (vf32 *)src;
  vf32 *p_dest   = (vf32 *)dest;

  for (u32 r=0; r<nr; r+=2)
  {
    // initialize queue
    a47 = p_source[0];

    for (u32 c=0; c<nc; c+=4)
    {
      // advance queue
      a03 = a47;

      // read 2 qwords from source
      a47 = p_source[1];
      b03 = p_source[qstride];

      // shuffle
      a0a2a2a4 = spu_shuffle(a03, a47, (vu8)shuf_ACCa);
      b0a0b2a2 = spu_shuffle(a03, b03, (vu8)shuf_aAcC);

      // average
      avg = half * (a0a2a2a4 + b0a0b2a2);

      // interleave
      out0 = spu_sel(a03, avg, mask_0F0F);
      out1 = spu_sel(a03, avg, mask_F0F0);

      // store
      p_dest[0]       = out0;
      p_dest[qstride] = out1;

      // advance pointers
      p_source++;
      p_dest++;
    }

    p_dest += qstride;
  }
}
#endif


//--------------------------------------------------------------------------------------------------------------------------

#if 0
void R2O_Geomorph(f32 *dest, f32 *full_morph, f32 *src, u32 nc, u32 nr, f32 z0, f32 z1)
{
  // check alignments
  IG_ASSERT(((u32)dest & 0xF) == 0);
  IG_ASSERT(((u32)src  & 0xF) == 0);
  IG_ASSERT((nc        & 0x3) == 0);

  f32  scale    = 1.0f / (z1-z0);
  vf32 blend0   = spu_splats((spu_extract(origin_camera,2) - z0) * scale);
  vf32 dblend_c = spu_splats(spu_extract(dvc_camera,2) * scale);
  vf32 dblend_r = spu_splats(spu_extract(dvr_camera,2) * scale);
  vf32 dblend_y = spu_splats(spu_extract(*(vf32 *)&g_pViewData->m_world_to_camera_matrix.m_v1,2) * scale);

  blend0 = blend0 + (vf32){0,1,2,3} * dblend_c;
  dblend_c = dblend_c * spu_splats(4.0f);

  vf32 y, dy, blend, y_morphed, y_full;
  vu32 blend_vu32;

  u32  qstride = nc >> 2;

  for (u32 c=0; c<nc; c+=4)
  {
    vf32 *p_src  = (vf32 *)&src[c];
    vf32 *p_dest = (vf32 *)&dest[c];
    vf32 *p_full = (vf32 *)&full_morph[c];

    vf32 grid_blend = blend0;

    for (u32 r=0; r<nr; r++)
    {
      // read 1 qword from source array, 1 from dest array
      y  = *p_src;
      dy = *p_dest;

      // compute blend factors
      blend = grid_blend + y * dblend_y;

      // clamp blend factors
      blend_vu32 = spu_convtu(blend, 32);
      blend = spu_convtf(blend_vu32, 32);

      // morph
      y_morphed = y + blend * dy;
      *p_dest   = y_morphed;

      // full morph
      y_full  = y + dy;
      *p_full = y_full;

      // step pointers
      p_src  += qstride;
      p_dest += qstride;
      p_full += qstride;

      // step grid blends
      grid_blend += dblend_r;
    }

    blend0 += dblend_c;
  }
}
#endif


//--------------------------------------------------------------------------------------------------------------------------

#if 0
void R2O_ReplicateAmbient(f32 *dest, u32 dc, u32 dr, f32 *amb, i32 c_amb, i32 r_amb, f32 amplitude)
{
  // check alignments
  IG_ASSERT(((u32)dest  & 0xF) == 0);
  IG_ASSERT((dc         & 0x3) == 0);
  IG_ASSERT((c_amb      & 0x3) == 0);

  u32 qdc = dc >> 2;
  vf32 *p_amb = (vf32 *)amb;
  vf32 qamplitude = spu_splats(amplitude);

  // loop over destination columns
  for (u32 c=0; c<dc; c+=4)
  {
    i32 idx_amb = (r_amb*32+c_amb) >> 2;

    vf32 *p_dest    = (vf32 *)&dest[c];

    // loop over destination rows
    for (u32 r=0; r<dr; r++)
    {
      // copy 1 qword, with scaling
      *p_dest = p_amb[idx_amb] * qamplitude;

      // step pointer and index
      p_dest += qdc;
      idx_amb = (idx_amb + 8) & 0xFF;
    }

    c_amb = (c_amb + 4) & 31;
  }
}
#endif

//--------------------------------------------------------------------------------------------------------------------------

#if 0
u32 R2O_GenerateFans(u16 fans[], u16 indices[], u32 dc, u32 dr, u32 stride)
{
  u32  u0, u1, u2, u3, u4;
  vu16 L0, L1, L2, L3, L4;
  vu16 R0, R1, R2, R3, R4;
  vu8  shufE, shufO;

  vu16 *p_dest = (vu16 *)fans;

  vu16 indices0, indices1, indices2, indices3, indices4;

  vu16 fan0, fan1;
  vu16 cmp0, cmp1;
  vu16 center;
  u32  dest_inc;

  for (u32 r=0; r<dr; r+=4)
  {
    u0 = (u32)indices + 2*(r+0)*stride;
    u1 = (u32)indices + 2*(r+1)*stride;
    u2 = (u32)indices + 2*(r+2)*stride;
    u3 = (u32)indices + 2*(r+3)*stride;
    u4 = (u32)indices + 2*(r+4)*stride;

    for (u32 c=0; c<dc; c+=4)
    {
      // get 5 rows of indices
      L0 = ((vu16 *)u0)[0];
      L1 = ((vu16 *)u1)[0];
      L2 = ((vu16 *)u2)[0];
      L3 = ((vu16 *)u3)[0];
      L4 = ((vu16 *)u4)[0];

      R0 = ((vu16 *)u0)[1];
      R1 = ((vu16 *)u1)[1];
      R2 = ((vu16 *)u2)[1];
      R3 = ((vu16 *)u3)[1];
      R4 = ((vu16 *)u4)[1];

      shufE = (u0 & 8) ? (vu8)shuf_EFGHabcd : (vu8)shuf_ABCDEFGH;
      shufO = (u1 & 8) ? (vu8)shuf_EFGHabcd : (vu8)shuf_ABCDEFGH;

      indices0 = spu_shuffle(L0, R0, shufE);
      indices1 = spu_shuffle(L1, R1, shufO);
      indices2 = spu_shuffle(L2, R2, shufE);
      indices3 = spu_shuffle(L3, R3, shufO);
      indices4 = spu_shuffle(L4, R4, shufE);

      u0 += 8;
      u1 += 8;
      u2 += 8;
      u3 += 8;
      u4 += 8;


      /*
      
      a--b--c--d--e
      |\ | /|\ | /|
      | \|/ | \|/ |
      f--g--h--i--j
      | /|\ | /|\ |
      |/ | \|/ | \|
      k--l--m--n--o
      |\ | /|\ | /|
      | \|/ | \|/ |
      p--q--r--s--t
      | /|\ | /|\ |
      |/ | \|/ | \|
      u--v--w--x--y
      
      */

      // ------------------------------------------------------------------------
      // top left fan
      fan0   = spu_shuffle(indices0, indices1, (vu8)shuf_b000aABC);   // g000fabc
      fan0   = spu_shuffle(fan0,     indices2, (vu8)shuf_AcbaEFGH);   // gmlkfabc
      fan1   = spu_shuffle(indices1, indices2, (vu8)shuf_CcXXXXXX);   // hm111111
                                                                            
      center = spu_shuffle(indices1, indices1, (vu8)shuf_BBBBBBBB);   // gggggggg
                                                                            
      cmp0   = spu_cmpeq(fan0, (vu16)mask_FFFF);                      // ????????
      cmp1   = spu_cmpeq(fan1, (vu16)mask_F000);                      // ??000000
                                                                            
      fan0   = spu_sel(fan0, center, cmp0);                           // ????????
      fan1   = spu_sel(fan1, center, cmp1);                           // ??111111

      p_dest[0] = fan0;
      p_dest[1] = fan1;
      
      dest_inc = ~spu_extract(cmp0, 0) & 2;
      p_dest += dest_inc;
      
      // ------------------------------------------------------------------------
      // top right fan
      fan0   = spu_shuffle(indices0, indices1, (vu8)shuf_d0cCDEe0);   // i0hcdej0
      fan0   = spu_shuffle(fan0,     indices2, (vu8)shuf_AcCDEFGe);   // imhcdejo
      fan1   = spu_shuffle(indices1, indices2, (vu8)shuf_dcXXXXXX);   // nm111111

      center = spu_shuffle(indices1, indices1, (vu8)shuf_DDDDDDDD);   // iiiiiiii

      cmp0   = spu_cmpeq(fan0, (vu16)mask_FFFF);                      // ????????
      cmp1   = spu_cmpeq(fan1, (vu16)mask_F000);                      // ??000000

      fan0   = spu_sel(fan0, center, cmp0);                           // ????????
      fan1   = spu_sel(fan1, center, cmp1);                           // ??111111

      p_dest[0] = fan0;
      p_dest[1] = fan1;
      
      dest_inc = ~spu_extract(cmp0, 0) & 2;
      p_dest += dest_inc;

      // ------------------------------------------------------------------------
      // bottom right fan
      fan0   = spu_shuffle(indices2, indices3, (vu8)shuf_dCDEe000);   // smnot000
      fan0   = spu_shuffle(fan0,     indices4, (vu8)shuf_ABCDEedc);   // smnotyxw
      fan1   = spu_shuffle(indices2, indices3, (vu8)shuf_cCXXXXXX);   // rm111111
                                                                    
      center = spu_shuffle(indices3, indices3, (vu8)shuf_DDDDDDDD);   // ssssssss
                                                                            
      cmp0   = spu_cmpeq(fan0, (vu16)mask_FFFF);                      // ????????
      cmp1   = spu_cmpeq(fan1, (vu16)mask_F000);                      // ??000000
                                                                            
      fan0   = spu_sel(fan0, center, cmp0);                           // ????????
      fan1   = spu_sel(fan1, center, cmp1);                           // ??111111

      p_dest[0] = fan0;
      p_dest[1] = fan1;
      
      dest_inc = ~spu_extract(cmp0, 0) & 2;
      p_dest += dest_inc;

      // ------------------------------------------------------------------------
      // bottom left fan
      fan0   = spu_shuffle(indices2, indices3, (vu8)shuf_bCc000aA);   // qmr000pk
      fan0   = spu_shuffle(fan0,     indices4, (vu8)shuf_ABCcbaGH);   // qmrwvupk
      fan1   = spu_shuffle(indices2, indices3, (vu8)shuf_BCXXXXXX);   // lm111111
                                                                    
      center = spu_shuffle(indices3, indices3, (vu8)shuf_BBBBBBBB);   // qqqqqqqq
                                                                            
      cmp0   = spu_cmpeq(fan0, (vu16)mask_FFFF);                      // ????????
      cmp1   = spu_cmpeq(fan1, (vu16)mask_F000);                      // ??000000
                                                                            
      fan0   = spu_sel(fan0, center, cmp0);                           // ????????
      fan1   = spu_sel(fan1, center, cmp1);                           // ??111111

      p_dest[0] = fan0;
      p_dest[1] = fan1;
      
      dest_inc = ~spu_extract(cmp0, 0) & 2;
      p_dest += dest_inc;

      // ------------------------------------------------------------------------
    }
  }

  return (u16 *)p_dest - fans;
}
#endif


//--------------------------------------------------------------------------------------------------------------------------

#if 0
#if 0

u32 R2O_AssignIndices(u16 indices[], u8 outcodes[], u32 cnt)
{
  u32 idx=0;

  for (u32 i=0; i<cnt; i++)
  {
    u32 render_bit = ((outcodes[i] & 0xC0) == 0x40);
    indices[i] = render_bit ? idx : 0xFFFF;
    idx += render_bit;
  }

  return idx;
}

#else

u32 R2O_AssignIndices(u16 indices[], u8 outcodes[], u32 cnt)
{
  IG_ASSERT((cnt           & 0xF)==0);
  IG_ASSERT(((u32)indices  & 0xF)==0);
  IG_ASSERT(((u32)outcodes & 0xF)==0);

  vu16 *p_idx8 = (vu16 *)indices;
  vu8  *p_o16  = (vu8  *)outcodes;
  vu16 idx_assign = spu_splats((u16)0xFFFF);

  vu8  o16, cmp, ones, shifted, acc;
  vu16 acc_lo, acc_hi, cmp_lo, cmp_hi, idx_out;

  for (u32 i=0; i<cnt; i+=16)
  {
    o16 = p_o16[0];

    o16  = spu_and(o16, 0xC0);
    cmp  = spu_cmpeq(o16, 0x40);                                  // eg cmp     = (0,F,F,0,0,0,F,0,0,F,F,F,0,F,0,0)
    ones = spu_and(cmp, 0x01);                                    // eg ones    = (0,1,1,0,0,0,1,0,0,1,1,1,0,1,0,0)

    shifted = spu_rlmaskqwbyte(ones, -1);                         // eg shifted = (0,0,1,1,0,0,0,1,0,0,1,1,1,0,1,0)
    acc = (vu8)spu_add((vu16)ones, (vu16)shifted);                // eg acc     = (0,1,2,1,0,0,1,1,0,1,2,2,1,1,1,0)

    shifted = spu_rlmaskqwbyte(acc, -2);                          // eg shifted = (0,0,0,1,2,1,0,0,1,1,0,1,2,2,1,1)
    acc = (vu8)spu_add((vu16)acc, (vu16)shifted);                 // eg acc     = (0,1,2,2,2,1,1,1,1,2,2,3,3,3,2,1)

    shifted = spu_rlmaskqwbyte(acc, -4);                          // eg shifted = (0,0,0,0,0,1,2,2,2,1,1,1,1,2,2,3)
    acc = (vu8)spu_add((vu16)acc, (vu16)shifted);                 // eg acc     = (0,1,2,2,2,2,3,3,3,3,3,4,4,5,4,4)

    shifted = spu_rlmaskqwbyte(acc, -8);                          // eg shifted = (0,0,0,0,0,0,0,0,0,1,2,2,2,2,3,3)
    acc = (vu8)spu_add((vu16)acc, (vu16)shifted);                 // eg acc     = (0,1,2,2,2,2,3,3,3,4,5,6,6,7,7,7)



    acc_lo = (vu16)spu_shuffle(acc, acc, (vu8)shuf_0A0B0C0D0E0F0G0H);  // eg acc_lo  = (0,1,2,2,2,2,3,3)
    acc_hi = (vu16)spu_shuffle(acc, acc, (vu8)shuf_0I0J0K0L0M0N0O0P);  // eg acc_hi  = (3,4,5,6,6,7,7,7)

    acc_lo  = spu_add(acc_lo, idx_assign);                        // eg acc_lo  = (4,5,6,6,6,6,7,7)
    acc_hi  = spu_add(acc_hi, idx_assign);                        // eg acc_hi  = (7,8,9,A,A,B,B,B)

    cmp_lo = (vu16)spu_shuffle(cmp, cmp, (vu8)shuf_AABBCCDDEEFFGGHH);  // eg cmp_lo  = (0,F,F,0,0,0,F,0
    cmp_hi = (vu16)spu_shuffle(cmp, cmp, (vu8)shuf_IIJJKKLLMMNNOOPP);  // eg cmp_hi  = (0,F,F,F,0,F,0,0)

    idx_out = spu_sel((vu16)mask_FFFF, acc_lo, cmp_lo);           // eg idx_out = (X,5,6,X,X,X,7,X)
    p_idx8[0] = idx_out;

    idx_out = spu_sel((vu16)mask_FFFF, acc_hi, cmp_hi);           // eg idx_out = (X,8,9,A,X,B,X,X)
    p_idx8[1] = idx_out;



    idx_assign = spu_shuffle(acc_hi, acc_hi, (vu8)shuf_HHHHHHHH);



    p_o16  += 1;
    p_idx8 += 2;
  }

  return spu_extract(idx_assign, 0) + 1;
}

#endif
#endif


//--------------------------------------------------------------------------------------------------------------------------

#if 0
#if 0
void R2O_FlagVertsRender(u8 outcodes[], i32 nc, i32 nr)
{
  for (i32 r=0; r<nr; r++)
  {
    for (i32 c=0; c<nc; c++)
    {
      u32 i = r*nc + c;
      u8 o_or = outcodes[i] | outcodes[i+1] | outcodes[i+nc] | outcodes[i+nc+1];

      if (o_or & 0x40)
      {
        outcodes[i] |= 0x40;
      }
      else
      {
        outcodes[i] &= ~0x40;
      }
    }
  }
}
#else
void R2O_FlagVertsRender(u8 outcodes[], i32 nc, i32 nr)
{
  vu8 *p_o = (vu8 *)outcodes;
  vu8 oTL, oTR, oBL, oBR, oT1, oB0, oB1, orL, or1, o_or, o_out, shuf0, shuf1;

  u32 cnt    = nc*nr;
  u32 stride = nc>>4;

  shuf0 = (nc&8) ? (vu8)shuf_IJKLMNOPabcdefgh : (vu8)shuf_ABCDEFGHIJKLMNOP;
  shuf1 = (nc&8) ? (vu8)shuf_JKLMNOPabcdefghi : (vu8)shuf_BCDEFGHIJKLMNOPa;

  oTR = p_o[0];
  oBR = p_o[stride];

  for (u32 i=0; i<cnt; i+=16)
  {
    oTL = oTR;
    oBL = oBR;
    oTR = p_o[1];
    oBR = p_o[stride+1];

    oT1 = spu_shuffle(oTL, oTR, (vu8)shuf_BCDEFGHIJKLMNOPa);

    oB0 = spu_shuffle(oBL, oBR, shuf0);
    oB1 = spu_shuffle(oBL, oBR, shuf1);

    orL = spu_or(oTL, oB0);
    or1 = spu_or(oT1, oB1);

    o_or = spu_or(orL, or1);

    o_out = spu_sel(oTL, o_or, vu8_0x40);

    p_o[0] = o_out;
    p_o++;
  }
}
#endif
#endif


//--------------------------------------------------------------------------------------------------------------------------


#if 0
#if 0
void R2O_FlagVertsKeep(u8 outcodes[], i32 nc, i32 nr)
{
  for (i32 r=0; r<nr; r++)
  {
    for (i32 c=0; c<nc; c++)
    {
      u32 i = r*nc + c;
      u8 o_or = outcodes[i] | outcodes[i+1] | outcodes[i+nc] | outcodes[i+nc+1];

      if ((o_or & 0x20) == 0)
      {
        outcodes[i] |= 0x80;
      }
    }
  }
}
#else
void R2O_FlagVertsKeep(u8 outcodes[], i32 nc, i32 nr)
{
  vu8 *p_o = (vu8 *)outcodes;
  vu8 oTL, oTR, oBL, oBR, oT1, oB0, oB1, orL, or1, o_or, or_masked, cmp, cull, o_out, shuf0, shuf1;

  u32 cnt    = nc*nr;
  u32 stride = nc>>4;

  shuf0 = (nc&8) ? (vu8)shuf_IJKLMNOPabcdefgh : (vu8)shuf_ABCDEFGHIJKLMNOP;
  shuf1 = (nc&8) ? (vu8)shuf_JKLMNOPabcdefghi : (vu8)shuf_BCDEFGHIJKLMNOPa;

  oTR = p_o[0];
  oBR = p_o[stride];

  for (u32 i=0; i<cnt; i+=16)
  {
    oTL = oTR;
    oTR = p_o[1];
    oBL = oBR;
    oBR = p_o[stride+1];

    oT1 = spu_shuffle(oTL, oTR, (vu8)shuf_BCDEFGHIJKLMNOPa);

    oB0 = spu_shuffle(oBL, oBR, shuf0);
    oB1 = spu_shuffle(oBL, oBR, shuf1);

    orL = spu_or(oTL, oB0);
    or1 = spu_or(oT1, oB1);

    o_or = spu_or(orL, or1);

    or_masked = spu_and(o_or, 0x20);
    cmp = spu_cmpeq(or_masked, 0x00);

    cull = spu_and(cmp, 0x80);
    o_out = spu_or(oTL, cull);

    p_o[0] = o_out;
    p_o++;
  }
}
#endif
#endif


//--------------------------------------------------------------------------------------------------------------------------


#if 0
#if 0
void R2O_FlagQuadsRenderKeep(u8 outcodes[], i32 nc, i32 nr)
{
  u32 cnt = nc*nr;
  
  for (u32 i=0; i<cnt; i++)
  {
    u8 o_and = outcodes[i] & outcodes[i-1] & outcodes[i-nc] & outcodes[i-nc-1];

    // mark the quad 'render' if it's all beyond the near plane, isn't all outside any of planes TBLR,
    // and isn't incident to a culled vertex
    if ((o_and & 0x9F) == 0x01)
    {
      outcodes[i] |= 0x40;
    }

    // mark the quad 'keep' if it isn't all beyond the near plane, isn't all outside any of planes TBLR,
    // and isn't incident to a culled vertex
    if ((o_and & 0x9F) == 0x00)
    {
      outcodes[i] |= 0x20;
    }
  }
}
#else
void R2O_FlagQuadsRenderKeep(u8 outcodes[], i32 nc, i32 nr)
{
  vu8 *p_o = (vu8 *)outcodes;
  vu8 oTL, oTR, oBL, oBR, oT0, oT1, oB0, andR, and0, o_and, and_masked, cmp, bit, o_out, shuf0, shuf1;

  u32 cnt = nc*nr;
  u32 stride = nc>>4;

  shuf0 = (nc&8) ? (vu8)shuf_HIJKLMNOPabcdefg : (vu8)shuf_Pabcdefghijklmno;
  shuf1 = (nc&8) ? (vu8)shuf_IJKLMNOPabcdefgh : (vu8)shuf_abcdefghijklmnop;

  oTR = p_o[-stride-1];
  oBR = p_o[-1];

  for (u32 i=0; i<cnt; i+=16)
  {
    oTL = oTR;
    oBL = oBR;
    oTR = p_o[-stride];
    oBR = p_o[0];

    oT0 = spu_shuffle(oTL, oTR, shuf0);
    oT1 = spu_shuffle(oTL, oTR, shuf1);

    oB0 = spu_shuffle(oBL, oBR, (vu8)shuf_Pabcdefghijklmno);

    and0 = spu_and(oT0, oB0);
    andR = spu_and(oT1, oBR);

    o_and = spu_and(and0, andR);

    and_masked = spu_and(o_and, 0x9F);

    cmp = spu_cmpeq(and_masked, 0x01);
    bit = spu_and(cmp, 0x40);
    o_out = spu_or(oBR, bit);

    cmp = spu_cmpeq(and_masked, 0x00);
    bit = spu_and(cmp, 0x20);
    o_out = spu_or(o_out, bit);

    p_o[0] = o_out;
    p_o++;
  }
}
#endif
#endif



//--------------------------------------------------------------------------------------------------------------------------

#if 0
void R2O_GenerateVertexDerivs(i16 *dest, f32 *src, f32 step, vf32 basis_col, vf32 basis_row, i32 nc, i32 nr)
{
  // check alignments
  IG_ASSERT(((u32)src   & 0xF) == 0);
  IG_ASSERT(((u32)dest  & 0xF) == 0);
  IG_ASSERT((nc         & 0x3) == 0);

  f32 vert_deriv_scale = (0.2f / g_WaterObject.m_amplitude);
  vf32 scale = spu_splats((0.5f / step) * vert_deriv_scale);    // why is there a 0.5 in there?
  basis_col *= scale;
  basis_row *= scale;

  vf32 basis_col_x = spu_shuffle(basis_col, basis_col, (vu8)shuf_AAAA);
  vf32 basis_col_z = spu_shuffle(basis_col, basis_col, (vu8)shuf_CCCC);
  vf32 basis_row_x = spu_shuffle(basis_row, basis_row, (vu8)shuf_AAAA);
  vf32 basis_row_z = spu_shuffle(basis_row, basis_row, (vu8)shuf_CCCC);
  
  vf32 a03, bCF, bF2, b03=spu_splats(0.0f), b14, b47, c03;
  vf32 diff_x, diff_z, dy_dx, dy_dz;
  vi32 dy_dx_int, dy_dz_int;
  vi16 out;

  u32 qstride = nc >> 2;

  for (i32 c=0; c<nc; c+=4)
  {
    vi16 *p_dest    = (vi16 *)&dest[c<<1];
    vf32 *p_src     = (vf32 *)&src[c];

    // initialize queue
    c03 = p_src[0];

    // loop over destination rows
    for (i32 r=0; r<nr; r++)
    {
      // advance queue
      a03 = b03;
      b03 = c03;

      // load 3 qwords
      bCF = p_src[-1];
      b47 = p_src[1];
      c03 = p_src[qstride];

      // shifts
      bF2 = spu_shuffle(bCF, b03, (vu8)shuf_Dabc);
      b14 = spu_shuffle(b03, b47, (vu8)shuf_BCDa);

      // form differences
      diff_x = b14-bF2;
      diff_z = c03-a03;

      // compute derivs
      dy_dx = (diff_x * basis_col_x) + (diff_z * basis_row_x);
      dy_dz = (diff_x * basis_col_z) + (diff_z * basis_row_z);

      // convert to i16
      dy_dx_int = spu_convts(dy_dx, 0);
      dy_dz_int = spu_convts(dy_dz, 0);

      // interleave
      out = spu_shuffle(*(vi16 *)&dy_dx_int, *(vi16 *)&dy_dz_int, (vu8)shuf_AaCcEeGg);

      // store
      *p_dest = out;

      // step pointers
      p_src  += qstride;
      p_dest += qstride;
    }
  }
}
#endif


//--------------------------------------------------------------------------------------------------------------------------
/*

              a4 a5 a6 a7 a8 a9 aA aB
  b0 b1 b2 b3 b4 b5 b6 b7 b8 b9 bA bB bC bD bE bF
              c4 c5 c6 c7 c8 c9 cA cB
              
*/

#if 0
void R2O_GenerateMapDerivs(i8 derivs[], f32 heights[], i32 cols, i32 rows, u32 src_stride, u32 dest_stride,
                           u32 wrap, f32 inv_step, vf32 basis_col, vf32 basis_row)
{
  IG_ASSERT((cols & 7) == 0);

  f32 map_deriv_scale  = 128.0f;
  vf32 scale = spu_splats(inv_step * map_deriv_scale);
  basis_col *= scale;
  basis_row *= scale;

  vf32 basis_col_x = spu_shuffle(basis_col, basis_col, (vu8)shuf_AAAA);
  vf32 basis_col_z = spu_shuffle(basis_col, basis_col, (vu8)shuf_CCCC);
  vf32 basis_row_x = spu_shuffle(basis_row, basis_row, (vu8)shuf_AAAA);
  vf32 basis_row_z = spu_shuffle(basis_row, basis_row, (vu8)shuf_CCCC);
  
  vf32 a47, a8B, b03, b36, b47, b58, b7A, b8B, b9C, bCF, c47, c8B;
  vf32 diff_x_L, diff_x_R, diff_z_L, diff_z_R;
  vf32 dy_dx_L, dy_dx_R, dy_dz_L, dy_dz_R;
  vi32 dy_dx_L_int, dy_dx_R_int, dy_dz_L_int, dy_dz_R_int;
  vi16 dy_dx, dy_dz;
  vi8  results;
  i32  idx_src;

  vf32 *p_srcN, *p_src0, *p_src1, *p_src2;
  vi8  *p_dest;
  
  i32 src_qstride  = src_stride >>2;
  i32 dest_qstride = dest_stride>>4;

  // loop over source cols
  for (i32 c=0; c<cols; c+=8)
  {
    // set pointers
    p_dest = (vi8  *)&derivs[c<<1];
    p_src0 = (vf32 *)&heights[c];
    p_src1 = p_src0 + 1;
    p_srcN = p_src0 + ((wrap && (cols-c == cols)) ? src_qstride :  0) - 1;
    p_src2 = p_src0 - ((wrap && (cols-c == 8   )) ? src_qstride :  0) + 2;
    idx_src = 0;

    // initialize queue
    b47 = p_src0[idx_src - src_qstride]; // don't need to wrap here because the calling code
    b8B = p_src1[idx_src - src_qstride]; //     has duplicated top & bottom rows
    c47 = p_src0[idx_src];
    c8B = p_src1[idx_src];

    // loop over source rows
    for (i32 r=0; r<rows; r++)
    {
      // advance queue
      a47 = b47;
      a8B = b8B;
      b47 = c47;
      b8B = c8B;

      // load 2 qwords from current line
      b03 = p_srcN[idx_src];
      bCF = p_src2[idx_src];

      // step source index
      // (we don't need to wrap here, because the calling code adds top & bottom row duplicates at the bottom & top resp.)
      idx_src = idx_src + src_qstride;

      // load 2 qwords from next line
      c47 = p_src0[idx_src];
      c8B = p_src1[idx_src];

      // shuffles
      b36 = (vf32)si_shufb((qword)b03, (qword)b47, shuf_Dabc);
      b7A = (vf32)si_shufb((qword)b47, (qword)b8B, shuf_Dabc);
      b58 = (vf32)si_shufb((qword)b47, (qword)b8B, shuf_BCDa);
      b9C = (vf32)si_shufb((qword)b8B, (qword)bCF, shuf_BCDa);

      // diffs
      diff_x_L = b58 - b36;
      diff_x_R = b9C - b7A;
      diff_z_L = c47 - a47;
      diff_z_R = c8B - a8B;

      // transform to world coords
      dy_dx_L = (diff_x_L * basis_col_x) + (diff_z_L * basis_row_x);
      dy_dz_L = (diff_x_L * basis_col_z) + (diff_z_L * basis_row_z);
      dy_dx_R = (diff_x_R * basis_col_x) + (diff_z_R * basis_row_x);
      dy_dz_R = (diff_x_R * basis_col_z) + (diff_z_R * basis_row_z);

      // x-derivs
      dy_dx_L_int = (vi32)si_cflts((qword)dy_dx_L, 24);
      dy_dx_R_int = (vi32)si_cflts((qword)dy_dx_R, 24);

      // z-derivs
      dy_dz_L_int = (vi32)si_cflts((qword)dy_dz_L, 24);
      dy_dz_R_int = (vi32)si_cflts((qword)dy_dz_R, 24);

      // interleave
      dy_dx   = (vi16)si_shufb((qword)dy_dx_L_int, (qword)dy_dx_R_int, shuf_ACEGaceg);
      dy_dz   = (vi16)si_shufb((qword)dy_dz_L_int, (qword)dy_dz_R_int, shuf_ACEGaceg);
      results = (vi8 )si_shufb((qword)dy_dx,       (qword)dy_dz,       shuf_AaCcEeGgIiKkMmOo);

      // store results
      *p_dest = results;

      // step dest pointer
      p_dest += dest_qstride;
    }
  }
}
#endif

// ------------------------------------------------------------------------------------------------------------------------

#if 0
void R2O_InterpTileBilinear(f32 *dest, f32 *src)
{
  vf32 a03,a14,a47,k03,h03,v03;
  vf32 *p_src0, *p_src1, *p_dest;
  u32  idx0, idx1, addr_diff;

  p_src0  = (vf32 *)&src[28];
  p_src1  = (vf32 *)&src[0];
  addr_diff = (u32)dest-(u32)src;

  for (i32 c=0; c<32; c+=4)
  {
    // set dest from p_src0
    p_dest = (vf32 *)((u32)p_src0 + addr_diff);
  
    // initialize queue using bottom row
    idx0 = 31*8;
    idx1 = 0;
    a03 = p_src0[idx0];
    a47 = p_src1[idx0];
    a14 = spu_shuffle(a03, a47, (vu8)shuf_BCDa);
    h03 = half * (a03 + a14);

    // loop over rows
    for (i32 r=0; r<32; r++)
    {
      // get next 2 qwords
      a03 = p_src0[idx1];
      a47 = p_src1[idx1];

      // shift left by 1
      a14 = spu_shuffle(a03, a47, (vu8)shuf_BCDa);

      // interpolate horizontally
      k03 = h03;
      h03 = half * (a03 + a14);

      // interpolate vertically
      v03 = half * (h03 + k03);

      // store result
      p_dest[idx0] = v03;

      // step index
      idx0 = idx1;
      idx1 += 8;
    }
    
    // set pointers for next col
    p_src0 = p_src1;
    p_src1++;
  }
}
#endif


//------------------------------------------------------------------------------------------------------------------------

/*


              /\
             /w0\ dest
            /    \
           /Z0  W0\
          /       /\
         /y1  z1 /w1\
       +/------ /----\-----------------+
       /X1  Y1 /Z1  W1\ ..  ..  ..  .. | source
       \      /       /\               |
       |\ x2 /y2  z2 /w2\ ..  ..  ..  ..
       | \  /       /    \             |
       |..\/X2  Y2 /Z2  W2\ ..  ..  .. |
       |   \      /       /\           |
       |  ..\ x3 /y3  z3 /w3\ ..  ..  ..
       |     \  /       /    \         |
       |..  ..\/X3  Y3 /Z3  W3\ ..  .. |
       |       \      /       /        |
       |  ..  ..\ x4 /y4  z4 /w4  ..  ..
       |         \  /       /          |
       |..  ..  ..\/X4  Y4 /Z4  W4  .. |
       |           \      /            |
       |  ..  ..  ..\ x5 /y5  z5  w5  ..
       +-------------\--/--------------+
                      \/



*/



#if 0
void R2O_RotateTileClockwise(f32 *dest, f32 *src_even, f32 *src_odd)
{
  vf32 L2, R2, L3, R3, L4, R4;
  vf32 Q0, Q1, Q2, Q3, Q4;
  vf32 l3, r3, l4, r4, l5, r5;
  vf32 q0, q1, q2, q3, q4, q5;
  vf32 Q10, Q21, Q32, Q43;
  vf32 q10, q21, q32, q43, q54;
  vf32 out0, out1, out2, out3, out4, out5, out6, out7;
  u32 idx_src0, idx_src1, idx_dest;

  // set source & dest pointers
  vf32 *p_src_even_0x00 = (vf32 *)src_even + 0x00;
  vf32 *p_src_even_0x08 = (vf32 *)src_even + 0x08;
  vf32 *p_src_even_0x10 = (vf32 *)src_even + 0x10;
  vf32 *p_src_even_0x18 = (vf32 *)src_even + 0x18;

  vf32 *p_src_odd__0x00 = (vf32 *)src_odd  + 0x00;
  vf32 *p_src_odd__0x08 = (vf32 *)src_odd  + 0x08;
  vf32 *p_src_odd__0x10 = (vf32 *)src_odd  + 0x10;
  vf32 *p_src_odd__0x18 = (vf32 *)src_odd  + 0x18;

  vf32 *p_dest_0x00 = (vf32 *)dest + 0x00;
  vf32 *p_dest_0x08 = (vf32 *)dest + 0x08;
  vf32 *p_dest_0x10 = (vf32 *)dest + 0x10;
  vf32 *p_dest_0x18 = (vf32 *)dest + 0x18;

  // loop start pos over source cols
  for (i32 c=0; c<32; c+=4)
  {
    // initialize source indices
    idx_src0 = c>>2;
    idx_src1 = (idx_src0 + 1) & 7;

    // initialize dest index
    idx_dest = (c>>2) + (c<<3);

    // initialize queue
    q4 = p_src_odd__0x10[idx_src0 + 0xE0];
    Q4 = p_src_even_0x18[idx_src0 + 0xE0];
    q5 = p_src_odd__0x18[idx_src0 + 0xE0];
    q4 = spu_shuffle(q4, q4, (vu8)shuf_CDab);
    Q4 = spu_shuffle(Q4, Q4, (vu8)shuf_Dabc);
    q5 = spu_shuffle(q5, q5, (vu8)shuf_Dabc);

    // loop over dest rows
    for (i32 r=0; r<64; r+=8)
    {
      // advance queue
      q0 = q4;
      Q0 = Q4;
      q1 = q5;

      // load qwords
      Q1 = p_src_even_0x00[idx_src0];
      q2 = p_src_odd__0x00[idx_src0];
      L2 = p_src_even_0x08[idx_src0];
      R2 = p_src_even_0x08[idx_src1];
      l3 = p_src_odd__0x08[idx_src0];
      r3 = p_src_odd__0x08[idx_src1];
      L3 = p_src_even_0x10[idx_src0];
      R3 = p_src_even_0x10[idx_src1];
      l4 = p_src_odd__0x10[idx_src0];
      r4 = p_src_odd__0x10[idx_src1];
      L4 = p_src_even_0x18[idx_src0];
      R4 = p_src_even_0x18[idx_src1];
      l5 = p_src_odd__0x18[idx_src0];
      r5 = p_src_odd__0x18[idx_src1];

      // shuffle
      Q2 = spu_shuffle(L2, R2, (vu8)shuf_BCDa);
      q3 = spu_shuffle(l3, r3, (vu8)shuf_BCDa);
      Q3 = spu_shuffle(L3, R3, (vu8)shuf_CDab);
      q4 = spu_shuffle(l4, r4, (vu8)shuf_CDab);
      Q4 = spu_shuffle(L4, R4, (vu8)shuf_Dabc);
      q5 = spu_shuffle(l5, r5, (vu8)shuf_Dabc);

      // select
      Q10 = spu_sel(Q1, Q0, mask_00FF);
      Q21 = spu_sel(Q2, Q1, mask_00FF);
      Q32 = spu_sel(Q3, Q2, mask_00FF);
      Q43 = spu_sel(Q4, Q3, mask_00FF);

      q10 = spu_sel(q1, q0, mask_00FF);
      q21 = spu_sel(q2, q1, mask_00FF);
      q32 = spu_sel(q3, q2, mask_00FF);
      q43 = spu_sel(q4, q3, mask_00FF);
      q54 = spu_sel(q5, q4, mask_00FF);

      out0 = spu_sel(Q10, q10, mask_0F0F);
      out1 = spu_sel(q21, Q10, mask_0F0F);
      out2 = spu_sel(Q21, q21, mask_0F0F);
      out3 = spu_sel(q32, Q21, mask_0F0F);
      out4 = spu_sel(Q32, q32, mask_0F0F);
      out5 = spu_sel(q43, Q32, mask_0F0F);
      out6 = spu_sel(Q43, q43, mask_0F0F);
      out7 = spu_sel(q54, Q43, mask_0F0F);

      // output first 4x4 results
      p_dest_0x00[idx_dest] = out0;
      p_dest_0x08[idx_dest] = out1;
      p_dest_0x10[idx_dest] = out2;
      p_dest_0x18[idx_dest] = out3;

      // step 4 rows down in dest
      idx_dest = (idx_dest + 0x20) & 0x1FF;

      // output second 4x4 results
      p_dest_0x00[idx_dest] = out4;
      p_dest_0x08[idx_dest] = out5;
      p_dest_0x10[idx_dest] = out6;
      p_dest_0x18[idx_dest] = out7;

      // step 4 rows down in dest
      idx_dest = (idx_dest + 0x20) & 0x1FF;

      // step 4 rows down and 4 cols (1 qword) right in source
      idx_src0 = (idx_src1 + 0x20) & 0xFF;
      idx_src1 = (idx_src0 & 0xF8)  | ((idx_src0 + 1) & 7);
    }
  }
}
#endif

//------------------------------------------------------------------------------------------------------------------------


/*




                              /\
                             /X0\
                            /    \ dest
                           /x0  y0\
                          /\       \
                         /X1\ Y1  Z1\
                        /    \       \
                       /x1  y1\ z1  w1\
     +----------------/\-------\-----+/
     |..  ..  ..  .. /X2\ Y2  Z2\ W2 /(O)origin
     |              /    \       \  /| 
     |  ..  ..  .. /x2  y2\ z2  w2\/..
     |            /\       \      /  | 
     |..  ..  .. /X3\ Y3  Z3\ W3 /.. | 
     |          /    \       \  /    | 
     |  ..  .. /x3  y3\ z3  w3\/..  ..
     |         \       \      /      |
     |..  ..  X4\ Y4  Z4\ W4 /..  .. |
     |           \       \  /        |
     |  ..  x4  y4\ z4  w4\/..  ..  ..
     |             \      /          |
     |..  X5  Y5  Z5\ W5 /..  ..  .. |
     |               \  /            | source
     |  x5  y5  z5  w5\/..  ..  ..  ..
     +-------------------------------+

*/


#if 0
void R2O_RotateTileAntiClockwise(f32 *dest, f32 *src_even, f32 *src_odd)
{
  vf32 L3, R3, L4, R4, L5, R5;
  vf32 Q0, Q1, Q2, Q3, Q4, Q5;
  vf32 l2, r2, l3, r3, l4, r4;
  vf32 q0, q1, q2, q3, q4, q5;
  vf32 Q01, Q12, Q23, Q34, Q45;
  vf32 q01, q12, q23, q34;
  vf32 out0, out1, out2, out3, out4, out5, out6, out7;
  u32 idx_src0, idx_src1, idx_dest;

  // set source & dest pointers
  vf32 *p_src_even_0x00 = (vf32 *)src_even + 0x00;
  vf32 *p_src_even_0x08 = (vf32 *)src_even + 0x08;
  vf32 *p_src_even_0x10 = (vf32 *)src_even + 0x10;
  vf32 *p_src_even_0x18 = (vf32 *)src_even + 0x18;

  vf32 *p_src_odd__0x00 = (vf32 *)src_odd  + 0x00;
  vf32 *p_src_odd__0x08 = (vf32 *)src_odd  + 0x08;
  vf32 *p_src_odd__0x10 = (vf32 *)src_odd  + 0x10;
  vf32 *p_src_odd__0x18 = (vf32 *)src_odd  + 0x18;

  vf32 *p_dest_0x00 = (vf32 *)dest + 0x00;
  vf32 *p_dest_0x08 = (vf32 *)dest + 0x08;
  vf32 *p_dest_0x10 = (vf32 *)dest + 0x10;
  vf32 *p_dest_0x18 = (vf32 *)dest + 0x18;

  // loop start pos over source cols, right to left
  for (i32 c=28; c>=0; c-=4)
  {
    // initialize source indices
    idx_src1 = c>>2;
    idx_src0 = (idx_src1 - 1) & 7;

    // initialize dest index
    idx_dest = (c>>2) + ((60-c)<<3);

    // initialize queue
    Q4 = p_src_even_0x10[idx_src1 + 0xE0];
    q4 = p_src_odd__0x10[idx_src1 + 0xE0];
    Q5 = p_src_even_0x18[idx_src1 + 0xE0];
    q5 = p_src_odd__0x18[idx_src1 + 0xE0];
    Q4 = spu_shuffle(Q4, Q4, (vu8)shuf_CDab);
    q4 = spu_shuffle(q4, q4, (vu8)shuf_BCDa);
    Q5 = spu_shuffle(Q5, Q5, (vu8)shuf_BCDa);

    // loop over dest rows
    for (i32 r=0; r<64; r+=8)
    {
      // advance queue
      Q0 = Q4;
      q0 = q4;
      Q1 = Q5;
      q1 = q5;

      // load qwords
      Q2 = p_src_even_0x00[idx_src1];
      l2 = p_src_odd__0x00[idx_src0];
      r2 = p_src_odd__0x00[idx_src1];
      L3 = p_src_even_0x08[idx_src0];
      R3 = p_src_even_0x08[idx_src1];
      l3 = p_src_odd__0x08[idx_src0];
      r3 = p_src_odd__0x08[idx_src1];
      L4 = p_src_even_0x10[idx_src0];
      R4 = p_src_even_0x10[idx_src1];
      l4 = p_src_odd__0x10[idx_src0];
      r4 = p_src_odd__0x10[idx_src1];
      L5 = p_src_even_0x18[idx_src0];
      R5 = p_src_even_0x18[idx_src1];
      q5 = p_src_odd__0x18[idx_src0];

      // shuffle
      q2 = spu_shuffle(l2, r2, (vu8)shuf_Dabc);
      Q3 = spu_shuffle(L3, R3, (vu8)shuf_Dabc);
      q3 = spu_shuffle(l3, r3, (vu8)shuf_CDab);
      Q4 = spu_shuffle(L4, R4, (vu8)shuf_CDab);
      q4 = spu_shuffle(l4, r4, (vu8)shuf_BCDa);
      Q5 = spu_shuffle(L5, R5, (vu8)shuf_BCDa);

      // select
      Q01 = spu_sel(Q0, Q1, mask_00FF);
      Q12 = spu_sel(Q1, Q2, mask_00FF);
      Q23 = spu_sel(Q2, Q3, mask_00FF);
      Q34 = spu_sel(Q3, Q4, mask_00FF);
      Q45 = spu_sel(Q4, Q5, mask_00FF);

      q01 = spu_sel(q0, q1, mask_00FF);
      q12 = spu_sel(q1, q2, mask_00FF);
      q23 = spu_sel(q2, q3, mask_00FF);
      q34 = spu_sel(q3, q4, mask_00FF);

      out0 = spu_sel(Q01, q01, mask_0F0F);
      out1 = spu_sel(q01, Q12, mask_0F0F);
      out2 = spu_sel(Q12, q12, mask_0F0F);
      out3 = spu_sel(q12, Q23, mask_0F0F);
      out4 = spu_sel(Q23, q23, mask_0F0F);
      out5 = spu_sel(q23, Q34, mask_0F0F);
      out6 = spu_sel(Q34, q34, mask_0F0F);
      out7 = spu_sel(q34, Q45, mask_0F0F);

      // output first 4x4 results
      p_dest_0x00[idx_dest] = out0;
      p_dest_0x08[idx_dest] = out1;
      p_dest_0x10[idx_dest] = out2;
      p_dest_0x18[idx_dest] = out3;

      // step 4 rows down in dest
      idx_dest = (idx_dest + 0x20) & 0x1FF;

      p_dest_0x00[idx_dest] = out4;
      p_dest_0x08[idx_dest] = out5;
      p_dest_0x10[idx_dest] = out6;
      p_dest_0x18[idx_dest] = out7;

      // step 4 rows down in dest
      idx_dest = (idx_dest + 0x20) & 0x1FF;

      // step 4 rows down and 4 cols (1 qword) left in source
      idx_src1 = (idx_src0 + 0x20) & 0xFF;
      idx_src0 = (idx_src1 & 0xF8)  | ((idx_src1 - 1) & 7);

    }
  }
}
#endif


//------------------------------------------------------------------------------------------------------------------------

#if 0
void R2O_RotateTile(f32 *dest, f32 *src_even, f32 *src_odd, u32 b_anticlockwise)
{
  if (b_anticlockwise)
  {
    R2O_RotateTileAntiClockwise(dest, src_even, src_odd);
  }
  else
  {
    R2O_RotateTileClockwise(dest, src_even, src_odd);
  }
}
#endif

//------------------------------------------------------------------------------------------------------------------------



// Make 1 pass over the low lod map, which is either a 32x32 fft tile or a 64x64 intermediate-step map.
// The mid lod map is a 32x32 fft tile, rotated 45 degrees in the appropriate direction and bilinearly interpolated
// yielding a 32x64 map.
// The high lod map is a straight 32x32 fft tile.
// Process row-pairs from bottom to top, so we can initialize the queue with the wrapped values.

#if 0
void R2O_AccumulateNormalMaps(f32 *dest, f32 *lo, f32 scale_lo, f32 *mid, f32 scale_mid, f32 *hi, f32 scale_hi,
                              i32 nc_dest, i32 nr_dest)
{
  vf32 *p_lo, *p_mid, *p_hi, *p_dest;
  vf32 lo0, lo1, lo2, lo3, mid0, mid1, mid2, mid3, hi0, hi1, hi2, hi3, out0, out1, out2, out3;
  vf32 loL, loR, loLa, loRa, loLb, loRb, lo0_next, lo1_next;
  u32 ofs_lo, idx_mid, idx_hi;

  u32 qstride_lo   = nc_dest >> 3;
  u32 qstride_dest = nc_dest >> 2;

  vf32 qscale_lo  = spu_splats(scale_lo);
  vf32 qscale_mid = spu_splats(scale_mid);
  vf32 qscale_hi  = spu_splats(scale_hi);

  // fix mid/hi ptrs at start of array
  p_hi  = (vf32 *)hi;
  p_mid = (vf32 *)mid;

  // loop over destination columns
  for (i32 c=0; c<nc_dest; c+=8)
  {
    // init lo ptr/ofs temporarily on the top row (because we need to init the queue with the wrapped values)
    p_lo = (vf32 *)&lo[(c>>1)];
    ofs_lo = (c+8==nc_dest) ? 1-qstride_lo : 1;

    // init mid idx 2 rows up from bottom or from half-way down
    idx_mid = ((c>>2) & 7) + 30*8 + ((c&32) ? 0 : 32*8);

    // init hi idx on penultimate row
    idx_hi  = ((c>>2) & 7) + 30*8;

    // init dest ptr on the penultimate row
    p_dest = (vf32 *)&dest[c + nc_dest*(nr_dest-2)];

    // initialize input queue for low lod
    loL = p_lo[0];
    loR = p_lo[ofs_lo];

    loLa = spu_shuffle(loL, loR, (vu8)shuf_AABB);
    loLb = spu_shuffle(loL, loR, (vu8)shuf_ABBC);

    loRa = spu_shuffle(loL, loR, (vu8)shuf_CCDD);
    loRb = spu_shuffle(loL, loR, (vu8)shuf_CDDa);

    lo0_next = half * (loLa + loLb);
    lo1_next = half * (loRa + loRb);

    // wrap the lo ptr to the bottom row
    p_lo = (vf32 *)&lo[(c>>1) + (nc_dest>>1)*((nr_dest>>1)-1)];

    // loop over destination row pairs in reverse order
    for (i32 r=nr_dest-2; r>=0; r-=2)
    {
      // low lod
      loL = p_lo[0];                                // (xL,yL,zL,wL)
      loR = p_lo[ofs_lo];                           // (xR,yR,zR,wR)

      loLa = spu_shuffle(loL, loR, (vu8)shuf_AABB); // (xL,xL,yL,yL)
      loLb = spu_shuffle(loL, loR, (vu8)shuf_ABBC); // (xL,yL,yL,zL)

      loRa = spu_shuffle(loL, loR, (vu8)shuf_CCDD); // (zL,zL,wL,wL)
      loRb = spu_shuffle(loL, loR, (vu8)shuf_CDDa); // (zL,wL,wL,xR))

      lo0 = half * (loLa + loLb);                   // (xL, (xL+yL)/2, yL, (yL+zL)/2)
      lo1 = half * (loRa + loRb);                   // (zL, (zL+wL)/2, wL, (wL+xR)/2)

      lo2 = half * (lo0 + lo0_next);                // ((xL+xL')/2, (xL+yL+xL'+yL')/4, (yL+yL')/2, (yL+zL+yL'+zL')/4)
      lo3 = half * (lo1 + lo1_next);                // ((zL+zL')/2, (zL+wL+zL'+wL')/4, (wL+wL')/2, (wL+xR+wL'+xR')/4)
      
      lo0_next = lo0;
      lo1_next = lo1;

      // mid lod
      mid0 = p_mid[idx_mid + 0];
      mid1 = p_mid[idx_mid + 1];
      mid2 = p_mid[idx_mid + 8];
      mid3 = p_mid[idx_mid + 9];

      // high lod
      hi0 = p_hi[idx_hi + 0];
      hi1 = p_hi[idx_hi + 1];
      hi2 = p_hi[idx_hi + 8];
      hi3 = p_hi[idx_hi + 9];

      // accumulate
      out0 = (qscale_lo) * lo0 + (qscale_mid) * mid0 + (qscale_hi) * hi0;
      out1 = (qscale_lo) * lo1 + (qscale_mid) * mid1 + (qscale_hi) * hi1;
      out2 = (qscale_lo) * lo2 + (qscale_mid) * mid2 + (qscale_hi) * hi2;
      out3 = (qscale_lo) * lo3 + (qscale_mid) * mid3 + (qscale_hi) * hi3;

      // store
      p_dest[0]              = out0;
      p_dest[1]              = out1;
      p_dest[qstride_dest]   = out2;
      p_dest[qstride_dest+1] = out3;

      // step pointers
      p_lo    -= qstride_lo;      // up 1 row
      idx_mid -= 16;              // up 2 rows
      idx_mid &= 0x1FF;           // stay within 8K (512 qwords)
      idx_hi  -= 16;              // up 2 rows
      idx_hi  &= 0xFF;            // stay within 4K (256 qwords)
      p_dest  -= qstride_dest<<1; // up 2 rows
    }
  }
}
#endif

