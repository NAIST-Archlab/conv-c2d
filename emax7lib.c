
/* EMAX7 library                        */
/*         Copyright (C) 2013- by NAIST */
/*          Primary writer: Y.Nakashima */
/*                 nakashim@is.naist.jp */
#include "emax7lib.h"

/*******************************************************************************/
/******************************** EMAX7 NCLIB (no conv-c2c)*********************/
/*******************************************************************************/
#ifndef NO_EMAX7LIB_BODY
void /*__attribute__((always_inline))*/
cex(Uint op_cx, Ull *ex, Ull c3, Ull c2, Ull c1, Ull c0, Ushort pattern)
{
  Uint index1, index0;

  switch (op_cx) {
  case OP_NOP:
    if (ex)
      *ex = 3; /* for bsim */
    break;
  case OP_CEXE:
    index1 = ((c3>>32&1)<<3)|((c2>>32&1)<<2)|((c1>>32&1)<<1)|(c0>>32&1);
    index0 = ((c3    &1)<<3)|((c2    &1)<<2)|((c1    &1)<<1)|(c0    &1);
    *ex = 0;
    if (pattern>>index1&1) *ex |= 2;
    if (pattern>>index0&1) *ex |= 1;
    break;
  default:
    printf("emax7lib: cex: undefined op_cx=%d\n", op_cx);
    break;
  }  
}

void /*__attribute__((always_inline))*/
ex4(Uint op_ex1, Ull *d, Ull *r1, Uint exp1, Ull *r2, Uint exp2, Ull *r3, Uint exp3, Uint op_ex2, Ull *r4, Uint op_ex3, Ull *r5)
{
  switch (op_ex1) {
  case OP_SFMA: /* 3in 8bit*32 stochastic r1+r2*r3 -> 8bit */
    exe(op_ex1, (d+0), (Ull)r1, exp1, *(r2+0), exp2, *(r3+0), exp3, OP_NOP, (Ull)r4, OP_NOP, (Ull)r5);
    exe(op_ex1, (d+0), *(d+0),  exp1, *(r2+1), exp2, *(r3+1), exp3, OP_NOP, (Ull)r4, OP_NOP, (Ull)r5);
    exe(op_ex1, (d+0), *(d+0),  exp1, *(r2+2), exp2, *(r3+2), exp3, OP_NOP, (Ull)r4, OP_NOP, (Ull)r5);
    exe(op_ex1, (d+0), *(d+0),  exp1, *(r2+3), exp2, *(r3+3), exp3, OP_NOP, (Ull)r4, OP_NOP, (Ull)r5);
    break;
  case OP_NOP:
  case OP_CFMA: /* 3in [idx|32bit]*2 (idx2==idx3)?r1+r2*r3:r1 */
  case OP_FMA:  /* 3in 32bit*2 floating-point r1+r2*r3 */
  case OP_FMS:  /* 3in 32bit*2 floating-point r1-r2*r3 */
  case OP_FML:  /* 2in 32bit*2 floating-point r1*r2 */
  case OP_FAD:  /* 2in 32bit*2 floating-point r1+r2 */
  case OP_FML3: /* 3in 32bit*2 floating-point r1*r2[idx:r3]) */
  case OP_ADD3: /* 3in 32bit*2 fixed-point r1+(r2+r3) */
  case OP_SUB3: /* 3in 32bit*2 fixed-point r1-(r2+r3) */
  case OP_ADD:  /* 2in 32bit*2 fixed-point r1+r2 */
  case OP_SUB:  /* 2in 32bit*2 fixed-point r1-r2 */
    exe(op_ex1, (d+0), *(r1+0), exp1, *(r2+0), exp2, *(r3+0), exp3, OP_NOP, 0LL, OP_NOP, 0LL);
    exe(op_ex1, (d+1), *(r1+1), exp1, *(r2+1), exp2, *(r3+1), exp3, OP_NOP, 0LL, OP_NOP, 0LL);
    exe(op_ex1, (d+2), *(r1+2), exp1, *(r2+2), exp2, *(r3+2), exp3, OP_NOP, 0LL, OP_NOP, 0LL);
    exe(op_ex1, (d+3), *(r1+3), exp1, *(r2+3), exp2, *(r3+3), exp3, OP_NOP, 0LL, OP_NOP, 0LL);
    break;
  default:
    printf("emax7lib: ex4: undefined op_ex1=%d\n", op_ex1);
    break;
  }

  switch (op_ex2) {
  case OP_NOP:
    break;
  default:
    printf("emax7lib: ex4: illegal op_ex2=%d\n", op_ex2);
    break;
  }

  switch (op_ex3) {
  case OP_NOP:
    break;
  default:
    printf("emax7lib: ex4: illegal op_ex3=%d\n", op_ex3);
    break;
  }
}

int convf32tou7(Uchar *out, float in)
{
  //  convf32tou7 e=126     0.992 -> s0111111  0111111111111111111111111111111111111111111111111111111111111111
  //  convf32tou7 e=126 f=0 0.500 -> s0100000  0000000000000000000000000000000011111111111111111111111111111111
  //  convf32tou7 e=125 f=0 0.250 -> s0010000  0000000000000000000000000000000000000000000000001111111111111111
  //  convf32tou7 e=124 f=0 0.125 -> s0001000  0000000000000000000000000000000000000000000000000000000011111111
  //  convf32tou7 e=123 f=0 0.062 -> s0000100  0000000000000000000000000000000000000000000000000000000000001111
  //  convf32tou7 e=122 f=0 0.031 -> s0000010  0000000000000000000000000000000000000000000000000000000000000011
  //  convf32tou7 e=121 f=0 0.016 -> s0000001  0000000000000000000000000000000000000000000000000000000000000001
  //                        0.000 -> s0000000  0000000000000000000000000000000000000000000000000000000000000000
  f32bit in_f32;
  wu7bit out_u7;

  *(float*)&in_f32 = in;

  out_u7.s = in_f32.s;
  out_u7.b = 0;

  in = abs(in);
  if  (in >= 1.0) out_u7.e = 63;    /* �軻����6bitɽ��(-1.0+1.0) */
  else            out_u7.e = in*64; /* number of 1 */    

  *out = *(Uchar*)&out_u7;
//printf("%7.4f -> %02.2x\n", *(float*)&in_f32, *out);
}

int convf32tou8(Uchar *out, float in)
{
  f32bit in_f32;
  wu8bit out_u8;

  *(float*)&in_f32 = in;

  out_u8.s = in_f32.s;

  in = abs(in);
  if  (in >= 2.0) out_u8.e = 127;   /* �軻����6bitɽ��(-1.0+1.0) */
  else            out_u8.e = in*64; /* number of 1 */    

  *out = *(Uchar*)&out_u8;
//printf("%7.4f -> %02.2x\n", *(float*)&in_f32, *out);
}

int convu8tof32(float *out, Uchar in)
{
  wu8bit in_u8;
  f32bit out_f32;

  *(Uchar*)&in_u8 = in;
  *(float*)&out_f32 = (float)in_u8.e/64; /* 6bitɽ��(-2.0+2.0) */
  out_f32.s = in_u8.s;
  *out = *(float*)&out_f32;

//printf("%02.2x -> %7.4f\n", in, *out);
}

Ull urand(int no)
{
  static Ull urand_seed[8]
    = {0xc3c3c3c3a5a5a5a5LL, 0x123456789abcdef0LL, 0xe1e1e1e1d4d4d4d4LL, 0x8888777766665555LL,
       0x8787878796969696LL, 0xfedcba9876543210LL, 0x5a5a5a5a3c3c3c3cLL, 0xbbbbccccddddeeeeLL};
  Ull retval = urand_seed[no];

//urand_seed = urand_seed * 1103515245LL + 12345LL;

  urand_seed[no] ^= (urand_seed[no]<<29);
  urand_seed[no] ^= (urand_seed[no]>>27);
  urand_seed[no] ^= (urand_seed[no]<<37);
  return (retval);
}

Ull shfl(Ull in, Ull r)
{
  int i;
  for (i=0; i<32; i++) {
    if (r&(1LL<<(i+16)))
      in = (in&~(1LL<<(i+32)|1LL<<i)) | (in>>i&1)<<(i+32) | (in>>(i+32)&1)<<i;
  }
  for (i=0; i<48; i++) {
    if (r&(1LL<<(i+8)))
      in = (in&~(1LL<<(i+16)|1LL<<i)) | (in>>i&1)<<(i+16) | (in>>(i+16)&1)<<i;
  }
  for (i=0; i<56; i++) {
    if (r&(1LL<<(i+4)))
      in = (in&~(1LL<<(i+ 8)|1LL<<i)) | (in>>i&1)<<(i+ 8) | (in>>(i+ 8)&1)<<i;
  }
  for (i=0; i<60; i++) {
    if (r&(1LL<<(i+2)))
      in = (in&~(1LL<<(i+ 4)|1LL<<i)) | (in>>i&1)<<(i+ 4) | (in>>(i+ 4)&1)<<i;
  }
  for (i=0; i<62; i++) {
    if (r&(1LL<<(i+1)))
      in = (in&~(1LL<<(i+ 2)|1LL<<i)) | (in>>i&1)<<(i+ 2) | (in>>(i+ 2)&1)<<i;
  }
  for (i=0; i<63; i++) {
    if (r&(1LL<<(i+0)))
      in = (in&~(1LL<<(i+ 1)|1LL<<i)) | (in>>i&1)<<(i+ 1) | (in>>(i+ 1)&1)<<i;
  }
  return(in);
}

void x11_softu64_dist(float, float);
int softu64(int stage, Ull *o1, Ull *o2, Ull *o3, Ull r1, Ull r2, Ull r3, Ull r4) /* o <- s1 + s2 * s3 */
     /* stage:1 stage_2 in EXEC:  r2*r3 64bit*2  -> *o1 32bit*8 b mult     */
     /* stage:2 stage_3 in EXEC:  *o1,r4 32bit*8 -> *o2 8bit+8bit count up */
     /* stage:3 stage_4 in EXEC:  r1 + *o2��     -> *o3 8bit               */
{
  int i, j;
  Ull u[8];
  Ull ss[8];
  Ull s2[8], s3[8];
  int pc, nc; /* number of 1 */
  int os, oc;

//#define SPU_DATA_BITS 31
//#define SPU_DATA_DIST 2
//#define SPU_COUT_BITS 31
#define SPU_DATA_BITS 15
#define SPU_DATA_DIST 4
#define SPU_COUT_BITS 12

  switch (stage) {
  case 1: /* stage2 */
    for (i=0; i<8; i++) /* s2 * s3 -> ad2 */
      u[i] = urand(i);
    for (i=0; i<8; i++) { /* s2 * s3 -> ad2 */
      ss[i] = (r2>>(i*8+7))&1 ^ (r3>>(i*8+7))&1;
  int s2e   = (r2>>(i*8))&0x7f; s2e = s2e<SPU_DATA_BITS?s2e:SPU_DATA_BITS;
  int s3e   = (r3>>(i*8))&0x7f; s3e = s3e<SPU_DATA_BITS?s3e:SPU_DATA_BITS;
#if 0
      s2[i] = (Ull)0x7fffffffffffffffLL>>(63-s2e); //�軻��6bit*6bit->6bit
      s3[i] = (Ull)0x7fffffffffffffffLL>>(63-s3e); //�軻��6bit*6bit->6bit
      // ����64bit��ǥ���åե�
      s2[i] = shfl(s2[i], u[2]);
      s3[i] = shfl(s3[i], u[3]);
#else
      // �����SPU_DATA_WIDTH bit��˻���.�����ͤ�6bit���ۤȤ��15�ʲ��Ǥ��뤳�Ȥ�����(63�᤯�ʤ�ȸ������Ф�Ϥ�)
      s2[i] = 0LL;
      s3[i] = 0LL;
      for (j=0; j<SPU_COUT_BITS; j++) {
	int k = j * SPU_DATA_DIST; /* SPU_DATA_BITS=15�ʤ�4bit�� */
	s2[i] |= ((u[(i+0)%8]>>k&SPU_DATA_BITS)<=s2e)<<j;
	s3[i] |= ((u[(i+1)%8]>>k&SPU_DATA_BITS)<=s3e)<<j;
      }
      //printf("%08.8x_%08.8x %08.8x_%08.8x %d:%08.8x %d:%08.8x\n", (Uint)(u2>>32), (Uint)u2, (Uint)(u3>>32), (Uint)u3, s2e, (Uint)s2[i], s3e, (Uint)s3[i]);
#endif
      // s2*s3 �����Ǥ�stochastic�軻
      o1[i] = s2[i] & s3[i];                         // 1*1=1�ˤʤ� �ºݤϾ��SPU_DATA_BITS�Τ�AND
      o1[i] = ss[i]<<63|(o1[i]&0x7fffffffffffffffLL);// stage2�ν��Ϥ�(��Ƭ���bit|SPU_DATA_BITS bit) * 8
    }
    break;
  case 2: /* stage3 */
    pc = 0;
    nc = 0;
    // ����/������롼�פ��Ȥˡ�������ʬ�򥹥ʥåץ���å�
    for (j=0; j<SPU_COUT_BITS; j++) {
      for (i=0; i<8; i++) { /* s2 * s3 -> ad2 */
	if (!(o1[i]>>63)) pc += (o1[i] & (1LL<<j))!=0;
	else              nc += (o1[i] & (1LL<<j))!=0;
      }
    }
    pc = pc>>r4; // r4=3 for MNIST/CIFAR10
    nc = nc>>r4; // r4=2 for test021
    *o2 = (Ull)(pc&0xffff)<<32 | (Ull)(nc&0xffff);
    break;
  case 3: /* stage4 */
    pc = *o2>>32&0xffff; /* high */
    nc = *o2    &0xffff; /* low */
    // s1�򤵤�˲û�
    if (!(r1&0x80)) pc += (r1&0x7f); /* merge pos s1 s1.e�Ϻ���7bit */
    else            nc += (r1&0x7f); /* merge neg s1 s1.e�Ϻ���7bit */
    // ����������βû�(s1:7bit + s2*s3:6bit->7bit)
    if (pc >= nc) {
      os = 0x00; /* pos */
      oc = pc-nc; /* # of 1 */
    }
    else {
      os = 0x80; /* neg */
      oc = nc-pc; /* # of 1 */
    }
    if (oc >= 128) oc = 127;
    *o3 = os|oc;
#if !defined(ARMSIML) && defined(TRACE_SPIKE)
    if (enable_x11) {
      int i;
      Uchar r2_u8;
      Uchar r3_u8;
      float r1_f32;
      float r2_f32;
      float r3_f32;
      float o3_f32;
      convu8tof32(&o3_f32, *(Uchar*)o3);   /* for graph */
      convu8tof32(&r1_f32, *(Uchar*)&r1); /* for graph */
      for (i=0; i<8; i++) { /* s2 * s3 -> ad2 */
	r2_u8 = r2>>(i*8)&0xff;
	r3_u8 = r3>>(i*8)&0xff;
	convu8tof32(&r2_f32, r2_u8); /* for graph */
	convu8tof32(&r3_f32, r3_u8); /* for graph */
	r1_f32 += r2_f32*r3_f32;
      }
      x11_softu64_dist(r1_f32, o3_f32);
    }
#endif
    break;
  }

  return (0);
}

int  /*__attribute__((always_inline))*/
exe(Uint op_ex1, Ull *d, Ull s1, Uint exp1, Ull s2, Uint exp2, Ull s3, Uint exp3, Uint op_ex2, Ull r4, Uint op_ex3, Ull r5)
{
  /* return 0:normal, 1:OP_WHILE breaks */
  union { Uint i; float f; } f3, f2, f1, f0;
  Ull r1, r2, r3;
  Ull t3, t2, t1, t0;
  Ull ro00, ro01, ro02, ro10, ro11, ro12;
  Ull c1, c0;
  Ull ex1_outd;
  Ull ex1_outd_sfma[8];
  Ull ex2_outd;
  int retval = 0;
  float convi4f32[16] = {-8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

  r1 = exm(s1, exp1);
  r2 = exm(s2, exp2);
  r3 = exm(s3, exp3);

  switch (op_ex1) {
  case OP_NOP:
    ex1_outd = r1;
    break;
  case OP_WHILE: /* emax7nc��lib�Ȥ��Ƥϻ��Ѥ���,bsim/emax7.c��siml�˻��� */
    t0 = (r1&0x00000000ffffffffLL)+(r2&0x00000000ffffffffLL);
    t0 &= 0x00000000ffffffffLL;
    ex1_outd = t0;
    if (t0==0) retval = 1;
    break;
  case OP_FOR: /* emax7nc��lib�Ȥ��Ƥϻ��Ѥ���,bsim/emax7.c��siml�˻��� */
    t0 = (r1&0x00000000ffffffffLL)+(r2&0x00000000ffffffffLL);
    t0 &= 0x00000000ffffffffLL;
    ex1_outd = t0;
    if (t0==0) retval = 1;
    break;
  case OP_SFMA: /* 3in 8bit*32 stochastic r1+r2*r3 -> 8bit */
    softu64(1, ex1_outd_sfma, NULL, NULL, r1, r2, r3, r4);
    break;
  case OP_CFMA: /* 3in [idx|32bit]*2 (idx2==idx3)?r1+r2*r3:r1 */
    f1.i = (Uint)(r1);
    f2.i = (Uint)(r2>>32);
    f3.i = (Uint)(r3>>32);
    if (f2.i != -1 && f2.i == f3.i) {
      f2.i = (Uint)(r2);
      f3.i = (Uint)(r3);
      f0.f = f1.f + (f2.f * f3.f);
    }
    else {
      f0.f = f1.f;
    }
    t0 = f0.i;
    ex1_outd = t0;
    break;
  case OP_FMA: /* 3in 32bit*2 floating-point r1+r2*r3 */
  case OP_FMS: /* 3in 32bit*2 floating-point r1-r2*r3 */
    /* *(double*)&ex1_outd = *(double*)&r1 + (*(double*)&r2 * *(double*)&r3);*/
    f1.i = (Uint)(r1>>32);
    f2.i = (Uint)(r2>>32)^(op_ex1==OP_FMA?0:0x80000000);
    f3.i = (Uint)(r3>>32);
    f0.f = f1.f + (f2.f * f3.f);
    t2 = f0.i;
    f1.i = (Uint)(r1);
    f2.i = (Uint)(r2)^(op_ex1==OP_FMA?0:0x80000000);
    f3.i = (Uint)(r3);
    f0.f = f1.f + (f2.f * f3.f);
    t0 = f0.i;
    ex1_outd = (t2<<32)|(t0);
    break;
  case OP_FML: /* 2in 32bit*2 floating-point r1*r2 */
    /* *(double*)&ex1_outd = *(double*)&r1 * *(double*)&r2;*/
    f1.i = (Uint)(r1>>32);
    f2.i = (Uint)(r2>>32);
    f0.f = f1.f * f2.f;
    t2 = f0.i;
    f1.i = (Uint)(r1);
    f2.i = (Uint)(r2);
    f0.f = f1.f * f2.f;
    t0 = f0.i;
    ex1_outd = (t2<<32)|(t0);
    break;
  case OP_FAD: /* 2in 32bit*2 floating-point r1+r2 */
    /* *(double*)&ex1_outd = *(double*)&r1 + *(double*)&r2;*/
    f1.i = (Uint)(r1>>32);
    f2.i = (Uint)(r2>>32);
    f0.f = f1.f + f2.f;
    t2 = f0.i;
    f1.i = (Uint)(r1);
    f2.i = (Uint)(r2);
    f0.f = f1.f + f2.f;
    t0 = f0.i;
    ex1_outd = (t2<<32)|(t0);
    break;
  case OP_FML3: /* 3in 32bit*2 floating-point r1*r2[idx:r3] */
    /* *(double*)&ex1_outd = *(double*)&r1 * *(int*)&r2[idx];*/
    f1.i = (Uint)(r1>>32);
    f2.i = (Uint)(r2>>32);
    f3.i = (Uint)(r3>>32);
    f0.f = f1.f * convi4f32[(f2.i>>((f3.i&7)*4))&0xf];
    t2 = f0.i;
    f1.i = (Uint)(r1);
    f2.i = (Uint)(r2);
    f3.i = (Uint)(r3);
    f0.f = f1.f * convi4f32[(f2.i>>((f3.i&7)*4))&0xf];
    t0 = f0.i;
    ex1_outd = (t2<<32)|(t0);
    break;
  case OP_ADD3: /* 3in 32bit*2 integer add s1+(s2+s3) */
    t2 = (r1>>32&0x00000000ffffffffLL)+((r2>>32&0x00000000ffffffffLL)+(r3>>32&0x00000000ffffffffLL));
    t2 &= 0x00000000ffffffffLL;
    t0 = (r1    &0x00000000ffffffffLL)+((r2    &0x00000000ffffffffLL)+(r3    &0x00000000ffffffffLL));
    t0 &= 0x00000000ffffffffLL;
    ex1_outd = (t2<<32)|(t0);
    break;
  case OP_SUB3: /* 3in 32bit*2 integer subtract s1-(s2+s3) */
    t2 = (r1>>32&0x00000000ffffffffLL)-((r2>>32&0x00000000ffffffffLL)+(r3>>32&0x00000000ffffffffLL));
    t2 &= 0x00000000ffffffffLL;
    t0 = (r1    &0x00000000ffffffffLL)-((r2    &0x00000000ffffffffLL)+(r3    &0x00000000ffffffffLL));
    t0 &= 0x00000000ffffffffLL;
    ex1_outd = (t2<<32)|(t0);
    break;
  case OP_ADD: /* 2in 32bit*2 integer add s1+s2 */
    t2 = (r1>>32&0x00000000ffffffffLL)+(r2>>32&0x00000000ffffffffLL);
    t2 &= 0x00000000ffffffffLL;
    t0 = (r1    &0x00000000ffffffffLL)+(r2    &0x00000000ffffffffLL);
    t0 &= 0x00000000ffffffffLL;
    ex1_outd = (t2<<32)|(t0);
    break;
  case OP_SUB: /* 2in 32bit*2 integer subtract s1-s2 */
    t2 = (r1>>32&0x00000000ffffffffLL)-(r2>>32&0x00000000ffffffffLL);
    t2 &= 0x00000000ffffffffLL;
    t0 = (r1    &0x00000000ffffffffLL)-(r2    &0x00000000ffffffffLL);
    t0 &= 0x00000000ffffffffLL;
    ex1_outd = (t2<<32)|(t0);
    break;
  case OP_CMP_EQ: /* 2in 32bit*2 compare and set 1*2bit-CC */
    c1 = (r1>>32&0x00000000ffffffffLL) == (r2>>32&0x00000000ffffffffLL);
    c0 = (r1    &0x00000000ffffffffLL) == (r2    &0x00000000ffffffffLL);
    ex1_outd = (c1<<32)|c0;
    break;
  case OP_CMP_NE: /* 2in 32bit*2 compare and set 1*2bit-CC */
    c1 = (r1>>32&0x00000000ffffffffLL) != (r2>>32&0x00000000ffffffffLL);
    c0 = (r1    &0x00000000ffffffffLL) != (r2    &0x00000000ffffffffLL);
    ex1_outd = (c1<<32)|c0;
    break;
  case OP_CMP_LT: /* 2in 32bit*2 compare and set 1*2bit-CC */
    c1 = (r1>>32&0x00000000ffffffffLL) < (r2>>32&0x00000000ffffffffLL);
    c0 = (r1    &0x00000000ffffffffLL) < (r2    &0x00000000ffffffffLL);
    ex1_outd = (c1<<32)|c0;
    break;
  case OP_CMP_LE: /* 2in 32bit*2 compare and set 1*2bit-CC */
    c1 = (r1>>32&0x00000000ffffffffLL) <= (r2>>32&0x00000000ffffffffLL);
    c0 = (r1    &0x00000000ffffffffLL) <= (r2    &0x00000000ffffffffLL);
    ex1_outd = (c1<<32)|c0;
    break;
  case OP_CMP_GT: /* 2in 32bit*2 compare and set 1*2bit-CC */
    c1 = (r1>>32&0x00000000ffffffffLL) > (r2>>32&0x00000000ffffffffLL);
    c0 = (r1    &0x00000000ffffffffLL) > (r2    &0x00000000ffffffffLL);
    ex1_outd = (c1<<32)|c0;
    break;
  case OP_CMP_GE: /* 2in 32bit*2 compare and set 1*2bit-CC */
    c1 = (r1>>32&0x00000000ffffffffLL) >= (r2>>32&0x00000000ffffffffLL);
    c0 = (r1    &0x00000000ffffffffLL) >= (r2    &0x00000000ffffffffLL);
    ex1_outd = (c1<<32)|c0;
    break;
  case OP_CMOV: /* 2in 32bit*2 word-wise(w1/w0) conditional move */
    c1 = r1>>32&1;
    c0 = r1    &1;
    t2 = c1 ? (r2&0xffffffff00000000LL) : (r3&0xffffffff00000000LL);
    t0 = c0 ? (r2&0x00000000ffffffffLL) : (r3&0x00000000ffffffffLL);
    ex1_outd = t2 | t0;
    break;
  case OP_MAUH3: /* 3in 16bit*4 r1.pos+(r2.pos+r3.pos) */
    t3 = (r1>>48&0x000000000000ffffLL)+((r2>>48&0x000000000000ffffLL)+(r3>>48&0x000000000000ffffLL));
    if (t3 > 0x000000000000ffffLL) t3 = 0x000000000000ffffLL;
    t2 = (r1>>32&0x000000000000ffffLL)+((r2>>32&0x000000000000ffffLL)+(r3>>32&0x000000000000ffffLL));
    if (t2 > 0x000000000000ffffLL) t2 = 0x000000000000ffffLL;
    t1 = (r1>>16&0x000000000000ffffLL)+((r2>>16&0x000000000000ffffLL)+(r3>>16&0x000000000000ffffLL));
    if (t1 > 0x000000000000ffffLL) t1 = 0x000000000000ffffLL;
    t0 = (r1    &0x000000000000ffffLL)+((r2    &0x000000000000ffffLL)+(r3    &0x000000000000ffffLL));
    if (t0 > 0x000000000000ffffLL) t0 = 0x000000000000ffffLL;
    ex1_outd = (t3<<48)|(t2<<32)|(t1<<16)|(t0);
    break;
  case OP_MAUH: /* 2in 16bit*4 r1.pos+r2.pos */
    t3 = (r1>>48&0x000000000000ffffLL)+(r2>>48&0x000000000000ffffLL);
    if (t3 > 0x000000000000ffffLL) t3 = 0x000000000000ffffLL;
    t2 = (r1>>32&0x000000000000ffffLL)+(r2>>32&0x000000000000ffffLL);
    if (t2 > 0x000000000000ffffLL) t2 = 0x000000000000ffffLL;
    t1 = (r1>>16&0x000000000000ffffLL)+(r2>>16&0x000000000000ffffLL);
    if (t1 > 0x000000000000ffffLL) t1 = 0x000000000000ffffLL;
    t0 = (r1    &0x000000000000ffffLL)+(r2    &0x000000000000ffffLL);
    if (t0 > 0x000000000000ffffLL) t0 = 0x000000000000ffffLL;
    ex1_outd = (t3<<48)|(t2<<32)|(t1<<16)|(t0);
    break;
  case OP_MSUH3: /* 3in 16bit*4 r1.pos-(r2.pos+r3.pos) */
    t3 = (r1>>48&0x000000000000ffffLL)-((r2>>48&0x000000000000ffffLL)+(r3>>48&0x000000000000ffffLL));
    if (t3 > 0x000000000000ffffLL) t3 = 0x0000000000000000LL;
    t2 = (r1>>32&0x000000000000ffffLL)-((r2>>32&0x000000000000ffffLL)+(r3>>32&0x000000000000ffffLL));
    if (t2 > 0x000000000000ffffLL) t2 = 0x0000000000000000LL;
    t1 = (r1>>16&0x000000000000ffffLL)-((r2>>16&0x000000000000ffffLL)+(r3>>16&0x000000000000ffffLL));
    if (t1 > 0x000000000000ffffLL) t1 = 0x0000000000000000LL;
    t0 = (r1    &0x000000000000ffffLL)-((r2    &0x000000000000ffffLL)+(r3    &0x000000000000ffffLL));
    if (t0 > 0x000000000000ffffLL) t0 = 0x0000000000000000LL;
    ex1_outd = (t3<<48)|(t2<<32)|(t1<<16)|(t0);
    break;
  case OP_MSUH: /* 2in 16bit*4 r1.pos-r2.pos */
    t3 = (r1>>48&0x000000000000ffffLL)-(r2>>48&0x000000000000ffffLL);
    if (t3 > 0x000000000000ffffLL) t3 = 0x0000000000000000LL;
    t2 = (r1>>32&0x000000000000ffffLL)-(r2>>32&0x000000000000ffffLL);
    if (t2 > 0x000000000000ffffLL) t2 = 0x0000000000000000LL;
    t1 = (r1>>16&0x000000000000ffffLL)-(r2>>16&0x000000000000ffffLL);
    if (t1 > 0x000000000000ffffLL) t1 = 0x0000000000000000LL;
    t0 = (r1    &0x000000000000ffffLL)-(r2    &0x000000000000ffffLL);
    if (t0 > 0x000000000000ffffLL) t0 = 0x0000000000000000LL;
    ex1_outd = (t3<<48)|(t2<<32)|(t1<<16)|(t0);
    break;
  case OP_MLUH: /* (11bit*4)*9bit r1.pos*r2.pos */
    t3 = (r1>>48&0x00000000000007ffLL)*(r2&0x00000000000001ffLL);
    if (t3 > 0x000000000000ffffLL) t3 = 0x000000000000ffffLL;
    t2 = (r1>>32&0x00000000000007ffLL)*(r2&0x00000000000001ffLL);
    if (t2 > 0x000000000000ffffLL) t2 = 0x000000000000ffffLL;
    t1 = (r1>>16&0x00000000000007ffLL)*(r2&0x00000000000001ffLL);
    if (t1 > 0x000000000000ffffLL) t1 = 0x000000000000ffffLL;
    t0 = (r1    &0x00000000000007ffLL)*(r2&0x00000000000001ffLL);
    if (t0 > 0x000000000000ffffLL) t0 = 0x000000000000ffffLL;
    ex1_outd = (t3<<48)|(t2<<32)|(t1<<16)|(t0);
    break;
  case OP_MMRG: /* 3in 8bit*2 r1.b4|r2.b4|r3.b4|0->w1, r1.b0|r2.b0|r3.b0|0->w0 */
    ex1_outd = ((r1&0x000000ff00000000LL)<<24) | ((r2&0x000000ff00000000LL)<<16) | ((r3&0x000000ff00000000LL)<<8)
             | ((r1&0x00000000000000ffLL)<<24) | ((r2&0x00000000000000ffLL)<<16) | ((r3&0x00000000000000ffLL)<<8);
    break;
  case OP_MSSAD: /* 2in 16bit*4 8bit*8 r1.h3+df(r2.b7,r3.b7)+df(r2.b6,r3.b6)->d.h3
                                       r1.h2+df(r2.b5,r3.b5)+df(r2.b4,r3.b4)->d.h2
                                       r1.h1+df(r2.b3,r3.b3)+df(r2.b2,r3.b2)->d.h1
                                       r1.h0+df(r2.b1,r3.b1)+df(r2.b0,r3.b0)->d.h0 */
    t3 = (r1>>48&0x000000000000ffffLL) + ad(r2>>56&0x00000000000000ffLL, r3>>56&0x00000000000000ffLL) + ad(r2>>48&0x00000000000000ffLL, r3>>48&0x00000000000000ffLL);
    if (t3 > 0x000000000000ffffLL) t3 = 0x000000000000ffffLL;
    t2 = (r1>>32&0x000000000000ffffLL) + ad(r2>>40&0x00000000000000ffLL, r3>>40&0x00000000000000ffLL) + ad(r2>>32&0x00000000000000ffLL, r3>>32&0x00000000000000ffLL);
    if (t2 > 0x000000000000ffffLL) t2 = 0x000000000000ffffLL;
    t1 = (r1>>16&0x000000000000ffffLL) + ad(r2>>24&0x00000000000000ffLL, r3>>24&0x00000000000000ffLL) + ad(r2>>16&0x00000000000000ffLL, r3>>16&0x00000000000000ffLL);
    if (t1 > 0x000000000000ffffLL) t1 = 0x000000000000ffffLL;
    t0 = (r1    &0x000000000000ffffLL) + ad(r2>> 8&0x00000000000000ffLL, r3>> 8&0x00000000000000ffLL) + ad(r2    &0x00000000000000ffLL, r3    &0x00000000000000ffLL);
    if (t0 > 0x000000000000ffffLL) t0 = 0x000000000000ffffLL;
    ex1_outd = (t3<<48)|(t2<<32)|(t1<<16)|(t0);
    break;
  case OP_MSAD: /* 2in 16bit*4 8bit*8 df(r1.b7,r2.b7)+df(r1.b6,r2.b6)->d.h3
                                      df(r1.b5,r2.b5)+df(r1.b4,r2.b4)->d.h2
                                      df(r1.b3,r2.b3)+df(r1.b2,r2.b2)->d.h1
                                      df(r1.b1,r2.b1)+df(r1.b0,r2.b0)->d.h0 */
    t3 = ad(r1>>56&0x00000000000000ffLL, r2>>56&0x00000000000000ffLL) + ad(r1>>48&0x00000000000000ffLL, r2>>48&0x00000000000000ffLL);
    if (t3 > 0x000000000000ffffLL) t3 = 0x000000000000ffffLL;
    t2 = ad(r1>>40&0x00000000000000ffLL, r2>>40&0x00000000000000ffLL) + ad(r1>>32&0x00000000000000ffLL, r2>>32&0x00000000000000ffLL);
    if (t2 > 0x000000000000ffffLL) t2 = 0x000000000000ffffLL;
    t1 = ad(r1>>24&0x00000000000000ffLL, r2>>24&0x00000000000000ffLL) + ad(r1>>16&0x00000000000000ffLL, r2>>16&0x00000000000000ffLL);
    if (t1 > 0x000000000000ffffLL) t1 = 0x000000000000ffffLL;
    t0 = ad(r1>> 8&0x00000000000000ffLL, r2>> 8&0x00000000000000ffLL) + ad(r1    &0x00000000000000ffLL, r2    &0x00000000000000ffLL);
    if (t0 > 0x000000000000ffffLL) t0 = 0x000000000000ffffLL;
    ex1_outd = (t3<<48)|(t2<<32)|(t1<<16)|(t0);
    break;
  case OP_MINL3: /* 3in 16bit*4 (r3.h3<r3.h2)?r1.h3|r3.h3:r2.h3|r3.h2->d.w1
                                (r3.h1<r3.h0)?r1.h1|r3.h1:r2.h1|r3.h0->d.w0 */
    t3 = r3>>48&0x000000000000ffffLL;
    t2 = r3>>32&0x000000000000ffffLL;
    t1 = r3>>16&0x000000000000ffffLL;
    t0 = r3    &0x000000000000ffffLL;
    if (t3<t2) t2 = (r1&0xffff000000000000LL)|(r3>>16&0x0000ffff00000000LL);
    else       t2 = (r2&0xffff000000000000LL)|(r3    &0x0000ffff00000000LL);
    if (t1<t0) t0 = (r1&0x00000000ffff0000LL)|(r3>>16&0x000000000000ffffLL);
    else       t0 = (r2&0x00000000ffff0000LL)|(r3    &0x000000000000ffffLL);
    ex1_outd = t2 | t0;
    break;
  case OP_MINL: /* 2in 16bit*4 (r1.h2<r2.h2)?r1.w1:r2.w1->d.w1
	                       (r1.h0<r2.h0)?r1.w0:r2.w0->d.w0 */
    if ((r1&0x0000ffff00000000LL)<(r2&0x0000ffff00000000LL)) t2 = r1&0xffffffff00000000LL;
    else                                                     t2 = r2&0xffffffff00000000LL;
    if ((r1&0x000000000000ffffLL)<(r2&0x000000000000ffffLL)) t0 = r1&0x00000000ffffffffLL;
    else                                                     t0 = r2&0x00000000ffffffffLL;
    ex1_outd = t2 | t0;
   break;
  case OP_MH2BW: /* 2in 16bit*4 r1.b6|r1.b4|r2.b6|r2.b4|r1.b2|r1.b0|r2.b2|r2.b0 */
    ex1_outd = (((r1>>48&0x000000000000ff00LL) ? 255 : (r1>>48&0x00000000000000ffLL))<<56)
             | (((r1>>32&0x000000000000ff00LL) ? 255 : (r1>>32&0x00000000000000ffLL))<<48)
             | (((r2>>48&0x000000000000ff00LL) ? 255 : (r2>>48&0x00000000000000ffLL))<<40)
             | (((r2>>32&0x000000000000ff00LL) ? 255 : (r2>>32&0x00000000000000ffLL))<<32)
             | (((r1>>16&0x000000000000ff00LL) ? 255 : (r1>>16&0x00000000000000ffLL))<<24)
             | (((r1    &0x000000000000ff00LL) ? 255 : (r1    &0x00000000000000ffLL))<<16)
             | (((r2>>16&0x000000000000ff00LL) ? 255 : (r2>>16&0x00000000000000ffLL))<< 8)
             | (((r2    &0x000000000000ff00LL) ? 255 : (r2    &0x00000000000000ffLL))    );
    break;
  case OP_MCAS: /* 2in 16bit*2 (r1.h2<r2.h2)?0:0xff->d.b1
                               (r1.h0<r2.h0)?0:0xff->d.b0 */
    t2 = ((r1&0x0000ffff00000000LL)<(r2&0x0000ffff00000000LL))?0:0x000000ff00000000LL;
    t0 = ((r1&0x000000000000ffffLL)<(r2&0x000000000000ffffLL))?0:0x00000000000000ffLL;
    ex1_outd = t2 | t0;
    break;
  case OP_MMID3: /* 3in 8bit*8 bytewise compare and output middle */
    t1 = ((r1&0xff00000000000000LL)<(r2&0xff00000000000000LL)?(r1&0xff00000000000000LL):(r2&0xff00000000000000LL))
       | ((r1&0x00ff000000000000LL)<(r2&0x00ff000000000000LL)?(r1&0x00ff000000000000LL):(r2&0x00ff000000000000LL))
       | ((r1&0x0000ff0000000000LL)<(r2&0x0000ff0000000000LL)?(r1&0x0000ff0000000000LL):(r2&0x0000ff0000000000LL))
       | ((r1&0x000000ff00000000LL)<(r2&0x000000ff00000000LL)?(r1&0x000000ff00000000LL):(r2&0x000000ff00000000LL))
       | ((r1&0x00000000ff000000LL)<(r2&0x00000000ff000000LL)?(r1&0x00000000ff000000LL):(r2&0x00000000ff000000LL))
       | ((r1&0x0000000000ff0000LL)<(r2&0x0000000000ff0000LL)?(r1&0x0000000000ff0000LL):(r2&0x0000000000ff0000LL))
       | ((r1&0x000000000000ff00LL)<(r2&0x000000000000ff00LL)?(r1&0x000000000000ff00LL):(r2&0x000000000000ff00LL))
       | ((r1&0x00000000000000ffLL)<(r2&0x00000000000000ffLL)?(r1&0x00000000000000ffLL):(r2&0x00000000000000ffLL));
    t2 = ((r1&0xff00000000000000LL)>(r2&0xff00000000000000LL)?(r1&0xff00000000000000LL):(r2&0xff00000000000000LL))
       | ((r1&0x00ff000000000000LL)>(r2&0x00ff000000000000LL)?(r1&0x00ff000000000000LL):(r2&0x00ff000000000000LL))
       | ((r1&0x0000ff0000000000LL)>(r2&0x0000ff0000000000LL)?(r1&0x0000ff0000000000LL):(r2&0x0000ff0000000000LL))
       | ((r1&0x000000ff00000000LL)>(r2&0x000000ff00000000LL)?(r1&0x000000ff00000000LL):(r2&0x000000ff00000000LL))
       | ((r1&0x00000000ff000000LL)>(r2&0x00000000ff000000LL)?(r1&0x00000000ff000000LL):(r2&0x00000000ff000000LL))
       | ((r1&0x0000000000ff0000LL)>(r2&0x0000000000ff0000LL)?(r1&0x0000000000ff0000LL):(r2&0x0000000000ff0000LL))
       | ((r1&0x000000000000ff00LL)>(r2&0x000000000000ff00LL)?(r1&0x000000000000ff00LL):(r2&0x000000000000ff00LL))
       | ((r1&0x00000000000000ffLL)>(r2&0x00000000000000ffLL)?(r1&0x00000000000000ffLL):(r2&0x00000000000000ffLL));
    ex1_outd = ((r3&0xff00000000000000LL)<(t1&0xff00000000000000LL)?(t1&0xff00000000000000LL):((r3&0xff00000000000000LL)<(t2&0xff00000000000000LL)?(r3&0xff00000000000000LL):(t2&0xff00000000000000LL)))
             | ((r3&0x00ff000000000000LL)<(t1&0x00ff000000000000LL)?(t1&0x00ff000000000000LL):((r3&0x00ff000000000000LL)<(t2&0x00ff000000000000LL)?(r3&0x00ff000000000000LL):(t2&0x00ff000000000000LL)))
             | ((r3&0x0000ff0000000000LL)<(t1&0x0000ff0000000000LL)?(t1&0x0000ff0000000000LL):((r3&0x0000ff0000000000LL)<(t2&0x0000ff0000000000LL)?(r3&0x0000ff0000000000LL):(t2&0x0000ff0000000000LL)))
             | ((r3&0x000000ff00000000LL)<(t1&0x000000ff00000000LL)?(t1&0x000000ff00000000LL):((r3&0x000000ff00000000LL)<(t2&0x000000ff00000000LL)?(r3&0x000000ff00000000LL):(t2&0x000000ff00000000LL)))
             | ((r3&0x00000000ff000000LL)<(t1&0x00000000ff000000LL)?(t1&0x00000000ff000000LL):((r3&0x00000000ff000000LL)<(t2&0x00000000ff000000LL)?(r3&0x00000000ff000000LL):(t2&0x00000000ff000000LL)))
             | ((r3&0x0000000000ff0000LL)<(t1&0x0000000000ff0000LL)?(t1&0x0000000000ff0000LL):((r3&0x0000000000ff0000LL)<(t2&0x0000000000ff0000LL)?(r3&0x0000000000ff0000LL):(t2&0x0000000000ff0000LL)))
             | ((r3&0x000000000000ff00LL)<(t1&0x000000000000ff00LL)?(t1&0x000000000000ff00LL):((r3&0x000000000000ff00LL)<(t2&0x000000000000ff00LL)?(r3&0x000000000000ff00LL):(t2&0x000000000000ff00LL)))
             | ((r3&0x00000000000000ffLL)<(t1&0x00000000000000ffLL)?(t1&0x00000000000000ffLL):((r3&0x00000000000000ffLL)<(t2&0x00000000000000ffLL)?(r3&0x00000000000000ffLL):(t2&0x00000000000000ffLL)));
    break;
  case OP_MMAX3: /* 3in 8bit*8 bytewise compare and output maximum */
    t1 = ((r1&0xff00000000000000LL)>(r2&0xff00000000000000LL)?(r1&0xff00000000000000LL):(r2&0xff00000000000000LL))
       | ((r1&0x00ff000000000000LL)>(r2&0x00ff000000000000LL)?(r1&0x00ff000000000000LL):(r2&0x00ff000000000000LL))
       | ((r1&0x0000ff0000000000LL)>(r2&0x0000ff0000000000LL)?(r1&0x0000ff0000000000LL):(r2&0x0000ff0000000000LL))
       | ((r1&0x000000ff00000000LL)>(r2&0x000000ff00000000LL)?(r1&0x000000ff00000000LL):(r2&0x000000ff00000000LL))
       | ((r1&0x00000000ff000000LL)>(r2&0x00000000ff000000LL)?(r1&0x00000000ff000000LL):(r2&0x00000000ff000000LL))
       | ((r1&0x0000000000ff0000LL)>(r2&0x0000000000ff0000LL)?(r1&0x0000000000ff0000LL):(r2&0x0000000000ff0000LL))
       | ((r1&0x000000000000ff00LL)>(r2&0x000000000000ff00LL)?(r1&0x000000000000ff00LL):(r2&0x000000000000ff00LL))
       | ((r1&0x00000000000000ffLL)>(r2&0x00000000000000ffLL)?(r1&0x00000000000000ffLL):(r2&0x00000000000000ffLL));
    ex1_outd = ((t1&0xff00000000000000LL)>(r3&0xff00000000000000LL)?(t1&0xff00000000000000LL):(r3&0xff00000000000000LL))
             | ((t1&0x00ff000000000000LL)>(r3&0x00ff000000000000LL)?(t1&0x00ff000000000000LL):(r3&0x00ff000000000000LL))
             | ((t1&0x0000ff0000000000LL)>(r3&0x0000ff0000000000LL)?(t1&0x0000ff0000000000LL):(r3&0x0000ff0000000000LL))
             | ((t1&0x000000ff00000000LL)>(r3&0x000000ff00000000LL)?(t1&0x000000ff00000000LL):(r3&0x000000ff00000000LL))
             | ((t1&0x00000000ff000000LL)>(r3&0x00000000ff000000LL)?(t1&0x00000000ff000000LL):(r3&0x00000000ff000000LL))
             | ((t1&0x0000000000ff0000LL)>(r3&0x0000000000ff0000LL)?(t1&0x0000000000ff0000LL):(r3&0x0000000000ff0000LL))
             | ((t1&0x000000000000ff00LL)>(r3&0x000000000000ff00LL)?(t1&0x000000000000ff00LL):(r3&0x000000000000ff00LL))
             | ((t1&0x00000000000000ffLL)>(r3&0x00000000000000ffLL)?(t1&0x00000000000000ffLL):(r3&0x00000000000000ffLL));
    break;
  case OP_MMIN3: /* 3in 8bit*8 bytewise compare and output minimum */
    t1 = ((r1&0xff00000000000000LL)<(r2&0xff00000000000000LL)?(r1&0xff00000000000000LL):(r2&0xff00000000000000LL))
       | ((r1&0x00ff000000000000LL)<(r2&0x00ff000000000000LL)?(r1&0x00ff000000000000LL):(r2&0x00ff000000000000LL))
       | ((r1&0x0000ff0000000000LL)<(r2&0x0000ff0000000000LL)?(r1&0x0000ff0000000000LL):(r2&0x0000ff0000000000LL))
       | ((r1&0x000000ff00000000LL)<(r2&0x000000ff00000000LL)?(r1&0x000000ff00000000LL):(r2&0x000000ff00000000LL))
       | ((r1&0x00000000ff000000LL)<(r2&0x00000000ff000000LL)?(r1&0x00000000ff000000LL):(r2&0x00000000ff000000LL))
       | ((r1&0x0000000000ff0000LL)<(r2&0x0000000000ff0000LL)?(r1&0x0000000000ff0000LL):(r2&0x0000000000ff0000LL))
       | ((r1&0x000000000000ff00LL)<(r2&0x000000000000ff00LL)?(r1&0x000000000000ff00LL):(r2&0x000000000000ff00LL))
       | ((r1&0x00000000000000ffLL)<(r2&0x00000000000000ffLL)?(r1&0x00000000000000ffLL):(r2&0x00000000000000ffLL));
    ex1_outd = ((t1&0xff00000000000000LL)<(r3&0xff00000000000000LL)?(t1&0xff00000000000000LL):(r3&0xff00000000000000LL))
             | ((t1&0x00ff000000000000LL)<(r3&0x00ff000000000000LL)?(t1&0x00ff000000000000LL):(r3&0x00ff000000000000LL))
             | ((t1&0x0000ff0000000000LL)<(r3&0x0000ff0000000000LL)?(t1&0x0000ff0000000000LL):(r3&0x0000ff0000000000LL))
             | ((t1&0x000000ff00000000LL)<(r3&0x000000ff00000000LL)?(t1&0x000000ff00000000LL):(r3&0x000000ff00000000LL))
             | ((t1&0x00000000ff000000LL)<(r3&0x00000000ff000000LL)?(t1&0x00000000ff000000LL):(r3&0x00000000ff000000LL))
             | ((t1&0x0000000000ff0000LL)<(r3&0x0000000000ff0000LL)?(t1&0x0000000000ff0000LL):(r3&0x0000000000ff0000LL))
             | ((t1&0x000000000000ff00LL)<(r3&0x000000000000ff00LL)?(t1&0x000000000000ff00LL):(r3&0x000000000000ff00LL))
             | ((t1&0x00000000000000ffLL)<(r3&0x00000000000000ffLL)?(t1&0x00000000000000ffLL):(r3&0x00000000000000ffLL));
    break;
  case OP_MMAX: /* 2in 8bit*8 bytewise compare and output maximum */
    ex1_outd = ((r1&0xff00000000000000LL)>(r2&0xff00000000000000LL)?(r1&0xff00000000000000LL):(r2&0xff00000000000000LL))
             | ((r1&0x00ff000000000000LL)>(r2&0x00ff000000000000LL)?(r1&0x00ff000000000000LL):(r2&0x00ff000000000000LL))
             | ((r1&0x0000ff0000000000LL)>(r2&0x0000ff0000000000LL)?(r1&0x0000ff0000000000LL):(r2&0x0000ff0000000000LL))
             | ((r1&0x000000ff00000000LL)>(r2&0x000000ff00000000LL)?(r1&0x000000ff00000000LL):(r2&0x000000ff00000000LL))
             | ((r1&0x00000000ff000000LL)>(r2&0x00000000ff000000LL)?(r1&0x00000000ff000000LL):(r2&0x00000000ff000000LL))
             | ((r1&0x0000000000ff0000LL)>(r2&0x0000000000ff0000LL)?(r1&0x0000000000ff0000LL):(r2&0x0000000000ff0000LL))
             | ((r1&0x000000000000ff00LL)>(r2&0x000000000000ff00LL)?(r1&0x000000000000ff00LL):(r2&0x000000000000ff00LL))
             | ((r1&0x00000000000000ffLL)>(r2&0x00000000000000ffLL)?(r1&0x00000000000000ffLL):(r2&0x00000000000000ffLL));
    break;
  case OP_MMIN: /* 2in 8bit*8 bytewise compare and output minimum */
    ex1_outd = ((r1&0xff00000000000000LL)<(r2&0xff00000000000000LL)?(r1&0xff00000000000000LL):(r2&0xff00000000000000LL))
             | ((r1&0x00ff000000000000LL)<(r2&0x00ff000000000000LL)?(r1&0x00ff000000000000LL):(r2&0x00ff000000000000LL))
             | ((r1&0x0000ff0000000000LL)<(r2&0x0000ff0000000000LL)?(r1&0x0000ff0000000000LL):(r2&0x0000ff0000000000LL))
             | ((r1&0x000000ff00000000LL)<(r2&0x000000ff00000000LL)?(r1&0x000000ff00000000LL):(r2&0x000000ff00000000LL))
             | ((r1&0x00000000ff000000LL)<(r2&0x00000000ff000000LL)?(r1&0x00000000ff000000LL):(r2&0x00000000ff000000LL))
             | ((r1&0x0000000000ff0000LL)<(r2&0x0000000000ff0000LL)?(r1&0x0000000000ff0000LL):(r2&0x0000000000ff0000LL))
             | ((r1&0x000000000000ff00LL)<(r2&0x000000000000ff00LL)?(r1&0x000000000000ff00LL):(r2&0x000000000000ff00LL))
             | ((r1&0x00000000000000ffLL)<(r2&0x00000000000000ffLL)?(r1&0x00000000000000ffLL):(r2&0x00000000000000ffLL));
    break;
  case OP_MAJ: /* (((x) & (y))^((x) & (z))^((y) & (z))) */
    ex1_outd = (r1&0xffffffff00000000LL) | (((r1 & r2)^(r1 & r3)^(r2 & r3))&0xffffffffLL);
    break;
  case OP_CH: /*  (((x) & (y))^(~(x) & (z))) */
    ex1_outd = (r1&0xffffffff00000000LL) | (((r1 & r2)^(~r1 & r3))&0xffffffffLL);
    break;
  default:
    printf("emax7lib: exe: undefined op_ex1=%d\n", op_ex1);
    break;
  }

  switch (op_ex2) {
  case OP_NOP:
    if (op_ex1 == OP_SFMA)
      softu64(2, ex1_outd_sfma, &ex2_outd, NULL, r1, r2, r3, r4);
    else
      ex2_outd = ex1_outd;
    break;
  case OP_AND: /* 2in 64bit logical and s1&s2 */
    ex2_outd = ex1_outd & r4;
    break;
  case OP_OR: /* 2in 64bit logical or s1|s2 */
    ex2_outd = ex1_outd | r4;
    break;
  case OP_XOR: /* 2in 64bit logical xor s1^s2 */
    ex2_outd = ex1_outd ^ r4;
    break;
  case OP_SUMHH: /* 1in 16bit*4 & s1.h3+s1.h2->d.h3, s1.h1+s1.h0->d.h1 */
    t3 = ex1_outd>>48&0x000000000000ffffLL;
    t2 = ex1_outd>>32&0x000000000000ffffLL;
    t1 = ex1_outd>>16&0x000000000000ffffLL;
    t0 = ex1_outd    &0x000000000000ffffLL;
    t3 += t2;
    if (t3 > 0x000000000000ffffLL) t3 = 0x000000000000ffffLL;
    t1 += t0;
    if (t1 > 0x000000000000ffffLL) t1 = 0x000000000000ffffLL;
    ex2_outd = (t3<<48)|(t1<<16);
    break;
  case OP_SUMHL: /* 1in 16bit*4 & s1.h3+s1.h2->d.h2, s1.h1+s1.h0->d.h0 */
    t3 = ex1_outd>>48&0x000000000000ffffLL;
    t2 = ex1_outd>>32&0x000000000000ffffLL;
    t1 = ex1_outd>>16&0x000000000000ffffLL;
    t0 = ex1_outd    &0x000000000000ffffLL;
    t2 += t3;
    if (t2 > 0x000000000000ffffLL) t2 = 0x000000000000ffffLL;
    t0 += t1;
    if (t0 > 0x000000000000ffffLL) t0 = 0x000000000000ffffLL;
    ex2_outd = (t2<<32)|(t0);
    break;
//case OP_WSWAP: /* 32bit 2in swap and mask words */
//  ex2_outd = ((ex1_outd<<32)|(ex1_outd>>32)) & r4;
//  break;
  case OP_ROTS: /* hi-32bit #define ROTRIGHT (((a) >> (b)) | ((a) << (32-(b)))) */
    t2 = ex1_outd & 0xffffffff00000000LL;
    ro10 = r4>>32 & 0xff;
    ro11 = r4>>40 & 0xff;
    ro12 = r4>>48 & 0xff;
    t0 = ex1_outd & 0x00000000ffffffffLL;
    ro00 = r4     & 0xff;
    ro01 = r4>> 8 & 0xff;
    ro02 = r4>>16 & 0xff;
    ex2_outd = (((t2>>ro12|t2<<(32-ro12))^(t2>>ro11|t2<<(32-ro11))^(t2>>ro10|t2<<(32-ro10)))&0xffffffff00000000LL)
              |(((t0>>ro02|t0<<(32-ro02))^(t0>>ro01|t0<<(32-ro01))^(t0>>ro00|t0<<(32-ro00)))&0x00000000ffffffffLL);
    break;
  default:
    printf("emax7lib: exe: undefined op_ex2=%d\n", op_ex2);
    break;
  }

  switch (op_ex3) {
  case OP_NOP:
    if (op_ex1 == OP_SFMA)
      softu64(3, NULL, &ex2_outd, d, r1, r2, r3, r4);
    else
      if (d) *d = ex2_outd;
    break;
  case OP_SLL: /* 2in 32bit*2 32bit logical shift to left */
    t1 = (Ull)(ex2_outd    &0xffffffff00000000LL)<<r5;
    t0 = (Ull)(ex2_outd<<r5&0x00000000ffffffffLL);
    if (d) *d = t1 | t0;
    break;
  case OP_SRL: /* 2in 32bit*2 32bit logical shift to right */
    t1 = (Ull)(ex2_outd>>r5&0xffffffff00000000LL);
    t0 = (Ull)(ex2_outd    &0x00000000ffffffffLL)>>r5;
    if (d) *d = t1 | t0;
    break;
  case OP_SRAA: /* 2in 32bit*2 32bit arith shift to right (bit63,31 is ext.) */
    t1 = (Sll)(ex2_outd    )>>r5&0xffffffff00000000LL;
    t0 = (Sll)(ex2_outd<<32)>>r5&0xffffffff00000000LL;
    if (d) *d = t1 | (t0>>32);
    break;
  case OP_SRAB: /* 2in 32bit*2 32bit arith shift to right (bit55,23 is ext.) */
    t1 = (Sll)(ex2_outd<< 8)>>(r5+8)&0xffffffff00000000LL;
    t0 = (Sll)(ex2_outd<<40)>>(r5+8)&0xffffffff00000000LL;
    if (d) *d = t1 | (t0>>32);
    break;
//case OP_SRAC: /* 2in 32bit*2 32bit arith shift to right (bit47,15 is ext.) */
//  t1 = (Sll)(ex2_outd<<16)>>(r5+16)&0xffffffff00000000LL;
//  t0 = (Sll)(ex2_outd<<48)>>(r5+16)&0xffffffff00000000LL;
//  if (d) *d = t1 | (t0>>32);
//  break;
//case OP_SRAD: /* 2in 32bit*2 32bit arith shift to right (bit39,7 is ext.) */
//  t1 = (Sll)(ex2_outd<<24)>>(r5+24)&0xffffffff00000000LL;
//  t0 = (Sll)(ex2_outd<<56)>>(r5+24)&0xffffffff00000000LL;
//  if (d) *d = t1 | (t0>>32);
//  break;
  case OP_SRLM: /* 2in 16bit*4 16bit arith shift to right */
    t3 = (Ull)(ex2_outd    )>>r5&0xffff000000000000LL;
    t2 = (Ull)(ex2_outd<<16)>>r5&0xffff000000000000LL;
    t1 = (Ull)(ex2_outd<<32)>>r5&0xffff000000000000LL;
    t0 = (Ull)(ex2_outd<<48)>>r5&0xffff000000000000LL;
    if (d) *d = t3 | (t2>>16) | (t1>>32) | (t0>>48);
    break;
  default:
    printf("emax7lib: exe: undefined op_ex3=%d\n", op_ex3);
    break;
  }

  return (retval);
}

void /*__attribute__((always_inline))*/
mex(Uint op_mex2, Uchar **d2, Uchar *base2, Ull ofs2, Uint op_mex1, Uchar **d1, Uchar *base1, Ull ofs1, Ull limit, Ull s2, Ull s1)
{
  /* limit:  0, 8, 16, .... 4096, 8192, 16384, 32768     */
  /* encode: 0, 1, 2,  3,   10    11    12     13 (4bit) */
  Uint limit2 = limit*2;
  Uint ss2 = s2>>32;
  Uint ss1 = s1>>32;

  switch (op_mex1) {
  case OP_NOP:
    *d1 = base1;
    break;
  case OP_ALWAYS: /* base++ �б� */
    *d1 = base1 + ofs1;
    break;
  case OP_CMPA_GE:
    //d1��ʬ(ss1)��ffffffff�ʤ����. base1==limit2�ʤ����. base2==limit �ʤ�����
    if (!limit) /* sparse matrix */
      *d1 = base1 + ((ss1!=0xffffffff && ss2>=ss1) ? ofs1:0);
    else { /* merge sort */
      if ((base2==limit && base1+ofs1==limit2)||(base2+ofs2==limit && base1==limit2))
	*d1 = limit;
      else
	*d1 = base1 + (base1!=limit2 && ((base2!=limit  && ss2>=ss1)||base2==limit ) ? ofs1:0);
    }
    break;
  default:
    printf("emax7lib: mex: undefined op_mex1=%d\n", op_mex1);
    break;
  }  

  switch (op_mex2) {
  case OP_NOP:
    *d2 = base2;
    break;
  case OP_ALWAYS: /* base++ �б� */
    *d2 = base2 + ofs2;
    break;
  case OP_CMPA_LE:
    //d2��ʬ(ss2)��ffffffff�ʤ����. base2==limit �ʤ����. base1==limit2�ʤ�����
    if (!limit) /* sparse matrix */
      *d2 = base2 + ((ss2!=0xffffffff && ss2<=ss1) ? ofs2:0);
    else { /* merge sort */
      if ((base2==limit && base1+ofs1==limit2)||(base2+ofs2==limit && base1==limit2))
	*d2 = 0;
      else
	*d2 = base2 + (base2!=limit  && ((base1!=limit2 && ss2<=ss1)||base1==limit2) ? ofs2:0);
    }
    break;
  default:
    printf("emax7lib: mex: undefined op_mex2=%d\n", op_mex2);
    break;
  }  
}

void /*__attribute__((always_inline))*/
eag(Ull *adr, Ull base, Ull ofs)
{
  *adr = base + ofs;
}

void /*__attribute__((always_inline))*/
mop(Uint op_mm, Ull ex, Ull *d, Ull base, Ull offset, Uchar msk, Ull top, Uint len, Uint blk, Uchar force, Ull ptop, Uint plen)
{
  Ull adr, ofs;

  eag(&adr, base, eam(offset, msk));
  mmp(op_mm, ex, d, adr, top, len, blk);
}

void /*__attribute__((always_inline))*/
mo4(Uint op_mm, Ull ex, Ull *d, Ull base, Ull offset, Uchar msk, Ull top, Uint len, Uint blk, Uchar force, Ull ptop, Uint plen)
{
  Ull adr, ofs;

  eag(&adr, base, eam(offset, msk));
  mmp(op_mm, ex, d, adr, top, len, blk);
}

void /*__attribute__((always_inline))*/
mmp(Uint op_mm, Ull ex, Ull *d, Ull adr, Ull top, Uint len, Uint blk)
{
  Ull c1, c0, load64;

#if defined(__i386)
  adr &= (1LL<<32)-1;
  top &= (1LL<<32)-1;
#endif  

  if (!((op_mm==OP_LDRQ && blk) || op_mm==OP_LDDMQ || op_mm==OP_TR) && (!adr || !top)) return; /* NULL skip DMA */

#define CHECK_MMP_MARGIN 12
  if (!((op_mm==OP_LDRQ && blk) || op_mm==OP_LDDMQ || op_mm==OP_TR) && ex && len && (adr < top || adr >= top+len*sizeof(Uint)+CHECK_MMP_MARGIN)) {
    printf("mmp: adr=%08.8x_%08.8x out of range (top=%08.8x_%08.8x len=%dB)\n", (Uint)(adr>>32), (Uint)adr, (Uint)(top>>32), (Uint)top, len*sizeof(Uint));
    fflush(stdout);
  }

  c1 = ex>>1&1;
  c0 = ex   &1;

  switch (op_mm) {
  case OP_NOP:
    break;

    /* MOP */
  case OP_LDR: /* 64bit lmm LMM is preloaded, random-access */
    load64 = *(Ull*)(adr&~7LL);
    if ((adr&7) == 0)
      *d = load64;
    else if (!emax7_unaligned_load_valid) { /* BR[][][1] */
      emax7_unaligned_load_valid = 1;
      emax7_unaligned_load_high = load64;
      *d = load64 >> (adr&7)*8;
    }
    else { /* BR[][][0] */
      emax7_unaligned_load_valid = 0; 
      *d = emax7_unaligned_load_high << (8-(adr&7))*8 | load64 >> (adr&7)*8;
    }
    break;
  case OP_LDWR: /* u32bit lmm LMM is preloaded, random-access */
    *d = (Ull)*(Uint*)(adr&~3LL);
    break;
//case OP_LDHR: /* u16bit lmm LMM is preloaded, random-access */
//  *d = (Ull)(Uint)*(Ushort*)(adr&~1LL);
//  break;
  case OP_LDBR: /* u8bit lmm LMM is preloaded, random-access */
    *d = (Ull)(Uint)*(Uchar*)adr;
    break;
  case OP_STR: /* 64bit lmm LMM is drained. random-access */
    if (c1) *((Uint*)(adr&~7LL)+1) = *d>>32;
    if (c0) *((Uint*)(adr&~7LL)  ) = *d;
    break;
  case OP_STWR: /* 32bit lmm LMM is drained. random-access */
    if (c0) *(Uint*)(adr&~3LL) = *d;
    break;
//case OP_STHR: /* 16bit lmm LMM is drained. random-access */
//  if (c0) *(Ushort*)(adr&~1LL) = *d;
//  break;
  case OP_STBR: /* 8bit lmm LMM is drained. random-access */
    if (c0) *(Uchar*)adr = *d;
    break;

    /* MO4 */
  case OP_LDRQ: /* 64bit*4 lmm LMM is preloaded, random-access */
    switch (blk) {
    case 0: /* normal */
      /* adr=0,32,64,... */
      *(d+0) = *((Ull*)(adr&~31LL)+0);
      *(d+1) = *((Ull*)(adr&~31LL)+1);
      *(d+2) = *((Ull*)(adr&~31LL)+2);
      *(d+3) = *((Ull*)(adr&~31LL)+3);
      break;
    case 1: /* block_size=16-members */
      /* adr=0,32,64,... memadr = mem(top + (adr/32/16*ptr)) + (adr/32/&15)*4 */
      *(d+0) = *(*(Ull**)(top + (adr/32/16*sizeof(Ull*))) + (adr/32&15)*4 + 0);
      *(d+1) = *(*(Ull**)(top + (adr/32/16*sizeof(Ull*))) + (adr/32&15)*4 + 1);
      *(d+2) = *(*(Ull**)(top + (adr/32/16*sizeof(Ull*))) + (adr/32&15)*4 + 2);
      *(d+3) = *(*(Ull**)(top + (adr/32/16*sizeof(Ull*))) + (adr/32&15)*4 + 3);
      break;
    case 2: /* block_size=32-members */
      /* adr=0,32,64,... memadr = mem(top + (adr/32/32*ptr)) + (adr/32/&31)*4 */
      *(d+0) = *(*(Ull**)(top + (adr/32/32*sizeof(Ull*))) + (adr/32&31)*4 + 0);
      *(d+1) = *(*(Ull**)(top + (adr/32/32*sizeof(Ull*))) + (adr/32&31)*4 + 1);
      *(d+2) = *(*(Ull**)(top + (adr/32/32*sizeof(Ull*))) + (adr/32&31)*4 + 2);
      *(d+3) = *(*(Ull**)(top + (adr/32/32*sizeof(Ull*))) + (adr/32&31)*4 + 3);
      break;
    default:/* block_size=64-members */
      /* adr=0,32,64,... memadr = mem(top + (adr/32/64*ptr)) + (adr/32/&63)*4 */
      *(d+0) = *(*(Ull**)(top + (adr/32/64*sizeof(Ull*))) + (adr/32&63)*4 + 0);
      *(d+1) = *(*(Ull**)(top + (adr/32/64*sizeof(Ull*))) + (adr/32&63)*4 + 1);
      *(d+2) = *(*(Ull**)(top + (adr/32/64*sizeof(Ull*))) + (adr/32&63)*4 + 2);
      *(d+3) = *(*(Ull**)(top + (adr/32/64*sizeof(Ull*))) + (adr/32&63)*4 + 3);
      break;
    }
    break;
  case OP_LDDMQ: /* 64bit*4 mem Direct access to MM */
    if (c0) {
      *(d+0) = *((Ull*)(adr&~31LL)+0);
      *(d+1) = *((Ull*)(adr&~31LL)+1);
      *(d+2) = *((Ull*)(adr&~31LL)+2);
      *(d+3) = *((Ull*)(adr&~31LL)+3);
    }
    break;
  case OP_STRQ: /* 64bit*4 lmm LMM is drained. random-access */
    *((Ull*)(adr&~31LL)+0) = *(d+0);
    *((Ull*)(adr&~31LL)+1) = *(d+1);
    *((Ull*)(adr&~31LL)+2) = *(d+2);
    *((Ull*)(adr&~31LL)+3) = *(d+3);
    break;
  case OP_TR: /* 64bit*4 exec Send transaction */
    /* addr��transaction()����˻��� */
    if (c0) {
      Ull (*trans)() = top;
      (trans)(*(d+0), *(d+1), *(d+2), *(d+3));
    }
    break;
  default:
    printf("emax7lib: mmp: undefined op_mm=%d\n", op_mm);
    break;
  }
}
#endif