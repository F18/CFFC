/*
 * Copyright (c) 2003, 2007 Matteo Frigo
 * Copyright (c) 2003, 2007 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Mon Mar 19 08:01:47 CDT 2007 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2c -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -twiddle-log3 -precompute-twiddles -n 16 -dif -name hc2cb2_16 -include hc2cb.h */

/*
 * This function contains 196 FP additions, 134 FP multiplications,
 * (or, 104 additions, 42 multiplications, 92 fused multiply/add),
 * 112 stack variables, and 64 memory accesses
 */
#include "hc2cb.h"

static void hc2cb2_16(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     INT m;
     for (m = mb, W = W + ((mb - 1) * 8); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 8, MAKE_VOLATILE_STRIDE(rs)) {
	  E Tv, TB, TF, Ty, T1J, T1O, T1N, T1K;
	  {
	       E Tw, T2z, T2C, Tx, T3f, T3l, T2F, T3r, Tz;
	       Tv = W[0];
	       Tw = W[2];
	       T2z = W[6];
	       T2C = W[7];
	       TB = W[4];
	       Tx = Tv * Tw;
	       T3f = Tv * T2z;
	       T3l = Tv * T2C;
	       T2F = Tv * TB;
	       T3r = Tw * TB;
	       TF = W[5];
	       Ty = W[1];
	       Tz = W[3];
	       {
		    E T2G, T3z, T3m, T3g, T3L, T3s, T1V, TA, T3w, T3Q, T30, T3C, TE, T1X, T1D;
		    E TG, T1G, T1o, T2p, T1Y, T2u, T2c, T1Z, TL, T1t, T2d, T3n, T35, T3R, T3F;
		    E T1w, T20, T3M, Tf, T3h, T2L, T2e, TW, T2Q, T36, T3I, T3N, T2V, T37, T1d;
		    E Tu, T3S, T18, T1z, T1i, T24, T2g, T27, T2h;
		    {
			 E T2K, TQ, TV, T2H;
			 {
			      E TH, T3, T32, T1s, T1p, T6, T33, TK, TM, Ta, TS, T2J, TP, TR, Td;
			      E TT, TI, TJ;
			      {
				   E T1q, T1r, T4, T5;
				   {
					E T1, T1n, TC, T2b, T1W, T2, T3v, T2Z, TD;
					T1 = Rp[0];
					T3v = Tw * TF;
					T2Z = Tv * TF;
					T2G = FNMS(Ty, TF, T2F);
					T3z = FMA(Ty, TF, T2F);
					T3m = FNMS(Ty, T2z, T3l);
					T3g = FMA(Ty, T2C, T3f);
					T3L = FNMS(Tz, TF, T3r);
					T3s = FMA(Tz, TF, T3r);
					T1V = FMA(Ty, Tz, Tx);
					TA = FNMS(Ty, Tz, Tx);
					TD = Tv * Tz;
					T3w = FNMS(Tz, TB, T3v);
					T3Q = FMA(Tz, TB, T3v);
					T30 = FMA(Ty, TB, T2Z);
					T3C = FNMS(Ty, TB, T2Z);
					T1n = TA * TF;
					TC = TA * TB;
					T2b = T1V * TF;
					T1W = T1V * TB;
					TE = FMA(Ty, Tw, TD);
					T1X = FNMS(Ty, Tw, TD);
					T2 = Rm[WS(rs, 7)];
					T1q = Ip[0];
					T1D = FMA(TE, TF, TC);
					TG = FNMS(TE, TF, TC);
					T1G = FNMS(TE, TB, T1n);
					T1o = FMA(TE, TB, T1n);
					T2p = FMA(T1X, TF, T1W);
					T1Y = FNMS(T1X, TF, T1W);
					T2u = FNMS(T1X, TB, T2b);
					T2c = FMA(T1X, TB, T2b);
					TH = T1 - T2;
					T3 = T1 + T2;
					T1r = Im[WS(rs, 7)];
				   }
				   T4 = Rp[WS(rs, 4)];
				   T5 = Rm[WS(rs, 3)];
				   TI = Ip[WS(rs, 4)];
				   T32 = T1q - T1r;
				   T1s = T1q + T1r;
				   T1p = T4 - T5;
				   T6 = T4 + T5;
				   TJ = Im[WS(rs, 3)];
			      }
			      {
				   E TN, TO, T8, T9, Tb, Tc;
				   T8 = Rp[WS(rs, 2)];
				   T9 = Rm[WS(rs, 5)];
				   TN = Ip[WS(rs, 2)];
				   T33 = TI - TJ;
				   TK = TI + TJ;
				   TM = T8 - T9;
				   Ta = T8 + T9;
				   TO = Im[WS(rs, 5)];
				   Tb = Rm[WS(rs, 1)];
				   Tc = Rp[WS(rs, 6)];
				   TS = Ip[WS(rs, 6)];
				   T2J = TN - TO;
				   TP = TN + TO;
				   TR = Tb - Tc;
				   Td = Tb + Tc;
				   TT = Im[WS(rs, 1)];
			      }
			      {
				   E T2I, TU, Te, T31, T34, T3D;
				   T1Z = TH + TK;
				   TL = TH - TK;
				   T1t = T1p + T1s;
				   T2d = T1s - T1p;
				   T2I = TS - TT;
				   TU = TS + TT;
				   Te = Ta + Td;
				   T31 = Ta - Td;
				   T34 = T32 - T33;
				   T3D = T32 + T33;
				   {
					E T1u, T1v, T3E, T7;
					T3E = T2J + T2I;
					T2K = T2I - T2J;
					TQ = TM - TP;
					T1u = TM + TP;
					T3n = T34 - T31;
					T35 = T31 + T34;
					TV = TR - TU;
					T1v = TR + TU;
					T3R = T3D - T3E;
					T3F = T3D + T3E;
					T2H = T3 - T6;
					T7 = T3 + T6;
					T1w = T1u - T1v;
					T20 = T1u + T1v;
					T3M = T7 - Te;
					Tf = T7 + Te;
				   }
			      }
			 }
			 {
			      E T1e, Ti, T2N, T1c, T19, Tl, T2O, T1h, Tq, T13, Tp, T2S, T11, Tr, T14;
			      E T15;
			      {
				   E Tj, Tk, T1f, T1g;
				   {
					E Tg, Th, T1a, T1b;
					Tg = Rp[WS(rs, 1)];
					T3h = T2H - T2K;
					T2L = T2H + T2K;
					T2e = TQ - TV;
					TW = TQ + TV;
					Th = Rm[WS(rs, 6)];
					T1a = Ip[WS(rs, 1)];
					T1b = Im[WS(rs, 6)];
					Tj = Rp[WS(rs, 5)];
					T1e = Tg - Th;
					Ti = Tg + Th;
					T2N = T1a - T1b;
					T1c = T1a + T1b;
					Tk = Rm[WS(rs, 2)];
					T1f = Ip[WS(rs, 5)];
					T1g = Im[WS(rs, 2)];
				   }
				   {
					E Tn, To, TZ, T10;
					Tn = Rm[0];
					T19 = Tj - Tk;
					Tl = Tj + Tk;
					T2O = T1f - T1g;
					T1h = T1f + T1g;
					To = Rp[WS(rs, 7)];
					TZ = Ip[WS(rs, 7)];
					T10 = Im[0];
					Tq = Rp[WS(rs, 3)];
					T13 = Tn - To;
					Tp = Tn + To;
					T2S = TZ - T10;
					T11 = TZ + T10;
					Tr = Rm[WS(rs, 4)];
					T14 = Ip[WS(rs, 3)];
					T15 = Im[WS(rs, 4)];
				   }
			      }
			      {
				   E TY, T16, Tm, Tt;
				   {
					E T2P, T3G, Ts, T2M, T3H, T2U, T2T, T2R;
					T2P = T2N - T2O;
					T3G = T2N + T2O;
					TY = Tq - Tr;
					Ts = Tq + Tr;
					T2T = T14 - T15;
					T16 = T14 + T15;
					T2M = Ti - Tl;
					Tm = Ti + Tl;
					T3H = T2S + T2T;
					T2U = T2S - T2T;
					Tt = Tp + Ts;
					T2R = Tp - Ts;
					T2Q = T2M - T2P;
					T36 = T2M + T2P;
					T3I = T3G + T3H;
					T3N = T3H - T3G;
					T2V = T2R + T2U;
					T37 = T2U - T2R;
				   }
				   {
					E T25, T26, T22, T23, T12, T17;
					T12 = TY - T11;
					T25 = TY + T11;
					T26 = T13 + T16;
					T17 = T13 - T16;
					T22 = T1c - T19;
					T1d = T19 + T1c;
					Tu = Tm + Tt;
					T3S = Tm - Tt;
					T18 = FNMS(KP414213562, T17, T12);
					T1z = FMA(KP414213562, T12, T17);
					T1i = T1e - T1h;
					T23 = T1e + T1h;
					T24 = FNMS(KP414213562, T23, T22);
					T2g = FMA(KP414213562, T22, T23);
					T27 = FNMS(KP414213562, T26, T25);
					T2h = FMA(KP414213562, T25, T26);
				   }
			      }
			 }
		    }
		    {
			 E T1j, T1y, T3V, T3X, T3W, T38, T3i, T3o, T2W, T3K, T3B, T3A;
			 Rp[0] = Tf + Tu;
			 T3A = Tf - Tu;
			 T1j = FMA(KP414213562, T1i, T1d);
			 T1y = FNMS(KP414213562, T1d, T1i);
			 T3K = T3C * T3A;
			 T3B = T3z * T3A;
			 {
			      E T3O, T3T, T3J, T3P, T3U;
			      T3O = T3M - T3N;
			      T3V = T3M + T3N;
			      T3X = T3S + T3R;
			      T3T = T3R - T3S;
			      Rm[0] = T3F + T3I;
			      T3J = T3F - T3I;
			      T3P = T3L * T3O;
			      T3U = T3L * T3T;
			      T3W = TA * T3V;
			      Rp[WS(rs, 4)] = FNMS(T3C, T3J, T3B);
			      Rm[WS(rs, 4)] = FMA(T3z, T3J, T3K);
			      Rp[WS(rs, 6)] = FNMS(T3Q, T3T, T3P);
			      Rm[WS(rs, 6)] = FMA(T3Q, T3O, T3U);
			      T38 = T36 + T37;
			      T3i = T37 - T36;
			      T3o = T2Q - T2V;
			      T2W = T2Q + T2V;
			 }
			 {
			      E T2q, T21, T28, T2w, T2v, T2f, T2i, T2r;
			      {
				   E T2Y, T3a, T3c, T3d, T39, T3e, T3b, T2X, T3Y;
				   Rp[WS(rs, 2)] = FNMS(TE, T3X, T3W);
				   T3Y = TA * T3X;
				   {
					E T3t, T3j, T3x, T3p;
					T3t = FMA(KP707106781, T3i, T3h);
					T3j = FNMS(KP707106781, T3i, T3h);
					T3x = FMA(KP707106781, T3o, T3n);
					T3p = FNMS(KP707106781, T3o, T3n);
					Rm[WS(rs, 2)] = FMA(TE, T3V, T3Y);
					{
					     E T3u, T3k, T3y, T3q;
					     T3u = T3s * T3t;
					     T3k = T3g * T3j;
					     T3y = T3s * T3x;
					     T3q = T3g * T3p;
					     Rp[WS(rs, 3)] = FNMS(T3w, T3x, T3u);
					     Rp[WS(rs, 7)] = FNMS(T3m, T3p, T3k);
					     Rm[WS(rs, 3)] = FMA(T3w, T3t, T3y);
					     Rm[WS(rs, 7)] = FMA(T3m, T3j, T3q);
					     T3b = FMA(KP707106781, T2W, T2L);
					     T2X = FNMS(KP707106781, T2W, T2L);
					}
				   }
				   T2Y = T2G * T2X;
				   T3a = T30 * T2X;
				   T3c = T1V * T3b;
				   T3d = FMA(KP707106781, T38, T35);
				   T39 = FNMS(KP707106781, T38, T35);
				   T3e = T1X * T3b;
				   T2q = FMA(KP707106781, T20, T1Z);
				   T21 = FNMS(KP707106781, T20, T1Z);
				   Rp[WS(rs, 1)] = FNMS(T1X, T3d, T3c);
				   Rm[WS(rs, 5)] = FMA(T2G, T39, T3a);
				   Rp[WS(rs, 5)] = FNMS(T30, T39, T2Y);
				   Rm[WS(rs, 1)] = FMA(T1V, T3d, T3e);
				   T28 = T24 + T27;
				   T2w = T27 - T24;
				   T2v = FNMS(KP707106781, T2e, T2d);
				   T2f = FMA(KP707106781, T2e, T2d);
				   T2i = T2g - T2h;
				   T2r = T2g + T2h;
			      }
			      {
				   E TX, T1k, T1x, T1A;
				   T1J = FMA(KP707106781, TW, TL);
				   TX = FNMS(KP707106781, TW, TL);
				   {
					E T2l, T29, T2n, T2j;
					T2l = FNMS(KP923879532, T28, T21);
					T29 = FMA(KP923879532, T28, T21);
					T2n = FMA(KP923879532, T2i, T2f);
					T2j = FNMS(KP923879532, T2i, T2f);
					{
					     E T2o, T2m, T2k, T2a;
					     T2o = Tz * T2l;
					     T2m = Tw * T2l;
					     T2k = T2c * T29;
					     T2a = T1Y * T29;
					     Im[WS(rs, 1)] = FMA(Tw, T2n, T2o);
					     Ip[WS(rs, 1)] = FNMS(Tz, T2n, T2m);
					     Im[WS(rs, 5)] = FMA(T1Y, T2j, T2k);
					     Ip[WS(rs, 5)] = FNMS(T2c, T2j, T2a);
					     T1k = T18 - T1j;
					     T1O = T1j + T18;
					}
				   }
				   T1N = FMA(KP707106781, T1w, T1t);
				   T1x = FNMS(KP707106781, T1w, T1t);
				   T1A = T1y - T1z;
				   T1K = T1y + T1z;
				   {
					E T1E, T1l, T1H, T1B;
					T1E = FMA(KP923879532, T1k, TX);
					T1l = FNMS(KP923879532, T1k, TX);
					T1H = FMA(KP923879532, T1A, T1x);
					T1B = FNMS(KP923879532, T1A, T1x);
					{
					     E T1I, T1F, T1C, T1m;
					     T1I = T1G * T1E;
					     T1F = T1D * T1E;
					     T1C = T1o * T1l;
					     T1m = TG * T1l;
					     Im[WS(rs, 2)] = FMA(T1D, T1H, T1I);
					     Ip[WS(rs, 2)] = FNMS(T1G, T1H, T1F);
					     Im[WS(rs, 6)] = FMA(TG, T1B, T1C);
					     Ip[WS(rs, 6)] = FNMS(T1o, T1B, T1m);
					}
				   }
				   {
					E T2A, T2s, T2D, T2x;
					T2A = FMA(KP923879532, T2r, T2q);
					T2s = FNMS(KP923879532, T2r, T2q);
					T2D = FNMS(KP923879532, T2w, T2v);
					T2x = FMA(KP923879532, T2w, T2v);
					{
					     E T2B, T2t, T2E, T2y;
					     T2B = T2z * T2A;
					     T2t = T2p * T2s;
					     T2E = T2z * T2D;
					     T2y = T2p * T2x;
					     Ip[WS(rs, 7)] = FNMS(T2C, T2D, T2B);
					     Ip[WS(rs, 3)] = FNMS(T2u, T2x, T2t);
					     Im[WS(rs, 7)] = FMA(T2C, T2A, T2E);
					     Im[WS(rs, 3)] = FMA(T2u, T2s, T2y);
					}
				   }
			      }
			 }
		    }
	       }
	  }
	  {
	       E T1L, T1R, T1P, T1T;
	       T1L = FNMS(KP923879532, T1K, T1J);
	       T1R = FMA(KP923879532, T1K, T1J);
	       T1P = FNMS(KP923879532, T1O, T1N);
	       T1T = FMA(KP923879532, T1O, T1N);
	       {
		    E T1S, T1M, T1U, T1Q;
		    T1S = Tv * T1R;
		    T1M = TB * T1L;
		    T1U = Tv * T1T;
		    T1Q = TB * T1P;
		    Ip[0] = FNMS(Ty, T1T, T1S);
		    Ip[WS(rs, 4)] = FNMS(TF, T1P, T1M);
		    Im[0] = FMA(Ty, T1R, T1U);
		    Im[WS(rs, 4)] = FMA(TF, T1L, T1Q);
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_CEXP, 1, 1},
     {TW_CEXP, 1, 3},
     {TW_CEXP, 1, 9},
     {TW_CEXP, 1, 15},
     {TW_NEXT, 1, 0}
};

static const hc2c_desc desc = { 16, "hc2cb2_16", twinstr, &GENUS, {104, 42, 92, 0} };

void X(codelet_hc2cb2_16) (planner *p) {
     X(khc2c_register) (p, hc2cb2_16, &desc, HC2C_VIA_RDFT);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2c -compact -variables 4 -pipeline-latency 4 -sign 1 -twiddle-log3 -precompute-twiddles -n 16 -dif -name hc2cb2_16 -include hc2cb.h */

/*
 * This function contains 196 FP additions, 108 FP multiplications,
 * (or, 156 additions, 68 multiplications, 40 fused multiply/add),
 * 80 stack variables, and 64 memory accesses
 */
#include "hc2cb.h"

static void hc2cb2_16(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT m;
     for (m = mb, W = W + ((mb - 1) * 8); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 8, MAKE_VOLATILE_STRIDE(rs)) {
	  E Tv, Ty, T1l, T1n, T1p, T1t, T27, T25, Tz, Tw, TB, T21, T1P, T1H, T1X;
	  E T17, T1L, T1N, T1v, T1w, T1x, T1B, T2F, T2T, T2b, T2R, T3j, T3x, T35, T3t;
	  {
	       E TA, T1J, T15, T1G, Tx, T1K, T16, T1F;
	       {
		    E T1m, T1s, T1o, T1r;
		    Tv = W[0];
		    Ty = W[1];
		    T1l = W[2];
		    T1n = W[3];
		    T1m = Tv * T1l;
		    T1s = Ty * T1l;
		    T1o = Ty * T1n;
		    T1r = Tv * T1n;
		    T1p = T1m + T1o;
		    T1t = T1r - T1s;
		    T27 = T1r + T1s;
		    T25 = T1m - T1o;
		    Tz = W[5];
		    TA = Ty * Tz;
		    T1J = T1l * Tz;
		    T15 = Tv * Tz;
		    T1G = T1n * Tz;
		    Tw = W[4];
		    Tx = Tv * Tw;
		    T1K = T1n * Tw;
		    T16 = Ty * Tw;
		    T1F = T1l * Tw;
	       }
	       TB = Tx - TA;
	       T21 = T1J + T1K;
	       T1P = T15 - T16;
	       T1H = T1F + T1G;
	       T1X = T1F - T1G;
	       T17 = T15 + T16;
	       T1L = T1J - T1K;
	       T1N = Tx + TA;
	       T1v = W[6];
	       T1w = W[7];
	       T1x = FMA(Tv, T1v, Ty * T1w);
	       T1B = FNMS(Ty, T1v, Tv * T1w);
	       {
		    E T2D, T2E, T29, T2a;
		    T2D = T25 * Tz;
		    T2E = T27 * Tw;
		    T2F = T2D + T2E;
		    T2T = T2D - T2E;
		    T29 = T25 * Tw;
		    T2a = T27 * Tz;
		    T2b = T29 - T2a;
		    T2R = T29 + T2a;
	       }
	       {
		    E T3h, T3i, T33, T34;
		    T3h = T1p * Tz;
		    T3i = T1t * Tw;
		    T3j = T3h + T3i;
		    T3x = T3h - T3i;
		    T33 = T1p * Tw;
		    T34 = T1t * Tz;
		    T35 = T33 - T34;
		    T3t = T33 + T34;
	       }
	  }
	  {
	       E T7, T36, T3k, TC, T1f, T2e, T2I, T1Q, Te, TJ, T1R, T18, T2L, T37, T2l;
	       E T3l, Tm, T1T, TT, T1h, T2A, T2N, T3b, T3n, Tt, T1U, T12, T1i, T2t, T2O;
	       E T3e, T3o;
	       {
		    E T3, T2c, T1b, T2H, T6, T2G, T1e, T2d;
		    {
			 E T1, T2, T19, T1a;
			 T1 = Rp[0];
			 T2 = Rm[WS(rs, 7)];
			 T3 = T1 + T2;
			 T2c = T1 - T2;
			 T19 = Ip[0];
			 T1a = Im[WS(rs, 7)];
			 T1b = T19 - T1a;
			 T2H = T19 + T1a;
		    }
		    {
			 E T4, T5, T1c, T1d;
			 T4 = Rp[WS(rs, 4)];
			 T5 = Rm[WS(rs, 3)];
			 T6 = T4 + T5;
			 T2G = T4 - T5;
			 T1c = Ip[WS(rs, 4)];
			 T1d = Im[WS(rs, 3)];
			 T1e = T1c - T1d;
			 T2d = T1c + T1d;
		    }
		    T7 = T3 + T6;
		    T36 = T2c + T2d;
		    T3k = T2H - T2G;
		    TC = T3 - T6;
		    T1f = T1b - T1e;
		    T2e = T2c - T2d;
		    T2I = T2G + T2H;
		    T1Q = T1b + T1e;
	       }
	       {
		    E Ta, T2f, TI, T2g, Td, T2i, TF, T2j;
		    {
			 E T8, T9, TG, TH;
			 T8 = Rp[WS(rs, 2)];
			 T9 = Rm[WS(rs, 5)];
			 Ta = T8 + T9;
			 T2f = T8 - T9;
			 TG = Ip[WS(rs, 2)];
			 TH = Im[WS(rs, 5)];
			 TI = TG - TH;
			 T2g = TG + TH;
		    }
		    {
			 E Tb, Tc, TD, TE;
			 Tb = Rm[WS(rs, 1)];
			 Tc = Rp[WS(rs, 6)];
			 Td = Tb + Tc;
			 T2i = Tb - Tc;
			 TD = Ip[WS(rs, 6)];
			 TE = Im[WS(rs, 1)];
			 TF = TD - TE;
			 T2j = TD + TE;
		    }
		    Te = Ta + Td;
		    TJ = TF - TI;
		    T1R = TI + TF;
		    T18 = Ta - Td;
		    {
			 E T2J, T2K, T2h, T2k;
			 T2J = T2f + T2g;
			 T2K = T2i + T2j;
			 T2L = KP707106781 * (T2J - T2K);
			 T37 = KP707106781 * (T2J + T2K);
			 T2h = T2f - T2g;
			 T2k = T2i - T2j;
			 T2l = KP707106781 * (T2h + T2k);
			 T3l = KP707106781 * (T2h - T2k);
		    }
	       }
	       {
		    E Ti, T2x, TO, T2v, Tl, T2u, TR, T2y, TL, TS;
		    {
			 E Tg, Th, TM, TN;
			 Tg = Rp[WS(rs, 1)];
			 Th = Rm[WS(rs, 6)];
			 Ti = Tg + Th;
			 T2x = Tg - Th;
			 TM = Ip[WS(rs, 1)];
			 TN = Im[WS(rs, 6)];
			 TO = TM - TN;
			 T2v = TM + TN;
		    }
		    {
			 E Tj, Tk, TP, TQ;
			 Tj = Rp[WS(rs, 5)];
			 Tk = Rm[WS(rs, 2)];
			 Tl = Tj + Tk;
			 T2u = Tj - Tk;
			 TP = Ip[WS(rs, 5)];
			 TQ = Im[WS(rs, 2)];
			 TR = TP - TQ;
			 T2y = TP + TQ;
		    }
		    Tm = Ti + Tl;
		    T1T = TO + TR;
		    TL = Ti - Tl;
		    TS = TO - TR;
		    TT = TL - TS;
		    T1h = TL + TS;
		    {
			 E T2w, T2z, T39, T3a;
			 T2w = T2u + T2v;
			 T2z = T2x - T2y;
			 T2A = FMA(KP923879532, T2w, KP382683432 * T2z);
			 T2N = FNMS(KP382683432, T2w, KP923879532 * T2z);
			 T39 = T2x + T2y;
			 T3a = T2v - T2u;
			 T3b = FNMS(KP923879532, T3a, KP382683432 * T39);
			 T3n = FMA(KP382683432, T3a, KP923879532 * T39);
		    }
	       }
	       {
		    E Tp, T2q, TX, T2o, Ts, T2n, T10, T2r, TU, T11;
		    {
			 E Tn, To, TV, TW;
			 Tn = Rm[0];
			 To = Rp[WS(rs, 7)];
			 Tp = Tn + To;
			 T2q = Tn - To;
			 TV = Ip[WS(rs, 7)];
			 TW = Im[0];
			 TX = TV - TW;
			 T2o = TV + TW;
		    }
		    {
			 E Tq, Tr, TY, TZ;
			 Tq = Rp[WS(rs, 3)];
			 Tr = Rm[WS(rs, 4)];
			 Ts = Tq + Tr;
			 T2n = Tq - Tr;
			 TY = Ip[WS(rs, 3)];
			 TZ = Im[WS(rs, 4)];
			 T10 = TY - TZ;
			 T2r = TY + TZ;
		    }
		    Tt = Tp + Ts;
		    T1U = TX + T10;
		    TU = Tp - Ts;
		    T11 = TX - T10;
		    T12 = TU + T11;
		    T1i = T11 - TU;
		    {
			 E T2p, T2s, T3c, T3d;
			 T2p = T2n - T2o;
			 T2s = T2q - T2r;
			 T2t = FNMS(KP382683432, T2s, KP923879532 * T2p);
			 T2O = FMA(KP382683432, T2p, KP923879532 * T2s);
			 T3c = T2q + T2r;
			 T3d = T2n + T2o;
			 T3e = FNMS(KP923879532, T3d, KP382683432 * T3c);
			 T3o = FMA(KP382683432, T3d, KP923879532 * T3c);
		    }
	       }
	       {
		    E Tf, Tu, T1O, T1S, T1V, T1W;
		    Tf = T7 + Te;
		    Tu = Tm + Tt;
		    T1O = Tf - Tu;
		    T1S = T1Q + T1R;
		    T1V = T1T + T1U;
		    T1W = T1S - T1V;
		    Rp[0] = Tf + Tu;
		    Rm[0] = T1S + T1V;
		    Rp[WS(rs, 4)] = FNMS(T1P, T1W, T1N * T1O);
		    Rm[WS(rs, 4)] = FMA(T1P, T1O, T1N * T1W);
	       }
	       {
		    E T3g, T3r, T3q, T3s;
		    {
			 E T38, T3f, T3m, T3p;
			 T38 = T36 - T37;
			 T3f = T3b + T3e;
			 T3g = T38 - T3f;
			 T3r = T38 + T3f;
			 T3m = T3k + T3l;
			 T3p = T3n - T3o;
			 T3q = T3m - T3p;
			 T3s = T3m + T3p;
		    }
		    Ip[WS(rs, 5)] = FNMS(T3j, T3q, T35 * T3g);
		    Im[WS(rs, 5)] = FMA(T3j, T3g, T35 * T3q);
		    Ip[WS(rs, 1)] = FNMS(T1n, T3s, T1l * T3r);
		    Im[WS(rs, 1)] = FMA(T1n, T3r, T1l * T3s);
	       }
	       {
		    E T3w, T3B, T3A, T3C;
		    {
			 E T3u, T3v, T3y, T3z;
			 T3u = T36 + T37;
			 T3v = T3n + T3o;
			 T3w = T3u - T3v;
			 T3B = T3u + T3v;
			 T3y = T3k - T3l;
			 T3z = T3b - T3e;
			 T3A = T3y + T3z;
			 T3C = T3y - T3z;
		    }
		    Ip[WS(rs, 3)] = FNMS(T3x, T3A, T3t * T3w);
		    Im[WS(rs, 3)] = FMA(T3t, T3A, T3x * T3w);
		    Ip[WS(rs, 7)] = FNMS(T1w, T3C, T1v * T3B);
		    Im[WS(rs, 7)] = FMA(T1v, T3C, T1w * T3B);
	       }
	       {
		    E T14, T1q, T1k, T1u;
		    {
			 E TK, T13, T1g, T1j;
			 TK = TC + TJ;
			 T13 = KP707106781 * (TT + T12);
			 T14 = TK - T13;
			 T1q = TK + T13;
			 T1g = T18 + T1f;
			 T1j = KP707106781 * (T1h + T1i);
			 T1k = T1g - T1j;
			 T1u = T1g + T1j;
		    }
		    Rp[WS(rs, 5)] = FNMS(T17, T1k, TB * T14);
		    Rm[WS(rs, 5)] = FMA(T17, T14, TB * T1k);
		    Rp[WS(rs, 1)] = FNMS(T1t, T1u, T1p * T1q);
		    Rm[WS(rs, 1)] = FMA(T1t, T1q, T1p * T1u);
	       }
	       {
		    E T1A, T1I, T1E, T1M;
		    {
			 E T1y, T1z, T1C, T1D;
			 T1y = TC - TJ;
			 T1z = KP707106781 * (T1i - T1h);
			 T1A = T1y - T1z;
			 T1I = T1y + T1z;
			 T1C = T1f - T18;
			 T1D = KP707106781 * (TT - T12);
			 T1E = T1C - T1D;
			 T1M = T1C + T1D;
		    }
		    Rp[WS(rs, 7)] = FNMS(T1B, T1E, T1x * T1A);
		    Rm[WS(rs, 7)] = FMA(T1x, T1E, T1B * T1A);
		    Rp[WS(rs, 3)] = FNMS(T1L, T1M, T1H * T1I);
		    Rm[WS(rs, 3)] = FMA(T1H, T1M, T1L * T1I);
	       }
	       {
		    E T2C, T2S, T2Q, T2U;
		    {
			 E T2m, T2B, T2M, T2P;
			 T2m = T2e - T2l;
			 T2B = T2t - T2A;
			 T2C = T2m - T2B;
			 T2S = T2m + T2B;
			 T2M = T2I - T2L;
			 T2P = T2N - T2O;
			 T2Q = T2M - T2P;
			 T2U = T2M + T2P;
		    }
		    Ip[WS(rs, 6)] = FNMS(T2F, T2Q, T2b * T2C);
		    Im[WS(rs, 6)] = FMA(T2F, T2C, T2b * T2Q);
		    Ip[WS(rs, 2)] = FNMS(T2T, T2U, T2R * T2S);
		    Im[WS(rs, 2)] = FMA(T2T, T2S, T2R * T2U);
	       }
	       {
		    E T2X, T31, T30, T32;
		    {
			 E T2V, T2W, T2Y, T2Z;
			 T2V = T2e + T2l;
			 T2W = T2N + T2O;
			 T2X = T2V - T2W;
			 T31 = T2V + T2W;
			 T2Y = T2I + T2L;
			 T2Z = T2A + T2t;
			 T30 = T2Y - T2Z;
			 T32 = T2Y + T2Z;
		    }
		    Ip[WS(rs, 4)] = FNMS(Tz, T30, Tw * T2X);
		    Im[WS(rs, 4)] = FMA(Tw, T30, Tz * T2X);
		    Ip[0] = FNMS(Ty, T32, Tv * T31);
		    Im[0] = FMA(Tv, T32, Ty * T31);
	       }
	       {
		    E T20, T26, T24, T28;
		    {
			 E T1Y, T1Z, T22, T23;
			 T1Y = T7 - Te;
			 T1Z = T1U - T1T;
			 T20 = T1Y - T1Z;
			 T26 = T1Y + T1Z;
			 T22 = T1Q - T1R;
			 T23 = Tm - Tt;
			 T24 = T22 - T23;
			 T28 = T23 + T22;
		    }
		    Rp[WS(rs, 6)] = FNMS(T21, T24, T1X * T20);
		    Rm[WS(rs, 6)] = FMA(T1X, T24, T21 * T20);
		    Rp[WS(rs, 2)] = FNMS(T27, T28, T25 * T26);
		    Rm[WS(rs, 2)] = FMA(T25, T28, T27 * T26);
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_CEXP, 1, 1},
     {TW_CEXP, 1, 3},
     {TW_CEXP, 1, 9},
     {TW_CEXP, 1, 15},
     {TW_NEXT, 1, 0}
};

static const hc2c_desc desc = { 16, "hc2cb2_16", twinstr, &GENUS, {156, 68, 40, 0} };

void X(codelet_hc2cb2_16) (planner *p) {
     X(khc2c_register) (p, hc2cb2_16, &desc, HC2C_VIA_RDFT);
}
#endif				/* HAVE_FMA */
