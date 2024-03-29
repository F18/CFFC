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
/* Generated on Mon Mar 19 08:00:17 CDT 2007 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2c -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 12 -dif -name hc2cb_12 -include hc2cb.h */

/*
 * This function contains 118 FP additions, 68 FP multiplications,
 * (or, 72 additions, 22 multiplications, 46 fused multiply/add),
 * 64 stack variables, and 48 memory accesses
 */
#include "hc2cb.h"

static void hc2cb_12(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT m;
     for (m = mb, W = W + ((mb - 1) * 22); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 22, MAKE_VOLATILE_STRIDE(rs)) {
	  E T1U, T1X, T1W, T1Y, T1V;
	  {
	       E T18, T20, T21, T1b, T2a, T1s, T29, T1p, TO, T11, To, Tb, Tg, T23, T1f;
	       E Tl, Ty, Tt, T1i, T24, T1z, T2d, T1w, T2c;
	       {
		    E T5, Ta, TN, TI;
		    {
			 E T1, TE, T6, TM, T7, T1o, T4, T17, TH, T8, TJ, TK;
			 T1 = Rp[0];
			 TE = Ip[0];
			 T6 = Rm[WS(rs, 5)];
			 TM = Im[WS(rs, 5)];
			 {
			      E T2, T3, TF, TG;
			      T2 = Rp[WS(rs, 4)];
			      T3 = Rm[WS(rs, 3)];
			      TF = Ip[WS(rs, 4)];
			      TG = Im[WS(rs, 3)];
			      T7 = Rm[WS(rs, 1)];
			      T1o = T2 - T3;
			      T4 = T2 + T3;
			      T17 = TF + TG;
			      TH = TF - TG;
			      T8 = Rp[WS(rs, 2)];
			      TJ = Ip[WS(rs, 2)];
			      TK = Im[WS(rs, 1)];
			 }
			 {
			      E T1r, T1a, T19, T1q, T9, TL, T16, T1n;
			      T5 = T1 + T4;
			      T16 = FNMS(KP500000000, T4, T1);
			      T1r = T7 - T8;
			      T9 = T7 + T8;
			      T1a = TJ + TK;
			      TL = TJ - TK;
			      T18 = FNMS(KP866025403, T17, T16);
			      T20 = FMA(KP866025403, T17, T16);
			      T19 = FNMS(KP500000000, T9, T6);
			      Ta = T6 + T9;
			      TN = TL - TM;
			      T1q = FMA(KP500000000, TL, TM);
			      T1n = FNMS(KP500000000, TH, TE);
			      TI = TE + TH;
			      T21 = FNMS(KP866025403, T1a, T19);
			      T1b = FMA(KP866025403, T1a, T19);
			      T2a = FMA(KP866025403, T1r, T1q);
			      T1s = FNMS(KP866025403, T1r, T1q);
			      T29 = FNMS(KP866025403, T1o, T1n);
			      T1p = FMA(KP866025403, T1o, T1n);
			 }
		    }
		    {
			 E Tc, Tp, Th, Tx, Ti, Tf, T1v, Ts, T1e, Tj, Tu, Tv;
			 Tc = Rp[WS(rs, 3)];
			 TO = TI - TN;
			 T11 = TI + TN;
			 Tp = Ip[WS(rs, 3)];
			 To = T5 - Ta;
			 Tb = T5 + Ta;
			 Th = Rm[WS(rs, 2)];
			 Tx = Im[WS(rs, 2)];
			 {
			      E Td, Te, Tq, Tr;
			      Td = Rm[WS(rs, 4)];
			      Te = Rm[0];
			      Tq = Im[WS(rs, 4)];
			      Tr = Im[0];
			      Ti = Rp[WS(rs, 1)];
			      Tf = Td + Te;
			      T1v = Td - Te;
			      Ts = Tq + Tr;
			      T1e = Tq - Tr;
			      Tj = Rp[WS(rs, 5)];
			      Tu = Ip[WS(rs, 1)];
			      Tv = Ip[WS(rs, 5)];
			 }
			 {
			      E T1y, T1h, T1g, T1x, Tk, Tw, T1d, T1u;
			      T1d = FNMS(KP500000000, Tf, Tc);
			      Tg = Tc + Tf;
			      Tk = Ti + Tj;
			      T1y = Ti - Tj;
			      Tw = Tu + Tv;
			      T1h = Tv - Tu;
			      T23 = FNMS(KP866025403, T1e, T1d);
			      T1f = FMA(KP866025403, T1e, T1d);
			      Tl = Th + Tk;
			      T1g = FNMS(KP500000000, Tk, Th);
			      T1x = FMA(KP500000000, Tw, Tx);
			      Ty = Tw - Tx;
			      Tt = Tp - Ts;
			      T1u = FMA(KP500000000, Ts, Tp);
			      T1i = FMA(KP866025403, T1h, T1g);
			      T24 = FNMS(KP866025403, T1h, T1g);
			      T1z = FNMS(KP866025403, T1y, T1x);
			      T2d = FMA(KP866025403, T1y, T1x);
			      T1w = FMA(KP866025403, T1v, T1u);
			      T2c = FNMS(KP866025403, T1v, T1u);
			 }
		    }
	       }
	       {
		    E TY, T13, TX, T10;
		    {
			 E Tn, T12, TC, Tm, TD, TS, TA, Tz;
			 Tn = W[16];
			 T12 = Tt + Ty;
			 Tz = Tt - Ty;
			 TC = W[17];
			 Tm = Tg + Tl;
			 TD = Tg - Tl;
			 TS = To + Tz;
			 TA = To - Tz;
			 {
			      E TV, TU, TW, TT;
			      {
				   E TQ, TR, TP, TB;
				   TV = TO - TD;
				   TP = TD + TO;
				   Rp[0] = Tb + Tm;
				   TB = Tn * TA;
				   TQ = Tn * TP;
				   TR = W[4];
				   Ip[WS(rs, 4)] = FNMS(TC, TP, TB);
				   TU = W[5];
				   Im[WS(rs, 4)] = FMA(TC, TA, TQ);
				   TW = TR * TV;
				   TT = TR * TS;
			      }
			      Im[WS(rs, 1)] = FMA(TU, TS, TW);
			      Ip[WS(rs, 1)] = FNMS(TU, TV, TT);
			      TY = Tb - Tm;
			      T13 = T11 - T12;
			      TX = W[10];
			      T10 = W[11];
			      Rm[0] = T11 + T12;
			 }
		    }
		    {
			 E T1K, T1Q, T1P, T1L, T2o, T2u, T2t, T2p;
			 {
			      E T1E, T1D, T1H, T1F, T1G, T1t, T1k, T1A;
			      {
				   E T1c, TZ, T14, T1j;
				   T1K = T18 - T1b;
				   T1c = T18 + T1b;
				   TZ = TX * TY;
				   T14 = T10 * TY;
				   T1j = T1f + T1i;
				   T1Q = T1f - T1i;
				   T1P = T1p + T1s;
				   T1t = T1p - T1s;
				   Rp[WS(rs, 3)] = FNMS(T10, T13, TZ);
				   Rm[WS(rs, 3)] = FMA(TX, T13, T14);
				   T1E = T1c + T1j;
				   T1k = T1c - T1j;
				   T1A = T1w - T1z;
				   T1L = T1w + T1z;
			      }
			      {
				   E T15, T1m, T1B, T1l, T1C;
				   T15 = W[18];
				   T1m = W[19];
				   T1D = W[6];
				   T1H = T1t + T1A;
				   T1B = T1t - T1A;
				   T1l = T15 * T1k;
				   T1C = T1m * T1k;
				   T1F = T1D * T1E;
				   T1G = W[7];
				   Rp[WS(rs, 5)] = FNMS(T1m, T1B, T1l);
				   Rm[WS(rs, 5)] = FMA(T15, T1B, T1C);
			      }
			      {
				   E T26, T2i, T2l, T2f, T1Z, T28;
				   {
					E T22, T1I, T25, T2b, T2e;
					T22 = T20 + T21;
					T2o = T20 - T21;
					Rp[WS(rs, 2)] = FNMS(T1G, T1H, T1F);
					T1I = T1G * T1E;
					T2u = T23 - T24;
					T25 = T23 + T24;
					T2b = T29 - T2a;
					T2t = T29 + T2a;
					T2p = T2c + T2d;
					T2e = T2c - T2d;
					Rm[WS(rs, 2)] = FMA(T1D, T1H, T1I);
					T26 = T22 - T25;
					T2i = T22 + T25;
					T2l = T2b + T2e;
					T2f = T2b - T2e;
				   }
				   T1Z = W[2];
				   T28 = W[3];
				   {
					E T2h, T2k, T27, T2g, T2j, T2m;
					T2h = W[14];
					T2k = W[15];
					T27 = T1Z * T26;
					T2g = T28 * T26;
					T2j = T2h * T2i;
					T2m = T2k * T2i;
					Rp[WS(rs, 1)] = FNMS(T28, T2f, T27);
					Rm[WS(rs, 1)] = FMA(T1Z, T2f, T2g);
					Rp[WS(rs, 4)] = FNMS(T2k, T2l, T2j);
					Rm[WS(rs, 4)] = FMA(T2h, T2l, T2m);
				   }
			      }
			 }
			 {
			      E T2y, T2B, T2A, T2C, T2z;
			      {
				   E T2n, T2q, T2v, T2s, T2r, T2x, T2w;
				   T2n = W[8];
				   T2y = T2o + T2p;
				   T2q = T2o - T2p;
				   T2B = T2t - T2u;
				   T2v = T2t + T2u;
				   T2s = W[9];
				   T2r = T2n * T2q;
				   T2x = W[20];
				   T2w = T2n * T2v;
				   T2A = W[21];
				   Ip[WS(rs, 2)] = FNMS(T2s, T2v, T2r);
				   T2C = T2x * T2B;
				   T2z = T2x * T2y;
				   Im[WS(rs, 2)] = FMA(T2s, T2q, T2w);
			      }
			      Im[WS(rs, 5)] = FMA(T2A, T2y, T2C);
			      Ip[WS(rs, 5)] = FNMS(T2A, T2B, T2z);
			      {
				   E T1J, T1M, T1R, T1O, T1N, T1T, T1S;
				   T1J = W[0];
				   T1U = T1K + T1L;
				   T1M = T1K - T1L;
				   T1X = T1P - T1Q;
				   T1R = T1P + T1Q;
				   T1O = W[1];
				   T1N = T1J * T1M;
				   T1T = W[12];
				   T1S = T1J * T1R;
				   T1W = W[13];
				   Ip[0] = FNMS(T1O, T1R, T1N);
				   T1Y = T1T * T1X;
				   T1V = T1T * T1U;
				   Im[0] = FMA(T1O, T1M, T1S);
			      }
			 }
		    }
	       }
	  }
	  Im[WS(rs, 3)] = FMA(T1W, T1U, T1Y);
	  Ip[WS(rs, 3)] = FNMS(T1W, T1X, T1V);
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 12},
     {TW_NEXT, 1, 0}
};

static const hc2c_desc desc = { 12, "hc2cb_12", twinstr, &GENUS, {72, 22, 46, 0} };

void X(codelet_hc2cb_12) (planner *p) {
     X(khc2c_register) (p, hc2cb_12, &desc, HC2C_VIA_RDFT);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2c -compact -variables 4 -pipeline-latency 4 -sign 1 -n 12 -dif -name hc2cb_12 -include hc2cb.h */

/*
 * This function contains 118 FP additions, 60 FP multiplications,
 * (or, 88 additions, 30 multiplications, 30 fused multiply/add),
 * 39 stack variables, and 48 memory accesses
 */
#include "hc2cb.h"

static void hc2cb_12(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     INT m;
     for (m = mb, W = W + ((mb - 1) * 22); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 22, MAKE_VOLATILE_STRIDE(rs)) {
	  E T5, TH, T12, T1M, T1i, T1U, Tl, Ty, T1c, T1Y, T1s, T1Q, Ta, TM, T15;
	  E T1N, T1l, T1V, Tg, Tt, T19, T1X, T1p, T1P;
	  {
	       E T1, TD, T4, T1g, TG, T11, T10, T1h;
	       T1 = Rp[0];
	       TD = Ip[0];
	       {
		    E T2, T3, TE, TF;
		    T2 = Rp[WS(rs, 4)];
		    T3 = Rm[WS(rs, 3)];
		    T4 = T2 + T3;
		    T1g = KP866025403 * (T2 - T3);
		    TE = Ip[WS(rs, 4)];
		    TF = Im[WS(rs, 3)];
		    TG = TE - TF;
		    T11 = KP866025403 * (TE + TF);
	       }
	       T5 = T1 + T4;
	       TH = TD + TG;
	       T10 = FNMS(KP500000000, T4, T1);
	       T12 = T10 - T11;
	       T1M = T10 + T11;
	       T1h = FNMS(KP500000000, TG, TD);
	       T1i = T1g + T1h;
	       T1U = T1h - T1g;
	  }
	  {
	       E Th, Tx, Tk, T1a, Tw, T1r, T1b, T1q;
	       Th = Rm[WS(rs, 2)];
	       Tx = Im[WS(rs, 2)];
	       {
		    E Ti, Tj, Tu, Tv;
		    Ti = Rp[WS(rs, 1)];
		    Tj = Rp[WS(rs, 5)];
		    Tk = Ti + Tj;
		    T1a = KP866025403 * (Ti - Tj);
		    Tu = Ip[WS(rs, 1)];
		    Tv = Ip[WS(rs, 5)];
		    Tw = Tu + Tv;
		    T1r = KP866025403 * (Tv - Tu);
	       }
	       Tl = Th + Tk;
	       Ty = Tw - Tx;
	       T1b = FMA(KP500000000, Tw, Tx);
	       T1c = T1a - T1b;
	       T1Y = T1a + T1b;
	       T1q = FNMS(KP500000000, Tk, Th);
	       T1s = T1q + T1r;
	       T1Q = T1q - T1r;
	  }
	  {
	       E T6, TL, T9, T1j, TK, T14, T13, T1k;
	       T6 = Rm[WS(rs, 5)];
	       TL = Im[WS(rs, 5)];
	       {
		    E T7, T8, TI, TJ;
		    T7 = Rm[WS(rs, 1)];
		    T8 = Rp[WS(rs, 2)];
		    T9 = T7 + T8;
		    T1j = KP866025403 * (T7 - T8);
		    TI = Ip[WS(rs, 2)];
		    TJ = Im[WS(rs, 1)];
		    TK = TI - TJ;
		    T14 = KP866025403 * (TI + TJ);
	       }
	       Ta = T6 + T9;
	       TM = TK - TL;
	       T13 = FNMS(KP500000000, T9, T6);
	       T15 = T13 + T14;
	       T1N = T13 - T14;
	       T1k = FMA(KP500000000, TK, TL);
	       T1l = T1j - T1k;
	       T1V = T1j + T1k;
	  }
	  {
	       E Tc, Tp, Tf, T17, Ts, T1o, T18, T1n;
	       Tc = Rp[WS(rs, 3)];
	       Tp = Ip[WS(rs, 3)];
	       {
		    E Td, Te, Tq, Tr;
		    Td = Rm[WS(rs, 4)];
		    Te = Rm[0];
		    Tf = Td + Te;
		    T17 = KP866025403 * (Td - Te);
		    Tq = Im[WS(rs, 4)];
		    Tr = Im[0];
		    Ts = Tq + Tr;
		    T1o = KP866025403 * (Tq - Tr);
	       }
	       Tg = Tc + Tf;
	       Tt = Tp - Ts;
	       T18 = FMA(KP500000000, Ts, Tp);
	       T19 = T17 + T18;
	       T1X = T18 - T17;
	       T1n = FNMS(KP500000000, Tf, Tc);
	       T1p = T1n + T1o;
	       T1P = T1n - T1o;
	  }
	  {
	       E Tb, Tm, TU, TW, TX, TY, TT, TV;
	       Tb = T5 + Ta;
	       Tm = Tg + Tl;
	       TU = Tb - Tm;
	       TW = TH + TM;
	       TX = Tt + Ty;
	       TY = TW - TX;
	       Rp[0] = Tb + Tm;
	       Rm[0] = TW + TX;
	       TT = W[10];
	       TV = W[11];
	       Rp[WS(rs, 3)] = FNMS(TV, TY, TT * TU);
	       Rm[WS(rs, 3)] = FMA(TV, TU, TT * TY);
	  }
	  {
	       E TA, TQ, TO, TS;
	       {
		    E To, Tz, TC, TN;
		    To = T5 - Ta;
		    Tz = Tt - Ty;
		    TA = To - Tz;
		    TQ = To + Tz;
		    TC = Tg - Tl;
		    TN = TH - TM;
		    TO = TC + TN;
		    TS = TN - TC;
	       }
	       {
		    E Tn, TB, TP, TR;
		    Tn = W[16];
		    TB = W[17];
		    Ip[WS(rs, 4)] = FNMS(TB, TO, Tn * TA);
		    Im[WS(rs, 4)] = FMA(Tn, TO, TB * TA);
		    TP = W[4];
		    TR = W[5];
		    Ip[WS(rs, 1)] = FNMS(TR, TS, TP * TQ);
		    Im[WS(rs, 1)] = FMA(TP, TS, TR * TQ);
	       }
	  }
	  {
	       E T28, T2e, T2c, T2g;
	       {
		    E T26, T27, T2a, T2b;
		    T26 = T1M - T1N;
		    T27 = T1X + T1Y;
		    T28 = T26 - T27;
		    T2e = T26 + T27;
		    T2a = T1U + T1V;
		    T2b = T1P - T1Q;
		    T2c = T2a + T2b;
		    T2g = T2a - T2b;
	       }
	       {
		    E T25, T29, T2d, T2f;
		    T25 = W[8];
		    T29 = W[9];
		    Ip[WS(rs, 2)] = FNMS(T29, T2c, T25 * T28);
		    Im[WS(rs, 2)] = FMA(T25, T2c, T29 * T28);
		    T2d = W[20];
		    T2f = W[21];
		    Ip[WS(rs, 5)] = FNMS(T2f, T2g, T2d * T2e);
		    Im[WS(rs, 5)] = FMA(T2d, T2g, T2f * T2e);
	       }
	  }
	  {
	       E T1S, T22, T20, T24;
	       {
		    E T1O, T1R, T1W, T1Z;
		    T1O = T1M + T1N;
		    T1R = T1P + T1Q;
		    T1S = T1O - T1R;
		    T22 = T1O + T1R;
		    T1W = T1U - T1V;
		    T1Z = T1X - T1Y;
		    T20 = T1W - T1Z;
		    T24 = T1W + T1Z;
	       }
	       {
		    E T1L, T1T, T21, T23;
		    T1L = W[2];
		    T1T = W[3];
		    Rp[WS(rs, 1)] = FNMS(T1T, T20, T1L * T1S);
		    Rm[WS(rs, 1)] = FMA(T1T, T1S, T1L * T20);
		    T21 = W[14];
		    T23 = W[15];
		    Rp[WS(rs, 4)] = FNMS(T23, T24, T21 * T22);
		    Rm[WS(rs, 4)] = FMA(T23, T22, T21 * T24);
	       }
	  }
	  {
	       E T1C, T1I, T1G, T1K;
	       {
		    E T1A, T1B, T1E, T1F;
		    T1A = T12 + T15;
		    T1B = T1p + T1s;
		    T1C = T1A - T1B;
		    T1I = T1A + T1B;
		    T1E = T1i + T1l;
		    T1F = T19 + T1c;
		    T1G = T1E - T1F;
		    T1K = T1E + T1F;
	       }
	       {
		    E T1z, T1D, T1H, T1J;
		    T1z = W[18];
		    T1D = W[19];
		    Rp[WS(rs, 5)] = FNMS(T1D, T1G, T1z * T1C);
		    Rm[WS(rs, 5)] = FMA(T1D, T1C, T1z * T1G);
		    T1H = W[6];
		    T1J = W[7];
		    Rp[WS(rs, 2)] = FNMS(T1J, T1K, T1H * T1I);
		    Rm[WS(rs, 2)] = FMA(T1J, T1I, T1H * T1K);
	       }
	  }
	  {
	       E T1e, T1w, T1u, T1y;
	       {
		    E T16, T1d, T1m, T1t;
		    T16 = T12 - T15;
		    T1d = T19 - T1c;
		    T1e = T16 - T1d;
		    T1w = T16 + T1d;
		    T1m = T1i - T1l;
		    T1t = T1p - T1s;
		    T1u = T1m + T1t;
		    T1y = T1m - T1t;
	       }
	       {
		    E TZ, T1f, T1v, T1x;
		    TZ = W[0];
		    T1f = W[1];
		    Ip[0] = FNMS(T1f, T1u, TZ * T1e);
		    Im[0] = FMA(TZ, T1u, T1f * T1e);
		    T1v = W[12];
		    T1x = W[13];
		    Ip[WS(rs, 3)] = FNMS(T1x, T1y, T1v * T1w);
		    Im[WS(rs, 3)] = FMA(T1v, T1y, T1x * T1w);
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 12},
     {TW_NEXT, 1, 0}
};

static const hc2c_desc desc = { 12, "hc2cb_12", twinstr, &GENUS, {88, 30, 30, 0} };

void X(codelet_hc2cb_12) (planner *p) {
     X(khc2c_register) (p, hc2cb_12, &desc, HC2C_VIA_RDFT);
}
#endif				/* HAVE_FMA */
