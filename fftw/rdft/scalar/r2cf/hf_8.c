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
/* Generated on Mon Mar 19 08:02:29 CDT 2007 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 8 -dit -name hf_8 -include hf.h */

/*
 * This function contains 66 FP additions, 36 FP multiplications,
 * (or, 44 additions, 14 multiplications, 22 fused multiply/add),
 * 61 stack variables, and 32 memory accesses
 */
#include "hf.h"

static void hf_8(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT m;
     for (m = mb, W = W + ((mb - 1) * 14); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 14, MAKE_VOLATILE_STRIDE(rs)) {
	  E T1f, T1g, T1e, Tm, T1q, T1o, T1p, TN, T1h, T1i;
	  {
	       E T1, T1m, T1l, T7, TS, Tk, TQ, Te, To, Tr, T17, TM, T12, Tu, TW;
	       E Tp, Tx, Tt, Tq, Tw;
	       {
		    E T3, T6, T2, T5;
		    T1 = cr[0];
		    T1m = ci[0];
		    T3 = cr[WS(rs, 4)];
		    T6 = ci[WS(rs, 4)];
		    T2 = W[6];
		    T5 = W[7];
		    {
			 E Ta, Td, T9, Tc;
			 {
			      E Tg, Tj, Ti, TR, Th, T1k, T4, Tf;
			      Tg = cr[WS(rs, 6)];
			      Tj = ci[WS(rs, 6)];
			      T1k = T2 * T6;
			      T4 = T2 * T3;
			      Tf = W[10];
			      Ti = W[11];
			      T1l = FNMS(T5, T3, T1k);
			      T7 = FMA(T5, T6, T4);
			      TR = Tf * Tj;
			      Th = Tf * Tg;
			      Ta = cr[WS(rs, 2)];
			      Td = ci[WS(rs, 2)];
			      TS = FNMS(Ti, Tg, TR);
			      Tk = FMA(Ti, Tj, Th);
			      T9 = W[2];
			      Tc = W[3];
			 }
			 {
			      E TB, TE, TH, T13, TC, TK, TG, TD, TJ, TP, Tb, TA, Tn;
			      TB = cr[WS(rs, 7)];
			      TE = ci[WS(rs, 7)];
			      TP = T9 * Td;
			      Tb = T9 * Ta;
			      TA = W[12];
			      TH = cr[WS(rs, 3)];
			      TQ = FNMS(Tc, Ta, TP);
			      Te = FMA(Tc, Td, Tb);
			      T13 = TA * TE;
			      TC = TA * TB;
			      TK = ci[WS(rs, 3)];
			      TG = W[4];
			      TD = W[13];
			      TJ = W[5];
			      {
				   E T14, TF, T16, TL, T15, TI;
				   To = cr[WS(rs, 1)];
				   T15 = TG * TK;
				   TI = TG * TH;
				   T14 = FNMS(TD, TB, T13);
				   TF = FMA(TD, TE, TC);
				   T16 = FNMS(TJ, TH, T15);
				   TL = FMA(TJ, TK, TI);
				   Tr = ci[WS(rs, 1)];
				   Tn = W[0];
				   T17 = T14 - T16;
				   T1f = T14 + T16;
				   TM = TF + TL;
				   T12 = TF - TL;
			      }
			      Tu = cr[WS(rs, 5)];
			      TW = Tn * Tr;
			      Tp = Tn * To;
			      Tx = ci[WS(rs, 5)];
			      Tt = W[8];
			      Tq = W[1];
			      Tw = W[9];
			 }
		    }
	       }
	       {
		    E T8, T1j, Tl, Tz, T1a, TU, T1n, T1b, T1c, T1v, T1t, T1u, T19, T1w, T1d;
		    {
			 E T1r, T10, TV, T1s, T11, T18;
			 {
			      E TO, TX, Ts, TZ, Ty, TT, TY, Tv;
			      T8 = T1 + T7;
			      TO = T1 - T7;
			      TY = Tt * Tx;
			      Tv = Tt * Tu;
			      TX = FNMS(Tq, To, TW);
			      Ts = FMA(Tq, Tr, Tp);
			      TZ = FNMS(Tw, Tu, TY);
			      Ty = FMA(Tw, Tx, Tv);
			      TT = TQ - TS;
			      T1j = TQ + TS;
			      Tl = Te + Tk;
			      T1r = Te - Tk;
			      T10 = TX - TZ;
			      T1g = TX + TZ;
			      Tz = Ts + Ty;
			      TV = Ts - Ty;
			      T1a = TO - TT;
			      TU = TO + TT;
			      T1s = T1m - T1l;
			      T1n = T1l + T1m;
			 }
			 T1b = TV - T10;
			 T11 = TV + T10;
			 T18 = T12 - T17;
			 T1c = T12 + T17;
			 T1v = T1s - T1r;
			 T1t = T1r + T1s;
			 T1u = T18 - T11;
			 T19 = T11 + T18;
		    }
		    ci[WS(rs, 4)] = FMA(KP707106781, T1u, T1t);
		    cr[WS(rs, 7)] = FMS(KP707106781, T1u, T1t);
		    cr[WS(rs, 1)] = FMA(KP707106781, T19, TU);
		    ci[WS(rs, 2)] = FNMS(KP707106781, T19, TU);
		    T1w = T1c - T1b;
		    T1d = T1b + T1c;
		    ci[WS(rs, 6)] = FMA(KP707106781, T1w, T1v);
		    cr[WS(rs, 5)] = FMS(KP707106781, T1w, T1v);
		    ci[0] = FMA(KP707106781, T1d, T1a);
		    cr[WS(rs, 3)] = FNMS(KP707106781, T1d, T1a);
		    T1e = T8 - Tl;
		    Tm = T8 + Tl;
		    T1q = T1n - T1j;
		    T1o = T1j + T1n;
		    T1p = TM - Tz;
		    TN = Tz + TM;
	       }
	  }
	  ci[WS(rs, 5)] = T1p + T1q;
	  cr[WS(rs, 6)] = T1p - T1q;
	  cr[0] = Tm + TN;
	  ci[WS(rs, 3)] = Tm - TN;
	  T1h = T1f - T1g;
	  T1i = T1g + T1f;
	  ci[WS(rs, 7)] = T1i + T1o;
	  cr[WS(rs, 4)] = T1i - T1o;
	  ci[WS(rs, 1)] = T1e + T1h;
	  cr[WS(rs, 2)] = T1e - T1h;
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 8},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 8, "hf_8", twinstr, &GENUS, {44, 14, 22, 0} };

void X(codelet_hf_8) (planner *p) {
     X(khc2hc_register) (p, hf_8, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2hc -compact -variables 4 -pipeline-latency 4 -n 8 -dit -name hf_8 -include hf.h */

/*
 * This function contains 66 FP additions, 32 FP multiplications,
 * (or, 52 additions, 18 multiplications, 14 fused multiply/add),
 * 28 stack variables, and 32 memory accesses
 */
#include "hf.h"

static void hf_8(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT m;
     for (m = mb, W = W + ((mb - 1) * 14); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 14, MAKE_VOLATILE_STRIDE(rs)) {
	  E T7, T1f, TH, T19, TF, T12, TR, TU, Ti, T1e, TK, T16, Tu, T13, TM;
	  E TP;
	  {
	       E T1, T18, T6, T17;
	       T1 = cr[0];
	       T18 = ci[0];
	       {
		    E T3, T5, T2, T4;
		    T3 = cr[WS(rs, 4)];
		    T5 = ci[WS(rs, 4)];
		    T2 = W[6];
		    T4 = W[7];
		    T6 = FMA(T2, T3, T4 * T5);
		    T17 = FNMS(T4, T3, T2 * T5);
	       }
	       T7 = T1 + T6;
	       T1f = T18 - T17;
	       TH = T1 - T6;
	       T19 = T17 + T18;
	  }
	  {
	       E Tz, TS, TE, TT;
	       {
		    E Tw, Ty, Tv, Tx;
		    Tw = cr[WS(rs, 7)];
		    Ty = ci[WS(rs, 7)];
		    Tv = W[12];
		    Tx = W[13];
		    Tz = FMA(Tv, Tw, Tx * Ty);
		    TS = FNMS(Tx, Tw, Tv * Ty);
	       }
	       {
		    E TB, TD, TA, TC;
		    TB = cr[WS(rs, 3)];
		    TD = ci[WS(rs, 3)];
		    TA = W[4];
		    TC = W[5];
		    TE = FMA(TA, TB, TC * TD);
		    TT = FNMS(TC, TB, TA * TD);
	       }
	       TF = Tz + TE;
	       T12 = TS + TT;
	       TR = Tz - TE;
	       TU = TS - TT;
	  }
	  {
	       E Tc, TI, Th, TJ;
	       {
		    E T9, Tb, T8, Ta;
		    T9 = cr[WS(rs, 2)];
		    Tb = ci[WS(rs, 2)];
		    T8 = W[2];
		    Ta = W[3];
		    Tc = FMA(T8, T9, Ta * Tb);
		    TI = FNMS(Ta, T9, T8 * Tb);
	       }
	       {
		    E Te, Tg, Td, Tf;
		    Te = cr[WS(rs, 6)];
		    Tg = ci[WS(rs, 6)];
		    Td = W[10];
		    Tf = W[11];
		    Th = FMA(Td, Te, Tf * Tg);
		    TJ = FNMS(Tf, Te, Td * Tg);
	       }
	       Ti = Tc + Th;
	       T1e = Tc - Th;
	       TK = TI - TJ;
	       T16 = TI + TJ;
	  }
	  {
	       E To, TN, Tt, TO;
	       {
		    E Tl, Tn, Tk, Tm;
		    Tl = cr[WS(rs, 1)];
		    Tn = ci[WS(rs, 1)];
		    Tk = W[0];
		    Tm = W[1];
		    To = FMA(Tk, Tl, Tm * Tn);
		    TN = FNMS(Tm, Tl, Tk * Tn);
	       }
	       {
		    E Tq, Ts, Tp, Tr;
		    Tq = cr[WS(rs, 5)];
		    Ts = ci[WS(rs, 5)];
		    Tp = W[8];
		    Tr = W[9];
		    Tt = FMA(Tp, Tq, Tr * Ts);
		    TO = FNMS(Tr, Tq, Tp * Ts);
	       }
	       Tu = To + Tt;
	       T13 = TN + TO;
	       TM = To - Tt;
	       TP = TN - TO;
	  }
	  {
	       E Tj, TG, T1b, T1c;
	       Tj = T7 + Ti;
	       TG = Tu + TF;
	       ci[WS(rs, 3)] = Tj - TG;
	       cr[0] = Tj + TG;
	       T1b = TF - Tu;
	       T1c = T19 - T16;
	       cr[WS(rs, 6)] = T1b - T1c;
	       ci[WS(rs, 5)] = T1b + T1c;
	       {
		    E TX, T1i, T10, T1h, TY, TZ;
		    TX = TH - TK;
		    T1i = T1f - T1e;
		    TY = TM - TP;
		    TZ = TR + TU;
		    T10 = KP707106781 * (TY + TZ);
		    T1h = KP707106781 * (TZ - TY);
		    cr[WS(rs, 3)] = TX - T10;
		    ci[WS(rs, 6)] = T1h + T1i;
		    ci[0] = TX + T10;
		    cr[WS(rs, 5)] = T1h - T1i;
	       }
	  }
	  {
	       E T15, T1a, T11, T14;
	       T15 = T13 + T12;
	       T1a = T16 + T19;
	       cr[WS(rs, 4)] = T15 - T1a;
	       ci[WS(rs, 7)] = T15 + T1a;
	       T11 = T7 - Ti;
	       T14 = T12 - T13;
	       cr[WS(rs, 2)] = T11 - T14;
	       ci[WS(rs, 1)] = T11 + T14;
	       {
		    E TL, T1g, TW, T1d, TQ, TV;
		    TL = TH + TK;
		    T1g = T1e + T1f;
		    TQ = TM + TP;
		    TV = TR - TU;
		    TW = KP707106781 * (TQ + TV);
		    T1d = KP707106781 * (TV - TQ);
		    ci[WS(rs, 2)] = TL - TW;
		    ci[WS(rs, 4)] = T1d + T1g;
		    cr[WS(rs, 1)] = TL + TW;
		    cr[WS(rs, 7)] = T1d - T1g;
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 8},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 8, "hf_8", twinstr, &GENUS, {52, 18, 14, 0} };

void X(codelet_hf_8) (planner *p) {
     X(khc2hc_register) (p, hf_8, &desc);
}
#endif				/* HAVE_FMA */
