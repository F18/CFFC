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
/* Generated on Mon Mar 19 07:44:03 CDT 2007 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_notw_c -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -n 15 -name n1fv_15 -include n1f.h */

/*
 * This function contains 78 FP additions, 49 FP multiplications,
 * (or, 36 additions, 7 multiplications, 42 fused multiply/add),
 * 78 stack variables, and 30 memory accesses
 */
#include "n1f.h"

static void n1fv_15(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP823639103, +0.823639103546331925877420039278190003029660514);
     DVK(KP910592997, +0.910592997310029334643087372129977886038870291);
     DVK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DVK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DVK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DVK(KP618033988, +0.618033988749894848204586834365638117720309180);
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     const R *xi;
     R *xo;
     xi = ri;
     xo = ro;
     for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  V Tb, TX, TM, TQ, Th, TB, T5, Ti, Ta, TC, TN, Te, TG, Tq, Tj;
	  V T1, T2, T3;
	  T1 = LD(&(xi[0]), ivs, &(xi[0]));
	  T2 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
	  T3 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
	  {
	       V T6, T7, T8, Tm, Tn, To;
	       T6 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
	       T7 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
	       T8 = LD(&(xi[WS(is, 13)]), ivs, &(xi[WS(is, 1)]));
	       Tm = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
	       Tn = LD(&(xi[WS(is, 14)]), ivs, &(xi[0]));
	       To = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
	       {
		    V T4, Tc, T9, Td, Tp;
		    Tb = LD(&(xi[WS(is, 12)]), ivs, &(xi[0]));
		    T4 = VADD(T2, T3);
		    TX = VSUB(T3, T2);
		    Tc = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    TM = VSUB(T8, T7);
		    T9 = VADD(T7, T8);
		    Td = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
		    Tp = VADD(Tn, To);
		    TQ = VSUB(To, Tn);
		    Th = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
		    TB = VFNMS(LDK(KP500000000), T4, T1);
		    T5 = VADD(T1, T4);
		    Ti = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
		    Ta = VADD(T6, T9);
		    TC = VFNMS(LDK(KP500000000), T9, T6);
		    TN = VSUB(Td, Tc);
		    Te = VADD(Tc, Td);
		    TG = VFNMS(LDK(KP500000000), Tp, Tm);
		    Tq = VADD(Tm, Tp);
		    Tj = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
	       }
	  }
	  {
	       V TY, TO, Tf, TD, TP, Tk;
	       TY = VADD(TM, TN);
	       TO = VSUB(TM, TN);
	       Tf = VADD(Tb, Te);
	       TD = VFNMS(LDK(KP500000000), Te, Tb);
	       TP = VSUB(Tj, Ti);
	       Tk = VADD(Ti, Tj);
	       {
		    V Tx, Tg, TE, TU, TZ, TR, Tl, TF;
		    Tx = VSUB(Ta, Tf);
		    Tg = VADD(Ta, Tf);
		    TE = VADD(TC, TD);
		    TU = VSUB(TC, TD);
		    TZ = VADD(TP, TQ);
		    TR = VSUB(TP, TQ);
		    Tl = VADD(Th, Tk);
		    TF = VFNMS(LDK(KP500000000), Tk, Th);
		    {
			 V T12, T10, T18, TS, Tw, Tr, TH, TV, T11, T1g;
			 T12 = VSUB(TY, TZ);
			 T10 = VADD(TY, TZ);
			 T18 = VFNMS(LDK(KP618033988), TO, TR);
			 TS = VFMA(LDK(KP618033988), TR, TO);
			 Tw = VSUB(Tl, Tq);
			 Tr = VADD(Tl, Tq);
			 TH = VADD(TF, TG);
			 TV = VSUB(TF, TG);
			 T11 = VFNMS(LDK(KP250000000), T10, TX);
			 T1g = VMUL(LDK(KP866025403), VADD(TX, T10));
			 {
			      V TA, Ty, Tu, TK, TI, T1a, TW, T1b, T13, Tt, Ts, TJ, T1f;
			      TA = VMUL(LDK(KP951056516), VFMA(LDK(KP618033988), Tw, Tx));
			      Ty = VMUL(LDK(KP951056516), VFNMS(LDK(KP618033988), Tx, Tw));
			      Ts = VADD(Tg, Tr);
			      Tu = VSUB(Tg, Tr);
			      TK = VSUB(TE, TH);
			      TI = VADD(TE, TH);
			      T1a = VFNMS(LDK(KP618033988), TU, TV);
			      TW = VFMA(LDK(KP618033988), TV, TU);
			      T1b = VFNMS(LDK(KP559016994), T12, T11);
			      T13 = VFMA(LDK(KP559016994), T12, T11);
			      ST(&(xo[0]), VADD(T5, Ts), ovs, &(xo[0]));
			      Tt = VFNMS(LDK(KP250000000), Ts, T5);
			      TJ = VFNMS(LDK(KP250000000), TI, TB);
			      T1f = VADD(TB, TI);
			      {
				   V T1c, T1e, T16, T14, Tv, Tz, T17, TL;
				   T1c = VMUL(LDK(KP951056516), VFNMS(LDK(KP910592997), T1b, T1a));
				   T1e = VMUL(LDK(KP951056516), VFMA(LDK(KP910592997), T1b, T1a));
				   T16 = VMUL(LDK(KP951056516), VFMA(LDK(KP910592997), T13, TW));
				   T14 = VMUL(LDK(KP951056516), VFNMS(LDK(KP910592997), T13, TW));
				   Tv = VFNMS(LDK(KP559016994), Tu, Tt);
				   Tz = VFMA(LDK(KP559016994), Tu, Tt);
				   T17 = VFNMS(LDK(KP559016994), TK, TJ);
				   TL = VFMA(LDK(KP559016994), TK, TJ);
				   ST(&(xo[WS(os, 10)]), VFMAI(T1g, T1f), ovs, &(xo[0]));
				   ST(&(xo[WS(os, 5)]), VFNMSI(T1g, T1f), ovs, &(xo[WS(os, 1)]));
				   {
					V T19, T1d, T15, TT;
					ST(&(xo[WS(os, 12)]), VFMAI(Ty, Tv), ovs, &(xo[0]));
					ST(&(xo[WS(os, 3)]), VFNMSI(Ty, Tv), ovs, &(xo[WS(os, 1)]));
					ST(&(xo[WS(os, 9)]), VFMAI(TA, Tz), ovs, &(xo[WS(os, 1)]));
					ST(&(xo[WS(os, 6)]), VFNMSI(TA, Tz), ovs, &(xo[0]));
					T19 = VFMA(LDK(KP823639103), T18, T17);
					T1d = VFNMS(LDK(KP823639103), T18, T17);
					T15 = VFNMS(LDK(KP823639103), TS, TL);
					TT = VFMA(LDK(KP823639103), TS, TL);
					ST(&(xo[WS(os, 2)]), VFMAI(T1c, T19), ovs, &(xo[0]));
					ST(&(xo[WS(os, 13)]), VFNMSI(T1c, T19), ovs, &(xo[WS(os, 1)]));
					ST(&(xo[WS(os, 7)]), VFMAI(T1e, T1d), ovs, &(xo[WS(os, 1)]));
					ST(&(xo[WS(os, 8)]), VFNMSI(T1e, T1d), ovs, &(xo[0]));
					ST(&(xo[WS(os, 4)]), VFMAI(T16, T15), ovs, &(xo[0]));
					ST(&(xo[WS(os, 11)]), VFNMSI(T16, T15), ovs, &(xo[WS(os, 1)]));
					ST(&(xo[WS(os, 14)]), VFMAI(T14, TT), ovs, &(xo[0]));
					ST(&(xo[WS(os, 1)]), VFNMSI(T14, TT), ovs, &(xo[WS(os, 1)]));
				   }
			      }
			 }
		    }
	       }
	  }
     }
}

static const kdft_desc desc = { 15, "n1fv_15", {36, 7, 42, 0}, &GENUS, 0, 0, 0, 0 };
void X(codelet_n1fv_15) (planner *p) {
     X(kdft_register) (p, n1fv_15, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_notw_c -simd -compact -variables 4 -pipeline-latency 8 -n 15 -name n1fv_15 -include n1f.h */

/*
 * This function contains 78 FP additions, 25 FP multiplications,
 * (or, 64 additions, 11 multiplications, 14 fused multiply/add),
 * 55 stack variables, and 30 memory accesses
 */
#include "n1f.h"

static void n1fv_15(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP216506350, +0.216506350946109661690930792688234045867850657);
     DVK(KP509036960, +0.509036960455127183450980863393907648510733164);
     DVK(KP823639103, +0.823639103546331925877420039278190003029660514);
     DVK(KP587785252, +0.587785252292473129168705954639072768597652438);
     DVK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DVK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DVK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DVK(KP484122918, +0.484122918275927110647408174972799951354115213);
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     const R *xi;
     R *xo;
     xi = ri;
     xo = ro;
     for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  V T5, T10, TB, TO, TU, TV, TR, Ta, Tf, Tg, Tl, Tq, Tr, TE, TH;
	  V TI, TZ, T11, T1f, T1g;
	  {
	       V T1, T2, T3, T4;
	       T1 = LD(&(xi[0]), ivs, &(xi[0]));
	       T2 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
	       T3 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
	       T4 = VADD(T2, T3);
	       T5 = VADD(T1, T4);
	       T10 = VSUB(T3, T2);
	       TB = VFNMS(LDK(KP500000000), T4, T1);
	  }
	  {
	       V T6, T9, TC, TP, Tm, Tp, TG, TN, Tb, Te, TD, TQ, Th, Tk, TF;
	       V TM, TX, TY;
	       {
		    V T7, T8, Tn, To;
		    T6 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
		    T7 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
		    T8 = LD(&(xi[WS(is, 13)]), ivs, &(xi[WS(is, 1)]));
		    T9 = VADD(T7, T8);
		    TC = VFNMS(LDK(KP500000000), T9, T6);
		    TP = VSUB(T8, T7);
		    Tm = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
		    Tn = LD(&(xi[WS(is, 14)]), ivs, &(xi[0]));
		    To = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
		    Tp = VADD(Tn, To);
		    TG = VFNMS(LDK(KP500000000), Tp, Tm);
		    TN = VSUB(To, Tn);
	       }
	       {
		    V Tc, Td, Ti, Tj;
		    Tb = LD(&(xi[WS(is, 12)]), ivs, &(xi[0]));
		    Tc = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    Td = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
		    Te = VADD(Tc, Td);
		    TD = VFNMS(LDK(KP500000000), Te, Tb);
		    TQ = VSUB(Td, Tc);
		    Th = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
		    Ti = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
		    Tj = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    Tk = VADD(Ti, Tj);
		    TF = VFNMS(LDK(KP500000000), Tk, Th);
		    TM = VSUB(Tj, Ti);
	       }
	       TO = VSUB(TM, TN);
	       TU = VSUB(TF, TG);
	       TV = VSUB(TC, TD);
	       TR = VSUB(TP, TQ);
	       Ta = VADD(T6, T9);
	       Tf = VADD(Tb, Te);
	       Tg = VADD(Ta, Tf);
	       Tl = VADD(Th, Tk);
	       Tq = VADD(Tm, Tp);
	       Tr = VADD(Tl, Tq);
	       TE = VADD(TC, TD);
	       TH = VADD(TF, TG);
	       TI = VADD(TE, TH);
	       TX = VADD(TP, TQ);
	       TY = VADD(TM, TN);
	       TZ = VMUL(LDK(KP484122918), VSUB(TX, TY));
	       T11 = VADD(TX, TY);
	  }
	  T1f = VADD(TB, TI);
	  T1g = VBYI(VMUL(LDK(KP866025403), VADD(T10, T11)));
	  ST(&(xo[WS(os, 5)]), VSUB(T1f, T1g), ovs, &(xo[WS(os, 1)]));
	  ST(&(xo[WS(os, 10)]), VADD(T1f, T1g), ovs, &(xo[0]));
	  {
	       V Tu, Ts, Tt, Ty, TA, Tw, Tx, Tz, Tv;
	       Tu = VMUL(LDK(KP559016994), VSUB(Tg, Tr));
	       Ts = VADD(Tg, Tr);
	       Tt = VFNMS(LDK(KP250000000), Ts, T5);
	       Tw = VSUB(Tl, Tq);
	       Tx = VSUB(Ta, Tf);
	       Ty = VBYI(VFNMS(LDK(KP587785252), Tx, VMUL(LDK(KP951056516), Tw)));
	       TA = VBYI(VFMA(LDK(KP951056516), Tx, VMUL(LDK(KP587785252), Tw)));
	       ST(&(xo[0]), VADD(T5, Ts), ovs, &(xo[0]));
	       Tz = VADD(Tu, Tt);
	       ST(&(xo[WS(os, 6)]), VSUB(Tz, TA), ovs, &(xo[0]));
	       ST(&(xo[WS(os, 9)]), VADD(TA, Tz), ovs, &(xo[WS(os, 1)]));
	       Tv = VSUB(Tt, Tu);
	       ST(&(xo[WS(os, 3)]), VSUB(Tv, Ty), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[WS(os, 12)]), VADD(Ty, Tv), ovs, &(xo[0]));
	  }
	  {
	       V TS, TW, T1b, T18, T13, T1a, TL, T17, T12, TJ, TK;
	       TS = VFNMS(LDK(KP509036960), TR, VMUL(LDK(KP823639103), TO));
	       TW = VFNMS(LDK(KP587785252), TV, VMUL(LDK(KP951056516), TU));
	       T1b = VFMA(LDK(KP951056516), TV, VMUL(LDK(KP587785252), TU));
	       T18 = VFMA(LDK(KP823639103), TR, VMUL(LDK(KP509036960), TO));
	       T12 = VFNMS(LDK(KP216506350), T11, VMUL(LDK(KP866025403), T10));
	       T13 = VSUB(TZ, T12);
	       T1a = VADD(TZ, T12);
	       TJ = VFNMS(LDK(KP250000000), TI, TB);
	       TK = VMUL(LDK(KP559016994), VSUB(TE, TH));
	       TL = VSUB(TJ, TK);
	       T17 = VADD(TK, TJ);
	       {
		    V TT, T14, T1d, T1e;
		    TT = VSUB(TL, TS);
		    T14 = VBYI(VSUB(TW, T13));
		    ST(&(xo[WS(os, 8)]), VSUB(TT, T14), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 7)]), VADD(TT, T14), ovs, &(xo[WS(os, 1)]));
		    T1d = VSUB(T17, T18);
		    T1e = VBYI(VADD(T1b, T1a));
		    ST(&(xo[WS(os, 11)]), VSUB(T1d, T1e), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 4)]), VADD(T1d, T1e), ovs, &(xo[0]));
	       }
	       {
		    V T15, T16, T19, T1c;
		    T15 = VADD(TL, TS);
		    T16 = VBYI(VADD(TW, T13));
		    ST(&(xo[WS(os, 13)]), VSUB(T15, T16), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 2)]), VADD(T15, T16), ovs, &(xo[0]));
		    T19 = VADD(T17, T18);
		    T1c = VBYI(VSUB(T1a, T1b));
		    ST(&(xo[WS(os, 14)]), VSUB(T19, T1c), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 1)]), VADD(T19, T1c), ovs, &(xo[WS(os, 1)]));
	       }
	  }
     }
}

static const kdft_desc desc = { 15, "n1fv_15", {64, 11, 14, 0}, &GENUS, 0, 0, 0, 0 };
void X(codelet_n1fv_15) (planner *p) {
     X(kdft_register) (p, n1fv_15, &desc);
}

#endif				/* HAVE_FMA */
