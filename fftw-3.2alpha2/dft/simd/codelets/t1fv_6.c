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
/* Generated on Mon Mar 19 07:48:20 CDT 2007 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle_c -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -n 6 -name t1fv_6 -include t1f.h */

/*
 * This function contains 23 FP additions, 18 FP multiplications,
 * (or, 17 additions, 12 multiplications, 6 fused multiply/add),
 * 27 stack variables, and 12 memory accesses
 */
#include "t1f.h"

static void t1fv_6(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
     INT m;
     R *x;
     x = ri;
     for (m = mb, W = W + (mb * ((TWVL / VL) * 10)); m < me; m = m + VL, x = x + (VL * ms), W = W + (TWVL * 10), MAKE_VOLATILE_STRIDE(rs)) {
	  V T1, T2, Ta, Tc, T5, T7;
	  T1 = LD(&(x[0]), ms, &(x[0]));
	  T2 = LD(&(x[WS(rs, 3)]), ms, &(x[WS(rs, 1)]));
	  Ta = LD(&(x[WS(rs, 4)]), ms, &(x[0]));
	  Tc = LD(&(x[WS(rs, 1)]), ms, &(x[WS(rs, 1)]));
	  T5 = LD(&(x[WS(rs, 2)]), ms, &(x[0]));
	  T7 = LD(&(x[WS(rs, 5)]), ms, &(x[WS(rs, 1)]));
	  {
	       V T3, Tb, Td, T6, T8;
	       T3 = BYTWJ(&(W[TWVL * 4]), T2);
	       Tb = BYTWJ(&(W[TWVL * 6]), Ta);
	       Td = BYTWJ(&(W[0]), Tc);
	       T6 = BYTWJ(&(W[TWVL * 2]), T5);
	       T8 = BYTWJ(&(W[TWVL * 8]), T7);
	       {
		    V Ti, T4, Tk, Te, Tj, T9;
		    Ti = VADD(T1, T3);
		    T4 = VSUB(T1, T3);
		    Tk = VADD(Tb, Td);
		    Te = VSUB(Tb, Td);
		    Tj = VADD(T6, T8);
		    T9 = VSUB(T6, T8);
		    {
			 V Tl, Tn, Tf, Th, Tm, Tg;
			 Tl = VADD(Tj, Tk);
			 Tn = VMUL(LDK(KP866025403), VSUB(Tk, Tj));
			 Tf = VADD(T9, Te);
			 Th = VMUL(LDK(KP866025403), VSUB(Te, T9));
			 ST(&(x[0]), VADD(Ti, Tl), ms, &(x[0]));
			 Tm = VFNMS(LDK(KP500000000), Tl, Ti);
			 ST(&(x[WS(rs, 3)]), VADD(T4, Tf), ms, &(x[WS(rs, 1)]));
			 Tg = VFNMS(LDK(KP500000000), Tf, T4);
			 ST(&(x[WS(rs, 2)]), VFNMSI(Tn, Tm), ms, &(x[0]));
			 ST(&(x[WS(rs, 4)]), VFMAI(Tn, Tm), ms, &(x[0]));
			 ST(&(x[WS(rs, 5)]), VFNMSI(Th, Tg), ms, &(x[WS(rs, 1)]));
			 ST(&(x[WS(rs, 1)]), VFMAI(Th, Tg), ms, &(x[WS(rs, 1)]));
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     VTW(0, 1),
     VTW(0, 2),
     VTW(0, 3),
     VTW(0, 4),
     VTW(0, 5),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 6, "t1fv_6", twinstr, &GENUS, {17, 12, 6, 0}, 0, 0, 0 };

void X(codelet_t1fv_6) (planner *p) {
     X(kdft_dit_register) (p, t1fv_6, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle_c -simd -compact -variables 4 -pipeline-latency 8 -n 6 -name t1fv_6 -include t1f.h */

/*
 * This function contains 23 FP additions, 14 FP multiplications,
 * (or, 21 additions, 12 multiplications, 2 fused multiply/add),
 * 19 stack variables, and 12 memory accesses
 */
#include "t1f.h"

static void t1fv_6(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
     INT m;
     R *x;
     x = ri;
     for (m = mb, W = W + (mb * ((TWVL / VL) * 10)); m < me; m = m + VL, x = x + (VL * ms), W = W + (TWVL * 10), MAKE_VOLATILE_STRIDE(rs)) {
	  V T4, Ti, Te, Tk, T9, Tj, T1, T3, T2;
	  T1 = LD(&(x[0]), ms, &(x[0]));
	  T2 = LD(&(x[WS(rs, 3)]), ms, &(x[WS(rs, 1)]));
	  T3 = BYTWJ(&(W[TWVL * 4]), T2);
	  T4 = VSUB(T1, T3);
	  Ti = VADD(T1, T3);
	  {
	       V Tb, Td, Ta, Tc;
	       Ta = LD(&(x[WS(rs, 4)]), ms, &(x[0]));
	       Tb = BYTWJ(&(W[TWVL * 6]), Ta);
	       Tc = LD(&(x[WS(rs, 1)]), ms, &(x[WS(rs, 1)]));
	       Td = BYTWJ(&(W[0]), Tc);
	       Te = VSUB(Tb, Td);
	       Tk = VADD(Tb, Td);
	  }
	  {
	       V T6, T8, T5, T7;
	       T5 = LD(&(x[WS(rs, 2)]), ms, &(x[0]));
	       T6 = BYTWJ(&(W[TWVL * 2]), T5);
	       T7 = LD(&(x[WS(rs, 5)]), ms, &(x[WS(rs, 1)]));
	       T8 = BYTWJ(&(W[TWVL * 8]), T7);
	       T9 = VSUB(T6, T8);
	       Tj = VADD(T6, T8);
	  }
	  {
	       V Th, Tf, Tg, Tn, Tl, Tm;
	       Th = VBYI(VMUL(LDK(KP866025403), VSUB(Te, T9)));
	       Tf = VADD(T9, Te);
	       Tg = VFNMS(LDK(KP500000000), Tf, T4);
	       ST(&(x[WS(rs, 3)]), VADD(T4, Tf), ms, &(x[WS(rs, 1)]));
	       ST(&(x[WS(rs, 1)]), VADD(Tg, Th), ms, &(x[WS(rs, 1)]));
	       ST(&(x[WS(rs, 5)]), VSUB(Tg, Th), ms, &(x[WS(rs, 1)]));
	       Tn = VBYI(VMUL(LDK(KP866025403), VSUB(Tk, Tj)));
	       Tl = VADD(Tj, Tk);
	       Tm = VFNMS(LDK(KP500000000), Tl, Ti);
	       ST(&(x[0]), VADD(Ti, Tl), ms, &(x[0]));
	       ST(&(x[WS(rs, 4)]), VADD(Tm, Tn), ms, &(x[0]));
	       ST(&(x[WS(rs, 2)]), VSUB(Tm, Tn), ms, &(x[0]));
	  }
     }
}

static const tw_instr twinstr[] = {
     VTW(0, 1),
     VTW(0, 2),
     VTW(0, 3),
     VTW(0, 4),
     VTW(0, 5),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 6, "t1fv_6", twinstr, &GENUS, {21, 12, 2, 0}, 0, 0, 0 };

void X(codelet_t1fv_6) (planner *p) {
     X(kdft_dit_register) (p, t1fv_6, &desc);
}
#endif				/* HAVE_FMA */
