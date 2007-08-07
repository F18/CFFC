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
/* Generated on Mon Mar 19 07:46:23 CDT 2007 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_notw_c -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -sign 1 -n 8 -name n2bv_8 -with-ostride 2 -include n2b.h -store-multiple 2 */

/*
 * This function contains 26 FP additions, 10 FP multiplications,
 * (or, 16 additions, 0 multiplications, 10 fused multiply/add),
 * 38 stack variables, and 20 memory accesses
 */
#include "n2b.h"

static void n2bv_8(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     const R *xi;
     R *xo;
     xi = ii;
     xo = io;
     for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  V T1, T2, Tc, Td, T4, T5, T7, T8;
	  T1 = LD(&(xi[0]), ivs, &(xi[0]));
	  T2 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
	  Tc = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
	  Td = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
	  T4 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
	  T5 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
	  T7 = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
	  T8 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
	  {
	       V T3, Tj, Te, Tk, T6, Tm, T9, Tn, Tp, Tl;
	       T3 = VSUB(T1, T2);
	       Tj = VADD(T1, T2);
	       Te = VSUB(Tc, Td);
	       Tk = VADD(Tc, Td);
	       T6 = VSUB(T4, T5);
	       Tm = VADD(T4, T5);
	       T9 = VSUB(T7, T8);
	       Tn = VADD(T7, T8);
	       Tp = VADD(Tj, Tk);
	       Tl = VSUB(Tj, Tk);
	       {
		    V Tq, To, Ta, Tf;
		    Tq = VADD(Tm, Tn);
		    To = VSUB(Tm, Tn);
		    Ta = VADD(T6, T9);
		    Tf = VSUB(T6, T9);
		    {
			 V Tr, Ts, Tt, Tu, Tg, Ti, Tb, Th;
			 Tr = VFMAI(To, Tl);
			 STM2(&(xo[4]), Tr, ovs, &(xo[0]));
			 Ts = VFNMSI(To, Tl);
			 STM2(&(xo[12]), Ts, ovs, &(xo[0]));
			 Tt = VADD(Tp, Tq);
			 STM2(&(xo[0]), Tt, ovs, &(xo[0]));
			 Tu = VSUB(Tp, Tq);
			 STM2(&(xo[8]), Tu, ovs, &(xo[0]));
			 Tg = VFNMS(LDK(KP707106781), Tf, Te);
			 Ti = VFMA(LDK(KP707106781), Tf, Te);
			 Tb = VFNMS(LDK(KP707106781), Ta, T3);
			 Th = VFMA(LDK(KP707106781), Ta, T3);
			 {
			      V Tv, Tw, Tx, Ty;
			      Tv = VFNMSI(Ti, Th);
			      STM2(&(xo[14]), Tv, ovs, &(xo[2]));
			      STN2(&(xo[12]), Ts, Tv, ovs);
			      Tw = VFMAI(Ti, Th);
			      STM2(&(xo[2]), Tw, ovs, &(xo[2]));
			      STN2(&(xo[0]), Tt, Tw, ovs);
			      Tx = VFMAI(Tg, Tb);
			      STM2(&(xo[10]), Tx, ovs, &(xo[2]));
			      STN2(&(xo[8]), Tu, Tx, ovs);
			      Ty = VFNMSI(Tg, Tb);
			      STM2(&(xo[6]), Ty, ovs, &(xo[2]));
			      STN2(&(xo[4]), Tr, Ty, ovs);
			 }
		    }
	       }
	  }
     }
}

static const kdft_desc desc = { 8, "n2bv_8", {16, 0, 10, 0}, &GENUS, 0, 2, 0, 0 };
void X(codelet_n2bv_8) (planner *p) {
     X(kdft_register) (p, n2bv_8, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_notw_c -simd -compact -variables 4 -pipeline-latency 8 -sign 1 -n 8 -name n2bv_8 -with-ostride 2 -include n2b.h -store-multiple 2 */

/*
 * This function contains 26 FP additions, 2 FP multiplications,
 * (or, 26 additions, 2 multiplications, 0 fused multiply/add),
 * 24 stack variables, and 20 memory accesses
 */
#include "n2b.h"

static void n2bv_8(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     const R *xi;
     R *xo;
     xi = ii;
     xo = io;
     for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  V Ta, Tk, Te, Tj, T7, Tn, Tf, Tm, Tr, Tu;
	  {
	       V T8, T9, Tc, Td;
	       T8 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
	       T9 = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
	       Ta = VSUB(T8, T9);
	       Tk = VADD(T8, T9);
	       Tc = LD(&(xi[0]), ivs, &(xi[0]));
	       Td = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
	       Te = VSUB(Tc, Td);
	       Tj = VADD(Tc, Td);
	       {
		    V T1, T2, T3, T4, T5, T6;
		    T1 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    T2 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
		    T3 = VSUB(T1, T2);
		    T4 = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
		    T5 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
		    T6 = VSUB(T4, T5);
		    T7 = VMUL(LDK(KP707106781), VSUB(T3, T6));
		    Tn = VADD(T4, T5);
		    Tf = VMUL(LDK(KP707106781), VADD(T3, T6));
		    Tm = VADD(T1, T2);
	       }
	  }
	  {
	       V Ts, Tb, Tg, Tp, Tq, Tt;
	       Tb = VBYI(VSUB(T7, Ta));
	       Tg = VSUB(Te, Tf);
	       Tr = VADD(Tb, Tg);
	       STM2(&(xo[6]), Tr, ovs, &(xo[2]));
	       Ts = VSUB(Tg, Tb);
	       STM2(&(xo[10]), Ts, ovs, &(xo[2]));
	       Tp = VADD(Tj, Tk);
	       Tq = VADD(Tm, Tn);
	       Tt = VSUB(Tp, Tq);
	       STM2(&(xo[8]), Tt, ovs, &(xo[0]));
	       STN2(&(xo[8]), Tt, Ts, ovs);
	       Tu = VADD(Tp, Tq);
	       STM2(&(xo[0]), Tu, ovs, &(xo[0]));
	  }
	  {
	       V Tw, Th, Ti, Tv;
	       Th = VBYI(VADD(Ta, T7));
	       Ti = VADD(Te, Tf);
	       Tv = VADD(Th, Ti);
	       STM2(&(xo[2]), Tv, ovs, &(xo[2]));
	       STN2(&(xo[0]), Tu, Tv, ovs);
	       Tw = VSUB(Ti, Th);
	       STM2(&(xo[14]), Tw, ovs, &(xo[2]));
	       {
		    V Tl, To, Tx, Ty;
		    Tl = VSUB(Tj, Tk);
		    To = VBYI(VSUB(Tm, Tn));
		    Tx = VSUB(Tl, To);
		    STM2(&(xo[12]), Tx, ovs, &(xo[0]));
		    STN2(&(xo[12]), Tx, Tw, ovs);
		    Ty = VADD(Tl, To);
		    STM2(&(xo[4]), Ty, ovs, &(xo[0]));
		    STN2(&(xo[4]), Ty, Tr, ovs);
	       }
	  }
     }
}

static const kdft_desc desc = { 8, "n2bv_8", {26, 2, 0, 0}, &GENUS, 0, 2, 0, 0 };
void X(codelet_n2bv_8) (planner *p) {
     X(kdft_register) (p, n2bv_8, &desc);
}

#endif				/* HAVE_FMA */
