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
/* Generated on Mon Mar 19 07:44:17 CDT 2007 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_notw_c -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -sign 1 -n 4 -name n1bv_4 -include n1b.h */

/*
 * This function contains 8 FP additions, 2 FP multiplications,
 * (or, 6 additions, 0 multiplications, 2 fused multiply/add),
 * 11 stack variables, and 8 memory accesses
 */
#include "n1b.h"

static void n1bv_4(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     INT i;
     const R *xi;
     R *xo;
     xi = ii;
     xo = io;
     for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  V T1, T2, T4, T5;
	  T1 = LD(&(xi[0]), ivs, &(xi[0]));
	  T2 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
	  T4 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
	  T5 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
	  {
	       V T3, T7, T6, T8;
	       T3 = VSUB(T1, T2);
	       T7 = VADD(T1, T2);
	       T6 = VSUB(T4, T5);
	       T8 = VADD(T4, T5);
	       ST(&(xo[WS(os, 2)]), VSUB(T7, T8), ovs, &(xo[0]));
	       ST(&(xo[0]), VADD(T7, T8), ovs, &(xo[0]));
	       ST(&(xo[WS(os, 1)]), VFMAI(T6, T3), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[WS(os, 3)]), VFNMSI(T6, T3), ovs, &(xo[WS(os, 1)]));
	  }
     }
}

static const kdft_desc desc = { 4, "n1bv_4", {6, 0, 2, 0}, &GENUS, 0, 0, 0, 0 };
void X(codelet_n1bv_4) (planner *p) {
     X(kdft_register) (p, n1bv_4, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_notw_c -simd -compact -variables 4 -pipeline-latency 8 -sign 1 -n 4 -name n1bv_4 -include n1b.h */

/*
 * This function contains 8 FP additions, 0 FP multiplications,
 * (or, 8 additions, 0 multiplications, 0 fused multiply/add),
 * 11 stack variables, and 8 memory accesses
 */
#include "n1b.h"

static void n1bv_4(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     INT i;
     const R *xi;
     R *xo;
     xi = ii;
     xo = io;
     for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  V T3, T7, T6, T8;
	  {
	       V T1, T2, T4, T5;
	       T1 = LD(&(xi[0]), ivs, &(xi[0]));
	       T2 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
	       T3 = VSUB(T1, T2);
	       T7 = VADD(T1, T2);
	       T4 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
	       T5 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
	       T6 = VBYI(VSUB(T4, T5));
	       T8 = VADD(T4, T5);
	  }
	  ST(&(xo[WS(os, 3)]), VSUB(T3, T6), ovs, &(xo[WS(os, 1)]));
	  ST(&(xo[0]), VADD(T7, T8), ovs, &(xo[0]));
	  ST(&(xo[WS(os, 1)]), VADD(T3, T6), ovs, &(xo[WS(os, 1)]));
	  ST(&(xo[WS(os, 2)]), VSUB(T7, T8), ovs, &(xo[0]));
     }
}

static const kdft_desc desc = { 4, "n1bv_4", {8, 0, 0, 0}, &GENUS, 0, 0, 0, 0 };
void X(codelet_n1bv_4) (planner *p) {
     X(kdft_register) (p, n1bv_4, &desc);
}

#endif				/* HAVE_FMA */
