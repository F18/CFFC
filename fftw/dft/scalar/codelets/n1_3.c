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
/* Generated on Mon Mar 19 07:34:40 CDT 2007 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_notw -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 3 -name n1_3 -include n.h */

/*
 * This function contains 12 FP additions, 6 FP multiplications,
 * (or, 6 additions, 0 multiplications, 6 fused multiply/add),
 * 15 stack variables, and 12 memory accesses
 */
#include "n.h"

static void n1_3(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  E T1, T9, T2, T3, T6, T7;
	  T1 = ri[0];
	  T9 = ii[0];
	  T2 = ri[WS(is, 1)];
	  T3 = ri[WS(is, 2)];
	  T6 = ii[WS(is, 1)];
	  T7 = ii[WS(is, 2)];
	  {
	       E T4, Tc, T8, Ta, T5, Tb;
	       T4 = T2 + T3;
	       Tc = T3 - T2;
	       T8 = T6 - T7;
	       Ta = T6 + T7;
	       T5 = FNMS(KP500000000, T4, T1);
	       ro[0] = T1 + T4;
	       Tb = FNMS(KP500000000, Ta, T9);
	       io[0] = T9 + Ta;
	       ro[WS(os, 1)] = FMA(KP866025403, T8, T5);
	       ro[WS(os, 2)] = FNMS(KP866025403, T8, T5);
	       io[WS(os, 2)] = FNMS(KP866025403, Tc, Tb);
	       io[WS(os, 1)] = FMA(KP866025403, Tc, Tb);
	  }
     }
}

static const kdft_desc desc = { 3, "n1_3", {6, 0, 6, 0}, &GENUS, 0, 0, 0, 0 };
void X(codelet_n1_3) (planner *p) {
     X(kdft_register) (p, n1_3, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_notw -compact -variables 4 -pipeline-latency 4 -n 3 -name n1_3 -include n.h */

/*
 * This function contains 12 FP additions, 4 FP multiplications,
 * (or, 10 additions, 2 multiplications, 2 fused multiply/add),
 * 15 stack variables, and 12 memory accesses
 */
#include "n.h"

static void n1_3(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  E T1, Ta, T4, T9, T8, Tb, T5, Tc;
	  T1 = ri[0];
	  Ta = ii[0];
	  {
	       E T2, T3, T6, T7;
	       T2 = ri[WS(is, 1)];
	       T3 = ri[WS(is, 2)];
	       T4 = T2 + T3;
	       T9 = KP866025403 * (T3 - T2);
	       T6 = ii[WS(is, 1)];
	       T7 = ii[WS(is, 2)];
	       T8 = KP866025403 * (T6 - T7);
	       Tb = T6 + T7;
	  }
	  ro[0] = T1 + T4;
	  io[0] = Ta + Tb;
	  T5 = FNMS(KP500000000, T4, T1);
	  ro[WS(os, 2)] = T5 - T8;
	  ro[WS(os, 1)] = T5 + T8;
	  Tc = FNMS(KP500000000, Tb, Ta);
	  io[WS(os, 1)] = T9 + Tc;
	  io[WS(os, 2)] = Tc - T9;
     }
}

static const kdft_desc desc = { 3, "n1_3", {10, 2, 2, 0}, &GENUS, 0, 0, 0, 0 };
void X(codelet_n1_3) (planner *p) {
     X(kdft_register) (p, n1_3, &desc);
}

#endif				/* HAVE_FMA */
