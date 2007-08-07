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
/* Generated on Mon Mar 19 07:59:37 CDT 2007 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_r2cb -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 6 -name r2cbIII_6 -dft-III -include r2cbIII.h */

/*
 * This function contains 12 FP additions, 8 FP multiplications,
 * (or, 6 additions, 2 multiplications, 6 fused multiply/add),
 * 15 stack variables, and 12 memory accesses
 */
#include "r2cbIII.h"

static void r2cbIII_6(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP1_732050807, +1.732050807568877293527446341505872366942805254);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, R0 = R0 + ovs, R1 = R1 + ovs, Cr = Cr + ivs, Ci = Ci + ivs, MAKE_VOLATILE_STRIDE(rs), MAKE_VOLATILE_STRIDE(csr), MAKE_VOLATILE_STRIDE(csi)) {
	  E T1, T8, T2, T3, T5, T6;
	  T1 = Cr[WS(csr, 1)];
	  T8 = Ci[WS(csi, 1)];
	  T2 = Cr[WS(csr, 2)];
	  T3 = Cr[0];
	  T5 = Ci[WS(csi, 2)];
	  T6 = Ci[0];
	  {
	       E T4, Ta, T7, Tc, Tb, T9;
	       T4 = T2 + T3;
	       Ta = T2 - T3;
	       T7 = T5 + T6;
	       Tc = T5 - T6;
	       Tb = FNMS(KP2_000000000, T1, T4);
	       R0[0] = KP2_000000000 * (T1 + T4);
	       T9 = FMA(KP2_000000000, T8, T7);
	       R1[WS(rs, 1)] = KP2_000000000 * (T8 - T7);
	       R0[WS(rs, 2)] = FMS(KP1_732050807, Tc, Tb);
	       R0[WS(rs, 1)] = FMA(KP1_732050807, Tc, Tb);
	       R1[WS(rs, 2)] = FMS(KP1_732050807, Ta, T9);
	       R1[0] = -(FMA(KP1_732050807, Ta, T9));
	  }
     }
}

static const kr2c_desc desc = { 6, "r2cbIII_6", {6, 2, 6, 0}, &GENUS };

void X(codelet_r2cbIII_6) (planner *p) {
     X(kr2c_register) (p, r2cbIII_6, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_r2cb -compact -variables 4 -pipeline-latency 4 -sign 1 -n 6 -name r2cbIII_6 -dft-III -include r2cbIII.h */

/*
 * This function contains 12 FP additions, 6 FP multiplications,
 * (or, 10 additions, 4 multiplications, 2 fused multiply/add),
 * 15 stack variables, and 12 memory accesses
 */
#include "r2cbIII.h"

static void r2cbIII_6(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     DK(KP1_732050807, +1.732050807568877293527446341505872366942805254);
     INT i;
     for (i = v; i > 0; i = i - 1, R0 = R0 + ovs, R1 = R1 + ovs, Cr = Cr + ivs, Ci = Ci + ivs, MAKE_VOLATILE_STRIDE(rs), MAKE_VOLATILE_STRIDE(csr), MAKE_VOLATILE_STRIDE(csi)) {
	  E T1, T6, T4, T5, T9, Tb, Ta, Tc;
	  T1 = Cr[WS(csr, 1)];
	  T6 = Ci[WS(csi, 1)];
	  {
	       E T2, T3, T7, T8;
	       T2 = Cr[WS(csr, 2)];
	       T3 = Cr[0];
	       T4 = T2 + T3;
	       T5 = KP1_732050807 * (T2 - T3);
	       T7 = Ci[WS(csi, 2)];
	       T8 = Ci[0];
	       T9 = T7 + T8;
	       Tb = KP1_732050807 * (T7 - T8);
	  }
	  R0[0] = KP2_000000000 * (T1 + T4);
	  R1[WS(rs, 1)] = KP2_000000000 * (T6 - T9);
	  Ta = FMA(KP2_000000000, T6, T9);
	  R1[0] = -(T5 + Ta);
	  R1[WS(rs, 2)] = T5 - Ta;
	  Tc = FMS(KP2_000000000, T1, T4);
	  R0[WS(rs, 1)] = Tb - Tc;
	  R0[WS(rs, 2)] = Tc + Tb;
     }
}

static const kr2c_desc desc = { 6, "r2cbIII_6", {10, 4, 2, 0}, &GENUS };

void X(codelet_r2cbIII_6) (planner *p) {
     X(kr2c_register) (p, r2cbIII_6, &desc);
}

#endif				/* HAVE_FMA */
