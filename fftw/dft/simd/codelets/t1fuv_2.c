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
/* Generated on Mon Mar 19 07:52:20 CDT 2007 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle_c -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -n 2 -name t1fuv_2 -include t1fu.h */

/*
 * This function contains 3 FP additions, 2 FP multiplications,
 * (or, 3 additions, 2 multiplications, 0 fused multiply/add),
 * 5 stack variables, and 4 memory accesses
 */
#include "t1fu.h"

static void t1fuv_2(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     INT m;
     R *x;
     x = ri;
     for (m = mb, W = W + (mb * ((TWVL / VL) * 2)); m < me; m = m + VL, x = x + (VL * ms), W = W + (TWVL * 2), MAKE_VOLATILE_STRIDE(rs)) {
	  V T1, T2, T3;
	  T1 = LD(&(x[0]), ms, &(x[0]));
	  T2 = LD(&(x[WS(rs, 1)]), ms, &(x[WS(rs, 1)]));
	  T3 = BYTWJ(&(W[0]), T2);
	  ST(&(x[0]), VADD(T1, T3), ms, &(x[0]));
	  ST(&(x[WS(rs, 1)]), VSUB(T1, T3), ms, &(x[WS(rs, 1)]));
     }
}

static const tw_instr twinstr[] = {
     VTW(0, 1),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 2, "t1fuv_2", twinstr, &GENUS, {3, 2, 0, 0}, 0, 0, 0 };

void X(codelet_t1fuv_2) (planner *p) {
     X(kdft_dit_register) (p, t1fuv_2, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle_c -simd -compact -variables 4 -pipeline-latency 8 -n 2 -name t1fuv_2 -include t1fu.h */

/*
 * This function contains 3 FP additions, 2 FP multiplications,
 * (or, 3 additions, 2 multiplications, 0 fused multiply/add),
 * 5 stack variables, and 4 memory accesses
 */
#include "t1fu.h"

static void t1fuv_2(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     INT m;
     R *x;
     x = ri;
     for (m = mb, W = W + (mb * ((TWVL / VL) * 2)); m < me; m = m + VL, x = x + (VL * ms), W = W + (TWVL * 2), MAKE_VOLATILE_STRIDE(rs)) {
	  V T1, T3, T2;
	  T1 = LD(&(x[0]), ms, &(x[0]));
	  T2 = LD(&(x[WS(rs, 1)]), ms, &(x[WS(rs, 1)]));
	  T3 = BYTWJ(&(W[0]), T2);
	  ST(&(x[WS(rs, 1)]), VSUB(T1, T3), ms, &(x[WS(rs, 1)]));
	  ST(&(x[0]), VADD(T1, T3), ms, &(x[0]));
     }
}

static const tw_instr twinstr[] = {
     VTW(0, 1),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 2, "t1fuv_2", twinstr, &GENUS, {3, 2, 0, 0}, 0, 0, 0 };

void X(codelet_t1fuv_2) (planner *p) {
     X(kdft_dit_register) (p, t1fuv_2, &desc);
}
#endif				/* HAVE_FMA */
