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
/* Generated on Mon Mar 19 07:35:31 CDT 2007 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 4 -name t1_4 -include t.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 16 additions, 6 multiplications, 6 fused multiply/add),
 * 31 stack variables, and 16 memory accesses
 */
#include "t.h"

static void t1_4(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     INT m;
     for (m = mb, W = W + (mb * 6); m < me; m = m + 1, ri = ri + ms, ii = ii + ms, W = W + 6, MAKE_VOLATILE_STRIDE(rs)) {
	  E To, Te, Tm, T8, Tw, Tx, Tq, Tk;
	  {
	       E T1, Tv, Tu, T7, Tg, Tj, Tf, Ti, Tp, Th;
	       T1 = ri[0];
	       Tv = ii[0];
	       {
		    E T3, T6, T2, T5;
		    T3 = ri[WS(rs, 2)];
		    T6 = ii[WS(rs, 2)];
		    T2 = W[2];
		    T5 = W[3];
		    {
			 E Ta, Td, Tc, Tn, Tb, Tt, T4, T9;
			 Ta = ri[WS(rs, 1)];
			 Td = ii[WS(rs, 1)];
			 Tt = T2 * T6;
			 T4 = T2 * T3;
			 T9 = W[0];
			 Tc = W[1];
			 Tu = FNMS(T5, T3, Tt);
			 T7 = FMA(T5, T6, T4);
			 Tn = T9 * Td;
			 Tb = T9 * Ta;
			 Tg = ri[WS(rs, 3)];
			 Tj = ii[WS(rs, 3)];
			 To = FNMS(Tc, Ta, Tn);
			 Te = FMA(Tc, Td, Tb);
			 Tf = W[4];
			 Ti = W[5];
		    }
	       }
	       Tm = T1 - T7;
	       T8 = T1 + T7;
	       Tw = Tu + Tv;
	       Tx = Tv - Tu;
	       Tp = Tf * Tj;
	       Th = Tf * Tg;
	       Tq = FNMS(Ti, Tg, Tp);
	       Tk = FMA(Ti, Tj, Th);
	  }
	  {
	       E Ts, Tr, Tl, Ty;
	       Ts = To + Tq;
	       Tr = To - Tq;
	       Tl = Te + Tk;
	       Ty = Te - Tk;
	       ri[WS(rs, 1)] = Tm + Tr;
	       ri[WS(rs, 3)] = Tm - Tr;
	       ii[WS(rs, 2)] = Tw - Ts;
	       ii[0] = Ts + Tw;
	       ii[WS(rs, 3)] = Ty + Tx;
	       ii[WS(rs, 1)] = Tx - Ty;
	       ri[0] = T8 + Tl;
	       ri[WS(rs, 2)] = T8 - Tl;
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 4},
     {TW_NEXT, 1, 0}
};

static const ct_desc desc = { 4, "t1_4", twinstr, &GENUS, {16, 6, 6, 0}, 0, 0, 0 };

void X(codelet_t1_4) (planner *p) {
     X(kdft_dit_register) (p, t1_4, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle -compact -variables 4 -pipeline-latency 4 -n 4 -name t1_4 -include t.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 16 additions, 6 multiplications, 6 fused multiply/add),
 * 13 stack variables, and 16 memory accesses
 */
#include "t.h"

static void t1_4(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     INT m;
     for (m = mb, W = W + (mb * 6); m < me; m = m + 1, ri = ri + ms, ii = ii + ms, W = W + 6, MAKE_VOLATILE_STRIDE(rs)) {
	  E T1, Tp, T6, To, Tc, Tk, Th, Tl;
	  T1 = ri[0];
	  Tp = ii[0];
	  {
	       E T3, T5, T2, T4;
	       T3 = ri[WS(rs, 2)];
	       T5 = ii[WS(rs, 2)];
	       T2 = W[2];
	       T4 = W[3];
	       T6 = FMA(T2, T3, T4 * T5);
	       To = FNMS(T4, T3, T2 * T5);
	  }
	  {
	       E T9, Tb, T8, Ta;
	       T9 = ri[WS(rs, 1)];
	       Tb = ii[WS(rs, 1)];
	       T8 = W[0];
	       Ta = W[1];
	       Tc = FMA(T8, T9, Ta * Tb);
	       Tk = FNMS(Ta, T9, T8 * Tb);
	  }
	  {
	       E Te, Tg, Td, Tf;
	       Te = ri[WS(rs, 3)];
	       Tg = ii[WS(rs, 3)];
	       Td = W[4];
	       Tf = W[5];
	       Th = FMA(Td, Te, Tf * Tg);
	       Tl = FNMS(Tf, Te, Td * Tg);
	  }
	  {
	       E T7, Ti, Tn, Tq;
	       T7 = T1 + T6;
	       Ti = Tc + Th;
	       ri[WS(rs, 2)] = T7 - Ti;
	       ri[0] = T7 + Ti;
	       Tn = Tk + Tl;
	       Tq = To + Tp;
	       ii[0] = Tn + Tq;
	       ii[WS(rs, 2)] = Tq - Tn;
	  }
	  {
	       E Tj, Tm, Tr, Ts;
	       Tj = T1 - T6;
	       Tm = Tk - Tl;
	       ri[WS(rs, 3)] = Tj - Tm;
	       ri[WS(rs, 1)] = Tj + Tm;
	       Tr = Tp - To;
	       Ts = Tc - Th;
	       ii[WS(rs, 1)] = Tr - Ts;
	       ii[WS(rs, 3)] = Ts + Tr;
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 4},
     {TW_NEXT, 1, 0}
};

static const ct_desc desc = { 4, "t1_4", twinstr, &GENUS, {16, 6, 6, 0}, 0, 0, 0 };

void X(codelet_t1_4) (planner *p) {
     X(kdft_dit_register) (p, t1_4, &desc);
}
#endif				/* HAVE_FMA */
