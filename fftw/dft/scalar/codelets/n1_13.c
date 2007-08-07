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
/* Generated on Mon Mar 19 07:34:43 CDT 2007 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_notw -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 13 -name n1_13 -include n.h */

/*
 * This function contains 176 FP additions, 114 FP multiplications,
 * (or, 62 additions, 0 multiplications, 114 fused multiply/add),
 * 87 stack variables, and 52 memory accesses
 */
#include "n.h"

static void n1_13(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP875502302, +0.875502302409147941146295545768755143177842006);
     DK(KP520028571, +0.520028571888864619117130500499232802493238139);
     DK(KP575140729, +0.575140729474003121368385547455453388461001608);
     DK(KP600477271, +0.600477271932665282925769253334763009352012849);
     DK(KP300462606, +0.300462606288665774426601772289207995520941381);
     DK(KP516520780, +0.516520780623489722840901288569017135705033622);
     DK(KP968287244, +0.968287244361984016049539446938120421179794516);
     DK(KP503537032, +0.503537032863766627246873853868466977093348562);
     DK(KP251768516, +0.251768516431883313623436926934233488546674281);
     DK(KP581704778, +0.581704778510515730456870384989698884939833902);
     DK(KP859542535, +0.859542535098774820163672132761689612766401925);
     DK(KP083333333, +0.083333333333333333333333333333333333333333333);
     DK(KP957805992, +0.957805992594665126462521754605754580515587217);
     DK(KP522026385, +0.522026385161275033714027226654165028300441940);
     DK(KP853480001, +0.853480001859823990758994934970528322872359049);
     DK(KP769338817, +0.769338817572980603471413688209101117038278899);
     DK(KP612264650, +0.612264650376756543746494474777125408779395514);
     DK(KP038632954, +0.038632954644348171955506895830342264440241080);
     DK(KP302775637, +0.302775637731994646559610633735247973125648287);
     DK(KP514918778, +0.514918778086315755491789696138117261566051239);
     DK(KP686558370, +0.686558370781754340655719594850823015421401653);
     DK(KP226109445, +0.226109445035782405468510155372505010481906348);
     DK(KP301479260, +0.301479260047709873958013540496673347309208464);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  E T2B, T2H, T2I, T2G;
	  {
	       E T1, T1P, T2n, T2o, To, TH, T2h, T2k, TE, TB, TF, Tw, T2j, T2c, T1m;
	       E T1W, T1X, T1c, T19, T1j, T12, T1f, T21, T24, T27, T1U;
	       T1 = ri[0];
	       T1P = ii[0];
	       {
		    E T2b, Tv, Ts, T2a;
		    {
			 E T2d, Tf, Tq, Ty, Tb, Tr, T6, Tx, Ti, Tt, Tu, Tl;
			 {
			      E T7, T8, T9, Td, Te;
			      Td = ri[WS(is, 8)];
			      Te = ri[WS(is, 5)];
			      T7 = ri[WS(is, 12)];
			      T8 = ri[WS(is, 10)];
			      T9 = ri[WS(is, 4)];
			      T2d = Td - Te;
			      Tf = Td + Te;
			      {
				   E T2, Ta, T3, T4;
				   T2 = ri[WS(is, 1)];
				   Ta = T8 + T9;
				   Tq = T8 - T9;
				   T3 = ri[WS(is, 3)];
				   T4 = ri[WS(is, 9)];
				   {
					E Tg, T5, Th, Tj, Tk;
					Tg = ri[WS(is, 11)];
					Ty = FMS(KP500000000, Ta, T7);
					Tb = T7 + Ta;
					Tr = T4 - T3;
					T5 = T3 + T4;
					Th = ri[WS(is, 6)];
					Tj = ri[WS(is, 7)];
					Tk = ri[WS(is, 2)];
					T6 = T2 + T5;
					Tx = FNMS(KP500000000, T5, T2);
					Ti = Tg + Th;
					Tt = Tg - Th;
					Tu = Tj - Tk;
					Tl = Tj + Tk;
				   }
			      }
			 }
			 {
			      E Tc, Tm, T2e, T2g;
			      Tc = T6 + Tb;
			      T2n = T6 - Tb;
			      T2b = Ti - Tl;
			      Tm = Ti + Tl;
			      T2e = Tt + Tu;
			      Tv = Tt - Tu;
			      Ts = Tq - Tr;
			      T2g = Tr + Tq;
			      {
				   E Tz, TA, Tn, T2f;
				   Tz = Tx - Ty;
				   T2a = Tx + Ty;
				   TA = FNMS(KP500000000, Tm, Tf);
				   Tn = Tf + Tm;
				   T2f = FNMS(KP500000000, T2e, T2d);
				   T2o = T2d + T2e;
				   To = Tc + Tn;
				   TH = Tc - Tn;
				   T2h = FMA(KP866025403, T2g, T2f);
				   T2k = FNMS(KP866025403, T2g, T2f);
				   TE = Tz - TA;
				   TB = Tz + TA;
			      }
			 }
		    }
		    {
			 E T1R, TM, T10, T18, T1l, TX, T1k, T15, TP, T1a, T1b, TS;
			 {
			      E T16, TY, TZ, TK, TL;
			      TK = ii[WS(is, 8)];
			      TF = Ts - Tv;
			      Tw = Ts + Tv;
			      T2j = FNMS(KP866025403, T2b, T2a);
			      T2c = FMA(KP866025403, T2b, T2a);
			      TL = ii[WS(is, 5)];
			      T16 = ii[WS(is, 12)];
			      TY = ii[WS(is, 10)];
			      TZ = ii[WS(is, 4)];
			      T1R = TK + TL;
			      TM = TK - TL;
			      {
				   E T13, T17, TV, TW;
				   T13 = ii[WS(is, 1)];
				   T17 = TY + TZ;
				   T10 = TY - TZ;
				   TV = ii[WS(is, 9)];
				   TW = ii[WS(is, 3)];
				   {
					E TN, T14, TO, TQ, TR;
					TN = ii[WS(is, 11)];
					T18 = FMS(KP500000000, T17, T16);
					T1l = T16 + T17;
					TX = TV - TW;
					T14 = TW + TV;
					TO = ii[WS(is, 6)];
					TQ = ii[WS(is, 7)];
					TR = ii[WS(is, 2)];
					T1k = T13 + T14;
					T15 = FNMS(KP500000000, T14, T13);
					TP = TN - TO;
					T1a = TN + TO;
					T1b = TQ + TR;
					TS = TQ - TR;
				   }
			      }
			 }
			 {
			      E T1Q, T11, TT, T1S;
			      T1Q = T1k + T1l;
			      T1m = T1k - T1l;
			      T11 = TX + T10;
			      T1W = T10 - TX;
			      T1X = TP - TS;
			      TT = TP + TS;
			      T1S = T1a + T1b;
			      T1c = T1a - T1b;
			      {
				   E T1Z, TU, T1T, T20;
				   T19 = T15 + T18;
				   T1Z = T15 - T18;
				   T1j = TM + TT;
				   TU = FNMS(KP500000000, TT, TM);
				   T1T = T1R + T1S;
				   T20 = FNMS(KP500000000, T1S, T1R);
				   T12 = FMA(KP866025403, T11, TU);
				   T1f = FNMS(KP866025403, T11, TU);
				   T21 = T1Z + T20;
				   T24 = T1Z - T20;
				   T27 = T1Q - T1T;
				   T1U = T1Q + T1T;
			      }
			 }
		    }
	       }
	       {
		    E T1g, T1d, T25, T1Y;
		    ro[0] = T1 + To;
		    T1g = FNMS(KP866025403, T1c, T19);
		    T1d = FMA(KP866025403, T1c, T19);
		    T25 = T1W - T1X;
		    T1Y = T1W + T1X;
		    io[0] = T1P + T1U;
		    {
			 E T1C, T1B, T1F, T1K;
			 {
			      E TC, T1J, T1z, T1w, T1I, T1O, Tp, T1E, T1q, TI, T1o, T1s;
			      {
				   E TG, T1n, T1G, T1u, T1e, T1h, T1v, T1x, T1y, T1H, T1i;
				   TC = FMA(KP301479260, TB, Tw);
				   T1x = FNMS(KP226109445, Tw, TB);
				   T1y = FMA(KP686558370, TE, TF);
				   TG = FNMS(KP514918778, TF, TE);
				   T1n = FNMS(KP302775637, T1m, T1j);
				   T1G = FMA(KP302775637, T1j, T1m);
				   T1u = FNMS(KP038632954, T12, T1d);
				   T1e = FMA(KP038632954, T1d, T12);
				   T1h = FMA(KP612264650, T1g, T1f);
				   T1v = FNMS(KP612264650, T1f, T1g);
				   T1J = FMA(KP769338817, T1y, T1x);
				   T1z = FNMS(KP769338817, T1y, T1x);
				   T1H = FNMS(KP853480001, T1v, T1u);
				   T1w = FMA(KP853480001, T1v, T1u);
				   T1I = FNMS(KP522026385, T1H, T1G);
				   T1O = FMA(KP957805992, T1G, T1H);
				   Tp = FNMS(KP083333333, To, T1);
				   T1E = FMA(KP853480001, T1h, T1e);
				   T1i = FNMS(KP853480001, T1h, T1e);
				   T1q = FNMS(KP859542535, TG, TH);
				   TI = FMA(KP581704778, TH, TG);
				   T1o = FMA(KP957805992, T1n, T1i);
				   T1s = FNMS(KP522026385, T1i, T1n);
			      }
			      {
				   E T1A, T1D, T1t, T1L, T1M;
				   {
					E T1p, TD, TJ, T1N, T1r;
					T1p = FNMS(KP251768516, TC, Tp);
					TD = FMA(KP503537032, TC, Tp);
					T1C = FNMS(KP968287244, T1z, T1w);
					T1A = FMA(KP968287244, T1z, T1w);
					TJ = FMA(KP516520780, TI, TD);
					T1N = FNMS(KP516520780, TI, TD);
					T1D = FNMS(KP300462606, T1q, T1p);
					T1r = FMA(KP300462606, T1q, T1p);
					ro[WS(os, 8)] = FNMS(KP600477271, T1O, T1N);
					ro[WS(os, 12)] = FMA(KP600477271, T1o, TJ);
					ro[WS(os, 1)] = FNMS(KP600477271, T1o, TJ);
					T1t = FNMS(KP575140729, T1s, T1r);
					T1B = FMA(KP575140729, T1s, T1r);
					ro[WS(os, 5)] = FMA(KP600477271, T1O, T1N);
				   }
				   T1L = FNMS(KP520028571, T1E, T1D);
				   T1F = FMA(KP520028571, T1E, T1D);
				   T1K = FMA(KP875502302, T1J, T1I);
				   T1M = FNMS(KP875502302, T1J, T1I);
				   ro[WS(os, 3)] = FMA(KP520028571, T1A, T1t);
				   ro[WS(os, 9)] = FNMS(KP520028571, T1A, T1t);
				   ro[WS(os, 6)] = FMA(KP575140729, T1M, T1L);
				   ro[WS(os, 11)] = FNMS(KP575140729, T1M, T1L);
			      }
			 }
			 {
			      E T22, T2F, T2N, T2K, T2w, T2A, T1V, T2C, T28, T2y, T2M, T2q;
			      {
				   E T26, T2v, T2p, T2i, T2s, T2t, T2l, T2D, T2E, T2u, T2m;
				   T2D = FNMS(KP226109445, T1Y, T21);
				   T22 = FMA(KP301479260, T21, T1Y);
				   ro[WS(os, 2)] = FMA(KP575140729, T1K, T1F);
				   ro[WS(os, 7)] = FNMS(KP575140729, T1K, T1F);
				   ro[WS(os, 4)] = FMA(KP520028571, T1C, T1B);
				   ro[WS(os, 10)] = FNMS(KP520028571, T1C, T1B);
				   T26 = FNMS(KP514918778, T25, T24);
				   T2E = FMA(KP686558370, T24, T25);
				   T2v = FNMS(KP302775637, T2n, T2o);
				   T2p = FMA(KP302775637, T2o, T2n);
				   T2i = FNMS(KP038632954, T2h, T2c);
				   T2s = FMA(KP038632954, T2c, T2h);
				   T2t = FMA(KP612264650, T2j, T2k);
				   T2l = FNMS(KP612264650, T2k, T2j);
				   T2F = FNMS(KP769338817, T2E, T2D);
				   T2N = FMA(KP769338817, T2E, T2D);
				   T2K = FMA(KP853480001, T2t, T2s);
				   T2u = FNMS(KP853480001, T2t, T2s);
				   T2w = FMA(KP957805992, T2v, T2u);
				   T2A = FNMS(KP522026385, T2u, T2v);
				   T1V = FNMS(KP083333333, T1U, T1P);
				   T2m = FNMS(KP853480001, T2l, T2i);
				   T2C = FMA(KP853480001, T2l, T2i);
				   T28 = FMA(KP581704778, T27, T26);
				   T2y = FNMS(KP859542535, T26, T27);
				   T2M = FNMS(KP522026385, T2m, T2p);
				   T2q = FMA(KP957805992, T2p, T2m);
			      }
			      {
				   E T2O, T2Q, T2z, T2P, T2L;
				   {
					E T23, T2x, T2r, T29, T2J;
					T23 = FMA(KP503537032, T22, T1V);
					T2x = FNMS(KP251768516, T22, T1V);
					T2O = FNMS(KP875502302, T2N, T2M);
					T2Q = FMA(KP875502302, T2N, T2M);
					T2r = FMA(KP516520780, T28, T23);
					T29 = FNMS(KP516520780, T28, T23);
					T2z = FMA(KP300462606, T2y, T2x);
					T2J = FNMS(KP300462606, T2y, T2x);
					io[WS(os, 12)] = FNMS(KP600477271, T2w, T2r);
					io[WS(os, 1)] = FMA(KP600477271, T2w, T2r);
					io[WS(os, 8)] = FMA(KP600477271, T2q, T29);
					io[WS(os, 5)] = FNMS(KP600477271, T2q, T29);
					T2P = FMA(KP520028571, T2K, T2J);
					T2L = FNMS(KP520028571, T2K, T2J);
				   }
				   T2B = FMA(KP575140729, T2A, T2z);
				   T2H = FNMS(KP575140729, T2A, T2z);
				   io[WS(os, 11)] = FMA(KP575140729, T2Q, T2P);
				   io[WS(os, 6)] = FNMS(KP575140729, T2Q, T2P);
				   io[WS(os, 7)] = FMA(KP575140729, T2O, T2L);
				   io[WS(os, 2)] = FNMS(KP575140729, T2O, T2L);
				   T2I = FMA(KP968287244, T2F, T2C);
				   T2G = FNMS(KP968287244, T2F, T2C);
			      }
			 }
		    }
	       }
	  }
	  io[WS(os, 10)] = FMA(KP520028571, T2I, T2H);
	  io[WS(os, 4)] = FNMS(KP520028571, T2I, T2H);
	  io[WS(os, 9)] = FMA(KP520028571, T2G, T2B);
	  io[WS(os, 3)] = FNMS(KP520028571, T2G, T2B);
     }
}

static const kdft_desc desc = { 13, "n1_13", {62, 0, 114, 0}, &GENUS, 0, 0, 0, 0 };
void X(codelet_n1_13) (planner *p) {
     X(kdft_register) (p, n1_13, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_notw -compact -variables 4 -pipeline-latency 4 -n 13 -name n1_13 -include n.h */

/*
 * This function contains 176 FP additions, 68 FP multiplications,
 * (or, 138 additions, 30 multiplications, 38 fused multiply/add),
 * 71 stack variables, and 52 memory accesses
 */
#include "n.h"

static void n1_13(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     DK(KP083333333, +0.083333333333333333333333333333333333333333333);
     DK(KP251768516, +0.251768516431883313623436926934233488546674281);
     DK(KP075902986, +0.075902986037193865983102897245103540356428373);
     DK(KP132983124, +0.132983124607418643793760531921092974399165133);
     DK(KP258260390, +0.258260390311744861420450644284508567852516811);
     DK(KP1_732050807, +1.732050807568877293527446341505872366942805254);
     DK(KP300238635, +0.300238635966332641462884626667381504676006424);
     DK(KP011599105, +0.011599105605768290721655456654083252189827041);
     DK(KP156891391, +0.156891391051584611046832726756003269660212636);
     DK(KP256247671, +0.256247671582936600958684654061725059144125175);
     DK(KP174138601, +0.174138601152135905005660794929264742616964676);
     DK(KP575140729, +0.575140729474003121368385547455453388461001608);
     DK(KP503537032, +0.503537032863766627246873853868466977093348562);
     DK(KP113854479, +0.113854479055790798974654345867655310534642560);
     DK(KP265966249, +0.265966249214837287587521063842185948798330267);
     DK(KP387390585, +0.387390585467617292130675966426762851778775217);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP300462606, +0.300462606288665774426601772289207995520941381);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  E T1, T1q, Tt, Tu, To, T22, T20, T24, TF, TH, TA, TI, T1X, T25, T2a;
	  E T2d, T18, T1n, T2k, T2n, T1l, T1r, T1f, T1o, T2h, T2m;
	  T1 = ri[0];
	  T1q = ii[0];
	  {
	       E Tf, Tp, Tb, TC, Tx, T6, TB, Tw, Ti, Tq, Tl, Tr, Tm, Ts, Td;
	       E Te, Tc, Tn;
	       Td = ri[WS(is, 8)];
	       Te = ri[WS(is, 5)];
	       Tf = Td + Te;
	       Tp = Td - Te;
	       {
		    E T7, T8, T9, Ta;
		    T7 = ri[WS(is, 12)];
		    T8 = ri[WS(is, 10)];
		    T9 = ri[WS(is, 4)];
		    Ta = T8 + T9;
		    Tb = T7 + Ta;
		    TC = T8 - T9;
		    Tx = FNMS(KP500000000, Ta, T7);
	       }
	       {
		    E T2, T3, T4, T5;
		    T2 = ri[WS(is, 1)];
		    T3 = ri[WS(is, 3)];
		    T4 = ri[WS(is, 9)];
		    T5 = T3 + T4;
		    T6 = T2 + T5;
		    TB = T3 - T4;
		    Tw = FNMS(KP500000000, T5, T2);
	       }
	       {
		    E Tg, Th, Tj, Tk;
		    Tg = ri[WS(is, 11)];
		    Th = ri[WS(is, 6)];
		    Ti = Tg + Th;
		    Tq = Tg - Th;
		    Tj = ri[WS(is, 7)];
		    Tk = ri[WS(is, 2)];
		    Tl = Tj + Tk;
		    Tr = Tj - Tk;
	       }
	       Tm = Ti + Tl;
	       Ts = Tq + Tr;
	       Tt = Tp + Ts;
	       Tu = T6 - Tb;
	       Tc = T6 + Tb;
	       Tn = Tf + Tm;
	       To = Tc + Tn;
	       T22 = KP300462606 * (Tc - Tn);
	       {
		    E T1Y, T1Z, TD, TE;
		    T1Y = TB + TC;
		    T1Z = Tq - Tr;
		    T20 = T1Y - T1Z;
		    T24 = T1Y + T1Z;
		    TD = KP866025403 * (TB - TC);
		    TE = FNMS(KP500000000, Ts, Tp);
		    TF = TD - TE;
		    TH = TD + TE;
	       }
	       {
		    E Ty, Tz, T1V, T1W;
		    Ty = Tw - Tx;
		    Tz = KP866025403 * (Ti - Tl);
		    TA = Ty + Tz;
		    TI = Ty - Tz;
		    T1V = Tw + Tx;
		    T1W = FNMS(KP500000000, Tm, Tf);
		    T1X = T1V - T1W;
		    T25 = T1V + T1W;
	       }
	  }
	  {
	       E TZ, T2b, TV, T1i, T1a, TQ, T1h, T19, T12, T1d, T15, T1c, T16, T2c, TX;
	       E TY, TW, T17;
	       TX = ii[WS(is, 8)];
	       TY = ii[WS(is, 5)];
	       TZ = TX + TY;
	       T2b = TX - TY;
	       {
		    E TR, TS, TT, TU;
		    TR = ii[WS(is, 12)];
		    TS = ii[WS(is, 10)];
		    TT = ii[WS(is, 4)];
		    TU = TS + TT;
		    TV = FNMS(KP500000000, TU, TR);
		    T1i = TR + TU;
		    T1a = TS - TT;
	       }
	       {
		    E TM, TN, TO, TP;
		    TM = ii[WS(is, 1)];
		    TN = ii[WS(is, 3)];
		    TO = ii[WS(is, 9)];
		    TP = TN + TO;
		    TQ = FNMS(KP500000000, TP, TM);
		    T1h = TM + TP;
		    T19 = TN - TO;
	       }
	       {
		    E T10, T11, T13, T14;
		    T10 = ii[WS(is, 11)];
		    T11 = ii[WS(is, 6)];
		    T12 = T10 + T11;
		    T1d = T10 - T11;
		    T13 = ii[WS(is, 7)];
		    T14 = ii[WS(is, 2)];
		    T15 = T13 + T14;
		    T1c = T13 - T14;
	       }
	       T16 = T12 + T15;
	       T2c = T1d + T1c;
	       T2a = T1h - T1i;
	       T2d = T2b + T2c;
	       TW = TQ + TV;
	       T17 = FNMS(KP500000000, T16, TZ);
	       T18 = TW - T17;
	       T1n = TW + T17;
	       {
		    E T2i, T2j, T1j, T1k;
		    T2i = TQ - TV;
		    T2j = KP866025403 * (T15 - T12);
		    T2k = T2i + T2j;
		    T2n = T2i - T2j;
		    T1j = T1h + T1i;
		    T1k = TZ + T16;
		    T1l = KP300462606 * (T1j - T1k);
		    T1r = T1j + T1k;
	       }
	       {
		    E T1b, T1e, T2f, T2g;
		    T1b = T19 + T1a;
		    T1e = T1c - T1d;
		    T1f = T1b + T1e;
		    T1o = T1e - T1b;
		    T2f = FNMS(KP500000000, T2c, T2b);
		    T2g = KP866025403 * (T1a - T19);
		    T2h = T2f - T2g;
		    T2m = T2g + T2f;
	       }
	  }
	  ro[0] = T1 + To;
	  io[0] = T1q + T1r;
	  {
	       E T1D, T1N, T1y, T1x, T1E, T1O, Tv, TK, T1J, T1Q, T1m, T1R, T1t, T1I, TG;
	       E TJ;
	       {
		    E T1B, T1C, T1v, T1w;
		    T1B = FMA(KP387390585, T1f, KP265966249 * T18);
		    T1C = FMA(KP113854479, T1o, KP503537032 * T1n);
		    T1D = T1B + T1C;
		    T1N = T1C - T1B;
		    T1y = FMA(KP575140729, Tu, KP174138601 * Tt);
		    T1v = FNMS(KP156891391, TH, KP256247671 * TI);
		    T1w = FMA(KP011599105, TF, KP300238635 * TA);
		    T1x = T1v - T1w;
		    T1E = T1y + T1x;
		    T1O = KP1_732050807 * (T1v + T1w);
	       }
	       Tv = FNMS(KP174138601, Tu, KP575140729 * Tt);
	       TG = FNMS(KP300238635, TF, KP011599105 * TA);
	       TJ = FMA(KP256247671, TH, KP156891391 * TI);
	       TK = TG - TJ;
	       T1J = KP1_732050807 * (TJ + TG);
	       T1Q = Tv - TK;
	       {
		    E T1g, T1H, T1p, T1s, T1G;
		    T1g = FNMS(KP132983124, T1f, KP258260390 * T18);
		    T1H = T1l - T1g;
		    T1p = FNMS(KP251768516, T1o, KP075902986 * T1n);
		    T1s = FNMS(KP083333333, T1r, T1q);
		    T1G = T1s - T1p;
		    T1m = FMA(KP2_000000000, T1g, T1l);
		    T1R = T1H + T1G;
		    T1t = FMA(KP2_000000000, T1p, T1s);
		    T1I = T1G - T1H;
	       }
	       {
		    E TL, T1u, T1P, T1S;
		    TL = FMA(KP2_000000000, TK, Tv);
		    T1u = T1m + T1t;
		    io[WS(os, 1)] = TL + T1u;
		    io[WS(os, 12)] = T1u - TL;
		    {
			 E T1z, T1A, T1T, T1U;
			 T1z = FMS(KP2_000000000, T1x, T1y);
			 T1A = T1t - T1m;
			 io[WS(os, 5)] = T1z + T1A;
			 io[WS(os, 8)] = T1A - T1z;
			 T1T = T1R - T1Q;
			 T1U = T1O + T1N;
			 io[WS(os, 4)] = T1T - T1U;
			 io[WS(os, 10)] = T1U + T1T;
		    }
		    T1P = T1N - T1O;
		    T1S = T1Q + T1R;
		    io[WS(os, 3)] = T1P + T1S;
		    io[WS(os, 9)] = T1S - T1P;
		    {
			 E T1L, T1M, T1F, T1K;
			 T1L = T1J + T1I;
			 T1M = T1E + T1D;
			 io[WS(os, 6)] = T1L - T1M;
			 io[WS(os, 11)] = T1M + T1L;
			 T1F = T1D - T1E;
			 T1K = T1I - T1J;
			 io[WS(os, 2)] = T1F + T1K;
			 io[WS(os, 7)] = T1K - T1F;
		    }
	       }
	  }
	  {
	       E T2y, T2I, T2J, T2K, T2B, T2L, T2e, T2p, T2u, T2G, T23, T2F, T28, T2t, T2l;
	       E T2o;
	       {
		    E T2w, T2x, T2z, T2A;
		    T2w = FMA(KP387390585, T20, KP265966249 * T1X);
		    T2x = FNMS(KP503537032, T25, KP113854479 * T24);
		    T2y = T2w + T2x;
		    T2I = T2w - T2x;
		    T2J = FMA(KP575140729, T2a, KP174138601 * T2d);
		    T2z = FNMS(KP300238635, T2n, KP011599105 * T2m);
		    T2A = FNMS(KP156891391, T2h, KP256247671 * T2k);
		    T2K = T2z + T2A;
		    T2B = KP1_732050807 * (T2z - T2A);
		    T2L = T2J + T2K;
	       }
	       T2e = FNMS(KP575140729, T2d, KP174138601 * T2a);
	       T2l = FMA(KP256247671, T2h, KP156891391 * T2k);
	       T2o = FMA(KP300238635, T2m, KP011599105 * T2n);
	       T2p = T2l - T2o;
	       T2u = T2e - T2p;
	       T2G = KP1_732050807 * (T2o + T2l);
	       {
		    E T21, T2r, T26, T27, T2s;
		    T21 = FNMS(KP132983124, T20, KP258260390 * T1X);
		    T2r = T22 - T21;
		    T26 = FMA(KP251768516, T24, KP075902986 * T25);
		    T27 = FNMS(KP083333333, To, T1);
		    T2s = T27 - T26;
		    T23 = FMA(KP2_000000000, T21, T22);
		    T2F = T2s - T2r;
		    T28 = FMA(KP2_000000000, T26, T27);
		    T2t = T2r + T2s;
	       }
	       {
		    E T29, T2q, T2N, T2O;
		    T29 = T23 + T28;
		    T2q = FMA(KP2_000000000, T2p, T2e);
		    ro[WS(os, 12)] = T29 - T2q;
		    ro[WS(os, 1)] = T29 + T2q;
		    {
			 E T2v, T2C, T2P, T2Q;
			 T2v = T2t - T2u;
			 T2C = T2y - T2B;
			 ro[WS(os, 10)] = T2v - T2C;
			 ro[WS(os, 4)] = T2v + T2C;
			 T2P = T28 - T23;
			 T2Q = FMS(KP2_000000000, T2K, T2J);
			 ro[WS(os, 5)] = T2P - T2Q;
			 ro[WS(os, 8)] = T2P + T2Q;
		    }
		    T2N = T2F - T2G;
		    T2O = T2L - T2I;
		    ro[WS(os, 11)] = T2N - T2O;
		    ro[WS(os, 6)] = T2N + T2O;
		    {
			 E T2H, T2M, T2D, T2E;
			 T2H = T2F + T2G;
			 T2M = T2I + T2L;
			 ro[WS(os, 7)] = T2H - T2M;
			 ro[WS(os, 2)] = T2H + T2M;
			 T2D = T2t + T2u;
			 T2E = T2y + T2B;
			 ro[WS(os, 3)] = T2D - T2E;
			 ro[WS(os, 9)] = T2D + T2E;
		    }
	       }
	  }
     }
}

static const kdft_desc desc = { 13, "n1_13", {138, 30, 38, 0}, &GENUS, 0, 0, 0, 0 };
void X(codelet_n1_13) (planner *p) {
     X(kdft_register) (p, n1_13, &desc);
}

#endif				/* HAVE_FMA */
