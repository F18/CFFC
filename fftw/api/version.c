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


#include "api.h"

const char X(cc)[] = FFTW_CC;
const char X(codelet_optim)[] = CODELET_OPTIM;

const char X(version)[] = PACKAGE "-" PACKAGE_VERSION

#if HAVE_FMA
   "-fma"
#endif

#if HAVE_SSE
   "-sse"
#endif

#if HAVE_SSE2
   "-sse2"
#endif

#if HAVE_ALTIVEC
   "-altivec"
#endif

;
