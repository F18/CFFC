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
#include "rdft.h"

X(plan) XGURU(dft_r2c)(int rank, const IODIM *dims,
		       int howmany_rank,
		       const IODIM *howmany_dims,
		       R *in, C *out, unsigned flags)
{
     R *ro, *io;

     if (!GURU_KOSHERP(rank, dims, howmany_rank, howmany_dims)) return 0;

     EXTRACT_REIM(FFT_SIGN, out, &ro, &io);

     return X(mkapiplan)(
	  0, flags,
	  X(mkproblem_rdft2_d_3pointers)(
	       MKTENSOR_IODIMS(rank, dims, 1, 2),
	       MKTENSOR_IODIMS(howmany_rank, howmany_dims, 1, 2),
	       TAINT_UNALIGNED(in, flags),
	       TAINT_UNALIGNED(ro, flags),
	       TAINT_UNALIGNED(io, flags), R2HC));
}
