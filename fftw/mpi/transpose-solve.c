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

#include "mpi-transpose.h"

/* use the apply() operation for MPI_TRANSPOSE problems */
void XM(transpose_solve)(const plan *ego_, const problem *p_)
{
     const plan_mpi_transpose *ego = (const plan_mpi_transpose *) ego_;
     const problem_mpi_transpose *p = (const problem_mpi_transpose *) p_;
     ego->apply(ego_, UNTAINT(p->I), UNTAINT(p->O));
}
