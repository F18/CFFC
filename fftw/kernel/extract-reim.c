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

#include "ifftw.h"

/* decompose complex pointer into real and imaginary parts.
   Flip real and imaginary if there the sign does not match
   FFTW's idea of what the sign should be */

void X(extract_reim)(int sign, R *c, R **r, R **i)
{
     if (sign == FFT_SIGN) {
          *r = c + 0;
          *i = c + 1;
     } else {
          *r = c + 1;
          *i = c + 0;
     }
}
