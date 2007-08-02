/************ ICEMCFDctype.c *******************************
  Added to allow ICEMCFD libraries to be used with   
  gcc 3.x compilers and glibc 2.3.  
************************************************************/

#ifdef _GNU_GCC_296
// do nothing
#else 

// assumming using ICC or GCC V3.x
#ifdef _ICEMCFD_VERSION

/*
 Static ctype data for IFC-7.1 / RedHat-9 workaround.
 This is taken from the GLIBC source code.
 Hack by Joe Krahn <krahn@niehs.nih.gov>

 To use, compile with gcc (ifc should work as well):

   gcc -c ctype.c

 Next, include the resulting ctype.o when creating executables.
 Do this automatically using a ~/.ifcrc file with a line like this:
  
   -Wl,/some/path/ctype.o

 That's all.

 My .ifcrc also turns off the typically excessive warnings:
   -cm -w90 -w95 -Wl,/home/krahn/Prog/ctype/ctype.o
 
 To IFC maintainters: it would be nice to flag on/off specific
 warning types, like -woff=22,33,44

*/

/*
  file: ctype.c

  Derived from C-ctype.c and ctype-info.c in the GLIBC 2.3 source.

  To compile: (either gcc or icc is OK)
      icc -c ctype.c

  To use, just add the object file ctype.o at link stage. Example:
      ifc -o hello hello.f ctype.o

*/

/* Copyright (C) 1995-1999, 2000, 2001, 2002 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Ulrich Drepper <drepper@cygnus.com>, 1995.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

#ifdef HAVE_CONFIG_H
#  include<config.h>
#endif

#include<errno.h>

// THIS BREAKS errno FYI !!! ( only important if your are depending on it)
#ifdef errno
# undef errno
#endif
int errno=0;

#include <stdint.h>

/* This table's entries are taken from POSIX.2 Table 2-6
   ``LC_CTYPE Category Definition in the POSIX Locale''.

   The `_nl_C_LC_CTYPE_width' array is a GNU extension.

   In the `_nl_C_LC_CTYPE_class' array the value for EOF (== -1)
   is set to always return 0 and the conversion arrays return EOF.  */

const char _nl_C_LC_CTYPE_class[768] =
  /* 0x80 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x86 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x8c */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x92 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x98 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x9e */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xa4 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xaa */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xb0 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xb6 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xbc */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xc2 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xc8 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xce */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xd4 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xda */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xe0 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xe6 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xec */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xf2 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xf8 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xfe */ "\000\000" "\000\000" "\002\000" "\002\000" "\002\000" "\002\000"
  /* 0x04 */ "\002\000" "\002\000" "\002\000" "\002\000" "\002\000" "\003\040"
  /* 0x0a */ "\002\040" "\002\040" "\002\040" "\002\040" "\002\000" "\002\000"
  /* 0x10 */ "\002\000" "\002\000" "\002\000" "\002\000" "\002\000" "\002\000"
  /* 0x16 */ "\002\000" "\002\000" "\002\000" "\002\000" "\002\000" "\002\000"
  /* 0x1c */ "\002\000" "\002\000" "\002\000" "\002\000" "\001\140" "\004\300"
  /* 0x22 */ "\004\300" "\004\300" "\004\300" "\004\300" "\004\300" "\004\300"
  /* 0x28 */ "\004\300" "\004\300" "\004\300" "\004\300" "\004\300" "\004\300"
  /* 0x2e */ "\004\300" "\004\300" "\010\330" "\010\330" "\010\330" "\010\330"
  /* 0x34 */ "\010\330" "\010\330" "\010\330" "\010\330" "\010\330" "\010\330"
  /* 0x3a */ "\004\300" "\004\300" "\004\300" "\004\300" "\004\300" "\004\300"
  /* 0x40 */ "\004\300" "\010\325" "\010\325" "\010\325" "\010\325" "\010\325"
  /* 0x46 */ "\010\325" "\010\305" "\010\305" "\010\305" "\010\305" "\010\305"
  /* 0x4c */ "\010\305" "\010\305" "\010\305" "\010\305" "\010\305" "\010\305"
  /* 0x52 */ "\010\305" "\010\305" "\010\305" "\010\305" "\010\305" "\010\305"
  /* 0x58 */ "\010\305" "\010\305" "\010\305" "\004\300" "\004\300" "\004\300"
  /* 0x5e */ "\004\300" "\004\300" "\004\300" "\010\326" "\010\326" "\010\326"
  /* 0x64 */ "\010\326" "\010\326" "\010\326" "\010\306" "\010\306" "\010\306"
  /* 0x6a */ "\010\306" "\010\306" "\010\306" "\010\306" "\010\306" "\010\306"
  /* 0x70 */ "\010\306" "\010\306" "\010\306" "\010\306" "\010\306" "\010\306"
  /* 0x76 */ "\010\306" "\010\306" "\010\306" "\010\306" "\010\306" "\004\300"
  /* 0x7c */ "\004\300" "\004\300" "\004\300" "\002\000" "\000\000" "\000\000"
  /* 0x82 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x88 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x8e */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x94 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x9a */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xa0 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xa6 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xac */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xb2 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xb8 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xbe */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xc4 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xca */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xd0 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xd6 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xdc */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xe2 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xe8 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xee */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xf4 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xfa */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
;
const char _nl_C_LC_CTYPE_class32[1024] =
  /* 0x00 */ "\000\000\002\000" "\000\000\002\000" "\000\000\002\000"
  /* 0x03 */ "\000\000\002\000" "\000\000\002\000" "\000\000\002\000"
  /* 0x06 */ "\000\000\002\000" "\000\000\002\000" "\000\000\002\000"
  /* 0x09 */ "\000\000\003\040" "\000\000\002\040" "\000\000\002\040"
  /* 0x0c */ "\000\000\002\040" "\000\000\002\040" "\000\000\002\000"
  /* 0x0f */ "\000\000\002\000" "\000\000\002\000" "\000\000\002\000"
  /* 0x12 */ "\000\000\002\000" "\000\000\002\000" "\000\000\002\000"
  /* 0x15 */ "\000\000\002\000" "\000\000\002\000" "\000\000\002\000"
  /* 0x18 */ "\000\000\002\000" "\000\000\002\000" "\000\000\002\000"
  /* 0x1b */ "\000\000\002\000" "\000\000\002\000" "\000\000\002\000"
  /* 0x1e */ "\000\000\002\000" "\000\000\002\000" "\000\000\001\140"
  /* 0x21 */ "\000\000\004\300" "\000\000\004\300" "\000\000\004\300"
  /* 0x24 */ "\000\000\004\300" "\000\000\004\300" "\000\000\004\300"
  /* 0x27 */ "\000\000\004\300" "\000\000\004\300" "\000\000\004\300"
  /* 0x2a */ "\000\000\004\300" "\000\000\004\300" "\000\000\004\300"
  /* 0x2d */ "\000\000\004\300" "\000\000\004\300" "\000\000\004\300"
  /* 0x30 */ "\000\000\010\330" "\000\000\010\330" "\000\000\010\330"
  /* 0x33 */ "\000\000\010\330" "\000\000\010\330" "\000\000\010\330"
  /* 0x36 */ "\000\000\010\330" "\000\000\010\330" "\000\000\010\330"
  /* 0x39 */ "\000\000\010\330" "\000\000\004\300" "\000\000\004\300"
  /* 0x3c */ "\000\000\004\300" "\000\000\004\300" "\000\000\004\300"
  /* 0x3f */ "\000\000\004\300" "\000\000\004\300" "\000\000\010\325"
  /* 0x42 */ "\000\000\010\325" "\000\000\010\325" "\000\000\010\325"
  /* 0x45 */ "\000\000\010\325" "\000\000\010\325" "\000\000\010\305"
  /* 0x48 */ "\000\000\010\305" "\000\000\010\305" "\000\000\010\305"
  /* 0x4b */ "\000\000\010\305" "\000\000\010\305" "\000\000\010\305"
  /* 0x4e */ "\000\000\010\305" "\000\000\010\305" "\000\000\010\305"
  /* 0x51 */ "\000\000\010\305" "\000\000\010\305" "\000\000\010\305"
  /* 0x54 */ "\000\000\010\305" "\000\000\010\305" "\000\000\010\305"
  /* 0x57 */ "\000\000\010\305" "\000\000\010\305" "\000\000\010\305"
  /* 0x5a */ "\000\000\010\305" "\000\000\004\300" "\000\000\004\300"
  /* 0x5d */ "\000\000\004\300" "\000\000\004\300" "\000\000\004\300"
  /* 0x60 */ "\000\000\004\300" "\000\000\010\326" "\000\000\010\326"
  /* 0x63 */ "\000\000\010\326" "\000\000\010\326" "\000\000\010\326"
  /* 0x66 */ "\000\000\010\326" "\000\000\010\306" "\000\000\010\306"
  /* 0x69 */ "\000\000\010\306" "\000\000\010\306" "\000\000\010\306"
  /* 0x6c */ "\000\000\010\306" "\000\000\010\306" "\000\000\010\306"
  /* 0x6f */ "\000\000\010\306" "\000\000\010\306" "\000\000\010\306"
  /* 0x72 */ "\000\000\010\306" "\000\000\010\306" "\000\000\010\306"
  /* 0x75 */ "\000\000\010\306" "\000\000\010\306" "\000\000\010\306"
  /* 0x78 */ "\000\000\010\306" "\000\000\010\306" "\000\000\010\306"
  /* 0x7b */ "\000\000\004\300" "\000\000\004\300" "\000\000\004\300"
  /* 0x7e */ "\000\000\004\300" "\000\000\002\000" "\000\000\000\000"
  /* 0x81 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0x84 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0x87 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0x8a */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0x8d */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0x90 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0x93 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0x96 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0x99 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0x9c */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0x9f */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xa2 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xa5 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xa8 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xab */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xae */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xb1 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xb4 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xb7 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xba */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xbd */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xc0 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xc3 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xc6 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xc9 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xcc */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xcf */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xd2 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xd5 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xd8 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xdb */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xde */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xe1 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xe4 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xe7 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xea */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xed */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xf0 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xf3 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xf6 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xf9 */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xfc */ "\000\000\000\000" "\000\000\000\000" "\000\000\000\000"
  /* 0xff */ "\000\000\000\000"
;
const uint32_t _nl_C_LC_CTYPE_toupper[384] =
{
  /* 0x80 */ 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
  /* 0x88 */ 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
  /* 0x90 */ 0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97,
  /* 0x98 */ 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
  /* 0xa0 */ 0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
  /* 0xa8 */ 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
  /* 0xb0 */ 0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7,
  /* 0xb8 */ 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf,
  /* 0xc0 */ 0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7,
  /* 0xc8 */ 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf,
  /* 0xd0 */ 0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7,
  /* 0xd8 */ 0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf,
  /* 0xe0 */ 0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7,
  /* 0xe8 */ 0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
  /* 0xf0 */ 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7,
  /* 0xf8 */ 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xffffffff,
  /* 0x00 */ 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
  /* 0x08 */ 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
  /* 0x10 */ 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
  /* 0x18 */ 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
  /* 0x20 */ 0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
  /* 0x28 */ 0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
  /* 0x30 */ 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
  /* 0x38 */ 0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f,
  /* 0x40 */ 0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47,
  /* 0x48 */ 0x48, 0x49, 0x4a, 0x4b, 0x4c, 0x4d, 0x4e, 0x4f,
  /* 0x50 */ 0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57,
  /* 0x58 */ 0x58, 0x59, 0x5a, 0x5b, 0x5c, 0x5d, 0x5e, 0x5f,
  /* 0x60 */ 0x60, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47,
  /* 0x68 */ 0x48, 0x49, 0x4a, 0x4b, 0x4c, 0x4d, 0x4e, 0x4f,
  /* 0x70 */ 0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57,
  /* 0x78 */ 0x58, 0x59, 0x5a, 0x7b, 0x7c, 0x7d, 0x7e, 0x7f,
  /* 0x80 */ 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
  /* 0x88 */ 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
  /* 0x90 */ 0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97,
  /* 0x98 */ 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
  /* 0xa0 */ 0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
  /* 0xa8 */ 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
  /* 0xb0 */ 0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7,
  /* 0xb8 */ 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf,
  /* 0xc0 */ 0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7,
  /* 0xc8 */ 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf,
  /* 0xd0 */ 0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7,
  /* 0xd8 */ 0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf,
  /* 0xe0 */ 0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7,
  /* 0xe8 */ 0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
  /* 0xf0 */ 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7,
  /* 0xf8 */ 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff
};
const uint32_t _nl_C_LC_CTYPE_tolower[384] =
{
  /* 0x80 */ 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
  /* 0x88 */ 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
  /* 0x90 */ 0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97,
  /* 0x98 */ 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
  /* 0xa0 */ 0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
  /* 0xa8 */ 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
  /* 0xb0 */ 0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7,
  /* 0xb8 */ 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf,
  /* 0xc0 */ 0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7,
  /* 0xc8 */ 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf,
  /* 0xd0 */ 0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7,
  /* 0xd8 */ 0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf,
  /* 0xe0 */ 0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7,
  /* 0xe8 */ 0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
  /* 0xf0 */ 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7,
  /* 0xf8 */ 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xffffffff,
  /* 0x00 */ 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
  /* 0x08 */ 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
  /* 0x10 */ 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
  /* 0x18 */ 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
  /* 0x20 */ 0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
  /* 0x28 */ 0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
  /* 0x30 */ 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
  /* 0x38 */ 0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f,
  /* 0x40 */ 0x40, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67,
  /* 0x48 */ 0x68, 0x69, 0x6a, 0x6b, 0x6c, 0x6d, 0x6e, 0x6f,
  /* 0x50 */ 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77,
  /* 0x58 */ 0x78, 0x79, 0x7a, 0x5b, 0x5c, 0x5d, 0x5e, 0x5f,
  /* 0x60 */ 0x60, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67,
  /* 0x68 */ 0x68, 0x69, 0x6a, 0x6b, 0x6c, 0x6d, 0x6e, 0x6f,
  /* 0x70 */ 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77,
  /* 0x78 */ 0x78, 0x79, 0x7a, 0x7b, 0x7c, 0x7d, 0x7e, 0x7f,
  /* 0x80 */ 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
  /* 0x88 */ 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
  /* 0x90 */ 0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97,
  /* 0x98 */ 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
  /* 0xa0 */ 0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
  /* 0xa8 */ 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
  /* 0xb0 */ 0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7,
  /* 0xb8 */ 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf,
  /* 0xc0 */ 0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7,
  /* 0xc8 */ 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf,
  /* 0xd0 */ 0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7,
  /* 0xd8 */ 0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf,
  /* 0xe0 */ 0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7,
  /* 0xe8 */ 0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
  /* 0xf0 */ 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7,
  /* 0xf8 */ 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff
};

#define STRUCT_CTYPE_CLASS(p, q) \
  struct                                                                      \
    {                                                                         \
      uint32_t isctype_data[8];                                               \
      uint32_t header[5];                                                     \
      uint32_t level1[1];                                                     \
      uint32_t level2[1 << q];                                                \
      uint32_t level3[1 << p];                                                \
    }

const STRUCT_CTYPE_CLASS(1, 1) _nl_C_LC_CTYPE_class_upper =
{
  { 0x00000000, 0x00000000, 0x07fffffe, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 7, 1, 6, 1, 1 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 0, 8 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0x07fffffe, 0x00000000 }
};
const STRUCT_CTYPE_CLASS(1, 1) _nl_C_LC_CTYPE_class_lower =
{
  { 0x00000000, 0x00000000, 0x00000000, 0x07fffffe,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 7, 1, 6, 1, 1 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 0, 8 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0x00000000, 0x07fffffe }
};
const STRUCT_CTYPE_CLASS(1, 1) _nl_C_LC_CTYPE_class_alpha =
{
  { 0x00000000, 0x00000000, 0x07fffffe, 0x07fffffe,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 7, 1, 6, 1, 1 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 0, 8 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0x07fffffe, 0x07fffffe }
};
const STRUCT_CTYPE_CLASS(1, 0) _nl_C_LC_CTYPE_class_digit =
{
  { 0x00000000, 0x03ff0000, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 6, 1, 6, 0, 1 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 7 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0x00000000, 0x03ff0000 }
};
const STRUCT_CTYPE_CLASS(2, 0) _nl_C_LC_CTYPE_class_xdigit =
{
  { 0x00000000, 0x03ff0000, 0x0000007e, 0x0000007e,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 7, 1, 7, 0, 3 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 7 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0x00000000, 0x03ff0000, 0x0000007e, 0x0000007e }
};
const STRUCT_CTYPE_CLASS(1, 0) _nl_C_LC_CTYPE_class_space =
{
  { 0x00003e00, 0x00000001, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 6, 1, 6, 0, 1 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 7 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0x00003e00, 0x00000001 }
};
const STRUCT_CTYPE_CLASS(2, 0) _nl_C_LC_CTYPE_class_print =
{
  { 0x00000000, 0xffffffff, 0xffffffff, 0x7fffffff,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 7, 1, 7, 0, 3 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 7 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0x00000000, 0xffffffff, 0xffffffff, 0x7fffffff }
};
const STRUCT_CTYPE_CLASS(2, 0) _nl_C_LC_CTYPE_class_graph =
{
  { 0x00000000, 0xfffffffe, 0xffffffff, 0x7fffffff,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 7, 1, 7, 0, 3 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 7 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0x00000000, 0xfffffffe, 0xffffffff, 0x7fffffff }
};
const STRUCT_CTYPE_CLASS(1, 0) _nl_C_LC_CTYPE_class_blank =
{
  { 0x00000200, 0x00000001, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 6, 1, 6, 0, 1 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 7 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0x00000200, 0x00000001 }
};
const STRUCT_CTYPE_CLASS(2, 0) _nl_C_LC_CTYPE_class_cntrl =
{
  { 0xffffffff, 0x00000000, 0x00000000, 0x80000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 7, 1, 7, 0, 3 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 7 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0xffffffff, 0x00000000, 0x00000000, 0x80000000 }
};
const STRUCT_CTYPE_CLASS(2, 0) _nl_C_LC_CTYPE_class_punct =
{
  { 0x00000000, 0xfc00fffe, 0xf8000001, 0x78000001,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 7, 1, 7, 0, 3 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 7 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0x00000000, 0xfc00fffe, 0xf8000001, 0x78000001 }
};
const STRUCT_CTYPE_CLASS(2, 0) _nl_C_LC_CTYPE_class_alnum =
{
  { 0x00000000, 0x03ff0000, 0x07fffffe, 0x07fffffe,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  },
  { 7, 1, 7, 0, 3 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 7 * sizeof (uint32_t) },
  /* 3rd-level table */
  { 0x00000000, 0x03ff0000, 0x07fffffe, 0x07fffffe }
};

const struct
{
  uint32_t header[5];
  uint32_t level1[1];
  uint32_t level2[4];
  int32_t level3[32];
}
_nl_C_LC_CTYPE_map_toupper =
{
  { 7, 1, 5, 3, 31 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 0, 0, 0, 10 * sizeof (uint32_t) },
  /* 3rd-level table */
  {
    0x00000000, 0xffffffe0, 0xffffffe0, 0xffffffe0,
    0xffffffe0, 0xffffffe0, 0xffffffe0, 0xffffffe0,
    0xffffffe0, 0xffffffe0, 0xffffffe0, 0xffffffe0,
    0xffffffe0, 0xffffffe0, 0xffffffe0, 0xffffffe0,
    0xffffffe0, 0xffffffe0, 0xffffffe0, 0xffffffe0,
    0xffffffe0, 0xffffffe0, 0xffffffe0, 0xffffffe0,
    0xffffffe0, 0xffffffe0, 0xffffffe0, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  }
},
_nl_C_LC_CTYPE_map_tolower =
{
  { 7, 1, 5, 3, 31 },
  /* 1st-level table */
  { 6 * sizeof (uint32_t) },
  /* 2nd-level table */
  { 0, 0, 10 * sizeof (uint32_t), 0 },
  /* 3rd-level table */
  {
    0x00000000, 0x00000020, 0x00000020, 0x00000020,
    0x00000020, 0x00000020, 0x00000020, 0x00000020,
    0x00000020, 0x00000020, 0x00000020, 0x00000020,
    0x00000020, 0x00000020, 0x00000020, 0x00000020,
    0x00000020, 0x00000020, 0x00000020, 0x00000020,
    0x00000020, 0x00000020, 0x00000020, 0x00000020,
    0x00000020, 0x00000020, 0x00000020, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000
  }
};

#define b(t,x,o) (((const t *) _nl_C_LC_CTYPE_##x) + o)

const uint16_t *__ctype_b = b (uint16_t, class, 128);
const uint32_t *__ctype32_b = b (uint32_t, class32, 0);
const int32_t *__ctype_tolower = b (int32_t, tolower, 128);
const int32_t *__ctype_toupper = b (int32_t, toupper, 128);
const uint32_t *__ctype32_tolower = b (uint32_t, tolower, 128);
const uint32_t *__ctype32_toupper = b (uint32_t, toupper, 128);

#endif // end of _ICEMCFD_VERSION 
#endif // end of _GCC_VERSION
