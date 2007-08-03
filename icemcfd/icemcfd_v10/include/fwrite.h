#ifndef FWRITE_H
#define FWRITE_H

#if defined(_WIN32) || defined(linux) || defined(__alpha) || defined(__i386)

#include <stdio.h>
#include "util/fortran.h"

/* these are defined on machines who's sex is opposite that of sgi/sun/ibm */

EXTERN_C int fwrite_ushort(const unsigned short *ptr, 
		size_t size, size_t nitems, FILE *stream);
EXTERN_C int fwrite_short(const short *ptr, size_t size, size_t nitems, 
		FILE *stream);
EXTERN_C int fread_short(short *ptr, size_t size, size_t nitems, 
		FILE *stream);
EXTERN_C int fread_ushort(unsigned short *ptr, size_t size, size_t nitems, 
		FILE *stream);
EXTERN_C int fwrite_int(const int *ptr, size_t size, size_t nitems,
		FILE *stream);
EXTERN_C int fwrite_uint(const unsigned int *ptr, 
			size_t size, size_t nitems, FILE *stream);
EXTERN_C int fwrite_float(const float *ptr, 
			size_t size, size_t nitems, FILE *stream);
EXTERN_C int fread_int(int *ptr, size_t size, size_t nitems, FILE *stream);
EXTERN_C int fread_uint(unsigned int *ptr, size_t size, 
					size_t nitems, FILE *stream);
EXTERN_C int fread_float(float *ptr, size_t size, size_t nitems, FILE *stream);
EXTERN_C int fwrite_double(const double *ptr, size_t size, size_t nitems, 
				FILE *stream);
EXTERN_C int fread_double(double *ptr, size_t size, size_t nitems, 
				FILE *stream);

#if defined(_WIN32) && !defined(_WINSOCKAPI_) && !defined(_WINSOCK2API_)

EXTERN_C int htonl(int x);
EXTERN_C int ntohl(int x);
EXTERN_C short htons(short x);
EXTERN_C short ntohs(short x);

#endif /* _WIN32 */

#else /* _WIN32 linux etc */

/* otherwise they are simply fwrite's */

#define fwrite_short(a, b, c, d)	fwrite(a, b, c, d)
#define fwrite_ushort(a, b, c, d)	fwrite(a, b, c, d)
#define fread_ushort(a, b, c, d) 	fread(a, b, c, d)
#define fread_short(a, b, c, d) 	fread(a, b, c, d)
#define fwrite_float(a, b, c, d)	fwrite(a, b, c, d)
#define fwrite_int(a, b, c, d)		fwrite(a, b, c, d)
#define fwrite_uint(a, b, c, d)		fwrite(a, b, c, d)
#define fread_int(a, b, c, d) 		fread(a, b, c, d)
#define fread_uint(a, b, c, d) 		fread(a, b, c, d)
#define fread_float(a, b, c, d)		fread(a, b, c, d)
#define fwrite_double(a, b, c, d) 	fwrite(a, b, c, d)
#define fread_double(a, b, c, d) 	fread(a, b, c, d)

#endif /* _WIN32 linux etc */

#endif /* FWRITE_H */
