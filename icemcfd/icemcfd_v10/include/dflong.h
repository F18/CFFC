#ifndef DFLONG_H
#define DFLONG_H

#if defined(__mips)
typedef unsigned long long dflong;
#if __mips > 2
#define df_fseek(s, o, w)	fseek64((s), (o), (w))
#else
#define df_fseek(s, o, w)       fseek((s), (o), (w))
#endif

#endif

#if defined(__hpux) && !defined(__LP64__)
/* hp10 */
typedef unsigned long long dflong;
#define df_fseek(s, o, w)       fseek((s), (o), (w))
#endif



#if defined(_WIN32) || defined(__sun)

/* windows */
typedef unsigned long dflong;
#define df_fseek(s, o, w)	fseek((s), (o), (w))

#endif

#ifdef __linux
typedef unsigned long long dflong;
#define df_fseek(s, o, w)       fseek((s), (o), (w))
#endif

#ifndef df_fseek

#if defined(__hpux) && defined(__LP64__) && !defined(_OFF_T)
typedef unsigned long off_t;
#endif
#if defined(__linux) || defined(__MACH__)
#include <sys/types.h>
#endif

/*
 *  this is most os's 
 */

typedef off_t dflong;
#define df_fseek(s, o, w)	fseeko((s), (o), (w))

#endif


#endif

