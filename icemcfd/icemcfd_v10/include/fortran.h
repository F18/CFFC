
#ifndef	FORTRAN_H
#define	FORTRAN_H

/**************************************************************/
/* Include file used to interface back and forth with FORTRAN */
/**************************************************************/

/* L4NG is the type of a FORTRAN INTEGER, it must be at least 4 bytes and   */
/*	some systems is 8 bytes.  It is normally declared to be a long,     */
/*	but one system is a long long.  A variable should not be declared   */
/*	to be long but instead L4NG.                                        */
/* UL4NG is the same as above, but is the unsigned version of L4NG          */

/* FTNPTR is a variable used to hold either a pointer or an index of an     */
/*        offset of a pointer                                               */

/* FMPREFIX is a prefix needed for prototyping or declaring subroutines.    */
/*	    It is only currently needed for the Windows NT system.          */
/* FCNAME is used to access FORTRAN commons from C.  To access COMMON /ABC/ */
/*	  #define Abc	FCNAME(abc,ABC)                                     */
/*	  struct {} Abc;		                                    */
/* FMNAME is used when you prototype a FORTRAN routine for calling from C   */
/*        (which you HAVE to do), or when you create the subroutine         */
/*	  INTEGER FMNAME(abc,ABC) (INTEGER ival, REAL rval);                */
/*	  void FMNAME(abc,ABC) (INTEGER ival, REAL rval)                    */
/*	  { *ival = 1; *rval = 2.0; return; }                               */
/* FMCALL is used when you call a FORTRAN subroutine from C                 */
/*	  VINTEGER ival;                                                    */
/*	  VREAL rval;                                                       */
/*	  FMCALL(abc,ABC) (&ival, &rval);                                   */

/* FNC_CPTR is used to call a function                                      */
/* FNC_PPTR is used in receiving a function from FORTRAN                    */
/* FNC_FPTR is used as an argument to FORTRAN                               */
/*	    void FMNAME(callme,CALLME) (void);                              */
/*	    #define CallMe FMCALL(callme,CALLME)                            */
/*	    INTEGER FMNAME(abc,ABC) (FNC_PPTR(func), INTEGER ival)          */
/*	    { (FNC_CPTR(func)) (ival, FNC_FPTR(CallMe)); }                  */

/* STR_PSTR is used in receiving arguments from FORTRAN                     */
/* STR_PLEN is used in end of arguments from FORTRAN (no comma before it)   */
/* STR_PTR  is used to get the address of a string from FORTRAN             */
/* STR_LEN  is used to get the length of a string from FORTRAN              */
/*	    INTEGER FMNAME(abc,ABC) (STR_PSTR(str), INTEGER ival            */
/*				     STR_PLEN(str))                         */
/*	    { char *pointer = STR_PTR(str); int length = STR_LEN(str); }    */

/* STR_SDEF is used to define a local string to be passed to FORTRAN        */
/* STR_SPTR is used to set the address of a local string to FORTRAN         */
/* STR_SLEN is used to set the length of a local string to FORTRAN          */
/* STR_FSTR is used as an argument to FORTRAN                               */
/* STR_FLEN is used at end of arguments when passing a string to FORTRAN    */
/*	    (no comma before it)					    */
/*									    */
/*	    char local_str[128];					    */
/*	    STR_SDEF(abc);						    */
/*	    VINTEGER ival;						    */
/*									    */
/*          STR_SPTR(abc, local_str);					    */
/*	    STR_SLEN(abc, sizeof (local_str));				    */
/*	    FMCALL(doit,DOIT) (STR_FSTR(abc), &ival STR_FLEN(abc));	    */

#ifdef	_AIX
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#	define FCNAME(lname, uname) lname 	/* Fortran commons  */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) lname  	/* Fortran Module Name */
#	define FMCALL(lname, uname) lname 	/* Fortran Module Call */
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#	define STR_PSTR(str)	char *str
#	define STR_PLEN(str)	, int Len ## str
#	define STR_PTR(str)	str
#	define STR_LEN(str)	Len ## str
#	define STR_SDEF(str)	char *str; int Len ## str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_SLEN(str, len) Len ## str = (int)len
#	define STR_FSTR(str)	(char *)str
#	define STR_FLEN(str)	, (int)Len ## str
#endif /* AIX */

#ifdef	__convex__
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#	define FCNAME(lname, uname) _ ## lname ## _   /* Fortran commons  */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) lname ## _	      /* Fortran Module Name */
#	define FMCALL(lname, uname) lname ## _	      /* Fortran Module Call */
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#	define STR_PSTR(str)	char *str
#	define STR_PLEN(str)	, int Len ## str
#	define STR_PTR(str)	str
#	define STR_LEN(str)	Len ## str
#	define STR_SDEF(str)	char *str; int Len ## str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_SLEN(str, len) Len ## str = (int)len
#	define STR_FSTR(str)	(char *)str
#	define STR_FLEN(str)	, (int)Len ## str
#endif

#ifdef	_CRAY
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
/* #	include "/usr/include/fortran.h" */
#       include "/opt/ctl/craylibs/craylibs/include/fortran.h"
#	define FCNAME(lname, uname) uname /* Fortran commons */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) uname /* Fortran Module Name */
#	define FMCALL(lname, uname) uname /* Fortran Module Call */
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#	define STR_PSTR(str)	_fcd str
#	define STR_PLEN(str)
#	define STR_PTR(str)	_fcdtocp (str)
#	define STR_LEN(str)	_fcdlen (str)
#	define STR_SDEF(str)	_fcd str
#	define STR_SPTR(str, ptr) str = _cptofcd (ptr, _fcdlen (str))
#	define STR_SLEN(str, len) str = _cptofcd (_fcdtocp (str), len)
#	define STR_FSTR(str)	(_fcd)str
#	define STR_FLEN(str)
#endif

#ifdef  __hp9000s300
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#       define FCNAME(lname, uname) lname       /* Fortran commons  */
#       define FMPREFIX                         /* Fortran Module Prefix */
#       define FMNAME(lname, uname) lname       /* Fortran Module Name */
#       define FMCALL(lname, uname) lname       /* Fortran Module Call */
#       define FNC_CPTR(func)   * func
#       define FNC_PPTR(func)   * func
#       define FNC_FPTR(func)   func
#       define STR_PSTR(str)    char *str
#       define STR_PLEN(str)    , int Len ## str
#       define STR_PTR(str)     str
#       define STR_LEN(str)     Len ## str
#       define STR_SDEF(str)    char *str; int Len ## str
#       define STR_SPTR(str, ptr) str = (char *)ptr
#       define STR_SLEN(str, len) Len ## str = (int)len
#       define STR_FSTR(str)    (char *)str
#       define STR_FLEN(str)    , (int)Len ## str
#endif

#ifdef  __hp9000s700
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#       define FCNAME(lname, uname) lname       /* Fortran commons  */
#       define FMPREFIX                         /* Fortran Module Prefix */
#       define FMNAME(lname, uname) lname       /* Fortran Module Name */
#       define FMCALL(lname, uname) lname       /* Fortran Module Call */
#       define FNC_CPTR(func)   * func
#       define FNC_PPTR(func)   * func
#       define FNC_FPTR(func)   func
#       define STR_PSTR(str)    char *str
#       define STR_PLEN(str)    , int Len ## str
#       define STR_PTR(str)     str
#       define STR_LEN(str)     Len ## str
#       define STR_SDEF(str)    char *str; int Len ## str
#       define STR_SPTR(str, ptr) str = (char *)ptr
#       define STR_SLEN(str, len) Len ## str = (int)len
#       define STR_FSTR(str)    (char *)str
#       define STR_FLEN(str)    , (int)Len ## str
#else
#ifdef  __hp9000s800
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
/*  HP-UX 8.0
#       define FCNAME(lname, uname) lname
#       define FMPREFIX
#       define FMNAME(lname, uname) lname
#       define FMCALL(lname, uname) lname
#       define FNC_CPTR(func)   ** func
#       define FNC_PPTR(func)   ** func
#       define FNC_FPTR(func)   &(func)
#       define STR_PSTR(str)    char *str, int Len ## str
#       define STR_PLEN(str)
#       define STR_PTR(str)     str
#       define STR_LEN(str)     Len ## str
#       define STR_SDEF(str)    char *str; int Len ## str
#       define STR_SPTR(str, ptr) str = (char *)ptr
#       define STR_SLEN(str, len) Len ## str = (int)len
#       define STR_FSTR(str)    (char *)str, (int)Len ## str
#       define STR_FLEN(str)
*/
#       define FCNAME(lname, uname) lname       /* Fortran commons  */
#       define FMPREFIX                         /* Fortran Module Prefix */
#       define FMNAME(lname, uname) lname       /* Fortran Module Name */
#       define FMCALL(lname, uname) lname       /* Fortran Module Call */
#       define FNC_CPTR(func)   * func
#       define FNC_PPTR(func)   * func
#       define FNC_FPTR(func)   func
#       define STR_PSTR(str)    char *str
#       define STR_PLEN(str)    , int Len ## str
#       define STR_PTR(str)     str
#       define STR_LEN(str)     Len ## str
#       define STR_SDEF(str)    char *str; int Len ## str
#       define STR_SPTR(str, ptr) str = (char *)ptr
#       define STR_SLEN(str, len) Len ## str = (int)len
#       define STR_FSTR(str)    (char *)str
#       define STR_FLEN(str)    , (int)Len ## str
#endif
#endif

#ifdef  __linux
#     define L4NG long
#     define UL4NG unsigned long
#     define FTNPTR L4NG
#     define FCNAME(lname, uname) lname ## _  /* Fortran commons  */
#     define FMPREFIX                         /* Fortran Module Prefix */
#     define FMNAME(lname, uname) lname ## _  /* Fortran Module Name */
#     define FMCALL(lname, uname) lname ## _  /* Fortran Module Call */
#     define FNC_CPTR(func)   * func
#     define FNC_PPTR(func)   * func
#     define FNC_FPTR(func)   func
#     define STR_PSTR(str)    char *str
#     define STR_PLEN(str)    , int Len ## str
#     define STR_PTR(str)     str
#     define STR_LEN(str)     Len ## str
#     define STR_SDEF(str)    char *str; int Len ## str
#     define STR_SPTR(str, ptr) str = (char *)ptr
#     define STR_SLEN(str, len) Len ## str = (int)len
#     define STR_FSTR(str)    (char *)str
#     define STR_FLEN(str)    , (int)Len ## str
#endif

#ifdef	__sgi
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#	define FCNAME(lname, uname) lname ## _	/* Fortran commons  */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) lname ## _	/* Fortran Module Name */
#	define FMCALL(lname, uname) lname ## _	/* Fortran Module Call */
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#	define STR_PSTR(str)	char *str
#	define STR_PLEN(str)	, int Len ## str
#	define STR_PTR(str)	str
#	define STR_LEN(str)	Len ## str
#	define STR_SDEF(str)	char *str; int Len ## str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_SLEN(str, len) Len ## str = (int)len
#	define STR_FSTR(str)	(char *)str
#	define STR_FLEN(str)	, (int)Len ## str
#endif

#ifdef	sparc
#       define L4NG long
#       define UL4NG unsigned long
#       define FMPREFIX                         /* Fortran Module Prefix */
#       define FCNAME(lname, uname) lname ## _  /* Fortran commons  */
#       define FMNAME(lname, uname) lname ## _  /* Fortran Module Name */
#       define FMCALL(lname, uname) lname ## _  /* Fortran Module Call */
#       define STR_PLEN(str)    , int Len ## str
#       define STR_LEN(str)     Len ## str
#       define STR_SDEF(str)    char *str; int Len ## str
#       define STR_SLEN(str, len) Len ## str = (int)len
#       define STR_FLEN(str)    , (int)Len ## str
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#	define STR_PSTR(str)	char *str
#	define STR_PTR(str)	str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_FSTR(str)	(char *)str
#endif

#ifdef	__ultrix
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#	define FCNAME(lname, uname) lname ## _	/* Fortran commons  */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) lname ## _	/* Fortran Module Name */
#	define FMCALL(lname, uname) lname ## _	/* Fortran Module Call */
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#	define STR_PSTR(str)	char *str
#	define STR_PLEN(str)	, int Len ## str
#	define STR_PTR(str)	str
#	define STR_LEN(str)	Len ## str
#	define STR_SDEF(str)	char *str; int Len ## str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_SLEN(str, len) Len ## str = (int)len
#	define STR_FSTR(str)	(char *)str
#	define STR_FLEN(str)	, (int)Len ## str
#endif

#ifdef  VMS
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#       include <descrip.h>
#       define FCNAME(lname, uname) uname /* Fortran commons */
#	define FMPREFIX			  /* Fortran Module Prefix */
#       define FMNAME(lname, uname) uname /* Fortran Module Name */
#       define FMCALL(lname, uname) uname /* Fortran Module Call */
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#       define STR_PSTR(str)    struct dsc$descriptor_s *str
#       define STR_PLEN(str)
#       define STR_PTR(str)     str->dsc$a_pointer
#       define STR_LEN(str)     str->dsc$w_length
#       define STR_SDEF(str)    struct dsc$descriptor_s str ## Int = \
				  {0, DSC$K_DTYPE_T, DSC$K_CLASS_S, NULL}; \
	struct dsc$descriptor_s *str = &str ## Int
#       define STR_SPTR(str, ptr) str->dsc$a_pointer = ptr
#       define STR_SLEN(str, len) str->dsc$w_length = len
#       define STR_FSTR(str)    str
#       define STR_FLEN(str)
#endif

#ifdef	hitachi
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#	define FCNAME(lname, uname) lname 	/* Fortran commons  */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) lname 	/* Fortran Module Name */
#	define FMCALL(lname, uname) lname 	/* Fortran Module Call */
#	define FNC_CPTR(func)	** func
#	define FNC_PPTR(func)	** func
#	define FNC_FPTR(func)   &(func)
#	define STR_PSTR(str)	char *str
#	define STR_PLEN(str)	, int Len ## str
#	define STR_PTR(str)	str
#	define STR_LEN(str)	Len ## str
#	define STR_SDEF(str)	char *str; int Len ## str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_SLEN(str, len) Len ## str = (int)len
#	define STR_FSTR(str)	(char *)str
#	define STR_FLEN(str)	, (int)Len ## str
#endif

#if   	defined(__osf__) && !defined(VMS) && !defined(_WIN32)
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR long
#	define FCNAME(lname, uname) lname ## _	/* Fortran commons  */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) lname ## _	/* Fortran Module Name */
#	define FMCALL(lname, uname) lname ## _	/* Fortran Module Call */
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#	define STR_PSTR(str)	char *str
#	define STR_PLEN(str)	, int Len ## str
#	define STR_PTR(str)	str
#	define STR_LEN(str)	Len ## str
#	define STR_SDEF(str)	char *str; int Len ## str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_SLEN(str, len) Len ## str = (int)len
#	define STR_FSTR(str)	(char *)str
#	define STR_FLEN(str)	, (int)Len ## str
#endif

#ifdef	__ksr__    /* Kendall Square */
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#	define FCNAME(lname, uname) lname ## _	/* Fortran commons  */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) lname ## _	/* Fortran Module Name */
#	define FMCALL(lname, uname) lname ## _	/* Fortran Module Call */
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#	define STR_PSTR(str)	char *str
#	define STR_PLEN(str)	, int Len ## str
#	define STR_PTR(str)	str
#	define STR_LEN(str)	Len ## str
#	define STR_SDEF(str)	char *str; int Len ## str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_SLEN(str, len) Len ## str = (int)len
#	define STR_FSTR(str)	(char *)str
#	define STR_FLEN(str)	, (int)Len ## str
#endif

#ifdef	nec_ews
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#	define FCNAME(lname, uname) lname ## _	/* Fortran commons  */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) lname ## _	/* Fortran Module Name */
#	define FMCALL(lname, uname) lname ## _	/* Fortran Module Call */
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#	define STR_PSTR(str)	char *str
#	define STR_PLEN(str)	, int Len ## str
#	define STR_PTR(str)	str
#	define STR_LEN(str)	Len ## str
#	define STR_SDEF(str)	char *str; int Len ## str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_SLEN(str, len) Len ## str = (int)len
#	define STR_FSTR(str)	(char *)str
#	define STR_FLEN(str)	, (int)Len ## str
#endif

#ifdef	SX   /* NEC Supercomputer */
#ifndef mode_t
# define mode_t int
#endif
#ifndef pid_t
# define pid_t int
#endif
#       define L4NG long long
#       define UL4NG unsigned long long
#       define FTNPTR L4NG
#	define FCNAME(lname, uname) lname ## _	/* Fortran commons  */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) lname ## _	/* Fortran Module Name */
#	define FMCALL(lname, uname) lname ## _	/* Fortran Module Call */
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#	define STR_PSTR(str)	char *str
#	define STR_PLEN(str)	, int Len ## str
#	define STR_PTR(str)	str
#	define STR_LEN(str)	Len ## str
#	define STR_SDEF(str)	char *str; int Len ## str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_SLEN(str, len) Len ## str = (int)len
#	define STR_FSTR(str)	(char *)str
#	define STR_FLEN(str)	, (int)Len ## str
#endif

#ifdef	cdc
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#	define FCNAME(lname, uname) lname ## _	/* Fortran commons  */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) lname ## _	/* Fortran Module Name */
#	define FMCALL(lname, uname) lname ## _	/* Fortran Module Call */
#	define FNC_CPTR(func)	* func
#	define FNC_PPTR(func)	* func
#	define FNC_FPTR(func)   func
#	define STR_PSTR(str)	char *str
#	define STR_PLEN(str)	, int Len ## str
#	define STR_PTR(str)	str
#	define STR_LEN(str)	Len ## str
#	define STR_SDEF(str)	char *str; int Len ## str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_SLEN(str, len) Len ## str = (int)len
#	define STR_FSTR(str)	(char *)str
#	define STR_FLEN(str)	, (int)Len ## str
#endif

#ifdef	__uxp__   /* FUJITSU */
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#	define FCNAME(lname, uname) lname ## _	/* Fortran commons  */
#	define FMPREFIX				/* Fortran Module Prefix */
#	define FMNAME(lname, uname) lname ## _	/* Fortran Module Name */
#	define FMCALL(lname, uname) lname ## _	/* Fortran Module Call */
#	define FNC_CPTR(func)	** func
#	define FNC_PPTR(func)	** func
#	define FNC_FPTR(func)   &(func)
#	define STR_PSTR(str)	char *str
#	define STR_PLEN(str)	, int Len ## str
#	define STR_PTR(str)	str
#	define STR_LEN(str)	Len ## str
#	define STR_SDEF(str)	char *str; int Len ## str
#	define STR_SPTR(str, ptr) str = (char *)ptr
#	define STR_SLEN(str, len) Len ## str = (int)len
#	define STR_FSTR(str)	(char *)str
#	define STR_FLEN(str)	, (int)Len ## str
#endif

#ifdef  _WIN32   /* WINDOW_NT */
#       define L4NG long
#       define UL4NG unsigned long
#       define FTNPTR L4NG
#       define FCNAME(lname, uname) uname           /* Fortran Common Name */
#	define FMPREFIX	__stdcall		    /* Fortran Module Prefix */
#if defined(__NUTC__)
#       define FMNAME(lname, uname) lname  /* Fortran Module Name */
#else
#       define FMNAME(lname, uname) FMPREFIX uname  /* Fortran Module Name */
#endif
#       define FMCALL(lname, uname) uname           /* Fortran Module Call */
#       define FNC_CPTR(func)   * func
#       define FNC_PPTR(func)   FMPREFIX * func
#       define FNC_FPTR(func)   func
#       define STR_PSTR(str)    char *str, int Len ## str
#       define STR_PLEN(str)
#       define STR_PTR(str)     str
#       define STR_LEN(str)     Len ## str
#       define STR_SDEF(str)    char *str; int Len ## str
#       define STR_SPTR(str, ptr) str = (char *)ptr
#       define STR_SLEN(str, len) Len ## str = (int)len
#       define STR_FSTR(str)    (char *)str, (int)Len ## str
#       define STR_FLEN(str)
#endif /* _WIN32 */

#ifndef FMNAME2
#define FMNAME2(x,y) FMNAME(x,y)
#endif
#ifndef FMCALL2
#define FMCALL2(x,y) FMCALL(x,y)
#endif

#ifndef FMNAME
	*** You need to add information about this system in __FILE__ ***
#endif


/*************/
/* Datatypes */
/*************/

#define  L8NG    L4NG
#define  UL8NG   UL4NG

typedef	char		VCHARACTER;
#ifdef INTEGER_IS_LONG
typedef L4NG		VINTEGER;
#else
typedef int		VINTEGER;
#endif
typedef double		VREAL;
typedef float		VFLOAT;

#if 1

typedef VCHARACTER	* CHARACTER;
typedef	VINTEGER	* INTEGER;
typedef VREAL		* REAL;
#if !defined(_WIN32) || !defined(_WINDEF_)
typedef VFLOAT		* PFLOAT;
#endif

#else

typedef VCHARACTER	*const CHARACTER;
typedef	VINTEGER	*const INTEGER;
typedef VREAL		*const REAL;
#ifndef _WIN32
typedef VFLOAT		*const PFLOAT;
#endif

#endif

/**********/
/* Macros */
/**********/

/* PTR_TO_INT is a macro to convert from a pointer to an int */

#define PTR_TO_INT(ptr) (int)((FTNPTR)(ptr))

/**************************/
/* Fortran Boolean Values */
/**************************/
#define FORTRAN_FALSE  (VINTEGER)(0x00000000)

#if defined(_CRAY) || defined(CRAY_OS)
# define FORTRAN_TRUE   (VINTEGER)(-1)
#endif

#if defined(__hp9000s800) || defined(HP800_OS)
# define FORTRAN_TRUE   (VINTEGER)(1)
#endif

#ifndef FORTRAN_TRUE
# define FORTRAN_TRUE   (VINTEGER)(0xFFFFFFFF)
#endif


/**************/
/* Prototypes */
/**************/

/*****************************************************************************/
/*  tocstr - convert a fortran character string into a fortran integer array */
/*	     that has a trailing null character to look like a c-string	     */
/*									     */
/*  character	str			(in) Fortran chacter string	     */
/*  integer	icstr(*)		(out) Fortran integer array	     */
/*									     */
/*  notes:								     */
/*    1) Trailing blanks are removed, and leading/trailing '"' are also.     */
/*    2) To keep the trailing blanks, quote them '"'			     */
/*    3) It is a little faster to call this routine without any trailing     */
/*       blanks; ex. call tocstr (abc[1:notblk], icstr)                      */
/*****************************************************************************/

void FMNAME(tocstr,TOCSTR) (
STR_PSTR(str),				/* (in)  Fortran character string */
CHARACTER icstr				/* (out) Fortran integer array */
STR_PLEN(str)				/* (in)  Compiler passed len of str */
			   );

/*****************************************************************************/
/*  frcstr - convert a fortran integer array into a fortran character string */
/*									     */
/*  integer	icstr(*)		(in) Fortran integer array	     */
/*  character	str			(out) Fortran chacter string	     */
/*****************************************************************************/

void FMNAME(frcstr,FRCSTR) (
CHARACTER icstr,			/* (in)  Fortran integer array */
STR_PSTR(str)				/* (out) Fortran character string */
STR_PLEN(str)				/* (in)  Compiler passed len of str */
			   );

/*************************************************************************/
/*  fstr_to_cstr - convert a fortran character string into a c character */
/*		   string						 */
/*************************************************************************/

void fstr_to_cstr (
char *str,				/* (in) Pointer to character string */
int   ilen,				/* (in) Max length of str */
char *icstr				/* (out) C character string */
		  );

/*************************************************************************/
/*  cstr_to_fstr - convert a c character string into a fortran character */
/*		   string						 */
/*************************************************************************/

void cstr_to_fstr (
char *icstr,				/* (in) C character string */
int   ilen,				/* (in) Max length of str */
char *str				/* (out) Pointer to character string */
		  );


/******************************************************************************
  The following are convenience macros that help the C programmer call
  FORTRAN routines and FORTRAN callable C routines that require strings . . .
******************************************************************************/
#define CSTR(fString) cString_ ## fString

/******************************************************************************
  Use CSTR_DEFINE() to define both C and FORTRAN variables.

     fString  -  root name for FORTRAN variables
     cSize    -  dimension for C char cString[cSize]
******************************************************************************/
#define CSTR_DEFINE( fString, cSize )                   \
   char CSTR(fString)[ ((cSize)+1) ];                   \
   STR_SDEF( fString );                                 \
   STR_SPTR( fString, CSTR(fString) );                  \
   STR_SLEN( fString, cSize )

/******************************************************************************
  Use CSTR_ASSOC() when you already have a C string to pass to FORTRAN.

     fString  -  root name for FORTRAN variables
     cString  -  name of C (char *)
     cSize    -  strlen( cString )  when  (char *)cString
                 sizeof( cString )  when  (char *)cString  or  char cString[x]
******************************************************************************/
#define CSTR_ASSOC( fString, cString, cSize )           \
   STR_SDEF( fString );                                 \
   STR_SPTR( fString, cString );                        \
   STR_SLEN( fString, cSize )

#define CSTR_ASSOC_STRLEN( fString, cString )           \
        CSTR_ASSOC( fString, cString, strlen(cString) )

#define CSTR_ASSOC_SIZEOF( fString, cString )           \
        CSTR_ASSOC( fString, cString, sizeof(cString) )


#ifndef EXTERN_C
/*
 * Note:    __cplusplus is what's specified in the ARM sec. 16.3c
 */
#if defined(_LANGUAGE_C_PLUS_PLUS) || defined(__cplusplus)
#define EXTERN_C        extern "C"
#else
#define EXTERN_C        extern
#endif
#endif

#endif	/* FORTRAN_H */
