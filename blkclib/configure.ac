#                                               -*- Autoconf -*-
#--------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
#--------------------------------------------------------------------

AC_PREREQ([2.59])
AC_INIT([Block Connectivity], [1.0], [sguzik@utias.utoronto.ca], [blkc.tar.bz2])
AC_CONFIG_SRCDIR([src/BlkC_block.h])

#--------------------------------------------------------------------
# System
#--------------------------------------------------------------------

UNAME=`uname -m`

#--------------------------------------------------------------------
# Library name
#--------------------------------------------------------------------

AC_SUBST([LIBNAME],[blkc])

#--------------------------------------------------------------------
# Source directories
#--------------------------------------------------------------------

srcdir="src"

#--------------------------------------------------------------------
# Get the C++ compiler
#--------------------------------------------------------------------

# Known compilers
if test -z "$CXX"; then
   case "$UNAME" in
      alpha) CXX=cxx;;
      ia64)  CXX=icpc;;
   esac
fi

# Test for compiler - if not found, use g++
if test "$CXX" = "cxx" || test "$CXX" = "icpc"; then
   AC_CHECK_PROG([CXX], [$CXX], [g++])
fi
if test -z "$CXX" || test "$CXX" = "g++"; then
   AC_PROG_CXX([$CXX])
else
   OBJEXT=o
fi

# Check for debug
AC_MSG_CHECKING([if debugging is enabled])
AC_ARG_ENABLE(debug,
   AC_HELP_STRING([--enable-debug], [build debugging version (default=no)]),,
   enableval=no)
if test "$enableval" = "yes"; then
   DEBUG=yes
else
   DEBUG=no
fi
AC_MSG_RESULT($DEBUG)

# Set debug, optimization, and regular flags for known compilers
case "$CXX" in
   g++)  DEBUGFLAGS="-O -g"
         OPTFLAGS=-O3
         REGFLAGS=""
         DEPFLAGS=-MM ;;
   cxx)  DEBUGFLAGS=-g
         OPTFLAGS=-fast
         REGFLAGS="-assume accuracy_sensitive -std strict_ansi"
         DEPFLAGS=-M
         LDFLAGS=-lm
         REPOSOBJ="cxx_repository/*.o" ;;
   icpc) DEBUGFLAGS="-O -g"
         OPTFLAGS=-O3
         REGFLAGS="-mp -w"
         DEPFLAGS=-MM ;;
esac

# Set compilation flags
if test "$DEBUG" = "yes"; then
   CXXFLAGS="$DEBUGFLAGS $REGFLAGS"
else
   CXXFLAGS="$OPTFLAGS $REGFLAGS"
fi

#--------------------------------------------------------------------
# Program characteristics
#--------------------------------------------------------------------

BLKCDEF=""

# Size of size_t
AC_CHECK_SIZEOF([size_t])
if test $ac_cv_sizeof_size_t != 4; then
  if test $ac_cv_sizeof_size_t != 8; then
    AC_MSG_WARN([Unable to detect size of size_t - assuming 64 bit for FNV hashing.])
  fi
  BLKCDEF="-DHAVE_64BIT_SIZE_T $BLKCDEF"
fi

# Check for hash map
AC_MSG_CHECKING([for type of map])
AC_ARG_ENABLE(hashmap,
   AC_HELP_STRING([--enable-hashmap], [use hash map in place of maps (default=no)]),,
   enableval=no)
if test "$enableval" = "yes"; then
   BLKCDEF="-DHAVE_HASH_MAP $BLKCDEF"
   hashmap="hash map"
else
   hashmap="std::map"
fi
AC_MSG_RESULT($hashmap)
if test "$hashmap" = "hash map"; then
   AC_MSG_NOTICE([Set hash map in src/hashmap.h])
fi

#--------------------------------------------------------------------
# Standard checks for required programs
#--------------------------------------------------------------------

AC_PROG_INSTALL

#--------------------------------------------------------------------
# Additional programs
#--------------------------------------------------------------------

# Library archiver
AC_ARG_VAR([AR], [Library archiver])
AC_CHECK_PROG([AR], [ar], [ar rcs])

# Stream editor
AC_ARG_VAR([SED], [Stream editor])
AC_CHECK_PROG([SED], [sed], [sed], [])

#--------------------------------------------------------------------
# Extensions
#--------------------------------------------------------------------

AC_SUBST([LIBEXT],[a])

#--------------------------------------------------------------------
# Install directories
#--------------------------------------------------------------------

if test -z $prefix || test $prefix = NONE; then
  prefix=/usr/local
fi
exec_prefix=$prefix

AC_MSG_CHECKING([where to install library])
eval "INSTLIBDIR=$libdir"
AC_MSG_RESULT([$INSTLIBDIR])

AC_MSG_CHECKING([where to install headers])
eval "INSTINCLUDEDIR=$includedir"
AC_MSG_RESULT([$INSTINCLUDEDIR])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

AC_SUBST([DEPFLAGS])
AC_SUBST([REPOSOBJ])
AC_SUBST([BLKCDEF])
AC_SUBST([INSTLIBDIR])
AC_SUBST([INSTINCLUDEDIR])

AC_CONFIG_FILES([make.defs])
AC_OUTPUT
