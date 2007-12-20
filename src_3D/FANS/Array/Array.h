#ifndef _ARRAY_INCLUDED
#define _ARRAY_INCLUDED

/*******************************************************************************
 *
 * Include this file to use the arrays even if only one type is desired
 *
 ******************************************************************************/

//--MPI support

#ifdef USE_MPI
#include "mpi.h"
#endif

//--Set debugging flags

#ifdef _ARRAY_DEBUG
#define _DYNAMIC_ARRAY_DEBUG
#define _STATIC_ARRAY_DEBUG
#define _REFERENCE_ARRAY_DEBUG
#endif

#ifdef _ARRAY_CHECK_BOUNDS
#define _DYNAMIC_ARRAY_CHECK_BOUNDS
#define _STATIC_ARRAY_CHECK_BOUNDS
#define _REFERENCE_ARRAY_CHECK_BOUNDS
#endif

#ifdef _DYNAMIC_ARRAY_DEBUG
#define _DYNAMIC_ARRAY_CHECK_ALLOCATED
#define _DYNAMIC_ARRAY_CHECK_BOUNDS
#endif

#ifdef _STATIC_ARRAY_DEBUG
#define _STATIC_ARRAY_CHECK_BOUNDS
#endif

#ifdef _REFERENCE_ARRAY_DEBUG
#define _REFERENCE_ARRAY_CHECK_LENGTH
#define _REFERENCE_ARRAY_CHECK_SET
#define _REFERENCE_ARRAY_CHECK_BOUNDS
#endif

//--Helper includes

#include "Math_Array1D.h"
#include "Array_traits.h"
#include "Array_policies.h"

//--Array classes

#include "Array_reference.h"
#include "Array_dynamic.h"
#include "Array_static.h"

#endif  // _ARRAY_INCLUDED
