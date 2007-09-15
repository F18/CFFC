#ifndef _TYPEDEFINITION_CFRS_INCLUDED
#define _TYPEDEFINITION_CFRS_INCLUDED

#ifndef _Cell1D_INCLUDED
#include "Grid/Grid1D/Cell1D.h"
#endif //_Cell1D_INCLUDED

#ifndef _CELL2D_INCLUDED
#include "Grid/Grid2D/Cell2D.h"
#endif //_CELL2D_INCLUDED

#ifndef _CELL3D_INCLUDED
#include "Grid/Grid3D/Cell3D.h"
#endif //_CELL3D_INCLUDED

#ifndef _GAUSSIAN2D_STATE_INCLUDED
#include "CFD/Gaussian2DState.h"
#endif

#ifndef _EULER2D_STATE_INCLUDED
#include "CFD/Euler2DState.h"
#endif

#include "../../../src_2D/Utilities/TypeDefinition.h"

typedef enum {FALSE=0, TRUE=1} Bool;

// Geometry Traits;
template<SpaceType SpaceDimension>
class GeometryTraits;

template<>
class GeometryTraits<OneD> {
 public:
  typedef double NodeType;
  typedef double VectorType;
};

template<>
class GeometryTraits<TwoD> {
 public:
  typedef Node2D NodeType;
  typedef Vector2D VectorType;
};

template<>
class GeometryTraits<ThreeD> {
 public:
  typedef Node3D NodeType;
  typedef Vector3D VectorType;
};

template<>
class SolutionParameters<Gaussian2D_pState> {
 public:
  static const int NUM_OF_VARIABLES = 8; 
};

template<>
class SolutionParameters<Euler2D_pState> {
 public:
  static const int NUM_OF_VARIABLES = 4; 
};

#endif
