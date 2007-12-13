/**********************************************************************
 * Interface2D.h: Header file defining 2D interfaces.                 *
 **********************************************************************/

#ifndef _INTERFACE2D_INCLUDED
#define _INTERFACE2D_INCLUDED

// Include required C++ libraries.

#include <cassert>
#include <cstdlib>

using namespace std;

// Include math macro file.

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

// Include Spline2D header.

#ifndef _SPLINE2D_INCLUDED
#include "../Math/Spline2D.h"
#endif // _SPLINE2D_INCLUDED

// Include Vector2D header.

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif // _VECTOR2D_INCLUDED

// Include CFD header.

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

// Include Polygon header.

#ifndef _POLYGON_INCLUDED
#include "../Math/Polygon.h"
#endif // _POLYGON_INCLUDED

// Interface flag definitions.

#define INTERFACE_LINE                            500
#define INTERFACE_CIRCLE                          501
#define INTERFACE_ELLIPSE                         502
#define INTERFACE_SQUARE                          503
#define INTERFACE_RECTANGLE                       504
#define INTERFACE_STRAIGHT                        505
#define INTERFACE_ROCKET_PROPELLANT_GRAIN         506
#define INTERFACE_NACA0012_AEROFOIL               507
#define INTERFACE_NACA0015_AEROFOIL               508
#define INTERFACE_NASA_ROTOR_37                   509
#define INTERFACE_NASA_ROTOR_67                   510
#define INTERFACE_ZALESAK                         511
#define INTERFACE_RINGLEB                         512
#define INTERFACE_BUMP                            513
#define INTERFACE_FLAT_PLATE                      514
#define INTERFACE_USER_SPECIFIED                  515
#define INTERFACE_UNION                           516
#define INTERFACE_STAR                            517
#define INTERFACE_RESTART                         599

#define INTERFACE_BC_REFLECTION               BC_REFLECTION
#define INTERFACE_BC_ABSORPTION               BC_ABSORPTION
#define INTERFACE_BC_BURNING_SURFACE          BC_BURNING_SURFACE
#define INTERFACE_BC_WALL_VISCOUS_HEATFLUX    BC_WALL_VISCOUS_HEATFLUX
#define INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL  BC_WALL_VISCOUS_ISOTHERMAL
#define INTERFACE_BC_RINGLEB                  BC_RINGLEB_FLOW
#define INTERFACE_BC_ADIABATIC_WALL           BC_ADIABATIC_WALL
#define INTERFACE_BC_FLAT_PLATE                 99990      
#define INTERFACE_BC_DETERMINE                  99991
#define INTERFACE_BC_MIXED                      99992

#define INTERFACE_MOTION_STATIONARY               700
#define INTERFACE_MOTION_CONSTANT                 701
#define INTERFACE_MOTION_EXPAND                   702
#define INTERFACE_MOTION_UNIFORM                  703
#define INTERFACE_MOTION_TRANSLATE                704
#define INTERFACE_MOTION_STRETCH                  705
#define INTERFACE_MOTION_ROTATE                   706
#define INTERFACE_MOTION_BURNING_SURFACE          707
#define INTERFACE_MOTION_MOMENTUM_TRANSFER        708
#define INTERFACE_MOTION_MIXED                    709
#define INTERFACE_MOTION_COMPUTE                  710 // Not necessary?
#define INTERFACE_MOTION_LEVELSET                 720
#define INTERFACE_MOTION_LEVELSET_STATIONARY      721
#define INTERFACE_MOTION_LEVELSET_EXPAND          722
#define INTERFACE_MOTION_LEVELSET_STRETCH         723
#define INTERFACE_MOTION_LEVELSET_BULKFLOW        724

#define INTERFACE_BULKFLOWFIELD_NONE              800
#define INTERFACE_BULKFLOWFIELD_UNIFORM           801
#define INTERFACE_BULKFLOWFIELD_SWIRL             802

#define SOLUTION_REDISTRIBUTION_NONE              900
#define SOLUTION_REDISTRIBUTION_INJECTION         901
#define SOLUTION_REDISTIBUTION_WEIGHTED_INJECTION 902
#define SOLUTION_REDISTIBUTION_NEIGHBOUR_AVERAGE  903

/*!
 * Class: Interface2D
 *
 * @brief Class definition of the embedded boundaries (interfaces) used 
 * in the level set code or the embedded boundary code.
 *
 * \verbatim
 * Member functions
 *   Type        -- Interface type.
 *   Length      -- Interface characteristic length.
 *   BC_Type     -- Interface boundary condition type.
 *   Motion      -- Interface motion type.
 *   Speed       -- Interface characteristic speed.
 *   Xref        -- Interface reference point.
 *   Spline      -- Interface spline.
 *   Xmin        -- Interface bounding box minimum point.
 *   Xmax        -- Interface bounding box maximum point.
 *   F_Type      -- Velocity function type at interface nodes.
 *   F           -- Velocity function at interface nodes.
 *   allocate    -- Allocate memory for interface.
 *   deallocate  -- Deallocate memory for interface.
 *
 * Member operators
 *      I -- an interface.
 *
 * I = I;
 * cout << I; (output function)
 * cin  >> I; (input function)
 *\endverbatim
 */
class Interface2D{
 private:
 public:
  int             Type; //!< Interface type.
  double       Length1; //!< Interface characteristic length.
  double       Length2; //!< Interface characteristic length.
  int          BC_Type; //!< Interface boundary condition type.
  int           Motion; //!< Interface motion type.
  Vector2D       Speed; //!< Interface characteristic speed.
  Vector2D        Xref; //!< Interface reference point.
  Spline2D      Spline; //!< Interface spline.
  Vector2D        Xmin; //!< Interface bounding box minimum point.
  Vector2D        Xmax; //!< Interface bounding box maximum point.
  int          *F_Type; //!< Velocity function type defined at interface nodes.
  Vector2D          *F; //!< Velocity function defined at interface nodes.

  //! Creation constructor.
  Interface2D(void) {
    Type    = INTERFACE_LINE;
    Length1 = ZERO;
    Length2 = ZERO;
    BC_Type = INTERFACE_BC_REFLECTION;
    Motion  = INTERFACE_MOTION_STATIONARY;
    Speed   = Vector2D_ZERO;
    Xref    = Vector2D_ZERO;
    Spline.deallocate();
    Xmin    = Vector2D_ZERO;
    Xmax    = Vector2D_ZERO;
    F_Type  = NULL;
    F       = NULL;
  }

  //! Copy constructor.
  Interface2D(const Interface2D &I) {
    deallocate();
    Type    = I.Type;
    Length1 = I.Length1;
    Length2 = I.Length2;
    BC_Type = I.BC_Type;
    Motion  = I.Motion;
    Speed   = I.Speed;
    Xref    = I.Xref;
    Spline  = I.Spline;
    Xmin    = I.Xmin;
    Xmax    = I.Xmax;
    Initialize_Velocity_Function();
    for (int n = 0; n < Spline.np; n++) {
      F_Type[n] = I.F_Type[n];
      F[n]      = I.F[n];
    }
  }
  
  //! Destructor.
  ~Interface2D(void) {
    Spline.deallocate();
    if (F_Type != NULL) { delete []F_Type; F_Type = NULL; }
    if (F != NULL) { delete []F; F = NULL; }
  }

  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int npts) {
    deallocate();
    Spline.allocate(npts);
    Initialize_Velocity_Function();
  }

  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void) {
    Spline.deallocate();
    if (F_Type != NULL) { delete []F_Type; F_Type = NULL; }
    if (F != NULL) { delete []F; F = NULL; }
  }

  //! Copy interface.
  void Copy(Interface2D &I);

  //! Copy interface.
  void Copy(const Interface2D &I);

  //! Initialize the front speed array.
  void Initialize_Velocity_Function(void);

  //! Initialize the front speed array.
  void Initialize_Velocity_Function(const int &ftype);

  //! Set the interface velocity function type according to the motion type.
  void Set_Velocity_Function_Type(void);

  //! Set the interface velocity function according to the motion type.
  void Set_Velocity_Function(const double &val);

  //! Reset the interface velocity function into the normal direction.
  void Set_Normal_Velocity_Function(void);

  //! Set the interface velocity function to zero.
  void Zero_Velocity_Function(void);

  //! Sort the interface into counter-clockwise order.
  void Sort(void);

  //! Translate the interface.
  void Translate(const Vector2D &X);  

  //! Scale the interface.
  void Scale(const double &Scaling_Factor);

  //! Rotate the interface.
  void Rotate(const double &Angle);

  //! Calculate the area of the interface.
  double Area(void);

  //! Calculate the centroid of the interface.
  void Centroid(void);

  //! Determine the bounding box of the interface.
  void BoundingBox(void);

  //! Determine an extended bounding box of the interface.
  void BoundingBox(const double &val);

  //! Determine if the given point is in the bounding box of the interface.
  int Point_In_Bounding_Box(const Vector2D &X);

  //! Determine if the given bounding box intersects with the current bounding box.
  int Bounding_Box_Intersection(const Vector2D &Ymin, const Vector2D &Ymax);

  //! Determine if the given point is inside the interface.
  int Point_In_Interface(const Vector2D &Xt);

  //! Return the boundary condition at a location X on the interface spline.
  int Determine_Interface_BC_Type(const Vector2D &X);

  //! Return the velocity at a location X on the interface spline.
  Vector2D Determine_Interface_Velocity(const Vector2D &X,
					const double &time);

  //! Return normal velocity at the specified interface point.
  double Fn(const int &np, const Vector2D &norm_dir) const;
  double Fn(const int &np) const;

  //! Set the normal velocity at the specified interface point.
  void Fn(const int &np, const double &fn);

  //! Return interface normal at the specified interface point.
  Vector2D normal(const int &np) const;

  //! Return the maximum speed.
  double max_speed(void);

  //! Return the temperature at a location X on the interface spline.
  double Determine_Interface_Temperature(const Vector2D &X);

  //! Create Zalesak's disk interface spline.
  void Zalesak(const Vector2D &Origin,
	       const double &Radius);

  //! Create multi-point star interface spline.
  void Star(const Vector2D &Origin,
	    const double &Radius,
	    const int &Num_Ext_Pts);

  //! Create Ringleb's flow interface spline.
  void Ringleb(const double &Inner_Streamline_Number,
	       const double &Outer_Streamline_Number,
	       const double &Isotach_Line);

  //! Create the bump channel flow interface spline.
  void Bump_Channel_Flow(void);

  //! Create the flat plate interface spline.
  void Flat_Plate(void);

  //! Create from solution block.
  template <class Quad_Grid_Block> void Create_from_Solution_Block(Quad_Grid_Block &Grid);

  //! Determine the union of the interface with the given interface.
  int Interface_Union(Interface2D &I2);

  //! Determine the intersection of the interface with the given interface.
  void Interface_Intersection(Interface2D &I2);

  int Check_Interface_Intersection(Interface2D &I2);

  //! Assignment Operator.
  Interface2D& operator =(const Interface2D &I);

  //@{ Input-output operator.
  friend ostream &operator << (ostream &out_file, const Interface2D &I);
  friend istream &operator >> (istream &in_file, Interface2D &I);
  //@}

};

/*!
 * Class Interface2D_List
 *
 * @brief List of 2D interfaces.
 *
 * \verbatim
 * Member functions
 *   Ni          -- Number of interfaces in the list.
 *   Interface   -- List (array) of interfaces.
 *   relations   -- List of associated interface components.
 *   allocate    -- Allocate memory for the list.
 *   deallocate  -- Deallocate memory for the list.
 *
 * Member operators
 *      I -- an interface list.
 *
 * I = I;
 * \endverbatim
 */
class Interface2D_List {
 private:
 public:
  int                     Ni; //!< Number of interfaces in the list.
  Interface2D     *Interface; //!< List of interfaces.
  LinkedList<int> *relations; //!< List of associated interface components.

  //@{ @name Constructors
  //! Creation constructor.
  Interface2D_List(void) {
    Ni = 0;
    Interface = NULL;
    relations = NULL;
  }
  //! Copy constructor.
  Interface2D_List(const Interface2D_List &List) {
    deallocate();
    allocate(List.Ni);
    //for (int n = 0; n < Ni; n++) Interface[n] = List.Interface[n];
    for (int n = 0; n < Ni; n++) Interface[n].Copy(List.Interface[n]);
    for (int n = 0; n < Ni; n++) relations[n] = List.relations[n];
  }
  //! Destructor.
  ~Interface2D_List(void) {
    if (Interface != NULL) {
      for (int n = 0; n < Ni; n++) Interface[n].deallocate();
      delete []Interface; Interface = NULL;
    }
    if (relations != NULL) {
      for (int n = 0; n < Ni; n++) relations[n].deallocate();
      delete []relations; relations = NULL;
    }
    Ni = 0;
  }
  //@}

  //@{ Allocation and deallocation functions.
  //! Allocation function.
  void allocate(const int n) {
    deallocate();
    Ni = n;
    Interface = new Interface2D[Ni];
    relations = new LinkedList<int>[Ni];
  }
  //! Deallocation function.
  void deallocate(void) {
    if (Interface != NULL) {
      for (int n = 0; n < Ni; n++) Interface[n].deallocate();
      delete[]Interface; Interface = NULL;
    }
    if (relations != NULL) {
      for (int n = 0; n < Ni; n++) relations[n].deallocate();
      delete []relations; relations = NULL;
    }
    Ni = 0;
  }
  //@}

  //@{ Interface manipulation functions.

  //! Copy interface list.
  void Copy(Interface2D_List &List);

  //! Copy interface list.
  void Copy(const Interface2D_List &List);

  //! Construct interface union list from component list.
  int Construct_Union_List(Interface2D_List &Component_List);

  //@}

  //! Return the maximum interface speed.
  double max_speed(void);

  //@{ @name Index operator.
  Interface2D &operator[](int index) {
    assert(index >= 1 && index <= Ni);
    return Interface[index-1];
  }  
  const Interface2D &operator[](int index) const {
    assert(index >= 1 && index <= Ni);
    return Interface[index-1];
  }
  //@}

  //@{ Assignment operator.
  Interface2D_List& operator =(const Interface2D_List &List) {
    if (Ni != List.Ni) { deallocate(); allocate(List.Ni); }
    for (int n = 0; n < Ni; n++) Interface[n] = List.Interface[n];
    for (int n = 0; n < Ni; n++) relations[n] = List.relations[n];
    return *this;
  }
  //@}

  //@{ Input-output operator.
  friend ostream &operator << (ostream &out_file, const Interface2D_List &List);
  friend istream &operator >> (istream &in_file, Interface2D_List &List);
  //@}

};

/**********************************************************************
 * Interface2D::Define member functions.                              *
 **********************************************************************/

/**********************************************************************
 * Interface2D::Copy -- Copy interface.                               *
 **********************************************************************/
inline void Interface2D::Copy(Interface2D &I) {
  deallocate();
  Type = I.Type;
  Length1 = I.Length1;
  Length2 = I.Length2;
  BC_Type = I.BC_Type;
  Motion = I.Motion;
  Speed = I.Speed;
  Xref = I.Xref;
  Xmin = I.Xmin;
  Xmax = I.Xmax;
  Copy_Spline(Spline,I.Spline);
  Initialize_Velocity_Function();
  for (int np = 0; np < Spline.np; np++) {
    F_Type[np] = I.F_Type[np];
    F[np] = I.F[np];
  }
}

inline void Interface2D::Copy(const Interface2D &I) {
  deallocate();
  Type = I.Type;
  Length1 = I.Length1;
  Length2 = I.Length2;
  BC_Type = I.BC_Type;
  Motion = I.Motion;
  Speed = I.Speed;
  Xref = I.Xref;
  Xmin = I.Xmin;
  Xmax = I.Xmax;
  Copy_Spline(Spline,I.Spline);
  Initialize_Velocity_Function();
  for (int np = 0; np < Spline.np; np++) {
    F_Type[np] = I.F_Type[np];
    F[np] = I.F[np];
  }
}

/**********************************************************************
 * Interface2D::Initialize_Velocity_Function -- Initialize the        *
 *                                              interface velocity    *
 *                                              type and function.    *
 **********************************************************************/
inline void Interface2D::Initialize_Velocity_Function(void) {
  if (F_Type != NULL) { delete []F_Type; F_Type = NULL; }
  if (F != NULL) { delete []F; F = NULL; }
  F_Type = new int[Spline.np];
  F = new Vector2D[Spline.np];
  for (int np = 0; np < Spline.np; np++) {
    F_Type[np] = INTERFACE_MOTION_STATIONARY;
    F[np] = Vector2D_ZERO;
  }
}

/**********************************************************************
 * Interface2D::Initialize_Velocity_Function -- Initialize the        *
 *                                              interface velocity    *
 *                                              type and function.    *
 **********************************************************************/
inline void Interface2D::Initialize_Velocity_Function(const int &ftype) {
  if (F_Type != NULL) { delete []F_Type; F_Type = NULL; }
  if (F != NULL) { delete []F; F = NULL; }
  F_Type = new int[Spline.np];
  F = new Vector2D[Spline.np];
  for (int np = 0; np < Spline.np; np++) {
    F_Type[np] = ftype;
    F[np] = Vector2D_ZERO;
  }
}

/**********************************************************************
 * Interface2D::Set_Velocity_Function_Type -- Set the interface       *
 *                                            velocity function type  *
 *                                            according to the motion *
 *                                            type.                   *
 **********************************************************************/
inline void Interface2D::Set_Velocity_Function_Type(void) {
  for (int np = 0; np < Spline.np; np++) F_Type[np] = Motion;
}

/**********************************************************************
 * Interface2D::Set_Velocity_Function -- Set the interface velocity   *
 *                                       function according to the    *
 *                                       motion type.                 *
 **********************************************************************/
inline void Interface2D::Set_Velocity_Function(const double &val) {
  switch(Motion) {
  case INTERFACE_MOTION_STATIONARY :
    // The embedded boundary is stationary.
    break;
  case INTERFACE_MOTION_CONSTANT :
  case INTERFACE_MOTION_EXPAND :
    // The embedded boundary has a constant speed in the radial direction.
    for (int np = 0; np < Spline.np; np++)
      F[np] = Speed.x*(Spline.Xp[np] - Xref)/abs(Spline.Xp[np] - Xref);
    break;
  case INTERFACE_MOTION_UNIFORM :
  case INTERFACE_MOTION_TRANSLATE :
    // The embedded boundary has a constant velocity.
    for (int np = 0; np < Spline.np; np++) F[np] = Speed;
    break;
  case INTERFACE_MOTION_STRETCH :
    if (Type == INTERFACE_LINE) {
    } else if (Type == INTERFACE_CIRCLE) {
      for (int np = 0; np < Spline.np; np++)
	F[np].x = Speed.x*(ONE + fabs(Spline.Xp[np].y/abs(Spline.Xp[np])));
    } else if (Type == INTERFACE_SQUARE) {
      for (int np = 0; np < Spline.np; np++) {
	if (Spline.Xp[np].x <= ZERO) {
	  F[np].x = ONE;
	} else if (Spline.Xp[np].x > ZERO) {
	  F[np].x = TWO;
	}
      }
    }
    break;
  case INTERFACE_MOTION_ROTATE :
    // The interface has a constant rotational velocity.
    switch(Type) {
    case INTERFACE_CIRCLE :
      for (int np = 0; np < Spline.np; np++) {
	F[np] = Vector2D((TWO*PI*Speed.x/Speed.y)*cos(TWO*PI*val/Speed.y),ZERO);
      }
      break;
    case INTERFACE_ELLIPSE :
      break;
    case INTERFACE_SQUARE :
    case INTERFACE_USER_SPECIFIED :
    case INTERFACE_ZALESAK :
      for (int np = 0; np < Spline.np; np++) {
	F[np].x = -Speed.x*Spline.Xp[np].y;
	F[np].y =  Speed.x*Spline.Xp[np].x;
      }
      break;
    case INTERFACE_NACA0012_AEROFOIL :
      for (int np = 0; np < Spline.np; np++) {
// 	F[np].x =  (2.51/(2.0*PI*62.5))*cos(2.0*PI*62.5*(val))*(Spline.Xp[np].y - Xref.y);
// 	F[np].y = -(2.51/(2.0*PI*62.5))*cos(2.0*PI*62.5*(val))*(Spline.Xp[np].x - Xref.x);
	F[np].x =  (2.0*PI*62.5*(PI*2.51/180.0))*cos(2.0*PI*62.5*(val))*(Spline.Xp[np].y - Xref.y);
	F[np].y = -(2.0*PI*62.5*(PI*2.51/180.0))*cos(2.0*PI*62.5*(val))*(Spline.Xp[np].x - Xref.x);
      }
      break;
    case INTERFACE_NACA0015_AEROFOIL :
      for (int np = 0; np < Spline.np; np++) {
	F[np].x =  (2.0*PI*62.5*(PI*2.51/180.0))*cos(2.0*PI*62.5*(val))*(Spline.Xp[np].y - Xref.y);
	F[np].y = -(2.0*PI*62.5*(PI*2.51/180.0))*cos(2.0*PI*62.5*(val))*(Spline.Xp[np].x - Xref.x);
      }
      break;
    };
    break;
  case INTERFACE_MOTION_BURNING_SURFACE :
  case INTERFACE_MOTION_MOMENTUM_TRANSFER :
  case INTERFACE_MOTION_MIXED :
    // The embedded boundary velocity must be externally computed.
    break;
  case INTERFACE_MOTION_LEVELSET_STATIONARY :
    // The embedded boundary is stationary.
    break;
  case INTERFACE_MOTION_LEVELSET_EXPAND :
    // The embedded boundary has a constant speed in the radial direction.
    for (int np = 0; np < Spline.np; np++)
      F[np] = Speed;//*(Spline.Xp[np] - Xref)/abs(Spline.Xp[np] - Xref);
    break;
  case INTERFACE_MOTION_LEVELSET_STRETCH :
    if (Type == INTERFACE_LINE) {
    } else if (Type == INTERFACE_CIRCLE) {
      //for (int np = 0; np < Spline.np; np++)
      //F[np].x = Speed.x*(ONE + fabs(Spline.Xp[np].y/abs(Spline.Xp[np])));
      for (int np = 0; np < Spline.np; np++) {
 	F[np].x = Speed.x*Spline.Xp[np].x/Length1;
      }
    } else if (Type == INTERFACE_SQUARE) {
      for (int np = 0; np < Spline.np; np++) {
	if (Spline.Xp[np].x <= ZERO) {
	  F[np].x = ONE;
	} else if (Spline.Xp[np].x > ZERO) {
	  F[np].x = TWO;
	}
      }
    }
    break;
  case INTERFACE_MOTION_LEVELSET_BULKFLOW :
    // The embedded boundary velocity must be externally computed.
    break;
  case INTERFACE_MOTION_LEVELSET :
    // The embedded boundary velocity is extracted from the level set
    // solution.
    break;
  };
}

/**********************************************************************
 * Interface2D::Set_Normal_Velocity_Function -- Reset the interface   *
 *                                              velocity function     *
 *                                              into the normal       *
 *                                              direction.            *
 **********************************************************************/
inline void Interface2D::Set_Normal_Velocity_Function(void) {
  double fn;
  for (int np = 0; np < Spline.np; np++) {
    fn = abs(F[np].x);
    Fn(np,fn);
  }
}

/**********************************************************************
 * Interface2D::Zero_Velocity_Function -- Set the interface velocity  *
 *                                        function to zero.           *
 **********************************************************************/
inline void Interface2D::Zero_Velocity_Function(void) {
  for (int np = 0; np < Spline.np; np++) F[np] = Vector2D_ZERO;
}

/**********************************************************************
 * Interface2D::Sort -- Sort the interface into counter-clockwise     *
 *                      order.                                        *
 **********************************************************************/
inline void Interface2D::Sort(void) {
  int sum;
  double xprdct;
  Vector2D Xp, Fp;
  Interface2D Itemp;
  if (Spline.np > 2) {
    sum = 0;
    for (int n = 0; n < Spline.np-2; n++) {
      xprdct = (Spline.Xp[n]-Spline.Xp[n+1])^(Spline.Xp[n+1]-Spline.Xp[n+2]);
      if (xprdct < ZERO) sum--;
      else sum++;
    }
    if (sum == 2-Spline.np) {
      for (int np = 0; np < Spline.np/2; np++) {
	Xp = Spline.Xp[np];
	Spline.Xp[np] = Spline.Xp[Spline.np-np-1];
	Spline.Xp[Spline.np-np-1] = Xp;
	Fp = F[np];
	F[np] = F[Spline.np-np-1];
	F[Spline.np-np-1] = Fp;
      }
//       Itemp.Copy(*this);
//       for (int np = 0; np < Spline.np; np++) {
// 	Spline.Xp[np] = Itemp.Spline.Xp[Spline.np-np-1];
// 	Spline.tp[np] = Itemp.Spline.tp[Spline.np-np-1];
// 	Spline.bc[np] = Itemp.Spline.bc[Spline.np-np-1];
//  	F[np] = Itemp.F[Spline.np-np-1];
//       }
      Spline.pathlength();
    }
  }
}

/**********************************************************************
 * Interface2D::Translate -- Translates the interface.                *
 **********************************************************************/
inline void Interface2D::Translate(const Vector2D &X) {
  Translate_Spline(Spline,X);
  Xref += X;
  BoundingBox();
}

/**********************************************************************
 * Interface2D::Scale -- Scales the interface.                        *
 **********************************************************************/
inline void Interface2D::Scale(const double &Scaling_Factor) {
  Translate_Spline(Spline,-Xref);
  Scale_Spline(Spline,Scaling_Factor);
  Translate_Spline(Spline,Xref);
  BoundingBox();
}

/**********************************************************************
 * Interface2D::Rotate -- Rotates the interface.                      *
 **********************************************************************/
inline void Interface2D::Rotate(const double &Angle) {
  Translate_Spline(Spline,-Xref);
  Rotate_Spline(Spline,Angle);
  Translate_Spline(Spline,Xref);
  BoundingBox();
}

/**********************************************************************
 * Interface2D::Create_from_Solution_Block -- Creates an interface    *
 *                                            from the domain of a    *  
 *                                            solution block.         *
 **********************************************************************/
template <class Quad_Grid_Block>
inline void Interface2D::Create_from_Solution_Block(Quad_Grid_Block &Grid) {
  int n, npts;
  if (Spline.np != 0) deallocate();
  // Set spline type.
  Spline.settype(SPLINE2D_LINEAR);
  // Allocate memory for the interface spline.
  Spline.allocate(2*Grid.NNi + 2*Grid.NNj - 3);
  for (int i = 0; i < Grid.NNi; i++) {
    Spline.Xp[i] = Grid.Node[i][0].X;
    Spline.tp[i] = SPLINE2D_POINT_NORMAL;
    Spline.bc[i] = BC_NONE;
  }
  npts = Grid.NNi-1;
  for (int j = 1; j < Grid.NNj; j++) {
    Spline.Xp[npts+j] = Grid.Node[Grid.NNi-1][j].X;
    Spline.tp[npts+j] = SPLINE2D_POINT_NORMAL;
    Spline.bc[npts+j] = BC_NONE;
  }
  npts += Grid.NNj-1;
  for (int i = Grid.NNi-1; i >= 0; i--) {
    n = Grid.NNi-1 - i;
    Spline.Xp[npts+n] = Grid.Node[i][Grid.NNj-1].X;
    Spline.tp[npts+n] = SPLINE2D_POINT_NORMAL;
    Spline.bc[npts+n] = BC_NONE;
  }
  npts += Grid.NNi-1;
  for (int j = Grid.NNj-1; j >= 0; j--) {
    n = Grid.NNj-1 - j;
    Spline.Xp[npts+n] = Grid.Node[0][j].X;
    Spline.tp[npts+n] = SPLINE2D_POINT_NORMAL;
    Spline.bc[npts+n] = BC_NONE;
  }
  // Calculate the spline pathlengths.
  Spline.pathlength();
  // Initialize the velocity function.
  Initialize_Velocity_Function();
  // Determine the bounding box.
  BoundingBox();
  // Determine the reeference point/centroid.
  Centroid();
}

/**********************************************************************
 * Interface2D::Assignment operator.                                  *
 **********************************************************************/
inline Interface2D& Interface2D::operator =(const Interface2D &I) {
  deallocate();
  Type = I.Type;
  Length1 = I.Length1;
  Length2 = I.Length2;
  BC_Type = I.BC_Type;
  Motion = I.Motion;
  Speed = I.Speed;
  Xref = I.Xref;
  Xmin = I.Xmin;
  Xmax = I.Xmax;
  Copy_Spline(Spline,I.Spline);
  Initialize_Velocity_Function();
  for (int np = 0; np < Spline.np; np++) {
    F_Type[np] = I.F_Type[np];
    F[np] = I.F[np];
  }
  return *this;
}

/**********************************************************************
 * Interface2D -- Output operator.                                    *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const Interface2D &I) {
  out_file << "Interface Type: " << I.Type << endl;
  out_file << "Interface characteristic length: " << I.Length1 << endl;
  out_file << "Interface characteristic length: " << I.Length2 << endl;
  out_file << "Interface boundary condition type: " << I.BC_Type << endl;
  out_file << "Interface motion type: " << I.Motion << endl;
  out_file << "Interface characteristic speed: " << I.Speed << endl;
  out_file << "Interface reference point: " << I.Xref << endl;
  out_file << "Interface bounding box minimum point: " << I.Xmin << endl;
  out_file << "Interface bounding box maximum point: " << I.Xmax << endl;
  out_file << "Interface spline, number of points: " << I.Spline.np << endl;
  out_file << "Interface spline, type: " << I.Spline.type << endl;
  out_file << "Spline pt  |  Spline tp  |  Spline bc" << endl;
  out_file << I.Spline;
  out_file << "At spline points:  F_Type  |  F" << endl;
  for (int np = 0; np < I.Spline.np; np++)
    out_file << I.F_Type[np] << I.F[np] << endl;
  return out_file;
}

/**********************************************************************
 * Interface2D -- Input operator.                                     *
 **********************************************************************/
inline istream &operator >> (istream &in_file, Interface2D &I) {
  int nsp, stype, f;
  in_file.setf(ios::skipws);
  in_file >> I.Type;
  in_file >> I.Length1;
  in_file >> I.Length2;
  in_file >> I.BC_Type;
  in_file >> I.Motion;
  in_file.unsetf(ios::skipws);
  in_file >> I.Speed;
  in_file >> I.Xref;
  in_file >> I.Xmin;
  in_file >> I.Xmax;
  in_file.setf(ios::skipws);
  in_file >> nsp;
  in_file.unsetf(ios::skipws);
  I.Spline.allocate(nsp);
  in_file.setf(ios::skipws);
  in_file >> stype;
  in_file.unsetf(ios::skipws);
  I.Spline.settype(stype);
  in_file >> I.Spline;
  I.Spline.pathlength();
  I.Initialize_Velocity_Function(I.Motion);
  in_file.setf(ios::skipws);
  for (int np = 0; np < nsp; np++) in_file >> I.F_Type[np] >> I.F[np].x >> I.F[np].y;
  in_file.unsetf(ios::skipws);
  I.Set_Velocity_Function(ZERO);
  return in_file;
}

/**********************************************************************
 * Interface2D_List::Define member functions.                         *
 **********************************************************************/

/**********************************************************************
 * Interface2D_List::Copy -- Copy interface list.                     *
 **********************************************************************/
inline void Interface2D_List::Copy(Interface2D_List &List) {
  allocate(List.Ni);
  for (int n = 0; n < Ni; n++) Interface[n].Copy(List.Interface[n]);
  for (int n = 0; n < Ni; n++) relations[n] = List.relations[n];
}

inline void Interface2D_List::Copy(const Interface2D_List &List) {
  allocate(List.Ni);
  for (int n = 0; n < Ni; n++) Interface[n].Copy(List.Interface[n]);
  for (int n = 0; n < Ni; n++) relations[n] = List.relations[n];
}

/**********************************************************************
 * Interface2D_List::max_speed -- Return the maximum interface speed. *
 **********************************************************************/
inline double Interface2D_List::max_speed(void) {
//   double vmax, v; vmax = ZERO;
//   for (int n = 0; n < Ni; n++) {
//     v = Interface[n].max_speed();
//     if (v > vmax) vmax = v;
//   }
//   return vmax;
  double vmax = ZERO;
  for (int n = 0; n < Ni; n++) vmax = max(ZERO,fabs(Interface[n].max_speed()));
  return vmax;
}

/**********************************************************************
 * Interface2D_List -- Output operator.                               *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const Interface2D_List &List) {
  out_file << List.Ni << endl;
  for (int n = 1; n <= List.Ni; n++) {
    out_file << "Ni=" << n << " of " << List.Ni << endl;
    out_file << List[n];
    List.relations[n-1].write(out_file);
  }
  return out_file;
}

/**********************************************************************
 * Interface2D_List -- Input operator.                                *
 **********************************************************************/
inline istream &operator >> (istream &in_file, Interface2D_List &List) {
  int num;
  if (List.Interface != NULL) List.deallocate();
  in_file.setf(ios::skipws);
  in_file >> num;
  in_file.unsetf(ios::skipws);
  if (num) {
    List.allocate(num);
    for (int n = 1; n <= List.Ni; n++) {
      in_file >> List[n];
      List.relations[n-1].read(in_file);
    }
  }
  return in_file;
}

/**********************************************************************
 * Interface2D -- External subroutines.                               *
 **********************************************************************/

extern void Broadcast_Interface(Interface2D &I);

#ifdef _MPI_VERSION
extern void Broadcast_Interface(Interface2D &I,
				MPI::Intracomm &Communicator, 
				const int Source_CPU);
#endif

extern double Intersection_Area(Interface2D &I1, Interface2D &I2);

/**********************************************************************
 * Interface2D_List -- External subroutines.                          *
 **********************************************************************/

extern void Broadcast_Interface_List(Interface2D_List &List);

#ifdef _MPI_VERSION
extern void Broadcast_Interface_List(Interface2D_List &List,
				     MPI::Intracomm &Communicator, 
				     const int Source_CPU);
#endif

#endif // _INTERFACE2D_INCLUDED

