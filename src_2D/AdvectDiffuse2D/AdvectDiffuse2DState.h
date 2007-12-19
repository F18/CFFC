/*!\file New_AdvectDiffuse2DState.h
  \brief Temporary header file defining 2D Advection Diffusion Equation Solution State Class. */

#ifndef _NEW_ADVECTDIFFUSE2D_STATE_INCLUDED
#define _NEW_ADVECTDIFFUSE2D_STATE_INCLUDED

/* Include required C++ libraries. */
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

/* Using std namespace functions */
using namespace std;

/* Include CFFC header files */
#include "../Math/Math.h" /* Include math macro header files. */
#include "../CFD/CFD.h"   /* Include CFD header files. */
#include "../Math/Vector2D.h" /* Include vector header files. */
#include "AdvectDiffuse2DExactSolutions.h" /* Include 2D advection diffusion exact solutions header file */
#include "AdvectDiffuse2DParameterFields.h" /* Include 2D advection diffusion parameter fields */

/* Define the classes. */

#define	NUM_VAR_ADVECTDIFFUSE2D    1

/*!
 * \class AdvectDiffuse2D_State_New
 *
 * @brief Solution state class definition for the 2D advection-diffusion
 *        equation.
 *
 * Solution state class definition for the 2D advection-diffusion
 * equation.
 *
 * \verbatim
 * Member variables
 *     u        -- Return solution.
 *     V        -- Return advection velocity field
 *     k        -- Return diffusion coefficient field
 * SourceTerm   -- Return source term field
 *
 * Member functions
 *    Fa        -- Return advective flux.
 *    Fd        -- Return diffusive flux.
 *
 * Member operators
 *      U -- a solution state
 *      c -- a scalar (double)
 *
 * U = U;
 * U = U + U;
 * U = U - U;
 * U = U * U;
 * U = c * U;
 * U = U * c;
 * U = U / c;
 * U = +U;
 * U = -U;
 * U += U;
 * U -= U;
 * U *= a;
 * U /= a;
 * U = U ^ U; (a useful product)
 * U == U;
 * U != U;
 * U <= U;
 * U >= U;
 * cout << U; (output function)
 * cin  >> U; (input function)
 * \endverbatim
 */
class AdvectDiffuse2D_State_New{
public:
  //! @name Defined datatypes
  //@{
  typedef VelocityFields::VelocityFieldType AdvectionVelocityType;
  typedef DiffusionFields::NonlinearDiffusionFieldType DiffusionFieldType;
  //@}

  //! @name Solution state variables and associated constants:
  //@{
  double  u;                            //!< Solution.
  static AdvectionVelocityType  V;      //!< Advection velocity as function of 'x' and 'y'.
  static DiffusionFieldType k;		//!< Non-linear diffusion coefficient as function of 'x', 'y' and 'u'
  static SourceTermFields *SourceTerm;  //!< Non-linear source term as function of 'x', 'y' and 'u' or only '\f$\bar{u}\f$'
  //@}

  //! @name Creation and copy constructors.
  //@{
  //! Default constructor.
  AdvectDiffuse2D_State_New(void);

  //! Value Constructor
  AdvectDiffuse2D_State_New(const double &uu);

  //! Copy constructor.
  AdvectDiffuse2D_State_New(const AdvectDiffuse2D_State_New &U){ u = U.u; }
  //@}

  //! Destructor.
  ~AdvectDiffuse2D_State_New(void){}

  //! @name Useful operators.
  //@{
  //! Return the number of variables.
  int NumVar(void) { return NUM_VAR_ADVECTDIFFUSE2D; }

  //! Copy operator.
  void Copy(const AdvectDiffuse2D_State_New &U) { u = U.u; }

  //! Vacuum/zero operator.
  void Vacuum(void) { u = ZERO; }

  //! Set the pointers to the fields (i.e. velocity, diffusion and source)
  static void Set_Pointers_To_Fields(void);
  //@}

  //! @name Functions required for multigrid.
  //@{
  //! Copy variables solved by multigrid only.
  void Copy_Multigrid_State_Variables(const AdvectDiffuse2D_State_New &Ufine) { Copy(Ufine); }

  //! Zero variables not-solved by multigrid.
  void Zero_Non_Multigrid_State_Variables(void){};

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const{};
  //@}

  //! @name Compute upwind state
  //@{
  friend AdvectDiffuse2D_State_New RiemannState(const AdvectDiffuse2D_State_New &Ul,
						const AdvectDiffuse2D_State_New &Ur,
						const double &x, const double &y,
						const Vector2D & norm_dir);
  //@}

  //! @name Advective Flux.
  //@{
  Vector2D Fa(const double &x, const double &y) const;
  friend Vector2D Fa(const AdvectDiffuse2D_State_New &U, const double &x, const double &y);
  AdvectDiffuse2D_State_New Fa(const double &x, const double &y, const Vector2D & norm_dir) const;
  friend AdvectDiffuse2D_State_New Fa(const AdvectDiffuse2D_State_New &Ul,
				      const AdvectDiffuse2D_State_New &Ur,
				      const double &x, const double &y,
				      const Vector2D & norm_dir);
  friend AdvectDiffuse2D_State_New Fa(const AdvectDiffuse2D_State_New &Ul,
				      const AdvectDiffuse2D_State_New &Ur,
				      const Vector2D & ComputationPoint,
				      const Vector2D & norm_dir);
  //@}

  //! @name Computation of diffusive flux.
  //@{
  //! Calculate flux based on the position (x,y) and the gradient of the scalar
  Vector2D Fd(const double &dudx, const double &dudy,
	      const double &x, const double &y) const;
  friend Vector2D Fd(const AdvectDiffuse2D_State_New &U,
		     const double &dudx, const double &dudy,
		     const double &x, const double &y);
  Vector2D Fd(const Vector2D &GradU, 
	      const double &x, const double &y) const;
  friend Vector2D Fd(const AdvectDiffuse2D_State_New &U, 
		     const Vector2D &GradU,
		     const double &x, const double &y);

  //! Calculate flux in the normal direction norm_dir
  AdvectDiffuse2D_State_New Fd(const double &dudx, const double &dudy,
			       const double &x, const double &y,
			       const Vector2D & norm_dir) const{
    return AdvectDiffuse2D_State_New(Fd(dudx,dudy,x,y)*norm_dir);
  }

  friend AdvectDiffuse2D_State_New Fd(const AdvectDiffuse2D_State_New &U,
				      const double &dudx, const double &dudy,
				      const double &x, const double &y,
				      const Vector2D & norm_dir){
    return U.Fd(dudx,dudy,x,y,norm_dir);
  }

  //! Calculate the interface flux in the normal direction norm_dir and use the mean average to estimate k(x,y,u)
  friend AdvectDiffuse2D_State_New Fd(const AdvectDiffuse2D_State_New &Ul,
				      const AdvectDiffuse2D_State_New &Ur,
				      const double &dudx, const double &dudy,
				      const double &x, const double &y,
				      const Vector2D & norm_dir) {
    return (0.5*(Ul+Ur)).Fd(dudx,dudy,x,y,norm_dir);
  }

  //! Calculate flux using the projection of the gradient on the normal direction \f$GradU_n = \nabla u \cdot \vec{n} \f$
  AdvectDiffuse2D_State_New Fd(const double &GradU_n,
			       const double &x, const double &y) const;

  //! Calculate flux using a provided state, the normal gradient projection and the position of interest.
  friend AdvectDiffuse2D_State_New Fd(const AdvectDiffuse2D_State_New &U,
				      const double &GradU_n,
				      const double &x, const double &y){
    return U.Fd(GradU_n,x,y);
  }

  //! Calculate the interface flux at the position of interest in the normal direction,
  //  based on the interface states and the scalar gradient
  friend AdvectDiffuse2D_State_New Fd(const AdvectDiffuse2D_State_New &Ul,
				      const AdvectDiffuse2D_State_New &Ur,
				      const Vector2D & GradU,
				      const double &x, const double &y,
				      const Vector2D & norm_dir);

  friend AdvectDiffuse2D_State_New Fd(const AdvectDiffuse2D_State_New &U,
				      const Vector2D & GradU,
				      const Vector2D & ComputationPoint,
				      const Vector2D & norm_dir);

  friend AdvectDiffuse2D_State_New Fd(const AdvectDiffuse2D_State_New &Ul,
				      const AdvectDiffuse2D_State_New &Ur,
				      const Vector2D & GradU,
				      const Vector2D & ComputationPoint,
				      const Vector2D & norm_dir);
  //@}

  //! @name Regular source term.
  //@{
  //! \brief Evaluate integral source field
  double s(void) const;
  friend double s(const AdvectDiffuse2D_State_New &U);
  //! \brief Evaluate pointwise source field
  double s(const double &x, const double &y, const double &u) const;
  //! \brief Evaluate the source field at a given location
  double source(const double &x, const double &y, const AdvectDiffuse2D_State_New &U);
  //@}

  //! @name Axisymmetric source term.
  //@{
  double s_axi(const Vector2D &X) const;
  friend double s_axi(const AdvectDiffuse2D_State_New &U, const Vector2D &X);
  //@}

  /*! @name Assignment operator. */
  //@{
  AdvectDiffuse2D_State_New& operator= (const AdvectDiffuse2D_State_New &rhs);
  //@}

  /*! @name Index operators */
  //@{
  double &operator[](const int & index);
  const double &operator[](const int & index) const;
  //@}

  /*! Absolute value */
  friend AdvectDiffuse2D_State_New fabs(const AdvectDiffuse2D_State_New &U){ return AdvectDiffuse2D_State_New(fabs(U.u));}

  //! @name Binary arithmetic operators.
  //@{
  AdvectDiffuse2D_State_New operator +(const AdvectDiffuse2D_State_New &U) const;
  AdvectDiffuse2D_State_New operator -(const AdvectDiffuse2D_State_New &U) const;
  AdvectDiffuse2D_State_New operator *(const AdvectDiffuse2D_State_New &U) const;
  AdvectDiffuse2D_State_New operator *(const double &a) const;
  friend AdvectDiffuse2D_State_New operator *(const double &a, const AdvectDiffuse2D_State_New &U){ return U*a;}
  AdvectDiffuse2D_State_New operator /(const double &a){ return *this * (1.0/a); }
  AdvectDiffuse2D_State_New operator ^(const AdvectDiffuse2D_State_New &U) const;
  //@}

  //! @name Unary arithmetic operators.
  //@{
  friend AdvectDiffuse2D_State_New operator +(const AdvectDiffuse2D_State_New &U);
  friend AdvectDiffuse2D_State_New operator -(const AdvectDiffuse2D_State_New &U);
  //@}

  //! @name Shortcut arithmetic operators.
  //@{
  AdvectDiffuse2D_State_New & operator +=(const AdvectDiffuse2D_State_New &U);
  AdvectDiffuse2D_State_New & operator -=(const AdvectDiffuse2D_State_New &U);
  AdvectDiffuse2D_State_New & operator *=(const AdvectDiffuse2D_State_New &U);
  AdvectDiffuse2D_State_New & operator *=(const double &a);
  AdvectDiffuse2D_State_New & operator /=(const double &a);
  //@}

  //! @name Relational operators.
  //@{
  friend bool operator ==(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs);
  friend bool operator !=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs);
  friend bool operator >=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs);
  friend bool operator <=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs);
  //@}

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file, const AdvectDiffuse2D_State_New &U);
  friend istream &operator >> (istream &in_file,  AdvectDiffuse2D_State_New &U);
  //@}
};

// Default constructor.
inline AdvectDiffuse2D_State_New::AdvectDiffuse2D_State_New(void): u(0.0){
}

// Value Constructor
inline AdvectDiffuse2D_State_New::AdvectDiffuse2D_State_New(const double &uu): u(uu){
}


inline void AdvectDiffuse2D_State_New::Set_Pointers_To_Fields(void){
  // Set pointer to the velocity field
  VelocityFields::Connect_Pointer_To_Flow_Field(V);

  // Set pointer to the diffusion field
  DiffusionFields::Connect_Pointer_To_Diffusion_Field(k);

  // Get access to the SourceTermField
  SourceTerm = &SourceTermFields::getInstance();
}

/*!
 * Compute the Riemann state
 * 
 * \param [in] x x-coordinate location
 * \param [in] y y-coordinate location
 * \return the solution to Riemann problem in the 'norm_dir'-direction (defined from left to right)
 */
inline AdvectDiffuse2D_State_New RiemannState(const AdvectDiffuse2D_State_New &Ul,
					      const AdvectDiffuse2D_State_New &Ur,
					      const double &x, const double &y,
					      const Vector2D & norm_dir){
  return (AdvectDiffuse2D_State_New::V(x,y)*norm_dir > 0) ? Ul : Ur;
}

/*!
 * Compute the advective flux \f$Fa=\vec{V}(x,y) \, u\f$
 * 
 * \param [in] x x-coordinate location
 * \param [in] y y-coordinate location
 * \return advective flux at a given location
 */
inline Vector2D AdvectDiffuse2D_State_New::Fa(const double &x, const double &y) const {
  return V(x,y)*u;
}

/*!
 * Compute the advective flux \f$Fa=\vec{V}(x,y) \, u\f$
 *
 * \param [in] x x-coordinate location
 * \param [in] y y-coordinate location
 * \param [in] U solution state used to compute the advective flux
 * \return advective flux at a given location based on the solution state U
 */
inline Vector2D Fa(const AdvectDiffuse2D_State_New &U,
		   const double &x, const double &y) {
  return U.Fa(x,y);
}

/*!
 * Compute the advective flux \f$Fa=(\vec{V}(x,y) \cdot \vec{n} )\, u\f$
 *
 * \param [in] x        x-coordinate location
 * \param [in] y        y-coordinate location
 * \param [in] norm_dir the direction used to project the flux onto
 * \return the projection of the advective flux onto 
 * the 'norm-dir' direction at a given location based on the solution state U
 */
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::Fa(const double &x, const double &y,
							       const Vector2D & norm_dir) const {
  return AdvectDiffuse2D_State_New((V(x,y)*norm_dir)*u);
}

/*!
 * Compute the upwind advective flux \f$Fa=(\vec{V}(x,y) \cdot \vec{n})\, u\f$ 
 * at an interface between 2 states.
 *
 * \param [in] Ul       the state to the left of the interface
 * \param [in] Ur       the state to the right of the interface
 * \param [in] x        x-coordinate location
 * \param [in] y        y-coordinate location
 * \param [in] norm_dir the direction used to project the flux onto
 * \return the projection of the upwind advective flux onto the 'norm-dir' direction 
 * at a given location based on the solution state chosen by the upwinding process
 */
inline AdvectDiffuse2D_State_New Fa(const AdvectDiffuse2D_State_New &Ul,
				    const AdvectDiffuse2D_State_New &Ur,
				    const double &x, const double &y,
				    const Vector2D & norm_dir) {
  return (Ul.V(x,y)*norm_dir > 0) ? Ul.Fa(x,y,norm_dir) : Ur.Fa(x,y,norm_dir);
}

/*!
 * Compute the upwind advective flux \f$Fa=(\vec{V}(x,y) \cdot \vec{n})\, u\f$ 
 * at an interface between 2 states.
 *
 * \param [in] Ul       the state to the left of the interface
 * \param [in] Ur       the state to the right of the interface
 * \param [in] ComputationPoint the Cartesian coordinates of the location where the flux is computed
 * \param [in] norm_dir         the direction used to project the flux onto
 * \return the projection of the upwind advective flux onto the 'norm-dir' direction 
 * at a given location based on the solution state chosen by the upwinding process
 */
inline AdvectDiffuse2D_State_New Fa(const AdvectDiffuse2D_State_New &Ul,
				    const AdvectDiffuse2D_State_New &Ur,
				    const Vector2D & ComputationPoint,
				    const Vector2D & norm_dir){
  return Fa(Ul,Ur,ComputationPoint.x,ComputationPoint.y,norm_dir);
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y,u) \, \nabla u\f$, where 'u' is the solution scalar.
 * 
 * \param [in] x    x-coordinate location
 * \param [in] y    y-coordinate location
 * \param [in] dudx solution gradient with respect to x-direction
 * \param [in] dudy solution gradient with respect to y-direction
 * \return diffusive flux at a given location
 */
inline Vector2D AdvectDiffuse2D_State_New::Fd(const double &dudx, const double &dudy,
					      const double &x, const double &y) const {
  return Vector2D(-k(x,y,u)*dudx ,-k(x,y,u)*dudy);
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y,u) \, \nabla u\f$, where 'u' is the solution scalar.
 * 
 * \param [in] x     x-coordinate location
 * \param [in] y     y-coordinate location
 * \param [in] GradU solution gradient vector
 * \return diffusive flux at a given location
 */
inline Vector2D AdvectDiffuse2D_State_New::Fd(const Vector2D &GradU,
					      const double &x, const double &y) const {
  return Fd(GradU.x,GradU.y,x,y);
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y,u) \, \nabla u\f$.
 * 
 * \param [in] x    x-coordinate location
 * \param [in] y    y-coordinate location
 * \param [in] dudx solution gradient with respect to x-direction
 * \param [in] dudy solution gradient with respect to y-direction
 * \param [in] U    solution state used to compute the diffusive flux
 * \return diffusive flux at a given location
 */
inline Vector2D Fd(const AdvectDiffuse2D_State_New &U,
		   const double &dudx, const double &dudy,
		   const double &x, const double &y) {
  return U.Fd(dudx,dudy,x,y);
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y,u) \, \nabla u\f$, where 'u' is the solution scalar.
 * 
 * \param [in] x     x-coordinate location
 * \param [in] y     y-coordinate location
 * \param [in] GradU solution gradient vector
 * \param [in] U     solution state used to compute the diffusive flux
 * \return diffusive flux at a given location
 */
inline Vector2D Fd(const AdvectDiffuse2D_State_New &U,
		   const Vector2D &GradU,
		   const double &x, const double &y) {
  return U.Fd(GradU.x,GradU.y,x,y);
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y,u) ( \nabla u \cdot \vec{n})\f$, where 'u' is the solution scalar.
 * 
 * \param [in] x        x-coordinate location
 * \param [in] y        y-coordinate location
 * \param [in] GradU_n  the projection of the solution gradient onto the normal direction
 *
 * \return diffusive flux at a given location which has the diffusion coefficient estimated based on the state scalar
 */
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::Fd(const double &GradU_n,
							       const double &x, const double &y) const {
  return AdvectDiffuse2D_State_New(-k(x,y,u)*GradU_n);
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y,u) ( \nabla u \cdot \vec{n})\f$ at the interface between 2 states. \n
 * 
 * The diffusion coefficient is estimated based on the average of the two interface-solution states. 
 *
 * \param [in] Ul         the state to the left of the interface
 * \param [in] Ur         the state to the right of the interface
 * \param [in] GradU      solution gradient vector
 * \param [in] x          x-coordinate location
 * \param [in] y          y-coordinate location
 * \param [in] norm_dir   the direction used to project the flux onto
 *
 * \return the projection of the diffusive flux onto the 'norm-dir' direction  at a given location 
 *         based on the solution gradient and the mean average state at the interface.
 */
inline AdvectDiffuse2D_State_New Fd(const AdvectDiffuse2D_State_New &Ul,
				    const AdvectDiffuse2D_State_New &Ur,
				    const Vector2D & GradU,
				    const double &x, const double &y,
				    const Vector2D & norm_dir){
  return Fd(0.5*(Ul+Ur), GradU*norm_dir, x, y);  
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y,u) ( \nabla u \cdot \vec{n})\f$ at the interface between 2 states. \n
 * 
 * \param [in] U          the state used for computing the diffusion coefficient
 * \param [in] GradU      solution gradient vector
 * \param [in] ComputationPoint the Cartesian coordinates of the location where the flux is computed
 * \param [in] norm_dir   the direction used to project the flux onto
 *
 * \return the projection of the diffusive flux onto the 'norm-dir' direction  at a given location 
 *         based on the solution gradient and the state U at the interface.
 */
inline AdvectDiffuse2D_State_New Fd(const AdvectDiffuse2D_State_New &U,
				    const Vector2D & GradU,
				    const Vector2D & ComputationPoint,
				    const Vector2D & norm_dir){
  return Fd(U,GradU*norm_dir,ComputationPoint.x,ComputationPoint.y);  
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y,u) ( \nabla u \cdot \vec{n})\f$ at the interface between 2 states. \n
 * 
 * \param [in] Ul         the state to the left of the interface
 * \param [in] Ur         the state to the right of the interface
 * \param [in] GradU      solution gradient vector
 * \param [in] ComputationPoint the Cartesian coordinates of the location where the flux is computed
 * \param [in] norm_dir   the direction used to project the flux onto
 *
 * \return the projection of the diffusive flux onto the 'norm-dir' direction  at a given location 
 *         based on the solution gradient and the mean average state at the interface.
 */
inline AdvectDiffuse2D_State_New Fd(const AdvectDiffuse2D_State_New &Ul,
				    const AdvectDiffuse2D_State_New &Ur,
				    const Vector2D & GradU,
				    const Vector2D & ComputationPoint,
				    const Vector2D & norm_dir){
  return Fd(Ul,Ur,GradU,ComputationPoint.x,ComputationPoint.y,norm_dir);  
}

/*!
 * Integral source term (i.e function of the average solution)
 */
inline double AdvectDiffuse2D_State_New::s(void) const{
  return SourceTerm->operator()(u);
}

/*!
 * Integral source term (i.e function of the average solution)
 */
inline double s(const AdvectDiffuse2D_State_New &U){
  return U.s();
}

/*!
 * Pointwise defined source term (i.e function of the Cartesian coordinates and the solution at the point)
 */
inline double AdvectDiffuse2D_State_New::s(const double &x, const double &y, const double &u) const{
  return SourceTerm->operator()(x,y,u);
}

/*!
 * Return the value of the source field at a given location.
 * This subroutine checks if the source field is pointwise or integral.
 */
inline double source(const double &x, const double &y, const AdvectDiffuse2D_State_New &U) {
  if (U.SourceTerm->FieldRequireIntegration()){
    return U.SourceTerm->operator()(x,y,U.u); // compute the pointwise value
  } else {
    return U.SourceTerm->operator()(U.u);   // compute the integral value
  }
}

/*!
 * Axisymmetric source term.
 */
inline double AdvectDiffuse2D_State_New::s_axi(const Vector2D &X) const {
    return (ZERO);
}

/*!
 * Axisymmetric source term.
 */
inline double s_axi(const AdvectDiffuse2D_State_New &U, const Vector2D &X) {
  return (ZERO);
}

/*!
 * Assignment operator
 */
inline AdvectDiffuse2D_State_New& AdvectDiffuse2D_State_New::operator= (const AdvectDiffuse2D_State_New &rhs){
  if(this == &rhs) return *this;
  u = rhs.u;
  return *this;
}

//! index operator 
inline double & AdvectDiffuse2D_State_New::operator[](const int & index){
  assert( index >= 1 && index <= NUM_VAR_ADVECTDIFFUSE2D );
  switch(index) {
  case 1 :
    return (u);
  default:
    return (u);
  };
}

//! constant index operator 
inline const double & AdvectDiffuse2D_State_New::operator[](const int & index) const {
  assert( index >= 1 && index <= NUM_VAR_ADVECTDIFFUSE2D );
  switch(index) {
  case 1 :
    return (u);
  default:
    return (u);
  };
}

//! summation between states
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::operator +(const AdvectDiffuse2D_State_New &U) const {
  return AdvectDiffuse2D_State_New(u + U.u);
}

//! subtraction between states
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::operator -(const AdvectDiffuse2D_State_New &U) const {
  return AdvectDiffuse2D_State_New(u - U.u);
}

//! multiplication between states
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::operator *(const AdvectDiffuse2D_State_New &U) const {
  return AdvectDiffuse2D_State_New(u * U.u);
}

//! multiplication with scalar
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::operator *(const double &a) const{
  return AdvectDiffuse2D_State_New(a * u);
}

//! A useful solution state product operator.
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::operator ^(const AdvectDiffuse2D_State_New &U) const{
  return AdvectDiffuse2D_State_New(u * U.u);
}

//! unary arithmetic + operator
inline AdvectDiffuse2D_State_New operator +(const AdvectDiffuse2D_State_New &U){
  return U;
}

//! unary arithmetic - operator
inline AdvectDiffuse2D_State_New operator -(const AdvectDiffuse2D_State_New &U){
  return AdvectDiffuse2D_State_New(-U.u);
}

//! shortcut summation
inline AdvectDiffuse2D_State_New & AdvectDiffuse2D_State_New::operator +=(const AdvectDiffuse2D_State_New &U){
  u += U.u;
  return *this;
}

//! shortcut subtraction
inline AdvectDiffuse2D_State_New & AdvectDiffuse2D_State_New::operator -=(const AdvectDiffuse2D_State_New &U){
  u -= U.u;
  return *this;
}

//! shortcut multiplication
inline AdvectDiffuse2D_State_New & AdvectDiffuse2D_State_New::operator *=(const AdvectDiffuse2D_State_New &U){
  u *= U.u;
  return *this;
}

//! shortcut multiplication with scalar
inline AdvectDiffuse2D_State_New & AdvectDiffuse2D_State_New::operator *=(const double &a){
  u *= a;
  return *this;
}

//! shortcut division with scalar
inline AdvectDiffuse2D_State_New & AdvectDiffuse2D_State_New::operator /=(const double &a){
  u /= a;
  return *this;
}

//! equal operator 
inline bool operator ==(const AdvectDiffuse2D_State_New &lhs, const AdvectDiffuse2D_State_New &rhs) {
  return (lhs.u == rhs.u);
}

//! not-equal operator
inline bool operator !=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs){
  return !(lhs == rhs);
}

//! greater or equal operator
inline bool operator >=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs){
  return (lhs.u >= rhs.u);
}

//! less or equal operator 
inline bool operator <=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs){
  return (lhs.u <= rhs.u);
}

//! output operator
inline ostream &operator << (ostream &out_file, const AdvectDiffuse2D_State_New &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.u ;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, AdvectDiffuse2D_State_New &U) {
  in_file.setf(ios::skipws);
  in_file >> U.u ;
  in_file.unsetf(ios::skipws);
  return (in_file);
}


#endif /* _DVECTDIFFUSE2D_STATE_INCLUDED  */
