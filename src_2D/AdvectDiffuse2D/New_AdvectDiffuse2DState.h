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
 *
 * Member functions
 *    Fd        -- Return diffusive flux.
 *     F        -- Return advective flux.
 *     s        -- Return regular source term.
 *     s_axi    -- Return axisymmetric source term.
 *    F_diff    -- Evaluates diffusive flux.
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
 * U == U;
 * U != U;
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
  double  u;   // Solution.
  static AdvectionVelocityType  V;      // Advection velocity as function of 'x' and 'y'.
  static DiffusionFieldType k;		// Non-linear diffusion coefficient as function of 'x', 'y' and 'u'
  //@}

  //! @name Creation and copy constructors.
  //@{
  //! Default constructor.
  AdvectDiffuse2D_State_New(void): u(0.0){ }

  //! Value Constructor
  AdvectDiffuse2D_State_New(const double &Val): u(Val){ }

  //! Copy constructor.
  AdvectDiffuse2D_State_New(const AdvectDiffuse2D_State_New &U){ }

  //! Assignment constructor.
  AdvectDiffuse2D_State_New(const double &uu,
			const Vector2D &VV,
			const double &kappa,
			const double &tau,
			const Vector2D &FF) {

  }

  //! Assignment constructor.
  AdvectDiffuse2D_State_New(const double &uu,
			const double &Vx,
			const double &Vy,
			const double &kappa,
			const double &tau) {

  }
    
  //! Assignment constructor.
  AdvectDiffuse2D_State_New(const double &uu,
			const double &Vx,
			const double &Vy,
			const double &kappa,
			const double &tau,
			const double &Fdx,
			const double &Fdy) {

  }

  /* Destructor. */
  // ~AdvectDiffuse2D_State_New(void);
  // Use automatically generated destructor.
  //@}

  //! @name Useful operators.
  //@{
  //! Vacuum/zero operator.
  //@}

  //! @name Advective Flux.
  //@{
  //@}

  //! @name Regular source term.
  //@{
  //@}

  //! @name Axisymmetric source term.
  //@{
  //@}

  //! @name Evaluates diffusive flux.
  //@{
  //@}

  /* @name Assignment operator. */
  //@{
  // AdvectDiffuse2D_State_New operator = (const AdvectDiffuse2D_State_New &W);
  // Use automatically generated assignment operator.

  //! @name Binary arithmetic operators.
  //@{
  //@}

  //! @name Unary arithmetic operators.
  //@{
  //@}

  //! @name Shortcut arithmetic operators.
  //@{
  //@}

  //! @name Relational operators.
  //@{
  //@}

  //! @name Input-output operators.
  //@{
  //@}

};



#endif /* _DVECTDIFFUSE2D_STATE_INCLUDED  */
