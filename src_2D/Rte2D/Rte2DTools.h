/******************* Rte2DTools.h ************************************

  Defines various functions for Rte2D that do not really belong 
  anywhere.

***********************************************************************/
#ifndef _RTE2D_TOOLS_INCLUDED
#define _RTE2D_TOOLS_INCLUDED

/*******************************************************
 * CFFC includes                                       *
 *******************************************************/  
#include "Rte2DInput.h"
#include "Rte2DState.h"


/*******************************************************
 * Main function to setup the Rte2D_State static       *
 * parameters                                          *
 *******************************************************/  
void SetupStateStatic( const Rte2D_Input_Parameters &IP );

/*******************************************************
 * Depending on the input parameters, initilize the    *
 * non-solution parameters of state U accordingly.     *
 *******************************************************/  
void SetInitialValues( Rte2D_State &U, 
		       const Rte2D_Input_Parameters &IP );

/********************************************************
 * Exact solution functions                             *
 *                                                      *
 * These functions are used to compute exact solutions  *
 * for simple cases and output the results.             *
 ********************************************************/
void CylindricalEnclosure( const double gas_temp,
			   const double c,
			   const double tau,
			   const double rpos,
			   const double zpos,
			   double &G, 
			   double &qr, 
			   double &qz );

void RectangularEnclosure( const double gas_temp,
			   const double kappa,
			   const double left,
			   const double right,
			   const double bot,
			   const double top,
			   const double xpos,
			   const double ypos,
			   double &G, 
			   double &qx, 
			   double &qy );


#endif // _RTE2D_TOOLS_INCLUDED
