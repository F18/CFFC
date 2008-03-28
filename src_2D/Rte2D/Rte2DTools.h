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


/********************************************************
 * STRUCTS REQUIRED FOR INTEGRATION                     *
 ********************************************************/

// Struct needed to integrate the exact solution for radiative
// heat transfer in a cylindrical enclosure.
struct exact_cyl_param {
  double z;       // the non-dimensional axial distance
  double r;       // the non-dimensional radial distance 
  double c;       // the half length divided by the outer radius
  double kappa;   // non-dimensional absorbsion coefficient
  int term_flag;  // a flag for which parameter we are computing
  int coord_flag; // a flag for which coordinate this is for
};

// Struct needed to integrate the exact solution for radiative
// heat transfer in a rectangular enclosure.
struct exact_rect_param {
  double x;       // the dimensional x location
  double y;       // the dimensional y location 
  double a1;      // west wall location
  double a2;      // east wall location
  double b1;      // south wall location
  double b2;      // north wall location
  double kappa;   // absorbsion coefficient
  int term_flag;  // a flag for which parameter we are computing
  int coord_flag; // a flag for which coordinate this is for
};


/********************************************************
 * Spacial Differencing Schemes                         *
 ********************************************************/
// Balsara's Multidimensional Skew-Upwind High Resolution scheme
// See: D. Balsara, "Fast and accurate discrete ordinates methods
// for multidimensional radiative heat transfer," in Journal       
// of Quantitative Spectroscopy, v69, 2001.    
// Note: this scheme uses the deffered correction method
void GM_Scheme( const double &dx,     // x-dir avg step size
		const double &dy,     // y-dir avg step size
		const double &vx,     // x-dir convective speed
		const double &vy,     // y-dir convective speed
		const double &S_id_j, // x-dir downstream Source
		const double &S_i_jd, // y-dir downstream Source
		const double &I_id_j, // x-dir downstream intensity
		const double &I_i_jd, // y-dir downstream intensity
		const double &I_iu_j, // x-dir upstream intensity
		const double &I_i_ju, // y-dir upstream intensity
		const double &I_i_j,  // cell intensity
		double &Ix_out,       // outgoing x-dir face intensity
		double &Iy_out );     // outgoing y-dir face intensity


/********************************************************
 * Exact solution functions                             *
 *                                                      *
 * These functions are used to compute exact solutions  *
 * for simple cases and output the results.             *
 ********************************************************/
double func_exact_cyl(int ndim, double *x, void *params);
double func_exact_rect(int ndim, double *x, void *params);

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
