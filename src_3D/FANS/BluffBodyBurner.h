/**********************************************************************
 Flowfield data for bluff-body stabilised turbulent non-reacting flows
 Burner bluff-body diameter, D_b=50mm, (R_b=25mm), jet diameter=3.6mm
 **********************************************************************/

#ifndef _BLUFFBODYBURNER_INCLUDED
#define _BLUFFBODYBURNER_INCLUDED

/* Include required C++ libraries. */

#include <cassert>
#include <cstdlib>

using namespace std;

/* Include math macro file. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

/* Include Vector2D header. */

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

/* Include Spline2D header. */

#ifndef _SPLINE2D_INCLUDED
#include "../Math/Spline2D.h"
#endif //_SPLINE2D_INCLUDED

/*!
 * Class: NonreactiveBluffBodyBurner
 *
 * @brief Normalized initial conditions for the bluff-body stabilised
 * turbulent non-reacting flows.
 *
 *
 * \verbatim
 * Private member functions
 *      Ns -- Number of stations.
 *      U  -- Array of u-velocity splines at each station.
 *      V  -- Array of v-velocity splines at each station.
 *      x  -- Array of station x-coordinates (axial direction).
 * Public member functions
 *      allocate      -- Allocate memory and initialize arrays.
 *      deallocate    -- Deallocate memory for splines and arrays.
 *      interpolation -- Interpolation function.
 * \endverbatim
 */
class NonreactiveVelocityField {
/* Velocity of coflowing air is 20 m/s. */
/* Velocity of coflowing fuel (in this case "air" ) is 61 m/s . */
 private:
  int       Ns; //!< Number of stations.
  Spline2D  *U; //!< Array of u-velocity splines at each station.
  Spline2D  *V; //!< Array of v-velocity splines at each station.
  double    *x; //!< Array of station x-coordinates (axial direction).
 
 public:
  //! Creation constructor.
  NonreactiveVelocityField(void) { allocate(); }

  //! Destructor.
  ~NonreactiveVelocityField(void) { deallocate(); }

  //! Allocate memory and initialize arrays.
  void allocate(void);

  //! Deallocate memory for splines and arrays.
  void deallocate(void);

  //! Interpolation function.
  Vector2D interpolation(const Vector2D &Xt);

};

//! NonreactiveVectorField::allocate -- Allocate memory and
//!                                          initialize arrays.
inline void NonreactiveVelocityField::allocate(void) {
/* Velocity of coflowing air is 20 m/s. */
/* Velocity of coflowing fuel (in this case "air" ) is 61 m/s . */
// Note the original data in the website is the commented section
// however, the y coordinat is different in this computation frame from it in experimental
// frame , so "-" goes to "+" and "+" goes to "-" for vertical velocity
  // Set the number of stations.
  Ns = 13;
  // Allocate the arrays.
  U = new Spline2D[Ns];
  V = new Spline2D[Ns];
  x = new double[Ns];
 // STATION 1
  x[0] = 0.06;
  U[0].allocate(24);                              V[0].allocate(24);
  U[0].settype(SPLINE2D_LINEAR);               V[0].settype(SPLINE2D_LINEAR);
  U[0].Xp[ 0].y = 0.00;  U[0].Xp[ 0].x = 75.180;  V[0].Xp[ 0].x = - 0.10;  V[0].Xp[ 0].y = 0.00;
  U[0].Xp[ 1].y = 0.02;  U[0].Xp[ 1].x = 74.400;  V[0].Xp[ 1].x = - 0.01;  V[0].Xp[ 1].y = 0.02;
  U[0].Xp[ 2].y = 0.04;  U[0].Xp[ 2].x = 70.990;  V[0].Xp[ 2].x = 0.06;  V[0].Xp[ 2].y = 0.04;
  U[0].Xp[ 3].y = 0.06;  U[0].Xp[ 3].x = 53.030;  V[0].Xp[ 3].x = - 0.14;  V[0].Xp[ 3].y = 0.06;
  U[0].Xp[ 4].y = 0.08;  U[0].Xp[ 4].x = 27.300;  V[0].Xp[ 4].x = 2.24;  V[0].Xp[ 4].y = 0.08;
  U[0].Xp[ 5].y = 0.10;  U[0].Xp[ 5].x =  5.040;  V[0].Xp[ 5].x = 0.62;  V[0].Xp[ 5].y = 0.10;
  U[0].Xp[ 6].y = 0.12;  U[0].Xp[ 6].x = -1.440;  V[0].Xp[ 6].x = - 1.58;  V[0].Xp[ 6].y = 0.12;
  U[0].Xp[ 7].y = 0.16;  U[0].Xp[ 7].x = -1.210;  V[0].Xp[ 7].x = - 1.21;  V[0].Xp[ 7].y = 0.16;
  U[0].Xp[ 8].y = 0.20;  U[0].Xp[ 8].x = -1.200;  V[0].Xp[ 8].x = - 0.76;  V[0].Xp[ 8].y = 0.20;
  U[0].Xp[ 9].y = 0.24;  U[0].Xp[ 9].x = -1.030;  V[0].Xp[ 9].x = - 0.36;  V[0].Xp[ 9].y = 0.24;
  U[0].Xp[10].y = 0.32;  U[0].Xp[10].x = -0.820;  V[0].Xp[10].x = 0.61;  V[0].Xp[10].y = 0.32;
  U[0].Xp[11].y = 0.40;  U[0].Xp[11].x = -0.870;  V[0].Xp[11].x = 0.96;  V[0].Xp[11].y = 0.40;
  U[0].Xp[12].y = 0.48;  U[0].Xp[12].x = -0.680;  V[0].Xp[12].x = 1.57;  V[0].Xp[12].y = 0.48;
  U[0].Xp[13].y = 0.56;  U[0].Xp[13].x = -0.550;  V[0].Xp[13].x = 1.98;  V[0].Xp[13].y = 0.56;
  U[0].Xp[14].y = 0.64;  U[0].Xp[14].x = -0.380;  V[0].Xp[14].x = 2.27;  V[0].Xp[14].y = 0.64;
  U[0].Xp[15].y = 0.72;  U[0].Xp[15].x = -0.020;  V[0].Xp[15].x = 2.38;  V[0].Xp[15].y = 0.72;
  U[0].Xp[16].y = 0.80;  U[0].Xp[16].x =  0.200;  V[0].Xp[16].x = 2.22;  V[0].Xp[16].y = 0.80;
  U[0].Xp[17].y = 0.88;  U[0].Xp[17].x =  0.440;  V[0].Xp[17].x = 1.99;  V[0].Xp[17].y = 0.88;
  U[0].Xp[18].y = 0.92;  U[0].Xp[18].x =  0.480;  V[0].Xp[18].x = 1.74;  V[0].Xp[18].y = 0.92;
  U[0].Xp[19].y = 0.96;  U[0].Xp[19].x =  0.390;  V[0].Xp[19].x = 1.19;  V[0].Xp[19].y = 0.96;
  U[0].Xp[20].y = 1.00;  U[0].Xp[20].x =  4.740;  V[0].Xp[20].x = 0.50;  V[0].Xp[20].y = 1.00;
  U[0].Xp[21].y = 1.04;  U[0].Xp[21].x = 15.460;  V[0].Xp[21].x = - 0.79;  V[0].Xp[21].y = 1.04;
  U[0].Xp[22].y = 1.12;  U[0].Xp[22].x = 17.280;  V[0].Xp[22].x = - 0.78;  V[0].Xp[22].y = 1.12;
  U[0].Xp[23].y = 1.20;  U[0].Xp[23].x = 20.330;  V[0].Xp[23].x = - 0.74;  V[0].Xp[23].y = 1.20;
  U[0].pathlength();                              V[0].pathlength();
  
  // STATION 2
  x[1] = 0.20;
  U[1].allocate(24);                              V[1].allocate(24);
  U[1].settype(SPLINE2D_LINEAR);               V[1].settype(SPLINE2D_LINEAR);
  U[1].Xp[ 0].y = 0.00;  U[1].Xp[ 0].x = 74.830;  V[1].Xp[ 0].x = 0.06;  V[1].Xp[ 0].y = 0.00;
  U[1].Xp[ 1].y = 0.02;  U[1].Xp[ 1].x = 69.120;  V[1].Xp[ 1].x = - 0.53;  V[1].Xp[ 1].y = 0.02;
  U[1].Xp[ 2].y = 0.04;  U[1].Xp[ 2].x = 63.520;  V[1].Xp[ 2].x = - 0.21;  V[1].Xp[ 2].y = 0.04;
  U[1].Xp[ 3].y = 0.06;  U[1].Xp[ 3].x = 43.950;  V[1].Xp[ 3].x = 1.34;  V[1].Xp[ 3].y = 0.06;
  U[1].Xp[ 4].y = 0.08;  U[1].Xp[ 4].x = 33.000;  V[1].Xp[ 4].x = 1.81;  V[1].Xp[ 4].y = 0.08;
  U[1].Xp[ 5].y = 0.10;  U[1].Xp[ 5].x = 15.810;  V[1].Xp[ 5].x = 1.48;  V[1].Xp[ 5].y = 0.10;
  U[1].Xp[ 6].y = 0.12;  U[1].Xp[ 6].x =  9.660;  V[1].Xp[ 6].x = 1.27;  V[1].Xp[ 6].y = 0.12;
  U[1].Xp[ 7].y = 0.16;  U[1].Xp[ 7].x =  0.780;  V[1].Xp[ 7].x = - 0.22;  V[1].Xp[ 7].y = 0.16;
  U[1].Xp[ 8].y = 0.20;  U[1].Xp[ 8].x = -3.090;  V[1].Xp[ 8].x =  -1.01;  V[1].Xp[ 8].y = 0.20;
  U[1].Xp[ 9].y = 0.24;  U[1].Xp[ 9].x = -3.580;  V[1].Xp[ 9].x =  -0.62;  V[1].Xp[ 9].y = 0.24;
  U[1].Xp[10].y = 0.32;  U[1].Xp[10].x = -3.700;  V[1].Xp[10].x = 0.04;  V[1].Xp[10].y = 0.32;
  U[1].Xp[11].y = 0.40;  U[1].Xp[11].x = -3.690;  V[1].Xp[11].x = 0.50;  V[1].Xp[11].y = 0.40;
  U[1].Xp[12].y = 0.48;  U[1].Xp[12].x = -3.320;  V[1].Xp[12].x = 0.98;  V[1].Xp[12].y = 0.48;
  U[1].Xp[13].y = 0.56;  U[1].Xp[13].x = -2.740;  V[1].Xp[13].x = 1.29;  V[1].Xp[13].y = 0.56;
  U[1].Xp[14].y = 0.64;  U[1].Xp[14].x = -1.850;  V[1].Xp[14].x = 1.44;  V[1].Xp[14].y = 0.64;
  U[1].Xp[15].y = 0.72;  U[1].Xp[15].x = -0.990;  V[1].Xp[15].x = 1.43;  V[1].Xp[15].y = 0.72;
  U[1].Xp[16].y = 0.80;  U[1].Xp[16].x =  0.070;  V[1].Xp[16].x = 1.26;  V[1].Xp[16].y = 0.80;
  U[1].Xp[17].y = 0.88;  U[1].Xp[17].x =  0.990;  V[1].Xp[17].x = 0.81;  V[1].Xp[17].y = 0.88;
  U[1].Xp[18].y = 0.92;  U[1].Xp[18].x =  2.250;  V[1].Xp[18].x = 0.22;  V[1].Xp[18].y = 0.92;
  U[1].Xp[19].y = 0.96;  U[1].Xp[19].x =  5.300;  V[1].Xp[19].x = - 0.32;  V[1].Xp[19].y = 0.96;
  U[1].Xp[20].y = 1.00;  U[1].Xp[20].x =  9.950;  V[1].Xp[20].x = - 0.77;  V[1].Xp[20].y = 1.00;
  U[1].Xp[21].y = 1.04;  U[1].Xp[21].x = 14.980;  V[1].Xp[21].x = - 0.85;  V[1].Xp[21].y = 1.04;
  U[1].Xp[22].y = 1.12;  U[1].Xp[22].x = 17.610;  V[1].Xp[22].x = - 1.06;  V[1].Xp[22].y = 1.12;
  U[1].Xp[23].y = 1.20;  U[1].Xp[23].x = 19.710;  V[1].Xp[23].x = - 1.05;  V[1].Xp[23].y = 1.20;
  U[1].pathlength();                              V[1].pathlength();

  // STATION 3
  x[2] = 0.40;
  U[2].allocate(23);                              V[2].allocate(23);
  U[2].settype(SPLINE2D_LINEAR);               V[2].settype(SPLINE2D_LINEAR);
  U[2].Xp[ 0].y = 0.00;  U[2].Xp[ 0].x = 65.610;  V[2].Xp[ 0].x =  -0.16;  V[2].Xp[ 0].y = 0.00;
  U[2].Xp[ 1].y = 0.02;  U[2].Xp[ 1].x = 60.850;  V[2].Xp[ 1].x = 0.56;  V[2].Xp[ 1].y = 0.02;
  U[2].Xp[ 2].y = 0.04;  U[2].Xp[ 2].x = 52.750;  V[2].Xp[ 2].x = 1.29;  V[2].Xp[ 2].y = 0.04;
  U[2].Xp[ 3].y = 0.06;  U[2].Xp[ 3].x = 42.690;  V[2].Xp[ 3].x = 2.01;  V[2].Xp[ 3].y = 0.06;
  U[2].Xp[ 4].y = 0.08;  U[2].Xp[ 4].x = 32.550;  V[2].Xp[ 4].x = 1.98;  V[2].Xp[ 4].y = 0.08;
  U[2].Xp[ 5].y = 0.12;  U[2].Xp[ 5].x = 16.770;  V[2].Xp[ 5].x = 1.92;  V[2].Xp[ 5].y = 0.12;
  U[2].Xp[ 6].y = 0.16;  U[2].Xp[ 6].x =  7.320;  V[2].Xp[ 6].x = 1.33;  V[2].Xp[ 6].y = 0.16;
  U[2].Xp[ 7].y = 0.20;  U[2].Xp[ 7].x = -0.010;  V[2].Xp[ 7].x = 0.01;  V[2].Xp[ 7].y = 0.20;
  U[2].Xp[ 8].y = 0.24;  U[2].Xp[ 8].x = -3.850;  V[2].Xp[ 8].x = - 0.25;  V[2].Xp[ 8].y = 0.24;
  U[2].Xp[ 9].y = 0.32;  U[2].Xp[ 9].x = -5.550;  V[2].Xp[ 9].x = - 0.09;  V[2].Xp[ 9].y = 0.32;
  U[2].Xp[10].y = 0.40;  U[2].Xp[10].x = -5.810;  V[2].Xp[10].x = 0.54;  V[2].Xp[10].y = 0.40;
  U[2].Xp[11].y = 0.48;  U[2].Xp[11].x = -5.400;  V[2].Xp[11].x = 0.96;  V[2].Xp[11].y = 0.48;
  U[2].Xp[12].y = 0.56;  U[2].Xp[12].x = -4.570;  V[2].Xp[12].x = 1.28;  V[2].Xp[12].y = 0.56;
  U[2].Xp[13].y = 0.64;  U[2].Xp[13].x = -3.070;  V[2].Xp[13].x = 1.31;  V[2].Xp[13].y = 0.64;
  U[2].Xp[14].y = 0.72;  U[2].Xp[14].x = -1.150;  V[2].Xp[14].x = 1.22;  V[2].Xp[14].y = 0.72;
  U[2].Xp[15].y = 0.80;  U[2].Xp[15].x =  1.620;  V[2].Xp[15].x = 0.49;  V[2].Xp[15].y = 0.80;
  U[2].Xp[16].y = 0.88;  U[2].Xp[16].x =  4.650;  V[2].Xp[16].x = - 0.25;  V[2].Xp[16].y = 0.88;
  U[2].Xp[17].y = 0.92;  U[2].Xp[17].x =  7.270;  V[2].Xp[17].x = - 0.81;  V[2].Xp[17].y = 0.92;
  U[2].Xp[18].y = 0.96;  U[2].Xp[18].x = 10.050;  V[2].Xp[18].x = - 1.06;  V[2].Xp[18].y = 0.96;
  U[2].Xp[19].y = 1.00;  U[2].Xp[19].x = 13.030;  V[2].Xp[19].x = - 1.18;  V[2].Xp[19].y = 1.00;
  U[2].Xp[20].y = 1.04;  U[2].Xp[20].x = 16.160;  V[2].Xp[20].x = - 1.88;  V[2].Xp[20].y = 1.04;
  U[2].Xp[21].y = 1.12;  U[2].Xp[21].x = 19.210;  V[2].Xp[21].x = - 2.13;  V[2].Xp[21].y = 1.12;
  U[2].Xp[22].y = 1.20;  U[2].Xp[22].x = 19.680;  V[2].Xp[22].x = - 1.70;  V[2].Xp[22].y = 1.20;
  U[2].pathlength();                              V[2].pathlength();

  // STATION 4
  x[3] = 0.60;
  U[3].allocate(23);                             V[3].allocate(23);
  U[3].settype(SPLINE2D_LINEAR);              V[3].settype(SPLINE2D_LINEAR);
  U[3].Xp[ 0].x = 47.10;  U[3].Xp[ 0].y = 0.00;  V[3].Xp[ 0].x = 0.39;  V[3].Xp[ 0].y = 0.00;
  U[3].Xp[ 1].x = 45.13;  U[3].Xp[ 1].y = 0.02;  V[3].Xp[ 1].x = 1.54;  V[3].Xp[ 1].y = 0.02;
  U[3].Xp[ 2].x = 42.38;  U[3].Xp[ 2].y = 0.04;  V[3].Xp[ 2].x = 1.92;  V[3].Xp[ 2].y = 0.04;
  U[3].Xp[ 3].x = 38.22;  U[3].Xp[ 3].y = 0.06;	 V[3].Xp[ 3].x = 2.21;  V[3].Xp[ 3].y = 0.06;
  U[3].Xp[ 4].x = 32.19;  U[3].Xp[ 4].y = 0.08;	 V[3].Xp[ 4].x = 2.36;  V[3].Xp[ 4].y = 0.08;
  U[3].Xp[ 5].x = 21.03;  U[3].Xp[ 5].y = 0.12;	 V[3].Xp[ 5].x = 2.78;  V[3].Xp[ 5].y = 0.12;
  U[3].Xp[ 6].x = 12.85;  U[3].Xp[ 6].y = 0.16;	 V[3].Xp[ 6].x = 2.42;  V[3].Xp[ 6].y = 0.16;
  U[3].Xp[ 7].x =  6.98;  U[3].Xp[ 7].y = 0.20;  V[3].Xp[ 7].x = 2.41;  V[3].Xp[ 7].y = 0.20;
  U[3].Xp[ 8].x =  2.37;  U[3].Xp[ 8].y = 0.24;	 V[3].Xp[ 8].x = 1.96;  V[3].Xp[ 8].y = 0.24;
  U[3].Xp[ 9].x = -4.83;  U[3].Xp[ 9].y = 0.32;	 V[3].Xp[ 9].x = 1.02;  V[3].Xp[ 9].y = 0.32;
  U[3].Xp[10].x = -6.67;  U[3].Xp[10].y = 0.40;  V[3].Xp[10].x = 0.75;  V[3].Xp[10].y = 0.40;
  U[3].Xp[11].x = -6.51;  U[3].Xp[11].y = 0.48;	 V[3].Xp[11].x = 0.48;  V[3].Xp[11].y = 0.48;
  U[3].Xp[12].x = -5.09;  U[3].Xp[12].y = 0.56;	 V[3].Xp[12].x = - 0.02;  V[3].Xp[12].y = 0.56;
  U[3].Xp[13].x = -2.22;  U[3].Xp[13].y = 0.64;	 V[3].Xp[13].x = - 0.75;  V[3].Xp[13].y = 0.64;
  U[3].Xp[14].x =  1.80;  U[3].Xp[14].y = 0.72;	 V[3].Xp[14].x = - 1.44;  V[3].Xp[14].y = 0.72;
  U[3].Xp[15].x =  6.17;  U[3].Xp[15].y = 0.80;  V[3].Xp[15].x = - 2.53;  V[3].Xp[15].y = 0.80;
  U[3].Xp[16].x = 10.25;  U[3].Xp[16].y = 0.88;	 V[3].Xp[16].x = - 3.15;  V[3].Xp[16].y = 0.88;
  U[3].Xp[17].x = 11.99;  U[3].Xp[17].y = 0.92;	 V[3].Xp[17].x = - 3.11;  V[3].Xp[17].y = 0.92;
  U[3].Xp[18].x = 14.76;  U[3].Xp[18].y = 0.96;	 V[3].Xp[18].x = - 3.54;  V[3].Xp[18].y = 0.96;
  U[3].Xp[19].x = 16.42;  U[3].Xp[19].y = 1.00;  V[3].Xp[19].x = - 3.60;  V[3].Xp[19].y = 1.00;
  U[3].Xp[20].x = 18.14;  U[3].Xp[20].y = 1.04;	 V[3].Xp[20].x = - 3.61;  V[3].Xp[20].y = 1.04;
  U[3].Xp[21].x = 19.14;  U[3].Xp[21].y = 1.12;	 V[3].Xp[21].x = - 3.14;  V[3].Xp[21].y = 1.12;
  U[3].Xp[22].x = 19.75;  U[3].Xp[22].y = 1.20;  V[3].Xp[22].x = - 2.85;  V[3].Xp[22].y = 1.20;
  U[3].pathlength();                             V[3].pathlength();
  // STATION 5
  x[4] = 0.80;
  U[4].allocate(23);                             V[4].allocate(23);
  U[4].settype(SPLINE2D_LINEAR);              V[4].settype(SPLINE2D_LINEAR);
  U[4].Xp[ 0].x = 35.38;  U[4].Xp[ 0].y = 0.00;  V[4].Xp[ 0].x = 0.51;  V[4].Xp[ 0].y = 0.00;
  U[4].Xp[ 1].x = 34.13;  U[4].Xp[ 1].y = 0.02;  V[4].Xp[ 1].x = 1.54;  V[4].Xp[ 1].y = 0.02;
  U[4].Xp[ 2].x = 32.41;  U[4].Xp[ 2].y = 0.04;  V[4].Xp[ 2].x = 2.02;  V[4].Xp[ 2].y = 0.04;
  U[4].Xp[ 3].x = 29.78;  U[4].Xp[ 3].y = 0.06;  V[4].Xp[ 3].x = 2.46;  V[4].Xp[ 3].y = 0.06;
  U[4].Xp[ 4].x = 26.56;  U[4].Xp[ 4].y = 0.08;  V[4].Xp[ 4].x = 2.89;  V[4].Xp[ 4].y = 0.08;
  U[4].Xp[ 5].x = 19.81;  U[4].Xp[ 5].y = 0.12;  V[4].Xp[ 5].x = 3.34;  V[4].Xp[ 5].y = 0.12;
  U[4].Xp[ 6].x = 14.71;  U[4].Xp[ 6].y = 0.16;  V[4].Xp[ 6].x = 3.30;  V[4].Xp[ 6].y = 0.16;
  U[4].Xp[ 7].x = 10.24;  U[4].Xp[ 7].y = 0.20;  V[4].Xp[ 7].x = 2.98;  V[4].Xp[ 7].y = 0.20;
  U[4].Xp[ 8].x =  6.53;  U[4].Xp[ 8].y = 0.24;  V[4].Xp[ 8].x = 3.43;  V[4].Xp[ 8].y = 0.24;
  U[4].Xp[ 9].x = -1.28;  U[4].Xp[ 9].y = 0.32;  V[4].Xp[ 9].x = 1.82;  V[4].Xp[ 9].y = 0.32;
  U[4].Xp[10].x = -4.11;  U[4].Xp[10].y = 0.40;  V[4].Xp[10].x = 0.43;  V[4].Xp[10].y = 0.40;
  U[4].Xp[11].x = -3.45;  U[4].Xp[11].y = 0.48;  V[4].Xp[11].x = - 1.11;  V[4].Xp[11].y = 0.48;
  U[4].Xp[12].x = -0.75;  U[4].Xp[12].y = 0.56;  V[4].Xp[12].x = - 2.63;  V[4].Xp[12].y = 0.56;
  U[4].Xp[13].x =  3.83;  U[4].Xp[13].y = 0.64;  V[4].Xp[13].x = - 3.96;  V[4].Xp[13].y = 0.64;
  U[4].Xp[14].x =  7.54;  U[4].Xp[14].y = 0.72;  V[4].Xp[14].x = - 4.47;  V[4].Xp[14].y = 0.72;
  U[4].Xp[15].x = 10.46;  U[4].Xp[15].y = 0.80;  V[4].Xp[15].x = - 4.52;  V[4].Xp[15].y = 0.80;
  U[4].Xp[16].x = 14.05;  U[4].Xp[16].y = 0.88;  V[4].Xp[16].x = - 4.63;  V[4].Xp[16].y = 0.88;
  U[4].Xp[17].x = 15.61;  U[4].Xp[17].y = 0.92;  V[4].Xp[17].x = - 4.47;  V[4].Xp[17].y = 0.92;
  U[4].Xp[18].x = 16.69;  U[4].Xp[18].y = 0.96;  V[4].Xp[18].x = - 4.29;  V[4].Xp[18].y = 0.96;
  U[4].Xp[19].x = 17.40;  U[4].Xp[19].y = 1.00;  V[4].Xp[19].x = - 4.16;  V[4].Xp[19].y = 1.00;
  U[4].Xp[20].x = 18.01;  U[4].Xp[20].y = 1.04;  V[4].Xp[20].x = - 3.97;  V[4].Xp[20].y = 1.04;
  U[4].Xp[21].x = 18.73;  U[4].Xp[21].y = 1.12;  V[4].Xp[21].x = - 3.60;  V[4].Xp[21].y = 1.12;
  U[4].Xp[22].x = 19.23;  U[4].Xp[22].y = 1.20;  V[4].Xp[22].x = - 3.18;  V[4].Xp[22].y = 1.20;
  U[4].pathlength();                             V[4].pathlength();

/*   // STATION 6 */
  x[5] = 1.0;
  U[5].allocate(23);                             V[5].allocate(23);
  U[5].settype(SPLINE2D_LINEAR);              V[5].settype(SPLINE2D_LINEAR);
  U[5].Xp[ 0].x = 24.07;  U[5].Xp[ 0].y = 0.00;  V[5].Xp[ 0].x = 0.21;  V[5].Xp[ 0].y = 0.00;
  U[5].Xp[ 1].x = 23.91;  U[5].Xp[ 1].y = 0.02;  V[5].Xp[ 1].x = 0.98;  V[5].Xp[ 1].y = 0.02;
  U[5].Xp[ 2].x = 23.31;  U[5].Xp[ 2].y = 0.04;  V[5].Xp[ 2].x = 1.08;  V[5].Xp[ 2].y = 0.04;
  U[5].Xp[ 3].x = 22.48;  U[5].Xp[ 3].y = 0.06;  V[5].Xp[ 3].x = 1.69;  V[5].Xp[ 3].y = 0.06;
  U[5].Xp[ 4].x = 21.33;  U[5].Xp[ 4].y = 0.08;  V[5].Xp[ 4].x = 2.10;  V[5].Xp[ 4].y = 0.08;
  U[5].Xp[ 5].x = 18.18;  U[5].Xp[ 5].y = 0.12;  V[5].Xp[ 5].x = 2.93;  V[5].Xp[ 5].y = 0.12;
  U[5].Xp[ 6].x = 15.29;  U[5].Xp[ 6].y = 0.16;  V[5].Xp[ 6].x = 2.96;  V[5].Xp[ 6].y = 0.16;
  U[5].Xp[ 7].x = 12.01;  U[5].Xp[ 7].y = 0.20;  V[5].Xp[ 7].x = 3.26;  V[5].Xp[ 7].y = 0.20;
  U[5].Xp[ 8].x =  9.07;  U[5].Xp[ 8].y = 0.24;  V[5].Xp[ 8].x = 3.11;  V[5].Xp[ 8].y = 0.24;
  U[5].Xp[ 9].x =  4.10;  U[5].Xp[ 9].y = 0.32;  V[5].Xp[ 9].x = 2.39;  V[5].Xp[ 9].y = 0.32;
  U[5].Xp[10].x =  1.30;  U[5].Xp[10].y = 0.40;  V[5].Xp[10].x = 0.38;  V[5].Xp[10].y = 0.40;
  U[5].Xp[11].x =  2.12;  U[5].Xp[11].y = 0.46;  V[5].Xp[11].x = - 2.02;  V[5].Xp[11].y = 0.48;
  U[5].Xp[12].x =  5.52;  U[5].Xp[12].y = 0.56;  V[5].Xp[12].x = - 4.20;  V[5].Xp[12].y = 0.56;
  U[5].Xp[13].x =  7.95;  U[5].Xp[13].y = 0.64;  V[5].Xp[13].x = - 4.46;  V[5].Xp[13].y = 0.64;
  U[5].Xp[14].x = 10.85;  U[5].Xp[14].y = 0.72;  V[5].Xp[14].x = - 4.49;  V[5].Xp[14].y = 0.72;
  U[5].Xp[15].x = 13.09;  U[5].Xp[15].y = 0.80;  V[5].Xp[15].x = - 4.22;  V[5].Xp[15].y = 0.80;
  U[5].Xp[16].x = 15.39;  U[5].Xp[16].y = 0.88;  V[5].Xp[16].x = - 3.93;  V[5].Xp[16].y = 0.88;
  U[5].Xp[17].x = 16.39;  U[5].Xp[17].y = 0.92;  V[5].Xp[17].x = - 3.85;  V[5].Xp[17].y = 0.92;
  U[5].Xp[18].x = 16.72;  U[5].Xp[18].y = 0.96;  V[5].Xp[18].x = - 3.65;  V[5].Xp[18].y = 0.96;
  U[5].Xp[19].x = 17.34;  U[5].Xp[19].y = 1.00;  V[5].Xp[19].x = - 3.52;  V[5].Xp[19].y = 1.00;
  U[5].Xp[20].x = 17.77;  U[5].Xp[20].y = 1.04;  V[5].Xp[20].x = - 3.40;  V[5].Xp[20].y = 1.04;
  U[5].Xp[21].x = 18.56;  U[5].Xp[21].y = 1.12;  V[5].Xp[21].x = - 3.13;  V[5].Xp[21].y = 1.12;
  U[5].Xp[22].x = 19.06;  U[5].Xp[22].y = 1.20;  V[5].Xp[22].x = - 2.89;  V[5].Xp[22].y = 1.20;
  U[5].pathlength();                             V[5].pathlength();

  // STATION 7
  x[6] = 1.20;
  U[6].allocate(23);                             V[6].allocate(23);
  U[6].settype(SPLINE2D_LINEAR);              V[6].settype(SPLINE2D_LINEAR);
  U[6].Xp[ 0].x = 17.79;  U[6].Xp[ 0].y = 0.00;  V[6].Xp[ 0].x = 0.22;  V[6].Xp[ 0].y = 0.00;
  U[6].Xp[ 1].x = 17.63;  U[6].Xp[ 1].y = 0.02;  V[6].Xp[ 1].x = 0.81;  V[6].Xp[ 1].y = 0.02;
  U[6].Xp[ 2].x = 17.23;  U[6].Xp[ 2].y = 0.04;  V[6].Xp[ 2].x = 1.28;  V[6].Xp[ 2].y = 0.04;
  U[6].Xp[ 3].x = 17.02;  U[6].Xp[ 3].y = 0.06;  V[6].Xp[ 3].x = 1.47;  V[6].Xp[ 3].y = 0.06;
  U[6].Xp[ 4].x = 16.27;  U[6].Xp[ 4].y = 0.08;  V[6].Xp[ 4].x = 1.49;  V[6].Xp[ 4].y = 0.08;
  U[6].Xp[ 5].x = 15.26;  U[6].Xp[ 5].y = 0.12;  V[6].Xp[ 5].x = 2.03;  V[6].Xp[ 5].y = 0.12;
  U[6].Xp[ 6].x = 14.10;  U[6].Xp[ 6].y = 0.16;  V[6].Xp[ 6].x = 2.16;  V[6].Xp[ 6].y = 0.16;
  U[6].Xp[ 7].x = 12.19;  U[6].Xp[ 7].y = 0.20;  V[6].Xp[ 7].x = 2.44;  V[6].Xp[ 7].y = 0.20;
  U[6].Xp[ 8].x = 10.65;  U[6].Xp[ 8].y = 0.24;  V[6].Xp[ 8].x = 2.13;  V[6].Xp[ 8].y = 0.24;
  U[6].Xp[ 9].x =  7.70;  U[6].Xp[ 9].y = 0.32;  V[6].Xp[ 9].x = 1.73;  V[6].Xp[ 9].y = 0.32;
  U[6].Xp[10].x =  6.12;  U[6].Xp[10].y = 0.40;  V[6].Xp[10].x = 0.95;  V[6].Xp[10].y = 0.40;
  U[6].Xp[11].x =  7.22;  U[6].Xp[11].y = 0.48;  V[6].Xp[11].x = - 2.93;  V[6].Xp[11].y = 0.48;
  U[6].Xp[12].x =  8.99;  U[6].Xp[12].y = 0.56;  V[6].Xp[12].x = - 3.25;  V[6].Xp[12].y = 0.56;
  U[6].Xp[13].x = 10.82;  U[6].Xp[13].y = 0.64;  V[6].Xp[13].x = - 3.36;  V[6].Xp[13].y = 0.64;
  U[6].Xp[14].x = 13.14;  U[6].Xp[14].y = 0.72;  V[6].Xp[14].x = - 3.24;  V[6].Xp[14].y = 0.72;
  U[6].Xp[15].x = 14.78;  U[6].Xp[15].y = 0.80;  V[6].Xp[15].x = - 3.01;  V[6].Xp[15].y = 0.80;
  U[6].Xp[16].x = 15.92;  U[6].Xp[16].y = 0.88;  V[6].Xp[16].x = - 2.85;  V[6].Xp[16].y = 0.88;
  U[6].Xp[17].x = 16.59;  U[6].Xp[17].y = 0.92;  V[6].Xp[17].x = - 2.68;  V[6].Xp[17].y = 0.92;
  U[6].Xp[18].x = 16.78;  U[6].Xp[18].y = 0.96;  V[6].Xp[18].x = - 2.65;  V[6].Xp[18].y = 0.96;
  U[6].Xp[19].x = 17.27;  U[6].Xp[19].y = 1.00;  V[6].Xp[19].x = - 2.60;  V[6].Xp[19].y = 1.00;
  U[6].Xp[20].x = 17.65;  U[6].Xp[20].y = 1.04;  V[6].Xp[20].x = - 2.57;  V[6].Xp[20].y = 1.04;
  U[6].Xp[21].x = 17.94;  U[6].Xp[21].y = 1.12;  V[6].Xp[21].x = - 2.38;  V[6].Xp[21].y = 1.12;
  U[6].Xp[22].x = 18.67;  U[6].Xp[22].y = 1.20;  V[6].Xp[22].x = - 2.24;  V[6].Xp[22].y = 1.20;
  U[6].pathlength();                             V[6].pathlength();

  // STATION 8
  x[7] = 1.40;
  U[7].allocate(23);                             V[7].allocate(23);
  U[7].settype(SPLINE2D_LINEAR);              V[7].settype(SPLINE2D_LINEAR);
  U[7].Xp[ 0].x = 13.71;  U[7].Xp[ 0].y = 0.00;  V[7].Xp[ 0].x = - 0.10;  V[7].Xp[ 0].y = 0.00;
  U[7].Xp[ 1].x = 13.62;  U[7].Xp[ 1].y = 0.02;  V[7].Xp[ 1].x = 0.03;  V[7].Xp[ 1].y = 0.02;
  U[7].Xp[ 2].x = 13.72;  U[7].Xp[ 2].y = 0.04;  V[7].Xp[ 2].x = 0.25;  V[7].Xp[ 2].y = 0.04;
  U[7].Xp[ 3].x = 13.66;  U[7].Xp[ 3].y = 0.06;  V[7].Xp[ 3].x = 0.35;  V[7].Xp[ 3].y = 0.06;
  U[7].Xp[ 4].x = 13.46;  U[7].Xp[ 4].y = 0.08;  V[7].Xp[ 4].x = 0.55;  V[7].Xp[ 4].y = 0.08;
  U[7].Xp[ 5].x = 13.08;  U[7].Xp[ 5].y = 0.12;  V[7].Xp[ 5].x = 0.80;  V[7].Xp[ 5].y = 0.12;
  U[7].Xp[ 6].x = 12.28;  U[7].Xp[ 6].y = 0.16;  V[7].Xp[ 6].x = 0.87;  V[7].Xp[ 6].y = 0.16;
  U[7].Xp[ 7].x = 11.80;  U[7].Xp[ 7].y = 0.20;  V[7].Xp[ 7].x = 1.22;  V[7].Xp[ 7].y = 0.20;
  U[7].Xp[ 8].x = 11.11;  U[7].Xp[ 8].y = 0.24;  V[7].Xp[ 8].x = 1.42;  V[7].Xp[ 8].y = 0.24;
  U[7].Xp[ 9].x =  9.62;  U[7].Xp[ 9].y = 0.32;  V[7].Xp[ 9].x = 1.33;  V[7].Xp[ 9].y = 0.32;
  U[7].Xp[10].x =  8.55;  U[7].Xp[10].y = 0.40;  V[7].Xp[10].x = - 0.31;  V[7].Xp[10].y = 0.40;
  U[7].Xp[11].x =  9.09;  U[7].Xp[11].y = 0.48;  V[7].Xp[11].x = - 2.10;  V[7].Xp[11].y = 0.48;
  U[7].Xp[12].x =  9.98;  U[7].Xp[12].y = 0.56;  V[7].Xp[12].x = - 1.82;  V[7].Xp[12].y = 0.56;
  U[7].Xp[13].x = 11.58;  U[7].Xp[13].y = 0.64;  V[7].Xp[13].x = - 1.88;  V[7].Xp[13].y = 0.64;
  U[7].Xp[14].x = 13.75;  U[7].Xp[14].y = 0.72;  V[7].Xp[14].x = - 2.23;  V[7].Xp[14].y = 0.72;
  U[7].Xp[15].x = 15.28;  U[7].Xp[15].y = 0.80;  V[7].Xp[15].x = - 2.04;  V[7].Xp[15].y = 0.80;
  U[7].Xp[16].x = 16.59;  U[7].Xp[16].y = 0.88;  V[7].Xp[16].x = - 1.95;  V[7].Xp[16].y = 0.88;
  U[7].Xp[17].x = 17.22;  U[7].Xp[17].y = 0.92;  V[7].Xp[17].x = - 1.85;  V[7].Xp[17].y = 0.92;
  U[7].Xp[18].x = 17.63;  U[7].Xp[18].y = 0.96;  V[7].Xp[18].x = - 1.87;  V[7].Xp[18].y = 0.96;
  U[7].Xp[19].x = 17.85;  U[7].Xp[19].y = 1.00;  V[7].Xp[19].x = - 1.81;  V[7].Xp[19].y = 1.00;
  U[7].Xp[20].x = 18.14;  U[7].Xp[20].y = 1.04;  V[7].Xp[20].x = - 1.80;  V[7].Xp[20].y = 1.04;
  U[7].Xp[21].x = 18.70;  U[7].Xp[21].y = 1.12;  V[7].Xp[21].x = - 1.66;  V[7].Xp[21].y = 1.12;
  U[7].Xp[22].x = 19.21;  U[7].Xp[22].y = 1.20;  V[7].Xp[22].x = - 1.60;  V[7].Xp[22].y = 1.20;
  U[7].pathlength();                             V[7].pathlength();

  // STATION 9
  x[8] = 1.8;
  U[8].allocate(23);                             V[8].allocate(23);
  U[8].settype(SPLINE2D_LINEAR);              V[8].settype(SPLINE2D_LINEAR);
  U[8].Xp[ 0].x = 11.91;  U[8].Xp[ 0].y = 0.00;  V[8].Xp[ 0].x = 0.02;  V[8].Xp[ 0].y = 0.00;
  U[8].Xp[ 1].x = 11.95;  U[8].Xp[ 1].y = 0.02;  V[8].Xp[ 1].x = 0.15;  V[8].Xp[ 1].y = 0.02;
  U[8].Xp[ 2].x = 11.86;  U[8].Xp[ 2].y = 0.04;  V[8].Xp[ 2].x = 0.17;  V[8].Xp[ 2].y = 0.04;
  U[8].Xp[ 3].x = 11.76;  U[8].Xp[ 3].y = 0.06;  V[8].Xp[ 3].x = 0.13;  V[8].Xp[ 3].y = 0.06;
  U[8].Xp[ 4].x = 11.76;  U[8].Xp[ 4].y = 0.08;  V[8].Xp[ 4].x = 0.27;  V[8].Xp[ 4].y = 0.08;
  U[8].Xp[ 5].x = 11.79;  U[8].Xp[ 5].y = 0.12;  V[8].Xp[ 5].x = 0.32;  V[8].Xp[ 5].y = 0.12;
  U[8].Xp[ 6].x = 11.70;  U[8].Xp[ 6].y = 0.16;  V[8].Xp[ 6].x = 0.45;  V[8].Xp[ 6].y = 0.16;
  U[8].Xp[ 7].x = 11.58;  U[8].Xp[ 7].y = 0.20;  V[8].Xp[ 7].x = 0.53;  V[8].Xp[ 7].y = 0.20;
  U[8].Xp[ 8].x = 11.49;  U[8].Xp[ 8].y = 0.24;  V[8].Xp[ 8].x = 0.57;  V[8].Xp[ 8].y = 0.24;
  U[8].Xp[ 9].x = 11.36;  U[8].Xp[ 9].y = 0.32;  V[8].Xp[ 9].x = 0.55;  V[8].Xp[ 9].y = 0.32;
  U[8].Xp[10].x = 11.52;  U[8].Xp[10].y = 0.40;  V[8].Xp[10].x = - 0.42;  V[8].Xp[10].y = 0.40;
  U[8].Xp[11].x = 12.26;  U[8].Xp[11].y = 0.48;  V[8].Xp[11].x = - 0.82;  V[8].Xp[11].y = 0.48;
  U[8].Xp[12].x = 13.19;  U[8].Xp[12].y = 0.56;  V[8].Xp[12].x = - 1.07;  V[8].Xp[12].y = 0.56;
  U[8].Xp[13].x = 14.04;  U[8].Xp[13].y = 0.64;  V[8].Xp[13].x = - 0.89;  V[8].Xp[13].y = 0.64;
  U[8].Xp[14].x = 15.40;  U[8].Xp[14].y = 0.72;  V[8].Xp[14].x = - 0.97;  V[8].Xp[14].y = 0.72;
  U[8].Xp[15].x = 16.30;  U[8].Xp[15].y = 0.80;  V[8].Xp[15].x = - 0.89;  V[8].Xp[15].y = 0.80;
  U[8].Xp[16].x = 17.20;  U[8].Xp[16].y = 0.88;  V[8].Xp[16].x = - 0.84;  V[8].Xp[16].y = 0.88;
  U[8].Xp[17].x = 17.53;  U[8].Xp[17].y = 0.92;  V[8].Xp[17].x = - 0.70;  V[8].Xp[17].y = 0.92;
  U[8].Xp[18].x = 17.86;  U[8].Xp[18].y = 0.96;  V[8].Xp[18].x = - 0.64;  V[8].Xp[18].y = 0.96;
  U[8].Xp[19].x = 18.32;  U[8].Xp[19].y = 1.00;  V[8].Xp[19].x = - 0.88;  V[8].Xp[19].y = 1.00;
  U[8].Xp[20].x = 18.49;  U[8].Xp[20].y = 1.04;  V[8].Xp[20].x = - 0.72;  V[8].Xp[20].y = 1.04;
  U[8].Xp[21].x = 18.91;  U[8].Xp[21].y = 1.12;  V[8].Xp[21].x = - 0.75;  V[8].Xp[21].y = 1.12;
  U[8].Xp[22].x = 19.15;  U[8].Xp[22].y = 1.20;  V[8].Xp[22].x = - 0.71;  V[8].Xp[22].y = 1.20;
  U[8].pathlength();                             V[8].pathlength();

  // STATION 10
  x[ 9] = 2.40;
  U[ 9].allocate(23);                              V[ 9].allocate(23);
  U[ 9].settype(SPLINE2D_LINEAR);               V[ 9].settype(SPLINE2D_LINEAR);
  U[ 9].Xp[ 0].x = 13.05;  U[ 9].Xp[ 0].y = 0.00;  V[ 9].Xp[ 0].x = 0.14;  V[ 9].Xp[ 0].y = 0.00;
  U[ 9].Xp[ 1].x = 13.08;  U[ 9].Xp[ 1].y = 0.02;  V[ 9].Xp[ 1].x = 0.15;  V[ 9].Xp[ 1].y = 0.02;
  U[ 9].Xp[ 2].x = 12.98;  U[ 9].Xp[ 2].y = 0.04;  V[ 9].Xp[ 2].x = 0.20;  V[ 9].Xp[ 2].y = 0.04;
  U[ 9].Xp[ 3].x = 13.06;  U[ 9].Xp[ 3].y = 0.06;  V[ 9].Xp[ 3].x = 0.16;  V[ 9].Xp[ 3].y = 0.06;
  U[ 9].Xp[ 4].x = 12.99;  U[ 9].Xp[ 4].y = 0.08;  V[ 9].Xp[ 4].x = 0.18;  V[ 9].Xp[ 4].y = 0.08;
  U[ 9].Xp[ 5].x = 13.03;  U[ 9].Xp[ 5].y = 0.12;  V[ 9].Xp[ 5].x = 0.29;  V[ 9].Xp[ 5].y = 0.12;
  U[ 9].Xp[ 6].x = 13.01;  U[ 9].Xp[ 6].y = 0.16;  V[ 9].Xp[ 6].x = 0.43;  V[ 9].Xp[ 6].y = 0.16;
  U[ 9].Xp[ 7].x = 13.02;  U[ 9].Xp[ 7].y = 0.20;  V[ 9].Xp[ 7].x = 0.40;  V[ 9].Xp[ 7].y = 0.20;
  U[ 9].Xp[ 8].x = 13.04;  U[ 9].Xp[ 8].y = 0.24;  V[ 9].Xp[ 8].x = 0.54;  V[ 9].Xp[ 8].y = 0.24;
  U[ 9].Xp[ 9].x = 13.20;  U[ 9].Xp[ 9].y = 0.32;  V[ 9].Xp[ 9].x = 0.69;  V[ 9].Xp[ 9].y = 0.32;
  U[ 9].Xp[10].x = 13.34;  U[ 9].Xp[10].y = 0.40;  V[ 9].Xp[10].x = 0.86;  V[ 9].Xp[10].y = 0.40;
  U[ 9].Xp[11].x = 13.59;  U[ 9].Xp[11].y = 0.48;  V[ 9].Xp[11].x = 0.92;  V[ 9].Xp[11].y = 0.48;
  U[ 9].Xp[12].x = 14.34;  U[ 9].Xp[12].y = 0.56;  V[ 9].Xp[12].x = 0.25;  V[ 9].Xp[12].y = 0.56;
  U[ 9].Xp[13].x = 15.37;  U[ 9].Xp[13].y = 0.64;  V[ 9].Xp[13].x = - 0.43;  V[ 9].Xp[13].y = 0.64;
  U[ 9].Xp[14].x = 16.07;  U[ 9].Xp[14].y = 0.72;  V[ 9].Xp[14].x = - 0.36;  V[ 9].Xp[14].y = 0.72;
  U[ 9].Xp[15].x = 16.86;  U[ 9].Xp[15].y = 0.80;  V[ 9].Xp[15].x = - 0.42;  V[ 9].Xp[15].y = 0.80;
  U[ 9].Xp[16].x = 17.28;  U[ 9].Xp[16].y = 0.88;  V[ 9].Xp[16].x = - 0.32;  V[ 9].Xp[16].y = 0.88;
  U[ 9].Xp[17].x = 17.62;  U[ 9].Xp[17].y = 0.92;  V[ 9].Xp[17].x = - 0.34;  V[ 9].Xp[17].y = 0.92;
  U[ 9].Xp[18].x = 18.00;  U[ 9].Xp[18].y = 0.96;  V[ 9].Xp[18].x = - 0.30;  V[ 9].Xp[18].y = 0.96;
  U[ 9].Xp[19].x = 18.22;  U[ 9].Xp[19].y = 1.00;  V[ 9].Xp[19].x = - 0.24;  V[ 9].Xp[19].y = 1.00;
  U[ 9].Xp[20].x = 18.33;  U[ 9].Xp[20].y = 1.04;  V[ 9].Xp[20].x = - 0.30;  V[ 9].Xp[20].y = 1.04;
  U[ 9].Xp[21].x = 18.77;  U[ 9].Xp[21].y = 1.12;  V[ 9].Xp[21].x = - 0.18;  V[ 9].Xp[21].y = 1.12;
  U[ 9].Xp[22].x = 19.17;  U[ 9].Xp[22].y = 1.20;  V[ 9].Xp[22].x = - 0.16;  V[ 9].Xp[22].y = 1.20;
  U[ 9].pathlength();                              V[ 9].pathlength();

  // STATION 11
  x[10] = 3.40;
  U[10].allocate(23);                              V[10].allocate(23);
  U[10].settype(SPLINE2D_LINEAR);               V[10].settype(SPLINE2D_LINEAR);
  U[10].Xp[ 0].x = 14.49;  U[10].Xp[ 0].y = 0.00;  V[10].Xp[ 0].x =  -0.74;  V[10].Xp[ 0].y = 0.00;
  U[10].Xp[ 1].x = 14.42;  U[10].Xp[ 1].y = 0.02;  V[10].Xp[ 1].x =  -0.31;  V[10].Xp[ 1].y = 0.02;
  U[10].Xp[ 2].x = 14.58;  U[10].Xp[ 2].y = 0.04;  V[10].Xp[ 2].x =  -0.73;  V[10].Xp[ 2].y = 0.04;
  U[10].Xp[ 3].x = 14.52;  U[10].Xp[ 3].y = 0.06;  V[10].Xp[ 3].x =  -0.59;  V[10].Xp[ 3].y = 0.06;
  U[10].Xp[ 4].x = 14.58;  U[10].Xp[ 4].y = 0.08;  V[10].Xp[ 4].x =  -0.43;  V[10].Xp[ 4].y = 0.08;
  U[10].Xp[ 5].x = 14.74;  U[10].Xp[ 5].y = 0.12;  V[10].Xp[ 5].x =  -0.71;  V[10].Xp[ 5].y = 0.12;
  U[10].Xp[ 6].x = 14.83;  U[10].Xp[ 6].y = 0.16;  V[10].Xp[ 6].x =  -0.69;  V[10].Xp[ 6].y = 0.16;
  U[10].Xp[ 7].x = 14.70;  U[10].Xp[ 7].y = 0.20;  V[10].Xp[ 7].x =  -0.50;  V[10].Xp[ 7].y = 0.20;
  U[10].Xp[ 8].x = 14.84;  U[10].Xp[ 8].y = 0.24;  V[10].Xp[ 8].x =  -0.53;  V[10].Xp[ 8].y = 0.24;
  U[10].Xp[ 9].x = 14.86;  U[10].Xp[ 9].y = 0.32;  V[10].Xp[ 9].x =  -0.27;  V[10].Xp[ 9].y = 0.32;
  U[10].Xp[10].x = 15.02;  U[10].Xp[10].y = 0.40;  V[10].Xp[10].x =  -0.29;  V[10].Xp[10].y = 0.40;
  U[10].Xp[11].x = 15.04;  U[10].Xp[11].y = 0.48;  V[10].Xp[11].x = 0.04;  V[10].Xp[11].y = 0.48;
  U[10].Xp[12].x = 15.39;  U[10].Xp[12].y = 0.56;  V[10].Xp[12].x = - 0.13;  V[10].Xp[12].y = 0.56;
  U[10].Xp[13].x = 15.88;  U[10].Xp[13].y = 0.64;  V[10].Xp[13].x = - 0.11;  V[10].Xp[13].y = 0.64;
  U[10].Xp[14].x = 16.40;  U[10].Xp[14].y = 0.72;  V[10].Xp[14].x = - 0.16;  V[10].Xp[14].y = 0.72;
  U[10].Xp[15].x = 16.98;  U[10].Xp[15].y = 0.80;  V[10].Xp[15].x = - 0.11;  V[10].Xp[15].y = 0.80;
  U[10].Xp[16].x = 17.41;  U[10].Xp[16].y = 0.88;  V[10].Xp[16].x = - 0.05;  V[10].Xp[16].y = 0.88;
  U[10].Xp[17].x = 17.75;  U[10].Xp[17].y = 0.92;  V[10].Xp[17].x = - 0.05;  V[10].Xp[17].y = 0.92;
  U[10].Xp[18].x = 17.90;  U[10].Xp[18].y = 0.96;  V[10].Xp[18].x = - 0.02;  V[10].Xp[18].y = 0.96;
  U[10].Xp[19].x = 18.23;  U[10].Xp[19].y = 1.00;  V[10].Xp[19].x = - 0.03;  V[10].Xp[19].y = 1.00;
  U[10].Xp[20].x = 18.42;  U[10].Xp[20].y = 1.04;  V[10].Xp[20].x = - 0.01;  V[10].Xp[20].y = 1.04;
  U[10].Xp[21].x = 18.80;  U[10].Xp[21].y = 1.12;  V[10].Xp[21].x = 0.01;  V[10].Xp[21].y = 1.12;
  U[10].Xp[22].x = 18.99;  U[10].Xp[22].y = 1.20;  V[10].Xp[22].x = 0.01;  V[10].Xp[22].y = 1.20;
  U[10].pathlength();                              V[10].pathlength();

  // STATION 12
  x[11] = 4.40;
  U[11].allocate(20);                              V[11].allocate(20);
  U[11].settype(SPLINE2D_LINEAR);               V[11].settype(SPLINE2D_LINEAR);
  U[11].Xp[ 0].x = 15.44;  U[11].Xp[ 0].y = 0.00;  V[11].Xp[ 0].x = - 0.08;  V[11].Xp[ 0].y = 0.00;
  U[11].Xp[ 1].x = 15.35;  U[11].Xp[ 1].y = 0.04;  V[11].Xp[ 1].x = 0.05;  V[11].Xp[ 1].y = 0.04;
  U[11].Xp[ 2].x = 15.35;  U[11].Xp[ 2].y = 0.08;  V[11].Xp[ 2].x = 0.07;  V[11].Xp[ 2].y = 0.08;
  U[11].Xp[ 3].x = 15.36;  U[11].Xp[ 3].y = 0.12;  V[11].Xp[ 3].x = - 0.15;  V[11].Xp[ 3].y = 0.12;
  U[11].Xp[ 4].x = 15.52;  U[11].Xp[ 4].y = 0.16;  V[11].Xp[ 4].x = - 0.27;  V[11].Xp[ 4].y = 0.16;
  U[11].Xp[ 5].x = 15.74;  U[11].Xp[ 5].y = 0.24;  V[11].Xp[ 5].x = - 0.50;  V[11].Xp[ 5].y = 0.24;
  U[11].Xp[ 6].x = 15.83;  U[11].Xp[ 6].y = 0.32;  V[11].Xp[ 6].x = - 0.46;  V[11].Xp[ 6].y = 0.32;
  U[11].Xp[ 7].x = 15.90;  U[11].Xp[ 7].y = 0.40;  V[11].Xp[ 7].x = - 0.33;  V[11].Xp[ 7].y = 0.40;
  U[11].Xp[ 8].x = 16.09;  U[11].Xp[ 8].y = 0.48;  V[11].Xp[ 8].x = - 0.34;  V[11].Xp[ 8].y = 0.48;
  U[11].Xp[ 9].x = 16.18;  U[11].Xp[ 9].y = 0.56;  V[11].Xp[ 9].x = - 0.37;  V[11].Xp[ 9].y = 0.56;
  U[11].Xp[10].x = 16.42;  U[11].Xp[10].y = 0.64;  V[11].Xp[10].x = - 0.23;  V[11].Xp[10].y = 0.64;
  U[11].Xp[11].x = 16.77;  U[11].Xp[11].y = 0.72;  V[11].Xp[11].x = - 0.25;  V[11].Xp[11].y = 0.72;
  U[11].Xp[12].x = 17.02;  U[11].Xp[12].y = 0.80;  V[11].Xp[12].x = - 0.21;  V[11].Xp[12].y = 0.80;
  U[11].Xp[13].x = 17.34;  U[11].Xp[13].y = 0.88;  V[11].Xp[13].x = - 0.12;  V[11].Xp[13].y = 0.88;
  U[11].Xp[14].x = 17.55;  U[11].Xp[14].y = 0.92;  V[11].Xp[14].x = - 0.16;  V[11].Xp[14].y = 0.92;
  U[11].Xp[15].x = 17.96;  U[11].Xp[15].y = 0.96;  V[11].Xp[15].x = - 0.11;  V[11].Xp[15].y = 0.96;
  U[11].Xp[16].x = 17.97;  U[11].Xp[16].y = 1.00;  V[11].Xp[16].x = - 0.12;  V[11].Xp[16].y = 1.00;
  U[11].Xp[17].x = 18.21;  U[11].Xp[17].y = 1.04;  V[11].Xp[17].x = - 0.06;  V[11].Xp[17].y = 1.04;
  U[11].Xp[18].x = 18.57;  U[11].Xp[18].y = 1.12;  V[11].Xp[18].x = - 0.02;  V[11].Xp[18].y = 1.12;
  U[11].Xp[19].x = 18.84;  U[11].Xp[19].y = 1.20;  V[11].Xp[19].x = 0.06;  V[11].Xp[19].y = 1.20;
  U[11].pathlength();                              V[11].pathlength();

  // STATION 13
  x[12] = 5.20;
  U[12].allocate(20);                              V[12].allocate(20);
  U[12].settype(SPLINE2D_LINEAR);               V[12].settype(SPLINE2D_LINEAR);
  U[12].Xp[ 0].x = 15.76;  U[12].Xp[ 0].y = 0.00;  V[12].Xp[ 0].x = 0.07;  V[12].Xp[ 0].y = 0.00;
  U[12].Xp[ 1].x = 15.83;  U[12].Xp[ 1].y = 0.04;  V[12].Xp[ 1].x = 0.01;  V[12].Xp[ 1].y = 0.04;
  U[12].Xp[ 2].x = 15.79;  U[12].Xp[ 2].y = 0.08;  V[12].Xp[ 2].x = - 0.44;  V[12].Xp[ 2].y = 0.08;
  U[12].Xp[ 3].x = 15.90;  U[12].Xp[ 3].y = 0.12;  V[12].Xp[ 3].x = - 0.25;  V[12].Xp[ 3].y = 0.12;
  U[12].Xp[ 4].x = 15.88;  U[12].Xp[ 4].y = 0.16;  V[12].Xp[ 4].x = - 0.29;  V[12].Xp[ 4].y = 0.16;
  U[12].Xp[ 5].x = 15.99;  U[12].Xp[ 5].y = 0.24;  V[12].Xp[ 5].x = - 0.44;  V[12].Xp[ 5].y = 0.24;
  U[12].Xp[ 6].x = 16.04;  U[12].Xp[ 6].y = 0.32;  V[12].Xp[ 6].x = - 0.39;  V[12].Xp[ 6].y = 0.32;
  U[12].Xp[ 7].x = 16.15;  U[12].Xp[ 7].y = 0.40;  V[12].Xp[ 7].x = - 0.34;  V[12].Xp[ 7].y = 0.40;
  U[12].Xp[ 8].x = 16.26;  U[12].Xp[ 8].y = 0.48;  V[12].Xp[ 8].x = - 0.24;  V[12].Xp[ 8].y = 0.48;
  U[12].Xp[ 9].x = 16.39;  U[12].Xp[ 9].y = 0.56;  V[12].Xp[ 9].x = - 0.23;  V[12].Xp[ 9].y = 0.56;
  U[12].Xp[10].x = 16.65;  U[12].Xp[10].y = 0.64;  V[12].Xp[10].x = - 0.26;  V[12].Xp[10].y = 0.64;
  U[12].Xp[11].x = 16.84;  U[12].Xp[11].y = 0.72;  V[12].Xp[11].x = - 0.17;  V[12].Xp[11].y = 0.72;
  U[12].Xp[12].x = 16.99;  U[12].Xp[12].y = 0.80;  V[12].Xp[12].x = - 0.08;  V[12].Xp[12].y = 0.80;
  U[12].Xp[13].x = 17.55;  U[12].Xp[13].y = 0.88;  V[12].Xp[13].x = - 0.15;  V[12].Xp[13].y = 0.88;
  U[12].Xp[14].x = 17.31;  U[12].Xp[14].y = 0.92;  V[12].Xp[14].x = - 0.06;  V[12].Xp[14].y = 0.92;
  U[12].Xp[15].x = 17.66;  U[12].Xp[15].y = 0.96;  V[12].Xp[15].x = - 0.14;  V[12].Xp[15].y = 0.96;
  U[12].Xp[16].x = 17.79;  U[12].Xp[16].y = 1.00;  V[12].Xp[16].x = - 0.06;  V[12].Xp[16].y = 1.00;
  U[12].Xp[17].x = 18.08;  U[12].Xp[17].y = 1.04;  V[12].Xp[17].x = - 0.08;  V[12].Xp[17].y = 1.04;
  U[12].Xp[18].x = 18.39;  U[12].Xp[18].y = 1.12;  V[12].Xp[18].x = - 0.10;  V[12].Xp[18].y = 1.12;
  U[12].Xp[19].x = 18.56;  U[12].Xp[19].y = 1.20;  V[12].Xp[19].x = 0.07;  V[12].Xp[19].y = 1.20;
  U[12].pathlength();                              V[12].pathlength();

 }

//! NonreactiveVelocityField::deallocate -- Deallocate memory for
//!                                            splines and arrays.
inline void NonreactiveVelocityField::deallocate(void) {
  // Deallocate splines.
  for (int n = 0; n < Ns; n++) {
    U[n].deallocate();
    V[n].deallocate();
  }
  // Deallocate arrays.
  delete []U; U = NULL;
  delete []V; V = NULL;
  delete []x; x = NULL;
}

/* Private member functions
 *      Ns -- Number of stations.
 *      M  -- Array of mixture fraction splines at each station.
 *      x  -- Array of station x-coordinates (axial direction).
 * Public member functions
 *      allocate      -- Allocate memory and initialize arrays.
 *      deallocate    -- Deallocate memory for splines and arrays.
 *      interpolation -- Interpolation function.
 * \endverbatim
 */
class NonreactiveScalarField {
 private:
   int       Ns; //!< Number of stations.
   Spline2D  *M; //!< Array of mixture fraction splines at each station.
   double    *x; //!< Array of station x-coordinates (axial direction).
  
 public:
  //! Creation constructor.
  NonreactiveScalarField() { allocate(); }

  //! Destructor.
  ~NonreactiveScalarField(void) { deallocate(); }

  //! Allocate memory and initialize arrays.
  void allocate(void);

  //! Deallocate memory for splines and arrays.
  void deallocate(void);

  //! Interpolation function.
  Vector2D interpolation(const Vector2D &Xt);

};

//! NonreactiveVelocityField::allocate -- Allocate memory and
//!                                          initialize arrays.
inline void NonreactiveScalarField::allocate(void) {
/* Velocity of coflowing air is 20 m/s. */
/* Velocity of coflowing fuel (in this case "CH4" ) is 50 m/s . */
 // mixing field ---  mixture fraction (CH4)
  // Set the number of stations.
  Ns = 54;
  // Allocate the arrays.
  M = new Spline2D[Ns];
  x = new double[Ns];
/* //station 1 */
/*   x[0] = 20.0; */
/*   M[0].allocate(44); */
/*   M[0].settype(SPLINE2D_LINEAR); */
/*   M[0].Xp[0].x =0.829	;	M[0].Xp[0].y =0	; */
/*   M[0].Xp[1].x =0.687	;	M[0].Xp[1].y =0.6	; */
/*   M[0].Xp[2].x =0.631	;	M[0].Xp[2].y =1.2	; */
/*   M[0].Xp[3].x =0.465	;	M[0].Xp[3].y =1.8	; */
/*   M[0].Xp[4].x =0.458	;	M[0].Xp[4].y =2.4	; */
/*   M[0].Xp[5].x =0.344	;	M[0].Xp[5].y =3	; */
/*   M[0].Xp[6].x =0.259	;	M[0].Xp[6].y =3.6	; */
/*   M[0].Xp[7].x =0.265	;	M[0].Xp[7].y =4.2	; */
/*   M[0].Xp[8].x =0.239	;	M[0].Xp[8].y =4.8	; */
/*   M[0].Xp[9].x =0.218	;	M[0].Xp[9].y =5.4	; */
/*   M[0].Xp[10].x =0.16	;	M[0].Xp[10].y =6	; */
/*   M[0].Xp[11].x =0.165	;	M[0].Xp[11].y =6.6	; */
/*   M[0].Xp[12].x =0.179	;	M[0].Xp[12].y =7.2	; */
/*   M[0].Xp[13].x =0.162	;	M[0].Xp[13].y =7.8	; */
/*   M[0].Xp[14].x =0.139	;	M[0].Xp[14].y =8.4	; */
/*   M[0].Xp[15].x =0.135	;	M[0].Xp[15].y =9	; */
/*   M[0].Xp[16].x =0.146	;	M[0].Xp[16].y =9.6	; */
/*   M[0].Xp[17].x =0.155	;	M[0].Xp[17].y =10.2	; */
/*   M[0].Xp[18].x =0.139	;	M[0].Xp[18].y =10.8	; */
/*   M[0].Xp[19].x =0.127	;	M[0].Xp[19].y =11.4	; */
/*   M[0].Xp[20].x =0.155	;	M[0].Xp[20].y =12	; */
/*   M[0].Xp[21].x =0.128	;	M[0].Xp[21].y =12.6	; */
/*   M[0].Xp[22].x =0.13	;	M[0].Xp[22].y =13.2	; */
/*   M[0].Xp[23].x =0.115	;	M[0].Xp[23].y =13.8	; */
/*   M[0].Xp[24].x =0.15	;	M[0].Xp[24].y =14.4	; */
/*   M[0].Xp[25].x =0.101	;	M[0].Xp[25].y =15	; */
/*   M[0].Xp[26].x =0.139	;	M[0].Xp[26].y =15.6	; */
/*   M[0].Xp[27].x =0.111	;	M[0].Xp[27].y =16.2	; */
/*   M[0].Xp[28].x =0.108	;	M[0].Xp[28].y =16.8	; */
/*   M[0].Xp[29].x =0.121	;	M[0].Xp[29].y =17.4	; */
/*   M[0].Xp[30].x =0.111	;	M[0].Xp[30].y =18	; */
/*   M[0].Xp[31].x =0.124	;	M[0].Xp[31].y =18.6	; */
/*   M[0].Xp[32].x =0.108	;	M[0].Xp[32].y =19.2	; */
/*   M[0].Xp[33].x =0.093	;	M[0].Xp[33].y =19.8	; */
/*   M[0].Xp[34].x =0.1	;	M[0].Xp[34].y =20.4	; */
/*   M[0].Xp[35].x =0.076	;	M[0].Xp[35].y =21	; */
/*   M[0].Xp[36].x =0.065	;	M[0].Xp[36].y =21.6	; */
/*   M[0].Xp[37].x =0.049	;	M[0].Xp[37].y =22.2	; */
/*   M[0].Xp[38].x =0.012	;	M[0].Xp[38].y =22.8	; */
/*   M[0].Xp[39].x =0.007	;	M[0].Xp[39].y =23.4	; */
/*   M[0].Xp[40].x =0.008	;	M[0].Xp[40].y =24	; */
/*   M[0].Xp[41].x =0.004	;	M[0].Xp[41].y =24.6	; */
/*   M[0].Xp[42].x =0.007	;	M[0].Xp[42].y =25.2	; */
/*   M[0].Xp[43].x =0.002	;	M[0].Xp[43].y =25.8	; */
/*   M[0].pathlength(); */
/*   //Station 2 */
/*   x[1] = 20.80; */
/*   M[1].allocate(44); */
/*   M[1].settype(SPLINE2D_LINEAR); */
/*   M[1].Xp[0].x=0.733	;	M[1].Xp[0].y =0	; */
/*   M[1].Xp[1].x=0.668	;	M[1].Xp[1].y =0.6	; */
/*   M[1].Xp[2].x=0.564	;	M[1].Xp[2].y =1.2	; */
/*   M[1].Xp[3].x=0.524	;	M[1].Xp[3].y =1.8	; */
/*   M[1].Xp[4].x=0.385	;	M[1].Xp[4].y =2.4	; */
/*   M[1].Xp[5].x=0.398	;	M[1].Xp[5].y =3	; */
/*   M[1].Xp[6].x=0.314	;	M[1].Xp[6].y =3.6	; */
/*   M[1].Xp[7].x=0.286	;	M[1].Xp[7].y =4.2	; */
/*   M[1].Xp[8].x=0.226	;	M[1].Xp[8].y =4.8	; */
/*   M[1].Xp[9].x=0.208	;	M[1].Xp[9].y =5.4	; */
/*   M[1].Xp[10].x=0.191	;	M[1].Xp[10].y =6	; */
/*   M[1].Xp[11].x=0.176	;	M[1].Xp[11].y =6.6	; */
/*   M[1].Xp[12].x=0.171	;	M[1].Xp[12].y =7.2	; */
/*   M[1].Xp[13].x=0.167	;	M[1].Xp[13].y =7.8	; */
/*   M[1].Xp[14].x=0.169	;	M[1].Xp[14].y =8.4	; */
/*   M[1].Xp[15].x=0.151	;	M[1].Xp[15].y =9	; */
/*   M[1].Xp[16].x=0.166	;	M[1].Xp[16].y =9.6	; */
/*   M[1].Xp[17].x=0.151	;	M[1].Xp[17].y =10.2	; */
/*   M[1].Xp[18].x=0.141	;	M[1].Xp[18].y =10.8	; */
/*   M[1].Xp[19].x=0.128	;	M[1].Xp[19].y =11.4	; */
/*   M[1].Xp[20].x=0.12	;	M[1].Xp[20].y =12	; */
/*   M[1].Xp[21].x=0.131	;	M[1].Xp[21].y =12.6	; */
/*   M[1].Xp[22].x=0.15	;	M[1].Xp[22].y =13.2	; */
/*   M[1].Xp[23].x=0.123	;	M[1].Xp[23].y =13.8	; */
/*  M[1].Xp[24].x=0.133	;	M[1].Xp[24].y =14.4	; */
/*  M[1].Xp[25].x=0.142	;	M[1].Xp[25].y =15	; */
/*  M[1].Xp[26].x=0.113	;	M[1].Xp[26].y =15.6	; */
/*  M[1].Xp[27].x=0.118	;	M[1].Xp[27].y =16.2	; */
/*  M[1].Xp[28].x=0.123	;	M[1].Xp[28].y =16.8	; */
/*  M[1].Xp[29].x=0.119	;	M[1].Xp[29].y =17.4	; */
/*  M[1].Xp[30].x=0.116	;	M[1].Xp[30].y =18	; */
/*  M[1].Xp[31].x=0.128	;	M[1].Xp[31].y =18.6	; */
/*  M[1].Xp[32].x=0.118	;	M[1].Xp[32].y =19.2	; */
/*  M[1].Xp[33].x=0.125	;	M[1].Xp[33].y =19.8	; */
/*  M[1].Xp[34].x=0.091	;	M[1].Xp[34].y =20.4	; */
/*  M[1].Xp[35].x=0.084	;	M[1].Xp[35].y =21	; */
/*  M[1].Xp[36].x=0.067	;	M[1].Xp[36].y =21.6	; */
/*  M[1].Xp[37].x=0.038	;	M[1].Xp[37].y =22.2	; */
/*  M[1].Xp[38].x=0.027	;	M[1].Xp[38].y =22.8	; */
/*  M[1].Xp[39].x=0.006	;	M[1].Xp[39].y =23.4	; */
/*  M[1].Xp[40].x=0.005	;	M[1].Xp[40].y =24	; */
/*  M[1].Xp[41].x=0.001	;	M[1].Xp[41].y =24.6	; */
/*  M[1].Xp[42].x=0.002	;	M[1].Xp[42].y =25.2	; */
/*  M[1].Xp[43].x=0.001	;	M[1].Xp[43].y =25.8	; */
/*  M[1].pathlength(); */
/*  //station 3 */
/*  x[2] = 21.50; */
/*  M[2].allocate(44); */
/*  M[2].settype(SPLINE2D_LINEAR); */
/*  M[2].Xp[0].x = 	0.708	;	M[2].Xp[0].y =0	; */
/*  M[2].Xp[1].x = 	0.697	;	M[2].Xp[1].y =0.6	; */
/*  M[2].Xp[2].x = 	0.629	;	M[2].Xp[2].y =1.2	; */
/*  M[2].Xp[3].x = 	0.495	;	M[2].Xp[3].y =1.8	; */
/*  M[2].Xp[4].x = 	0.419	;	M[2].Xp[4].y =2.4	; */
/*  M[2].Xp[5].x = 	0.353	;	M[2].Xp[5].y =3	; */
/*  M[2].Xp[6].x = 	0.322	;	M[2].Xp[6].y =3.6	; */
/*  M[2].Xp[7].x = 	0.257	;	M[2].Xp[7].y =4.2	; */
/*  M[2].Xp[8].x = 	0.243	;	M[2].Xp[8].y =4.8	; */
/*  M[2].Xp[9].x = 	0.208	;	M[2].Xp[9].y = 5.4	; */
/*  M[2].Xp[10].x=        	0.208	;	M[2].Xp[10].y =6	; */
/*  M[2].Xp[11].x = 	0.164	;	M[2].Xp[11].y =6.6	; */
/*  M[2].Xp[12].x = 	0.189	;	M[2].Xp[12].y =7.2	; */
/*  M[2].Xp[13].x = 	0.166	;	M[2].Xp[13].y =7.8	; */
/*  M[2].Xp[14].x = 	0.156	;	M[2].Xp[14].y =8.4	; */
/*  M[2].Xp[15].x = 	0.139	;	M[2].Xp[15].y =9	; */
/*  M[2].Xp[16].x = 	0.14	;	M[2].Xp[16].y =9.6	; */
/*  M[2].Xp[17].x = 	0.134	;	M[2].Xp[17].y =10.2	; */
/*  M[2].Xp[18].x = 	0.125	;	M[2].Xp[18].y =10.8	; */
/*  M[2].Xp[19].x = 	0.145	;	M[2].Xp[19].y =11.4	; */
/*  M[2].Xp[20].x = 	0.145	;	M[2].Xp[20].y =12	; */
/*  M[2].Xp[21].x = 	0.144	;	M[2].Xp[21].y =12.6	; */
/*  M[2].Xp[22].x = 	0.149	;	M[2].Xp[22].y =13.2	; */
/*  M[2].Xp[23].x = 	0.129	;	M[2].Xp[23].y =13.8	; */
/*  M[2].Xp[24].x = 	0.147	;	M[2].Xp[24].y =14.4	; */
/*  M[2].Xp[25].x = 	0.131	;	M[2].Xp[25].y =15	; */
/*  M[2].Xp[26].x = 	0.141	;	M[2].Xp[26].y =15.6	; */
/*  M[2].Xp[27].x = 	0.129	;	M[2].Xp[27].y =16.2	; */
/*  M[2].Xp[28].x = 	0.142	;	M[2].Xp[28].y =16.8	; */
/*  M[2].Xp[29].x = 	0.129	;	M[2].Xp[29].y =17.4	; */
/*  M[2].Xp[30].x = 	0.118	;	M[2].Xp[30].y =18	; */
/*  M[2].Xp[31].x = 	0.106	;	M[2].Xp[31].y =18.6	; */
/*  M[2].Xp[32].x = 	0.112	;	M[2].Xp[32].y =19.2	; */
/*  M[2].Xp[33].x = 	0.11	;	M[2].Xp[33].y =19.8	; */
/*  M[2].Xp[34].x = 	0.103	;	M[2].Xp[34].y =20.4	; */
/*  M[2].Xp[35].x = 	0.087	;	M[2].Xp[35].y =21	; */
/*  M[2].Xp[36].x = 	0.058	;	M[2].Xp[36].y =21.6	; */
/*  M[2].Xp[37].x = 	0.036	;	M[2].Xp[37].y =22.2	; */
/*  M[2].Xp[38].x = 	0.024	;	M[2].Xp[38].y =22.8	; */
/*  M[2].Xp[39].x = 	0.01	;	M[2].Xp[39].y =23.4	; */
/*  M[2].Xp[40].x = 	0.004	;	M[2].Xp[40].y =24	; */
/*  M[2].Xp[41].x = 	0.001	;	M[2].Xp[41].y =24.6	; */
/*  M[2].Xp[42].x = 	0.002	;	M[2].Xp[42].y =25.2	; */
/*  M[2].Xp[43].x = 	0.002	;	M[2].Xp[43].y =25.8	; */
/*  M[2].pathlength(); */
/*  //station 4 */
/*  x[3] = 22.30; */
/*  M[3].allocate(44); */
/*  M[3].settype(SPLINE2D_LINEAR); */
/*  M[3].Xp[0].x = 	0.747	;	M[3].Xp[0].y =0	; */
/*  M[3].Xp[1].x = 	0.699	;	M[3].Xp[1].y =0.6	; */
/*  M[3].Xp[2].x = 	0.584	;	M[3].Xp[2].y =1.2	; */
/*  M[3].Xp[3].x = 	0.511	;	M[3].Xp[3].y =1.8	; */
/*  M[3].Xp[4].x = 	0.452	;	M[3].Xp[4].y =2.4	; */
/*  M[3].Xp[5].x = 	0.405	;	M[3].Xp[5].y =3	; */
/*  M[3].Xp[6].x = 	0.323	;	M[3].Xp[6].y =3.6	; */
/*  M[3].Xp[7].x = 	0.273	;	M[3].Xp[7].y =4.2	; */
/*  M[3].Xp[8].x = 	0.212	;	M[3].Xp[8].y =4.8	; */
/*  M[3].Xp[9].x = 	0.188	;	M[3].Xp[9].y =5.4	; */
/*  M[3].Xp[10].x = 	0.186	;	M[3].Xp[10].y =6	; */
/*  M[3].Xp[11].x = 	0.192	;	M[3].Xp[11].y =6.6	; */
/*  M[3].Xp[12].x = 	0.176	;	M[3].Xp[12].y =7.2	; */
/*  M[3].Xp[13].x = 	0.17	;	M[3].Xp[13].y =7.8	; */
/*  M[3].Xp[14].x = 	0.167	;	M[3].Xp[14].y =8.4	; */
/*  M[3].Xp[15].x = 	0.141	;	M[3].Xp[15].y =9	; */
/*  M[3].Xp[16].x = 	0.147	;	M[3].Xp[16].y =9.6	; */
/*  M[3].Xp[17].x = 	0.147	;	M[3].Xp[17].y =10.2	; */
/*  M[3].Xp[18].x = 	0.135	;	M[3].Xp[18].y =10.8	; */
/*  M[3].Xp[19].x = 	0.14	;	M[3].Xp[19].y =11.4	; */
/*  M[3].Xp[20].x = 	0.136	;	M[3].Xp[20].y =12	; */
/*  M[3].Xp[21].x = 	0.14	;	M[3].Xp[21].y =12.6	; */
/*  M[3].Xp[22].x = 	0.135	;	M[3].Xp[22].y =13.2	; */
/*  M[3].Xp[23].x = 	0.12	;	M[3].Xp[23].y =13.8	; */
/*  M[3].Xp[24].x = 	0.128	;	M[3].Xp[24].y =14.4	; */
/*  M[3].Xp[25].x = 	0.122	;	M[3].Xp[25].y =15	; */
/*  M[3].Xp[26].x = 	0.13	;	M[3].Xp[26].y =15.6	; */
/*  M[3].Xp[27].x = 	0.118	;	M[3].Xp[27].y =16.2	; */
/*  M[3].Xp[28].x = 	0.127	;	M[3].Xp[28].y =16.8	; */
/*  M[3].Xp[29].x = 	0.116	;	M[3].Xp[29].y =17.4	; */
/*  M[3].Xp[30].x = 	0.11	;	M[3].Xp[30].y =18	; */
/*  M[3].Xp[31].x = 	0.116	;	M[3].Xp[31].y =18.6	; */
/*  M[3].Xp[32].x = 	0.108	;	M[3].Xp[32].y =19.2	; */
/*  M[3].Xp[33].x = 	0.105	;	M[3].Xp[33].y =19.8	; */
/*  M[3].Xp[34].x = 	0.093	;	M[3].Xp[34].y =20.4	; */
/*  M[3].Xp[35].x = 	0.078	;	M[3].Xp[35].y =21	; */
/*  M[3].Xp[36].x = 	0.058	;	M[3].Xp[36].y =21.6	; */
/*  M[3].Xp[37].x = 	0.031	;	M[3].Xp[37].y =22.2	; */
/*  M[3].Xp[38].x = 	0.022	;	M[3].Xp[38].y =22.8	; */
/*  M[3].Xp[39].x = 	0.014	;	M[3].Xp[39].y =23.4	; */
/*  M[3].Xp[40].x = 	0.004	;	M[3].Xp[40].y =24	; */
/*  M[3].Xp[41].x = 	0.002	;	M[3].Xp[41].y =24.6	; */
/*  M[3].Xp[42].x = 	0.004	;	M[3].Xp[42].y =25.2	; */
/*  M[3].Xp[43].x = 	0.006	;	M[3].Xp[43].y =25.8	; */
/*  M[3].pathlength(); */
 
/*  //station 5 */
/*  x[4] = 23.00; */
/*  M[4].allocate(44); */
/*  M[4].settype(SPLINE2D_LINEAR); */
/*  M[4].Xp[0].x = 	0.733	;	M[4].Xp[0].y =0	; */
/*  M[4].Xp[1].x = 	0.672	;	M[4].Xp[1].y =0.6	; */
/*  M[4].Xp[2].x = 	0.58	;	M[4].Xp[2].y =1.2	; */
/*  M[4].Xp[3].x = 	0.555	;	M[4].Xp[3].y =1.8	; */
/*  M[4].Xp[4].x = 	0.462	;	M[4].Xp[4].y =2.4	; */
/*  M[4].Xp[5].x = 	0.402	;	M[4].Xp[5].y =3	; */
/*  M[4].Xp[6].x = 	0.321	;	M[4].Xp[6].y =3.6	; */
/*  M[4].Xp[7].x = 	0.274	;	M[4].Xp[7].y =4.2	; */
/*  M[4].Xp[8].x = 	0.235	;	M[4].Xp[8].y =4.8	; */
/*  M[4].Xp[9].x = 	0.2	;	M[4].Xp[9].y =5.4	; */
/*  M[4].Xp[10].x = 	0.188	;	M[4].Xp[10].y =6	; */
/*  M[4].Xp[11].x = 	0.196	;	M[4].Xp[11].y =6.6	; */
/*  M[4].Xp[12].x = 	0.185	;	M[4].Xp[12].y =7.2	; */
/*  M[4].Xp[13].x = 	0.164	;	M[4].Xp[13].y =7.8	; */
/*  M[4].Xp[14].x = 	0.164	;	M[4].Xp[14].y =8.4	; */
/*  M[4].Xp[15].x = 	0.152	;	M[4].Xp[15].y =9	; */
/*  M[4].Xp[16].x = 	0.156	;	M[4].Xp[16].y =9.6	; */
/*  M[4].Xp[17].x = 	0.145	;	M[4].Xp[17].y =10.2	; */
/*  M[4].Xp[18].x = 	0.144	;	M[4].Xp[18].y =10.8	; */
/*  M[4].Xp[19].x = 	0.143	;	M[4].Xp[19].y =11.4	; */
/*  M[4].Xp[20].x = 	0.137	;	M[4].Xp[20].y =12	; */
/*  M[4].Xp[21].x = 	0.146	;	M[4].Xp[21].y =12.6	; */
/*  M[4].Xp[22].x = 	0.133	;	M[4].Xp[22].y =13.2	; */
/*  M[4].Xp[23].x = 	0.137	;	M[4].Xp[23].y =13.8	; */
/*  M[4].Xp[24].x = 	0.138	;	M[4].Xp[24].y =14.4	; */
/*  M[4].Xp[25].x = 	0.127	;	M[4].Xp[25].y =15	; */
/*  M[4].Xp[26].x = 	0.127	;	M[4].Xp[26].y =15.6	; */
/*  M[4].Xp[27].x = 	0.124	;	M[4].Xp[27].y =16.2	; */
/*  M[4].Xp[28].x = 	0.124	;	M[4].Xp[28].y =16.8	; */
/*  M[4].Xp[29].x = 	0.121	;	M[4].Xp[29].y =17.4	; */
/*  M[4].Xp[30].x = 	0.123	;	M[4].Xp[30].y =18	; */
/*  M[4].Xp[31].x = 	0.119	;	M[4].Xp[31].y =18.6	; */
/*  M[4].Xp[32].x = 	0.113	;	M[4].Xp[32].y =19.2	; */
/*  M[4].Xp[33].x = 	0.099	;	M[4].Xp[33].y =19.8	; */
/*  M[4].Xp[34].x = 	0.087	;	M[4].Xp[34].y =20.4	; */
/*  M[4].Xp[35].x = 	0.08	;	M[4].Xp[35].y =21	; */
/*  M[4].Xp[36].x = 	0.042	;	M[4].Xp[36].y =21.6	; */
/*  M[4].Xp[37].x = 	0.036	;	M[4].Xp[37].y =22.2	; */
/*  M[4].Xp[38].x = 	0.019	;	M[4].Xp[38].y =22.8	; */
/*  M[4].Xp[39].x = 	0.016	;	M[4].Xp[39].y =23.4	; */
/*  M[4].Xp[40].x = 	0.005	;	M[4].Xp[40].y =24	; */
/*  M[4].Xp[41].x = 	0.002	;	M[4].Xp[41].y =24.6	; */
/*  M[4].Xp[42].x = 	0.001	;	M[4].Xp[42].y =25.2	; */
/*  M[4].Xp[43].x = 	0	;	M[4].Xp[43].y =25.8	; */
/*  M[4].pathlength(); */
/*  //station 6 */
/*  x[5] = 23.80; */
/*  M[5].allocate(44); */
/*  M[5].settype(SPLINE2D_LINEAR); */
/*  M[5].Xp[0].x = 	0.696	;	M[5].Xp[0].y =0	; */
/*  M[5].Xp[1].x = 	0.687	;	M[5].Xp[1].y =0.6	; */
/*  M[5].Xp[2].x = 	0.568	;	M[5].Xp[2].y =1.2	; */
/*  M[5].Xp[3].x = 	0.554	;	M[5].Xp[3].y =1.8	; */
/*  M[5].Xp[4].x = 	0.467	;	M[5].Xp[4].y =2.4	; */
/*  M[5].Xp[5].x = 	0.387	;	M[5].Xp[5].y =3	; */
/*  M[5].Xp[6].x = 	0.321	;	M[5].Xp[6].y =3.6	; */
/*  M[5].Xp[7].x = 	0.272	;	M[5].Xp[7].y =4.2	; */
/*  M[5].Xp[8].x = 	0.248	;	M[5].Xp[8].y =4.8	; */
/*  M[5].Xp[9].x = 	0.217	;	M[5].Xp[9].y =5.4	; */
/*  M[5].Xp[10].x = 	0.189	;	M[5].Xp[10].y =6	; */
/*  M[5].Xp[11].x = 	0.203	;	M[5].Xp[11].y =6.6	; */
/*  M[5].Xp[12].x = 	0.175	;	M[5].Xp[12].y =7.2	; */
/*  M[5].Xp[13].x = 	0.172	;	M[5].Xp[13].y =7.8	; */
/*  M[5].Xp[14].x = 	0.166	;	M[5].Xp[14].y =8.4	; */
/*  M[5].Xp[15].x = 	0.151	;	M[5].Xp[15].y =9	; */
/*  M[5].Xp[16].x = 	0.157	;	M[5].Xp[16].y =9.6	; */
/*  M[5].Xp[17].x = 	0.14	;	M[5].Xp[17].y =10.2	; */
/*  M[5].Xp[18].x = 	0.14	;	M[5].Xp[18].y =10.8	; */
/*   M[5].Xp[19].x = 	0.145	;	 M[5].Xp[19].y =11.4	; */
/*  M[5].Xp[20].x = 	0.154	;	 M[5].Xp[20].y =12	; */
/*  M[5].Xp[21].x = 	0.142	;	 M[5].Xp[21].y =12.6	; */
/*  M[5].Xp[22].x = 	0.133	;	 M[5].Xp[22].y =13.2	; */
/*  M[5].Xp[23].x = 	0.136	;	 M[5].Xp[23].y =13.8	; */
/*  M[5].Xp[24].x = 	0.141	;	 M[5].Xp[24].y =14.4	; */
/*  M[5].Xp[25].x = 	0.12	;	 M[5].Xp[25].y =15	; */
/*  M[5].Xp[26].x = 	0.118	;	 M[5].Xp[26].y =15.6	; */
/*  M[5].Xp[27].x = 	0.115	;	 M[5].Xp[27].y =16.2	; */
/*  M[5].Xp[28].x = 	0.116	;	 M[5].Xp[28].y =16.8	; */
/*  M[5].Xp[29].x = 	0.124	;	 M[5].Xp[29].y =17.4	; */
/*  M[5].Xp[30].x = 	0.118	;	 M[5].Xp[30].y =18	; */
/*  M[5].Xp[31].x = 	0.114	;	 M[5].Xp[31].y =18.6	; */
/*  M[5].Xp[32].x = 	0.103	;	 M[5].Xp[32].y =19.2	; */
/*  M[5].Xp[33].x = 	0.1	;	 M[5].Xp[33].y =19.8	; */
/*  M[5].Xp[34].x = 	0.093	;	 M[5].Xp[34].y =20.4	; */
/*  M[5].Xp[35].x = 	0.077	;	 M[5].Xp[35].y =21	; */
/*  M[5].Xp[36].x = 	0.051	;	 M[5].Xp[36].y =21.6	; */
/*  M[5].Xp[37].x = 	0.032	;	 M[5].Xp[37].y =22.2	; */
/*  M[5].Xp[38].x = 	0.014	;	 M[5].Xp[38].y =22.8	; */
/*  M[5].Xp[39].x = 	0.01	;	 M[5].Xp[39].y =23.4	; */
/*  M[5].Xp[40].x = 	0.001	;	 M[5].Xp[40].y =24	; */
/*  M[5].Xp[41].x = 	0.002	;	 M[5].Xp[41].y =24.6	; */
/*  M[5].Xp[42].x = 	0.001	;	 M[5].Xp[42].y =25.2	; */
/*  M[5].Xp[43].x = 	0.001	;	 M[5].Xp[43].y =25.8	; */
/*  M[5].pathlength(); */
/* //station 7 */
/*  x[6] = 24.50; */
/*  M[6].allocate(44); */
/*  M[6].settype(SPLINE2D_LINEAR); */
/*  M[6].Xp[0].x = 	0.731	;	 M[6].Xp[0].y =0	; */
/*  M[6].Xp[1].x = 	0.691	;	 M[6].Xp[1].y =0.6	; */
/*  M[6].Xp[2].x = 	0.584	;	 M[6].Xp[2].y =1.2	; */
/*  M[6].Xp[3].x = 	0.533	;	 M[6].Xp[3].y =1.8	; */
/*  M[6].Xp[4].x = 	0.47	;	 M[6].Xp[4].y =2.4	; */
/*  M[6].Xp[5].x = 	0.366	;	 M[6].Xp[5].y =3	; */
/*  M[6].Xp[6].x = 	0.317	;	 M[6].Xp[6].y =3.6	; */
/*  M[6].Xp[7].x = 	0.302	;	 M[6].Xp[7].y =4.2	; */
/*  M[6].Xp[8].x = 	0.25	;	 M[6].Xp[8].y =4.8	; */
/*  M[6].Xp[9].x = 	0.222	;	 M[6].Xp[9].y =5.4	; */
/*  M[6].Xp[10].x = 	0.203	;	 M[6].Xp[10].y =6	; */
/*  M[6].Xp[11].x = 	0.187	;	 M[6].Xp[11].y =6.6	; */
/*  M[6].Xp[12].x = 	0.176	;	 M[6].Xp[12].y =7.2	; */
/*  M[6].Xp[13].x = 	0.174	;	 M[6].Xp[13].y =7.8	; */
/*  M[6].Xp[14].x = 	0.163	;	 M[6].Xp[14].y =8.4	; */
/*  M[6].Xp[15].x = 	0.165	;	 M[6].Xp[15].y =9	; */
/*  M[6].Xp[16].x = 	0.146	;	 M[6].Xp[16].y =9.6	; */
/*  M[6].Xp[17].x = 	0.153	;	 M[6].Xp[17].y =10.2	; */
/*  M[6].Xp[18].x = 	0.14	;	 M[6].Xp[18].y =10.8	; */
/*  M[6].Xp[19].x = 	0.141	;	 M[6].Xp[19].y =11.4	; */
/*  M[6].Xp[20].x = 	0.135	;	 M[6].Xp[20].y =12	; */
/*  M[6].Xp[21].x = 	0.144	;	 M[6].Xp[21].y =12.6	; */
/*  M[6].Xp[22].x = 	0.142	;	 M[6].Xp[22].y =13.2	; */
/*  M[6].Xp[23].x = 	0.125	;	 M[6].Xp[23].y =13.8	; */
/*  M[6].Xp[24].x = 	0.136	;	 M[6].Xp[24].y =14.4	; */
/*  M[6].Xp[25].x = 	0.129	;	 M[6].Xp[25].y =15	; */
/*  M[6].Xp[26].x = 	0.124	;	 M[6].Xp[26].y =15.6	; */
/*  M[6].Xp[27].x = 	0.13	;	 M[6].Xp[27].y =16.2	; */
/*  M[6].Xp[28].x = 	0.116	;	 M[6].Xp[28].y =16.8	; */
/*  M[6].Xp[29].x = 	0.124	;	 M[6].Xp[29].y =17.4	; */
/*  M[6].Xp[30].x = 	0.122	;	 M[6].Xp[30].y =18	; */
/*  M[6].Xp[31].x = 	0.118	;	 M[6].Xp[31].y =18.6	; */
/*  M[6].Xp[32].x = 	0.098	;	 M[6].Xp[32].y =19.2	; */
/*  M[6].Xp[33].x = 	0.097	;	 M[6].Xp[33].y =19.8	; */
/*  M[6].Xp[34].x = 	0.081	;	 M[6].Xp[34].y =20.4	; */
/*  M[6].Xp[35].x = 	0.056	;	 M[6].Xp[35].y =21	; */
/*  M[6].Xp[36].x = 	0.039	;	 M[6].Xp[36].y =21.6	; */
/*  M[6].Xp[37].x = 	0.024	;	 M[6].Xp[37].y =22.2	; */
/*  M[6].Xp[38].x = 	0.015	;	 M[6].Xp[38].y =22.8	; */
/*  M[6].Xp[39].x = 	0.008	;	 M[6].Xp[39].y =23.4	; */
/*  M[6].Xp[40].x = 	0.003	;	 M[6].Xp[40].y =24	; */
/*  M[6].Xp[41].x = 	0	;	 M[6].Xp[41].y =24.6	; */
/*  M[6].Xp[42].x = 	0.001	;	 M[6].Xp[42].y =25.2	; */
/*  M[6].Xp[43].x = 	0	;	 M[6].Xp[43].y =25.8	; */
/*  M[6].pathlength(); */

/* //station 8 */
/*  x[7] = 25.30; */
/*  M[7].allocate(44); */
/*  M[7].settype(SPLINE2D_LINEAR); */
/*  M[7].Xp[0].x = 	0.724	;	 M[7].Xp[0].y =0	; */
/*  M[7].Xp[1].x = 	0.664	;	 M[7].Xp[1].y =0.6	; */
/*  M[7].Xp[2].x = 	0.641	;	 M[7].Xp[2].y =1.2	; */
/*  M[7].Xp[3].x = 	0.575	;	 M[7].Xp[3].y =1.8	; */
/*  M[7].Xp[4].x = 	0.473	;	 M[7].Xp[4].y =2.4	; */
/*  M[7].Xp[5].x = 	0.385	;	 M[7].Xp[5].y =3	; */
/*  M[7].Xp[6].x = 	0.327	;	 M[7].Xp[6].y =3.6	; */
/*  M[7].Xp[7].x = 	0.284	;	 M[7].Xp[7].y =4.2	; */
/*  M[7].Xp[8].x = 	0.266	;	 M[7].Xp[8].y =4.8	; */
/*  M[7].Xp[9].x = 	0.245	;	 M[7].Xp[9].y =5.4	; */
/*  M[7].Xp[10].x = 	0.223	;	 M[7].Xp[10].y =6	; */
/*  M[7].Xp[11].x = 	0.188	;	 M[7].Xp[11].y =6.6	; */
/*  M[7].Xp[12].x = 	0.184	;	 M[7].Xp[12].y =7.2	; */
/*  M[7].Xp[13].x = 	0.172	;	 M[7].Xp[13].y =7.8	; */
/*  M[7].Xp[14].x = 	0.165	;	 M[7].Xp[14].y =8.4	; */
/*  M[7].Xp[15].x = 	0.18	;	 M[7].Xp[15].y =9	; */
/*  M[7].Xp[16].x = 	0.157	;	 M[7].Xp[16].y =9.6	; */
/*  M[7].Xp[17].x = 	0.161	;	 M[7].Xp[17].y =10.2	; */
/*  M[7].Xp[18].x = 	0.149	;	 M[7].Xp[18].y =10.8	; */
/*  M[7].Xp[19].x = 	0.149	;	 M[7].Xp[19].y =11.4	; */
/*  M[7].Xp[20].x = 	0.147	;	 M[7].Xp[20].y =12	; */
/*  M[7].Xp[21].x = 	0.134	;	 M[7].Xp[21].y =12.6	; */
/*  M[7].Xp[22].x = 	0.142	;	 M[7].Xp[22].y =13.2	; */
/*  M[7].Xp[23].x = 	0.138	;	 M[7].Xp[23].y =13.8	; */
/*  M[7].Xp[24].x = 	0.132	;	 M[7].Xp[24].y =14.4	; */
/*  M[7].Xp[25].x = 	0.123	;	 M[7].Xp[25].y =15	; */
/*  M[7].Xp[26].x = 	0.127	;	 M[7].Xp[26].y =15.6	; */
/*  M[7].Xp[27].x = 	0.133	;	 M[7].Xp[27].y =16.2	; */
/*  M[7].Xp[28].x = 	0.128	;	 M[7].Xp[28].y =16.8	; */
/*  M[7].Xp[29].x = 	0.127	;	 M[7].Xp[29].y =17.4	; */
/*  M[7].Xp[30].x = 	0.117	;	 M[7].Xp[30].y =18	; */
/*  M[7].Xp[31].x = 	0.111	;	 M[7].Xp[31].y =18.6	; */
/*  M[7].Xp[32].x = 	0.11	;	 M[7].Xp[32].y =19.2	; */
/*  M[7].Xp[33].x = 	0.085	;	 M[7].Xp[33].y =19.8	; */
/*  M[7].Xp[34].x = 	0.072	;	 M[7].Xp[34].y =20.4	; */
/*  M[7].Xp[35].x = 	0.058	;	 M[7].Xp[35].y =21	; */
/*  M[7].Xp[36].x = 	0.042	;	 M[7].Xp[36].y =21.6	; */
/*  M[7].Xp[37].x = 	0.029	;	 M[7].Xp[37].y =22.2	; */
/*  M[7].Xp[38].x = 	0.019	;	 M[7].Xp[38].y =22.8	; */
/*  M[7].Xp[39].x = 	0.009	;	 M[7].Xp[39].y =23.4	; */
/*  M[7].Xp[40].x = 	0.007	;	 M[7].Xp[40].y =24	; */
/*  M[7].Xp[41].x = 	0.001	;	 M[7].Xp[41].y =24.6	; */
/*  M[7].Xp[42].x = 	0.004	;	 M[7].Xp[42].y =25.2	; */
/*  M[7].Xp[43].x = 	0.002	;	 M[7].Xp[43].y =25.8	; */

/*  M[7].pathlength(); */
/* //station 9 */
/*  x[8] = 26.00; */
/*  M[8].allocate(44); */
/*  M[8].settype(SPLINE2D_LINEAR); */
/*  M[8].Xp[0].x = 	0.68	;	 M[8].Xp[0].y =0	; */
/*  M[8].Xp[1].x = 	0.652	;	 M[8].Xp[1].y =0.6	; */
/*  M[8].Xp[2].x = 	0.607	;	 M[8].Xp[2].y =1.2	; */
/*  M[8].Xp[3].x = 	0.561	;	 M[8].Xp[3].y =1.8	; */
/*  M[8].Xp[4].x = 	0.444	;	 M[8].Xp[4].y =2.4	; */
/*  M[8].Xp[5].x = 	0.371	;	 M[8].Xp[5].y =3	; */
/*  M[8].Xp[6].x = 	0.334	;	 M[8].Xp[6].y =3.6	; */
/*  M[8].Xp[7].x = 	0.286	;	 M[8].Xp[7].y =4.2	; */
/*  M[8].Xp[8].x = 	0.246	;	 M[8].Xp[8].y =4.8	; */
/*  M[8].Xp[9].x = 	0.224	;	 M[8].Xp[9].y =5.4	; */
/*  M[8].Xp[10].x = 	0.204	;	 M[8].Xp[10].y =6	; */
/*  M[8].Xp[11].x = 	0.187	;	 M[8].Xp[11].y =6.6	; */
/*  M[8].Xp[12].x = 	0.185	;	 M[8].Xp[12].y =7.2	; */
/*  M[8].Xp[13].x = 	0.192	;	 M[8].Xp[13].y =7.8	; */
/*  M[8].Xp[14].x = 	0.16	;	 M[8].Xp[14].y =8.4	; */
/*  M[8].Xp[15].x = 	0.158	;	 M[8].Xp[15].y =9	; */
/*  M[8].Xp[16].x = 	0.135	;	 M[8].Xp[16].y =9.6	; */
/*  M[8].Xp[17].x = 	0.141	;	 M[8].Xp[17].y =10.2	; */
/*  M[8].Xp[18].x = 	0.143	;	 M[8].Xp[18].y =10.8	; */
/*  M[8].Xp[19].x = 	0.153	;	 M[8].Xp[19].y =11.4	; */
/*  M[8].Xp[20].x = 	0.138	;	 M[8].Xp[20].y =12	; */
/*  M[8].Xp[21].x = 	0.128	;	 M[8].Xp[21].y =12.6	; */
/*  M[8].Xp[22].x = 	0.132	;	 M[8].Xp[22].y =13.2	; */
/*  M[8].Xp[23].x = 	0.123	;	 M[8].Xp[23].y =13.8	; */
/*  M[8].Xp[24].x = 	0.128	;	 M[8].Xp[24].y =14.4	; */
/*  M[8].Xp[25].x = 	0.135	;	 M[8].Xp[25].y =15	; */
/*  M[8].Xp[26].x = 	0.124	;	 M[8].Xp[26].y =15.6	; */
/*  M[8].Xp[27].x = 	0.125	;	 M[8].Xp[27].y =16.2	; */
/*  M[8].Xp[28].x = 	0.121	;	 M[8].Xp[28].y =16.8	; */
/*  M[8].Xp[29].x = 	0.115	;	 M[8].Xp[29].y =17.4	; */
/*  M[8].Xp[30].x = 	0.11	;	 M[8].Xp[30].y =18	; */
/*  M[8].Xp[31].x = 	0.106	;	 M[8].Xp[31].y =18.6	; */
/*  M[8].Xp[32].x = 	0.097	;	 M[8].Xp[32].y =19.2	; */
/*  M[8].Xp[33].x = 	0.085	;	 M[8].Xp[33].y =19.8	; */
/*  M[8].Xp[34].x = 	0.075	;	 M[8].Xp[34].y =20.4	; */
/*  M[8].Xp[35].x = 	0.058	;	 M[8].Xp[35].y =21	; */
/*  M[8].Xp[36].x = 	0.034	;	 M[8].Xp[36].y =21.6	; */
/*  M[8].Xp[37].x = 	0.025	;	 M[8].Xp[37].y =22.2	; */
/*  M[8].Xp[38].x = 	0.016	;	 M[8].Xp[38].y =22.8	; */
/*  M[8].Xp[39].x = 	0.012	;	 M[8].Xp[39].y =23.4	; */
/*  M[8].Xp[40].x = 	0.003	;	 M[8].Xp[40].y =24	; */
/*  M[8].Xp[41].x = 	0.002	;	 M[8].Xp[41].y =24.6	; */
/*  M[8].Xp[42].x = 	0.001	;	 M[8].Xp[42].y =25.2	; */
/*  M[8].Xp[43].x = 	0.001	;	 M[8].Xp[43].y =25.8	; */
/*  M[8].pathlength(); */
/* //station 10 */
/*  x[9] = 26.800; */
/*  M[9].allocate(44); */
/*  M[9].settype(SPLINE2D_LINEAR); */
/*  M[9].Xp[0].x = 	0.652	;	 M[9].Xp[0].y =0	; */
/*  M[9].Xp[1].x = 	0.631	;	 M[9].Xp[1].y =0.6	; */
/*  M[9].Xp[2].x = 	0.631	;	 M[9].Xp[2].y =1.2	; */
/*  M[9].Xp[3].x = 	0.58	;	 M[9].Xp[3].y =1.8	; */
/*  M[9].Xp[4].x = 	0.481	;	 M[9].Xp[4].y =2.4	; */
/*  M[9].Xp[5].x = 	0.377	;	 M[9].Xp[5].y =3	; */
/*  M[9].Xp[6].x = 	0.326	;	 M[9].Xp[6].y =3.6	; */
/*  M[9].Xp[7].x = 	0.283	;	 M[9].Xp[7].y =4.2	; */
/*  M[9].Xp[8].x = 	0.256	;	 M[9].Xp[8].y =4.8	; */
/*  M[9].Xp[9].x = 	0.225	;	 M[9].Xp[9].y =5.4	; */
/*  M[9].Xp[10].x = 	0.196	;	 M[9].Xp[10].y =6	; */
/*  M[9].Xp[11].x = 	0.193	;	 M[9].Xp[11].y =6.6	; */
/*  M[9].Xp[12].x = 	0.177	;	 M[9].Xp[12].y =7.2	; */
/*  M[9].Xp[13].x = 	0.179	;	 M[9].Xp[13].y =7.8	; */
/*  M[9].Xp[14].x = 	0.173	;	 M[9].Xp[14].y =8.4	; */
/*  M[9].Xp[15].x = 	0.162	;	 M[9].Xp[15].y =9	; */
/*  M[9].Xp[16].x = 	0.146	;	 M[9].Xp[16].y =9.6	; */
/*  M[9].Xp[17].x = 	0.144	;	 M[9].Xp[17].y =10.2	; */
/*  M[9].Xp[18].x = 	0.138	;	 M[9].Xp[18].y =10.8	; */
/*  M[9].Xp[19].x = 	0.14	;	 M[9].Xp[19].y =11.4	; */
/*  M[9].Xp[20].x = 	0.136	;	 M[9].Xp[20].y =12	; */
/*  M[9].Xp[21].x = 	0.118	;	 M[9].Xp[21].y =12.6	; */
/*  M[9].Xp[22].x = 	0.138	;	 M[9].Xp[22].y =13.2	; */
/*  M[9].Xp[23].x = 	0.124	;	 M[9].Xp[23].y =13.8	; */
/*  M[9].Xp[24].x = 	0.126	;	 M[9].Xp[24].y =14.4	; */
/*  M[9].Xp[25].x = 	0.121	;	 M[9].Xp[25].y =15	; */
/*  M[9].Xp[26].x = 	0.125	;	 M[9].Xp[26].y =15.6	; */
/*  M[9].Xp[27].x = 	0.119	;	 M[9].Xp[27].y =16.2	; */
/*  M[9].Xp[28].x = 	0.112	;	 M[9].Xp[28].y =16.8	; */
/*  M[9].Xp[29].x = 	0.116	;	 M[9].Xp[29].y =17.4	; */
/*  M[9].Xp[30].x = 	0.113	;	 M[9].Xp[30].y =18	; */
/*  M[9].Xp[31].x = 	0.096	;	 M[9].Xp[31].y =18.6	; */
/*  M[9].Xp[32].x = 	0.09	;	 M[9].Xp[32].y =19.2	; */
/*  M[9].Xp[33].x = 	0.087	;	 M[9].Xp[33].y =19.8	; */
/*  M[9].Xp[34].x = 	0.069	;	 M[9].Xp[34].y =20.4	; */
/*  M[9].Xp[35].x = 	0.048	;	 M[9].Xp[35].y =21	; */
/*  M[9].Xp[36].x = 	0.034	;	 M[9].Xp[36].y =21.6	; */
/*  M[9].Xp[37].x = 	0.023	;	 M[9].Xp[37].y =22.2	; */
/*  M[9].Xp[38].x = 	0.014	;	 M[9].Xp[38].y =22.8	; */
/*  M[9].Xp[39].x = 	0.009	;	 M[9].Xp[39].y =23.4	; */
/*  M[9].Xp[40].x = 	0.003	;	 M[9].Xp[40].y =24	; */
/*  M[9].Xp[41].x = 	0.002	;	 M[9].Xp[41].y =24.6	; */
/*  M[9].Xp[42].x = 	0.002	;	 M[9].Xp[42].y =25.2	; */
/*  M[9].Xp[43].x = 	0.002	;	 M[9].Xp[43].y =25.8	; */
/*  M[9].pathlength(); */
/* //station 11 */
/*  x[10] = 27.50; */
/*  M[10].allocate(44); */
/*  M[10].settype(SPLINE2D_LINEAR); */
/*  M[10].Xp[0].x = 	0.699	;	 M[10].Xp[0].y =0	; */
/*  M[10].Xp[1].x = 	0.683	;	 M[10].Xp[1].y =0.6	; */
/*  M[10].Xp[2].x = 	0.62	;	 M[10].Xp[2].y =1.2	; */
/*  M[10].Xp[3].x = 	0.545	;	 M[10].Xp[3].y =1.8	; */
/*  M[10].Xp[4].x = 	0.445	;	 M[10].Xp[4].y =2.4	; */
/*  M[10].Xp[5].x = 	0.381	;	 M[10].Xp[5].y =3	; */
/*  M[10].Xp[6].x = 	0.326	;	 M[10].Xp[6].y =3.6	; */
/*  M[10].Xp[7].x = 	0.283	;	 M[10].Xp[7].y =4.2	; */
/*  M[10].Xp[8].x = 	0.248	;	 M[10].Xp[8].y =4.8	; */
/*  M[10].Xp[9].x = 	0.222	;	 M[10].Xp[9].y =5.4	; */
/*  M[10].Xp[10].x = 	0.222	;	 M[10].Xp[10].y =6	; */
/*  M[10].Xp[11].x = 	0.208	;	 M[10].Xp[11].y =6.6	; */
/*  M[10].Xp[12].x = 	0.191	;	 M[10].Xp[12].y =7.2	; */
/*  M[10].Xp[13].x = 	0.189	;	 M[10].Xp[13].y =7.8	; */
/*  M[10].Xp[14].x = 	0.177	;	 M[10].Xp[14].y =8.4	; */
/*  M[10].Xp[15].x = 	0.165	;	 M[10].Xp[15].y =9	; */
/*  M[10].Xp[16].x = 	0.155	;	 M[10].Xp[16].y =9.6	; */
/*  M[10].Xp[17].x = 	0.15	;	 M[10].Xp[17].y =10.2	; */
/*  M[10].Xp[18].x = 	0.142	;	 M[10].Xp[18].y =10.8	; */
/*  M[10].Xp[19].x = 	0.151	;	 M[10].Xp[19].y =11.4	; */
/*  M[10].Xp[20].x = 	0.126	;	 M[10].Xp[20].y =12	; */
/*  M[10].Xp[21].x = 	0.123	;	 M[10].Xp[21].y =12.6	; */
/*  M[10].Xp[22].x = 	0.123	;	 M[10].Xp[22].y =13.2	; */
/*  M[10].Xp[23].x = 	0.132	;	 M[10].Xp[23].y =13.8	; */
/*  M[10].Xp[24].x = 	0.132	;	 M[10].Xp[24].y =14.4	; */
/*  M[10].Xp[25].x = 	0.117	;	 M[10].Xp[25].y =15	; */
/*  M[10].Xp[26].x = 	0.118	;	 M[10].Xp[26].y =15.6	; */
/*  M[10].Xp[27].x = 	0.109	;	 M[10].Xp[27].y =16.2	; */
/*  M[10].Xp[28].x = 	0.114	;	 M[10].Xp[28].y =16.8	; */
/*  M[10].Xp[29].x = 	0.113	;	 M[10].Xp[29].y =17.4	; */
/*  M[10].Xp[30].x = 	0.106	;	 M[10].Xp[30].y =18	; */
/*  M[10].Xp[31].x = 	0.098	;	 M[10].Xp[31].y =18.6	; */
/*  M[10].Xp[32].x = 	0.101	;	 M[10].Xp[32].y =19.2	; */
/*  M[10].Xp[33].x = 	0.077	;	 M[10].Xp[33].y =19.8	; */
/*  M[10].Xp[34].x = 	0.068	;	 M[10].Xp[34].y =20.4	; */
/*  M[10].Xp[35].x = 	0.055	;	 M[10].Xp[35].y =21	; */
/*  M[10].Xp[36].x = 	0.043	;	 M[10].Xp[36].y =21.6	; */
/*  M[10].Xp[37].x = 	0.026	;	 M[10].Xp[37].y =22.2	; */
/*  M[10].Xp[38].x = 	0.012	;	 M[10].Xp[38].y =22.8	; */
/*  M[10].Xp[39].x = 	0.01	;	 M[10].Xp[39].y =23.4	; */
/*  M[10].Xp[40].x = 	0.005	;	 M[10].Xp[40].y =24	; */
/*  M[10].Xp[41].x = 	0.001	;	 M[10].Xp[41].y =24.6	; */
/*  M[10].Xp[42].x = 	0.003	;	 M[10].Xp[42].y =25.2	; */
/*  M[10].Xp[43].x = 	0.002	;	 M[10].Xp[43].y =25.8	; */

/*  M[10].pathlength(); */


/* //station 12 */
/*  x[11] = 28.30; */
/*  M[11].allocate(44); */
/*  M[11].settype(SPLINE2D_LINEAR); */
/*  M[11].Xp[0].x = 	0.626	;	 M[11].Xp[0].y =0	; */
/*  M[11].Xp[1].x = 	0.619	;	 M[11].Xp[1].y =0.6	; */
/*  M[11].Xp[2].x = 	0.578	;	 M[11].Xp[2].y =1.2	; */
/*  M[11].Xp[3].x = 	0.506	;	 M[11].Xp[3].y =1.8	; */
/*  M[11].Xp[4].x = 	0.418	;	 M[11].Xp[4].y =2.4	; */
/*  M[11].Xp[5].x = 	0.371	;	 M[11].Xp[5].y =3	; */
/*  M[11].Xp[6].x = 	0.318	;	 M[11].Xp[6].y =3.6	; */
/*  M[11].Xp[7].x = 	0.274	;	 M[11].Xp[7].y =4.2	; */
/*  M[11].Xp[8].x = 	0.242	;	 M[11].Xp[8].y =4.8	; */
/*  M[11].Xp[9].x = 	0.23	;	 M[11].Xp[9].y =5.4	; */
/*  M[11].Xp[10].x = 	0.209	;	 M[11].Xp[10].y =6	; */
/*  M[11].Xp[11].x = 	0.198	;	 M[11].Xp[11].y =6.6	; */
/*  M[11].Xp[12].x = 	0.193	;	 M[11].Xp[12].y =7.2	; */
/*  M[11].Xp[13].x = 	0.186	;	 M[11].Xp[13].y =7.8	; */
/*  M[11].Xp[14].x = 	0.184	;	 M[11].Xp[14].y =8.4	; */
/*  M[11].Xp[15].x = 	0.17	;	 M[11].Xp[15].y =9	; */
/*  M[11].Xp[16].x = 	0.159	;	 M[11].Xp[16].y =9.6	; */
/*  M[11].Xp[17].x = 	0.15	;	 M[11].Xp[17].y =10.2	; */
/*  M[11].Xp[18].x = 	0.14	;	 M[11].Xp[18].y =10.8	; */
/*  M[11].Xp[19].x = 	0.149	;	 M[11].Xp[19].y =11.4	; */
/*  M[11].Xp[20].x = 	0.123	;	 M[11].Xp[20].y =12	; */
/*  M[11].Xp[21].x = 	0.111	;	 M[11].Xp[21].y =12.6	; */
/*  M[11].Xp[22].x = 	0.12	;	 M[11].Xp[22].y =13.2	; */
/*  M[11].Xp[23].x = 	0.125	;	 M[11].Xp[23].y =13.8	; */
/*  M[11].Xp[24].x = 	0.118	;	 M[11].Xp[24].y =14.4	; */
/*  M[11].Xp[25].x = 	0.114	;	 M[11].Xp[25].y =15	; */
/*  M[11].Xp[26].x = 	0.116	;	 M[11].Xp[26].y =15.6	; */
/*  M[11].Xp[27].x = 	0.118	;	 M[11].Xp[27].y =16.2	; */
/*  M[11].Xp[28].x = 	0.107	;	 M[11].Xp[28].y =16.8	; */
/*  M[11].Xp[29].x = 	0.108	;	 M[11].Xp[29].y =17.4	; */
/*  M[11].Xp[30].x = 	0.104	;	 M[11].Xp[30].y =18	; */
/*  M[11].Xp[31].x = 	0.088	;	 M[11].Xp[31].y =18.6	; */
/*  M[11].Xp[32].x = 	0.096	;	 M[11].Xp[32].y =19.2	; */
/*  M[11].Xp[33].x = 	0.071	;	 M[11].Xp[33].y =19.8	; */
/*  M[11].Xp[34].x = 	0.061	;	 M[11].Xp[34].y =20.4	; */
/*  M[11].Xp[35].x = 	0.044	;	 M[11].Xp[35].y =21	; */
/*  M[11].Xp[36].x = 	0.036	;	 M[11].Xp[36].y =21.6	; */
/*  M[11].Xp[37].x = 	0.027	;	 M[11].Xp[37].y =22.2	; */
/*  M[11].Xp[38].x = 	0.018	;	 M[11].Xp[38].y =22.8	; */
/*  M[11].Xp[39].x = 	0.011	;	 M[11].Xp[39].y =23.4	; */
/*  M[11].Xp[40].x = 	0.004	;	 M[11].Xp[40].y =24	; */
/*  M[11].Xp[41].x = 	0.004	;	 M[11].Xp[41].y =24.6	; */
/*  M[11].Xp[42].x = 	0.002	;	 M[11].Xp[42].y =25.2	; */
/*  M[11].Xp[43].x = 	0.002	;	 M[11].Xp[43].y =25.8	; */

/*  M[11].pathlength(); */

/* //station 13 */
/*  x[12] = 29.00; */
/*  M[12].allocate(44); */
/*  M[12].settype(SPLINE2D_LINEAR); */
/*  M[12].Xp[0].x = 	0.603	;	 M[12].Xp[0].y =0	; */
/*  M[12].Xp[1].x = 	0.584	;	 M[12].Xp[1].y =0.6	; */
/*  M[12].Xp[2].x = 	0.553	;	 M[12].Xp[2].y =1.2	; */
/*  M[12].Xp[3].x = 	0.501	;	 M[12].Xp[3].y =1.8	; */
/*  M[12].Xp[4].x = 	0.388	;	 M[12].Xp[4].y =2.4	; */
/*  M[12].Xp[5].x = 	0.36	;	 M[12].Xp[5].y =3	; */
/*  M[12].Xp[6].x = 	0.3	;	 M[12].Xp[6].y =3.6	; */
/*  M[12].Xp[7].x = 	0.283	;	 M[12].Xp[7].y =4.2	; */
/*  M[12].Xp[8].x = 	0.248	;	 M[12].Xp[8].y =4.8	; */
/*  M[12].Xp[9].x = 	0.242	;	 M[12].Xp[9].y =5.4	; */
/*  M[12].Xp[10].x = 	0.21	;	 M[12].Xp[10].y =6	; */
/*  M[12].Xp[11].x = 	0.195	;	 M[12].Xp[11].y =6.6	; */
/*  M[12].Xp[12].x = 	0.199	;	 M[12].Xp[12].y =7.2	; */
/*  M[12].Xp[13].x = 	0.179	;	 M[12].Xp[13].y =7.8	; */
/*  M[12].Xp[14].x = 	0.174	;	 M[12].Xp[14].y =8.4	; */
/*  M[12].Xp[15].x = 	0.168	;	 M[12].Xp[15].y =9	; */
/*  M[12].Xp[16].x = 	0.151	;	 M[12].Xp[16].y =9.6	; */
/*  M[12].Xp[17].x = 	0.142	;	 M[12].Xp[17].y =10.2	; */
/*  M[12].Xp[18].x = 	0.138	;	 M[12].Xp[18].y =10.8	; */
/*  M[12].Xp[19].x = 	0.133	;	 M[12].Xp[19].y =11.4	; */
/*  M[12].Xp[20].x = 	0.132	;	 M[12].Xp[20].y =12	; */
/*  M[12].Xp[21].x = 	0.122	;	 M[12].Xp[21].y =12.6	; */
/*  M[12].Xp[22].x = 	0.114	;	 M[12].Xp[22].y =13.2	; */
/*  M[12].Xp[23].x = 	0.115	;	 M[12].Xp[23].y =13.8	; */
/*  M[12].Xp[24].x = 	0.11	;	 M[12].Xp[24].y =14.4	; */
/*  M[12].Xp[25].x = 	0.114	;	 M[12].Xp[25].y =15	; */
/*  M[12].Xp[26].x = 	0.112	;	 M[12].Xp[26].y =15.6	; */
/*  M[12].Xp[27].x = 	0.109	;	 M[12].Xp[27].y =16.2	; */
/*  M[12].Xp[28].x = 	0.11	;	 M[12].Xp[28].y =16.8	; */
/*  M[12].Xp[29].x = 	0.106	;	 M[12].Xp[29].y =17.4	; */
/*  M[12].Xp[30].x = 	0.103	;	 M[12].Xp[30].y =18	; */
/*  M[12].Xp[31].x = 	0.082	;	 M[12].Xp[31].y =18.6	; */
/*  M[12].Xp[32].x = 	0.075	;	 M[12].Xp[32].y =19.2	; */
/*  M[12].Xp[33].x = 	0.071	;	 M[12].Xp[33].y =19.8	; */
/*  M[12].Xp[34].x = 	0.052	;	 M[12].Xp[34].y =20.4	; */
/*  M[12].Xp[35].x = 	0.043	;	 M[12].Xp[35].y =21	; */
/*  M[12].Xp[36].x = 	0.028	;	 M[12].Xp[36].y =21.6	; */
/*  M[12].Xp[37].x = 	0.016	;	 M[12].Xp[37].y =22.2	; */
/*  M[12].Xp[38].x = 	0.012	;	 M[12].Xp[38].y =22.8	; */
/*  M[12].Xp[39].x = 	0.007	;	 M[12].Xp[39].y =23.4	; */
/*  M[12].Xp[40].x = 	0.004	;	 M[12].Xp[40].y =24	; */
/*  M[12].Xp[41].x = 	0.002	;	 M[12].Xp[41].y =24.6	; */
/*  M[12].Xp[42].x = 	0.001	;	 M[12].Xp[42].y =25.2	; */
/*  M[12].Xp[43].x = 	0	;	 M[12].Xp[43].y =25.8	; */

/*  M[12].pathlength(); */

/*  cout.flush(); */
/* //station 14 */
/*  x[13] = 29.80; */
/*  M[13].allocate(44); */
/*  M[13].settype(SPLINE2D_LINEAR); */
/*  M[13].Xp[0].x = 	0.631	;	 M[13].Xp[0].y =0	; */
/*  M[13].Xp[1].x = 	0.579	;	 M[13].Xp[1].y =0.6	; */
/*  M[13].Xp[2].x = 	0.579	;	 M[13].Xp[2].y =1.2	; */
/*  M[13].Xp[3].x = 	0.476	;	 M[13].Xp[3].y =1.8	; */
/*  M[13].Xp[4].x = 	0.42	;	 M[13].Xp[4].y =2.4	; */
/*  M[13].Xp[5].x = 	0.361	;	 M[13].Xp[5].y =3	; */
/*  M[13].Xp[6].x = 	0.295	;	 M[13].Xp[6].y =3.6	; */
/*  M[13].Xp[7].x = 	0.288	;	 M[13].Xp[7].y =4.2	; */
/*  M[13].Xp[8].x = 	0.25	;	 M[13].Xp[8].y =4.8	; */
/*  M[13].Xp[9].x = 	0.254	;	 M[13].Xp[9].y =5.4	; */
/*  M[13].Xp[10].x = 	0.218	;	 M[13].Xp[10].y =6	; */
/*  M[13].Xp[11].x = 	0.205	;	 M[13].Xp[11].y =6.6	; */
/*  M[13].Xp[12].x = 	0.189	;	 M[13].Xp[12].y =7.2	; */
/*  M[13].Xp[13].x = 	0.189	;	 M[13].Xp[13].y =7.8	; */
/*  M[13].Xp[14].x = 	0.181	;	 M[13].Xp[14].y =8.4	; */
/*  M[13].Xp[15].x = 	0.183	;	 M[13].Xp[15].y =9	; */
/*  M[13].Xp[16].x = 	0.134	;	 M[13].Xp[16].y =9.6	; */
/*  M[13].Xp[17].x = 	0.155	;	 M[13].Xp[17].y =10.2	; */
/*  M[13].Xp[18].x = 	0.135	;	 M[13].Xp[18].y =10.8	; */
/*  M[13].Xp[19].x = 	0.127	;	 M[13].Xp[19].y =11.4	; */
/*  M[13].Xp[20].x = 	0.127	;	 M[13].Xp[20].y =12	; */
/*  M[13].Xp[21].x = 	0.119	;	 M[13].Xp[21].y =12.6	; */
/*  M[13].Xp[22].x = 	0.115	;	 M[13].Xp[22].y =13.2	; */
/*  M[13].Xp[23].x = 	0.119	;	 M[13].Xp[23].y =13.8	; */
/*  M[13].Xp[24].x = 	0.113	;	 M[13].Xp[24].y =14.4	; */
/*  M[13].Xp[25].x = 	0.121	;	 M[13].Xp[25].y =15	; */
/*  M[13].Xp[26].x = 	0.113	;	 M[13].Xp[26].y =15.6	; */
/*  M[13].Xp[27].x = 	0.104	;	 M[13].Xp[27].y =16.2	; */
/*  M[13].Xp[28].x = 	0.103	;	 M[13].Xp[28].y =16.8	; */
/*  M[13].Xp[29].x = 	0.106	;	 M[13].Xp[29].y =17.4	; */
/*  M[13].Xp[30].x = 	0.088	;	 M[13].Xp[30].y =18	; */
/*  M[13].Xp[31].x = 	0.089	;	 M[13].Xp[31].y =18.6	; */
/*  M[13].Xp[32].x = 	0.084	;	 M[13].Xp[32].y =19.2	; */
/*  M[13].Xp[33].x = 	0.069	;	 M[13].Xp[33].y =19.8	; */
/*  M[13].Xp[34].x = 	0.055	;	 M[13].Xp[34].y =20.4	; */
/*  M[13].Xp[35].x = 	0.037	;	 M[13].Xp[35].y =21	; */
/*  M[13].Xp[36].x = 	0.029	;	 M[13].Xp[36].y =21.6	; */
/*  M[13].Xp[37].x = 	0.023	;	 M[13].Xp[37].y =22.2	; */
/*  M[13].Xp[38].x = 	0.009	;	 M[13].Xp[38].y =22.8	; */
/*  M[13].Xp[39].x = 	0.005	;	 M[13].Xp[39].y =23.4	; */
/*  M[13].Xp[40].x = 	0.003	;	 M[13].Xp[40].y =24	; */
/*  M[13].Xp[41].x = 	0.001	;	 M[13].Xp[41].y =24.6	; */
/*  M[13].Xp[42].x = 	0.002	;	 M[13].Xp[42].y =25.2	; */
/*  M[13].Xp[43].x = 	0.002	;	 M[13].Xp[43].y =25.8	; */
/*  M[13].pathlength(); */
 
/*  cout.flush(); */
/* //station 15 */
/*  x[14] = 30.50; */
/*  M[14].allocate(44); */
/*  M[14].settype(SPLINE2D_LINEAR); */
/*  M[14].Xp[0].x = 	0.597	;	 M[14].Xp[0].y =0	; */
/*  M[14].Xp[1].x = 	0.541	;	 M[14].Xp[1].y =0.6	; */
/*  M[14].Xp[2].x = 	0.53	;	 M[14].Xp[2].y =1.2	; */
/*  M[14].Xp[3].x = 	0.483	;	 M[14].Xp[3].y =1.8	; */
/*  M[14].Xp[4].x = 	0.42	;	 M[14].Xp[4].y =2.4	; */
/*  M[14].Xp[5].x = 	0.378	;	 M[14].Xp[5].y =3	; */
/*  M[14].Xp[6].x = 	0.313	;	 M[14].Xp[6].y =3.6	; */
/*  M[14].Xp[7].x = 	0.271	;	 M[14].Xp[7].y =4.2	; */
/*  M[14].Xp[8].x = 	0.265	;	 M[14].Xp[8].y =4.8	; */
/*  M[14].Xp[9].x = 	0.238	;	 M[14].Xp[9].y =5.4	; */
/*  M[14].Xp[1].x = 	0.235	;	 M[14].Xp[10].y =6	; */
/*  M[14].Xp[11].x = 	0.213	;	 M[14].Xp[11].y =6.6	; */
/*  M[14].Xp[12].x = 	0.184	;	 M[14].Xp[12].y =7.2	; */
/*  M[14].Xp[13].x = 	0.194	;	 M[14].Xp[13].y =7.8	; */
/*  M[14].Xp[14].x = 	0.172	;	 M[14].Xp[14].y =8.4	; */
/*  M[14].Xp[15].x = 	0.154	;	 M[14].Xp[15].y =9	; */
/*  M[14].Xp[16].x = 	0.141	;	 M[14].Xp[16].y =9.6	; */
/*  M[14].Xp[17].x = 	0.134	;	 M[14].Xp[17].y =10.2	; */
/*  M[14].Xp[18].x = 	0.137	;	 M[14].Xp[18].y =10.8	; */
/*  M[14].Xp[19].x = 	0.119	;	 M[14].Xp[19].y =11.4	; */
/*  M[14].Xp[20].x = 	0.12	;	 M[14].Xp[20].y =12	; */
/*  M[14].Xp[21].x = 	0.123	;	 M[14].Xp[21].y =12.6	; */
/*  M[14].Xp[22].x = 	0.12	;	 M[14].Xp[22].y =13.2	; */
/*  M[14].Xp[23].x = 	0.104	;	 M[14].Xp[23].y =13.8	; */
/*  M[14].Xp[24].x = 	0.102	;	 M[14].Xp[24].y =14.4	; */
/*  M[14].Xp[25].x = 	0.102	;	 M[14].Xp[25].y =15	; */
/*  M[14].Xp[26].x = 	0.121	;	 M[14].Xp[26].y =15.6	; */
/*  M[14].Xp[27].x = 	0.101	;	 M[14].Xp[27].y =16.2	; */
/*  M[14].Xp[28].x = 	0.098	;	 M[14].Xp[28].y =16.8	; */
/*  M[14].Xp[29].x = 	0.096	;	 M[14].Xp[29].y =17.4	; */
/*  M[14].Xp[30].x = 	0.088	;	 M[14].Xp[30].y =18	; */
/*  M[14].Xp[31].x = 	0.079	;	 M[14].Xp[31].y =18.6	; */
/*  M[14].Xp[32].x = 	0.082	;	 M[14].Xp[32].y =19.2	; */
/*  M[14].Xp[33].x = 	0.063	;	 M[14].Xp[33].y =19.8	; */
/*  M[14].Xp[34].x = 	0.047	;	 M[14].Xp[34].y =20.4	; */
/*  M[14].Xp[35].x = 	0.028	;	 M[14].Xp[35].y =21	; */
/*  M[14].Xp[36].x = 	0.02	;	 M[14].Xp[36].y =21.6	; */
/*  M[14].Xp[37].x = 	0.015	;	 M[14].Xp[37].y =22.2	; */
/*  M[14].Xp[38].x = 	0.008	;	 M[14].Xp[38].y =22.8	; */
/*  M[14].Xp[39].x = 	0.004	;	 M[14].Xp[39].y =23.4	; */
/*  M[14].Xp[40].x = 	0.001	;	 M[14].Xp[40].y =24	; */
/*  M[14].Xp[41].x = 	0.004	;	 M[14].Xp[41].y =24.6	; */
/*  M[14].Xp[42].x = 	0.01	;	 M[14].Xp[42].y =25.2	; */
/*  M[14].Xp[43].x = 	0.004	;	 M[14].Xp[43].y =25.8	; */
/*  M[14].pathlength(); */
 
/*  cout.flush(); */
/* //station 16 */
/*  x[15] = 31.30; */
/*  M[15].allocate(44); */
/*  M[15].settype(SPLINE2D_LINEAR); */
/*  M[15].Xp[0].x = 	0.577	;	 M[15].Xp[0].y =0	; */
/*  M[15].Xp[1].x = 	0.553	;	 M[15].Xp[1].y =0.6	; */
/*  M[15].Xp[2].x = 	0.545	;	 M[15].Xp[2].y =1.2	; */
/*  M[15].Xp[3].x = 	0.479	;	 M[15].Xp[3].y =1.8	; */
/*  M[15].Xp[4].x = 	0.414	;	 M[15].Xp[4].y =2.4	; */
/*  M[15].Xp[5].x = 	0.333	;	 M[15].Xp[5].y =3	; */
/*  M[15].Xp[6].x = 	0.326	;	 M[15].Xp[6].y =3.6	; */
/*  M[15].Xp[7].x = 	0.27	;	 M[15].Xp[7].y =4.2	; */
/*  M[15].Xp[8].x = 	0.254	;	 M[15].Xp[8].y =4.8	; */
/*  M[15].Xp[9].x = 	0.248	;	 M[15].Xp[9].y =5.4	; */
/*  M[15].Xp[10].x = 	0.237	;	 M[15].Xp[10].y =6	; */
/*  M[15].Xp[11].x = 	0.187	;	 M[15].Xp[11].y =6.6	; */
/*  M[15].Xp[12].x = 	0.196	;	 M[15].Xp[12].y =7.2	; */
/*  M[15].Xp[13].x = 	0.207	;	 M[15].Xp[13].y =7.8	; */
/*  M[15].Xp[14].x = 	0.195	;	 M[15].Xp[14].y =8.4	; */
/*  M[15].Xp[15].x = 	0.162	;	 M[15].Xp[15].y =9	; */
/*  M[15].Xp[16].x = 	0.156	;	 M[15].Xp[16].y =9.6	; */
/*  M[15].Xp[17].x = 	0.137	;	 M[15].Xp[17].y =10.2	; */
/*  M[15].Xp[18].x = 	0.138	;	 M[15].Xp[18].y =10.8	; */
/*  M[15].Xp[19].x = 	0.122	;	 M[15].Xp[19].y =11.4	; */
/*  M[15].Xp[20].x = 	0.121	;	 M[15].Xp[20].y =12	; */
/*  M[15].Xp[21].x = 	0.135	;	 M[15].Xp[21].y =12.6	; */
/*  M[15].Xp[22].x = 	0.127	;	 M[15].Xp[22].y =13.2	; */
/*  M[15].Xp[23].x = 	0.114	;	 M[15].Xp[23].y =13.8	; */
/*  M[15].Xp[24].x = 	0.105	;	 M[15].Xp[24].y =14.4	; */
/*  M[15].Xp[25].x = 	0.114	;	 M[15].Xp[25].y =15	; */
/*  M[15].Xp[26].x = 	0.105	;	 M[15].Xp[26].y =15.6	; */
/*  M[15].Xp[27].x = 	0.096	;	 M[15].Xp[27].y =16.2	; */
/*  M[15].Xp[28].x = 	0.105	;	 M[15].Xp[28].y =16.8	; */
/*  M[15].Xp[29].x = 	0.106	;	 M[15].Xp[29].y =17.4	; */
/*  M[15].Xp[30].x = 	0.091	;	 M[15].Xp[30].y =18	; */
/*  M[15].Xp[31].x = 	0.089	;	 M[15].Xp[31].y =18.6	; */
/*  M[15].Xp[32].x = 	0.083	;	 M[15].Xp[32].y =19.2	; */
/*  M[15].Xp[33].x = 	0.065	;	 M[15].Xp[33].y =19.8	; */
/*  M[15].Xp[34].x = 	0.05	;	 M[15].Xp[34].y =20.4	; */
/*  M[15].Xp[35].x = 	0.036	;	 M[15].Xp[35].y =21	; */
/*  M[15].Xp[36].x = 	0.017	;	 M[15].Xp[36].y =21.6	; */
/*  M[15].Xp[37].x = 	0.021	;	 M[15].Xp[37].y =22.2	; */
/*  M[15].Xp[38].x = 	0.006	;	 M[15].Xp[38].y =22.8	; */
/*  M[15].Xp[39].x = 	0.004	;	 M[15].Xp[39].y =23.4	; */
/*  M[15].Xp[40].x = 	0.004	;	 M[15].Xp[40].y =24	; */
/*  M[15].Xp[41].x = 	0.002	;	 M[15].Xp[41].y =24.6	; */
/*  M[15].Xp[42].x = 	0.002	;	 M[15].Xp[42].y =25.2	; */
/*  M[15].Xp[43].x = 	0.004	;	 M[15].Xp[43].y =25.8	; */
/*  M[15].pathlength(); */

/*  cout.flush(); */
/* //station 17 */
/*  x[16] = 32.00; */
/*  M[16].allocate(44); */
/*  M[16].settype(SPLINE2D_LINEAR); */
/*  M[16].Xp[0].x = 	0.571	;	 M[16].Xp[0].y =0	; */
/*  M[16].Xp[1].x = 	0.528	;	 M[16].Xp[1].y =0.6	; */
/*  M[16].Xp[2].x = 	0.505	;	 M[16].Xp[2].y =1.2	; */
/*  M[16].Xp[3].x = 	0.524	;	 M[16].Xp[3].y =1.8	; */
/*  M[16].Xp[4].x = 	0.412	;	 M[16].Xp[4].y =2.4	; */
/*  M[16].Xp[5].x = 	0.381	;	 M[16].Xp[5].y =3	; */
/*  M[16].Xp[6].x = 	0.315	;	 M[16].Xp[6].y =3.6	; */
/*  M[16].Xp[7].x = 	0.321	;	 M[16].Xp[7].y =4.2	; */
/*  M[16].Xp[8].x = 	0.277	;	 M[16].Xp[8].y =4.8	; */
/*  M[16].Xp[9].x = 	0.248	;	 M[16].Xp[9].y =5.4	; */
/*  M[16].Xp[10].x = 	0.217	;	 M[16].Xp[10].y =6	; */
/*  M[16].Xp[11].x = 	0.215	;	 M[16].Xp[11].y =6.6	; */
/*  M[16].Xp[12].x = 	0.22	;	 M[16].Xp[12].y =7.2	; */
/*  M[16].Xp[13].x = 	0.188	;	 M[16].Xp[13].y =7.8	; */
/*  M[16].Xp[14].x = 	0.175	;	 M[16].Xp[14].y =8.4	; */
/*  M[16].Xp[15].x = 	0.165	;	 M[16].Xp[15].y =9	; */
/*  M[16].Xp[16].x = 	0.151	;	 M[16].Xp[16].y =9.6	; */
/*  M[16].Xp[17].x = 	0.156	;	 M[16].Xp[17].y =10.2	; */
/*  M[16].Xp[18].x = 	0.136	;	 M[16].Xp[18].y =10.8	; */
/*  M[16].Xp[19].x = 	0.121	;	 M[16].Xp[19].y =11.4	; */
/*  M[16].Xp[20].x = 	0.144	;	 M[16].Xp[20].y =12	; */
/*  M[16].Xp[21].x = 	0.123	;	 M[16].Xp[21].y =12.6	; */
/*  M[16].Xp[22].x = 	0.114	;	 M[16].Xp[22].y =13.2	; */
/*  M[16].Xp[23].x = 	0.113	;	 M[16].Xp[23].y =13.8	; */
/*  M[16].Xp[24].x = 	0.109	;	 M[16].Xp[24].y =14.4	; */
/*  M[16].Xp[25].x = 	0.084	;	 M[16].Xp[25].y =15	; */
/*  M[16].Xp[26].x = 	0.083	;	 M[16].Xp[26].y =15.6	; */
/*  M[16].Xp[27].x = 	0.098	;	 M[16].Xp[27].y =16.2	; */
/*  M[16].Xp[28].x = 	0.092	;	 M[16].Xp[28].y =16.8	; */
/*  M[16].Xp[29].x = 	0.09	;	 M[16].Xp[29].y =17.4	; */
/*  M[16].Xp[30].x = 	0.096	;	 M[16].Xp[30].y =18	; */
/*  M[16].Xp[31].x = 	0.074	;	 M[16].Xp[31].y =18.6	; */
/*  M[16].Xp[32].x = 	0.081	;	 M[16].Xp[32].y =19.2	; */
/*  M[16].Xp[33].x = 	0.067	;	 M[16].Xp[33].y =19.8	; */
/*  M[16].Xp[34].x = 	0.049	;	 M[16].Xp[34].y =20.4	; */
/*  M[16].Xp[35].x = 	0.031	;	 M[16].Xp[35].y =21	; */
/*  M[16].Xp[36].x = 	0.029	;	 M[16].Xp[36].y =21.6	; */
/*  M[16].Xp[37].x = 	0.018	;	 M[16].Xp[37].y =22.2	; */
/*  M[16].Xp[38].x = 	0.008	;	 M[16].Xp[38].y =22.8	; */
/*  M[16].Xp[39].x = 	0.007	;	 M[16].Xp[39].y =23.4	; */
/*  M[16].Xp[40].x = 	0.007	;	 M[16].Xp[40].y =24	; */
/*  M[16].Xp[41].x = 	0.008	;	 M[16].Xp[41].y =24.6	; */
/*  M[16].Xp[42].x = 	0.005	;	 M[16].Xp[42].y =25.2	; */
/*  M[16].Xp[43].x = 	0.004	;	 M[16].Xp[43].y =25.8	; */
/*  M[16].pathlength(); */

/* //station 18 */
/*  x[17] = 32.80; */
/*  M[17].allocate(44); */
/*  M[17].settype(SPLINE2D_LINEAR); */
/*  M[17].Xp[0].x = 	0.545	;	 M[17].Xp[0].y =0	; */
/*  M[17].Xp[1].x = 	0.532	;	 M[17].Xp[1].y =0.6	; */
/*  M[17].Xp[2].x = 	0.475	;	 M[17].Xp[2].y =1.2	; */
/*  M[17].Xp[3].x = 	0.479	;	 M[17].Xp[3].y =1.8	; */
/*  M[17].Xp[4].x = 	0.399	;	 M[17].Xp[4].y =2.4	; */
/*  M[17].Xp[5].x = 	0.372	;	 M[17].Xp[5].y =3	; */
/*  M[17].Xp[6].x = 	0.329	;	 M[17].Xp[6].y =3.6	; */
/*  M[17].Xp[7].x = 	0.303	;	 M[17].Xp[7].y =4.2	; */
/*  M[17].Xp[8].x = 	0.26	;	 M[17].Xp[8].y =4.8	; */
/*  M[17].Xp[9].x = 	0.207	;	 M[17].Xp[9].y =5.4	; */
/*  M[17].Xp[10].x = 	0.22	;	 M[17].Xp[10].y =6	; */
/*  M[17].Xp[11].x = 	0.222	;	 M[17].Xp[11].y =6.6	; */
/*  M[17].Xp[12].x = 	0.18	;	 M[17].Xp[12].y =7.2	; */
/*  M[17].Xp[13].x = 	0.191	;	 M[17].Xp[13].y =7.8	; */
/*  M[17].Xp[14].x = 	0.159	;	 M[17].Xp[14].y =8.4	; */
/*  M[17].Xp[15].x = 	0.157	;	 M[17].Xp[15].y =9	; */
/*  M[17].Xp[16].x = 	0.152	;	 M[17].Xp[16].y =9.6	; */
/*  M[17].Xp[17].x = 	0.147	;	 M[17].Xp[17].y =10.2	; */
/*  M[17].Xp[18].x = 	0.145	;	 M[17].Xp[18].y =10.8	; */
/*  M[17].Xp[19].x = 	0.123	;	 M[17].Xp[19].y =11.4	; */
/*  M[17].Xp[20].x = 	0.121	;	 M[17].Xp[20].y =12	; */
/*  M[17].Xp[21].x = 	0.111	;	 M[17].Xp[21].y =12.6	; */
/*  M[17].Xp[22].x = 	0.107	;	 M[17].Xp[22].y =13.2	; */
/*  M[17].Xp[23].x = 	0.1	;	 M[17].Xp[23].y =13.8	; */
/*  M[17].Xp[24].x = 	0.085	;	 M[17].Xp[24].y =14.4	; */
/*  M[17].Xp[25].x = 	0.098	;	 M[17].Xp[25].y =15	; */
/*  M[17].Xp[26].x = 	0.083	;	 M[17].Xp[26].y =15.6	; */
/*  M[17].Xp[27].x = 	0.097	;	 M[17].Xp[27].y =16.2	; */
/*  M[17].Xp[28].x = 	0.097	;	 M[17].Xp[28].y =16.8	; */
/*  M[17].Xp[29].x = 	0.074	;	 M[17].Xp[29].y =17.4	; */
/*  M[17].Xp[30].x = 	0.074	;	 M[17].Xp[30].y =18	; */
/*  M[17].Xp[31].x = 	0.094	;	 M[17].Xp[31].y =18.6	; */
/*  M[17].Xp[32].x = 	0.081	;	 M[17].Xp[32].y =19.2	; */
/*  M[17].Xp[33].x = 	0.054	;	 M[17].Xp[33].y =19.8	; */
/*  M[17].Xp[34].x = 	0.054	;	 M[17].Xp[34].y =20.4	; */
/*  M[17].Xp[35].x = 	0.032	;	 M[17].Xp[35].y =21	; */
/*  M[17].Xp[36].x = 	0.026	;	 M[17].Xp[36].y =21.6	; */
/*  M[17].Xp[37].x = 	0.013	;	 M[17].Xp[37].y =22.2	; */
/*  M[17].Xp[38].x = 	0.01	;	 M[17].Xp[38].y =22.8	; */
/*  M[17].Xp[39].x = 	0.005	;	 M[17].Xp[39].y =23.4	; */
/*  M[17].Xp[40].x = 	0.012	;	 M[17].Xp[40].y =24	; */
/*  M[17].Xp[41].x = 	0.007	;	 M[17].Xp[41].y =24.6	; */
/*  M[17].Xp[42].x = 	0.005	;	 M[17].Xp[42].y =25.2	; */
/*  M[17].Xp[43].x = 	0.01	;	 M[17].Xp[43].y =25.8	; */
/*  M[17].pathlength(); */
/* //station 19 */
/*  x[18] = 40.0; */
/*  M[18].allocate(44); */
/*  M[18].settype(SPLINE2D_LINEAR); */
/*  M[18].Xp[0].x = 	0.474	;	 M[18].Xp[0].y =0	; */
/*  M[18].Xp[1].x = 	0.407	;	 M[18].Xp[1].y =0.6	; */
/*  M[18].Xp[2].x = 	0.405	;	 M[18].Xp[2].y =1.2	; */
/*  M[18].Xp[3].x = 	0.442	;	 M[18].Xp[3].y =1.8	; */
/*  M[18].Xp[4].x = 	0.418	;	 M[18].Xp[4].y =2.4	; */
/*  M[18].Xp[5].x = 	0.354	;	 M[18].Xp[5].y =3	; */
/*  M[18].Xp[6].x = 	0.381	;	 M[18].Xp[6].y =3.6	; */
/*  M[18].Xp[7].x = 	0.342	;	 M[18].Xp[7].y =4.2	; */
/*  M[18].Xp[8].x = 	0.291	;	 M[18].Xp[8].y =4.8	; */
/*  M[18].Xp[9].x = 	0.291	;	 M[18].Xp[9].y =5.4	; */
/*  M[18].Xp[10].x = 	0.248	;	 M[18].Xp[10].y =6	; */
/*  M[18].Xp[11].x = 	0.209	;	 M[18].Xp[11].y =6.6	; */
/*  M[18].Xp[12].x = 	0.178	;	 M[18].Xp[12].y =7.2	; */
/*  M[18].Xp[13].x = 	0.163	;	 M[18].Xp[13].y =7.8	; */
/*  M[18].Xp[14].x = 	0.15	;	 M[18].Xp[14].y =8.4	; */
/*  M[18].Xp[15].x = 	0.133	;	 M[18].Xp[15].y =9	; */
/*  M[18].Xp[16].x = 	0.106	;	 M[18].Xp[16].y =9.6	; */
/*  M[18].Xp[17].x = 	0.12	;	 M[18].Xp[17].y =10.2	; */
/*  M[18].Xp[18].x = 	0.09	;	 M[18].Xp[18].y =10.8	; */
/*  M[18].Xp[19].x = 	0.075	;	 M[18].Xp[19].y =11.4	; */
/*  M[18].Xp[20].x = 	0.079	;	 M[18].Xp[20].y =12	; */
/*  M[18].Xp[21].x = 	0.064	;	 M[18].Xp[21].y =12.6	; */
/*  M[18].Xp[22].x = 	0.06	;	 M[18].Xp[22].y =13.2	; */
/*  M[18].Xp[23].x = 	0.048	;	 M[18].Xp[23].y =13.8	; */
/*  M[18].Xp[24].x = 	0.074	;	 M[18].Xp[24].y =14.4	; */
/*  M[18].Xp[25].x = 	0.061	;	 M[18].Xp[25].y =15	; */
/*  M[18].Xp[26].x = 	0.045	;	 M[18].Xp[26].y =15.6	; */
/*  M[18].Xp[27].x = 	0.06	;	 M[18].Xp[27].y =16.2	; */
/*  M[18].Xp[28].x = 	0.053	;	 M[18].Xp[28].y =16.8	; */
/*  M[18].Xp[29].x = 	0.041	;	 M[18].Xp[29].y =17.4	; */
/*  M[18].Xp[30].x = 	0.035	;	 M[18].Xp[30].y =18	; */
/*  M[18].Xp[31].x = 	0.04	;	 M[18].Xp[31].y =18.6	; */
/*  M[18].Xp[32].x = 	0.036	;	 M[18].Xp[32].y =19.2	; */
/*  M[18].Xp[33].x = 	0.029	;	 M[18].Xp[33].y =19.8	; */
/*  M[18].Xp[34].x = 	0.025	;	 M[18].Xp[34].y =20.4	; */
/*  M[18].Xp[35].x = 	0.022	;	 M[18].Xp[35].y =21	; */
/*  M[18].Xp[36].x = 	0.027	;	 M[18].Xp[36].y =21.6	; */
/*  M[18].Xp[37].x = 	0.01	;	 M[18].Xp[37].y =22.2	; */
/*  M[18].Xp[38].x = 	0.006	;	 M[18].Xp[38].y =22.8	; */
/*  M[18].Xp[39].x = 	0.008	;	 M[18].Xp[39].y =23.4	; */
/*  M[18].Xp[40].x = 	0.004	;	 M[18].Xp[40].y =24	; */
/*  M[18].Xp[41].x = 	0.007	;	 M[18].Xp[41].y =24.6	; */
/*  M[18].Xp[42].x = 	0.01	;	 M[18].Xp[42].y =25.2	; */
/*  M[18].Xp[43].x = 	0.013	;	 M[18].Xp[43].y =25.8	; */
/*  M[18].pathlength(); */

/* //station 20 */
/*  x[19] = 40.8; */
/*  M[19].allocate(44); */
/*  M[19].settype(SPLINE2D_LINEAR); */
/*  M[19].Xp[0].x = 	0.429	;	 M[19].Xp[0].y =0	; */
/*  M[19].Xp[1].x = 	0.489	;	 M[19].Xp[1].y =0.6	; */
/*  M[19].Xp[2].x = 	0.437	;	 M[19].Xp[2].y =1.2	; */
/*  M[19].Xp[3].x = 	0.459	;	 M[19].Xp[3].y =1.8	; */
/*  M[19].Xp[4].x = 	0.427	;	 M[19].Xp[4].y =2.4	; */
/*  M[19].Xp[5].x = 	0.414	;	 M[19].Xp[5].y =3	; */
/*  M[19].Xp[6].x = 	0.389	;	 M[19].Xp[6].y =3.6	; */
/*  M[19].Xp[7].x = 	0.341	;	 M[19].Xp[7].y =4.2	; */
/*  M[19].Xp[8].x = 	0.336	;	 M[19].Xp[8].y =4.8	; */
/*  M[19].Xp[9].x = 	0.283	;	 M[19].Xp[9].y =5.4	; */
/*  M[19].Xp[10].x = 	0.231	;	 M[19].Xp[10].y =6	; */
/*  M[19].Xp[11].x = 	0.201	;	 M[19].Xp[11].y =6.6	; */
/*  M[19].Xp[12].x = 	0.187	;	 M[19].Xp[12].y =7.2	; */
/*  M[19].Xp[13].x = 	0.168	;	 M[19].Xp[13].y =7.8	; */
/*  M[19].Xp[14].x = 	0.172	;	 M[19].Xp[14].y =8.4	; */
/*  M[19].Xp[15].x = 	0.135	;	 M[19].Xp[15].y =9	; */
/*  M[19].Xp[16].x = 	0.123	;	 M[19].Xp[16].y =9.6	; */
/*  M[19].Xp[17].x = 	0.107	;	 M[19].Xp[17].y =10.2	; */
/*  M[19].Xp[18].x = 	0.095	;	 M[19].Xp[18].y =10.8	; */
/*  M[19].Xp[19].x = 	0.078	;	 M[19].Xp[19].y =11.4	; */
/*  M[19].Xp[20].x = 	0.072	;	 M[19].Xp[20].y =12	; */
/*  M[19].Xp[21].x = 	0.065	;	 M[19].Xp[21].y =12.6	; */
/*  M[19].Xp[22].x = 	0.062	;	 M[19].Xp[22].y =13.2	; */
/*  M[19].Xp[23].x = 	0.052	;	 M[19].Xp[23].y =13.8	; */
/*  M[19].Xp[24].x = 	0.057	;	 M[19].Xp[24].y =14.4	; */
/*  M[19].Xp[25].x = 	0.054	;	 M[19].Xp[25].y =15	; */
/*  M[19].Xp[26].x = 	0.063	;	 M[19].Xp[26].y =15.6	; */
/*  M[19].Xp[27].x = 	0.045	;	 M[19].Xp[27].y =16.2	; */
/*  M[19].Xp[28].x = 	0.037	;	 M[19].Xp[28].y =16.8	; */
/*  M[19].Xp[29].x = 	0.048	;	 M[19].Xp[29].y =17.4	; */
/*  M[19].Xp[30].x = 	0.048	;	 M[19].Xp[30].y =18	; */
/*  M[19].Xp[31].x = 	0.03	;	 M[19].Xp[31].y =18.6	; */
/*  M[19].Xp[32].x = 	0.031	;	 M[19].Xp[32].y =19.2	; */
/*  M[19].Xp[33].x = 	0.038	;	 M[19].Xp[33].y =19.8	; */
/*  M[19].Xp[34].x = 	0.03	;	 M[19].Xp[34].y =20.4	; */
/*  M[19].Xp[35].x = 	0.023	;	 M[19].Xp[35].y =21	; */
/*  M[19].Xp[36].x = 	0.015	;	 M[19].Xp[36].y =21.6	; */
/*  M[19].Xp[37].x = 	0.017	;	 M[19].Xp[37].y =22.2	; */
/*  M[19].Xp[38].x = 	0.015	;	 M[19].Xp[38].y =22.8	; */
/*  M[19].Xp[39].x = 	0.009	;	 M[19].Xp[39].y =23.4	; */
/*  M[19].Xp[40].x = 	0.004	;	 M[19].Xp[40].y =24	; */
/*  M[19].Xp[41].x = 	0.004	;	 M[19].Xp[41].y =24.6	; */
/*  M[19].Xp[42].x = 	0.008	;	 M[19].Xp[42].y =25.2	; */
/*  M[19].Xp[43].x = 	0.002	;	 M[19].Xp[43].y =25.8	; */
/*  M[19].pathlength(); */
/* //station 21 */
/*  x[20] = 41.5; */
/*  M[20].allocate(44); */
/*  M[20].settype(SPLINE2D_LINEAR); */
/*  M[20].Xp[0].x = 	0.457	;	 M[20].Xp[0].y =0	; */
/*  M[20].Xp[1].x = 	0.458	;	 M[20].Xp[1].y =0.6	; */
/*  M[20].Xp[2].x = 	0.421	;	 M[20].Xp[2].y =1.2	; */
/*  M[20].Xp[3].x = 	0.454	;	 M[20].Xp[3].y =1.8	; */
/*  M[20].Xp[4].x = 	0.416	;	 M[20].Xp[4].y =2.4	; */
/*  M[20].Xp[5].x = 	0.386	;	 M[20].Xp[5].y =3	; */
/*  M[20].Xp[6].x = 	0.383	;	 M[20].Xp[6].y =3.6	; */
/*  M[20].Xp[7].x = 	0.364	;	 M[20].Xp[7].y =4.2	; */
/*  M[20].Xp[8].x = 	0.327	;	 M[20].Xp[8].y =4.8	; */
/*  M[20].Xp[9].x = 	0.292	;	 M[20].Xp[9].y =5.4	; */
/*  M[20].Xp[10].x = 	0.248	;	 M[20].Xp[10].y =6	; */
/*  M[20].Xp[11].x = 	0.231	;	 M[20].Xp[11].y =6.6	; */
/*  M[20].Xp[12].x = 	0.221	;	 M[20].Xp[12].y =7.2	; */
/*  M[20].Xp[13].x = 	0.16	;	 M[20].Xp[13].y =7.8	; */
/*  M[20].Xp[14].x = 	0.148	;	 M[20].Xp[14].y =8.4	; */
/*  M[20].Xp[15].x = 	0.122	;	 M[20].Xp[15].y =9	; */
/*  M[20].Xp[16].x = 	0.125	;	 M[20].Xp[16].y =9.6	; */
/*  M[20].Xp[17].x = 	0.1	;	 M[20].Xp[17].y =10.2	; */
/*  M[20].Xp[18].x = 	0.085	;	 M[20].Xp[18].y =10.8	; */
/*  M[20].Xp[19].x = 	0.084	;	 M[20].Xp[19].y =11.4	; */
/*  M[20].Xp[20].x = 	0.071	;	 M[20].Xp[20].y =12	; */
/*  M[20].Xp[21].x = 	0.065	;	 M[20].Xp[21].y =12.6	; */
/*  M[20].Xp[22].x = 	0.063	;	 M[20].Xp[22].y =13.2	; */
/*  M[20].Xp[23].x = 	0.051	;	 M[20].Xp[23].y =13.8	; */
/*  M[20].Xp[24].x = 	0.053	;	 M[20].Xp[24].y =14.4	; */
/*  M[20].Xp[25].x = 	0.049	;	 M[20].Xp[25].y =15	; */
/*  M[20].Xp[26].x = 	0.058	;	 M[20].Xp[26].y =15.6	; */
/*  M[20].Xp[27].x = 	0.046	;	 M[20].Xp[27].y =16.2	; */
/*  M[20].Xp[28].x = 	0.042	;	 M[20].Xp[28].y =16.8	; */
/*  M[20].Xp[29].x = 	0.042	;	 M[20].Xp[29].y =17.4	; */
/*  M[20].Xp[30].x = 	0.031	;	 M[20].Xp[30].y =18	; */
/*  M[20].Xp[31].x = 	0.035	;	 M[20].Xp[31].y =18.6	; */
/*  M[20].Xp[32].x = 	0.03	;	 M[20].Xp[32].y =19.2	; */
/*  M[20].Xp[33].x = 	0.021	;	 M[20].Xp[33].y =19.8	; */
/*  M[20].Xp[34].x = 	0.023	;	 M[20].Xp[34].y =20.4	; */
/*  M[20].Xp[35].x = 	0.022	;	 M[20].Xp[35].y =21	; */
/*  M[20].Xp[36].x = 	0.015	;	 M[20].Xp[36].y =21.6	; */
/*  M[20].Xp[37].x = 	0.01	;	 M[20].Xp[37].y =22.2	; */
/*  M[20].Xp[38].x = 	0.007	;	 M[20].Xp[38].y =22.8	; */
/*  M[20].Xp[39].x = 	0.003	;	 M[20].Xp[39].y =23.4	; */
/*  M[20].Xp[40].x = 	0.007	;	 M[20].Xp[40].y =24	; */
/*  M[20].Xp[41].x = 	0.006	;	 M[20].Xp[41].y =24.6	; */
/*  M[20].Xp[42].x = 	0.005	;	 M[20].Xp[42].y =25.2	; */
/*  M[20].Xp[43].x = 	0.003	;	 M[20].Xp[43].y =25.8	; */
/*  M[20].pathlength(); */
/* //station 22 */
/*  x[21] = 42.3; */
/*  M[21].allocate(44); */
/*  M[21].settype(SPLINE2D_LINEAR); */
/*  M[21].Xp[0].x = 	0.445	;	 M[21].Xp[0].y =0	; */
/*  M[21].Xp[1].x = 	0.466	;	 M[21].Xp[1].y =0.6	; */
/*  M[21].Xp[2].x = 	0.452	;	 M[21].Xp[2].y =1.2	; */
/*  M[21].Xp[3].x = 	0.432	;	 M[21].Xp[3].y =1.8	; */
/*  M[21].Xp[4].x = 	0.444	;	 M[21].Xp[4].y =2.4	; */
/*  M[21].Xp[5].x = 	0.386	;	 M[21].Xp[5].y =3	; */
/*  M[21].Xp[6].x = 	0.35	;	 M[21].Xp[6].y =3.6	; */
/*  M[21].Xp[7].x = 	0.337	;	 M[21].Xp[7].y =4.2	; */
/*  M[21].Xp[8].x = 	0.335	;	 M[21].Xp[8].y =4.8	; */
/*  M[21].Xp[9].x = 	0.293	;	 M[21].Xp[9].y =5.4	; */
/*  M[21].Xp[10].x = 	0.262	;	 M[21].Xp[10].y =6	; */
/*  M[21].Xp[11].x = 	0.244	;	 M[21].Xp[11].y =6.6	; */
/*  M[21].Xp[12].x = 	0.198	;	 M[21].Xp[12].y =7.2	; */
/*  M[21].Xp[13].x = 	0.189	;	 M[21].Xp[13].y =7.8	; */
/*  M[21].Xp[14].x = 	0.137	;	 M[21].Xp[14].y =8.4	; */
/*  M[21].Xp[15].x = 	0.134	;	 M[21].Xp[15].y =9	; */
/*  M[21].Xp[16].x = 	0.117	;	 M[21].Xp[16].y =9.6	; */
/*  M[21].Xp[17].x = 	0.112	;	 M[21].Xp[17].y =10.2	; */
/*  M[21].Xp[18].x = 	0.079	;	 M[21].Xp[18].y =10.8	; */
/*  M[21].Xp[19].x = 	0.068	;	 M[21].Xp[19].y =11.4	; */
/*  M[21].Xp[20].x = 	0.06	;	 M[21].Xp[20].y =12	; */
/*  M[21].Xp[21].x = 	0.06	;	 M[21].Xp[21].y =12.6	; */
/*  M[21].Xp[22].x = 	0.062	;	 M[21].Xp[22].y =13.2	; */
/*  M[21].Xp[23].x = 	0.061	;	 M[21].Xp[23].y =13.8	; */
/*  M[21].Xp[24].x = 	0.062	;	 M[21].Xp[24].y =14.4	; */
/*  M[21].Xp[25].x = 	0.048	;	 M[21].Xp[25].y =15	; */
/*  M[21].Xp[26].x = 	0.052	;	 M[21].Xp[26].y =15.6	; */
/*  M[21].Xp[27].x = 	0.054	;	 M[21].Xp[27].y =16.2	; */
/*  M[21].Xp[28].x = 	0.042	;	 M[21].Xp[28].y =16.8	; */
/*  M[21].Xp[29].x = 	0.038	;	 M[21].Xp[29].y =17.4	; */
/*  M[21].Xp[30].x = 	0.031	;	 M[21].Xp[30].y =18	; */
/*  M[21].Xp[31].x = 	0.03	;	 M[21].Xp[31].y =18.6	; */
/*  M[21].Xp[32].x = 	0.027	;	 M[21].Xp[32].y =19.2	; */
/*  M[21].Xp[33].x = 	0.026	;	 M[21].Xp[33].y =19.8	; */
/*  M[21].Xp[34].x = 	0.021	;	 M[21].Xp[34].y =20.4	; */
/*  M[21].Xp[35].x = 	0.019	;	 M[21].Xp[35].y =21	; */
/*  M[21].Xp[36].x = 	0.012	;	 M[21].Xp[36].y =21.6	; */
/*  M[21].Xp[37].x = 	0.012	;	 M[21].Xp[37].y =22.2	; */
/*  M[21].Xp[38].x = 	0.007	;	 M[21].Xp[38].y =22.8	; */
/*  M[21].Xp[39].x = 	0.004	;	 M[21].Xp[39].y =23.4	; */
/*  M[21].Xp[40].x = 	0.003	;	 M[21].Xp[40].y =24	; */
/*  M[21].Xp[41].x = 	0.004	;	 M[21].Xp[41].y =24.6	; */
/*  M[21].Xp[42].x = 	0.002	;	 M[21].Xp[42].y =25.2	; */
/*  M[21].Xp[43].x = 	0.001	;	 M[21].Xp[43].y =25.8	; */
/*  M[21].pathlength(); */
/* //station 23 */
/*  x[22] = 43.0; */
/*  M[22].allocate(44); */
/*  M[22].settype(SPLINE2D_LINEAR); */
/*  M[22].Xp[0].x = 	0.438	;	 M[22].Xp[0].y =0	; */
/*  M[22].Xp[1].x = 	0.466	;	 M[22].Xp[1].y =0.6	; */
/*  M[22].Xp[2].x = 	0.463	;	 M[22].Xp[2].y =1.2	; */
/*  M[22].Xp[3].x = 	0.445	;	 M[22].Xp[3].y =1.8	; */
/*  M[22].Xp[4].x = 	0.429	;	 M[22].Xp[4].y =2.4	; */
/*  M[22].Xp[5].x = 	0.41	;	 M[22].Xp[5].y =3	; */
/*  M[22].Xp[6].x = 	0.406	;	 M[22].Xp[6].y =3.6	; */
/*  M[22].Xp[7].x = 	0.356	;	 M[22].Xp[7].y =4.2	; */
/*  M[22].Xp[8].x = 	0.318	;	 M[22].Xp[8].y =4.8	; */
/*  M[22].Xp[9].x = 	0.307	;	 M[22].Xp[9].y =5.4	; */
/*  M[22].Xp[10].x = 	0.29	;	 M[22].Xp[10].y =6	; */
/*  M[22].Xp[11].x = 	0.236	;	 M[22].Xp[11].y =6.6	; */
/*  M[22].Xp[12].x = 	0.205	;	 M[22].Xp[12].y =7.2	; */
/*  M[22].Xp[13].x = 	0.17	;	 M[22].Xp[13].y =7.8	; */
/*  M[22].Xp[14].x = 	0.155	;	 M[22].Xp[14].y =8.4	; */
/*  M[22].Xp[15].x = 	0.117	;	 M[22].Xp[15].y =9	; */
/*  M[22].Xp[16].x = 	0.093	;	 M[22].Xp[16].y =9.6	; */
/*  M[22].Xp[17].x = 	0.103	;	 M[22].Xp[17].y =10.2	; */
/*  M[22].Xp[18].x = 	0.094	;	 M[22].Xp[18].y =10.8	; */
/*  M[22].Xp[19].x = 	0.084	;	 M[22].Xp[19].y =11.4	; */
/*  M[22].Xp[20].x = 	0.074	;	 M[22].Xp[20].y =12	; */
/*  M[22].Xp[21].x = 	0.07	;	 M[22].Xp[21].y =12.6	; */
/*  M[22].Xp[22].x = 	0.061	;	 M[22].Xp[22].y =13.2	; */
/*  M[22].Xp[23].x = 	0.06	;	 M[22].Xp[23].y =13.8	; */
/*  M[22].Xp[24].x = 	0.057	;	 M[22].Xp[24].y =14.4	; */
/*  M[22].Xp[25].x = 	0.048	;	 M[22].Xp[25].y =15	; */
/*  M[22].Xp[26].x = 	0.059	;	 M[22].Xp[26].y =15.6	; */
/*  M[22].Xp[27].x = 	0.047	;	 M[22].Xp[27].y =16.2	; */
/*  M[22].Xp[28].x = 	0.043	;	 M[22].Xp[28].y =16.8	; */
/*  M[22].Xp[29].x = 	0.038	;	 M[22].Xp[29].y =17.4	; */
/*  M[22].Xp[30].x = 	0.03	;	 M[22].Xp[30].y =18	; */
/*  M[22].Xp[31].x = 	0.028	;	 M[22].Xp[31].y =18.6	; */
/*  M[22].Xp[32].x = 	0.026	;	 M[22].Xp[32].y =19.2	; */
/*  M[22].Xp[33].x = 	0.028	;	 M[22].Xp[33].y =19.8	; */
/*  M[22].Xp[34].x = 	0.021	;	 M[22].Xp[34].y =20.4	; */
/*  M[22].Xp[35].x = 	0.014	;	 M[22].Xp[35].y =21	; */
/*  M[22].Xp[36].x = 	0.013	;	 M[22].Xp[36].y =21.6	; */
/*  M[22].Xp[37].x = 	0.006	;	 M[22].Xp[37].y =22.2	; */
/*  M[22].Xp[38].x = 	0.007	;	 M[22].Xp[38].y =22.8	; */
/*  M[22].Xp[39].x = 	0.003	;	 M[22].Xp[39].y =23.4	; */
/*  M[22].Xp[40].x = 	0.005	;	 M[22].Xp[40].y =24	; */
/*  M[22].Xp[41].x = 	0.004	;	 M[22].Xp[41].y =24.6	; */
/*  M[22].Xp[42].x = 	0	;	 M[22].Xp[42].y =25.2	; */
/*  M[22].Xp[43].x = 	0.002	;	 M[22].Xp[43].y =25.8	; */
/*  M[22].pathlength(); */
/* //station 24 */
/*  x[23] = 43.8; */
/*  M[23].allocate(44); */
/*  M[23].settype(SPLINE2D_LINEAR); */
/*  M[23].Xp[0].x = 	0.463	;	 M[23].Xp[0].y =0	; */
/*  M[23].Xp[1].x = 	0.461	;	 M[23].Xp[1].y =0.6	; */
/*  M[23].Xp[2].x = 	0.471	;	 M[23].Xp[2].y =1.2	; */
/*  M[23].Xp[3].x = 	0.475	;	 M[23].Xp[3].y =1.8	; */
/*  M[23].Xp[4].x = 	0.409	;	 M[23].Xp[4].y =2.4	; */
/*  M[23].Xp[5].x = 	0.391	;	 M[23].Xp[5].y =3	; */
/*  M[23].Xp[6].x = 	0.381	;	 M[23].Xp[6].y =3.6	; */
/*  M[23].Xp[7].x = 	0.346	;	 M[23].Xp[7].y =4.2	; */
/*  M[23].Xp[8].x = 	0.31	;	 M[23].Xp[8].y =4.8	; */
/*  M[23].Xp[9].x = 	0.298	;	 M[23].Xp[9].y =5.4	; */
/*  M[23].Xp[10].x = 	0.281	;	 M[23].Xp[10].y =6	; */
/*  M[23].Xp[11].x = 	0.238	;	 M[23].Xp[11].y =6.6	; */
/*  M[23].Xp[12].x = 	0.225	;	 M[23].Xp[12].y =7.2	; */
/*  M[23].Xp[13].x = 	0.18	;	 M[23].Xp[13].y =7.8	; */
/*  M[23].Xp[14].x = 	0.163	;	 M[23].Xp[14].y =8.4	; */
/*  M[23].Xp[15].x = 	0.14	;	 M[23].Xp[15].y =9	; */
/*  M[23].Xp[16].x = 	0.109	;	 M[23].Xp[16].y =9.6	; */
/*  M[23].Xp[17].x = 	0.1	;	 M[23].Xp[17].y =10.2	; */
/*  M[23].Xp[18].x = 	0.091	;	 M[23].Xp[18].y =10.8	; */
/*  M[23].Xp[19].x = 	0.077	;	 M[23].Xp[19].y =11.4	; */
/*  M[23].Xp[20].x = 	0.071	;	 M[23].Xp[20].y =12	; */
/*  M[23].Xp[21].x = 	0.064	;	 M[23].Xp[21].y =12.6	; */
/*  M[23].Xp[22].x = 	0.06	;	 M[23].Xp[22].y =13.2	; */
/*  M[23].Xp[23].x = 	0.052	;	 M[23].Xp[23].y =13.8	; */
/*  M[23].Xp[24].x = 	0.05	;	 M[23].Xp[24].y =14.4	; */
/*  M[23].Xp[25].x = 	0.044	;	 M[23].Xp[25].y =15	; */
/*  M[23].Xp[26].x = 	0.043	;	 M[23].Xp[26].y =15.6	; */
/*  M[23].Xp[27].x = 	0.042	;	 M[23].Xp[27].y =16.2	; */
/*  M[23].Xp[28].x = 	0.042	;	 M[23].Xp[28].y =16.8	; */
/*  M[23].Xp[29].x = 	0.039	;	 M[23].Xp[29].y =17.4	; */
/*  M[23].Xp[30].x = 	0.027	;	 M[23].Xp[30].y =18	; */
/*  M[23].Xp[31].x = 	0.032	;	 M[23].Xp[31].y =18.6	; */
/*  M[23].Xp[32].x = 	0.02	;	 M[23].Xp[32].y =19.2	; */
/*  M[23].Xp[33].x = 	0.023	;	 M[23].Xp[33].y =19.8	; */
/*  M[23].Xp[34].x = 	0.017	;	 M[23].Xp[34].y =20.4	; */
/*  M[23].Xp[35].x = 	0.018	;	 M[23].Xp[35].y =21	; */
/*  M[23].Xp[36].x = 	0.009	;	 M[23].Xp[36].y =21.6	; */
/*  M[23].Xp[37].x = 	0.002	;	 M[23].Xp[37].y =22.2	; */
/*  M[23].Xp[38].x = 	0.004	;	 M[23].Xp[38].y =22.8	; */
/*  M[23].Xp[39].x = 	0.005	;	 M[23].Xp[39].y =23.4	; */
/*  M[23].Xp[40].x = 	0.003	;	 M[23].Xp[40].y =24	; */
/*  M[23].Xp[41].x = 	0.001	;	 M[23].Xp[41].y =24.6	; */
/*  M[23].Xp[42].x = 	0.002	;	 M[23].Xp[42].y =25.2	; */
/*  M[23].Xp[43].x = 	0.003	;	 M[23].Xp[43].y =25.8	; */
/*  M[23].pathlength(); */
/* //station 25 */
/*   x[24] = 44.5; */
/*  M[24].allocate(44); */
/*  M[24].settype(SPLINE2D_LINEAR); */
/*  M[24].Xp[0].x = 	0.405	;	 M[24].Xp[0].y =0	; */
/*  M[24].Xp[1].x = 	0.425	;	 M[24].Xp[1].y =0.6	; */
/*  M[24].Xp[2].x = 	0.441	;	 M[24].Xp[2].y =1.2	; */
/*  M[24].Xp[3].x = 	0.452	;	 M[24].Xp[3].y =1.8	; */
/*  M[24].Xp[4].x = 	0.438	;	 M[24].Xp[4].y =2.4	; */
/*  M[24].Xp[5].x = 	0.368	;	 M[24].Xp[5].y =3	; */
/*  M[24].Xp[6].x = 	0.351	;	 M[24].Xp[6].y =3.6	; */
/*  M[24].Xp[7].x = 	0.344	;	 M[24].Xp[7].y =4.2	; */
/*  M[24].Xp[8].x = 	0.328	;	 M[24].Xp[8].y =4.8	; */
/*  M[24].Xp[9].x = 	0.316	;	 M[24].Xp[9].y =5.4	; */
/*  M[24].Xp[10].x = 	0.293	;	 M[24].Xp[10].y =6	; */
/*  M[24].Xp[11].x = 	0.254	;	 M[24].Xp[11].y =6.6	; */
/*  M[24].Xp[12].x = 	0.212	;	 M[24].Xp[12].y =7.2	; */
/*  M[24].Xp[13].x = 	0.202	;	 M[24].Xp[13].y =7.8	; */
/*  M[24].Xp[14].x = 	0.18	;	 M[24].Xp[14].y =8.4	; */
/*  M[24].Xp[15].x = 	0.143	;	 M[24].Xp[15].y =9	; */
/*  M[24].Xp[16].x = 	0.105	;	 M[24].Xp[16].y =9.6	; */
/*  M[24].Xp[17].x = 	0.109	;	 M[24].Xp[17].y =10.2	; */
/*  M[24].Xp[18].x = 	0.092	;	 M[24].Xp[18].y =10.8	; */
/*  M[24].Xp[19].x = 	0.073	;	 M[24].Xp[19].y =11.4	; */
/*  M[24].Xp[20].x = 	0.076	;	 M[24].Xp[20].y =12	; */
/*  M[24].Xp[21].x = 	0.053	;	 M[24].Xp[21].y =12.6	; */
/*  M[24].Xp[22].x = 	0.056	;	 M[24].Xp[22].y =13.2	; */
/*  M[24].Xp[23].x = 	0.048	;	 M[24].Xp[23].y =13.8	; */
/*  M[24].Xp[24].x = 	0.05	;	 M[24].Xp[24].y =14.4	; */
/*  M[24].Xp[25].x = 	0.048	;	 M[24].Xp[25].y =15	; */
/*  M[24].Xp[26].x = 	0.043	;	 M[24].Xp[26].y =15.6	; */
/*  M[24].Xp[27].x = 	0.041	;	 M[24].Xp[27].y =16.2	; */
/*  M[24].Xp[28].x = 	0.039	;	 M[24].Xp[28].y =16.8	; */
/*  M[24].Xp[29].x = 	0.037	;	 M[24].Xp[29].y =17.4	; */
/*  M[24].Xp[30].x = 	0.035	;	 M[24].Xp[30].y =18	; */
/*  M[24].Xp[31].x = 	0.029	;	 M[24].Xp[31].y =18.6	; */
/*  M[24].Xp[32].x = 	0.025	;	 M[24].Xp[32].y =19.2	; */
/*  M[24].Xp[33].x = 	0.023	;	 M[24].Xp[33].y =19.8	; */
/*  M[24].Xp[34].x = 	0.015	;	 M[24].Xp[34].y =20.4	; */
/*  M[24].Xp[35].x = 	0.013	;	 M[24].Xp[35].y =21	; */
/*  M[24].Xp[36].x = 	0.009	;	 M[24].Xp[36].y =21.6	; */
/*  M[24].Xp[37].x = 	0.009	;	 M[24].Xp[37].y =22.2	; */
/*  M[24].Xp[38].x = 	0.004	;	 M[24].Xp[38].y =22.8	; */
/*  M[24].Xp[39].x = 	0.002	;	 M[24].Xp[39].y =23.4	; */
/*  M[24].Xp[40].x = 	0.002	;	 M[24].Xp[40].y =24	; */
/*  M[24].Xp[41].x = 	0.002	;	 M[24].Xp[41].y =24.6	; */
/*  M[24].Xp[42].x = 	0.002	;	 M[24].Xp[42].y =25.2	; */
/*  M[24].Xp[43].x = 	0.005	;	 M[24].Xp[43].y =25.8	; */
/*  M[24].pathlength(); */

/* //station 26 */
/*  x[25] = 45.3; */
/*  M[25].allocate(44); */
/*  M[25].settype(SPLINE2D_LINEAR); */
/*  M[25].Xp[0].x = 	0.398	;	 M[25].Xp[0].y =0	; */
/*  M[25].Xp[1].x = 	0.415	;	 M[25].Xp[1].y =0.6	; */
/*  M[25].Xp[2].x = 	0.457	;	 M[25].Xp[2].y =1.2	; */
/*  M[25].Xp[3].x = 	0.461	;	 M[25].Xp[3].y =1.8	; */
/*  M[25].Xp[4].x = 	0.408	;	 M[25].Xp[4].y =2.4	; */
/*  M[25].Xp[5].x = 	0.385	;	 M[25].Xp[5].y =3	; */
/*  M[25].Xp[6].x = 	0.381	;	 M[25].Xp[6].y =3.6	; */
/*  M[25].Xp[7].x = 	0.341	;	 M[25].Xp[7].y =4.2	; */
/*  M[25].Xp[8].x = 	0.338	;	 M[25].Xp[8].y =4.8	; */
/*  M[25].Xp[9].x = 	0.319	;	 M[25].Xp[9].y =5.4	; */
/*  M[25].Xp[10].x = 	0.283	;	 M[25].Xp[10].y =6	; */
/*  M[25].Xp[11].x = 	0.247	;	 M[25].Xp[11].y =6.6	; */
/*  M[25].Xp[12].x = 	0.194	;	 M[25].Xp[12].y =7.2	; */
/*  M[25].Xp[13].x = 	0.185	;	 M[25].Xp[13].y =7.8	; */
/*  M[25].Xp[14].x = 	0.154	;	 M[25].Xp[14].y =8.4	; */
/*  M[25].Xp[15].x = 	0.143	;	 M[25].Xp[15].y =9	; */
/*  M[25].Xp[16].x = 	0.126	;	 M[25].Xp[16].y =9.6	; */
/*  M[25].Xp[17].x = 	0.106	;	 M[25].Xp[17].y =10.2	; */
/*  M[25].Xp[18].x = 	0.092	;	 M[25].Xp[18].y =10.8	; */
/*  M[25].Xp[19].x = 	0.082	;	 M[25].Xp[19].y =11.4	; */
/*  M[25].Xp[20].x = 	0.061	;	 M[25].Xp[20].y =12	; */
/*  M[25].Xp[21].x = 	0.057	;	 M[25].Xp[21].y =12.6	; */
/*  M[25].Xp[22].x = 	0.055	;	 M[25].Xp[22].y =13.2	; */
/*  M[25].Xp[23].x = 	0.046	;	 M[25].Xp[23].y =13.8	; */
/*  M[25].Xp[24].x = 	0.052	;	 M[25].Xp[24].y =14.4	; */
/*  M[25].Xp[25].x = 	0.047	;	 M[25].Xp[25].y =15	; */
/*  M[25].Xp[26].x = 	0.044	;	 M[25].Xp[26].y =15.6	; */
/*  M[25].Xp[27].x = 	0.052	;	 M[25].Xp[27].y =16.2	; */
/*  M[25].Xp[28].x = 	0.043	;	 M[25].Xp[28].y =16.8	; */
/*  M[25].Xp[29].x = 	0.045	;	 M[25].Xp[29].y =17.4	; */
/*  M[25].Xp[30].x = 	0.033	;	 M[25].Xp[30].y =18	; */
/*  M[25].Xp[31].x = 	0.027	;	 M[25].Xp[31].y =18.6	; */
/*  M[25].Xp[32].x = 	0.025	;	 M[25].Xp[32].y =19.2	; */
/*  M[25].Xp[33].x = 	0.022	;	 M[25].Xp[33].y =19.8	; */
/*  M[25].Xp[34].x = 	0.016	;	 M[25].Xp[34].y =20.4	; */
/*  M[25].Xp[35].x = 	0.015	;	 M[25].Xp[35].y =21	; */
/*  M[25].Xp[36].x = 	0.011	;	 M[25].Xp[36].y =21.6	; */
/*  M[25].Xp[37].x = 	0.006	;	 M[25].Xp[37].y =22.2	; */
/*  M[25].Xp[38].x = 	0.006	;	 M[25].Xp[38].y =22.8	; */
/*  M[25].Xp[39].x = 	0.006	;	 M[25].Xp[39].y =23.4	; */
/*  M[25].Xp[40].x = 	0.003	;	 M[25].Xp[40].y =24	; */
/*  M[25].Xp[41].x = 	0.003	;	 M[25].Xp[41].y =24.6	; */
/*  M[25].Xp[42].x = 	0.006	;	 M[25].Xp[42].y =25.2	; */
/*  M[25].Xp[43].x = 	0.004	;	 M[25].Xp[43].y =25.8	; */
/*  M[25].pathlength(); */
/* //station 27 */
/*  x[26] = 46.0; */
/*  M[26].allocate(44); */
/*  M[26].settype(SPLINE2D_LINEAR); */
/*  M[26].Xp[0].x = 	0.431	;	 M[26].Xp[0].y =0	; */
/*  M[26].Xp[1].x = 	0.409	;	 M[26].Xp[1].y =0.6	; */
/*  M[26].Xp[2].x = 	0.416	;	 M[26].Xp[2].y =1.2	; */
/*  M[26].Xp[3].x = 	0.421	;	 M[26].Xp[3].y =1.8	; */
/*  M[26].Xp[4].x = 	0.405	;	 M[26].Xp[4].y =2.4	; */
/*  M[26].Xp[5].x = 	0.385	;	 M[26].Xp[5].y =3	; */
/*  M[26].Xp[6].x = 	0.376	;	 M[26].Xp[6].y =3.6	; */
/*  M[26].Xp[7].x = 	0.352	;	 M[26].Xp[7].y =4.2	; */
/*  M[26].Xp[8].x = 	0.337	;	 M[26].Xp[8].y =4.8	; */
/*  M[26].Xp[9].x = 	0.307	;	 M[26].Xp[9].y =5.4	; */
/*  M[26].Xp[10].x = 	0.272	;	 M[26].Xp[10].y =6	; */
/*  M[26].Xp[11].x = 	0.238	;	 M[26].Xp[11].y =6.6	; */
/*  M[26].Xp[12].x = 	0.187	;	 M[26].Xp[12].y =7.2	; */
/*  M[26].Xp[13].x = 	0.182	;	 M[26].Xp[13].y =7.8	; */
/*  M[26].Xp[14].x = 	0.164	;	 M[26].Xp[14].y =8.4	; */
/*  M[26].Xp[15].x = 	0.136	;	 M[26].Xp[15].y =9	; */
/*  M[26].Xp[16].x = 	0.116	;	 M[26].Xp[16].y =9.6	; */
/*  M[26].Xp[17].x = 	0.113	;	 M[26].Xp[17].y =10.2	; */
/*  M[26].Xp[18].x = 	0.094	;	 M[26].Xp[18].y =10.8	; */
/*  M[26].Xp[19].x = 	0.074	;	 M[26].Xp[19].y =11.4	; */
/*  M[26].Xp[20].x = 	0.071	;	 M[26].Xp[20].y =12	; */
/*  M[26].Xp[21].x = 	0.064	;	 M[26].Xp[21].y =12.6	; */
/*  M[26].Xp[22].x = 	0.053	;	 M[26].Xp[22].y =13.2	; */
/*  M[26].Xp[23].x = 	0.058	;	 M[26].Xp[23].y =13.8	; */
/*  M[26].Xp[24].x = 	0.053	;	 M[26].Xp[24].y =14.4	; */
/*  M[26].Xp[25].x = 	0.047	;	 M[26].Xp[25].y =15	; */
/*  M[26].Xp[26].x = 	0.044	;	 M[26].Xp[26].y =15.6	; */
/*  M[26].Xp[27].x = 	0.036	;	 M[26].Xp[27].y =16.2	; */
/*  M[26].Xp[28].x = 	0.038	;	 M[26].Xp[28].y =16.8	; */
/*  M[26].Xp[29].x = 	0.033	;	 M[26].Xp[29].y =17.4	; */
/*  M[26].Xp[30].x = 	0.034	;	 M[26].Xp[30].y =18	; */
/*  M[26].Xp[31].x = 	0.027	;	 M[26].Xp[31].y =18.6	; */
/*  M[26].Xp[32].x = 	0.02	;	 M[26].Xp[32].y =19.2	; */
/*  M[26].Xp[33].x = 	0.023	;	 M[26].Xp[33].y =19.8	; */
/*  M[26].Xp[34].x = 	0.014	;	 M[26].Xp[34].y =20.4	; */
/*  M[26].Xp[35].x = 	0.016	;	 M[26].Xp[35].y =21	; */
/*  M[26].Xp[36].x = 	0.01	;	 M[26].Xp[36].y =21.6	; */
/*  M[26].Xp[37].x = 	0.009	;	 M[26].Xp[37].y =22.2	; */
/*  M[26].Xp[38].x = 	0.004	;	 M[26].Xp[38].y =22.8	; */
/*  M[26].Xp[39].x = 	0.004	;	 M[26].Xp[39].y =23.4	; */
/*  M[26].Xp[40].x = 	0.006	;	 M[26].Xp[40].y =24	; */
/*  M[26].Xp[41].x = 	0.007	;	 M[26].Xp[41].y =24.6	; */
/*  M[26].Xp[42].x = 	0.004	;	 M[26].Xp[42].y =25.2	; */
/*  M[26].Xp[43].x = 	0.006	;	 M[26].Xp[43].y =25.8	; */

/*  M[26].pathlength(); */
/* //station 28 */
/*  x[27] = 46.8; */
/*  M[27].allocate(44); */
/*  M[27].settype(SPLINE2D_LINEAR); */
/*  M[27].Xp[0].x = 	0.425	;	 M[27].Xp[0].y =0	; */
/*  M[27].Xp[1].x = 	0.411	;	 M[27].Xp[1].y =0.6	; */
/*  M[27].Xp[2].x = 	0.41	;	 M[27].Xp[2].y =1.2	; */
/*  M[27].Xp[3].x = 	0.366	;	 M[27].Xp[3].y =1.8	; */
/*  M[27].Xp[4].x = 	0.381	;	 M[27].Xp[4].y =2.4	; */
/*  M[27].Xp[5].x = 	0.364	;	 M[27].Xp[5].y =3	; */
/*  M[27].Xp[6].x = 	0.355	;	 M[27].Xp[6].y =3.6	; */
/*  M[27].Xp[7].x = 	0.336	;	 M[27].Xp[7].y =4.2	; */
/*  M[27].Xp[8].x = 	0.316	;	 M[27].Xp[8].y =4.8	; */
/*  M[27].Xp[9].x = 	0.298	;	 M[27].Xp[9].y =5.4	; */
/*  M[27].Xp[10].x = 	0.248	;	 M[27].Xp[10].y =6	; */
/*  M[27].Xp[11].x = 	0.229	;	 M[27].Xp[11].y =6.6	; */
/*  M[27].Xp[12].x = 	0.193	;	 M[27].Xp[12].y =7.2	; */
/*  M[27].Xp[13].x = 	0.181	;	 M[27].Xp[13].y =7.8	; */
/*  M[27].Xp[14].x = 	0.155	;	 M[27].Xp[14].y =8.4	; */
/*  M[27].Xp[15].x = 	0.123	;	 M[27].Xp[15].y =9	; */
/*  M[27].Xp[16].x = 	0.106	;	 M[27].Xp[16].y =9.6	; */
/*  M[27].Xp[17].x = 	0.089	;	 M[27].Xp[17].y =10.2	; */
/*  M[27].Xp[18].x = 	0.084	;	 M[27].Xp[18].y =10.8	; */
/*  M[27].Xp[19].x = 	0.074	;	 M[27].Xp[19].y =11.4	; */
/*  M[27].Xp[20].x = 	0.061	;	 M[27].Xp[20].y =12	; */
/*  M[27].Xp[21].x = 	0.062	;	 M[27].Xp[21].y =12.6	; */
/*  M[27].Xp[22].x = 	0.059	;	 M[27].Xp[22].y =13.2	; */
/*  M[27].Xp[23].x = 	0.046	;	 M[27].Xp[23].y =13.8	; */
/*  M[27].Xp[24].x = 	0.049	;	 M[27].Xp[24].y =14.4	; */
/*  M[27].Xp[25].x = 	0.05	;	 M[27].Xp[25].y =15	; */
/*  M[27].Xp[26].x = 	0.042	;	 M[27].Xp[26].y =15.6	; */
/*  M[27].Xp[27].x = 	0.039	;	 M[27].Xp[27].y =16.2	; */
/*  M[27].Xp[28].x = 	0.032	;	 M[27].Xp[28].y =16.8	; */
/*  M[27].Xp[29].x = 	0.031	;	 M[27].Xp[29].y =17.4	; */
/*  M[27].Xp[30].x = 	0.029	;	 M[27].Xp[30].y =18	; */
/*  M[27].Xp[31].x = 	0.024	;	 M[27].Xp[31].y =18.6	; */
/*  M[27].Xp[32].x = 	0.017	;	 M[27].Xp[32].y =19.2	; */
/*  M[27].Xp[33].x = 	0.019	;	 M[27].Xp[33].y =19.8	; */
/*  M[27].Xp[34].x = 	0.013	;	 M[27].Xp[34].y =20.4	; */
/*  M[27].Xp[35].x = 	0.013	;	 M[27].Xp[35].y =21	; */
/*  M[27].Xp[36].x = 	0.011	;	 M[27].Xp[36].y =21.6	; */
/*  M[27].Xp[37].x = 	0.005	;	 M[27].Xp[37].y =22.2	; */
/*  M[27].Xp[38].x = 	0.002	;	 M[27].Xp[38].y =22.8	; */
/*  M[27].Xp[39].x = 	0.01	;	 M[27].Xp[39].y =23.4	; */
/*  M[27].Xp[40].x = 	0.005	;	 M[27].Xp[40].y =24	; */
/*  M[27].Xp[41].x = 	0.002	;	 M[27].Xp[41].y =24.6	; */
/*  M[27].Xp[42].x = 	0.007	;	 M[27].Xp[42].y =25.2	; */
/*  M[27].Xp[43].x = 	0.002	;	 M[27].Xp[43].y =25.8	; */

/*  M[27].pathlength(); */
/* //station 29 */
/*  x[28] = 47.5; */
/*  M[28].allocate(44); */
/*  M[28].settype(SPLINE2D_LINEAR); */
/*  M[28].Xp[0].x = 	0.408	;	 M[28].Xp[0].y =0	; */
/*  M[28].Xp[1].x = 	0.372	;	 M[28].Xp[1].y =0.6	; */
/*  M[28].Xp[2].x = 	0.385	;	 M[28].Xp[2].y =1.2	; */
/*  M[28].Xp[3].x = 	0.379	;	 M[28].Xp[3].y =1.8	; */
/*  M[28].Xp[4].x = 	0.371	;	 M[28].Xp[4].y =2.4	; */
/*  M[28].Xp[5].x = 	0.346	;	 M[28].Xp[5].y =3	; */
/*  M[28].Xp[6].x = 	0.34	;	 M[28].Xp[6].y =3.6	; */
/*  M[28].Xp[7].x = 	0.302	;	 M[28].Xp[7].y =4.2	; */
/*  M[28].Xp[8].x = 	0.324	;	 M[28].Xp[8].y =4.8	; */
/*  M[28].Xp[9].x = 	0.263	;	 M[28].Xp[9].y =5.4	; */
/*  M[28].Xp[10].x = 	0.238	;	 M[28].Xp[10].y =6	; */
/*  M[28].Xp[11].x = 	0.191	;	 M[28].Xp[11].y =6.6	; */
/*  M[28].Xp[12].x = 	0.178	;	 M[28].Xp[12].y =7.2	; */
/*  M[28].Xp[13].x = 	0.181	;	 M[28].Xp[13].y =7.8	; */
/*  M[28].Xp[14].x = 	0.155	;	 M[28].Xp[14].y =8.4	; */
/*  M[28].Xp[15].x = 	0.128	;	 M[28].Xp[15].y =9	; */
/*  M[28].Xp[16].x = 	0.102	;	 M[28].Xp[16].y =9.6	; */
/*  M[28].Xp[17].x = 	0.089	;	 M[28].Xp[17].y =10.2	; */
/*  M[28].Xp[18].x = 	0.078	;	 M[28].Xp[18].y =10.8	; */
/*  M[28].Xp[19].x = 	0.08	;	 M[28].Xp[19].y =11.4	; */
/*  M[28].Xp[20].x = 	0.059	;	 M[28].Xp[20].y =12	; */
/*  M[28].Xp[21].x = 	0.056	;	 M[28].Xp[21].y =12.6	; */
/*  M[28].Xp[22].x = 	0.042	;	 M[28].Xp[22].y =13.2	; */
/*  M[28].Xp[23].x = 	0.044	;	 M[28].Xp[23].y =13.8	; */
/*  M[28].Xp[24].x = 	0.042	;	 M[28].Xp[24].y =14.4	; */
/*  M[28].Xp[25].x = 	0.048	;	 M[28].Xp[25].y =15	; */
/*  M[28].Xp[26].x = 	0.044	;	 M[28].Xp[26].y =15.6	; */
/*  M[28].Xp[27].x = 	0.031	;	 M[28].Xp[27].y =16.2	; */
/*  M[28].Xp[28].x = 	0.032	;	 M[28].Xp[28].y =16.8	; */
/*  M[28].Xp[29].x = 	0.03	;	 M[28].Xp[29].y =17.4	; */
/*  M[28].Xp[30].x = 	0.026	;	 M[28].Xp[30].y =18	; */
/*  M[28].Xp[31].x = 	0.022	;	 M[28].Xp[31].y =18.6	; */
/*  M[28].Xp[32].x = 	0.021	;	 M[28].Xp[32].y =19.2	; */
/*  M[28].Xp[33].x = 	0.016	;	 M[28].Xp[33].y =19.8	; */
/*  M[28].Xp[34].x = 	0.01	;	 M[28].Xp[34].y =20.4	; */
/*  M[28].Xp[35].x = 	0.013	;	 M[28].Xp[35].y =21	; */
/*  M[28].Xp[36].x = 	0.008	;	 M[28].Xp[36].y =21.6	; */
/*  M[28].Xp[37].x = 	0.006	;	 M[28].Xp[37].y =22.2	; */
/*  M[28].Xp[38].x = 	0.004	;	 M[28].Xp[38].y =22.8	; */
/*  M[28].Xp[39].x = 	0.003	;	 M[28].Xp[39].y =23.4	; */
/*  M[28].Xp[40].x = 	0.003	;	 M[28].Xp[40].y =24	; */
/*  M[28].Xp[41].x = 	0.003	;	 M[28].Xp[41].y =24.6	; */
/*  M[28].Xp[42].x = 	0.002	;	 M[28].Xp[42].y =25.2	; */
/*  M[28].Xp[43].x = 	0.006	;	 M[28].Xp[43].y =25.8	; */

/*  M[28].pathlength(); */
/* //station 30 */
/*  x[29] = 48.3; */
/*  M[29].allocate(44); */
/*  M[29].settype(SPLINE2D_LINEAR); */
/*  M[29].Xp[0].x = 	0.39	;	 M[29].Xp[0].y =0	; */
/*  M[29].Xp[1].x = 	0.38	;	 M[29].Xp[1].y =0.6	; */
/*  M[29].Xp[2].x = 	0.379	;	 M[29].Xp[2].y =1.2	; */
/*  M[29].Xp[3].x = 	0.382	;	 M[29].Xp[3].y =1.8	; */
/*  M[29].Xp[4].x = 	0.347	;	 M[29].Xp[4].y =2.4	; */
/*  M[29].Xp[5].x = 	0.335	;	 M[29].Xp[5].y =3	; */
/*  M[29].Xp[6].x = 	0.325	;	 M[29].Xp[6].y =3.6	; */
/*  M[29].Xp[7].x = 	0.296	;        M[29].Xp[7].y =4.2	; */
/*  M[29].Xp[8].x = 	0.279	;	 M[29].Xp[8].y =4.8	; */
/*  M[29].Xp[9].x = 	0.27	;	 M[29].Xp[9].y =5.4	; */
/*  M[29].Xp[10].x = 	0.256	;	 M[29].Xp[10].y =6	; */
/*  M[29].Xp[11].x = 	0.19	;	 M[29].Xp[11].y =6.6	; */
/*  M[29].Xp[12].x = 	0.182	;	 M[29].Xp[12].y =7.2	; */
/*  M[29].Xp[13].x = 	0.16	;	 M[29].Xp[13].y =7.8	; */
/*  M[29].Xp[14].x = 	0.158	;	 M[29].Xp[14].y =8.4	; */
/*  M[29].Xp[15].x = 	0.116	;	 M[29].Xp[15].y =9	; */
/*  M[29].Xp[16].x = 	0.092	;	 M[29].Xp[16].y =9.6	; */
/*  M[29].Xp[17].x = 	0.092	;	 M[29].Xp[17].y =10.2	; */
/*  M[29].Xp[18].x = 	0.069	;	 M[29].Xp[18].y =10.8	; */
/*  M[29].Xp[19].x = 	0.067	;	 M[29].Xp[19].y =11.4	; */
/*  M[29].Xp[20].x = 	0.063	;	 M[29].Xp[20].y =12	; */
/*  M[29].Xp[21].x = 	0.058	;	 M[29].Xp[21].y =12.6	; */
/*  M[29].Xp[22].x = 	0.05	;	 M[29].Xp[22].y =13.2	; */
/*  M[29].Xp[23].x = 	0.046	;	 M[29].Xp[23].y =13.8	; */
/*  M[29].Xp[24].x = 	0.045	;	 M[29].Xp[24].y =14.4	; */
/*  M[29].Xp[25].x = 	0.043	;	 M[29].Xp[25].y =15	; */
/*  M[29].Xp[26].x = 	0.045	;	 M[29].Xp[26].y =15.6	; */
/*  M[29].Xp[27].x = 	0.031	;	 M[29].Xp[27].y =16.2	; */
/*  M[29].Xp[28].x = 	0.034	;	 M[29].Xp[28].y =16.8	; */
/*  M[29].Xp[29].x = 	0.027	;	 M[29].Xp[29].y =17.4	; */
/*  M[29].Xp[30].x = 	0.022	;	 M[29].Xp[30].y =18	; */
/*  M[29].Xp[31].x = 	0.022	;	 M[29].Xp[31].y =18.6	; */
/*  M[29].Xp[32].x = 	0.015	;	 M[29].Xp[32].y =19.2	; */
/*  M[29].Xp[33].x = 	0.017	;	 M[29].Xp[33].y =19.8	; */
/*  M[29].Xp[34].x = 	0.011	;	 M[29].Xp[34].y =20.4	; */
/*  M[29].Xp[35].x = 	0.014	;	 M[29].Xp[35].y =21	; */
/*  M[29].Xp[36].x = 	0.006	;	 M[29].Xp[36].y =21.6	; */
/*  M[29].Xp[37].x = 	0.009	;	 M[29].Xp[37].y =22.2	; */
/*  M[29].Xp[38].x = 	0.009	;	 M[29].Xp[38].y =22.8	; */
/*  M[29].Xp[39].x = 	0.004	;	 M[29].Xp[39].y =23.4	; */
/*  M[29].Xp[40].x = 	0.003	;	 M[29].Xp[40].y =24	; */
/*  M[29].Xp[41].x = 	0.005	;	 M[29].Xp[41].y =24.6	; */
/*  M[29].Xp[42].x = 	0.006	;	 M[29].Xp[42].y =25.2	; */
/*  M[29].Xp[43].x = 	0.004	;	 M[29].Xp[43].y =25.8	; */
/*  M[29].pathlength(); */

/* //station 31 */
/*  x[30] = 49.0; */
/*  M[30].allocate(44); */
/*  M[30].settype(SPLINE2D_LINEAR); */
/*  M[30].Xp[0].x = 	0.406	;	 M[30].Xp[0].y =0	; */
/*  M[30].Xp[1].x = 	0.387	;	 M[30].Xp[1].y =0.6	; */
/*  M[30].Xp[2].x = 	0.366	;	 M[30].Xp[2].y =1.2	; */
/*  M[30].Xp[3].x = 	0.354	;	 M[30].Xp[3].y =1.8	; */
/*  M[30].Xp[4].x = 	0.337	;	 M[30].Xp[4].y =2.4	; */
/*  M[30].Xp[5].x = 	0.301	;	 M[30].Xp[5].y =3	; */
/*  M[30].Xp[6].x = 	0.288	;	 M[30].Xp[6].y =3.6	; */
/*  M[30].Xp[7].x = 	0.279	;	 M[30].Xp[7].y =4.2	; */
/*  M[30].Xp[8].x = 	0.272	;	 M[30].Xp[8].y =4.8	; */
/*  M[30].Xp[9].x = 	0.244	;	 M[30].Xp[9].y =5.4	; */
/*  M[30].Xp[10].x = 	0.206	;	 M[30].Xp[10].y =6	; */
/*  M[30].Xp[11].x = 	0.2	;	 M[30].Xp[11].y =6.6	; */
/*  M[30].Xp[12].x = 	0.18	;	 M[30].Xp[12].y =7.2	; */
/*  M[30].Xp[13].x = 	0.145	;	 M[30].Xp[13].y =7.8	; */
/*  M[30].Xp[14].x = 	0.132	;	 M[30].Xp[14].y =8.4	; */
/*  M[30].Xp[15].x = 	0.103	;	 M[30].Xp[15].y =9	; */
/*  M[30].Xp[16].x = 	0.095	;	 M[30].Xp[16].y =9.6	; */
/*  M[30].Xp[17].x = 	0.085	;	 M[30].Xp[17].y =10.2	; */
/*  M[30].Xp[18].x = 	0.073	;	 M[30].Xp[18].y =10.8	; */
/*  M[30].Xp[19].x = 	0.067	;	 M[30].Xp[19].y =11.4	; */
/*  M[30].Xp[20].x = 	0.051	;	 M[30].Xp[20].y =12	; */
/*  M[30].Xp[21].x = 	0.047	;	 M[30].Xp[21].y =12.6	; */
/*  M[30].Xp[22].x = 	0.053	;	 M[30].Xp[22].y =13.2	; */
/*  M[30].Xp[23].x = 	0.038	;	 M[30].Xp[23].y =13.8	; */
/*  M[30].Xp[24].x = 	0.044	;	 M[30].Xp[24].y =14.4	; */
/*  M[30].Xp[25].x = 	0.04	;	 M[30].Xp[25].y =15	; */
/*  M[30].Xp[26].x = 	0.037	;	 M[30].Xp[26].y =15.6	; */
/*  M[30].Xp[27].x = 	0.03	;	 M[30].Xp[27].y =16.2	; */
/*  M[30].Xp[28].x = 	0.022	;	 M[30].Xp[28].y =16.8	; */
/*  M[30].Xp[29].x = 	0.027	;	 M[30].Xp[29].y =17.4	; */
/*  M[30].Xp[30].x = 	0.028	;	 M[30].Xp[30].y =18	; */
/*  M[30].Xp[31].x = 	0.023	;	 M[30].Xp[31].y =18.6	; */
/*  M[30].Xp[32].x = 	0.024	;	 M[30].Xp[32].y =19.2	; */
/*  M[30].Xp[33].x = 	0.018	;	 M[30].Xp[33].y =19.8	; */
/*  M[30].Xp[34].x = 	0.01	;	 M[30].Xp[34].y =20.4	; */
/*  M[30].Xp[35].x = 	0.006	;	 M[30].Xp[35].y =21	; */
/*  M[30].Xp[36].x = 	0.008	;	 M[30].Xp[36].y =21.6	; */
/*  M[30].Xp[37].x = 	0.004	;	 M[30].Xp[37].y =22.2	; */
/*  M[30].Xp[38].x = 	0.002	;	 M[30].Xp[38].y =22.8	; */
/*  M[30].Xp[39].x = 	0.003	;	 M[30].Xp[39].y =23.4	; */
/*  M[30].Xp[40].x = 	0.002	;	 M[30].Xp[40].y =24	; */
/*  M[30].Xp[41].x = 	0.002	;	 M[30].Xp[41].y =24.6	; */
/*  M[30].Xp[42].x = 	0.001	;	 M[30].Xp[42].y =25.2	; */
/*  M[30].Xp[43].x = 	0.008	;	 M[30].Xp[43].y =25.8	; */
/*  M[30].pathlength(); */
/* //station 32 */
/*  x[31] = 49.8; */
/*  M[31].allocate(44); */
/*  M[31].settype(SPLINE2D_LINEAR); */
/*  M[31].Xp[0].x = 	0.37	;	 M[31].Xp[0].y =0	; */
/*  M[31].Xp[1].x = 	0.38	;	 M[31].Xp[1].y =0.6	; */
/*  M[31].Xp[2].x = 	0.377	;	 M[31].Xp[2].y =1.2	; */
/*  M[31].Xp[3].x = 	0.36	;	 M[31].Xp[3].y =1.8	; */
/*  M[31].Xp[4].x = 	0.344	;	 M[31].Xp[4].y =2.4	; */
/*  M[31].Xp[5].x = 	0.317	;	 M[31].Xp[5].y =3	; */
/*  M[31].Xp[6].x = 	0.304	;	 M[31].Xp[6].y =3.6	; */
/*  M[31].Xp[7].x = 	0.269	;	 M[31].Xp[7].y =4.2	; */
/*  M[31].Xp[8].x = 	0.257	;	 M[31].Xp[8].y =4.8	; */
/*  M[31].Xp[9].x = 	0.25	;	 M[31].Xp[9].y =5.4	; */
/*  M[31].Xp[10].x = 	0.2	;	 M[31].Xp[10].y =6	; */
/*  M[31].Xp[11].x = 	0.184	;	 M[31].Xp[11].y =6.6	; */
/*  M[31].Xp[12].x = 	0.163	;	 M[31].Xp[12].y =7.2	; */
/*  M[31].Xp[13].x = 	0.145	;	 M[31].Xp[13].y =7.8	; */
/*  M[31].Xp[14].x = 	0.123	;	 M[31].Xp[14].y =8.4	; */
/*  M[31].Xp[15].x = 	0.099	;	 M[31].Xp[15].y =9	; */
/*  M[31].Xp[16].x = 	0.095	;	 M[31].Xp[16].y =9.6	; */
/*  M[31].Xp[17].x = 	0.083	;	 M[31].Xp[17].y =10.2	; */
/*  M[31].Xp[18].x = 	0.075	;	 M[31].Xp[18].y =10.8	; */
/*  M[31].Xp[19].x = 	0.065	;	 M[31].Xp[19].y =11.4	; */
/*  M[31].Xp[20].x = 	0.059	;	 M[31].Xp[20].y =12	; */
/*  M[31].Xp[21].x = 	0.053	;	 M[31].Xp[21].y =12.6	; */
/*  M[31].Xp[22].x = 	0.046	;	 M[31].Xp[22].y =13.2	; */
/*  M[31].Xp[23].x = 	0.046	;	 M[31].Xp[23].y =13.8	; */
/*  M[31].Xp[24].x = 	0.042	;	 M[31].Xp[24].y =14.4	; */
/*  M[31].Xp[25].x = 	0.04	;	 M[31].Xp[25].y =15	; */
/*  M[31].Xp[26].x = 	0.043	;	 M[31].Xp[26].y =15.6	; */
/*  M[31].Xp[27].x = 	0.036	;	 M[31].Xp[27].y =16.2	; */
/*  M[31].Xp[28].x = 	0.038	;	 M[31].Xp[28].y =16.8	; */
/*  M[31].Xp[29].x = 	0.027	;	 M[31].Xp[29].y =17.4	; */
/*  M[31].Xp[30].x = 	0.026	;	 M[31].Xp[30].y =18	; */
/*  M[31].Xp[31].x = 	0.014	;	 M[31].Xp[31].y =18.6	; */
/*  M[31].Xp[32].x = 	0.016	;	 M[31].Xp[32].y =19.2	; */
/*  M[31].Xp[33].x = 	0.017	;	 M[31].Xp[33].y =19.8	; */
/*  M[31].Xp[34].x = 	0.009	;	 M[31].Xp[34].y =20.4	; */
/*  M[31].Xp[35].x = 	0.009	;	 M[31].Xp[35].y =21	; */
/*  M[31].Xp[36].x = 	0.003	;	 M[31].Xp[36].y =21.6	; */
/*  M[31].Xp[37].x = 	0.004	;	 M[31].Xp[37].y =22.2	; */
/*  M[31].Xp[38].x = 	0.004	;	 M[31].Xp[38].y =22.8	; */
/*  M[31].Xp[39].x = 	0.009	;	 M[31].Xp[39].y =23.4	; */
/*  M[31].Xp[40].x = 	0.002	;	 M[31].Xp[40].y =24	; */
/*  M[31].Xp[41].x = 	0.003	;	 M[31].Xp[41].y =24.6	; */
/*  M[31].Xp[42].x = 	0.003	;	 M[31].Xp[42].y =25.2	; */
/*  M[31].Xp[43].x = 	0.001	;	 M[31].Xp[43].y =25.8	; */
/*  M[31].pathlength(); */
/* //station 33 */
/*  x[32] = 50.5; */
/*  M[32].allocate(44); */
/*  M[32].settype(SPLINE2D_LINEAR); */
/*  M[32].Xp[0].x = 	0.392	;	 M[32].Xp[0].y =0	; */
/*  M[32].Xp[1].x = 	0.366	;	 M[32].Xp[1].y =0.6	; */
/*  M[32].Xp[2].x = 	0.369	;	 M[32].Xp[2].y =1.2	; */
/*  M[32].Xp[3].x = 	0.353	;	 M[32].Xp[3].y =1.8	; */
/*  M[32].Xp[4].x = 	0.326	;	 M[32].Xp[4].y =2.4	; */
/*  M[32].Xp[5].x = 	0.292	;	 M[32].Xp[5].y =3	; */
/*  M[32].Xp[6].x = 	0.291	;	 M[32].Xp[6].y =3.6	; */
/*  M[32].Xp[7].x = 	0.285	;	 M[32].Xp[7].y =4.2	; */
/*  M[32].Xp[8].x = 	0.269	;	 M[32].Xp[8].y =4.8	; */
/*  M[32].Xp[9].x = 	0.246	;	 M[32].Xp[9].y =5.4	; */
/*  M[32].Xp[10].x = 	0.23	;	 M[32].Xp[10].y =6	; */
/*  M[32].Xp[11].x = 	0.212	;	 M[32].Xp[11].y =6.6	; */
/*  M[32].Xp[12].x = 	0.187	;	 M[32].Xp[12].y =7.2	; */
/*  M[32].Xp[13].x = 	0.149	;	 M[32].Xp[13].y =7.8	; */
/*  M[32].Xp[14].x = 	0.128	;	 M[32].Xp[14].y =8.4	; */
/*  M[32].Xp[15].x = 	0.111	;	 M[32].Xp[15].y =9	; */
/*  M[32].Xp[16].x = 	0.088	;	 M[32].Xp[16].y =9.6	; */
/*  M[32].Xp[17].x = 	0.08	;	 M[32].Xp[17].y =10.2	; */
/*  M[32].Xp[18].x = 	0.072	;	 M[32].Xp[18].y =10.8	; */
/*  M[32].Xp[19].x = 	0.065	;	 M[32].Xp[19].y =11.4	; */
/*  M[32].Xp[20].x = 	0.068	;	 M[32].Xp[20].y =12	; */
/*  M[32].Xp[21].x = 	0.049	;	 M[32].Xp[21].y =12.6	; */
/*  M[32].Xp[22].x = 	0.047	;	 M[32].Xp[22].y =13.2	; */
/*  M[32].Xp[23].x = 	0.045	;	 M[32].Xp[23].y =13.8	; */
/*  M[32].Xp[24].x = 	0.044	;	 M[32].Xp[24].y =14.4	; */
/*  M[32].Xp[25].x = 	0.045	;	 M[32].Xp[25].y =15	; */
/*  M[32].Xp[26].x = 	0.039	;	 M[32].Xp[26].y =15.6	; */
/*  M[32].Xp[27].x = 	0.037	;	 M[32].Xp[27].y =16.2	; */
/*  M[32].Xp[28].x = 	0.028	;	 M[32].Xp[28].y =16.8	; */
/*  M[32].Xp[29].x = 	0.039	;	 M[32].Xp[29].y =17.4	; */
/*  M[32].Xp[30].x = 	0.025	;	 M[32].Xp[30].y =18	; */
/*  M[32].Xp[31].x = 	0.019	;	 M[32].Xp[31].y =18.6	; */
/*  M[32].Xp[32].x = 	0.017	;	 M[32].Xp[32].y =19.2	; */
/*  M[32].Xp[33].x = 	0.012	;	 M[32].Xp[33].y =19.8	; */
/*  M[32].Xp[34].x = 	0.016	;	 M[32].Xp[34].y =20.4	; */
/*  M[32].Xp[35].x = 	0.009	;	 M[32].Xp[35].y =21	; */
/*  M[32].Xp[36].x = 	0.005	;	 M[32].Xp[36].y =21.6	; */
/*  M[32].Xp[37].x = 	0.009	;	 M[32].Xp[37].y =22.2	; */
/*  M[32].Xp[38].x = 	0.005	;	 M[32].Xp[38].y =22.8	; */
/*  M[32].Xp[39].x = 	0.009	;	 M[32].Xp[39].y =23.4	; */
/*  M[32].Xp[40].x = 	0.004	;	 M[32].Xp[40].y =24	; */
/*  M[32].Xp[41].x = 	0.006	;	 M[32].Xp[41].y =24.6	; */
/*  M[32].Xp[42].x = 	0.002	;	 M[32].Xp[42].y =25.2	; */
/*  M[32].Xp[43].x = 	0.008	;	 M[32].Xp[43].y =25.8	; */
/*  M[32].pathlength(); */
/* //station 34 */
/*  x[33] = 51.3; */
/*  M[33].allocate(44); */
/*  M[33].settype(SPLINE2D_LINEAR); */
/*  M[33].Xp[0].x = 	0.347	;	 M[33].Xp[0].y =0	; */
/*  M[33].Xp[1].x = 	0.383	;	 M[33].Xp[1].y =0.6	; */
/*  M[33].Xp[2].x = 	0.366	;	 M[33].Xp[2].y =1.2	; */
/*  M[33].Xp[3].x = 	0.344	;	 M[33].Xp[3].y =1.8	; */
/*  M[33].Xp[4].x = 	0.331	;	 M[33].Xp[4].y =2.4	; */
/*  M[33].Xp[5].x = 	0.311	;	 M[33].Xp[5].y =3	; */
/*  M[33].Xp[6].x = 	0.313	;	 M[33].Xp[6].y =3.6	; */
/*  M[33].Xp[7].x = 	0.291	;	 M[33].Xp[7].y =4.2	; */
/*  M[33].Xp[8].x = 	0.246	;	 M[33].Xp[8].y =4.8	; */
/*  M[33].Xp[9].x = 	0.25	;	 M[33].Xp[9].y =5.4	; */
/*  M[33].Xp[10].x = 	0.213	;	 M[33].Xp[10].y =6	; */
/*  M[33].Xp[11].x = 	0.174	;	 M[33].Xp[11].y =6.6	; */
/*  M[33].Xp[12].x = 	0.166	;	 M[33].Xp[12].y =7.2	; */
/*  M[33].Xp[13].x = 	0.14	;	 M[33].Xp[13].y =7.8	; */
/*  M[33].Xp[14].x = 	0.1	;	 M[33].Xp[14].y =8.4	; */
/*  M[33].Xp[15].x = 	0.099	;	 M[33].Xp[15].y =9	; */
/*  M[33].Xp[16].x = 	0.076	;	 M[33].Xp[16].y =9.6	; */
/*  M[33].Xp[17].x = 	0.071	;	 M[33].Xp[17].y =10.2	; */
/*  M[33].Xp[18].x = 	0.076	;	 M[33].Xp[18].y =10.8	; */
/*  M[33].Xp[19].x = 	0.065	;	 M[33].Xp[19].y =11.4	; */
/*  M[33].Xp[20].x = 	0.066	;	 M[33].Xp[20].y =12	; */
/*  M[33].Xp[21].x = 	0.054	;	 M[33].Xp[21].y =12.6	; */
/*  M[33].Xp[22].x = 	0.043	;	 M[33].Xp[22].y =13.2	; */
/*  M[33].Xp[23].x = 	0.046	;	 M[33].Xp[23].y =13.8	; */
/*  M[33].Xp[24].x = 	0.041	;	 M[33].Xp[24].y =14.4	; */
/*  M[33].Xp[25].x = 	0.047	;	 M[33].Xp[25].y =15	; */
/*  M[33].Xp[26].x = 	0.03	;	 M[33].Xp[26].y =15.6	; */
/*  M[33].Xp[27].x = 	0.035	;	 M[33].Xp[27].y =16.2	; */
/*  M[33].Xp[28].x = 	0.031	;	 M[33].Xp[28].y =16.8	; */
/*  M[33].Xp[29].x = 	0.024	;	 M[33].Xp[29].y =17.4	; */
/*  M[33].Xp[30].x = 	0.026	;	 M[33].Xp[30].y =18	; */
/*  M[33].Xp[31].x = 	0.023	;	 M[33].Xp[31].y =18.6	; */
/*  M[33].Xp[32].x = 	0.011	;	 M[33].Xp[32].y =19.2	; */
/*  M[33].Xp[33].x = 	0.018	;	 M[33].Xp[33].y =19.8	; */
/*  M[33].Xp[34].x = 	0.013	;	 M[33].Xp[34].y =20.4	; */
/*  M[33].Xp[35].x = 	0.012	;	 M[33].Xp[35].y =21	; */
/*  M[33].Xp[36].x = 	0.011	;	 M[33].Xp[36].y =21.6	; */
/*  M[33].Xp[37].x = 	0.006	;	 M[33].Xp[37].y =22.2	; */
/*  M[33].Xp[38].x = 	0.003	;	 M[33].Xp[38].y =22.8	; */
/*  M[33].Xp[39].x = 	0.006	;	 M[33].Xp[39].y =23.4	; */
/*  M[33].Xp[40].x = 	0.006	;	 M[33].Xp[40].y =24	; */
/*  M[33].Xp[41].x = 	0.006	;	 M[33].Xp[41].y =24.6	; */
/*  M[33].Xp[42].x = 	0.005	;	 M[33].Xp[42].y =25.2	; */
/*  M[33].Xp[43].x = 	0.005	;	 M[33].Xp[43].y =25.8	; */
/*  M[33].pathlength(); */
/* //station 35 */
/*  x[34] = 52.0; */
/*  M[34].allocate(44); */
/*  M[34].settype(SPLINE2D_LINEAR); */
/*  M[34].Xp[0].x = 	0.309	;	 M[34].Xp[0].y =0	; */
/*  M[34].Xp[1].x = 	0.324	;	 M[34].Xp[1].y =0.6	; */
/*  M[34].Xp[2].x = 	0.325	;	 M[34].Xp[2].y =1.2	; */
/*  M[34].Xp[3].x = 	0.296	;	 M[34].Xp[3].y =1.8	; */
/*  M[34].Xp[4].x = 	0.277	;	 M[34].Xp[4].y =2.4	; */
/*  M[34].Xp[5].x = 	0.278	;	 M[34].Xp[5].y =3	; */
/*  M[34].Xp[6].x = 	0.295	;	 M[34].Xp[6].y =3.6	; */
/*  M[34].Xp[7].x = 	0.245	;	 M[34].Xp[7].y =4.2	; */
/*  M[34].Xp[8].x = 	0.243	;	 M[34].Xp[8].y =4.8	; */
/*  M[34].Xp[9].x = 	0.213	;	 M[34].Xp[9].y =5.4	; */
/*  M[34].Xp[10].x = 	0.181	;	 M[34].Xp[10].y =6	; */
/*  M[34].Xp[11].x = 	0.185	;	 M[34].Xp[11].y =6.6	; */
/*  M[34].Xp[12].x = 	0.16	;	 M[34].Xp[12].y =7.2	; */
/*  M[34].Xp[13].x = 	0.128	;	 M[34].Xp[13].y =7.8	; */
/*  M[34].Xp[14].x = 	0.127	;	 M[34].Xp[14].y =8.4	; */
/*  M[34].Xp[15].x = 	0.101	;	 M[34].Xp[15].y =9	; */
/*  M[34].Xp[16].x = 	0.066	;	 M[34].Xp[16].y =9.6	; */
/*  M[34].Xp[17].x = 	0.069	;	 M[34].Xp[17].y =10.2	; */
/*  M[34].Xp[18].x = 	0.057	;	 M[34].Xp[18].y =10.8	; */
/*  M[34].Xp[19].x = 	0.055	;	 M[34].Xp[19].y =11.4	; */
/*  M[34].Xp[20].x = 	0.055	;	 M[34].Xp[20].y =12	; */
/*  M[34].Xp[21].x = 	0.06	;	 M[34].Xp[21].y =12.6	; */
/*  M[34].Xp[22].x = 	0.042	;	 M[34].Xp[22].y =13.2	; */
/*  M[34].Xp[23].x = 	0.037	;	 M[34].Xp[23].y =13.8	; */
/*  M[34].Xp[24].x = 	0.046	;	 M[34].Xp[24].y =14.4	; */
/*  M[34].Xp[25].x = 	0.033	;	 M[34].Xp[25].y =15	; */
/*  M[34].Xp[26].x = 	0.038	;	 M[34].Xp[26].y =15.6	; */
/*  M[34].Xp[27].x = 	0.027	;	 M[34].Xp[27].y =16.2	; */
/*  M[34].Xp[28].x = 	0.032	;	 M[34].Xp[28].y =16.8	; */
/*  M[34].Xp[29].x = 	0.015	;	 M[34].Xp[29].y =17.4	; */
/*  M[34].Xp[30].x = 	0.026	;	 M[34].Xp[30].y =18	; */
/*  M[34].Xp[31].x = 	0.015	;	 M[34].Xp[31].y =18.6	; */
/*  M[34].Xp[32].x = 	0.015	;	 M[34].Xp[32].y =19.2	; */
/*  M[34].Xp[33].x = 	0.009	;	 M[34].Xp[33].y =19.8	; */
/*  M[34].Xp[34].x = 	0.014	;	 M[34].Xp[34].y =20.4	; */
/*  M[34].Xp[35].x = 	0.011	;	 M[34].Xp[35].y =21	; */
/*  M[34].Xp[36].x = 	0.01	;	 M[34].Xp[36].y =21.6	; */
/*  M[34].Xp[37].x = 	0.006	;	 M[34].Xp[37].y =22.2	; */
/*  M[34].Xp[38].x = 	0.002	;	 M[34].Xp[38].y =22.8	; */
/*  M[34].Xp[39].x = 	0.004	;	 M[34].Xp[39].y =23.4	; */
/*  M[34].Xp[40].x = 	0.005	;	 M[34].Xp[40].y =24	; */
/*  M[34].Xp[41].x = 	0.006	;	 M[34].Xp[41].y =24.6	; */
/*  M[34].Xp[42].x = 	0.005	;	 M[34].Xp[42].y =25.2	; */
/*  M[34].Xp[43].x = 	0.008	;	 M[34].Xp[43].y =25.8	; */
/*  M[34].pathlength(); */
/* //station 36 */
/*  x[35] = 52.8; */
/*  M[35].allocate(44); */
/*  M[35].settype(SPLINE2D_LINEAR); */
/*  M[35].Xp[0].x = 	0.271	;	 M[35].Xp[0].y =0	; */
/*  M[35].Xp[1].x = 	0.303	;	 M[35].Xp[1].y =0.6	; */
/*  M[35].Xp[2].x = 	0.308	;	 M[35].Xp[2].y =1.2	; */
/*  M[35].Xp[3].x = 	0.28	;	 M[35].Xp[3].y =1.8	; */
/*  M[35].Xp[4].x = 	0.28	;	 M[35].Xp[4].y =2.4	; */
/*  M[35].Xp[5].x = 	0.282	;	 M[35].Xp[5].y =3	; */
/*  M[35].Xp[6].x = 	0.23	;	 M[35].Xp[6].y =3.6	; */
/*  M[35].Xp[7].x = 	0.227	;	 M[35].Xp[7].y =4.2	; */
/*  M[35].Xp[8].x = 	0.236	;	 M[35].Xp[8].y =4.8	; */
/*  M[35].Xp[9].x = 	0.203	;	 M[35].Xp[9].y =5.4	; */
/*  M[35].Xp[10].x = 	0.153	;	 M[35].Xp[10].y =6	; */
/*  M[35].Xp[11].x = 	0.159	;	 M[35].Xp[11].y =6.6	; */
/*  M[35].Xp[12].x = 	0.115	;	 M[35].Xp[12].y =7.2	; */
/*  M[35].Xp[13].x = 	0.11	;	 M[35].Xp[13].y =7.8	; */
/*  M[35].Xp[14].x = 	0.105	;	 M[35].Xp[14].y =8.4	; */
/*  M[35].Xp[15].x = 	0.101	;	 M[35].Xp[15].y =9	; */
/*  M[35].Xp[16].x = 	0.072	;	 M[35].Xp[16].y =9.6	; */
/*  M[35].Xp[17].x = 	0.066	;	 M[35].Xp[17].y =10.2	; */
/*  M[35].Xp[18].x = 	0.068	;	 M[35].Xp[18].y =10.8	; */
/*  M[35].Xp[19].x = 	0.073	;	 M[35].Xp[19].y =11.4	; */
/*  M[35].Xp[20].x = 	0.055	;	 M[35].Xp[20].y =12	; */
/*  M[35].Xp[21].x = 	0.057	;	 M[35].Xp[21].y =12.6	; */
/*  M[35].Xp[22].x = 	0.055	;	 M[35].Xp[22].y =13.2	; */
/*  M[35].Xp[23].x = 	0.042	;	 M[35].Xp[23].y =13.8	; */
/*  M[35].Xp[24].x = 	0.039	;	 M[35].Xp[24].y =14.4	; */
/*  M[35].Xp[25].x = 	0.048	;	 M[35].Xp[25].y =15	; */
/*  M[35].Xp[26].x = 	0.024	;	 M[35].Xp[26].y =15.6	; */
/*  M[35].Xp[27].x = 	0.035	;	 M[35].Xp[27].y =16.2	; */
/*  M[35].Xp[28].x = 	0.019	;	 M[35].Xp[28].y =16.8	; */
/*  M[35].Xp[29].x = 	0.027	;	 M[35].Xp[29].y =17.4	; */
/*  M[35].Xp[30].x = 	0.021	;	 M[35].Xp[30].y =18	; */
/*  M[35].Xp[31].x = 	0.027	;	 M[35].Xp[31].y =18.6	; */
/*  M[35].Xp[32].x = 	0.012	;	 M[35].Xp[32].y =19.2	; */
/*  M[35].Xp[33].x = 	0.019	;	 M[35].Xp[33].y =19.8	; */
/*  M[35].Xp[34].x = 	0.01	;	 M[35].Xp[34].y =20.4	; */
/*  M[35].Xp[35].x = 	0.007	;	 M[35].Xp[35].y =21	; */
/*  M[35].Xp[36].x = 	0.007	;	 M[35].Xp[36].y =21.6	; */
/*  M[35].Xp[37].x = 	0.004	;	 M[35].Xp[37].y =22.2	; */
/*  M[35].Xp[38].x = 	0.007	;	 M[35].Xp[38].y =22.8	; */
/*  M[35].Xp[39].x = 	0.012	;	 M[35].Xp[39].y =23.4	; */
/*  M[35].Xp[40].x = 	0.007	;	 M[35].Xp[40].y =24	; */
/*  M[35].Xp[41].x = 	0.007	;	 M[35].Xp[41].y =24.6	; */
/*  M[35].Xp[42].x = 	0.005	;	 M[35].Xp[42].y =25.2	; */
/*  M[35].Xp[43].x = 	0.002	;	 M[35].Xp[43].y =25.8	; */

/*  M[35].pathlength(); */
/* //station 37 */
/*  x[36] = 55.0; */
/*  M[36].allocate(44); */
/*  M[36].settype(SPLINE2D_LINEAR); */
/*  M[36].Xp[0].x = 	0.222	;	 M[36].Xp[0].y =0.2	; */
/*  M[36].Xp[1].x = 	0.216	;	 M[36].Xp[1].y =0.8	; */
/*  M[36].Xp[2].x = 	0.217	;	 M[36].Xp[2].y =1.4	; */
/*  M[36].Xp[3].x = 	0.195	;	 M[36].Xp[3].y =2	; */
/*  M[36].Xp[4].x = 	0.169	;	 M[36].Xp[4].y =2.6	; */
/*  M[36].Xp[5].x = 	0.157	;	 M[36].Xp[5].y =3.2	; */
/*  M[36].Xp[6].x = 	0.16	;	 M[36].Xp[6].y =3.8	; */
/*  M[36].Xp[7].x = 	0.161	;	 M[36].Xp[7].y =4.4	; */
/*  M[36].Xp[8].x = 	0.135	;	 M[36].Xp[8].y =5	; */
/*  M[36].Xp[9].x = 	0.128	;	 M[36].Xp[9].y =5.6	; */
/*  M[36].Xp[10].x = 	0.112	;	 M[36].Xp[10].y =6.2	; */
/*  M[36].Xp[11].x = 	0.106	;	 M[36].Xp[11].y =6.8	; */
/*  M[36].Xp[12].x = 	0.075	;	 M[36].Xp[12].y =7.4	; */
/*  M[36].Xp[13].x = 	0.075	;	 M[36].Xp[13].y =8	; */
/*  M[36].Xp[14].x = 	0.055	;	 M[36].Xp[14].y =8.6	; */
/*  M[36].Xp[15].x = 	0.069	;	 M[36].Xp[15].y =9.2	; */
/*  M[36].Xp[16].x = 	0.051	;	 M[36].Xp[16].y =9.8	; */
/*  M[36].Xp[17].x = 	0.029	;	 M[36].Xp[17].y =10.4	; */
/*  M[36].Xp[18].x = 	0.054	;	 M[36].Xp[18].y =11	; */
/*  M[36].Xp[19].x = 	0.055	;	 M[36].Xp[19].y =11.6	; */
/*  M[36].Xp[20].x = 	0.037	;	 M[36].Xp[20].y =12.2	; */
/*  M[36].Xp[21].x = 	0.036	;	 M[36].Xp[21].y =12.8	; */
/*  M[36].Xp[22].x = 	0.031	;	 M[36].Xp[22].y =13.4	; */
/*  M[36].Xp[23].x = 	0.03	;	 M[36].Xp[23].y =14	; */
/*  M[36].Xp[24].x = 	0.036	;	 M[36].Xp[24].y =14.6	; */
/*  M[36].Xp[25].x = 	0.033	;	 M[36].Xp[25].y =15.2	; */
/*  M[36].Xp[26].x = 	0.034	;	 M[36].Xp[26].y =15.8	; */
/*  M[36].Xp[27].x = 	0.023	;	 M[36].Xp[27].y =16.4	; */
/*  M[36].Xp[28].x = 	0.017	;	 M[36].Xp[28].y =17	; */
/*  M[36].Xp[29].x = 	0.025	;	 M[36].Xp[29].y =17.6	; */
/*  M[36].Xp[30].x = 	0.017	;	 M[36].Xp[30].y =18.2	; */
/*  M[36].Xp[31].x = 	0.016	;	 M[36].Xp[31].y =18.8	; */
/*  M[36].Xp[32].x = 	0.011	;	 M[36].Xp[32].y =19.4	; */
/*  M[36].Xp[33].x = 	0.009	;	 M[36].Xp[33].y =20	; */
/*  M[36].Xp[34].x = 	0.007	;	 M[36].Xp[34].y =20.6	; */
/*  M[36].Xp[35].x = 	0.006	;	 M[36].Xp[35].y =21.2	; */
/*  M[36].Xp[36].x = 	0.004	;	 M[36].Xp[36].y =21.8	; */
/*  M[36].Xp[37].x = 	0.009	;	 M[36].Xp[37].y =22.4	; */
/*  M[36].Xp[38].x = 	0.006	;	 M[36].Xp[38].y =23	; */
/*  M[36].Xp[39].x = 	0.002	;	 M[36].Xp[39].y =23.6	; */
/*  M[36].Xp[40].x = 	0.009	;	 M[36].Xp[40].y =24.2	; */
/*  M[36].Xp[41].x = 	0.01	;	 M[36].Xp[41].y =24.8	; */
/*  M[36].Xp[42].x = 	0.008	;	 M[36].Xp[42].y =25.4	; */
/*  M[36].Xp[43].x = 	0.005	;	 M[36].Xp[43].y =26	; */
/*  M[36].pathlength(); */
/* //station 38 */
/*  x[37] = 55.8; */
/*  M[37].allocate(44); */
/*  M[37].settype(SPLINE2D_LINEAR); */
/*  M[37].Xp[0].x = 	0.196	;	 M[37].Xp[0].y =0.2	; */
/*  M[37].Xp[1].x = 	0.201	;	 M[37].Xp[1].y =0.8	; */
/*  M[37].Xp[2].x = 	0.184	;	 M[37].Xp[2].y =1.4	; */
/*  M[37].Xp[3].x = 	0.197	;	 M[37].Xp[3].y =2	; */
/*  M[37].Xp[4].x = 	0.172	;	 M[37].Xp[4].y =2.6	; */
/*  M[37].Xp[5].x = 	0.136	;	 M[37].Xp[5].y =3.2	; */
/*  M[37].Xp[6].x = 	0.134	;	 M[37].Xp[6].y =3.8	; */
/*  M[37].Xp[7].x = 	0.146	;	 M[37].Xp[7].y =4.4	; */
/*  M[37].Xp[8].x = 	0.129	;	 M[37].Xp[8].y =5	; */
/*  M[37].Xp[9].x = 	0.103	;	 M[37].Xp[9].y =5.6	; */
/*  M[37].Xp[10].x = 	0.089	;	 M[37].Xp[10].y =6.2	; */
/*  M[37].Xp[11].x = 	0.101	;	 M[37].Xp[11].y =6.8	; */
/*  M[37].Xp[12].x = 	0.092	;	 M[37].Xp[12].y =7.4	; */
/*  M[37].Xp[13].x = 	0.064	;	 M[37].Xp[13].y =8	; */
/*  M[37].Xp[14].x = 	0.045	;	 M[37].Xp[14].y =8.6	; */
/*  M[37].Xp[15].x = 	0.046	;	 M[37].Xp[15].y =9.2	; */
/*  M[37].Xp[16].x = 	0.052	;	 M[37].Xp[16].y =9.8	; */
/*  M[37].Xp[17].x = 	0.041	;	 M[37].Xp[17].y =10.4	; */
/*  M[37].Xp[18].x = 	0.042	;	 M[37].Xp[18].y =11	; */
/*  M[37].Xp[19].x = 	0.039	;	 M[37].Xp[19].y =11.6	; */
/*  M[37].Xp[20].x = 	0.019	;	 M[37].Xp[20].y =12.2	; */
/*  M[37].Xp[21].x = 	0.022	;	 M[37].Xp[21].y =12.8	; */
/*  M[37].Xp[22].x = 	0.034	;	 M[37].Xp[22].y =13.4	; */
/*  M[37].Xp[23].x = 	0.037	;	 M[37].Xp[23].y =14	; */
/*  M[37].Xp[24].x = 	0.027	;	 M[37].Xp[24].y =14.6	; */
/*  M[37].Xp[25].x = 	0.028	;	 M[37].Xp[25].y =15.2	; */
/*  M[37].Xp[26].x = 	0.021	;	 M[37].Xp[26].y =15.8	; */
/*  M[37].Xp[27].x = 	0.021	;	 M[37].Xp[27].y =16.4	; */
/*  M[37].Xp[28].x = 	0.028	;	 M[37].Xp[28].y =17	; */
/*  M[37].Xp[29].x = 	0.018	;	 M[37].Xp[29].y =17.6	; */
/*  M[37].Xp[30].x = 	0.018	;	 M[37].Xp[30].y =18.2	; */
/*  M[37].Xp[31].x = 	0.01	;	 M[37].Xp[31].y =18.8	; */
/*  M[37].Xp[32].x = 	0.015	;	 M[37].Xp[32].y =19.4	; */
/*  M[37].Xp[33].x = 	0.012	;	 M[37].Xp[33].y =20	; */
/*  M[37].Xp[34].x = 	0.007	;	 M[37].Xp[34].y =20.6	; */
/*  M[37].Xp[35].x = 	0.008	;	 M[37].Xp[35].y =21.2	; */
/*  M[37].Xp[36].x = 	0.005	;	 M[37].Xp[36].y =21.8	; */
/*  M[37].Xp[37].x = 	0.006	;	 M[37].Xp[37].y =22.4	; */
/*  M[37].Xp[38].x = 	0.002	;	 M[37].Xp[38].y =23	; */
/*  M[37].Xp[39].x = 	0.002	;	 M[37].Xp[39].y =23.6	; */
/*  M[37].Xp[40].x = 	0.006	;	 M[37].Xp[40].y =24.2	; */
/*  M[37].Xp[41].x = 	0.002	;	 M[37].Xp[41].y =24.8	; */
/*  M[37].Xp[42].x = 	0.003	;	 M[37].Xp[42].y =25.4	; */
/*  M[37].Xp[43].x = 	0.001	;	 M[37].Xp[43].y =26	; */
/*  M[37].pathlength(); */
/* //station 39 */
/*  x[38] = 56.5; */
/*  M[38].allocate(44); */
/*  M[38].settype(SPLINE2D_LINEAR); */
/*  M[38].Xp[0].x = 	0.181	;	 M[38].Xp[0].y =0.2	; */
/*  M[38].Xp[1].x = 	0.196	;	 M[38].Xp[1].y =0.8	; */
/*  M[38].Xp[2].x = 	0.182	;	 M[38].Xp[2].y =1.4	; */
/*  M[38].Xp[3].x = 	0.184	;	 M[38].Xp[3].y =2	; */
/*  M[38].Xp[4].x = 	0.186	;	 M[38].Xp[4].y =2.6	; */
/*  M[38].Xp[5].x = 	0.17	;	 M[38].Xp[5].y =3.2	; */
/*  M[38].Xp[6].x = 	0.141	;	 M[38].Xp[6].y =3.8	; */
/*  M[38].Xp[7].x = 	0.135	;	 M[38].Xp[7].y =4.4	; */
/*  M[38].Xp[8].x = 	0.127	;	 M[38].Xp[8].y =5	; */
/*  M[38].Xp[9].x = 	0.129	;	 M[38].Xp[9].y =5.6	; */
/*  M[38].Xp[10].x = 	0.083	;	 M[38].Xp[10].y =6.2	; */
/*  M[38].Xp[11].x = 	0.09	;	 M[38].Xp[11].y =6.8	; */
/*  M[38].Xp[12].x = 	0.088	;	 M[38].Xp[12].y =7.4	; */
/*  M[38].Xp[13].x = 	0.071	;	 M[38].Xp[13].y =8	; */
/*  M[38].Xp[14].x = 	0.057	;	 M[38].Xp[14].y =8.6	; */
/*  M[38].Xp[15].x = 	0.055	;	 M[38].Xp[15].y =9.2	; */
/*  M[38].Xp[16].x = 	0.057	;	 M[38].Xp[16].y =9.8	; */
/*  M[38].Xp[17].x = 	0.055	;	 M[38].Xp[17].y =10.4	; */
/*  M[38].Xp[18].x = 	0.035	;	 M[38].Xp[18].y =11	; */
/*  M[38].Xp[19].x = 	0.035	;	 M[38].Xp[19].y =11.6	; */
/*  M[38].Xp[20].x = 	0.035	;	 M[38].Xp[20].y =12.2	; */
/*  M[38].Xp[21].x = 	0.036	;	 M[38].Xp[21].y =12.8	; */
/*  M[38].Xp[22].x = 	0.037	;	 M[38].Xp[22].y =13.4	; */
/*  M[38].Xp[23].x = 	0.024	;	 M[38].Xp[23].y =14	; */
/*  M[38].Xp[24].x = 	0.031	;	 M[38].Xp[24].y =14.6	; */
/*  M[38].Xp[25].x = 	0.025	;	 M[38].Xp[25].y =15.2	; */
/*  M[38].Xp[26].x = 	0.025	;	 M[38].Xp[26].y =15.8	; */
/*  M[38].Xp[27].x = 	0.021	;	 M[38].Xp[27].y =16.4	; */
/*  M[38].Xp[28].x = 	0.017	;	 M[38].Xp[28].y =17	; */
/*  M[38].Xp[29].x = 	0.022	;	 M[38].Xp[29].y =17.6	; */
/*  M[38].Xp[30].x = 	0.022	;	 M[38].Xp[30].y =18.2	; */
/*  M[38].Xp[31].x = 	0.021	;	 M[38].Xp[31].y =18.8	; */
/*  M[38].Xp[32].x = 	0.01	;	 M[38].Xp[32].y =19.4	; */
/*  M[38].Xp[33].x = 	0.012	;	 M[38].Xp[33].y =20	; */
/*  M[38].Xp[34].x = 	0.011	;	 M[38].Xp[34].y =20.6	; */
/*  M[38].Xp[35].x = 	0.003	;	 M[38].Xp[35].y =21.2	; */
/*  M[38].Xp[36].x = 	0.004	;	 M[38].Xp[36].y =21.8	; */
/*  M[38].Xp[37].x = 	0.008	;	 M[38].Xp[37].y =22.4	; */
/*  M[38].Xp[38].x = 	0.002	;	 M[38].Xp[38].y =23	; */
/*  M[38].Xp[39].x = 	0.002	;	 M[38].Xp[39].y =23.6	; */
/*  M[38].Xp[40].x = 	0.003	;	 M[38].Xp[40].y =24.2	; */
/*  M[38].Xp[41].x = 	0.003	;	 M[38].Xp[41].y =24.8	; */
/*  M[38].Xp[42].x = 	0.001	;	 M[38].Xp[42].y =25.4	; */
/*  M[38].Xp[43].x = 	0.005	;	 M[38].Xp[43].y =26	; */
/*  M[38].pathlength(); */
/* //station 40 */
/*  x[39] = 57.3; */
/*  M[39].allocate(44); */
/*  M[39].settype(SPLINE2D_LINEAR); */
/*  M[39].Xp[0].x = 	0.177	;	 M[39].Xp[0].y =0.2	; */
/*  M[39].Xp[1].x = 	0.194	;	 M[39].Xp[1].y =0.8	; */
/*  M[39].Xp[2].x = 	0.165	;	 M[39].Xp[2].y =1.4	; */
/*  M[39].Xp[3].x = 	0.166	;	 M[39].Xp[3].y =2	; */
/*  M[39].Xp[4].x = 	0.154	;	 M[39].Xp[4].y =2.6	; */
/*  M[39].Xp[5].x = 	0.166	;	 M[39].Xp[5].y =3.2	; */
/*  M[39].Xp[6].x = 	0.13	;	 M[39].Xp[6].y =3.8	; */
/*  M[39].Xp[7].x = 	0.131	;	 M[39].Xp[7].y =4.4	; */
/*  M[39].Xp[8].x = 	0.127	;	 M[39].Xp[8].y =5	; */
/*  M[39].Xp[9].x = 	0.111	;	 M[39].Xp[9].y =5.6	; */
/*  M[39].Xp[10].x = 	0.082	;	 M[39].Xp[10].y =6.2	; */
/*  M[39].Xp[11].x = 	0.08	;	 M[39].Xp[11].y =6.8	; */
/*  M[39].Xp[12].x = 	0.07	;	 M[39].Xp[12].y =7.4	; */
/*  M[39].Xp[13].x = 	0.082	;	 M[39].Xp[13].y =8	; */
/*  M[39].Xp[14].x = 	0.063	;	 M[39].Xp[14].y =8.6	; */
/*  M[39].Xp[15].x = 	0.045	;	 M[39].Xp[15].y =9.2	; */
/*  M[39].Xp[16].x = 	0.042	;	 M[39].Xp[16].y =9.8	; */
/*  M[39].Xp[17].x = 	0.064	;	 M[39].Xp[17].y =10.4	; */
/*  M[39].Xp[18].x = 	0.045	;	 M[39].Xp[18].y =11	; */
/*  M[39].Xp[19].x = 	0.045	;	 M[39].Xp[19].y =11.6	; */
/*  M[39].Xp[20].x = 	0.042	;	 M[39].Xp[20].y =12.2	; */
/*  M[39].Xp[21].x = 	0.022	;	 M[39].Xp[21].y =12.8	; */
/*  M[39].Xp[22].x = 	0.029	;	 M[39].Xp[22].y =13.4	; */
/*  M[39].Xp[23].x = 	0.034	;	 M[39].Xp[23].y =14	; */
/*  M[39].Xp[24].x = 	0.031	;	 M[39].Xp[24].y =14.6	; */
/*  M[39].Xp[25].x = 	0.029	;	 M[39].Xp[25].y =15.2	; */
/*  M[39].Xp[26].x = 	0.033	;	 M[39].Xp[26].y =15.8	; */
/*  M[39].Xp[27].x = 	0.02	;	 M[39].Xp[27].y =16.4	; */
/*  M[39].Xp[28].x = 	0.025	;	 M[39].Xp[28].y =17	; */
/*  M[39].Xp[29].x = 	0.023	;	 M[39].Xp[29].y =17.6	; */
/*  M[39].Xp[30].x = 	0.02	;	 M[39].Xp[30].y =18.2	; */
/*  M[39].Xp[31].x = 	0.012	;	 M[39].Xp[31].y =18.8	; */
/*  M[39].Xp[32].x = 	0.014	;	 M[39].Xp[32].y =19.4	; */
/*  M[39].Xp[33].x = 	0.007	;	 M[39].Xp[33].y =20	; */
/*  M[39].Xp[34].x = 	0.004	;	 M[39].Xp[34].y =20.6	; */
/*  M[39].Xp[35].x = 	0.006	;	 M[39].Xp[35].y =21.2	; */
/*  M[39].Xp[36].x = 	0.004	;	 M[39].Xp[36].y =21.8	; */
/*  M[39].Xp[37].x = 	0.002	;	 M[39].Xp[37].y =22.4	; */
/*  M[39].Xp[38].x = 	0.004	;	 M[39].Xp[38].y =23	; */
/*  M[39].Xp[39].x = 	0.008	;	 M[39].Xp[39].y =23.6	; */
/*  M[39].Xp[40].x = 	0.009	;	 M[39].Xp[40].y =24.2	; */
/*  M[39].Xp[41].x = 	0.005	;	 M[39].Xp[41].y =24.8	; */
/*  M[39].Xp[42].x = 	0.01	;	 M[39].Xp[42].y =25.4	; */
/*  M[39].Xp[43].x = 	0.004	;	 M[39].Xp[43].y =26	; */
/*  M[39].pathlength(); */
/* //station 41 */
/*  x[40] = 58.0; */
/*  M[40].allocate(44); */
/*  M[40].settype(SPLINE2D_LINEAR); */
/*  M[40].Xp[0].x = 	0.169	;	 M[40].Xp[0].y =0.2	; */
/*  M[40].Xp[1].x = 	0.167	;	 M[40].Xp[1].y =0.8	; */
/*  M[40].Xp[2].x = 	0.163	;	 M[40].Xp[2].y =1.4	; */
/*  M[40].Xp[3].x = 	0.161	;	 M[40].Xp[3].y =2	; */
/*  M[40].Xp[4].x = 	0.125	;	 M[40].Xp[4].y =2.6	; */
/*  M[40].Xp[5].x = 	0.125	;	 M[40].Xp[5].y =3.2	; */
/*  M[40].Xp[6].x = 	0.14	;	 M[40].Xp[6].y =3.8	; */
/*  M[40].Xp[7].x = 	0.124	;	 M[40].Xp[7].y =4.4	; */
/*  M[40].Xp[8].x = 	0.11	;	 M[40].Xp[8].y =5	; */
/*  M[40].Xp[9].x = 	0.083	;	 M[40].Xp[9].y =5.6	; */
/*  M[40].Xp[10].x = 	0.089	;	 M[40].Xp[10].y =6.2	; */
/*  M[40].Xp[11].x = 	0.082	;	 M[40].Xp[11].y =6.8	; */
/*  M[40].Xp[12].x = 	0.07	;	 M[40].Xp[12].y =7.4	; */
/*  M[40].Xp[13].x = 	0.064	;	 M[40].Xp[13].y =8	; */
/*  M[40].Xp[14].x = 	0.064	;	 M[40].Xp[14].y =8.6	; */
/*  M[40].Xp[15].x = 	0.048	;	 M[40].Xp[15].y =9.2	; */
/*  M[40].Xp[16].x = 	0.053	;	 M[40].Xp[16].y =9.8	; */
/*  M[40].Xp[17].x = 	0.053	;	 M[40].Xp[17].y =10.4	; */
/*  M[40].Xp[18].x = 	0.058	;	 M[40].Xp[18].y =11	; */
/*  M[40].Xp[19].x = 	0.039	;	 M[40].Xp[19].y =11.6	; */
/*  M[40].Xp[20].x = 	0.045	;	 M[40].Xp[20].y =12.2	; */
/*  M[40].Xp[21].x = 	0.03	;	 M[40].Xp[21].y =12.8	; */
/*  M[40].Xp[22].x = 	0.041	;	 M[40].Xp[22].y =13.4	; */
/*  M[40].Xp[23].x = 	0.027	;	 M[40].Xp[23].y =14	; */
/*  M[40].Xp[24].x = 	0.029	;	 M[40].Xp[24].y =14.6	; */
/*  M[40].Xp[25].x = 	0.025	;	 M[40].Xp[25].y =15.2	; */
/*  M[40].Xp[26].x = 	0.026	;	 M[40].Xp[26].y =15.8	; */
/*  M[40].Xp[27].x = 	0.021	;	 M[40].Xp[27].y =16.4	; */
/*  M[40].Xp[28].x = 	0.022	;	 M[40].Xp[28].y =17	; */
/*  M[40].Xp[29].x = 	0.014	;	 M[40].Xp[29].y =17.6	; */
/*  M[40].Xp[30].x = 	0.022	;	 M[40].Xp[30].y =18.2	; */
/*  M[40].Xp[31].x = 	0.02	;	 M[40].Xp[31].y =18.8	; */
/*  M[40].Xp[32].x = 	0.013	;	 M[40].Xp[32].y =19.4	; */
/*  M[40].Xp[33].x = 	0.009	;	 M[40].Xp[33].y =20	; */
/*  M[40].Xp[34].x = 	0.01	;	 M[40].Xp[34].y =20.6	; */
/*  M[40].Xp[35].x = 	0.007	;	 M[40].Xp[35].y =21.2	; */
/*  M[40].Xp[36].x = 	0.009	;	 M[40].Xp[36].y =21.8	; */
/*  M[40].Xp[37].x = 	0.003	;	 M[40].Xp[37].y =22.4	; */
/*  M[40].Xp[38].x = 	0.002	;	 M[40].Xp[38].y =23	; */
/*  M[40].Xp[39].x = 	0.002	;	 M[40].Xp[39].y =23.6	; */
/*  M[40].Xp[40].x = 	0.004	;	 M[40].Xp[40].y =24.2	; */
/*  M[40].Xp[41].x = 	0.002	;	 M[40].Xp[41].y =24.8	; */
/*  M[40].Xp[42].x = 	0.004	;	 M[40].Xp[42].y =25.4	; */
/*  M[40].Xp[43].x = 	0.008	;	 M[40].Xp[43].y =26	; */

/*  M[40].pathlength(); */
/* //station 42 */
/*  x[41] = 58.8; */
/*  M[41].allocate(44); */
/*  M[41].settype(SPLINE2D_LINEAR); */
/*  M[41].Xp[0].x = 	0.169	;	 M[41].Xp[0].y =0.2	; */
/*  M[41].Xp[1].x = 	0.139	;	 M[41].Xp[1].y =0.8	; */
/*  M[41].Xp[2].x = 	0.159	;	 M[41].Xp[2].y =1.4	; */
/*  M[41].Xp[3].x = 	0.139	;	 M[41].Xp[3].y =2	; */
/*  M[41].Xp[4].x = 	0.136	;	 M[41].Xp[4].y =2.6	; */
/*  M[41].Xp[5].x = 	0.144	;	 M[41].Xp[5].y =3.2	; */
/*  M[41].Xp[6].x = 	0.132	;	 M[41].Xp[6].y =3.8	; */
/*  M[41].Xp[7].x = 	0.121	;	 M[41].Xp[7].y =4.4	; */
/*  M[41].Xp[8].x = 	0.083	;	 M[41].Xp[8].y =5	; */
/*  M[41].Xp[9].x = 	0.077	;	 M[41].Xp[9].y =5.6	; */
/*  M[41].Xp[10].x = 	0.088	;	 M[41].Xp[10].y =6.2	; */
/*  M[41].Xp[11].x = 	0.065	;	 M[41].Xp[11].y =6.8	; */
/*  M[41].Xp[12].x = 	0.059	;	 M[41].Xp[12].y =7.4	; */
/*  M[41].Xp[13].x = 	0.059	;	 M[41].Xp[13].y =8	; */
/*  M[41].Xp[14].x = 	0.059	;	 M[41].Xp[14].y =8.6	; */
/*  M[41].Xp[15].x = 	0.044	;	 M[41].Xp[15].y =9.2	; */
/*  M[41].Xp[16].x = 	0.05	;	 M[41].Xp[16].y =9.8	; */
/*  M[41].Xp[17].x = 	0.053	;	 M[41].Xp[17].y =10.4	; */
/*  M[41].Xp[18].x = 	0.044	;	 M[41].Xp[18].y =11	; */
/*  M[41].Xp[19].x = 	0.032	;	 M[41].Xp[19].y =11.6	; */
/*  M[41].Xp[20].x = 	0.04	;	 M[41].Xp[20].y =12.2	; */
/*  M[41].Xp[21].x = 	0.038	;	 M[41].Xp[21].y =12.8	; */
/*  M[41].Xp[22].x = 	0.039	;	 M[41].Xp[22].y =13.4	; */
/*  M[41].Xp[23].x = 	0.032	;	 M[41].Xp[23].y =14	; */
/*  M[41].Xp[24].x = 	0.026	;	 M[41].Xp[24].y =14.6	; */
/*  M[41].Xp[25].x = 	0.029	;	 M[41].Xp[25].y =15.2	; */
/*  M[41].Xp[26].x = 	0.019	;	 M[41].Xp[26].y =15.8	; */
/*  M[41].Xp[27].x = 	0.018	;	 M[41].Xp[27].y =16.4	; */
/*  M[41].Xp[28].x = 	0.018	;	 M[41].Xp[28].y =17	; */
/*  M[41].Xp[29].x = 	0.014	;	 M[41].Xp[29].y =17.6	; */
/*  M[41].Xp[30].x = 	0.018	;	 M[41].Xp[30].y =18.2	; */
/*  M[41].Xp[31].x = 	0.014	;	 M[41].Xp[31].y =18.8	; */
/*  M[41].Xp[32].x = 	0.012	;	 M[41].Xp[32].y =19.4	; */
/*  M[41].Xp[33].x = 	0.015	;	 M[41].Xp[33].y =20	; */
/*  M[41].Xp[34].x = 	0.005	;	 M[41].Xp[34].y =20.6	; */
/*  M[41].Xp[35].x = 	0.004	;	 M[41].Xp[35].y =21.2	; */
/*  M[41].Xp[36].x = 	0.003	;	 M[41].Xp[36].y =21.8	; */
/*  M[41].Xp[37].x = 	0.004	;	 M[41].Xp[37].y =22.4	; */
/*  M[41].Xp[38].x = 	0.001	;	 M[41].Xp[38].y =23	; */
/*  M[41].Xp[39].x = 	0.004	;	 M[41].Xp[39].y =23.6	; */
/*  M[41].Xp[40].x = 	0.006	;	 M[41].Xp[40].y =24.2	; */
/*  M[41].Xp[41].x = 	0.006	;	 M[41].Xp[41].y =24.8	; */
/*  M[41].Xp[42].x = 	0.006	;	 M[41].Xp[42].y =25.4	; */
/*  M[41].Xp[43].x = 	0	;	 M[41].Xp[43].y =26	; */
/*  M[41].pathlength(); */
/* //station 43 */
/*  x[42] = 59.5; */
/*  M[42].allocate(44); */
/*  M[42].settype(SPLINE2D_LINEAR); */
/*  M[42].Xp[0].x = 	0.154	;	 M[42].Xp[0].y =0.2	; */
/*  M[42].Xp[1].x = 	0.137	;	 M[42].Xp[1].y =0.8	; */
/*  M[42].Xp[2].x = 	0.126	;	 M[42].Xp[2].y =1.4	; */
/*  M[42].Xp[3].x = 	0.121	;	 M[42].Xp[3].y =2	; */
/*  M[42].Xp[4].x = 	0.15	;	 M[42].Xp[4].y =2.6	; */
/*  M[42].Xp[5].x = 	0.131	;	 M[42].Xp[5].y =3.2	; */
/*  M[42].Xp[6].x = 	0.121	;	 M[42].Xp[6].y =3.8	; */
/*  M[42].Xp[7].x = 	0.112	;	 M[42].Xp[7].y =4.4	; */
/*  M[42].Xp[8].x = 	0.1	;	 M[42].Xp[8].y =5	; */
/*  M[42].Xp[9].x = 	0.085	;	 M[42].Xp[9].y =5.6	; */
/*  M[42].Xp[10].x = 	0.078	;	 M[42].Xp[10].y =6.2	; */
/*  M[42].Xp[11].x = 	0.061	;	 M[42].Xp[11].y =6.8	; */
/*  M[42].Xp[12].x = 	0.05	;	 M[42].Xp[12].y =7.4	; */
/*  M[42].Xp[13].x = 	0.046	;	 M[42].Xp[13].y =8	; */
/*  M[42].Xp[14].x = 	0.054	;	 M[42].Xp[14].y =8.6	; */
/*  M[42].Xp[15].x = 	0.051	;	 M[42].Xp[15].y =9.2	; */
/*  M[42].Xp[16].x = 	0.046	;	 M[42].Xp[16].y =9.8	; */
/*  M[42].Xp[17].x = 	0.051	;	 M[42].Xp[17].y =10.4	; */
/*  M[42].Xp[18].x = 	0.038	;	 M[42].Xp[18].y =11	; */
/*  M[42].Xp[19].x = 	0.029	;	 M[42].Xp[19].y =11.6	; */
/*  M[42].Xp[20].x = 	0.035	;	 M[42].Xp[20].y =12.2	; */
/*  M[42].Xp[21].x = 	0.039	;	 M[42].Xp[21].y =12.8	; */
/*  M[42].Xp[22].x = 	0.027	;	 M[42].Xp[22].y =13.4	; */
/*  M[42].Xp[23].x = 	0.028	;	 M[42].Xp[23].y =14	; */
/*  M[42].Xp[24].x = 	0.025	;	 M[42].Xp[24].y =14.6	; */
/*  M[42].Xp[25].x = 	0.024	;	 M[42].Xp[25].y =15.2	; */
/*  M[42].Xp[26].x = 	0.023	;	 M[42].Xp[26].y =15.8	; */
/*  M[42].Xp[27].x = 	0.026	;	 M[42].Xp[27].y =16.4	; */
/*  M[42].Xp[28].x = 	0.019	;	 M[42].Xp[28].y =17	; */
/*  M[42].Xp[29].x = 	0.015	;	 M[42].Xp[29].y =17.6	; */
/*  M[42].Xp[30].x = 	0.014	;	 M[42].Xp[30].y =18.2	; */
/*  M[42].Xp[31].x = 	0.014	;	 M[42].Xp[31].y =18.8	; */
/*  M[42].Xp[32].x = 	0.01	;	 M[42].Xp[32].y =19.4	; */
/*  M[42].Xp[33].x = 	0.009	;	 M[42].Xp[33].y =20	; */
/*  M[42].Xp[34].x = 	0.001	;	 M[42].Xp[34].y =20.6	; */
/*  M[42].Xp[35].x = 	0.002	;	 M[42].Xp[35].y =21.2	; */
/*  M[42].Xp[36].x = 	0.003	;	 M[42].Xp[36].y =21.8	; */
/*  M[42].Xp[37].x = 	0.002	;	 M[42].Xp[37].y =22.4	; */
/*  M[42].Xp[38].x = 	0.001	;	 M[42].Xp[38].y =23	; */
/*  M[42].Xp[39].x = 	0.002	;	 M[42].Xp[39].y =23.6	; */
/*  M[42].Xp[40].x = 	0.001	;	 M[42].Xp[40].y =24.2	; */
/*  M[42].Xp[41].x = 	0.006	;	 M[42].Xp[41].y =24.8	; */
/*  M[42].Xp[42].x = 	0.004	;	 M[42].Xp[42].y =25.4	; */
/*  M[42].Xp[43].x = 	0.003	;	 M[42].Xp[43].y =26	; */
/*  M[42].pathlength(); */
/* //station 44 */
/*  x[43] = 60.3; */
/*  M[43].allocate(44); */
/*  M[43].settype(SPLINE2D_LINEAR); */
/*  M[43].Xp[0].x = 	0.12	;	 M[43].Xp[0].y =0.2	; */
/*  M[43].Xp[1].x = 	0.135	;	 M[43].Xp[1].y =0.8	; */
/*  M[43].Xp[2].x = 	0.124	;	 M[43].Xp[2].y =1.4	; */
/*  M[43].Xp[3].x = 	0.128	;	 M[43].Xp[3].y =2	; */
/*  M[43].Xp[4].x = 	0.131	;	 M[43].Xp[4].y =2.6	; */
/*  M[43].Xp[5].x = 	0.132	;	 M[43].Xp[5].y =3.2	; */
/*  M[43].Xp[6].x = 	0.102	;	 M[43].Xp[6].y =3.8	; */
/*  M[43].Xp[7].x = 	0.105	;	 M[43].Xp[7].y =4.4	; */
/*  M[43].Xp[8].x = 	0.084	;	 M[43].Xp[8].y =5	; */
/*  M[43].Xp[9].x = 	0.077	;	 M[43].Xp[9].y =5.6	; */
/*  M[43].Xp[10].x = 	0.082	;	 M[43].Xp[10].y =6.2	; */
/*  M[43].Xp[11].x = 	0.062	;	 M[43].Xp[11].y =6.8	; */
/*  M[43].Xp[12].x = 	0.06	;	 M[43].Xp[12].y =7.4	; */
/*  M[43].Xp[13].x = 	0.043	;	 M[43].Xp[13].y =8	; */
/*  M[43].Xp[14].x = 	0.057	;	 M[43].Xp[14].y =8.6	; */
/*  M[43].Xp[15].x = 	0.053	;	 M[43].Xp[15].y =9.2	; */
/*  M[43].Xp[16].x = 	0.039	;	 M[43].Xp[16].y =9.8	; */
/*  M[43].Xp[17].x = 	0.043	;	 M[43].Xp[17].y =10.4	; */
/*  M[43].Xp[18].x = 	0.042	;	 M[43].Xp[18].y =11	; */
/*  M[43].Xp[19].x = 	0.039	;	 M[43].Xp[19].y =11.6	; */
/*  M[43].Xp[20].x = 	0.037	;	 M[43].Xp[20].y =12.2	; */
/*  M[43].Xp[21].x = 	0.034	;	 M[43].Xp[21].y =12.8	; */
/*  M[43].Xp[22].x = 	0.03	;	 M[43].Xp[22].y =13.4	; */
/*  M[43].Xp[23].x = 	0.023	;	 M[43].Xp[23].y =14	; */
/*  M[43].Xp[24].x = 	0.026	;	 M[43].Xp[24].y =14.6	; */
/*  M[43].Xp[25].x = 	0.02	;	 M[43].Xp[25].y =15.2	; */
/*  M[43].Xp[26].x = 	0.021	;	 M[43].Xp[26].y =15.8	; */
/*  M[43].Xp[27].x = 	0.022	;	 M[43].Xp[27].y =16.4	; */
/*  M[43].Xp[28].x = 	0.017	;	 M[43].Xp[28].y =17	; */
/*  M[43].Xp[29].x = 	0.012	;	 M[43].Xp[29].y =17.6	; */
/*  M[43].Xp[30].x = 	0.011	;	 M[43].Xp[30].y =18.2	; */
/*  M[43].Xp[31].x = 	0.013	;	 M[43].Xp[31].y =18.8	; */
/*  M[43].Xp[32].x = 	0.007	;	 M[43].Xp[32].y =19.4	; */
/*  M[43].Xp[33].x = 	0.007	;	 M[43].Xp[33].y =20	; */
/*  M[43].Xp[34].x = 	0.004	;	 M[43].Xp[34].y =20.6	; */
/*  M[43].Xp[35].x = 	0.007	;	 M[43].Xp[35].y =21.2	; */
/*  M[43].Xp[36].x = 	0.005	;	 M[43].Xp[36].y =21.8	; */
/*  M[43].Xp[37].x = 	0.008	;	 M[43].Xp[37].y =22.4	; */
/*  M[43].Xp[38].x = 	0.003	;	 M[43].Xp[38].y =23	; */
/*  M[43].Xp[39].x = 	0.004	;	 M[43].Xp[39].y =23.6	; */
/*  M[43].Xp[40].x = 	0.004	;	 M[43].Xp[40].y =24.2	; */
/*  M[43].Xp[41].x = 	0.004	;	 M[43].Xp[41].y =24.8	; */
/*  M[43].Xp[42].x = 	0	;	 M[43].Xp[42].y =25.4	; */
/*  M[43].Xp[43].x = 	0.004	;	 M[43].Xp[43].y =26	; */
/*  M[43].pathlength(); */
/* //station 45 */
/*  x[44] = 61.0; */
/*  M[44].allocate(44); */
/*  M[44].settype(SPLINE2D_LINEAR); */
/*  M[44].Xp[0].x = 	0.117	;	 M[44].Xp[0].y =0.2	; */
/*  M[44].Xp[1].x = 	0.128	;	 M[44].Xp[1].y =0.8	; */
/*  M[44].Xp[2].x = 	0.124	;	 M[44].Xp[2].y =1.4	; */
/*  M[44].Xp[3].x = 	0.134	;	 M[44].Xp[3].y =2	; */
/*  M[44].Xp[4].x = 	0.109	;	 M[44].Xp[4].y =2.6	; */
/*  M[44].Xp[5].x = 	0.111	;	 M[44].Xp[5].y =3.2	; */
/*  M[44].Xp[6].x = 	0.108	;	 M[44].Xp[6].y =3.8	; */
/*  M[44].Xp[7].x = 	0.08	;	 M[44].Xp[7].y =4.4	; */
/*  M[44].Xp[8].x = 	0.083	;	 M[44].Xp[8].y =5	; */
/*  M[44].Xp[9].x = 	0.079	;	 M[44].Xp[9].y =5.6	; */
/*  M[44].Xp[10].x = 	0.072	;	 M[44].Xp[10].y =6.2	; */
/*  M[44].Xp[11].x = 	0.064	;	 M[44].Xp[11].y =6.8	; */
/*  M[44].Xp[12].x = 	0.06	;	 M[44].Xp[12].y =7.4	; */
/*  M[44].Xp[13].x = 	0.05	;	 M[44].Xp[13].y =8	; */
/*  M[44].Xp[14].x = 	0.055	;	 M[44].Xp[14].y =8.6	; */
/*  M[44].Xp[15].x = 	0.045	;	 M[44].Xp[15].y =9.2	; */
/*  M[44].Xp[16].x = 	0.047	;	 M[44].Xp[16].y =9.8	; */
/*  M[44].Xp[17].x = 	0.04	;	 M[44].Xp[17].y =10.4	; */
/*  M[44].Xp[18].x = 	0.032	;	 M[44].Xp[18].y =11	; */
/*  M[44].Xp[19].x = 	0.046	;	 M[44].Xp[19].y =11.6	; */
/*  M[44].Xp[20].x = 	0.038	;	 M[44].Xp[20].y =12.2	; */
/*  M[44].Xp[21].x = 	0.03	;	 M[44].Xp[21].y =12.8	; */
/*  M[44].Xp[22].x = 	0.029	;	 M[44].Xp[22].y =13.4	; */
/*  M[44].Xp[23].x = 	0.026	;	 M[44].Xp[23].y =14	; */
/*  M[44].Xp[24].x = 	0.019	;	 M[44].Xp[24].y =14.6	; */
/*  M[44].Xp[25].x = 	0.02	;	 M[44].Xp[25].y =15.2	; */
/*  M[44].Xp[26].x = 	0.022	;	 M[44].Xp[26].y =15.8	; */
/*  M[44].Xp[27].x = 	0.016	;	 M[44].Xp[27].y =16.4	; */
/*  M[44].Xp[28].x = 	0.013	;	 M[44].Xp[28].y =17	; */
/*  M[44].Xp[29].x = 	0.013	;	 M[44].Xp[29].y =17.6	; */
/*  M[44].Xp[30].x = 	0.01	;	 M[44].Xp[30].y =18.2	; */
/*  M[44].Xp[31].x = 	0.011	;	 M[44].Xp[31].y =18.8	; */
/*  M[44].Xp[32].x = 	0.011	;	 M[44].Xp[32].y =19.4	; */
/*  M[44].Xp[33].x = 	0.006	;	 M[44].Xp[33].y =20	; */
/*  M[44].Xp[34].x = 	0.005	;	 M[44].Xp[34].y =20.6	; */
/*  M[44].Xp[35].x = 	0.005	;	 M[44].Xp[35].y =21.2	; */
/*  M[44].Xp[36].x = 	0.004	;	 M[44].Xp[36].y =21.8	; */
/*  M[44].Xp[37].x = 	0.004	;	 M[44].Xp[37].y =22.4	; */
/*  M[44].Xp[38].x = 	0	;	 M[44].Xp[38].y =23	; */
/*  M[44].Xp[39].x = 	0	;	 M[44].Xp[39].y =23.6	; */
/*  M[44].Xp[40].x = 	0.001	;	 M[44].Xp[40].y =24.2	; */
/*  M[44].Xp[41].x = 	0.001	;	 M[44].Xp[41].y =24.8	; */
/*  M[44].Xp[42].x = 	0.002	;	 M[44].Xp[42].y =25.4	; */
/*  M[44].Xp[43].x = 	0	;	 M[44].Xp[43].y =26	; */
/*  M[44].pathlength(); */
/* //station 46 */
/*  x[45] = 61.8; */
/*  M[45].allocate(44); */
/*  M[45].settype(SPLINE2D_LINEAR); */
/*  M[45].Xp[0].x = 	0.089	;	 M[45].Xp[0].y =0.2	; */
/*  M[45].Xp[1].x = 	0.108	;	 M[45].Xp[1].y =0.8	; */
/*  M[45].Xp[2].x = 	0.108	;	 M[45].Xp[2].y =1.4	; */
/*  M[45].Xp[3].x = 	0.116	;	 M[45].Xp[3].y =2	; */
/*  M[45].Xp[4].x = 	0.13	;	 M[45].Xp[4].y =2.6	; */
/*  M[45].Xp[5].x = 	0.1	;	 M[45].Xp[5].y =3.2	; */
/*  M[45].Xp[6].x = 	0.083	;	 M[45].Xp[6].y =3.8	; */
/*  M[45].Xp[7].x = 	0.079	;	 M[45].Xp[7].y =4.4	; */
/*  M[45].Xp[8].x = 	0.085	;	 M[45].Xp[8].y =5	; */
/*  M[45].Xp[9].x = 	0.073	;	 M[45].Xp[9].y =5.6	; */
/*  M[45].Xp[10].x = 	0.06	;	 M[45].Xp[10].y =6.2	; */
/*  M[45].Xp[11].x = 	0.057	;	 M[45].Xp[11].y =6.8	; */
/*  M[45].Xp[12].x = 	0.054	;	 M[45].Xp[12].y =7.4	; */
/*  M[45].Xp[13].x = 	0.059	;	 M[45].Xp[13].y =8	; */
/*  M[45].Xp[14].x = 	0.047	;	 M[45].Xp[14].y =8.6	; */
/*  M[45].Xp[15].x = 	0.044	;	 M[45].Xp[15].y =9.2	; */
/*  M[45].Xp[16].x = 	0.035	;	 M[45].Xp[16].y =9.8	; */
/*  M[45].Xp[17].x = 	0.042	;	 M[45].Xp[17].y =10.4	; */
/*  M[45].Xp[18].x = 	0.035	;	 M[45].Xp[18].y =11	; */
/*  M[45].Xp[19].x = 	0.031	;	 M[45].Xp[19].y =11.6	; */
/*  M[45].Xp[20].x = 	0.034	;	 M[45].Xp[20].y =12.2	; */
/*  M[45].Xp[21].x = 	0.027	;	 M[45].Xp[21].y =12.8	; */
/*  M[45].Xp[22].x = 	0.028	;	 M[45].Xp[22].y =13.4	; */
/*  M[45].Xp[23].x = 	0.03	;	 M[45].Xp[23].y =14	; */
/*  M[45].Xp[24].x = 	0.022	;	 M[45].Xp[24].y =14.6	; */
/*  M[45].Xp[25].x = 	0.026	;	 M[45].Xp[25].y =15.2	; */
/*  M[45].Xp[26].x = 	0.015	;	 M[45].Xp[26].y =15.8	; */
/*  M[45].Xp[27].x = 	0.022	;	 M[45].Xp[27].y =16.4	; */
/*  M[45].Xp[28].x = 	0.012	;	 M[45].Xp[28].y =17	; */
/*  M[45].Xp[29].x = 	0.01	;	 M[45].Xp[29].y =17.6	; */
/*  M[45].Xp[30].x = 	0.009	;	 M[45].Xp[30].y =18.2	; */
/*  M[45].Xp[31].x = 	0.006	;	 M[45].Xp[31].y =18.8	; */
/*  M[45].Xp[32].x = 	0.009	;	 M[45].Xp[32].y =19.4	; */
/*  M[45].Xp[33].x = 	0.004	;	 M[45].Xp[33].y =20	; */
/*  M[45].Xp[34].x = 	0.006	;	 M[45].Xp[34].y =20.6	; */
/*  M[45].Xp[35].x = 	0.003	;	 M[45].Xp[35].y =21.2	; */
/*  M[45].Xp[36].x = 	0.004	;	 M[45].Xp[36].y =21.8	; */
/*  M[45].Xp[37].x = 	0.003	;	 M[45].Xp[37].y =22.4	; */
/*  M[45].Xp[38].x = 	0.002	;	 M[45].Xp[38].y =23	; */
/*  M[45].Xp[39].x = 	0.003	;	 M[45].Xp[39].y =23.6	; */
/*  M[45].Xp[40].x = 	0.002	;	 M[45].Xp[40].y =24.2	; */
/*  M[45].Xp[41].x = 	0.001	;	 M[45].Xp[41].y =24.8	; */
/*  M[45].Xp[42].x = 	0.003	;	 M[45].Xp[42].y =25.4	; */
/*  M[45].Xp[43].x = 	0	;	 M[45].Xp[43].y =26	; */
/*  M[45].pathlength(); */
/* //station 47 */
/*  x[46] = 62.5; */
/*  M[46].allocate(44); */
/*  M[46].settype(SPLINE2D_LINEAR); */
/*  M[46].Xp[0].x = 	0.107	;	 M[46].Xp[0].y =0.2	; */
/*  M[46].Xp[1].x = 	0.098	;	 M[46].Xp[1].y =0.8	; */
/*  M[46].Xp[2].x = 	0.102	;	 M[46].Xp[2].y =1.4	; */
/*  M[46].Xp[3].x = 	0.112	;	 M[46].Xp[3].y =2	; */
/*  M[46].Xp[4].x = 	0.097	;	 M[46].Xp[4].y =2.6	; */
/*  M[46].Xp[5].x = 	0.087	;	 M[46].Xp[5].y =3.2	; */
/*  M[46].Xp[6].x = 	0.099	;	 M[46].Xp[6].y =3.8	; */
/*  M[46].Xp[7].x = 	0.087	;	 M[46].Xp[7].y =4.4	; */
/*  M[46].Xp[8].x = 	0.082	;	 M[46].Xp[8].y =5	; */
/*  M[46].Xp[9].x = 	0.066	;	 M[46].Xp[9].y =5.6	; */
/*  M[46].Xp[10].x = 	0.071	;	 M[46].Xp[10].y =6.2	; */
/*  M[46].Xp[11].x = 	0.066	;	 M[46].Xp[11].y =6.8	; */
/*  M[46].Xp[12].x = 	0.077	;	 M[46].Xp[12].y =7.4	; */
/*  M[46].Xp[13].x = 	0.056	;	 M[46].Xp[13].y =8	; */
/*  M[46].Xp[14].x = 	0.05	;	 M[46].Xp[14].y =8.6	; */
/*  M[46].Xp[15].x = 	0.047	;	 M[46].Xp[15].y =9.2	; */
/*  M[46].Xp[16].x = 	0.048	;	 M[46].Xp[16].y =9.8	; */
/*  M[46].Xp[17].x = 	0.046	;	 M[46].Xp[17].y =10.4	; */
/*  M[46].Xp[18].x = 	0.039	;	 M[46].Xp[18].y =11	; */
/*  M[46].Xp[19].x = 	0.045	;	 M[46].Xp[19].y =11.6	; */
/*  M[46].Xp[20].x = 	0.043	;	 M[46].Xp[20].y =12.2	; */
/*  M[46].Xp[21].x = 	0.029	;	 M[46].Xp[21].y =12.8	; */
/*  M[46].Xp[22].x = 	0.032	;	 M[46].Xp[22].y =13.4	; */
/*  M[46].Xp[23].x = 	0.033	;	 M[46].Xp[23].y =14	; */
/*  M[46].Xp[24].x = 	0.03	;	 M[46].Xp[24].y =14.6	; */
/*  M[46].Xp[25].x = 	0.027	;	 M[46].Xp[25].y =15.2	; */
/*  M[46].Xp[26].x = 	0.022	;	 M[46].Xp[26].y =15.8	; */
/*  M[46].Xp[27].x = 	0.022	;	 M[46].Xp[27].y =16.4	; */
/*  M[46].Xp[28].x = 	0.016	;	 M[46].Xp[28].y =17	; */
/*  M[46].Xp[29].x = 	0.007	;	 M[46].Xp[29].y =17.6	; */
/*  M[46].Xp[30].x = 	0.009	;	 M[46].Xp[30].y =18.2	; */
/*  M[46].Xp[31].x = 	0.009	;	 M[46].Xp[31].y =18.8	; */
/*  M[46].Xp[32].x = 	0.004	;	 M[46].Xp[32].y =19.4	; */
/*  M[46].Xp[33].x = 	0.005	;	 M[46].Xp[33].y =20	; */
/*  M[46].Xp[34].x = 	0.007	;	 M[46].Xp[34].y =20.6	; */
/*  M[46].Xp[35].x = 	0.006	;	 M[46].Xp[35].y =21.2	; */
/*  M[46].Xp[36].x = 	0.007	;	 M[46].Xp[36].y =21.8	; */
/*  M[46].Xp[37].x = 	0.004	;	 M[46].Xp[37].y =22.4	; */
/*  M[46].Xp[38].x = 	0.007	;	 M[46].Xp[38].y =23	; */
/*  M[46].Xp[39].x = 	0.002	;	 M[46].Xp[39].y =23.6	; */
/*  M[46].Xp[40].x = 	0.003	;	 M[46].Xp[40].y =24.2	; */
/*  M[46].Xp[41].x = 	0.004	;	 M[46].Xp[41].y =24.8	; */
/*  M[46].Xp[42].x = 	0.004	;	 M[46].Xp[42].y =25.4	; */
/*  M[46].Xp[43].x = 	0.002	;	 M[46].Xp[43].y =26	; */
/*  M[46].pathlength(); */
/* //station 48 */
/*  x[47] = 63.3; */
/*  M[47].allocate(44); */
/*  M[47].settype(SPLINE2D_LINEAR); */
/*  M[47].Xp[0].x = 	0.094	;	 M[47].Xp[0].y =0.2	; */
/*  M[47].Xp[1].x = 	0.091	;	 M[47].Xp[1].y =0.8	; */
/*  M[47].Xp[2].x = 	0.081	;	 M[47].Xp[2].y =1.4	; */
/*  M[47].Xp[3].x = 	0.088	;	 M[47].Xp[3].y =2	; */
/*  M[47].Xp[4].x = 	0.079	;	 M[47].Xp[4].y =2.6	; */
/*  M[47].Xp[5].x = 	0.089	;	 M[47].Xp[5].y =3.2	; */
/*  M[47].Xp[6].x = 	0.075	;	 M[47].Xp[6].y =3.8	; */
/*  M[47].Xp[7].x = 	0.067	;	 M[47].Xp[7].y =4.4	; */
/*  M[47].Xp[8].x = 	0.071	;	 M[47].Xp[8].y =5	; */
/*  M[47].Xp[9].x = 	0.067	;	 M[47].Xp[9].y =5.6	; */
/*  M[47].Xp[10].x = 	0.066	;	 M[47].Xp[10].y =6.2	; */
/*  M[47].Xp[11].x = 	0.064	;	 M[47].Xp[11].y =6.8	; */
/*  M[47].Xp[12].x = 	0.069	;	 M[47].Xp[12].y =7.4	; */
/*  M[47].Xp[13].x = 	0.062	;	 M[47].Xp[13].y =8	; */
/*  M[47].Xp[14].x = 	0.051	;	 M[47].Xp[14].y =8.6	; */
/*  M[47].Xp[15].x = 	0.042	;	 M[47].Xp[15].y =9.2	; */
/*  M[47].Xp[16].x = 	0.045	;	 M[47].Xp[16].y =9.8	; */
/*  M[47].Xp[17].x = 	0.046	;	 M[47].Xp[17].y =10.4	; */
/*  M[47].Xp[18].x = 	0.04	;	 M[47].Xp[18].y =11	; */
/*  M[47].Xp[19].x = 	0.045	;	 M[47].Xp[19].y =11.6	; */
/*  M[47].Xp[20].x = 	0.037	;	 M[47].Xp[20].y =12.2	; */
/*  M[47].Xp[21].x = 	0.027	;	 M[47].Xp[21].y =12.8	; */
/*  M[47].Xp[22].x = 	0.031	;	 M[47].Xp[22].y =13.4	; */
/*  M[47].Xp[23].x = 	0.026	;	 M[47].Xp[23].y =14	; */
/*  M[47].Xp[24].x = 	0.023	;	 M[47].Xp[24].y =14.6	; */
/*  M[47].Xp[25].x = 	0.025	;	 M[47].Xp[25].y =15.2	; */
/*  M[47].Xp[26].x = 	0.023	;	 M[47].Xp[26].y =15.8	; */
/*  M[47].Xp[27].x = 	0.016	;	 M[47].Xp[27].y =16.4	; */
/*  M[47].Xp[28].x = 	0.016	;	 M[47].Xp[28].y =17	; */
/*  M[47].Xp[29].x = 	0.012	;	 M[47].Xp[29].y =17.6	; */
/*  M[47].Xp[30].x = 	0.01	;	 M[47].Xp[30].y =18.2	; */
/*  M[47].Xp[31].x = 	0.01	;	 M[47].Xp[31].y =18.8	; */
/*  M[47].Xp[32].x = 	0.009	;	 M[47].Xp[32].y =19.4	; */
/*  M[47].Xp[33].x = 	0.005	;	 M[47].Xp[33].y =20	; */
/*  M[47].Xp[34].x = 	0.003	;	 M[47].Xp[34].y =20.6	; */
/*  M[47].Xp[35].x = 	0.005	;	 M[47].Xp[35].y =21.2	; */
/*  M[47].Xp[36].x = 	0.006	;	 M[47].Xp[36].y =21.8	; */
/*  M[47].Xp[37].x = 	0.004	;	 M[47].Xp[37].y =22.4	; */
/*  M[47].Xp[38].x = 	0.004	;	 M[47].Xp[38].y =23	; */
/*  M[47].Xp[39].x = 	0.006	;	 M[47].Xp[39].y =23.6	; */
/*  M[47].Xp[40].x = 	0.004	;	 M[47].Xp[40].y =24.2	; */
/*  M[47].Xp[41].x = 	0.006	;	 M[47].Xp[41].y =24.8	; */
/*  M[47].Xp[42].x = 	0.002	;	 M[47].Xp[42].y =25.4	; */
/*  M[47].Xp[43].x = 	0.004	;	 M[47].Xp[43].y =26	; */
/*  M[47].pathlength(); */
/* //station 49 */
/*  x[48] = 64.0; */
/*  M[48].allocate(44); */
/*  M[48].settype(SPLINE2D_LINEAR); */
/*  M[48].Xp[0].x = 	0.073	;	 M[48].Xp[0].y =0.2	; */
/*  M[48].Xp[1].x = 	0.074	;	 M[48].Xp[1].y =0.8	; */
/*  M[48].Xp[2].x = 	0.094	;	 M[48].Xp[2].y =1.4	; */
/*  M[48].Xp[3].x = 	0.074	;	 M[48].Xp[3].y =2	; */
/*  M[48].Xp[4].x = 	0.072	;	 M[48].Xp[4].y =2.6	; */
/*  M[48].Xp[5].x = 	0.088	;	 M[48].Xp[5].y =3.2	; */
/*  M[48].Xp[6].x = 	0.076	;	 M[48].Xp[6].y =3.8	; */
/*  M[48].Xp[7].x = 	0.07	;	 M[48].Xp[7].y =4.4	; */
/*  M[48].Xp[8].x = 	0.08	;	 M[48].Xp[8].y =5	; */
/*  M[48].Xp[9].x = 	0.075	;	 M[48].Xp[9].y =5.6	; */
/*  M[48].Xp[10].x = 	0.067	;	 M[48].Xp[10].y =6.2	; */
/*  M[48].Xp[11].x = 	0.056	;	 M[48].Xp[11].y =6.8	; */
/*  M[48].Xp[12].x = 	0.053	;	 M[48].Xp[12].y =7.4	; */
/*  M[48].Xp[13].x = 	0.051	;	 M[48].Xp[13].y =8	; */
/*  M[48].Xp[14].x = 	0.048	;	 M[48].Xp[14].y =8.6	; */
/*  M[48].Xp[15].x = 	0.038	;	 M[48].Xp[15].y =9.2	; */
/*  M[48].Xp[16].x = 	0.038	;	 M[48].Xp[16].y =9.8	; */
/*  M[48].Xp[17].x = 	0.036	;	 M[48].Xp[17].y =10.4	; */
/*  M[48].Xp[18].x = 	0.032	;	 M[48].Xp[18].y =11	; */
/*  M[48].Xp[19].x = 	0.032	;	 M[48].Xp[19].y =11.6	; */
/*  M[48].Xp[20].x = 	0.026	;	 M[48].Xp[20].y =12.2	; */
/*  M[48].Xp[21].x = 	0.024	;	 M[48].Xp[21].y =12.8	; */
/*  M[48].Xp[22].x = 	0.023	;	 M[48].Xp[22].y =13.4	; */
/*  M[48].Xp[23].x = 	0.029	;	 M[48].Xp[23].y =14	; */
/*  M[48].Xp[24].x = 	0.029	;	 M[48].Xp[24].y =14.6	; */
/*  M[48].Xp[25].x = 	0.022	;	 M[48].Xp[25].y =15.2	; */
/*  M[48].Xp[26].x = 	0.024	;	 M[48].Xp[26].y =15.8	; */
/*  M[48].Xp[27].x = 	0.021	;	 M[48].Xp[27].y =16.4	; */
/*  M[48].Xp[28].x = 	0.014	;	 M[48].Xp[28].y =17	; */
/*  M[48].Xp[29].x = 	0.013	;	 M[48].Xp[29].y =17.6	; */
/*  M[48].Xp[30].x = 	0.015	;	 M[48].Xp[30].y =18.2	; */
/*  M[48].Xp[31].x = 	0.007	;	 M[48].Xp[31].y =18.8	; */
/*  M[48].Xp[32].x = 	0.007	;	 M[48].Xp[32].y =19.4	; */
/*  M[48].Xp[33].x = 	0.003	;	 M[48].Xp[33].y =20	; */
/*  M[48].Xp[34].x = 	0.007	;	 M[48].Xp[34].y =20.6	; */
/*  M[48].Xp[35].x = 	0.005	;	 M[48].Xp[35].y =21.2	; */
/*  M[48].Xp[36].x = 	0.004	;	 M[48].Xp[36].y =21.8	; */
/*  M[48].Xp[37].x = 	0.002	;	 M[48].Xp[37].y =22.4	; */
/*  M[48].Xp[38].x = 	0.002	;	 M[48].Xp[38].y =23	; */
/*  M[48].Xp[39].x = 	0.003	;	 M[48].Xp[39].y =23.6	; */
/*  M[48].Xp[40].x = 	0.002	;	 M[48].Xp[40].y =24.2	; */
/*  M[48].Xp[41].x = 	0	;	 M[48].Xp[41].y =24.8	; */
/*  M[48].Xp[42].x = 	0.003	;	 M[48].Xp[42].y =25.4	; */
/*  M[48].Xp[43].x = 	0.001	;	 M[48].Xp[43].y =26	; */
/*  M[48].pathlength(); */
/* //station 50 */
/*  x[49] = 64.8; */
/*  M[49].allocate(44); */
/*  M[49].settype(SPLINE2D_LINEAR); */
/*  M[49].Xp[0].x = 	0.072	;	 M[49].Xp[0].y =0.2	; */
/*  M[49].Xp[1].x = 	0.06	;	 M[49].Xp[1].y =0.8	; */
/*  M[49].Xp[2].x = 	0.072	;	 M[49].Xp[2].y =1.4	; */
/*  M[49].Xp[3].x = 	0.056	;	 M[49].Xp[3].y =2	; */
/*  M[49].Xp[4].x = 	0.07	;	 M[49].Xp[4].y =2.6	; */
/*  M[49].Xp[5].x = 	0.075	;	 M[49].Xp[5].y =3.2	; */
/*  M[49].Xp[6].x = 	0.067	;	 M[49].Xp[6].y =3.8	; */
/*  M[49].Xp[7].x = 	0.069	;	 M[49].Xp[7].y =4.4	; */
/*  M[49].Xp[8].x = 	0.065	;	 M[49].Xp[8].y =5	; */
/*  M[49].Xp[9].x = 	0.054	;	 M[49].Xp[9].y =5.6	; */
/*  M[49].Xp[10].x = 	0.054	;	 M[49].Xp[10].y =6.2	; */
/*  M[49].Xp[11].x = 	0.052	;	 M[49].Xp[11].y =6.8	; */
/*  M[49].Xp[12].x = 	0.049	;	 M[49].Xp[12].y =7.4	; */
/*  M[49].Xp[13].x = 	0.051	;	 M[49].Xp[13].y =8	; */
/*  M[49].Xp[14].x = 	0.051	;	 M[49].Xp[14].y =8.6	; */
/*  M[49].Xp[15].x = 	0.05	;	 M[49].Xp[15].y =9.2	; */
/*  M[49].Xp[16].x = 	0.039	;	 M[49].Xp[16].y =9.8	; */
/*  M[49].Xp[17].x = 	0.038	;	 M[49].Xp[17].y =10.4	; */
/*  M[49].Xp[18].x = 	0.031	;	 M[49].Xp[18].y =11	; */
/*  M[49].Xp[19].x = 	0.03	;	 M[49].Xp[19].y =11.6	; */
/*  M[49].Xp[20].x = 	0.024	;	 M[49].Xp[20].y =12.2	; */
/*  M[49].Xp[21].x = 	0.023	;	 M[49].Xp[21].y =12.8	; */
/*  M[49].Xp[22].x = 	0.027	;	 M[49].Xp[22].y =13.4	; */
/*  M[49].Xp[23].x = 	0.029	;	 M[49].Xp[23].y =14	; */
/*  M[49].Xp[24].x = 	0.02	;	 M[49].Xp[24].y =14.6	; */
/*  M[49].Xp[25].x = 	0.02	;	 M[49].Xp[25].y =15.2	; */
/*  M[49].Xp[26].x = 	0.016	;	 M[49].Xp[26].y =15.8	; */
/*  M[49].Xp[27].x = 	0.016	;	 M[49].Xp[27].y =16.4	; */
/*  M[49].Xp[28].x = 	0.016	;	 M[49].Xp[28].y =17	; */
/*  M[49].Xp[29].x = 	0.011	;	 M[49].Xp[29].y =17.6	; */
/*  M[49].Xp[30].x = 	0.014	;	 M[49].Xp[30].y =18.2	; */
/*  M[49].Xp[31].x = 	0.007	;	 M[49].Xp[31].y =18.8	; */
/*  M[49].Xp[32].x = 	0.009	;	 M[49].Xp[32].y =19.4	; */
/*  M[49].Xp[33].x = 	0.012	;	 M[49].Xp[33].y =20	; */
/*  M[49].Xp[34].x = 	0.003	;	 M[49].Xp[34].y =20.6	; */
/*  M[49].Xp[35].x = 	0.001	;	 M[49].Xp[35].y =21.2	; */
/*  M[49].Xp[36].x = 	0.005	;	 M[49].Xp[36].y =21.8	; */
/*  M[49].Xp[37].x = 	0.003	;	 M[49].Xp[37].y =22.4	; */
/*  M[49].Xp[38].x = 	0.006	;	 M[49].Xp[38].y =23	; */
/*  M[49].Xp[39].x = 	0.001	;	 M[49].Xp[39].y =23.6	; */
/*  M[49].Xp[40].x = 	0.002	;	 M[49].Xp[40].y =24.2	; */
/*  M[49].Xp[41].x = 	0.001	;	 M[49].Xp[41].y =24.8	; */
/*  M[49].Xp[42].x = 	0.001	;	 M[49].Xp[42].y =25.4	; */
/*  M[49].Xp[43].x = 	0.001	;	 M[49].Xp[43].y =26	; */
/*  M[49].pathlength(); */
/* //station 51 */
/*  x[50] = 65.5; */
/*  M[50].allocate(44); */
/*  M[50].settype(SPLINE2D_LINEAR); */
/*  M[50].Xp[0].x = 	0.087	;	 M[50].Xp[0].y =0.2	; */
/*  M[50].Xp[1].x = 	0.085	;	 M[50].Xp[1].y =0.8	; */
/*  M[50].Xp[2].x = 	0.069	;	 M[50].Xp[2].y =1.4	; */
/*  M[50].Xp[3].x = 	0.053	;	 M[50].Xp[3].y =2	; */
/*  M[50].Xp[4].x = 	0.076	;	 M[50].Xp[4].y =2.6	; */
/*  M[50].Xp[5].x = 	0.079	;	 M[50].Xp[5].y =3.2	; */
/*  M[50].Xp[6].x = 	0.074	;	 M[50].Xp[6].y =3.8	; */
/*  M[50].Xp[7].x = 	0.082	;	 M[50].Xp[7].y =4.4	; */
/*  M[50].Xp[8].x = 	0.084	;	 M[50].Xp[8].y =5	; */
/*  M[50].Xp[9].x = 	0.066	;	 M[50].Xp[9].y =5.6	; */
/*  M[50].Xp[10].x = 	0.068	;	 M[50].Xp[10].y =6.2	; */
/*  M[50].Xp[11].x = 	0.065	;	 M[50].Xp[11].y =6.8	; */
/*  M[50].Xp[12].x = 	0.047	;	 M[50].Xp[12].y =7.4	; */
/*  M[50].Xp[13].x = 	0.053	;	 M[50].Xp[13].y =8	; */
/*  M[50].Xp[14].x = 	0.032	;	 M[50].Xp[14].y =8.6	; */
/*  M[50].Xp[15].x = 	0.045	;	 M[50].Xp[15].y =9.2	; */
/*  M[50].Xp[16].x = 	0.039	;	 M[50].Xp[16].y =9.8	; */
/*  M[50].Xp[17].x = 	0.036	;	 M[50].Xp[17].y =10.4	; */
/*  M[50].Xp[18].x = 	0.038	;	 M[50].Xp[18].y =11	; */
/*  M[50].Xp[19].x = 	0.033	;	 M[50].Xp[19].y =11.6	; */
/*  M[50].Xp[20].x = 	0.018	;	 M[50].Xp[20].y =12.2	; */
/*  M[50].Xp[21].x = 	0.02	;	 M[50].Xp[21].y =12.8	; */
/*  M[50].Xp[22].x = 	0.027	;	 M[50].Xp[22].y =13.4	; */
/*  M[50].Xp[23].x = 	0.025	;	 M[50].Xp[23].y =14	; */
/*  M[50].Xp[24].x = 	0.028	;	 M[50].Xp[24].y =14.6	; */
/*  M[50].Xp[25].x = 	0.018	;	 M[50].Xp[25].y =15.2	; */
/*  M[50].Xp[26].x = 	0.031	;	 M[50].Xp[26].y =15.8	; */
/*  M[50].Xp[27].x = 	0.021	;	 M[50].Xp[27].y =16.4	; */
/*  M[50].Xp[28].x = 	0.021	;	 M[50].Xp[28].y =17	; */
/*  M[50].Xp[29].x = 	0.013	;	 M[50].Xp[29].y =17.6	; */
/*  M[50].Xp[30].x = 	0.012	;	 M[50].Xp[30].y =18.2	; */
/*  M[50].Xp[31].x = 	0.004	;	 M[50].Xp[31].y =18.8	; */
/*  M[50].Xp[32].x = 	0.009	;	 M[50].Xp[32].y =19.4	; */
/*  M[50].Xp[33].x = 	0.008	;	 M[50].Xp[33].y =20	; */
/*  M[50].Xp[34].x = 	0.004	;	 M[50].Xp[34].y =20.6	; */
/*  M[50].Xp[35].x = 	0.004	;	 M[50].Xp[35].y =21.2	; */
/*  M[50].Xp[36].x = 	0.004	;	 M[50].Xp[36].y =21.8	; */
/*  M[50].Xp[37].x = 	0.004	;	 M[50].Xp[37].y =22.4	; */
/*  M[50].Xp[38].x = 	0.009	;	 M[50].Xp[38].y =23	; */
/*  M[50].Xp[39].x = 	0.002	;	 M[50].Xp[39].y =23.6	; */
/*  M[50].Xp[40].x = 	0.002	;	 M[50].Xp[40].y =24.2	; */
/*  M[50].Xp[41].x = 	0.005	;	 M[50].Xp[41].y =24.8	; */
/*  M[50].Xp[42].x = 	0.002	;	 M[50].Xp[42].y =25.4	; */
/*  M[50].Xp[43].x = 	0.003	;	 M[50].Xp[43].y =26	; */
/*  M[50].pathlength(); */

/*  x[51] = 66.3; */
/*  M[51].allocate(44); */
/*  M[51].settype(SPLINE2D_LINEAR); */
/*  M[51].Xp[0].x = 	0.08	;	 M[51].Xp[0].y =0.2	; */
/*  M[51].Xp[1].x = 	0.068	;	 M[51].Xp[1].y =0.8	; */
/*  M[51].Xp[2].x = 	0.059	;	 M[51].Xp[2].y =1.4	; */
/*  M[51].Xp[3].x = 	0.056	;	 M[51].Xp[3].y =2	; */
/*  M[51].Xp[4].x = 	0.089	;	 M[51].Xp[4].y =2.6	; */
/*  M[51].Xp[5].x = 	0.09	;	 M[51].Xp[5].y =3.2	; */
/*  M[51].Xp[6].x = 	0.09	;	 M[51].Xp[6].y =3.8	; */
/*  M[51].Xp[7].x = 	0.085	;	 M[51].Xp[7].y =4.4	; */
/*  M[51].Xp[8].x = 	0.066	;	 M[51].Xp[8].y =5	; */
/*  M[51].Xp[9].x = 	0.063	;	 M[51].Xp[9].y =5.6	; */
/*  M[51].Xp[10].x = 	0.065	;	 M[51].Xp[10].y =6.2	; */
/*  M[51].Xp[11].x = 	0.072	;	 M[51].Xp[11].y =6.8	; */
/*  M[51].Xp[12].x = 	0.038	;	 M[51].Xp[12].y =7.4	; */
/*  M[51].Xp[13].x = 	0.06	;	 M[51].Xp[13].y =8	; */
/*  M[51].Xp[14].x = 	0.037	;	 M[51].Xp[14].y =8.6	; */
/*  M[51].Xp[15].x = 	0.051	;	 M[51].Xp[15].y =9.2	; */
/*  M[51].Xp[16].x = 	0.032	;	 M[51].Xp[16].y =9.8	; */
/*  M[51].Xp[17].x = 	0.035	;	 M[51].Xp[17].y =10.4	; */
/*  M[51].Xp[18].x = 	0.031	;	 M[51].Xp[18].y =11	; */
/*  M[51].Xp[19].x = 	0.025	;	 M[51].Xp[19].y =11.6	; */
/*  M[51].Xp[20].x = 	0.029	;	 M[51].Xp[20].y =12.2	; */
/*  M[51].Xp[21].x = 	0.02	;	 M[51].Xp[21].y =12.8	; */
/*  M[51].Xp[22].x = 	0.016	;	 M[51].Xp[22].y =13.4	; */
/*  M[51].Xp[23].x = 	0.02	;	 M[51].Xp[23].y =14	; */
/*  M[51].Xp[24].x = 	0.029	;	 M[51].Xp[24].y =14.6	; */
/*  M[51].Xp[25].x = 	0.022	;	 M[51].Xp[25].y =15.2	; */
/*  M[51].Xp[26].x = 	0.024	;	 M[51].Xp[26].y =15.8	; */
/*  M[51].Xp[27].x = 	0.016	;	 M[51].Xp[27].y =16.4	; */
/*  M[51].Xp[28].x = 	0.017	;	 M[51].Xp[28].y =17	; */
/*  M[51].Xp[29].x = 	0.017	;	 M[51].Xp[29].y =17.6	; */
/*  M[51].Xp[30].x = 	0.015	;	 M[51].Xp[30].y =18.2	; */
/*  M[51].Xp[31].x = 	0.005	;	 M[51].Xp[31].y =18.8	; */
/*  M[51].Xp[32].x = 	0.007	;	 M[51].Xp[32].y =19.4	; */
/*  M[51].Xp[33].x = 	0.006	;	 M[51].Xp[33].y =20	; */
/*  M[51].Xp[34].x = 	0.008	;	 M[51].Xp[34].y =20.6	; */
/*  M[51].Xp[35].x = 	0.006	;	 M[51].Xp[35].y =21.2	; */
/*  M[51].Xp[36].x = 	0.008	;	 M[51].Xp[36].y =21.8	; */
/*  M[51].Xp[37].x = 	0.005	;	 M[51].Xp[37].y =22.4	; */
/*  M[51].Xp[38].x = 	0.007	;	 M[51].Xp[38].y =23	; */
/*  M[51].Xp[39].x = 	0.004	;	 M[51].Xp[39].y =23.6	; */
/*  M[51].Xp[40].x = 	0.002	;	 M[51].Xp[40].y =24.2	; */
/*  M[51].Xp[41].x = 	0.007	;	 M[51].Xp[41].y =24.8	; */
/*  M[51].Xp[42].x = 	0.002	;	 M[51].Xp[42].y =25.4	; */
/*  M[51].Xp[43].x = 	0.004	;	 M[51].Xp[43].y =26	; */
/*  M[51].pathlength(); */
/* //station 53 */
/*  x[52] = 67.0; */
/*  M[52].allocate(44); */
/*  M[52].settype(SPLINE2D_LINEAR); */
/*  M[52].Xp[0].x = 	0.068	;	 M[52].Xp[0].y =0.2	; */
/*  M[52].Xp[1].x = 	0.055	;	 M[52].Xp[1].y =0.8	; */
/*  M[52].Xp[2].x = 	0.058	;	 M[52].Xp[2].y =1.4	; */
/*  M[52].Xp[3].x = 	0.06	;	 M[52].Xp[3].y =2	; */
/*  M[52].Xp[4].x = 	0.076	;	 M[52].Xp[4].y =2.6	; */
/*  M[52].Xp[5].x = 	0.067	;	 M[52].Xp[5].y =3.2	; */
/*  M[52].Xp[6].x = 	0.073	;	 M[52].Xp[6].y =3.8	; */
/*  M[52].Xp[7].x = 	0.07	;	 M[52].Xp[7].y =4.4	; */
/*  M[52].Xp[8].x = 	0.06	;	 M[52].Xp[8].y =5	; */
/*  M[52].Xp[9].x = 	0.067	;	 M[52].Xp[9].y =5.6	; */
/*  M[52].Xp[10].x = 	0.042	;	 M[52].Xp[10].y =6.2	; */
/*  M[52].Xp[11].x = 	0.057	;	 M[52].Xp[11].y =6.8	; */
/*  M[52].Xp[12].x = 	0.052	;	 M[52].Xp[12].y =7.4	; */
/*  M[52].Xp[13].x = 	0.049	;	 M[52].Xp[13].y =8	; */
/*  M[52].Xp[14].x = 	0.047	;	 M[52].Xp[14].y =8.6	; */
/*  M[52].Xp[15].x = 	0.034	;	 M[52].Xp[15].y =9.2	; */
/*  M[52].Xp[16].x = 	0.038	;	 M[52].Xp[16].y =9.8	; */
/*  M[52].Xp[17].x = 	0.039	;	 M[52].Xp[17].y =10.4	; */
/*  M[52].Xp[18].x = 	0.031	;	 M[52].Xp[18].y =11	; */
/*  M[52].Xp[19].x = 	0.037	;	 M[52].Xp[19].y =11.6	; */
/*  M[52].Xp[20].x = 	0.022	;	 M[52].Xp[20].y =12.2	; */
/*  M[52].Xp[21].x = 	0.025	;	 M[52].Xp[21].y =12.8	; */
/*  M[52].Xp[22].x = 	0.018	;	 M[52].Xp[22].y =13.4	; */
/*  M[52].Xp[23].x = 	0.013	;	 M[52].Xp[23].y =14	; */
/*  M[52].Xp[24].x = 	0.015	;	 M[52].Xp[24].y =14.6	; */
/*  M[52].Xp[25].x = 	0.031	;	 M[52].Xp[25].y =15.2	; */
/*  M[52].Xp[26].x = 	0.017	;	 M[52].Xp[26].y =15.8	; */
/*  M[52].Xp[27].x = 	0.02	;	 M[52].Xp[27].y =16.4	; */
/*  M[52].Xp[28].x = 	0.011	;	 M[52].Xp[28].y =17	; */
/*  M[52].Xp[29].x = 	0.009	;	 M[52].Xp[29].y =17.6	; */
/*  M[52].Xp[30].x = 	0.012	;	 M[52].Xp[30].y =18.2	; */
/*  M[52].Xp[31].x = 	0.006	;	 M[52].Xp[31].y =18.8	; */
/*  M[52].Xp[32].x = 	0.008	;	 M[52].Xp[32].y =19.4	; */
/*  M[52].Xp[33].x = 	0.01	;	 M[52].Xp[33].y =20	; */
/*  M[52].Xp[34].x = 	0.002	;	 M[52].Xp[34].y =20.6	; */
/*  M[52].Xp[35].x = 	0.007	;	 M[52].Xp[35].y =21.2	; */
/*  M[52].Xp[36].x = 	0.003	;	 M[52].Xp[36].y =21.8	; */
/*  M[52].Xp[37].x = 	0.006	;	 M[52].Xp[37].y =22.4	; */
/*  M[52].Xp[38].x = 	0.007	;	 M[52].Xp[38].y =23	; */
/*  M[52].Xp[39].x = 	0.006	;	 M[52].Xp[39].y =23.6	; */
/*  M[52].Xp[40].x = 	0.001	;	 M[52].Xp[40].y =24.2	; */
/*  M[52].Xp[41].x = 	0.004	;	 M[52].Xp[41].y =24.8	; */
/*  M[52].Xp[42].x = 	0.004	;	 M[52].Xp[42].y =25.4	; */
/*  M[52].Xp[43].x = 	0.003	;	 M[52].Xp[43].y =26	; */
/*  M[52].pathlength(); */
/* //station 54 */
/*  x[53] = 67.8; */
/*  M[53].allocate(44); */
/*  M[53].settype(SPLINE2D_LINEAR); */
/*  M[53].Xp[0].x = 	0.079	;	 M[53].Xp[0].y =0.2	; */
/*  M[53].Xp[1].x = 	0.053	;	 M[53].Xp[1].y =0.8	; */
/*  M[53].Xp[2].x = 	0.052	;	 M[53].Xp[2].y =1.4	; */
/*  M[53].Xp[3].x = 	0.069	;	 M[53].Xp[3].y =2	; */
/*  M[53].Xp[4].x = 	0.074	;	 M[53].Xp[4].y =2.6	; */
/*  M[53].Xp[5].x = 	0.075	;	 M[53].Xp[5].y =3.2	; */
/*  M[53].Xp[6].x = 	0.066	;	 M[53].Xp[6].y =3.8	; */
/*  M[53].Xp[7].x = 	0.067	;	 M[53].Xp[7].y =4.4	; */
/*  M[53].Xp[8].x = 	0.075	;	 M[53].Xp[8].y =5	; */
/*  M[53].Xp[9].x = 	0.064	;	 M[53].Xp[9].y =5.6	; */
/*  M[53].Xp[10].x = 	0.049	;	 M[53].Xp[10].y =6.2	; */
/*  M[53].Xp[11].x = 	0.044	;	 M[53].Xp[11].y =6.8	; */
/*  M[53].Xp[12].x = 	0.045	;	 M[53].Xp[12].y =7.4	; */
/*  M[53].Xp[13].x = 	0.064	;	 M[53].Xp[13].y =8	; */
/*  M[53].Xp[14].x = 	0.037	;	 M[53].Xp[14].y =8.6	; */
/*  M[53].Xp[15].x = 	0.047	;	 M[53].Xp[15].y =9.2	; */
/*  M[53].Xp[16].x = 	0.042	;	 M[53].Xp[16].y =9.8	; */
/*  M[53].Xp[17].x = 	0.026	;	 M[53].Xp[17].y =10.4	; */
/*  M[53].Xp[18].x = 	0.036	;	 M[53].Xp[18].y =11	; */
/*  M[53].Xp[19].x = 	0.038	;	 M[53].Xp[19].y =11.6	; */
/*  M[53].Xp[20].x = 	0.034	;	 M[53].Xp[20].y =12.2	; */
/*  M[53].Xp[21].x = 	0.033	;	 M[53].Xp[21].y =12.8	; */
/*  M[53].Xp[22].x = 	0.019	;	 M[53].Xp[22].y =13.4	; */
/*  M[53].Xp[23].x = 	0.035	;	 M[53].Xp[23].y =14	; */
/*  M[53].Xp[24].x = 	0.022	;	 M[53].Xp[24].y =14.6	; */
/*  M[53].Xp[25].x = 	0.022	;	 M[53].Xp[25].y =15.2	; */
/*  M[53].Xp[26].x = 	0.031	;	 M[53].Xp[26].y =15.8	; */
/*  M[53].Xp[27].x = 	0.018	;	 M[53].Xp[27].y =16.4	; */
/*  M[53].Xp[28].x = 	0.026	;	 M[53].Xp[28].y =17	; */
/*  M[53].Xp[29].x = 	0.013	;	 M[53].Xp[29].y =17.6	; */
/*  M[53].Xp[30].x = 	0.019	;	 M[53].Xp[30].y =18.2	; */
/*  M[53].Xp[31].x = 	0.018	;	 M[53].Xp[31].y =18.8	; */
/*  M[53].Xp[32].x = 	0.004	;	 M[53].Xp[32].y =19.4	; */
/*  M[53].Xp[33].x = 	0.012	;	 M[53].Xp[33].y =20	; */
/*  M[53].Xp[34].x = 	0.012	;	 M[53].Xp[34].y =20.6	; */
/*  M[53].Xp[35].x = 	0.005	;	 M[53].Xp[35].y =21.2	; */
/*  M[53].Xp[36].x = 	0.005	;	 M[53].Xp[36].y =21.8	; */
/*  M[53].Xp[37].x = 	0.009	;	 M[53].Xp[37].y =22.4	; */
/*  M[53].Xp[38].x = 	0.008	;	 M[53].Xp[38].y =23	; */
/*  M[53].Xp[39].x = 	0.021	;	 M[53].Xp[39].y =23.6	; */
/*  M[53].Xp[40].x = 	0.013	;	 M[53].Xp[40].y =24.2	; */
/*  M[53].Xp[41].x = 	0.007	;	 M[53].Xp[41].y =24.8	; */
/*  M[53].Xp[42].x = 	0.007	;	 M[53].Xp[42].y =25.4	; */
/*  M[53].Xp[43].x = 	0.005	;	 M[53].Xp[43].y =26	; */

/*  M[53].pathlength(); */


/* Velocity of coflowing air is 20 m/s. */
/* Velocity of coflowing fuel (in this case "C2H4" ) is 50 m/s . */
 // mixing field ---  mixture fraction (C2H4)
  // Set the number of stations.
  x[0] = 20.0;
 M[0].allocate(44);                              
 M[0].settype(SPLINE2D_LINEAR);                                                                           
 M[ 0] . Xp[ 0] .y =    0.40;             M[ 0] . Xp[ 0] .x=  0.738;          
 M[ 0] . Xp[ 1] .y =    1.00;             M[ 0] . Xp[ 1] .x=  0.761;          
 M[ 0] . Xp[ 2] .y =    1.60;             M[ 0] . Xp[ 2] .x=  0.684;          
 M[ 0] . Xp[ 3] .y =    2.20;             M[ 0] . Xp[ 3] .x=  0.581;          
 M[ 0] . Xp[ 4] .y =    2.80;             M[ 0] . Xp[ 4] .x=  0.524;          
 M[ 0] . Xp[ 5] .y =    3.40;             M[ 0] . Xp[ 5] .x=  0.423;          
 M[ 0] . Xp[ 6] .y =    4.00;             M[ 0] . Xp[ 6] .x=  0.305;          
 M[ 0] . Xp[ 7] .y =    4.60;             M[ 0] . Xp[ 7] .x=  0.288;          
 M[ 0] . Xp[ 8] .y =    5.20;             M[ 0] . Xp[ 8] .x=  0.246;          
 M[ 0] . Xp[ 9] .y =    5.80;             M[ 0] . Xp[ 9] .x=  0.210;          
 M[ 0] . Xp[10] .y =    6.40;             M[ 0] . Xp[10] .x=  0.186;          
 M[ 0] . Xp[11] .y =    7.00;             M[ 0] . Xp[11] .x=  0.174;          
 M[ 0] . Xp[12] .y =    7.60;             M[ 0] . Xp[12] .x=  0.171;          
 M[ 0] . Xp[13] .y =    8.20;             M[ 0] . Xp[13] .x=  0.180;          
 M[ 0] . Xp[14] .y =    8.80;             M[ 0] . Xp[14] .x=  0.152;          
 M[ 0] . Xp[15] .y =    9.40;             M[ 0] . Xp[15] .x=  0.153;          
 M[ 0] . Xp[16] .y =   10.00;             M[ 0] . Xp[16] .x=  0.150;          
 M[ 0] . Xp[17] .y =   10.60;             M[ 0] . Xp[17] .x=  0.140;          
 M[ 0] . Xp[18] .y =   11.20;             M[ 0] . Xp[18] .x=  0.151;          
 M[ 0] . Xp[19] .y =   11.80;             M[ 0] . Xp[19] .x=  0.154;          
 M[ 0] . Xp[20] .y =   12.40;             M[ 0] . Xp[20] .x=  0.151;          
 M[ 0] . Xp[21] .y =   13.00;             M[ 0] . Xp[21] .x=  0.138;          
 M[ 0] . Xp[22] .y =   13.60;             M[ 0] . Xp[22] .x=  0.135;          
 M[ 0] . Xp[23] .y =   14.20;             M[ 0] . Xp[23] .x=  0.138;          
 M[ 0] . Xp[24] .y =   14.80;             M[ 0] . Xp[24] .x=  0.138;          
 M[ 0] . Xp[25] .y =   15.40;             M[ 0] . Xp[25] .x=  0.136;          
 M[ 0] . Xp[26] .y =   16.00;             M[ 0] . Xp[26] .x=  0.146;          
 M[ 0] . Xp[27] .y =   16.60;             M[ 0] . Xp[27] .x=  0.126;          
 M[ 0] . Xp[28] .y =   17.20;             M[ 0] . Xp[28] .x=  0.130;          
 M[ 0] . Xp[29] .y =   17.80;             M[ 0] . Xp[29] .x=  0.115;          
 M[ 0] . Xp[30] .y =   18.40;             M[ 0] . Xp[30] .x=  0.128;          
 M[ 0] . Xp[31] .y =   19.00;             M[ 0] . Xp[31] .x=  0.132;          
 M[ 0] . Xp[32] .y =   19.60;             M[ 0] . Xp[32] .x=  0.131;          
 M[ 0] . Xp[33] .y =   20.20;             M[ 0] . Xp[33] .x=  0.115;          
 M[ 0] . Xp[34] .y =   20.80;             M[ 0] . Xp[34] .x=  0.107;          
 M[ 0] . Xp[35] .y =   21.40;             M[ 0] . Xp[35] .x=  0.086;          
 M[ 0] . Xp[36] .y =   22.00;             M[ 0] . Xp[36] .x=  0.075;          
 M[ 0] . Xp[37] .y =   22.60;             M[ 0] . Xp[37] .x=  0.044;          
 M[ 0] . Xp[38] .y =   23.20;             M[ 0] . Xp[38] .x=  0.026;          
 M[ 0] . Xp[39] .y =   23.80;             M[ 0] . Xp[39] .x=  0.014;          
 M[ 0] . Xp[40] .y =   24.40;             M[ 0] . Xp[40] .x=  0.007;          
 M[ 0] . Xp[41] .y =   25.00;             M[ 0] . Xp[41] .x=  0.004;          
 M[ 0] . Xp[42] .y =   25.60;             M[ 0] . Xp[42] .x=  0.001;          
 M[ 0] . Xp[43] .y =   26.20;             M[ 0] . Xp[43] .x=  0.006;          
 M[0].pathlength();
//Station 2
 x[1] = 20.80;
 M[1].allocate(44);                              
 M[1].settype(SPLINE2D_LINEAR);                                                                             
                                                                      
 M[ 1] . Xp[ 0] .y =    0.40;             M[ 1] . Xp[ 0] .x=  0.823;          
 M[ 1] . Xp[ 1] .y =    1.00;             M[ 1] . Xp[ 1] .x=  0.713;          
 M[ 1] . Xp[ 2] .y =    1.60;             M[ 1] . Xp[ 2] .x=  0.653;          
 M[ 1] . Xp[ 3] .y =    2.20;             M[ 1] . Xp[ 3] .x=  0.606;          
 M[ 1] . Xp[ 4] .y =    2.80;             M[ 1] . Xp[ 4] .x=  0.498;          
 M[ 1] . Xp[ 5] .y =    3.40;             M[ 1] . Xp[ 5] .x=  0.397;          
 M[ 1] . Xp[ 6] .y =    4.00;             M[ 1] . Xp[ 6] .x=  0.371;          
 M[ 1] . Xp[ 7] .y =    4.60;             M[ 1] . Xp[ 7] .x=  0.298;          
 M[ 1] . Xp[ 8] .y =    5.20;             M[ 1] . Xp[ 8] .x=  0.279;          
 M[ 1] . Xp[ 9] .y =    5.80;             M[ 1] . Xp[ 9] .x=  0.194;          
 M[ 1] . Xp[10] .y =    6.40;             M[ 1] . Xp[10] .x=  0.197;          
 M[ 1] . Xp[11] .y =    7.00;             M[ 1] . Xp[11] .x=  0.167;          
 M[ 1] . Xp[12] .y =    7.60;             M[ 1] . Xp[12] .x=  0.160;          
 M[ 1] . Xp[13] .y =    8.20;             M[ 1] . Xp[13] .x=  0.155;          
 M[ 1] . Xp[14] .y =    8.80;             M[ 1] . Xp[14] .x=  0.149;          
 M[ 1] . Xp[15] .y =    9.40;             M[ 1] . Xp[15] .x=  0.150;          
 M[ 1] . Xp[16] .y =   10.00;             M[ 1] . Xp[16] .x=  0.145;          
 M[ 1] . Xp[17] .y =   10.60;             M[ 1] . Xp[17] .x=  0.160;          
 M[ 1] . Xp[18] .y =   11.20;             M[ 1] . Xp[18] .x=  0.162;          
 M[ 1] . Xp[19] .y =   11.80;             M[ 1] . Xp[19] .x=  0.156;          
 M[ 1] . Xp[20] .y =   12.40;             M[ 1] . Xp[20] .x=  0.143;          
 M[ 1] . Xp[21] .y =   13.00;             M[ 1] . Xp[21] .x=  0.140;          
 M[ 1] . Xp[22] .y =   13.60;             M[ 1] . Xp[22] .x=  0.146;          
 M[ 1] . Xp[23] .y =   14.20;             M[ 1] . Xp[23] .x=  0.141;          
 M[ 1] . Xp[24] .y =   14.80;             M[ 1] . Xp[24] .x=  0.137;          
 M[ 1] . Xp[25] .y =   15.40;             M[ 1] . Xp[25] .x=  0.133;          
 M[ 1] . Xp[26] .y =   16.00;             M[ 1] . Xp[26] .x=  0.135;          
 M[ 1] . Xp[27] .y =   16.60;             M[ 1] . Xp[27] .x=  0.141;          
 M[ 1] . Xp[28] .y =   17.20;             M[ 1] . Xp[28] .x=  0.130;          
 M[ 1] . Xp[29] .y =   17.80;             M[ 1] . Xp[29] .x=  0.125;          
 M[ 1] . Xp[30] .y =   18.40;             M[ 1] . Xp[30] .x=  0.136;          
 M[ 1] . Xp[31] .y =   19.00;             M[ 1] . Xp[31] .x=  0.118;          
 M[ 1] . Xp[32] .y =   19.60;             M[ 1] . Xp[32] .x=  0.121;          
 M[ 1] . Xp[33] .y =   20.20;             M[ 1] . Xp[33] .x=  0.115;          
 M[ 1] . Xp[34] .y =   20.80;             M[ 1] . Xp[34] .x=  0.105;          
 M[ 1] . Xp[35] .y =   21.40;             M[ 1] . Xp[35] .x=  0.090;          
 M[ 1] . Xp[36] .y =   22.00;             M[ 1] . Xp[36] .x=  0.066;          
 M[ 1] . Xp[37] .y =   22.60;             M[ 1] . Xp[37] .x=  0.052;          
 M[ 1] . Xp[38] .y =   23.20;             M[ 1] . Xp[38] .x=  0.024;          
 M[ 1] . Xp[39] .y =   23.80;             M[ 1] . Xp[39] .x=  0.019;          
 M[ 1] . Xp[40] .y =   24.40;             M[ 1] . Xp[40] .x=  0.003;          
 M[ 1] . Xp[41] .y =   25.00;             M[ 1] . Xp[41] .x=  0.001;          
 M[ 1] . Xp[42] .y =   25.60;             M[ 1] . Xp[42] .x=  0.001;          
 M[ 1] . Xp[43] .y =   26.20;             M[ 1] . Xp[43] .x=  0.001;          
   M[1].pathlength();   
//station 3
 x[2] = 21.50;
  M[2].allocate(44);                              
  M[2].settype(SPLINE2D_LINEAR);                                                                             
                                                                       
 M[ 2] . Xp[ 0] .y =    0.40;             M[ 2] . Xp[ 0] .x=  0.819;          
 M[ 2] . Xp[ 1] .y =    1.00;             M[ 2] . Xp[ 1] .x=  0.717;          
 M[ 2] . Xp[ 2] .y =    1.60;             M[ 2] . Xp[ 2] .x=  0.685;          
 M[ 2] . Xp[ 3] .y =    2.20;             M[ 2] . Xp[ 3] .x=  0.600;          
 M[ 2] . Xp[ 4] .y =    2.80;             M[ 2] . Xp[ 4] .x=  0.490;          
 M[ 2] . Xp[ 5] .y =    3.40;             M[ 2] . Xp[ 5] .x=  0.457;          
 M[ 2] . Xp[ 6] .y =    4.00;             M[ 2] . Xp[ 6] .x=  0.367;          
 M[ 2] . Xp[ 7] .y =    4.60;             M[ 2] . Xp[ 7] .x=  0.317;          
 M[ 2] . Xp[ 8] .y =    5.20;             M[ 2] . Xp[ 8] .x=  0.241;          
 M[ 2] . Xp[ 9] .y =    5.80;             M[ 2] . Xp[ 9] .x=  0.216;          
 M[ 2] . Xp[10] .y =    6.40;             M[ 2] . Xp[10] .x=  0.197;          
 M[ 2] . Xp[11] .y =    7.00;             M[ 2] . Xp[11] .x=  0.175;          
 M[ 2] . Xp[12] .y =    7.60;             M[ 2] . Xp[12] .x=  0.184;          
 M[ 2] . Xp[13] .y =    8.20;             M[ 2] . Xp[13] .x=  0.176;          
 M[ 2] . Xp[14] .y =    8.80;             M[ 2] . Xp[14] .x=  0.172;          
 M[ 2] . Xp[15] .y =    9.40;             M[ 2] . Xp[15] .x=  0.174;          
 M[ 2] . Xp[16] .y =   10.00;             M[ 2] . Xp[16] .x=  0.165;          
 M[ 2] . Xp[17] .y =   10.60;             M[ 2] . Xp[17] .x=  0.164;          
 M[ 2] . Xp[18] .y =   11.20;             M[ 2] . Xp[18] .x=  0.157;          
 M[ 2] . Xp[19] .y =   11.80;             M[ 2] . Xp[19] .x=  0.157;          
 M[ 2] . Xp[20] .y =   12.40;             M[ 2] . Xp[20] .x=  0.156;          
 M[ 2] . Xp[21] .y =   13.00;             M[ 2] . Xp[21] .x=  0.146;          
 M[ 2] . Xp[22] .y =   13.60;             M[ 2] . Xp[22] .x=  0.141;          
 M[ 2] . Xp[23] .y =   14.20;             M[ 2] . Xp[23] .x=  0.141;          
 M[ 2] . Xp[24] .y =   14.80;             M[ 2] . Xp[24] .x=  0.155;          
 M[ 2] . Xp[25] .y =   15.40;             M[ 2] . Xp[25] .x=  0.139;          
 M[ 2] . Xp[26] .y =   16.00;             M[ 2] . Xp[26] .x=  0.142;          
 M[ 2] . Xp[27] .y =   16.60;             M[ 2] . Xp[27] .x=  0.134;          
 M[ 2] . Xp[28] .y =   17.20;             M[ 2] . Xp[28] .x=  0.133;          
 M[ 2] . Xp[29] .y =   17.80;             M[ 2] . Xp[29] .x=  0.130;          
 M[ 2] . Xp[30] .y =   18.40;             M[ 2] . Xp[30] .x=  0.121;          
 M[ 2] . Xp[31] .y =   19.00;             M[ 2] . Xp[31] .x=  0.128;          
 M[ 2] . Xp[32] .y =   19.60;             M[ 2] . Xp[32] .x=  0.121;          
 M[ 2] . Xp[33] .y =   20.20;             M[ 2] . Xp[33] .x=  0.111;          
 M[ 2] . Xp[34] .y =   20.80;             M[ 2] . Xp[34] .x=  0.110;          
 M[ 2] . Xp[35] .y =   21.40;             M[ 2] . Xp[35] .x=  0.090;          
 M[ 2] . Xp[36] .y =   22.00;             M[ 2] . Xp[36] .x=  0.070;          
 M[ 2] . Xp[37] .y =   22.60;             M[ 2] . Xp[37] .x=  0.038;          
 M[ 2] . Xp[38] .y =   23.20;             M[ 2] . Xp[38] .x=  0.026;          
 M[ 2] . Xp[39] .y =   23.80;             M[ 2] . Xp[39] .x=  0.016;          
 M[ 2] . Xp[40] .y =   24.40;             M[ 2] . Xp[40] .x=  0.006;          
 M[ 2] . Xp[41] .y =   25.00;             M[ 2] . Xp[41] .x=  0.002;          
 M[ 2] . Xp[42] .y =   25.60;             M[ 2] . Xp[42] .x=  0.001;          
 M[ 2] . Xp[43] .y =   26.20;             M[ 2] . Xp[43] .x=  0.001;          
   M[2].pathlength();                                                                             
//station 4
 x[3] = 22.30;
  M[3].allocate(44);                              
  M[3].settype(SPLINE2D_LINEAR);                                                                           
 M[ 3] . Xp[ 0] .y =    0.40;             M[ 3] . Xp[ 0] .x=  0.845;          
 M[ 3] . Xp[ 1] .y =    1.00;             M[ 3] . Xp[ 1] .x=  0.781;          
 M[ 3] . Xp[ 2] .y =    1.60;             M[ 3] . Xp[ 2] .x=  0.697;          
 M[ 3] . Xp[ 3] .y =    2.20;             M[ 3] . Xp[ 3] .x=  0.611;          
 M[ 3] . Xp[ 4] .y =    2.80;             M[ 3] . Xp[ 4] .x=  0.539;          
 M[ 3] . Xp[ 5] .y =    3.40;             M[ 3] . Xp[ 5] .x=  0.479;          
 M[ 3] . Xp[ 6] .y =    4.00;             M[ 3] . Xp[ 6] .x=  0.381;          
 M[ 3] . Xp[ 7] .y =    4.60;             M[ 3] . Xp[ 7] .x=  0.307;          
 M[ 3] . Xp[ 8] .y =    5.20;             M[ 3] . Xp[ 8] .x=  0.247;          
 M[ 3] . Xp[ 9] .y =    5.80;             M[ 3] . Xp[ 9] .x=  0.225;          
 M[ 3] . Xp[10] .y =    6.40;             M[ 3] . Xp[10] .x=  0.221;          
 M[ 3] . Xp[11] .y =    7.00;             M[ 3] . Xp[11] .x=  0.185;          
 M[ 3] . Xp[12] .y =    7.60;             M[ 3] . Xp[12] .x=  0.177;          
 M[ 3] . Xp[13] .y =    8.20;             M[ 3] . Xp[13] .x=  0.186;          
 M[ 3] . Xp[14] .y =    8.80;             M[ 3] . Xp[14] .x=  0.172;          
 M[ 3] . Xp[15] .y =    9.40;             M[ 3] . Xp[15] .x=  0.173;          
 M[ 3] . Xp[16] .y =   10.00;             M[ 3] . Xp[16] .x=  0.166;          
 M[ 3] . Xp[17] .y =   10.60;             M[ 3] . Xp[17] .x=  0.171;          
 M[ 3] . Xp[18] .y =   11.20;             M[ 3] . Xp[18] .x=  0.163;          
 M[ 3] . Xp[19] .y =   11.80;             M[ 3] . Xp[19] .x=  0.175;          
 M[ 3] . Xp[20] .y =   12.40;             M[ 3] . Xp[20] .x=  0.148;          
 M[ 3] . Xp[21] .y =   13.00;             M[ 3] . Xp[21] .x=  0.155;          
 M[ 3] . Xp[22] .y =   13.60;             M[ 3] . Xp[22] .x=  0.166;          
 M[ 3] . Xp[23] .y =   14.20;             M[ 3] . Xp[23] .x=  0.149;          
 M[ 3] . Xp[24] .y =   14.80;             M[ 3] . Xp[24] .x=  0.145;          
 M[ 3] . Xp[25] .y =   15.40;             M[ 3] . Xp[25] .x=  0.142;          
 M[ 3] . Xp[26] .y =   16.00;             M[ 3] . Xp[26] .x=  0.146;          
 M[ 3] . Xp[27] .y =   16.60;             M[ 3] . Xp[27] .x=  0.139;          
 M[ 3] . Xp[28] .y =   17.20;             M[ 3] . Xp[28] .x=  0.141;          
 M[ 3] . Xp[29] .y =   17.80;             M[ 3] . Xp[29] .x=  0.129;          
 M[ 3] . Xp[30] .y =   18.40;             M[ 3] . Xp[30] .x=  0.135;          
 M[ 3] . Xp[31] .y =   19.00;             M[ 3] . Xp[31] .x=  0.132;          
 M[ 3] . Xp[32] .y =   19.60;             M[ 3] . Xp[32] .x=  0.127;          
 M[ 3] . Xp[33] .y =   20.20;             M[ 3] . Xp[33] .x=  0.117;          
 M[ 3] . Xp[34] .y =   20.80;             M[ 3] . Xp[34] .x=  0.099;          
 M[ 3] . Xp[35] .y =   21.40;             M[ 3] . Xp[35] .x=  0.091;          
 M[ 3] . Xp[36] .y =   22.00;             M[ 3] . Xp[36] .x=  0.062;          
 M[ 3] . Xp[37] .y =   22.60;             M[ 3] . Xp[37] .x=  0.037;          
 M[ 3] . Xp[38] .y =   23.20;             M[ 3] . Xp[38] .x=  0.030;          
 M[ 3] . Xp[39] .y =   23.80;             M[ 3] . Xp[39] .x=  0.017;          
 M[ 3] . Xp[40] .y =   24.40;             M[ 3] . Xp[40] .x=  0.006;          
 M[ 3] . Xp[41] .y =   25.00;             M[ 3] . Xp[41] .x=  0.002;          
 M[ 3] . Xp[42] .y =   25.60;             M[ 3] . Xp[42] .x=  0.002;          
 M[ 3] . Xp[43] .y =   26.20;             M[ 3] . Xp[43] .x=  0.001;          
   M[3].pathlength();
 //station 5
 x[4] = 23.00;
  M[4].allocate(44);                              
  M[4].settype(SPLINE2D_LINEAR);                                                                               
                                                                          
 M[ 4] . Xp[ 0] .y =    0.40;             M[ 4] . Xp[ 0] .x=  0.830;          
 M[ 4] . Xp[ 1] .y =    1.00;             M[ 4] . Xp[ 1] .x=  0.764;          
 M[ 4] . Xp[ 2] .y =    1.60;             M[ 4] . Xp[ 2] .x=  0.714;          
 M[ 4] . Xp[ 3] .y =    2.20;             M[ 4] . Xp[ 3] .x=  0.641;          
 M[ 4] . Xp[ 4] .y =    2.80;             M[ 4] . Xp[ 4] .x=  0.568;          
 M[ 4] . Xp[ 5] .y =    3.40;             M[ 4] . Xp[ 5] .x=  0.500;          
 M[ 4] . Xp[ 6] .y =    4.00;             M[ 4] . Xp[ 6] .x=  0.408;          
 M[ 4] . Xp[ 7] .y =    4.60;             M[ 4] . Xp[ 7] .x=  0.314;          
 M[ 4] . Xp[ 8] .y =    5.20;             M[ 4] . Xp[ 8] .x=  0.281;          
 M[ 4] . Xp[ 9] .y =    5.80;             M[ 4] . Xp[ 9] .x=  0.228;          
 M[ 4] . Xp[10] .y =    6.40;             M[ 4] . Xp[10] .x=  0.221;          
 M[ 4] . Xp[11] .y =    7.00;             M[ 4] . Xp[11] .x=  0.199;          
 M[ 4] . Xp[12] .y =    7.60;             M[ 4] . Xp[12] .x=  0.183;          
 M[ 4] . Xp[13] .y =    8.20;             M[ 4] . Xp[13] .x=  0.181;          
 M[ 4] . Xp[14] .y =    8.80;             M[ 4] . Xp[14] .x=  0.181;          
 M[ 4] . Xp[15] .y =    9.40;             M[ 4] . Xp[15] .x=  0.169;          
 M[ 4] . Xp[16] .y =   10.00;             M[ 4] . Xp[16] .x=  0.154;          
 M[ 4] . Xp[17] .y =   10.60;             M[ 4] . Xp[17] .x=  0.166;          
 M[ 4] . Xp[18] .y =   11.20;             M[ 4] . Xp[18] .x=  0.181;          
 M[ 4] . Xp[19] .y =   11.80;             M[ 4] . Xp[19] .x=  0.164;          
 M[ 4] . Xp[20] .y =   12.40;             M[ 4] . Xp[20] .x=  0.155;          
 M[ 4] . Xp[21] .y =   13.00;             M[ 4] . Xp[21] .x=  0.152;          
 M[ 4] . Xp[22] .y =   13.60;             M[ 4] . Xp[22] .x=  0.150;          
 M[ 4] . Xp[23] .y =   14.20;             M[ 4] . Xp[23] .x=  0.143;          
 M[ 4] . Xp[24] .y =   14.80;             M[ 4] . Xp[24] .x=  0.146;          
 M[ 4] . Xp[25] .y =   15.40;             M[ 4] . Xp[25] .x=  0.144;          
 M[ 4] . Xp[26] .y =   16.00;             M[ 4] . Xp[26] .x=  0.136;          
 M[ 4] . Xp[27] .y =   16.60;             M[ 4] . Xp[27] .x=  0.140;          
 M[ 4] . Xp[28] .y =   17.20;             M[ 4] . Xp[28] .x=  0.147;          
 M[ 4] . Xp[29] .y =   17.80;             M[ 4] . Xp[29] .x=  0.140;          
 M[ 4] . Xp[30] .y =   18.40;             M[ 4] . Xp[30] .x=  0.134;          
 M[ 4] . Xp[31] .y =   19.00;             M[ 4] . Xp[31] .x=  0.131;          
 M[ 4] . Xp[32] .y =   19.60;             M[ 4] . Xp[32] .x=  0.131;          
 M[ 4] . Xp[33] .y =   20.20;             M[ 4] . Xp[33] .x=  0.116;          
 M[ 4] . Xp[34] .y =   20.80;             M[ 4] . Xp[34] .x=  0.101;          
 M[ 4] . Xp[35] .y =   21.40;             M[ 4] . Xp[35] .x=  0.083;          
 M[ 4] . Xp[36] .y =   22.00;             M[ 4] . Xp[36] .x=  0.061;          
 M[ 4] . Xp[37] .y =   22.60;             M[ 4] . Xp[37] .x=  0.046;          
 M[ 4] . Xp[38] .y =   23.20;             M[ 4] . Xp[38] .x=  0.024;          
 M[ 4] . Xp[39] .y =   23.80;             M[ 4] . Xp[39] .x=  0.019;          
 M[ 4] . Xp[40] .y =   24.40;             M[ 4] . Xp[40] .x=  0.007;          
 M[ 4] . Xp[41] .y =   25.00;             M[ 4] . Xp[41] .x=  0.005;          
 M[ 4] . Xp[42] .y =   25.60;             M[ 4] . Xp[42] .x=  0.001;          
 M[ 4] . Xp[43] .y =   26.20;             M[ 4] . Xp[43] .x=  0.000;          
  //station 6
   M[4].pathlength();
 x[5] = 23.80;
  M[5].allocate(44);                              
  M[5].settype(SPLINE2D_LINEAR);                                                                                                                                           
 M[ 5] . Xp[ 0] .y =    0.40;             M[ 5] . Xp[ 0] .x=  0.844;          
 M[ 5] . Xp[ 1] .y =    1.00;             M[ 5] . Xp[ 1] .x=  0.792;          
 M[ 5] . Xp[ 2] .y =    1.60;             M[ 5] . Xp[ 2] .x=  0.735;          
 M[ 5] . Xp[ 3] .y =    2.20;             M[ 5] . Xp[ 3] .x=  0.628;          
 M[ 5] . Xp[ 4] .y =    2.80;             M[ 5] . Xp[ 4] .x=  0.574;          
 M[ 5] . Xp[ 5] .y =    3.40;             M[ 5] . Xp[ 5] .x=  0.505;          
 M[ 5] . Xp[ 6] .y =    4.00;             M[ 5] . Xp[ 6] .x=  0.399;          
 M[ 5] . Xp[ 7] .y =    4.60;             M[ 5] . Xp[ 7] .x=  0.329;          
 M[ 5] . Xp[ 8] .y =    5.20;             M[ 5] . Xp[ 8] .x=  0.316;          
 M[ 5] . Xp[ 9] .y =    5.80;             M[ 5] . Xp[ 9] .x=  0.260;          
 M[ 5] . Xp[10] .y =    6.40;             M[ 5] . Xp[10] .x=  0.213;          
 M[ 5] . Xp[11] .y =    7.00;             M[ 5] . Xp[11] .x=  0.201;          
 M[ 5] . Xp[12] .y =    7.60;             M[ 5] . Xp[12] .x=  0.208;          
 M[ 5] . Xp[13] .y =    8.20;             M[ 5] . Xp[13] .x=  0.190;          
 M[ 5] . Xp[14] .y =    8.80;             M[ 5] . Xp[14] .x=  0.184;          
 M[ 5] . Xp[15] .y =    9.40;             M[ 5] . Xp[15] .x=  0.179;          
 M[ 5] . Xp[16] .y =   10.00;             M[ 5] . Xp[16] .x=  0.172;          
 M[ 5] . Xp[17] .y =   10.60;             M[ 5] . Xp[17] .x=  0.174;          
 M[ 5] . Xp[18] .y =   11.20;             M[ 5] . Xp[18] .x=  0.180;          
 M[ 5] . Xp[19] .y =   11.80;             M[ 5] . Xp[19] .x=  0.167;          
 M[ 5] . Xp[20] .y =   12.40;             M[ 5] . Xp[20] .x=  0.159;          
 M[ 5] . Xp[21] .y =   13.00;             M[ 5] . Xp[21] .x=  0.155;          
 M[ 5] . Xp[22] .y =   13.60;             M[ 5] . Xp[22] .x=  0.156;          
 M[ 5] . Xp[23] .y =   14.20;             M[ 5] . Xp[23] .x=  0.149;          
 M[ 5] . Xp[24] .y =   14.80;             M[ 5] . Xp[24] .x=  0.146;          
 M[ 5] . Xp[25] .y =   15.40;             M[ 5] . Xp[25] .x=  0.143;          
 M[ 5] . Xp[26] .y =   16.00;             M[ 5] . Xp[26] .x=  0.151;          
 M[ 5] . Xp[27] .y =   16.60;             M[ 5] . Xp[27] .x=  0.148;          
 M[ 5] . Xp[28] .y =   17.20;             M[ 5] . Xp[28] .x=  0.141;          
 M[ 5] . Xp[29] .y =   17.80;             M[ 5] . Xp[29] .x=  0.140;          
 M[ 5] . Xp[30] .y =   18.40;             M[ 5] . Xp[30] .x=  0.140;          
 M[ 5] . Xp[31] .y =   19.00;             M[ 5] . Xp[31] .x=  0.136;          
 M[ 5] . Xp[32] .y =   19.60;             M[ 5] . Xp[32] .x=  0.132;          
 M[ 5] . Xp[33] .y =   20.20;             M[ 5] . Xp[33] .x=  0.123;          
 M[ 5] . Xp[34] .y =   20.80;             M[ 5] . Xp[34] .x=  0.098;          
 M[ 5] . Xp[35] .y =   21.40;             M[ 5] . Xp[35] .x=  0.084;          
 M[ 5] . Xp[36] .y =   22.00;             M[ 5] . Xp[36] .x=  0.062;          
 M[ 5] . Xp[37] .y =   22.60;             M[ 5] . Xp[37] .x=  0.051;          
 M[ 5] . Xp[38] .y =   23.20;             M[ 5] . Xp[38] .x=  0.029;          
 M[ 5] . Xp[39] .y =   23.80;             M[ 5] . Xp[39] .x=  0.021;          
 M[ 5] . Xp[40] .y =   24.40;             M[ 5] . Xp[40] .x=  0.010;          
 M[ 5] . Xp[41] .y =   25.00;             M[ 5] . Xp[41] .x=  0.004;          
 M[ 5] . Xp[42] .y =   25.60;             M[ 5] . Xp[42] .x=  0.003;          
 M[ 5] . Xp[43] .y =   26.20;             M[ 5] . Xp[43] .x=  0.002;          
   M[5].pathlength();
//station 7
   x[6] = 24.50;
 M[6].allocate(44);                              
 M[6].settype(SPLINE2D_LINEAR);                                                                                                                                                         
 M[ 6] . Xp[ 0] .y =    0.40;             M[ 6] . Xp[ 0] .x=  0.818;          
 M[ 6] . Xp[ 1] .y =    1.00;             M[ 6] . Xp[ 1] .x=  0.763;          
 M[ 6] . Xp[ 2] .y =    1.60;             M[ 6] . Xp[ 2] .x=  0.708;          
 M[ 6] . Xp[ 3] .y =    2.20;             M[ 6] . Xp[ 3] .x=  0.635;          
 M[ 6] . Xp[ 4] .y =    2.80;             M[ 6] . Xp[ 4] .x=  0.582;          
 M[ 6] . Xp[ 5] .y =    3.40;             M[ 6] . Xp[ 5] .x=  0.463;          
 M[ 6] . Xp[ 6] .y =    4.00;             M[ 6] . Xp[ 6] .x=  0.411;          
 M[ 6] . Xp[ 7] .y =    4.60;             M[ 6] . Xp[ 7] .x=  0.342;          
 M[ 6] . Xp[ 8] .y =    5.20;             M[ 6] . Xp[ 8] .x=  0.291;          
 M[ 6] . Xp[ 9] .y =    5.80;             M[ 6] . Xp[ 9] .x=  0.253;          
 M[ 6] . Xp[10] .y =    6.40;             M[ 6] . Xp[10] .x=  0.217;          
 M[ 6] . Xp[11] .y =    7.00;             M[ 6] . Xp[11] .x=  0.204;          
 M[ 6] . Xp[12] .y =    7.60;             M[ 6] . Xp[12] .x=  0.199;          
 M[ 6] . Xp[13] .y =    8.20;             M[ 6] . Xp[13] .x=  0.186;          
 M[ 6] . Xp[14] .y =    8.80;             M[ 6] . Xp[14] .x=  0.179;          
 M[ 6] . Xp[15] .y =    9.40;             M[ 6] . Xp[15] .x=  0.185;          
 M[ 6] . Xp[16] .y =   10.00;             M[ 6] . Xp[16] .x=  0.180;          
 M[ 6] . Xp[17] .y =   10.60;             M[ 6] . Xp[17] .x=  0.178;          
 M[ 6] . Xp[18] .y =   11.20;             M[ 6] . Xp[18] .x=  0.172;          
 M[ 6] . Xp[19] .y =   11.80;             M[ 6] . Xp[19] .x=  0.163;          
 M[ 6] . Xp[20] .y =   12.40;             M[ 6] . Xp[20] .x=  0.161;          
 M[ 6] . Xp[21] .y =   13.00;             M[ 6] . Xp[21] .x=  0.148;          
 M[ 6] . Xp[22] .y =   13.60;             M[ 6] . Xp[22] .x=  0.149;          
 M[ 6] . Xp[23] .y =   14.20;             M[ 6] . Xp[23] .x=  0.140;          
 M[ 6] . Xp[24] .y =   14.80;             M[ 6] . Xp[24] .x=  0.136;          
 M[ 6] . Xp[25] .y =   15.40;             M[ 6] . Xp[25] .x=  0.151;          
 M[ 6] . Xp[26] .y =   16.00;             M[ 6] . Xp[26] .x=  0.146;          
 M[ 6] . Xp[27] .y =   16.60;             M[ 6] . Xp[27] .x=  0.144;          
 M[ 6] . Xp[28] .y =   17.20;             M[ 6] . Xp[28] .x=  0.140;          
 M[ 6] . Xp[29] .y =   17.80;             M[ 6] . Xp[29] .x=  0.142;          
 M[ 6] . Xp[30] .y =   18.40;             M[ 6] . Xp[30] .x=  0.128;          
 M[ 6] . Xp[31] .y =   19.00;             M[ 6] . Xp[31] .x=  0.133;          
 M[ 6] . Xp[32] .y =   19.60;             M[ 6] . Xp[32] .x=  0.123;          
 M[ 6] . Xp[33] .y =   20.20;             M[ 6] . Xp[33] .x=  0.107;          
 M[ 6] . Xp[34] .y =   20.80;             M[ 6] . Xp[34] .x=  0.105;          
 M[ 6] . Xp[35] .y =   21.40;             M[ 6] . Xp[35] .x=  0.081;          
 M[ 6] . Xp[36] .y =   22.00;             M[ 6] . Xp[36] .x=  0.062;          
 M[ 6] . Xp[37] .y =   22.60;             M[ 6] . Xp[37] .x=  0.044;          
 M[ 6] . Xp[38] .y =   23.20;             M[ 6] . Xp[38] .x=  0.028;          
 M[ 6] . Xp[39] .y =   23.80;             M[ 6] . Xp[39] .x=  0.015;          
 M[ 6] . Xp[40] .y =   24.40;             M[ 6] . Xp[40] .x=  0.008;          
 M[ 6] . Xp[41] .y =   25.00;             M[ 6] . Xp[41] .x=  0.004;          
 M[ 6] . Xp[42] .y =   25.60;             M[ 6] . Xp[42] .x=  0.002;          
 M[ 6] . Xp[43] .y =   26.20;             M[ 6] . Xp[43] .x=  0.002;          
   M[6].pathlength();
//station 8
 x[7] = 25.30;
  M[7].allocate(44);                              
  M[7].settype(SPLINE2D_LINEAR);                                                               
                                                                       
 M[ 7] . Xp[ 0] .y =    0.40;             M[ 7] . Xp[ 0] .x=  0.843;          
 M[ 7] . Xp[ 1] .y =    1.00;             M[ 7] . Xp[ 1] .x=  0.793;          
 M[ 7] . Xp[ 2] .y =    1.60;             M[ 7] . Xp[ 2] .x=  0.725;          
 M[ 7] . Xp[ 3] .y =    2.20;             M[ 7] . Xp[ 3] .x=  0.647;          
 M[ 7] . Xp[ 4] .y =    2.80;             M[ 7] . Xp[ 4] .x=  0.557;          
 M[ 7] . Xp[ 5] .y =    3.40;             M[ 7] . Xp[ 5] .x=  0.483;          
 M[ 7] . Xp[ 6] .y =    4.00;             M[ 7] . Xp[ 6] .x=  0.397;          
 M[ 7] . Xp[ 7] .y =    4.60;             M[ 7] . Xp[ 7] .x=  0.347;          
 M[ 7] . Xp[ 8] .y =    5.20;             M[ 7] . Xp[ 8] .x=  0.294;          
 M[ 7] . Xp[ 9] .y =    5.80;             M[ 7] . Xp[ 9] .x=  0.273;          
 M[ 7] . Xp[10] .y =    6.40;             M[ 7] . Xp[10] .x=  0.226;          
 M[ 7] . Xp[11] .y =    7.00;             M[ 7] . Xp[11] .x=  0.211;          
 M[ 7] . Xp[12] .y =    7.60;             M[ 7] . Xp[12] .x=  0.198;          
 M[ 7] . Xp[13] .y =    8.20;             M[ 7] . Xp[13] .x=  0.190;          
 M[ 7] . Xp[14] .y =    8.80;             M[ 7] . Xp[14] .x=  0.188;          
 M[ 7] . Xp[15] .y =    9.40;             M[ 7] . Xp[15] .x=  0.192;          
 M[ 7] . Xp[16] .y =   10.00;             M[ 7] . Xp[16] .x=  0.185;          
 M[ 7] . Xp[17] .y =   10.60;             M[ 7] . Xp[17] .x=  0.178;          
 M[ 7] . Xp[18] .y =   11.20;             M[ 7] . Xp[18] .x=  0.167;          
 M[ 7] . Xp[19] .y =   11.80;             M[ 7] . Xp[19] .x=  0.174;          
 M[ 7] . Xp[20] .y =   12.40;             M[ 7] . Xp[20] .x=  0.167;          
 M[ 7] . Xp[21] .y =   13.00;             M[ 7] . Xp[21] .x=  0.146;          
 M[ 7] . Xp[22] .y =   13.60;             M[ 7] . Xp[22] .x=  0.155;          
 M[ 7] . Xp[23] .y =   14.20;             M[ 7] . Xp[23] .x=  0.134;          
 M[ 7] . Xp[24] .y =   14.80;             M[ 7] . Xp[24] .x=  0.146;          
 M[ 7] . Xp[25] .y =   15.40;             M[ 7] . Xp[25] .x=  0.152;          
 M[ 7] . Xp[26] .y =   16.00;             M[ 7] . Xp[26] .x=  0.152;          
 M[ 7] . Xp[27] .y =   16.60;             M[ 7] . Xp[27] .x=  0.146;          
 M[ 7] . Xp[28] .y =   17.20;             M[ 7] . Xp[28] .x=  0.136;          
 M[ 7] . Xp[29] .y =   17.80;             M[ 7] . Xp[29] .x=  0.140;          
 M[ 7] . Xp[30] .y =   18.40;             M[ 7] . Xp[30] .x=  0.125;          
 M[ 7] . Xp[31] .y =   19.00;             M[ 7] . Xp[31] .x=  0.128;          
 M[ 7] . Xp[32] .y =   19.60;             M[ 7] . Xp[32] .x=  0.130;          
 M[ 7] . Xp[33] .y =   20.20;             M[ 7] . Xp[33] .x=  0.108;          
 M[ 7] . Xp[34] .y =   20.80;             M[ 7] . Xp[34] .x=  0.089;          
 M[ 7] . Xp[35] .y =   21.40;             M[ 7] . Xp[35] .x=  0.078;          
 M[ 7] . Xp[36] .y =   22.00;             M[ 7] . Xp[36] .x=  0.064;          
 M[ 7] . Xp[37] .y =   22.60;             M[ 7] . Xp[37] .x=  0.044;          
 M[ 7] . Xp[38] .y =   23.20;             M[ 7] . Xp[38] .x=  0.032;          
 M[ 7] . Xp[39] .y =   23.80;             M[ 7] . Xp[39] .x=  0.014;          
 M[ 7] . Xp[40] .y =   24.40;             M[ 7] . Xp[40] .x=  0.006;          
 M[ 7] . Xp[41] .y =   25.00;             M[ 7] . Xp[41] .x=  0.003;          
 M[ 7] . Xp[42] .y =   25.60;             M[ 7] . Xp[42] .x=  0.001;          
 M[ 7] . Xp[43] .y =   26.20;             M[ 7] . Xp[43] .x=  0.001;          
   M[7].pathlength();
//station 9
 x[8] = 26.00;
  M[8].allocate(44);                              
  M[8].settype(SPLINE2D_LINEAR); 
                                                                                                                                                  
 M[ 8] . Xp[ 0] .y =    0.40;             M[ 8] . Xp[ 0] .x=  0.803;          
 M[ 8] . Xp[ 1] .y =    1.00;             M[ 8] . Xp[ 1] .x=  0.733;          
 M[ 8] . Xp[ 2] .y =    1.60;             M[ 8] . Xp[ 2] .x=  0.715;          
 M[ 8] . Xp[ 3] .y =    2.20;             M[ 8] . Xp[ 3] .x=  0.643;          
 M[ 8] . Xp[ 4] .y =    2.80;             M[ 8] . Xp[ 4] .x=  0.564;          
 M[ 8] . Xp[ 5] .y =    3.40;             M[ 8] . Xp[ 5] .x=  0.473;          
 M[ 8] . Xp[ 6] .y =    4.00;             M[ 8] . Xp[ 6] .x=  0.389;          
 M[ 8] . Xp[ 7] .y =    4.60;             M[ 8] . Xp[ 7] .x=  0.337;          
 M[ 8] . Xp[ 8] .y =    5.20;             M[ 8] . Xp[ 8] .x=  0.299;          
 M[ 8] . Xp[ 9] .y =    5.80;             M[ 8] . Xp[ 9] .x=  0.247;          
 M[ 8] . Xp[10] .y =    6.40;             M[ 8] . Xp[10] .x=  0.224;          
 M[ 8] . Xp[11] .y =    7.00;             M[ 8] . Xp[11] .x=  0.212;          
 M[ 8] . Xp[12] .y =    7.60;             M[ 8] . Xp[12] .x=  0.196;          
 M[ 8] . Xp[13] .y =    8.20;             M[ 8] . Xp[13] .x=  0.186;          
 M[ 8] . Xp[14] .y =    8.80;             M[ 8] . Xp[14] .x=  0.188;          
 M[ 8] . Xp[15] .y =    9.40;             M[ 8] . Xp[15] .x=  0.170;          
 M[ 8] . Xp[16] .y =   10.00;             M[ 8] . Xp[16] .x=  0.169;          
 M[ 8] . Xp[17] .y =   10.60;             M[ 8] . Xp[17] .x=  0.181;          
 M[ 8] . Xp[18] .y =   11.20;             M[ 8] . Xp[18] .x=  0.172;          
 M[ 8] . Xp[19] .y =   11.80;             M[ 8] . Xp[19] .x=  0.156;          
 M[ 8] . Xp[20] .y =   12.40;             M[ 8] . Xp[20] .x=  0.153;          
 M[ 8] . Xp[21] .y =   13.00;             M[ 8] . Xp[21] .x=  0.146;          
 M[ 8] . Xp[22] .y =   13.60;             M[ 8] . Xp[22] .x=  0.146;          
 M[ 8] . Xp[23] .y =   14.20;             M[ 8] . Xp[23] .x=  0.142;          
 M[ 8] . Xp[24] .y =   14.80;             M[ 8] . Xp[24] .x=  0.137;          
 M[ 8] . Xp[25] .y =   15.40;             M[ 8] . Xp[25] .x=  0.141;          
 M[ 8] . Xp[26] .y =   16.00;             M[ 8] . Xp[26] .x=  0.139;          
 M[ 8] . Xp[27] .y =   16.60;             M[ 8] . Xp[27] .x=  0.137;          
 M[ 8] . Xp[28] .y =   17.20;             M[ 8] . Xp[28] .x=  0.139;          
 M[ 8] . Xp[29] .y =   17.80;             M[ 8] . Xp[29] .x=  0.142;          
 M[ 8] . Xp[30] .y =   18.40;             M[ 8] . Xp[30] .x=  0.120;          
 M[ 8] . Xp[31] .y =   19.00;             M[ 8] . Xp[31] .x=  0.117;          
 M[ 8] . Xp[32] .y =   19.60;             M[ 8] . Xp[32] .x=  0.116;          
 M[ 8] . Xp[33] .y =   20.20;             M[ 8] . Xp[33] .x=  0.102;          
 M[ 8] . Xp[34] .y =   20.80;             M[ 8] . Xp[34] .x=  0.083;          
 M[ 8] . Xp[35] .y =   21.40;             M[ 8] . Xp[35] .x=  0.072;          
 M[ 8] . Xp[36] .y =   22.00;             M[ 8] . Xp[36] .x=  0.063;          
 M[ 8] . Xp[37] .y =   22.60;             M[ 8] . Xp[37] .x=  0.042;          
 M[ 8] . Xp[38] .y =   23.20;             M[ 8] . Xp[38] .x=  0.029;          
 M[ 8] . Xp[39] .y =   23.80;             M[ 8] . Xp[39] .x=  0.016;          
 M[ 8] . Xp[40] .y =   24.40;             M[ 8] . Xp[40] .x=  0.005;          
 M[ 8] . Xp[41] .y =   25.00;             M[ 8] . Xp[41] .x=  0.004;          
 M[ 8] . Xp[42] .y =   25.60;             M[ 8] . Xp[42] .x=  0.001;          
 M[ 8] . Xp[43] .y =   26.20;             M[ 8] . Xp[43] .x=  0.002;          
   M[8].pathlength();
//station 10
 x[9] = 26.800;
  M[9].allocate(44);                              
  M[9].settype(SPLINE2D_LINEAR);                                                                                                                                                     
 M[ 9] . Xp[ 0] .y =    0.40;             M[ 9] . Xp[ 0] .x=  0.802;          
 M[ 9] . Xp[ 1] .y =    1.00;             M[ 9] . Xp[ 1] .x=  0.775;          
 M[ 9] . Xp[ 2] .y =    1.60;             M[ 9] . Xp[ 2] .x=  0.721;          
 M[ 9] . Xp[ 3] .y =    2.20;             M[ 9] . Xp[ 3] .x=  0.635;          
 M[ 9] . Xp[ 4] .y =    2.80;             M[ 9] . Xp[ 4] .x=  0.528;          
 M[ 9] . Xp[ 5] .y =    3.40;             M[ 9] . Xp[ 5] .x=  0.453;          
 M[ 9] . Xp[ 6] .y =    4.00;             M[ 9] . Xp[ 6] .x=  0.389;          
 M[ 9] . Xp[ 7] .y =    4.60;             M[ 9] . Xp[ 7] .x=  0.341;          
 M[ 9] . Xp[ 8] .y =    5.20;             M[ 9] . Xp[ 8] .x=  0.302;          
 M[ 9] . Xp[ 9] .y =    5.80;             M[ 9] . Xp[ 9] .x=  0.279;          
 M[ 9] . Xp[10] .y =    6.40;             M[ 9] . Xp[10] .x=  0.230;          
 M[ 9] . Xp[11] .y =    7.00;             M[ 9] . Xp[11] .x=  0.222;          
 M[ 9] . Xp[12] .y =    7.60;             M[ 9] . Xp[12] .x=  0.201;          
 M[ 9] . Xp[13] .y =    8.20;             M[ 9] . Xp[13] .x=  0.193;          
 M[ 9] . Xp[14] .y =    8.80;             M[ 9] . Xp[14] .x=  0.181;          
 M[ 9] . Xp[15] .y =    9.40;             M[ 9] . Xp[15] .x=  0.192;          
 M[ 9] . Xp[16] .y =   10.00;             M[ 9] . Xp[16] .x=  0.179;          
 M[ 9] . Xp[17] .y =   10.60;             M[ 9] . Xp[17] .x=  0.179;          
 M[ 9] . Xp[18] .y =   11.20;             M[ 9] . Xp[18] .x=  0.164;          
 M[ 9] . Xp[19] .y =   11.80;             M[ 9] . Xp[19] .x=  0.160;          
 M[ 9] . Xp[20] .y =   12.40;             M[ 9] . Xp[20] .x=  0.145;          
 M[ 9] . Xp[21] .y =   13.00;             M[ 9] . Xp[21] .x=  0.147;          
 M[ 9] . Xp[22] .y =   13.60;             M[ 9] . Xp[22] .x=  0.140;          
 M[ 9] . Xp[23] .y =   14.20;             M[ 9] . Xp[23] .x=  0.139;          
 M[ 9] . Xp[24] .y =   14.80;             M[ 9] . Xp[24] .x=  0.135;          
 M[ 9] . Xp[25] .y =   15.40;             M[ 9] . Xp[25] .x=  0.135;          
 M[ 9] . Xp[26] .y =   16.00;             M[ 9] . Xp[26] .x=  0.134;          
 M[ 9] . Xp[27] .y =   16.60;             M[ 9] . Xp[27] .x=  0.130;          
 M[ 9] . Xp[28] .y =   17.20;             M[ 9] . Xp[28] .x=  0.134;          
 M[ 9] . Xp[29] .y =   17.80;             M[ 9] . Xp[29] .x=  0.119;          
 M[ 9] . Xp[30] .y =   18.40;             M[ 9] . Xp[30] .x=  0.123;          
 M[ 9] . Xp[31] .y =   19.00;             M[ 9] . Xp[31] .x=  0.109;          
 M[ 9] . Xp[32] .y =   19.60;             M[ 9] . Xp[32] .x=  0.112;          
 M[ 9] . Xp[33] .y =   20.20;             M[ 9] . Xp[33] .x=  0.094;          
 M[ 9] . Xp[34] .y =   20.80;             M[ 9] . Xp[34] .x=  0.077;          
 M[ 9] . Xp[35] .y =   21.40;             M[ 9] . Xp[35] .x=  0.070;          
 M[ 9] . Xp[36] .y =   22.00;             M[ 9] . Xp[36] .x=  0.049;          
 M[ 9] . Xp[37] .y =   22.60;             M[ 9] . Xp[37] .x=  0.035;          
 M[ 9] . Xp[38] .y =   23.20;             M[ 9] . Xp[38] .x=  0.026;          
 M[ 9] . Xp[39] .y =   23.80;             M[ 9] . Xp[39] .x=  0.016;          
 M[ 9] . Xp[40] .y =   24.40;             M[ 9] . Xp[40] .x=  0.005;          
 M[ 9] . Xp[41] .y =   25.00;             M[ 9] . Xp[41] .x=  0.004;          
 M[ 9] . Xp[42] .y =   25.60;             M[ 9] . Xp[42] .x=  0.001;          
 M[ 9] . Xp[43] .y =   26.20;             M[ 9] . Xp[43] .x=  0.002;          
   M[9].pathlength();
//station 11
 x[10] = 27.50;
  M[10].allocate(44);                              
  M[10].settype(SPLINE2D_LINEAR);                                                                               
                                                                         
 M[10] . Xp[ 0] .y =    0.40;             M[10] . Xp[ 0] .x=  0.768;          
 M[10] . Xp[ 1] .y =    1.00;             M[10] . Xp[ 1] .x=  0.773;          
 M[10] . Xp[ 2] .y =    1.60;             M[10] . Xp[ 2] .x=  0.720;          
 M[10] . Xp[ 3] .y =    2.20;             M[10] . Xp[ 3] .x=  0.589;          
 M[10] . Xp[ 4] .y =    2.80;             M[10] . Xp[ 4] .x=  0.526;          
 M[10] . Xp[ 5] .y =    3.40;             M[10] . Xp[ 5] .x=  0.469;          
 M[10] . Xp[ 6] .y =    4.00;             M[10] . Xp[ 6] .x=  0.384;          
 M[10] . Xp[ 7] .y =    4.60;             M[10] . Xp[ 7] .x=  0.349;          
 M[10] . Xp[ 8] .y =    5.20;             M[10] . Xp[ 8] .x=  0.318;          
 M[10] . Xp[ 9] .y =    5.80;             M[10] . Xp[ 9] .x=  0.272;          
 M[10] . Xp[10] .y =    6.40;             M[10] . Xp[10] .x=  0.241;          
 M[10] . Xp[11] .y =    7.00;             M[10] . Xp[11] .x=  0.213;          
 M[10] . Xp[12] .y =    7.60;             M[10] . Xp[12] .x=  0.208;          
 M[10] . Xp[13] .y =    8.20;             M[10] . Xp[13] .x=  0.191;          
 M[10] . Xp[14] .y =    8.80;             M[10] . Xp[14] .x=  0.185;          
 M[10] . Xp[15] .y =    9.40;             M[10] . Xp[15] .x=  0.183;          
 M[10] . Xp[16] .y =   10.00;             M[10] . Xp[16] .x=  0.192;          
 M[10] . Xp[17] .y =   10.60;             M[10] . Xp[17] .x=  0.179;          
 M[10] . Xp[18] .y =   11.20;             M[10] . Xp[18] .x=  0.170;          
 M[10] . Xp[19] .y =   11.80;             M[10] . Xp[19] .x=  0.154;          
 M[10] . Xp[20] .y =   12.40;             M[10] . Xp[20] .x=  0.163;          
 M[10] . Xp[21] .y =   13.00;             M[10] . Xp[21] .x=  0.151;          
 M[10] . Xp[22] .y =   13.60;             M[10] . Xp[22] .x=  0.148;          
 M[10] . Xp[23] .y =   14.20;             M[10] . Xp[23] .x=  0.140;          
 M[10] . Xp[24] .y =   14.80;             M[10] . Xp[24] .x=  0.138;          
 M[10] . Xp[25] .y =   15.40;             M[10] . Xp[25] .x=  0.124;          
 M[10] . Xp[26] .y =   16.00;             M[10] . Xp[26] .x=  0.125;          
 M[10] . Xp[27] .y =   16.60;             M[10] . Xp[27] .x=  0.124;          
 M[10] . Xp[28] .y =   17.20;             M[10] . Xp[28] .x=  0.130;          
 M[10] . Xp[29] .y =   17.80;             M[10] . Xp[29] .x=  0.122;          
 M[10] . Xp[30] .y =   18.40;             M[10] . Xp[30] .x=  0.119;          
 M[10] . Xp[31] .y =   19.00;             M[10] . Xp[31] .x=  0.107;          
 M[10] . Xp[32] .y =   19.60;             M[10] . Xp[32] .x=  0.103;          
 M[10] . Xp[33] .y =   20.20;             M[10] . Xp[33] .x=  0.091;          
 M[10] . Xp[34] .y =   20.80;             M[10] . Xp[34] .x=  0.077;          
 M[10] . Xp[35] .y =   21.40;             M[10] . Xp[35] .x=  0.062;          
 M[10] . Xp[36] .y =   22.00;             M[10] . Xp[36] .x=  0.043;          
 M[10] . Xp[37] .y =   22.60;             M[10] . Xp[37] .x=  0.029;          
 M[10] . Xp[38] .y =   23.20;             M[10] . Xp[38] .x=  0.023;          
 M[10] . Xp[39] .y =   23.80;             M[10] . Xp[39] .x=  0.014;          
 M[10] . Xp[40] .y =   24.40;             M[10] . Xp[40] .x=  0.003;          
 M[10] . Xp[41] .y =   25.00;             M[10] . Xp[41] .x=  0.003;          
 M[10] . Xp[42] .y =   25.60;             M[10] . Xp[42] .x=  0.002;          
 M[10] . Xp[43] .y =   26.20;             M[10] . Xp[43] .x=  0.001;          
  M[10].pathlength(); 


//station 12
 x[11] = 28.30;
  M[11].allocate(44);                              
  M[11].settype(SPLINE2D_LINEAR);                                                                              
                                                                          
 M[11] . Xp[ 0] .y =    0.40;             M[11] . Xp[ 0] .x=  0.707;          
 M[11] . Xp[ 1] .y =    1.00;             M[11] . Xp[ 1] .x=  0.684;          
 M[11] . Xp[ 2] .y =    1.60;             M[11] . Xp[ 2] .x=  0.634;          
 M[11] . Xp[ 3] .y =    2.20;             M[11] . Xp[ 3] .x=  0.596;          
 M[11] . Xp[ 4] .y =    2.80;             M[11] . Xp[ 4] .x=  0.541;          
 M[11] . Xp[ 5] .y =    3.40;             M[11] . Xp[ 5] .x=  0.444;          
 M[11] . Xp[ 6] .y =    4.00;             M[11] . Xp[ 6] .x=  0.374;          
 M[11] . Xp[ 7] .y =    4.60;             M[11] . Xp[ 7] .x=  0.326;          
 M[11] . Xp[ 8] .y =    5.20;             M[11] . Xp[ 8] .x=  0.303;          
 M[11] . Xp[ 9] .y =    5.80;             M[11] . Xp[ 9] .x=  0.281;          
 M[11] . Xp[10] .y =    6.40;             M[11] . Xp[10] .x=  0.242;          
 M[11] . Xp[11] .y =    7.00;             M[11] . Xp[11] .x=  0.234;          
 M[11] . Xp[12] .y =    7.60;             M[11] . Xp[12] .x=  0.193;          
 M[11] . Xp[13] .y =    8.20;             M[11] . Xp[13] .x=  0.192;          
 M[11] . Xp[14] .y =    8.80;             M[11] . Xp[14] .x=  0.182;          
 M[11] . Xp[15] .y =    9.40;             M[11] . Xp[15] .x=  0.179;          
 M[11] . Xp[16] .y =   10.00;             M[11] . Xp[16] .x=  0.175;          
 M[11] . Xp[17] .y =   10.60;             M[11] . Xp[17] .x=  0.176;          
 M[11] . Xp[18] .y =   11.20;             M[11] . Xp[18] .x=  0.167;          
 M[11] . Xp[19] .y =   11.80;             M[11] . Xp[19] .x=  0.152;          
 M[11] . Xp[20] .y =   12.40;             M[11] . Xp[20] .x=  0.153;          
 M[11] . Xp[21] .y =   13.00;             M[11] . Xp[21] .x=  0.143;          
 M[11] . Xp[22] .y =   13.60;             M[11] . Xp[22] .x=  0.137;          
 M[11] . Xp[23] .y =   14.20;             M[11] . Xp[23] .x=  0.131;          
 M[11] . Xp[24] .y =   14.80;             M[11] . Xp[24] .x=  0.132;          
 M[11] . Xp[25] .y =   15.40;             M[11] . Xp[25] .x=  0.133;          
 M[11] . Xp[26] .y =   16.00;             M[11] . Xp[26] .x=  0.130;          
 M[11] . Xp[27] .y =   16.60;             M[11] . Xp[27] .x=  0.129;          
 M[11] . Xp[28] .y =   17.20;             M[11] . Xp[28] .x=  0.125;          
 M[11] . Xp[29] .y =   17.80;             M[11] . Xp[29] .x=  0.120;          
 M[11] . Xp[30] .y =   18.40;             M[11] . Xp[30] .x=  0.105;          
 M[11] . Xp[31] .y =   19.00;             M[11] . Xp[31] .x=  0.099;          
 M[11] . Xp[32] .y =   19.60;             M[11] . Xp[32] .x=  0.088;          
 M[11] . Xp[33] .y =   20.20;             M[11] . Xp[33] .x=  0.087;          
 M[11] . Xp[34] .y =   20.80;             M[11] . Xp[34] .x=  0.067;          
 M[11] . Xp[35] .y =   21.40;             M[11] . Xp[35] .x=  0.054;          
 M[11] . Xp[36] .y =   22.00;             M[11] . Xp[36] .x=  0.047;          
 M[11] . Xp[37] .y =   22.60;             M[11] . Xp[37] .x=  0.029;          
 M[11] . Xp[38] .y =   23.20;             M[11] . Xp[38] .x=  0.022;          
 M[11] . Xp[39] .y =   23.80;             M[11] . Xp[39] .x=  0.015;          
 M[11] . Xp[40] .y =   24.40;             M[11] . Xp[40] .x=  0.003;          
 M[11] . Xp[41] .y =   25.00;             M[11] . Xp[41] .x=  0.002;          
 M[11] . Xp[42] .y =   25.60;             M[11] . Xp[42] .x=  0.000;          
 M[11] . Xp[43] .y =   26.20;             M[11] . Xp[43] .x=  0.000;          
  M[11].pathlength(); 

//station 13
 x[12] = 29.00;
  M[12].allocate(44);                              
  M[12].settype(SPLINE2D_LINEAR);                                                                                
                                                                         
 M[12] . Xp[ 0] .y =    0.40;             M[12] . Xp[ 0] .x=  0.737;          
 M[12] . Xp[ 1] .y =    1.00;             M[12] . Xp[ 1] .x=  0.709;          
 M[12] . Xp[ 2] .y =    1.60;             M[12] . Xp[ 2] .x=  0.681;          
 M[12] . Xp[ 3] .y =    2.20;             M[12] . Xp[ 3] .x=  0.590;          
 M[12] . Xp[ 4] .y =    2.80;             M[12] . Xp[ 4] .x=  0.531;          
 M[12] . Xp[ 5] .y =    3.40;             M[12] . Xp[ 5] .x=  0.454;          
 M[12] . Xp[ 6] .y =    4.00;             M[12] . Xp[ 6] .x=  0.365;          
 M[12] . Xp[ 7] .y =    4.60;             M[12] . Xp[ 7] .x=  0.352;          
 M[12] . Xp[ 8] .y =    5.20;             M[12] . Xp[ 8] .x=  0.316;          
 M[12] . Xp[ 9] .y =    5.80;             M[12] . Xp[ 9] .x=  0.283;          
 M[12] . Xp[10] .y =    6.40;             M[12] . Xp[10] .x=  0.255;          
 M[12] . Xp[11] .y =    7.00;             M[12] . Xp[11] .x=  0.236;          
 M[12] . Xp[12] .y =    7.60;             M[12] . Xp[12] .x=  0.212;          
 M[12] . Xp[13] .y =    8.20;             M[12] . Xp[13] .x=  0.199;          
 M[12] . Xp[14] .y =    8.80;             M[12] . Xp[14] .x=  0.191;          
 M[12] . Xp[15] .y =    9.40;             M[12] . Xp[15] .x=  0.188;          
 M[12] . Xp[16] .y =   10.00;             M[12] . Xp[16] .x=  0.184;          
 M[12] . Xp[17] .y =   10.60;             M[12] . Xp[17] .x=  0.170;          
 M[12] . Xp[18] .y =   11.20;             M[12] . Xp[18] .x=  0.158;          
 M[12] . Xp[19] .y =   11.80;             M[12] . Xp[19] .x=  0.146;          
 M[12] . Xp[20] .y =   12.40;             M[12] . Xp[20] .x=  0.151;          
 M[12] . Xp[21] .y =   13.00;             M[12] . Xp[21] .x=  0.140;          
 M[12] . Xp[22] .y =   13.60;             M[12] . Xp[22] .x=  0.132;          
 M[12] . Xp[23] .y =   14.20;             M[12] . Xp[23] .x=  0.133;          
 M[12] . Xp[24] .y =   14.80;             M[12] . Xp[24] .x=  0.131;          
 M[12] . Xp[25] .y =   15.40;             M[12] . Xp[25] .x=  0.139;          
 M[12] . Xp[26] .y =   16.00;             M[12] . Xp[26] .x=  0.135;          
 M[12] . Xp[27] .y =   16.60;             M[12] . Xp[27] .x=  0.139;          
 M[12] . Xp[28] .y =   17.20;             M[12] . Xp[28] .x=  0.121;          
 M[12] . Xp[29] .y =   17.80;             M[12] . Xp[29] .x=  0.124;          
 M[12] . Xp[30] .y =   18.40;             M[12] . Xp[30] .x=  0.108;          
 M[12] . Xp[31] .y =   19.00;             M[12] . Xp[31] .x=  0.119;          
 M[12] . Xp[32] .y =   19.60;             M[12] . Xp[32] .x=  0.093;          
 M[12] . Xp[33] .y =   20.20;             M[12] . Xp[33] .x=  0.089;          
 M[12] . Xp[34] .y =   20.80;             M[12] . Xp[34] .x=  0.076;          
 M[12] . Xp[35] .y =   21.40;             M[12] . Xp[35] .x=  0.058;          
 M[12] . Xp[36] .y =   22.00;             M[12] . Xp[36] .x=  0.052;          
 M[12] . Xp[37] .y =   22.60;             M[12] . Xp[37] .x=  0.030;          
 M[12] . Xp[38] .y =   23.20;             M[12] . Xp[38] .x=  0.018;          
 M[12] . Xp[39] .y =   23.80;             M[12] . Xp[39] .x=  0.014;          
 M[12] . Xp[40] .y =   24.40;             M[12] . Xp[40] .x=  0.005;          
 M[12] . Xp[41] .y =   25.00;             M[12] . Xp[41] .x=  0.003;          
 M[12] . Xp[42] .y =   25.60;             M[12] . Xp[42] .x=  0.000;          
 M[12] . Xp[43] .y =   26.20;             M[12] . Xp[43] .x=  0.003;          
  M[12].pathlength(); 

 cout.flush();
//station 14
 x[13] = 29.80;
  M[13].allocate(44);
  M[13].settype(SPLINE2D_LINEAR);
                                                                                                                                                       
 M[13] . Xp[ 0] .y =    0.40;             M[13] . Xp[ 0] .x=  0.706;          
 M[13] . Xp[ 1] .y =    1.00;             M[13] . Xp[ 1] .x=  0.701;          
 M[13] . Xp[ 2] .y =    1.60;             M[13] . Xp[ 2] .x=  0.622;          
 M[13] . Xp[ 3] .y =    2.20;             M[13] . Xp[ 3] .x=  0.560;          
 M[13] . Xp[ 4] .y =    2.80;             M[13] . Xp[ 4] .x=  0.502;          
 M[13] . Xp[ 5] .y =    3.40;             M[13] . Xp[ 5] .x=  0.429;          
 M[13] . Xp[ 6] .y =    4.00;             M[13] . Xp[ 6] .x=  0.377;          
 M[13] . Xp[ 7] .y =    4.60;             M[13] . Xp[ 7] .x=  0.346;          
 M[13] . Xp[ 8] .y =    5.20;             M[13] . Xp[ 8] .x=  0.316;          
 M[13] . Xp[ 9] .y =    5.80;             M[13] . Xp[ 9] .x=  0.283;          
 M[13] . Xp[10] .y =    6.40;             M[13] . Xp[10] .x=  0.251;          
 M[13] . Xp[11] .y =    7.00;             M[13] . Xp[11] .x=  0.239;          
 M[13] . Xp[12] .y =    7.60;             M[13] . Xp[12] .x=  0.213;          
 M[13] . Xp[13] .y =    8.20;             M[13] . Xp[13] .x=  0.203;          
 M[13] . Xp[14] .y =    8.80;             M[13] . Xp[14] .x=  0.190;          
 M[13] . Xp[15] .y =    9.40;             M[13] . Xp[15] .x=  0.167;          
 M[13] . Xp[16] .y =   10.00;             M[13] . Xp[16] .x=  0.179;          
 M[13] . Xp[17] .y =   10.60;             M[13] . Xp[17] .x=  0.168;          
 M[13] . Xp[18] .y =   11.20;             M[13] . Xp[18] .x=  0.159;          
 M[13] . Xp[19] .y =   11.80;             M[13] . Xp[19] .x=  0.146;          
 M[13] . Xp[20] .y =   12.40;             M[13] . Xp[20] .x=  0.150;          
 M[13] . Xp[21] .y =   13.00;             M[13] . Xp[21] .x=  0.130;          
 M[13] . Xp[22] .y =   13.60;             M[13] . Xp[22] .x=  0.139;          
 M[13] . Xp[23] .y =   14.20;             M[13] . Xp[23] .x=  0.120;          
 M[13] . Xp[24] .y =   14.80;             M[13] . Xp[24] .x=  0.124;          
 M[13] . Xp[25] .y =   15.40;             M[13] . Xp[25] .x=  0.128;          
 M[13] . Xp[26] .y =   16.00;             M[13] . Xp[26] .x=  0.124;          
 M[13] . Xp[27] .y =   16.60;             M[13] . Xp[27] .x=  0.128;          
 M[13] . Xp[28] .y =   17.20;             M[13] . Xp[28] .x=  0.117;          
 M[13] . Xp[29] .y =   17.80;             M[13] . Xp[29] .x=  0.107;          
 M[13] . Xp[30] .y =   18.40;             M[13] . Xp[30] .x=  0.106;          
 M[13] . Xp[31] .y =   19.00;             M[13] . Xp[31] .x=  0.103;          
 M[13] . Xp[32] .y =   19.60;             M[13] . Xp[32] .x=  0.097;          
 M[13] . Xp[33] .y =   20.20;             M[13] . Xp[33] .x=  0.082;          
 M[13] . Xp[34] .y =   20.80;             M[13] . Xp[34] .x=  0.061;          
 M[13] . Xp[35] .y =   21.40;             M[13] . Xp[35] .x=  0.055;          
 M[13] . Xp[36] .y =   22.00;             M[13] . Xp[36] .x=  0.042;          
 M[13] . Xp[37] .y =   22.60;             M[13] . Xp[37] .x=  0.026;          
 M[13] . Xp[38] .y =   23.20;             M[13] . Xp[38] .x=  0.017;          
 M[13] . Xp[39] .y =   23.80;             M[13] . Xp[39] .x=  0.009;          
 M[13] . Xp[40] .y =   24.40;             M[13] . Xp[40] .x=  0.003;          
 M[13] . Xp[41] .y =   25.00;             M[13] . Xp[41] .x=  0.002;          
 M[13] . Xp[42] .y =   25.60;             M[13] . Xp[42] .x=  0.003;          
 M[13] . Xp[43] .y =   26.20;             M[13] . Xp[43] .x=  0.003;          
 M[13].pathlength();
 
//station 15
 x[14] = 30.50;
  M[14].allocate(44);
  M[14].settype(SPLINE2D_LINEAR);                                                                          
                                                                          
 M[14] . Xp[ 0] .y =    0.40;             M[14] . Xp[ 0] .x=  0.711;          
 M[14] . Xp[ 1] .y =    1.00;             M[14] . Xp[ 1] .x=  0.660;          
 M[14] . Xp[ 2] .y =    1.60;             M[14] . Xp[ 2] .x=  0.633;          
 M[14] . Xp[ 3] .y =    2.20;             M[14] . Xp[ 3] .x=  0.566;          
 M[14] . Xp[ 4] .y =    2.80;             M[14] . Xp[ 4] .x=  0.480;          
 M[14] . Xp[ 5] .y =    3.40;             M[14] . Xp[ 5] .x=  0.428;          
 M[14] . Xp[ 6] .y =    4.00;             M[14] . Xp[ 6] .x=  0.393;          
 M[14] . Xp[ 7] .y =    4.60;             M[14] . Xp[ 7] .x=  0.334;          
 M[14] . Xp[ 8] .y =    5.20;             M[14] . Xp[ 8] .x=  0.319;          
 M[14] . Xp[ 9] .y =    5.80;             M[14] . Xp[ 9] .x=  0.289;          
 M[14] . Xp[10] .y =    6.40;             M[14] . Xp[10] .x=  0.271;          
 M[14] . Xp[11] .y =    7.00;             M[14] . Xp[11] .x=  0.243;          
 M[14] . Xp[12] .y =    7.60;             M[14] . Xp[12] .x=  0.231;          
 M[14] . Xp[13] .y =    8.20;             M[14] . Xp[13] .x=  0.206;          
 M[14] . Xp[14] .y =    8.80;             M[14] . Xp[14] .x=  0.182;          
 M[14] . Xp[15] .y =    9.40;             M[14] . Xp[15] .x=  0.176;          
 M[14] . Xp[16] .y =   10.00;             M[14] . Xp[16] .x=  0.169;          
 M[14] . Xp[17] .y =   10.60;             M[14] . Xp[17] .x=  0.157;          
 M[14] . Xp[18] .y =   11.20;             M[14] . Xp[18] .x=  0.166;          
 M[14] . Xp[19] .y =   11.80;             M[14] . Xp[19] .x=  0.156;          
 M[14] . Xp[20] .y =   12.40;             M[14] . Xp[20] .x=  0.136;          
 M[14] . Xp[21] .y =   13.00;             M[14] . Xp[21] .x=  0.147;          
 M[14] . Xp[22] .y =   13.60;             M[14] . Xp[22] .x=  0.132;          
 M[14] . Xp[23] .y =   14.20;             M[14] . Xp[23] .x=  0.135;          
 M[14] . Xp[24] .y =   14.80;             M[14] . Xp[24] .x=  0.126;          
 M[14] . Xp[25] .y =   15.40;             M[14] . Xp[25] .x=  0.134;          
 M[14] . Xp[26] .y =   16.00;             M[14] . Xp[26] .x=  0.117;          
 M[14] . Xp[27] .y =   16.60;             M[14] . Xp[27] .x=  0.123;          
 M[14] . Xp[28] .y =   17.20;             M[14] . Xp[28] .x=  0.123;          
 M[14] . Xp[29] .y =   17.80;             M[14] . Xp[29] .x=  0.116;          
 M[14] . Xp[30] .y =   18.40;             M[14] . Xp[30] .x=  0.106;          
 M[14] . Xp[31] .y =   19.00;             M[14] . Xp[31] .x=  0.092;          
 M[14] . Xp[32] .y =   19.60;             M[14] . Xp[32] .x=  0.090;          
 M[14] . Xp[33] .y =   20.20;             M[14] . Xp[33] .x=  0.076;          
 M[14] . Xp[34] .y =   20.80;             M[14] . Xp[34] .x=  0.060;          
 M[14] . Xp[35] .y =   21.40;             M[14] . Xp[35] .x=  0.046;          
 M[14] . Xp[36] .y =   22.00;             M[14] . Xp[36] .x=  0.034;          
 M[14] . Xp[37] .y =   22.60;             M[14] . Xp[37] .x=  0.018;          
 M[14] . Xp[38] .y =   23.20;             M[14] . Xp[38] .x=  0.017;          
 M[14] . Xp[39] .y =   23.80;             M[14] . Xp[39] .x=  0.009;          
 M[14] . Xp[40] .y =   24.40;             M[14] . Xp[40] .x=  0.006;          
 M[14] . Xp[41] .y =   25.00;             M[14] . Xp[41] .x=  0.004;          
 M[14] . Xp[42] .y =   25.60;             M[14] . Xp[42] .x=  0.001;          
 M[14] . Xp[43] .y =   26.20;             M[14] . Xp[43] .x=  0.004;          
 M[14].pathlength();
//station 16
 x[15] = 31.30;
  M[15].allocate(44);
  M[15].settype(SPLINE2D_LINEAR);                                                                             
                                                                          
 M[15] . Xp[ 0] .y =    0.40;             M[15] . Xp[ 0] .x=  0.653;          
 M[15] . Xp[ 1] .y =    1.00;             M[15] . Xp[ 1] .x=  0.691;          
 M[15] . Xp[ 2] .y =    1.60;             M[15] . Xp[ 2] .x=  0.588;          
 M[15] . Xp[ 3] .y =    2.20;             M[15] . Xp[ 3] .x=  0.536;          
 M[15] . Xp[ 4] .y =    2.80;             M[15] . Xp[ 4] .x=  0.493;          
 M[15] . Xp[ 5] .y =    3.40;             M[15] . Xp[ 5] .x=  0.440;          
 M[15] . Xp[ 6] .y =    4.00;             M[15] . Xp[ 6] .x=  0.380;          
 M[15] . Xp[ 7] .y =    4.60;             M[15] . Xp[ 7] .x=  0.324;          
 M[15] . Xp[ 8] .y =    5.20;             M[15] . Xp[ 8] .x=  0.296;          
 M[15] . Xp[ 9] .y =    5.80;             M[15] . Xp[ 9] .x=  0.262;          
 M[15] . Xp[10] .y =    6.40;             M[15] . Xp[10] .x=  0.260;          
 M[15] . Xp[11] .y =    7.00;             M[15] . Xp[11] .x=  0.252;          
 M[15] . Xp[12] .y =    7.60;             M[15] . Xp[12] .x=  0.229;          
 M[15] . Xp[13] .y =    8.20;             M[15] . Xp[13] .x=  0.198;          
 M[15] . Xp[14] .y =    8.80;             M[15] . Xp[14] .x=  0.208;          
 M[15] . Xp[15] .y =    9.40;             M[15] . Xp[15] .x=  0.186;          
 M[15] . Xp[16] .y =   10.00;             M[15] . Xp[16] .x=  0.170;          
 M[15] . Xp[17] .y =   10.60;             M[15] . Xp[17] .x=  0.168;          
 M[15] . Xp[18] .y =   11.20;             M[15] . Xp[18] .x=  0.171;          
 M[15] . Xp[19] .y =   11.80;             M[15] . Xp[19] .x=  0.155;          
 M[15] . Xp[20] .y =   12.40;             M[15] . Xp[20] .x=  0.148;          
 M[15] . Xp[21] .y =   13.00;             M[15] . Xp[21] .x=  0.131;          
 M[15] . Xp[22] .y =   13.60;             M[15] . Xp[22] .x=  0.135;          
 M[15] . Xp[23] .y =   14.20;             M[15] . Xp[23] .x=  0.136;          
 M[15] . Xp[24] .y =   14.80;             M[15] . Xp[24] .x=  0.128;          
 M[15] . Xp[25] .y =   15.40;             M[15] . Xp[25] .x=  0.126;          
 M[15] . Xp[26] .y =   16.00;             M[15] . Xp[26] .x=  0.116;          
 M[15] . Xp[27] .y =   16.60;             M[15] . Xp[27] .x=  0.118;          
 M[15] . Xp[28] .y =   17.20;             M[15] . Xp[28] .x=  0.115;          
 M[15] . Xp[29] .y =   17.80;             M[15] . Xp[29] .x=  0.104;          
 M[15] . Xp[30] .y =   18.40;             M[15] . Xp[30] .x=  0.097;          
 M[15] . Xp[31] .y =   19.00;             M[15] . Xp[31] .x=  0.078;          
 M[15] . Xp[32] .y =   19.60;             M[15] . Xp[32] .x=  0.076;          
 M[15] . Xp[33] .y =   20.20;             M[15] . Xp[33] .x=  0.072;          
 M[15] . Xp[34] .y =   20.80;             M[15] . Xp[34] .x=  0.054;          
 M[15] . Xp[35] .y =   21.40;             M[15] . Xp[35] .x=  0.042;          
 M[15] . Xp[36] .y =   22.00;             M[15] . Xp[36] .x=  0.037;          
 M[15] . Xp[37] .y =   22.60;             M[15] . Xp[37] .x=  0.019;          
 M[15] . Xp[38] .y =   23.20;             M[15] . Xp[38] .x=  0.012;          
 M[15] . Xp[39] .y =   23.80;             M[15] . Xp[39] .x=  0.012;          
 M[15] . Xp[40] .y =   24.40;             M[15] . Xp[40] .x=  0.005;          
 M[15] . Xp[41] .y =   25.00;             M[15] . Xp[41] .x=  0.003;          
 M[15] . Xp[42] .y =   25.60;             M[15] . Xp[42] .x=  0.005;          
 M[15] . Xp[43] .y =   26.20;             M[15] . Xp[43] .x=  0.002;          
  M[15].pathlength();

 cout.flush();
//station 17
 x[16] = 32.00;
  M[16].allocate(44);
  M[16].settype(SPLINE2D_LINEAR);                                                                             
                                                                          
 M[16] . Xp[ 0] .y =    0.40;             M[16] . Xp[ 0] .x=  0.650;          
 M[16] . Xp[ 1] .y =    1.00;             M[16] . Xp[ 1] .x=  0.641;          
 M[16] . Xp[ 2] .y =    1.60;             M[16] . Xp[ 2] .x=  0.580;          
 M[16] . Xp[ 3] .y =    2.20;             M[16] . Xp[ 3] .x=  0.543;          
 M[16] . Xp[ 4] .y =    2.80;             M[16] . Xp[ 4] .x=  0.528;          
 M[16] . Xp[ 5] .y =    3.40;             M[16] . Xp[ 5] .x=  0.438;          
 M[16] . Xp[ 6] .y =    4.00;             M[16] . Xp[ 6] .x=  0.370;          
 M[16] . Xp[ 7] .y =    4.60;             M[16] . Xp[ 7] .x=  0.359;          
 M[16] . Xp[ 8] .y =    5.20;             M[16] . Xp[ 8] .x=  0.289;          
 M[16] . Xp[ 9] .y =    5.80;             M[16] . Xp[ 9] .x=  0.316;          
 M[16] . Xp[10] .y =    6.40;             M[16] . Xp[10] .x=  0.264;          
 M[16] . Xp[11] .y =    7.00;             M[16] . Xp[11] .x=  0.251;          
 M[16] . Xp[12] .y =    7.60;             M[16] . Xp[12] .x=  0.231;          
 M[16] . Xp[13] .y =    8.20;             M[16] . Xp[13] .x=  0.211;          
 M[16] . Xp[14] .y =    8.80;             M[16] . Xp[14] .x=  0.201;          
 M[16] . Xp[15] .y =    9.40;             M[16] . Xp[15] .x=  0.177;          
 M[16] . Xp[16] .y =   10.00;             M[16] . Xp[16] .x=  0.169;          
 M[16] . Xp[17] .y =   10.60;             M[16] . Xp[17] .x=  0.168;          
 M[16] . Xp[18] .y =   11.20;             M[16] . Xp[18] .x=  0.168;          
 M[16] . Xp[19] .y =   11.80;             M[16] . Xp[19] .x=  0.137;          
 M[16] . Xp[20] .y =   12.40;             M[16] . Xp[20] .x=  0.142;          
 M[16] . Xp[21] .y =   13.00;             M[16] . Xp[21] .x=  0.133;          
 M[16] . Xp[22] .y =   13.60;             M[16] . Xp[22] .x=  0.123;          
 M[16] . Xp[23] .y =   14.20;             M[16] . Xp[23] .x=  0.120;          
 M[16] . Xp[24] .y =   14.80;             M[16] . Xp[24] .x=  0.113;          
 M[16] . Xp[25] .y =   15.40;             M[16] . Xp[25] .x=  0.125;          
 M[16] . Xp[26] .y =   16.00;             M[16] . Xp[26] .x=  0.107;          
 M[16] . Xp[27] .y =   16.60;             M[16] . Xp[27] .x=  0.112;          
 M[16] . Xp[28] .y =   17.20;             M[16] . Xp[28] .x=  0.114;          
 M[16] . Xp[29] .y =   17.80;             M[16] . Xp[29] .x=  0.104;          
 M[16] . Xp[30] .y =   18.40;             M[16] . Xp[30] .x=  0.096;          
 M[16] . Xp[31] .y =   19.00;             M[16] . Xp[31] .x=  0.083;          
 M[16] . Xp[32] .y =   19.60;             M[16] . Xp[32] .x=  0.078;          
 M[16] . Xp[33] .y =   20.20;             M[16] . Xp[33] .x=  0.069;          
 M[16] . Xp[34] .y =   20.80;             M[16] . Xp[34] .x=  0.052;          
 M[16] . Xp[35] .y =   21.40;             M[16] . Xp[35] .x=  0.040;          
 M[16] . Xp[36] .y =   22.00;             M[16] . Xp[36] .x=  0.026;          
 M[16] . Xp[37] .y =   22.60;             M[16] . Xp[37] .x=  0.023;          
 M[16] . Xp[38] .y =   23.20;             M[16] . Xp[38] .x=  0.013;          
 M[16] . Xp[39] .y =   23.80;             M[16] . Xp[39] .x=  0.009;          
 M[16] . Xp[40] .y =   24.40;             M[16] . Xp[40] .x=  0.004;          
 M[16] . Xp[41] .y =   25.00;             M[16] . Xp[41] .x=  0.010;          
 M[16] . Xp[42] .y =   25.60;             M[16] . Xp[42] .x=  0.005;          
 M[16] . Xp[43] .y =   26.20;             M[16] . Xp[43] .x=  0.005;          
 M[16].pathlength();

//station 18
 x[17] = 32.80;
  M[17].allocate(44);
  M[17].settype(SPLINE2D_LINEAR);                                                                            
                                                                          
 M[17] . Xp[ 0] .y =    0.40;             M[17] . Xp[ 0] .x=  0.617;          
 M[17] . Xp[ 1] .y =    1.00;             M[17] . Xp[ 1] .x=  0.590;          
 M[17] . Xp[ 2] .y =    1.60;             M[17] . Xp[ 2] .x=  0.595;          
 M[17] . Xp[ 3] .y =    2.20;             M[17] . Xp[ 3] .x=  0.546;          
 M[17] . Xp[ 4] .y =    2.80;             M[17] . Xp[ 4] .x=  0.462;          
 M[17] . Xp[ 5] .y =    3.40;             M[17] . Xp[ 5] .x=  0.400;          
 M[17] . Xp[ 6] .y =    4.00;             M[17] . Xp[ 6] .x=  0.340;          
 M[17] . Xp[ 7] .y =    4.60;             M[17] . Xp[ 7] .x=  0.321;          
 M[17] . Xp[ 8] .y =    5.20;             M[17] . Xp[ 8] .x=  0.329;          
 M[17] . Xp[ 9] .y =    5.80;             M[17] . Xp[ 9] .x=  0.296;          
 M[17] . Xp[10] .y =    6.40;             M[17] . Xp[10] .x=  0.263;          
 M[17] . Xp[11] .y =    7.00;             M[17] . Xp[11] .x=  0.241;          
 M[17] . Xp[12] .y =    7.60;             M[17] . Xp[12] .x=  0.208;          
 M[17] . Xp[13] .y =    8.20;             M[17] . Xp[13] .x=  0.206;          
 M[17] . Xp[14] .y =    8.80;             M[17] . Xp[14] .x=  0.183;          
 M[17] . Xp[15] .y =    9.40;             M[17] . Xp[15] .x=  0.158;          
 M[17] . Xp[16] .y =   10.00;             M[17] . Xp[16] .x=  0.172;          
 M[17] . Xp[17] .y =   10.60;             M[17] . Xp[17] .x=  0.156;          
 M[17] . Xp[18] .y =   11.20;             M[17] . Xp[18] .x=  0.125;          
 M[17] . Xp[19] .y =   11.80;             M[17] . Xp[19] .x=  0.141;          
 M[17] . Xp[20] .y =   12.40;             M[17] . Xp[20] .x=  0.119;          
 M[17] . Xp[21] .y =   13.00;             M[17] . Xp[21] .x=  0.118;          
 M[17] . Xp[22] .y =   13.60;             M[17] . Xp[22] .x=  0.111;          
 M[17] . Xp[23] .y =   14.20;             M[17] . Xp[23] .x=  0.116;          
 M[17] . Xp[24] .y =   14.80;             M[17] . Xp[24] .x=  0.119;          
 M[17] . Xp[25] .y =   15.40;             M[17] . Xp[25] .x=  0.108;          
 M[17] . Xp[26] .y =   16.00;             M[17] . Xp[26] .x=  0.116;          
 M[17] . Xp[27] .y =   16.60;             M[17] . Xp[27] .x=  0.098;          
 M[17] . Xp[28] .y =   17.20;             M[17] . Xp[28] .x=  0.104;          
 M[17] . Xp[29] .y =   17.80;             M[17] . Xp[29] .x=  0.094;          
 M[17] . Xp[30] .y =   18.40;             M[17] . Xp[30] .x=  0.089;          
 M[17] . Xp[31] .y =   19.00;             M[17] . Xp[31] .x=  0.079;          
 M[17] . Xp[32] .y =   19.60;             M[17] . Xp[32] .x=  0.060;          
 M[17] . Xp[33] .y =   20.20;             M[17] . Xp[33] .x=  0.064;          
 M[17] . Xp[34] .y =   20.80;             M[17] . Xp[34] .x=  0.042;          
 M[17] . Xp[35] .y =   21.40;             M[17] . Xp[35] .x=  0.037;          
 M[17] . Xp[36] .y =   22.00;             M[17] . Xp[36] .x=  0.022;          
 M[17] . Xp[37] .y =   22.60;             M[17] . Xp[37] .x=  0.019;          
 M[17] . Xp[38] .y =   23.20;             M[17] . Xp[38] .x=  0.012;          
 M[17] . Xp[39] .y =   23.80;             M[17] . Xp[39] .x=  0.003;          
 M[17] . Xp[40] .y =   24.40;             M[17] . Xp[40] .x=  0.000;          
 M[17] . Xp[41] .y =   25.00;             M[17] . Xp[41] .x=  0.010;          
 M[17] . Xp[42] .y =   25.60;             M[17] . Xp[42] .x=  0.007;          
 M[17] . Xp[43] .y =   26.20;             M[17] . Xp[43] .x=  0.002;          
  M[17].pathlength();
//station 19
 x[18] = 40.0;
  M[18].allocate(43);
  M[18].settype(SPLINE2D_LINEAR);                                                                             
                                                                        
 M[18] . Xp[ 0] .y =    0.20;             M[18] . Xp[ 0] .x=  0.540;          
 M[18] . Xp[ 1] .y =    0.80;             M[18] . Xp[ 1] .x=  0.499;          
 M[18] . Xp[ 2] .y =    1.40;             M[18] . Xp[ 2] .x=  0.459;          
 M[18] . Xp[ 3] .y =    2.00;             M[18] . Xp[ 3] .x=  0.459;          
 M[18] . Xp[ 4] .y =    2.60;             M[18] . Xp[ 4] .x=  0.417;          
 M[18] . Xp[ 5] .y =    3.20;             M[18] . Xp[ 5] .x=  0.412;          
 M[18] . Xp[ 6] .y =    3.80;             M[18] . Xp[ 6] .x=  0.339;          
 M[18] . Xp[ 7] .y =    4.40;             M[18] . Xp[ 7] .x=  0.330;          
 M[18] . Xp[ 8] .y =    5.00;             M[18] . Xp[ 8] .x=  0.319;          
 M[18] . Xp[ 9] .y =    5.60;             M[18] . Xp[ 9] .x=  0.296;          
 M[18] . Xp[10] .y =    6.20;             M[18] . Xp[10] .x=  0.261;          
 M[18] . Xp[11] .y =    6.80;             M[18] . Xp[11] .x=  0.219;          
 M[18] . Xp[12] .y =    7.40;             M[18] . Xp[12] .x=  0.220;          
 M[18] . Xp[13] .y =    8.00;             M[18] . Xp[13] .x=  0.171;          
 M[18] . Xp[14] .y =    8.60;             M[18] . Xp[14] .x=  0.157;          
 M[18] . Xp[15] .y =    9.20;             M[18] . Xp[15] .x=  0.125;          
 M[18] . Xp[16] .y =    9.80;             M[18] . Xp[16] .x=  0.132;          
 M[18] . Xp[17] .y =   10.40;             M[18] . Xp[17] .x=  0.109;          
 M[18] . Xp[18] .y =   11.00;             M[18] . Xp[18] .x=  0.105;          
 M[18] . Xp[19] .y =   11.60;             M[18] . Xp[19] .x=  0.097;          
 M[18] . Xp[20] .y =   12.20;             M[18] . Xp[20] .x=  0.077;          
 M[18] . Xp[21] .y =   12.80;             M[18] . Xp[21] .x=  0.092;          
 M[18] . Xp[22] .y =   13.40;             M[18] . Xp[22] .x=  0.071;          
 M[18] . Xp[23] .y =   14.00;             M[18] . Xp[23] .x=  0.073;          
 M[18] . Xp[24] .y =   14.60;             M[18] . Xp[24] .x=  0.071;          
 M[18] . Xp[25] .y =   15.20;             M[18] . Xp[25] .x=  0.059;          
 M[18] . Xp[26] .y =   15.80;             M[18] . Xp[26] .x=  0.066;          
 M[18] . Xp[27] .y =   16.40;             M[18] . Xp[27] .x=  0.054;          
 M[18] . Xp[28] .y =   17.00;             M[18] . Xp[28] .x=  0.052;          
 M[18] . Xp[29] .y =   17.60;             M[18] . Xp[29] .x=  0.044;          
 M[18] . Xp[30] .y =   18.20;             M[18] . Xp[30] .x=  0.047;          
 M[18] . Xp[31] .y =   18.80;             M[18] . Xp[31] .x=  0.042;          
 M[18] . Xp[32] .y =   19.40;             M[18] . Xp[32] .x=  0.025;          
 M[18] . Xp[33] .y =   20.00;             M[18] . Xp[33] .x=  0.021;          
 M[18] . Xp[34] .y =   20.60;             M[18] . Xp[34] .x=  0.016;          
 M[18] . Xp[35] .y =   21.20;             M[18] . Xp[35] .x=  0.014;          
 M[18] . Xp[36] .y =   21.80;             M[18] . Xp[36] .x=  0.012;          
 M[18] . Xp[37] .y =   22.40;             M[18] . Xp[37] .x=  0.011;          
 M[18] . Xp[38] .y =   23.00;             M[18] . Xp[38] .x=  0.004;          
 M[18] . Xp[39] .y =   23.60;             M[18] . Xp[39] .x=  0.004;          
 M[18] . Xp[40] .y =   24.20;             M[18] . Xp[40] .x=  0.002;          
 M[18] . Xp[41] .y =   24.80;             M[18] . Xp[41] .x=  0.009;          
 M[18] . Xp[42] .y =   25.40;             M[18] . Xp[42] .x=  0.006;          
  M[18].pathlength();

//station 20
 x[19] = 40.8;
  M[19].allocate(43);
  M[19].settype(SPLINE2D_LINEAR);                                                                             
                                                                         
 M[19] . Xp[ 0] .y =    0.20;             M[19] . Xp[ 0] .x=  0.559;          
 M[19] . Xp[ 1] .y =    0.80;             M[19] . Xp[ 1] .x=  0.522;          
 M[19] . Xp[ 2] .y =    1.40;             M[19] . Xp[ 2] .x=  0.489;          
 M[19] . Xp[ 3] .y =    2.00;             M[19] . Xp[ 3] .x=  0.450;          
 M[19] . Xp[ 4] .y =    2.60;             M[19] . Xp[ 4] .x=  0.431;          
 M[19] . Xp[ 5] .y =    3.20;             M[19] . Xp[ 5] .x=  0.392;          
 M[19] . Xp[ 6] .y =    3.80;             M[19] . Xp[ 6] .x=  0.373;          
 M[19] . Xp[ 7] .y =    4.40;             M[19] . Xp[ 7] .x=  0.329;          
 M[19] . Xp[ 8] .y =    5.00;             M[19] . Xp[ 8] .x=  0.318;          
 M[19] . Xp[ 9] .y =    5.60;             M[19] . Xp[ 9] .x=  0.295;          
 M[19] . Xp[10] .y =    6.20;             M[19] . Xp[10] .x=  0.254;          
 M[19] . Xp[11] .y =    6.80;             M[19] . Xp[11] .x=  0.230;          
 M[19] . Xp[12] .y =    7.40;             M[19] . Xp[12] .x=  0.202;          
 M[19] . Xp[13] .y =    8.00;             M[19] . Xp[13] .x=  0.190;          
 M[19] . Xp[14] .y =    8.60;             M[19] . Xp[14] .x=  0.160;          
 M[19] . Xp[15] .y =    9.20;             M[19] . Xp[15] .x=  0.150;          
 M[19] . Xp[16] .y =    9.80;             M[19] . Xp[16] .x=  0.133;          
 M[19] . Xp[17] .y =   10.40;             M[19] . Xp[17] .x=  0.105;          
 M[19] . Xp[18] .y =   11.00;             M[19] . Xp[18] .x=  0.096;          
 M[19] . Xp[19] .y =   11.60;             M[19] . Xp[19] .x=  0.095;          
 M[19] . Xp[20] .y =   12.20;             M[19] . Xp[20] .x=  0.078;          
 M[19] . Xp[21] .y =   12.80;             M[19] . Xp[21] .x=  0.075;          
 M[19] . Xp[22] .y =   13.40;             M[19] . Xp[22] .x=  0.078;          
 M[19] . Xp[23] .y =   14.00;             M[19] . Xp[23] .x=  0.073;          
 M[19] . Xp[24] .y =   14.60;             M[19] . Xp[24] .x=  0.060;          
 M[19] . Xp[25] .y =   15.20;             M[19] . Xp[25] .x=  0.053;          
 M[19] . Xp[26] .y =   15.80;             M[19] . Xp[26] .x=  0.060;          
 M[19] . Xp[27] .y =   16.40;             M[19] . Xp[27] .x=  0.054;          
 M[19] . Xp[28] .y =   17.00;             M[19] . Xp[28] .x=  0.051;          
 M[19] . Xp[29] .y =   17.60;             M[19] . Xp[29] .x=  0.052;          
 M[19] . Xp[30] .y =   18.20;             M[19] . Xp[30] .x=  0.035;          
 M[19] . Xp[31] .y =   18.80;             M[19] . Xp[31] .x=  0.034;          
 M[19] . Xp[32] .y =   19.40;             M[19] . Xp[32] .x=  0.032;          
 M[19] . Xp[33] .y =   20.00;             M[19] . Xp[33] .x=  0.022;          
 M[19] . Xp[34] .y =   20.60;             M[19] . Xp[34] .x=  0.014;          
 M[19] . Xp[35] .y =   21.20;             M[19] . Xp[35] .x=  0.018;          
 M[19] . Xp[36] .y =   21.80;             M[19] . Xp[36] .x=  0.013;          
 M[19] . Xp[37] .y =   22.40;             M[19] . Xp[37] .x=  0.008;          
 M[19] . Xp[38] .y =   23.00;             M[19] . Xp[38] .x=  0.003;          
 M[19] . Xp[39] .y =   23.60;             M[19] . Xp[39] .x=  0.002;          
 M[19] . Xp[40] .y =   24.20;             M[19] . Xp[40] .x=  0.001;          
 M[19] . Xp[41] .y =   24.80;             M[19] . Xp[41] .x=  0.003;          
 M[19] . Xp[42] .y =   25.40;             M[19] . Xp[42] .x=  0.003;          
  M[19].pathlength();
//station 21
 x[20] = 41.5;
  M[20].allocate(43);
  M[20].settype(SPLINE2D_LINEAR);                                                                                                                                                      
 M[20] . Xp[ 0] .y =    0.20;             M[20] . Xp[ 0] .x=  0.555;          
 M[20] . Xp[ 1] .y =    0.80;             M[20] . Xp[ 1] .x=  0.562;          
 M[20] . Xp[ 2] .y =    1.40;             M[20] . Xp[ 2] .x=  0.499;          
 M[20] . Xp[ 3] .y =    2.00;             M[20] . Xp[ 3] .x=  0.487;          
 M[20] . Xp[ 4] .y =    2.60;             M[20] . Xp[ 4] .x=  0.434;          
 M[20] . Xp[ 5] .y =    3.20;             M[20] . Xp[ 5] .x=  0.413;          
 M[20] . Xp[ 6] .y =    3.80;             M[20] . Xp[ 6] .x=  0.379;          
 M[20] . Xp[ 7] .y =    4.40;             M[20] . Xp[ 7] .x=  0.339;          
 M[20] . Xp[ 8] .y =    5.00;             M[20] . Xp[ 8] .x=  0.335;          
 M[20] . Xp[ 9] .y =    5.60;             M[20] . Xp[ 9] .x=  0.288;          
 M[20] . Xp[10] .y =    6.20;             M[20] . Xp[10] .x=  0.275;          
 M[20] . Xp[11] .y =    6.80;             M[20] . Xp[11] .x=  0.240;          
 M[20] . Xp[12] .y =    7.40;             M[20] . Xp[12] .x=  0.208;          
 M[20] . Xp[13] .y =    8.00;             M[20] . Xp[13] .x=  0.193;          
 M[20] . Xp[14] .y =    8.60;             M[20] . Xp[14] .x=  0.159;          
 M[20] . Xp[15] .y =    9.20;             M[20] . Xp[15] .x=  0.151;          
 M[20] . Xp[16] .y =    9.80;             M[20] . Xp[16] .x=  0.128;          
 M[20] . Xp[17] .y =   10.40;             M[20] . Xp[17] .x=  0.114;          
 M[20] . Xp[18] .y =   11.00;             M[20] . Xp[18] .x=  0.098;          
 M[20] . Xp[19] .y =   11.60;             M[20] . Xp[19] .x=  0.096;          
 M[20] . Xp[20] .y =   12.20;             M[20] . Xp[20] .x=  0.079;          
 M[20] . Xp[21] .y =   12.80;             M[20] . Xp[21] .x=  0.087;          
 M[20] . Xp[22] .y =   13.40;             M[20] . Xp[22] .x=  0.070;          
 M[20] . Xp[23] .y =   14.00;             M[20] . Xp[23] .x=  0.063;          
 M[20] . Xp[24] .y =   14.60;             M[20] . Xp[24] .x=  0.055;          
 M[20] . Xp[25] .y =   15.20;             M[20] . Xp[25] .x=  0.058;          
 M[20] . Xp[26] .y =   15.80;             M[20] . Xp[26] .x=  0.056;          
 M[20] . Xp[27] .y =   16.40;             M[20] . Xp[27] .x=  0.056;          
 M[20] . Xp[28] .y =   17.00;             M[20] . Xp[28] .x=  0.051;          
 M[20] . Xp[29] .y =   17.60;             M[20] . Xp[29] .x=  0.045;          
 M[20] . Xp[30] .y =   18.20;             M[20] . Xp[30] .x=  0.038;          
 M[20] . Xp[31] .y =   18.80;             M[20] . Xp[31] .x=  0.033;          
 M[20] . Xp[32] .y =   19.40;             M[20] . Xp[32] .x=  0.030;          
 M[20] . Xp[33] .y =   20.00;             M[20] . Xp[33] .x=  0.021;          
 M[20] . Xp[34] .y =   20.60;             M[20] . Xp[34] .x=  0.016;          
 M[20] . Xp[35] .y =   21.20;             M[20] . Xp[35] .x=  0.009;          
 M[20] . Xp[36] .y =   21.80;             M[20] . Xp[36] .x=  0.008;          
 M[20] . Xp[37] .y =   22.40;             M[20] . Xp[37] .x=  0.004;          
 M[20] . Xp[38] .y =   23.00;             M[20] . Xp[38] .x=  0.002;          
 M[20] . Xp[39] .y =   23.60;             M[20] . Xp[39] .x=  0.001;          
 M[20] . Xp[40] .y =   24.20;             M[20] . Xp[40] .x=  0.005;          
 M[20] . Xp[41] .y =   24.80;             M[20] . Xp[41] .x=  0.004;          
 M[20] . Xp[42] .y =   25.40;             M[20] . Xp[42] .x=  0.005;          
 M[20].pathlength();
//station 22
 x[21] = 42.3;
  M[21].allocate(43);
  M[21].settype(SPLINE2D_LINEAR);                                                                           
                                                                          
 M[21] . Xp[ 0] .y =    0.20;             M[21] . Xp[ 0] .x=  0.602;          
 M[21] . Xp[ 1] .y =    0.80;             M[21] . Xp[ 1] .x=  0.514;          
 M[21] . Xp[ 2] .y =    1.40;             M[21] . Xp[ 2] .x=  0.522;          
 M[21] . Xp[ 3] .y =    2.00;             M[21] . Xp[ 3] .x=  0.462;          
 M[21] . Xp[ 4] .y =    2.60;             M[21] . Xp[ 4] .x=  0.439;          
 M[21] . Xp[ 5] .y =    3.20;             M[21] . Xp[ 5] .x=  0.410;          
 M[21] . Xp[ 6] .y =    3.80;             M[21] . Xp[ 6] .x=  0.409;          
 M[21] . Xp[ 7] .y =    4.40;             M[21] . Xp[ 7] .x=  0.362;          
 M[21] . Xp[ 8] .y =    5.00;             M[21] . Xp[ 8] .x=  0.358;          
 M[21] . Xp[ 9] .y =    5.60;             M[21] . Xp[ 9] .x=  0.299;          
 M[21] . Xp[10] .y =    6.20;             M[21] . Xp[10] .x=  0.278;          
 M[21] . Xp[11] .y =    6.80;             M[21] . Xp[11] .x=  0.245;          
 M[21] . Xp[12] .y =    7.40;             M[21] . Xp[12] .x=  0.210;          
 M[21] . Xp[13] .y =    8.00;             M[21] . Xp[13] .x=  0.189;          
 M[21] . Xp[14] .y =    8.60;             M[21] . Xp[14] .x=  0.173;          
 M[21] . Xp[15] .y =    9.20;             M[21] . Xp[15] .x=  0.145;          
 M[21] . Xp[16] .y =    9.80;             M[21] . Xp[16] .x=  0.140;          
 M[21] . Xp[17] .y =   10.40;             M[21] . Xp[17] .x=  0.113;          
 M[21] . Xp[18] .y =   11.00;             M[21] . Xp[18] .x=  0.079;          
 M[21] . Xp[19] .y =   11.60;             M[21] . Xp[19] .x=  0.080;          
 M[21] . Xp[20] .y =   12.20;             M[21] . Xp[20] .x=  0.077;          
 M[21] . Xp[21] .y =   12.80;             M[21] . Xp[21] .x=  0.080;          
 M[21] . Xp[22] .y =   13.40;             M[21] . Xp[22] .x=  0.076;          
 M[21] . Xp[23] .y =   14.00;             M[21] . Xp[23] .x=  0.067;          
 M[21] . Xp[24] .y =   14.60;             M[21] . Xp[24] .x=  0.059;          
 M[21] . Xp[25] .y =   15.20;             M[21] . Xp[25] .x=  0.050;          
 M[21] . Xp[26] .y =   15.80;             M[21] . Xp[26] .x=  0.051;          
 M[21] . Xp[27] .y =   16.40;             M[21] . Xp[27] .x=  0.051;          
 M[21] . Xp[28] .y =   17.00;             M[21] . Xp[28] .x=  0.050;          
 M[21] . Xp[29] .y =   17.60;             M[21] . Xp[29] .x=  0.042;          
 M[21] . Xp[30] .y =   18.20;             M[21] . Xp[30] .x=  0.031;          
 M[21] . Xp[31] .y =   18.80;             M[21] . Xp[31] .x=  0.026;          
 M[21] . Xp[32] .y =   19.40;             M[21] . Xp[32] .x=  0.026;          
 M[21] . Xp[33] .y =   20.00;             M[21] . Xp[33] .x=  0.023;          
 M[21] . Xp[34] .y =   20.60;             M[21] . Xp[34] .x=  0.013;          
 M[21] . Xp[35] .y =   21.20;             M[21] . Xp[35] .x=  0.007;          
 M[21] . Xp[36] .y =   21.80;             M[21] . Xp[36] .x=  0.006;          
 M[21] . Xp[37] .y =   22.40;             M[21] . Xp[37] .x=  0.000;          
 M[21] . Xp[38] .y =   23.00;             M[21] . Xp[38] .x=  0.003;          
 M[21] . Xp[39] .y =   23.60;             M[21] . Xp[39] .x=  0.003;          
 M[21] . Xp[40] .y =   24.20;             M[21] . Xp[40] .x=  0.004;          
 M[21] . Xp[41] .y =   24.80;             M[21] . Xp[41] .x=  0.004;          
 M[21] . Xp[42] .y =   25.40;             M[21] . Xp[42] .x=  0.005;          
  M[21].pathlength();
//station 23
 x[22] = 43.0;
  M[22].allocate(43);
  M[22].settype(SPLINE2D_LINEAR);                                                                         
                                                                     
 M[22] . Xp[ 0] .y =    0.20;             M[22] . Xp[ 0] .x=  0.567;          
 M[22] . Xp[ 1] .y =    0.80;             M[22] . Xp[ 1] .x=  0.546;          
 M[22] . Xp[ 2] .y =    1.40;             M[22] . Xp[ 2] .x=  0.513;          
 M[22] . Xp[ 3] .y =    2.00;             M[22] . Xp[ 3] .x=  0.493;          
 M[22] . Xp[ 4] .y =    2.60;             M[22] . Xp[ 4] .x=  0.427;          
 M[22] . Xp[ 5] .y =    3.20;             M[22] . Xp[ 5] .x=  0.428;          
 M[22] . Xp[ 6] .y =    3.80;             M[22] . Xp[ 6] .x=  0.392;          
 M[22] . Xp[ 7] .y =    4.40;             M[22] . Xp[ 7] .x=  0.394;          
 M[22] . Xp[ 8] .y =    5.00;             M[22] . Xp[ 8] .x=  0.359;          
 M[22] . Xp[ 9] .y =    5.60;             M[22] . Xp[ 9] .x=  0.325;          
 M[22] . Xp[10] .y =    6.20;             M[22] . Xp[10] .x=  0.285;          
 M[22] . Xp[11] .y =    6.80;             M[22] . Xp[11] .x=  0.234;          
 M[22] . Xp[12] .y =    7.40;             M[22] . Xp[12] .x=  0.236;          
 M[22] . Xp[13] .y =    8.00;             M[22] . Xp[13] .x=  0.197;          
 M[22] . Xp[14] .y =    8.60;             M[22] . Xp[14] .x=  0.179;          
 M[22] . Xp[15] .y =    9.20;             M[22] . Xp[15] .x=  0.147;          
 M[22] . Xp[16] .y =    9.80;             M[22] . Xp[16] .x=  0.136;          
 M[22] . Xp[17] .y =   10.40;             M[22] . Xp[17] .x=  0.114;          
 M[22] . Xp[18] .y =   11.00;             M[22] . Xp[18] .x=  0.091;          
 M[22] . Xp[19] .y =   11.60;             M[22] . Xp[19] .x=  0.088;          
 M[22] . Xp[20] .y =   12.20;             M[22] . Xp[20] .x=  0.091;          
 M[22] . Xp[21] .y =   12.80;             M[22] . Xp[21] .x=  0.089;          
 M[22] . Xp[22] .y =   13.40;             M[22] . Xp[22] .x=  0.072;          
 M[22] . Xp[23] .y =   14.00;             M[22] . Xp[23] .x=  0.062;          
 M[22] . Xp[24] .y =   14.60;             M[22] . Xp[24] .x=  0.059;          
 M[22] . Xp[25] .y =   15.20;             M[22] . Xp[25] .x=  0.066;          
 M[22] . Xp[26] .y =   15.80;             M[22] . Xp[26] .x=  0.052;          
 M[22] . Xp[27] .y =   16.40;             M[22] . Xp[27] .x=  0.058;          
 M[22] . Xp[28] .y =   17.00;             M[22] . Xp[28] .x=  0.058;          
 M[22] . Xp[29] .y =   17.60;             M[22] . Xp[29] .x=  0.043;          
 M[22] . Xp[30] .y =   18.20;             M[22] . Xp[30] .x=  0.041;          
 M[22] . Xp[31] .y =   18.80;             M[22] . Xp[31] .x=  0.030;          
 M[22] . Xp[32] .y =   19.40;             M[22] . Xp[32] .x=  0.023;          
 M[22] . Xp[33] .y =   20.00;             M[22] . Xp[33] .x=  0.020;          
 M[22] . Xp[34] .y =   20.60;             M[22] . Xp[34] .x=  0.014;          
 M[22] . Xp[35] .y =   21.20;             M[22] . Xp[35] .x=  0.009;          
 M[22] . Xp[36] .y =   21.80;             M[22] . Xp[36] .x=  0.003;          
 M[22] . Xp[37] .y =   22.40;             M[22] . Xp[37] .x=  0.004;          
 M[22] . Xp[38] .y =   23.00;             M[22] . Xp[38] .x=  0.004;          
 M[22] . Xp[39] .y =   23.60;             M[22] . Xp[39] .x=  0.004;          
 M[22] . Xp[40] .y =   24.20;             M[22] . Xp[40] .x=  0.002;          
 M[22] . Xp[41] .y =   24.80;             M[22] . Xp[41] .x=  0.005;          
 M[22] . Xp[42] .y =   25.40;             M[22] . Xp[42] .x=  0.003;          
  M[22].pathlength();
//station 24
 x[23] = 43.8;
  M[23].allocate(43);
  M[23].settype(SPLINE2D_LINEAR);                                                                                                                                                       
 M[23] . Xp[ 0] .y =    0.20;             M[23] . Xp[ 0] .x=  0.537;          
 M[23] . Xp[ 1] .y =    0.80;             M[23] . Xp[ 1] .x=  0.536;          
 M[23] . Xp[ 2] .y =    1.40;             M[23] . Xp[ 2] .x=  0.496;          
 M[23] . Xp[ 3] .y =    2.00;             M[23] . Xp[ 3] .x=  0.483;          
 M[23] . Xp[ 4] .y =    2.60;             M[23] . Xp[ 4] .x=  0.459;          
 M[23] . Xp[ 5] .y =    3.20;             M[23] . Xp[ 5] .x=  0.410;          
 M[23] . Xp[ 6] .y =    3.80;             M[23] . Xp[ 6] .x=  0.398;          
 M[23] . Xp[ 7] .y =    4.40;             M[23] . Xp[ 7] .x=  0.380;          
 M[23] . Xp[ 8] .y =    5.00;             M[23] . Xp[ 8] .x=  0.359;          
 M[23] . Xp[ 9] .y =    5.60;             M[23] . Xp[ 9] .x=  0.334;          
 M[23] . Xp[10] .y =    6.20;             M[23] . Xp[10] .x=  0.277;          
 M[23] . Xp[11] .y =    6.80;             M[23] . Xp[11] .x=  0.239;          
 M[23] . Xp[12] .y =    7.40;             M[23] . Xp[12] .x=  0.211;          
 M[23] . Xp[13] .y =    8.00;             M[23] . Xp[13] .x=  0.185;          
 M[23] . Xp[14] .y =    8.60;             M[23] . Xp[14] .x=  0.159;          
 M[23] . Xp[15] .y =    9.20;             M[23] . Xp[15] .x=  0.142;          
 M[23] . Xp[16] .y =    9.80;             M[23] . Xp[16] .x=  0.129;          
 M[23] . Xp[17] .y =   10.40;             M[23] . Xp[17] .x=  0.124;          
 M[23] . Xp[18] .y =   11.00;             M[23] . Xp[18] .x=  0.099;          
 M[23] . Xp[19] .y =   11.60;             M[23] . Xp[19] .x=  0.092;          
 M[23] . Xp[20] .y =   12.20;             M[23] . Xp[20] .x=  0.082;          
 M[23] . Xp[21] .y =   12.80;             M[23] . Xp[21] .x=  0.076;          
 M[23] . Xp[22] .y =   13.40;             M[23] . Xp[22] .x=  0.068;          
 M[23] . Xp[23] .y =   14.00;             M[23] . Xp[23] .x=  0.067;          
 M[23] . Xp[24] .y =   14.60;             M[23] . Xp[24] .x=  0.062;          
 M[23] . Xp[25] .y =   15.20;             M[23] . Xp[25] .x=  0.055;          
 M[23] . Xp[26] .y =   15.80;             M[23] . Xp[26] .x=  0.049;          
 M[23] . Xp[27] .y =   16.40;             M[23] . Xp[27] .x=  0.047;          
 M[23] . Xp[28] .y =   17.00;             M[23] . Xp[28] .x=  0.050;          
 M[23] . Xp[29] .y =   17.60;             M[23] . Xp[29] .x=  0.037;          
 M[23] . Xp[30] .y =   18.20;             M[23] . Xp[30] .x=  0.030;          
 M[23] . Xp[31] .y =   18.80;             M[23] . Xp[31] .x=  0.031;          
 M[23] . Xp[32] .y =   19.40;             M[23] . Xp[32] .x=  0.025;          
 M[23] . Xp[33] .y =   20.00;             M[23] . Xp[33] .x=  0.022;          
 M[23] . Xp[34] .y =   20.60;             M[23] . Xp[34] .x=  0.012;          
 M[23] . Xp[35] .y =   21.20;             M[23] . Xp[35] .x=  0.007;          
 M[23] . Xp[36] .y =   21.80;             M[23] . Xp[36] .x=  0.004;          
 M[23] . Xp[37] .y =   22.40;             M[23] . Xp[37] .x=  0.004;          
 M[23] . Xp[38] .y =   23.00;             M[23] . Xp[38] .x=  0.004;          
 M[23] . Xp[39] .y =   23.60;             M[23] . Xp[39] .x=  0.001;          
 M[23] . Xp[40] .y =   24.20;             M[23] . Xp[40] .x=  0.007;          
 M[23] . Xp[41] .y =   24.80;             M[23] . Xp[41] .x=  0.004;          
 M[23] . Xp[42] .y =   25.40;             M[23] . Xp[42] .x=  0.003;          
   M[23].pathlength();
//station 25
  x[24] = 44.5;
  M[24].allocate(43);
  M[24].settype(SPLINE2D_LINEAR);                                                                                                                                                   
 M[24] . Xp[ 0] .y =    0.20;             M[24] . Xp[ 0] .x=  0.514;          
 M[24] . Xp[ 1] .y =    0.80;             M[24] . Xp[ 1] .x=  0.513;          
 M[24] . Xp[ 2] .y =    1.40;             M[24] . Xp[ 2] .x=  0.493;          
 M[24] . Xp[ 3] .y =    2.00;             M[24] . Xp[ 3] .x=  0.483;          
 M[24] . Xp[ 4] .y =    2.60;             M[24] . Xp[ 4] .x=  0.449;          
 M[24] . Xp[ 5] .y =    3.20;             M[24] . Xp[ 5] .x=  0.398;          
 M[24] . Xp[ 6] .y =    3.80;             M[24] . Xp[ 6] .x=  0.387;          
 M[24] . Xp[ 7] .y =    4.40;             M[24] . Xp[ 7] .x=  0.357;          
 M[24] . Xp[ 8] .y =    5.00;             M[24] . Xp[ 8] .x=  0.353;          
 M[24] . Xp[ 9] .y =    5.60;             M[24] . Xp[ 9] .x=  0.323;          
 M[24] . Xp[10] .y =    6.20;             M[24] . Xp[10] .x=  0.282;          
 M[24] . Xp[11] .y =    6.80;             M[24] . Xp[11] .x=  0.217;          
 M[24] . Xp[12] .y =    7.40;             M[24] . Xp[12] .x=  0.205;          
 M[24] . Xp[13] .y =    8.00;             M[24] . Xp[13] .x=  0.180;          
 M[24] . Xp[14] .y =    8.60;             M[24] . Xp[14] .x=  0.147;          
 M[24] . Xp[15] .y =    9.20;             M[24] . Xp[15] .x=  0.132;          
 M[24] . Xp[16] .y =    9.80;             M[24] . Xp[16] .x=  0.124;          
 M[24] . Xp[17] .y =   10.40;             M[24] . Xp[17] .x=  0.116;          
 M[24] . Xp[18] .y =   11.00;             M[24] . Xp[18] .x=  0.101;          
 M[24] . Xp[19] .y =   11.60;             M[24] . Xp[19] .x=  0.097;          
 M[24] . Xp[20] .y =   12.20;             M[24] . Xp[20] .x=  0.080;          
 M[24] . Xp[21] .y =   12.80;             M[24] . Xp[21] .x=  0.074;          
 M[24] . Xp[22] .y =   13.40;             M[24] . Xp[22] .x=  0.069;          
 M[24] . Xp[23] .y =   14.00;             M[24] . Xp[23] .x=  0.058;          
 M[24] . Xp[24] .y =   14.60;             M[24] . Xp[24] .x=  0.052;          
 M[24] . Xp[25] .y =   15.20;             M[24] . Xp[25] .x=  0.058;          
 M[24] . Xp[26] .y =   15.80;             M[24] . Xp[26] .x=  0.046;          
 M[24] . Xp[27] .y =   16.40;             M[24] . Xp[27] .x=  0.048;          
 M[24] . Xp[28] .y =   17.00;             M[24] . Xp[28] .x=  0.041;          
 M[24] . Xp[29] .y =   17.60;             M[24] . Xp[29] .x=  0.037;          
 M[24] . Xp[30] .y =   18.20;             M[24] . Xp[30] .x=  0.035;          
 M[24] . Xp[31] .y =   18.80;             M[24] . Xp[31] .x=  0.025;          
 M[24] . Xp[32] .y =   19.40;             M[24] . Xp[32] .x=  0.021;          
 M[24] . Xp[33] .y =   20.00;             M[24] . Xp[33] .x=  0.020;          
 M[24] . Xp[34] .y =   20.60;             M[24] . Xp[34] .x=  0.010;          
 M[24] . Xp[35] .y =   21.20;             M[24] . Xp[35] .x=  0.006;          
 M[24] . Xp[36] .y =   21.80;             M[24] . Xp[36] .x=  0.004;          
 M[24] . Xp[37] .y =   22.40;             M[24] . Xp[37] .x=  0.005;          
 M[24] . Xp[38] .y =   23.00;             M[24] . Xp[38] .x=  0.003;          
 M[24] . Xp[39] .y =   23.60;             M[24] . Xp[39] .x=  0.002;          
 M[24] . Xp[40] .y =   24.20;             M[24] . Xp[40] .x=  0.004;          
 M[24] . Xp[41] .y =   24.80;             M[24] . Xp[41] .x=  0.002;          
 M[24] . Xp[42] .y =   25.40;             M[24] . Xp[42] .x=  0.000;          
  M[24].pathlength();

//station 26
 x[25] = 45.3;
  M[25].allocate(43);
  M[25].settype(SPLINE2D_LINEAR);                                                                             
                                                                          
 M[25] . Xp[ 0] .y =    0.20;             M[25] . Xp[ 0] .x=  0.481;          
 M[25] . Xp[ 1] .y =    0.80;             M[25] . Xp[ 1] .x=  0.481;          
 M[25] . Xp[ 2] .y =    1.40;             M[25] . Xp[ 2] .x=  0.481;          
 M[25] . Xp[ 3] .y =    2.00;             M[25] . Xp[ 3] .x=  0.473;          
 M[25] . Xp[ 4] .y =    2.60;             M[25] . Xp[ 4] .x=  0.423;          
 M[25] . Xp[ 5] .y =    3.20;             M[25] . Xp[ 5] .x=  0.390;          
 M[25] . Xp[ 6] .y =    3.80;             M[25] . Xp[ 6] .x=  0.389;          
 M[25] . Xp[ 7] .y =    4.40;             M[25] . Xp[ 7] .x=  0.341;          
 M[25] . Xp[ 8] .y =    5.00;             M[25] . Xp[ 8] .x=  0.320;          
 M[25] . Xp[ 9] .y =    5.60;             M[25] . Xp[ 9] .x=  0.295;          
 M[25] . Xp[10] .y =    6.20;             M[25] . Xp[10] .x=  0.260;          
 M[25] . Xp[11] .y =    6.80;             M[25] . Xp[11] .x=  0.210;          
 M[25] . Xp[12] .y =    7.40;             M[25] . Xp[12] .x=  0.197;          
 M[25] . Xp[13] .y =    8.00;             M[25] . Xp[13] .x=  0.165;          
 M[25] . Xp[14] .y =    8.60;             M[25] . Xp[14] .x=  0.141;          
 M[25] . Xp[15] .y =    9.20;             M[25] . Xp[15] .x=  0.123;          
 M[25] . Xp[16] .y =    9.80;             M[25] . Xp[16] .x=  0.112;          
 M[25] . Xp[17] .y =   10.40;             M[25] . Xp[17] .x=  0.107;          
 M[25] . Xp[18] .y =   11.00;             M[25] . Xp[18] .x=  0.091;          
 M[25] . Xp[19] .y =   11.60;             M[25] . Xp[19] .x=  0.080;          
 M[25] . Xp[20] .y =   12.20;             M[25] . Xp[20] .x=  0.080;          
 M[25] . Xp[21] .y =   12.80;             M[25] . Xp[21] .x=  0.063;          
 M[25] . Xp[22] .y =   13.40;             M[25] . Xp[22] .x=  0.061;          
 M[25] . Xp[23] .y =   14.00;             M[25] . Xp[23] .x=  0.063;          
 M[25] . Xp[24] .y =   14.60;             M[25] . Xp[24] .x=  0.052;          
 M[25] . Xp[25] .y =   15.20;             M[25] . Xp[25] .x=  0.047;          
 M[25] . Xp[26] .y =   15.80;             M[25] . Xp[26] .x=  0.049;          
 M[25] . Xp[27] .y =   16.40;             M[25] . Xp[27] .x=  0.046;          
 M[25] . Xp[28] .y =   17.00;             M[25] . Xp[28] .x=  0.040;          
 M[25] . Xp[29] .y =   17.60;             M[25] . Xp[29] .x=  0.037;          
 M[25] . Xp[30] .y =   18.20;             M[25] . Xp[30] .x=  0.031;          
 M[25] . Xp[31] .y =   18.80;             M[25] . Xp[31] .x=  0.024;          
 M[25] . Xp[32] .y =   19.40;             M[25] . Xp[32] .x=  0.021;          
 M[25] . Xp[33] .y =   20.00;             M[25] . Xp[33] .x=  0.019;          
 M[25] . Xp[34] .y =   20.60;             M[25] . Xp[34] .x=  0.012;          
 M[25] . Xp[35] .y =   21.20;             M[25] . Xp[35] .x=  0.007;          
 M[25] . Xp[36] .y =   21.80;             M[25] . Xp[36] .x=  0.004;          
 M[25] . Xp[37] .y =   22.40;             M[25] . Xp[37] .x=  0.003;          
 M[25] . Xp[38] .y =   23.00;             M[25] . Xp[38] .x=  0.003;          
 M[25] . Xp[39] .y =   23.60;             M[25] . Xp[39] .x=  0.003;          
 M[25] . Xp[40] .y =   24.20;             M[25] . Xp[40] .x=  0.002;          
 M[25] . Xp[41] .y =   24.80;             M[25] . Xp[41] .x=  0.001;          
 M[25] . Xp[42] .y =   25.40;             M[25] . Xp[42] .x=  0.002;          
  M[25].pathlength();
//station 27
 x[26] = 46.0;
  M[26].allocate(43);
  M[26].settype(SPLINE2D_LINEAR);                                                                             
                                                                         
 M[26] . Xp[ 0] .y =    0.20;             M[26] . Xp[ 0] .x=  0.520;          
 M[26] . Xp[ 1] .y =    0.80;             M[26] . Xp[ 1] .x=  0.493;          
 M[26] . Xp[ 2] .y =    1.40;             M[26] . Xp[ 2] .x=  0.474;          
 M[26] . Xp[ 3] .y =    2.00;             M[26] . Xp[ 3] .x=  0.435;          
 M[26] . Xp[ 4] .y =    2.60;             M[26] . Xp[ 4] .x=  0.408;          
 M[26] . Xp[ 5] .y =    3.20;             M[26] . Xp[ 5] .x=  0.391;          
 M[26] . Xp[ 6] .y =    3.80;             M[26] . Xp[ 6] .x=  0.382;          
 M[26] . Xp[ 7] .y =    4.40;             M[26] . Xp[ 7] .x=  0.343;          
 M[26] . Xp[ 8] .y =    5.00;             M[26] . Xp[ 8] .x=  0.320;          
 M[26] . Xp[ 9] .y =    5.60;             M[26] . Xp[ 9] .x=  0.299;          
 M[26] . Xp[10] .y =    6.20;             M[26] . Xp[10] .x=  0.260;          
 M[26] . Xp[11] .y =    6.80;             M[26] . Xp[11] .x=  0.202;          
 M[26] . Xp[12] .y =    7.40;             M[26] . Xp[12] .x=  0.181;          
 M[26] . Xp[13] .y =    8.00;             M[26] . Xp[13] .x=  0.164;          
 M[26] . Xp[14] .y =    8.60;             M[26] . Xp[14] .x=  0.133;          
 M[26] . Xp[15] .y =    9.20;             M[26] . Xp[15] .x=  0.114;          
 M[26] . Xp[16] .y =    9.80;             M[26] . Xp[16] .x=  0.125;          
 M[26] . Xp[17] .y =   10.40;             M[26] . Xp[17] .x=  0.112;          
 M[26] . Xp[18] .y =   11.00;             M[26] . Xp[18] .x=  0.090;          
 M[26] . Xp[19] .y =   11.60;             M[26] . Xp[19] .x=  0.082;          
 M[26] . Xp[20] .y =   12.20;             M[26] . Xp[20] .x=  0.079;          
 M[26] . Xp[21] .y =   12.80;             M[26] . Xp[21] .x=  0.074;          
 M[26] . Xp[22] .y =   13.40;             M[26] . Xp[22] .x=  0.064;          
 M[26] . Xp[23] .y =   14.00;             M[26] . Xp[23] .x=  0.057;          
 M[26] . Xp[24] .y =   14.60;             M[26] . Xp[24] .x=  0.054;          
 M[26] . Xp[25] .y =   15.20;             M[26] . Xp[25] .x=  0.043;          
 M[26] . Xp[26] .y =   15.80;             M[26] . Xp[26] .x=  0.048;          
 M[26] . Xp[27] .y =   16.40;             M[26] . Xp[27] .x=  0.042;          
 M[26] . Xp[28] .y =   17.00;             M[26] . Xp[28] .x=  0.040;          
 M[26] . Xp[29] .y =   17.60;             M[26] . Xp[29] .x=  0.040;          
 M[26] . Xp[30] .y =   18.20;             M[26] . Xp[30] .x=  0.035;          
 M[26] . Xp[31] .y =   18.80;             M[26] . Xp[31] .x=  0.025;          
 M[26] . Xp[32] .y =   19.40;             M[26] . Xp[32] .x=  0.021;          
 M[26] . Xp[33] .y =   20.00;             M[26] . Xp[33] .x=  0.017;          
 M[26] . Xp[34] .y =   20.60;             M[26] . Xp[34] .x=  0.011;          
 M[26] . Xp[35] .y =   21.20;             M[26] . Xp[35] .x=  0.003;          
 M[26] . Xp[36] .y =   21.80;             M[26] . Xp[36] .x=  0.005;          
 M[26] . Xp[37] .y =   22.40;             M[26] . Xp[37] .x=  0.004;          
 M[26] . Xp[38] .y =   23.00;             M[26] . Xp[38] .x=  0.004;          
 M[26] . Xp[39] .y =   23.60;             M[26] . Xp[39] .x=  0.003;          
 M[26] . Xp[40] .y =   24.20;             M[26] . Xp[40] .x=  0.004;          
 M[26] . Xp[41] .y =   24.80;             M[26] . Xp[41] .x=  0.002;          
 M[26] . Xp[42] .y =   25.40;             M[26] . Xp[42] .x=  0.008;          
  M[26].pathlength();
//station 28
 x[27] = 46.8;
  M[27].allocate(43);
  M[27].settype(SPLINE2D_LINEAR);                                                                             
                                                                         
 M[27] . Xp[ 0] .y =    0.20;             M[27] . Xp[ 0] .x=  0.525;          
 M[27] . Xp[ 1] .y =    0.80;             M[27] . Xp[ 1] .x=  0.488;          
 M[27] . Xp[ 2] .y =    1.40;             M[27] . Xp[ 2] .x=  0.476;          
 M[27] . Xp[ 3] .y =    2.00;             M[27] . Xp[ 3] .x=  0.447;          
 M[27] . Xp[ 4] .y =    2.60;             M[27] . Xp[ 4] .x=  0.412;          
 M[27] . Xp[ 5] .y =    3.20;             M[27] . Xp[ 5] .x=  0.384;          
 M[27] . Xp[ 6] .y =    3.80;             M[27] . Xp[ 6] .x=  0.357;          
 M[27] . Xp[ 7] .y =    4.40;             M[27] . Xp[ 7] .x=  0.343;          
 M[27] . Xp[ 8] .y =    5.00;             M[27] . Xp[ 8] .x=  0.338;          
 M[27] . Xp[ 9] .y =    5.60;             M[27] . Xp[ 9] .x=  0.275;          
 M[27] . Xp[10] .y =    6.20;             M[27] . Xp[10] .x=  0.244;          
 M[27] . Xp[11] .y =    6.80;             M[27] . Xp[11] .x=  0.192;          
 M[27] . Xp[12] .y =    7.40;             M[27] . Xp[12] .x=  0.157;          
 M[27] . Xp[13] .y =    8.00;             M[27] . Xp[13] .x=  0.153;          
 M[27] . Xp[14] .y =    8.60;             M[27] . Xp[14] .x=  0.149;          
 M[27] . Xp[15] .y =    9.20;             M[27] . Xp[15] .x=  0.129;          
 M[27] . Xp[16] .y =    9.80;             M[27] . Xp[16] .x=  0.122;          
 M[27] . Xp[17] .y =   10.40;             M[27] . Xp[17] .x=  0.100;          
 M[27] . Xp[18] .y =   11.00;             M[27] . Xp[18] .x=  0.081;          
 M[27] . Xp[19] .y =   11.60;             M[27] . Xp[19] .x=  0.074;          
 M[27] . Xp[20] .y =   12.20;             M[27] . Xp[20] .x=  0.077;          
 M[27] . Xp[21] .y =   12.80;             M[27] . Xp[21] .x=  0.075;          
 M[27] . Xp[22] .y =   13.40;             M[27] . Xp[22] .x=  0.064;          
 M[27] . Xp[23] .y =   14.00;             M[27] . Xp[23] .x=  0.063;          
 M[27] . Xp[24] .y =   14.60;             M[27] . Xp[24] .x=  0.052;          
 M[27] . Xp[25] .y =   15.20;             M[27] . Xp[25] .x=  0.041;          
 M[27] . Xp[26] .y =   15.80;             M[27] . Xp[26] .x=  0.042;          
 M[27] . Xp[27] .y =   16.40;             M[27] . Xp[27] .x=  0.039;          
 M[27] . Xp[28] .y =   17.00;             M[27] . Xp[28] .x=  0.037;          
 M[27] . Xp[29] .y =   17.60;             M[27] . Xp[29] .x=  0.034;          
 M[27] . Xp[30] .y =   18.20;             M[27] . Xp[30] .x=  0.027;          
 M[27] . Xp[31] .y =   18.80;             M[27] . Xp[31] .x=  0.024;          
 M[27] . Xp[32] .y =   19.40;             M[27] . Xp[32] .x=  0.016;          
 M[27] . Xp[33] .y =   20.00;             M[27] . Xp[33] .x=  0.016;          
 M[27] . Xp[34] .y =   20.60;             M[27] . Xp[34] .x=  0.012;          
 M[27] . Xp[35] .y =   21.20;             M[27] . Xp[35] .x=  0.007;          
 M[27] . Xp[36] .y =   21.80;             M[27] . Xp[36] .x=  0.004;          
 M[27] . Xp[37] .y =   22.40;             M[27] . Xp[37] .x=  0.005;          
 M[27] . Xp[38] .y =   23.00;             M[27] . Xp[38] .x=  0.004;          
 M[27] . Xp[39] .y =   23.60;             M[27] . Xp[39] .x=  0.002;          
 M[27] . Xp[40] .y =   24.20;             M[27] . Xp[40] .x=  0.003;          
 M[27] . Xp[41] .y =   24.80;             M[27] . Xp[41] .x=  0.003;          
 M[27] . Xp[42] .y =   25.40;             M[27] . Xp[42] .x=  0.002;          
  M[27].pathlength();
//station 29
 x[28] = 47.5;
  M[28].allocate(43);
  M[28].settype(SPLINE2D_LINEAR);                                                                             
                                                                        
 M[28] . Xp[ 0] .y =    0.20;             M[28] . Xp[ 0] .x=  0.480;          
 M[28] . Xp[ 1] .y =    0.80;             M[28] . Xp[ 1] .x=  0.455;          
 M[28] . Xp[ 2] .y =    1.40;             M[28] . Xp[ 2] .x=  0.442;          
 M[28] . Xp[ 3] .y =    2.00;             M[28] . Xp[ 3] .x=  0.444;          
 M[28] . Xp[ 4] .y =    2.60;             M[28] . Xp[ 4] .x=  0.417;          
 M[28] . Xp[ 5] .y =    3.20;             M[28] . Xp[ 5] .x=  0.403;          
 M[28] . Xp[ 6] .y =    3.80;             M[28] . Xp[ 6] .x=  0.381;          
 M[28] . Xp[ 7] .y =    4.40;             M[28] . Xp[ 7] .x=  0.350;          
 M[28] . Xp[ 8] .y =    5.00;             M[28] . Xp[ 8] .x=  0.308;          
 M[28] . Xp[ 9] .y =    5.60;             M[28] . Xp[ 9] .x=  0.280;          
 M[28] . Xp[10] .y =    6.20;             M[28] . Xp[10] .x=  0.235;          
 M[28] . Xp[11] .y =    6.80;             M[28] . Xp[11] .x=  0.205;          
 M[28] . Xp[12] .y =    7.40;             M[28] . Xp[12] .x=  0.173;          
 M[28] . Xp[13] .y =    8.00;             M[28] . Xp[13] .x=  0.160;          
 M[28] . Xp[14] .y =    8.60;             M[28] . Xp[14] .x=  0.146;          
 M[28] . Xp[15] .y =    9.20;             M[28] . Xp[15] .x=  0.121;          
 M[28] . Xp[16] .y =    9.80;             M[28] . Xp[16] .x=  0.109;          
 M[28] . Xp[17] .y =   10.40;             M[28] . Xp[17] .x=  0.096;          
 M[28] . Xp[18] .y =   11.00;             M[28] . Xp[18] .x=  0.080;          
 M[28] . Xp[19] .y =   11.60;             M[28] . Xp[19] .x=  0.070;          
 M[28] . Xp[20] .y =   12.20;             M[28] . Xp[20] .x=  0.063;          
 M[28] . Xp[21] .y =   12.80;             M[28] . Xp[21] .x=  0.068;          
 M[28] . Xp[22] .y =   13.40;             M[28] . Xp[22] .x=  0.065;          
 M[28] . Xp[23] .y =   14.00;             M[28] . Xp[23] .x=  0.054;          
 M[28] . Xp[24] .y =   14.60;             M[28] . Xp[24] .x=  0.050;          
 M[28] . Xp[25] .y =   15.20;             M[28] . Xp[25] .x=  0.043;          
 M[28] . Xp[26] .y =   15.80;             M[28] . Xp[26] .x=  0.040;          
 M[28] . Xp[27] .y =   16.40;             M[28] . Xp[27] .x=  0.028;          
 M[28] . Xp[28] .y =   17.00;             M[28] . Xp[28] .x=  0.031;          
 M[28] . Xp[29] .y =   17.60;             M[28] . Xp[29] .x=  0.030;          
 M[28] . Xp[30] .y =   18.20;             M[28] . Xp[30] .x=  0.026;          
 M[28] . Xp[31] .y =   18.80;             M[28] . Xp[31] .x=  0.019;          
 M[28] . Xp[32] .y =   19.40;             M[28] . Xp[32] .x=  0.014;          
 M[28] . Xp[33] .y =   20.00;             M[28] . Xp[33] .x=  0.010;          
 M[28] . Xp[34] .y =   20.60;             M[28] . Xp[34] .x=  0.005;          
 M[28] . Xp[35] .y =   21.20;             M[28] . Xp[35] .x=  0.003;          
 M[28] . Xp[36] .y =   21.80;             M[28] . Xp[36] .x=  0.004;          
 M[28] . Xp[37] .y =   22.40;             M[28] . Xp[37] .x=  0.002;          
 M[28] . Xp[38] .y =   23.00;             M[28] . Xp[38] .x=  0.004;          
 M[28] . Xp[39] .y =   23.60;             M[28] . Xp[39] .x=  0.002;          
 M[28] . Xp[40] .y =   24.20;             M[28] . Xp[40] .x=  0.001;          
 M[28] . Xp[41] .y =   24.80;             M[28] . Xp[41] .x=  0.002;          
 M[28] . Xp[42] .y =   25.40;             M[28] . Xp[42] .x=  0.002;          
   M[28].pathlength();
//station 30
 x[29] = 48.3;
  M[29].allocate(43);
  M[29].settype(SPLINE2D_LINEAR);                                                                                                                                                     
 M[29] . Xp[ 0] .y =    0.20;             M[29] . Xp[ 0] .x=  0.435;          
 M[29] . Xp[ 1] .y =    0.80;             M[29] . Xp[ 1] .x=  0.457;          
 M[29] . Xp[ 2] .y =    1.40;             M[29] . Xp[ 2] .x=  0.421;          
 M[29] . Xp[ 3] .y =    2.00;             M[29] . Xp[ 3] .x=  0.414;          
 M[29] . Xp[ 4] .y =    2.60;             M[29] . Xp[ 4] .x=  0.393;          
 M[29] . Xp[ 5] .y =    3.20;             M[29] . Xp[ 5] .x=  0.390;          
 M[29] . Xp[ 6] .y =    3.80;             M[29] . Xp[ 6] .x=  0.380;          
 M[29] . Xp[ 7] .y =    4.40;             M[29] . Xp[ 7] .x=  0.323;          
 M[29] . Xp[ 8] .y =    5.00;             M[29] . Xp[ 8] .x=  0.314;          
 M[29] . Xp[ 9] .y =    5.60;             M[29] . Xp[ 9] .x=  0.262;          
 M[29] . Xp[10] .y =    6.20;             M[29] . Xp[10] .x=  0.251;          
 M[29] . Xp[11] .y =    6.80;             M[29] . Xp[11] .x=  0.226;          
 M[29] . Xp[12] .y =    7.40;             M[29] . Xp[12] .x=  0.192;          
 M[29] . Xp[13] .y =    8.00;             M[29] . Xp[13] .x=  0.171;          
 M[29] . Xp[14] .y =    8.60;             M[29] . Xp[14] .x=  0.160;          
 M[29] . Xp[15] .y =    9.20;             M[29] . Xp[15] .x=  0.120;          
 M[29] . Xp[16] .y =    9.80;             M[29] . Xp[16] .x=  0.111;          
 M[29] . Xp[17] .y =   10.40;             M[29] . Xp[17] .x=  0.103;          
 M[29] . Xp[18] .y =   11.00;             M[29] . Xp[18] .x=  0.089;          
 M[29] . Xp[19] .y =   11.60;             M[29] . Xp[19] .x=  0.070;          
 M[29] . Xp[20] .y =   12.20;             M[29] . Xp[20] .x=  0.066;          
 M[29] . Xp[21] .y =   12.80;             M[29] . Xp[21] .x=  0.057;          
 M[29] . Xp[22] .y =   13.40;             M[29] . Xp[22] .x=  0.060;          
 M[29] . Xp[23] .y =   14.00;             M[29] . Xp[23] .x=  0.057;          
 M[29] . Xp[24] .y =   14.60;             M[29] . Xp[24] .x=  0.046;          
 M[29] . Xp[25] .y =   15.20;             M[29] . Xp[25] .x=  0.044;          
 M[29] . Xp[26] .y =   15.80;             M[29] . Xp[26] .x=  0.040;          
 M[29] . Xp[27] .y =   16.40;             M[29] . Xp[27] .x=  0.028;          
 M[29] . Xp[28] .y =   17.00;             M[29] . Xp[28] .x=  0.031;          
 M[29] . Xp[29] .y =   17.60;             M[29] . Xp[29] .x=  0.032;          
 M[29] . Xp[30] .y =   18.20;             M[29] . Xp[30] .x=  0.020;          
 M[29] . Xp[31] .y =   18.80;             M[29] . Xp[31] .x=  0.021;          
 M[29] . Xp[32] .y =   19.40;             M[29] . Xp[32] .x=  0.013;          
 M[29] . Xp[33] .y =   20.00;             M[29] . Xp[33] .x=  0.006;          
 M[29] . Xp[34] .y =   20.60;             M[29] . Xp[34] .x=  0.008;          
 M[29] . Xp[35] .y =   21.20;             M[29] . Xp[35] .x=  0.006;          
 M[29] . Xp[36] .y =   21.80;             M[29] . Xp[36] .x=  0.003;          
 M[29] . Xp[37] .y =   22.40;             M[29] . Xp[37] .x=  0.002;          
 M[29] . Xp[38] .y =   23.00;             M[29] . Xp[38] .x=  0.003;          
 M[29] . Xp[39] .y =   23.60;             M[29] . Xp[39] .x=  0.001;          
 M[29] . Xp[40] .y =   24.20;             M[29] . Xp[40] .x=  0.002;          
 M[29] . Xp[41] .y =   24.80;             M[29] . Xp[41] .x=  0.004;          
 M[29] . Xp[42] .y =   25.40;             M[29] . Xp[42] .x=  0.003;          
  M[29].pathlength();

//station 31
 x[30] = 49.0;
  M[30].allocate(43);
  M[30].settype(SPLINE2D_LINEAR);                                                                             
                                                                         
 M[30] . Xp[ 0] .y =    0.20;             M[30] . Xp[ 0] .x=  0.444;          
 M[30] . Xp[ 1] .y =    0.80;             M[30] . Xp[ 1] .x=  0.449;          
 M[30] . Xp[ 2] .y =    1.40;             M[30] . Xp[ 2] .x=  0.413;          
 M[30] . Xp[ 3] .y =    2.00;             M[30] . Xp[ 3] .x=  0.391;          
 M[30] . Xp[ 4] .y =    2.60;             M[30] . Xp[ 4] .x=  0.409;          
 M[30] . Xp[ 5] .y =    3.20;             M[30] . Xp[ 5] .x=  0.376;          
 M[30] . Xp[ 6] .y =    3.80;             M[30] . Xp[ 6] .x=  0.349;          
 M[30] . Xp[ 7] .y =    4.40;             M[30] . Xp[ 7] .x=  0.333;          
 M[30] . Xp[ 8] .y =    5.00;             M[30] . Xp[ 8] .x=  0.313;          
 M[30] . Xp[ 9] .y =    5.60;             M[30] . Xp[ 9] .x=  0.289;          
 M[30] . Xp[10] .y =    6.20;             M[30] . Xp[10] .x=  0.247;          
 M[30] . Xp[11] .y =    6.80;             M[30] . Xp[11] .x=  0.221;          
 M[30] . Xp[12] .y =    7.40;             M[30] . Xp[12] .x=  0.202;          
 M[30] . Xp[13] .y =    8.00;             M[30] . Xp[13] .x=  0.173;          
 M[30] . Xp[14] .y =    8.60;             M[30] . Xp[14] .x=  0.140;          
 M[30] . Xp[15] .y =    9.20;             M[30] . Xp[15] .x=  0.128;          
 M[30] . Xp[16] .y =    9.80;             M[30] . Xp[16] .x=  0.101;          
 M[30] . Xp[17] .y =   10.40;             M[30] . Xp[17] .x=  0.098;          
 M[30] . Xp[18] .y =   11.00;             M[30] . Xp[18] .x=  0.075;          
 M[30] . Xp[19] .y =   11.60;             M[30] . Xp[19] .x=  0.070;          
 M[30] . Xp[20] .y =   12.20;             M[30] . Xp[20] .x=  0.065;          
 M[30] . Xp[21] .y =   12.80;             M[30] . Xp[21] .x=  0.052;          
 M[30] . Xp[22] .y =   13.40;             M[30] . Xp[22] .x=  0.062;          
 M[30] . Xp[23] .y =   14.00;             M[30] . Xp[23] .x=  0.053;          
 M[30] . Xp[24] .y =   14.60;             M[30] . Xp[24] .x=  0.050;          
 M[30] . Xp[25] .y =   15.20;             M[30] . Xp[25] .x=  0.043;          
 M[30] . Xp[26] .y =   15.80;             M[30] . Xp[26] .x=  0.034;          
 M[30] . Xp[27] .y =   16.40;             M[30] . Xp[27] .x=  0.029;          
 M[30] . Xp[28] .y =   17.00;             M[30] . Xp[28] .x=  0.031;          
 M[30] . Xp[29] .y =   17.60;             M[30] . Xp[29] .x=  0.027;          
 M[30] . Xp[30] .y =   18.20;             M[30] . Xp[30] .x=  0.026;          
 M[30] . Xp[31] .y =   18.80;             M[30] . Xp[31] .x=  0.020;          
 M[30] . Xp[32] .y =   19.40;             M[30] . Xp[32] .x=  0.017;          
 M[30] . Xp[33] .y =   20.00;             M[30] . Xp[33] .x=  0.010;          
 M[30] . Xp[34] .y =   20.60;             M[30] . Xp[34] .x=  0.005;          
 M[30] . Xp[35] .y =   21.20;             M[30] . Xp[35] .x=  0.003;          
 M[30] . Xp[36] .y =   21.80;             M[30] . Xp[36] .x=  0.002;          
 M[30] . Xp[37] .y =   22.40;             M[30] . Xp[37] .x=  0.004;          
 M[30] . Xp[38] .y =   23.00;             M[30] . Xp[38] .x=  0.001;          
 M[30] . Xp[39] .y =   23.60;             M[30] . Xp[39] .x=  0.001;          
 M[30] . Xp[40] .y =   24.20;             M[30] . Xp[40] .x=  0.002;          
 M[30] . Xp[41] .y =   24.80;             M[30] . Xp[41] .x=  0.004;          
 M[30] . Xp[42] .y =   25.40;             M[30] . Xp[42] .x=  0.004;          
   M[30].pathlength();
//station 32
 x[31] = 49.8;
  M[31].allocate(43);
  M[31].settype(SPLINE2D_LINEAR);                                                                            
                                                                         
 M[31] . Xp[ 0] .y =    0.20;             M[31] . Xp[ 0] .x=  0.446;          
 M[31] . Xp[ 1] .y =    0.80;             M[31] . Xp[ 1] .x=  0.434;          
 M[31] . Xp[ 2] .y =    1.40;             M[31] . Xp[ 2] .x=  0.387;          
 M[31] . Xp[ 3] .y =    2.00;             M[31] . Xp[ 3] .x=  0.386;          
 M[31] . Xp[ 4] .y =    2.60;             M[31] . Xp[ 4] .x=  0.373;          
 M[31] . Xp[ 5] .y =    3.20;             M[31] . Xp[ 5] .x=  0.345;          
 M[31] . Xp[ 6] .y =    3.80;             M[31] . Xp[ 6] .x=  0.352;          
 M[31] . Xp[ 7] .y =    4.40;             M[31] . Xp[ 7] .x=  0.324;          
 M[31] . Xp[ 8] .y =    5.00;             M[31] . Xp[ 8] .x=  0.316;          
 M[31] . Xp[ 9] .y =    5.60;             M[31] . Xp[ 9] .x=  0.277;          
 M[31] . Xp[10] .y =    6.20;             M[31] . Xp[10] .x=  0.241;          
 M[31] . Xp[11] .y =    6.80;             M[31] . Xp[11] .x=  0.218;          
 M[31] . Xp[12] .y =    7.40;             M[31] . Xp[12] .x=  0.197;          
 M[31] . Xp[13] .y =    8.00;             M[31] . Xp[13] .x=  0.180;          
 M[31] . Xp[14] .y =    8.60;             M[31] . Xp[14] .x=  0.135;          
 M[31] . Xp[15] .y =    9.20;             M[31] . Xp[15] .x=  0.110;          
 M[31] . Xp[16] .y =    9.80;             M[31] . Xp[16] .x=  0.095;          
 M[31] . Xp[17] .y =   10.40;             M[31] . Xp[17] .x=  0.093;          
 M[31] . Xp[18] .y =   11.00;             M[31] . Xp[18] .x=  0.076;          
 M[31] . Xp[19] .y =   11.60;             M[31] . Xp[19] .x=  0.068;          
 M[31] . Xp[20] .y =   12.20;             M[31] . Xp[20] .x=  0.065;          
 M[31] . Xp[21] .y =   12.80;             M[31] . Xp[21] .x=  0.061;          
 M[31] . Xp[22] .y =   13.40;             M[31] . Xp[22] .x=  0.057;          
 M[31] . Xp[23] .y =   14.00;             M[31] . Xp[23] .x=  0.055;          
 M[31] . Xp[24] .y =   14.60;             M[31] . Xp[24] .x=  0.047;          
 M[31] . Xp[25] .y =   15.20;             M[31] . Xp[25] .x=  0.039;          
 M[31] . Xp[26] .y =   15.80;             M[31] . Xp[26] .x=  0.036;          
 M[31] . Xp[27] .y =   16.40;             M[31] . Xp[27] .x=  0.029;          
 M[31] . Xp[28] .y =   17.00;             M[31] . Xp[28] .x=  0.023;          
 M[31] . Xp[29] .y =   17.60;             M[31] . Xp[29] .x=  0.022;          
 M[31] . Xp[30] .y =   18.20;             M[31] . Xp[30] .x=  0.022;          
 M[31] . Xp[31] .y =   18.80;             M[31] . Xp[31] .x=  0.015;          
 M[31] . Xp[32] .y =   19.40;             M[31] . Xp[32] .x=  0.013;          
 M[31] . Xp[33] .y =   20.00;             M[31] . Xp[33] .x=  0.009;          
 M[31] . Xp[34] .y =   20.60;             M[31] . Xp[34] .x=  0.005;          
 M[31] . Xp[35] .y =   21.20;             M[31] . Xp[35] .x=  0.005;          
 M[31] . Xp[36] .y =   21.80;             M[31] . Xp[36] .x=  0.003;          
 M[31] . Xp[37] .y =   22.40;             M[31] . Xp[37] .x=  0.003;          
 M[31] . Xp[38] .y =   23.00;             M[31] . Xp[38] .x=  0.002;          
 M[31] . Xp[39] .y =   23.60;             M[31] . Xp[39] .x=  0.001;          
 M[31] . Xp[40] .y =   24.20;             M[31] . Xp[40] .x=  0.000;          
 M[31] . Xp[41] .y =   24.80;             M[31] . Xp[41] .x=  0.002;          
 M[31] . Xp[42] .y =   25.40;             M[31] . Xp[42] .x=  0.001;          
  M[31].pathlength();
//station 33
 x[32] = 50.5;
  M[32].allocate(43);
  M[32].settype(SPLINE2D_LINEAR);                                                                             
                                                                         
 M[32] . Xp[ 0] .y =    0.20;             M[32] . Xp[ 0] .x=  0.447;          
 M[32] . Xp[ 1] .y =    0.80;             M[32] . Xp[ 1] .x=  0.420;          
 M[32] . Xp[ 2] .y =    1.40;             M[32] . Xp[ 2] .x=  0.414;          
 M[32] . Xp[ 3] .y =    2.00;             M[32] . Xp[ 3] .x=  0.413;          
 M[32] . Xp[ 4] .y =    2.60;             M[32] . Xp[ 4] .x=  0.390;          
 M[32] . Xp[ 5] .y =    3.20;             M[32] . Xp[ 5] .x=  0.363;          
 M[32] . Xp[ 6] .y =    3.80;             M[32] . Xp[ 6] .x=  0.343;          
 M[32] . Xp[ 7] .y =    4.40;             M[32] . Xp[ 7] .x=  0.334;          
 M[32] . Xp[ 8] .y =    5.00;             M[32] . Xp[ 8] .x=  0.310;          
 M[32] . Xp[ 9] .y =    5.60;             M[32] . Xp[ 9] .x=  0.282;          
 M[32] . Xp[10] .y =    6.20;             M[32] . Xp[10] .x=  0.248;          
 M[32] . Xp[11] .y =    6.80;             M[32] . Xp[11] .x=  0.224;          
 M[32] . Xp[12] .y =    7.40;             M[32] . Xp[12] .x=  0.201;          
 M[32] . Xp[13] .y =    8.00;             M[32] . Xp[13] .x=  0.173;          
 M[32] . Xp[14] .y =    8.60;             M[32] . Xp[14] .x=  0.135;          
 M[32] . Xp[15] .y =    9.20;             M[32] . Xp[15] .x=  0.110;          
 M[32] . Xp[16] .y =    9.80;             M[32] . Xp[16] .x=  0.100;          
 M[32] . Xp[17] .y =   10.40;             M[32] . Xp[17] .x=  0.083;          
 M[32] . Xp[18] .y =   11.00;             M[32] . Xp[18] .x=  0.081;          
 M[32] . Xp[19] .y =   11.60;             M[32] . Xp[19] .x=  0.077;          
 M[32] . Xp[20] .y =   12.20;             M[32] . Xp[20] .x=  0.062;          
 M[32] . Xp[21] .y =   12.80;             M[32] . Xp[21] .x=  0.058;          
 M[32] . Xp[22] .y =   13.40;             M[32] . Xp[22] .x=  0.064;          
 M[32] . Xp[23] .y =   14.00;             M[32] . Xp[23] .x=  0.053;          
 M[32] . Xp[24] .y =   14.60;             M[32] . Xp[24] .x=  0.048;          
 M[32] . Xp[25] .y =   15.20;             M[32] . Xp[25] .x=  0.036;          
 M[32] . Xp[26] .y =   15.80;             M[32] . Xp[26] .x=  0.039;          
 M[32] . Xp[27] .y =   16.40;             M[32] . Xp[27] .x=  0.032;          
 M[32] . Xp[28] .y =   17.00;             M[32] . Xp[28] .x=  0.023;          
 M[32] . Xp[29] .y =   17.60;             M[32] . Xp[29] .x=  0.029;          
 M[32] . Xp[30] .y =   18.20;             M[32] . Xp[30] .x=  0.025;          
 M[32] . Xp[31] .y =   18.80;             M[32] . Xp[31] .x=  0.014;          
 M[32] . Xp[32] .y =   19.40;             M[32] . Xp[32] .x=  0.013;          
 M[32] . Xp[33] .y =   20.00;             M[32] . Xp[33] .x=  0.008;          
 M[32] . Xp[34] .y =   20.60;             M[32] . Xp[34] .x=  0.006;          
 M[32] . Xp[35] .y =   21.20;             M[32] . Xp[35] .x=  0.003;          
 M[32] . Xp[36] .y =   21.80;             M[32] . Xp[36] .x=  0.002;          
 M[32] . Xp[37] .y =   22.40;             M[32] . Xp[37] .x=  0.002;          
 M[32] . Xp[38] .y =   23.00;             M[32] . Xp[38] .x=  0.004;          
 M[32] . Xp[39] .y =   23.60;             M[32] . Xp[39] .x=  0.006;          
 M[32] . Xp[40] .y =   24.20;             M[32] . Xp[40] .x=  0.002;          
 M[32] . Xp[41] .y =   24.80;             M[32] . Xp[41] .x=  0.004;          
 M[32] . Xp[42] .y =   25.40;             M[32] . Xp[42] .x=  0.002;          
  M[32].pathlength();
//station 34
 x[33] = 51.3;
  M[33].allocate(43);
  M[33].settype(SPLINE2D_LINEAR);                                                                             
                                                                         
 M[33] . Xp[ 0] .y =    0.20;             M[33] . Xp[ 0] .x=  0.404;          
 M[33] . Xp[ 1] .y =    0.80;             M[33] . Xp[ 1] .x=  0.419;          
 M[33] . Xp[ 2] .y =    1.40;             M[33] . Xp[ 2] .x=  0.430;          
 M[33] . Xp[ 3] .y =    2.00;             M[33] . Xp[ 3] .x=  0.393;          
 M[33] . Xp[ 4] .y =    2.60;             M[33] . Xp[ 4] .x=  0.364;          
 M[33] . Xp[ 5] .y =    3.20;             M[33] . Xp[ 5] .x=  0.343;          
 M[33] . Xp[ 6] .y =    3.80;             M[33] . Xp[ 6] .x=  0.318;          
 M[33] . Xp[ 7] .y =    4.40;             M[33] . Xp[ 7] .x=  0.308;          
 M[33] . Xp[ 8] .y =    5.00;             M[33] . Xp[ 8] .x=  0.297;          
 M[33] . Xp[ 9] .y =    5.60;             M[33] . Xp[ 9] .x=  0.266;          
 M[33] . Xp[10] .y =    6.20;             M[33] . Xp[10] .x=  0.226;          
 M[33] . Xp[11] .y =    6.80;             M[33] . Xp[11] .x=  0.210;          
 M[33] . Xp[12] .y =    7.40;             M[33] . Xp[12] .x=  0.173;          
 M[33] . Xp[13] .y =    8.00;             M[33] . Xp[13] .x=  0.155;          
 M[33] . Xp[14] .y =    8.60;             M[33] . Xp[14] .x=  0.123;          
 M[33] . Xp[15] .y =    9.20;             M[33] . Xp[15] .x=  0.114;          
 M[33] . Xp[16] .y =    9.80;             M[33] . Xp[16] .x=  0.105;          
 M[33] . Xp[17] .y =   10.40;             M[33] . Xp[17] .x=  0.075;          
 M[33] . Xp[18] .y =   11.00;             M[33] . Xp[18] .x=  0.079;          
 M[33] . Xp[19] .y =   11.60;             M[33] . Xp[19] .x=  0.064;          
 M[33] . Xp[20] .y =   12.20;             M[33] . Xp[20] .x=  0.057;          
 M[33] . Xp[21] .y =   12.80;             M[33] . Xp[21] .x=  0.051;          
 M[33] . Xp[22] .y =   13.40;             M[33] . Xp[22] .x=  0.059;          
 M[33] . Xp[23] .y =   14.00;             M[33] . Xp[23] .x=  0.048;          
 M[33] . Xp[24] .y =   14.60;             M[33] . Xp[24] .x=  0.048;          
 M[33] . Xp[25] .y =   15.20;             M[33] . Xp[25] .x=  0.038;          
 M[33] . Xp[26] .y =   15.80;             M[33] . Xp[26] .x=  0.038;          
 M[33] . Xp[27] .y =   16.40;             M[33] . Xp[27] .x=  0.032;          
 M[33] . Xp[28] .y =   17.00;             M[33] . Xp[28] .x=  0.023;          
 M[33] . Xp[29] .y =   17.60;             M[33] . Xp[29] .x=  0.030;          
 M[33] . Xp[30] .y =   18.20;             M[33] . Xp[30] .x=  0.025;          
 M[33] . Xp[31] .y =   18.80;             M[33] . Xp[31] .x=  0.023;          
 M[33] . Xp[32] .y =   19.40;             M[33] . Xp[32] .x=  0.013;          
 M[33] . Xp[33] .y =   20.00;             M[33] . Xp[33] .x=  0.009;          
 M[33] . Xp[34] .y =   20.60;             M[33] . Xp[34] .x=  0.009;          
 M[33] . Xp[35] .y =   21.20;             M[33] . Xp[35] .x=  0.004;          
 M[33] . Xp[36] .y =   21.80;             M[33] . Xp[36] .x=  0.004;          
 M[33] . Xp[37] .y =   22.40;             M[33] . Xp[37] .x=  0.004;          
 M[33] . Xp[38] .y =   23.00;             M[33] . Xp[38] .x=  0.004;          
 M[33] . Xp[39] .y =   23.60;             M[33] . Xp[39] .x=  0.002;          
 M[33] . Xp[40] .y =   24.20;             M[33] . Xp[40] .x=  0.009;          
 M[33] . Xp[41] .y =   24.80;             M[33] . Xp[41] .x=  0.003;          
 M[33] . Xp[42] .y =   25.40;             M[33] . Xp[42] .x=  0.007;          
   M[33].pathlength();
//station 35
 x[34] = 52.0;
  M[34].allocate(43);
  M[34].settype(SPLINE2D_LINEAR);                                                                            
                                                                         
 M[34] . Xp[ 0] .y =    0.20;             M[34] . Xp[ 0] .x=  0.397;          
 M[34] . Xp[ 1] .y =    0.80;             M[34] . Xp[ 1] .x=  0.413;          
 M[34] . Xp[ 2] .y =    1.40;             M[34] . Xp[ 2] .x=  0.393;          
 M[34] . Xp[ 3] .y =    2.00;             M[34] . Xp[ 3] .x=  0.393;          
 M[34] . Xp[ 4] .y =    2.60;             M[34] . Xp[ 4] .x=  0.359;          
 M[34] . Xp[ 5] .y =    3.20;             M[34] . Xp[ 5] .x=  0.323;          
 M[34] . Xp[ 6] .y =    3.80;             M[34] . Xp[ 6] .x=  0.340;          
 M[34] . Xp[ 7] .y =    4.40;             M[34] . Xp[ 7] .x=  0.292;          
 M[34] . Xp[ 8] .y =    5.00;             M[34] . Xp[ 8] .x=  0.283;          
 M[34] . Xp[ 9] .y =    5.60;             M[34] . Xp[ 9] .x=  0.241;          
 M[34] . Xp[10] .y =    6.20;             M[34] . Xp[10] .x=  0.216;          
 M[34] . Xp[11] .y =    6.80;             M[34] . Xp[11] .x=  0.194;          
 M[34] . Xp[12] .y =    7.40;             M[34] . Xp[12] .x=  0.151;          
 M[34] . Xp[13] .y =    8.00;             M[34] . Xp[13] .x=  0.117;          
 M[34] . Xp[14] .y =    8.60;             M[34] . Xp[14] .x=  0.116;          
 M[34] . Xp[15] .y =    9.20;             M[34] . Xp[15] .x=  0.106;          
 M[34] . Xp[16] .y =    9.80;             M[34] . Xp[16] .x=  0.096;          
 M[34] . Xp[17] .y =   10.40;             M[34] . Xp[17] .x=  0.073;          
 M[34] . Xp[18] .y =   11.00;             M[34] . Xp[18] .x=  0.066;          
 M[34] . Xp[19] .y =   11.60;             M[34] . Xp[19] .x=  0.064;          
 M[34] . Xp[20] .y =   12.20;             M[34] . Xp[20] .x=  0.050;          
 M[34] . Xp[21] .y =   12.80;             M[34] . Xp[21] .x=  0.038;          
 M[34] . Xp[22] .y =   13.40;             M[34] . Xp[22] .x=  0.050;          
 M[34] . Xp[23] .y =   14.00;             M[34] . Xp[23] .x=  0.039;          
 M[34] . Xp[24] .y =   14.60;             M[34] . Xp[24] .x=  0.030;          
 M[34] . Xp[25] .y =   15.20;             M[34] . Xp[25] .x=  0.035;          
 M[34] . Xp[26] .y =   15.80;             M[34] . Xp[26] .x=  0.029;          
 M[34] . Xp[27] .y =   16.40;             M[34] . Xp[27] .x=  0.029;          
 M[34] . Xp[28] .y =   17.00;             M[34] . Xp[28] .x=  0.025;          
 M[34] . Xp[29] .y =   17.60;             M[34] . Xp[29] .x=  0.022;          
 M[34] . Xp[30] .y =   18.20;             M[34] . Xp[30] .x=  0.014;          
 M[34] . Xp[31] .y =   18.80;             M[34] . Xp[31] .x=  0.015;          
 M[34] . Xp[32] .y =   19.40;             M[34] . Xp[32] .x=  0.014;          
 M[34] . Xp[33] .y =   20.00;             M[34] . Xp[33] .x=  0.007;          
 M[34] . Xp[34] .y =   20.60;             M[34] . Xp[34] .x=  0.003;          
 M[34] . Xp[35] .y =   21.20;             M[34] . Xp[35] .x=  0.003;          
 M[34] . Xp[36] .y =   21.80;             M[34] . Xp[36] .x=  0.008;          
 M[34] . Xp[37] .y =   22.40;             M[34] . Xp[37] .x=  0.002;          
 M[34] . Xp[38] .y =   23.00;             M[34] . Xp[38] .x=  0.004;          
 M[34] . Xp[39] .y =   23.60;             M[34] . Xp[39] .x=  0.006;          
 M[34] . Xp[40] .y =   24.20;             M[34] . Xp[40] .x=  0.005;          
 M[34] . Xp[41] .y =   24.80;             M[34] . Xp[41] .x=  0.003;          
 M[34] . Xp[42] .y =   25.40;             M[34] . Xp[42] .x=  0.004;          
   M[34].pathlength();
//station 36
 x[35] = 52.8;
  M[35].allocate(43);
  M[35].settype(SPLINE2D_LINEAR);                                                                                                                                                       
 M[35] . Xp[ 0] .y =    0.20;             M[35] . Xp[ 0] .x=  0.407;          
 M[35] . Xp[ 1] .y =    0.80;             M[35] . Xp[ 1] .x=  0.402;          
 M[35] . Xp[ 2] .y =    1.40;             M[35] . Xp[ 2] .x=  0.349;          
 M[35] . Xp[ 3] .y =    2.00;             M[35] . Xp[ 3] .x=  0.350;          
 M[35] . Xp[ 4] .y =    2.60;             M[35] . Xp[ 4] .x=  0.343;          
 M[35] . Xp[ 5] .y =    3.20;             M[35] . Xp[ 5] .x=  0.320;          
 M[35] . Xp[ 6] .y =    3.80;             M[35] . Xp[ 6] .x=  0.304;          
 M[35] . Xp[ 7] .y =    4.40;             M[35] . Xp[ 7] .x=  0.300;          
 M[35] . Xp[ 8] .y =    5.00;             M[35] . Xp[ 8] .x=  0.273;          
 M[35] . Xp[ 9] .y =    5.60;             M[35] . Xp[ 9] .x=  0.250;          
 M[35] . Xp[10] .y =    6.20;             M[35] . Xp[10] .x=  0.200;          
 M[35] . Xp[11] .y =    6.80;             M[35] . Xp[11] .x=  0.147;          
 M[35] . Xp[12] .y =    7.40;             M[35] . Xp[12] .x=  0.145;          
 M[35] . Xp[13] .y =    8.00;             M[35] . Xp[13] .x=  0.128;          
 M[35] . Xp[14] .y =    8.60;             M[35] . Xp[14] .x=  0.128;          
 M[35] . Xp[15] .y =    9.20;             M[35] . Xp[15] .x=  0.103;          
 M[35] . Xp[16] .y =    9.80;             M[35] . Xp[16] .x=  0.082;          
 M[35] . Xp[17] .y =   10.40;             M[35] . Xp[17] .x=  0.070;          
 M[35] . Xp[18] .y =   11.00;             M[35] . Xp[18] .x=  0.055;          
 M[35] . Xp[19] .y =   11.60;             M[35] . Xp[19] .x=  0.047;          
 M[35] . Xp[20] .y =   12.20;             M[35] . Xp[20] .x=  0.040;          
 M[35] . Xp[21] .y =   12.80;             M[35] . Xp[21] .x=  0.041;          
 M[35] . Xp[22] .y =   13.40;             M[35] . Xp[22] .x=  0.037;          
 M[35] . Xp[23] .y =   14.00;             M[35] . Xp[23] .x=  0.041;          
 M[35] . Xp[24] .y =   14.60;             M[35] . Xp[24] .x=  0.043;          
 M[35] . Xp[25] .y =   15.20;             M[35] . Xp[25] .x=  0.038;          
 M[35] . Xp[26] .y =   15.80;             M[35] . Xp[26] .x=  0.025;          
 M[35] . Xp[27] .y =   16.40;             M[35] . Xp[27] .x=  0.022;          
 M[35] . Xp[28] .y =   17.00;             M[35] . Xp[28] .x=  0.021;          
 M[35] . Xp[29] .y =   17.60;             M[35] . Xp[29] .x=  0.026;          
 M[35] . Xp[30] .y =   18.20;             M[35] . Xp[30] .x=  0.013;          
 M[35] . Xp[31] .y =   18.80;             M[35] . Xp[31] .x=  0.015;          
 M[35] . Xp[32] .y =   19.40;             M[35] . Xp[32] .x=  0.013;          
 M[35] . Xp[33] .y =   20.00;             M[35] . Xp[33] .x=  0.010;          
 M[35] . Xp[34] .y =   20.60;             M[35] . Xp[34] .x=  0.005;          
 M[35] . Xp[35] .y =   21.20;             M[35] . Xp[35] .x=  0.004;          
 M[35] . Xp[36] .y =   21.80;             M[35] . Xp[36] .x=  0.003;          
 M[35] . Xp[37] .y =   22.40;             M[35] . Xp[37] .x=  0.004;          
 M[35] . Xp[38] .y =   23.00;             M[35] . Xp[38] .x=  0.007;          
 M[35] . Xp[39] .y =   23.60;             M[35] . Xp[39] .x=  0.003;          
 M[35] . Xp[40] .y =   24.20;             M[35] . Xp[40] .x=  0.006;          
 M[35] . Xp[41] .y =   24.80;             M[35] . Xp[41] .x=  0.002;          
 M[35] . Xp[42] .y =   25.40;             M[35] . Xp[42] .x=  0.012;          
  M[35].pathlength();
//station 37
 x[36] = 55.0;
  M[36].allocate(43);
  M[36].settype(SPLINE2D_LINEAR);                                                                             
                                                                      
 M[36] . Xp[ 0] .y =    0.00;             M[36] . Xp[ 0] .x=  0.359;          
 M[36] . Xp[ 1] .y =    0.60;             M[36] . Xp[ 1] .x=  0.319;          
 M[36] . Xp[ 2] .y =    1.20;             M[36] . Xp[ 2] .x=  0.327;          
 M[36] . Xp[ 3] .y =    1.80;             M[36] . Xp[ 3] .x=  0.334;          
 M[36] . Xp[ 4] .y =    2.40;             M[36] . Xp[ 4] .x=  0.305;          
 M[36] . Xp[ 5] .y =    3.00;             M[36] . Xp[ 5] .x=  0.263;          
 M[36] . Xp[ 6] .y =    3.60;             M[36] . Xp[ 6] .x=  0.235;          
 M[36] . Xp[ 7] .y =    4.20;             M[36] . Xp[ 7] .x=  0.228;          
 M[36] . Xp[ 8] .y =    4.80;             M[36] . Xp[ 8] .x=  0.204;          
 M[36] . Xp[ 9] .y =    5.40;             M[36] . Xp[ 9] .x=  0.200;          
 M[36] . Xp[10] .y =    6.00;             M[36] . Xp[10] .x=  0.154;          
 M[36] . Xp[11] .y =    6.60;             M[36] . Xp[11] .x=  0.145;          
 M[36] . Xp[12] .y =    7.20;             M[36] . Xp[12] .x=  0.132;          
 M[36] . Xp[13] .y =    7.80;             M[36] . Xp[13] .x=  0.127;          
 M[36] . Xp[14] .y =    8.40;             M[36] . Xp[14] .x=  0.111;          
 M[36] . Xp[15] .y =    9.00;             M[36] . Xp[15] .x=  0.088;          
 M[36] . Xp[16] .y =    9.60;             M[36] . Xp[16] .x=  0.089;          
 M[36] . Xp[17] .y =   10.20;             M[36] . Xp[17] .x=  0.069;          
 M[36] . Xp[18] .y =   10.80;             M[36] . Xp[18] .x=  0.048;          
 M[36] . Xp[19] .y =   11.40;             M[36] . Xp[19] .x=  0.055;          
 M[36] . Xp[20] .y =   12.00;             M[36] . Xp[20] .x=  0.041;          
 M[36] . Xp[21] .y =   12.60;             M[36] . Xp[21] .x=  0.040;          
 M[36] . Xp[22] .y =   13.20;             M[36] . Xp[22] .x=  0.033;          
 M[36] . Xp[23] .y =   13.80;             M[36] . Xp[23] .x=  0.029;          
 M[36] . Xp[24] .y =   14.40;             M[36] . Xp[24] .x=  0.037;          
 M[36] . Xp[25] .y =   15.00;             M[36] . Xp[25] .x=  0.020;          
 M[36] . Xp[26] .y =   15.60;             M[36] . Xp[26] .x=  0.023;          
 M[36] . Xp[27] .y =   16.20;             M[36] . Xp[27] .x=  0.017;          
 M[36] . Xp[28] .y =   16.80;             M[36] . Xp[28] .x=  0.016;          
 M[36] . Xp[29] .y =   17.40;             M[36] . Xp[29] .x=  0.011;          
 M[36] . Xp[30] .y =   18.00;             M[36] . Xp[30] .x=  0.007;          
 M[36] . Xp[31] .y =   18.60;             M[36] . Xp[31] .x=  0.009;          
 M[36] . Xp[32] .y =   19.20;             M[36] . Xp[32] .x=  0.013;          
 M[36] . Xp[33] .y =   19.80;             M[36] . Xp[33] .x=  0.007;          
 M[36] . Xp[34] .y =   20.40;             M[36] . Xp[34] .x=  0.008;          
 M[36] . Xp[35] .y =   21.00;             M[36] . Xp[35] .x=  0.005;          
 M[36] . Xp[36] .y =   21.60;             M[36] . Xp[36] .x=  0.004;          
 M[36] . Xp[37] .y =   22.20;             M[36] . Xp[37] .x=  0.003;          
 M[36] . Xp[38] .y =   22.80;             M[36] . Xp[38] .x=  0.003;          
 M[36] . Xp[39] .y =   23.40;             M[36] . Xp[39] .x=  0.003;          
 M[36] . Xp[40] .y =   24.00;             M[36] . Xp[40] .x=  0.001;          
 M[36] . Xp[41] .y =   24.60;             M[36] . Xp[41] .x=  0.001;          
 M[36] . Xp[42] .y =   25.20;             M[36] . Xp[42] .x=  0.003;          
   M[36].pathlength();
//station 38
 x[37] = 55.8;
  M[37].allocate(43);
  M[37].settype(SPLINE2D_LINEAR);                                                                            
                                                                        
 M[37] . Xp[ 0] .y =    0.00;             M[37] . Xp[ 0] .x=  0.300;          
 M[37] . Xp[ 1] .y =    0.60;             M[37] . Xp[ 1] .x=  0.339;          
 M[37] . Xp[ 2] .y =    1.20;             M[37] . Xp[ 2] .x=  0.304;          
 M[37] . Xp[ 3] .y =    1.80;             M[37] . Xp[ 3] .x=  0.326;          
 M[37] . Xp[ 4] .y =    2.40;             M[37] . Xp[ 4] .x=  0.293;          
 M[37] . Xp[ 5] .y =    3.00;             M[37] . Xp[ 5] .x=  0.278;          
 M[37] . Xp[ 6] .y =    3.60;             M[37] . Xp[ 6] .x=  0.247;          
 M[37] . Xp[ 7] .y =    4.20;             M[37] . Xp[ 7] .x=  0.240;          
 M[37] . Xp[ 8] .y =    4.80;             M[37] . Xp[ 8] .x=  0.196;          
 M[37] . Xp[ 9] .y =    5.40;             M[37] . Xp[ 9] .x=  0.162;          
 M[37] . Xp[10] .y =    6.00;             M[37] . Xp[10] .x=  0.153;          
 M[37] . Xp[11] .y =    6.60;             M[37] . Xp[11] .x=  0.138;          
 M[37] . Xp[12] .y =    7.20;             M[37] . Xp[12] .x=  0.124;          
 M[37] . Xp[13] .y =    7.80;             M[37] . Xp[13] .x=  0.127;          
 M[37] . Xp[14] .y =    8.40;             M[37] . Xp[14] .x=  0.102;          
 M[37] . Xp[15] .y =    9.00;             M[37] . Xp[15] .x=  0.090;          
 M[37] . Xp[16] .y =    9.60;             M[37] . Xp[16] .x=  0.082;          
 M[37] . Xp[17] .y =   10.20;             M[37] . Xp[17] .x=  0.069;          
 M[37] . Xp[18] .y =   10.80;             M[37] . Xp[18] .x=  0.056;          
 M[37] . Xp[19] .y =   11.40;             M[37] . Xp[19] .x=  0.046;          
 M[37] . Xp[20] .y =   12.00;             M[37] . Xp[20] .x=  0.043;          
 M[37] . Xp[21] .y =   12.60;             M[37] . Xp[21] .x=  0.035;          
 M[37] . Xp[22] .y =   13.20;             M[37] . Xp[22] .x=  0.028;          
 M[37] . Xp[23] .y =   13.80;             M[37] . Xp[23] .x=  0.030;          
 M[37] . Xp[24] .y =   14.40;             M[37] . Xp[24] .x=  0.024;          
 M[37] . Xp[25] .y =   15.00;             M[37] . Xp[25] .x=  0.025;          
 M[37] . Xp[26] .y =   15.60;             M[37] . Xp[26] .x=  0.021;          
 M[37] . Xp[27] .y =   16.20;             M[37] . Xp[27] .x=  0.015;          
 M[37] . Xp[28] .y =   16.80;             M[37] . Xp[28] .x=  0.010;          
 M[37] . Xp[29] .y =   17.40;             M[37] . Xp[29] .x=  0.011;          
 M[37] . Xp[30] .y =   18.00;             M[37] . Xp[30] .x=  0.006;          
 M[37] . Xp[31] .y =   18.60;             M[37] . Xp[31] .x=  0.006;          
 M[37] . Xp[32] .y =   19.20;             M[37] . Xp[32] .x=  0.006;          
 M[37] . Xp[33] .y =   19.80;             M[37] . Xp[33] .x=  0.004;          
 M[37] . Xp[34] .y =   20.40;             M[37] . Xp[34] .x=  0.003;          
 M[37] . Xp[35] .y =   21.00;             M[37] . Xp[35] .x=  0.002;          
 M[37] . Xp[36] .y =   21.60;             M[37] . Xp[36] .x=  0.001;          
 M[37] . Xp[37] .y =   22.20;             M[37] . Xp[37] .x=  0.003;          
 M[37] . Xp[38] .y =   22.80;             M[37] . Xp[38] .x=  0.001;          
 M[37] . Xp[39] .y =   23.40;             M[37] . Xp[39] .x=  0.001;          
 M[37] . Xp[40] .y =   24.00;             M[37] . Xp[40] .x=  0.003;          
 M[37] . Xp[41] .y =   24.60;             M[37] . Xp[41] .x=  0.002;          
 M[37] . Xp[42] .y =   25.20;             M[37] . Xp[42] .x=  0.001;          
    M[37].pathlength();
//station 39
 x[38] = 56.5;
  M[38].allocate(43);
  M[38].settype(SPLINE2D_LINEAR);                                                                           
                                                                         
 M[38] . Xp[ 0] .y =    0.00;             M[38] . Xp[ 0] .x=  0.352;          
 M[38] . Xp[ 1] .y =    0.60;             M[38] . Xp[ 1] .x=  0.314;          
 M[38] . Xp[ 2] .y =    1.20;             M[38] . Xp[ 2] .x=  0.331;          
 M[38] . Xp[ 3] .y =    1.80;             M[38] . Xp[ 3] .x=  0.323;          
 M[38] . Xp[ 4] .y =    2.40;             M[38] . Xp[ 4] .x=  0.308;          
 M[38] . Xp[ 5] .y =    3.00;             M[38] . Xp[ 5] .x=  0.293;          
 M[38] . Xp[ 6] .y =    3.60;             M[38] . Xp[ 6] .x=  0.260;          
 M[38] . Xp[ 7] .y =    4.20;             M[38] . Xp[ 7] .x=  0.231;          
 M[38] . Xp[ 8] .y =    4.80;             M[38] . Xp[ 8] .x=  0.196;          
 M[38] . Xp[ 9] .y =    5.40;             M[38] . Xp[ 9] .x=  0.155;          
 M[38] . Xp[10] .y =    6.00;             M[38] . Xp[10] .x=  0.175;          
 M[38] . Xp[11] .y =    6.60;             M[38] . Xp[11] .x=  0.148;          
 M[38] . Xp[12] .y =    7.20;             M[38] . Xp[12] .x=  0.140;          
 M[38] . Xp[13] .y =    7.80;             M[38] . Xp[13] .x=  0.120;          
 M[38] . Xp[14] .y =    8.40;             M[38] . Xp[14] .x=  0.109;          
 M[38] . Xp[15] .y =    9.00;             M[38] . Xp[15] .x=  0.082;          
 M[38] . Xp[16] .y =    9.60;             M[38] . Xp[16] .x=  0.078;          
 M[38] . Xp[17] .y =   10.20;             M[38] . Xp[17] .x=  0.074;          
 M[38] . Xp[18] .y =   10.80;             M[38] . Xp[18] .x=  0.055;          
 M[38] . Xp[19] .y =   11.40;             M[38] . Xp[19] .x=  0.054;          
 M[38] . Xp[20] .y =   12.00;             M[38] . Xp[20] .x=  0.045;          
 M[38] . Xp[21] .y =   12.60;             M[38] . Xp[21] .x=  0.035;          
 M[38] . Xp[22] .y =   13.20;             M[38] . Xp[22] .x=  0.038;          
 M[38] . Xp[23] .y =   13.80;             M[38] . Xp[23] .x=  0.034;          
 M[38] . Xp[24] .y =   14.40;             M[38] . Xp[24] .x=  0.028;          
 M[38] . Xp[25] .y =   15.00;             M[38] . Xp[25] .x=  0.029;          
 M[38] . Xp[26] .y =   15.60;             M[38] . Xp[26] .x=  0.022;          
 M[38] . Xp[27] .y =   16.20;             M[38] . Xp[27] .x=  0.019;          
 M[38] . Xp[28] .y =   16.80;             M[38] . Xp[28] .x=  0.019;          
 M[38] . Xp[29] .y =   17.40;             M[38] . Xp[29] .x=  0.017;          
 M[38] . Xp[30] .y =   18.00;             M[38] . Xp[30] .x=  0.012;          
 M[38] . Xp[31] .y =   18.60;             M[38] . Xp[31] .x=  0.009;          
 M[38] . Xp[32] .y =   19.20;             M[38] . Xp[32] .x=  0.008;          
 M[38] . Xp[33] .y =   19.80;             M[38] . Xp[33] .x=  0.007;          
 M[38] . Xp[34] .y =   20.40;             M[38] . Xp[34] .x=  0.004;          
 M[38] . Xp[35] .y =   21.00;             M[38] . Xp[35] .x=  0.004;          
 M[38] . Xp[36] .y =   21.60;             M[38] . Xp[36] .x=  0.004;          
 M[38] . Xp[37] .y =   22.20;             M[38] . Xp[37] .x=  0.003;          
 M[38] . Xp[38] .y =   22.80;             M[38] . Xp[38] .x=  0.002;          
 M[38] . Xp[39] .y =   23.40;             M[38] . Xp[39] .x=  0.005;          
 M[38] . Xp[40] .y =   24.00;             M[38] . Xp[40] .x=  0.004;          
 M[38] . Xp[41] .y =   24.60;             M[38] . Xp[41] .x=  0.005;          
 M[38] . Xp[42] .y =   25.20;             M[38] . Xp[42] .x=  0.004;          
  M[38].pathlength();
//station 40
 x[39] = 57.3;
  M[39].allocate(43);
  M[39].settype(SPLINE2D_LINEAR);                                                                             
                                                                         
 M[39] . Xp[ 0] .y =    0.00;             M[39] . Xp[ 0] .x=  0.323;          
 M[39] . Xp[ 1] .y =    0.60;             M[39] . Xp[ 1] .x=  0.315;          
 M[39] . Xp[ 2] .y =    1.20;             M[39] . Xp[ 2] .x=  0.347;          
 M[39] . Xp[ 3] .y =    1.80;             M[39] . Xp[ 3] .x=  0.328;          
 M[39] . Xp[ 4] .y =    2.40;             M[39] . Xp[ 4] .x=  0.301;          
 M[39] . Xp[ 5] .y =    3.00;             M[39] . Xp[ 5] .x=  0.290;          
 M[39] . Xp[ 6] .y =    3.60;             M[39] . Xp[ 6] .x=  0.279;          
 M[39] . Xp[ 7] .y =    4.20;             M[39] . Xp[ 7] .x=  0.219;          
 M[39] . Xp[ 8] .y =    4.80;             M[39] . Xp[ 8] .x=  0.197;          
 M[39] . Xp[ 9] .y =    5.40;             M[39] . Xp[ 9] .x=  0.172;          
 M[39] . Xp[10] .y =    6.00;             M[39] . Xp[10] .x=  0.167;          
 M[39] . Xp[11] .y =    6.60;             M[39] . Xp[11] .x=  0.142;          
 M[39] . Xp[12] .y =    7.20;             M[39] . Xp[12] .x=  0.131;          
 M[39] . Xp[13] .y =    7.80;             M[39] . Xp[13] .x=  0.107;          
 M[39] . Xp[14] .y =    8.40;             M[39] . Xp[14] .x=  0.091;          
 M[39] . Xp[15] .y =    9.00;             M[39] . Xp[15] .x=  0.086;          
 M[39] . Xp[16] .y =    9.60;             M[39] . Xp[16] .x=  0.085;          
 M[39] . Xp[17] .y =   10.20;             M[39] . Xp[17] .x=  0.076;          
 M[39] . Xp[18] .y =   10.80;             M[39] . Xp[18] .x=  0.054;          
 M[39] . Xp[19] .y =   11.40;             M[39] . Xp[19] .x=  0.046;          
 M[39] . Xp[20] .y =   12.00;             M[39] . Xp[20] .x=  0.040;          
 M[39] . Xp[21] .y =   12.60;             M[39] . Xp[21] .x=  0.036;          
 M[39] . Xp[22] .y =   13.20;             M[39] . Xp[22] .x=  0.027;          
 M[39] . Xp[23] .y =   13.80;             M[39] . Xp[23] .x=  0.026;          
 M[39] . Xp[24] .y =   14.40;             M[39] . Xp[24] .x=  0.027;          
 M[39] . Xp[25] .y =   15.00;             M[39] . Xp[25] .x=  0.021;          
 M[39] . Xp[26] .y =   15.60;             M[39] . Xp[26] .x=  0.021;          
 M[39] . Xp[27] .y =   16.20;             M[39] . Xp[27] .x=  0.015;          
 M[39] . Xp[28] .y =   16.80;             M[39] . Xp[28] .x=  0.018;          
 M[39] . Xp[29] .y =   17.40;             M[39] . Xp[29] .x=  0.011;          
 M[39] . Xp[30] .y =   18.00;             M[39] . Xp[30] .x=  0.011;          
 M[39] . Xp[31] .y =   18.60;             M[39] . Xp[31] .x=  0.009;          
 M[39] . Xp[32] .y =   19.20;             M[39] . Xp[32] .x=  0.005;          
 M[39] . Xp[33] .y =   19.80;             M[39] . Xp[33] .x=  0.003;          
 M[39] . Xp[34] .y =   20.40;             M[39] . Xp[34] .x=  0.002;          
 M[39] . Xp[35] .y =   21.00;             M[39] . Xp[35] .x=  0.002;          
 M[39] . Xp[36] .y =   21.60;             M[39] . Xp[36] .x=  0.001;          
 M[39] . Xp[37] .y =   22.20;             M[39] . Xp[37] .x=  0.004;          
 M[39] . Xp[38] .y =   22.80;             M[39] . Xp[38] .x=  0.001;          
 M[39] . Xp[39] .y =   23.40;             M[39] . Xp[39] .x=  0.002;          
 M[39] . Xp[40] .y =   24.00;             M[39] . Xp[40] .x=  0.003;          
 M[39] . Xp[41] .y =   24.60;             M[39] . Xp[41] .x=  0.004;          
 M[39] . Xp[42] .y =   25.20;             M[39] . Xp[42] .x=  0.003;          
   M[39].pathlength();
//station 41
 x[40] = 58.0;
  M[40].allocate(43);
  M[40].settype(SPLINE2D_LINEAR);                                                                            
                                                                         
 M[40] . Xp[ 0] .y =    0.00;             M[40] . Xp[ 0] .x=  0.328;          
 M[40] . Xp[ 1] .y =    0.60;             M[40] . Xp[ 1] .x=  0.340;          
 M[40] . Xp[ 2] .y =    1.20;             M[40] . Xp[ 2] .x=  0.308;          
 M[40] . Xp[ 3] .y =    1.80;             M[40] . Xp[ 3] .x=  0.301;          
 M[40] . Xp[ 4] .y =    2.40;             M[40] . Xp[ 4] .x=  0.300;          
 M[40] . Xp[ 5] .y =    3.00;             M[40] . Xp[ 5] .x=  0.287;          
 M[40] . Xp[ 6] .y =    3.60;             M[40] . Xp[ 6] .x=  0.271;          
 M[40] . Xp[ 7] .y =    4.20;             M[40] . Xp[ 7] .x=  0.201;          
 M[40] . Xp[ 8] .y =    4.80;             M[40] . Xp[ 8] .x=  0.194;          
 M[40] . Xp[ 9] .y =    5.40;             M[40] . Xp[ 9] .x=  0.180;          
 M[40] . Xp[10] .y =    6.00;             M[40] . Xp[10] .x=  0.172;          
 M[40] . Xp[11] .y =    6.60;             M[40] . Xp[11] .x=  0.145;          
 M[40] . Xp[12] .y =    7.20;             M[40] . Xp[12] .x=  0.142;          
 M[40] . Xp[13] .y =    7.80;             M[40] . Xp[13] .x=  0.121;          
 M[40] . Xp[14] .y =    8.40;             M[40] . Xp[14] .x=  0.107;          
 M[40] . Xp[15] .y =    9.00;             M[40] . Xp[15] .x=  0.082;          
 M[40] . Xp[16] .y =    9.60;             M[40] . Xp[16] .x=  0.075;          
 M[40] . Xp[17] .y =   10.20;             M[40] . Xp[17] .x=  0.062;          
 M[40] . Xp[18] .y =   10.80;             M[40] . Xp[18] .x=  0.055;          
 M[40] . Xp[19] .y =   11.40;             M[40] . Xp[19] .x=  0.045;          
 M[40] . Xp[20] .y =   12.00;             M[40] . Xp[20] .x=  0.039;          
 M[40] . Xp[21] .y =   12.60;             M[40] . Xp[21] .x=  0.032;          
 M[40] . Xp[22] .y =   13.20;             M[40] . Xp[22] .x=  0.031;          
 M[40] . Xp[23] .y =   13.80;             M[40] . Xp[23] .x=  0.028;          
 M[40] . Xp[24] .y =   14.40;             M[40] . Xp[24] .x=  0.019;          
 M[40] . Xp[25] .y =   15.00;             M[40] . Xp[25] .x=  0.020;          
 M[40] . Xp[26] .y =   15.60;             M[40] . Xp[26] .x=  0.021;          
 M[40] . Xp[27] .y =   16.20;             M[40] . Xp[27] .x=  0.016;          
 M[40] . Xp[28] .y =   16.80;             M[40] . Xp[28] .x=  0.016;          
 M[40] . Xp[29] .y =   17.40;             M[40] . Xp[29] .x=  0.011;          
 M[40] . Xp[30] .y =   18.00;             M[40] . Xp[30] .x=  0.010;          
 M[40] . Xp[31] .y =   18.60;             M[40] . Xp[31] .x=  0.008;          
 M[40] . Xp[32] .y =   19.20;             M[40] . Xp[32] .x=  0.006;          
 M[40] . Xp[33] .y =   19.80;             M[40] . Xp[33] .x=  0.004;          
 M[40] . Xp[34] .y =   20.40;             M[40] . Xp[34] .x=  0.002;          
 M[40] . Xp[35] .y =   21.00;             M[40] . Xp[35] .x=  0.004;          
 M[40] . Xp[36] .y =   21.60;             M[40] . Xp[36] .x=  0.001;          
 M[40] . Xp[37] .y =   22.20;             M[40] . Xp[37] .x=  0.001;          
 M[40] . Xp[38] .y =   22.80;             M[40] . Xp[38] .x=  0.002;          
 M[40] . Xp[39] .y =   23.40;             M[40] . Xp[39] .x=  0.000;          
 M[40] . Xp[40] .y =   24.00;             M[40] . Xp[40] .x=  0.002;          
 M[40] . Xp[41] .y =   24.60;             M[40] . Xp[41] .x=  0.003;          
 M[40] . Xp[42] .y =   25.20;             M[40] . Xp[42] .x=  0.001;          

  M[40].pathlength();
//station 42
 x[41] = 58.8;
  M[41].allocate(43);
  M[41].settype(SPLINE2D_LINEAR);                                                                             
                                                                         
 M[41] . Xp[ 0] .y =    0.00;             M[41] . Xp[ 0] .x=  0.293;          
 M[41] . Xp[ 1] .y =    0.60;             M[41] . Xp[ 1] .x=  0.311;          
 M[41] . Xp[ 2] .y =    1.20;             M[41] . Xp[ 2] .x=  0.315;          
 M[41] . Xp[ 3] .y =    1.80;             M[41] . Xp[ 3] .x=  0.301;          
 M[41] . Xp[ 4] .y =    2.40;             M[41] . Xp[ 4] .x=  0.287;          
 M[41] . Xp[ 5] .y =    3.00;             M[41] . Xp[ 5] .x=  0.280;          
 M[41] . Xp[ 6] .y =    3.60;             M[41] . Xp[ 6] .x=  0.235;          
 M[41] . Xp[ 7] .y =    4.20;             M[41] . Xp[ 7] .x=  0.211;          
 M[41] . Xp[ 8] .y =    4.80;             M[41] . Xp[ 8] .x=  0.204;          
 M[41] . Xp[ 9] .y =    5.40;             M[41] . Xp[ 9] .x=  0.173;          
 M[41] . Xp[10] .y =    6.00;             M[41] . Xp[10] .x=  0.156;          
 M[41] . Xp[11] .y =    6.60;             M[41] . Xp[11] .x=  0.139;          
 M[41] . Xp[12] .y =    7.20;             M[41] . Xp[12] .x=  0.128;          
 M[41] . Xp[13] .y =    7.80;             M[41] . Xp[13] .x=  0.117;          
 M[41] . Xp[14] .y =    8.40;             M[41] . Xp[14] .x=  0.098;          
 M[41] . Xp[15] .y =    9.00;             M[41] . Xp[15] .x=  0.090;          
 M[41] . Xp[16] .y =    9.60;             M[41] . Xp[16] .x=  0.075;          
 M[41] . Xp[17] .y =   10.20;             M[41] . Xp[17] .x=  0.051;          
 M[41] . Xp[18] .y =   10.80;             M[41] . Xp[18] .x=  0.053;          
 M[41] . Xp[19] .y =   11.40;             M[41] . Xp[19] .x=  0.042;          
 M[41] . Xp[20] .y =   12.00;             M[41] . Xp[20] .x=  0.034;          
 M[41] . Xp[21] .y =   12.60;             M[41] . Xp[21] .x=  0.033;          
 M[41] . Xp[22] .y =   13.20;             M[41] . Xp[22] .x=  0.031;          
 M[41] . Xp[23] .y =   13.80;             M[41] . Xp[23] .x=  0.026;          
 M[41] . Xp[24] .y =   14.40;             M[41] . Xp[24] .x=  0.023;          
 M[41] . Xp[25] .y =   15.00;             M[41] . Xp[25] .x=  0.020;          
 M[41] . Xp[26] .y =   15.60;             M[41] . Xp[26] .x=  0.019;          
 M[41] . Xp[27] .y =   16.20;             M[41] . Xp[27] .x=  0.017;          
 M[41] . Xp[28] .y =   16.80;             M[41] . Xp[28] .x=  0.012;          
 M[41] . Xp[29] .y =   17.40;             M[41] . Xp[29] .x=  0.014;          
 M[41] . Xp[30] .y =   18.00;             M[41] . Xp[30] .x=  0.007;          
 M[41] . Xp[31] .y =   18.60;             M[41] . Xp[31] .x=  0.005;          
 M[41] . Xp[32] .y =   19.20;             M[41] . Xp[32] .x=  0.004;          
 M[41] . Xp[33] .y =   19.80;             M[41] . Xp[33] .x=  0.003;          
 M[41] . Xp[34] .y =   20.40;             M[41] . Xp[34] .x=  0.003;          
 M[41] . Xp[35] .y =   21.00;             M[41] . Xp[35] .x=  0.000;          
 M[41] . Xp[36] .y =   21.60;             M[41] . Xp[36] .x=  0.004;          
 M[41] . Xp[37] .y =   22.20;             M[41] . Xp[37] .x=  0.002;          
 M[41] . Xp[38] .y =   22.80;             M[41] . Xp[38] .x=  0.001;          
 M[41] . Xp[39] .y =   23.40;             M[41] . Xp[39] .x=  0.001;          
 M[41] . Xp[40] .y =   24.00;             M[41] . Xp[40] .x=  0.001;          
 M[41] . Xp[41] .y =   24.60;             M[41] . Xp[41] .x=  0.002;          
 M[41] . Xp[42] .y =   25.20;             M[41] . Xp[42] .x=  0.003;          
   M[41].pathlength();
//station 43
 x[42] = 59.5;
  M[42].allocate(43);
  M[42].settype(SPLINE2D_LINEAR);                                                                            
                                                                         
 M[42] . Xp[ 0] .y =    0.00;             M[42] . Xp[ 0] .x=  0.301;          
 M[42] . Xp[ 1] .y =    0.60;             M[42] . Xp[ 1] .x=  0.302;          
 M[42] . Xp[ 2] .y =    1.20;             M[42] . Xp[ 2] .x=  0.299;          
 M[42] . Xp[ 3] .y =    1.80;             M[42] . Xp[ 3] .x=  0.303;          
 M[42] . Xp[ 4] .y =    2.40;             M[42] . Xp[ 4] .x=  0.275;          
 M[42] . Xp[ 5] .y =    3.00;             M[42] . Xp[ 5] .x=  0.285;          
 M[42] . Xp[ 6] .y =    3.60;             M[42] . Xp[ 6] .x=  0.251;          
 M[42] . Xp[ 7] .y =    4.20;             M[42] . Xp[ 7] .x=  0.219;          
 M[42] . Xp[ 8] .y =    4.80;             M[42] . Xp[ 8] .x=  0.197;          
 M[42] . Xp[ 9] .y =    5.40;             M[42] . Xp[ 9] .x=  0.186;          
 M[42] . Xp[10] .y =    6.00;             M[42] . Xp[10] .x=  0.169;          
 M[42] . Xp[11] .y =    6.60;             M[42] . Xp[11] .x=  0.149;          
 M[42] . Xp[12] .y =    7.20;             M[42] . Xp[12] .x=  0.138;          
 M[42] . Xp[13] .y =    7.80;             M[42] . Xp[13] .x=  0.111;          
 M[42] . Xp[14] .y =    8.40;             M[42] . Xp[14] .x=  0.094;          
 M[42] . Xp[15] .y =    9.00;             M[42] . Xp[15] .x=  0.089;          
 M[42] . Xp[16] .y =    9.60;             M[42] . Xp[16] .x=  0.074;          
 M[42] . Xp[17] .y =   10.20;             M[42] . Xp[17] .x=  0.056;          
 M[42] . Xp[18] .y =   10.80;             M[42] . Xp[18] .x=  0.052;          
 M[42] . Xp[19] .y =   11.40;             M[42] . Xp[19] .x=  0.040;          
 M[42] . Xp[20] .y =   12.00;             M[42] . Xp[20] .x=  0.041;          
 M[42] . Xp[21] .y =   12.60;             M[42] . Xp[21] .x=  0.042;          
 M[42] . Xp[22] .y =   13.20;             M[42] . Xp[22] .x=  0.031;          
 M[42] . Xp[23] .y =   13.80;             M[42] . Xp[23] .x=  0.028;          
 M[42] . Xp[24] .y =   14.40;             M[42] . Xp[24] .x=  0.023;          
 M[42] . Xp[25] .y =   15.00;             M[42] . Xp[25] .x=  0.023;          
 M[42] . Xp[26] .y =   15.60;             M[42] . Xp[26] .x=  0.020;          
 M[42] . Xp[27] .y =   16.20;             M[42] . Xp[27] .x=  0.016;          
 M[42] . Xp[28] .y =   16.80;             M[42] . Xp[28] .x=  0.013;          
 M[42] . Xp[29] .y =   17.40;             M[42] . Xp[29] .x=  0.013;          
 M[42] . Xp[30] .y =   18.00;             M[42] . Xp[30] .x=  0.010;          
 M[42] . Xp[31] .y =   18.60;             M[42] . Xp[31] .x=  0.009;          
 M[42] . Xp[32] .y =   19.20;             M[42] . Xp[32] .x=  0.005;          
 M[42] . Xp[33] .y =   19.80;             M[42] . Xp[33] .x=  0.004;          
 M[42] . Xp[34] .y =   20.40;             M[42] . Xp[34] .x=  0.003;          
 M[42] . Xp[35] .y =   21.00;             M[42] . Xp[35] .x=  0.003;          
 M[42] . Xp[36] .y =   21.60;             M[42] . Xp[36] .x=  0.002;          
 M[42] . Xp[37] .y =   22.20;             M[42] . Xp[37] .x=  0.003;          
 M[42] . Xp[38] .y =   22.80;             M[42] . Xp[38] .x=  0.001;          
 M[42] . Xp[39] .y =   23.40;             M[42] . Xp[39] .x=  0.001;          
 M[42] . Xp[40] .y =   24.00;             M[42] . Xp[40] .x=  0.006;          
 M[42] . Xp[41] .y =   24.60;             M[42] . Xp[41] .x=  0.002;          
 M[42] . Xp[42] .y =   25.20;             M[42] . Xp[42] .x=  0.003;          
  M[42].pathlength();
//station 44
 x[43] = 60.3;
  M[43].allocate(43);
  M[43].settype(SPLINE2D_LINEAR);                                                                            
                                                                      
 M[43] . Xp[ 0] .y =    0.00;             M[43] . Xp[ 0] .x=  0.292;          
 M[43] . Xp[ 1] .y =    0.60;             M[43] . Xp[ 1] .x=  0.306;          
 M[43] . Xp[ 2] .y =    1.20;             M[43] . Xp[ 2] .x=  0.300;          
 M[43] . Xp[ 3] .y =    1.80;             M[43] . Xp[ 3] .x=  0.288;          
 M[43] . Xp[ 4] .y =    2.40;             M[43] . Xp[ 4] .x=  0.251;          
 M[43] . Xp[ 5] .y =    3.00;             M[43] . Xp[ 5] .x=  0.269;          
 M[43] . Xp[ 6] .y =    3.60;             M[43] . Xp[ 6] .x=  0.226;          
 M[43] . Xp[ 7] .y =    4.20;             M[43] . Xp[ 7] .x=  0.216;          
 M[43] . Xp[ 8] .y =    4.80;             M[43] . Xp[ 8] .x=  0.198;          
 M[43] . Xp[ 9] .y =    5.40;             M[43] . Xp[ 9] .x=  0.184;          
 M[43] . Xp[10] .y =    6.00;             M[43] . Xp[10] .x=  0.155;          
 M[43] . Xp[11] .y =    6.60;             M[43] . Xp[11] .x=  0.150;          
 M[43] . Xp[12] .y =    7.20;             M[43] . Xp[12] .x=  0.138;          
 M[43] . Xp[13] .y =    7.80;             M[43] . Xp[13] .x=  0.122;          
 M[43] . Xp[14] .y =    8.40;             M[43] . Xp[14] .x=  0.096;          
 M[43] . Xp[15] .y =    9.00;             M[43] . Xp[15] .x=  0.082;          
 M[43] . Xp[16] .y =    9.60;             M[43] . Xp[16] .x=  0.073;          
 M[43] . Xp[17] .y =   10.20;             M[43] . Xp[17] .x=  0.055;          
 M[43] . Xp[18] .y =   10.80;             M[43] . Xp[18] .x=  0.053;          
 M[43] . Xp[19] .y =   11.40;             M[43] . Xp[19] .x=  0.049;          
 M[43] . Xp[20] .y =   12.00;             M[43] . Xp[20] .x=  0.044;          
 M[43] . Xp[21] .y =   12.60;             M[43] . Xp[21] .x=  0.040;          
 M[43] . Xp[22] .y =   13.20;             M[43] . Xp[22] .x=  0.038;          
 M[43] . Xp[23] .y =   13.80;             M[43] . Xp[23] .x=  0.033;          
 M[43] . Xp[24] .y =   14.40;             M[43] . Xp[24] .x=  0.023;          
 M[43] . Xp[25] .y =   15.00;             M[43] . Xp[25] .x=  0.023;          
 M[43] . Xp[26] .y =   15.60;             M[43] . Xp[26] .x=  0.024;          
 M[43] . Xp[27] .y =   16.20;             M[43] . Xp[27] .x=  0.015;          
 M[43] . Xp[28] .y =   16.80;             M[43] . Xp[28] .x=  0.015;          
 M[43] . Xp[29] .y =   17.40;             M[43] . Xp[29] .x=  0.012;          
 M[43] . Xp[30] .y =   18.00;             M[43] . Xp[30] .x=  0.009;          
 M[43] . Xp[31] .y =   18.60;             M[43] . Xp[31] .x=  0.010;          
 M[43] . Xp[32] .y =   19.20;             M[43] . Xp[32] .x=  0.006;          
 M[43] . Xp[33] .y =   19.80;             M[43] . Xp[33] .x=  0.004;          
 M[43] . Xp[34] .y =   20.40;             M[43] . Xp[34] .x=  0.005;          
 M[43] . Xp[35] .y =   21.00;             M[43] . Xp[35] .x=  0.003;          
 M[43] . Xp[36] .y =   21.60;             M[43] . Xp[36] .x=  0.002;          
 M[43] . Xp[37] .y =   22.20;             M[43] . Xp[37] .x=  0.002;          
 M[43] . Xp[38] .y =   22.80;             M[43] . Xp[38] .x=  0.003;          
 M[43] . Xp[39] .y =   23.40;             M[43] . Xp[39] .x=  0.001;          
 M[43] . Xp[40] .y =   24.00;             M[43] . Xp[40] .x=  0.003;          
 M[43] . Xp[41] .y =   24.60;             M[43] . Xp[41] .x=  0.003;          
 M[43] . Xp[42] .y =   25.20;             M[43] . Xp[42] .x=  0.003;          
  M[43].pathlength();
//station 45
 x[44] = 61.0;
  M[44].allocate(43);
  M[44].settype(SPLINE2D_LINEAR);                                                                             
                                                                          
 M[44] . Xp[ 0] .y =    0.00;             M[44] . Xp[ 0] .x=  0.296;          
 M[44] . Xp[ 1] .y =    0.60;             M[44] . Xp[ 1] .x=  0.301;          
 M[44] . Xp[ 2] .y =    1.20;             M[44] . Xp[ 2] .x=  0.272;          
 M[44] . Xp[ 3] .y =    1.80;             M[44] . Xp[ 3] .x=  0.255;          
 M[44] . Xp[ 4] .y =    2.40;             M[44] . Xp[ 4] .x=  0.248;          
 M[44] . Xp[ 5] .y =    3.00;             M[44] . Xp[ 5] .x=  0.243;          
 M[44] . Xp[ 6] .y =    3.60;             M[44] . Xp[ 6] .x=  0.228;          
 M[44] . Xp[ 7] .y =    4.20;             M[44] . Xp[ 7] .x=  0.205;          
 M[44] . Xp[ 8] .y =    4.80;             M[44] . Xp[ 8] .x=  0.186;          
 M[44] . Xp[ 9] .y =    5.40;             M[44] . Xp[ 9] .x=  0.160;          
 M[44] . Xp[10] .y =    6.00;             M[44] . Xp[10] .x=  0.147;          
 M[44] . Xp[11] .y =    6.60;             M[44] . Xp[11] .x=  0.132;          
 M[44] . Xp[12] .y =    7.20;             M[44] . Xp[12] .x=  0.121;          
 M[44] . Xp[13] .y =    7.80;             M[44] . Xp[13] .x=  0.097;          
 M[44] . Xp[14] .y =    8.40;             M[44] . Xp[14] .x=  0.087;          
 M[44] . Xp[15] .y =    9.00;             M[44] . Xp[15] .x=  0.078;          
 M[44] . Xp[16] .y =    9.60;             M[44] . Xp[16] .x=  0.056;          
 M[44] . Xp[17] .y =   10.20;             M[44] . Xp[17] .x=  0.054;          
 M[44] . Xp[18] .y =   10.80;             M[44] . Xp[18] .x=  0.054;          
 M[44] . Xp[19] .y =   11.40;             M[44] . Xp[19] .x=  0.057;          
 M[44] . Xp[20] .y =   12.00;             M[44] . Xp[20] .x=  0.042;          
 M[44] . Xp[21] .y =   12.60;             M[44] . Xp[21] .x=  0.034;          
 M[44] . Xp[22] .y =   13.20;             M[44] . Xp[22] .x=  0.026;          
 M[44] . Xp[23] .y =   13.80;             M[44] . Xp[23] .x=  0.026;          
 M[44] . Xp[24] .y =   14.40;             M[44] . Xp[24] .x=  0.023;          
 M[44] . Xp[25] .y =   15.00;             M[44] . Xp[25] .x=  0.020;          
 M[44] . Xp[26] .y =   15.60;             M[44] . Xp[26] .x=  0.020;          
 M[44] . Xp[27] .y =   16.20;             M[44] . Xp[27] .x=  0.013;          
 M[44] . Xp[28] .y =   16.80;             M[44] . Xp[28] .x=  0.013;          
 M[44] . Xp[29] .y =   17.40;             M[44] . Xp[29] .x=  0.009;          
 M[44] . Xp[30] .y =   18.00;             M[44] . Xp[30] .x=  0.008;          
 M[44] . Xp[31] .y =   18.60;             M[44] . Xp[31] .x=  0.007;          
 M[44] . Xp[32] .y =   19.20;             M[44] . Xp[32] .x=  0.003;          
 M[44] . Xp[33] .y =   19.80;             M[44] . Xp[33] .x=  0.004;          
 M[44] . Xp[34] .y =   20.40;             M[44] . Xp[34] .x=  0.004;          
 M[44] . Xp[35] .y =   21.00;             M[44] . Xp[35] .x=  0.002;          
 M[44] . Xp[36] .y =   21.60;             M[44] . Xp[36] .x=  0.001;          
 M[44] . Xp[37] .y =   22.20;             M[44] . Xp[37] .x=  0.002;          
 M[44] . Xp[38] .y =   22.80;             M[44] . Xp[38] .x=  0.002;          
 M[44] . Xp[39] .y =   23.40;             M[44] . Xp[39] .x=  0.003;          
 M[44] . Xp[40] .y =   24.00;             M[44] . Xp[40] .x=  0.003;          
 M[44] . Xp[41] .y =   24.60;             M[44] . Xp[41] .x=  0.002;          
 M[44] . Xp[42] .y =   25.20;             M[44] . Xp[42] .x=  0.003;          
  M[44].pathlength();
//station 46
 x[45] = 61.8;
  M[45].allocate(43);
  M[45].settype(SPLINE2D_LINEAR);                                                                             
                                                                          
 M[45] . Xp[ 0] .y =    0.00;             M[45] . Xp[ 0] .x=  0.296;          
 M[45] . Xp[ 1] .y =    0.60;             M[45] . Xp[ 1] .x=  0.278;          
 M[45] . Xp[ 2] .y =    1.20;             M[45] . Xp[ 2] .x=  0.262;          
 M[45] . Xp[ 3] .y =    1.80;             M[45] . Xp[ 3] .x=  0.238;          
 M[45] . Xp[ 4] .y =    2.40;             M[45] . Xp[ 4] .x=  0.244;          
 M[45] . Xp[ 5] .y =    3.00;             M[45] . Xp[ 5] .x=  0.228;          
 M[45] . Xp[ 6] .y =    3.60;             M[45] . Xp[ 6] .x=  0.211;          
 M[45] . Xp[ 7] .y =    4.20;             M[45] . Xp[ 7] .x=  0.197;          
 M[45] . Xp[ 8] .y =    4.80;             M[45] . Xp[ 8] .x=  0.166;          
 M[45] . Xp[ 9] .y =    5.40;             M[45] . Xp[ 9] .x=  0.149;          
 M[45] . Xp[10] .y =    6.00;             M[45] . Xp[10] .x=  0.129;          
 M[45] . Xp[11] .y =    6.60;             M[45] . Xp[11] .x=  0.127;          
 M[45] . Xp[12] .y =    7.20;             M[45] . Xp[12] .x=  0.116;          
 M[45] . Xp[13] .y =    7.80;             M[45] . Xp[13] .x=  0.100;          
 M[45] . Xp[14] .y =    8.40;             M[45] . Xp[14] .x=  0.084;          
 M[45] . Xp[15] .y =    9.00;             M[45] . Xp[15] .x=  0.075;          
 M[45] . Xp[16] .y =    9.60;             M[45] . Xp[16] .x=  0.058;          
 M[45] . Xp[17] .y =   10.20;             M[45] . Xp[17] .x=  0.058;          
 M[45] . Xp[18] .y =   10.80;             M[45] . Xp[18] .x=  0.047;          
 M[45] . Xp[19] .y =   11.40;             M[45] . Xp[19] .x=  0.044;          
 M[45] . Xp[20] .y =   12.00;             M[45] . Xp[20] .x=  0.037;          
 M[45] . Xp[21] .y =   12.60;             M[45] . Xp[21] .x=  0.028;          
 M[45] . Xp[22] .y =   13.20;             M[45] . Xp[22] .x=  0.029;          
 M[45] . Xp[23] .y =   13.80;             M[45] . Xp[23] .x=  0.026;          
 M[45] . Xp[24] .y =   14.40;             M[45] . Xp[24] .x=  0.019;          
 M[45] . Xp[25] .y =   15.00;             M[45] . Xp[25] .x=  0.016;          
 M[45] . Xp[26] .y =   15.60;             M[45] . Xp[26] .x=  0.017;          
 M[45] . Xp[27] .y =   16.20;             M[45] . Xp[27] .x=  0.011;          
 M[45] . Xp[28] .y =   16.80;             M[45] . Xp[28] .x=  0.011;          
 M[45] . Xp[29] .y =   17.40;             M[45] . Xp[29] .x=  0.008;          
 M[45] . Xp[30] .y =   18.00;             M[45] . Xp[30] .x=  0.007;          
 M[45] . Xp[31] .y =   18.60;             M[45] . Xp[31] .x=  0.006;          
 M[45] . Xp[32] .y =   19.20;             M[45] . Xp[32] .x=  0.004;          
 M[45] . Xp[33] .y =   19.80;             M[45] . Xp[33] .x=  0.003;          
 M[45] . Xp[34] .y =   20.40;             M[45] . Xp[34] .x=  0.002;          
 M[45] . Xp[35] .y =   21.00;             M[45] . Xp[35] .x=  0.001;          
 M[45] . Xp[36] .y =   21.60;             M[45] . Xp[36] .x=  0.002;          
 M[45] . Xp[37] .y =   22.20;             M[45] . Xp[37] .x=  0.001;          
 M[45] . Xp[38] .y =   22.80;             M[45] . Xp[38] .x=  0.003;          
 M[45] . Xp[39] .y =   23.40;             M[45] . Xp[39] .x=  0.001;          
 M[45] . Xp[40] .y =   24.00;             M[45] . Xp[40] .x=  0.001;          
 M[45] . Xp[41] .y =   24.60;             M[45] . Xp[41] .x=  0.003;          
 M[45] . Xp[42] .y =   25.20;             M[45] . Xp[42] .x=  0.002;          
  M[45].pathlength();
//station 47
 x[46] = 62.5;
  M[46].allocate(43);
  M[46].settype(SPLINE2D_LINEAR);                                                                             
                                                                         
 M[46] . Xp[ 0] .y =    0.00;             M[46] . Xp[ 0] .x=  0.271;          
 M[46] . Xp[ 1] .y =    0.60;             M[46] . Xp[ 1] .x=  0.281;          
 M[46] . Xp[ 2] .y =    1.20;             M[46] . Xp[ 2] .x=  0.253;          
 M[46] . Xp[ 3] .y =    1.80;             M[46] . Xp[ 3] .x=  0.227;          
 M[46] . Xp[ 4] .y =    2.40;             M[46] . Xp[ 4] .x=  0.208;          
 M[46] . Xp[ 5] .y =    3.00;             M[46] . Xp[ 5] .x=  0.208;          
 M[46] . Xp[ 6] .y =    3.60;             M[46] . Xp[ 6] .x=  0.185;          
 M[46] . Xp[ 7] .y =    4.20;             M[46] . Xp[ 7] .x=  0.167;          
 M[46] . Xp[ 8] .y =    4.80;             M[46] . Xp[ 8] .x=  0.161;          
 M[46] . Xp[ 9] .y =    5.40;             M[46] . Xp[ 9] .x=  0.139;          
 M[46] . Xp[10] .y =    6.00;             M[46] . Xp[10] .x=  0.136;          
 M[46] . Xp[11] .y =    6.60;             M[46] . Xp[11] .x=  0.123;          
 M[46] . Xp[12] .y =    7.20;             M[46] . Xp[12] .x=  0.110;          
 M[46] . Xp[13] .y =    7.80;             M[46] . Xp[13] .x=  0.099;          
 M[46] . Xp[14] .y =    8.40;             M[46] . Xp[14] .x=  0.100;          
 M[46] . Xp[15] .y =    9.00;             M[46] . Xp[15] .x=  0.075;          
 M[46] . Xp[16] .y =    9.60;             M[46] . Xp[16] .x=  0.061;          
 M[46] . Xp[17] .y =   10.20;             M[46] . Xp[17] .x=  0.056;          
 M[46] . Xp[18] .y =   10.80;             M[46] . Xp[18] .x=  0.050;          
 M[46] . Xp[19] .y =   11.40;             M[46] . Xp[19] .x=  0.036;          
 M[46] . Xp[20] .y =   12.00;             M[46] . Xp[20] .x=  0.032;          
 M[46] . Xp[21] .y =   12.60;             M[46] . Xp[21] .x=  0.030;          
 M[46] . Xp[22] .y =   13.20;             M[46] . Xp[22] .x=  0.029;          
 M[46] . Xp[23] .y =   13.80;             M[46] . Xp[23] .x=  0.026;          
 M[46] . Xp[24] .y =   14.40;             M[46] . Xp[24] .x=  0.024;          
 M[46] . Xp[25] .y =   15.00;             M[46] . Xp[25] .x=  0.018;          
 M[46] . Xp[26] .y =   15.60;             M[46] . Xp[26] .x=  0.015;          
 M[46] . Xp[27] .y =   16.20;             M[46] . Xp[27] .x=  0.012;          
 M[46] . Xp[28] .y =   16.80;             M[46] . Xp[28] .x=  0.012;          
 M[46] . Xp[29] .y =   17.40;             M[46] . Xp[29] .x=  0.011;          
 M[46] . Xp[30] .y =   18.00;             M[46] . Xp[30] .x=  0.007;          
 M[46] . Xp[31] .y =   18.60;             M[46] . Xp[31] .x=  0.004;          
 M[46] . Xp[32] .y =   19.20;             M[46] . Xp[32] .x=  0.004;          
 M[46] . Xp[33] .y =   19.80;             M[46] . Xp[33] .x=  0.004;          
 M[46] . Xp[34] .y =   20.40;             M[46] . Xp[34] .x=  0.002;          
 M[46] . Xp[35] .y =   21.00;             M[46] . Xp[35] .x=  0.003;          
 M[46] . Xp[36] .y =   21.60;             M[46] . Xp[36] .x=  0.004;          
 M[46] . Xp[37] .y =   22.20;             M[46] . Xp[37] .x=  0.001;          
 M[46] . Xp[38] .y =   22.80;             M[46] . Xp[38] .x=  0.002;          
 M[46] . Xp[39] .y =   23.40;             M[46] . Xp[39] .x=  0.002;          
 M[46] . Xp[40] .y =   24.00;             M[46] . Xp[40] .x=  0.002;          
 M[46] . Xp[41] .y =   24.60;             M[46] . Xp[41] .x=  0.003;          
 M[46] . Xp[42] .y =   25.20;             M[46] . Xp[42] .x=  0.005;          
   M[46].pathlength();
//station 48
 x[47] = 63.3;
  M[47].allocate(43);
  M[47].settype(SPLINE2D_LINEAR);                                                                            
                                                                      
 M[47] . Xp[ 0] .y =    0.00;             M[47] . Xp[ 0] .x=  0.270;          
 M[47] . Xp[ 1] .y =    0.60;             M[47] . Xp[ 1] .x=  0.251;          
 M[47] . Xp[ 2] .y =    1.20;             M[47] . Xp[ 2] .x=  0.248;          
 M[47] . Xp[ 3] .y =    1.80;             M[47] . Xp[ 3] .x=  0.200;          
 M[47] . Xp[ 4] .y =    2.40;             M[47] . Xp[ 4] .x=  0.183;          
 M[47] . Xp[ 5] .y =    3.00;             M[47] . Xp[ 5] .x=  0.179;          
 M[47] . Xp[ 6] .y =    3.60;             M[47] . Xp[ 6] .x=  0.170;          
 M[47] . Xp[ 7] .y =    4.20;             M[47] . Xp[ 7] .x=  0.160;          
 M[47] . Xp[ 8] .y =    4.80;             M[47] . Xp[ 8] .x=  0.139;          
 M[47] . Xp[ 9] .y =    5.40;             M[47] . Xp[ 9] .x=  0.140;          
 M[47] . Xp[10] .y =    6.00;             M[47] . Xp[10] .x=  0.141;          
 M[47] . Xp[11] .y =    6.60;             M[47] . Xp[11] .x=  0.121;          
 M[47] . Xp[12] .y =    7.20;             M[47] . Xp[12] .x=  0.097;          
 M[47] . Xp[13] .y =    7.80;             M[47] . Xp[13] .x=  0.092;          
 M[47] . Xp[14] .y =    8.40;             M[47] . Xp[14] .x=  0.082;          
 M[47] . Xp[15] .y =    9.00;             M[47] . Xp[15] .x=  0.073;          
 M[47] . Xp[16] .y =    9.60;             M[47] . Xp[16] .x=  0.066;          
 M[47] . Xp[17] .y =   10.20;             M[47] . Xp[17] .x=  0.057;          
 M[47] . Xp[18] .y =   10.80;             M[47] . Xp[18] .x=  0.043;          
 M[47] . Xp[19] .y =   11.40;             M[47] . Xp[19] .x=  0.037;          
 M[47] . Xp[20] .y =   12.00;             M[47] . Xp[20] .x=  0.027;          
 M[47] . Xp[21] .y =   12.60;             M[47] . Xp[21] .x=  0.027;          
 M[47] . Xp[22] .y =   13.20;             M[47] . Xp[22] .x=  0.026;          
 M[47] . Xp[23] .y =   13.80;             M[47] . Xp[23] .x=  0.025;          
 M[47] . Xp[24] .y =   14.40;             M[47] . Xp[24] .x=  0.024;          
 M[47] . Xp[25] .y =   15.00;             M[47] . Xp[25] .x=  0.016;          
 M[47] . Xp[26] .y =   15.60;             M[47] . Xp[26] .x=  0.015;          
 M[47] . Xp[27] .y =   16.20;             M[47] . Xp[27] .x=  0.014;          
 M[47] . Xp[28] .y =   16.80;             M[47] . Xp[28] .x=  0.006;          
 M[47] . Xp[29] .y =   17.40;             M[47] . Xp[29] .x=  0.010;          
 M[47] . Xp[30] .y =   18.00;             M[47] . Xp[30] .x=  0.004;          
 M[47] . Xp[31] .y =   18.60;             M[47] . Xp[31] .x=  0.004;          
 M[47] . Xp[32] .y =   19.20;             M[47] . Xp[32] .x=  0.003;          
 M[47] . Xp[33] .y =   19.80;             M[47] . Xp[33] .x=  0.002;          
 M[47] . Xp[34] .y =   20.40;             M[47] . Xp[34] .x=  0.000;          
 M[47] . Xp[35] .y =   21.00;             M[47] . Xp[35] .x=  0.003;          
 M[47] . Xp[36] .y =   21.60;             M[47] . Xp[36] .x=  0.001;          
 M[47] . Xp[37] .y =   22.20;             M[47] . Xp[37] .x=  0.003;          
 M[47] . Xp[38] .y =   22.80;             M[47] . Xp[38] .x=  0.002;          
 M[47] . Xp[39] .y =   23.40;             M[47] . Xp[39] .x=  0.001;          
 M[47] . Xp[40] .y =   24.00;             M[47] . Xp[40] .x=  0.001;          
 M[47] . Xp[41] .y =   24.60;             M[47] . Xp[41] .x=  0.002;          
 M[47] . Xp[42] .y =   25.20;             M[47] . Xp[42] .x=  0.003;          
  M[47].pathlength();
//station 49
 x[48] = 64.0;
  M[48].allocate(43);
  M[48].settype(SPLINE2D_LINEAR);                                                                            
                                                                          
 M[48] . Xp[ 0] .y =    0.00;             M[48] . Xp[ 0] .x=  0.243;          
 M[48] . Xp[ 1] .y =    0.60;             M[48] . Xp[ 1] .x=  0.229;          
 M[48] . Xp[ 2] .y =    1.20;             M[48] . Xp[ 2] .x=  0.212;          
 M[48] . Xp[ 3] .y =    1.80;             M[48] . Xp[ 3] .x=  0.189;          
 M[48] . Xp[ 4] .y =    2.40;             M[48] . Xp[ 4] .x=  0.161;          
 M[48] . Xp[ 5] .y =    3.00;             M[48] . Xp[ 5] .x=  0.159;          
 M[48] . Xp[ 6] .y =    3.60;             M[48] . Xp[ 6] .x=  0.155;          
 M[48] . Xp[ 7] .y =    4.20;             M[48] . Xp[ 7] .x=  0.137;          
 M[48] . Xp[ 8] .y =    4.80;             M[48] . Xp[ 8] .x=  0.131;          
 M[48] . Xp[ 9] .y =    5.40;             M[48] . Xp[ 9] .x=  0.141;          
 M[48] . Xp[10] .y =    6.00;             M[48] . Xp[10] .x=  0.129;          
 M[48] . Xp[11] .y =    6.60;             M[48] . Xp[11] .x=  0.115;          
 M[48] . Xp[12] .y =    7.20;             M[48] . Xp[12] .x=  0.097;          
 M[48] . Xp[13] .y =    7.80;             M[48] . Xp[13] .x=  0.093;          
 M[48] . Xp[14] .y =    8.40;             M[48] . Xp[14] .x=  0.071;          
 M[48] . Xp[15] .y =    9.00;             M[48] . Xp[15] .x=  0.067;          
 M[48] . Xp[16] .y =    9.60;             M[48] . Xp[16] .x=  0.054;          
 M[48] . Xp[17] .y =   10.20;             M[48] . Xp[17] .x=  0.058;          
 M[48] . Xp[18] .y =   10.80;             M[48] . Xp[18] .x=  0.044;          
 M[48] . Xp[19] .y =   11.40;             M[48] . Xp[19] .x=  0.036;          
 M[48] . Xp[20] .y =   12.00;             M[48] . Xp[20] .x=  0.031;          
 M[48] . Xp[21] .y =   12.60;             M[48] . Xp[21] .x=  0.025;          
 M[48] . Xp[22] .y =   13.20;             M[48] . Xp[22] .x=  0.022;          
 M[48] . Xp[23] .y =   13.80;             M[48] . Xp[23] .x=  0.023;          
 M[48] . Xp[24] .y =   14.40;             M[48] . Xp[24] .x=  0.025;          
 M[48] . Xp[25] .y =   15.00;             M[48] . Xp[25] .x=  0.019;          
 M[48] . Xp[26] .y =   15.60;             M[48] . Xp[26] .x=  0.016;          
 M[48] . Xp[27] .y =   16.20;             M[48] . Xp[27] .x=  0.012;          
 M[48] . Xp[28] .y =   16.80;             M[48] . Xp[28] .x=  0.011;          
 M[48] . Xp[29] .y =   17.40;             M[48] . Xp[29] .x=  0.007;          
 M[48] . Xp[30] .y =   18.00;             M[48] . Xp[30] .x=  0.004;          
 M[48] . Xp[31] .y =   18.60;             M[48] . Xp[31] .x=  0.005;          
 M[48] . Xp[32] .y =   19.20;             M[48] . Xp[32] .x=  0.003;          
 M[48] . Xp[33] .y =   19.80;             M[48] . Xp[33] .x=  0.003;          
 M[48] . Xp[34] .y =   20.40;             M[48] . Xp[34] .x=  0.004;          
 M[48] . Xp[35] .y =   21.00;             M[48] . Xp[35] .x=  0.002;          
 M[48] . Xp[36] .y =   21.60;             M[48] . Xp[36] .x=  0.002;          
 M[48] . Xp[37] .y =   22.20;             M[48] . Xp[37] .x=  0.003;          
 M[48] . Xp[38] .y =   22.80;             M[48] . Xp[38] .x=  0.002;          
 M[48] . Xp[39] .y =   23.40;             M[48] . Xp[39] .x=  0.000;          
 M[48] . Xp[40] .y =   24.00;             M[48] . Xp[40] .x=  0.002;          
 M[48] . Xp[41] .y =   24.60;             M[48] . Xp[41] .x=  0.001;          
 M[48] . Xp[42] .y =   25.20;             M[48] . Xp[42] .x=  0.004;          
  M[48].pathlength();
//station 50
 x[49] = 64.8;
  M[49].allocate(43);
  M[49].settype(SPLINE2D_LINEAR);                                                                             
                                                                         
 M[49] . Xp[ 0] .y =    0.00;             M[49] . Xp[ 0] .x=  0.228;          
 M[49] . Xp[ 1] .y =    0.60;             M[49] . Xp[ 1] .x=  0.222;          
 M[49] . Xp[ 2] .y =    1.20;             M[49] . Xp[ 2] .x=  0.196;          
 M[49] . Xp[ 3] .y =    1.80;             M[49] . Xp[ 3] .x=  0.189;          
 M[49] . Xp[ 4] .y =    2.40;             M[49] . Xp[ 4] .x=  0.171;          
 M[49] . Xp[ 5] .y =    3.00;             M[49] . Xp[ 5] .x=  0.148;          
 M[49] . Xp[ 6] .y =    3.60;             M[49] . Xp[ 6] .x=  0.133;          
 M[49] . Xp[ 7] .y =    4.20;             M[49] . Xp[ 7] .x=  0.132;          
 M[49] . Xp[ 8] .y =    4.80;             M[49] . Xp[ 8] .x=  0.134;          
 M[49] . Xp[ 9] .y =    5.40;             M[49] . Xp[ 9] .x=  0.122;          
 M[49] . Xp[10] .y =    6.00;             M[49] . Xp[10] .x=  0.111;          
 M[49] . Xp[11] .y =    6.60;             M[49] . Xp[11] .x=  0.107;          
 M[49] . Xp[12] .y =    7.20;             M[49] . Xp[12] .x=  0.096;          
 M[49] . Xp[13] .y =    7.80;             M[49] . Xp[13] .x=  0.080;          
 M[49] . Xp[14] .y =    8.40;             M[49] . Xp[14] .x=  0.076;          
 M[49] . Xp[15] .y =    9.00;             M[49] . Xp[15] .x=  0.055;          
 M[49] . Xp[16] .y =    9.60;             M[49] . Xp[16] .x=  0.052;          
 M[49] . Xp[17] .y =   10.20;             M[49] . Xp[17] .x=  0.046;          
 M[49] . Xp[18] .y =   10.80;             M[49] . Xp[18] .x=  0.044;          
 M[49] . Xp[19] .y =   11.40;             M[49] . Xp[19] .x=  0.035;          
 M[49] . Xp[20] .y =   12.00;             M[49] . Xp[20] .x=  0.034;          
 M[49] . Xp[21] .y =   12.60;             M[49] . Xp[21] .x=  0.023;          
 M[49] . Xp[22] .y =   13.20;             M[49] . Xp[22] .x=  0.022;          
 M[49] . Xp[23] .y =   13.80;             M[49] . Xp[23] .x=  0.018;          
 M[49] . Xp[24] .y =   14.40;             M[49] . Xp[24] .x=  0.024;          
 M[49] . Xp[25] .y =   15.00;             M[49] . Xp[25] .x=  0.022;          
 M[49] . Xp[26] .y =   15.60;             M[49] . Xp[26] .x=  0.014;          
 M[49] . Xp[27] .y =   16.20;             M[49] . Xp[27] .x=  0.014;          
 M[49] . Xp[28] .y =   16.80;             M[49] . Xp[28] .x=  0.010;          
 M[49] . Xp[29] .y =   17.40;             M[49] . Xp[29] .x=  0.006;          
 M[49] . Xp[30] .y =   18.00;             M[49] . Xp[30] .x=  0.005;          
 M[49] . Xp[31] .y =   18.60;             M[49] . Xp[31] .x=  0.005;          
 M[49] . Xp[32] .y =   19.20;             M[49] . Xp[32] .x=  0.002;          
 M[49] . Xp[33] .y =   19.80;             M[49] . Xp[33] .x=  0.001;          
 M[49] . Xp[34] .y =   20.40;             M[49] . Xp[34] .x=  0.002;          
 M[49] . Xp[35] .y =   21.00;             M[49] . Xp[35] .x=  0.000;          
 M[49] . Xp[36] .y =   21.60;             M[49] . Xp[36] .x=  0.003;          
 M[49] . Xp[37] .y =   22.20;             M[49] . Xp[37] .x=  0.000;          
 M[49] . Xp[38] .y =   22.80;             M[49] . Xp[38] .x=  0.002;          
 M[49] . Xp[39] .y =   23.40;             M[49] . Xp[39] .x=  0.002;          
 M[49] . Xp[40] .y =   24.00;             M[49] . Xp[40] .x=  0.002;          
 M[49] . Xp[41] .y =   24.60;             M[49] . Xp[41] .x=  0.000;          
 M[49] . Xp[42] .y =   25.20;             M[49] . Xp[42] .x=  0.005;          
  M[49].pathlength();
//station 51
 x[50] = 65.5;
  M[50].allocate(43);
  M[50].settype(SPLINE2D_LINEAR);                                                                             
                                                                        
 M[50] . Xp[ 0] .y =    0.00;             M[50] . Xp[ 0] .x=  0.224;          
 M[50] . Xp[ 1] .y =    0.60;             M[50] . Xp[ 1] .x=  0.215;          
 M[50] . Xp[ 2] .y =    1.20;             M[50] . Xp[ 2] .x=  0.186;          
 M[50] . Xp[ 3] .y =    1.80;             M[50] . Xp[ 3] .x=  0.186;          
 M[50] . Xp[ 4] .y =    2.40;             M[50] . Xp[ 4] .x=  0.174;          
 M[50] . Xp[ 5] .y =    3.00;             M[50] . Xp[ 5] .x=  0.142;          
 M[50] . Xp[ 6] .y =    3.60;             M[50] . Xp[ 6] .x=  0.139;          
 M[50] . Xp[ 7] .y =    4.20;             M[50] . Xp[ 7] .x=  0.135;          
 M[50] . Xp[ 8] .y =    4.80;             M[50] . Xp[ 8] .x=  0.129;          
 M[50] . Xp[ 9] .y =    5.40;             M[50] . Xp[ 9] .x=  0.124;          
 M[50] . Xp[10] .y =    6.00;             M[50] . Xp[10] .x=  0.109;          
 M[50] . Xp[11] .y =    6.60;             M[50] . Xp[11] .x=  0.102;          
 M[50] . Xp[12] .y =    7.20;             M[50] . Xp[12] .x=  0.087;          
 M[50] . Xp[13] .y =    7.80;             M[50] . Xp[13] .x=  0.084;          
 M[50] . Xp[14] .y =    8.40;             M[50] . Xp[14] .x=  0.060;          
 M[50] . Xp[15] .y =    9.00;             M[50] . Xp[15] .x=  0.055;          
 M[50] . Xp[16] .y =    9.60;             M[50] . Xp[16] .x=  0.053;          
 M[50] . Xp[17] .y =   10.20;             M[50] . Xp[17] .x=  0.045;          
 M[50] . Xp[18] .y =   10.80;             M[50] . Xp[18] .x=  0.034;          
 M[50] . Xp[19] .y =   11.40;             M[50] . Xp[19] .x=  0.033;          
 M[50] . Xp[20] .y =   12.00;             M[50] . Xp[20] .x=  0.034;          
 M[50] . Xp[21] .y =   12.60;             M[50] . Xp[21] .x=  0.024;          
 M[50] . Xp[22] .y =   13.20;             M[50] . Xp[22] .x=  0.025;          
 M[50] . Xp[23] .y =   13.80;             M[50] . Xp[23] .x=  0.025;          
 M[50] . Xp[24] .y =   14.40;             M[50] . Xp[24] .x=  0.018;          
 M[50] . Xp[25] .y =   15.00;             M[50] . Xp[25] .x=  0.020;          
 M[50] . Xp[26] .y =   15.60;             M[50] . Xp[26] .x=  0.016;          
 M[50] . Xp[27] .y =   16.20;             M[50] . Xp[27] .x=  0.017;          
 M[50] . Xp[28] .y =   16.80;             M[50] . Xp[28] .x=  0.008;          
 M[50] . Xp[29] .y =   17.40;             M[50] . Xp[29] .x=  0.010;          
 M[50] . Xp[30] .y =   18.00;             M[50] . Xp[30] .x=  0.012;          
 M[50] . Xp[31] .y =   18.60;             M[50] . Xp[31] .x=  0.004;          
 M[50] . Xp[32] .y =   19.20;             M[50] . Xp[32] .x=  0.004;          
 M[50] . Xp[33] .y =   19.80;             M[50] . Xp[33] .x=  0.004;          
 M[50] . Xp[34] .y =   20.40;             M[50] . Xp[34] .x=  0.004;          
 M[50] . Xp[35] .y =   21.00;             M[50] . Xp[35] .x=  0.001;          
 M[50] . Xp[36] .y =   21.60;             M[50] . Xp[36] .x=  0.002;          
 M[50] . Xp[37] .y =   22.20;             M[50] . Xp[37] .x=  0.004;          
 M[50] . Xp[38] .y =   22.80;             M[50] . Xp[38] .x=  0.002;          
 M[50] . Xp[39] .y =   23.40;             M[50] . Xp[39] .x=  0.005;          
 M[50] . Xp[40] .y =   24.00;             M[50] . Xp[40] .x=  0.004;          
 M[50] . Xp[41] .y =   24.60;             M[50] . Xp[41] .x=  0.004;          
 M[50] . Xp[42] .y =   25.20;             M[50] . Xp[42] .x=  0.002;          
   M[50].pathlength();

 x[51] = 66.3;
  M[51].allocate(43);
  M[51].settype(SPLINE2D_LINEAR);                                                                            
                                                                          
 M[51] . Xp[ 0] .y =    0.00;             M[51] . Xp[ 0] .x=  0.222;          
 M[51] . Xp[ 1] .y =    0.60;             M[51] . Xp[ 1] .x=  0.191;          
 M[51] . Xp[ 2] .y =    1.20;             M[51] . Xp[ 2] .x=  0.186;          
 M[51] . Xp[ 3] .y =    1.80;             M[51] . Xp[ 3] .x=  0.172;          
 M[51] . Xp[ 4] .y =    2.40;             M[51] . Xp[ 4] .x=  0.153;          
 M[51] . Xp[ 5] .y =    3.00;             M[51] . Xp[ 5] .x=  0.136;          
 M[51] . Xp[ 6] .y =    3.60;             M[51] . Xp[ 6] .x=  0.126;          
 M[51] . Xp[ 7] .y =    4.20;             M[51] . Xp[ 7] .x=  0.129;          
 M[51] . Xp[ 8] .y =    4.80;             M[51] . Xp[ 8] .x=  0.114;          
 M[51] . Xp[ 9] .y =    5.40;             M[51] . Xp[ 9] .x=  0.107;          
 M[51] . Xp[10] .y =    6.00;             M[51] . Xp[10] .x=  0.100;          
 M[51] . Xp[11] .y =    6.60;             M[51] . Xp[11] .x=  0.093;          
 M[51] . Xp[12] .y =    7.20;             M[51] . Xp[12] .x=  0.086;          
 M[51] . Xp[13] .y =    7.80;             M[51] . Xp[13] .x=  0.082;          
 M[51] . Xp[14] .y =    8.40;             M[51] . Xp[14] .x=  0.056;          
 M[51] . Xp[15] .y =    9.00;             M[51] . Xp[15] .x=  0.056;          
 M[51] . Xp[16] .y =    9.60;             M[51] . Xp[16] .x=  0.049;          
 M[51] . Xp[17] .y =   10.20;             M[51] . Xp[17] .x=  0.037;          
 M[51] . Xp[18] .y =   10.80;             M[51] . Xp[18] .x=  0.036;          
 M[51] . Xp[19] .y =   11.40;             M[51] . Xp[19] .x=  0.026;          
 M[51] . Xp[20] .y =   12.00;             M[51] . Xp[20] .x=  0.029;          
 M[51] . Xp[21] .y =   12.60;             M[51] . Xp[21] .x=  0.023;          
 M[51] . Xp[22] .y =   13.20;             M[51] . Xp[22] .x=  0.024;          
 M[51] . Xp[23] .y =   13.80;             M[51] . Xp[23] .x=  0.023;          
 M[51] . Xp[24] .y =   14.40;             M[51] . Xp[24] .x=  0.019;          
 M[51] . Xp[25] .y =   15.00;             M[51] . Xp[25] .x=  0.021;          
 M[51] . Xp[26] .y =   15.60;             M[51] . Xp[26] .x=  0.019;          
 M[51] . Xp[27] .y =   16.20;             M[51] . Xp[27] .x=  0.019;          
 M[51] . Xp[28] .y =   16.80;             M[51] . Xp[28] .x=  0.012;          
 M[51] . Xp[29] .y =   17.40;             M[51] . Xp[29] .x=  0.008;          
 M[51] . Xp[30] .y =   18.00;             M[51] . Xp[30] .x=  0.011;          
 M[51] . Xp[31] .y =   18.60;             M[51] . Xp[31] .x=  0.008;          
 M[51] . Xp[32] .y =   19.20;             M[51] . Xp[32] .x=  0.004;          
 M[51] . Xp[33] .y =   19.80;             M[51] . Xp[33] .x=  0.003;          
 M[51] . Xp[34] .y =   20.40;             M[51] . Xp[34] .x=  0.002;          
 M[51] . Xp[35] .y =   21.00;             M[51] . Xp[35] .x=  0.002;          
 M[51] . Xp[36] .y =   21.60;             M[51] . Xp[36] .x=  0.001;          
 M[51] . Xp[37] .y =   22.20;             M[51] . Xp[37] .x=  0.005;          
 M[51] . Xp[38] .y =   22.80;             M[51] . Xp[38] .x=  0.002;          
 M[51] . Xp[39] .y =   23.40;             M[51] . Xp[39] .x=  0.003;          
 M[51] . Xp[40] .y =   24.00;             M[51] . Xp[40] .x=  0.001;          
 M[51] . Xp[41] .y =   24.60;             M[51] . Xp[41] .x=  0.004;          
 M[51] . Xp[42] .y =   25.20;             M[51] . Xp[42] .x=  0.004;          
 M[51].pathlength();
//station 53
 x[52] = 67.0;
  M[52].allocate(43);
  M[52].settype(SPLINE2D_LINEAR);                                                                             
                                                                          
 M[52] . Xp[ 0] .y =    0.00;             M[52] . Xp[ 0] .x=  0.182;          
 M[52] . Xp[ 1] .y =    0.60;             M[52] . Xp[ 1] .x=  0.182;          
 M[52] . Xp[ 2] .y =    1.20;             M[52] . Xp[ 2] .x=  0.173;          
 M[52] . Xp[ 3] .y =    1.80;             M[52] . Xp[ 3] .x=  0.177;          
 M[52] . Xp[ 4] .y =    2.40;             M[52] . Xp[ 4] .x=  0.142;          
 M[52] . Xp[ 5] .y =    3.00;             M[52] . Xp[ 5] .x=  0.147;          
 M[52] . Xp[ 6] .y =    3.60;             M[52] . Xp[ 6] .x=  0.128;          
 M[52] . Xp[ 7] .y =    4.20;             M[52] . Xp[ 7] .x=  0.123;          
 M[52] . Xp[ 8] .y =    4.80;             M[52] . Xp[ 8] .x=  0.114;          
 M[52] . Xp[ 9] .y =    5.40;             M[52] . Xp[ 9] .x=  0.094;          
 M[52] . Xp[10] .y =    6.00;             M[52] . Xp[10] .x=  0.099;          
 M[52] . Xp[11] .y =    6.60;             M[52] . Xp[11] .x=  0.086;          
 M[52] . Xp[12] .y =    7.20;             M[52] . Xp[12] .x=  0.082;          
 M[52] . Xp[13] .y =    7.80;             M[52] . Xp[13] .x=  0.071;          
 M[52] . Xp[14] .y =    8.40;             M[52] . Xp[14] .x=  0.055;          
 M[52] . Xp[15] .y =    9.00;             M[52] . Xp[15] .x=  0.055;          
 M[52] . Xp[16] .y =    9.60;             M[52] . Xp[16] .x=  0.047;          
 M[52] . Xp[17] .y =   10.20;             M[52] . Xp[17] .x=  0.043;          
 M[52] . Xp[18] .y =   10.80;             M[52] . Xp[18] .x=  0.037;          
 M[52] . Xp[19] .y =   11.40;             M[52] . Xp[19] .x=  0.030;          
 M[52] . Xp[20] .y =   12.00;             M[52] . Xp[20] .x=  0.027;          
 M[52] . Xp[21] .y =   12.60;             M[52] . Xp[21] .x=  0.028;          
 M[52] . Xp[22] .y =   13.20;             M[52] . Xp[22] .x=  0.025;          
 M[52] . Xp[23] .y =   13.80;             M[52] . Xp[23] .x=  0.027;          
 M[52] . Xp[24] .y =   14.40;             M[52] . Xp[24] .x=  0.017;          
 M[52] . Xp[25] .y =   15.00;             M[52] . Xp[25] .x=  0.018;          
 M[52] . Xp[26] .y =   15.60;             M[52] . Xp[26] .x=  0.017;          
 M[52] . Xp[27] .y =   16.20;             M[52] . Xp[27] .x=  0.016;          
 M[52] . Xp[28] .y =   16.80;             M[52] . Xp[28] .x=  0.016;          
 M[52] . Xp[29] .y =   17.40;             M[52] . Xp[29] .x=  0.009;          
 M[52] . Xp[30] .y =   18.00;             M[52] . Xp[30] .x=  0.010;          
 M[52] . Xp[31] .y =   18.60;             M[52] . Xp[31] .x=  0.007;          
 M[52] . Xp[32] .y =   19.20;             M[52] . Xp[32] .x=  0.001;          
 M[52] . Xp[33] .y =   19.80;             M[52] . Xp[33] .x=  0.006;          
 M[52] . Xp[34] .y =   20.40;             M[52] . Xp[34] .x=  0.005;          
 M[52] . Xp[35] .y =   21.00;             M[52] . Xp[35] .x=  0.002;          
 M[52] . Xp[36] .y =   21.60;             M[52] . Xp[36] .x=  0.005;          
 M[52] . Xp[37] .y =   22.20;             M[52] . Xp[37] .x=  0.005;          
 M[52] . Xp[38] .y =   22.80;             M[52] . Xp[38] .x=  0.003;          
 M[52] . Xp[39] .y =   23.40;             M[52] . Xp[39] .x=  0.009;          
 M[52] . Xp[40] .y =   24.00;             M[52] . Xp[40] .x=  0.008;          
 M[52] . Xp[41] .y =   24.60;             M[52] . Xp[41] .x=  0.002;          
 M[52] . Xp[42] .y =   25.20;             M[52] . Xp[42] .x=  0.002;          
 M[52].pathlength();
//station 54
 x[53] = 67.8;
  M[53].allocate(43);
  M[53].settype(SPLINE2D_LINEAR);                                                                            
                                                                         
 M[53] . Xp[ 0] .y =    0.00;             M[53] . Xp[ 0] .x=  0.152;          
 M[53] . Xp[ 1] .y =    0.60;             M[53] . Xp[ 1] .x=  0.161;          
 M[53] . Xp[ 2] .y =    1.20;             M[53] . Xp[ 2] .x=  0.168;          
 M[53] . Xp[ 3] .y =    1.80;             M[53] . Xp[ 3] .x=  0.155;          
 M[53] . Xp[ 4] .y =    2.40;             M[53] . Xp[ 4] .x=  0.148;          
 M[53] . Xp[ 5] .y =    3.00;             M[53] . Xp[ 5] .x=  0.137;          
 M[53] . Xp[ 6] .y =    3.60;             M[53] . Xp[ 6] .x=  0.132;          
 M[53] . Xp[ 7] .y =    4.20;             M[53] . Xp[ 7] .x=  0.108;          
 M[53] . Xp[ 8] .y =    4.80;             M[53] . Xp[ 8] .x=  0.108;          
 M[53] . Xp[ 9] .y =    5.40;             M[53] . Xp[ 9] .x=  0.089;          
 M[53] . Xp[10] .y =    6.00;             M[53] . Xp[10] .x=  0.096;          
 M[53] . Xp[11] .y =    6.60;             M[53] . Xp[11] .x=  0.073;          
 M[53] . Xp[12] .y =    7.20;             M[53] . Xp[12] .x=  0.078;          
 M[53] . Xp[13] .y =    7.80;             M[53] . Xp[13] .x=  0.071;          
 M[53] . Xp[14] .y =    8.40;             M[53] . Xp[14] .x=  0.056;          
 M[53] . Xp[15] .y =    9.00;             M[53] . Xp[15] .x=  0.041;          
 M[53] . Xp[16] .y =    9.60;             M[53] . Xp[16] .x=  0.053;          
 M[53] . Xp[17] .y =   10.20;             M[53] . Xp[17] .x=  0.044;          
 M[53] . Xp[18] .y =   10.80;             M[53] . Xp[18] .x=  0.023;          
 M[53] . Xp[19] .y =   11.40;             M[53] . Xp[19] .x=  0.037;          
 M[53] . Xp[20] .y =   12.00;             M[53] . Xp[20] .x=  0.025;          
 M[53] . Xp[21] .y =   12.60;             M[53] . Xp[21] .x=  0.027;          
 M[53] . Xp[22] .y =   13.20;             M[53] . Xp[22] .x=  0.020;          
 M[53] . Xp[23] .y =   13.80;             M[53] . Xp[23] .x=  0.026;          
 M[53] . Xp[24] .y =   14.40;             M[53] . Xp[24] .x=  0.022;          
 M[53] . Xp[25] .y =   15.00;             M[53] . Xp[25] .x=  0.012;          
 M[53] . Xp[26] .y =   15.60;             M[53] . Xp[26] .x=  0.011;          
 M[53] . Xp[27] .y =   16.20;             M[53] . Xp[27] .x=  0.009;          
 M[53] . Xp[28] .y =   16.80;             M[53] . Xp[28] .x=  0.009;          
 M[53] . Xp[29] .y =   17.40;             M[53] . Xp[29] .x=  0.006;          
 M[53] . Xp[30] .y =   18.00;             M[53] . Xp[30] .x=  0.007;          
 M[53] . Xp[31] .y =   18.60;             M[53] . Xp[31] .x=  0.011;          
 M[53] . Xp[32] .y =   19.20;             M[53] . Xp[32] .x=  0.006;          
 M[53] . Xp[33] .y =   19.80;             M[53] . Xp[33] .x=  0.003;          
 M[53] . Xp[34] .y =   20.40;             M[53] . Xp[34] .x=  0.002;          
 M[53] . Xp[35] .y =   21.00;             M[53] . Xp[35] .x=  0.001;          
 M[53] . Xp[36] .y =   21.60;             M[53] . Xp[36] .x=  0.005;          
 M[53] . Xp[37] .y =   22.20;             M[53] . Xp[37] .x=  0.005;          
 M[53] . Xp[38] .y =   22.80;             M[53] . Xp[38] .x=  0.003;          
 M[53] . Xp[39] .y =   23.40;             M[53] . Xp[39] .x=  0.004;          
 M[53] . Xp[40] .y =   24.00;             M[53] . Xp[40] .x=  0.005;          
 M[53] . Xp[41] .y =   24.60;             M[53] . Xp[41] .x=  0.003;          
 M[53] . Xp[42] .y =   25.20;             M[53] . Xp[42] .x=  0.005;          
   M[53].pathlength();
}

//! NonreactiveVelocityField::deallocate -- Deallocate memory for
//!                                            splines and arrays.
inline void NonreactiveScalarField::deallocate(void) {
  // Deallocate splines.
  for (int n = 0; n < Ns; n++) {
    M[n].deallocate();
  }
  // Deallocate arrays.
  delete []M; M = NULL;
  delete []x; x = NULL;
}

/*!
 * Class: ReactiveBluffBodyBurner
 *
 * @brief Normalized initial conditions for the bluff-body stabilised
 * turbulent non-reacting flows.
 *
 * \verbatim
 * Private member functions
 *      Ns -- Number of stations.
 *      U  -- Array of u-velocity splines at each station.
 *      V  -- Array of v-velocity splines at each station.
 *      x  -- Array of station x-coordinates (axial direction).
 * Public member functions
 *      allocate      -- Allocate memory and initialize arrays.
 *      deallocate    -- Deallocate memory for splines and arrays.
 *      interpolation -- Interpolation function.
 * \endverbatim
 */
class ReactiveVelocityField_Fuel_CH4H2 {
/* Coflow air + Fuel (CH4/H2)  */
 private:
  int       Ns; //!< Number of stations.
  Spline2D  *U; //!< Array of u-velocity splines at each station.
  Spline2D  *V; //!< Array of v-velocity splines at each station.
  double    *x; //!< Array of station x-coordinates (axial direction).
  
 public:
  //! Creation constructor.
  ReactiveVelocityField_Fuel_CH4H2(void) { allocate(); }

  //! Destructor.
  ~ReactiveVelocityField_Fuel_CH4H2(void) { deallocate(); }

  //! Allocate memory and initialize arrays.
  void allocate(void);

  //! Deallocate memory for splines and arrays.
  void deallocate(void);

  //! Interpolation function.
  Vector2D interpolation(const Vector2D &Xt);

};

//! ReactiveVelocityField::allocate -- Allocate memory and
//!                                          initialize arrays.
inline void ReactiveVelocityField_Fuel_CH4H2::allocate(void) {
/* Velocity of coflowing air is 40 m/s. */
/* Velocity of coflowing fuel (in this case "CH4/H2" ) is 118 m/s . */
  // Set the number of stations.
  Ns = 6 ;
  // Allocate the arrays.
  U = new Spline2D[Ns];
  V = new Spline2D[Ns];
  x = new double[Ns];
 
  // STATION 1
  x[0] = 0.1;
  U[0].allocate(21);                              V[0].allocate(21);
  U[0].settype(SPLINE2D_LINEAR);               V[0].settype(SPLINE2D_LINEAR);
  U[0].Xp[0].x=	142.41;	U[0].Xp[0].y=	0;	V[0].Xp[0].x=	1.01;	V[0].Xp[0].y=	0;
  U[0].Xp[1].x=	142.95;	U[0].Xp[1].y=	0.02;	V[0].Xp[1].x=	1.08;	V[0].Xp[1].y=	0.02;
  U[0].Xp[2].x=	141.34;	U[0].Xp[2].y=	0.06;	V[0].Xp[2].x=	1.66;	V[0].Xp[2].y=	0.06;
  U[0].Xp[3].x=	49.34;	U[0].Xp[3].y=	0.08;	V[0].Xp[3].x=	4.68;	V[0].Xp[3].y=	0.08;
  U[0].Xp[4].x=	10.51;	U[0].Xp[4].y=	0.12;	V[0].Xp[4].x=	-1.92;	V[0].Xp[4].y=	0.12;
  U[0].Xp[5].x=	-1.03;	U[0].Xp[5].y=	0.16;	V[0].Xp[5].x=	-1.36;	V[0].Xp[5].y=	0.16;
  U[0].Xp[6].x=	-1.53;	U[0].Xp[6].y=	0.2;	V[0].Xp[6].x=-1.35;	V[0].Xp[6].y=	0.2;
  U[0].Xp[7].x=-1.16;	U[0].Xp[7].y=	0.24;	V[0].Xp[7].x=-0.72;	V[0].Xp[7].y=	0.24;
  U[0].Xp[8].x=0.6;	U[0].Xp[8].y=	0.32;	V[0].Xp[8].x=0.17;	V[0].Xp[8].y=	0.32;
  U[0].Xp[9].x=0.63;	U[0].Xp[9].y=	0.4;	V[0].Xp[9].x=1.17;	V[0].Xp[9].y=	0.4;
  U[0].Xp[10].x=1.27;	U[0].Xp[10].y=	0.48;	V[0].Xp[10].x=2.68;	V[0].Xp[10].y=	0.48;
  U[0].Xp[11].x=1.75;	U[0].Xp[11].y=	0.56;	V[0].Xp[11].x=3.44;	V[0].Xp[11].y=	0.56;
  U[0].Xp[12].x=2.25;	U[0].Xp[12].y=	0.64;	V[0].Xp[12].x=3.63;	V[0].Xp[12].y=	0.64;
  U[0].Xp[13].x=2.33;	U[0].Xp[13].y=	0.72;	V[0].Xp[13].x=3.8;	V[0].Xp[13].y=	0.72;
  U[0].Xp[14].x=2.74;	U[0].Xp[14].y=	0.8;	V[0].Xp[14].x=3.82;	V[0].Xp[14].y=	0.8;
  U[0].Xp[15].x=2.87;	U[0].Xp[15].y=	0.92;	V[0].Xp[15].x=2.99;	V[0].Xp[15].y=	0.92;
  U[0].Xp[16].x=2.82;	U[0].Xp[16].y=	0.96;	V[0].Xp[16].x=2.16;	V[0].Xp[16].y=	0.96;
  U[0].Xp[17].x=4.85;	U[0].Xp[17].y=	1;	V[0].Xp[17].x=0.84;	V[0].Xp[17].y=	1;
  U[0].Xp[18].x=6.81;	U[0].Xp[18].y=	1.02;	V[0].Xp[18].x=0.72;	V[0].Xp[18].y=	1.02;
  U[0].Xp[19].x=36.44;	U[0].Xp[19].y=	1.08;	V[0].Xp[19].x=-0.86;	V[0].Xp[19].y=	1.08;
  U[0].Xp[20].x=38.72;	U[0].Xp[20].y=	1.2;	V[0].Xp[20].x=-1.37;	V[0].Xp[20].y=	1.2;
  U[0].pathlength(); V[0].pathlength();

  // STATION 2
  x[1] = 0.26;
  U[1].allocate(22);                              V[1].allocate(22);
  U[1].settype(SPLINE2D_LINEAR);               V[1].settype(SPLINE2D_LINEAR);
 
  U[1].Xp[0].x=136.84;	U[1].Xp[0].y=0;	        V[1].Xp[0].x=0.53;	V[1].Xp[0].y=	0;
  U[1].Xp[1].x=134.48;	U[1].Xp[1].y=0.04;	V[1].Xp[1].x=0.53;	V[1].Xp[1].y=	0.04;
  U[1].Xp[2].x=49.55;	U[1].Xp[2].y=0.08;	V[1].Xp[2].x=0.84;	V[1].Xp[2].y=	0.08;
  U[1].Xp[3].x=17.13;	U[1].Xp[3].y=0.12;	V[1].Xp[3].x=-0.46;	V[1].Xp[3].y=	0.12;
  U[1].Xp[4].x=1.42;	U[1].Xp[4].y=0.16;	V[1].Xp[4].x=-1.15;	V[1].Xp[4].y=	0.16;
  U[1].Xp[5].x=-1.93;	U[1].Xp[5].y=0.2;	V[1].Xp[5].x=-2.11;	V[1].Xp[5].y=	0.2;
  U[1].Xp[6].x=-2.01;	U[1].Xp[6].y=0.24;	V[1].Xp[6].x=-1.91;	V[1].Xp[6].y=	0.24;
  U[1].Xp[7].x=-3.81;	U[1].Xp[7].y=0.28;	V[1].Xp[7].x=-0.73;	V[1].Xp[7].y=	0.28;
  U[1].Xp[8].x=-3.95;	U[1].Xp[8].y=0.32;	V[1].Xp[8].x=-0.48;	V[1].Xp[8].y=	0.32;
  U[1].Xp[9].x=-5;	U[1].Xp[9].y=0.4;	V[1].Xp[9].x=0.1;	V[1].Xp[9].y=	0.4;
  U[1].Xp[10].x=-5.48;	U[1].Xp[10].y=0.48;	V[1].Xp[10].x=-0.03;	V[1].Xp[10].y=	0.48;
  U[1].Xp[11].x=-4.4;	U[1].Xp[11].y=0.56;	V[1].Xp[11].x=1.29;	V[1].Xp[11].y=	0.56;
  U[1].Xp[12].x=-3.45;	U[1].Xp[12].y=0.6;	V[1].Xp[12].x=1.41;	V[1].Xp[12].y=	0.6;
  U[1].Xp[13].x=-2.52;	U[1].Xp[13].y=0.68;	V[1].Xp[13].x=0.92;	V[1].Xp[13].y=	0.68;
  U[1].Xp[14].x=-1.08;	U[1].Xp[14].y=0.76;	V[1].Xp[14].x=0.94;	V[1].Xp[14].y=	0.76;
  U[1].Xp[15].x=-0.53;	U[1].Xp[15].y=0.8;	V[1].Xp[15].x=0.84;	V[1].Xp[15].y=	0.8;
  U[1].Xp[16].x=0.59;	U[1].Xp[16].y=0.88;	V[1].Xp[16].x=0.79;	V[1].Xp[16].y=	0.88;
  U[1].Xp[17].x=1.3;	U[1].Xp[17].y=0.92;	V[1].Xp[17].x=0.69;	V[1].Xp[17].y=	0.92;
  U[1].Xp[18].x=5.47;	U[1].Xp[18].y=0.96;	V[1].Xp[18].x=0.01;	V[1].Xp[18].y=	0.96;
  U[1].Xp[19].x=13.98;	U[1].Xp[19].y=1;	V[1].Xp[19].x=-0.2;	V[1].Xp[19].y=	1;
  U[1].Xp[20].x=37.2;	U[1].Xp[20].y=1.08;	V[1].Xp[20].x=-1.37;	V[1].Xp[20].y=	1.08;
  U[1].Xp[21].x=38.54;	U[1].Xp[21].y=1.2;	V[1].Xp[21].x=-1.48;	V[1].Xp[21].y=	1.2;

  U[1].pathlength(); V[1].pathlength();

  // STATION 3
  x[2] = 0.6;
  U[2].allocate(21);                              V[2].allocate(21);
  U[2].settype(SPLINE2D_LINEAR);               V[2].settype(SPLINE2D_LINEAR);
  U[2].Xp[0].x=113.51;	U[2].Xp[0].y=	0.04;	V[2].Xp[0].x=1.42;	V[2].Xp[0].y=	0.04;
  U[2].Xp[1].x=96.08;	U[2].Xp[1].y=	0.08;	V[2].Xp[1].x=2.2;	V[2].Xp[1].y=	0.08;
  U[2].Xp[2].x=28.22;	U[2].Xp[2].y=	0.12;	V[2].Xp[2].x=2.72;	V[2].Xp[2].y=	0.12;
  U[2].Xp[3].x=9.04;	U[2].Xp[3].y=	0.16;	V[2].Xp[3].x=2.17;	V[2].Xp[3].y=	0.16;
  U[2].Xp[4].x=1.35;	U[2].Xp[4].y=	0.2;	V[2].Xp[4].x=2.19;	V[2].Xp[4].y=	0.2;
  U[2].Xp[5].x=-1.62;	U[2].Xp[5].y=	0.24;	V[2].Xp[5].x=2.53;	V[2].Xp[5].y=	0.24;
  U[2].Xp[6].x=-5.14;	U[2].Xp[6].y=	0.32;	V[2].Xp[6].x=0.68;	V[2].Xp[6].y=	0.32;
  U[2].Xp[7].x=-7.02;	U[2].Xp[7].y=	0.36;	V[2].Xp[7].x=-0.25;	V[2].Xp[7].y=	0.36;
  U[2].Xp[8].x=-9.76;	U[2].Xp[8].y=	0.4;	V[2].Xp[8].x=-0.09;	V[2].Xp[8].y=	0.4;
  U[2].Xp[9].x=-11.37;	U[2].Xp[9].y=	0.48;	V[2].Xp[9].x=-0.07;	V[2].Xp[9].y=	0.48;
  U[2].Xp[10].x=-11.34;	U[2].Xp[10].y=	0.52;	V[2].Xp[10].x=0.67;	V[2].Xp[10].y=	0.52;
  U[2].Xp[11].x=-9.86;	U[2].Xp[11].y=	0.6;	V[2].Xp[11].x=1.4;	V[2].Xp[11].y=	0.6;
  U[2].Xp[12].x=-8.36;	U[2].Xp[12].y=	0.68;	V[2].Xp[12].x=2.07;	V[2].Xp[12].y=	0.68;
  U[2].Xp[13].x=-7.27;	U[2].Xp[13].y=	0.72;	V[2].Xp[13].x=2.08;	V[2].Xp[13].y=	0.72;
  U[2].Xp[14].x=-2.5;	U[2].Xp[14].y=	0.8;	V[2].Xp[14].x=0.33;	V[2].Xp[14].y=	0.8;
  U[2].Xp[15].x=0.15;	U[2].Xp[15].y=	0.88;	V[2].Xp[15].x=-0.07;	V[2].Xp[15].y=	0.88;
  U[2].Xp[16].x=11.18;	U[2].Xp[16].y=	0.92;	V[2].Xp[16].x=-1.11;	V[2].Xp[16].y=	0.92;
  U[2].Xp[17].x=11.98;	U[2].Xp[17].y=	0.96;	V[2].Xp[17].x=-1.09;	V[2].Xp[17].y=	0.96;
  U[2].Xp[18].x=26.16;	U[2].Xp[18].y=	1;	V[2].Xp[18].x=-1.11;	V[2].Xp[18].y=	1;
  U[2].Xp[19].x=33.43;	U[2].Xp[19].y=	1.08;	V[2].Xp[19].x=-1.16;	V[2].Xp[19].y=	1.08;
  U[2].Xp[20].x=40.74;	U[2].Xp[20].y=	1.2;	V[2].Xp[20].x=-1.81;	V[2].Xp[20].y=	1.2;

  U[2].pathlength(); V[2].pathlength();

 // STATION 4
  x[3] = 0.9;
  U[3].allocate(22);                              V[3].allocate(22);
  U[3].settype(SPLINE2D_LINEAR);               V[3].settype(SPLINE2D_LINEAR);
  U[3].Xp[0].x=63.54;	U[3].Xp[0].y=	0.04;	V[3].Xp[0].x=2;	V[3].Xp[0].y=	0.04;
  U[3].Xp[1].x=63.15;	U[3].Xp[1].y=	0.08;	V[3].Xp[1].x=4.7;	V[3].Xp[1].y=	0.08;
  U[3].Xp[2].x=56.65;	U[3].Xp[2].y=	0.12;	V[3].Xp[2].x=3.9;	V[3].Xp[2].y=	0.12;
  U[3].Xp[3].x=25.85;	U[3].Xp[3].y=	0.24;	V[3].Xp[3].x=2.47;	V[3].Xp[3].y=	0.24;
  U[3].Xp[4].x=-4.89;	U[3].Xp[4].y=	0.32;	V[3].Xp[4].x=1.87;	V[3].Xp[4].y=	0.32;
  U[3].Xp[5].x=-8.44;	U[3].Xp[5].y=	0.4;	V[3].Xp[5].x=0.17;	V[3].Xp[5].y=	0.4;
  U[3].Xp[6].x=-10.4;	U[3].Xp[6].y=	0.44;	V[3].Xp[6].x=-0.79;	V[3].Xp[6].y=	0.44;
  U[3].Xp[7].x=-11.49;	U[3].Xp[7].y=	0.48;	V[3].Xp[7].x=-0.42;	V[3].Xp[7].y=	0.48;
  U[3].Xp[8].x=-12.21;	U[3].Xp[8].y=	0.52;	V[3].Xp[8].x=-0.74;	V[3].Xp[8].y=	0.52;
  U[3].Xp[9].x=-11.88;	U[3].Xp[9].y=	0.56;	V[3].Xp[9].x=-0.57;	V[3].Xp[9].y=	0.56;
  U[3].Xp[10].x=-10.97;	U[3].Xp[10].y=	0.6;	V[3].Xp[10].x=-0.18;	V[3].Xp[10].y=	0.6;
  U[3].Xp[11].x=-7.87;	U[3].Xp[11].y=	0.68;	V[3].Xp[11].x=1.04;	V[3].Xp[11].y=	0.68;
  U[3].Xp[12].x=-6.02;	U[3].Xp[12].y=	0.72;	V[3].Xp[12].x=1.01;	V[3].Xp[12].y=	0.72;
  U[3].Xp[13].x=-3.74;	U[3].Xp[13].y=	0.76;	V[3].Xp[13].x=-1.47;	V[3].Xp[13].y=	0.76;
  U[3].Xp[14].x=2.81;	U[3].Xp[14].y=	0.8;	V[3].Xp[14].x=-2.05;	V[3].Xp[14].y=	0.8;
  U[3].Xp[15].x=14.33;	U[3].Xp[15].y=	0.88;	V[3].Xp[15].x=-2.42;	V[3].Xp[15].y=	0.88;
  U[3].Xp[16].x=27.73;	U[3].Xp[16].y=	0.92;	V[3].Xp[16].x=-2.98;	V[3].Xp[16].y=	0.92;
  U[3].Xp[17].x=28.87;	U[3].Xp[17].y=	0.96;	V[3].Xp[17].x=-2.84;	V[3].Xp[17].y=	0.96;
  U[3].Xp[18].x=31.19;	U[3].Xp[18].y=	1;	V[3].Xp[18].x=-2.27;	V[3].Xp[18].y=	1;
  U[3].Xp[19].x=31.27;	U[3].Xp[19].y=	1;	V[3].Xp[19].x=-2.38;	V[3].Xp[19].y=	1;
  U[3].Xp[20].x=37;	U[3].Xp[20].y=	1.08;	V[3].Xp[20].x=-2.62;	V[3].Xp[20].y=	1.08;
  U[3].Xp[21].x=40.07;	U[3].Xp[21].y=	1.2;	V[3].Xp[21].x=-2.22;	V[3].Xp[21].y=	1.2;
  U[3].pathlength(); V[3].pathlength();
  // STATION 5
  x[4] = 1.3;
  U[4].allocate(24);                              V[4].allocate(24);
  U[4].settype(SPLINE2D_LINEAR);               V[4].settype(SPLINE2D_LINEAR);
  U[4].Xp[0].x=56.5;	U[4].Xp[0].y=	0;	V[4].Xp[0].x=0.42;	V[4].Xp[0].y=	0;
  U[4].Xp[1].x=58.9;	U[4].Xp[1].y=	0.04;	V[4].Xp[1].x=1.73;	V[4].Xp[1].y=	0.04;
  U[4].Xp[2].x=54.9;	U[4].Xp[2].y=	0.08;	V[4].Xp[2].x=2.09;	V[4].Xp[2].y=	0.08;
  U[4].Xp[3].x=49.86;	U[4].Xp[3].y=	0.12;	V[4].Xp[3].x=2.34;	V[4].Xp[3].y=	0.12;
  U[4].Xp[4].x=45.71;	U[4].Xp[4].y=	0.16;	V[4].Xp[4].x=2.54;	V[4].Xp[4].y=	0.16;
  U[4].Xp[5].x=35.68;	U[4].Xp[5].y=	0.2;	V[4].Xp[5].x=2.79;	V[4].Xp[5].y=	0.2;
  U[4].Xp[6].x=23.8;	U[4].Xp[6].y=	0.24;	V[4].Xp[6].x=2.93;	V[4].Xp[6].y=	0.24;
  U[4].Xp[7].x=9.18;	U[4].Xp[7].y=	0.28;	V[4].Xp[7].x=3.48;	V[4].Xp[7].y=	0.28;
  U[4].Xp[8].x=2.49;	U[4].Xp[8].y=	0.32;	V[4].Xp[8].x=2.78;	V[4].Xp[8].y=	0.32;
  U[4].Xp[9].x=-6.43;	U[4].Xp[9].y=	0.4;	V[4].Xp[9].x=3.25;	V[4].Xp[9].y=	0.4;
  U[4].Xp[10].x=-6.86;	U[4].Xp[10].y=	0.44;	V[4].Xp[10].x=2.84;	V[4].Xp[10].y=	0.44;
  U[4].Xp[11].x=-8.98;	U[4].Xp[11].y=	0.52;	V[4].Xp[11].x=2.25;	V[4].Xp[11].y=	0.52;
  U[4].Xp[12].x=-9.32;	U[4].Xp[12].y=	0.56;	V[4].Xp[12].x=2.91;	V[4].Xp[12].y=	0.56;
  U[4].Xp[13].x=-8.37;	U[4].Xp[13].y=	0.6;	V[4].Xp[13].x=2.2;	V[4].Xp[13].y=	0.6;
  U[4].Xp[14].x=-6.83;	U[4].Xp[14].y=	0.64;	V[4].Xp[14].x=2.44;	V[4].Xp[14].y=	0.64;
  U[4].Xp[15].x=-3.73;	U[4].Xp[15].y=	0.68;	V[4].Xp[15].x=-0.79;	V[4].Xp[15].y=	0.68;
  U[4].Xp[16].x=2.35;	U[4].Xp[16].y=	0.72;	V[4].Xp[16].x=-1.76;	V[4].Xp[16].y=	0.72;
  U[4].Xp[17].x=11.63;	U[4].Xp[17].y=	0.8;	V[4].Xp[17].x=-2.66;	V[4].Xp[17].y=	0.8;
  U[4].Xp[18].x=21.64;	U[4].Xp[18].y=	0.84;	V[4].Xp[18].x=-3.69;	V[4].Xp[18].y=	0.84;
  U[4].Xp[19].x=23.59;	U[4].Xp[19].y=	0.88;	V[4].Xp[19].x=-3.29;	V[4].Xp[19].y=	0.88;
  U[4].Xp[20].x=32.44;	U[4].Xp[20].y=	0.96;	V[4].Xp[20].x=-3.48;	V[4].Xp[20].y=	0.96;
  U[4].Xp[21].x=34.44;	U[4].Xp[21].y=	1;	V[4].Xp[21].x=-3.12;	V[4].Xp[21].y=	1;
  U[4].Xp[22].x=38.23;	U[4].Xp[22].y=	1.08;	V[4].Xp[22].x=-2.78;	V[4].Xp[22].y=	1.08;
  U[4].Xp[23].x=38.77;	U[4].Xp[23].y=	1.2;	V[4].Xp[23].x=-2.85;	V[4].Xp[23].y=	1.2;
  U[4].pathlength(); V[4].pathlength();
  // STATION 6
  x[5] = 1.8;
  U[5].allocate(20);                              V[5].allocate(20);
  U[5].settype(SPLINE2D_LINEAR);               V[5].settype(SPLINE2D_LINEAR);
  
  U[5].Xp[0].x=32.8;	U[5].Xp[0].y=	0;	V[5].Xp[0].x=-0.72;	V[5].Xp[0].y=	0;
  U[5].Xp[1].x=32.13;	U[5].Xp[1].y=	0.04;	V[5].Xp[1].x=0.02;	V[5].Xp[1].y=	0.04;
  U[5].Xp[2].x=31.88;	U[5].Xp[2].y=	0.08;	V[5].Xp[2].x=1.44;	V[5].Xp[2].y=	0.08;
  U[5].Xp[3].x=13.61;	U[5].Xp[3].y=	0.12;	V[5].Xp[3].x=2.01;	V[5].Xp[3].y=	0.12;
  U[5].Xp[4].x=11.31;	U[5].Xp[4].y=	0.16;	V[5].Xp[4].x=3.98;	V[5].Xp[4].y=	0.16;
  U[5].Xp[5].x=8.4;	U[5].Xp[5].y=	0.2;	V[5].Xp[5].x=4.05;	V[5].Xp[5].y=	0.2;
  U[5].Xp[6].x=2.75;	U[5].Xp[6].y=	0.24;	V[5].Xp[6].x=4.19;	V[5].Xp[6].y=	0.24;
  U[5].Xp[7].x=3.03;	U[5].Xp[7].y=	0.28;	V[5].Xp[7].x=3.78;	V[5].Xp[7].y=	0.28;
  U[5].Xp[8].x=2.85;	U[5].Xp[8].y=	0.32;	V[5].Xp[8].x=4.45;	V[5].Xp[8].y=	0.32;
  U[5].Xp[9].x=2.32;	U[5].Xp[9].y=	0.4;	V[5].Xp[9].x=4.5;	V[5].Xp[9].y=	0.4;
  U[5].Xp[10].x=2.93;	U[5].Xp[10].y=	0.48;	V[5].Xp[10].x=3.6;	V[5].Xp[10].y=	0.48;
  U[5].Xp[11].x=3.29;	U[5].Xp[11].y=	0.52;	V[5].Xp[11].x=2.66;	V[5].Xp[11].y=	0.52;
  U[5].Xp[12].x=3.24;	U[5].Xp[12].y=	0.56;	V[5].Xp[12].x=1.33;	V[5].Xp[12].y=	0.56;
  U[5].Xp[13].x=4.53;	U[5].Xp[13].y=	0.64;	V[5].Xp[13].x=0.3;	V[5].Xp[13].y=	0.64;
  U[5].Xp[14].x=11.81;	U[5].Xp[14].y=	0.72;	V[5].Xp[14].x=-2.75;	V[5].Xp[14].y=	0.72;
  U[5].Xp[15].x=13.4;	U[5].Xp[15].y=	0.8;	V[5].Xp[15].x=-2.57;	V[5].Xp[15].y=	0.8;
  U[5].Xp[16].x=19.52;	U[5].Xp[16].y=	0.88;	V[5].Xp[16].x=-2.43;	V[5].Xp[16].y=	0.88;
  U[5].Xp[17].x=31.21;	U[5].Xp[17].y=	0.92;	V[5].Xp[17].x=-2;	V[5].Xp[17].y=	0.92;
  U[5].Xp[18].x=38.01;	U[5].Xp[18].y=	1;	V[5].Xp[18].x=-2.77;	V[5].Xp[18].y=	1;
  U[5].Xp[19].x=39.1;	U[5].Xp[19].y=	1.12;	V[5].Xp[19].x=-2.72;	V[5].Xp[19].y=	1.12;
  U[5].pathlength(); V[5].pathlength();

 
}

//! ReactiveVelocityField::deallocate -- Deallocate memory for
//!                                            splines and arrays.
inline void ReactiveVelocityField_Fuel_CH4H2::deallocate(void) {
  // Deallocate splines.
  for (int n = 0; n < Ns; n++) {
    U[n].deallocate();
    V[n].deallocate();
  }
  // Deallocate arrays.
  delete []U; U = NULL;
  delete []V; V = NULL;
  delete []x; x = NULL;
}
/* Private member functions
 *      Ns -- Number of stations.
 *      M  -- Array of mixture fraction splines at each station.
 *      x  -- Array of station x-coordinates (axial direction).
 * Public member functions
 *      allocate      -- Allocate memory and initialize arrays.
 *      deallocate    -- Deallocate memory for splines and arrays.
 *      interpolation -- Interpolation function.
 * \endverbatim
 */
class ReactiveScalarField_Fuel_CH4H2 {
/* Coflow air + Fuel (CH4/H2)  */
/* Velocity of coflowing air is 40 m/s. */
/* Velocity of coflowing fuel (in this case "CH4/H2" ) is 118 m/s . */
 private:
  int       Ns; //!< Number of stations.
  Spline2D  *M; //!< Array of mixture fraction splines at each station.
  Spline2D  *T; //!< Array of mixture temperature splines at each station.
  Spline2D  *O2; //!< Array of species fraction splines at each station.
  Spline2D  *N2; //!< Array of species fraction splines at each station.
  Spline2D  *H2; //!< Array of species fraction splines at each station.
  Spline2D  *H2O; //!< Array of species fraction splines at each station.
  Spline2D  *CO; //!< Array of species fraction splines at each station.
  Spline2D  *CO2; //!< Array of species fraction splines at each station.
  Spline2D  *HC; //!< Array of species fraction splines at each station.
  Spline2D  *OH; //!< Array of species fraction splines at each station.
 
  double    *x; //!< Array of station x-coordinates (axial direction).
  
 public:
  //! Creation constructor.
  ReactiveScalarField_Fuel_CH4H2(void) { allocate(); }
  
  //! Destructor.
  ~ReactiveScalarField_Fuel_CH4H2(void) { deallocate(); }
  
  //! Allocate memory and initialize arrays.
  void allocate(void);
  
  //! Deallocate memory for splines and arrays.
  void deallocate(void);

  //! Interpolation function.
 
  Vector2D interpolation(const Vector2D &Xt, double []);

};

//! NonreactiveBluffBodeuBurner::allocate -- Allocate memory and
//!                                          initialize arrays.
inline void ReactiveScalarField_Fuel_CH4H2::allocate(void) {
  // Set the number of stations.
  Ns = 6;
  // Allocate the arrays.
  M = new Spline2D[Ns];
  T = new Spline2D[Ns];
  O2 = new Spline2D[Ns];
  N2 = new Spline2D[Ns];
  H2 = new Spline2D[Ns];
  H2O= new Spline2D[Ns];
  CO= new Spline2D[Ns];
  CO2 = new Spline2D[Ns];
  HC = new Spline2D[Ns];
  OH = new Spline2D[Ns];

 /*  C = new Spline2D *[9]; //nine species */
/*   for(int i; i<9; i++){ */
/*     C[i] = new Spline2D[Ns]; */
/*   } */

  x = new double[Ns];
  //Station 1
  x[0 ] = 0.3;
  T[0].allocate(9);
  T[0].settype(SPLINE2D_LINEAR);
  O2[0].allocate(9);
  O2[0].settype(SPLINE2D_LINEAR);
  N2[0].allocate(9);
  N2[0].settype(SPLINE2D_LINEAR);
  H2[0].allocate(9);
  H2[0].settype(SPLINE2D_LINEAR);
  H2O[0].allocate(9);
  H2O[0].settype(SPLINE2D_LINEAR);
  CO[0].allocate(9);
  CO[0].settype(SPLINE2D_LINEAR);
  CO2[0].allocate(9);
  CO2[0].settype(SPLINE2D_LINEAR);
  HC[0].allocate(9);
  HC[0].settype(SPLINE2D_LINEAR);
  OH[0].allocate(9);
  OH[0].settype(SPLINE2D_LINEAR);
  M[0].allocate(9);
  M[0].settype(SPLINE2D_LINEAR);
  
  T[0].Xp[0].x=	323;	T[0].Xp[0].y=	0.1;
  T[0].Xp[1].x=	365;	T[0].Xp[1].y=	1.1;
  T[0].Xp[2].x=	563;	T[0].Xp[2].y=	2.1;
  T[0].Xp[3].x=	1364;	T[0].Xp[3].y=	4.1;
  T[0].Xp[4].x=	1620;	T[0].Xp[4].y=	9.1;
  T[0].Xp[5].x=	1641;	T[0].Xp[5].y=	15.1;
  T[0].Xp[6].x=	1613;	T[0].Xp[6].y=	23;
  T[0].Xp[7].x=	1470;	T[0].Xp[7].y=	24;
  T[0].Xp[8].x=	797;	T[0].Xp[8].y=	25;
  //O2 fraction
  O2[0].Xp[0].x=2.04;	O2[0].Xp[0].y=0.1;
  O2[0].Xp[1].x=5.89;	O2[0].Xp[1].y=1.1;
  O2[0].Xp[2].x=20.48;	O2[0].Xp[2].y=2.1;
  O2[0].Xp[3].x=55.75;	O2[0].Xp[3].y=4.1;
  O2[0].Xp[4].x=62.4;	O2[0].Xp[4].y=9.1;
  O2[0].Xp[5].x=63.25;	O2[0].Xp[5].y=15.1;
  O2[0].Xp[6].x=63.92;	O2[0].Xp[6].y=23;
  O2[0].Xp[7].x=65.34;	O2[0].Xp[7].y=24;
  O2[0].Xp[8].x=73.2;	O2[0].Xp[8].y=25;

  N2[0].Xp[0].x=9.2595 ;N2[0].Xp[0].y=0.1;
  N2[0].Xp[1].x=9.0388;	N2[0].Xp[1].y=1.1;
  N2[0].Xp[2].x=7.0184;	N2[0].Xp[2].y=2.1;
  N2[0].Xp[3].x=2.3256;	N2[0].Xp[3].y=4.1;
  N2[0].Xp[4].x=1.245;	N2[0].Xp[4].y=9.1;
  N2[0].Xp[5].x=1.0738;	N2[0].Xp[5].y=15.1;
  N2[0].Xp[6].x=0.8403;	N2[0].Xp[6].y=23;
  N2[0].Xp[7].x=0.6679;	N2[0].Xp[7].y=24;
  N2[0].Xp[8].x=0.2964;	N2[0].Xp[8].y=25;
  CO[0].Xp[0].x=0.231;	CO[0].Xp[0].y=0.1;
  CO[0].Xp[1].x=0.84;	CO[0].Xp[1].y=1.1;
  CO[0].Xp[2].x=2.55;	CO[0].Xp[2].y=2.1;
  CO[0].Xp[3].x=5.666;	CO[0].Xp[3].y=4.1;
  CO[0].Xp[4].x=6.633;	CO[0].Xp[4].y=9.1;
  CO[0].Xp[5].x=7.155;	CO[0].Xp[5].y=15.1;
  CO[0].Xp[6].x=7.516;	CO[0].Xp[6].y=23;
  CO[0].Xp[7].x=6.136;	CO[0].Xp[7].y=24;
  CO[0].Xp[8].x=2.089;	CO[0].Xp[8].y=25;


  H2[0].Xp[0].x=0.353;	H2[0].Xp[0].y=0.1;
  H2[0].Xp[1].x=1.093;	H2[0].Xp[1].y=1.1;
  H2[0].Xp[2].x=4.658;	H2[0].Xp[2].y=2.1;
  H2[0].Xp[3].x=12.846;	H2[0].Xp[3].y=4.1;
  H2[0].Xp[4].x=14.803;	H2[0].Xp[4].y=9.1;
  H2[0].Xp[5].x=15.302;	H2[0].Xp[5].y=15.1;
  H2[0].Xp[6].x=15.002;	H2[0].Xp[6].y=23;
  H2[0].Xp[7].x=11.748;	H2[0].Xp[7].y=24;
  H2[0].Xp[8].x=4.806;	H2[0].Xp[8].y=25;
 
  H2O[0].Xp[0].x=0.081;	H2O[0].Xp[0].y=	0.1;
  H2O[0].Xp[1].x=0.402;	H2O[0].Xp[1].y=	1.1;
  H2O[0].Xp[2].x=1.677;	H2O[0].Xp[2].y=	2.1;
  H2O[0].Xp[3].x=4.583;	H2O[0].Xp[3].y=	4.1;
  H2O[0].Xp[4].x=5.227;	H2O[0].Xp[4].y=	9.1;
  H2O[0].Xp[5].x=5.195;	H2O[0].Xp[5].y=	15.1;
  H2O[0].Xp[6].x=5.327;	H2O[0].Xp[6].y=	23;
  H2O[0].Xp[7].x=4.307;	H2O[0].Xp[7].y=	24;
  H2O[0].Xp[8].x=1.601;	H2O[0].Xp[8].y=	25;
  CO2[0].Xp[0].x=86.616;CO2[0].Xp[0].y=0.1;
  CO2[0].Xp[1].x=81.633;CO2[0].Xp[1].y=1.1;
  CO2[0].Xp[2].x=62.46;	CO2[0].Xp[2].y=2.1;
  CO2[0].Xp[3].x=16.4;	CO2[0].Xp[3].y=4.1;
  CO2[0].Xp[4].x=6.617;	CO2[0].Xp[4].y=9.1;
  CO2[0].Xp[5].x=4.699;	CO2[0].Xp[5].y=15.1;
  CO2[0].Xp[6].x=3.309;	CO2[0].Xp[6].y=23;
  CO2[0].Xp[7].x=2.781;	CO2[0].Xp[7].y=24;
  CO2[0].Xp[8].x=1.399;	CO2[0].Xp[8].y=25;
  HC[0].Xp[0].x=0;	HC[0].Xp[0].y=0.1;
  HC[0].Xp[1].x=0;	HC[0].Xp[1].y=1.1;
  HC[0].Xp[2].x=0;	HC[0].Xp[2].y=2.1;
  HC[0].Xp[3].x=0.0057;	HC[0].Xp[3].y=4.1;
  HC[0].Xp[4].x=0.0165;	HC[0].Xp[4].y=9.1;
  HC[0].Xp[5].x=0.0208;	HC[0].Xp[5].y=15.1;
  HC[0].Xp[6].x=0.0229;	HC[0].Xp[6].y=23;
  HC[0].Xp[7].x=0.0392;	HC[0].Xp[7].y=24;
  HC[0].Xp[8].x=0.0261;	HC[0].Xp[8].y=25;
  OH[0].Xp[0].x=0.0018;	OH[0].Xp[0].y=0.1;
  OH[0].Xp[1].x=0.013;	OH[0].Xp[1].y=1.1;
  OH[0].Xp[2].x=0.0526;	OH[0].Xp[2].y=2.1;
  OH[0].Xp[3].x=0.2781;	OH[0].Xp[3].y=4.1;
  OH[0].Xp[4].x=0.4119;	OH[0].Xp[4].y=9.1;
  OH[0].Xp[5].x=0.4292;	OH[0].Xp[5].y=15.1;
  OH[0].Xp[6].x=0.4137;	OH[0].Xp[6].y=23;
  OH[0].Xp[7].x=0.3561;	OH[0].Xp[7].y=24;
  OH[0].Xp[8].x=0.1037;	OH[0].Xp[8].y=25;
  M[0].Xp[0].x=1.129;	M[0].Xp[0].y=0.1;
  M[0].Xp[1].x=1.179;	M[0].Xp[1].y=1.1;
  M[0].Xp[2].x=1.272;	M[0].Xp[2].y=2.1;
  M[0].Xp[3].x=1.175;	M[0].Xp[3].y=4.1;
  M[0].Xp[4].x=1.084;	M[0].Xp[4].y=9.1;
  M[0].Xp[5].x=1.046;	M[0].Xp[5].y=15.1;
  M[0].Xp[6].x=1.022;	M[0].Xp[6].y=23;
  M[0].Xp[7].x=1.047;	M[0].Xp[7].y=24;
  M[0].Xp[8].x=1.025;	M[0].Xp[8].y=25;
  T[0].pathlength();
  O2[0].pathlength();
  N2[0].pathlength();
  H2[0].pathlength();
  H2O[0].pathlength();
  CO[0].pathlength();
  CO2[0].pathlength();
  HC[0].pathlength();
  OH[0].pathlength();
  M[0].pathlength();

  //Station 2
  x[1 ] = 0.6;
  T[1].allocate(13);
  T[1].settype(SPLINE2D_LINEAR);
  O2[1].allocate(13);
  O2[1].settype(SPLINE2D_LINEAR);
  N2[1].allocate(13);
  N2[1].settype(SPLINE2D_LINEAR);
  H2[1].allocate(13);
  H2[1].settype(SPLINE2D_LINEAR);
  H2O[1].allocate(13);
  H2O[1].settype(SPLINE2D_LINEAR);
  CO[1].allocate(13);
  CO[1].settype(SPLINE2D_LINEAR);
  CO2[1].allocate(13);
  CO2[1].settype(SPLINE2D_LINEAR);
  HC[1].allocate(13);
  HC[1].settype(SPLINE2D_LINEAR);
  OH[1].allocate(13);
  OH[1].settype(SPLINE2D_LINEAR);
  M[1].allocate(13);
  M[1].settype(SPLINE2D_LINEAR);
  
  T[1].Xp[0].x=471;	T[1].Xp[0].y=0.1;
  T[1].Xp[1].x=514;	T[1].Xp[1].y=1.1;
  T[1].Xp[2].x=637;	T[1].Xp[2].y=2.1;
  T[1].Xp[3].x=808;	T[1].Xp[3].y=3.1;
  T[1].Xp[4].x=1126;	T[1].Xp[4].y=5.1;
  T[1].Xp[5].x=1502;	T[1].Xp[5].y=9.1;
  T[1].Xp[6].x=1653;	T[1].Xp[6].y=13.1;
  T[1].Xp[7].x=1756;	T[1].Xp[7].y=17;
  T[1].Xp[8].x=1819;	T[1].Xp[8].y=21;
  T[1].Xp[9].x=1561;	T[1].Xp[9].y=23;
  T[1].Xp[10].x=1056;	T[1].Xp[10].y=24;
  T[1].Xp[11].x=618;	T[1].Xp[11].y=25.1;
  T[1].Xp[12].x=414;	T[1].Xp[12].y=26;
  O2[1].Xp[0].x=13.49;	O2[1].Xp[0].y=0.1;
  O2[1].Xp[1].x=16.38;	O2[1].Xp[1].y=1.1;
  O2[1].Xp[2].x=24.03;	O2[1].Xp[2].y=2.1;
  O2[1].Xp[3].x=32.94;	O2[1].Xp[3].y=3.1;
  O2[1].Xp[4].x=47.44;	O2[1].Xp[4].y=5.1;
  O2[1].Xp[5].x=60.12;	O2[1].Xp[5].y=9.1;
  O2[1].Xp[6].x=61.42;	O2[1].Xp[6].y=13.1;
  O2[1].Xp[7].x=61.46;	O2[1].Xp[7].y=17;
  O2[1].Xp[8].x=61.3;	O2[1].Xp[8].y=21;
  O2[1].Xp[9].x=69.23;	O2[1].Xp[9].y=23;
  O2[1].Xp[10].x=73.37;	O2[1].Xp[10].y=24;
  O2[1].Xp[11].x=75.46;	O2[1].Xp[11].y=25.1;
  O2[1].Xp[12].x=75.97;	O2[1].Xp[12].y=26;

  N2[1].Xp[0].x=7.9166;	N2[1].Xp[0].y=0.1;
  N2[1].Xp[1].x=7.5744;	N2[1].Xp[1].y=1.1;
  N2[1].Xp[2].x=6.4644;	N2[1].Xp[2].y=2.1;
  N2[1].Xp[3].x=5.161;	N2[1].Xp[3].y=3.1;
  N2[1].Xp[4].x=3.181;	N2[1].Xp[4].y=5.1;
  N2[1].Xp[5].x=1.4646;	N2[1].Xp[5].y=9.1;
  N2[1].Xp[6].x=0.9603;	N2[1].Xp[6].y=13.1;
  N2[1].Xp[7].x=0.6598;	N2[1].Xp[7].y=17;
  N2[1].Xp[8].x=0.4609;	N2[1].Xp[8].y=21;
  N2[1].Xp[9].x=0.2489;	N2[1].Xp[9].y=23;
  N2[1].Xp[10].x=0.1552;N2[1].Xp[10].y=24;
  N2[1].Xp[11].x=0.0998;N2[1].Xp[11].y=25.1;
  N2[1].Xp[12].x=0.0463;N2[1].Xp[12].y=26;
  H2[1].Xp[0].x=3.052;	H2[1].Xp[0].y=0.1;
  H2[1].Xp[1].x=3.721;	H2[1].Xp[1].y=1.1;
  H2[1].Xp[2].x=5.529;	H2[1].Xp[2].y=2.1;
  H2[1].Xp[3].x=7.629;	H2[1].Xp[3].y=3.1;
  H2[1].Xp[4].x=10.788;	H2[1].Xp[4].y=5.1;
  H2[1].Xp[5].x=13.633;	H2[1].Xp[5].y=9.1;
  H2[1].Xp[6].x=14.422;	H2[1].Xp[6].y=13.1;
  H2[1].Xp[7].x=14.345;	H2[1].Xp[7].y=17;
  H2[1].Xp[8].x=13.623;	H2[1].Xp[8].y=21;
  H2[1].Xp[9].x=9.998;	H2[1].Xp[9].y=23;
  H2[1].Xp[10].x=6.162;	H2[1].Xp[10].y=24;
  H2[1].Xp[11].x=3.155;	H2[1].Xp[11].y=25.1;
  H2[1].Xp[12].x=1.762;	H2[1].Xp[12].y=26;

  H2O[1].Xp[0].x=1.091;	H2O[1].Xp[0].y=0.1;
  H2O[1].Xp[1].x=1.317;	H2O[1].Xp[1].y=1.1;
  H2O[1].Xp[2].x=1.937;	H2O[1].Xp[2].y=2.1;
  H2O[1].Xp[3].x=2.577;	H2O[1].Xp[3].y=3.1;
  H2O[1].Xp[4].x=3.656;	H2O[1].Xp[4].y=5.1;
  H2O[1].Xp[5].x=4.527;	H2O[1].Xp[5].y=9.1;
  H2O[1].Xp[6].x=4.835;	H2O[1].Xp[6].y=13.1;
  H2O[1].Xp[7].x=5.054;	H2O[1].Xp[7].y=17;
  H2O[1].Xp[8].x=4.865;	H2O[1].Xp[8].y=21;
  H2O[1].Xp[9].x=2.696;	H2O[1].Xp[9].y=23;
  H2O[1].Xp[10].x=1.461;	H2O[1].Xp[10].y=24;
  H2O[1].Xp[11].x=0.735;	H2O[1].Xp[11].y=25.1;
  H2O[1].Xp[12].x=0.312;	H2O[1].Xp[12].y=26;
  CO[1].Xp[0].x=1.521;	CO[1].Xp[0].y=0.1;
  CO[1].Xp[1].x=1.982;	CO[1].Xp[1].y=1.1;
  CO[1].Xp[2].x=2.84;	CO[1].Xp[2].y=2.1;
  CO[1].Xp[3].x=3.657;	CO[1].Xp[3].y=3.1;
  CO[1].Xp[4].x=4.787;	CO[1].Xp[4].y=5.1;
  CO[1].Xp[5].x=5.902;	CO[1].Xp[5].y=9.1;
  CO[1].Xp[6].x=6.758;	CO[1].Xp[6].y=13.1;
  CO[1].Xp[7].x=7.305;	CO[1].Xp[7].y=17;
  CO[1].Xp[8].x=7.436;	CO[1].Xp[8].y=21;
  CO[1].Xp[9].x=5.834;	CO[1].Xp[9].y=23;
  CO[1].Xp[10].x=3.212;	CO[1].Xp[10].y=24;
  CO[1].Xp[11].x=1.079;	CO[1].Xp[11].y=25.1;
  CO[1].Xp[12].x=0.306;	CO[1].Xp[12].y=26;
  CO2[1].Xp[0].x=71.701;CO2[1].Xp[0].y=0.1;
  CO2[1].Xp[1].x=67.867;CO2[1].Xp[1].y=1.1;
  CO2[1].Xp[2].x=57.891;CO2[1].Xp[2].y=2.1;
  CO2[1].Xp[3].x=46.326;CO2[1].Xp[3].y=3.1;
  CO2[1].Xp[4].x=27.731;CO2[1].Xp[4].y=5.1;
  CO2[1].Xp[5].x=10.54;	CO2[1].Xp[5].y=9.1;
  CO2[1].Xp[6].x=5.466;	CO2[1].Xp[6].y=13.1;
  CO2[1].Xp[7].x=2.293;	CO2[1].Xp[7].y=17;
  CO2[1].Xp[8].x=0.881;	CO2[1].Xp[8].y=21;
  CO2[1].Xp[9].x=0.686;	CO2[1].Xp[9].y=23;
  CO2[1].Xp[10].x=0.667;CO2[1].Xp[10].y=24;
  CO2[1].Xp[11].x=0.431;CO2[1].Xp[11].y=25.1;
  CO2[1].Xp[12].x=0.21;	CO2[1].Xp[12].y=26;
  HC[1].Xp[0].x=0;      HC[1].Xp[0].y=0.1;
  HC[1].Xp[1].x=0;	HC[1].Xp[1].y=1.1;
  HC[1].Xp[2].x=0;	HC[1].Xp[2].y=2.1;
  HC[1].Xp[3].x=0;	HC[1].Xp[3].y=3.1;
  HC[1].Xp[4].x=0.0005;	HC[1].Xp[4].y=5.1;
  HC[1].Xp[5].x=0.009;	HC[1].Xp[5].y=9.1;
  HC[1].Xp[6].x=0.0159;	HC[1].Xp[6].y=13.1;
  HC[1].Xp[7].x=0.0209;	HC[1].Xp[7].y=17;
  HC[1].Xp[8].x=0.0511;	HC[1].Xp[8].y=21;
  HC[1].Xp[9].x=0.1618;	HC[1].Xp[9].y=23;
  HC[1].Xp[10].x=0.0979;HC[1].Xp[10].y=24;
  HC[1].Xp[11].x=0.0365;HC[1].Xp[11].y=25.1;
  HC[1].Xp[12].x=0.0132;HC[1].Xp[12].y=26;
  OH[1].Xp[0].x=0.0321;	OH[1].Xp[0].y=0.1;
  OH[1].Xp[1].x=0.0386;	OH[1].Xp[1].y=1.1;
  OH[1].Xp[2].x=0.057;	OH[1].Xp[2].y=2.1;
  OH[1].Xp[3].x=0.0966;	OH[1].Xp[3].y=3.1;
  OH[1].Xp[4].x=0.1842;	OH[1].Xp[4].y=5.1;
  OH[1].Xp[5].x=0.338;	OH[1].Xp[5].y=9.1;
  OH[1].Xp[6].x=0.4302;	OH[1].Xp[6].y=13.1;
  OH[1].Xp[7].x=0.494;	OH[1].Xp[7].y=17;
  OH[1].Xp[8].x=0.5009;	OH[1].Xp[8].y=21;
  OH[1].Xp[9].x=0.3516;	OH[1].Xp[9].y=23;
  OH[1].Xp[10].x=0.1631;OH[1].Xp[10].y=24;
  OH[1].Xp[11].x=0.0467;OH[1].Xp[11].y=25.1;
  OH[1].Xp[12].x=0.0137;OH[1].Xp[12].y=26;
  M[1].Xp[0].x=1.198;	M[1].Xp[0].y=0.1;
  M[1].Xp[1].x=1.221;	M[1].Xp[1].y=1.1;
  M[1].Xp[2].x=1.259;	M[1].Xp[2].y=2.1;
  M[1].Xp[3].x=1.274;	M[1].Xp[3].y=3.1;
  M[1].Xp[4].x=1.222;	M[1].Xp[4].y=5.1;
  M[1].Xp[5].x=1.083;	M[1].Xp[5].y=9.1;
  M[1].Xp[6].x=1.063;	M[1].Xp[6].y=13.1;
  M[1].Xp[7].x=1.081;	M[1].Xp[7].y=17;
  M[1].Xp[8].x=1.117;	M[1].Xp[8].y=21;
  M[1].Xp[9].x=1.061;	M[1].Xp[9].y=23;
  M[1].Xp[10].x=1.025;	M[1].Xp[10].y=24;
  M[1].Xp[11].x=1.018;	M[1].Xp[11].y=25.1;
  M[1].Xp[12].x=1.016;	M[1].Xp[12].y=26;
  T[1].pathlength();
  O2[1].pathlength();
  N2[1].pathlength();
  H2[1].pathlength();
  H2O[1].pathlength();
  CO[1].pathlength();
  CO2[1].pathlength();
  HC[1].pathlength();
  OH[1].pathlength();
  M[1].pathlength();

  //Station 3
  x[2 ] = 0.9;
  T[2].allocate(12);
  T[2].settype(SPLINE2D_LINEAR);
  O2[2].allocate(12);
  O2[2].settype(SPLINE2D_LINEAR);
  N2[2].allocate(12);
  N2[2].settype(SPLINE2D_LINEAR);
  H2[2].allocate(12);
  H2[2].settype(SPLINE2D_LINEAR);
  H2O[2].allocate(12);
  H2O[2].settype(SPLINE2D_LINEAR);
  CO[2].allocate(12);
  CO[2].settype(SPLINE2D_LINEAR);
  CO2[2].allocate(12);
  CO2[2].settype(SPLINE2D_LINEAR);
  HC[2].allocate(12);
  HC[2].settype(SPLINE2D_LINEAR);
  OH[2].allocate(12);
  OH[2].settype(SPLINE2D_LINEAR);
  M[2].allocate(12);
  M[2].settype(SPLINE2D_LINEAR);
  		
  T[2].Xp[0].x=646;	T[2].Xp[0].y=0.1;
  T[2].Xp[1].x=667;	T[2].Xp[1].y=1.1;
  T[2].Xp[2].x=730;	T[2].Xp[2].y=2.1;
  T[2].Xp[3].x=934;	T[2].Xp[3].y=4.1;
  T[2].Xp[4].x=1281;	T[2].Xp[4].y=8;
  T[2].Xp[5].x=1594;	T[2].Xp[5].y=12.1;
  T[2].Xp[6].x=1851;	T[2].Xp[6].y=16.1;
  T[2].Xp[7].x=1926;	T[2].Xp[7].y=20;
  T[2].Xp[8].x=1584;	T[2].Xp[8].y=22;
  T[2].Xp[9].x=780;	T[2].Xp[9].y=24;
  T[2].Xp[10].x=513;	T[2].Xp[10].y=25;
  T[2].Xp[11].x=382;	T[2].Xp[11].y=26;
  O2[2].Xp[0].x=24.2;	O2[2].Xp[0].y=0.1;
  O2[2].Xp[1].x=25.32;	O2[2].Xp[1].y=1.1;
  O2[2].Xp[2].x=28.58;	O2[2].Xp[2].y=2.1;
  O2[2].Xp[3].x=38.47;	O2[2].Xp[3].y=4.1;
  O2[2].Xp[4].x=52.68;	O2[2].Xp[4].y=8;
  O2[2].Xp[5].x=59.69;	O2[2].Xp[5].y=12.1;
  O2[2].Xp[6].x=61.96;	O2[2].Xp[6].y=16.1;
  O2[2].Xp[7].x=66.69;	O2[2].Xp[7].y=20;
  O2[2].Xp[8].x=72.26;	O2[2].Xp[8].y=22;
  O2[2].Xp[9].x=75.61;	O2[2].Xp[9].y=24;
  O2[2].Xp[10].x=76.08;	O2[2].Xp[10].y=25;
  O2[2].Xp[11].x=76.13;	O2[2].Xp[11].y=26;
  N2[2].Xp[0].x=6.4158;	N2[2].Xp[0].y=0.1;
  N2[2].Xp[1].x=6.2656;	N2[2].Xp[1].y=1.1;
  N2[2].Xp[2].x=5.7965;	N2[2].Xp[2].y=2.1;
  N2[2].Xp[3].x=4.376;	N2[2].Xp[3].y=4.1;
  N2[2].Xp[4].x=2.4834;	N2[2].Xp[4].y=8;
  N2[2].Xp[5].x=1.2971;	N2[2].Xp[5].y=12.1;
  N2[2].Xp[6].x=0.6185;	N2[2].Xp[6].y=16.1;
  N2[2].Xp[7].x=0.2504;	N2[2].Xp[7].y=20;
  N2[2].Xp[8].x=0.0805;	N2[2].Xp[8].y=22;
  N2[2].Xp[9].x=0.06;	N2[2].Xp[9].y=24;
  N2[2].Xp[10].x=0.0379;N2[2].Xp[10].y=25;
  N2[2].Xp[11].x=0.0273;N2[2].Xp[11].y=26;
  H2[2].Xp[0].x=5.544;	H2[2].Xp[0].y=0.1;
  H2[2].Xp[1].x=5.822;	H2[2].Xp[1].y=1.1;
  H2[2].Xp[2].x=6.561;	H2[2].Xp[2].y=2.1;
  H2[2].Xp[3].x=8.837;	H2[2].Xp[3].y=4.1;
  H2[2].Xp[4].x=12.026;	H2[2].Xp[4].y=8;
  H2[2].Xp[5].x=13.824;	H2[2].Xp[5].y=12.1;
  H2[2].Xp[6].x=14.004;	H2[2].Xp[6].y=16.1;
  H2[2].Xp[7].x=12.677;	H2[2].Xp[7].y=20;
  H2[2].Xp[8].x=9.355;	H2[2].Xp[8].y=22;
  H2[2].Xp[9].x=3.961;	H2[2].Xp[9].y=24;
  H2[2].Xp[10].x=2.311;	H2[2].Xp[10].y=25;
  H2[2].Xp[11].x=1.469;	H2[2].Xp[11].y=26;
  H2O[2].Xp[0].x=2.001;	H2O[2].Xp[0].y=0.1;
  H2O[2].Xp[1].x=2.099;	H2O[2].Xp[1].y=1.1;
  H2O[2].Xp[2].x=2.264;	H2O[2].Xp[2].y=2.1;
  H2O[2].Xp[3].x=3.073;	H2O[2].Xp[3].y=4.1;
  H2O[2].Xp[4].x=3.995;	H2O[2].Xp[4].y=8;
  H2O[2].Xp[5].x=4.339;	H2O[2].Xp[5].y=12.1;
  H2O[2].Xp[6].x=4.11;	H2O[2].Xp[6].y=16.1;
  H2O[2].Xp[7].x=2.784;	H2O[2].Xp[7].y=20;
  H2O[2].Xp[8].x=1.179;	H2O[2].Xp[8].y=22;
  H2O[2].Xp[9].x=0.572;	H2O[2].Xp[9].y=24;
  H2O[2].Xp[10].x=0.376;H2O[2].Xp[10].y=25;
  H2O[2].Xp[11].x=0.202;H2O[2].Xp[11].y=26;
  CO[2].Xp[0].x=2.785;	CO[2].Xp[0].y=0.1;
  CO[2].Xp[1].x=2.841;	CO[2].Xp[1].y=1.1;
  CO[2].Xp[2].x=3.117;	CO[2].Xp[2].y=2.1;
  CO[2].Xp[3].x=3.91;	CO[2].Xp[3].y=4.1;
  CO[2].Xp[4].x=4.71;	CO[2].Xp[4].y=8;
  CO[2].Xp[5].x=5.967;	CO[2].Xp[5].y=12.1;
  CO[2].Xp[6].x=7.019;	CO[2].Xp[6].y=16.1;
  CO[2].Xp[7].x=7.671;	CO[2].Xp[7].y=20;
  CO[2].Xp[8].x=5.966;	CO[2].Xp[8].y=22;
  CO[2].Xp[9].x=1.802;	CO[2].Xp[9].y=24;
  CO[2].Xp[10].x=0.628;	CO[2].Xp[10].y=25;
  CO[2].Xp[11].x=0.178;	CO[2].Xp[11].y=26;

  CO2[2].Xp[0].x=57.786;CO2[2].Xp[0].y=0.1;
  CO2[2].Xp[1].x=56.361;CO2[2].Xp[1].y=1.1;
  CO2[2].Xp[2].x=52.265;CO2[2].Xp[2].y=2.1;
  CO2[2].Xp[3].x=39.543;CO2[2].Xp[3].y=4.1;
  CO2[2].Xp[4].x=21.413;CO2[2].Xp[4].y=8;
  CO2[2].Xp[5].x=9.639;	CO2[2].Xp[5].y=12.1;
  CO2[2].Xp[6].x=2.94;	CO2[2].Xp[6].y=16.1;
  CO2[2].Xp[7].x=0.429;	CO2[2].Xp[7].y=20;
  CO2[2].Xp[8].x=0.121;	CO2[2].Xp[8].y=22;
  CO2[2].Xp[9].x=0.209;	CO2[2].Xp[9].y=24;
  CO2[2].Xp[10].x=0.164;CO2[2].Xp[10].y=25;
  CO2[2].Xp[11].x=0.102;CO2[2].Xp[11].y=26;

  HC[2].Xp[0].x=0;	HC[2].Xp[0].y=0.1;
  HC[2].Xp[1].x=0;	HC[2].Xp[1].y=1.1;
  HC[2].Xp[2].x=0;	HC[2].Xp[2].y=2.1;
  HC[2].Xp[3].x=0;	HC[2].Xp[3].y=4.1;
  HC[2].Xp[4].x=0.0024;	HC[2].Xp[4].y=8;
  HC[2].Xp[5].x=0.0208;	HC[2].Xp[5].y=12.1;
  HC[2].Xp[6].x=0.069;	HC[2].Xp[6].y=16.1;
  HC[2].Xp[7].x=0.1903;	HC[2].Xp[7].y=20;
  HC[2].Xp[8].x=0.1866;	HC[2].Xp[8].y=22;
  HC[2].Xp[9].x=0.0435;	HC[2].Xp[9].y=24;
  HC[2].Xp[10].x=0.0164;	HC[2].Xp[10].y=25;
  HC[2].Xp[11].x=0.0164;	HC[2].Xp[11].y=26;

  OH[2].Xp[0].x=0.0539;	OH[2].Xp[0].y=0.1;
  OH[2].Xp[1].x=0.057;	OH[2].Xp[1].y=1.1;
  OH[2].Xp[2].x=0.0677;	OH[2].Xp[2].y=2.1;
  OH[2].Xp[3].x=0.1138;	OH[2].Xp[3].y=4.1;
  OH[2].Xp[4].x=0.2193;	OH[2].Xp[4].y=8;
  OH[2].Xp[5].x=0.3688;	OH[2].Xp[5].y=12.1;
  OH[2].Xp[6].x=0.5077;	OH[2].Xp[6].y=16.1;
  OH[2].Xp[7].x=0.5285;	OH[2].Xp[7].y=20;
  OH[2].Xp[8].x=0.3316;	OH[2].Xp[8].y=22;
  OH[2].Xp[9].x=0.0731;	OH[2].Xp[9].y=24;
  OH[2].Xp[10].x=0.0265;	OH[2].Xp[10].y=25;
  OH[2].Xp[11].x=0.0073;	OH[2].Xp[11].y=26;

  M[2].Xp[0].x=1.256;	M[2].Xp[0].y=0.1;
  M[2].Xp[1].x=1.264;	M[2].Xp[1].y=1.1;
  M[2].Xp[2].x=1.277;	M[2].Xp[2].y=2.1;
  M[2].Xp[3].x=1.27;	M[2].Xp[3].y=4.1;
  M[2].Xp[4].x=1.174;	M[2].Xp[4].y=8;
  M[2].Xp[5].x=1.1;	M[2].Xp[5].y=12.1;
  M[2].Xp[6].x=1.113;	M[2].Xp[6].y=16.1;
  M[2].Xp[7].x=1.087;	M[2].Xp[7].y=20;
  M[2].Xp[8].x=1.038;	M[2].Xp[8].y=22;
  M[2].Xp[9].x=1.018;	M[2].Xp[9].y=24;
  M[2].Xp[10].x=1.019;	M[2].Xp[10].y=25;
  M[2].Xp[11].x=1.021;	M[2].Xp[11].y=26;
  T[2].pathlength();
  O2[2].pathlength();
  N2[2].pathlength();
  H2[2].pathlength();
  H2O[2].pathlength();
  CO[2].pathlength();
  CO2[2].pathlength();
  HC[2].pathlength();
  OH[2].pathlength();
  M[2].pathlength();

  //Station 4
  x[3 ] = 1.3;
  T[3].allocate(10);
  T[3].settype(SPLINE2D_LINEAR);
  O2[3].allocate(10);
  O2[3].settype(SPLINE2D_LINEAR);
  N2[3].allocate(10);
  N2[3].settype(SPLINE2D_LINEAR);
  H2[3].allocate(10);
  H2[3].settype(SPLINE2D_LINEAR);
  H2O[3].allocate(10);
  H2O[3].settype(SPLINE2D_LINEAR);
  CO[3].allocate(10);
  CO[3].settype(SPLINE2D_LINEAR);
  CO2[3].allocate(10);
  CO2[3].settype(SPLINE2D_LINEAR);
  HC[3].allocate(10);
  HC[3].settype(SPLINE2D_LINEAR);
  OH[3].allocate(10);
  OH[3].settype(SPLINE2D_LINEAR);
  M[3].allocate(10);
  M[3].settype(SPLINE2D_LINEAR);
  				
  T[3].Xp[0].x=833;	T[3].Xp[0].y=0.1;
  T[3].Xp[1].x=826;	T[3].Xp[1].y=1.1;
  T[3].Xp[2].x=854;	T[3].Xp[2].y=2.1;
  T[3].Xp[3].x=1132;	T[3].Xp[3].y=8.1;
  T[3].Xp[4].x=1420;	T[3].Xp[4].y=12.1;
  T[3].Xp[5].x=1744;	T[3].Xp[5].y=16.1;
  T[3].Xp[6].x=1407;	T[3].Xp[6].y=20;
  T[3].Xp[7].x=892;	T[3].Xp[7].y=22;
  T[3].Xp[8].x=498;	T[3].Xp[8].y=24;
  T[3].Xp[9].x=338;	T[3].Xp[9].y=26;
 
  O2[3].Xp[0].x=34.06;	O2[3].Xp[0].y=0.1;
  O2[3].Xp[1].x=33.75;	O2[3].Xp[1].y=1.1;
  O2[3].Xp[2].x=35.19;	O2[3].Xp[2].y=2.1;
  O2[3].Xp[3].x=47.28;	O2[3].Xp[3].y=8.1;
  O2[3].Xp[4].x=56.47;	O2[3].Xp[4].y=12.1;
  O2[3].Xp[5].x=62.96;	O2[3].Xp[5].y=16.1;
  O2[3].Xp[6].x=73.13;	O2[3].Xp[6].y=20;
  O2[3].Xp[7].x=75.6;	O2[3].Xp[7].y=22;
  O2[3].Xp[8].x=76.22;	O2[3].Xp[8].y=24;
  O2[3].Xp[9].x=76.14;	O2[3].Xp[9].y=26;
  N2[3].Xp[0].x=5.0038;	N2[3].Xp[0].y=0.1;
  N2[3].Xp[1].x=5.0297;	N2[3].Xp[1].y=1.1;
  N2[3].Xp[2].x=4.8426;	N2[3].Xp[2].y=2.1;
  N2[3].Xp[3].x=3.185;	N2[3].Xp[3].y=8.1;
  N2[3].Xp[4].x=1.8474;	N2[3].Xp[4].y=12.1;
  N2[3].Xp[5].x=0.6229;	N2[3].Xp[5].y=16.1;
  N2[3].Xp[6].x=0.0879;	N2[3].Xp[6].y=20;
  N2[3].Xp[7].x=0.035;	N2[3].Xp[7].y=22;
  N2[3].Xp[8].x=0.0284;	N2[3].Xp[8].y=24;
  N2[3].Xp[9].x=0.0104;	N2[3].Xp[9].y=26;
  H2[3].Xp[0].x=7.849;	H2[3].Xp[0].y=0.1;
  H2[3].Xp[1].x=7.796;	H2[3].Xp[1].y=1.1;
  H2[3].Xp[2].x=8.117;	H2[3].Xp[2].y=2.1;
  H2[3].Xp[3].x=10.875;	H2[3].Xp[3].y=8.1;
  H2[3].Xp[4].x=13.049;	H2[3].Xp[4].y=12.1;
  H2[3].Xp[5].x=13.063;	H2[3].Xp[5].y=16.1;
  H2[3].Xp[6].x=8.089;	H2[3].Xp[6].y=20;
  H2[3].Xp[7].x=4.596;	H2[3].Xp[7].y=22;
  H2[3].Xp[8].x=2.183;	H2[3].Xp[8].y=24;
  H2[3].Xp[9].x=1.152;	H2[3].Xp[9].y=26;
  H2O[3].Xp[0].x=2.408;	H2O[3].Xp[0].y=0.1;
  H2O[3].Xp[1].x=2.365;	H2O[3].Xp[1].y=1.1;
  H2O[3].Xp[2].x=2.441;	H2O[3].Xp[2].y=2.1;
  H2O[3].Xp[3].x=3.21;	H2O[3].Xp[3].y=8.1;
  H2O[3].Xp[4].x=3.733;	H2O[3].Xp[4].y=12.1;
  H2O[3].Xp[5].x=2.91;	H2O[3].Xp[5].y=16.1;
  H2O[3].Xp[6].x=0.722;	H2O[3].Xp[6].y=20;
  H2O[3].Xp[7].x=0.368;	H2O[3].Xp[7].y=22;
  H2O[3].Xp[8].x=0.218;	H2O[3].Xp[8].y=24;
  H2O[3].Xp[9].x=0.076;	H2O[3].Xp[9].y=26;
  CO[3].Xp[0].x=3.745;	CO[3].Xp[0].y=0.1;
  CO[3].Xp[1].x=3.639;	CO[3].Xp[1].y=1.1;
  CO[3].Xp[2].x=3.77;	CO[3].Xp[2].y=2.1;
  CO[3].Xp[3].x=4.407;	CO[3].Xp[3].y=8.1;
  CO[3].Xp[4].x=5.258;	CO[3].Xp[4].y=12.1;
  CO[3].Xp[5].x=6.701;	CO[3].Xp[5].y=16.1;
  CO[3].Xp[6].x=4.863;	CO[3].Xp[6].y=20;
  CO[3].Xp[7].x=2.326;	CO[3].Xp[7].y=22;
  CO[3].Xp[8].x=0.611;	CO[3].Xp[8].y=24;
  CO[3].Xp[9].x=0.083;	CO[3].Xp[9].y=26;
  CO2[3].Xp[0].x=44.877;CO2[3].Xp[0].y=0.1;
  CO2[3].Xp[1].x=45.444;CO2[3].Xp[1].y=1.1;
  CO2[3].Xp[2].x=43.566;CO2[3].Xp[2].y=2.1;
  CO2[3].Xp[3].x=28.271;CO2[3].Xp[3].y=8.1;
  CO2[3].Xp[4].x=15.479;CO2[3].Xp[4].y=12.1;
  CO2[3].Xp[5].x=4.019;	CO2[3].Xp[5].y=16.1;
  CO2[3].Xp[6].x=0.237;	CO2[3].Xp[6].y=20;
  CO2[3].Xp[7].x=0.105;	CO2[3].Xp[7].y=22;
  CO2[3].Xp[8].x=0.099;	CO2[3].Xp[8].y=24;
  CO2[3].Xp[9].x=0.043;	CO2[3].Xp[9].y=26;

  HC[3].Xp[0].x=0;	HC[3].Xp[0].y=0.1;
  HC[3].Xp[1].x=0;	HC[3].Xp[1].y=1.1;
  HC[3].Xp[2].x=0;	HC[3].Xp[2].y=2.1;
  HC[3].Xp[3].x=0.0008;	HC[3].Xp[3].y=8.1;
  HC[3].Xp[4].x=0.0241;	HC[3].Xp[4].y=12.1;
  HC[3].Xp[5].x=0.1274;	HC[3].Xp[5].y=16.1;
  HC[3].Xp[6].x=0.1205;	HC[3].Xp[6].y=20;
  HC[3].Xp[7].x=0.0389;	HC[3].Xp[7].y=22;
  HC[3].Xp[8].x=0.0096;	HC[3].Xp[8].y=24;
  HC[3].Xp[9].x=0.0123;	HC[3].Xp[9].y=26;
  OH[3].Xp[0].x=0.0891;	OH[3].Xp[0].y=0.1;
  OH[3].Xp[1].x=0.0886;	OH[3].Xp[1].y=1.1;
  OH[3].Xp[2].x=0.0944;	OH[3].Xp[2].y=2.1;
  OH[3].Xp[3].x=0.1605;	OH[3].Xp[3].y=8.1;
  OH[3].Xp[4].x=0.2592;	OH[3].Xp[4].y=12.1;
  OH[3].Xp[5].x=0.4059;	OH[3].Xp[5].y=16.1;
  OH[3].Xp[6].x=0.2576;	OH[3].Xp[6].y=20;
  OH[3].Xp[7].x=0.101;	OH[3].Xp[7].y=22;
  OH[3].Xp[8].x=0.0204;	OH[3].Xp[8].y=24;
  OH[3].Xp[9].x=0.0025;	OH[3].Xp[9].y=26;

  M[3].Xp[0].x=1.267;	M[3].Xp[0].y=0.1;
  M[3].Xp[1].x=1.27;	M[3].Xp[1].y=1.1;
  M[3].Xp[2].x=1.267;	M[3].Xp[2].y=2.1;
  M[3].Xp[3].x=1.211;	M[3].Xp[3].y=8.1;
  M[3].Xp[4].x=1.114;	M[3].Xp[4].y=12.1;
  M[3].Xp[5].x=1.092;	M[3].Xp[5].y=16.1;
  M[3].Xp[6].x=1.042;	M[3].Xp[6].y=20;
  M[3].Xp[7].x=1.021;	M[3].Xp[7].y=22;
  M[3].Xp[8].x=1.015;	M[3].Xp[8].y=24;
  M[3].Xp[9].x=1.019;	M[3].Xp[9].y=26;
  T[3].pathlength();
  O2[3].pathlength();
  N2[3].pathlength();
  H2[3].pathlength();
  H2O[3].pathlength();
  CO[3].pathlength();
  CO2[3].pathlength();
  HC[3].pathlength();
  OH[3].pathlength();
  M[3].pathlength();


  //Station 5
  x[4 ] = 1.8;
  T[4].allocate(8);
  T[4].settype(SPLINE2D_LINEAR);
  O2[4].allocate(8);
  O2[4].settype(SPLINE2D_LINEAR);
  N2[4].allocate(8);
  N2[4].settype(SPLINE2D_LINEAR);
  H2[4].allocate(8);
  H2[4].settype(SPLINE2D_LINEAR);
  H2O[4].allocate(8);
  H2O[4].settype(SPLINE2D_LINEAR);
  CO[4].allocate(8);
  CO[4].settype(SPLINE2D_LINEAR);
  CO2[4].allocate(8);
  CO2[4].settype(SPLINE2D_LINEAR);
  HC[4].allocate(8);
  HC[4].settype(SPLINE2D_LINEAR);
  OH[4].allocate(8);
  OH[4].settype(SPLINE2D_LINEAR);
  M[4].allocate(8);
  M[4].settype(SPLINE2D_LINEAR);
  			
  T[4].Xp[0].x=997;	T[4].Xp[0].y=0.1;
  T[4].Xp[1].x=994;	T[4].Xp[1].y=1.1;
  T[4].Xp[2].x=1023;	T[4].Xp[2].y=3.1;
  T[4].Xp[3].x=1069;	T[4].Xp[3].y=5.1;
  T[4].Xp[4].x=1203;	T[4].Xp[4].y=9.1;
  T[4].Xp[5].x=1506;	T[4].Xp[5].y=13.1;
  T[4].Xp[6].x=1337;	T[4].Xp[6].y=17.1;
  T[4].Xp[7].x=653;	T[4].Xp[7].y=21.1;

  O2[4].Xp[0].x=41.73;	O2[4].Xp[0].y=0.1;
  O2[4].Xp[1].x=41.27;	O2[4].Xp[1].y=1.1;
  O2[4].Xp[2].x=42.68;	O2[4].Xp[2].y=3.1;
  O2[4].Xp[3].x=44.67;	O2[4].Xp[3].y=5.1;
  O2[4].Xp[4].x=49.84;	O2[4].Xp[4].y=9.1;
  O2[4].Xp[5].x=59.25;	O2[4].Xp[5].y=13.1;
  O2[4].Xp[6].x=71.03;	O2[4].Xp[6].y=17.1;
  O2[4].Xp[7].x=76.18;	O2[4].Xp[7].y=21.1;
  N2[4].Xp[0].x=3.9343;	N2[4].Xp[0].y=0.1;
  N2[4].Xp[1].x=3.9615;	N2[4].Xp[1].y=1.1;
  N2[4].Xp[2].x=3.7888;	N2[4].Xp[2].y=3.1;
  N2[4].Xp[3].x=3.5239;	N2[4].Xp[3].y=5.1;
  N2[4].Xp[4].x=2.831;	N2[4].Xp[4].y=9.1;
  N2[4].Xp[5].x=1.3311;	N2[4].Xp[5].y=13.1;
  N2[4].Xp[6].x=0.24;	N2[4].Xp[6].y=17.1;
  N2[4].Xp[7].x=0.0254;	N2[4].Xp[7].y=21.1;
  H2[4].Xp[0].x=9.538;	H2[4].Xp[0].y=0.1;
  H2[4].Xp[1].x=9.538;	H2[4].Xp[1].y=1.1;
  H2[4].Xp[2].x=9.831;	H2[4].Xp[2].y=3.1;
  H2[4].Xp[3].x=10.266;	H2[4].Xp[3].y=5.1;
  H2[4].Xp[4].x=11.414;	H2[4].Xp[4].y=9.1;
  H2[4].Xp[5].x=12.59;	H2[4].Xp[5].y=13.1;
  H2[4].Xp[6].x=8.411;	H2[4].Xp[6].y=17.1;
  H2[4].Xp[7].x=3.056;	H2[4].Xp[7].y=21.1;
  H2O[4].Xp[0].x=2.782;	H2O[4].Xp[0].y=0.1;
  H2O[4].Xp[1].x=2.768;	H2O[4].Xp[1].y=1.1;
  H2O[4].Xp[2].x=2.857;	H2O[4].Xp[2].y=3.1;
  H2O[4].Xp[3].x=2.953;	H2O[4].Xp[3].y=5.1;
  H2O[4].Xp[4].x=3.167;	H2O[4].Xp[4].y=9.1;
  H2O[4].Xp[5].x=3.225;	H2O[4].Xp[5].y=13.1;
  H2O[4].Xp[6].x=1.293;	H2O[4].Xp[6].y=17.1;
  H2O[4].Xp[7].x=0.248;	H2O[4].Xp[7].y=21.1;
  CO[4].Xp[0].x=4.064;	CO[4].Xp[0].y=0.1;
  CO[4].Xp[1].x=4.062;	CO[4].Xp[1].y=1.1;
  CO[4].Xp[2].x=4.089;	CO[4].Xp[2].y=3.1;
  CO[4].Xp[3].x=4.138;	CO[4].Xp[3].y=5.1;
  CO[4].Xp[4].x=4.409;	CO[4].Xp[4].y=9.1;
  CO[4].Xp[5].x=5.422;	CO[4].Xp[5].y=13.1;
  CO[4].Xp[6].x=4.288;	CO[4].Xp[6].y=17.1;
  CO[4].Xp[7].x=1.196;	CO[4].Xp[7].y=21.1;
  CO2[4].Xp[0].x=35.396;	CO2[4].Xp[0].y=0.1;
  CO2[4].Xp[1].x=35.847;	CO2[4].Xp[1].y=1.1;
  CO2[4].Xp[2].x=34.086;	CO2[4].Xp[2].y=3.1;
  CO2[4].Xp[3].x=31.572;	CO2[4].Xp[3].y=5.1;
  CO2[4].Xp[4].x=24.985;	CO2[4].Xp[4].y=9.1;
  CO2[4].Xp[5].x=10.818;	CO2[4].Xp[5].y=13.1;
  CO2[4].Xp[6].x=1.338;	CO2[4].Xp[6].y=17.1;
  CO2[4].Xp[7].x=0.09;	CO2[4].Xp[7].y=21.1;

  HC[4].Xp[0].x=0;	HC[4].Xp[0].y=0.1;
  HC[4].Xp[1].x=0;	HC[4].Xp[1].y=1.1;
  HC[4].Xp[2].x=0;	HC[4].Xp[2].y=3.1;
  HC[4].Xp[3].x=0.0003;	HC[4].Xp[3].y=5.1;
  HC[4].Xp[4].x=0.0088;	HC[4].Xp[4].y=9.1;
  HC[4].Xp[5].x=0.0661;	HC[4].Xp[5].y=13.1;
  HC[4].Xp[6].x=0.0993;	HC[4].Xp[6].y=17.1;
  HC[4].Xp[7].x=0.0187;	HC[4].Xp[7].y=21.1;

  OH[4].Xp[0].x=0.1245;	OH[4].Xp[0].y=0.1;
  OH[4].Xp[1].x=0.1246;	OH[4].Xp[1].y=1.1;
  OH[4].Xp[2].x=0.1307;	OH[4].Xp[2].y=3.1;
  OH[4].Xp[3].x=0.141;	OH[4].Xp[3].y=5.1;
  OH[4].Xp[4].x=0.1751;	OH[4].Xp[4].y=9.1;
  OH[4].Xp[5].x=0.2805;	OH[4].Xp[5].y=13.1;
  OH[4].Xp[6].x=0.2266;	OH[4].Xp[6].y=17.1;
  OH[4].Xp[7].x=0.0494;	OH[4].Xp[7].y=21.1;

  M[4].Xp[0].x=1.244;	M[4].Xp[0].y=0.1;
  M[4].Xp[1].x=1.246;	M[4].Xp[1].y=1.1;
  M[4].Xp[2].x=1.236;	M[4].Xp[2].y=3.1;
  M[4].Xp[3].x=1.227;	M[4].Xp[3].y=5.1;
  M[4].Xp[4].x=1.182;	M[4].Xp[4].y=9.1;
  M[4].Xp[5].x=1.1;	M[4].Xp[5].y=13.1;
  M[4].Xp[6].x=1.044;	M[4].Xp[6].y=17.1;
  M[4].Xp[7].x=1.016;	M[4].Xp[7].y=21.1;

  T[4].pathlength();
  O2[4].pathlength();
  N2[4].pathlength();
  H2[4].pathlength();
  H2O[4].pathlength();
  CO[4].pathlength();
  CO2[4].pathlength();
  HC[4].pathlength();
  OH[4].pathlength();
  M[4].pathlength();
 
  //Station 6
  x[5 ] = 2.4;
  T[5].allocate(9);
  T[5].settype(SPLINE2D_LINEAR);
  O2[5].allocate(9);
  O2[5].settype(SPLINE2D_LINEAR);
  N2[5].allocate(9);
  N2[5].settype(SPLINE2D_LINEAR);
  H2[5].allocate(9);
  H2[5].settype(SPLINE2D_LINEAR);
  H2O[5].allocate(9);
  H2O[5].settype(SPLINE2D_LINEAR);
  CO[5].allocate(9);
  CO[5].settype(SPLINE2D_LINEAR);
  CO2[5].allocate(9);
  CO2[5].settype(SPLINE2D_LINEAR);
  HC[5].allocate(9);
  HC[5].settype(SPLINE2D_LINEAR);
  OH[5].allocate(9);
  OH[5].settype(SPLINE2D_LINEAR);
  M[5].allocate(9);
  M[5].settype(SPLINE2D_LINEAR);
  			
  T[5].Xp[0].x=1216;	T[5].Xp[0].y=0.1;
  T[5].Xp[1].x=1224;	T[5].Xp[1].y=1.1;
  T[5].Xp[2].x=1244;	T[5].Xp[2].y=3.1;
  T[5].Xp[3].x=1282;	T[5].Xp[3].y=5.1;
  T[5].Xp[4].x=1370;	T[5].Xp[4].y=7.1;
  T[5].Xp[5].x=1562;	T[5].Xp[5].y=11.1;
  T[5].Xp[6].x=1294;	T[5].Xp[6].y=15.1;
  T[5].Xp[7].x=798;	T[5].Xp[7].y=19.1;
  T[5].Xp[8].x=557;	T[5].Xp[8].y=21.1;
  O2[5].Xp[0].x=49.98;	O2[5].Xp[0].y=0.1;	
  O2[5].Xp[1].x=50.04;	O2[5].Xp[1].y=1.1;
  O2[5].Xp[2].x=50.91;	O2[5].Xp[2].y=3.1;
  O2[5].Xp[3].x=52;	O2[5].Xp[3].y=5.1;
  O2[5].Xp[4].x=54.41;	O2[5].Xp[4].y=7.1;
  O2[5].Xp[5].x=59.85;	O2[5].Xp[5].y=11.1;
  O2[5].Xp[6].x=70.17;	O2[5].Xp[6].y=15.1;
  O2[5].Xp[7].x=75.5;	O2[5].Xp[7].y=19.1;
  O2[5].Xp[8].x=76.39;	O2[5].Xp[8].y=21.1;

  N2[5].Xp[0].x=2.829;	N2[5].Xp[0].y=0.1;
  N2[5].Xp[1].x=2.8073;	N2[5].Xp[1].y=1.1;
  N2[5].Xp[2].x=2.7051;	N2[5].Xp[2].y=3.1;
  N2[5].Xp[3].x=2.5537;	N2[5].Xp[3].y=5.1;
  N2[5].Xp[4].x=2.1532;	N2[5].Xp[4].y=7.1;
  N2[5].Xp[5].x=1.2713;	N2[5].Xp[5].y=11.1;
  N2[5].Xp[6].x=0.3089;	N2[5].Xp[6].y=15.1;
  N2[5].Xp[7].x=0.0538;	N2[5].Xp[7].y=19.1;
  N2[5].Xp[8].x=0.0177;	N2[5].Xp[8].y=21.1;

  H2[5].Xp[0].x=11.405;	H2[5].Xp[0].y=0.1;
  H2[5].Xp[1].x=11.439;	H2[5].Xp[1].y=1.1;
  H2[5].Xp[2].x=11.579;	H2[5].Xp[2].y=3.1;
  H2[5].Xp[3].x=11.949;	H2[5].Xp[3].y=5.1;
  H2[5].Xp[4].x=12.359;	H2[5].Xp[4].y=7.1;
  H2[5].Xp[5].x=12.875;	H2[5].Xp[5].y=11.1;
  H2[5].Xp[6].x=8.149;	H2[5].Xp[6].y=15.1;
  H2[5].Xp[7].x=3.931;	H2[5].Xp[7].y=19.1;
  H2[5].Xp[8].x=2.336;	H2[5].Xp[8].y=21.1;
  H2O[5].Xp[0].x=3.255;	H2O[5].Xp[0].y=0.1;
  H2O[5].Xp[1].x=3.329;	H2O[5].Xp[1].y=1.1;
  H2O[5].Xp[2].x=3.299;	H2O[5].Xp[2].y=3.1;
  H2O[5].Xp[3].x=3.434;	H2O[5].Xp[3].y=5.1;
  H2O[5].Xp[4].x=3.497;	H2O[5].Xp[4].y=7.1;
  H2O[5].Xp[5].x=3.306;	H2O[5].Xp[5].y=11.1;
  H2O[5].Xp[6].x=1.418;	H2O[5].Xp[6].y=15.1;
  H2O[5].Xp[7].x=0.381;	H2O[5].Xp[7].y=19.1;
  H2O[5].Xp[8].x=0.171;	H2O[5].Xp[8].y=21.1;

  CO[5].Xp[0].x=4.215;	CO[5].Xp[0].y=0.1;
  CO[5].Xp[1].x=4.165;	CO[5].Xp[1].y=1.1;
  CO[5].Xp[2].x=4.24;	CO[5].Xp[2].y=3.1;
  CO[5].Xp[3].x=4.312;	CO[5].Xp[3].y=5.1;
  CO[5].Xp[4].x=4.654;	CO[5].Xp[4].y=7.1;
  CO[5].Xp[5].x=5.272;	CO[5].Xp[5].y=11.1;
  CO[5].Xp[6].x=3.921;	CO[5].Xp[6].y=15.1;
  CO[5].Xp[7].x=1.746;	CO[5].Xp[7].y=19.1;
  CO[5].Xp[8].x=0.781;	CO[5].Xp[8].y=21.1;

  CO2[5].Xp[0].x=24.77;	CO2[5].Xp[0].y=0.1;
  CO2[5].Xp[1].x=24.793;CO2[5].Xp[1].y=1.1;
  CO2[5].Xp[2].x=23.783;CO2[5].Xp[2].y=3.1;
  CO2[5].Xp[3].x=22.194;CO2[5].Xp[3].y=5.1;
  CO2[5].Xp[4].x=18.224;CO2[5].Xp[4].y=7.1;
  CO2[5].Xp[5].x=10.104;CO2[5].Xp[5].y=11.1;
  CO2[5].Xp[6].x=1.964;	CO2[5].Xp[6].y=15.1;
  CO2[5].Xp[7].x=0.258;	CO2[5].Xp[7].y=19.1;
  CO2[5].Xp[8].x=0.07;	CO2[5].Xp[8].y=21.1;

  HC[5].Xp[0].x=0.0009;	HC[5].Xp[0].y=0.1;
  HC[5].Xp[1].x=0.0009;	HC[5].Xp[1].y=1.1;
  HC[5].Xp[2].x=0.0024;	HC[5].Xp[2].y=3.1;
  HC[5].Xp[3].x=0.0059;	HC[5].Xp[3].y=5.1;
  HC[5].Xp[4].x=0.0181;	HC[5].Xp[4].y=7.1;
  HC[5].Xp[5].x=0.0629;	HC[5].Xp[5].y=11.1;
  HC[5].Xp[6].x=0.0907;	HC[5].Xp[6].y=15.1;
  HC[5].Xp[7].x=0.0389;	HC[5].Xp[7].y=19.1;
  HC[5].Xp[8].x=0.0182;	HC[5].Xp[8].y=21.1;

  OH[5].Xp[0].x=0.1683;	OH[5].Xp[0].y=0.1;
  OH[5].Xp[1].x=0.1719;	OH[5].Xp[1].y=1.1;
  OH[5].Xp[2].x=0.1786;	OH[5].Xp[2].y=3.1;
  OH[5].Xp[3].x=0.1855;	OH[5].Xp[3].y=5.1;
  OH[5].Xp[4].x=0.2228;	OH[5].Xp[4].y=7.1;
  OH[5].Xp[5].x=0.2949;	OH[5].Xp[5].y=11.1;
  OH[5].Xp[6].x=0.2082;	OH[5].Xp[6].y=15.1;
  OH[5].Xp[7].x=0.0787;	OH[5].Xp[7].y=19.1;
  OH[5].Xp[8].x=0.0305;	OH[5].Xp[8].y=21.1;

  M[5].Xp[0].x=1.181;	M[5].Xp[0].y=0.1;
  M[5].Xp[1].x=1.184;	M[5].Xp[1].y=1.1;
  M[5].Xp[2].x=1.173;	M[5].Xp[2].y=3.1;
  M[5].Xp[3].x=1.16;	M[5].Xp[3].y=5.1;
  M[5].Xp[4].x=1.136;	M[5].Xp[4].y=7.1;
  M[5].Xp[5].x=1.098;	M[5].Xp[5].y=11.1;
  M[5].Xp[6].x=1.056;	M[5].Xp[6].y=15.1;
  M[5].Xp[7].x=1.02;	M[5].Xp[7].y=19.1;
  M[5].Xp[8].x=1.014;	M[5].Xp[8].y=21.1; 

  T[5].pathlength();
  O2[5].pathlength();
  N2[5].pathlength();
  H2[5].pathlength();
  H2O[5].pathlength();
  CO[5].pathlength();
  CO2[5].pathlength();
  HC[5].pathlength();
  OH[5].pathlength();
  M[5].pathlength();

  
  
}

//! NonreactiveVelocityField::deallocate -- Deallocate memory for
//!                                            splines and arrays.
inline void ReactiveScalarField_Fuel_CH4H2::deallocate(void) {
  // Deallocate splines.
  for (int n = 0; n <Ns; n++) {
    M[n].deallocate();
    T[n].deallocate();
    O2[n].deallocate();
    N2[n].deallocate();
    H2[n].deallocate();
    H2O[n].deallocate();
    CO[n].deallocate();
    CO2[n].deallocate();
    HC[n].deallocate();
    OH[n].deallocate();
  
 
  }

  // Deallocate arrays.
  delete []M; M = NULL;
  delete []T; T = NULL;
  delete []N2; N2 = NULL;
  delete []H2; H2 = NULL;
  delete []O2; O2 = NULL;
  delete []H2O; H2O = NULL;
  delete []CO; CO = NULL;
  delete []CO2; CO2 = NULL;
  delete []HC; HC = NULL;
  delete []OH; OH = NULL;
  delete []x; x = NULL;
}

class ReactiveScalarField_Fuel_CH4 {
/* Coflow air + Fuel (CH4)  */
/*                          */
 private:
  int       Ns; //!< Number of stations.
  Spline2D  *M; //!< Array of mixture fraction splines at each station.
  Spline2D  *T; //!< Array of mixture temperature splines at each station.
  Spline2D  *O2; //!< Array of species fraction splines at each station.
  Spline2D  *N2; //!< Array of species fraction splines at each station.
  Spline2D  *H2; //!< Array of species fraction splines at each station.
  Spline2D  *H2O; //!< Array of species fraction splines at each station.
  Spline2D  *CO; //!< Array of species fraction splines at each station.
  Spline2D  *CO2; //!< Array of species fraction splines at each station.
  Spline2D  *HC; //!< Array of species fraction splines at each station.
  Spline2D  *OH; //!< Array of species fraction splines at each station.
 
  double    *x; //!< Array of station x-coordinates (axial direction).
  
 public:
  //! Creation constructor.
  ReactiveScalarField_Fuel_CH4(void) { allocate(); }
  
  //! Destructor.
  ~ReactiveScalarField_Fuel_CH4(void) { deallocate(); }
  
  //! Allocate memory and initialize arrays.
  void allocate(void);
  
  //! Deallocate memory for splines and arrays.
  void deallocate(void);

  //! Interpolation function.
 
  Vector2D interpolation(const Vector2D &Xt, double []);

};

//! NonreactiveBluffBodeuBurner::allocate -- Allocate memory and
//!                                          initialize arrays.
 inline void ReactiveScalarField_Fuel_CH4::allocate(void) {
/*   // Set the number of stations. */
//------------------------------------------------------------------//
/* Need at least 3 stations, otherwise the interpolation of spline  */
/* can not be doen, the running error message will be 
   segmentation error. For this fuel case, no more experimental data
   available, so currently not testing this case                       */
//-----------------------------------------------------------------//
//  Ns = 1;
  // Allocate the arrays.
 /*  M = new Spline2D[Ns]; */
/*   T = new Spline2D[Ns]; */
/*   O2 = new Spline2D[Ns]; */
/*   N2 = new Spline2D[Ns]; */
/*   H2 = new Spline2D[Ns]; */
/*   H2O= new Spline2D[Ns]; */
/*   CO= new Spline2D[Ns]; */
/*   CO2 = new Spline2D[Ns]; */
/*   HC = new Spline2D[Ns]; */
/*   OH = new Spline2D[Ns]; */

/*  /\*  C = new Spline2D *[9]; //nine species *\/ */
/* /\*   for(int i; i<9; i++){ *\/ */
/* /\*     C[i] = new Spline2D[Ns]; *\/ */
/* /\*   } *\/ */

/*   x = new double[Ns]; */
/*   //Station 1 */
/*   x[0 ] = 1.9; */
/*   T[0].allocate(2); */
/*   T[0].settype(SPLINE2D_LINEAR); */
/*   O2[0].allocate(2); */
/*   O2[0].settype(SPLINE2D_LINEAR); */
/*   N2[0].allocate(2); */
/*   N2[0].settype(SPLINE2D_LINEAR); */
/*   H2[0].allocate(2); */
/*   H2[0].settype(SPLINE2D_LINEAR); */
/*   H2O[0].allocate(2); */
/*   H2O[0].settype(SPLINE2D_LINEAR); */
/*   CO[0].allocate(2); */
/*   CO[0].settype(SPLINE2D_LINEAR); */
/*   CO2[0].allocate(2); */
/*   CO2[0].settype(SPLINE2D_LINEAR); */
/*   HC[0].allocate(2); */
/*   HC[0].settype(SPLINE2D_LINEAR); */
/*   OH[0].allocate(2); */
/*   OH[0].settype(SPLINE2D_LINEAR); */
/*   M[0].allocate(2); */
/*   M[0].settype(SPLINE2D_LINEAR); */
  
/*   T[0].Xp[0].x=	1120;	T[0].Xp[0].y=	10.0; */
/*   T[0].Xp[1].x=	840;	T[0].Xp[1].y=	13.0; */
/*   //O2 fraction */
/*   O2[0].Xp[0].x=11.416;	O2[0].Xp[0].y= 10.0; */
/*   O2[0].Xp[1].x=14.825;	O2[0].Xp[1].y= 13.0; */

/*   N2[0].Xp[0].x=71.67 ;N2[0].Xp[0].y=10.0; */
/*   N2[0].Xp[1].x=74.40;	N2[0].Xp[1].y=13.0; */

/*   H2[0].Xp[0].x= 0.1377;	H2[0].Xp[0].y=10.0; */
/*   H2[0].Xp[1].x=0.0743;	H2[0].Xp[1].y=13.0; */

/*   H2O[0].Xp[0].x=5.052;	H2O[0].Xp[0].y=	10.0; */
/*   H2O[0].Xp[1].x=3.238;	H2O[0].Xp[1].y=	13.0; */

/*   CO[0].Xp[0].x=1.281;	CO[5].Xp[0].y=10.0; */
/*   CO[0].Xp[1].x=0.765;	CO[5].Xp[1].y=13.0; */
  
/*   CO2[0].Xp[0].x=7.169;  CO2[0].Xp[0].y=10.0; */
/*   CO2[0].Xp[1].x=5.439; CO2[0].Xp[1].y=13.0; */

/*   HC[0].Xp[0].x=3.214;	HC[0].Xp[0].y=10.0; */
/*   HC[0].Xp[1].x=1.227;	HC[0].Xp[1].y=13.0; */

/*   OH[0].Xp[0].x=0.0588;	OH[0].Xp[0].y=10.0; */
/*   OH[0].Xp[1].x=0.0349;	OH[0].Xp[1].y=13.0; */
  
/*   M[0].Xp[0].x=1.0;	M[0].Xp[0].y=10.0; */
/*   M[0].Xp[1].x=1.0;	M[0].Xp[1].y=13.0; */
 
/*   T[0].pathlength(); */
/*   O2[0].pathlength(); */
/*   N2[0].pathlength(); */
/*   H2[0].pathlength(); */
/*   H2O[0].pathlength(); */
/*   CO[0].pathlength(); */
/*   CO2[0].pathlength(); */
/*   HC[0].pathlength(); */
/*   OH[0].pathlength(); */
/*   M[0].pathlength(); */
 

/*   //Station 1 */
/*   x[0 ] = 1.9; */
/*   T[0].allocate(2); */
/*   T[0].settype(SPLINE2D_LINEAR); */
/*   O2[0].allocate(2); */
/*   O2[0].settype(SPLINE2D_LINEAR); */
/*   N2[0].allocate(2); */
/*   N2[0].settype(SPLINE2D_LINEAR); */
/*   H2[0].allocate(2); */
/*   H2[0].settype(SPLINE2D_LINEAR); */
/*   H2O[0].allocate(2); */
/*   H2O[0].settype(SPLINE2D_LINEAR); */
/*   CO[0].allocate(2); */
/*   CO[0].settype(SPLINE2D_LINEAR); */
/*   CO2[0].allocate(2); */
/*   CO2[0].settype(SPLINE2D_LINEAR); */
/*   HC[0].allocate(2); */
/*   HC[0].settype(SPLINE2D_LINEAR); */
/*   OH[0].allocate(2); */
/*   OH[0].settype(SPLINE2D_LINEAR); */
/*   M[0].allocate(2); */
/*   M[0].settype(SPLINE2D_LINEAR); */
  // Set the number of stations.


}
//! ReactiveVelocityField::deallocate -- Deallocate memory for
//!                                            splines and arrays.
inline void ReactiveScalarField_Fuel_CH4::deallocate(void) {
  // Deallocate splines.
  for (int n = 0; n <Ns; n++) {
    M[n].deallocate();
    T[n].deallocate();
    O2[n].deallocate();
    N2[n].deallocate();
    H2[n].deallocate();
    H2O[n].deallocate();
    CO[n].deallocate();
    CO2[n].deallocate();
    HC[n].deallocate();
    OH[n].deallocate();
  
 
  }

  // Deallocate arrays.
  delete []M; M = NULL;
  delete []T; T = NULL;
  delete []N2; N2 = NULL;
  delete []H2; H2 = NULL;
  delete []O2; O2 = NULL;
  delete []H2O; H2O = NULL;
  delete []CO; CO = NULL;
  delete []CO2; CO2 = NULL;
  delete []HC; HC = NULL;
  delete []OH; OH = NULL;
  delete []x; x = NULL;
}

//! Routine: get U -- Conduct the interpolation of the data spline.
extern double getU(const Vector2D &Xt, const Spline2D &U);

#endif /* _BLUFFBODYBURNER_INCLUDED  */
