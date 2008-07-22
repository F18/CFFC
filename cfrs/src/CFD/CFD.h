/* CFD.h:  Header file defining CFD subroutines and macros. */

#ifndef _CFD_INCLUDED
#define _CFD_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>
#include <ctime>

using namespace std;

/* Include math macro file. */

#include "../../../src_2D/Math/Math.h"

/* Include vector2d headers */
#include "../../../src_2D/Math/Vector2D.h"

/* Include require header file */
#include "../../../src_2D/Utilities/Utilities.h"

/********************************************************
 * CFD -- CFDkit+caboodle library name and version.     *
 ********************************************************/

// Name.
inline string CFDkit_Name() {
  return ("CFDkit+caboodle");
}

// Version.
inline string CFDkit_Version() {
  return ("CFDkit+caboodle, Version 1.5, UTIAS CFD & Propulsion Group, 1999-2004.");
}

/**********************************************************************
 * CFD -- CFFC library name and version.                              *
 **********************************************************************/

// Name.
inline string CFFC_Name() {
  return ("CFFC");
}

// Version.
inline string CFFC_Version() {
  return ("CFFC, Version 0.01, UTIAS CFD & Propulsion Group, 1999-2007.");
}

/********************************************************
 * CFD -- Date and time.                                *
 ********************************************************/

inline char *Date_And_Time() {
  char *string_ptr; time_t epoch_time;
  epoch_time = time(NULL);
  string_ptr = ctime(&epoch_time);
  string_ptr[strlen(string_ptr)-1] = '\0';
  return (string_ptr);
}

/********************************************************
 * CFD -- General Purpose Command Codes.                *
 ********************************************************/

#define	OFF                           0
#define	ON                            1

/********************************************************
 * CFD -- Input parameter command codes.                *
 ********************************************************/

#define EXECUTE_CODE                        10000
#define TERMINATE_CODE                      10001
#define CONTINUE_CODE                       10002
#define WRITE_OUTPUT_CODE                   10003
#define WRITE_OUTPUT_CELLS_CODE             10004
#define WRITE_OUTPUT_QUASI3D_CODE           10005
#define WRITE_RESTART_CODE                  10006
#define WRITE_OUTPUT_GRID_CODE              10007
#define WRITE_GRID_DEFINITION_CODE          10008
#define WRITE_OUTPUT_GRID_NODES_CODE        10009
#define WRITE_OUTPUT_GRID_CELLS_CODE        10010
#define COMMENT_CODE                        10011
#define REFINE_GRID_CODE                    10012

#define WRITE_OUTPUT_POINT_FORMAT_CODE      10020
#define WRITE_OUTPUT_INTERFACE_NODES_CODE   10021
#define WRITE_OUTPUT_LEVEL_SET_CODE         10022
#define WRITE_OUTPUT_LEVEL_SET_CELLS_CODE   10023

#define WRITE_OUTPUT_RHS_CODE               10026
#define WRITE_OUTPUT_PERTURB_CODE           10028

#define WRITE_OUTPUT_RINGLEB_FLOW_CODE      10030

/**************************
Lucian's stuff
**************************/
#define WRITE_INTERFACE_NODES            10040
#define WRITE_OUTPUT_DERIVATIVES_CODE    10041
#define WRITE_OUTPUT_COMPARISON_CODE     10042
#define WRITE_OUTPUT_FUNCTION_GRAPH_CODE 10043
#define WRITE_L1_NORM                    10044
#define WRITE_L2_NORM                    10045
#define WRITE_NORM_ON_SCREEN             10046
#define WRITE_OUTPUT_ACCURACY_CODE       10047
#define ASSESS_ACCURACY_CODE             10048
#define WRITE_OUTPUT_FULL_GRID_NODES_CODE 10049
#define WRITE_UNIFIED_OUTPUT_CODE         10050
#define WRITE_RECONSTRUCTED_FUNCTIONS    10051
/**************************
Lucian's stuff
**************************/

#define	INVALID_INPUT_CODE                 -10000
#define	INVALID_INPUT_VALUE                -10001

/********************************************************
 * CFD -- I/O Types.                                    *
 ********************************************************/

#define	IO_GNUPLOT       1
#define	IO_TECPLOT       2
#define	IO_MATLAB        3
#define	IO_OCTAVE        4

/********************************************************
 * CFD -- Types of algebraic grid point distribution    *
 *        stretching functions.                         *
 ********************************************************/

#define STRETCHING_FCN_LINEAR            0
#define STRETCHING_FCN_MIN_CLUSTERING    1
#define STRETCHING_FCN_MAX_CLUSTERING    2
#define STRETCHING_FCN_MINMAX_CLUSTERING 3
#define STRETCHING_FCN_MIDPT_CLUSTERING  4
#define STRETCHING_FCN_SINE              5
#define STRETCHING_FCN_COSINE            6

/********************************************************
 * CFD -- Grid Types.                                   *
 ********************************************************/

#define GRID_CARTESIAN_UNIFORM                0
#define GRID_SQUARE                           1
#define GRID_RECTANGULAR_BOX                  2
#define GRID_FLAT_PLATE                       3
#define GRID_PIPE                             4
#define GRID_BLUNT_BODY                       5
#define GRID_ROCKET_MOTOR                     6
#define GRID_ROCKET_MOTOR_COLD_FLOW           7
#define GRID_CIRCULAR_CYLINDER                8
#define GRID_ELLIPSE                          9
#define GRID_NACA_AEROFOIL                   10
#define GRID_NASA_ROTOR_37                   11
#define GRID_NASA_ROTOR_67                   12
#define GRID_FREE_JET                        13
#define GRID_WEDGE                           14
#define GRID_UNSTEADY_BLUNT_BODY             15
#define GRID_RINGLEB_FLOW                    16
#define GRID_BUMP_CHANNEL_FLOW               17



/**************************
Lucian's stuff
**************************/
#define GRID_UNIFORM                       20
#define GRID_NONUNIFORM                    21
#define GRID_PREDEFINED1                   22
#define GRID_PREDEFINED2                   23

/**************************
Lucian's stuff
**************************/


#define GRID_ICEMCFD                       1000
#define GRID_READ_FROM_DEFINITION_FILE    10000
#define GRID_READ_FROM_GRID_DATA_FILE     10001

#define GRID_COUETTE                         50
#define GRID_1DFLAME                         51
#define GRID_LAMINAR_FLAME                   52

/**************************
3D stuff
**************************/

#define GRID_CUBE                           223
#define GRID_CHANNEL                         24
#define GRID_CHANNEL_ZDIR          GRID_CHANNEL
#define GRID_CHANNEL_XDIR                    25
#define GRID_CHANNEL_YDIR                    26
// #define GRID_COUETTE               see above
#define GRID_COUETTE_ZDIR          GRID_COUETTE
#define GRID_COUETTE_XDIR                    28
#define GRID_COUETTE_YDIR                    29
#define GRID_BLUFF_BODY_BURNER               30
#define GRID_PERIODIC_BOX                    31
#define GRID_PERIODIC_BOX_WITH_INFLOW        32
#define GRID_BUNSEN_BURNER                   33
#define GRID_BUNSEN_BOX                      34
#define GRID_BUNSEN_INFLOW                   35
#define GRID_TURBULENCE_BOX                  36


/********************************************************
 * Method for solving the least squares problem         *
 ********************************************************/

#define NORMAL_EQ                       0
#define LSH                             1

/********************************************************
 * Reconstruction -- Function Definition Types          *
 ********************************************************/

#define FUNCTION_DEFAULT                   1
#define FUNCTION_EXAMPLE_1                 2
#define FUNCTION_EXAMPLE_2                 3 
#define FUNCTION_EXAMPLE_3                 4 
#define FUNCTION_EXAMPLE_4                 5 
#define FUNCTION_EXAMPLE_5                 6 
#define FUNCTION_EXAMPLE_6                 7 
#define FUNCTION_EXAMPLE_7                 8 
#define FUNCTION_EXAMPLE_8                 9 
#define FUNCTION_EXAMPLE_9                 10 
#define FUNCTION_EXAMPLE_10                 11

#define ENO                                11
#define ENO_LS                             12
#define WENO                               13
#define DD_ENO                             14
#define SpectralDiff                       15
#define CENO                               16

/********************************************************
 * CFD -- Boundary Condition Types.                     *
 ********************************************************/

//----DEFAULT NAMES----//                          Core | Scope

//--Undefined (or interior)
#define	BC_NULL                             0   // Core | CGNS
#define BC_NONE                       BC_NULL   //      | Duplicate

//--Speciality interior
#define BC_INTERIOR                      1000   // Core | CFDKit

//--Collapsed blocks
#define BC_DEGENERATE                    2000   // Core | Header
#define BC_DEGENERATE_LINE               2001   //      | CGNS
#define BC_DEGENERATE_POINT              2002   //      | CGNS
#define BC_AXISYMMETRIC_WEDGE            2003   //      | CGNS

//--General?
#define BC_GENERAL                       3000   // Core | CGNS

//--Extrapolation
#define BC_EXTRAPOLATE                   4000   // Core | CGNS
#define	BC_CONSTANT_EXTRAPOLATION        4001   //      | CFDKit
#define	BC_LINEAR_EXTRAPOLATION          4002   //      | CFDKit

//--1st Derivative
#define BC_NEUMANN                       5000   // Core | CGNS

//--Periodic
#define BC_PERIODIC                      6000   // Core | CFDKit

//--Symmetry plane
#define BC_SYMMETRY                      7000   // Core | Header
#define BC_SYMMETRY_PLANE                7001   //      | CGNS
#define BC_REFLECTION       BC_SYMMETRY_PLANE   //      | Duplicate
#define BC_SYMMETRY_POLAR                7002   //      | CGNS

//--Far field
#define BC_FARFIELD                      8000   // Core | CGNS
#define BC_CHARACTERISTIC                8001   //      | CFDKit ?
#define BC_FIXED_PRESSURE                8002   //      | CFDKit

//--Inflow
#define BC_INFLOW                        9000   // Core | CGNS
#define BC_INFLOW_SUBSONIC               9001   //      | CGNS
#define BC_INFLOW_SUPERSONIC             9002   //      | CGNS
#define BC_TUNNEL_INFLOW                 9003   //      | CGNS
#define BC_FLAME_INFLOW                  9003   //      | CFDKit 

//--Outflow
#define BC_OUTFLOW                      10000   // Core | CGNS
#define BC_OPEN_END                BC_OUTFLOW   //      | Duplicate ?
#define BC_OUTFLOW_SUBSONIC             10001   //      | CGNS
#define BC_OUTFLOW_SUPERSONIC           10002   //      | CGNS
#define BC_TUNNEL_OUTFLOW               10003   //      | CGNS
#define BC_FLAME_OUTFLOW                10004   //      | CFDKit

//--Frozen
#define BC_ROBIN                        10999   //      | CFDKit
#define BC_DIRICHLET                    11000   // Core | CGNS
#define BC_FIXED                 BC_DIRICHLET   //      | Duplicate

//--Wall
//  Speciality wall
#define BC_ABSORPTION                   11800   //      | CFDKit
//  Moving wall
#define BC_MOVING_WALL                  11900   //      | CFDKit
#define BC_MOVING_WALL_ISOTHERMAL       11901   //      | CFDKit
#define BC_MOVING_WALL_HEATFLUX         11902   //      | CFDkit
#define BC_BURNING_SURFACE              11903   //      | CFDKit
//  Core
#define BC_WALL                         12000   // Core | CGNS
//  Invisid walls
#define BC_WALL_INVISCID                12001   //      | CGNS
//  Viscous walls
#define BC_WALL_VISCOUS                 12100   //      | CGNS
#define BC_NO_SLIP            BC_WALL_VISCOUS   //      | Duplicate
#define BC_WALL_VISCOUS_ISOTHERMAL      12101   //      | CGNS
#define BC_FIXED_TEMP_WALL BC_WALL_VISCOUS_ISOTHERMAL
                                                //      | Duplicate
#define BC_WALL_VISCOUS_HEATFLUX        12102   //      | CGNS
#define BC_ADIABATIC_WALL BC_WALL_VISCOUS_HEATFLUX
#define BC_FREE_SLIP                    12103   //      | ????                                

//----UNKNOWN OR WEIRD----//

#define LEFT_END_BOUNDARY               99901   // ?
#define RIGHT_END_BOUNDARY              99902   // ?
#define BC_FIXED_TEMP                   99903   // ? Dirichlet ?
#define BC_FIXED_HEATFLUX               99904   // ? Neumann ?
#define BC_RINGLEB_FLOW                 99905   // ?

/********************************************************
 * CFD - BC OPTIONS                                     *
 ********************************************************/
#define ADIABATIC_WALL                      1
#define FIXED_TEMPERATURE_WALL              2

/********************************************************
 * CFD -- Flow Types.                                   *
 ********************************************************/
#define FLOWTYPE_INVISCID                      0
#define FLOWTYPE_LAMINAR                       1
#define FLOWTYPE_TURBULENT_RANS_K_OMEGA        2
#define FLOWTYPE_TURBULENT_RANS_K_EPSILON      3
#define FLOWTYPE_TURBULENT_LES                 4
#define FLOWTYPE_TURBULENT_DES                 5
#define FLOWTYPE_TURBULENT_DNS                 6
#define FLOWTYPE_TURBULENT_DES_K_OMEGA        2
/********************************************************
 * CFD -- Initial Condition Types.                      *
 ********************************************************/

#define IC_RESTART                     -1

#define	IC_CONSTANT                     0
#define	IC_UNIFORM                      1
#define	IC_SOD                          2
#define	IC_SOD_XDIR                     2
#define	IC_SOD_YDIR                     3
#define	IC_GROTH                        4
#define	IC_GROTH_XDIR                   4
#define	IC_GROTH_YDIR                   5
#define	IC_EINFELDT                     6
#define	IC_EINFELDT_XDIR                6
#define	IC_EINFELDT_YDIR                7
#define IC_SHOCK_WAVE                   8
#define IC_SHOCK_WAVE_XDIR              8
#define IC_SHOCK_WAVE_YDIR              9
#define IC_CONTACT_SURFACE             10
#define IC_CONTACT_SURFACE_XDIR        10
#define IC_CONTACT_SURFACE_YDIR        11
#define IC_RAREFACTION_WAVE            12
#define IC_RAREFACTION_WAVE_XDIR       12
#define IC_RAREFACTION_WAVE_YDIR       13
#define	IC_SHOCK_BOX                   14
#define	IC_BRIO_WU                     15
#define	IC_HIGH_PRESSURE_RESERVOIR     16
#define	IC_LOW_PRESSURE_RESERVOIR      17
#define	IC_RIEMANN_XDIR                18
#define	IC_RIEMANN_YDIR                19
#define	IC_WEDGE_FLOW                  20
#define	IC_UNSTEADY_BLUNT_BODY         21
#define	IC_RINGLEB_FLOW                22
#define	IC_CYLINDRICAL_EXPLOSION       23
#define	IC_CYLINDRICAL_IMPLOSION       24

#define IC_DUSTY_COMPRESSION_XDIR      50
#define IC_DUSTY_COMPRESSION_YDIR      51
#define IC_DUSTY_EXPANSION_XDIR        52
#define IC_DUSTY_EXPANSION_YDIR        53
#define IC_DUSTY_CROSSING_JETS         54
#define IC_DUSTY_ELECTROSTATIC_CHANNEL 55
#define IC_DUSTY_SEALED_BOX            56

#define IC_GAS_MIX		       70
#define IC_CHEM_CORE_FLAME             71
#define IC_CHEM_INVERSE_FLAME          72
#define IC_CHEM_1DFLAME                73
#define IC_PRESSURE_GRADIENT_X         74
#define IC_PRESSURE_GRADIENT_Y         75

#define IC_VISCOUS_COUETTE                     80
#define IC_VISCOUS_COUETTE_ISOTHERMAL_WALL     81
#define IC_VISCOUS_COUETTE_PRESSURE_GRADIENT   82
#define IC_VISCOUS_POISEUILLE_FLOW             83
#define IC_VISCOUS_PIPE_FLOW                   83 // Duplicate of Poiseullie Flow
#define IC_VISCOUS_STOKES_FLOW                 84

#define IC_ROCKET_COLD_FLOW_STEADY_STATE       90
#define IC_ROCKET_COLD_FLOW_UNSTEADY           91

#define	IC_RIEMANN                    100
#define	IC_SQUARE_WAVE                101
#define	IC_SINX2_WAVE                 102
#define	IC_IMPULSIVE_ROD              103
#define	IC_SINUSOIDAL_ROD1            104
#define	IC_SINUSOIDAL_ROD4            105
#define	IC_RIEMANN_IVP_QX0            106
#define	IC_RIEMANN_IVP_T0             107
#define	IC_RIEMANN_IVP                108
#define IC_SQUARE_BOX_IVP             109
#define IC_CIRCULAR_BOX_IVP           110
#define IC_CIRCULAR_FLOW              111

#define IC_ELECTRIC_FIELD_UNIFORM     200
#define IC_ELECTRIC_FIELD_QUADRAPOLE  201
#define IC_ELECTRIC_FIELD_OCTAPOLE    202

/********************************************************
 * CFD -- Time Integration (Time-Stepping) Types.       *
 ********************************************************/

#define	TIME_STEPPING_EXPLICIT_EULER                     0
#define	TIME_STEPPING_IMPLICIT_EULER                     1
#define	TIME_STEPPING_IMPLICIT_TRAPEZOIDAL               2
#define	TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR       3
#define	TIME_STEPPING_SEMI_IMPLICIT_EULER                4
#define	TIME_STEPPING_SEMI_IMPLICIT_TRAPEZOIDAL          5
#define	TIME_STEPPING_SEMI_IMPLICIT_PREDICTOR_CORRECTOR  6
#define	TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING       7
#define	TIME_STEPPING_EXPLICIT_RUNGE_KUTTA               8

#define TIME_STEPPING_MULTIGRID                          9
#define TIME_STEPPING_NKS                               10

#define TIME_STEPPING_DUAL_TIME_STEPPING                11

#define	TIME_STEPPING_LAX_FRIEDRICHS                   100
#define	TIME_STEPPING_LAX_WENDROFF                     101
#define	TIME_STEPPING_MACCORMACK                       102
#define	TIME_STEPPING_HANCOCK                          103
#define	TIME_STEPPING_SIMPLE_EXPLICIT                  104
#define	TIME_STEPPING_SIMPLE_IMPLICIT                  105
#define	TIME_STEPPING_CRANK_NICOLSON                   106
#define	TIME_STEPPING_ADE                              107

/********************************************************
 * CFD -- Multigrid Cycle Types                         *
 ********************************************************/

#define MULTIGRID_V_CYCLE               1
#define MULTIGRID_W_CYCLE               2

/********************************************************
 * CFD -- Reconstruction Types.                         *
 ********************************************************/

#define	RECONSTRUCTION_MUSCL                     1
#define	RECONSTRUCTION_GREEN_GAUSS               2
#define	RECONSTRUCTION_LEAST_SQUARES             3
#define RECONSTRUCTION_CHARACTERISTIC            4
#define RECONSTRUCTION_DIAMOND_PATH              5
#define RECONSTRUCTION_INTERFACE_LEAST_SQUARES   6

#define VISCOUS_RECONSTRUCTION_ARITHMETIC       11
#define VISCOUS_RECONSTRUCTION_DIAMONDPATH      12
#define VISCOUS_RECONSTRUCTION_CARTESIAN        13
#define VISCOUS_RECONSTRUCTION_MEAN_GRADIENT    14

/********************************************************
 * CFD -- Limiter Types.                                *
 ********************************************************/

#define	LIMITER_ONE                    -1
#define	LIMITER_UNLIMITED              -1
#define	LIMITER_ZERO                    0
#define	LIMITER_MINMOD                  1
#define	LIMITER_UMIST                   2
#define	LIMITER_DOUBLE_MINMOD           3
#define	LIMITER_SUPERBEE                4
#define LIMITER_PHI                     5
#define LIMITER_VANLEER                 6
#define LIMITER_VANALBADA               7
#define LIMITER_BARTH_JESPERSEN         8
#define LIMITER_VENKATAKRISHNAN         9

/********************************************************
 * CFD -- Flux Function Types.                          *
 ********************************************************/

#define FLUX_FUNCTION_GODUNOV               1
#define	FLUX_FUNCTION_ROE                   2
#define	FLUX_FUNCTION_RUSANOV               3
#define	FLUX_FUNCTION_HLLE                  4
#define FLUX_FUNCTION_HLLL                  5
#define FLUX_FUNCTION_LINDE                 5 // Duplicate of HLLL
#define FLUX_FUNCTION_HLLC                  6
#define FLUX_FUNCTION_OSHER                 7
#define	FLUX_FUNCTION_ROE_PRECON_WS         8
#define	FLUX_FUNCTION_HLLE_PRECON_WS        9

#define DUSTY_FLUX_FUNCTION_SAUREL         11
#define DUSTY_FLUX_FUNCTION_EQUILIBRIUM    12

/********************************************************
 * CFD -- Local time-stepping/preconditioning.          *
 ********************************************************/

#define GLOBAL_TIME_STEPPING                          0
#define SCALAR_LOCAL_TIME_STEPPING                    1
#define MATRIX_LOCAL_TIME_STEPPING                    2
#define LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER    3
#define SEMI_IMPLICIT_LOCAL_TIME_STEPPING             4
#define SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER  5

/********************************************************
 * CFD -- Directions.                                   *
 ********************************************************/
 
#define CENTER                          0
#define NORTH                           1
#define SOUTH                           2
#define EAST                            3
#define WEST                            4
#define SNORTH                          5
#define SSOUTH                          6
#define SEAST                           7
#define SWEST                           8



/********************************************************
 * CFD -- Inline functions.                             *
 ********************************************************/

// minmod(x,y)
inline float  minmod(float x, float y)
   { return sgn(x)*max(float(ZERO),min(float(fabs(x)),sgn(x)*y)); }
inline double minmod(const double &x, const double &y)
   { return sgn(x)*max(ZERO,min(fabs(x),sgn(x)*y)); }

// minmod(x,y,z)
inline float  minmod(float x, float y, float z)
   { return sgn(x)*max(float(ZERO),min(float(fabs(x)),min(sgn(x)*y,sgn(x)*z))); }
inline double minmod(const double &x, const double &y, const double &z)
   { return sgn(x)*max(ZERO,min(fabs(x),min(sgn(x)*y,sgn(x)*z))); }

// minmod(x,y,z,w)
inline float  minmod(float x, float y, float z, float w)
   { return sgn(x)*max(float(ZERO),min(float(fabs(x)),min(sgn(x)*y,min(sgn(x)*z,sgn(x)*w)))); }
inline double minmod(const double &x, const double &y,
		     const double &z, const double &w)
   { return sgn(x)*max(ZERO,min(fabs(x),min(sgn(x)*y,min(sgn(x)*z,sgn(x)*w)))); }

// superbee(x,y)
inline float  superbee(float x, float y)
   { return sgn(x)*max(float(ZERO),max(min(float(TWO*fabs(x)),sgn(x)*y),
			           min(float(fabs(x)),float(TWO)*sgn(x)*y))); }
inline double superbee(const double &x, const double &y)
   { return sgn(x)*max(ZERO,max(min(TWO*fabs(x),sgn(x)*y),
			        min(fabs(x),TWO*sgn(x)*y))); }

// philimiter(x,y)
inline float  philimiter(float x, float y, float phi)
   { return sgn(x)*max(fabs(minmod(phi*x, y)), fabs(minmod(x, phi*y))); }
inline double philimiter(const double &x, const double &y,
			 const double &phi)
   { return sgn(x)*max(fabs(minmod(phi*x, y)), fabs(minmod(x, phi*y))); }

// vanleer(x,y)
inline float  vanleer(float x, float y)
   { return (fabs(x*y)+x*y)/(x+y+sgn(x+y)*sqr(TOLER)); }
inline double vanleer(const double &x, const double &y)
   { return (fabs(x*y)+x*y)/(x+y+sgn(x+y)*sqr(TOLER)); }

// vanalbada(x,y)
inline float  vanalbada(float x, float y, float epsi)
   { return (x*y+sqr(epsi))*(x+y)/(sqr(x)+sqr(y)+TWO*sqr(epsi)); }
inline double vanalbada(const double &x, const double &y,
		        const double &epsi)
   { return (x*y+sqr(epsi))*(x+y)/(sqr(x)+sqr(y)+TWO*sqr(epsi)); }

/********************************************************
 * CFD -- Define CFD structures and classes.            *
 ********************************************************/

/********************************************************
 * Class: CPUTime (CPU time for a process)              *
 *                                                      *
 * Member functions                                     *
 *      cpu_time0 -- Return CPU time used before last   *
 *                   update as a clock_t.               *
 *      cpu_time1 -- Return CPU time used at last       *
 *                   update as a clock_t.               *
 *      cput      -- Return total cput used by process. *
 *      reset     -- Resets CPU time counters.          *
 *      zero      -- Sets the CPU time to zero.         *
 *      update    -- Updates total CPU time for the     *
 *                   process.                           *
 *      sec       -- Returns total CPU time in seconds. *
 *      min       -- Returns total CPU time in minutes. *
 *      hrs       -- Returns total CPU time in hours.   *
 *                                                      *
 * Member operators                                     *
 *      T -- cpu time                                   *
 *                                                      *
 * T = T;                                               *
 * T = T + T;                                           *
 * T = T - T;                                           *
 * T = +T;                                              *
 * T = -T;                                              *
 * T += T;                                              *
 * T -= T;                                              *
 * T == T;                                              *
 * T != T;                                              *
 * cout << T; (output function)                         *
 * cin  >> T; (input function)                          *
 *                                                      *
 ********************************************************/
class CPUTime{
  private:
  public:
    clock_t cpu_time0,  // Track CPU time used so far as a clock_t.
            cpu_time1;  // To get the number of seconds used, 
                        // divide by CLOCKS_PER_SEC.
    double  cput;       // CPU time used by current process.
                        // Made public so can access them.

    /* Creation, copy, and assignment constructors. */
    CPUTime(void) {
       cpu_time0 = clock(); cpu_time1 = cpu_time0; 
       cput = ZERO;
    }

    CPUTime(const CPUTime &T) {
       cpu_time0 = T.cpu_time0; cpu_time1 = T.cpu_time1; 
       cput = T.cput; 
    }

    CPUTime(const double &cpu,
	    const clock_t &cpu0,
	    const clock_t &cpu1) {
       cpu_time0 = cpu0; cpu_time1 = cpu1; cput = cpu; 
    }

    /* Destructor. */
    // ~CPUTime(void);
    // Use automatically generated destructor.

    /* Reset CPU time counters. */
    void reset(void);

    /* Zero CPU time. */
    void zero(void);

    /* Update CPU time. */
    void update(void);

    /* CPU time in seconds. */
    double sec(void);
    double sec(void) const;

    /* CPU time in minutes. */
    double min(void);
    double min(void) const;

    /* CPU time in hours. */
    double hrs(void);
    double hrs(void) const;

    /* Assignment operator. */
    // CPUTime operator = (const CPUTime &T);
    // Use automatically generated assignment operator.

    /* Binary arithmetic operators. */
    friend CPUTime operator +(const CPUTime &T1, const CPUTime &T2);
    friend CPUTime operator -(const CPUTime &T1, const CPUTime &T2);
    
    /* Unary arithmetic operators. */
    friend CPUTime operator +(const CPUTime &T);
    friend CPUTime operator -(const CPUTime &T);

    /* Shortcut arithmetic operators. */
    friend CPUTime &operator +=(CPUTime &T1, const CPUTime &T2);
    friend CPUTime &operator -=(CPUTime &T1, const CPUTime &T2);
    
    /* Relational operators. */
    friend int operator ==(const CPUTime &T1, const CPUTime &T2);
    friend int operator !=(const CPUTime &T1, const CPUTime &T2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const CPUTime &T);
    friend istream &operator >> (istream &in_file,  CPUTime &T);

};

/********************************************************
 * CPUTime::reset -- Reset CPU time counters.           *
 ********************************************************/
inline void CPUTime::reset(void) {
    cpu_time0 = clock(); 
    cpu_time1 = cpu_time0; 
}

/********************************************************
 * CPUTime::zero -- Set total CPU time to zero.         *
 ********************************************************/
inline void CPUTime::zero(void) {
    cpu_time0 = clock(); 
    cpu_time1 = cpu_time0; 
    cput = ZERO;
}

/********************************************************
 * CPUTime::update -- Update total CPU time.            *
 ********************************************************/
inline void CPUTime::update(void) {
    cpu_time0 = cpu_time1;
    cpu_time1 = clock();
    cput = cput + double(cpu_time1-cpu_time0)/double(CLOCKS_PER_SEC);
}

/********************************************************
 * CPUTime::sec -- CPU time in seconds.                 *
 ********************************************************/
inline double CPUTime::sec(void) {
    return (cput);
}

inline double CPUTime::sec(void) const {
    return (cput);
}

/********************************************************
 * CPUTime::min -- CPU time in minutes.                 *
 ********************************************************/
inline double CPUTime::min(void) {
    return (cput/SIXTY);
}

inline double CPUTime::min(void) const {
    return (cput/SIXTY);
}

/********************************************************
 * CPUTime::hrs -- CPU time in hours.                   *
 ********************************************************/
inline double CPUTime::hrs(void) {
    return (cput/sqr(SIXTY));
}

inline double CPUTime::hrs(void) const {
    return (cput/sqr(SIXTY));
}

/********************************************************
 * CPUTime -- Binary arithmetic operators.              *
 ********************************************************/
inline CPUTime operator +(const CPUTime &T1, const CPUTime &T2) {
  return (CPUTime(T1.cput+T2.cput,T1.cpu_time0,T1.cpu_time1));
}

inline CPUTime operator -(const CPUTime &T1, const CPUTime &T2) {
  return (CPUTime(T1.cput-T2.cput,T1.cpu_time0,T1.cpu_time1));
}

/********************************************************
 * CPUTime -- Unary arithmetic operators.               *
 ********************************************************/
inline CPUTime operator +(const CPUTime &T) {
  return (CPUTime(T.cput,T.cpu_time0,T.cpu_time1));
}

inline CPUTime operator -(const CPUTime &T) {
  return (CPUTime(-T.cput,T.cpu_time0,T.cpu_time1));
}

/********************************************************
 * CPUTime -- Shortcut arithmetic operators.            *
 ********************************************************/
inline CPUTime &operator +=(CPUTime &T1, const CPUTime &T2) {
  T1.cput += T2.cput;
  T1.cpu_time0 = T1.cpu_time0;
  T1.cpu_time1 = T1.cpu_time1;
  return (T1);
}

inline CPUTime &operator -=(CPUTime &T1, const CPUTime &T2) {
  T1.cput -= T2.cput;
  T1.cpu_time0 = T1.cpu_time0;
  T1.cpu_time1 = T1.cpu_time1;
  return (T1);
}

/********************************************************
 * CPUTime -- Relational operators.                     *
 ********************************************************/
inline int operator ==(const CPUTime &T1, const CPUTime &T2) {
  return (T1.cput == T2.cput);
}

inline int operator !=(const CPUTime &T1, const CPUTime &T2) {
  return (T1.cput != T2.cput);
}

/********************************************************
 * CPUTime -- Input-output operators.                   *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const CPUTime &T) {
  out_file.setf(ios::scientific);
  out_file << " " << T.cput;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, CPUTime &T) {
  in_file.setf(ios::skipws);
  in_file >> T.cput;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

#define	INPUT_PARAMETER_LENGTH_CFD1D    80

/********************************************************
 * Class:  CFD1D_Input_Parameters                       *
 ********************************************************/
class CFD1D_Input_Parameters{
  private:
  public:
  // Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_CFD1D];

  // Input file stream:
  ifstream Input_File;

  // Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
      Maximum_Number_of_Time_Steps, N_Stage;
  double CFL_Number, Time_Max;

  // Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Reconstruction;

  // Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Limiter;

  // Flux_Function_Type and related input parameters
  char Flux_Function_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Flux_Function;

  // Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_ICs;

  // Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Grid;
  int Number_of_Cells, Number_of_Nodes;
  double X_Min, X_Max;

  // Diffusion coefficient, wave speed, and relaxation time:
  double Kappa, a, Tau;

  // Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_CFD1D];

  // Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_CFD1D];

  // Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_CFD1D];

  // Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Output_Format;

  // Input file line number:
  int Line_Number;

  /* Input-output operators. */

  friend ostream &operator << (ostream &out_file,
		               const CFD1D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       CFD1D_Input_Parameters &IP);

};

/*************************************************************
 * CFD1D_Input_Parameters -- Input-output operators.         *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const CFD1D_Input_Parameters &IP) {
    if (IP.i_Grid == GRID_CARTESIAN_UNIFORM) {
       out_file << "\n\n Solving PDE(s) on 1D uniform Cartesian mesh.";
    } else {
       out_file << "\n\n Solving PDE(s) on 1D mesh.";
    } /* endif */
    out_file << "\n  -> Input File Name: " 
             << IP.Input_File_Name;
    if (IP.Time_Accurate) { 
       out_file << "\n  -> Time Accurate (Unsteady) Solution";
    } else {
       out_file << "\n  -> Time Invariant (Steady-State) Solution";
    }
    out_file << "\n  -> Time Integration: " 
             << IP.Time_Integration_Type;
    out_file << "\n  -> Number of Stages in Multi-Stage Scheme: " 
             << IP.N_Stage;
    if (IP.Local_Time_Stepping) 
       out_file << "\n  -> Local Time Stepping";
    out_file << "\n  -> Reconstruction: " 
             << IP.Reconstruction_Type;
    out_file << "\n  -> Limiter: " 
             << IP.Limiter_Type;
    out_file << "\n  -> Flux Function: " 
             << IP.Flux_Function_Type;
    out_file << "\n  -> Initial Conditions: " 
             << IP.ICs_Type;
    switch(IP.i_ICs) {
      case IC_CONSTANT :
        break;
      case IC_UNIFORM :
        break;
      case IC_SOD :
        break;
      case IC_GROTH :
        break;
      case IC_EINFELDT :
        break;
      default:
        break;
    } /* endswitch */
    out_file << "\n  -> Diffusion Coefficient : " 
             << IP.Kappa;
    out_file << "\n  -> Wave Speed : " 
             << IP.a;
    out_file << "\n  -> Relaxation Time : " 
             << IP.Tau;
    out_file << "\n  -> Grid: " 
             << IP.Grid_Type;
    switch(IP.i_Grid) {
      case GRID_CARTESIAN_UNIFORM :
        out_file << "\n  -> Minimum X Location : " 
                 << IP.X_Min;
        out_file << "\n  -> Maximum X Location : " 
                 << IP.X_Max;
        out_file << "\n  -> Width of Solution Domain : " 
                 << IP.X_Max-IP.X_Min;
        break;
      default:
        out_file << "\n  -> Minimum X Location : " 
                 << IP.X_Min;
        out_file << "\n  -> Maximum X Location : " 
                 << IP.X_Max;
        out_file << "\n  -> Width of Solution Domain : " 
                 << IP.X_Max-IP.X_Min;
        break;
    } /* endswitch */
    out_file << "\n  -> Number of Cells: "
             << IP.Number_of_Cells;
    out_file << "\n  -> Number of Nodes: " 
             << IP.Number_of_Nodes;
    out_file << "\n  -> CFL Number: " 
             << IP.CFL_Number;
    out_file << "\n  -> Maximum Time : " 
             << IP.Time_Max;
    out_file << "\n  -> Maximum Numer of Time Steps: " 
             << IP.Maximum_Number_of_Time_Steps;
    out_file << "\n  -> Output File Name: " 
             << IP.Output_File_Name;
    out_file << "\n  -> Output Format: " 
             << IP.Output_Format_Type;
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     CFD1D_Input_Parameters &IP) {
    return (in_file);
}

/*************************************************************
 * CFD1D_Input_Parameters -- External subroutines.           *
 *************************************************************/
extern void Open_Input_File(CFD1D_Input_Parameters &IP);

extern void Close_Input_File(CFD1D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(CFD1D_Input_Parameters &IP);

extern void Get_Next_Input_Control_Parameter(CFD1D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(CFD1D_Input_Parameters &IP);

/********************************************************
 * CFD -- External subroutines.                         *
 ********************************************************/

extern void Output_Progress(const int Number_of_Time_Steps,
                            const double &Time,
                            const CPUTime &CPU_Time,
                            const double &Residual_L1_Norm,
                            const int First_Step,
                            const int Frequency);

extern void Output_Progress_L2norm(const int Number_of_Time_Steps,
				   const double &Time,
				   const CPUTime &CPU_Time,
				   const double &Residual_L2_Norm,
				   const int First_Step,
				   const int Frequency);

extern int Open_Progress_File(ofstream &Progress_File,
                              char *File_Name,
                              const int Append_to_File);

extern int Close_Progress_File(ofstream &Progress_File);

extern void Output_Progress_to_File(ostream &Progress_File,
                                    const int Number_of_Time_Steps,
                                    const double &Time,
                                    const CPUTime &CPU_Time,
                                    const double &Residual_L1_Norm,
                                    const double &Residual_L2_Norm,
                                    const double &Residual_Max_Norm);

extern double StretchingFcn(const double &xx,
                            const double &beta,
	      	            const double &tau,
                            const int i_stretch);

extern double HartenFixPos(const double &lambda_a,
                           const double &lambda_l,
	      	           const double &lambda_r);

extern double HartenFixNeg(const double &lambda_a,
                           const double &lambda_l,
	      	           const double &lambda_r);

extern double HartenFixAbs(const double &lambda_a,
                           const double &lambda_l,
	      	           const double &lambda_r);

extern double Limiter_BarthJespersen(double *uQuad,
                                     const double &u0,
                                     const double &u0Min,
	      	                     const double &u0Max,
			             const int nQuad);

extern double Limiter_Venkatakrishnan(double *uQuad,
                                      const double &u0,
                                      const double &u0Min,
	      	                      const double &u0Max,
			              const int nQuad);

extern double Limiter_VanLeer(double *uQuad,
                              const double &u0,
                              const double &u0Min,
	      	              const double &u0Max,
			      const int nQuad);

extern double Limiter_VanAlbada(double *uQuad,
                                const double &u0,
                                const double &u0Min,
	      	                const double &u0Max,
			        const int nQuad);

extern double Runge_Kutta(const int I_stage,
                          const int N_stage);

extern double MultiStage_Optimally_Smoothing(const int I_stage,
                                             const int N_stage,
                                             const int Limiter_Type);

/********************************************************
 * Routine: Bilinear_Interpolation                      *
 *                                                      *
 * Given 4 pts P1,P2,P3,P4 that form an arbitrary       *
 * quadrilateral, and their associated solution         *
 * variables U1,U2,U3,U4, and a point P0 known to lie   *
 * inside this quadrilateral, this fxn calculates the   *
 * interpolated solution variable U0, and returns an    *
 * integer error code if P0 falls outside of the quad.  *
 *                                                      *
 ********************************************************/
template <class T>
int Bilinear_Interpolation(const T U1, const Vector2D P1,
			   const T U2, const Vector2D P2,
			   const T U3, const Vector2D P3,
			   const T U4, const Vector2D P4,
			   const Vector2D P0, T &U0) {
  double a0,a1,a2,a3;
  double b0,b1,b2,b3;
  double zeta0, eta0;
  double a,b,c;

  a0 = P1.x;
  a1 = P4.x - P1.x;
  a2 = P2.x - P1.x;
  a3 = P1.x + P3.x - P2.x - P4.x;

  b0 = P1.y;
  b1 = P4.y - P1.y;
  b2 = P2.y - P1.y;
  b3 = P1.y + P3.y - P2.y - P4.y;

  a = a1*b3 - a3*b1;
  b = a0*b3 - b3*P0.x + a1*b2 - a2*b1 + a3*P0.y - a3*b0;
  c = a0*b2 - b2*P0.x + a2*P0.y - a2*b0;

  if (a > TOLER || a < -TOLER) {
    // Check for complex roots
    if (b*b-4*a*c < 0) return(1);
    zeta0 = (-b + sqrt(b*b - 4*a*c))/(2*a);
    // Check if zeta0 is within range, if not calculate conjugate
    if (zeta0 < 0 || zeta0 > 1) {
      zeta0 = (-b - sqrt(b*b - 4*a*c))/(2*a);
      // Check if zeta0 is within range
      if (zeta0 < 0 || zeta0 > 1) return(1);
    }
  }
  else if (b > TOLER || b < -TOLER) { // if (a == 0 AND b != 0)
    zeta0 = -c/b;
    // Check if zeta0 is within range
    if (zeta0 < 0 || zeta0 > 1) return(1);
  }
  else { // if (a == 0 AND b == 0)
    return(1);
  }

  // check for divide by zero error
  if (b2 + b3*zeta0 < TOLER && b2 + b3*zeta0 > -TOLER) return(1);

  eta0 = (P0.y - b0 - b1*zeta0)/(b2 + b3*zeta0);

  U0 = U1 + (U4 - U1)*zeta0 + (U2 - U1)*eta0 + (U1 + U3 - U2 - U4)*zeta0*eta0;
  return(0);
}

/**********************************************************************
 * Routine: Green_Gauss_Integration                                   *
 *                                                                    *
 * Given four pts X1, X2, X3, and X4 that form an arbitrary           *
 * quadrilateral, and their associated solution variables U1, U2, U3, *
 * and U4, this function calculates the gradient of the solution      *
 * variables at the centre of the quadrilateral by conducting a       *
 * Green-Gauss reconstruction on the specified path.                  *
 *                                                                    *
 **********************************************************************/
template <class T>
int Green_Gauss_Integration(const Vector2D X1, const T U1,
			    const Vector2D X2, const T U2,
			    const Vector2D X3, const T U3,
			    const Vector2D X4, const T U4,
			    T &dUdx, T &dUdy) {

  int error_flag;
  double l21, l32, l43, l14, A;
  Vector2D n21, n32, n43, n14;
  T U;

  // Determine the lengths and normals of the faces and the area
  // of the region of Green-Gauss integration.
  l21 = abs(X2-X1);  n21 = Vector2D((X2.y-X1.y),-(X2.x-X1.x))/abs(X2-X1);
  l32 = abs(X3-X2);  n32 = Vector2D((X3.y-X2.y),-(X3.x-X2.x))/abs(X3-X2);
  l43 = abs(X4-X3);  n43 = Vector2D((X4.y-X3.y),-(X4.x-X3.x))/abs(X4-X3);
  l14 = abs(X1-X4);  n14 = Vector2D((X1.y-X4.y),-(X1.x-X4.x))/abs(X1-X4);
  A = HALF*(((X2-X1)^(X4-X1)) + ((X3-X4)^(X3-X2)));
  // Compute Green-Gauss integration to determine the gradient of the 
  // solution state at the centre of the quadrilateral.
  U = HALF*(U1+U2)*l21;
  dUdx = U*n21.x;
  dUdy = U*n21.y;
  U = HALF*(U2+U3)*l32;
  dUdx += U*n32.x;
  dUdy += U*n32.y;
  U = HALF*(U3+U4)*l43;
  dUdx += U*n43.x;
  dUdy += U*n43.y;
  U = HALF*(U4+U1)*l14;
  dUdx += U*n14.x;
  dUdy += U*n14.y;
  dUdx = dUdx/A;
  dUdy = dUdy/A;

  // Green-Gauss integration successfully applied.
  return 0;

}

#endif /* _CFD_INCLUDED  */
