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

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

/* Include Vector2D header. */

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

/* Include Tensor2D header. */

#ifndef _TENSOR2D_INCLUDED
#include "../Math/Tensor2D.h"
#endif //_TENSOR2D_INCLUDED

/* Include Polygon header. */

#ifndef _POLYGON_INCLUDED
#include "../Math/Polygon.h"
#endif //_POLYGON_INCLUDED

/* Include Vector3D header. */

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif //_VECTOR3D_INCLUDED

/* Include Tensor3D header. */

#ifndef _TENSOR3D_INCLUDED
#include "../Math/Tensor3D.h"
#endif //_TENSOR3D_INCLUDED

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

// Cantera Version.
inline string Cantera_Version() {
#ifdef _CANTERA_VERSION
#ifdef _CANTERA_V70
  return ("Cantera Version, Built using Cantera Release 1.7, California Institute of Technology, 2001-2006.");
#else
  return ("Cantera Version, Built using Cantera, California Institute of Technology.");
#endif // _CANTERA_V70
#else
  return ("No Cantera Version.");
#endif //_CANTERA_VERSION
}

/**********************************************************************
 * CFD -- Date and time.                                              *
 **********************************************************************/

inline char *Date_And_Time() {
  char *string_ptr; time_t epoch_time;
  epoch_time = time(NULL);
  string_ptr = ctime(&epoch_time);
  string_ptr[strlen(string_ptr)-1] = '\0';
  return (string_ptr);
}

/**********************************************************************
 * CFD -- General Purpose Command Codes.                              *
 **********************************************************************/

#define	OFF                           0
#define	ON                            1

/**********************************************************************
 * CFD -- Input parameter command codes.                              *
 **********************************************************************/

#define EXECUTE_CODE                                     10000
#define TERMINATE_CODE                                   10001
#define CONTINUE_CODE                                    10002
#define WRITE_OUTPUT_CODE                                10003
#define WRITE_OUTPUT_CELLS_CODE                          10004
#define WRITE_OUTPUT_NODES_CODE                          10005
#define WRITE_OUTPUT_GRADIENTS_CODE                      10006
#define WRITE_OUTPUT_QUASI3D_CODE                        10007
#define WRITE_RESTART_CODE                               10008
#define READ_RESTART_CODE                                10009
#define WRITE_OUTPUT_GRID_CODE                           10010
#define WRITE_GRID_DEFINITION_CODE                       10011
#define WRITE_OUTPUT_GRID_NODES_CODE                     10012
#define WRITE_OUTPUT_GRID_CELLS_CODE                     10013
#define COMMENT_CODE                                     10014
#define REFINE_GRID_CODE                                 10015
#define BOUNDING_BOX_REFINE_GRID_CODE                    10016
#define MORTON_ORDERING_CODE                             10017

#define WRITE_OUTPUT_ELEMENTS_CODE                       10020
#define WRITE_OUTPUT_CELL_STATUS_CODE                    10022
#define WRITE_OUTPUT_INTERFACE_COMPONENT_LIST_CODE       10023
#define WRITE_OUTPUT_INTERFACE_UNION_LIST_CODE           10024

#define WRITE_OUTPUT_LEVEL_SET_CODE                      10030
#define WRITE_OUTPUT_LEVEL_SET_CELLS_CODE                10031
#define WRITE_OUTPUT_LEVEL_SET_INTERFACE_LIST_CODE       10032
#define WRITE_OUTPUT_LEVEL_SET_CIRCLE_CODE               10033
#define WRITE_OUTPUT_LEVEL_SET_ELLIPSE_CODE              10034
#define WRITE_OUTPUT_LEVEL_SET_ZALESAK_CODE              10035

#define WRITE_OUTPUT_RHS_CODE                            10040
#define WRITE_OUTPUT_PERTURB_CODE                        10041

#define WRITE_OUTPUT_RINGLEB_CODE                        10050
#define WRITE_OUTPUT_VISCOUS_CHANNEL_CODE                10051
#define WRITE_OUTPUT_VISCOUS_PIPE_CODE                   10052
#define WRITE_OUTPUT_TURBULENT_PIPE_CODE                 10053
#define WRITE_OUTPUT_FLAT_PLATE_CODE                     10054
#define WRITE_OUTPUT_DRIVEN_CAVITY_FLOW_CODE             10055
#define WRITE_OUTPUT_BACKWARD_FACING_STEP_CODE           10056
#define WRITE_OUTPUT_MIXING_LAYER_CODE                   10057

#define DETERMINE_MAXIMUM_SURFACE_PRESSURE_CODE          10070
#define DETERMINE_MACH_STEM_HEIGHT_CODE                  10071
#define WRITE_OUTPUT_WEDGE_SOLUTION_DISTRIBUTION_CODE    10072

#define WRITE_OUTPUT_ELLIPTIC_OPERATOR_ANALYSIS          10080

#define WRITE_OUTPUT_AERODYNAMIC_COEFFICIENTS_CODE       10090

#define	INVALID_INPUT_CODE                              -10000
#define	INVALID_INPUT_VALUE                             -10001

/**********************************************************************
 * CFD -- I/O Types.                                                  *
 **********************************************************************/

#define	IO_GNUPLOT       1
#define	IO_TECPLOT       2
#define	IO_MATLAB        3
#define	IO_OCTAVE        4

/**********************************************************************
 * CFD -- Types of algebraic grid point distribution stretching       *
 *        functions.                                                  *
 **********************************************************************/

#define STRETCHING_FCN_LINEAR            0
#define STRETCHING_FCN_MIN_CLUSTERING    1
#define STRETCHING_FCN_MAX_CLUSTERING    2
#define STRETCHING_FCN_MINMAX_CLUSTERING 3
#define STRETCHING_FCN_MIDPT_CLUSTERING  4
#define STRETCHING_FCN_SINE              5
#define STRETCHING_FCN_COSINE            6

/**********************************************************************
 * CFD -- Grid Types.                                                 *
 **********************************************************************/

#define GRID_CARTESIAN_UNIFORM                0
#define GRID_SQUARE                           1
#define GRID_RECTANGULAR_BOX                  2
#define GRID_FLAT_PLATE                       3
#define GRID_PIPE                             4
#define GRID_BLUNT_BODY                       5
#define GRID_ROCKET_MOTOR                     6
#define GRID_CIRCULAR_CYLINDER                7
#define GRID_ELLIPSE                          8
#define GRID_NACA_AEROFOIL                    9
#define GRID_NASA_ROTOR_37                   10
#define GRID_NASA_ROTOR_67                   11
#define GRID_FREE_JET                        12
#define GRID_WEDGE                           13
#define GRID_UNSTEADY_BLUNT_BODY             14
#define GRID_RINGLEB_FLOW                    15
#define GRID_BUMP_CHANNEL_FLOW               16
#define GRID_DRIVEN_CAVITY_FLOW              17
#define GRID_BACKWARD_FACING_STEP            18
#define GRID_MIXING_LAYER                    19 
#define GRID_NOZZLE                          20
#define GRID_NOZZLELESS_ROCKET_MOTOR         21
#define GRID_DESOLVATION_CHAMBER             22
#define GRID_CUBE                            23
#define GRID_CHANNEL                         24
#define GRID_CHANNEL_ZDIR          GRID_CHANNEL
#define GRID_CHANNEL_XDIR                    25
#define GRID_CHANNEL_YDIR                    26
#define GRID_COUETTE                         27
#define GRID_COUETTE_ZDIR          GRID_COUETTE
#define GRID_COUETTE_XDIR                    28
#define GRID_COUETTE_YDIR                    29
#define GRID_BLUFF_BODY_BURNER               30

#define GRID_ICEMCFD                       1000
#define GRID_READ_FROM_DEFINITION_FILE    10000
#define GRID_READ_FROM_GRID_DATA_FILE     10001

#define GRID_1DFLAME                         51
#define GRID_LAMINAR_FLAME                   52

/**********************************************************************
 * CFD -- Nozzle Types.                                               *
 **********************************************************************/

#define NOZZLE_CONICAL                     1770
#define NOZZLE_GOTTLIEB_FUNCTION           1771
#define NOZZLE_QUARTIC_FUNCTION            1772
#define NOZZLE_HYBRID_CONICAL_FUNCTION     1773

/**********************************************************************
 * CFD -- Boundary Condition Types.                                   *
 **********************************************************************/

//----DEFAULT NAMES----//                          Core | Scope

//--Undefined (or interior)
#define	BC_NULL                             0   // Core | CGNS
#define BC_NONE                       BC_NULL   //      | Duplicate

//--Speciality interior
#define BC_INTERIOR                      1000   // Core | CFFC

//--Collapsed blocks
#define BC_DEGENERATE                    2000   // Core | Header
#define BC_DEGENERATE_LINE               2001   //      | CGNS
#define BC_DEGENERATE_POINT              2002   //      | CGNS
#define BC_AXISYMMETRIC_WEDGE            2003   //      | CGNS

//--General?
#define BC_GENERAL                       3000   // Core | CGNS

//--Extrapolation
#define BC_EXTRAPOLATE                   4000   // Core | CGNS
#define	BC_CONSTANT_EXTRAPOLATION        4001   //      | CFFC
#define	BC_LINEAR_EXTRAPOLATION          4002   //      | CFFC

//--1st Derivative
#define BC_NEUMANN                       5000   // Core | CGNS

//--Periodic
#define BC_PERIODIC                      6000   // Core | CFFC

//--Symmetry plane
#define BC_SYMMETRY                      7000   // Core | Header
#define BC_SYMMETRY_PLANE                7001   //      | CGNS
#define BC_REFLECTION       BC_SYMMETRY_PLANE   //      | Duplicate
#define BC_SYMMETRY_POLAR                7002   //      | CGNS

//--Far field
#define BC_FARFIELD                      8000   // Core | CGNS
#define BC_CHARACTERISTIC                8001   //      | CFFC ?
#define BC_FIXED_PRESSURE                8002   //      | CFFC
#define BC_CHARACTERISTIC_VELOCITY       8003   //      | CFFC

//--Inflow
#define BC_INFLOW                        9000   // Core | CGNS
#define BC_INFLOW_SUBSONIC               9001   //      | CGNS
#define BC_INFLOW_SUPERSONIC             9002   //      | CGNS
#define BC_TUNNEL_INFLOW                 9003   //      | CGNS
#define BC_FLAME_INFLOW                  9003   //      | CFFC 
#define	BC_CHANNEL_INFLOW                9004   //      | CFFC
#define BC_INFLOW_EXTRAPOLATION          9005   //      | CFFC 

//--Outflow
#define BC_OUTFLOW                      10000   // Core | CGNS
#define BC_OPEN_END                BC_OUTFLOW   //      | Duplicate ?
#define BC_OUTFLOW_SUBSONIC             10001   //      | CGNS
#define BC_OUTFLOW_SUPERSONIC           10002   //      | CGNS
#define BC_TUNNEL_OUTFLOW               10003   //      | CGNS
#define BC_FLAME_OUTFLOW                10004   //      | CFFC
#define	BC_CHANNEL_OUTFLOW              10005   //      | CFFC
#define BC_OUTFLOW_EXTRAPOLATION        10006   //      | CFFC

//--Frozen
#define BC_ROBIN                        10999   //      | CFFC
#define BC_DIRICHLET                    11000   // Core | CGNS
#define BC_FIXED                 BC_DIRICHLET   //      | Duplicate

//--Wall
//  Speciality wall
#define BC_ABSORPTION                   11800   //      | CFFC
//  Moving wall
#define BC_MOVING_WALL                  11900   //      | CFFC
#define BC_MOVING_WALL_ISOTHERMAL       11901   //      | CFFC
#define BC_MOVING_WALL_HEATFLUX         11902   //      | CFFC
#define BC_BURNING_SURFACE              11903   //      | CFFC
#define BC_MASS_INJECTION               11904   //      | CFFC
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

/**********************************************************************
 * CFD - BC OPTIONS                                                   *
 **********************************************************************/

#define ADIABATIC_WALL                      1
#define FIXED_TEMPERATURE_WALL              2

/**********************************************************************
 * CFD -- Boundary Condition Data Types.                              *
 **********************************************************************/

#define BCS_SIZE                                        5
#define NUM_OF_DATA_TYPES   			        68
#define BC_DATA_UNDEFINED			        0
#define BC_DATA_POTENTIAL			        1
#define BC_DATA_STREAMFUNCTION			        2
#define BC_DATA_DENSITY			                3 
#define BC_DATA_PRESSURE			        4
#define BC_DATA_TEMPERATURE			        5
#define BC_DATA_ENERGYINTERNAL			        6
#define BC_DATA_ENTHALPY			        7
#define BC_DATA_ENTROPY			                8
#define BC_DATA_ENTROPYAPPROX			        9
#define BC_DATA_DENSITYSTAGNATION			10
#define BC_DATA_PRESSURESTAGNATION			11
#define BC_DATA_TEMPERATURESTAGNATION			12
#define BC_DATA_ENERGYSTAGNATION			13
#define BC_DATA_ENTHALPYSTAGNATION			14
#define BC_DATA_ENERGYSTAGNATIONDENSITY			15
#define BC_DATA_VELOCITYX			        16
#define BC_DATA_VELOCITYY			        17
#define BC_DATA_VELOCITYZ			        18
#define BC_DATA_VELOCITYR			        19
#define BC_DATA_VELOCITYTHETA			        20
#define BC_DATA_VELOCITYPHI			        21
#define BC_DATA_VELOCITYMAGNITUDE			22
#define BC_DATA_VELOCITYNORMAL			        23
#define BC_DATA_VELOCITYTANGENTIAL			24
#define BC_DATA_VELOCITYSOUND			        25
#define BC_DATA_VELOCITYSOUNDSTAGNATION			26
#define BC_DATA_MOMENTUMX			        27
#define BC_DATA_MOMENTUMY			        28
#define BC_DATA_MOMENTUMZ			        29
#define BC_DATA_MOMENTUMMAGNITUDE			30
#define BC_DATA_ENERGYKINETIC			        31
#define BC_DATA_PRESSUREDYNAMIC			        32
#define BC_DATA_VORTICITYX			        33
#define BC_DATA_VORTICITYY			        34
#define BC_DATA_VORTICITYZ			        35
#define BC_DATA_VORTICITYMAGNITUDE			36
#define BC_DATA_SKINFRICTIONX			        37
#define BC_DATA_SKINFRICTIONY			        38
#define BC_DATA_SKINFRICTIONZ			        39
#define BC_DATA_SKINFRICTIONMAGNITUDE			40
#define BC_DATA_VELOCITYANGLEX			        41
#define BC_DATA_VELOCITYANGLEY			        42
#define BC_DATA_VELOCITYANGLEZ			        43
#define BC_DATA_VELOCITYUNITVECTORX			44
#define BC_DATA_VELOCITYUNITVECTORY			45
#define BC_DATA_VELOCITYUNITVECTORZ			46
#define BC_DATA_MASSFLOW			        47
#define BC_DATA_VISCOSITYKINEMATIC			48
#define BC_DATA_VISCOSITYMOLECULAR			49
#define BC_DATA_VISCOSITYEDDYKINEMATIC			50
#define BC_DATA_VISCOSITYEDDY			        51 
#define BC_DATA_THERMALCONDUCTIVITY			52
#define BC_DATA_POWERLAWEXPONENT			53
#define BC_DATA_SUTHERLANDLAWCONSTANT			54
#define BC_DATA_TEMPERATUREREFERENCE			55
#define BC_DATA_VISCOSITYMOLECULARREFERENCE		56
#define BC_DATA_THERMALCONDUCTIVITYREFERENCE		57
#define BC_DATA_IDEALGASCONSTANT			58
#define BC_DATA_SPECIFICHEATPRESSURE			59
#define BC_DATA_SPECIFICHEATVOLUME			60
#define BC_DATA_REYNOLDSSTRESSXX			61
#define BC_DATA_REYNOLDSSTRESSXY			62
#define BC_DATA_REYNOLDSSTRESSXZ			63
#define BC_DATA_REYNOLDSSTRESSYY			64
#define BC_DATA_REYNOLDSSTRESSYZ			65
#define BC_DATA_REYNOLDSSTRESSZZ			66
#define BC_DATA_LENGTHREFERENCE			        67

/**********************************************************************
 * CFD -- Flow Types.                                                 *
 **********************************************************************/

#define FLOWTYPE_INVISCID                              0
#define FLOWTYPE_LAMINAR                               1
#define FLOWTYPE_TURBULENT_RANS_SPALART_ALLMARAS       2
#define FLOWTYPE_TURBULENT_RANS_K_EPSILON              3
#define FLOWTYPE_TURBULENT_RANS_WOLFSTEIN_K_EPSILON    4
#define FLOWTYPE_TURBULENT_RANS_TWO_LAYER_K_EPSILON    5
#define FLOWTYPE_TURBULENT_RANS_K_OMEGA                6
#define FLOWTYPE_TURBULENT_RANS_STRESS_OMEGA           7
#define FLOWTYPE_TURBULENT_LES                         8
#define FLOWTYPE_TURBULENT_DES                         9
#define FLOWTYPE_TURBULENT_DES_K_OMEGA                10
#define FLOWTYPE_TURBULENT_DNS                        11

/**********************************************************************
 * CFD -- Particle-phase formulation.                                 *
 **********************************************************************/

#define PARTICLE_PHASE_NONE                              0
#define PARTICLE_PHASE_EULERIAN_FORMULATION              1
#define PARTICLE_PHASE_LAGRANGIAN_FORMULATION            2

/**********************************************************************
 * CFD -- Turbulent flow boundary condition types.                    *
 **********************************************************************/

#define TURBULENT_BC_DIRECT_INTEGRATION           950
#define TURBULENT_BC_STANDARD_WALL_FUNCTION       951
#define TURBULENT_BC_TWO_LAYER_WALL_FUNCTION      952
#define TURBULENT_BC_THREE_LAYER_WALL_FUNCTION    953
#define TURBULENT_BC_AUTOMATIC_WALL_TREATMENT     954

/**********************************************************************
 * CFD -- Turbulent flow friction velocity evaluation types.          *
 **********************************************************************/

#define FRICTION_VELOCITY_LOCAL_SHEAR_STRESS      970
#define FRICTION_VELOCITY_WALL_SHEAR_STRESS       971
#define FRICTION_VELOCITY_ITERATIVE               972
#define FRICTION_VELOCITY_PIPE                    973

/**********************************************************************
 * CFD -- Initial Condition Types.                                    *
 **********************************************************************/

#define IC_RESTART                     -1

#define	IC_CONSTANT                     0
#define	IC_UNIFORM                      1
#define	IC_SOD                          2
#define	IC_SOD_XDIR                     2
#define	IC_SOD_YDIR                     3
#define IC_SOD_ZDIR                     4
#define	IC_SHOCK_BOX                    5
#define	IC_GROTH                        6
#define	IC_GROTH_XDIR                   6
#define	IC_GROTH_YDIR                   7
#define	IC_GROTH_ZDIR                   8
#define	IC_EINFELDT                     9
#define	IC_EINFELDT_XDIR                9
#define	IC_EINFELDT_YDIR               10
#define	IC_EINFELDT_ZDIR               11
#define IC_SHOCK_WAVE                  12
#define IC_SHOCK_WAVE_XDIR             12
#define IC_SHOCK_WAVE_YDIR             13
#define IC_SHOCK_WAVE_ZDIR             14
#define IC_CONTACT_SURFACE             15
#define IC_CONTACT_SURFACE_XDIR        15
#define IC_CONTACT_SURFACE_YDIR        16
#define IC_CONTACT_SURFACE_ZDIR        17
#define IC_RAREFACTION_WAVE            18
#define IC_RAREFACTION_WAVE_XDIR       18
#define IC_RAREFACTION_WAVE_YDIR       19
#define IC_RAREFACTION_WAVE_ZDIR       20
#define	IC_BRIO_WU                     21
#define	IC_HIGH_PRESSURE_RESERVOIR     22
#define	IC_LOW_PRESSURE_RESERVOIR      23
#define	IC_RIEMANN_XDIR                24
#define	IC_RIEMANN_YDIR                25
#define	IC_RIEMANN_ZDIR                26
#define	IC_WEDGE_FLOW                  27
#define	IC_UNSTEADY_BLUNT_BODY         28
#define	IC_RINGLEB_FLOW                29
#define	IC_CYLINDRICAL_EXPLOSION       30
#define	IC_CYLINDRICAL_IMPLOSION       31

#define IC_SQUARE_WAVE_XDIR            40
#define IC_SQUARE_WAVE_YDIR            41
#define IC_SQUARE_WAVE_ZDIR            42
#define IC_SINE_WAVE_XDIR              43
#define IC_SINE_WAVE_YDIR              44
#define IC_SINE_WAVE_ZDIR              45
#define IC_SINE2_WAVE_XDIR             46
#define IC_SINE2_WAVE_YDIR             47
#define IC_SINE2_WAVE_ZDIR             48

#define IC_COMPRESSION_XDIR            50
#define IC_COMPRESSION_YDIR            51
#define IC_COMPRESSION_ZDIR            52
#define IC_EXPANSION_XDIR              53
#define IC_EXPANSION_YDIR              54
#define IC_EXPANSION_ZDIR              55
#define IC_CROSSING_JETS               56
#define IC_ELECTROSTATIC_CHANNEL       57
#define IC_DESOLVATION_CHAMBER         58

#define IC_GAS_MIX		       70
#define IC_CHEM_CORE_FLAME             71
#define IC_CHEM_INVERSE_FLAME          72
#define IC_CHEM_1DFLAME                73
#define IC_PRESSURE_GRADIENT_X         74
#define IC_PRESSURE_GRADIENT_Y         75
#define IC_PRESSURE_GRADIENT_Z         76

#define IC_VISCOUS_CHANNEL_FLOW                80
#define IC_VISCOUS_COUETTE                     80 // Duplicate of IC_VISCOUS_CHANNEL_FLOW - TAKE THIS OUT!!!
#define IC_VISCOUS_COUETTE_ISOTHERMAL_WALL     80 // Duplicate of IC_VISCOUS_CHANNEL_FLOW - TAKE THIS OUT!!!
#define IC_VISCOUS_POISEUILLE_FLOW             80 // Duplicate of IC_VISCOUS_CHANNEL_FLOW - TAKE THIS OUT!!!
#define IC_VISCOUS_PIPE_FLOW                   81
#define IC_VISCOUS_FLAT_PLATE                  82
#define IC_VISCOUS_STOKES_FLOW                 83
#define IC_VISCOUS_DRIVEN_CAVITY_FLOW          84
#define IC_VISCOUS_BACKWARD_FACING_STEP        85
#define IC_VISCOUS_BRANCHED_DUCT               86
#define IC_MIXING_LAYER                        87

#define IC_TURBULENT_PIPE_FLOW                 90
#define IC_TURBULENT_CHANNEL_FLOW              91
#define IC_TURBULENT_COFLOW                    92 
#define IC_TURBULENT_DIFFUSION_FLAME           93 
#define IC_TURBULENT_DUMP_COMBUSTOR            94 
#define IC_CHANNEL_FLOW                        95
#define IC_FREE_JET_FLAME                      96
 
#define IC_VISCOUS_COUETTE_PRESSURE_GRADIENT   97
#define IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_X 97
#define IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Y 98
#define IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Z 99

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

#define IC_ELECTRIC_FIELD_UNIFORM            200
#define IC_ELECTRIC_FIELD_QUADRUPOLE         201
#define IC_ELECTRIC_FIELD_OCTAPOLE           202
#define IC_ELECTRIC_FIELD_DOUBLE_QUADRUPOLE  203
#define IC_ELECTRIC_FIELD_DOUBLE_OCTAPOLE    204

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
#define TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD     9
#define TIME_STEPPING_IMPLICIT_ADAMS_TYPE               10
#define TIME_STEPPING_IMPLICIT_LEES                     11
#define TIME_STEPPING_IMPLICIT_TWO_STEP_TRAPEZOIDAL     12
#define TIME_STEPPING_IMPLICIT_A_CONTRACTIVE            13

#define TIME_STEPPING_MULTIGRID                         20
#define TIME_STEPPING_NKS                               21
#define TIME_STEPPING_DUAL_TIME_STEPPING                22

#define	TIME_STEPPING_LAX_FRIEDRICHS                   100
#define	TIME_STEPPING_LAX_WENDROFF                     101
#define	TIME_STEPPING_MACCORMACK                       102
#define	TIME_STEPPING_HANCOCK                          103
#define	TIME_STEPPING_SIMPLE_EXPLICIT                  104
#define	TIME_STEPPING_SIMPLE_IMPLICIT                  105
#define	TIME_STEPPING_CRANK_NICOLSON                   106
#define	TIME_STEPPING_ADE                              107

/**********************************************************************
 * CFD -- Multigrid Cycle Types                                       *
 **********************************************************************/

#define MULTIGRID_V_CYCLE               1
#define MULTIGRID_W_CYCLE               2

/**********************************************************************
 * CFD -- Reconstruction Types.                                       *
 **********************************************************************/

#define	RECONSTRUCTION_MUSCL                                   1
#define	RECONSTRUCTION_GREEN_GAUSS                             2
#define	RECONSTRUCTION_LEAST_SQUARES                           3
#define	RECONSTRUCTION_LINEAR_LEAST_SQUARES                    3 // Duplicate of RECONSTRUCTION_LEAST_SQUARES
#define RECONSTRUCTION_QUADRATIC_LEAST_SQUARES                 4
#define RECONSTRUCTION_CHARACTERISTIC                          5
#define RECONSTRUCTION_DIAMOND_PATH                            6
#define RECONSTRUCTION_LINEAR_ESSENTIALLY_NON_OSCILLATORY      7
#define RECONSTRUCTION_QUADRATIC_ESSENTIALLY_NON_OSCILLATORY   8
#define RECONSTRUCTION_CUBIC_ESSENTIALLY_NON_OSCILLATORY       9
#define RECONSTRUCTION_WEIGHTED_ESSENTIALLY_NON_OSCILLATORY   10

#define VISCOUS_RECONSTRUCTION_CARTESIAN           21
#define VISCOUS_RECONSTRUCTION_DIAMOND_PATH        22
#define VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE  23
#define VISCOUS_RECONSTRUCTION_HYBRID              24

#define DIAMONDPATH_NONE                           0
#define DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX         1
#define DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL       2
#define DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX        3
#define DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL      4
#define DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION   5

/**********************************************************************
 * CFD -- Subcell reconstruction/interpolation types.                 *
 **********************************************************************/

#define SUBCELL_RECONSTRUCTION_LINEAR              4601
#define SUBCELL_RECONSTRUCTION_QUADRATIC           4602

/**********************************************************************
 * CFD -- Limiter Types.                                              *
 **********************************************************************/

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

/**********************************************************************
 * CFD -- Flux Function Types.                                        *
 **********************************************************************/

#define FLUX_FUNCTION_GODUNOV                          1
#define	FLUX_FUNCTION_ROE                              2
#define	FLUX_FUNCTION_RUSANOV                          3
#define	FLUX_FUNCTION_HLLE                             4
#define FLUX_FUNCTION_HLLL                             5
#define FLUX_FUNCTION_LINDE                            5 // Duplicate of HLLL
#define FLUX_FUNCTION_HLLC                             6
#define FLUX_FUNCTION_OSHER                            7
#define FLUX_FUNCTION_VANLEER                          8
#define FLUX_FUNCTION_AUSM                             9
#define FLUX_FUNCTION_AUSMplus                        10
#define	FLUX_FUNCTION_ROE_PRECON_WS                   11
#define	FLUX_FUNCTION_HLLE_PRECON_WS                  12

#define FLUX_FUNCTION_GODUNOV_MB                      21
#define FLUX_FUNCTION_ROE_MB                          22
#define FLUX_FUNCTION_HLLE_MB                         23
#define FLUX_FUNCTION_VANLEER_MB                      24

#define FLUX_FUNCTION_GODUNOV_WRS                     99

#define PARTICLE_PHASE_FLUX_FUNCTION_SAUREL           31
#define PARTICLE_PHASE_FLUX_FUNCTION_MULTIVELOCITY    32
#define PARTICLE_PHASE_FLUX_FUNCTION_SAUREL_MB        33
#define PARTICLE_PHASE_FLUX_FUNCTION_MULTIVELOCITY_MB 34

/**********************************************************************
 * CFD -- Local time-stepping/preconditioning.                        *
 **********************************************************************/

#define GLOBAL_TIME_STEPPING                           0
#define SCALAR_LOCAL_TIME_STEPPING                     1
#define MATRIX_LOCAL_TIME_STEPPING                     2
#define LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER     3
#define SEMI_IMPLICIT_LOCAL_TIME_STEPPING              4
#define SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER   5

#define MELSON_SCALAR_LOCAL_TIME_STEPPING             11

/**********************************************************************
 * CFD -- Eikonal equation constants.                                 *
 **********************************************************************/

#define LEVELSET_INITIAL_EXTENSION_GEOMETRIC         660
#define LEVELSET_INITIAL_EXTENSION_EXACT             661

#define EIKONAL_SCHEME_SUSSMAN                       670
#define EIKONAL_SCHEME_RUSSO_SMEREKA                 671
#define EIKONAL_SELECTION_GODUNOV                    680
#define EIKONAL_SELECTION_ORIGINAL                   681
#define EIKONAL_SIGN_FUNCTION_DISCRETE               691
#define EIKONAL_SIGN_FUNCTION_SMEARED                692
#define EIKONAL_SIGN_FUNCTION_DERIVATIVE             693

/**********************************************************************
 * CFD -- Directions.                                                 *
 **********************************************************************/

#define CENTER                          0
#define NORTH                           1
#define SOUTH                           2
#define EAST                            3
#define WEST                            4
#define SNORTH                          5
#define SSOUTH                          6
#define SEAST                           7
#define SWEST                           8
#define NORTH_EAST                     11
#define NORTH_WEST                     12
#define SOUTH_EAST                     13
#define SOUTH_WEST                     14

/********************************************************
 * CFD -- Turbulence energy spectrum.                   *
 ********************************************************/
 
#define LEE_REYNOLDS                    0
#define LAVAL_NAZARENKO                 1
#define VON_KARMAN_PAO                  2
#define HAWORTH_POINSOT                 3
#define CHASNOV                         4
#define BELL_DAY                        5

/**********************************************************************
 * CFD -- Inline functions.                                           *
 **********************************************************************/

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

// Inline functions for AUSMplusUP flux calculation.
// M+1
inline double Mplus_1(double M)
  { return 0.5*(M + fabs(M)); }

// M-1
inline double Mminus_1(double M)
  { return 0.5*(M - fabs(M)); }

// M+2
inline double Mplus_2(double M)
  { return 0.25*sqr(M + 1.0); }

// M-2
inline double Mminus_2(double M)
  { return -0.25*sqr(M - 1.0); }

/********************************************************
 * CFD -- Define CFD structures and classes.            *
 ********************************************************/

/*************************************************************
 * Class: BCinfoCGNS - Used by the Grid Classes to store     *
 *                 BC info                                   *
 * Member functions                                          *
 *      Type       -- Returns the type of boundary condition *
 *                    As per the CGNS list in CFD.h          *
 *      Data[]     -- Lists the BC data in the CGNS          *
 *                    Order in CFD.h                         *
 *************************************************************/
class BCinfoCGNS{
 private:
 public:
  int Type;
  int DataType[5];//**update
  double DataValue[5];//**update

  BCinfoCGNS(void){
    for (int i=0;i<5;i++){
      DataType[i]=BC_DATA_UNDEFINED; DataValue[i]=0.0;}
    Type = BC_NULL;
  }
};
inline ostream &operator << (ostream &out_file, 
                             const BCinfoCGNS &G){
  out_file<<G.Type;
  for (int i=0;i<BCS_SIZE;i++)
    out_file<<" "<<G.DataType[i]<<" "<<G.DataValue[i];
  out_file<<"\n";
  return (out_file);
}
inline istream &operator >> (istream &in_file, 
                             BCinfoCGNS &G){
  double d1;
  int i1;
  in_file>> G.Type;
  for (int i=0;i<BCS_SIZE;i++){
    in_file>>G.DataType[i]>>G.DataValue[i];
  }
  return (in_file);
}

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

extern void Output_Progress_L2norm(const int Number_of_Time_Steps,
				   const double &Time,
				   const CPUTime &CPU_Time,
				   const double &Residual_L2_Norm,
				   const int First_Step,
				   const int Frequency,
				   const int progress_character);

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

extern void A_Stable_Implicit_Method_Coefficients(double &theta,
						  double &xi,
						  double &phi,
						  const int Time_Integration_Scheme);

/**********************************************************************
* Routine: Trilinear_Interpolation                                   *
*                                                                    *
* Return the templated solution state at the specified node using    *
* the linear interpolation of Holmes and Connell (AIAA-1989-1932).   *
*                                                                    *
**********************************************************************/
template <class T>
T Trilinear_Interpolation(const Vector3D &X1, const T &U1,
             const Vector3D &X2, const T &U2,
             const Vector3D &X3, const T &U3,
             const Vector3D &X4, const T &U4,
             const Vector3D &X5, const T &U5,
             const Vector3D &X6, const T &U6,
             const Vector3D &X7, const T &U7,
             const Vector3D &X8, const T &U8,
             const Vector3D &X0) {

 double w1, w2, w3, w4, w5, w6, w7, w8, den;
 Vector3D lambda, R;
 Tensor3D I;

 // Determine weighting coefficients:
 R = X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 - EIGHT*X0;
 I.xx = sqr(X1.x - X0.x) + sqr(X2.x - X0.x) +
        sqr(X3.x - X0.x) + sqr(X4.x - X0.x) +
        sqr(X5.x - X0.x) + sqr(X6.x - X0.x) +
        sqr(X7.x - X0.x) + sqr(X8.x - X0.x);
 I.xy = (X1.x - X0.x)*(X1.y - X0.y) + (X2.x - X0.x)*(X2.y - X0.y) +
        (X3.x - X0.x)*(X3.y - X0.y) + (X4.x - X0.x)*(X4.y - X0.y) +
        (X5.x - X0.x)*(X5.y - X0.y) + (X6.x - X0.x)*(X6.y - X0.y) +
        (X7.x - X0.x)*(X7.y - X0.y) + (X8.x - X0.x)*(X8.y - X0.y);
 I.yy = sqr(X1.y - X0.y) + sqr(X2.y - X0.y) +
        sqr(X3.y - X0.y) + sqr(X4.y - X0.y) +
        sqr(X5.y - X0.y) + sqr(X6.y - X0.y) +
        sqr(X7.y - X0.y) + sqr(X8.y - X0.y);
 I.yz = (X1.y - X0.y)*(X1.z - X0.z) + (X2.y - X0.y)*(X2.z - X0.z) +
        (X3.y - X0.y)*(X3.z - X0.z) + (X4.y - X0.y)*(X4.z - X0.z) +
        (X5.y - X0.y)*(X5.z - X0.z) + (X6.y - X0.y)*(X6.z - X0.z) +
        (X7.y - X0.y)*(X7.z - X0.z) + (X8.y - X0.y)*(X8.z - X0.z);
 I.zz = sqr(X1.z - X0.z) + sqr(X2.z - X0.z) +
        sqr(X3.z - X0.z) + sqr(X4.z - X0.z) +
        sqr(X5.z - X0.z) + sqr(X6.z - X0.z) +
        sqr(X7.z - X0.z) + sqr(X8.z - X0.z);
 I.xz = (X1.x - X0.x)*(X1.z - X0.z) + (X2.x - X0.x)*(X2.z - X0.z) +
        (X3.x - X0.x)*(X3.z - X0.z) + (X4.x - X0.x)*(X4.z - X0.z) +
        (X5.x - X0.x)*(X5.z - X0.z) + (X6.x - X0.x)*(X6.z - X0.z) +
        (X7.x - X0.x)*(X7.z - X0.z) + (X8.x - X0.x)*(X8.z - X0.z);
 lambda.x = (I.xy*R.y + I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy) +
            (I.xz*R.z + I.zz*R.x)/(I.xx*I.zz - I.xz*I.xz);
 lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy) +
            (I.yz*R.z - I.zz*R.y)/(I.zz*I.yy - I.yz*I.yz);
 lambda.z = (I.xz*R.x + I.zz*R.z)/(I.xx*I.yy - I.xy*I.xy) +
            (I.yz*R.y - I.yy*R.z)/(I.zz*I.yy - I.yz*I.yz);
 den = (I.xx*I.yy*I.zz + 2.0*I.xy*I.xz*I.yz - I.xx*sqr(I.yz) - sqr(I.xz)*I.yy - sqr(I.xy)*I.zz);
 lambda.x = (R.x*(I.yy*I.zz - I.yz*I.yz) + R.y*(I.xz*I.yz - I.xy*I.zz) + R.z*(I.xy*I.yz - I.xz*I.yy))/den;
 lambda.y = (R.x*(I.xz*I.yz - I.xy*I.zz) + R.y*(I.xx*I.zz - I.xz*I.xz) + R.z*(I.xy*I.xz - I.xx*I.yz))/den;
 lambda.z = (R.x*(I.xy*I.yz - I.xz*I.yy) + R.y*(I.xy*I.xz - I.xx*I.yz) + R.z*(I.xx*I.yy - I.xy*I.xy))/den;

 // Determine the weights:
 w1 = ONE + lambda*(X1 - X0);
 w2 = ONE + lambda*(X2 - X0);
 w3 = ONE + lambda*(X3 - X0);
 w4 = ONE + lambda*(X4 - X0);
 w5 = ONE + lambda*(X5 - X0);
 w6 = ONE + lambda*(X6 - X0);
 w7 = ONE + lambda*(X7 - X0);
 w8 = ONE + lambda*(X8 - X0);

 // Return the interpolated state:
 return (w1*U1 + w2*U2 + w3*U3 + w4*U4 + w5*U5 + w6*U6 + w7*U7 + w8*U8)/
        (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8);

} 

/**********************************************************************
 * Routine: Bilinear_Interpolation                                    *
 *                                                                    *
 * This function returns the templated solution variable, U0, defined *
 * at apoint, P0, that is known to lie inside an arbitary             *
 * quadrilateral cell defined by four unique vectors, P1, P2, P3, and *
 * P4, and the associated templated solution variables (U1, U2, U3,   *
 * and U4) at these points using the bilinear interpolation routine   *
 * outlined by Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13   *
 * No. 3 1992).  This function returns an error code if the point P0  *
 * falls outside of the quadrilateral.                                *
 *                                                                    *
 **********************************************************************/
template <class T>
int Bilinear_Interpolation(const T U1, const Vector2D P1,
			   const T U2, const Vector2D P2,
			   const T U3, const Vector2D P3,
			   const T U4, const Vector2D P4,
			   const Vector2D P0, T &U0) {
  double a0, a1, a2, a3;
  double b0, b1, b2, b3;
  double zeta0, eta0;
  double a, b, c;

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
    if (b*b-4*a*c < 0) return 1;
    zeta0 = (-b + sqrt(b*b - 4*a*c))/(2*a);
    // Check if zeta0 is within range, if not calculate conjugate
    if (zeta0 < 0 || zeta0 > 1) {
      zeta0 = (-b - sqrt(b*b - 4*a*c))/(2*a);
      // Check if zeta0 is within range
      if (zeta0 < 0 || zeta0 > 1) return 1;
    }
  } else if (b > TOLER || b < -TOLER) { // if (a == 0 AND b != 0)
    zeta0 = -c/b;
    // Check if zeta0 is within range
    if (zeta0 < 0 || zeta0 > 1) return 1;
  } else { // if (a == 0 AND b == 0)
    return(1);
  }

  // check for divide by zero error
  if (b2 + b3*zeta0 < TOLER && b2 + b3*zeta0 > -TOLER) return(1);

  eta0 = (P0.y - b0 - b1*zeta0)/(b2 + b3*zeta0);

  U0 = U1 + (U4 - U1)*zeta0 + (U2 - U1)*eta0 + (U1 + U3 - U2 - U4)*zeta0*eta0;
  return 0;
}

template <class T>
int Bilinear_Interpolation_ZY(const T U1, const Vector2D P1,
			      const T U2, const Vector2D P2,
			      const T U3, const Vector2D P3,
			      const T U4, const Vector2D P4,
			      const Vector2D P0, T &U0) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
         eta1, zeta1, eta2, zeta2, eta, zeta;
  T A, B, C, D;
  x  = P0.x;
  y  = P0.y;
  ax = P1.x;
  bx = P2.x - P1.x; 
  cx = P4.x - P1.x; 
  dx = P3.x + P1.x - P2.x - P4.x;
  ay = P1.y;
  by = P2.y - P1.y; 
  cy = P4.y - P1.y; 
  dy = P3.y + P1.y - P2.y - P4.y;
  aa = bx*dy - dx*by;
  bb = dy*(ax-x) + bx*cy - cx*by + dx*(y-ay);
  cc = cy*(ax-x) + cx*(y-ay);
  if (fabs(aa) < TOLER*TOLER) {
    if (fabs(bb) >= TOLER*TOLER) {
      zeta1 = -cc/bb;
    } else { 
      zeta1 = -cc/sgn(bb)*(TOLER*TOLER); 
    } 
    if (fabs(cy+dy*zeta1) >= TOLER*TOLER) {
      eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1); 
    } else { 
      eta1 = HALF;
    } 
    zeta2 = zeta1;
    eta2  = eta1;
  } else {
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { 
      zeta1 = HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; 
    } else { zeta1 = -HALF*bb/aa;
    } 
    if (fabs(cy+dy*zeta1) < TOLER*TOLER) {
      eta1 = -ONE;
    } else {
      eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1);
    }
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) {
      zeta2 = HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa; 
    } else {
      zeta2 = -HALF*bb/aa;
    }
    if (fabs(cy+dy*zeta2) < TOLER*TOLER) { 
      eta2 = -ONE;
    } else {
      eta2 = (y-ay-by*zeta2)/(cy+dy*zeta2);
    }
  }
  if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
      eta1  > -TOLER && eta1  < ONE + TOLER) {
    zeta = zeta1;
    eta  = eta1;
  } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
	     eta2  > -TOLER && eta2  < ONE + TOLER) {
    zeta = zeta2;
    eta  = eta2;
  } else {
    zeta = HALF;
    eta  = HALF;
  }
  A = U1;
  B = U2 - U1; 
  C = U4 - U1;
  D = U3 + U1 - U2 - U4;
  U0 = A + B*zeta + C*eta + D*zeta*eta;
  return 0;
}

/**********************************************************************
 * Routine: Bilinear_Interpolation_HC                                 *
 *                                                                    *
 * Return the templated solution state at the specified node using    *
 * the bilinear interpolation of Holmes and Connell (AIAA-1989-1932). *
 *                                                                    *
 **********************************************************************/
template <class T>
int Bilinear_Interpolation_HC(const T U1, const Vector2D X1,
			      const T U2, const Vector2D X2,
			      const T U3, const Vector2D X3,
			      const T U4, const Vector2D X4,
			      const Vector2D X0, T &U0) {
  Polygon P;
  P.convert(X1,X4,X3,X2);
  if (!P.point_in_polygon(X0)) return 1;
  double w1, w2, w3, w4;
  Vector2D lambda, R;
  Tensor2D I;
  // Determine weighting coefficients:
  R = X1 + X2 + X3 + X4 - FOUR*X0;
  I.xx = sqr(X1.x - X0.x) + sqr(X2.x - X0.x) +
         sqr(X3.x - X0.x) + sqr(X4.x - X0.x);
  I.xy = (X1.x - X0.x)*(X1.y - X0.y) + (X2.x - X0.x)*(X2.y - X0.y) +
         (X3.x - X0.x)*(X3.y - X0.y) + (X4.x - X0.x)*(X4.y - X0.y);
  I.yy = sqr(X1.y - X0.y) + sqr(X2.y - X0.y) +
         sqr(X3.y - X0.y) + sqr(X4.y - X0.y);
  lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
  lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
  // Determine the weights:
  w1 = ONE + lambda*(X1 - X0);
  w2 = ONE + lambda*(X2 - X0);
  w3 = ONE + lambda*(X3 - X0);
  w4 = ONE + lambda*(X4 - X0);
  // Determine the interpolated state:
  U0 = (w1*U1 + w2*U2 + w3*U3+ w4*U4)/(w1 + w2 + w3 + w4);
  // Bilinear interpolation was successful.
  return 0;
}

/**********************************************************************
 * Routine: Linear_Interpolation                                      *
 *                                                                    *
 * Given 2 points X1 and X2, their associated solution states U1 and  *
 * U2, and a point X that is known to lie on the line defined by X1   *
 * and X2, this function calculates the solution state at X using     *
 * linear interpolation.                                              *
 *                                                                    *
 **********************************************************************/
template <class T>
T Linear_Interpolation(const Vector2D &X1, const T &U1,
		       const Vector2D &X2, const T &U2,
		       const Vector2D &X) {

  // Determine the linear interpolation weight.
  double s = abs(X-X1)/abs(X2-X1);

  // Compute the interpolation.
  return U1*(ONE - s) + U2*s;

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

  double A;
  Vector2D n21, n32, n43, n14;
  T U;
  // Determine the (not unit) normals of the faces and the area
  // of the region of Green-Gauss integration.
  n21 = Vector2D((X2.y-X1.y),-(X2.x-X1.x));
  n32 = Vector2D((X3.y-X2.y),-(X3.x-X2.x));
  n43 = Vector2D((X4.y-X3.y),-(X4.x-X3.x));
  n14 = Vector2D((X1.y-X4.y),-(X1.x-X4.x));
  A = HALF*(((X2-X1)^(X4-X1)) + ((X3-X4)^(X3-X2))); assert(A > ZERO);
  // Compute Green-Gauss integration to determine the gradient of the 
  // solution state at the centre of the quadrilateral.
  U = HALF*(U1+U2);
  dUdx = U*n21.x;
  dUdy = U*n21.y;
  U = HALF*(U2+U3);
  dUdx += U*n32.x;
  dUdy += U*n32.y;
  U = HALF*(U3+U4);
  dUdx += U*n43.x;
  dUdy += U*n43.y;
  U = HALF*(U4+U1);
  dUdx += U*n14.x;
  dUdy += U*n14.y;
  dUdx = dUdx/A;
  dUdy = dUdy/A;
  // Green-Gauss integration successfully applied.
  return 0;

}

#endif /* _CFD_INCLUDED  */
