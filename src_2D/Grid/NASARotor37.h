/*======================================================================*
 *                                                                      *
 *   CLASS NASARotor37 -                                                *
 *         Contains routines to generate geometry, CFD meshes           *
 *         and return EXPERIMENTALLY measured flow field                *
 *         for NASA Rotor 37                                            *
 *   ----------------------------------------------------------------   *
 *                         written by Tomas Dusatko                     *
 *                           email: dusatko@ecf.utoronto.ca             *
 *         University of Toronto Institute for Aerospace Studies        *
 *                            August 2001                               *
 *                                                                      *
 *     Body Geometry Data and flow data courtesy of A.J. Strazisar,     *
 *                       NASA Glenn Research Center                     * 
 *                                                                      *
 *                                                                      * 
 *======================================================================*/

#ifndef _NASA_ROTOR37_INCLUDED
#define _NASA_ROTOR37_INCLUDED

// Include C++ libraries
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include various CFFC header files
#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED
 
#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED
 
#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif //_VECTOR3D_INCLUDED

#ifndef _SPLINE2D_INCLUDED
#include "../Math/Spline2D.h"
#endif // _SPLINE2D_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

// Include the 2D quadrilateral multiblock grid header file.

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

// Define constants
#define PEAK_FLOW    1
#define STALL_FLOW   2

#define R37         37

#define R37_HUB_SHROUD_GEOMETRY_FILE           "rotor37hub-shroud.dat"
#define R37_BLADE_GEOMETRY_FILE                "rotor37geometry.dat"
#define UPSTREAM_DOWNSTREAM_FLOW_DATA_FILE     "rotorUp+DownstreamFlow.dat"

#define NUMBER_OF_POINTS_SHROUD_37      32
#define NUMBER_OF_POINTS_HUB_37         31

#define INITIAL_FLOW_DATA_POINTS_R37    20
#define DOWNSTREAM_FLOW_DATA_POINTS_R37 20

#define NUMBER_OF_BLADE_SURFACES_R37    13
#define POINTS_PER_SURFACE_R37          199
#define R37_LEADING_EDGE_PT             1
#define R37_TRAILING_EDGE_PT            100
#define R37GEOM_NUMPTS                  100

/*!
 * Class: NASARotor37
 *
 * @brief Contains routines to generate geometry, CFD meshes and return
 * the experimentally measured flow field for NASA Rotor 37
 *
 * \verbatim
 *  CONSTRUCTOR AND INITIALIZATION:
 *
 *    NASARotor37
 *    NASARotor37(int, char[]) - constructor for class rotor
 *                  USAGE: NASARotor37( <'PEAK_FLOW' or 'STALL_FLOW'>, $PATH TO DATA FILES$ )
 *    init(int, char[]) - constructor for class rotor
 *                  USAGE: init( <'PEAK_FLOW' or 'STALL_FLOW'>, $PATH TO DATA FILES$ )
 *    setgas (char[]) - set the gas type (DEFAULT: "AIR")
 * 
 *  MEMBER FUNCTIONS TO RETURN UPSTREAM FLOW BOUNDARY CONDITIONS:
 * 
 *    calcMassFlux(double pspan) -  returns the calcultated mass flux at any %span based on 
 *                                  a distribution derived from the mass flux of ROTOR67  
 *    outputTP_UpstreamFlowConditions (char[]) - Output all upstream flow conditions
 *                                               (at each %span) to a file in TECPLOT FORMAT  
 *    getT_up (double) - return the temperature at any percent span
 *    getP_up (double) - return the pressure at any percent span
 *    getA_up (double) - return the speed of sound at any percent span
 *    getRho_up (double) - return the density at any percent span
 *    getMach_up (double) - return the ABSOLUTE mach number at any percent span
 *    getV_z_up (double) - return the velocity component in the axial direction at any percent span
 *    getV_theta_up (double) - return the absolute velocity component in the theta direction at any percent span
 *    getVREL_theta_up (double) - return the relatvie velocity component in the theta direction at any percent span
 *    getMachREL_up (double) - return the RELATIVE mach number at any percent span
 *    getPstateABS_up (double) - return the ABSOLUTE pSTATE at any percent span
 *    getPstateREL_up (double) - return the RELATIVE pSTATE at any percent span
 *    getMassFlow_up (int) - return the integrated mass flow at the upstream station
 *
 *
 *  MEMBER FUNCTIONS TO RETURN DOWNSTREAM FLOW BOUNDARY CONDITIONS:
 * 
 *    outputTP_DownstreamFlowConditions (char[]) - Output all downstream flow conditions 
 *                                                 (at each %span) to a file in TECPLOT FORMAT  
 *    getT_down (double) - return the temperature at any percent span
 *    getP_down (double) - return the pressure at any percent span
 *    getA_down (double) - return the speed of sound at any percent span
 *    getRho_down (double) - return the density at any percent span
 *    getMach_down (double) - return the ABSOLUTE mach number at any percent span
 *    getTheta_down (double) - return the downstream flow angle at any percent span
 *    getV_z_down (double) - return the velocity component in the axial direction at any percent span
 *    getV_theta_down (double) - return the absolute velocity component in the theta direction at any percent span
 *    getVREL_theta_down (double) - return the relative velocity component in the theta direction at any percent span
 *    getMachREL_down (double) - return the RELATIVE mach number at any percent span
 *    getPstateABS_down (double) - return the ABSOLUTE pSTATE at any percent span
 *    getPstateREL_down (double) - return the RELATIVE pSTATE at any percent span
 *    getMassFlow_down (int) - return the integrated mass flow at the downstream station 
 *  
 *
 *  MEMBER FUNCTIONS TO RETURN GEOMETRY:
 *
 *    getBladeCS(double) - return a SPLINE2D of the cross section of the blade at any %span
 *    outputTP_Geometry(char[]) - create a 3D-data file in TECPLOT format for ROTOR 37 geometry
 *    hubR (double) - return the hub radius at any axial position
 *    shroudR (double) - return the shroud radius at any axial position 
 *    findPSpan(double r, double z) - Returns the %span at any axial position z given r
 *    findRadius(double ps, double z) - Return the radius at any %span and axial position z
 *    getCamberBlade(double, int) - return a spline of the camber line for the blade at any %span
 *    getCamberTrail(double, int, double) - return a spline of the camber line extended from trailing edge
 *    getCamberLead(double, int, double) - return a spline of the camber line extended from leading edge
 *    getCamber(double, int, double, double) - Return a spline of the *whole* camber line extending from 
 *                                             trailing edge to past leading edge
 *
 *
 *  MEMBER FUNCTIONS TO RETURN FLOW DATA THROUGH BLADE:
 *		
 *              NOT AVAILABLE  
 *
 *
 *  MEMBER FUNCTIONS THAT GENERATE SINGLE AND MULTI-BLOCK MESHES:
 *
 *    - generate a mesh for the cross-section of the blade at any %span - O GRID
 *    genMeshO(double pspan, double zMin, double zMax, int iMax, int jMax,
 *	       int StretchFcnI, double BetaI, double TauI,
 *	       int StretchFcnJ, double BetaJ, double TauJ,
 *	       int NbndType, int EbndType, int SbndType, int WbndType,
 *	       int OrthN, int OrthE, int OrthS, int OrthW, int n);
 *
 *    - generate a mesh for the cross-section of the blade at any %span - C GRID
 *    genMeshC(double pspan, double zMin, double zMax, int iMax, int jMax,
 *	       int StretchFcnI, double BetaI, double TauI,
 *	       int StretchFcnJ, double BetaJ, double TauJ,
 *	       int NbndType, int EbndType, int SbndType, int WbndType,
 *	       int OrthN, int OrthE, int OrthS, int OrthW, int n);
 *
 *    - generate a mesh for the cross-section of the blade at any %span - H GRID (multi-block)
 *    NOTE: Definition parameters for each block are automatically determined    
 *    genMeshH_3x2_AUTO(Grid2D_Quad_Block **mesh,
 *                      double pspan, double zMin, double zMax,  
 *		        int iMax, int jMax, int smooth)
 *
 *    NOTE: Must specify the definition parameters for each block 
 *    genMeshH_3x2(Grid2D_Quad_Block **mesh,
 *                 double pspan, double zMin, double zMax, int n, 
 *
 *	 	   //Parameters for block (0,0) - (BOTTOM LEFT CORNER)
 *	   	   int B00_iMax, int B00_jMax,
 *		   int B00_StretchFcnI, double B00_BetaI, double B00_TauI,
 *	  	   int B00_StretchFcnJ, double B00_BetaJ, double B00_TauJ,
 *		   int B00_NbndType, int B00_EbndType, int B00_SbndType, int B00_WbndType,
 *		   int B00_OrthN, int B00_OrthE, int B00_OrthS, int B00_OrthW,
 *		
 *		   //Parameters for block (1,0)
 *		   int B10_iMax, int B10_jMax,
 *		   int B10_StretchFcnI, double B10_BetaI, double B10_TauI,
 *		   int B10_StretchFcnJ, double B10_BetaJ, double B10_TauJ,
 *		   int B10_NbndType, int B10_EbndType, int B10_SbndType, int B10_WbndType,
 *	 	   int B10_OrthN, int B10_OrthE, int B10_OrthS, int B10_OrthW,
 *		
 *		   //Parameters for block (2,0)
 *		   int B20_iMax, int B20_jMax,
 *		   int B20_StretchFcnI, double B20_BetaI, double B20_TauI,
 *      	   int B20_StretchFcnJ, double B20_BetaJ, double B20_TauJ,
 *		   int B20_NbndType, int B20_EbndType, int B20_SbndType, int B20_WbndType,
 *		   int B20_OrthN, int B20_OrthE, int B20_OrthS, int B20_OrthW,
 *
 *		   //Parameters for block (0,1)
 *		   int B01_iMax, int B01_jMax,
 *		   int B01_StretchFcnI, double B01_BetaI, double B01_TauI,
 *		   int B01_StretchFcnJ, double B01_BetaJ, double B01_TauJ,
 *	           int B01_NbndType, int B01_EbndType, int B01_SbndType, int B01_WbndType,
 *		   int B01_OrthN, int B01_OrthE, int B01_OrthS, int B01_OrthW,
 *
 *		   //Parameters for block (1,1)
 *	           int B11_iMax, int B11_jMax,
 *		   int B11_StretchFcnI, double B11_BetaI, double B11_TauI,
 *		   int B11_StretchFcnJ, double B11_BetaJ, double B11_TauJ,
 *	           int B11_NbndType, int B11_EbndType, int B11_SbndType, int B11_WbndType,
 *		   int B11_OrthN, int B11_OrthE, int B11_OrthS, int B11_OrthW,
 *		
 *		   //Parameters for block (2,1)
 *		   int B21_iMax, int B21_jMax,
 *		   int B21_StretchFcnI, double B21_BetaI, double B21_TauI,
 *		   int B21_StretchFcnJ, double B21_BetaJ, double B21_TauJ,
 *		   int B21_NbndType, int B21_EbndType, int B21_SbndType, int B21_WbndType,
 *		   int B21_OrthN, int B21_OrthE, int B21_OrthS, int B21_OrthW)
 *
 *    outputTP_Mesh3D(int, int, int, double, double, char[], char[])
 *         - generate a 3D structured HEX mesh using the genMeshH routine and output in TECPLOT FORMAT
 *
 *    output_Tetin(char fname[], int layers)
 *         - creates tetin file for NASA ROTOR 37 geometry
 * \endverbatim
 */

class NASARotor37 {

private:
public:

  //@{ @name Structures
  
  //! Used to store blade geometry data
  struct bladePt {
    double z;  // axial pos
    double r;  // rad
    double t1; // theta to upper surface
    double t2; // theta to lower surface
  }; 
  
  //@}

  //@{ @name Variables
  
  char   datafile_path[256]; //!< Path of directory that hold all the data files
  int    rotor_type;         //!< Rotor type
  double rotor_mflow;        //!< Rotor mass flow rate
  int    num_blades;         //!< Number of blades on rotor
  double rpm;                //!< rotations per minute of rotor

  double R;           //!< Gas constant
  double gmma;        //!< Gamma - (ratio of Cp to Cv)
  double pRef, tRef;  //!< Reference pressure and temperature

  Spline2D hub, shroud; //!< Splines that represent the HUB and SHROUD

  /* The following variables store data about upstream flow conditions */
  
  Spline2D p_u;           //!< pressure
  Spline2D ps_u;          //!< measured static pressure
  Spline2D t_u;           //!< temperature
  Spline2D mach_u;        //!< Mach number
  double rHub_u, rShrd_u; //!< HUB and SHROUD radius at UPSTREAM MEASURMENT STATION
  double z_up;            //!< axial location of UPSTREAM MEASURMENT STATION

  /* The following variables store data about downstream flow conditions */
  
  Spline2D p_d;           //!< pressure
  Spline2D ps_d;          //!< measured static pressure
  Spline2D t_d;           //!< temperature
  Spline2D mach_d;        //!< Mach number
  Spline2D theta_d;       //!< flow angle
  double rHub_d, rShrd_d; //!< HUB and SHROUD radius at DOWNSTREAM MEASURMENT STATION
  double z_d;             //!< axial location of DOWNSTREAM MEASURMENT STATION

  /* The following variables store the blade geometry */

  bladePt bladeData[NUMBER_OF_BLADE_SURFACES_R37][R37GEOM_NUMPTS];                
  bladePt **bladeTmp;
  
  //@{ @name Public member functions

  //! Constructor for class NASARotor37
  NASARotor37(void) { 
  }

  //! Constructor for class NASARotor37
  NASARotor37(int flowtype, char *path) { 
    strcpy(datafile_path, path);

    // Set up constants
    setgas("AIR");
    
    rotor_type=R37;
    initializeHubShroud();

    initializeBoundaryFlowConditions(flowtype);  
    initializeBladeGeom();
  }

  //! Initialize rotor data.
  inline void init(int flowtype, char *path) { 
    strcpy(datafile_path, path);

    // Set up constants
    setgas("AIR");

    rotor_type=R37;
    initializeHubShroud();

    initializeBoundaryFlowConditions(flowtype);  
    initializeBladeGeom();
  }

  //! Set the type of gas used - DEFAULT used is 'AIR'
  void setgas(char *string_ptr); 

  //! Return the hub radius at any axial position Z
  inline double hubR (double z){
    if(z<hub.Xp[0].x)
      return(hub.Xp[0].y);
    else if(z>hub.Xp[hub.np-1].x)
      return(hub.Xp[hub.np-1].y);
    else
      return( getY(z, hub)[0].y );
  }  
  
  //! Return the shroud radius at any axial position Z
  inline double shroudR (double z){
    if(z<shroud.Xp[0].x)
      return(shroud.Xp[0].y);
    else if(z>shroud.Xp[shroud.np-1].x)
      return(shroud.Xp[shroud.np-1].y);
    else
      return( getY(z, shroud)[0].y );
  }  
  
  // Return the radius at any %span at any axial position z
  double findRadius(double ps, double z);

  //! Return the %span at any axial position z given r
  double findPSpan(double r, double z);  
  
  // **************** UPSTREAM STATION ******************

  //! Output all upstream flow conditions (at each %span) to a file in TECPLOT FORMAT 
  void outputTP_UpstreamFlowConditions(char fname[256]);

  //! Return the temperature at any percent span
  inline double getT_up(double pspan){
    if ( findRadius(pspan, z_up) > t_u.Xp[t_u.np-1].x ) { //avoids roundoff error
      return(t_u.Xp[t_u.np-1].y);
    } else if( findRadius(pspan, z_up) < t_u.Xp[0].x ) {
      return(t_u.Xp[0].y);
    } else {
      return(getY(findRadius(pspan, z_up), t_u)[0].y);
    } /* endif */
  }
  
  //! Return the static pressure at any percent span
  inline double getP_up(double pspan){
    if ( findRadius(pspan, z_up) > ps_u.Xp[ps_u.np-1].x ) { //avoids roundoff error
      return(ps_u.Xp[ps_u.np-1].y);
    } else if( findRadius(pspan, z_up) < ps_u.Xp[0].x ) {
      return(ps_u.Xp[0].y);
    } else {
      return(getY(findRadius(pspan, z_up), ps_u)[0].y);
    } /* endif */
  }
  
  //! Return the speed of sound at any percent span
  inline double getA_up(double pspan){
    return( sqrt(gmma*R*getT_up(pspan)));
  }

  //! Return the density at any percent span
  inline double getRho_up(double pspan){
    return( getP_up(pspan) / (R * getT_up(pspan) ));
  }

  //! Return the ABSOLUTE mach number at any percent span
  inline double getMach_up(double pspan){
    if ( findRadius(pspan, z_up) > mach_u.Xp[mach_u.np-1].x ) { //avoids roundoff error
      return(mach_u.Xp[mach_u.np-1].y);
    } else if( findRadius(pspan, z_up) < mach_u.Xp[0].x ) {
      return(mach_u.Xp[0].y);
    } else {
      return( getY(findRadius(pspan, z_up), mach_u)[0].y );
    } /* endif */
  }
 
  //! Return the velocity component in the axial direction at any percent span
  inline double getV_z_up(double pspan){
    return( getMach_up(pspan)*getA_up(pspan) );
  }
  
  //! Return the absolute velocity component in the theta direction at any percent span
  inline double getV_theta_up(double pspan){
    return( ZERO );
  }

  //! Return the relative velocity component in the theta direction at any percent span
  inline double getVREL_theta_up(double pspan){
    return( getV_theta_up(pspan) + rpm * findRadius(pspan, z_up) * M_PI / 30.0 );
  }

  //! Return the RELATIVE mach number at any percent span
  inline double getMachREL_up(double pspan){
    return(sqrt(pow(getV_z_up(pspan),2.0)+pow(getVREL_theta_up(pspan),2.0))/
           getA_up(pspan));
  }

  //! Return the ABSOLUTE pSTATE at any percent span
  template <class pState> void getPstateABS_up(pState &W, double ps){
    W[1] = getRho_up(ps);
    W[2] = getV_z_up(ps);
    W[3] = getV_theta_up(ps);
    W[4] = getP_up(ps);
  }

  //! Return the RELATIVE pSTATE at any percent span
  template <class pState> void getPstateREL_up(pState &W, double ps){
    W[1] = getRho_up(ps);
    W[2] = getV_z_up(ps);
    W[3] = getVREL_theta_up(ps);
    W[4] = getP_up(ps);
  }

  // **************** DOWNSTREAM STATION ******************

  //! Output all downstream boundary flow conditions (at each %span) to a file in TECPLOT FORMAT 
  void outputTP_DownstreamFlowConditions(char fname[256]);

  //! Return the temperature at any percent span
  inline double getT_down(double pspan){
    if ( findRadius(pspan, z_d) > t_d.Xp[t_d.np-1].x ) { //avoids roundoff error
      return(t_d.Xp[t_d.np-1].y);
    } else if( findRadius(pspan, z_d) < t_d.Xp[0].x ) {
      return(t_d.Xp[0].y);
    } else {
      return(getY(findRadius(pspan, z_d), t_d)[0].y);
    } /* endif */
  }
  
  //! Return the static pressure at any percent span
  inline double getP_down(double pspan){
    if ( findRadius(pspan, z_d) > ps_d.Xp[ps_d.np-1].x ) { //avoids roundoff error
      return(ps_d.Xp[ps_d.np-1].y);
    } else if( findRadius(pspan, z_d) < ps_d.Xp[0].x ) {
      return(ps_d.Xp[0].y);
    } else {
      return(getY(findRadius(pspan, z_d), ps_d)[0].y );
    } /* endif */
  }
  
  //! Return the speed of sound at any percent span
  inline double getA_down(double pspan){
    return( sqrt(gmma*R*getT_down(pspan)));
  }

  //! Return the density at any percent span
  inline double getRho_down(double pspan){
    return( getP_down(pspan) / (R * getT_down(pspan) ));
  }

  //! Return the ABSOLUTE mach number at any percent span
  inline double getMach_down(double pspan){
    if (findRadius(pspan, z_d)>mach_d.Xp[mach_d.np-1].x ) { //avoids roundoff error
      return(mach_d.Xp[mach_d.np-1].y);
    } else if( findRadius(pspan, z_d)<mach_d.Xp[0].x ) {
      return(mach_d.Xp[0].y);
    } else {
      return(getY(findRadius(pspan, z_d), mach_d)[0].y );
    } /* endif */
  }
 
  //! Return the flow angle at any percent span
  inline double getTheta_down(double pspan){
    return( ZERO );
  }

  //! Return the velocity component in the axial direction at any percent span
  inline double getV_z_down(double pspan){
    return( getMach_down(pspan)*getA_down(pspan)*cos(getTheta_down(pspan)) );
  }
  
  //! Return the absolute velocity component in the theta direction at any percent span
  inline double getV_theta_down(double pspan){
    return( -getMach_down(pspan)*getA_down(pspan)*sin(getTheta_down(pspan)) );
  }

  //! Return the relative velocity component in the theta direction at any percent span
  inline double getVREL_theta_down(double pspan){
    return( getV_theta_down(pspan) + rpm * findRadius(pspan, z_d) * M_PI / 30.0 );
  }

  //! Return the RELATIVE mach number at any percent span
  inline double getMachREL_down(double pspan){
    return(sqrt(pow(getV_z_down(pspan),2.0)+pow(getVREL_theta_down(pspan),2.0))/
           getA_down(pspan));
  }

  //! Return the ABSOLUTE pSTATE at any percent span
  template <class pState> void getPstateABS_down(pState &W, double ps){
    W[1] = getRho_down(ps);
    W[2] = getV_z_down(ps);
    W[3] = getV_theta_down(ps);
    W[4] = getP_down(ps);
  }

  //! Return the RELATIVE pSTATE at any percent span
  template <class pState> void getPstateREL_down(pState &W, double ps){
    W[1] = getRho_down(ps);
    W[2] = getV_z_down(ps);
    W[3] = getVREL_theta_down(ps);
    W[4] = getP_down(ps);
  }  

  //! Returns the calcultated mass flux at any %span based on a 
  //! distribution derived from the mass flux of ROTOR67  
  double calcMassFlux(double pspan);

  //! Return the integrated mass flow using splined data UPSTREAM
  double getMassFlow_up(int N);

  //! Return the integrated mass flow using splined data DOWNSTREAM
  double getMassFlow_down(int N);

  
  //! Return a spline of the cross section of one blade at any %span
  Spline2D getBladeCS(double pspan);

  //! Return a spline of the camber line for the blade at any %span
  //! n=number of points to use
  Spline2D getCamberBlade(double pspan, int n);

  //! Return a spline of the camber line extended from trailing edge
  //! n=number of points to use
  Spline2D getCamberTrail(double pspan, int n, double max);

  //! Return a spline of the camber line extended from leading edge
  //! n=number of points to use
  Spline2D getCamberLead(double pspan, int n, double min);

  //! Return a spline of the *whole* camber line extended from trailing edge to past leading edge
  //! n=number of points to use
  Spline2D getCamber(double pspan, int n, double min, double max);

  //! Create a 3D geometry file in TECPLOT format for ROTOR 37
  void outputTP_Geometry(char filename[]);

  //! Create a 3D-data file in ICEM-CFD (TETIN) format - ready for meshing
  //! Parameter l specifies the number of contour levels to generate around blade
  void output_Tetin(char fname[], int layers);
  
  //! Generate a mesh for the cross-section of the blade at any %span - O GRID
  Grid2D_Quad_Block genMeshO(double pspan, double zMin, double zMax,
			     const int Number_of_Cells_Idir,
			     const int Number_of_Cells_Jdir,
			     const int Number_of_Ghost_Cells,
			     int StretchFcnI, double BetaI, double TauI,
			     int StretchFcnJ, double BetaJ, double TauJ,
			     int NbndType, int EbndType, int SbndType, int WbndType,
			     int OrthN, int OrthE, int OrthS, int OrthW, int n);
  
  //! Generate a mesh for the cross-section of the blade at any %span - C GRID
  Grid2D_Quad_Block genMeshC(double pspan, double zMin, double zMax,
			     const int Number_of_Cells_Idir,
			     const int Number_of_Cells_Jdir,
			     const int Number_of_Ghost_Cells,
			     int StretchFcnI, double BetaI, double TauI,
			     int StretchFcnJ, double BetaJ, double TauJ,
			     int NbndType, int EbndType, int SbndType, int WbndType,
			     int OrthN, int OrthE, int OrthS, int OrthW, int n);

  //! Generate a mesh using genMeshH_3x2 with preset parameters.
  Grid2D_Quad_Block** genMeshH_3x2_AUTO(Grid2D_Quad_Block **mesh,
					double pspan, double zMin, double zMax,  
					const int Number_of_Cells_Idir,
					const int Number_of_Cells_Jdir,
					const int Number_of_Ghost_Cells,
					int smooth);

  //! Generate a H-grid for the cross-section of the blade at any %span.
  Grid2D_Quad_Block** genMeshH_3x2(Grid2D_Quad_Block **mesh,
                
                double pspan, double zMin, double zMax, int n, 
		
		// Parameters for block (0,0) - (BOTTOM LEFT CORNER)
		const int B00_Number_of_Cells_Idir, const int B00_Number_of_Cells_Jdir, const int B00_Number_of_Ghost_Cells,
		int B00_StretchFcnI, double B00_BetaI, double B00_TauI,
		int B00_StretchFcnJ, double B00_BetaJ, double B00_TauJ,
		int B00_NbndType, int B00_EbndType, int B00_SbndType, int B00_WbndType,
		int B00_OrthN, int B00_OrthE, int B00_OrthS, int B00_OrthW,
		
		// Parameters for block (1,0)
		const int B10_Number_of_Cells_Idir, const int B10_Number_of_Cells_Jdir, const int B10_Number_of_Ghost_Cells,
		int B10_StretchFcnI, double B10_BetaI, double B10_TauI,
		int B10_StretchFcnJ, double B10_BetaJ, double B10_TauJ,
		int B10_NbndType, int B10_EbndType, int B10_SbndType, int B10_WbndType,
		int B10_OrthN, int B10_OrthE, int B10_OrthS, int B10_OrthW,

		// Parameters for block (2,0)
		const int B20_Number_of_Cells_Idir, const int B20_Number_of_Cells_Jdir, const int B20_Number_of_Ghost_Cells,
		int B20_StretchFcnI, double B20_BetaI, double B20_TauI,
		int B20_StretchFcnJ, double B20_BetaJ, double B20_TauJ,
		int B20_NbndType, int B20_EbndType, int B20_SbndType, int B20_WbndType,
		int B20_OrthN, int B20_OrthE, int B20_OrthS, int B20_OrthW,

		// Parameters for block (0,1)
		const int B01_Number_of_Cells_Idir, const int B01_Number_of_Cells_Jdir, const int B01_Number_of_Ghost_Cells,
		int B01_StretchFcnI, double B01_BetaI, double B01_TauI,
		int B01_StretchFcnJ, double B01_BetaJ, double B01_TauJ,
		int B01_NbndType, int B01_EbndType, int B01_SbndType, int B01_WbndType,
		int B01_OrthN, int B01_OrthE, int B01_OrthS, int B01_OrthW,
		
		// Parameters for block (1,1)
		const int B11_Number_of_Cells_Idir, const int B11_Number_of_Cells_Jdir, const int B11_Number_of_Ghost_Cells,
		int B11_StretchFcnI, double B11_BetaI, double B11_TauI,
		int B11_StretchFcnJ, double B11_BetaJ, double B11_TauJ,
		int B11_NbndType, int B11_EbndType, int B11_SbndType, int B11_WbndType,
		int B11_OrthN, int B11_OrthE, int B11_OrthS, int B11_OrthW,
		
		// Parameters for block (2,1)
		const int B21_Number_of_Cells_Idir, const int B21_Number_of_Cells_Jdir, const int B21_Number_of_Ghost_Cells,
		int B21_StretchFcnI, double B21_BetaI, double B21_TauI,
		int B21_StretchFcnJ, double B21_BetaJ, double B21_TauJ,
		int B21_NbndType, int B21_EbndType, int B21_SbndType, int B21_WbndType,
		int B21_OrthN, int B21_OrthE, int B21_OrthS, int B21_OrthW);
  
  //! Generate a 3D structured HEX mesh using the genMeshH routine and output in TECPLOT FORMAT
  void outputTP_Mesh3D(const int Number_of_Cells_Idir, const int Number_of_Cells_Jdir, int kMax,
		       double zMin, double zMax,
		       char fnameTop[], char fnameBot[]);

  //@}

  //@{ @name Internally used member functions

  //! Utility member function: discard lines in file 'junkfile'
  inline void skipL(int i, ifstream &junkfile){
    char readtmp[1000];
    for(int j=1; j<=i; ++j){                                
      junkfile.getline(readtmp, sizeof(readtmp), '\n');  
    }  
  }

  //! Utility member function: converts a string to an integer
  void itostr(char *TargetString, unsigned theint);

  //! Member functions for HUB/SHROUD
  void initializeHubShroud (void);
  
  //! Member functions to return initial flow conditions.  
  void initializeBoundaryFlowConditions (int mflow); 
  double f(double);
  double r(double);
  double f_d(double r);
  double r_d(double r);  
  double calcMach(double p, double t, double A, double mflow); //calculate mach number using 
                                                               //Newton-Raphson method

  //! Member functions to return blade geometry
  void initializeBladeGeom(void);
  Vector3D lInt(const Vector3D &v1, const Vector3D &v2, const double psMin,
		const double psMax, const double pspan);
  void allocateTmpBlade(int surfaces, int points);
  void deallocateTmpBlade(int surfaces, int points);

  //! Member functions to return the camber line of a given cross-section
  double len(Vector2D v1, Vector2D v2);
  double findLT(double pspan, int &imin, int &imax);
  void rotateCS(Spline2D &cs, double rad);
  Vector2D rotate(Vector2D v, const double &rad);
  Spline2D transformCS(double pspan, double &m, int &imin, int &imax);
  Spline2D getCamberAndShift (double pspan, int n, double min, double max, double angle);
 
  //! Member functions to generate meshes  
  //! generates a mesh on a given block that already has the boundary splines defined
  void genMeshBlock(Grid2D_Quad_Block &block,
		    int n,
		    const int Number_of_Cells_Idir,
		    const int Number_of_Cells_Jdir,
		    const int Number_of_Ghost_Cells,
		    int node_init_procedure,
		    int  StretchFcnI, double  BetaI, double  TauI,
		    int  StretchFcnJ, double  BetaJ, double  TauJ,
		    int  NbndType, int  EbndType, int  SbndType, int  WbndType,
		    int  OrthN, int  OrthE, int  OrthS, int  OrthW) ;

  //@}

}; /* END OF CLASS NASARotor37 */

/**********************************************************************
 * Routine: initializeHubShroud                                       *
 *                                                                    *
 * Routine to initialize the hub and shroud splines.                  *
 *                                                                    *
 **********************************************************************/
inline void NASARotor37::initializeHubShroud (void) {

  int i;  
  char junk[100];
  char fname[256];
  ifstream rotorGeom;

  //open data file
  strcpy(fname, datafile_path);
  strcat(fname, R37_HUB_SHROUD_GEOMETRY_FILE);
  rotorGeom.open(fname, ios::in);
  if(rotorGeom.fail()) {
    cout << "\n ERROR: Unable to open hub/shroud data for NASA rotor 37.\n";
    cout.flush();
  } /* endif */
  assert(!rotorGeom.fail());

  rotorGeom.setf(ios::skipws);

  //Choose spline type and allocate memory for spline
  hub.settype(SPLINE2D_QUINTIC);
  hub.allocate(NUMBER_OF_POINTS_HUB_37);
  shroud.settype(SPLINE2D_QUINTIC);
  shroud.allocate(NUMBER_OF_POINTS_SHROUD_37);

  //Set point types for first and last (default=normal)
  hub.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  hub.tp[hub.np-1]=SPLINE2D_POINT_SHARP_CORNER;

  shroud.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  shroud.tp[shroud.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  
  for ( i=1; i<hub.np-1; ++i)
    hub.tp[i]=SPLINE2D_POINT_NORMAL;

  for ( i=1; i<shroud.np-1; ++i )
    shroud.tp[i]=SPLINE2D_POINT_NORMAL;

  for ( i=1; i<=5; ++i ) 
    rotorGeom.getline(junk, sizeof(junk), '\n');

  for ( i=0; i<=hub.np-1; ++i ) 
    rotorGeom >> hub.Xp[i].x;

  for ( i=1; i<=4; ++i ) 
    rotorGeom.getline(junk, sizeof(junk), '\n');
  
  for ( i=0; i<=hub.np-1; ++i )  
    rotorGeom >> hub.Xp[i].y;

  for ( i=1; i<=4; ++i ) 
    rotorGeom.getline(junk, sizeof(junk), '\n');

  for ( i=0; i<=shroud.np-1; ++i ) 
    rotorGeom >> shroud.Xp[i].x;

  for ( i=1; i<=4; ++i ) 
    rotorGeom.getline(junk, sizeof(junk), '\n');
  
  for ( i=0; i<=shroud.np-1; ++i )  
    rotorGeom >> shroud.Xp[i].y;

  rotorGeom.close();

 // Convert cm to m
  for ( i=0; i<shroud.np; ++i ) {
    shroud.Xp[i].x = shroud.Xp[i].x / 100;
    shroud.Xp[i].y = shroud.Xp[i].y / 100;
  }  

  // Convert cm to m
  for ( i=0; i<hub.np; ++i ) {
    hub.Xp[i].x = hub.Xp[i].x / 100;
    hub.Xp[i].y = hub.Xp[i].y / 100;
  }  

  hub.pathlength();
  shroud.pathlength();
}  

/**********************************************************************
 * Routine: calcMassFlux                                              *
 *                                                                    *
 * Returns the calcultated mass flux based on a distribution derived  *
 * from the mass flux of ROTOR67.                                     *
 *                                                                    *
 **********************************************************************/
inline double NASARotor37::calcMassFlux(double pspan) {
  if(rotor_mflow == 20.5114) //peak flow
    return( 627.9860043/M_PI * pow( sin(M_PI*(pspan/100)), 0.114059369488 ) );
  else //stall flow
    return( 592.7418919/M_PI * pow( sin(M_PI*(pspan/100)), 0.114059369488 ) );
}

/**********************************************************************
 * Routine: calcMach                                                  *
 *                                                                    *
 * Calculate mach number using Newton-Raphson method.                 *
 *                                                                    *
 **********************************************************************/
inline double NASARotor37::calcMach(double p, double t, double dA, double mflow) {
  int i;
  double A, B, C, m, mNew;

  m=0;
  mNew=0.8;

  A = mflow / (dA * p * sqrt(gmma/(t*R)) );
  B = 0.5*(gmma-1);
  C = 0.5*(gmma+1)/(gmma-1);
   
  while(fabs(m-mNew)>TOLER*TOLER) {
    m=mNew;
    mNew = m - (A*pow(1+B*m*m, C) - m) / (C*A*pow(1+B*m*m, C-1)*2*B*m - 1);
  }
  return(mNew);
}

/**********************************************************************
 * Routine: initializeBoundaryFlowConditions                          *
 **********************************************************************/
inline void NASARotor37::initializeBoundaryFlowConditions(int mflow) {

  ifstream flowData;
  double r[INITIAL_FLOW_DATA_POINTS_R37];
  double dA[INITIAL_FLOW_DATA_POINTS_R37];
  int i, nPts;
  double rTmp;
  char fname[256];

  strcpy(fname, datafile_path);
  strcat(fname, UPSTREAM_DOWNSTREAM_FLOW_DATA_FILE);
  flowData.open(fname, ios::in);
  if(flowData.fail()) {
    cout << "\n ERROR: Upstream flow data file for NASA rotor 37 not found.";
    cout.flush();
  } /* endif */
  assert(!flowData.fail());
      
  if (mflow != 1 && mflow != 2) {
    cout << "ERROR: Incorrect flow type for NASA rotor 37.";
    cout.flush();
  } /* endif */
  assert(mflow == 1 || mflow == 2);

  //set mass flow rate
  if (mflow == 1)
    rotor_mflow = 20.5114;
  else
    rotor_mflow = 19.36025;

  //Hub radius at AEROSTATION 1
  rHub_u = 17.5259e-2;

  //Shroud Radius at AEROSTATION 1
  rShrd_u = 25.6692e-2; 

  //Upstream and Downstream measurment station axial positions
  z_up=-4.19E-2;

  z_d=10.64E-2;

  //Set RPM
  rpm = 17188.7; 

  //Number of blades
  num_blades=36;

  //set reference temp and pressure
  pRef = 101325;
  tRef = 288.15;
	
  //SETUP splines to store initial pressure and temperature data
  p_u.settype(SPLINE2D_LINEAR);
  p_u.allocate(INITIAL_FLOW_DATA_POINTS_R37);
  ps_u.settype(SPLINE2D_LINEAR);
  ps_u.allocate(INITIAL_FLOW_DATA_POINTS_R37);
  t_u.settype(SPLINE2D_LINEAR);
  t_u.allocate(INITIAL_FLOW_DATA_POINTS_R37);
  mach_u.settype(SPLINE2D_LINEAR);
  mach_u.allocate(INITIAL_FLOW_DATA_POINTS_R37);
          
  p_u.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  p_u.tp[p_u.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  ps_u.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  ps_u.tp[ps_u.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  t_u.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  t_u.tp[t_u.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  mach_u.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mach_u.tp[mach_u.np-1]=SPLINE2D_POINT_SHARP_CORNER;

  for ( i=1; i<INITIAL_FLOW_DATA_POINTS_R37-1; ++i ) {
    t_u.tp[i]=SPLINE2D_POINT_NORMAL;
    p_u.tp[i]=SPLINE2D_POINT_NORMAL;
    ps_u.tp[i]=SPLINE2D_POINT_NORMAL;
    mach_u.tp[i]=SPLINE2D_POINT_NORMAL;
  }

  //add in TWO extra data points coresponding to the hub and shroud
  //static pressure values supplied in NASA data
      
  //set radius
  r[0]=rHub_u;
  r[INITIAL_FLOW_DATA_POINTS_R37-1]=rShrd_u;
  p_u.Xp[0].x=rHub_u;
  p_u.Xp[p_u.np-1].x=rShrd_u;
  ps_u.Xp[0].x=p_u.Xp[0].x;
  ps_u.Xp[ps_u.np-1].x=p_u.Xp[p_u.np-1].x;
  t_u.Xp[0].x=p_u.Xp[0].x;
  t_u.Xp[t_u.np-1].x=p_u.Xp[p_u.np-1].x;
  mach_u.Xp[0].x=p_u.Xp[0].x;
  mach_u.Xp[mach_u.np-1].x=p_u.Xp[p_u.np-1].x;

  if(mflow==1) { //peak
    ps_u.Xp[0].y = 0.9864 * pRef;           //HUB
    ps_u.Xp[ps_u.np-1].y = 0.9435 * pRef;      //SHROUD
    p_u.Xp[0].y = 0.9864 * pRef;            //HUB
    p_u.Xp[p_u.np-1].y = 0.9435 * pRef;       //SHROUD
  } else {      //stall
    ps_u.Xp[0].y = 0.9864 * pRef;           //HUB
    ps_u.Xp[ps_u.np-1].y = 0.9435 * pRef;      //SHROUD
    p_u.Xp[0].y = 0.9864 * pRef;            //HUB
    p_u.Xp[p_u.np-1].y = 0.9435 * pRef;       //SHROUD
  }

  //flow velocity goes to ZERO at hub and shroud
  t_u.Xp[0].y=tRef;
  t_u.Xp[t_u.np-1].y=tRef;
  mach_u.Xp[0].y=0;
  mach_u.Xp[mach_u.np-1].y=0;

  //read in %spans and convert to radius
  skipL(18, flowData);
  for ( i=INITIAL_FLOW_DATA_POINTS_R37-2; i>=1; --i ) {
    flowData >> r[i];
    r[i]/=100.0;
  }

  skipL(4, flowData);
  //read in stagnation pressure
  for ( i=INITIAL_FLOW_DATA_POINTS_R37-2; i>=1; --i ) {
    flowData >> p_u.Xp[i].y;
    p_u.Xp[i].y*=pRef;
    p_u.Xp[i].x=r[i];
    ps_u.Xp[i].x=r[i];
  }

  skipL(4, flowData);
  //read in stagnation TEMP
  for ( i=INITIAL_FLOW_DATA_POINTS_R37-2; i>=1; --i ) {
    flowData >> t_u.Xp[i].y;
    t_u.Xp[i].y*=tRef;
    t_u.Xp[i].x=r[i];
  }

  //calculate Area of annuluses
  for ( i=INITIAL_FLOW_DATA_POINTS_R37-2; i>=1; --i ) {
    dA[i]=M_PI*(r[i+1]+r[i])/2.0 - M_PI*(r[i-1]+r[i])/2.0;
  }
  
  //calculate discrete points for mach number distribution
  for ( i=INITIAL_FLOW_DATA_POINTS_R37-2; i>=1; --i ) {
    mach_u.Xp[i].x=r[i];
    mach_u.Xp[i].y=calcMach(p_u.Xp[i].y, t_u.Xp[i].y, dA[i], calcMassFlux(findPSpan(r[i], z_up))*dA[i]);
  }

  //convert stagnation to static temperatures and pressures
  for ( i=INITIAL_FLOW_DATA_POINTS_R37-2; i>=1; --i ) {
    t_u.Xp[i].y /= 1 + 0.5*(gmma-1)*mach_u.Xp[i].y*mach_u.Xp[i].y;
    ps_u.Xp[i].y = p_u.Xp[i].y / pow(1 + 0.5*(gmma-1)*mach_u.Xp[i].y*mach_u.Xp[i].y, gmma/(gmma-1));
  }

  //set boundary values
  ps_u.Xp[0].y=ps_u.Xp[1].y;
  ps_u.Xp[ps_u.np-1].y=ps_u.Xp[ps_u.np-2].y;


  //calculate pathlengths
  t_u.pathlength();
  ps_u.pathlength();
  p_u.pathlength();
  mach_u.pathlength();

  //read in ***************** DOWNSTREAM **************** data
  
  flowData.close();
  flowData.open(fname, ios::in);
  if(flowData.fail()) {
    cout << "\n ERROR: Downstream flow data file for NASA rotor 37 not found.";
    cout.flush();
  } /* endif */
  assert(!flowData.fail());

  //number of data pts
  nPts = DOWNSTREAM_FLOW_DATA_POINTS_R37;

  // Hub Radius at AEROSTATION 2
  rHub_d = 19.381E-2;

  // Shroud Radius at AEROSTATION 2
  rShrd_d = 23.823E-2;

  //Locate appropriate place in file
  if (mflow == 1) //peak
    skipL(163, flowData);
  else {          //stall
    skipL(131, flowData); 
  }	    

  //SETUP splines to store initial pressure and temperature data
  p_d.settype(SPLINE2D_LINEAR);
  p_d.allocate(DOWNSTREAM_FLOW_DATA_POINTS_R37);
  ps_d.settype(SPLINE2D_LINEAR);
  ps_d.allocate(DOWNSTREAM_FLOW_DATA_POINTS_R37);
  t_d.settype(SPLINE2D_LINEAR);
  t_d.allocate(DOWNSTREAM_FLOW_DATA_POINTS_R37);
  mach_d.settype(SPLINE2D_LINEAR);
  mach_d.allocate(7);
          
  p_d.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  p_d.tp[nPts-1]=SPLINE2D_POINT_SHARP_CORNER;
  ps_d.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  ps_d.tp[nPts-1]=SPLINE2D_POINT_SHARP_CORNER;
  t_d.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  t_d.tp[nPts-1]=SPLINE2D_POINT_SHARP_CORNER;
  mach_d.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mach_d.tp[mach_d.np-1]=SPLINE2D_POINT_SHARP_CORNER;
   
  for ( i=1; i<nPts-1; ++i ) {
    t_d.tp[i]=SPLINE2D_POINT_NORMAL;
    p_d.tp[i]=SPLINE2D_POINT_NORMAL;
    ps_d.tp[i]=SPLINE2D_POINT_NORMAL;
  }

  for ( i=1; i<mach_d.np-1; ++i )
    mach_d.tp[i]=SPLINE2D_POINT_NORMAL;

  //add in TWO extra data points coresponding to the hub and shroud
       
  //set radius
  p_d.Xp[0].x=rHub_d;
  p_d.Xp[nPts-1].x=rShrd_d;
  ps_d.Xp[0].x=p_d.Xp[0].x;
  ps_d.Xp[nPts-1].x=p_d.Xp[nPts-1].x;
  t_d.Xp[0].x=p_d.Xp[0].x;
  t_d.Xp[nPts-1].x=p_d.Xp[nPts-1].x;
  mach_d.Xp[0].x=p_d.Xp[0].x;
  mach_d.Xp[mach_d.np-1].x=p_d.Xp[nPts-1].x;

  if(mflow==1) { //peak
    ps_d.Xp[0].y =  2.11878 * pRef;           //HUB
    ps_d.Xp[nPts-1].y = 1.95925 * pRef;      //SHROUD
    p_d.Xp[0].y =  2.11878 * pRef;            //HUB
    p_d.Xp[nPts-1].y = 1.95925 * pRef;       //SHROUD
  } else {      //stall
    ps_d.Xp[0].y = 2.20102 * pRef;           //HUB
    ps_d.Xp[nPts-1].y = 2.05388 * pRef;      //SHROUD
    p_d.Xp[0].y = 2.20102 * pRef;            //HUB
    p_d.Xp[nPts-1].y = 2.05388 * pRef;       //SHROUD
  }

  if(mflow==1) { //peak
    t_d.Xp[0].y=1.25594*tRef;
    t_d.Xp[nPts-1].y=1.28355*tRef; 
  } else {      //stall
    t_d.Xp[0].y=1.27072*tRef;
    t_d.Xp[nPts-1].y=1.32758*tRef;
  }

  //flow velocity goes to ZERO at hub and shroud
  mach_d.Xp[0].y=0;
  mach_d.Xp[mach_d.np-1].y=0;

  //read in static pressure
  for ( i=nPts-2; i>=1 ; --i ) {
    flowData >> ps_d.Xp[i].y;
    ps_d.Xp[i].y *= pRef;
    p_d.Xp[i].y *= pRef;
  }

  skipL(2, flowData);

  //read in radius
  for ( i=nPts-2; i>=1 ; --i ) {
    flowData >> rTmp;
    rTmp=findRadius(rTmp, z_d);
    p_d.Xp[i].x=rTmp;
    ps_d.Xp[i].x=rTmp;
    t_d.Xp[i].x=rTmp;
  }

  //read in static temperature
  skipL(4, flowData);
  for ( i=nPts-2; i>=1 ; --i ) {
    flowData >> t_d.Xp[i].y;
    t_d.Xp[i].y *= tRef;
  }

  //read in relative mach number
  skipL(8, flowData);
  for ( i=1; i<6 ; ++i ) {
    flowData >> mach_d.Xp[i].y;
  }

  //read in radius
  for ( i=1; i<6 ; ++i ) {
    flowData >> mach_d.Xp[i].x;
    mach_d.Xp[i].x=findRadius(mach_d.Xp[i].x, z_d);
  }

// CAN'T BE DONE - DATA IS STRANGE?!
//   //calc absolute mach number
//    for ( i=1; i<6 ; ++i ) {
//      cout <<  mach_d.Xp[i].y << "\n";
//  cout <<  getA_down(findPSpan(mach_d.Xp[i].x, z_d)) << "\n";
//  cout <<  getV_theta_down(findPSpan(mach_d.Xp[i].x, z_d)) << "\n";

//      mach_d.Xp[i].y = sqrt(pow(getA_down(findPSpan(mach_d.Xp[i].x, z_d))*mach_d.Xp[i].y,2)-
//  			  pow(getV_theta_down(findPSpan(mach_d.Xp[i].x, z_d)),2))/getA_down(findPSpan(mach_d.Xp[i].x, z_d));
//      cout <<  mach_d.Xp[i].y << "\n";
//    }

  //calculate pathlengths
  p_d.pathlength();
  t_d.pathlength();
  ps_d.pathlength();
  mach_d.pathlength();  
  
}

/**********************************************************************
 * Routine: initializeBladeGeom                                       *
 **********************************************************************/
inline void NASARotor37::initializeBladeGeom(void) {
  int i, j;
  ifstream bladeDataFile;
  double skip;
  char fname[256];

  strcpy(fname,datafile_path);
  strcat(fname, R37_BLADE_GEOMETRY_FILE);
  bladeDataFile.open(fname, ios::in);
  if(bladeDataFile.fail()) {
    cout << "\n ERROR: Blade geometry data file for NASA rotor 37 not found.";
    cout.flush();
  } /* endif */
  assert(!bladeDataFile.fail());
  
  //skip
  skipL(5, bladeDataFile);

  //READ IN ALL GEOMETRY DATA - ONLY BLADE DATA
  for ( i=0; i<NUMBER_OF_BLADE_SURFACES_R37; ++i ) { 
    for ( j=0; j<R37GEOM_NUMPTS; ++j ) {
      //read Z
      bladeDataFile >> skip >> bladeData[i][j].z;
    
      //read in R
      bladeDataFile >> bladeData[i][j].r;
    
      //read in theta 1
      bladeDataFile >> bladeData[i][j].t1;
          
      //read in theta 2
      bladeDataFile >> bladeData[i][j].t2;
    }
    skipL(3, bladeDataFile);
  }
    
  //Convert to Z,R,theta1,theta2 format in proper units
  for ( i=0; i<NUMBER_OF_BLADE_SURFACES_R37; ++i)
    for ( j=0; j<R37GEOM_NUMPTS; ++j ) {
      bladeData[i][j].z/=100;
      bladeData[i][j].r/=100;
    }
 
}

/**********************************************************************
 * Routine: outputTP_UpstreamFlowConditions                           *
 *                                                                    *
 * Output all upstream flow conditions (at each %span) to a file in   *
 * TECPLOT FORMAT.                                                    *
 *                                                                    *
 **********************************************************************/
inline void NASARotor37::outputTP_UpstreamFlowConditions(char fname[256]) {

  int i;
  ofstream outf;

  outf.open(fname, ios::out);
  assert(!outf.fail());

  outf.precision(15);
  outf << "TITLE = \"DATA FOR ROTOR 37 @ UPSTREAM STATION\"\n"
       << "VARIABLES = \"%Span\"\n"
       << "\"r\"\n"
       << "\"z\"\n"
       << "\"V_z\"\n"
       << "\"V_theta\"\n"
       << "\"V_theta_rel\"\n"
       << "\"Vabs\"\n"
       << "\"Vrel\"\n"
       << "\"Mabs\"\n"
       << "\"Mrel\"\n"
       << "\"Temp\"\n"
       << "\"Pressure\"\n"
       << "\"Density\"\n"
       << "\"Sound Speed\"\n"
       << "\"Mass Flux\"\n"
       << "ZONE T = \"UPSTREAM FLOW DATA\"\n"
       << "I=201, J=1, K=1, F=POINT\n"
       << "DT=(SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE)\n";

  for ( i=0; i<=200; i=i+1) {
    outf << i/2.0 << "  " 
         << rHub_u + (double(i)/2.00)/(HUNDRED)*(rShrd_u-rHub_u) << " " 
         << z_up << " "
         << getV_z_up(i/2.0) << "  " 
         << getV_theta_up(i/2.0) << "  " 
         << getVREL_theta_up(i/2.0) << "  " 
         << getMach_up(i/2.0)*getA_up(i/2.0) << "  " 
         << getMachREL_up(i/2.0)*getA_up(i/2.0) << "  "
         << getMach_up(i/2.0) << "  " 
         << getMachREL_up(i/2.0) << "  " 
         << getT_up(i/2.0) << "  " 
         << getP_up(i/2.0) << "  " 
         << getRho_up(i/2.0) << "  " 
         << getA_up(i/2.0) << "  " 
         << getV_z_up(i/2.0)*getRho_up(i/2.0) << "\n";
  }

  outf.close();

}

/**********************************************************************
 * Routine: outputTP_DownstreamFlowConditions                         *
 *                                                                    *
 * Output all downstream flow conditions (at each %span) to a file in *
 * TECPLOT FORMAT.                                                    *
 *                                                                    *
 **********************************************************************/
inline void NASARotor37::outputTP_DownstreamFlowConditions(char fname[256]) {

  int i;
  ofstream outf;

  outf.open(fname, ios::out);
  assert(!outf.fail());

  outf.precision(15);
  outf << "TITLE = \"DATA FOR ROTOR 37 @ DOWNSTREAM STATION\"\n"
       << "VARIABLES = \"%Span\"\n"
       << "\"r\"\n"
       << "\"z\"\n"
       << "\"V_z\"\n"
       << "\"V_theta\"\n"
       << "\"V_theta_rel\"\n"
       << "\"Vabs\"\n"
       << "\"Vrel\"\n"
       << "\"Mabs\"\n"
       << "\"Mrel\"\n"
       << "\"Temp\"\n"
       << "\"Pressure\"\n"
       << "\"Density\"\n"
       << "\"Sound Speed\"\n"
       << "\"Mass Flux\"\n"
       << "ZONE T = \"DOWNSTREAM FLOW DATA\"\n"
       << "I=201, J=1, K=1, F=POINT\n"
       << "DT=(SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE, SINGLE)\n";

  for ( i=0; i<=200; i=i+1 ) {
    outf << i/2.0 << "  " 
         << rHub_d + (double(i)/2.00)/(HUNDRED)*(rShrd_d-rHub_d) << " " 
         << z_d << " "
         << getV_z_down(i/2.0) << "  " 
         << getV_theta_down(i/2.0) << "  " 
         << getVREL_theta_down(i/2.0) << "  " 
         << getMach_down(i/2.0)*getA_down(i/2.0) << "  " 
         << getMachREL_down(i/2.0)*getA_down(i/2.0) << "  " 
         << getMach_down(i/2.0) << "  " 
         << getMachREL_down(i/2.0) << "  " 
         << getT_down(i/2.0) << "  " 
         << getP_down(i/2.0) << "  " 
         << getRho_down(i/2.0) << "  " 
         << getA_down(i/2.0) << "  " 
         << getV_z_down(i/2.0)*getRho_down(i/2.0) << "\n"; 
  }

  outf.close();

}

/**********************************************************************
 * Routine: itostr                                                    *
 **********************************************************************/
inline void NASARotor37::itostr(char *TargetString, unsigned theint) {
  TargetString+=1;
  for (; theint > 0; theint /= 10 ) {
    *TargetString = (char)((int)'0' + (theint - (theint / 10) * 10));
    TargetString--;
  }      
  TargetString[3]='\0';
}

/**********************************************************************
 * Routine: f                                                         *
 *                                                                    *
 * Calculate total mass flow using Splined values for p and t         *
 * upstream and downstream.                                           *
 *                                                                    *
 **********************************************************************/
inline double NASARotor37::f(double r) {
  return( getY (r, mach_u)[0].y * getY( r, ps_u )[0].y*sqrt(gmma*R*getY(r,t_u)[0].y) / R / getY(r,t_u)[0].y );
}

/**********************************************************************
 * Routine: r                                                         *
 **********************************************************************/
inline double NASARotor37::r(double A) {
  return ( sqrt(A/M_PI + rHub_u*rHub_u));
}

/**********************************************************************
 * Routine: getMassFlow_up                                            *
 **********************************************************************/
inline double NASARotor37::getMassFlow_up(int N) {
  int i;
  double dA;
  double mdot = 0;
  
  N=N*2;

  dA=M_PI*(rShrd_u*rShrd_u-rHub_u*rHub_u)/N;

  for ( i=1; i<=N/2; ++i ) {    
    mdot+=( f( r((2*i-2)*dA) ) + 4* f( r((2*i-1)*dA) ) +  f( r((2*i)*dA) ) );
  }

  mdot*= dA/3.0;

  return (mdot);

}

/**********************************************************************
 * Routine: f_d                                                       *
 **********************************************************************/
inline double NASARotor37::f_d(double r) {
  return( getY (r, mach_d)[0].y * getY( r, ps_d )[0].y*sqrt(gmma*R*getY(r,t_d)[0].y) / R / getY(r,t_d)[0].y );
}

/**********************************************************************
 * Routine: r_d                                                       *
 **********************************************************************/
inline double NASARotor37::r_d(double A) {
  return ( sqrt(A/M_PI + rHub_d*rHub_d));
}

/**********************************************************************
 * Routine: getMassFlow_down                                          *
 **********************************************************************/
inline double NASARotor37::getMassFlow_down(int N) {
  int i;
  double dA;
  double mdot = 0;
  
  N=N*2;

  dA=M_PI*(rShrd_d*rShrd_d-rHub_d*rHub_d)/N;

  for ( i=1; i<=N/2; ++i ) {    
    mdot+=( f_d( r_d((2*i-2)*dA) ) + 4* f_d( r_d((2*i-1)*dA) ) +  f_d( r_d((2*i)*dA) ) );
  }

  mdot*= dA/3.0;

  return (mdot);

}

/**********************************************************************
 * Routine: findRadius                                                *
 *                                                                    *
 * Return the radius for any given percent span at any axial position *
 * z.                                                                 *
 *                                                                    *
 **********************************************************************/
inline double NASARotor37::findRadius(double ps, double z) {
  double r;

  r=hubR(z) + (shroudR(z)-hubR(z))*ps / 100;

  //prevent round-off error when using this function in conjuntion with
  //downstream and upstream flow functions
  if(z==z_up) //upstream
    if(ps==0 && r<rHub_u)
      r=rHub_u;
    else if(ps==100 && r>rShrd_u)
      r=rShrd_u;
  if(z==z_d) //downstream
    if(ps==0 && r<rHub_d)
      r=rHub_d;
    else if(ps==100 && r>rShrd_d)
      r=rShrd_d;
  return (r);
}

/**********************************************************************
 * Routine: setgas                                                    *
 **********************************************************************/
inline void NASARotor37::setgas(char *string_ptr) {
  if (strcmp(string_ptr, "AIR") == 0) {
    gmma = GAMMA_AIR;
    R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
  } else if (strcmp(string_ptr, "A") == 0) {
    gmma = GAMMA_A;
    R = R_UNIVERSAL/(MOLE_WT_A*MILLI);
  } else if (strcmp(string_ptr, "CO") == 0) {
    gmma = GAMMA_CO;
    R = R_UNIVERSAL/(MOLE_WT_CO*MILLI);
  } else if (strcmp(string_ptr, "CO2") == 0) {
    gmma = GAMMA_CO2;
    R = R_UNIVERSAL/(MOLE_WT_CO2*MILLI);
  } else if (strcmp(string_ptr, "CH4") == 0) {
    gmma = GAMMA_CH4;
    R = R_UNIVERSAL/(MOLE_WT_CH4*MILLI);
  } else if (strcmp(string_ptr, "H") == 0) {
    gmma = GAMMA_H;
    R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
  } else if (strcmp(string_ptr, "H2") == 0) {
    gmma = GAMMA_H2;
    R = R_UNIVERSAL/(MOLE_WT_H2*MILLI);
  } else if (strcmp(string_ptr, "HE") == 0) {
    gmma = GAMMA_HE;
    R = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
  } else if (strcmp(string_ptr, "H2O") == 0) {
    gmma = GAMMA_H2O;
    R = R_UNIVERSAL/(MOLE_WT_H2O*MILLI);
  } else if (strcmp(string_ptr, "N2") == 0) {
    gmma = GAMMA_N2;
    R = R_UNIVERSAL/(MOLE_WT_N2*MILLI);
  } else if (strcmp(string_ptr, "O") == 0) {
    gmma = GAMMA_O;
    R = R_UNIVERSAL/(MOLE_WT_O*MILLI);
  } else if (strcmp(string_ptr, "O2") == 0) {
    gmma = GAMMA_O2;
    R = R_UNIVERSAL/(MOLE_WT_O2*MILLI);
  } else if (strcmp(string_ptr, "e") == 0) {
    gmma = GAMMA_e;
    R = R_UNIVERSAL/(MOLE_WT_e*MILLI);
  } else {
    gmma = GAMMA_AIR;
    R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
  } /* endif */
}

/**********************************************************************
 * Routine: findPSpan                                                 *
 *                                                                    *
 * Find the %span given an r value at any axial position z.           *
 *                                                                    *
 **********************************************************************/
inline double NASARotor37::findPSpan(double r, double z) {
  double span;
 
  span=(r-hubR(z))/(shroudR(z)-hubR(z))*100;
  if(span<0)
    span=0;
  if(fabs(span)<0.0001)
    span=0;
  if(span>100)
    span=100;
  return(span); 
}

/**********************************************************************
 * Routine: getBladeCS                                                *
 *                                                                    *
 * Returns the blade cross section at any percent span.               *
 *                                                                    *
 **********************************************************************/
inline Spline2D NASARotor37::getBladeCS(double pspan) {
  int i, j, sPos;
  double r; //initial guess for r
  Spline2D cs;
  Vector3D v1,v2,vInt;
  bladePt vIntc;
  
  //init cross-section spline
  cs.allocate(R37GEOM_NUMPTS*2-1); //for every file data pt there exists two physical
  cs.settype(SPLINE2D_QUINTIC);


  //do for each data point
  for ( j=R37GEOM_NUMPTS-1; j>=0; --j ) {
  
    //find two surfaces in data file that bound pspan
    for ( i=0; i<NUMBER_OF_BLADE_SURFACES_R37-1; ++i ) {
      if(pspan >= findPSpan(bladeData[i][j].r, bladeData[i][j].z) && 
	 pspan < findPSpan(bladeData[i+1][j].r, bladeData[i+1][j].z) ) {
	sPos=i;
      }
    }

    //fringe case
    if(pspan>=findPSpan(bladeData[NUMBER_OF_BLADE_SURFACES_R37-1][j].r,
			bladeData[NUMBER_OF_BLADE_SURFACES_R37-1][j].z))
       sPos=NUMBER_OF_BLADE_SURFACES_R37-2; 

    //POINT #1

    //convert to cartesian
    v1.x=bladeData[sPos][j].r*cos(bladeData[sPos][j].t1);
    v1.y=bladeData[sPos][j].r*sin(bladeData[sPos][j].t1);
    v1.z=bladeData[sPos][j].z;
    v2.x=bladeData[sPos+1][j].r*cos(bladeData[sPos+1][j].t1);
    v2.y=bladeData[sPos+1][j].r*sin(bladeData[sPos+1][j].t1);
    v2.z=bladeData[sPos+1][j].z;

            
    //calculate the interpolated point for the given r
    vInt=lInt(v1, v2, findPSpan(bladeData[sPos][j].r, bladeData[sPos][j].z), 
	      findPSpan(bladeData[sPos+1][j].r, bladeData[sPos+1][j].z), pspan);

    //convert back to cylindrical
    vIntc.z=vInt.z;
    vIntc.r=sqrt(vInt.x*vInt.x+vInt.y*vInt.y);
    vIntc.t1=atan(vInt.y/vInt.x);
  
    
    //store in spline var
    cs.Xp[R37GEOM_NUMPTS-j-1].x=vIntc.z; //independent variable is z(axial pos)
    cs.Xp[R37GEOM_NUMPTS-j-1].y=vIntc.r*vIntc.t1; //dependent var is r*theta
    cs.tp[R37GEOM_NUMPTS-j-1]=SPLINE2D_POINT_NORMAL;

    //POINT #2
    
    //convert to cartesian
    v1.x=bladeData[sPos][j].r*cos(bladeData[sPos][j].t2);
    v1.y=bladeData[sPos][j].r*sin(bladeData[sPos][j].t2);
    v1.z=bladeData[sPos][j].z;
    v2.x=bladeData[sPos+1][j].r*cos(bladeData[sPos+1][j].t2);
    v2.y=bladeData[sPos+1][j].r*sin(bladeData[sPos+1][j].t2);
    v2.z=bladeData[sPos+1][j].z;
      
    //calculate the interpolated point for the given r
    vInt=lInt(v1, v2, findPSpan(bladeData[sPos][j].r, bladeData[sPos][j].z), 
	      findPSpan(bladeData[sPos+1][j].r, bladeData[sPos+1][j].z), pspan);

    //convert back to cylindrical
    vIntc.z=vInt.z;
    vIntc.r=sqrt(vInt.x*vInt.x+vInt.y*vInt.y);
    vIntc.t1=atan(vInt.y/vInt.x);

    //store in spline var
    cs.Xp[R37GEOM_NUMPTS + j -1].x= vIntc.z; //independent variable is z(axial pos)
    cs.Xp[R37GEOM_NUMPTS + j -1].y=vIntc.r*vIntc.t1; //dependent var is r*theta
    cs.tp[R37GEOM_NUMPTS + j -1]=SPLINE2D_POINT_NORMAL; 
    
  }

    cs.pathlength();
    return(cs);
}

/**********************************************************************
 * Routine: lInt                                                      *
 **********************************************************************/
inline Vector3D NASARotor37::lInt(const Vector3D &v1, const Vector3D &v2, const double psMin,
				  const double psMax, const double pspan) {
  double A,B,C, ax, ay, az, L, r, r1, r2, hr;
  double tInt, t1, t2;
  Vector3D vInt;

  if(pspan==0)
      return(v1);
  if(pspan==100)
      return(v2);

  r1=sqrt(v1.x*v1.x+v1.y*v1.y);
  r2=sqrt(v2.x*v2.x+v2.y*v2.y);
  
  
  L=sqrt((v2.x-v1.x)*(v2.x-v1.x)+(v2.y-v1.y)*(v2.y-v1.y)+(v2.z-v1.z)*(v2.z-v1.z));
  ax=(v2.x-v1.x)/L;
  ay=(v2.y-v1.y)/L;
  az=(v2.z-v1.z)/L;
  A=ax*ax+ay*ay;
  B=2*(ax*v1.x+ay*v1.y);

  //initial guess
  r = r1 + (pspan-psMin) / (psMax-psMin) * (r2-r1);
 
  while(0<1) {
 
    C=v1.x*v1.x+v1.y*v1.y-r*r;

    //find two roots
    t1= (-B+sqrt(B*B-4*A*C))/(2*A);
    t2= (-B-sqrt(B*B-4*A*C))/(2*A);

    //to compensate for numerical error
    if(fabs(t1)<0.000001 && t1<0)
      t1=0;
    if(fabs(t1-L)<0.000001 && t1>L)
      t1=L;
    if(fabs(t2)<0.000001 && t2<0)
      t2=0;
    if(fabs(t2-L)<0.000001 && t2>L)
      t2=L;

    if(t1>=0 && t1<=L)
      tInt=t1;
    else if(t2>=0 && t2<=L)
      tInt=t2;
    else {
      cout << "ERROR: INTERPOLATION ERROR!" << "\n";
      cout << "L:" << L << " t1:" << t1 << " t2:" << t2 << "\n";
      cout.flush();
      assert(1<0);
    }

    vInt.x=v1.x+tInt*ax;
    vInt.y=v1.y+tInt*ay;
    vInt.z=v1.z+tInt*az;

    if(fabs(findPSpan(sqrt(vInt.x*vInt.x + vInt.y*vInt.y), vInt.z)-pspan)<0.0000001)
      return(vInt);

    //re-calc r
    hr=hubR(vInt.z);
    r=pspan/100*(shroudR(vInt.z)-hr)+hr;
    
  }     

}

/**********************************************************************
 * Routine: allocateTmpBlade                                          *
 *                                                                    *
 * Allocate memory for an array to store the interpolated geometry of *
 * one blade.                                                         *
 *                                                                    *
 **********************************************************************/
inline void NASARotor37::allocateTmpBlade(int surfaces, int points) {
  int i;
  bladeTmp=new bladePt*[surfaces];
  for ( i=0; i<surfaces;++i )
    bladeTmp[i]=new bladePt[points];
}

/**********************************************************************
 * Routine: deallocateTmpBlade                                        *
 *                                                                    *
 * Deallocate memory for temporary variable 'bladeTmp'.               *
 *                                                                    *
 **********************************************************************/
inline void NASARotor37::deallocateTmpBlade(int surfaces, int points) {
  int i;
  for ( i=0; i<surfaces; ++i ) {
      delete []bladeTmp[i];
    bladeTmp[i]=NULL;
  }
  delete []bladeTmp;
  bladeTmp=NULL;
}

/**********************************************************************
 * Routine:  outputTP_Geometry                                        *
 **********************************************************************/
inline void NASARotor37::outputTP_Geometry(char filename[]) {
  int i, j, k, l, sPos, nS, nP;
  double r, theta_shift, z, dz, t, dt;
  Vector3D v1, v2;
  Vector3D cs[R37GEOM_NUMPTS*2-1]; //temporary array to help with organization of pts
  ofstream TPfile;

  TPfile.open(filename, ios::out);
  assert(!TPfile.fail());
 
  nS=21;  //number of cross-sections  per blade = every 5% span
  nP=R37GEOM_NUMPTS*2-2; //number of points per surface

  //allocate memory to hold a copy of the interpolated blade geometry - 'bladeTmp'
  allocateTmpBlade(nS, nP);
  
  cout.width(4);
  cout.fill(' ');
  cout << "\n Generating 3D TECPLOT geometry file for NASA Rotor 37."
       << "\n  Contour data span interval: " << 100/(nS-1) <<"% +-0.001.";

  //write file header
  TPfile << "TITLE = \"3D PLOT OF NASA ROTOR 37\"\n"
	 << "VARIABLES = \"x\"\n"
	 << "\"y\"\n"
	 << "\"z\"\n"
         << "\"R*theta\"\n"
         << "\"R\"\n"
         << "\"theta\"\n\n";
  
  //generate the first blade
  cout << "\n  Creating Blade 1 ";
      
  TPfile << "\nZONE T= \"Blade 1\"\n"
	 << "F=FEPOINT, N=" << nP*nS << " E=" << nP*(nS-1) << "\n"
	 << "ET=quadrilateral\n\n";

  theta_shift=0;
      
  //write one surface for each %span
  for ( k=0; k<=100; k=k+100/(nS-1) ) {

    cout << ".";
    cout.flush();

    //do for each data point
    for ( j=R37GEOM_NUMPTS-1; j>=0; --j ) {
	   
      //find two surfaces in data file that bound pspan
      for ( i=0; i<NUMBER_OF_BLADE_SURFACES_R37-1; ++i ) {
	if(k >= findPSpan(bladeData[i][j].r, bladeData[i][j].z) && 
	   k < findPSpan(bladeData[i+1][j].r, bladeData[i+1][j].z) ) {
	  sPos=i;
	}
      }

      //fringe case
      if(k>=findPSpan(bladeData[NUMBER_OF_BLADE_SURFACES_R37-1][j].r,
		      bladeData[NUMBER_OF_BLADE_SURFACES_R37-1][j].z))
	sPos=NUMBER_OF_BLADE_SURFACES_R37-2; 
         
      //POINT #1

      //convert to cartesian
      v1.x=bladeData[sPos][j].r*cos(bladeData[sPos][j].t1+theta_shift);
      v1.y=bladeData[sPos][j].r*sin(bladeData[sPos][j].t1+theta_shift);
      v1.z=bladeData[sPos][j].z;
      v2.x=bladeData[sPos+1][j].r*cos(bladeData[sPos+1][j].t1+theta_shift);
      v2.y=bladeData[sPos+1][j].r*sin(bladeData[sPos+1][j].t1+theta_shift);
      v2.z=bladeData[sPos+1][j].z;

         
      //calculate the interpolated point for the given r and store
      cs[R37GEOM_NUMPTS-j-1]=lInt(v1, v2,
		  findPSpan(bladeData[sPos][j].r, bladeData[sPos][j].z), 
		  findPSpan(bladeData[sPos+1][j].r, bladeData[sPos+1][j].z), 
		  k);

      //POINT #2

      //convert to cartesian
      v1.x=bladeData[sPos][j].r*cos(bladeData[sPos][j].t2+theta_shift);
      v1.y=bladeData[sPos][j].r*sin(bladeData[sPos][j].t2+theta_shift);
      v1.z=bladeData[sPos][j].z;
      v2.x=bladeData[sPos+1][j].r*cos(bladeData[sPos+1][j].t2+theta_shift);
      v2.y=bladeData[sPos+1][j].r*sin(bladeData[sPos+1][j].t2+theta_shift);
      v2.z=bladeData[sPos+1][j].z;
      

      //calculate the interpolated point for the given r and store
      cs[R37GEOM_NUMPTS + j -1]=lInt(v1, v2,
			     findPSpan(bladeData[sPos][j].r, bladeData[sPos][j].z), 
			     findPSpan(bladeData[sPos+1][j].r, bladeData[sPos+1][j].z), 
			     k);

    }
	
    //write the contour data to file and save cross-section data
    for ( j=0; j<nP; ++j ) {
      TPfile << cs[j].x << " "
	     << cs[j].y << " "
	     << cs[j].z << " "
	     << sqrt(cs[j].x*cs[j].x + cs[j].y*cs[j].y)*atan2(cs[j].y, cs[j].x) << " "
	     << sqrt(cs[j].x*cs[j].x + cs[j].y*cs[j].y) << " "
	     << atan2(cs[j].y, cs[j].x) << "\n";
      
      bladeTmp[k/(100/(nS-1))][j].z=cs[j].z;
      bladeTmp[k/(100/(nS-1))][j].r=sqrt(cs[j].x*cs[j].x + cs[j].y*cs[j].y);
      bladeTmp[k/(100/(nS-1))][j].t1=atan2(cs[j].y, cs[j].x);
    }

  }

  //write connectivity data for this blade
  for ( i=0; i<nS-1;++i ) {
    for ( j=1; j<=nP-1;++j )
      TPfile << i*nP+j   << ", "
	     << i*nP+j+1 << ", " 
	     << (i+1)*nP+j+1 << ", "
	     << (i+1)*nP+j << "\n";
	
    TPfile << i*nP+nP   << ", "
           << i*nP+1    << ", " 
           << (i+1)*nP+1 << ", "
           << (i+1)*nP+nP << "\n";
	
  }
  cout << " Done!";

  //write other blades 21 blades
  for ( l=1; l<num_blades; ++l ) {

    cout << "\n  Creating Blade " << l+1 << " ";

    TPfile << "\nZONE T= \"Blade " << l+1 << "\"\n"
	   << "F=FEPOINT, N=" << nP*nS << " E=" << nP*(nS-1) << "\n"
	   << "ET=quadrilateral\n\n";

    theta_shift=2*M_PI/num_blades*l;

    //write one surface for each %span
    for ( k=0; k<nS; ++k ) {

      cout << ".";
      cout.flush();

      //write the contour data to file 
      for ( j=0; j<nP; ++j ) {
	TPfile << bladeTmp[k][j].r*cos(bladeTmp[k][j].t1+theta_shift) << " "
	       << bladeTmp[k][j].r*sin(bladeTmp[k][j].t1+theta_shift) << " "
	       << bladeTmp[k][j].z << " "
	       << bladeTmp[k][j].r*(bladeTmp[k][j].t1+theta_shift) << " "		 
	       << bladeTmp[k][j].r << " "
	       << bladeTmp[k][j].t1+theta_shift << "\n";
      }
      
    }
    
    //write connectivity data for this blade
    for ( i=0; i<nS-1;++i ) {
      for ( j=1; j<=nP-1;++j )
	TPfile << i*nP+j   << ", "
	       << i*nP+j+1 << ", " 
	       << (i+1)*nP+j+1 << ", "
	       << (i+1)*nP+j << "\n";
	
      TPfile << i*nP+nP   << ", "
	     << i*nP+1    << ", " 
	     << (i+1)*nP+1 << ", "
	     << (i+1)*nP+nP << "\n";
    }
  
    cout << " Done!";
  }
    

  cout << "\n  Creating hub...";
  cout.flush();

  //Create new zone for hub
  TPfile << "\nZONE T= \"hub\"\n"
	 << "F=FEPOINT, N=" << nP*nS  << " E=" << nP*(nS-1) << "\n"
	 << "ET=quadrilateral\n"
	 << "D=(FECONNECT)\n\n";

  //write pts
  z=hub.Xp[0].x;
  dz=(hub.Xp[hub.np-1].x-hub.Xp[0].x)/(nS-1);
  for ( i=1; i<=nS; ++i ) {
    t=0;
    dt=2*M_PI/nP;
    for ( j=1; j<=nP; ++j ) {
      TPfile << hubR(z)*cos(t) << " "
	     << hubR(z)*sin(t) << " "
	     << z << " "
	     << hubR(z)*t << " "
	     << hubR(z) << " "
	     << t << "\n";
      t+=dt;
    }
    z+=dz;
  }

  cout << "\n  Creating shroud...";
  cout.flush();

  //Create new zone for shroud
  TPfile << "\nZONE T= \"shroud\"\n"
	 << "F=FEPOINT, N=" << nP*nS << " E=" << nP*(nS-1) << "\n"
	 << "ET=quadrilateral\n"
	 << "D=(FECONNECT)\n\n";
  
  //write pts
  z=shroud.Xp[0].x;
  dz=(shroud.Xp[shroud.np-1].x-shroud.Xp[0].x)/(nS-1);
  for ( i=1; i<=nS; ++i ) {
    t=0;
    dt=2*M_PI/nP;
    for ( j=1; j<=nP; ++j ) {
      TPfile << shroudR(z)*cos(t) << " "
	     << shroudR(z)*sin(t) << " "
	     << z << " "
	     << shroudR(z)*t << " "
	     << shroudR(z) << " "
	     << t << "\n";
      t+=dt;
    }
    z+=dz;
  }

  TPfile.close();

  cout << "\n 3D NASA Rotor 37 Rotor geometry file complete!";

  deallocateTmpBlade(nS, nP);

}

/**********************************************************************
 * Routine: getCamberBlade                                            *
 *                                                                    *
 * Return the camber line of a blade at any %span where n = number of *
 * points to generate on camber line.                                 *
 *                                                                    *
 **********************************************************************/
inline Spline2D NASARotor37::getCamberBlade(double pspan, int n) {
  int i, imin, imax;
  Spline2D camber, cs, csTop, csBot;
  double maxZ, minZ;
  double m, dz, z, z1, z2;
  double t, b;
 
  //get rotated cross section and index of leading and trailing edge
  cs=transformCS(pspan, m, imin, imax);
  cs.pathlength();
 
  camber.allocate(n+1);
  camber.settype(SPLINE2D_QUINTIC);
  
  //get min and max Z points of rotated airfoil
  minZ=cs.Xp[imax].x;
  maxZ=cs.Xp[imin].x;

  //split cs spline into two: one for the top of blade, one for the bottom
  csTop.allocate(imax-imin+1);
  csTop.settype(SPLINE2D_QUINTIC);
  
  assert(imax>imin); // make sure

  for ( i=imax; i>=imin; --i )  
    csTop.Xp[imax-i]=cs.Xp[i];
  
  //get bottom of spline
  csBot.allocate(cs.np-(imax-imin));
  csBot.settype(SPLINE2D_QUINTIC);

  for ( i=imax; i<cs.np; ++i )
    csBot.Xp[i-imax]=cs.Xp[i];
  
  for ( i=1; i<=imin; ++i )
    csBot.Xp[cs.np-imax+i-1]=cs.Xp[i];
  
  //set point types
  for ( i=1; i<csTop.np-1; ++i )
    csTop.tp[i]=SPLINE2D_POINT_NORMAL;
  for ( i=1; i<csBot.np-1; ++i )
    csBot.tp[i]=SPLINE2D_POINT_NORMAL;
  csTop.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  csTop.tp[csTop.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  csBot.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  csBot.tp[csBot.np-1]=SPLINE2D_POINT_SHARP_CORNER;

  csTop.pathlength();
  csBot.pathlength();

  dz=(maxZ-minZ)/n;
  z=minZ;

  for ( i=0; i<=n-1; ++i ) {
    b=getY(z, csBot)[0].y;
    t=getY(z, csTop)[0].y;
    camber.Xp[i].x=z;
    camber.Xp[i].y=(t + b)/2.0;
    z+=dz;
  }

  //large values of n cause round-off error so set last point manually
  camber.Xp[n]=cs.Xp[imin];
  

  rotateCS(camber, -m);
    
  //set point types for spline
  camber.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  camber.tp[n]=SPLINE2D_POINT_SHARP_CORNER;

  for ( i=1 ; i<n; ++i )
    camber.tp[i]=SPLINE2D_POINT_NORMAL;

  camber.pathlength();
  
  return(camber);
}

/**********************************************************************
 * Routine: getCamberTrail                                            *
 **********************************************************************/
inline Spline2D NASARotor37::getCamberTrail(double pspan, int n, double max) {
  int i;
  double m, dz, maxZ, minZ, z1, z2;
  Spline2D camber, camberT;

  camberT.allocate(n);
  camberT.settype(SPLINE2D_LINEAR);
  camber=getCamberBlade(pspan, 35);

  //extrapolate camber line beyond blade

  minZ=camber.Xp[0].x;
  maxZ=camber.Xp[camber.np-1].x;
  dz=(maxZ-minZ)/100.0;

  //find slope at trailing edge using second-order finite difference
  m=(3*getY(maxZ, camber)[0].y-4*getY(maxZ-2*dz, camber)[0].y+getY(maxZ-4*dz, camber)[0].y)/(4*dz);

  //extrapolate camber line from trailing edge
  z1=maxZ;
  z2=max;
  dz=(z2-z1)/(n-1);

  for ( i=0; i<n; ++i ) {
    camberT.Xp[i].x=z1+dz*i;
    camberT.Xp[i].y = m*(camberT.Xp[i].x-camber.Xp[camber.np-1].x)+camber.Xp[camber.np-1].y;
  }

  //set point types for spline
  camberT.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  camberT.tp[n-1]=SPLINE2D_POINT_SHARP_CORNER;

  for ( i=1; i<n-1; ++i )
    camberT.tp[i]=SPLINE2D_POINT_NORMAL;

  camberT.pathlength();
  
  return(camberT);
  
}

/**********************************************************************
 * Routine: getCamberLead                                             *
 *                                                                    *
 * Return the camber line extended from the leading edge.             *
 *                                                                    *
 **********************************************************************/
inline Spline2D NASARotor37::getCamberLead(double pspan, int n, double min) {
  int i;
  double m, dz, maxZ, minZ, z1, z2;
  Spline2D camber, camberL;

  camberL.allocate(n);
  camberL.settype(SPLINE2D_LINEAR);

  camber=getCamberBlade(pspan, 35);
  
  //extrapolate camber line beyond blade

  minZ=camber.Xp[0].x;
  maxZ=camber.Xp[camber.np-1].x;
  dz=(maxZ-minZ)/100.0;

  //find slope at trailing edge using second-order finite difference
  
  m=(-3*getY(minZ, camber)[0].y+4*getY(minZ+2*dz, camber)[0].y-getY(minZ+4*dz, camber)[0].y)/(4*dz);
 
  //extrapolate camber line from leading edge
  z2=minZ;
  z1=min;
  dz=(z2-z1)/(n-1);

  for ( i=0; i<n; ++i ) {
    camberL.Xp[i].x=z1+dz*(i);
    camberL.Xp[i].y = m*(camberL.Xp[i].x-camber.Xp[0].x)+camber.Xp[0].y;
  }

  //set point types for spline
  camberL.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  camberL.tp[n-1]=SPLINE2D_POINT_SHARP_CORNER;

  for ( i=1; i<n-1; ++i )
    camberL.tp[i]=SPLINE2D_POINT_NORMAL;

  camberL.pathlength();
  
  return(camberL);
  
}

/**********************************************************************
 * Routine: getCamber                                                 *
 *                                                                    *
 * Return the camber line extended from the trailing edge.            *
 *                                                                    *
 **********************************************************************/
inline Spline2D NASARotor37::getCamber(double pspan, int n, double min, double max) {
  int i;
  Spline2D blade, lead, trail, camber;

  assert(max>min);

  blade=getCamberBlade(pspan,n);
  trail=getCamberTrail(pspan,n,max);
  lead=getCamberLead(pspan,n,min);

  camber.allocate(3*n-1);
  camber.settype(SPLINE2D_QUINTIC);

  for ( i=0; i<n-1; ++i )
    camber.Xp[i]=lead.Xp[i];

  for ( i=n-1; i<2*n; ++i )
    camber.Xp[i]=blade.Xp[i-(n-1)];

  for ( i=2*n;  i<3*n-1; ++i )
    camber.Xp[i]=trail.Xp[i-(2*n)+1];

  //set point types for spline
  camber.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  camber.tp[n-1]=SPLINE2D_POINT_SHARP_CORNER;

  for ( i=1; i<camber.np-1; ++i )
    camber.tp[i]=SPLINE2D_POINT_NORMAL;

  camber.pathlength();
  
  return(camber);
  
}

/**********************************************************************
 * Routine: getCamberAndShift                                         *
 *                                                                    *
 * Return the full camber line and shift it by any angle.             *
 *                                                                    *
 **********************************************************************/
inline Spline2D NASARotor37::getCamberAndShift(double pspan, int n, double min, double max, double angle) {
  int i;
  Spline2D camber;
  double theta, r;

  camber=getCamber(pspan, n, min, max);

  for ( i=0; i<camber.np; ++i ) {
    r=findRadius(pspan, camber.Xp[i].x);
    theta=camber.Xp[i].y/r;
    theta+=angle;
    camber.Xp[i].y=theta*r;
  }

  return(camber);
}

/**********************************************************************
 * Routine: len                                                       *
 *                                                                    *
 * Returns the distance between two points.                           *
 *                                                                    *
 **********************************************************************/
inline double NASARotor37::len(Vector2D v1, Vector2D v2) {
  return(sqrt( (v1.x-v2.x)*(v1.x-v2.x)+(v1.y-v2.y)*(v1.y-v2.y) ));
}

/**********************************************************************
 * Routine: findLT                                                    *
 *                                                                    *
 * Finds leading and trailing edge of cross section (imax & imin) and *
 * returns the slope of the connecting chord.                         *
 *                                                                    *
 **********************************************************************/
inline double NASARotor37::findLT(double pspan, int &imin, int &imax) {
  Spline2D cs;
  int i, j, max, min;
  double Lmax=0;

  cs=getBladeCS(pspan);

  for ( i=0; i<30; ++i ) 
    for ( j=70; j<130; ++j ) {
      if(len(cs.Xp[i], cs.Xp[j]) > Lmax) {
	Lmax=len(cs.Xp[i], cs.Xp[j]);
	max=i;
	min=j;
      }
    }

  for ( i=170; i<cs.np; ++i ) 
    for ( j=70; j<130; ++j ) {
      if(len(cs.Xp[i], cs.Xp[j]) > Lmax) {
	Lmax=len(cs.Xp[i], cs.Xp[j]);
	max=i;
	min=j;
      }
    }
  
  imin=max;
  imax=min;
  return((cs.Xp[max].y-cs.Xp[min].y)/(cs.Xp[max].x-cs.Xp[min].x));

}

/**********************************************************************
 * Routine: rotateCS                                                  *
 *                                                                    *
 * Rotates a Spline2D by 'rad' radians around origin - positive is    *
 * counter clockwise.                                                 *
 *                                                                    *
 **********************************************************************/
inline void NASARotor37::rotateCS(Spline2D &cs, double rad) {
  int i;
  for ( i=0; i<cs.np; ++i )
    cs.Xp[i]=rotate(cs.Xp[i], rad); 
} 

/**********************************************************************
 * Routine: rotate                                                    *
 *                                                                    *
 * Rotate a vector 'v' by 'rad' radians around origin - positive is   *
 * counter-clockwise.                                                 *
 *                                                                    *
 **********************************************************************/
inline Vector2D NASARotor37::rotate(Vector2D v, const double &rad) {
  return(Vector2D(v.x*cos(rad)-v.y*sin(rad), v.x*sin(rad)+v.y*cos(rad)));
}

/**********************************************************************
 * Routine: transformCS                                               *
 *                                                                    *
 * Get a spline of the cross-section at pspan and rotate it so that   *
 * the chord connecting the leading and trailing edge is horizontal.  *
 *                                                                    *
 **********************************************************************/
inline Spline2D NASARotor37::transformCS(double pspan, double &m, int &imin, int &imax) {
  Spline2D cs;

  m=-atan(findLT(pspan, imin, imax));
  cs=getBladeCS(pspan);
  
  rotateCS(cs, m);
 
  return(cs);
}

/**********************************************************************
 * Routine: genMeshO                                                  *
 *                                                                    *
 * Generate a mesh for the cross-section of the blade at any %span.   *
 *                                                                    *
 **********************************************************************/
inline Grid2D_Quad_Block NASARotor37::genMeshO(double pspan, double zMin, double zMax,
					       const int Number_of_Cells_Idir,
					       const int Number_of_Cells_Jdir,
					       const int Number_of_Ghost_Cells,
					       int StretchFcnI, double BetaI, double TauI,
					       int StretchFcnJ, double BetaJ, double TauJ,
					       int NbndType, int EbndType, int SbndType, int WbndType,
					       int OrthN, int OrthE, int OrthS, int OrthW, int n) {

  Grid2D_Quad_Block Grid;
  Spline2D upperN, lowerN, cs, camberTrail;
  Spline2D Bnd_Spline_North, Bnd_Spline_South, Bnd_Spline_East, Bnd_Spline_West;
  Vector2D swap;  
  int i, iLead, iTrail;
   
  /* Create the splines defining the north, south,
     east, and west boundaries of the grid. */

  // **** NORTH ****
  //create north spline from two full camber lines and sides
  camberTrail=getCamberTrail(pspan, 70, zMax);
  Bnd_Spline_North.allocate(6*70);
  Bnd_Spline_North.settype(SPLINE2D_QUINTIC);
  upperN=getCamberAndShift(pspan, 70, zMin, zMax, M_PI/num_blades);
  lowerN=getCamberAndShift(pspan, 70, zMin, zMax, -M_PI/num_blades);

  //set first point equal to the last point on the camber line extension
  Bnd_Spline_North.Xp[0]=camberTrail.Xp[camberTrail.np-1];
  
  //read in bottom camber line (in reverse order)
  for ( i=0; i<lowerN.np; ++i )
    Bnd_Spline_North.Xp[1+i]=lowerN.Xp[lowerN.np-i-1];

  //read in upper camber line
  for ( i=0; i<upperN.np; ++i )
    Bnd_Spline_North.Xp[1+lowerN.np+i]=upperN.Xp[i];
  
  //set last point equal to first
  Bnd_Spline_North.Xp[6*70-1]=Bnd_Spline_North.Xp[0];
  
  //set point types
  Bnd_Spline_North.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  Bnd_Spline_North.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  for ( i=1; i<lowerN.np-1; ++i )
    Bnd_Spline_North.tp[1+i]=SPLINE2D_POINT_NORMAL;
  
  Bnd_Spline_North.tp[1+lowerN.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  Bnd_Spline_North.tp[1+lowerN.np]=SPLINE2D_POINT_SHARP_CORNER;

  for ( i=1; i<upperN.np-1; ++i )
    Bnd_Spline_North.tp[1+lowerN.np+i]=SPLINE2D_POINT_NORMAL;

  Bnd_Spline_North.tp[1+lowerN.np + upperN.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  Bnd_Spline_North.tp[1+lowerN.np + upperN.np]=SPLINE2D_POINT_SHARP_CORNER;

  //calculate pathlength
  Bnd_Spline_North.pathlength();

  // **** EAST ****
  
  Copy_Spline(Bnd_Spline_East, camberTrail);

  // **** WEST ****
  
  Copy_Spline(Bnd_Spline_West, camberTrail);

  // **** SOUTH ****
  
  findLT(pspan, iTrail, iLead);
  cs=getBladeCS(pspan);

  //re-organize pts so they ascend in the clockwise direction
  for ( i=0; i<R37GEOM_NUMPTS; ++i ) {
    swap=cs.Xp[2*R37GEOM_NUMPTS-2-i];
    cs.Xp[2*R37GEOM_NUMPTS-2-i]=cs.Xp[i];
    cs.Xp[i]=swap;
  }

  //init spline
  Bnd_Spline_South.allocate(cs.np);
  Bnd_Spline_South.settype(SPLINE2D_QUINTIC);

  //change iTrail accordingly
  iTrail=cs.np-1-iTrail;

  //read in cross-section data starting with the trailing edge pt
  for ( i=iTrail; i<cs.np; ++i )
    Bnd_Spline_South.Xp[i-iTrail]=cs.Xp[i];

  //read in the rest if the airfoil data
  for ( i=1; i<=iTrail; ++i )
    Bnd_Spline_South.Xp[(cs.np-iTrail)-1+i]=cs.Xp[i];

  //set end pt types
  Bnd_Spline_South.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  Bnd_Spline_South.tp[Bnd_Spline_South.np-1]=SPLINE2D_POINT_SHARP_CORNER;

  //set other pt types
  for ( i=1; i<Bnd_Spline_South.np-1; ++i )
    Bnd_Spline_South.tp[i]=SPLINE2D_POINT_NORMAL;

  //re-calculate pathlength
  Bnd_Spline_South.pathlength();

  //set north, south, east and west boundary types for each boundary point
  for ( i=0; i<Bnd_Spline_North.np; ++i )
    Bnd_Spline_North.bc[i]=NbndType;

  for ( i=0; i<Bnd_Spline_South.np; ++i )
    Bnd_Spline_South.bc[i]=SbndType;
  
  for ( i=0; i<Bnd_Spline_East.np; ++i )
    Bnd_Spline_East.bc[i]=EbndType;
  
  for ( i=0; i<Bnd_Spline_West.np; ++i )
    Bnd_Spline_West.bc[i]=WbndType;
  
  /* Create the 2D quadrilateral grid block. */

  Create_Quad_Block(Grid,
                    Bnd_Spline_North,
                    Bnd_Spline_South,
                    Bnd_Spline_East,
                    Bnd_Spline_West,
                    Number_of_Cells_Idir,
		    Number_of_Cells_Jdir,
		    Number_of_Ghost_Cells,
                    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                    StretchFcnI, 
                    BetaI, 
                    TauI,
		    StretchFcnJ, 
                    BetaJ, 
                    TauJ,
                    OrthN, 
                    OrthE, 
                    OrthS, 
                    OrthW);

  /* Smooth the 2D quadrilateral grid block. */

  Smooth_Quad_Block(Grid, n);
  
  /* Deallocate the memory for the boundary splines. */
  
  Bnd_Spline_North.deallocate();
  Bnd_Spline_South.deallocate();
  Bnd_Spline_East.deallocate();
  Bnd_Spline_West.deallocate();

  /* Return the grid. */

  return(Grid);

}

/**********************************************************************
 * Routine: genMeshC                                                  *
 **********************************************************************/
inline Grid2D_Quad_Block NASARotor37::genMeshC(double pspan, double zMin, double zMax,
					       const int Number_of_Cells_Idir,
					       const int Number_of_Cells_Jdir,
					       const int Number_of_Ghost_Cells,
					       int StretchFcnI, double BetaI, double TauI,
					       int StretchFcnJ, double BetaJ, double TauJ,
					       int NbndType, int EbndType, int SbndType, int WbndType,
					       int OrthN, int OrthE, int OrthS, int OrthW, 
					       int n) {

  Grid2D_Quad_Block Grid;
  Spline2D upperN, lowerN, camberTrail, cs;
  Spline2D Bnd_Spline_North, Bnd_Spline_South, Bnd_Spline_East, Bnd_Spline_West;
  Vector2D swap;
  int i, iTrail, iLead;
 
  /* Create the splines defining the north, south,
     east, and west boundaries of the grid. */

  // **** NORTH ****
  //create north spline from two full camber lines and sides
  camberTrail=getCamberTrail(pspan, 70, zMax); 
  Bnd_Spline_North.allocate(6*70-2);
  Bnd_Spline_North.settype(SPLINE2D_QUINTIC);
  upperN=getCamberAndShift(pspan, 70, zMin, zMax, M_PI/num_blades);
  lowerN=getCamberAndShift(pspan, 70, zMin, zMax, -M_PI/num_blades);
    
  //read in bottom camber line (in reverse order)
  for ( i=0; i<lowerN.np; ++i )
    Bnd_Spline_North.Xp[i]=lowerN.Xp[lowerN.np-i-1];

  //read in upper camber line
  for ( i=0; i<upperN.np; ++i )
    Bnd_Spline_North.Xp[lowerN.np+i]=upperN.Xp[i];
      
  //set point types
  Bnd_Spline_North.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<lowerN.np-1; ++i )
    Bnd_Spline_North.tp[i]=SPLINE2D_POINT_NORMAL;
  Bnd_Spline_North.tp[lowerN.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  Bnd_Spline_North.tp[lowerN.np]=SPLINE2D_POINT_SHARP_CORNER;

  for ( i=1; i<upperN.np-1; ++i )
    Bnd_Spline_North.tp[lowerN.np+i]=SPLINE2D_POINT_NORMAL;
  Bnd_Spline_North.tp[lowerN.np + upperN.np-1]=SPLINE2D_POINT_SHARP_CORNER;

  //calculate pathlength
  Bnd_Spline_North.pathlength();

   // **** EAST ****
  Bnd_Spline_East.allocate(2);
  Bnd_Spline_East.settype(SPLINE2D_LINEAR);

  //set end points
  Bnd_Spline_East.Xp[0]=camberTrail.Xp[camberTrail.np-1];
  Bnd_Spline_East.Xp[1]=Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
  
  //set pt types
  Bnd_Spline_East.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  Bnd_Spline_East.tp[1]=SPLINE2D_POINT_SHARP_CORNER;

  //calculate pathlength
  Bnd_Spline_East.pathlength();

  // **** WEST ****
  
  Bnd_Spline_West.allocate(2);
  Bnd_Spline_West.settype(SPLINE2D_LINEAR);

  //set end points
  Bnd_Spline_West.Xp[0]=camberTrail.Xp[camberTrail.np-1];
  Bnd_Spline_West.Xp[1]=Bnd_Spline_North.Xp[0];
  
  //set pt types
  Bnd_Spline_West.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  Bnd_Spline_West.tp[1]=SPLINE2D_POINT_SHARP_CORNER;

  //calculate pathlength
  Bnd_Spline_West.pathlength();

  // **** SOUTH ****
  
  //find the index for the trailing edge
  findLT(pspan, iTrail, iLead);

  //get cross section spline
  cs=getBladeCS(pspan);

  //re-organize pts so they ascend in the clockwise direction
  for ( i=0; i<R37GEOM_NUMPTS; ++i ) {
    swap=cs.Xp[2*R37GEOM_NUMPTS-2-i];
    cs.Xp[2*R37GEOM_NUMPTS-2-i]=cs.Xp[i];
    cs.Xp[i]=swap;
  }

  //change iTrail accordingly
  iTrail=cs.np-1-iTrail;

  //allocate space
  Bnd_Spline_South.allocate(cs.np+2*70-2);
  
  //set type
  Bnd_Spline_South.settype(SPLINE2D_LINEAR);

  //read in camber line in reverse order
  for ( i=0; i<camberTrail.np-1; ++i )
    Bnd_Spline_South.Xp[i]=camberTrail.Xp[camberTrail.np-1-i];

  //read in cross-section data starting with the trailing edge pt
  for ( i=iTrail; i<cs.np; ++i )
    Bnd_Spline_South.Xp[camberTrail.np-1+i-iTrail]=cs.Xp[i];

  //read in the rest if the airfoil data
  for ( i=1; i<=iTrail; ++i )
    Bnd_Spline_South.Xp[camberTrail.np-1+(cs.np-iTrail)-1+i]=cs.Xp[i];

  //read in the camber line
  for ( i=1; i<camberTrail.np; ++i )
    Bnd_Spline_South.Xp[camberTrail.np+cs.np+i-2]=camberTrail.Xp[i];
   
  //set end pt types
  Bnd_Spline_South.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  Bnd_Spline_South.tp[Bnd_Spline_South.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  Bnd_Spline_South.tp[69]=SPLINE2D_POINT_SHARP_CORNER;
  Bnd_Spline_South.tp[207]=SPLINE2D_POINT_SHARP_CORNER;

  //set other point types
  for ( i=1; i<Bnd_Spline_South.np-1; ++i )
    Bnd_Spline_South.tp[i]=SPLINE2D_POINT_NORMAL;

  //calculate pathlength
  Bnd_Spline_South.pathlength();
 
  //set north, south, east and west boundary types for each boundary point
  for ( i=0; i<Bnd_Spline_North.np; ++i )
    Bnd_Spline_North.bc[i]=NbndType;

  for ( i=0; i<Bnd_Spline_South.np; ++i )
    Bnd_Spline_South.bc[i]=SbndType;
  
  for ( i=0; i<Bnd_Spline_East.np; ++i )
    Bnd_Spline_East.bc[i]=EbndType;
  
  for ( i=0; i<Bnd_Spline_West.np; ++i )
    Bnd_Spline_West.bc[i]=WbndType;
  
  /* Create the 2D quadrilateral grid block. */

  Create_Quad_Block(Grid,
                    Bnd_Spline_North,
                    Bnd_Spline_South,
                    Bnd_Spline_East,
                    Bnd_Spline_West,
                    Number_of_Cells_Idir,
		    Number_of_Cells_Jdir,
		    Number_of_Ghost_Cells,
                    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                    StretchFcnI, 
                    BetaI, 
                    TauI,
		    StretchFcnJ, 
                    BetaJ, 
                    TauJ,
                    OrthN, 
                    OrthE, 
                    OrthS, 
                    OrthW);

  /* Smooth the 2D quadrilateral grid block. */

  Smooth_Quad_Block(Grid, n);

  /* Deallocate the memory for the boundary splines. */
  
  Bnd_Spline_North.deallocate();
  Bnd_Spline_South.deallocate();
  Bnd_Spline_East.deallocate();
  Bnd_Spline_West.deallocate();

  /* Return the grid. */

  return(Grid);

}

/**********************************************************************
 * Routine: genMeshBlock                                              *
 **********************************************************************/
inline void NASARotor37::genMeshBlock(Grid2D_Quad_Block &block,
				      int n,
				      const int Number_of_Cells_Idir,
				      const int Number_of_Cells_Jdir,
				      const int Number_of_Ghost_Cells,
				      int node_init_procedure,
				      int StretchFcnI, double  BetaI, double  TauI,
				      int StretchFcnJ, double  BetaJ, double  TauJ,
				      int NbndType, int  EbndType, int  SbndType, int WbndType,
				      int OrthN, int  OrthE, int  OrthS, int  OrthW) {

  int i;
  Spline2D Bnd_Spline_North, Bnd_Spline_South, Bnd_Spline_East, Bnd_Spline_West;

  //set north, south, east and west boundary types for each boundary point
  for ( i=0; i<block.BndNorthSpline.np; ++i )
    block.BndNorthSpline.bc[i]= NbndType;

  for ( i=0; i<block.BndSouthSpline.np; ++i )
    block.BndSouthSpline.bc[i]= SbndType;
  
  for ( i=0; i<block.BndEastSpline.np; ++i )
    block.BndEastSpline.bc[i]= EbndType;
  
  for ( i=0; i<block.BndWestSpline.np; ++i )
    block.BndWestSpline.bc[i]= WbndType;

  //copy boundary splines to temp variables
  Copy_Spline(Bnd_Spline_North, block.BndNorthSpline);
  Copy_Spline(Bnd_Spline_South, block.BndSouthSpline);
  Copy_Spline(Bnd_Spline_East, block.BndEastSpline);
  Copy_Spline(Bnd_Spline_West, block.BndWestSpline);
  
  /* Create the 2D quadrilateral grid block. */

  Create_Quad_Block(block,
                    Bnd_Spline_North,
                    Bnd_Spline_South,
                    Bnd_Spline_East,
                    Bnd_Spline_West,
                    Number_of_Cells_Idir,
		    Number_of_Cells_Jdir,
		    Number_of_Ghost_Cells,
                    node_init_procedure,
                    StretchFcnI, 
                    BetaI, 
                    TauI,
		    StretchFcnJ, 
                    BetaJ, 
                    TauJ,
                    OrthN, 
                    OrthE, 
                    OrthS, 
                    OrthW);

  /* Smooth the 2D quadrilateral grid block. */

  Smooth_Quad_Block(block, n);

  /* Deallocate the memory for the boundary splines. */
  
  Bnd_Spline_North.deallocate();
  Bnd_Spline_South.deallocate();
  Bnd_Spline_East.deallocate();
  Bnd_Spline_West.deallocate();

}

/**********************************************************************
 * Routine: genMeshH_3x2                                              *
 *                                                                    *
 * -> generate a mesh for the cross-section of the blade at any %span *
 * -> H GRID (multi-block)                                            *
 * -> PARAMETERS MUST BE SUPPLIED FOR EVERY BLOCK                     *
 *                                                                    *
 * __________________________________________________________________ *
 * |                /                              \                | *
 * |               /                                \               | *
 * | Block (0,1)  |           Block (1,1)            |  Block (2,1) | *
 * |               \                                /               | *
 * |                \                              /                | *
 * |----------------- <<<<<<<<<<<<BLADE>>>>>>>>>>> -----------------| *
 * |                /                              \                | *
 * |               /                                \               | *
 * | Block (0,0)  |           Block (1,0)            |  Block (2,0) | *
 * |               \                                /               | *
 * |                \                              /                | *
 * |________________________________________________________________| *
 *                                                                    *
 **********************************************************************/
inline Grid2D_Quad_Block** NASARotor37::genMeshH_3x2(

			Grid2D_Quad_Block **mesh,
                        double pspan, double zMin, double zMax, int n, 

                        //Parameters for block (0,0) - (BOTTOM LEFT CORNER)
			const int B00_Number_of_Cells_Idir, const int B00_Number_of_Cells_Jdir, const int B00_Number_of_Ghost_Cells,
			int B00_StretchFcnI, double B00_BetaI, double B00_TauI,
			int B00_StretchFcnJ, double B00_BetaJ, double B00_TauJ,
			int B00_NbndType, int B00_EbndType, int B00_SbndType, int B00_WbndType,
			int B00_OrthN, int B00_OrthE, int B00_OrthS, int B00_OrthW,
			
                        //Parameters for block (1,0)
			const int B10_Number_of_Cells_Idir, const int B10_Number_of_Cells_Jdir, const int B10_Number_of_Ghost_Cells,
			int B10_StretchFcnI, double B10_BetaI, double B10_TauI,
			int B10_StretchFcnJ, double B10_BetaJ, double B10_TauJ,
			int B10_NbndType, int B10_EbndType, int B10_SbndType, int B10_WbndType,
			int B10_OrthN, int B10_OrthE, int B10_OrthS, int B10_OrthW,

                        //Parameters for block (2,0)
			const int B20_Number_of_Cells_Idir, const int B20_Number_of_Cells_Jdir, const int B20_Number_of_Ghost_Cells,
			int B20_StretchFcnI, double B20_BetaI, double B20_TauI,
			int B20_StretchFcnJ, double B20_BetaJ, double B20_TauJ,
			int B20_NbndType, int B20_EbndType, int B20_SbndType, int B20_WbndType,
			int B20_OrthN, int B20_OrthE, int B20_OrthS, int B20_OrthW,

                        //Parameters for block (0,1)
			const int B01_Number_of_Cells_Idir, const int B01_Number_of_Cells_Jdir, const int B01_Number_of_Ghost_Cells,
			int B01_StretchFcnI, double B01_BetaI, double B01_TauI,
			int B01_StretchFcnJ, double B01_BetaJ, double B01_TauJ,
			int B01_NbndType, int B01_EbndType, int B01_SbndType, int B01_WbndType,
			int B01_OrthN, int B01_OrthE, int B01_OrthS, int B01_OrthW,
			
                        //Parameters for block (1,1)
			const int B11_Number_of_Cells_Idir, const int B11_Number_of_Cells_Jdir, const int B11_Number_of_Ghost_Cells,
			int B11_StretchFcnI, double B11_BetaI, double B11_TauI,
			int B11_StretchFcnJ, double B11_BetaJ, double B11_TauJ,
			int B11_NbndType, int B11_EbndType, int B11_SbndType, int B11_WbndType,
			int B11_OrthN, int B11_OrthE, int B11_OrthS, int B11_OrthW,

                        //Parameters for block (2,1)
			const int B21_Number_of_Cells_Idir, const int B21_Number_of_Cells_Jdir, const int B21_Number_of_Ghost_Cells,
			int B21_StretchFcnI, double B21_BetaI, double B21_TauI,
			int B21_StretchFcnJ, double B21_BetaJ, double B21_TauJ,
			int B21_NbndType, int B21_EbndType, int B21_SbndType, int B21_WbndType,
			int B21_OrthN, int B21_OrthE, int B21_OrthS, int B21_OrthW) {

  Spline2D upperB, lowerB, upperMiddleB, lowerMiddleB, camberTrail, camberLead, camberBlade, cs;
  int i, iTrail, iLead, pos, mIndex;
  Vector2D swap, leadV, trailV, x_temp;
  double zlU, zlL, zrU, zrL, mLead, mTrail, mLeadc, mTrailc, 
         mTop, mBot, m1, m2, z, dz, A, m, zLead, zTrail, dm, swapm;

  //*******************************************************************************
  //*******************************************************************************
  //*******************************  INITIALIZATION *******************************
  //*******************************************************************************
  //*******************************************************************************

  cout << "\n  Initializing Mesh Variables...";
  cout.flush();

  //Allocate memory for 3x2 multiblock grid.
  mesh = Allocate_Multi_Block_Grid(mesh, 3, 2);

  //find the index for the leading and trailing edge
  findLT(pspan, iTrail, iLead);

  //get cross section spline
  cs=getBladeCS(pspan);

  //re-organize pts so they ascend in the clockwise direction
  for ( i=0; i<R37GEOM_NUMPTS; ++i ) {
    swap=cs.Xp[2*R37GEOM_NUMPTS-2-i];
    cs.Xp[2*R37GEOM_NUMPTS-2-i]=cs.Xp[i];
    cs.Xp[i]=swap;
  }

  //change iTrail and iLead accordingly
  iTrail=cs.np-1-iTrail;
  iLead=cs.np-1-iLead;

  //save leading and trailing edge vectors
  leadV=cs.Xp[iLead];
  trailV=cs.Xp[iTrail];

  //get camber line splines - 40 pts
  camberTrail=getCamberTrail(pspan, 40, zMax); 
  camberLead=getCamberLead(pspan, 40, zMin);
  camberBlade=getCamberBlade(pspan, 40);
  
  //get upper, lower, upperMiddle, and lowerMiddle boundaries
  upperB=getCamberAndShift(pspan, 40, zMin, zMax, M_PI/num_blades);
  lowerB=getCamberAndShift(pspan, 40, zMin, zMax, -M_PI/num_blades);
  upperMiddleB=getCamberAndShift(pspan, 40, zMin, zMax, M_PI/(2*num_blades));
  lowerMiddleB=getCamberAndShift(pspan, 40, zMin, zMax, -M_PI/(2*num_blades));

  // Locate boundary corners to ensure that boundaries meeting at the leading
  // and trailing edges are separated by approximately 45deg
   
  //calculate slope of leading and trailing camber lines
  mLeadc=(camberLead.Xp[camberLead.np-1].y-camberLead.Xp[0].y)/
    (camberLead.Xp[camberLead.np-1].x-camberLead.Xp[0].x);

  mTrailc=(camberTrail.Xp[camberTrail.np-1].y-camberTrail.Xp[0].y)/
    (camberTrail.Xp[camberTrail.np-1].x-camberTrail.Xp[0].x);

  
  //calculate slope of tangent to cross-section at the leading and 
  //trailing edge using a central-difference
  z=cs.Xp[iLead].x;
  dz=(cs.Xp[iLead].x-cs.Xp[iLead+1].x)/10.0;
  mLead=(-getY(z+2*dz, cs)[0].y+8*getY(z+dz, cs)[0].y-8*getY(z-dz, cs)[0].y+
        getY(z-2*dz, cs)[0].y)/(12*dz);

  z=cs.Xp[iTrail].x;
  dz=(cs.Xp[iTrail].x-cs.Xp[iTrail-1].x)/10.0;
  mTrail=(-getY(z+2*dz, cs)[1].y+8*getY(z+dz, cs)[1].y-8*getY(z-dz, cs)[1].y+
         getY(z-2*dz, cs)[1].y)/(12*dz);

  //get slopes of lines approx 45deg to camber line at the LEADING EDGE
  A=(mLead+mLeadc)/(1-mLead*mLeadc);
  mTop=(-1+sqrt(1+A*A))/A;
  mBot=(-1-sqrt(1+A*A))/A;

  //find point on upperMiddle spline - leading edge
  //do search
  dm=1e15;
  mIndex=-1;
  for ( i=0; i<upperMiddleB.np; ++i ) {
    if (fabs(upperMiddleB.Xp[i].x-leadV.x) > TOLER) {
       m=(upperMiddleB.Xp[i].y-leadV.y)/(upperMiddleB.Xp[i].x-leadV.x);
      
       if( fabs(m-mTop)<dm) {
         dm=fabs(m-mTop);
         mIndex=i;
       } /* endif */
    } /* endif */
  }
  assert(mIndex!=-1);
  zlU=upperMiddleB.Xp[mIndex].x;  

  //find point on lowerMiddle spline - leading edge
  //do search
  dm=1e15;
  mIndex=-1;
  for ( i=0; i<lowerMiddleB.np; ++i ) {
    if (fabs(lowerMiddleB.Xp[i].x-leadV.x) > TOLER) {
       m=(lowerMiddleB.Xp[i].y-leadV.y)/(lowerMiddleB.Xp[i].x-leadV.x);
      
       if( fabs(m-mBot)<dm) {
         dm=fabs(m-mBot);
         mIndex=i;
       } /* endif */
    } /* endif */
  }
  assert(mIndex!=-1);
  zlL=lowerMiddleB.Xp[mIndex].x;  

  //get slopes of lines approx 45deg to camber line at the TRAILING EDGE
  A=(mTrail+mTrailc)/(1-mTrail*mTrailc);
  mBot=(-1+sqrt(1+A*A))/A;
  mTop=(-1-sqrt(1+A*A))/A;
  
  //find point on upperMiddle spline - trailing edge
  //do search
  dm=1e15;
  mIndex=-1;
  for ( i=0; i<upperMiddleB.np; ++i ) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
       m=(upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      
       if( fabs(m-mTop)<dm) {
         dm=fabs(m-mTop);
         mIndex=i;
       } /* endif */
    } /* endif */
  }
  assert(mIndex!=-1);
  zrU=upperMiddleB.Xp[mIndex].x;  

  dm=1e15;
  mIndex=-1;
  for ( i=0; i<upperMiddleB.np; ++i ) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
       m=(upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      
       if( fabs(m-mBot)<dm) {
         dm=fabs(m-mBot);
         mIndex=i;
       } /* endif */
    } /* endif */
  }
  assert(mIndex!=-1);
  
  //determine which to use
  if(upperMiddleB.Xp[mIndex].x>zrU) {
    zrU=upperMiddleB.Xp[mIndex].x;  
    swapm=mBot;
    mBot=mTop;
    mTop=swapm;
  };  

  //find point on lowerMiddle spline - trailing edge
  //do search
  dm=1e15;
  mIndex=-1;
  for ( i=0; i<lowerMiddleB.np; ++i ) {
    if (fabs(lowerMiddleB.Xp[i].x-trailV.x) > TOLER) {
       m=(lowerMiddleB.Xp[i].y-trailV.y)/(lowerMiddleB.Xp[i].x-trailV.x);
      
       if( fabs(m-mBot)<dm) {
         dm=fabs(m-mBot);
         mIndex=i;
       } /* endif */
    } /* endif */
  }
  assert(mIndex!=-1);
  zrL=lowerMiddleB.Xp[mIndex].x;  

  //*******************************************************************************
  //*******************************************************************************
  //****************************  GENERATE BOUNDARIES  ****************************
  //*******************************************************************************
  //*******************************************************************************  

  cout << "\n  Generating Boundaries...";
  cout.flush();

  // **********************
  // **** BLOCK (1,1) *****
  // **********************

  // ** NORTH SPLINE **
  //allocate space
  mesh[1][1].BndNorthSpline.allocate(70);

  //set type
  mesh[1][1].BndNorthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=(zlL+zlU)/2;
  zTrail=(zrL+zrU)/2;
  dz=(zTrail-zLead)/69;
  for ( i=0; i<70; ++i ) {
    mesh[1][1].BndNorthSpline.Xp[i]=getY(zLead+i*dz, upperB)[0];
  }

  //set point types
  mesh[1][1].BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][1].BndNorthSpline.tp[mesh[1][1].BndNorthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i< mesh[1][1].BndNorthSpline.np-1; ++i )
    mesh[1][1].BndNorthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh[1][1].BndNorthSpline.pathlength();

  // ** SOUTH SPLINE **
  //allocate space
  mesh[1][1].BndSouthSpline.allocate(iTrail-iLead+1);
  
  //set type
  mesh[1][1].BndSouthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  for ( i=iLead; i<=iTrail; ++i )
    mesh[1][1].BndSouthSpline.Xp[i-iLead]=cs.Xp[i];

  //set point types
  mesh[1][1].BndSouthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][1].BndSouthSpline.tp[mesh[1][1].BndSouthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh[1][1].BndSouthSpline.np-1; ++i )
    mesh[1][1].BndSouthSpline.tp[i]=SPLINE2D_POINT_NORMAL;
  
  //calc pathlength
  mesh[1][1].BndSouthSpline.pathlength();

  // ** EAST SPLINE **
  //allocate space
  mesh[1][1].BndEastSpline.allocate(5);
  
  //set type
  mesh[1][1].BndEastSpline.settype(SPLINE2D_CUBIC);

  //set points
  mesh[1][1].BndEastSpline.Xp[0]=cs.Xp[iTrail];
  x_temp = getY(zrU, upperMiddleB)[0];
  mesh[1][1].BndEastSpline.Xp[1]=cs.Xp[iTrail]+HALF*(x_temp-cs.Xp[iTrail]);
  mesh[1][1].BndEastSpline.Xp[2]=x_temp;
  mesh[1][1].BndEastSpline.Xp[3]=x_temp+
     HALF*(mesh[1][1].BndNorthSpline.Xp[mesh[1][1].BndNorthSpline.np-1]-x_temp);
  mesh[1][1].BndEastSpline.Xp[4]=mesh[1][1].BndNorthSpline.Xp[mesh[1][1].BndNorthSpline.np-1];
  
  //set point types
  mesh[1][1].BndEastSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][1].BndEastSpline.tp[1]=SPLINE2D_POINT_NORMAL;
  mesh[1][1].BndEastSpline.tp[2]=SPLINE2D_POINT_NORMAL;
  mesh[1][1].BndEastSpline.tp[3]=SPLINE2D_POINT_NORMAL;
  mesh[1][1].BndEastSpline.tp[4]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh[1][1].BndEastSpline.pathlength();

  // ** WEST SPLINE **
  //allocate space
  mesh[1][1].BndWestSpline.allocate(6);
  
  //set type
  mesh[1][1].BndWestSpline.settype(SPLINE2D_CUBIC);

  //set points
  mesh[1][1].BndWestSpline.Xp[0]=cs.Xp[iLead];
  x_temp = getY(zlU, upperMiddleB)[0];
  mesh[1][1].BndWestSpline.Xp[1]=cs.Xp[iLead]+0.3333*(x_temp-cs.Xp[iLead]);
  mesh[1][1].BndWestSpline.Xp[2]=HALF*(cs.Xp[iLead]+HALF*(x_temp-cs.Xp[iLead])+
                                       x_temp+0.20*(mesh[1][1].BndNorthSpline.Xp[0]-x_temp));
  mesh[1][1].BndWestSpline.Xp[3]=HALF*(cs.Xp[iLead]+0.6666*(x_temp-cs.Xp[iLead])+
                                       x_temp+0.3333*(mesh[1][1].BndNorthSpline.Xp[0]-x_temp));
  mesh[1][1].BndWestSpline.Xp[4]=x_temp+HALF*(mesh[1][1].BndNorthSpline.Xp[0]-x_temp);
  mesh[1][1].BndWestSpline.Xp[5]=mesh[1][1].BndNorthSpline.Xp[0];
  
  //set point types
  mesh[1][1].BndWestSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][1].BndWestSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][1].BndWestSpline.tp[2]=SPLINE2D_POINT_NORMAL;
  mesh[1][1].BndWestSpline.tp[3]=SPLINE2D_POINT_NORMAL;
  mesh[1][1].BndWestSpline.tp[4]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][1].BndWestSpline.tp[5]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh[1][1].BndWestSpline.pathlength();

  // **********************
  // **** BLOCK (1,0) *****
  // **********************
 
  // ** NORTH SPLINE **
  //allocate space
  mesh[1][0].BndNorthSpline.allocate(cs.np-iTrail+iLead);
  
  //set type
  mesh[1][0].BndNorthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  for ( i=iLead; i>=0; --i )
    mesh[1][0].BndNorthSpline.Xp[iLead-i]=cs.Xp[i];
  for ( i=cs.np-2; i>=iTrail; --i )
    mesh[1][0].BndNorthSpline.Xp[iLead+1+(cs.np-2)-i]=cs.Xp[i];

  //set point types
  mesh[1][0].BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][0].BndNorthSpline.tp[mesh[1][0].BndNorthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh[1][0].BndNorthSpline.np-1; ++i )
    mesh[1][0].BndNorthSpline.tp[i]=SPLINE2D_POINT_NORMAL;
  
  //calc pathlength
  mesh[1][0].BndNorthSpline.pathlength();  

  // ** SOUTH SPLINE **
  //allocate space
  mesh[1][0].BndSouthSpline.allocate(70);
  
  //set type
  mesh[1][0].BndSouthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=(zlL+zlU)/2;;
  zTrail=(zrL+zrU)/2;;
  dz=(zTrail-zLead)/69;
  for ( i=0; i<70; ++i ) {
    mesh[1][0].BndSouthSpline.Xp[i]=getY(zLead+i*dz, lowerB)[0];
  }

  //set point types
  mesh[1][0].BndSouthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][0].BndSouthSpline.tp[mesh[1][0].BndSouthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh[1][0].BndSouthSpline.np-1; ++i )
    mesh[1][0].BndSouthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh[1][0].BndSouthSpline.pathlength();

  // ** WEST SPLINE **
  //allocate space
  mesh[1][0].BndWestSpline.allocate(5);
  
  //set type
  mesh[1][0].BndWestSpline.settype(SPLINE2D_CUBIC);

  //set points
  mesh[1][0].BndWestSpline.Xp[0]=mesh[1][0].BndSouthSpline.Xp[0];
  x_temp = getY(zlL, lowerMiddleB)[0];
  mesh[1][0].BndWestSpline.Xp[1]=mesh[1][0].BndSouthSpline.Xp[0]+
                                 HALF*(x_temp-mesh[1][0].BndSouthSpline.Xp[0]);
  mesh[1][0].BndWestSpline.Xp[2]=x_temp;
  mesh[1][0].BndWestSpline.Xp[3]=x_temp+HALF*(cs.Xp[iLead]-x_temp);
  mesh[1][0].BndWestSpline.Xp[4]=cs.Xp[iLead];
  
  //set point types
  mesh[1][0].BndWestSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][0].BndWestSpline.tp[1]=SPLINE2D_POINT_NORMAL;
  mesh[1][0].BndWestSpline.tp[2]=SPLINE2D_POINT_NORMAL;
  mesh[1][0].BndWestSpline.tp[3]=SPLINE2D_POINT_NORMAL;
  mesh[1][0].BndWestSpline.tp[4]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh[1][0].BndWestSpline.pathlength();

  // ** EAST SPLINE **
  //allocate space
  mesh[1][0].BndEastSpline.allocate(6);
  
  //set type
  mesh[1][0].BndEastSpline.settype(SPLINE2D_CUBIC);

  //set points
  mesh[1][0].BndEastSpline.Xp[0]=mesh[1][0].BndSouthSpline.Xp[mesh[1][0].BndSouthSpline.np-1];
  x_temp = getY(zrL, lowerMiddleB)[0];
  mesh[1][0].BndEastSpline.Xp[1]=x_temp + 
     HALF*(mesh[1][0].BndSouthSpline.Xp[mesh[1][0].BndSouthSpline.np-1]-x_temp);
  mesh[1][0].BndEastSpline.Xp[2]=HALF*(cs.Xp[iTrail]+0.6666*(x_temp-cs.Xp[iTrail])+
     x_temp+0.3333*(mesh[1][0].BndSouthSpline.Xp[mesh[1][0].BndSouthSpline.np-1]-x_temp));
  mesh[1][0].BndEastSpline.Xp[3]=0.95*(cs.Xp[iTrail]+HALF*(x_temp-cs.Xp[iTrail]))+
     0.05*(x_temp+0.80*(mesh[1][0].BndSouthSpline.Xp[mesh[1][0].BndSouthSpline.np-1]-x_temp));
  mesh[1][0].BndEastSpline.Xp[4]=cs.Xp[iTrail]+0.3333*(x_temp-cs.Xp[iTrail]);
  mesh[1][0].BndEastSpline.Xp[5]=cs.Xp[iTrail];
   
  //set point types
  mesh[1][0].BndEastSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][0].BndEastSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][0].BndEastSpline.tp[2]=SPLINE2D_POINT_NORMAL;
  mesh[1][0].BndEastSpline.tp[3]=SPLINE2D_POINT_NORMAL;
  mesh[1][0].BndEastSpline.tp[4]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[1][0].BndEastSpline.tp[5]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh[1][0].BndEastSpline.pathlength();

  // **********************
  // **** BLOCK (0,0) *****
  // **********************
  
  // ** NORTH SPLINE **
  //allocate space
  mesh[0][0].BndNorthSpline.allocate(2);
  
  //set type
  mesh[0][0].BndNorthSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh[0][0].BndNorthSpline.Xp[0]=camberLead.Xp[0];
  mesh[0][0].BndNorthSpline.Xp[1]=cs.Xp[iLead];
  
  //set point types
  mesh[0][0].BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[0][0].BndNorthSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh[0][0].BndNorthSpline.pathlength();

  // ** SOUTH SPLINE ** 
  //allocate space
  mesh[0][0].BndSouthSpline.allocate(30);
  
  //set type
  mesh[0][0].BndSouthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=lowerB.Xp[0].x;
  zTrail=mesh[1][0].BndSouthSpline.Xp[0].x;
  dz=(zTrail-zLead)/29;
  for ( i=0; i<30; ++i ) {
    mesh[0][0].BndSouthSpline.Xp[i]=getY(zLead+i*dz, lowerB)[0];
  }
  mesh[0][0].BndSouthSpline.Xp[0]=lowerB.Xp[0];

  //set point types
  mesh[0][0].BndSouthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[0][0].BndSouthSpline.tp[mesh[0][0].BndSouthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh[0][0].BndSouthSpline.np-1; ++i )
    mesh[0][0].BndSouthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh[0][0].BndSouthSpline.pathlength();

  // ** EAST SPLINE **
  Copy_Spline(mesh[0][0].BndEastSpline, mesh[1][0].BndWestSpline);  

  // ** WEST SPLINE **
  //allocate space
  mesh[0][0].BndWestSpline.allocate(2);
  
  //set type
  mesh[0][0].BndWestSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh[0][0].BndWestSpline.Xp[0]=mesh[0][0].BndSouthSpline.Xp[0];
  mesh[0][0].BndWestSpline.Xp[1]=mesh[0][0].BndNorthSpline.Xp[0];
  
  //set point types
  mesh[0][0].BndWestSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[0][0].BndWestSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh[0][0].BndWestSpline.pathlength();

  // **********************
  // **** BLOCK (2,0) *****
  // **********************
  
  // ** NORTH SPLINE **
  //allocate space
  mesh[2][0].BndNorthSpline.allocate(2);
  
  //set type
  mesh[2][0].BndNorthSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh[2][0].BndNorthSpline.Xp[0]=cs.Xp[iTrail];
  mesh[2][0].BndNorthSpline.Xp[1]=camberTrail.Xp[camberTrail.np-1];
  
  //set point types
  mesh[2][0].BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[2][0].BndNorthSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh[2][0].BndNorthSpline.pathlength();

  // ** SOUTH SPLINE ** 
  //allocate space
  mesh[2][0].BndSouthSpline.allocate(30);
  
  //set type
  mesh[2][0].BndSouthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=(zrL+zrU)/2;;
  zTrail=lowerB.Xp[lowerB.np-1].x;
  dz=(zTrail-zLead)/29;

  for ( i=0; i<30; ++i ) {
    mesh[2][0].BndSouthSpline.Xp[i].x=zLead+i*dz;
    //prevent round-off error
    if(mesh[2][0].BndSouthSpline.Xp[i].x>zTrail) mesh[2][0].BndSouthSpline.Xp[i].x=zTrail;
    mesh[2][0].BndSouthSpline.Xp[i].y=getY(mesh[2][0].BndSouthSpline.Xp[i].x, lowerB)[0].y;
  }
  mesh[2][0].BndSouthSpline.Xp[mesh[2][0].BndSouthSpline.np-1]=lowerB.Xp[lowerB.np-1];

  //set point types
  mesh[2][0].BndSouthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[2][0].BndSouthSpline.tp[mesh[2][0].BndSouthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh[2][0].BndSouthSpline.np-1; ++i )
    mesh[2][0].BndSouthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh[2][0].BndSouthSpline.pathlength();

  // ** WEST SPLINE **
  Copy_Spline(mesh[2][0].BndWestSpline, mesh[1][0].BndEastSpline); 

  // ** EAST SPLINE **
  //allocate space
  mesh[2][0].BndEastSpline.allocate(2);
  
  //set type
  mesh[2][0].BndEastSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh[2][0].BndEastSpline.Xp[0]=mesh[2][0].BndSouthSpline.Xp[mesh[2][0].BndSouthSpline.np-1];
  mesh[2][0].BndEastSpline.Xp[1]=mesh[2][0].BndNorthSpline.Xp[mesh[2][0].BndNorthSpline.np-1];
  
  //set point types
  mesh[2][0].BndEastSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[2][0].BndEastSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh[2][0].BndEastSpline.pathlength();

  // **********************
  // **** BLOCK (0,1) *****
  // **********************
  
  // ** NORTH SPLINE **
  //allocate space
  mesh[0][1].BndNorthSpline.allocate(30);
  
  //set type
  mesh[0][1].BndNorthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=upperB.Xp[0].x;
  zTrail=(zlL+zlU)/2;
  dz=(zTrail-zLead)/29;
  for ( i=0; i<30; ++i ) {
    mesh[0][1].BndNorthSpline.Xp[i]=getY(zLead+i*dz, upperB)[0];
  }
  mesh[0][1].BndNorthSpline.Xp[0]=upperB.Xp[0];

  //set point types
  mesh[0][1].BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[0][1].BndNorthSpline.tp[mesh[0][1].BndNorthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh[0][1].BndNorthSpline.np-1; ++i )
    mesh[0][1].BndNorthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh[0][1].BndNorthSpline.pathlength();

  // ** SOUTH SPLINE **
  Copy_Spline(mesh[0][1].BndSouthSpline, mesh[0][0].BndNorthSpline);

  // ** EAST SPLINE **
  Copy_Spline(mesh[0][1].BndEastSpline, mesh[1][1].BndWestSpline);

  // ** WEST SPLINE **
  //allocate space
  mesh[0][1].BndWestSpline.allocate(2);
  
  //set type
  mesh[0][1].BndWestSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh[0][1].BndWestSpline.Xp[0]=mesh[0][1].BndSouthSpline.Xp[0];
  mesh[0][1].BndWestSpline.Xp[1]=mesh[0][1].BndNorthSpline.Xp[0];
  
  //set point types
  mesh[0][1].BndWestSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[0][1].BndWestSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh[0][1].BndWestSpline.pathlength();

  // **********************
  // **** BLOCK (2,1) *****
  // **********************
  
  // ** NORTH SPLINE **
  //allocate space
  mesh[2][1].BndNorthSpline.allocate(30);
  
  //set type
  mesh[2][1].BndNorthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=(zrL+zrU)/2;
  zTrail=upperB.Xp[upperB.np-1].x;
  dz=(zTrail-zLead)/29;

  for ( i=0; i<30; ++i ) {
    mesh[2][1].BndNorthSpline.Xp[i].x=zLead+i*dz;
    //prevent round-off error
    if(mesh[2][1].BndNorthSpline.Xp[i].x>zTrail)
      mesh[2][1].BndNorthSpline.Xp[i].x=zTrail;
    mesh[2][1].BndNorthSpline.Xp[i].y=getY(mesh[2][1].BndNorthSpline.Xp[i].x, upperB)[0].y;
  }
  mesh[2][1].BndNorthSpline.Xp[mesh[2][1].BndNorthSpline.np-1]=upperB.Xp[upperB.np-1];

  //set point types
  mesh[2][1].BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[2][1].BndNorthSpline.tp[mesh[2][1].BndNorthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh[2][1].BndNorthSpline.np-1; ++i )
    mesh[2][1].BndNorthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh[2][1].BndNorthSpline.pathlength();

  // ** SOUTH SPLINE **
  Copy_Spline(mesh[2][1].BndSouthSpline, mesh[2][0].BndNorthSpline);

  // ** WEST SPLINE **
  Copy_Spline(mesh[2][1].BndWestSpline, mesh[1][1].BndEastSpline);

  // ** EAST SPLINE **
  //allocate space
  mesh[2][1].BndEastSpline.allocate(2);
  
  //set type
  mesh[2][1].BndEastSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh[2][1].BndEastSpline.Xp[0]=mesh[2][1].BndSouthSpline.Xp[mesh[2][1].BndSouthSpline.np-1];
  mesh[2][1].BndEastSpline.Xp[1]=mesh[2][1].BndNorthSpline.Xp[mesh[2][1].BndNorthSpline.np-1];
  
  //set point types
  mesh[2][1].BndEastSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh[2][1].BndEastSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh[2][1].BndEastSpline.pathlength();

  //*******************************************************************************
  //*******************************************************************************
  //*******************************  GENERATE MESH ********************************
  //*******************************************************************************
  //*******************************************************************************

  cout << "\n  Generating Mesh - Block (0,0)...";
  cout.flush();

  genMeshBlock(mesh[0][0], n, B00_Number_of_Cells_Idir, B00_Number_of_Cells_Jdir, B00_Number_of_Ghost_Cells,
	       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
	       B00_StretchFcnI, B00_BetaI, B00_TauI,
	       B00_StretchFcnJ, B00_BetaJ, B00_TauJ,
	       B00_NbndType, B00_EbndType, B00_SbndType, B00_WbndType,
	       B00_OrthN, B00_OrthE, B00_OrthS, B00_OrthW);

  cout << "\n  Generating Mesh - Block (1,0)...";
  cout.flush();

  genMeshBlock(mesh[1][0], n, B10_Number_of_Cells_Idir, B10_Number_of_Cells_Jdir, B10_Number_of_Ghost_Cells,
	       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY,
	       B10_StretchFcnI, B10_BetaI, B10_TauI,
	       B10_StretchFcnJ, B10_BetaJ, B10_TauJ,
	       B10_NbndType, B10_EbndType, B10_SbndType, B10_WbndType,
	       B10_OrthN, B10_OrthE, B10_OrthS, B10_OrthW);

  cout << "\n  Generating Mesh - Block (2,0)...";
  cout.flush();

  genMeshBlock(mesh[2][0], n, B20_Number_of_Cells_Idir, B20_Number_of_Cells_Jdir, B20_Number_of_Ghost_Cells,
	       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
	       B20_StretchFcnI, B20_BetaI, B20_TauI,
	       B20_StretchFcnJ, B20_BetaJ, B20_TauJ,
	       B20_NbndType, B20_EbndType, B20_SbndType, B20_WbndType,
	       B20_OrthN, B20_OrthE, B20_OrthS, B20_OrthW);

  cout << "\n  Generating Mesh - Block (0,1)...";
  cout.flush();

  genMeshBlock(mesh[0][1], n, B01_Number_of_Cells_Idir, B01_Number_of_Cells_Jdir, B01_Number_of_Ghost_Cells,
               GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
	       B01_StretchFcnI, B01_BetaI, B01_TauI,
	       B01_StretchFcnJ, B01_BetaJ, B01_TauJ,
	       B01_NbndType, B01_EbndType, B01_SbndType, B01_WbndType,
	       B01_OrthN, B01_OrthE, B01_OrthS, B01_OrthW);       

  cout << "\n  Generating Mesh - Block (1,1)...";
  cout.flush();

  genMeshBlock(mesh[1][1], n, B11_Number_of_Cells_Idir, B11_Number_of_Cells_Jdir, B11_Number_of_Ghost_Cells,
	       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY,
	       B11_StretchFcnI, B11_BetaI, B11_TauI,
	       B11_StretchFcnJ, B11_BetaJ, B11_TauJ,
	       B11_NbndType, B11_EbndType, B11_SbndType, B11_WbndType,
	       B11_OrthN, B11_OrthE, B11_OrthS, B11_OrthW);

  cout << "\n  Generating Mesh - Block (2,1)...";
  cout.flush();

  genMeshBlock(mesh[2][1], n, B21_Number_of_Cells_Idir, B21_Number_of_Cells_Jdir, B21_Number_of_Ghost_Cells,
	       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
	       B21_StretchFcnI, B21_BetaI, B21_TauI,
	       B21_StretchFcnJ, B21_BetaJ, B21_TauJ,
	       B21_NbndType, B21_EbndType, B21_SbndType, B21_WbndType,
	       B21_OrthN, B21_OrthE, B21_OrthS, B21_OrthW);

  /* Return multiblock mesh. */
  return(mesh);

}

/**********************************************************************
 * Routine: genMeshH_3x2_AUTO                                         *
 *                                                                    *
 * Calls genMeshH_3x2 with preset parameters.                         *
 *                                                                    *
 **********************************************************************/
inline Grid2D_Quad_Block** NASARotor37::genMeshH_3x2_AUTO(Grid2D_Quad_Block **mesh,
							  double pspan, double zMin, double zMax,  
							  const int Number_of_Cells_Idir,
							  const int Number_of_Cells_Jdir,
							  const int Number_of_Ghost_Cells,
							  int n) {

  double stretch_i, stretch_j;

  cout << "\n  NASA Rotor 37 automatic 3x2 multi-block mesh generation...";

  //'hard-wire' stretching parameters
  if (Number_of_Cells_Idir/2 < 10) {
     stretch_i=1.05;
  } else if (Number_of_Cells_Idir/2 < 25) {
     stretch_i=1.01;
  } else if (Number_of_Cells_Idir/2 < 50) {
     stretch_i=1.005;
  } else if (Number_of_Cells_Idir/2 < 100) {
     stretch_i=1.0025;
  } else {
     stretch_i=1.001;
  } /* endif */
  if (Number_of_Cells_Jdir < 10) {
     stretch_j=1.05;
  } else if (Number_of_Cells_Jdir < 25) {
     stretch_j=1.01;
  } else if (Number_of_Cells_Jdir < 50) {
     stretch_j=1.005;
  } else if (Number_of_Cells_Jdir < 100) {
     stretch_j=1.0025;
  } else {
     stretch_j=1.001;
  } /* endif */

  //generate mesh
  mesh = genMeshH_3x2(mesh,

           pspan, zMin, zMax, n,
	   
           Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	   STRETCHING_FCN_MAX_CLUSTERING, stretch_i, 0,
	   STRETCHING_FCN_MAX_CLUSTERING, stretch_j, 0,
	   BC_NONE, BC_NONE, BC_NONE, BC_FIXED,
	   ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, NOT_ORTHOGONAL, 
	   
           Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	   STRETCHING_FCN_MINMAX_CLUSTERING, ONE+(stretch_i-ONE)/ONE, 0,
	   STRETCHING_FCN_MAX_CLUSTERING, stretch_j, 0,
	   BC_REFLECTION, BC_NONE, BC_NONE, BC_NONE,
	   ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, 
	   
           Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	   STRETCHING_FCN_MIN_CLUSTERING, stretch_i, 0,
	   STRETCHING_FCN_MAX_CLUSTERING, stretch_j, 0,
	   BC_NONE, BC_CHARACTERISTIC, BC_NONE, BC_NONE,
	   NOT_ORTHOGONAL, NOT_ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, 
	   
           Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	   STRETCHING_FCN_MAX_CLUSTERING, stretch_i, 0,
	   STRETCHING_FCN_MIN_CLUSTERING, stretch_j, 0,
	   BC_NONE, BC_NONE, BC_NONE, BC_FIXED,
	   ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, NOT_ORTHOGONAL, 
	   
           Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	   STRETCHING_FCN_MINMAX_CLUSTERING, ONE+(stretch_i-ONE)/ONE, 0,
	   STRETCHING_FCN_MIN_CLUSTERING, stretch_j, 0,
	   BC_NONE, BC_NONE, BC_REFLECTION, BC_NONE,
	   ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, 
	   
           Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	   STRETCHING_FCN_MIN_CLUSTERING, stretch_i, 0,
	   STRETCHING_FCN_MIN_CLUSTERING, stretch_j, 0,
	   BC_NONE, BC_CHARACTERISTIC, BC_NONE, BC_NONE,
	   ORTHOGONAL, NOT_ORTHOGONAL, ORTHOGONAL, ORTHOGONAL);   

  return(mesh);

}

/**********************************************************************
 * Routine: outputTP_Mesh3D                                           *
 *                                                                    *
 * Creates a structured 3D mesh about a complete blade and writes the *
 * mesh to a TECPLOT output file.                                     *
 *                                                                    *
 **********************************************************************/
inline void NASARotor37::outputTP_Mesh3D(const int Number_of_Cells_Idir,
					 const int Number_of_Cells_Jdir,
					 int kMax,
					 double zMin, double zMax,
					 char fnameTop[], char fnameBot[]) {

  Grid2D_Quad_Block **mesh;
  int i, j, k;
  double pspan;
  double r, theta;
  ofstream out_file_top, out_file_bot;

  out_file_top.open(fnameTop, ios::out);
  out_file_bot.open(fnameBot, ios::out);

  assert(!out_file_top.fail());
  assert(!out_file_bot.fail());

  cout << "\n NASA Rotor 37 3D structured mesh generation.";

  out_file_top << "TITLE = \"" << "NASA ROTOR 37 3D MULTI-BLOCK MESH: STRUCTURED FORMAT"
	       << "\"" << "\n"
	       << "VARIABLES = \"x\" \\ \n"
	       << "\"y\" \n"
	       << "\"z\" \n"
	       << "ZONE T = \"3D MESH - UPPER\"\n"
	       << "I = " << 3*Number_of_Cells_Idir+1 << " \\ \n"
	       << "J = " << Number_of_Cells_Jdir+1 << " \\ \n"
	       << "K = " << kMax+1 << " \\ \n"
	       << "F = POINT \n"
	       << "DT=(SINGLE, SINGLE, SINGLE)\n";

  out_file_bot << "\nTITLE = \"" << "NASA ROTOR 37 3D MULTI-BLOCK MESH: STRUCTURED FORMAT"
	       << "\"" << "\n"
	       << "VARIABLES = \"x\" \\ \n"
	       << "\"y\" \n"
	       << "\"z\" \n"
	       << "ZONE T = \"3D MESH - LOWER\"\n"
	       << "I = " << 3*Number_of_Cells_Idir+1 << " \\ \n"
	       << "J = " << Number_of_Cells_Jdir+1 << " \\ \n"
	       << "K = " << kMax+1 << " \\ \n"
	       << "F = POINT \n"
	       << "DT=(SINGLE, SINGLE, SINGLE)\n";
  
  for ( k = 0; k<=kMax; ++k ) {

    //calculate %span in fractional form
    pspan=double(k)/double(kMax);

    //apply stretching function
    pspan=StretchingFcn(pspan,
			1.01,
			0,
			STRETCHING_FCN_MINMAX_CLUSTERING);

    pspan=pspan*HUNDRED;
    //is pspan is too small problems occur - ??
    if(pspan<TOLER) pspan=ZERO;
    
    cout << setprecision(4);
    cout << "\n  Status of 3D mesh generation: " 
         << int((double(k-1)/double(kMax))*HUNDRED) 
         << "% completed, computing " 
         << pspan
         << "% span..." ;
    cout << setprecision(6);
    cout.flush();

    //generate grid for this %span
    mesh = genMeshH_3x2_AUTO(mesh,
                             pspan, zMin, zMax,  
                             Number_of_Cells_Idir, Number_of_Cells_Jdir, 2, 0);

    cout << "\n  Done.  Writing to file...";
    cout.flush();
           
    for ( j  = mesh[0][0].JNl ; j <= mesh[0][0].JNu ; ++j ) {
      for ( i = mesh[0][0].INl ; i <= mesh[0][0].INu ; ++i ) {
	 r=findRadius(pspan, mesh[0][0].Node[i][j].X.x);
	 theta=mesh[0][0].Node[i][j].X.y/r;
	 out_file_bot << r*cos(theta) << "\n"
		      << r*sin(theta) << "\n"
		      << mesh[0][0].Node[i][j].X.x << "\n";
	
      }
      for ( i = mesh[1][0].INl+1 ; i <= mesh[1][0].INu ; ++i ) {
	 r=findRadius(pspan, mesh[1][0].Node[i][j].X.x);
	 theta=mesh[1][0].Node[i][j].X.y/r;
	 out_file_bot << r*cos(theta) << "\n"
		      << r*sin(theta) << "\n"
		      << mesh[1][0].Node[i][j].X.x << "\n";
	 
      } 
      for ( i = mesh[2][0].INl+1 ; i <= mesh[2][0].INu ; ++i ) {
	 r=findRadius(pspan, mesh[2][0].Node[i][j].X.x);
	 theta=mesh[2][0].Node[i][j].X.y/r;
	 out_file_bot << r*cos(theta) << "\n"
		      << r*sin(theta) << "\n"
		      << mesh[2][0].Node[i][j].X.x << "\n"; 
      } 
     
    } 
   
    for ( j  = mesh[0][1].JNl ; j <= mesh[0][1].JNu ; ++j ) {
       for ( i = mesh[0][1].INl ; i <= mesh[0][1].INu ; ++i ) {
	  r=findRadius(pspan, mesh[0][1].Node[i][j].X.x);
	  theta=mesh[0][1].Node[i][j].X.y/r;
	  out_file_top << r*cos(theta) << "\n"
		       << r*sin(theta) << "\n"
		       << mesh[0][1].Node[i][j].X.x << "\n";
       } 
       for ( i = mesh[1][1].INl+1 ; i <= mesh[1][1].INu ; ++i ) {
	  r=findRadius(pspan, mesh[1][1].Node[i][j].X.x);
	  theta=mesh[1][1].Node[i][j].X.y/r;
	  out_file_top << r*cos(theta) << "\n"
		       << r*sin(theta) << "\n"
		       << mesh[1][1].Node[i][j].X.x << "\n";
       } 
       for ( i = mesh[2][1].INl+1 ; i <= mesh[2][1].INu ; ++i ) {
	  r=findRadius(pspan, mesh[2][1].Node[i][j].X.x);
	  theta=mesh[2][1].Node[i][j].X.y/r;
	  out_file_top << r*cos(theta) << "\n"
		       << r*sin(theta) << "\n"
	  	       << mesh[2][1].Node[i][j].X.x << "\n";
       }        
    }

    mesh = Deallocate_Multi_Block_Grid(mesh, 3, 2);

  } /* endfor */

  out_file_top.close();
  out_file_bot.close();

}

/**********************************************************************
 * Routine: output_Tetin                                              *
 *                                                                    *
 * Creates tetin file for NASA ROTOR 37 geometry.                     *
 *                                                                    *
 **********************************************************************/
inline void NASARotor37::output_Tetin(char fname[], int layers) {
  int i, j, l, n, counter;
  int *leadI, *trailI; //index for leading and trailign edge of each contour
  double pspan, r, theta, dtheta;
  double hub_shroud_size, blade_LT_size, blade_surf_size, blade_curve_size, bound_surf_size, 
         bound_curve_size, global_max_size;
  Vector3D v;   
  Vector3D **geom, **camberB, **hubB, **shroudB, **upstreamB, **downstreamB;
  Spline2D cs, camber, camber2;
  ofstream out_file;

  out_file.open(fname, ios::out);
  assert(!out_file.fail());

  leadI = new int[layers];
  trailI = new int[layers];

  geom = new Vector3D*[layers];
  camberB = new Vector3D*[layers];
  hubB = new Vector3D*[layers];
  shroudB = new Vector3D*[layers];
  upstreamB = new Vector3D*[layers];
  downstreamB = new Vector3D*[layers];
  for ( i = 0; i <= layers-1 ; ++i ) {
     geom[i] = new Vector3D[POINTS_PER_SURFACE_R37];
     camberB[i] = new Vector3D[40*3-1];
     hubB[i] = new Vector3D[40*3-1];
     shroudB[i] = new Vector3D[40*3-1];
     upstreamB[i] = new Vector3D[layers];
     downstreamB[i] = new Vector3D[layers];
  } /* endfor */
   
  counter=1;
  dtheta=M_PI/num_blades;
     
  //set tetra sizes
  hub_shroud_size=0.005; 
  blade_LT_size=0.0004;
  blade_surf_size=0.003; 
  blade_curve_size=0.001;
  bound_surf_size=0.05;
  bound_curve_size=0.01;
  global_max_size=0.05;

  //process blade data
  for ( i=0; i<layers; ++i ) {
    pspan= 100.0/(layers-1) * i;
    cs=getBladeCS(pspan);
    findLT(pspan, trailI[i], leadI[i]);
    for ( j=0; j<POINTS_PER_SURFACE_R37; ++j ) {
      r = findRadius(pspan, cs.Xp[j].x);
      theta = cs.Xp[j].y / r;
      geom[i][j].z = cs.Xp[j].x; 
      geom[i][j].x = r*cos(theta);
      geom[i][j].y = r*sin(theta);
    }
  }
 
  //create tetin file
  out_file << "// tetin file for NASA ROTOR 37 geometry - generated by Class NASARotor37\n\n";

  //define BLADE SURFACES
  for ( l=0; l<layers-1; ++l) {
    for ( n=0; n<POINTS_PER_SURFACE_R37-2; ++n) {
      //make mesh finer near leading and trailing edge
      if(n<=10 || (n>R37GEOM_NUMPTS-10 && n<R37GEOM_NUMPTS+10) || n>POINTS_PER_SURFACE_R37-10)
	out_file << "define_surface family BLADE name BLADE/" << counter++ << " tetra_size " 
		 << blade_LT_size << "\n";
      else
	out_file << "define_surface family BLADE name BLADE/" << counter++ << " tetra_size "
		 << blade_surf_size << "\n";
      out_file << "unstruct_mesh n_points 4 n_triangles 2\n"
	       << geom[l][n] << "\n"
	       << geom[l][n+1] << "\n"
	       << geom[l+1][n+1] << "\n"
	       << geom[l+1][n] << "\n"
	       << "0 1 2\n"
	       << "0 2 3\n";
    }
    //make mesh finer near leading and trailing edge
    out_file << "define_surface family BLADE name BLADE/" << counter++ << " tetra_size "
	     << blade_LT_size << "\n"
	     << "unstruct_mesh n_points 4 n_triangles 2\n"
	     << geom[l][POINTS_PER_SURFACE_R37-2] << "\n"
	     << geom[l][0] << "\n"
	     << geom[l+1][0] << "\n"
	     << geom[l+1][POINTS_PER_SURFACE_R37-2] << "\n"
	     << "0 1 2\n"
	     << "0 2 3\n";
  }

  //define BLADE curves 

  //bottom contour - divide into four
  out_file << "define_curve family CURVES name CURVES/" << 0 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << R37GEOM_NUMPTS/2-trailI[0]+1 << " n_edges " 
	   << R37GEOM_NUMPTS/2-trailI[0] << "\n";
  for ( n=trailI[0]; n<=R37GEOM_NUMPTS/2; ++n )
    out_file << geom[0][n] << "\n";
  for ( n=0; n<R37GEOM_NUMPTS/2-trailI[0]; ++n )
    out_file << n << " " << n+1 << "\n";

  out_file << "define_curve family CURVES name CURVES/" << 1 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << leadI[0]-R37GEOM_NUMPTS/2+1 << " n_edges " 
	   << leadI[0]-R37GEOM_NUMPTS/2 << "\n";
  for ( n=R37GEOM_NUMPTS/2; n<=leadI[0]; ++n )
    out_file << geom[0][n] << "\n";
  for ( n=0; n<leadI[0]-R37GEOM_NUMPTS/2; ++n )
    out_file << n << " " << n+1 << "\n";

  out_file << "define_curve family CURVES name CURVES/" << 2 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << POINTS_PER_SURFACE_R37*3/4-leadI[0]+1 << " n_edges " 
	   << POINTS_PER_SURFACE_R37*3/4-leadI[0] << "\n";
  for ( n=leadI[0]; n<=POINTS_PER_SURFACE_R37*3/4; ++n )
    out_file << geom[0][n] << "\n";
  for ( n=0; n<POINTS_PER_SURFACE_R37*3/4-leadI[0]; ++n )
    out_file << n << " " << n+1 << "\n";

  out_file << "define_curve family CURVES name CURVES/" << 3 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << (POINTS_PER_SURFACE_R37-1)-POINTS_PER_SURFACE_R37*3/4+trailI[0]+1 << " n_edges " 
	   << (POINTS_PER_SURFACE_R37-1)-POINTS_PER_SURFACE_R37*3/4+trailI[0] << "\n";
  for ( n=POINTS_PER_SURFACE_R37*3/4; n<=POINTS_PER_SURFACE_R37-1; ++n )
    out_file << geom[0][n] << "\n";
  for ( n=1; n<=trailI[0]; ++n )
    out_file << geom[0][n] << "\n";
  for ( n=0; n<POINTS_PER_SURFACE_R37-1-POINTS_PER_SURFACE_R37*3/4+trailI[0]; ++n )
    out_file << n << " " << n+1 << "\n";   

  // - top contour - divide into four curves
  out_file << "define_curve family CURVES name CURVES/" << 4 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << R37GEOM_NUMPTS/2-trailI[layers-1]+1 << " n_edges " 
	   << R37GEOM_NUMPTS/2-trailI[layers-1] << "\n";
  for ( n=trailI[layers-1]; n<=R37GEOM_NUMPTS/2; ++n )
    out_file << geom[layers-1][n] << "\n";
  for ( n=0; n<R37GEOM_NUMPTS/2-trailI[layers-1]; ++n )
    out_file << n << " " << n+1 << "\n";

  out_file << "define_curve family CURVES name CURVES/" << 5 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << leadI[layers-1]-R37GEOM_NUMPTS/2+1 << " n_edges " 
	   << leadI[layers-1]-R37GEOM_NUMPTS/2 << "\n";
  for ( n=R37GEOM_NUMPTS/2; n<=leadI[layers-1]; ++n )
    out_file << geom[layers-1][n] << "\n";
  for ( n=0; n<leadI[layers-1]-R37GEOM_NUMPTS/2; ++n )
    out_file << n << " " << n+1 << "\n";

  out_file << "define_curve family CURVES name CURVES/" << 6 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << POINTS_PER_SURFACE_R37*3/4-leadI[layers-1]+1 << " n_edges " 
	   << POINTS_PER_SURFACE_R37*3/4-leadI[layers-1] << "\n";
  for ( n=leadI[layers-1]; n<=POINTS_PER_SURFACE_R37*3/4; ++n )
    out_file << geom[layers-1][n] << "\n";
  for ( n=0; n<POINTS_PER_SURFACE_R37*3/4-leadI[layers-1]; ++n )
    out_file << n << " " << n+1 << "\n";

  out_file << "define_curve family CURVES name CURVES/" << 7 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << (POINTS_PER_SURFACE_R37-1)-POINTS_PER_SURFACE_R37*3/4+trailI[layers-1]+1 << " n_edges " 
	   << (POINTS_PER_SURFACE_R37-1)-POINTS_PER_SURFACE_R37*3/4+trailI[layers-1] << "\n";
  for( n=POINTS_PER_SURFACE_R37*3/4; n<=POINTS_PER_SURFACE_R37-1; ++n )
    out_file << geom[layers-1][n] << "\n";
  for ( n=1; n<=trailI[layers-1]; ++n)
    out_file << geom[layers-1][n] << "\n";
  for ( n=0; n<POINTS_PER_SURFACE_R37-1-POINTS_PER_SURFACE_R37*3/4+trailI[layers-1]; ++n)
    out_file << n << " " << n+1 << "\n";


  //leading and trailing edge
  out_file << "define_curve family CURVES name CURVES/" << 8 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( n=0; n<layers; ++n)
    out_file << geom[n][leadI[n]] << "\n";
  for ( n=0; n<layers-1; ++n)
    out_file << n << " " << n+1 << "\n";
  
  out_file << "define_curve family CURVES name CURVES/" << 9 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( n=0; n<layers; ++n)
    out_file << geom[n][trailI[n]] << "\n";
  for ( n=0; n<layers-1; ++n)
    out_file << n << " " << n+1 << "\n";

  out_file << "define_curve family CURVES name CURVES/" << 10 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( n=0; n<layers; ++n)
    out_file << geom[n][R37GEOM_NUMPTS/2] << "\n";
  for ( n=0; n<layers-1; ++n)
    out_file << n << " " << n+1 << "\n";  

  out_file << "define_curve family CURVES name CURVES/" << 11 << " tetra_size "
	   << blade_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( n=0; n<layers; ++n)
    out_file << geom[n][POINTS_PER_SURFACE_R37*3/4] << "\n";
  for ( n=0; n<layers-1; ++n)
    out_file << n << " " << n+1 << "\n";  

  //blade points
  out_file << "prescribed_point " << geom[0][trailI[0]] << " name POINTS/" << 0 
	   << " family POINTS\n";
  out_file << "prescribed_point " << geom[0][leadI[0]] << " name POINTS/" << 1
	   << " family POINTS\n";
  out_file << "prescribed_point " << geom[0][R37GEOM_NUMPTS/2] << " name POINTS/" << 2
	   << " family POINTS\n";
  out_file << "prescribed_point " << geom[0][POINTS_PER_SURFACE_R37*3/4] << " name POINTS/" << 3 
	   << " family POINTS\n";
  out_file << "prescribed_point " << geom[layers-1][trailI[layers-1]] << " name POINTS/" << 4 
	   << " family POINTS\n";
  out_file << "prescribed_point " << geom[layers-1][leadI[layers-1]] << " name POINTS/" << 5
	   << " family POINTS\n";
  out_file << "prescribed_point " << geom[layers-1][R37GEOM_NUMPTS/2] << " name POINTS/" << 6
	   << " family POINTS\n";
  out_file << "prescribed_point " << geom[layers-1][POINTS_PER_SURFACE_R37*3/4] << " name POINTS/" << 7
	   << " family POINTS\n";

  //process camber line data - right boundary
  for ( i=0; i<layers; ++i) {
    pspan= 100.0/(layers-1) * i;
    camber=getCamber(pspan, 40, z_up,z_d);
    for ( j=0; j<camber.np; ++j) {
      r = findRadius(pspan, camber.Xp[j].x);
      theta = camber.Xp[j].y / r;
      theta+=dtheta;
      camberB[i][j].z = camber.Xp[j].x; 
      camberB[i][j].x = r*cos(theta);
      camberB[i][j].y = r*sin(theta);
    }
  }

  counter=0;
  //define camber_A BOUNDARY surface
  for ( l=0; l<layers-1; ++l) {
    for ( n=0; n<camber.np-1; ++n) {
      out_file << "define_surface family CAMBER_A name CAMBER_A/" << counter++ << " tetra_size "
	       << bound_surf_size << "\n"
	       << "unstruct_mesh n_points 4 n_triangles 2\n"
	       << camberB[l][n] << "\n"
	       << camberB[l][n+1] << "\n"
	       << camberB[l+1][n+1] << "\n"
	       << camberB[l+1][n] << "\n"
	       << "0 1 2\n"
	       << "0 2 3\n";
    }
  }

  //define CAMBER_A BOUNDARY points
  out_file << "prescribed_point " << camberB[0][0] << " name POINTS/" << 8 
	   << " family POINTS\n";
  out_file << "prescribed_point " << camberB[0][camber.np-1] << " name POINTS/" << 9
	   << " family POINTS\n";
  out_file << "prescribed_point " << camberB[layers-1][0] << " name POINTS/" << 10
	   << " family POINTS\n";
  out_file << "prescribed_point " << camberB[layers-1][camber.np-1] << " name POINTS/" << 11 
	   << " family POINTS\n";

  //define CAMBER_A BOUNDARY curves 
  // - bottom-right curve
  out_file << "define_curve family CURVES name CURVES/" << 12 << " tetra_size "
	   << hub_shroud_size << "\n"
	   << "unstruct_curve n_points " << camber.np << " n_edges " 
	   << camber.np-1 << "\n";
  for ( n=0; n<camber.np; ++n)
    out_file << camberB[0][n] << "\n";
  for ( n=0; n<camber.np-1; ++n)
    out_file << n << " " << n+1 << "\n";
  // - top-right curve
  out_file << "define_curve family CURVES name CURVES/" << 13 << " tetra_size "
	   << hub_shroud_size << "\n"
	   << "unstruct_curve n_points " << camber.np << " n_edges " 
	   << camber.np-1 << "\n";
  for ( n=0; n<camber.np; ++n)
    out_file << camberB[layers-1][n] << "\n";
  for ( n=0; n<camber.np-1; ++n)
    out_file << n << " " << n+1 << "\n";
  // - upstream-right curve
  out_file << "define_curve family CURVES name CURVES/" << 14 << " tetra_size "
	   << bound_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( l=0; l<layers; ++l)
    out_file << camberB[l][0] << "\n";
  for ( l=0; l<layers-1; ++l)
    out_file << l << " " << l+1 << "\n";
  // - downstream-right curve
  out_file << "define_curve family CURVES name CURVES/" << 15 << " tetra_size "
	   << bound_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( l=0; l<layers; ++l)
    out_file << camberB[l][camber.np-1] << "\n";
  for ( l=0; l<layers-1; ++l)
    out_file << l << " " << l+1 << "\n";

  //surface curve 1
  out_file << "define_curve family CURVES name CURVES/" << 16 << " tetra_size "
	   << bound_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( l=0; l<layers; ++l)
    out_file << camberB[l][40] << "\n";
  for ( l=0; l<layers-1; ++l)
    out_file << l << " " << l+1 << "\n";


  //surface curve 2
  out_file << "define_curve family CURVES name CURVES/" << 17 << " tetra_size "
	   << bound_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( l=0; l<layers; ++l)
    out_file << camberB[l][80] << "\n";
  for ( l=0; l<layers-1; ++l)
    out_file << l << " " << l+1 << "\n";

  //surface curve 3
  out_file << "define_curve family CURVES name CURVES/" << 18 << " tetra_size "
	   << bound_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( l=0; l<layers; ++l)
    out_file << camberB[l][60] << "\n";
  for ( l=0; l<layers-1; ++l)
    out_file << l << " " << l+1 << "\n";

  //process camber line data - left boundary
  for ( i=0; i<layers; ++i) {
    for ( j=0; j<camber.np; ++j) {
      //rotate camber line about z axis
      v=camberB[i][j];
      camberB[i][j].x=v.x*cos(-2*dtheta)-v.y*sin(-2*dtheta);
      camberB[i][j].y=v.x*sin(-2*dtheta)+v.y*cos(-2*dtheta);
    }
  }

  counter=0;
  //define CAMBER_B BOUNDARY SURFACES
  for ( l=0; l<layers-1; ++l) {
    for ( n=0; n<camber.np-1; ++n) {
      out_file << "define_surface family CAMBER_B name CAMBER_B/" << counter++ << " tetra_size "
	       << bound_surf_size << "\n"
	       << "unstruct_mesh n_points 4 n_triangles 2\n"
	       << camberB[l][n] << "\n"
	       << camberB[l][n+1] << "\n"
	       << camberB[l+1][n+1] << "\n"
	       << camberB[l+1][n] << "\n"
	       << "0 1 2\n"
	       << "0 2 3\n";
    }    
  }

  //define CAMBER_B BOUNDARY points
  out_file << "prescribed_point " << camberB[0][0] << " name POINTS/" << 12 
	   << " family POINTS\n";
  out_file << "prescribed_point " << camberB[0][camber.np-1] << " name POINTS/" << 13
	   << " family POINTS\n";
  out_file << "prescribed_point " << camberB[layers-1][0] << " name POINTS/" << 14
	   << " family POINTS\n";
  out_file << "prescribed_point " << camberB[layers-1][camber.np-1] << " name POINTS/" << 15 
	   << " family POINTS\n";

  //define CAMBER_B BOUNDARY curves 
  // - bottom-left curve
  out_file << "define_curve family CURVES name CURVES/" << 19 << " tetra_size "
	   << hub_shroud_size << "\n"
	   << "unstruct_curve n_points " << camber.np << " n_edges " 
	   << camber.np-1 << "\n";
  for ( n=0; n<camber.np; ++n)
    out_file << camberB[0][n] << "\n";
  for ( n=0; n<camber.np-1; ++n)
    out_file << n << " " << n+1 << "\n";
  // - top-left curve
  out_file << "define_curve family CURVES name CURVES/" << 20 << " tetra_size "
	   << hub_shroud_size << "\n"
	   << "unstruct_curve n_points " << camber.np << " n_edges " 
	   << camber.np-1 << "\n";
  for ( n=0; n<camber.np; ++n)
    out_file << camberB[layers-1][n] << "\n";
  for ( n=0; n<camber.np-1; ++n)
    out_file << n << " " << n+1 << "\n";
  // - upstream-left curve
  out_file << "define_curve family CURVES name CURVES/" << 21 << " tetra_size "
	   << bound_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( l=0; l<layers; ++l)
    out_file << camberB[l][0] << "\n";
  for ( l=0; l<layers-1; ++l)
    out_file << l << " " << l+1 << "\n";
  // - downstream-left curve
  out_file << "define_curve family CURVES name CURVES/" << 22 << " tetra_size "
	   << bound_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( l=0; l<layers; ++l)
    out_file << camberB[l][camber.np-1] << "\n";
  for ( l=0; l<layers-1; ++l)
    out_file << l << " " << l+1 << "\n";

  //surface curve 1
  out_file << "define_curve family CURVES name CURVES/" << 23 << " tetra_size "
	   << bound_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( l=0; l<layers; ++l)
    out_file << camberB[l][40] << "\n";
  for ( l=0; l<layers-1; ++l)
    out_file << l << " " << l+1 << "\n";

  //surface curve 2
  out_file << "define_curve family CURVES name CURVES/" << 24 << " tetra_size "
	   << bound_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( l=0; l<layers; ++l)
    out_file << camberB[l][80] << "\n";
  for ( l=0; l<layers-1; ++l)
    out_file << l << " " << l+1 << "\n";

 //surface curve 3
  out_file << "define_curve family CURVES name CURVES/" << 25 << " tetra_size "
	   << bound_curve_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( l=0; l<layers; ++l)
    out_file << camberB[l][60] << "\n";
  for ( l=0; l<layers-1; ++l)
    out_file << l << " " << l+1 << "\n";

  //process hub data
  camber=getCamber(0, 40, z_up,z_d);
  for ( i=0; i<layers; ++i) {
    for ( j=0; j<camber.np; ++j) {
      r = findRadius(0, camber.Xp[j].x);
      theta = camber.Xp[j].y / r;
      theta=theta-dtheta+2*dtheta/(layers-1)*i;
      hubB[i][j].z = camber.Xp[j].x; 
      hubB[i][j].x = r*cos(theta);
      hubB[i][j].y = r*sin(theta);
    }
  }

  //define hub boundary
  counter=0;
  for ( l=0; l<layers-1; ++l) {
    for ( n=0; n<camber.np-1; ++n) {
      out_file << "define_surface family HUB name HUB/" << counter++ << " tetra_size "
	       << hub_shroud_size << "\n"
	       << "unstruct_mesh n_points 4 n_triangles 2\n"
	       << hubB[l][n] << "\n"
	       << hubB[l][n+1] << "\n"
	       << hubB[l+1][n+1] << "\n"
	       << hubB[l+1][n] << "\n"
	       << "0 1 2\n"
	       << "0 2 3\n";
    }
  }

  //process shroud data
  camber=getCamber(100, 40, z_up,z_d);
  for ( i=0; i<layers; ++i) {
    for ( j=0; j<camber.np; ++j) {
      r = findRadius(100, camber.Xp[j].x);
      theta = camber.Xp[j].y / r;
      theta=theta-dtheta+2*dtheta/(layers-1)*i;
      shroudB[i][j].z = camber.Xp[j].x; 
      shroudB[i][j].x = r*cos(theta);
      shroudB[i][j].y = r*sin(theta);
    }
  }

  //define shroud boundary surfaces
  counter=0;
  for ( l=0; l<layers-1; ++l) {
    for ( n=0; n<camber.np-1; ++n) {
      out_file << "define_surface family SHROUD name SHROUD/" << counter++ << " tetra_size "
	       << hub_shroud_size << "\n"
	       << "unstruct_mesh n_points 4 n_triangles 2\n"
	       << shroudB[l][n] << "\n"
	       << shroudB[l][n+1] << "\n"
	       << shroudB[l+1][n+1] << "\n"
	       << shroudB[l+1][n] << "\n"
	       << "0 1 2\n"
	       << "0 2 3\n";
    }
  }
   
  //create upstream boundary surface
  for ( i=0; i<layers; ++i) {
    for ( j=0; j<layers; ++j) {
      //rotate camber line about z axis
      v=camberB[i][0];
      upstreamB[i][j].x=v.x*cos(2*dtheta/(layers-1)*j)-v.y*sin(2*dtheta/(layers-1)*j);
      upstreamB[i][j].y=v.x*sin(2*dtheta/(layers-1)*j)+v.y*cos(2*dtheta/(layers-1)*j);
      upstreamB[i][j].z=v.z;
    }
  }
    
  //define upstream boundary
  counter=0;
  for ( l=0; l<layers-1; ++l) {
    for ( n=0; n<layers-1; ++n) {
      out_file << "define_surface family UPSTREAM name UPSTREAM/" << counter++ << " tetra_size "
	       << bound_surf_size << "\n"
	       << "unstruct_mesh n_points 4 n_triangles 2\n"
	       << upstreamB[l][n] << "\n"
	       << upstreamB[l][n+1] << "\n"
	       << upstreamB[l+1][n+1] << "\n"
	       << upstreamB[l+1][n] << "\n"
	       << "0 1 2\n"
	       << "0 2 3\n";
    }
  }

  //define UPSTREAM BOUNDARY curves 
  out_file << "define_curve family CURVES name CURVES/" << 26 << " tetra_size "
	   << hub_shroud_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( n=0; n<layers; ++n)
    out_file << upstreamB[0][n] << "\n";
  for ( n=0; n<layers-1; ++n)
    out_file << n << " " << n+1 << "\n";
  out_file << "define_curve family CURVES name CURVES/" << 27 << " tetra_size "
	   << hub_shroud_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( n=0; n<layers; ++n)
    out_file << upstreamB[layers-1][n] << "\n";
  for ( n=0; n<layers-1; ++n)
    out_file << n << " " << n+1 << "\n";

  //create downstream boundary surface
  for ( i=0; i<layers; ++i ) {
    for ( j = 0; j<layers; ++j ) {
      //rotate camber line about z axis
      v=camberB[i][40*3-2];
      downstreamB[i][j].x=v.x*cos(2*dtheta/(layers-1)*j)-v.y*sin(2*dtheta/(layers-1)*j);
      downstreamB[i][j].y=v.x*sin(2*dtheta/(layers-1)*j)+v.y*cos(2*dtheta/(layers-1)*j);
      downstreamB[i][j].z=v.z;
    }
  }
    
  //define downstream boundary
  counter=0;
  for ( l=0; l<layers-1; ++l) {
    for ( n=0; n<layers-1; ++n) {
      out_file << "define_surface family DOWNSTREAM name DOWNSTREAM/" << counter++ << " tetra_size "
	       << bound_surf_size << "\n"
	       << "unstruct_mesh n_points 4 n_triangles 2\n"
	       << downstreamB[l][n] << "\n"
	       << downstreamB[l][n+1] << "\n"
	       << downstreamB[l+1][n+1] << "\n"
	       << downstreamB[l+1][n] << "\n"
	       << "0 1 2\n"
	       << "0 2 3\n";
    }
  }

  //define DOWNSTREAM BOUNDARY curves 
  out_file << "define_curve family CURVES name CURVES/" << 28 << " tetra_size "
	   << hub_shroud_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( n=0; n<layers; ++n)
    out_file << downstreamB[0][n] << "\n";
  for ( n=0; n<layers-1; ++n)
    out_file << n << " " << n+1 << "\n";
  out_file << "define_curve family CURVES name CURVES/" << 29 << " tetra_size "
	   << hub_shroud_size << "\n"
	   << "unstruct_curve n_points " << layers << " n_edges " 
	   << layers-1 << "\n";
  for ( n=0; n<layers; ++n)
    out_file << downstreamB[layers-1][n] << "\n";
  for ( n=0; n<layers-1; ++n)
    out_file << n << " " << n+1 << "\n";

  //define material point
  out_file << "material_point " << upstreamB[layers/2][layers/2].x << " " << upstreamB[layers/2][layers/2].y 
	   << " " << upstreamB[layers/2][layers/2].z+0.005 
	   << " family VOLUME_MATERIAL name VOLUME_MATERIAL/0\n";

  //set affix
  out_file << "affix 0\n";

  //set model parameters
  out_file << "define_model " << global_max_size << " reference_size 1\n";

  out_file.close();

  delete []leadI; leadI = NULL;
  delete []trailI; trailI = NULL;

  for ( i = 0; i <= layers-1 ; ++i ) {
     delete []geom[i]; geom[i] = NULL;
     delete []camberB[i]; camberB[i] = NULL;
     delete []hubB[i]; hubB[i] = NULL;
     delete []shroudB[i]; shroudB[i] = NULL;
     delete []upstreamB[i]; upstreamB[i] = NULL;
     delete []downstreamB[i]; downstreamB[i] = NULL;
  } /* endfor */
  delete []geom; geom = NULL; 
  delete []camberB; camberB = NULL; 
  delete []hubB; hubB = NULL;
  delete []shroudB; shroudB = NULL;
  delete []downstreamB; downstreamB = NULL;

}

#endif // _NASA_ROTOR37_INCLUDED
