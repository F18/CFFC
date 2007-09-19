/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: Medium2DState.cc                                           **
 **                                                                  **
 ** Description: The radiation state class contains the properties   **
 **              of a gray, absorbing, emitting, anisotropically     **
 **              scattering medium.  This file defines the state     ** 
 **              class.                                              **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            04/03/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/
#include "Medium2DState.h"


/********************************************************
 * Static member initialization                         *
 ********************************************************/
int    Medium2D_State :: Nband          = 0;     // number of bands (& quad pts)
int    Medium2D_State :: NUM_VAR_MEDIUM2D = 0;   // total number of variables
SNBCK* Medium2D_State :: SNBCKdata      = NULL;  // SNBCK data object
int    Medium2D_State :: Absorb_Type    = MEDIUM2D_ABSORB_GRAY; // absorbsion model
bool   Medium2D_State :: Scatter_Iso    = true;  // isotropic or anisotropic scattering
Vector2D_Function<Medium2D_State>* Medium2D_State :: Field = NULL;   // functor pointer


/********************************************************
 * Compute values and initialize state.                 *
 ********************************************************/
void Medium2D_State :: SetInitialValues( const double &Pressure,
					 const double &Temperature,
					 const double &xco,
					 const double &xh2o,
					 const double &xco2,
					 const double &xo2,
					 const double &fsoot,
					 const double &AbsorptionCoef,
					 const double &ScatteringCoef)
{

  //------------------------------------------------
  // Absorbsion coefficient, Blackbody intensity 
  //------------------------------------------------
  // Use SNBCK
  if (Absorb_Type == MEDIUM2D_ABSORB_SNBCK) {
    SNBCKdata->CalculateAbsorb( Pressure/PRESSURE_STDATM, //[atm]
				Temperature,              //[K]
				xco,
				xh2o,
				xco2,
				xo2,
				fsoot,
				kappa );
    SNBCKdata->CalculatePlanck( Temperature, Ib );

  // Use Gray Gas (ie. constant)
  } else if (Absorb_Type == MEDIUM2D_ABSORB_GRAY) {
    for (int v=0; v<Nband; v++) {
      kappa[v] = AbsorptionCoef;
      Ib   [v] = BlackBody(Temperature);
    } // endfor

  // error
  } else{
    cerr << "Medium2D_State::SetInitialValues() - Invalid flag for Absorbsion model\n";
    exit(1);
  }

  //------------------------------------------------
  // Scattering coefficient 
  //------------------------------------------------
  // scattering coefficient always assumed gray
  for (int v=0; v<Nband; v++) { sigma[v] = ScatteringCoef; }

}



/********************************************************
 * Setup Static variables.                              *
 ********************************************************/
void Medium2D_State :: SetupStatic( const int &i_Absorb_Type,
				    const int &i_Scattering_Type,
				    const SNBCK_Input_Parameters &SNBCK_IP,
				    const char* PATH) {

  // deallocate to be sure
  DeallocateSNBCK();

  //------------------------------------------------
  // Absorbsion model 
  //------------------------------------------------
  // set the absorption type flag
  Absorb_Type = i_Absorb_Type;
  
  // GRAY
  if (Absorb_Type == MEDIUM2D_ABSORB_GRAY) {
    Nband = 1;

  // SNBCK
  } else if (Absorb_Type == MEDIUM2D_ABSORB_SNBCK) {
    AllocateSNBCK();
    SNBCKdata->Setup(SNBCK_IP, PATH);
    Nband = SNBCKdata->NumVar();

  // ERROR
  } else {
    cerr << "Medium2D_State::SetupState - Invalid flag for gas type\n";
    exit(-1);
  } // endif

  // set the isotropic scattering flag
  Scatter_Iso = (i_Scattering_Type == RTE2D_SCATTER_ISO);

  // set number of variables
  NUM_VAR_MEDIUM2D = 3*Nband;

}


/********************************************************
 * Set all scalar fields to the same function.          *
 ********************************************************/
// CONSTANT
void Medium2D_State :: SetConstantField(const Medium2D_State &M)
{
  // deallocate to be sure
  DeallocateField();

  // Create objects.  You will no longer own this dynamically created object,
  // once passed to the class 'Vector2D_SpecFunction', it's destructor
  // will take care of it.
  ConstantFunc<Medium2D_State>* funcObj = new ConstantFunc<Medium2D_State>(M);

  // assign all the fields
  Field = new Vector2D_SpecFunction< Medium2D_State, 
                                     ConstantFunc<Medium2D_State> >(funcObj);
}

// DISCONTINUOUS
void Medium2D_State :: SetDiscontinuousField( const Medium2D_State &inner, 
					      const Medium2D_State &outer, 
					      const Vector2D &x_SW, 
					      const Vector2D &x_NE )
{
  // deallocate to be sure
  DeallocateField();

  // Create objects.  You will no longer own this dynamically created object,
  // once passed to the class 'Vector2D_SpecFunction', it's destructor
  // will take care of it.
  DiscontinuousFunc<Medium2D_State>* funcObj = new DiscontinuousFunc<Medium2D_State>(inner,
										     outer,
										     x_SW,
										     x_NE);

  // assign all the fields
  Field = new Vector2D_SpecFunction< Medium2D_State, 
                                     DiscontinuousFunc<Medium2D_State> >(funcObj);

}

/********************************************************
 * Compute NEW values for the medium state dependent on *
 * relevant paremters.                                  *
 ********************************************************/
void Medium2D_State :: ComputeNewState( const double &Pressure,
					const double &Temperature,
					const double &xco,
					const double &xh2o,
					const double &xco2,
					const double &xo2,
					const double &fsoot )
{

  //------------------------------------------------
  // Absorbsion coefficient, Blackbody intensity 
  //------------------------------------------------
  // Use SNBCK
  if (Absorb_Type == MEDIUM2D_ABSORB_SNBCK) {
    SNBCKdata->CalculateAbsorb( Pressure/PRESSURE_STDATM, //[atm]
				Temperature,              //[K]
				xco,
				xh2o,
				xco2,
				xo2,
				fsoot,
				kappa );
    SNBCKdata->CalculatePlanck( Temperature, Ib );

  // Use Gray Gas (ie. constant)
  } else if (Absorb_Type == MEDIUM2D_ABSORB_GRAY) {
    for (int v=0; v<Nband; v++) {
      //kappa[v] = AbsorptionCoef; // <-- it is a treated as a specified constant
      Ib   [v] = BlackBody(Temperature);
    } // endfor

  // error
  } else{
    cerr << "Medium2D_State::SetInitialValues() - Invalid flag for Absorbsion model\n";
    exit(1);
  }

  //------------------------------------------------
  // Scattering coefficient 
  //------------------------------------------------
  // scattering coefficient always assumed gray
  // for (int v=0; v<Nband; v++) { sigma[v] = ScatteringCoef; } // <-- it is a treated specified constant

}
