/*!\file AccuracyAssessment2DMultiBlock.cc
  \brief Source file implementing/setting member functions/variables 
  prototyped in AccuracyAssessment2D_MultiBlock class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "AccuracyAssessment2DMultiBlock.h"
#include "../Math/Vector2D.h"

// ======== Member Variables =========

vector<double> AccuracyAssessment2D_MultiBlock::LNorms = vector<double>(3,0.0);
double AccuracyAssessment2D_MultiBlock::TotalDomainArea = 0.0;
bool AccuracyAssessment2D_MultiBlock::AccuracyAssessed_Flag = false;
bool AccuracyAssessment2D_MultiBlock::Title_Error_Norms = true;
bool AccuracyAssessment2D_MultiBlock::Verbose = false;
unsigned int AccuracyAssessment2D_MultiBlock::TotalCells = 0;
vector<double> AccuracyAssessment2D_MultiBlock::Lift(0);
vector<double> AccuracyAssessment2D_MultiBlock::Drag(0);
vector<double> AccuracyAssessment2D_MultiBlock::SolidBodyLength(0);
double AccuracyAssessment2D_MultiBlock::FreeStreamDensity = 0.0;
double AccuracyAssessment2D_MultiBlock::FreeStreamVelocity = 0.0;

// ======== Implementation of Member Functions =========

double AccuracyAssessment2D_MultiBlock::getLift(void){
  double TotalLift(0);

  // sum up all contributions from each solid body
  for (int sb=0; sb<Lift.size(); ++sb){
    TotalLift += Lift[sb];
  }
  
  return TotalLift;
}

double AccuracyAssessment2D_MultiBlock::getLift(const int & BodyID){
  if (BodyID >= 0 && BodyID < Lift.size()){
    return Lift[BodyID];
  } else {
    throw runtime_error("AccuracyAssessment2D_MultiBlock::Lift(BodyID) ERROR! Index out of range!");
  }
}

double AccuracyAssessment2D_MultiBlock::getLiftCoefficient(void){
  double TotalLift(0), TotalWettedSurface(0);

  // sum up all contributions from each solid body
  for (int sb=0; sb<Lift.size(); ++sb){
    TotalLift += Lift[sb];
    TotalWettedSurface += SolidBodyLength[sb];
  }
  
  return TotalLift/(HALF*FreeStreamDensity*FreeStreamVelocity*FreeStreamVelocity*TotalWettedSurface);
}

double AccuracyAssessment2D_MultiBlock::getLiftCoefficient(const int & BodyID){
  if (BodyID >= 0 && BodyID < Lift.size()){
    return Lift[BodyID]/(HALF*FreeStreamDensity*FreeStreamVelocity*FreeStreamVelocity*SolidBodyLength[BodyID]); 
  } else {
    throw runtime_error("AccuracyAssessment2D_MultiBlock::LiftCoefficient(BodyID) ERROR! Index out of range!");
  }
}

double AccuracyAssessment2D_MultiBlock::getDrag(void){
  double TotalDrag(0);

  // sum up all contributions from each solid body
  for (int sb=0; sb<Drag.size(); ++sb){
    TotalDrag += Drag[sb];
  }
  
  return TotalDrag;
}

double AccuracyAssessment2D_MultiBlock::getDrag(const int & BodyID){
  if (BodyID >= 0 && BodyID < Drag.size()){
    return Drag[BodyID];
  } else {
    throw runtime_error("AccuracyAssessment2D_MultiBlock::Drag(BodyID) ERROR! Index out of range!");
  }
}

double AccuracyAssessment2D_MultiBlock::getDragCoefficient(void){
  double TotalDrag(0), TotalWettedSurface(0);

  // sum up all contributions from each solid body
  for (int sb=0; sb<Drag.size(); ++sb){
    TotalDrag += Drag[sb];
    TotalWettedSurface += SolidBodyLength[sb];
  }
  
  return TotalDrag/(HALF*FreeStreamDensity*FreeStreamVelocity*FreeStreamVelocity*TotalWettedSurface);
}

double AccuracyAssessment2D_MultiBlock::getDragCoefficient(const int & BodyID){
  if (BodyID >= 0 && BodyID < Drag.size()){
    return Drag[BodyID]/(HALF*FreeStreamDensity*FreeStreamVelocity*FreeStreamVelocity*SolidBodyLength[BodyID]);
  } else {
    throw runtime_error("AccuracyAssessment2D_MultiBlock::DragCoefficient(BodyID) ERROR! Index out of range!");
  }
}

double AccuracyAssessment2D_MultiBlock::getWettedSurface(void){
  double TotalWettedSurface(0);

  // sum up all contributions from each solid body
  for (int sb=0; sb<Drag.size(); ++sb){
    TotalWettedSurface += SolidBodyLength[sb];
  }

  return TotalWettedSurface;
}

double AccuracyAssessment2D_MultiBlock::getWettedSurface(const int & BodyID){
  if (BodyID >= 0 && BodyID < SolidBodyLength.size()){
    return SolidBodyLength[BodyID];
  } else {
    throw runtime_error("AccuracyAssessment2D_MultiBlock::WettedSurface(BodyID) ERROR! Index out of range!");
  }
}

/*! 
 * Output the aerodynamics data (i.e. lift, drag, etc.)
 * to the required output stream.
 * Output data for the whole configuration and for each 
 * solid body (i.e part) individually.
 */
void AccuracyAssessment2D_MultiBlock::PrintAerodynamicsData(ostream & os){

  // Output data for the whole configuration
  os << endl
     << " ==================================================================== "
     << endl
     << " Aerodynamics data for the complete configuration:"   << endl
     << " (Normalize coefficients using HALF*rho*U^2*WettedSurface!!!)"   << endl
     << "   Lift  = "  << setprecision(10) << getLift()  << endl
     << "   Drag  = "  << setprecision(10) << getDrag()  << endl
     << "   Cl  = "  << setprecision(10) << getLiftCoefficient()  << endl
     << "   Cd  = "  << setprecision(10) << getDragCoefficient()  << endl
     << "   WettedSurface  = "  << setprecision(10) << getWettedSurface()  << endl
     << " ==================================================================== "
     << endl;

  // Output data for each solid body individually
  if (Lift.size() > 1){
    for (int nb = 0; nb < Lift.size(); ++nb){
      os << endl
	 << " ==================================================================== "
	 << endl
	 << " Aerodynamics data for part #" << nb + 1 << " :"   << endl
	 << "   Lift  = "  << setprecision(10) << getLift(nb)  << endl
	 << "   Drag  = "  << setprecision(10) << getDrag(nb)  << endl
	 << "   Cl  = "  << setprecision(10) << getLiftCoefficient(nb)  << endl
	 << "   Cd  = "  << setprecision(10) << getDragCoefficient(nb)  << endl
	 << "   WettedSurface  = "  << setprecision(10) << getWettedSurface(nb)  << endl
	 << " ==================================================================== "
	 << endl;
    }
  }
}
