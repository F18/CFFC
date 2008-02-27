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


vector<double> AccuracyAssessment2D_MultiBlock::LNorms = vector<double>(3,0.0);
double AccuracyAssessment2D_MultiBlock::TotalDomainArea = 0.0;
bool AccuracyAssessment2D_MultiBlock::AccuracyAssessed_Flag = false;
bool AccuracyAssessment2D_MultiBlock::Title_Error_Norms = true;
bool AccuracyAssessment2D_MultiBlock::Verbose = false;
unsigned int AccuracyAssessment2D_MultiBlock::TotalCells = 0;

