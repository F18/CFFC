/*!\file AccuracyAssessment1D.cc
  \brief Source file implementing/setting member functions/variables prototyped in AccuracyAssessment1D class. */

/* Include CFFC header files */
#include "AccuracyAssessment1D.h"

vector<double> AccuracyAssessment1D::LNorms = vector<double>(3);
bool AccuracyAssessment1D::AccuracyAssessed_Flag = false;
bool AccuracyAssessment1D::Title_Error_Norms = true;
