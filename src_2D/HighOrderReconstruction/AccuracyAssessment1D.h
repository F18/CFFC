/*!\file AccuracyAssessment1D.h
  \brief Header file implementing the templated AccuracyAssessment1D nothington class. */

#ifndef _ACCURACY_ASSESSMENT_1D_INCLUDED
#define _ACCURACY_ASSESSMENT_1D_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
using std::ostream;

/* Include CFFC header files */
// None

/*!
 * Class: AccuracyAssessment1D
 *
 * @brief Collection of subroutines that do post-processing operations related to accuracy assessment
 *
 */
class AccuracyAssessment1D{
  
public:
  template<class Soln_Block_Type, class Input_Parameters_Type>
  static void PrintErrorNorms(Soln_Block_Type * SolnBlk,
			      const Input_Parameters_Type & IP,
			      ostream & os);

  template<class Soln_Block_Type, class Input_Parameters_Type>
  static void OutputErrorNormsTecplot(Soln_Block_Type * SolnBlk,
				      const Input_Parameters_Type & IP);
  



private:
  // Private constructor
  AccuracyAssessment1D(void){};

  static vector<double> LNorms;
  static bool AccuracyAssessed_Flag;
  static bool Title_Error_Norms;
};

/*!
 * Print the error norms to the "os" output stream
 */
template<class Soln_Block_Type, class Input_Parameters_Type>
void AccuracyAssessment1D::PrintErrorNorms(Soln_Block_Type * SolnBlk,
					   const Input_Parameters_Type & IP,
					   ostream & os){
  
//       if(AccuracyAssessed_flag){
// 	std::cout << "\n The accuracy norms are:";
// 	PrintErrorNorms(LNorms,Input_Parameters);
//       }
//       else{
// 	if (Input_Parameters.ExactFunction == NULL){
// 	  // For problems without an exact solution
// 	  throw runtime_error("AssessSolutionAccuracy() ERROR! There is no exact solution specified for this problem.");
// 	} else {
// 	  // Analyze reconstruction to assess the accuracy
// 	  if (! batch_flag) std::cout << " \n  Assess the accuracy";
// 	  AssessSolutionAccuracy(Soln_ptr,LNorms,Input_Parameters);
// 	  if (! batch_flag) std::cout << " \n  Accuracy assessed\n";
// 	  // change flag
// 	  AccuracyAssessed_flag = true;
// 	  // print solution
// 	  if (! batch_flag) std::cout << "The accuracy norms of the reconstruction are:";
// 	  PrintErrorNorms(LNorms,Input_Parameters);
// 	}
//       }
 
  
}


/*!
 * Output the error norms to the output file in a tecplot format
 */
template<class Soln_Block_Type, class Input_Parameters_Type>
void AccuracyAssessment1D::OutputErrorNormsTecplot(Soln_Block_Type * SolnBlk,
						   const Input_Parameters_Type & IP){

//       if(AccuracyAssessed_flag){
// 	// output to file
// 	Output_Error_Norms_Tecplot(LNorms,Input_Parameters,Title_Error_Norms);
//       }
//       else {
// 	if (Input_Parameters.ExactFunction == NULL){
// 	  // For problems without an exact solution
// 	  throw runtime_error("AssessSolutionAccuracy() ERROR! There is no exact solution specified for this problem.");
// 	} else {
// 	  // Analyze reconstruction to assess the accuracy
// 	  if (! batch_flag) std::cout << " \n  Assess the accuracy:\n";
// 	  AssessSolutionAccuracy(Soln_ptr,LNorms,Input_Parameters);
// 	  if (! batch_flag) std::cout << " \n  Accuracy assessed\n";
// 	  // change flag
// 	  AccuracyAssessed_flag = true;
// 	  // output to file
// 	  Output_Error_Norms_Tecplot(LNorms,Input_Parameters,Title_Error_Norms);
// 	}
//       }

}



#endif //_ACCURACY_ASSESSMENT_1D_INCLUDED
