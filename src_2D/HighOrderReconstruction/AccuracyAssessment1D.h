/*!\file AccuracyAssessment1D.h
  \brief Header file implementing the templated AccuracyAssessment1D nothington class. */

#ifndef _ACCURACY_ASSESSMENT_1D_INCLUDED
#define _ACCURACY_ASSESSMENT_1D_INCLUDED

/* Include required C++ libraries. */
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

/* Using std namespace functions */
using std::ostream;
using std::vector;
using std::endl;
using std::runtime_error;

/* Include CFFC header files */
#include "../CFD/CFD.h"
#include "../Utilities/Utilities.h"
#include "../Math/NumericalLibrary.h"

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
  

  template<class Soln_Block_Type, class Input_Parameters_Type>
  static void AssessSolutionAccuracy(Soln_Block_Type * SolnBlk,
				     const Input_Parameters_Type & IP,
				     typename Soln_Block_Type::HighOrderType & 
				     (Soln_Block_Type::*AccessToHighOrderVar)(void) = 
				     &Soln_Block_Type::CellHighOrder,
				     double (Soln_Block_Type::*ComputeLowOrderSolutionAt)(const double &,
											  const unsigned) = 
				     &Soln_Block_Type::SolutionAtCoordinates_PWL);

  template<class Soln_Block_Type, class Function_Object_Type>
  static void ComputeSolutionErrorsHighOrder(Soln_Block_Type * SolnBlk,
					     Function_Object_Type FuncObj,
					     const unsigned parameter,
					     typename Soln_Block_Type::HighOrderType & 
					     (Soln_Block_Type::*AccessToHighOrderVar)(void) = 
					     &Soln_Block_Type::CellHighOrder);

  template<class Soln_Block_Type, class Function_Object_Type>
  static void ComputeSolutionErrorsLowOrder(Soln_Block_Type * SolnBlk,
					    Function_Object_Type FuncObj,
					    const unsigned parameter,
					    double (Soln_Block_Type::*ComputeLowOrderSolutionAt)(const double &,
												 const unsigned) = 
					    &Soln_Block_Type::SolutionAtCoordinates_PWL);
  
  static double & L1(void) {return LNorms[0]; }	//!< return the L1 component of the error-norm vector
  static double & L2(void) {return LNorms[1]; } //!< return the L2 component of the error-norm vector
  static double & LMax(void) {return LNorms[2]; } //!< returns the LMax component of the error-norm vector

  static void ResetForNewCalculation(void){ AccuracyAssessed_Flag = false; }

private:
  AccuracyAssessment1D(void);     //!< Private default constructor
  AccuracyAssessment1D(const AccuracyAssessment1D&); //!< Private copy constructor
  AccuracyAssessment1D& operator=(const AccuracyAssessment1D&); //!< Private assignment operator


  static vector<double> LNorms;	      //!< vector of error norms
  static bool AccuracyAssessed_Flag;  //!< internal flag used to avoid re-assessment of an already determined error
  static bool Title_Error_Norms;      //!< internal flag used to ensure that the output of the Tecplot title is done only once
  static bool Verbose;		      //!< internal flag controlling the screen output stream

  // Operating functions
  /*! @brief Integrate over the domain of a 1D cell */
  template<class Soln_Block_Type, class Function_Object_Type, class ReturnType>
  static ReturnType IntegrateOverTheCell(Soln_Block_Type *SolnBlk, const int iCell,
					 const Function_Object_Type FuncObj, 
					 const int & digits, ReturnType _dummy_param);

  /*! @brief Compute the L1 error norm using a low-order variable */
  template<class Soln_Block_Type, class Function_Object_Type>
  static double ComputeSolutionErrorL1(Soln_Block_Type *SolnBlk, const int iCell,
				       const Function_Object_Type FuncObj,
				       const unsigned parameter,
				       double (Soln_Block_Type::*ComputeLowOrderSolutionAt)(const double &,
											    const unsigned) = 
				       &Soln_Block_Type::SolutionAtCoordinates_PWL);

  /*! @brief Compute the L2 error norm using a low-order variable */
//using the computational domain SolnBlk and the ComputeLowOrderSolutionAt member function */
  template<class Soln_Block_Type, class Function_Object_Type>
  static double ComputeSolutionErrorL2(Soln_Block_Type *SolnBlk, const int iCell,
				       const Function_Object_Type FuncObj,
				       const unsigned parameter,
				       double (Soln_Block_Type::*ComputeLowOrderSolutionAt)(const double &,
											    const unsigned) = 
				       &Soln_Block_Type::SolutionAtCoordinates_PWL);
  
  template<class Input_Parameters_Type>
  static void SetVerboseFlag(const Input_Parameters_Type & IP);

  static void OutputProgress(const int & i);


};

/*!
 * Set the verbose flag based on the input parameters
 */
template<class Input_Parameters_Type> inline
void AccuracyAssessment1D::SetVerboseFlag(const Input_Parameters_Type & IP){
  if ( IP.Verbose() ){
    Verbose = true;
  } else {
    Verbose = false;
  }
}

/*!
 * Mark the finish of the error computation in a cell by printing a 
 * dot on the screen.
 */
inline void AccuracyAssessment1D::OutputProgress(const int & i){
  if (Verbose) {
    if (i - 55*(i/55) == 0 ){
      cout << "\n ." ;
    } else {
      cout << "." ;
    }
  }
}


/*!
 * Print the error norms to the "os" output stream
 */
template<class Soln_Block_Type, class Input_Parameters_Type>
void AccuracyAssessment1D::PrintErrorNorms(Soln_Block_Type * SolnBlk,
					   const Input_Parameters_Type & IP,
					   ostream & os){
  
  // Set verbose flag
  SetVerboseFlag(IP);

  if(AccuracyAssessed_Flag == false){
    AssessSolutionAccuracy(SolnBlk,IP);
  }

  // print solution
  if (Verbose){
    os << "\n The accuracy norms based on parameter " << IP.ErrorParameter << " are:\n";
    os << " L1 = " << L1() << std::endl
       << " L2 = " << L2() << std::endl
       << " LMax = " << LMax() << std::endl;
    os.flush();
  }
}


/*!
 * Output the error norms to the output file in a tecplot format
 */
template<class Soln_Block_Type, class Input_Parameters_Type>
void AccuracyAssessment1D::OutputErrorNormsTecplot(Soln_Block_Type * SolnBlk,
						   const Input_Parameters_Type & IP){

  int i;
  char prefix[256], extension[256], output_file_name[256];
  char order[30];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Set verbose flag
  SetVerboseFlag(IP);
  
  if(AccuracyAssessed_Flag == false){
    AssessSolutionAccuracy(SolnBlk,IP);
  }

  /* Determine prefix of output data file names. */
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > (int)strlen(IP.Output_File_Name) ) break;
  } /* endwhile */
  prefix[i] = '\0';

  /* Determine output data file name. */
  switch(IP.Space_Accuracy){
  case 1:
    strcpy(order,"1");
    break;
  case 2:
    strcpy(order,"2");
    break;
  case 3:
    strcpy(order,"3");
    break;
  case 4:
    strcpy(order,"4");
    break;
  case 5:
    strcpy(order,"5");
    break;
  case 6:
    strcpy(order,"6");
    break;
  case 7:
    strcpy(order,"7");
    break;
  case 8:
    strcpy(order,"8");
    break;
  default:
    strcpy(order,"0");
  }

  /* Determine output data file name. */
  strcpy(extension,"_AccuracyNorms_Order");
  strcat(extension, order);
  strcpy(output_file_name, prefix); 
  strcat(output_file_name, extension);
  strcat(output_file_name, ".dat");
  output_file_name_ptr = output_file_name;
  /* Open the output data file. */

  if (Title_Error_Norms){
    output_file.open(output_file_name_ptr, ios::trunc); 
    if (output_file.fail()) exit(1);
    // write the header
    output_file << "TITLE = \""	<< " Accuracy Norms for " << IP.ICs_Type
		<< " solved for Tmax=" << IP.Time_Max << ", on the domain: "
		<< IP.X_Min << " :  " << IP.X_Max
		<<"\""<< endl
		<< "VARIABLES = \"N<sub>Nodes</sub>\" \n "
		<< "\" <math>o</math><greek>E</greek><math>o</math><sub><greek>1</greek></sub> \" \n"
		<< "\" <math>o</math><greek>E</greek><math>o</math><sub><greek>2</greek></sub> \" \n"
		<< "\" <math>o</math><greek>E</greek><math>o</math><sub><math>%</math></sub> \" \n" ;

    output_file	<< "ZONE \n";
    // Make sure that the title is written only once}
    Title_Error_Norms = false;

  } else{
    output_file.open(output_file_name_ptr, ios::app);
    if (output_file.fail()) exit(1);
  } //endif(Title_Error_Norms)

  /* Write the solution data. */ 
  cout << "\n Writing error norms to output data file `"
       << output_file_name_ptr << "'.\n";
  output_file << setprecision(14);
  output_file << " " << IP.Number_of_Nodes
    /* Obs. The mesh size "DeltaX" is proportional with 1.0/#Nodes */
	      << "   " << L1()
	      << "   " << L2()
	      << "   " << LMax()
	      << "\n";

  output_file.unsetf(ios::scientific);
  output_file << setprecision(6);

  /* Close the output data file. */

  output_file.close();

  /* Writing of output data files complete. */
}

/*!
 * Compute all error norms.
 * Based on the reconstruction type the routine usess the high-order or the low order variables.
 */
template<class Soln_Block_Type, class Input_Parameters_Type>
void AccuracyAssessment1D::AssessSolutionAccuracy(Soln_Block_Type * SolnBlk,
						  const Input_Parameters_Type & IP,
						  typename Soln_Block_Type::HighOrderType & 
						  (Soln_Block_Type::*AccessToHighOrderVar)(void),
						  double (Soln_Block_Type::*ComputeLowOrderSolutionAt)(const double &,
												       const unsigned) ){

  double _dummy_var;

  // Set verbose flag
  SetVerboseFlag(IP);

  if (IP.ExactFunction == NULL){
    std::cerr << endl
	      << " ============================================================ " 
	      << endl
	      << " ERROR: The accuracy of the solution couldn't be determined!"
	      << endl
	      << " The exact solution is invalid!" << endl
	      << " Please check the pointer to the exact solution or don't use "
	      << endl
	      << "the \"Print_Accuracy\" command!"
	      << endl
	      << " ============================================================ " 
	      << endl;
    
    throw runtime_error("AssessSolutionAccuracy() ERROR! There is no exact solution specified for this problem.");
  }


  // Analyze reconstruction to assess the accuracy
  if ( Verbose) std::cout << " \n Assess solution accuracy:\n ";


  if (IP.i_Reconstruction != RECONSTRUCTION_HIGH_ORDER){
    // == Compute errors using ComputeLowOrderSolutionAt
    ComputeSolutionErrorsLowOrder(SolnBlk,
				  mapped_function(IP.ExactFunction,_dummy_var,IP.X_Min,IP.X_Max,
						  IP.X_ExactSolution_Min,IP.X_ExactSolution_Max),
				  IP.ErrorParameter,
				  ComputeLowOrderSolutionAt);
  } else {
    // == Compute errors using the high-order variable
    ComputeSolutionErrorsHighOrder(SolnBlk,
				   mapped_function(IP.ExactFunction,_dummy_var,IP.X_Min,IP.X_Max,
						   IP.X_ExactSolution_Min,IP.X_ExactSolution_Max),
				   IP.ErrorParameter,
				   AccessToHighOrderVar);
  }

  // Accuracy assessed
  AccuracyAssessed_Flag = true;  

  if ( Verbose ) std::cout << " \n Solution accuracy assessed.\n";

}


/*!
 * Compute all error norms using the high-order variable provided by 'AccessToHighOrderVar' pointer to member function.
 */
template<class Soln_Block_Type, class Function_Object_Type>
void AccuracyAssessment1D::ComputeSolutionErrorsHighOrder(Soln_Block_Type * SolnBlk,
							  Function_Object_Type FuncObj,
							  const unsigned parameter,
							  typename Soln_Block_Type::HighOrderType & 
							  (Soln_Block_Type::*AccessToHighOrderVar)(void) ){

  int ICl(SolnBlk[0].ICl), ICu(SolnBlk[0].ICu);
  double DomainLength(0.0);
  double l1_norm(0.0);

  // Reset the error norm values
  L1() = 0.0; L2() = 0.0; LMax() = 0.0;

  for (int i=ICl; i<=ICu; ++i){
    // Calculate l1_norm for each cell
    l1_norm = (SolnBlk[i].*AccessToHighOrderVar)().ComputeSolutionErrorL1(FuncObj,parameter);

    // Add the error contribution to L1 from the current cell
    L1() += l1_norm;

    // Add the error contribution to L2 from the current cell
    L2() += (SolnBlk[i].*AccessToHighOrderVar)().ComputeSolutionErrorL2(FuncObj,parameter);

    // Compute LMax norm
    LMax() = max(LMax(),l1_norm/(SolnBlk[i].*AccessToHighOrderVar)().CellDelta());

    // Add the domain size contribution of the current cell
    DomainLength += (SolnBlk[i].*AccessToHighOrderVar)().CellDelta();

    // Output progress 
    OutputProgress(i);
  }

  /* compute final expression for norms */
  L1() /= DomainLength;
  L2() = sqrt(L2()/DomainLength);
}

//! Integrate over the cell
template<class Soln_Block_Type, class Function_Object_Type, class ReturnType>
ReturnType AccuracyAssessment1D::IntegrateOverTheCell(Soln_Block_Type *SolnBlk, const int iCell,
						      const Function_Object_Type FuncObj, 
						      const int & digits, ReturnType _dummy_param){

  return AdaptiveGaussianQuadrature(FuncObj,
				    SolnBlk[iCell].CellCenter() - 0.5* SolnBlk[iCell].CellDelta(),
				    SolnBlk[iCell].CellCenter() + 0.5* SolnBlk[iCell].CellDelta(),
				    _dummy_param,digits);
}

/*! 
 * \param SolnBlk computational domain used to determine the error norm
 * \param FuncObj function object relative to which the error is evaluated
 * \param parameter the solution state variable used to evaluate the error norm
 * \param ComputeLowOrderSolutionAt Soln_Block_Type member function which provides the value of
 *        the reconstructed polynomial at a particular location.
 * \return the L1 error norm (i.e. integral of the error function over the whole domain divided by the length of the domain)
 */
template<class Soln_Block_Type, class Function_Object_Type>
double AccuracyAssessment1D::ComputeSolutionErrorL1(Soln_Block_Type *SolnBlk, const int iCell,
						    const Function_Object_Type FuncObj,
						    const unsigned parameter,
						    double (Soln_Block_Type::
							    *ComputeLowOrderSolutionAt)(const double &,
											const unsigned)){
  // Set the type of the returned value
  double _dummy_param(0.0);
  
  // Call the integration function
  return IntegrateOverTheCell(SolnBlk, iCell,
			      error_function(FuncObj,
					     wrapped_member_function_one_parameter(&SolnBlk[iCell],
										   ComputeLowOrderSolutionAt,
										   parameter,
										   _dummy_param),
					     _dummy_param),
			      10,_dummy_param);
}

/*! 
 * \param SolnBlk computational domain used to determine the error norm
 * \param FuncObj function object relative to which the error is evaluated
 * \param parameter the solution state variable used to evaluate the error norm
 * \param ComputeLowOrderSolutionAt Soln_Block_Type member function which provides the value of
 *        the reconstructed polynomial at a particular location.
 * \return the L2 error norm (i.e. square root of the integral of the squared error function over 
 *         the whole domain divided by the length of the domain)
 */
template<class Soln_Block_Type, class Function_Object_Type>
double AccuracyAssessment1D::ComputeSolutionErrorL2(Soln_Block_Type *SolnBlk, const int iCell,
						    const Function_Object_Type FuncObj,
						    const unsigned parameter,
						    double (Soln_Block_Type::*
							    ComputeLowOrderSolutionAt)(const double &,
										       const unsigned) ){
  // Set the type of the returned value
  double _dummy_param(0.0);
  
  // Call the integration function
  return IntegrateOverTheCell(SolnBlk, iCell,
			      square_error_function(FuncObj,
						    wrapped_member_function_one_parameter(&SolnBlk[iCell],
											  ComputeLowOrderSolutionAt,
											  parameter,
											  _dummy_param),
						    _dummy_param),
			      10,_dummy_param);
}


/*!
 * Compute all error norms using the member function ComputeLowOrderSolutionAt which evaluates the 
 * reconstructed polynomial at a particular location.
 */
template<class Soln_Block_Type, class Function_Object_Type>
void AccuracyAssessment1D::ComputeSolutionErrorsLowOrder(Soln_Block_Type * SolnBlk,
							 Function_Object_Type FuncObj,
							 const unsigned parameter,
							 double (Soln_Block_Type::
								 *ComputeLowOrderSolutionAt)(const double &,
											     const unsigned) ){

  int ICl(SolnBlk[0].ICl), ICu(SolnBlk[0].ICu);
  double DomainLength(0.0);
  double l1_norm(0.0);

  // Reset the error norm values
  L1() = 0.0; L2() = 0.0; LMax() = 0.0;

  for (int i=ICl; i<=ICu; ++i){
    // Calculate l1_norm for each cell
    l1_norm = ComputeSolutionErrorL1(SolnBlk,i,FuncObj,parameter,ComputeLowOrderSolutionAt);

    // Add the error contribution to L1 from the current cell
    L1() += l1_norm;

    // Add the error contribution to L2 from the current cell
    L2() += ComputeSolutionErrorL2(SolnBlk,i,FuncObj,parameter,ComputeLowOrderSolutionAt);

    // Compute LMax norm
    LMax() = max(LMax(),l1_norm/SolnBlk[i].CellDelta());

    // Add the domain size contribution of the current cell
    DomainLength += SolnBlk[i].CellDelta();

    // Output progress 
    OutputProgress(i);
  }

  /* compute final expression for norms */
  L1() /= DomainLength;
  L2() = sqrt(L2()/DomainLength);
}

#endif //_ACCURACY_ASSESSMENT_1D_INCLUDED
