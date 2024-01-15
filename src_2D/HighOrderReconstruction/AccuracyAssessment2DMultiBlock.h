/*!\file AccuracyAssessment2DMultiBlock.h
  \brief Header file defining 2D accuracy assessment class for multi-block solutions. */

#ifndef _ACCURACY_ASSESSMENT_2D_MULTIBLOCK_INCLUDED
#define _ACCURACY_ASSESSMENT_2D_MULTIBLOCK_INCLUDED

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
using std::cout;

/* Include CFFC header files */
#include "AccuracyAssessment2D.h"
#include "../Utilities/Utilities.h"
#include "../AMR/AdaptiveBlock.h"
#include "../MPI/MPI.h"
#include "AccuracyAssessment_ExecutionMode.h"

#include "../System/System_Linux.h"

/*!
 * \class AccuracyAssessment2D_MultiBlock
 *
 * @brief Collection of routines for assessing the solution accuracy of 2D problems
 *        solved on multi-block quadrilateral grids.
 *
 */
class AccuracyAssessment2D_MultiBlock{
  
public:
  /*! @brief Assess the solution accuracy of the problem */
  template<typename Quad_Soln_Block, typename Input_Parameters_Type>
  static int AssessSolutionAccuracy(Quad_Soln_Block *SolnBlk,
				    const AdaptiveBlock2D_List &Soln_Block_List,
				    const Input_Parameters_Type &IP);

  /*! @brief Print the error norms to the specified output stream */
  template<typename Quad_Soln_Block, typename Input_Parameters_Type>
  static int PrintErrorNorms(Quad_Soln_Block *SolnBlk,
			     const AdaptiveBlock2D_List &Soln_Block_List,
			     const Input_Parameters_Type &IP,
			     ostream & os);
  
  /*! @brief Output the error norms to a file suitable for tecplot plotting */
  template<typename Quad_Soln_Block, typename Input_Parameters_Type>
  static void OutputErrorNormsTecplot(Quad_Soln_Block *SolnBlk,
				      const Input_Parameters_Type &IP);

  /*! @brief Print the error norms to the output file specified in the input parameters */
  template<typename Quad_Soln_Block, typename Input_Parameters_Type>
  static int WriteErrorNormsToOutputFile(Quad_Soln_Block *SolnBlk,
					 const AdaptiveBlock2D_List &Soln_Block_List,
					 const Input_Parameters_Type &IP);

  /*! @brief Append the error norms to the output file specified in the input parameters */
  template<typename Quad_Soln_Block, typename Input_Parameters_Type>
  static int AppendErrorNormsToOutputFile(Quad_Soln_Block *SolnBlk,
					  const AdaptiveBlock2D_List &Soln_Block_List,
					  const Input_Parameters_Type &IP);

  /*! @brief Append the error norms to the output file specified in the input parameters in a format suitable for Tecplot */
  template<typename Quad_Soln_Block, typename Input_Parameters_Type>
  static int AppendErrorNormsToTecplotOutputFile(Quad_Soln_Block *SolnBlk,
						 const AdaptiveBlock2D_List &Soln_Block_List,
						 const Input_Parameters_Type &IP,
						 const int &Number_of_Time_Steps,
						 const double &Time);

  template<typename Input_Parameters_Type>
  static void PrintTecplotTitle(const Input_Parameters_Type &IP,
				const int &Number_of_Time_Steps,
				const double &Time,
				ostream & os);

  //! @name Access to the error data:
  //@{
  static double & L1(void) {return LNorms[0]; }                 //!< return the L1 component of the error-norm vector
  static double & L2(void) {return LNorms[1]; }                 //!< return the L2 component of the error-norm vector
  static double & LMax(void) {return LNorms[2]; }               //!< return the LMax component of the error-norm vector
  static double & TotalArea(void) {return TotalDomainArea; }    //!< return the total area of the block
  static unsigned int & TotalNumberOfCells(void) {return TotalCells; }  //!< return the total number of cells 
  //@}

  //! @name Access to aerodynamic data:
  //@{
  static double getLift(void);	            //!< return total lift force produced by the current configuration
  static double getLift(const int & BodyID); //!< return lift force produced by the solid body BodyID
  static double getLiftCoefficient(void);    //!< return the lift coefficient characteristic for the current configuration
  static double getLiftCoefficient(const int & BodyID); //!< return the lift coefficient for the solid body BodyID
  static double getDrag(void);               //!< return total drag force produced by the current configuration
  static double getDrag(const int & BodyID); //!< return drag force produced by the solid body BodyID
  static double getDragCoefficient(void);    //!< return the drag coefficient characteristic for the current configuration
  static double getDragCoefficient(const int & BodyID); //!< return the drag coefficient for the solid body BodyID
  static double getWettedSurface(void);		       //!< return the total surface contributing to generate aerodynamic forces
  static double getWettedSurface(const int & BodyID);   //!< return the wetted surface of the solid body BodyID

  //! @brief Output the aerodynamics data (i.e. lift, drag, etc.) to the required stream
  static void PrintAerodynamicsData(ostream & os);
  //@}

  //! Prepare class for an accuracy reassessment
  template<typename Quad_Soln_Block>
  static void ResetForNewCalculation(Quad_Soln_Block *SolnBlk,
				     const AdaptiveBlock2D_List &Soln_Block_List);

private:
  AccuracyAssessment2D_MultiBlock(void);     //!< Private default constructor
  AccuracyAssessment2D_MultiBlock(const AccuracyAssessment2D_MultiBlock&); //!< Private copy constructor
  AccuracyAssessment2D_MultiBlock& operator=(const AccuracyAssessment2D_MultiBlock&); //!< Private assignment operator

  static vector<double> LNorms;	      //!< vector of error norms
  static double TotalDomainArea;      //!< total area of the computational domain
  static bool AccuracyAssessed_Flag;  //!< internal flag used to avoid re-assessment of an already determined error
  static bool Title_Error_Norms;      //!< internal flag used to ensure that the output of the Tecplot title is done only once
  static bool Verbose;		      //!< internal flag controlling the screen output stream
  static unsigned int TotalCells;     //!< total number of interior cells used for assessing the accuracy

  //! @name Member variables/functions related to lift and drag calculation:
  //@{
  static vector<double> Lift;	         //!< vector of lift forces. One entry for each solid body
  static vector<double> Drag;		 //!< vector of drag forces. One entry for each solid body
  static vector<double> SolidBodyLength; //!< vector of "wetted surfaces"
  static double FreeStreamDensity;       //!< the free stream density
  static double FreeStreamVelocity;	 //!< the free stream velocity
  //@}

  template<typename Input_Parameters_Type>
  static void SetVerboseFlag(const Input_Parameters_Type & IP);

  template<typename Quad_Soln_Block, typename Input_Parameters_Type>
  static int AssessSolutionAccuracyBasedOnExactSolution(Quad_Soln_Block *SolnBlk,
							const AdaptiveBlock2D_List &Soln_Block_List,
							const Input_Parameters_Type &IP);
  template<typename Quad_Soln_Block, typename Input_Parameters_Type>
  static int AssessSolutionAccuracyBasedOnEntropyVariation(Quad_Soln_Block *SolnBlk,
							   const AdaptiveBlock2D_List &Soln_Block_List,
							   const Input_Parameters_Type &IP);
  template<typename Quad_Soln_Block, typename Input_Parameters_Type>
  static int AssessSolutionAccuracyBasedOnLiftAndDragCoefficients(Quad_Soln_Block *SolnBlk,
								  const AdaptiveBlock2D_List &Soln_Block_List,
								  const Input_Parameters_Type &IP);

};

/*!
 * Set the verbose flag based on the input parameters
 */
template<typename Input_Parameters_Type> inline
void AccuracyAssessment2D_MultiBlock::SetVerboseFlag(const Input_Parameters_Type & IP){
  if ( IP.Verbose() ){
    Verbose = true;
  } else {
    Verbose = false;
  }
}

/*!
 * Reset the assessment accuracy flags in all blocks
 */
template<typename Quad_Soln_Block> inline
void AccuracyAssessment2D_MultiBlock::ResetForNewCalculation(Quad_Soln_Block *SolnBlk,
							     const AdaptiveBlock2D_List &Soln_Block_List){
  // reset the internal flag of this class
  AccuracyAssessed_Flag = false;

  // reset the flag of each block
  for (int nb = 0; nb < Soln_Block_List.Nblk; ++nb) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      SolnBlk[nb].AssessAccuracy.ResetForNewCalculation();
    }
  }
}

/*!
 * Print the error norms to the specified output stream
 */
template<typename Quad_Soln_Block, typename Input_Parameters_Type>
int AccuracyAssessment2D_MultiBlock::PrintErrorNorms(Quad_Soln_Block *SolnBlk,
						     const AdaptiveBlock2D_List &Soln_Block_List,
						     const Input_Parameters_Type &IP,
						     ostream & os){
  int error_flag(0);

  // Set verbose flag
  SetVerboseFlag(IP);

  // Assess the accuracy if it hasn't been assessed
  if(AccuracyAssessed_Flag == false){
    error_flag = AssessSolutionAccuracy(SolnBlk,
					Soln_Block_List,
					IP);
  }

  if (error_flag){
    return error_flag;
  } else {
    if( CFFC_Primary_MPI_Processor() && (os != cout || Verbose) ){
      
      // Customize the output based on the method
      switch(AccuracyAssessment_Execution_Mode::Method()){

      case AccuracyAssessment_Execution_Mode::Based_On_Lift_And_Drag_Coefficients:
	// output lift and drag coefficients
	PrintAerodynamicsData(os);
	break;

      default:
	// output error norms to the os stream
	os << endl
	   << " ==================================================================== "
	   << endl
	   << " The error norms for the solved problem are:"   << endl
	   << "   #Cells   = "  << TotalCells                  << endl
	   << "   L1_Norm  = "  << setprecision(10) << L1()    << endl
	   << "   L2_Norm  = "  << setprecision(10) << L2()    << endl
	   << "   Max_Norm = "  << setprecision(10) << LMax()  << endl
	   << " ==================================================================== "
	   << endl;
      }	// endswitch

    } // endif

    return 0;
  } // endif
}

/*!
 * Write the error norms to an output file which has the 
 * base name specified in the input parameters.
 * This function overrides any previous entries in the file.
 */
template<typename Quad_Soln_Block, typename Input_Parameters_Type>
int AccuracyAssessment2D_MultiBlock::WriteErrorNormsToOutputFile(Quad_Soln_Block *SolnBlk,
								 const AdaptiveBlock2D_List &Soln_Block_List,
								 const Input_Parameters_Type &IP){

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  ofstream output_file;    

  /* Determine prefix of output data file names. */

  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } /* endwhile */
  prefix[i] = '\0';

  /* Determine output data file name. */
  strcpy(extension, "_ErrorNorms.dat");
  strcpy(output_file_name, prefix);
  strcat(output_file_name, extension);

  /* Open the output data file. */
  output_file.open(output_file_name, ios::out);
  if (output_file.fail()) return (1);

  /* Write the error norms to the output stream. */
  PrintErrorNorms(SolnBlk,
		  Soln_Block_List,
		  IP,
		  output_file);

  /* Close the output data file. */
  output_file.close();

  /* Writing of output data files complete.  Return zero value. */
  return(0);
}

/*!
 * Append the error norms to an output file which has the 
 * base name specified in the input parameters.
 */
template<typename Quad_Soln_Block, typename Input_Parameters_Type>
int AccuracyAssessment2D_MultiBlock::AppendErrorNormsToOutputFile(Quad_Soln_Block *SolnBlk,
								  const AdaptiveBlock2D_List &Soln_Block_List,
								  const Input_Parameters_Type &IP){

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  ofstream output_file;    

  /* Determine prefix of output data file names. */

  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } /* endwhile */
  prefix[i] = '\0';

  /* Determine output data file name. */
  strcpy(extension, "_ErrorNorms.dat");
  strcpy(output_file_name, prefix);
  strcat(output_file_name, extension);

  /* Open the output data file. */
  output_file.open(output_file_name, ios::app);
  if (output_file.fail()) return (1);

  /* Write the error norms to the output stream. */
  PrintErrorNorms(SolnBlk,
		  Soln_Block_List,
		  IP,
		  output_file);

  /* Close the output data file. */
  output_file.close();

  /* Writing of output data files complete.  Return zero value. */
  return(0);

}

/*!
 * Append the error norms to an output file which has the 
 * base name specified in the input parameters, in a format
 * suitable for plotting with Tecplot.
 */
template<typename Quad_Soln_Block, typename Input_Parameters_Type>
int AccuracyAssessment2D_MultiBlock::AppendErrorNormsToTecplotOutputFile(Quad_Soln_Block *SolnBlk,
									 const AdaptiveBlock2D_List &Soln_Block_List,
									 const Input_Parameters_Type &IP,
									 const int &Number_of_Time_Steps,
									 const double &Time){

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  ofstream output_file;    
  int error_flag(0);

  /* Determine prefix of output data file names on primary CPU. */

  if( CFFC_Primary_MPI_Processor() ){
    i = 0;
    while (1) {
      if (IP.Output_File_Name[i] == ' ' ||
	  IP.Output_File_Name[i] == '.') break;
      prefix[i]=IP.Output_File_Name[i];
      i = i + 1;
      if (i > strlen(IP.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';

    /* Determine output data file name. */
    strcpy(extension, "_ErrorNorms.dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);

    // Check if the file has been previously generated
    if (System::Check_If_File_Exists(output_file_name) ){
      // Don't write the header in case there is a restart run
      Title_Error_Norms = false;
    }

    /* Open the output data file. */
    output_file.open(output_file_name, ios::app);
    if (output_file.fail()) return (1);

    /* Print Tecplot title if required */
    PrintTecplotTitle(IP,
		      Number_of_Time_Steps,
		      Time,
		      output_file);
  }

  /* Assess solution errors */
  error_flag = AssessSolutionAccuracy(SolnBlk,
				      Soln_Block_List,
				      IP);
  if (error_flag){
    return error_flag;
  }

  /* Print error measurements to the output data file on primary CPU. */
  if( CFFC_Primary_MPI_Processor() ){
    // Customize the output based on the method
    switch(AccuracyAssessment_Execution_Mode::Method()){
      
    case AccuracyAssessment_Execution_Mode::Based_On_Lift_And_Drag_Coefficients:
      output_file << " " 
		  << Time <<" "
		  << Number_of_Time_Steps <<" "
		  << getLift()  <<" " 
		  << getDrag()  <<" " 
		  << getLiftCoefficient()  <<" "
		  << getDragCoefficient()  <<" "
		  << getWettedSurface() << "\n";
      break;
      
    default:
      // output error norms to the os stream
      output_file << " " 
		  << TotalCells  <<" " 
		  << L1()  <<" " 
		  << L2()  <<" "
		  << LMax()  <<"\n";
    } // endswitch


    /* Close the output data file. */
    output_file.close();
  }

  /* Writing of output data files complete.  Return zero value. */
  return(0);

}

/*!
 * Print the Tecplot file header based on the error assessment method,
 * to the provided stream.
 */
template<typename Input_Parameters_Type>
void AccuracyAssessment2D_MultiBlock::PrintTecplotTitle(const Input_Parameters_Type &IP,
							const int &Number_of_Time_Steps,
							const double &Time,
							ostream & os){

  if (Title_Error_Norms){
    
    // Customize the header based on the method
    switch(AccuracyAssessment_Execution_Mode::Method()){

    case AccuracyAssessment_Execution_Mode::Based_On_Lift_And_Drag_Coefficients:
      os << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Solver, Aerodynamic Measurements \"\n"
	 << "VARIABLES = \"t\" \\ \n"
	 << "\"Time Step\" \\ \n"
	 << "\"Lift\" \\ \n"
	 << "\"Drag\" \\ \n"
	 << "\"Cl\" \\ \n"
	 << "\"Cd\" \\ \n"
	 << "\"WettedSurface\" \\ \n"
	 << "ZONE \n";
      break;

    default:
      os << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Error Measurements \"\n"
	 << "VARIABLES = \"Cells\" \\ \n"
	 << "\"L1_Norm\" \\ \n"
	 << "\"L2_Norm\" \\ \n"
	 << "\"Max_Norm\" \\ \n"
	 << "ZONE \n";
    } // endswitch

    // Mark the fact that title has been printed
    Title_Error_Norms = false;
      
  } // endif

}

/*!
 * Assess the solution accuracy of the problem
 */
template<typename Quad_Soln_Block, typename Input_Parameters_Type>
int AccuracyAssessment2D_MultiBlock::AssessSolutionAccuracy(Quad_Soln_Block *SolnBlk,
							    const AdaptiveBlock2D_List &Soln_Block_List,
							    const Input_Parameters_Type &IP){

  int error_flag(0);		//< used to exit correctly from MPI jobs
  
  // Set verbose flag
  SetVerboseFlag(IP);

  // Write message to user
  if ( Verbose ){
    std::cout << "\n Assess solution accuracy based on ";
    std::cout.flush();
  }


  // Assess solution accuracy based on the required mode for the current problem
  switch(AccuracyAssessment_Execution_Mode::Method()){
  case AccuracyAssessment_Execution_Mode::Based_On_Exact_Solution:
    if ( Verbose ){
      std::cout << "exact solution\n"
		<< "   Wait patiently ... ";
      std::cout.flush();
    }
    
    error_flag =  AssessSolutionAccuracyBasedOnExactSolution(SolnBlk,Soln_Block_List,IP);
    break;
    
  case AccuracyAssessment_Execution_Mode::Based_On_Entropy_Variation:
    if ( Verbose ){
      std::cout << "entropy variation relative to reference value\n"
		<< "   Wait patiently ... ";
      std::cout.flush();
    }

    error_flag =  AssessSolutionAccuracyBasedOnEntropyVariation(SolnBlk,Soln_Block_List,IP);
    break;

  case AccuracyAssessment_Execution_Mode::Based_On_Lift_And_Drag_Coefficients:
    if ( Verbose ){
      std::cout << "calculation of lift and drag coefficients\n"
		<< "   Wait patiently ... ";
      std::cout.flush();
    }

    error_flag = AssessSolutionAccuracyBasedOnLiftAndDragCoefficients(SolnBlk,Soln_Block_List,IP);
    break;

  default: 
    // Compute the error norms based on an exact solution
    if ( Verbose ){
      std::cout << "exact solution\n"
		<< "   Wait patiently ... ";
      std::cout.flush();
    }

    error_flag =  AssessSolutionAccuracyBasedOnExactSolution(SolnBlk,Soln_Block_List,IP);
  }


  // Write message to user
  if (error_flag){
    return error_flag;
  } else {
    // Accuracy assessed
    AccuracyAssessed_Flag = true;  

    if ( Verbose ){
      std::cout << "\n Solution accuracy assessed.\n";
      std::cout.flush();
    }
    
    return 0;
  }
}


/*!
 * Assess the solution accuracy of the problem using the provided exact solution
 */
template<typename Quad_Soln_Block, typename Input_Parameters_Type>
int AccuracyAssessment2D_MultiBlock::AssessSolutionAccuracyBasedOnExactSolution(Quad_Soln_Block *SolnBlk,
										const AdaptiveBlock2D_List &Soln_Block_List,
										const Input_Parameters_Type &IP){
  double TotalCells_MPI;

  try {

    // reset errors values
    L1() = ZERO; L2() = ZERO; LMax() = ZERO;
    TotalDomainArea = ZERO;
    TotalCells = 0;
    
    // compute the errors on each CPU for 
    // all the blocks on that CPU
    for (int nb = 0; nb < Soln_Block_List.Nblk; ++nb) {
      if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	if (IP.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
	  // Use the first high-order object to assess the accuracy
	  SolnBlk[nb].AssessAccuracy.ComputeSolutionErrorsHighOrder(AccuracyAssessment_Execution_Mode::Assessment_Parameter(),
								    AccuracyAssessment_Execution_Mode::Exact_Digits());
	} else {
	  // Use the low-order reconstruction
	  SolnBlk[nb].AssessAccuracy.ComputeSolutionErrors(AccuracyAssessment_Execution_Mode::Assessment_Parameter(),
							   AccuracyAssessment_Execution_Mode::Exact_Digits());
	}
	
	L1() += SolnBlk[nb].AssessAccuracy.L1();
	L2() += SolnBlk[nb].AssessAccuracy.L2();
	LMax() = max(LMax(),SolnBlk[nb].AssessAccuracy.LMax());

	TotalDomainArea += SolnBlk[nb].AssessAccuracy.BlockArea();

	TotalCells += SolnBlk[nb].AssessAccuracy.UsedCells();
      }
    }

    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  
    /* Total error for L1 norm on all CPUs */
    L1() = CFFC_Summation_MPI(L1());

    /* Total error for L2 norm on all CPUs */
    L2() = CFFC_Summation_MPI(L2());
  
    /* Final LMax norm on all CPUs */
    LMax() = CFFC_Maximum_MPI(LMax());

    /* Total area of the computational domain */
    TotalDomainArea = CFFC_Summation_MPI(TotalDomainArea);

    /* Total number of used cells for accuracy assessment */
    TotalCells_MPI = TotalCells;
    TotalCells_MPI = CFFC_Summation_MPI(TotalCells_MPI);
    TotalCells = (int)TotalCells_MPI;

    // === Final L1 error norm ===
    L1() /= TotalDomainArea;

    // === Final L2 error norm ===
    L2() = sqrt(L2()/TotalDomainArea);
    
    return 0;
  }
  catch (const ArgumentNullException & ){
    std::cerr << endl
	      << " ============================================================ " 
	      << endl
	      << " ERROR: The accuracy of the solution couldn't be determined!" << endl
	      << " The exact solution/accuracy parameters for this flow "       << endl
	      << " calculation are not defined!"                                << endl
	      << " Please check that the exact solution is set or "             << endl
	      << " don't require accuracy assessment!"                          << endl
	      << " ============================================================ " 
	      << endl;
    
    return 1;
  }

}

/*!
 * Assess the solution accuracy of the problem based on entropy variation
 * relative to a reference state.
 *
 * \note For solution states to which this method doesn't apply,
 *       the methods ComputeSolutionEntropyErrorsHighOrder() and
 *       ComputeSolutionEntropyErrors() must be specialized and
 *       throw ArgumentNullException error.
 */
template<typename Quad_Soln_Block, typename Input_Parameters_Type>
int AccuracyAssessment2D_MultiBlock::AssessSolutionAccuracyBasedOnEntropyVariation(Quad_Soln_Block *SolnBlk,
										   const AdaptiveBlock2D_List &Soln_Block_List,
										   const Input_Parameters_Type &IP){
  double TotalCells_MPI;

  try {

    // reset errors values
    L1() = ZERO; L2() = ZERO; LMax() = ZERO;
    TotalDomainArea = ZERO;
    TotalCells = 0;
    
    // compute the errors on each CPU for 
    // all the blocks on that CPU
    for (int nb = 0; nb < Soln_Block_List.Nblk; ++nb) {
      if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	if (IP.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
	  // Use the first high-order object to assess the accuracy
	  SolnBlk[nb].AssessAccuracy.ComputeSolutionEntropyErrorsHighOrder(AccuracyAssessment_Execution_Mode::Exact_Digits(),
									   IP);
	} else {
	  // Use the low-order reconstruction
	  SolnBlk[nb].AssessAccuracy.ComputeSolutionEntropyErrors(AccuracyAssessment_Execution_Mode::Exact_Digits(),
								  IP);
	}
	
	L1() += SolnBlk[nb].AssessAccuracy.L1();
	L2() += SolnBlk[nb].AssessAccuracy.L2();
	LMax() = max(LMax(),SolnBlk[nb].AssessAccuracy.LMax());

	TotalDomainArea += SolnBlk[nb].AssessAccuracy.BlockArea();

	TotalCells += SolnBlk[nb].AssessAccuracy.UsedCells();
      }
    }

    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

    /* Total error for L1 norm on all CPUs */
    L1() = CFFC_Summation_MPI(L1());

    /* Total error for L2 norm on all CPUs */
    L2() = CFFC_Summation_MPI(L2());
  
    /* Final LMax norm on all CPUs */
    LMax() = CFFC_Maximum_MPI(LMax());

    /* Total area of the computational domain */
    TotalDomainArea = CFFC_Summation_MPI(TotalDomainArea);

    /* Total number of used cells for accuracy assessment */
    TotalCells_MPI = TotalCells;
    TotalCells_MPI = CFFC_Summation_MPI(TotalCells_MPI);
    TotalCells = (int)TotalCells_MPI;

    // === Final L1 error norm ===
    L1() /= TotalDomainArea;

    // === Final L2 error norm ===
    L2() = sqrt(L2()/TotalDomainArea);
    
    return 0;
  }
  catch (const ArgumentNullException & ){
    std::cerr << endl
	      << " ================================================================== " 
	      << endl
	      << " ERROR: The accuracy of the solution couldn't be determined!" << endl
	      << " The entropy variation criterion couldn't be used for this flow! " << endl
	      << " Please check the problem or don't require accuracy assessment!" << endl
	      << " ================================================================== " 
	      << endl;
    
    return 1;
  }

}


/*!
 * Assess the solution accuracy by calculating 
 * the lift and drag forces/coefficients.
 */
template<typename Quad_Soln_Block, typename Input_Parameters_Type>
int AccuracyAssessment2D_MultiBlock::
AssessSolutionAccuracyBasedOnLiftAndDragCoefficients(Quad_Soln_Block *SolnBlk,
						     const AdaptiveBlock2D_List &Soln_Block_List,
						     const Input_Parameters_Type &IP){

  try {
    
    // === Set variables for lift and drag calculation
    int NumberOfSolidBodies(Quad_Soln_Block::GridType::BndSplineType::NumberOfSolidBodies());
    Lift.clear(); 
    Lift.assign(NumberOfSolidBodies, 0.0);
    Drag.clear(); 
    Drag.assign(NumberOfSolidBodies, 0.0);
    SolidBodyLength.clear();
    SolidBodyLength.assign(NumberOfSolidBodies, 0.0);
    FreeStreamDensity = IP.FreeStreamDensity();
    FreeStreamVelocity = IP.FreeStreamVelocity();
    double alpha;			 // angle of attack
    alpha = TWO*PI*IP.Flow_Angle/360.00; // angle of attack in radians
    double cos_alpha(cos(alpha)), sin_alpha(sin(alpha)); // cosine and sine of the angle of attach
    double Fx, Fy;

    // Counters
    int nb;
    
    // Compute the contributions to aerodynamic forces (i.e. Fx and Fy)
    // on each solid body due to each block on the current CPU.
    // Lift variable stores Fx and Drag variable stores Fy
    for (nb = 0; nb < Soln_Block_List.Nblk; ++nb) {
      if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED && SolnBlk[nb].Grid.IsThereAnySolidBoundary() ) {
	
	if (IP.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
	  // Use the first high-order object to assess the accuracy
	  SolnBlk[nb].AssessAccuracy.addAerodynamicForcesHighOrder(Lift,Drag,SolidBodyLength);
	} else {
	  // Use the low-order reconstruction
	  SolnBlk[nb].AssessAccuracy.addAerodynamicForces(Lift,Drag,SolidBodyLength);
	}
	
      }	//endif
    } //endfor
   
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

    for (nb = 0; nb < NumberOfSolidBodies; ++nb){
      /* Total aerodynamic force Fx on all CPUs for each of the solid bodies */
      Lift[nb] = CFFC_Summation_MPI(Lift[nb]);
      
      /* Total aerodynamic force Fy on all CPUs for each of the solid bodies */
      Drag[nb] = CFFC_Summation_MPI(Drag[nb]);
            
      /* Total wetted surface on all CPUs for each of the solid bodies */
      SolidBodyLength[nb] = CFFC_Summation_MPI(SolidBodyLength[nb]);
    }

    // Calculate final lift and drag forces for each solid body
    for (nb = 0; nb < NumberOfSolidBodies; ++nb){
      // Current aerodynamic forces in 'x' and 'y' directions
      Fx = Lift[nb];
      Fy = Drag[nb];
      
      // Lift force for the current solid body
      Lift[nb] = Fy*cos_alpha - Fx*sin_alpha;

      // Drag force for the current solid body
      Drag[nb] = Fy*sin_alpha + Fx*cos_alpha;
    }

    return 0;
  }
  catch (const ArgumentNullException & ){
    std::cerr << endl
	      << " ================================================================== " 
	      << endl
	      << " ERROR: The lift and drag coefficients couldn't be determined!" << endl
	      << " Please check the problem or don't require accuracy assessment!" << endl
	      << " ================================================================== " 
	      << endl;
    
    return 1;
  }
}

#endif	// _ACCURACY_ASSESSMENT_2D_MULTIBLOCK_INCLUDED
