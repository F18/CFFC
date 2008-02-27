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

  //! @name Access to the error data:
  //@{
  static double & L1(void) {return LNorms[0]; }                 //!< return the L1 component of the error-norm vector
  static double & L2(void) {return LNorms[1]; }                 //!< return the L2 component of the error-norm vector
  static double & LMax(void) {return LNorms[2]; }               //!< return the LMax component of the error-norm vector
  static double & TotalArea(void) {return TotalDomainArea; }    //!< return the total area of the block
  static unsigned int & TotalNumberOfCells(void) {return TotalCells; }  //!< return the total number of cells 
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

  template<typename Input_Parameters_Type>
  static void SetVerboseFlag(const Input_Parameters_Type & IP);

  template<typename Quad_Soln_Block, typename Input_Parameters_Type>
  static int AssessSolutionAccuracyBasedOnExactSolution(Quad_Soln_Block *SolnBlk,
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
    // output error norms to the os stream
    if(Verbose){
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
    } // endif

    return 0;
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
  }


  // Assess solution accuracy based on the required mode for the current problem
  switch(IP.Accuracy_Assessment_Mode){
  case ACCURACY_ASSESSMENT_BASED_ON_EXACT_SOLUTION:
    if ( Verbose ){
      std::cout << "exact solution\n"
		<< "   Wait patiently ... ";
      std::cout.flush();
    }
    
    error_flag =  AssessSolutionAccuracyBasedOnExactSolution(SolnBlk,Soln_Block_List,IP);
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
	SolnBlk[nb].AssessAccuracy.ComputeSolutionErrors(IP.Accuracy_Assessment_Parameter,
							 IP.Accuracy_Assessment_Exact_Digits);
	
	L1() += SolnBlk[nb].AssessAccuracy.L1();
	L2() += SolnBlk[nb].AssessAccuracy.L2();
	LMax() = max(LMax(),SolnBlk[nb].AssessAccuracy.LMax());

	TotalDomainArea += SolnBlk[nb].AssessAccuracy.BlockArea();

	TotalCells += SolnBlk[nb].AssessAccuracy.UsedCells();
      }
    }

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

#endif	// _ACCURACY_ASSESSMENT_2D_MULTIBLOCK_INCLUDED
