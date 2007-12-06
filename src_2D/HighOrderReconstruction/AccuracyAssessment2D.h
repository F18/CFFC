/*!\file AccuracyAssessment2D.h
  \brief Header file defining 2D accuracy assessment class. */

#ifndef _ACCURACY_ASSESSMENT_2D_INCLUDED
#define _ACCURACY_ASSESSMENT_2D_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
// None


/*!
 * \class AccuracyAssessment2D
 *
 * @brief Collection of routines for assessing the solution accuracy of 2D problems
 *
 */
template<typename Quad_Soln_Block>
class AccuracyAssessment2D{
public:

  //! Constructor with solution block
  AccuracyAssessment2D(Quad_Soln_Block * AssociatedSolnBlock);

  //! Destructor
  ~AccuracyAssessment2D(void){ };

  //! Re-associate solution block pointer
  void AssociateSolutionBlock(Quad_Soln_Block * AssociatedSolnBlock){ SolnBlk = AssociatedSolnBlock; }

  //! Access to the solution block
  Quad_Soln_Block *getGrid(void) const { return SolnBlk; }

  template<class Input_Parameters_Type>
  void PrintErrorNorms(const Input_Parameters_Type & IP,
		       ostream & os);
  
  template<class Input_Parameters_Type>
  void OutputErrorNormsTecplot(const Input_Parameters_Type & IP);
  

  //   template<class Input_Parameters_Type>
  //   void AssessSolutionAccuracy(const Input_Parameters_Type & IP,
  // 			      typename Quad_Soln_Block::HighOrderType & 
  // 			      (Quad_Soln_Block::*AccessToHighOrderVar)(void) = 
  // 			      &Quad_Soln_Block::CellHighOrder,
  // 			      double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
  // 									   const Vector2D &,
  // 									   const unsigned int &) const = 
  // 			      &Quad_Soln_Block::PiecewiseLinearSolutionAtLocation);
  
  //   template<class Quad_Soln_Block, class Function_Object_Type>
  //   static void ComputeSolutionErrorsHighOrder(Quad_Soln_Block * SolnBlk,
  // 					     Function_Object_Type FuncObj,
  // 					     const unsigned parameter,
  // 					     typename Quad_Soln_Block::HighOrderType & 
  // 					     (Quad_Soln_Block::*AccessToHighOrderVar)(void) = 
  // 					     &Quad_Soln_Block::CellHighOrder);

  void ComputeSolutionErrors(const unsigned int &parameter,
			     double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
									  const Vector2D &,
									  const unsigned int &) const =
			     &Quad_Soln_Block::PiecewiseLinearSolutionAtLocation);
  //     LNorms[0] = fabs(SolnBlk->ExactSoln->Solution(0.0,0.0) - 
  // 		     (SolnBlk->*ComputeLowOrderSolutionAt)(0,0,Vector2D(0.0,0.0),parameter));
  
  void ComputeSolutionErrors(double (*FuncObj)(const double &, const double &),
			     const unsigned int &parameter,
			     double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
									  const Vector2D &,
									  const unsigned int &) const =
			     &Quad_Soln_Block::PiecewiseLinearSolutionAtLocation){ };
  
  double & L1(void) {return LNorms[0]; }     //!< return the L1 component of the error-norm vector
  double & L2(void) {return LNorms[1]; }     //!< return the L2 component of the error-norm vector
  double & LMax(void) {return LNorms[2]; }   //!< returns the LMax component of the error-norm vector
  
  void ResetForNewCalculation(void){ AccuracyAssessed_Flag = false; }

private:
  Quad_Soln_Block *SolnBlk;	//!< Pointer to the solution block associated to this object

  AccuracyAssessment2D(void);	//!< Private default constructor  

  vector<double> LNorms;	//!< vector of errors/error norms
  double TotalBlockAria;	//!< the total area of the block
  bool AccuracyAssessed_Flag;   //!< internal flag used to avoid re-assessment of an already determined error
  bool Title_Error_Norms;       //!< internal flag used to ensure that the output of the Tecplot title is done only once
  bool Verbose;		        //!< internal flag controlling the screen output stream
  
  // Operating functions
  /*! @brief Compute errors for L1 norm calculation using the ComputeSolutionAt solution block member function */
  double ComputeSolutionErrorL1(const int &iCell, const int &jCell,
				double (*FuncObj)(const double &, const double &),
 				const unsigned int &parameter,
 				double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
									     const Vector2D &,
 									     const unsigned int &) const = 
 				&Quad_Soln_Block::PiecewiseLinearSolutionAtLocation){ };

  /*! @brief Compute errors for L2 norm calculation using the ComputeSolutionAt solution block member function */
  double ComputeSolutionErrorL2(const int &iCell, const int &jCell,
				double (*FuncObj)(const double &, const double &),
 				const unsigned int &parameter,
 				double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
									     const Vector2D &,
 									     const unsigned int &) const = 
 				&Quad_Soln_Block::PiecewiseLinearSolutionAtLocation){ };
  
  template<class Input_Parameters_Type>
  void SetVerboseFlag(const Input_Parameters_Type & IP);

  void OutputProgress(const int & i);
};


/*!
 * Constructor the accuracy assessment object and set the link to the associated solution block
 */
template<typename Quad_Soln_Block> inline 
AccuracyAssessment2D<Quad_Soln_Block>::AccuracyAssessment2D(Quad_Soln_Block * AssociatedSolnBlock):
  AccuracyAssessed_Flag(false), Title_Error_Norms(true), Verbose(false) 
{
  LNorms = vector<double>(3,0.0);
  SolnBlk = AssociatedSolnBlock;
}

/*!
 * Set the verbose flag based on the input parameters
 */
template<typename Quad_Soln_Block>
template<typename Input_Parameters_Type> inline
void AccuracyAssessment2D<Quad_Soln_Block>::SetVerboseFlag(const Input_Parameters_Type & IP){
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
template<typename Quad_Soln_Block> inline 
void AccuracyAssessment2D<Quad_Soln_Block>::OutputProgress(const int & i){
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
template<typename Quad_Soln_Block>
template<typename Input_Parameters_Type> inline
void AccuracyAssessment2D<Quad_Soln_Block>::PrintErrorNorms(const Input_Parameters_Type & IP,
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
template<typename Quad_Soln_Block>
template<typename Input_Parameters_Type>
void AccuracyAssessment2D<Quad_Soln_Block>::OutputErrorNormsTecplot(const Input_Parameters_Type & IP){

}

/*!
 * Compute solution errors relative to the exact solution 
 * set in the Quad_Soln_Block.
 *
 * \param parameter the state variable which is used for computing the errors
 * \param ComputeSolutionAt Quad_Soln_Block member function which returns the solution 
 *                          at a given location (PointOfInterest) using the reconstruction 
 *                          of cell (i,j) for a specified parameter.
 */
template<typename Quad_Soln_Block>
void AccuracyAssessment2D<Quad_Soln_Block>::
ComputeSolutionErrors(const unsigned int &parameter,
		      double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
								   const Vector2D &,
								   const unsigned int &) const ){
  
  
}

#endif
