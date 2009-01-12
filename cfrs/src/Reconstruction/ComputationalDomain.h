/* ComputationalDomain.h: Header file defining the template computational domain*/

#ifndef _COMPUTATIONALDOMAIN_INCLUDED
#define _COMPUTATIONALDOMAIN_INCLUDED

/* Include defined header file. */

#include "include/ComputationalCell.h"
#include "Reconstruction1D/Reconstruct1DInput.h"
#include "Reconstruction2D/Reconstruct2DInput.h"
#include "Reconstruction3D/Reconstruct3DInput.h"
#include "Grid/Grid2D/QuadrilateralGrid.h"
#include "Grid/Grid3D/Grid3DHexaBlock.h"
#include "ReconstructionFunction.h"
#include "CENO_DataAnalysis.h"
#include "Common/TecplotAuxData.h"

/*************************************************************
               PRE-PROCESSING PARAMETERS                      
***************************************************************/
/* Define CENO_Padding if PWL reconstruction is required also in the cells adjacent to a non-smooth reconstruction */
//#define CENO_Padding

// If "LimiterBasedOnlyOnGaussPoints" is not define,
// then the values at the nodes of the cell are also used for
// determining the admissible value for the limiter
#define LimiterBasedOnlyOnGaussPoints


template<typename T> inline
void Print_Progress(const T & Counter, const T & Frequency){

  if ( (Counter - Frequency*(Counter/Frequency)) == 0 ){
    std::cout << ".";
    std::cout.flush();
  }
}

// Define the ComputationalDomain class
template< SpaceType SpaceDimension = TwoD, class GeometryType = Cell2D_Quad, class SolutionType = double>
class ComputationalDomain;

/************************************************
*     Friend Functions : ComputationalCell       *
************************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
bool operator==(const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& left,
		const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& right);

template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
bool operator!=(const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& left,
		const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& right);

template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
std::ostream& operator<< (std::ostream& os,
			  const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& Obj);


 /****************************************************************************
 * TEMPLATIZED CLASS: ComputationalDomain                                   *
 *                                                                          *
 * Container for ComputationalCells.                                        *
 * 
 * VARIABLES:                                                               *
 *                                                                       templated to SpaceDimension   *
 ***************************************************************************/

template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
class ComputationalDomain{
 public:
  typedef ComputationalCell<SpaceDimension,GeometryType,SolutionType> CompCellType;
  typedef typename CompCellType::GeometricIntegrals GeometricIntegrals;
#ifdef __Use_Iterator__	
  typedef typename CompCellType::DerivIterator DerivativesIteratorType;
  typedef typename CompCellType::GeomCoeffIterator GeomCoeffIteratorType;
#endif
  typedef ComputationalDomain<SpaceDimension,GeometryType,SolutionType> CompDomainType;

  static const int NumberOfParameters = CompCellType::NumberOfVariables;
  
private:
  
  vector<int> N_XYZ;		/* Total number of cells in the X, Y and Z direction */
  vector<int> IndexLow;		/* Start of the computational domain */
  vector<int> IndexUp;		/* End of the computational domain */
  int Nghost;			/* Number of ghost cells */
  SolutionType L1Norm;		/* L1 Norm -> accuracy paramter */
  SolutionType L2Norm;		/* L2 Norm -> accuracy paramter */
  SolutionType LMaxNorm;        /* L_Max_Norm -> accuracy paramter */
  vector<double> MaxDeltaSolutionOverDomain; /* That's the difference between
						MaxSolution and MinSolution over 
						the whole computational domain */
  double CharacteristicLength;	/* The characteristic length of the geometry */
  double CutoffKnob;

  CompCellType*** SolnPtr;      /* Solution variable */

 public:

  /* Static Variable */
  static HeaderData VarNames;

  /* Basic member functions */
  ComputationalDomain();
  ~ComputationalDomain(){ deallocate(); }
  // Copy constructor
  ComputationalDomain(const CompDomainType & rhs);
  // Assignment operator
  CompDomainType & operator=(const CompDomainType & rhs);

  /* Memory management functions */
  void allocate();
  void deallocate();

  void SetDomain(const int & NCx, const int & Nghost);
  void SetDomain(const int & NCx, const int & NCy, const int & Nghost);
  void SetDomain(const int & NCx, const int & NCy, const int & NCz, const int & Nghost);
  void SetDomain(const Reconstruct1D_Input_Parameters & IP);
  void SetDomain(Grid2D_Quad_Block & Grid, const Reconstruct2D_Input_Parameters & IP);
  void SetDomain(const Reconstruct2D_Input_Parameters & IP);
  void SetDomain(Grid3D_Hexa_Block & Grid, const Reconstruct3D_Input_Parameters & IP);
  void SetDomain(const Reconstruct3D_Input_Parameters & IP);

  /* Access functions */
  int & iStart(void) { return IndexLow[0];}
  const int & iStart(void) const {return IndexLow[0];}
  int & jStart(void) { return IndexLow[1];}
  const int & jStart(void) const { return IndexLow[1];}
  int & kStart(void) { return IndexLow[2];}
  const int & kStart(void) const { return IndexLow[2];}
  int & iEnd(void) { return IndexUp[0];}
  const int & iEnd(void) const { return IndexUp[0];}
  int & jEnd(void) { return IndexUp[1];}
  const int & jEnd(void) const { return IndexUp[1];}
  int & kEnd(void) { return IndexUp[2];}
  const int & kEnd(void) const { return IndexUp[2];}
  int iLastCell(void){ return N_XYZ[0]-1;} /* the index of the last cell in I-direction */
  int iLastCell(void) const { return N_XYZ[0]-1;}
  int jLastCell(void){ return N_XYZ[1]-1;} /* the index of the last cell in J-direction */
  int jLastCell(void) const { return N_XYZ[1]-1;}
  int kLastCell(void){ return N_XYZ[2]-1;} /* the index of the last cell in K-direction */
  int kLastCell(void) const { return N_XYZ[2]-1;}
  double CellL1Norm(void) const { return L1Norm; } 
  double CellL2Norm(void) const { return L2Norm; } 
  double CellLMaxNorm(void) const { return LMaxNorm; } 
  vector<double> & MaxDeltaSolutionDomain(void) { return MaxDeltaSolutionOverDomain; }
  double & MaxDeltaSolutionDomain(int parameter) { return MaxDeltaSolutionOverDomain[parameter]; }
  double & GeomCharacteristicLength(void) { return CharacteristicLength; }
  const double & KnobCutoff(void) { return CutoffKnob;}
  const int & NumberOfCellRings(void) const { return SolnPtr[0][0][0].CellRings(); }

  /* Operation functions */
  int NumberOfTaylorDerivatives() const {
    return SolnPtr[0][0][0].NumberOfTaylorDerivatives();
  }
  int NumberOfNodes(void) const;
  bool null(void) const;

  void BC_Reflection(void);

  void ReconstructSolution(const Reconstruct1D_Input_Parameters & IP);
  void ReconstructSolution(const Reconstruct2D_Input_Parameters & IP);
  void ReconstructZonalSolution(const Reconstruct2D_Input_Parameters & IP);
  void ReconstructSolution(const Reconstruct3D_Input_Parameters & IP);
  void ReconstructSolutionAndBoundary(const Reconstruct1D_Input_Parameters & IP);
  void ReconstructToLimitedPiecewiseLinear(int iCell, const int Limiter);
  void ReconstructToLimitedPiecewiseLinear(int iCell, int jCell, const int Limiter);
  void ReconstructToLimitedPiecewiseLinear(int iCell, int jCell, int kCell, const int Limiter);

  template<class InputParameters>
    void SetInitialData(const InputParameters & IP); /* template based on the class of the input parameters */

  void UpdateSubgridSolution(void);
  void SetMaxDeltaSolutionOverDomain(void);

  template<class InputParameters>
    void AssessReconstructionAccuracy(const InputParameters & IP);

  void ComputeMultipleCorrelationCoefficient(void);
  void ComputeMultipleCorrelationCoefficient2D(void);

  /* Output functions */
  void DefineHeader(const HeaderData & NameOfThePrintedVariables){ VarNames = NameOfThePrintedVariables; }
  void OutputMeshNodesTecplot(std::ofstream &output_file,const bool Title = true) const;
  void OutputMeshCellsTecplot(std::ofstream &output_file,const bool Title = true) const;
  void OutputNodesTecplot(std::ofstream &output_file,const bool Title = true) const;
  void OutputSolutionNodesTecplot(std::ofstream &output_file,const bool Title = true) const;
  void OutputFullSolutionNodesTecplot(std::ofstream &output_file,const bool Title = true) const;
  void OutputSolutionAtGaussPointsTecplot(std::ofstream &output_file,const bool Title = true) const;
  void OutputSolutionCellTecplot(std::ofstream &output_file,const bool Title = true) const;
  void OutputFullSolutionCellTecplot(std::ofstream &output_file,const bool Title = true) const;
  void OutputSolutionCellTecplotOneZone(std::ofstream &output_file,const bool Title = true) const;
  void OutputMultipleCorrelationCoefficient(std::ofstream &output_file,const bool Title = true) const;
  void OutputSmoothnessIndicatorFlag(std::ofstream &output_file,const bool Title = true) const;
  void OutputSmoothnessIndicatorAtNodesTecplot(std::ofstream &output_file, const bool Title = true) const;
  void OutputReconstructedSolutionTecplot(std::ofstream &output_file,const bool Title = true) const;
  void OutputPWC(std::ofstream &output_file, const bool Title = true) const;
  void PrintErrorNorms(void) const;

  /* Overloaded operators */
  CompCellType & operator( )(const int &i);
  CompCellType & operator( )(const int &i, const int &j);
  CompCellType & operator( )(const int &i, const int &j, const int &k);


  /* Friend functions */
  // == operator
  friend bool operator==<SpaceDimension,GeometryType,SolutionType>
    (const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& left,
     const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& right);

  // != operator
  friend bool operator!=<SpaceDimension,GeometryType,SolutionType>
    (const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& left,
     const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& right);
 
  // << operator
  friend std::ostream& operator<< <SpaceDimension,GeometryType,SolutionType>
    (std::ostream& os, const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& Obj);


};

/* Include implementation of member functions */
#include "ComputationalDomain_Implementations.h"

#endif
