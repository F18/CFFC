/* CompCell1D.h: Header file defining 1D computational cell types*/

#ifndef _COMPCELL1D_INCLUDED
#define _COMPCELL1D_INCLUDED

/* Include defined header file. */

#ifndef _NumericalLibrary_INCLUDED
#include "Math/NumericalLibrary.h"
#endif // _NumericalLibrary_INCLUDED

#ifndef _CELL1D_INCLUDED
#include "Grid/Grid1D/Cell1D.h"
#endif //_CELL1D_INCLUDED

#ifndef _TESTFUNCTIONS_INCLUDED
#include "TestFunctions/TestFunctions.h"
#endif // _TESTFUNCTIONS_INCLUDED

#ifndef _MATH_MACROS_INCLUDED
#include "Math/Math.h"
#endif // _MATH_MACROS_INCLUDED
 

using namespace std;

/* Define the classes */

/********************************************************
 * Class: CompCell1D_Uniform                            *
 *                                                      *
 * Member functions                                     *
 *   Geometry                                           *
 *     geom      -- Object of type Cell1D,              *
 *                  stores the geometry of the cell     *
 *                                                      *
 *   Subdomain                                          *
 *      In each cell are defined some points, in which  *
 *      flow parameters are computed.                   *
 *                                                      *
 *     *SubGrid  -- Pointer of type double,             *
 *                  points to the vector of the points  *
 *                  defining the subdomain              *
 *                                                      *
 *   Flow conditions                                    *
 *                                                      *
 * Truncated Taylor series expansion reconstruction for *
 *              order k                                 *
 *          k                                           *
 *     R = SUM alpha[p]*D[p]*(x - x_j)^p                *
 *         p=0                                          *
 *                                                      *
 *     U_cell    -- The average cell quantity           *
 *     *D        -- Pointer of type double, points to   *
 *                  the vector storing the coeff. D     *      
 *     RO        -- Reconstruction Order                *
 *                                                      *
 *   Flow conditions in the subdomain                   *
 *    *U_points  -- Pointer of type double,             *
 *                  points to the vector storing the    *
 *                  value of the flow parameters in     *
 *                  each point of the subdomain         *
 *    *U_deriv_num  -- Pointer of type double,          *
 *                  points to the vector storing the    *
 *                  representation of the numerical     *
 *                  function derivatives                *
 *    N          -- Number of subgrid points            *
 *    delta_sg   -- Delta_Sub_Grid                      *
 *                                                      *
 *    SetDynamicVariable -- Set the length of the       *
 *                          vectors allocated           *
 *                          dynamically and initial     *
 *                          values for them             *
 *    SetCellParameters                                 *
                   -- Set the all parameters that       *
 *                    define the cell                   *
 *    SetSubDomainSolution -- Computes the solution in  *
 *                          the points of the SubDomain *
 *                          using the truncated Taylor  *
 *                          series expansion            *
 *                                                      *
 *   FinalOrder  -- The final order obtained during     *
 *                  the reconstruction                  *
 *******************************************************/
class CompCell1D_Uniform{
 private:
 public:
  Cell1D_Uniform geom;
  double *SubGrid;         
  double U_cell;
  double *D;
  int RO;
  double *U_points;
  double *U_deriv; // must be initialized
  int N;
  double delta_sg;
                    // Made public sa can access them.
  int FinalOrder;

  /* Creation, copy, and assignment constructors. */
  //CompCell1D_Uniform()
  //Use the default constructor

  //~CompCell1D_Uniform();
  //Use the default destructor

  void allocate(int NbSubDomainPoints, int ro);
  void deallocate();
  void SetCellParameters(double loc, double delta, int nr = 3, int ro = 0);
  void SetCellParameters(double loc, double delta, double delta_sg, int ro);
  void SetDynamicVariable(int nr, int ro);
  void SetDynamicVariable(double delta, int ro);
  void SetSubDomainSolution();
};

/**************************************************************************
 * CompCell1D_Uniform::allocate -- Allocate memory.                       *
 **************************************************************************/

inline void CompCell1D_Uniform::allocate(int NbSubDomainPoints, int ro){

  N = NbSubDomainPoints;
  SubGrid = new double [N];
  U_points = new double [N];
  RO = ro;
  FinalOrder = RO;
  D = new double [RO+1];
}

/**************************************************************************
 * CompCell1D_Uniform::deallocate -- Deallocate memory.                 *
 **************************************************************************/

inline void CompCell1D_Uniform::deallocate(){
  
  delete [] SubGrid; SubGrid = NULL;
  delete [] U_points; U_points = NULL;
  delete [] D; D = NULL;
}

/********************************************************
 * Class: CompCell1D_NonUniform                         *
 *                                                      *
 * Member functions                                     *
 *   Geometry                                           *
 *     geom      -- Object of type Cell1D,              *
 *                  stores the geometry of the cell     *
 *                                                      *
 *   Subdomain                                          *
 *      In each cell are defined some points, in which  *
 *      flow parameters are computed.                   *
 *                                                      *
 *     *SubGrid  -- Pointer of type double,             *
 *                  points to the vector of the points  *
 *                  defining the subdomain              *
 *                                                      *
 *   Flow conditions                                    *
 * Truncated Taylor series expansion reconstruction for *
 *              order k                                 *
 *          k                                           *
 *     R = SUM alpha[p]*D[p]*(x - x_j)^p                *
 *         p=0                                          *
 *                                                      *
 *     U_cell    -- The average cell quantity           *
 *     *D        -- Pointer of type double, points to   *
 *                  the vector storing the coeff. D     *
 *     RO        -- Reconstruction Order                *      
 *                                                      *
 *   Flow conditions in the subdomain                   *
 *    *U_points  -- Pointer of type double,             *
 *                  points to the vector storing the    *
 *                  value of the flow parameters in     *
 *                  each point of the subdomain         *
 *    *U_deriv_num  -- Pointer of type double,          *
 *                  points to the vector storing the    *
 *                  representation of the function      *
 *                  numerical derivatives               *
 *    N          -- Number of subgrid points            *
 *    delta_sg   -- Delta_Sub_Grid                      *
 *    SetDynamicVariable -- Set the length of the       *
 *                          vectors allocated           *
 *                          dynamically and initial     *
 *                          values for them             *
 *    SetCellParameters                                 *
 *                  -- Set the all parameters that       *
 *                    define the cell                   *
 *    SetSubDomainSolution -- Computes the solution in  *
 *                          the points of the SubDomain *
 *                          using the truncated Taylor  *
 *                          series expansion            *
 *    SetDerivative_Representation                      *
 *                 -- Computes the deivative of any order *
 *                    for the SubDomain points            *
 *    ErrorRec = average reconstruction error           *
 *               ErrorRec = Integral(|F_numeric - F_analytic|)/ dx *
 *                                                      *
 *   FinalOrder  -- The final order obtained during     *
 *                  the reconstruction                  *
 *******************************************************/

class CompCell1D_NonUniform{
 private:
 public:
  Cell1D_NonUniform geom;
  double *SubGrid;         
  double U_cell;
  double *D;
  int RO;
  double *U_points;
  double *U_deriv_num;
  int N;
  double delta_sg;
  double ErrorRec;
  int FinalOrder;
                    // Made public sa can access them.

  /* Creation, copy, and assignment constructors. */
  //CompCell1D_NonUniform()
  //Use the default constructor

  //~CompCell1D_NonUniform();
  //Use the default destructor

  void allocate(int NbSubDomainPoints, int ro);
  void deallocate();
  void SetCellParameters(double loc, double delta, int nr = 3, int ro = 0);
  void SetCellParameters(double loc, double delta, double delta_sg, int ro = 0);
  void SetDynamicVariable(int nr, int ro);
  void SetDynamicVariable (double delta, int ro);
  void SetSubDomainSolution();
  void SetSubDomainSolution(CompCell1D_NonUniform &Cell);
  void SetDerivative_Representation(int order);
  double SolutionAtCoordinates(double x);
  void ComputeReconstructionError (FunctionType1D ExactFunction);
  friend  double ErrorNorm1D(double x);
};

/**************************************************************************
 * CompCell1D_Uniform::allocate -- Allocate memory.                       *
 **************************************************************************/

inline void CompCell1D_NonUniform::allocate(int NbSubDomainPoints, int ro){

  N = NbSubDomainPoints;
  SubGrid = new double [N];
  U_points = new double [N];
  U_deriv_num = new double [N];
  RO = ro;
  FinalOrder = RO;
  D = new double [RO+1];
  if ((SubGrid==NULL)||(U_points==NULL) || (U_deriv_num==NULL) || (D == NULL))
    cout << "\nNot enough memory!\n"; cout.flush();
}

/**************************************************************************
 * CompCell1D_Uniform::deallocate -- Deallocate memory.                 *
 **************************************************************************/

inline void CompCell1D_NonUniform::deallocate(){
  
  delete [] SubGrid; SubGrid=NULL;
  delete [] U_points; U_points=NULL;
  delete [] U_deriv_num; U_deriv_num=NULL;
  delete [] D; D=NULL;
}

#endif // _COMPCELL1D_INCLUDED
