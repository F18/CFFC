/*!\file test_test_ComputationalDomain.cc
  \brief Regression tests for templated class ComputationalDomain . */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */

/* Include CFFC header files */
#include "CFD/CFD.h"
#include "Reconstruction/ComputationalDomain.h"
#include "TestData.h"

/* Include files required for kExact Reconstruction */
#include <vector>
#include "../../../src_3D/Math/LinearSystems.h"
#include "../../../src_3D/Utilities/Utilities.h"
#include "include/TypeDefinition.h"
#include "Reconstruction/ReconstructionHelpers.h"

using namespace std;

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_ComputationalDomain : public TestData {

    // Local variables
  public:
    int error_flag;
    // Reconstruction3D input variables and parameters:
    Reconstruct3D_Input_Parameters Input_Parameters;
    int CurrentOrderOfReconstruction;
    int command_flag;

    // Solution variables.
    ComputationalDomain<ThreeD,Cell3D_Hexa,double> SolnBlkDouble;
    // Mesh variables
    Grid3D_Hexa_Block   Grid;

    // Constructor
    Data_ComputationalDomain();

    ~Data_ComputationalDomain(){};

  private:
    
  };


  // Constructor
  Data_ComputationalDomain::Data_ComputationalDomain(){
    /* Set the global path for this test suite. 
       It's a good practice to put it in the constructor of the data class in order to be set
       automatically for each individual test. Declare it relative to the /src_2D directory,
       otherwise the framework might not find the input and output files. */
    
    set_test_suite_path("Reconstruction/UnitTests/");
            
    CurrentOrderOfReconstruction = -1;
  }


  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_ComputationalDomain> ComputationalDomain_TestSuite;
  typedef ComputationalDomain_TestSuite::object ComputationalDomain_object;


  /******************************************************************************************************
   ******************************************************************************************************
   *                      ********    *******        ***       ********        ***                      *
   *                         **       **           **   **        **         **   **                    *
   *                         **       **          **              **        **                          *
   *                         **       *******      *****          **          *****                     *
   *                         **       **                **        **               **                   *
   *                         **       **            *    **       **          *   **                    *
   *                         **       *******        ****         **           ****                     *
   ******************************************************************************************************
   ******************************************************************************************************/


  /*********************************************
     TESTING OPERATIONS:
        -->  ensure("ConditionName", test_condition);
        -->  ensure_not("ConditionName", test_condition);
        -->  ensure_equals("ConditionName", expected_value, got_value);
                 // operators '<<' and '!=' must be defined!!!
        -->  ensure_distance("ConditionName", expected_value, got_value, tolerance);
	         // operators '<<', '>=' and '<=' must be defined!!!
        -->  fail("Message");

     Obs: "ConditionName" is optional
  */



  /* Test 1: Set up the Grid and output the nodal points and cell centers to seperate output files*/
  template<>
  template<>
  void ComputationalDomain_object::test<1>()
  {

    set_test_name("Grid check: nodal and cell center points");
    
    set_local_input_path("ComputationalDomain_Data/");
    set_local_output_path("ComputationalDomain_Data/");        
    Open_Input_File("reconstruct3D.in");

    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
                                                      input_file_name,
                                                      command_flag);
    // update the CurrentOrderOfReconstruction
    CurrentOrderOfReconstruction = Input_Parameters.RecOrder();
    
    // output the input parameters to the screen
    cout << Input_Parameters << "\n";
    cout.flush();

    std::cout << "\n Create Grid:\n";
    // Set Grid
    Grid.Create_Block(Input_Parameters.Box_Length,
		      Input_Parameters.Box_Width,
		      Input_Parameters.Box_Height,
		      ZERO,    // x-orig
		      ZERO,    // y-orig
		      ZERO,    // z-orig
		      ZERO,    // alpha
		      ZERO,    // beta
		      ZERO,    // gamma
		      BC_NONE, // top
		      BC_NONE, // bottom
		      BC_NONE, // north
		      BC_NONE, // south
		      BC_NONE, // west
		      BC_NONE, // east
		      Input_Parameters.NCells_Idir,
		      Input_Parameters.NCells_Jdir,
		      Input_Parameters.NCells_Kdir,
		      Input_Parameters.Nghost_Cells);


    std::cout << " Set the Computational Domain:\n";
    // Set the Computation Domain
    SolnBlkDouble.SetDomain(Grid,Input_Parameters);

    // Output the nodal points of the mesh to a file
    Output_Mesh_Nodes_Tecplot(SolnBlkDouble, Input_Parameters);

    // Output the cell center locations of the mesh to a file
    Output_Mesh_Cells_Tecplot(SolnBlkDouble, Input_Parameters);    

  }



  /* Test 2: Generate the coefficient matrix A for inspection. This is the matrix A which is associated with the (Ax=b) used in 
     the kExactReconstruction procedure */
  template<>
  template<>
  void ComputationalDomain_object::test<2>()
  {

    set_test_name("Coefficient Matrix");
    
    set_local_input_path("ComputationalDomain_Data/");
    set_local_output_path("ComputationalDomain_Data/");
    Open_Input_File("reconstruct3D.in");

    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
                                                      input_file_name,
                                                      command_flag);
    // update the CurrentOrderOfReconstruction
    CurrentOrderOfReconstruction = Input_Parameters.RecOrder();
    
    // output the input parameters to the screen
    cout << Input_Parameters << "\n";
    cout.flush();

    std::cout << "\n Create Grid:\n";
    // Set Grid
    Grid.Create_Block(Input_Parameters.Box_Length,
		      Input_Parameters.Box_Width,
		      Input_Parameters.Box_Height,
		      ZERO,    // x-orig
		      ZERO,    // y-orig
		      ZERO,    // z-orig
		      ZERO,    // alpha
		      ZERO,    // beta
		      ZERO,    // gamma
		      BC_NONE, // top
		      BC_NONE, // bottom
		      BC_NONE, // north
		      BC_NONE, // south
		      BC_NONE, // west
		      BC_NONE, // east
		      Input_Parameters.NCells_Idir,
		      Input_Parameters.NCells_Jdir,
		      Input_Parameters.NCells_Kdir,
		      Input_Parameters.Nghost_Cells);


    std::cout << " Set the Computational Domain:\n";
    // Set the Computation Domain
    SolnBlkDouble.SetDomain(Grid,Input_Parameters);

    // Output the nodal points of the mesh to a file
    Output_Mesh_Nodes_Tecplot(SolnBlkDouble, Input_Parameters);

    // Variables for Reconstruction
    // vector<int> i_index, j_index, k_index;
    int *i_index(NULL), *j_index(NULL), *k_index(NULL);
    int NumOfCellsInOneDirection, StencilSize;
    int i, j, k;


    /* Determine the number of cells in the stencil based on the number of rings */
    NumOfCellsInOneDirection = 2*SolnBlkDouble.SolnPtr[0][0][0].CellRings() + 1;
    StencilSize = NumOfCellsInOneDirection*NumOfCellsInOneDirection*NumOfCellsInOneDirection;
    i_index = new int [StencilSize];
    j_index = new int [StencilSize];
    k_index = new int [StencilSize];
    

    /* Loop Through each cell in the domain to make a stencil and reconstruct the solution
       Solve reconstruction for each cell in the domain plus 2 additional layers of cells for each boundary.
       These boundary cells are used for checking the goodness of fit of the first domain cell. */
    
    for (i=SolnBlkDouble.iStart()-2; i<=SolnBlkDouble.iEnd()+2; ++i){
      for (j=SolnBlkDouble.jStart()-2; j<=SolnBlkDouble.jEnd()+2; ++j){
        for (k=SolnBlkDouble.kStart()-2; k<=SolnBlkDouble.kEnd()+2; ++k){
          
          /* Make Stencil */
          MakeReconstructionStencil(SolnBlkDouble.SolnPtr[0][0][0].CellRings(),i,j,k,i_index,j_index,k_index);

          /***************************************************************************
           * kExact_Reconstruction for 3D                                             *
           *                                                                          *
           * This subroutine determines the coefficients of a Taylor series expansion *
           * which approximates the solution over the domain of the cell specified    *
           * by "i_index[0]", "j_index[0]", and "k_index[0]" indexes.                 *
           ***************************************************************************/

          static const int NumberOfParameters = 1;
    
          // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
          int ND(SolnBlkDouble.NumberOfTaylorDerivatives());	      /* number of Taylor expansion coefficients */
          DenseMatrix A(StencilSize-1,ND-1);                          /* the matrix which the linear system is solved for */
          DenseMatrix All_Delta_U (StencilSize-1,NumberOfParameters); /* matrix for storing U[neighbour]-U[cell] */
          ColumnVector GeomWeights(StencilSize);                      /* the column vector of the geometric weights */
          Vector3D* DeltaCellCenters;                                 /* array for storing the X-distance and Y-distance between  
                                                                      the cell center of neighbour cells and the one of i,j,k cell */

          int krank;                                                  /* the final rank of A matrix is returned here */
          int IndexSumZ, IndexSumY, IndexSumX, P1, P2, P3;
          double CombP1X, CombP2Y, CombP3Z;
          double PowDistanceXC, PowDistanceYC, PowDistanceZC;
          int cell, i_temp, parameter;
          double WeightsSum(0.0);
          double IntSum1(0.0),IntSum2(0.0);
    
          // Allocate memory
          DeltaCellCenters = new Vector3D [StencilSize];
    
          // *********  Assign the average solution to D00 ***********
          SolnBlkDouble(i_index[0],j_index[0],k_index[0]).CellDeriv(0,0,0) = SolnBlkDouble(i_index[0],j_index[0],k_index[0]).CellSolution();

          if (ND == 1){ // piecewise constant
            return;
          }

          // START:   Set the LHS and RHS of the linear system 
          // ***************************************************

          // Step1. Compute the normalized geometric weights
          for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

            /* Compute the X, Y, and Z component of the distance between
               the cell center of the neighbours and the reconstructed cell */
            DeltaCellCenters[cell] = SolnBlkDouble(i_index[cell],j_index[cell],k_index[cell]).CellCenter() - 
              SolnBlkDouble(i_index[0],j_index[0],k_index[0]).CellCenter();
            /* Compute the geometric weights and their sum (this is used for normalization)
               based on the distance to each control volume */
            GeomWeights(cell) = sqrt(DeltaCellCenters[cell].x*DeltaCellCenters[cell].x + 
                                     DeltaCellCenters[cell].y*DeltaCellCenters[cell].y +
                                     DeltaCellCenters[cell].z*DeltaCellCenters[cell].z);
            GeomWeights(cell) *= GeomWeights(cell);
            GeomWeights(cell) = 1.0/(1.0E-15 + GeomWeights(cell));
    
            WeightsSum += GeomWeights(cell);
          }

          // Step2. Set the approximate equations
          for (cell=1 ; cell<StencilSize; ++cell){ //for each cell in the stencil
    
            // compute the normalized geometric weight
            GeomWeights(cell) /= WeightsSum;

            // *** SET the matrix A of the linear system (LHS) ***
            /* compute for each derivative the corresponding entry in the matrix of the linear system */
            for (i_temp=1; i_temp<=SolnBlkDouble(i_index[0],j_index[0],k_index[0]).CellDeriv().LastElem(); ++i_temp){
              // build the row of the matrix
              P1 = SolnBlkDouble(i_index[0],j_index[0],k_index[0]).CellDeriv(i_temp).P1();  // identify P1
              P2 = SolnBlkDouble(i_index[0],j_index[0],k_index[0]).CellDeriv(i_temp).P2();  // identify P2
              P3 = SolnBlkDouble(i_index[0],j_index[0],k_index[0]).CellDeriv(i_temp).P3();  // identify P3

              A(cell-1,i_temp-1) = 0.0;  // set sumation variable to zero
              CombP3Z = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
              PowDistanceZC = 1.0;  // initialize PowDistanceZC
      
              // Compute geometric integral over the neighbour's domain      
              for (IndexSumZ = 0; IndexSumZ<=P3; ++IndexSumZ){
                CombP2Y = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
                PowDistanceYC = 1.0;  // initialize PowDistanceYC
                IntSum2 = 0.0;         // reset internal summation variable
                
                for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
                  CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
                  PowDistanceXC = 1.0; // initialize PowDistanceXC
                  IntSum1 = 0.0;        // reset internal summation variable
          
                  for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
                    IntSum1 += ( CombP1X*PowDistanceXC*
                                 SolnBlkDouble(i_index[cell],j_index[cell],k_index[cell]).CellGeomCoeff(P1-IndexSumX,P2-IndexSumY,P3-IndexSumZ) );
                    // update the binomial coefficients
                    CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
                    PowDistanceXC *= DeltaCellCenters[cell].x;      // Update PowDistanceXC
                  }//endfor

                  IntSum2 += CombP2Y*PowDistanceYC*IntSum1;
                  CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
                  PowDistanceYC *= DeltaCellCenters[cell].y;      // Update PowDistanceYC
                }//endfor

                A(cell-1,i_temp-1) += CombP3Z*PowDistanceZC*IntSum2;  // update the external sum
        
                CombP3Z = (P3-IndexSumZ)*CombP3Z/(IndexSumZ+1); // the index is still the old one => expression for "nC k+1"
                PowDistanceZC *= DeltaCellCenters[cell].z;      // Update PowDistanceYC
              }//endfor

              A(cell-1,i_temp-1) -= SolnBlkDouble(i_index[0],j_index[0],k_index[0]).CellGeomCoeff(i_temp,true,true,true);
#if 0
              // apply geometric weighting
              // A(cell-1,i_temp-1) *= GeomWeights(cell);
#endif
            } // end for (i_temp) -- Number of Derivatives

          } // end for (cell) -- Stencil Size
      
          // The following prints the A matrix for inspection. Only the cells along the diagonal of the grid are printed out here
          if (i==j && j==k) {
            std::cout <<"Mesh Row Number: " << i <<endl;
            Print_(A);
          }
          
          // Free memory
          delete [] DeltaCellCenters; DeltaCellCenters = NULL;
        } //end for k
      } //end for j
    }//end for i

    // Free memory
    delete [] i_index; i_index = NULL;
    delete [] j_index; j_index = NULL;
    delete [] k_index; k_index = NULL;
    
  }
}



// Test suite constructor
tut::ComputationalDomain_TestSuite ComputationalDomain_Test("Template Class:ComputationalDomain");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

