/* Rte2DQuadSolvers.cc:  2D Rte Equation 
                         Multi-Block Quadrilateral Mesh Solvers. */

#include "Rte2DQuadSolvers.h"


/********************************************************
 * Routine: Rte2DQuadSolver                             *
 *                                                      *
 * Computes solutions to 2D Rte equations on 2D         *
 * quadrilateral multi-block solution-adaptive mesh.    *
 *                                                      *
 ********************************************************/
int Rte2DQuadSolver(char *Input_File_Name_ptr,
		    int batch_flag) {

  // create a solver object
  Rte2DSolver SolverObj(Input_File_Name_ptr, batch_flag);

  // execute solver
  return SolverObj.StandAloneSolve();
}
