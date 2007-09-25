/*!\file HighOrder1D_CENO_Reconstruction.h
  \brief Subroutines for performing high-order CENO solution reconstruction in 1D solvers. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "ReconstructionSolvers1D.h"


// Compute pseudo-inverse for the all the cells required to performed CENO reconstruction
template<typename SolutionContainer, typename Input_Parameters_Type> inline
void Compute_CENO_PseudoInverse(SolutionContainer * SolnBlk,
				const Input_Parameters_Type & IP){

#ifdef _CENO_SPEED_EFFICIENT

  if(Input_Parameters.i_ReconstructionMethod == RECONSTRUCTION_CENO){

    if (IP.Verbose()){
      cout << "\n\n --> Preprocessing Step\n" << "---------------------------------------\n"
	   << "    > Compute the pseudo-inverse in \n       every computational cell.\n";
    }
    
    CENO1D_HighOrder_Reconstruction_PseudoInverse(Soln_ptr,Input_Parameters);


    if (IP.Verbose()){
      cout << "---------------------------------------\n" << " --> Preprocessing Step Done\n";
    }
  }

#endif //_CENO_SPEED_EFFICIENT

}


#ifdef _CENO_SPEED_EFFICIENT
template<typename SolutionContainer, typename Input_Parameters_Type> inline
void CENO1D_HighOrder_Reconstruction_PseudoInverse(SolutionContainer * SolnBlk,
						   const Input_Parameters_Type & IP){

  // Impose minimum number of ghost cells
  require(SolnBlk[0].Nghost >= 1 + 2*Soln[0].CellRings(),
	  "\nCENO1D_HighOrder_Reconstruction_PseudoInverse() ERROR: Insufficient number of ghost cells!");

  /* Determine the number of cells in the stencil based on the number of rings */
  int StencilSize(1 + 2*Soln[0].CellRings()); // CellsInOneDirection = 1 + 2*Soln[0].CellRings();
  int ExtraCells(1+Soln[0].CellRings()); // ghost cells in which solution needs to be reconstructed
  vector<int> i_index(StencilSize); 
  int cell;
  int ICl, ICu;		  //!< Mesh indexes

  /* Set the indexes for the mesh */
  ICl = Soln_ptr[0].Nghost;
  ICu = Soln_ptr[0].Nghost + IP.Number_of_Cells - 1;

  // Compute the pseudo-inverse of the LHS term in the k-exact reconstruction
  for (cell=ICl-ExtraCells; cell<=ICu+ExtraCells; ++cell){
    /* Make Stencil */
    MakeReconstructionStencil(Soln[0].CellRings(),cell,i_index);

    /* Determine the LHS term based on the stencil of
       the current cell and the geometry of the the mesh. */
    kExact_Reconstruction_LHS(Soln,i_index,cell);

    /* Compute the pseudo-inverse and override the LHS term */
    Soln[cell].LHS().pseudo_inverse_override();
  }

}
#endif //_CENO_SPEED_EFFICIENT
