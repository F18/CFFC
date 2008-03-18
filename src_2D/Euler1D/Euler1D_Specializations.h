/*! \file Euler1D_Specializations.h
  \brief Header file defining template specializations for 1D Euler subroutines. */

#ifndef _EULER_1D_SPECIALIZATIONS_INCLUDED
#define _EULER_1D_SPECIALIZATIONS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "Euler1D.h"


/*! 
 * Spcialization of SolutionClassTraits for Euler1D_pState.
 */
template<>
class SolutionClassTraits<Euler1D_pState> {
 public:
  typedef Euler1D_cState SOLN_cSTATE;
};


/*! 
 * Implementation of periodic boundary conditions
 */
inline void BC_Periodic (Euler1D_UniformMesh *SolnBlk,
			 const CFD1D_Input_Parameters &IP){

  int Nghost(SolnBlk[0].Nghost);
  int DoubleNghost(2*Nghost);

  for (int i=Nghost; i>0; --i){
    // ***** left end ******
    // primitive variables
    SolnBlk[i-1].W = SolnBlk[IP.Number_of_Cells + i - 1].W;
    // conserved variables
    SolnBlk[i-1].U = SolnBlk[IP.Number_of_Cells + i - 1].U;

    // ***** right end *****
    // primitive variables
    SolnBlk[IP.Number_of_Cells + DoubleNghost - i].W = SolnBlk[DoubleNghost - i].W;
    // conserved variables
    SolnBlk[IP.Number_of_Cells + DoubleNghost - i].U = SolnBlk[DoubleNghost - i].U;
  }
}

/*! 
 * Implementation of periodic boundary conditions.
 * This subroutine assigns the solution obtained at the cell centroids
 * with a high-order reconstruction to the corresponding ghost 
 * cells for periodic BC.
 */
inline void BC_Periodic_CentroidSol (Euler1D_UniformMesh *SolnBlk,
				     const CFD1D_Input_Parameters &IP,
				     Euler1D_UniformMesh::HighOrderType & 
				     (Euler1D_UniformMesh::*AccessToHighOrderVar)(void) =
				     &Euler1D_UniformMesh::CellHighOrder){

  int Nghost(SolnBlk[0].Nghost);
  int DoubleNghost(2*Nghost);

  for (int i=Nghost; i>0; --i){
    // ***** left end ******
    (SolnBlk[i-1].*AccessToHighOrderVar)().CellDeriv(0) = 
      (SolnBlk[IP.Number_of_Cells + i - 1].*AccessToHighOrderVar)().CellDeriv(0);
    
    // ***** right end *****
    (SolnBlk[IP.Number_of_Cells + DoubleNghost - i].*AccessToHighOrderVar)().CellDeriv(0) = 
      (SolnBlk[DoubleNghost - i].*AccessToHighOrderVar)().CellDeriv(0);
  }
}

/*! 
 * Implementation of wall boundary conditions.
 */
inline void BC_Wall (Euler1D_UniformMesh *SolnBlk,
		     const CFD1D_Input_Parameters &IP){

  int ICu;
  int Nghost(SolnBlk[0].Nghost);

  ICu = Nghost + IP.Number_of_Cells-1;

  for (int i=1; i<=Nghost; ++i){
    // **** left ghost cells ******
    // primitive variables
    SolnBlk[Nghost-i].W = SolnBlk[Nghost+i-1].W;
    // change the sign of the velocity
    SolnBlk[Nghost-i].W.v = -SolnBlk[Nghost-i].W.v;
    // update the conservative variables
    SolnBlk[Nghost-i].U = U(SolnBlk[Nghost-i].W);

    // **** right ghost cells ******
    // primitive variables
    SolnBlk[ICu+i].W = SolnBlk[ICu-i+1].W;
    // change the sign of the velocity
    SolnBlk[ICu+i].W.v = -SolnBlk[ICu+i].W.v;
    // update the conservative variables
    SolnBlk[ICu+i].U = U(SolnBlk[ICu+i].W);
  }
}

/*! 
 * Specialization of BCs template subroutine for Euler1D_UniformMesh.
 * Applies boundary conditions particular to Euler1D problems by setting
 * the average solution in the ghost cells.
 */
template<> inline
void BCs<Euler1D_UniformMesh>(Euler1D_UniformMesh *SolnBlk,
			      const CFD1D_Input_Parameters &IP,
			      Euler1D_UniformMesh::HighOrderType & 
			      (Euler1D_UniformMesh::*AccessToHighOrderVar)(void)){

  int ICl(SolnBlk[0].ICl);
  int ICu(SolnBlk[0].ICu);

  switch(IP.i_ICs){
    /* Apply periodic BCs for some specific problems */
  case IC_SIN_WAVE:
    BC_Periodic(SolnBlk,IP);
    break;
  case IC_JIANG_WAVE:
    BC_Periodic(SolnBlk,IP);
    break;
  case IC_DENSITY_STEP_WAVE:
    BC_Periodic(SolnBlk,IP);
    break;
  case IC_CONVECTION_OF_DIFFERENT_SHAPES:
    BC_Periodic(SolnBlk,IP);
    break;
    
    /* Apply reflective BCs for some specific problems */
  case IC_BLAST_WAVE_INTERACTION:
    BC_Wall(SolnBlk,IP);
    break;

    /* By default, constant extrapolation boundary
       conditions are applied at either end of the mesh. */
  default:    
    for (int i = 0; i<ICl; ++i){
      // left end
      SolnBlk[i].U = SolnBlk[ICl].U;
      SolnBlk[i].W = SolnBlk[ICl].W;
      // right end
      SolnBlk[ICu+i+1].U = SolnBlk[ICu].U;
      SolnBlk[ICu+i+1].W = SolnBlk[ICu].W;
    }
  }
}

/*! 
 * Specialization of 'CheckSolutionPositivity' template subroutine
 * for Euler1D_UniformMesh.
 * Checks if the updated solution in "iCell" is still physical 
 * (i.e. positive density and positive energy).
 * If negative values occur, an error message is output and a 
 * runtime error is thrown.
 */
template<> inline
void CheckSolutionPositivity<Euler1D_UniformMesh>(const Euler1D_UniformMesh *SolnBlk,
						  const int & iCell){
  if (SolnBlk[iCell].U.d   <= ZERO ||
      SolnBlk[iCell].U.E   <= ZERO ||
      SolnBlk[iCell].U.e() <= ZERO ) {
    cout << "\n " << CFFC_Name() << " ERROR: Negative Density and/or Energy: \n"
	 << " node = " << iCell 
	 << "\n U = " << SolnBlk[iCell].U 
	 << "\n dUdt = "  << SolnBlk[iCell].dUdt << "\n"
	 << " rho = " << SolnBlk[iCell].U.d << "\t v = " << SolnBlk[iCell].U.v()
	 << "\t p = " << SolnBlk[iCell].U.p() << "\t e = " << SolnBlk[iCell].U.e() << "\n";

    cout.flush();

    throw runtime_error("CFFC ERROR: Solution update error");
  }
}


/*! 
 * Specialization of 'GetRightAndLeftFluxStates' template subroutine
 * for Euler1D_UniformMesh.
 * Computes the left and right states for the Riemann problem solved
 * at the interface between "Cell" and "Cell+1" cells.
 */
template<> inline
void GetRightAndLeftFluxStates<Euler1D_UniformMesh>(Euler1D_pState &Wl,
						    Euler1D_pState &Wr,
						    const int &Cell,
						    Euler1D_UniformMesh *SolnBlk,
						    const CFD1D_Input_Parameters &IP,
						    Euler1D_UniformMesh::HighOrderType & 
						    (Euler1D_UniformMesh::*AccessToHighOrderVar)(void)){

  int ICl(SolnBlk[0].ICl), ICu(SolnBlk[0].ICu);
  bool Left_UnphysicalValuesDetected(false), Right_UnphysicalValuesDetected(false) ;

  Wl = (SolnBlk[Cell].*AccessToHighOrderVar)().right_state();
  Wr = (SolnBlk[Cell+1].*AccessToHighOrderVar)().left_state();

  if (IP.i_ReconstructionMethod == RECONSTRUCTION_ENO_CHARACTERISTIC){
    // get the primitive variables from the conserved ones.
    // The conserved variables are determined based on the characteristic ones.
    Wl = W(Wl.CharactVarToConservedVar(SolnBlk[Cell].CellSolutionConsVar()));
    Wr = W(Wr.CharactVarToConservedVar(SolnBlk[Cell+1].CellSolutionConsVar()));
  }// endif
  
  /* Apply the BCs before the flux evaluation (right at the interface) */
  // ***** Left boundary **********
  if (Cell == ICl-1){
    // extrapolation BC (by default)
    Wl = Wr;

    // wall BCs
    if (IP.i_ICs == IC_BLAST_WAVE_INTERACTION){
      // change velocity sign
      Wl.v = -Wl.v;      
    }

    // periodic BCs
    if ((IP.i_ICs == IC_SIN_WAVE) || (IP.i_ICs == IC_JIANG_WAVE) ||
	(IP.i_ICs == IC_DENSITY_STEP_WAVE) || (IP.i_ICs == IC_CONVECTION_OF_DIFFERENT_SHAPES)){
      Wl = (SolnBlk[ICu].*AccessToHighOrderVar)().right_state();
      if ( IP.i_ReconstructionMethod == RECONSTRUCTION_ENO_CHARACTERISTIC ){
	Wl = W(Wl.CharactVarToConservedVar(SolnBlk[ICu].CellSolutionConsVar()));
      }	// endif
    }
  } /* endif LeftBC*/

  // *****  Right boundary *********
  if (Cell == ICu){
    // extrapolation BC (by default)
    Wr = Wl;

    // wall BCs
    if (IP.i_ICs == IC_BLAST_WAVE_INTERACTION){
      // change velocity sign
      Wr.v = -Wr.v;      
    }

    // periodic BCs
    if ((IP.i_ICs == IC_SIN_WAVE) || (IP.i_ICs == IC_JIANG_WAVE) ||
	(IP.i_ICs == IC_DENSITY_STEP_WAVE) || (IP.i_ICs == IC_CONVECTION_OF_DIFFERENT_SHAPES)){
      Wr = (SolnBlk[ICl].*AccessToHighOrderVar)().left_state();
      if (IP.i_ReconstructionMethod == RECONSTRUCTION_ENO_CHARACTERISTIC){
	Wr = W(Wr.CharactVarToConservedVar(SolnBlk[ICl].CellSolutionConsVar()));
      }	// endif
    } // endif
  } /* endif RightBC*/

  // Check validity of the left interface state.
  if (Wl.d <= ZERO ||
      Wl.p <= ZERO ) {

    // mark detection of unphysical values
    Left_UnphysicalValuesDetected = true;

    // output info
    cout << "\n " << CFFC_Name() << " ERROR: Negative Density and/or Pressure at the left interface: \n"
    	 << " node = " << Cell << "\n Wl = " << Wl << "\n" << " Wr = " << Wr << "\n"; 
    if (CENO_Execution_Mode::USE_CENO_ALGORITHM) {
      cout << " ISrho = " << (SolnBlk[Cell].*AccessToHighOrderVar)().CellInadequateFit(1)
	   << "\t ValISrho = " << (SolnBlk[Cell].*AccessToHighOrderVar)().CellSmoothnessIndicator(1) << "\n"
	   << " ISp = " << (SolnBlk[Cell].*AccessToHighOrderVar)().CellInadequateFit(3)
	   << "\t ValISp = " << (SolnBlk[Cell].*AccessToHighOrderVar)().CellSmoothnessIndicator(3) << "\n";
    }
  }

  // Check validity of the right interface state.
  if (Wr.d <= ZERO ||
      Wr.p <= ZERO ) {

    // mark detection of unphysical values
    Right_UnphysicalValuesDetected = true;

    // output user info
    cout << "\n " << CFFC_Name() << " ERROR: Negative Density and/or Pressure at the right interface: \n"
    	 << " node = " << Cell << "\n Wl = " << Wl << "\n Wr = " << Wr << "\n"; 
    if (CENO_Execution_Mode::USE_CENO_ALGORITHM) {
      cout << " ISrho = " << (SolnBlk[Cell+1].*AccessToHighOrderVar)().CellInadequateFit(1)
	   << "\t ValISrho = " << (SolnBlk[Cell+1].*AccessToHighOrderVar)().CellSmoothnessIndicator(1) << "\n"
	   << " ISp = " << (SolnBlk[Cell+1].*AccessToHighOrderVar)().CellInadequateFit(3)
	   << "\t ValISp = " << (SolnBlk[Cell+1].*AccessToHighOrderVar)().CellSmoothnessIndicator(3) << "\n";
    }
  }

  if (Left_UnphysicalValuesDetected || Right_UnphysicalValuesDetected){

    // Force with piecewise constant solution if required.
    if ( (CENO_Execution_Mode::USE_CENO_ALGORITHM && CENO_Execution_Mode::FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE) ||
	 (ENO_Execution_Mode::USE_ENO_ALGORITHM   &&  ENO_Execution_Mode::FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE) ){
      Wl = SolnBlk[Cell].CellSolutionPrimVar();
      Wr = SolnBlk[Cell+1].CellSolutionPrimVar();

      // check positivity of the piecewise constant values
      if (Wl.d <= ZERO || Wl.p <= ZERO || Wr.d <= ZERO || Wr.p <= ZERO  ){
	// Throw an error
	throw runtime_error("Euler1D::GetRightAndLeftFluxStates() ERROR! Unphysical piecewise constant values encountered!");
      }

      // output user info
      cout << " New interface values equal to the average solution: \n" << " Wl = " << Wl << "\n Wr = " << Wr << "\n"; 

    } else {
      // Throw an error
      throw runtime_error("Euler1D::GetRightAndLeftFluxStates() ERROR! Unphysical values encountered!");
    }
  }//endif 

}

#endif
