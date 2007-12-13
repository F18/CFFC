/*! \file dUdt_HighOrder1D.h
  \brief Different time marching template subroutines for solving 1D problems with high-order spatial accuracy. */

#ifndef _DUDT_HIGHORDER_1D_
#define _DUDT_HIGHORDER_1D_

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
// None

// SolutionClass Traits;
template<class SOLN_pSTATE>
class SolutionClassTraits;

/******************************************************//**
 * Routine: dUdt_explicitEuler_upwind                   
 *                                                      
 * This routine updates the solution using a 1st-order  
 * explicit Euler time integration and HIGH-ORDER upwind
 * spatial discretization scheme in conjunction with    
 * either the Godunov, Roe, Rusanov, HLLE, Linde, or    
 * HLLC flux functions.                                 
 *                                                      
 ********************************************************/
template<class Soln_Block_Type>
int dUdt_explicitEuler_upwind(Soln_Block_Type *SolnBlk,
			      const CFD1D_Input_Parameters &IP,
			      double &dtMin,
			      const int Local_Time_Stepping,
			      typename Soln_Block_Type::HighOrderType & 
			      (Soln_Block_Type::*AccessToHighOrderVar)(void) =
			      &Soln_Block_Type::CellHighOrder) {

  typedef typename Soln_Block_Type::SOLN_pSTATE Soln_pState;
  typedef typename Soln_Block_Type::SOLN_cSTATE Soln_cState;

  int i;
  Soln_pState Wl, Wr;
  Soln_cState Flux, Soln_U_VACUUM(0.0);

  int ICl(SolnBlk[0].ICl), ICu(SolnBlk[0].ICu);
    
  /*************************************************************
   * Apply the BCs to update the ghost cell average solutions. *
   ************************************************************/
  BCs(SolnBlk,IP,AccessToHighOrderVar);

  // High-order Reconstruction in the whole computational domain
  HighOrderSolutionReconstructionOverDomain(SolnBlk,IP,AccessToHighOrderVar);
 
  /* Evaluate the time rate of change of the solution
     (i.e., the solution residuals) using the first-order
     upwind scheme with a variety of flux functions. */

  if (!Local_Time_Stepping) SolnBlk[ICl-1].dt = dtMin;
  SolnBlk[ICl-1].dUdt = Soln_U_VACUUM;

    
  for ( i = ICl-1; i <= ICu ; ++i ) {
    if ( !Local_Time_Stepping) SolnBlk[i+1].dt = dtMin;

    /* Evaluate the cell interface flux. */
    GetRightAndLeftFluxStates(Wl,Wr,i,SolnBlk,IP,AccessToHighOrderVar);

    SolnBlk[i+1].dUdt = Soln_U_VACUUM;

    /********* FLUX EVALUATION ***********/
    Flux = RiemannFlux(IP.i_Flux_Function,Wl,Wr);

    /* Evaluate cell-averaged solution changes. */
    SolnBlk[i].dUdt -= Flux/SolnBlk[i].X.dx;
    SolnBlk[i+1].dUdt += Flux/SolnBlk[i+1].X.dx;
  } /* endfor */

  /* Update both conserved and primitive solution
     variables using explicit Euler method. */
  for ( i = ICl ; i <= ICu ; ++i ) {
    if ( !Local_Time_Stepping ) SolnBlk[i].dt = dtMin;
	
    SolnBlk[i].U += (IP.CFL_Number*SolnBlk[i].dt)*SolnBlk[i].dUdt;
	
    /* Check for the positivity of the updated solution.
     * Specialize this function if checks different than "d" and "E" must be performed
     */
    CheckSolutionPositivity(SolnBlk,i);

    /* Update the primitive variable solution state. */
    SolnBlk[i].W = W(SolnBlk[i].U);
  } /* endfor */

  /* set dUdt to Zero at the boundaries for the next step */
  SolnBlk[ICl-1].dUdt = Soln_U_VACUUM;
  SolnBlk[ICu+1].dUdt = Soln_U_VACUUM;

  /* Solution successfully updated. */
  return (0);
}

/******************************************************//**
 * Routine: dUdt_2stage_HighOrder_upwind                
 *                                                      
 * This routine updates the solution using a two-stage  
 * second-order explicit time integration scheme        
 * and a high-order upwind spatial discretization 
 * scheme with either the Godunov, Roe,  
 * Rusanov, HLLE, Linde, or HLLC flux functions.        
 *                                                      
 ********************************************************/
template<class Soln_Block_Type>
int dUdt_2stage_HighOrder_upwind(Soln_Block_Type *SolnBlk,
				 const CFD1D_Input_Parameters &IP,
				 double &dtMin,
				 const int Local_Time_Stepping,
				 typename Soln_Block_Type::HighOrderType & 
				 (Soln_Block_Type::*AccessToHighOrderVar)(void) =
				 &Soln_Block_Type::CellHighOrder) {
  
  typedef typename Soln_Block_Type::SOLN_pSTATE Soln_pState;
  typedef typename Soln_Block_Type::SOLN_cSTATE Soln_cState;

  int i, n_stage;
  double omega;
  Soln_pState Wl, Wr;
  Soln_cState Flux, Soln_U_VACUUM(0.0);

  int ICl(SolnBlk[0].ICl), ICu(SolnBlk[0].ICu);
  
  /* Perform second-order two-stage semi-implicit update of solution
     varibles for new time level. */

  for ( n_stage = 1 ; n_stage <= 2 ; ++n_stage ){

    /* Evaluate the time step fraction for stage. */
    omega = ONE/double(n_stage);

    /**************************************************************************
     *      Apply the BC to update the average solution at boundaries.        *
     * This update is used in the reconstruction process of the current stage *
     *************************************************************************/
    BCs(SolnBlk,IP,AccessToHighOrderVar);

    // High-order Reconstruction in the whole computational domain
    HighOrderSolutionReconstructionOverDomain(SolnBlk,IP,AccessToHighOrderVar);

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a higher-order
       upwind scheme with a variety of flux functions. */

    if ( !Local_Time_Stepping && n_stage == 1 ) SolnBlk[ICl-1].dt = dtMin;
    if ( n_stage == 1 ){
      SolnBlk[ICl-1].Uo = SolnBlk[ICl-1].U;
    }
    SolnBlk[ICl-1].dUdt = Soln_U_VACUUM;

    for (i =ICl-1 ; i<=ICu ;++i){
      if ( !Local_Time_Stepping && n_stage == 1 ) SolnBlk[i+1].dt = dtMin;
      if ( n_stage == 1 ) {
	SolnBlk[i+1].Uo = SolnBlk[i+1].U;
	SolnBlk[i+1].dUdt = Soln_U_VACUUM;
      } else {
	SolnBlk[i+1].dUdt = HALF*SolnBlk[i+1].dUdt;
      } /* endif */ 

      /* Evaluate the cell interface flux. */
      GetRightAndLeftFluxStates(Wl,Wr,i,SolnBlk,IP,AccessToHighOrderVar);

      /********* FLUX EVALUATION ***********/
      Flux = RiemannFlux(IP.i_Flux_Function,Wl,Wr);

      /* Evaluate cell-averaged solution changes. */
      SolnBlk[i].dUdt -= (omega*IP.CFL_Number*SolnBlk[i].dt)*Flux/SolnBlk[i].X.dx;
      SolnBlk[i+1].dUdt += (omega*IP.CFL_Number*SolnBlk[i+1].dt)*Flux/SolnBlk[i+1].X.dx;

    } /* endfor */

    /* Update solution variables for this stage. */
    for (i = ICl ;i<= ICu ;++i) {
      SolnBlk[i].U = SolnBlk[i].Uo + SolnBlk[i].dUdt;

      /* Check for the positivity of the updated solution.
       * Specialize this function if checks different than "d" and "E" must be performed
       */
      CheckSolutionPositivity(SolnBlk,i);

      /* Update the primitive variable solution state. */
      SolnBlk[i].W = W(SolnBlk[i].U);
    } /* endfor */

    /* set dUdt to Zero at the boundaries for the next step */
    SolnBlk[ICl-1].dUdt = Soln_U_VACUUM;
    SolnBlk[ICu+1].dUdt = Soln_U_VACUUM;
    
  } /* endfor - stage */

  /* Solution successfully updated. */  
  return (0);   
}

/******************************************************//**
 * Routine: dUdt_4stage_HighOrder_upwind                
 *                                                      
 * This routine updates the solution using a four-stage 
 * fourth-order explicit time integration scheme        
 * and a high-order data-dependent upwind spatial       
 * discretization scheme with either the Godunov, Roe,  
 * Rusanov, HLLE, Linde, or HLLC flux functions.        
 *                                                      
 ********************************************************/
template<class Soln_Block_Type>
int dUdt_4stage_HighOrder_upwind(Soln_Block_Type *SolnBlk,
				 const CFD1D_Input_Parameters &IP,
				 double &dtMin,
				 const int Local_Time_Stepping,
				 typename Soln_Block_Type::HighOrderType & 
				 (Soln_Block_Type::*AccessToHighOrderVar)(void) =
				 &Soln_Block_Type::CellHighOrder){


  typedef typename Soln_Block_Type::SOLN_pSTATE Soln_pState;
  typedef typename Soln_Block_Type::SOLN_cSTATE Soln_cState;

  int i, n_stage;
  double omega, beta;
  Soln_pState Wl, Wr;
  Soln_cState Flux, Soln_U_VACUUM(0.0);

  int ICl(SolnBlk[0].ICl), ICu(SolnBlk[0].ICu);

  /* Perform fourth-order four-stage semi-implicit update of solution
     varibles for new time level. */
  
  for ( n_stage = 1 ; n_stage <= 4 ; ++n_stage ){
    /* Evaluate the time step fraction and the residual coefficient for stage. */

    switch(n_stage){
    case 1:
      omega = ONE/TWO;
      beta = ONE/THREE;
      break;
    case 2:
      omega = ONE/TWO;
      beta = TWO/THREE;
      break;
    case 3:
      omega = ONE;
      beta = ONE/THREE;
      break;
    case 4:
      omega = ONE;
      beta = ONE/SIX;
      break;
    }

    /**************************************************************************
     * Apply the BCs to update the ghost cell average solutions.              *
     * This update is used in the reconstruction process of the current stage *
     *************************************************************************/
    BCs(SolnBlk,IP,AccessToHighOrderVar);
    
    // High-order Reconstruction in the whole computational domain
    HighOrderSolutionReconstructionOverDomain(SolnBlk,IP,AccessToHighOrderVar);

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a higher-order
       data-dependent upwind scheme with a variety of flux functions. */
    
    if ( n_stage == 1 ){
      SolnBlk[ICl-1].dt = dtMin;
      SolnBlk[ICl-1].Uo = SolnBlk[ICl-1].U;
      SolnBlk[ICl-1].TotaldUdt = Soln_U_VACUUM;
      SolnBlk[ICl-1].dUdt = Soln_U_VACUUM;
    }

    for (i =ICl-1 ; i<=ICu ;++i){
      if ( n_stage == 1 ) {
	SolnBlk[i+1].dt = dtMin;
	SolnBlk[i+1].Uo = SolnBlk[i+1].U;
	// set the total solution residual to 0
	SolnBlk[i+1].TotaldUdt = Soln_U_VACUUM;
	// set the solution residual to 0
	SolnBlk[i+1].dUdt = Soln_U_VACUUM;
      } /* endif */

      /* Evaluate the cell interface flux. */
      GetRightAndLeftFluxStates(Wl,Wr,i,SolnBlk,IP,AccessToHighOrderVar);

      /********* FLUX EVALUATION ***********/	
      Flux = RiemannFlux(IP.i_Flux_Function,Wl,Wr);

      /* Evaluate cell-averaged solution changes. */
      SolnBlk[i].dUdt -= (omega*IP.CFL_Number*SolnBlk[i].dt)*Flux/SolnBlk[i].X.dx;
      SolnBlk[i+1].dUdt += (omega*IP.CFL_Number*SolnBlk[i+1].dt)*Flux/SolnBlk[i+1].X.dx;

    } /* endfor */

    /* Update solution variables for this stage. */
    for (i = ICl ;i<= ICu ;++i) {

      /* Update the total flux with that of the stage */
      SolnBlk[i].TotaldUdt += beta*SolnBlk[i].dUdt;

      if (n_stage==4)
	SolnBlk[i].U = SolnBlk[i].Uo + SolnBlk[i].TotaldUdt;
      else{
	SolnBlk[i].U = SolnBlk[i].Uo + SolnBlk[i].dUdt;
      }

      /* Check for the positivity of the updated solution.
       * Specialize this function if checks different than "d" and "E" must be performed
       */
      CheckSolutionPositivity(SolnBlk,i);
      
      /* set dUdt to Zero for the next step */
      SolnBlk[i].dUdt = Soln_U_VACUUM;

      /* Update the primitive variable solution state. */
      SolnBlk[i].W = W(SolnBlk[i].U);
    } /* endfor */

    /* set dUdt to Zero at the boundaries for the next step */
    SolnBlk[ICl-1].dUdt = Soln_U_VACUUM;
    SolnBlk[ICu+1].dUdt = Soln_U_VACUUM;

  } /* endfor - stage */

  /* Solution successfully updated. */

  return (0);   
}

/******************************************************//**
 * Routine: BCs
 *                                                      
 * This routine imposes the boundary conditions by setting 
 * the average solution values in the ghost cells.
 * Only constant extrapolation boundary conditions are 
 * imposed with this subroutine.
 * For more diverse BCs, a specialization of this routine
 * is required.                                             
 ********************************************************/
template<class Soln_Block_Type> inline
void BCs(Soln_Block_Type *SolnBlk,
	 const CFD1D_Input_Parameters &IP,
	 typename Soln_Block_Type::HighOrderType & 
	 (Soln_Block_Type::*AccessToHighOrderVar)(void) =
	 &Soln_Block_Type::CellHighOrder){

  /* By default, constant extrapolation boundary
     conditions are applied at either end of the mesh. */
  int ICl(SolnBlk[0].ICl);
  int ICu(SolnBlk[0].ICu);
  
  for (int i = 0; i<ICl; ++i){
    // left end
    SolnBlk[i].U = SolnBlk[ICl].U;
    SolnBlk[i].W = SolnBlk[ICl].W;
    // right end
    SolnBlk[ICu+i+1].U = SolnBlk[ICu].U;
    SolnBlk[ICu+i+1].W = SolnBlk[ICu].W;
  }
}

/******************************************************//**
 * Routine: RiemannFlux
 *                                                      
 * This routine calculates the conserved flux through an
 * interface of two adjacent cells by solving a Riemann
 * problem.
 * \param [in] Flux_Function the method used to calculate the flux
 * \param [in] Wl the solution to the left of the interface
 * \param [in] Wr the solution to the right of the interface
 * \return the conserved flux through the interface
 ********************************************************/
template<class SOLN_pSTATE> inline
typename SolutionClassTraits<SOLN_pSTATE>::SOLN_cSTATE RiemannFlux(const int & Flux_Function,
								   SOLN_pSTATE &Wl,
								   SOLN_pSTATE &Wr){
  switch(Flux_Function) {
  case FLUX_FUNCTION_GODUNOV :
    return FluxGodunov(Wl, Wr);
  case FLUX_FUNCTION_ROE :
    return FluxRoe(Wl, Wr);
  case FLUX_FUNCTION_RUSANOV :
    return FluxRusanov(Wl, Wr);
  case FLUX_FUNCTION_HLLE :
    return FluxHLLE(Wl, Wr);
  case FLUX_FUNCTION_LINDE :
    return FluxLinde(Wl, Wr);
  case FLUX_FUNCTION_HLLC :
    return FluxHLLC(Wl, Wr);
  case FLUX_FUNCTION_OSHER :
    return FluxOsher(Wl, Wr);
  default:
    throw runtime_error("RiemannFlux() ERROR: Unknown flux function type!");
  } /* endswitch */
}


/******************************************************//**
 * Routine: CheckSolutionPositivity
 *                                                      
 * This routine checks whether the solution of cell "iCell"
 * is still physical (i.e. require some variables to be positive)
 * The current implementation checks for density and energy, but
 * specializations can be provided if the required checks are
 * different than this general one.
 * 
 * \param [in] SolnBlk the solution block
 * \param [in] iCell the cell which is checked for positivity
 * \throw throws runtime_error if non-positivity is detected.
 * 
 ********************************************************/
template<class Soln_Block_Type> inline
void CheckSolutionPositivity(const Soln_Block_Type *SolnBlk,
			     const int & iCell){
  if (SolnBlk[iCell].U.d   <= ZERO ||
      SolnBlk[iCell].U.E   <= ZERO ||
      SolnBlk[iCell].U.e() <= ZERO ) {
    cout << "\n " << CFFC_Name() << " ERROR: Negative Density and/or Energy: \n"
	 << " node = " << iCell 
	 << "\n U = " << SolnBlk[iCell].U 
	 << "\n dUdt = "  << SolnBlk[iCell].dUdt << "\n";
    cout.flush();

    throw runtime_error("CFFC ERROR: Solution update error");
  }
}

/***************************************************************//**
 * Routine: GetRightAndLeftFluxStates                            
 *
 * Calculates the right and left solution states for an Riemann problem
 * solved at the interface between "Cell" and "Cell+1" cells 
 * using a high-order variable.
 *
 * \param [out] Wl the solution to the left of the interface
 * \param [out] Wr the solution to the right of the interface
 * \param [in] Cell the identification of the interface
 * \param [in] SolnBlk the solution block
 * \param [in] IP the input parameters
 * \param [in] AccessToHighOrderVar member function of SolnBlk that returns the high-order variable
 *                                  used for the current computation
 ***************************************************************************************************/
template<class Soln_Block_Type> inline
void GetRightAndLeftFluxStates(typename Soln_Block_Type::SOLN_pSTATE &Wl,
			       typename Soln_Block_Type::SOLN_pSTATE &Wr,
			       const int &Cell,
			       Soln_Block_Type *SolnBlk,
			       const CFD1D_Input_Parameters &IP,
			       typename Soln_Block_Type::HighOrderType & 
			       (Soln_Block_Type::*AccessToHighOrderVar)(void) =
			       &Soln_Block_Type::CellHighOrder){

  Wl = (SolnBlk[Cell].*AccessToHighOrderVar)().right_state();
  Wr = (SolnBlk[Cell+1].*AccessToHighOrderVar)().left_state();
}

#endif
