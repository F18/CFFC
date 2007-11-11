#ifndef _NKS_DTS_INCLUDED
#define _NKS_DTS_INCLUDED

// Generic DTS 
#include "DTS.h"

/******************** TEMPLATED FUNCTIONS **************************************************************/

/*! *******************************************************************
 * LHS_Time                                                           *
 *                                                                    *
 * Left Hand Side time componenet                                     *
 *        [ I*LHS_Time + J ] = -R(Un)                                 *
 *                                                                    *
 **********************************************************************/
template  <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline double LHS_Time(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &Input, double& d_tau, const double &DTS_dTime){
  
  // Dual Time Stepping 
  if (Input.NKS_IP.Dual_Time_Stepping){
    // Implicit Euler
    if (Input.NKS_IP.Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {
      return (ONE/d_tau + ONE/DTS_dTime);      
      //BDF2
    } else if (Input.NKS_IP.Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {
      return (ONE/d_tau + THREE/(TWO*DTS_dTime));     
    }
  }
   
  //Standard
  return ONE/(d_tau);  
}

/*! *************************************************************************
 *  Overloading dUdt_Residual_Evaluation function to include Dual Time 
 *  Stepping source term.
 * 
 ****************************************************************************/
template  <typename SOLN_pSTATE, typename SOLN_cSTATE> 
int dUdt_Residual_Evaluation_DTS(HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> *Solution_Data,
				 DTS_Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> *DTS_SolnBlk,
				 const int Block){
				 
  int error_flag = Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].dUdt_Residual_Evaluation(Solution_Data->Input);

 
//   error_flag = Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].dUdt_Multistage_Explicit(1,Solution_Data->Input);
//   cout<<"\n R2 "<<Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].dUdt[2][2][2][0]/
//     (Solution_Data->Input.CFL_Number*Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].dt[2][2][2]);

  // Add dual time stepping Source Term to Residual ie. dUdt[i][j][0]
  if ( Solution_Data->Input.NKS_IP.Dual_Time_Stepping) {                 //Do these line up with dUdt_Residual_Evaluation?????   
    for (int i = Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].ICl; i <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].ICu; i++) {   
      for (int j = Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].JCl; j <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].JCu; j++) {
	for (int k = Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].KCl; k <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].KCu; k++) {
	  	  
	  //Implicit Euler  R(U_n)* = R(U_n) + (U - Un)/dt
	  if ( Solution_Data->Input.NKS_IP.Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {		  
	    Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].dUdt[i][j][k][0] -= 
	      (Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].U[i][j][k] 
	       - DTS_SolnBlk[Block].Un[i][j][k])/(DTS_SolnBlk[Block].DTS_dTime); 
	    // BDF2         R(U_n)* = R(U_n) + (3U - 4Un +  Un-1)/2dt  
	  } else if (Solution_Data->Input.NKS_IP.Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {     
	    Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].dUdt[i][j][k][0] -= 
	      (THREE*Solution_Data->Local_Solution_Blocks.Soln_Blks[Block].U[i][j][k] - FOUR*DTS_SolnBlk[Block].Un[i][j][k]
	       + DTS_SolnBlk[Block].Unminus1[i][j][k])/(TWO*DTS_SolnBlk[Block].DTS_dTime);
	  } 
	} 	      
      }
    }
  }  
     
  return error_flag;
}

#endif // _NKS_DTS_INCLUDED

