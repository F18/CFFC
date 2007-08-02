/* FANS3DThermallyPerfectHexaMultiBlock.h:  Header file creating multiblock list. */

#ifndef _FANS3D_THERMALLYPERFECT_HEXA_MULTIBLOCK_INCLUDED
#define _FANS3D_THERMALLYPERFECT_HEXA_MULTIBLOCK_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_MULTIBLOCK_INCLUDED
#include "../HexaBlock/HexaMultiBlock.h"
#endif //_HEXA_MULTIBLOCK_INCLUDEDT

#ifndef _FANS3D_THERMALLYPERFECT_STATE_INCLUDED
#include "FANS3DYhermallyPerfectState.h"
#endif // FANS3D_THERMALLYPERFECT_STATE_INCLUDED   


/* template< >  */
/* int Hexa_MultiBlock<Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,   */
/*    FANS3D_ThermallyPerfect_KOmega_cState> >::Wall_Shear(void){  */
   
/*    for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){ */
/*       if(Block_Used[nblk]){ */
            
/*          int i, j, k; */
/*          for ( k =  Hexa_Block_List[nblk]->Grid->KCl- Hexa_Block_List[nblk]->Grid->Nghost ; k <=  Hexa_Block_List[nblk]->Grid->KCu+ Hexa_Block_List[nblk]->Grid->Nghost ; ++k ) */
/*             for ( j =  Hexa_Block_List[nblk]->Grid->JCl- Hexa_Block_List[nblk]->Grid->Nghost ; j <=  Hexa_Block_List[nblk]->Grid->JCu+ Hexa_Block_List[nblk]->Grid->Nghost ; ++j )  */
/*                for ( i =  Hexa_Block_List[nblk]->Grid->ICl- Hexa_Block_List[nblk]->Grid->Nghost ; i <=  Hexa_Block_List[nblk]->Grid->ICu+ Hexa_Block_List[nblk]->Grid->Nghost ; ++i ) { */
            
/*             // for ghost cells close to solid wall, the ywall is set to be zero, */
/*             // this is related to the application of wall function (see the logic in the wall fucntion code */
/*             // it may be not necessary, since the noslip boundary condition, tangiential velocity */
/*             // is zero, so the wall function will most likely to be applied correctly. */
/*             // however, it is not over-excessive to have these boundary checks. */
            
/*             // Check West boundary. */
/*             if (Hexa_Block_List[nblk]->Grid->BCtypeW[j][k] == BC_WALL_VISCOUS || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeW[j][k] == BC_NO_SLIP  || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeW[j][k] == BC_WALL_VISCOUS_ISOTHERMAL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeW[j][k] == BC_WALL_VISCOUS_HEATFLUX || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeW[j][k] == BC_ADIABATIC_WALL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeW[j][k] == BC_BURNING_SURFACE){ */
               
/*                Hexa_Block_List[nblk]->Wall[ Hexa_Block_List[nblk]->Grid->ICl-2][j][k].ywall = ZERO; */
/*                Hexa_Block_List[nblk]->Wall[ Hexa_Block_List[nblk]->Grid->ICl-1][j][k].ywall = ZERO; */
/*             } /\* endif *\/ */
            
/*             // Check East boundary. */
/*             if (Hexa_Block_List[nblk]->Grid->BCtypeE[j][k] == BC_WALL_VISCOUS || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeE[j][k] == BC_NO_SLIP  || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeE[j][k] == BC_WALL_VISCOUS_ISOTHERMAL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeE[j][k] == BC_WALL_VISCOUS_HEATFLUX || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeE[j][k] == BC_ADIABATIC_WALL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeE[j][k] == BC_BURNING_SURFACE) { */
/*                Hexa_Block_List[nblk]->Wall[ Hexa_Block_List[nblk]->Grid->ICu+1][j][k].ywall = ZERO; */
/*                Hexa_Block_List[nblk]->Wall[ Hexa_Block_List[nblk]->Grid->ICu+2][j][k].ywall = ZERO; */
               
/*             } /\* endif *\/ */
            
/*             // Check South boundary. */
/*             if (Hexa_Block_List[nblk]->Grid->BCtypeS[i][k] == BC_WALL_VISCOUS || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeS[i][k] == BC_NO_SLIP  || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeS[i][k] == BC_WALL_VISCOUS_ISOTHERMAL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeS[i][k] == BC_WALL_VISCOUS_HEATFLUX || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeS[i][k] == BC_ADIABATIC_WALL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeS[i][k] == BC_BURNING_SURFACE) { */
               
/*                Hexa_Block_List[nblk]->Wall[i][ Hexa_Block_List[nblk]->Grid->JCl-2][k].ywall = ZERO; */
/*                Hexa_Block_List[nblk]->Wall[i][ Hexa_Block_List[nblk]->Grid->JCl-1][k].ywall = ZERO; */
/*             } /\* endif *\/ */
               
               
/*             // Check North boundary. */
/*             if (Hexa_Block_List[nblk]->Grid->BCtypeN[i][k] == BC_WALL_VISCOUS || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeN[i][k] == BC_NO_SLIP  || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeN[i][k] == BC_WALL_VISCOUS_ISOTHERMAL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeN[i][k] == BC_WALL_VISCOUS_HEATFLUX || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeN[i][k] == BC_ADIABATIC_WALL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeN[i][k] == BC_BURNING_SURFACE) { */
               
/*                Hexa_Block_List[nblk]->Wall[i][ Hexa_Block_List[nblk]->Grid->JCu+1][k].ywall = ZERO; */
/*                Hexa_Block_List[nblk]->Wall[i][ Hexa_Block_List[nblk]->Grid->JCu+2][k].ywall = ZERO; */
                  
/*             } /\* endif *\/ */
               
/*             // Check Bottom boundary. */
            
/*             if (Hexa_Block_List[nblk]->Grid->BCtypeB[i][j] == BC_WALL_VISCOUS || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeB[i][j] == BC_NO_SLIP  || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeB[i][j] == BC_WALL_VISCOUS_ISOTHERMAL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeB[i][j] == BC_WALL_VISCOUS_HEATFLUX || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeB[i][j] == BC_ADIABATIC_WALL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeB[i][j] == BC_BURNING_SURFACE) { */
               
/*                Hexa_Block_List[nblk]->Wall[i][j][ Hexa_Block_List[nblk]->Grid->KCl-1].ywall = ZERO; */
/*                Hexa_Block_List[nblk]->Wall[i][j][ Hexa_Block_List[nblk]->Grid->KCl-2].ywall = ZERO; */
/*             } */
            
/*             // Check Top boundary. */
/*             if (Hexa_Block_List[nblk]->Grid->BCtypeT[i][j] == BC_WALL_VISCOUS || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeT[i][j] == BC_NO_SLIP  || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeT[i][j] == BC_WALL_VISCOUS_ISOTHERMAL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeT[i][j] == BC_WALL_VISCOUS_HEATFLUX || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeT[i][j] == BC_ADIABATIC_WALL || */
/*                 Hexa_Block_List[nblk]->Grid->BCtypeT[i][j] == BC_BURNING_SURFACE) { */
               
/*                Hexa_Block_List[nblk]->Wall[i][j][ Hexa_Block_List[nblk]->Grid->KCu+1].ywall = ZERO; */
/*                Hexa_Block_List[nblk]->Wall[i][j][ Hexa_Block_List[nblk]->Grid->KCu+2].ywall = ZERO; */
               
/*             } */
               
/*             Hexa_Block_List[nblk]->Wall[i][j][k].utau =  Hexa_Block_List[nblk]->Wall_Friction_Velocity(i, j, k); */
/*             Hexa_Block_List[nblk]->Wall[i][j][k].tauw =  Hexa_Block_List[nblk]->W[i][j][k].rho * sqr(Hexa_Block_List[nblk]->Wall[i][j][k].utau); */
/*             Hexa_Block_List[nblk]->Wall[i][j][k].yplus =  Hexa_Block_List[nblk]->Wall[i][j][k].utau*Hexa_Block_List[nblk]->Wall[i][j][k].ywall/ */
/*                ( Hexa_Block_List[nblk]->W[i][j][k].mu()/ Hexa_Block_List[nblk]->W[i][j][k].rho); */
            
                             
/*                } /\* endfor *\/ */

   
/*       } */
      
/*    } */
   
   
/*    return 0; */
      
  
    
/* } */



/* template<class HEXA_BLOCK> */
/* int Hexa_MultiBlock<HEXA_BLOCK<Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,  */
/*    FANS3D_ThermallyPerfect_KOmega_cState> > >::Create_Initial_Solution(void){ */
   
/*    if( Hexa_Block_List[0]->Flow_Type != FLOWTYPE_TURBULENT_RANS_K_OMEGA){ */
      
/*       return (0); */
      
/*    }else{ */
      
/*       int error_flag; */
/*       double y_wall_temp; */
/*       Vector3D X_wall_temp, n_wall_temp; */
/*       int BC_wall_temp; */
      
/*       y_wall_temp = std::numeric_limits<double>::max(); */
/*       for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){ */
/*          if(Block_Used[nblk]){ */
            
/*             // compute wall data ... ...  */
/*             for ( int k = Hexa_Block_List[nblk]->Grid->KCl- Hexa_Block_List[nblk]->Grid->Nghost ; k <=   Hexa_Block_List[nblk]->Grid->KCu+  Hexa_Block_List[nblk]->Grid->Nghost ; ++k ) */
/*                for ( int j = Hexa_Block_List[nblk]->Grid->JCl- Hexa_Block_List[nblk]->Grid->Nghost ; j <=   Hexa_Block_List[nblk]->Grid->JCu+  Hexa_Block_List[nblk]->Grid->Nghost ; ++j )  */
/*                   for ( int i = Hexa_Block_List[nblk]->Grid->ICl-  Hexa_Block_List[nblk]->Grid->Nghost ; i <=   Hexa_Block_List[nblk]->Grid->ICu+  Hexa_Block_List[nblk]->Grid->Nghost ; ++i ) { */
/*                      Hexa_Block_List[nblk]->Wall[i][j][k].ywall =  std::numeric_limits<double>::max(); */
/*                      error_flag = Wall_Distance(Hexa_Block_List[nblk], */
/*                                                 Hexa_Block_List[nblk]->Grid->Cell[i][j][k].Xc, */
/*                                                 X_wall_temp, n_wall_temp, */
/*                                                 y_wall_temp, BC_wall_temp); */
/*                      if (y_wall_temp < Hexa_Block_List[nblk]->Wall[i][j][k].ywall ) { */
/*                         Hexa_Block_List[nblk]->Wall[i][j][k].ywall = y_wall_temp; */
/*                         Hexa_Block_List[nblk]->Wall[i][j][k].Xwall = X_wall_temp; */
/*                         Hexa_Block_List[nblk]->Wall[i][j][k].nwall = n_wall_temp; */
/*                         Hexa_Block_List[nblk]->Wall[i][j][k].BCwall = BC_wall_temp; */
/*                      } */
                                     
/*                   } */
/*             // compute y+ for each cell */
/*             Wall_Shear(Hexa_Block_List[nblk]); */
            
/*          }// compute for used blocks */
/*       }//loop through all blocks */
/*    }//turbulent case */
   

/*    return 0; */
/* } */

#endif // _FANS_THERMALLYPERFECT_HEXA_MULTIBLOCK_INCLUDED
