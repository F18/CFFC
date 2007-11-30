
#ifndef  _NAVIERSTOKES3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED
#include "NavierStokes3DThermallyPerfectHexaBlock.h"
#endif //  _NAVIERSTOKES3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

template<>
int Hexa_Block<NavierStokes3D_ThermallyPerfect_pState, 
               NavierStokes3D_ThermallyPerfect_cState>::dUdt_Multistage_Explicit(
                  const int i_stage,
                  Input_Parameters<NavierStokes3D_ThermallyPerfect_pState, 
                  NavierStokes3D_ThermallyPerfect_cState> &IPs) {
    
   int i, j, k,  k_residual;
   double omega;
   Vector3D dX;
   
   NavierStokes3D_ThermallyPerfect_pState Wl, Wr;
   NavierStokes3D_ThermallyPerfect_cState Flux, Temp;

   NavierStokes3D_ThermallyPerfect_cState U_VACUUM;
   U_VACUUM.Vacuum();
   NavierStokes3D_ThermallyPerfect_pState W_VACUUM;
   W_VACUUM.Vacuum();

   /* Evaluate the solution residual for stage 
    *       i_stage of an N stage scheme. */
   /* Evaluate the time step fraction and residual storage location 
    *       for the stage. */
   switch(IPs.i_Time_Integration) {
   case TIME_STEPPING_EXPLICIT_EULER :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      
      break;
   case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      break;
   case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0; 
      if (IPs.N_Stage == 4) {
         if (i_stage == 4) {
            k_residual = 0;
         } else {
            k_residual = i_stage - 1;
         } /* endif */
      } /* endif */
      break;
   case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
      omega = MultiStage_Optimally_Smoothing(i_stage, 
                                             IPs.N_Stage,
                                             IPs.i_Limiter);
      k_residual = 0;
      break;
   default: 
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      break;
   } /* endswitch */
   
      
   /* Perform the linear reconstruction within each cell
      of the computational grid for this stage. */
   
   switch(IPs.i_Reconstruction) {
      
   case RECONSTRUCTION_LEAST_SQUARES :
      
      Linear_Reconstruction_LeastSquares(IPs.i_Limiter);

      break;
   default:
      Linear_Reconstruction_LeastSquares(IPs.i_Limiter);
      break;
   } /* endswitch */

    /* Evaluate the time rate of change of the solution
      (i.e., the solution residuals) using a second-order
      limited upwind scheme with a variety of flux functions. */
    
   // Add i-direction (zeta-direction) fluxes.
   for ( k  =  KCl-1 ; k <=  KCu+1 ; ++k ){
      for ( j  =  JCl-1 ; j <=  JCu+1 ; ++j ) {
         if ( i_stage == 1 ) {
             Uo[ ICl-1][j][k] =  U[ ICl-1][j][k];
             dUdt[ ICl-1][j][k][k_residual] = U_VACUUM;
      
             
         } else {
             dUdt[ ICl-1][j][k][k_residual] = U_VACUUM;
         } /* endif */
         
         
         for ( i =  ICl-1 ; i <=  ICu ; ++i ) {
            if ( i_stage == 1 ) {
               Uo[i+1][j][k] =  U[i+1][j][k];
               
               dUdt[i+1][j][k][k_residual] = U_VACUUM;
            } else if (( j >  JCl-1 && j <  JCu+1 ) 
                       && (k >  KCl-1 && k <  KCu+1) ){
               switch(IPs.i_Time_Integration) {

               case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
                  //dUdt[i+1][j][k][k_residual] = 
                  //  dUdt[i+1][j][k][k_residual];
                  break;
               case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
		  if (IPs.N_Stage == 2) {
                     //dUdt[i+1][j][k_residual] = 
                     // dUdt[i+1][j][k_residual];
                  } else if (IPs.N_Stage == 4 && i_stage == 4) {
                     dUdt[i+1][j][k][k_residual] = 
                        dUdt[i+1][j][k][0] + 
                        TWO*dUdt[i+1][j][k][1] +TWO*dUdt[i+1][j][k][2];
                  } else {
         	     dUdt[i+1][j][k][k_residual].Vacuum();
                  } /* endif */
                  break;
               case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
                  dUdt[i+1][j][k][k_residual] = U_VACUUM;
                  break;
               default:
                  dUdt[i+1][j][k][k_residual] = U_VACUUM;
                  break;
               } /* endswitch */
            } /* endif */
            if (( j >  JCl-1 && j <  JCu+1 ) && 
                ( k >  KCl-1 && k <  KCu+1 )) {
            
            
               
               /* Evaluate the cell interface i-direction fluxes. */
               if (i ==  ICl-1 &&  (Grid.BCtypeW[j][k] == BC_REFLECTION ||
                                    Grid.BCtypeW[j][k] == BC_NO_SLIP||
                                    Grid.BCtypeW[j][k] == BC_MOVING_WALL)  ) {
                
                  
                  
                  dX =  Grid.xfaceW(i+1, j, k)- Grid.Cell[i+1][j][k].Xc;
                  
                  Wr =  W[i+1][j][k] + 
                     ( phi[i+1][j][k]^dWdx[i+1][j][k])*dX.x +
                     ( phi[i+1][j][k]^dWdy[i+1][j][k])*dX.y +
                     ( phi[i+1][j][k]^dWdz[i+1][j][k])*dX.z;
                  
                  if ( Grid.BCtypeW[j][k] == BC_REFLECTION) {
                     Wl =  NavierStokes3D_ThermallyPerfect_pState::Reflect(Wr,  Grid.nfaceW(i+1, j, k));
                     
                     
                  }
                  if ( Grid.BCtypeW[j][k] == BC_NO_SLIP) {
                     Wl =  NavierStokes3D_ThermallyPerfect_pState::No_Slip(Wr, WoW[j][k], Grid.nfaceW(i+1, j, k), 
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
                     
                  }
                  if ( Grid.BCtypeW[j][k] == BC_MOVING_WALL) {
                     Wl =  NavierStokes3D_ThermallyPerfect_pState::Moving_Wall(Wr, WoW[j][k], Grid.nfaceW(i+1, j, k),
                                                                               IPs.Moving_Wall_Velocity,
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
                 
                  }
               } else if (i ==  ICu && 
                          ( Grid.BCtypeE[j][k] == BC_REFLECTION||
                            Grid.BCtypeE[j][k] == BC_NO_SLIP||
                            Grid.BCtypeE[j][k] == BC_MOVING_WALL)) {
                  
                  dX =  Grid.xfaceE(i, j, k)- Grid.Cell[i][j][k].Xc;
                  Wl =  W[i][j][k] + 
                     ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                     ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                     ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                  if ( Grid.BCtypeE[j][k] == BC_REFLECTION) {
                     Wr =  NavierStokes3D_ThermallyPerfect_pState::Reflect(Wl,  Grid.nfaceE(i, j, k));
                  } 
                  if ( Grid.BCtypeE[j][k] == BC_NO_SLIP) {
                     Wr =  NavierStokes3D_ThermallyPerfect_pState::No_Slip(Wl, WoE[j][k], Grid.nfaceE(i, j, k),  
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
                  } 
                  if ( Grid.BCtypeE[j][k] == BC_MOVING_WALL) {
                     Wr =  NavierStokes3D_ThermallyPerfect_pState::Moving_Wall(Wl, WoE[j][k], Grid.nfaceE(i, j, k),  
                                                                               IPs.Moving_Wall_Velocity,
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
                  } 
               } else {            
                  dX =  Grid.xfaceE(i, j, k)- Grid.Cell[i][j][k].Xc;
                  Wl =  W[i][j][k] + 
                     ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                     ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                     ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                  dX =  Grid.xfaceW(i+1, j, k)- Grid.Cell[i+1][j][k].Xc;
                  Wr =  W[i+1][j][k] + 
                     ( phi[i+1][j][k]^ dWdx[i+1][j][k])*dX.x +
                     ( phi[i+1][j][k]^ dWdy[i+1][j][k])*dX.y +
                     ( phi[i+1][j][k]^ dWdz[i+1][j][k])*dX.z;
               } /* endif */
               
                        
               switch(IPs.i_Flux_Function) {
               case FLUX_FUNCTION_HLLE :
                  
                  Flux =  NavierStokes3D_ThermallyPerfect_pState::FluxHLLE_n(Wl, Wr,  Grid.nfaceE(i, j, k));
                  
                
                  break;
               case FLUX_FUNCTION_ROE :
                      
                  Flux =  NavierStokes3D_ThermallyPerfect_pState::FluxRoe_n(Wl, Wr,  Grid.nfaceE(i, j, k));
                     
                  break;
                  
               } /* endswitch */
               
               // add viscous flux in i direction
               Flux -=  NavierStokes3D_ThermallyPerfect_pState::FluxViscous_n(Wl, Wr,W[i][j][k], W[i+1][j][k], 
                                     dWdx[i][j][k], dWdy[i][j][k],dWdz[i][j][k],
                                     dWdx[i+1][j][k], dWdy[i+1][j][k],dWdz[i+1][j][k],
                                     Grid.nfaceE(i, j, k), Grid.Voe(i, j, k),  
                                     Grid.delta_oe(i, j, k),
                                     Grid.volume(i, j, k), 
                                     Grid.volume(i+1, j, k));
               
               /* Evaluate cell-averaged solution changes. */
               
                dUdt[i][j][k][k_residual] -= 
                  (IPs.CFL_Number* dt[i][j][k])*
                  Flux* Grid.AfaceE(i, j, k)/
                   Grid.volume(i, j, k);
                 dUdt[i+1][j][k][k_residual] += 
                   (IPs.CFL_Number* dt[i+1][j][k])*
                   Flux* Grid.AfaceW(i+1, j, k)/
                    Grid.volume(i+1, j ,k);
                
                       
               /* Include source terms associated with the finite-rate chemistry and
                  turbulence/chemistry interactions */
               
                 if (W[i][j][k].React.reactset_flag != NO_REACTIONS) {
                    
                    dUdt[i][j][k][k_residual] += IPs.CFL_Number*dt[i][j][k]*W[i][j][k].Sw(
                       W[i][j][k].React.reactset_flag, IPs.i_Flow_Type);
                    
                 } /* endif */
                        
                /* Save west and east face boundary flux. */
                
                //    if (i ==  ICl-1) {
//                     FluxW[j] = -Flux* Grid.lfaceW(i+1, j);
//                 } else if (i ==  ICu) {
//                     FluxE[j] = Flux* Grid.lfaceE(i, j);
//                 } /* endif */ 
                
             } /* endif */
          } /* endfor */
          
          if (( j >  JCl-1 && j <  JCu+1 ) && 
              ( k >  KCl-1 && k <  KCu+1 ) ) {
              dUdt[ ICl-1][j][k][k_residual] = U_VACUUM;
              dUdt[ ICu+1][j][k][k_residual] = U_VACUUM;
          } /* endif */
      } /* endfor */
   } /* endfor */
   
   // Add j-direction (eta-direction) fluxes.
   for ( k =  KCl-1 ; k <=  KCu+1 ; ++k ) {
      for ( i =  ICl ; i <=  ICu ; ++i ) {
          for ( j  =  JCl-1 ; j <=  JCu ; ++j ) {
             
             /* Evaluate the cell interface j-direction fluxes. */
             if (j ==  JCl-1 && 
                 ( Grid.BCtypeS[i][k] == BC_REFLECTION ||
                   Grid.BCtypeS[i][k] == BC_NO_SLIP||
                   Grid.BCtypeS[i][k] == BC_MOVING_WALL)) {
                
                dX =  Grid.xfaceS(i, j+1, k)- Grid.Cell[i][j+1][k].Xc;
                Wr =  W[i][j+1][k] +
                   ( phi[i][j+1][k]^ dWdx[i][j+1][k])*dX.x +
                   ( phi[i][j+1][k]^ dWdy[i][j+1][k])*dX.y+
                   ( phi[i][j+1][k]^ dWdz[i][j+1][k])*dX.z;
                if ( Grid.BCtypeS[i][k] == BC_REFLECTION) {
                   Wl =  NavierStokes3D_ThermallyPerfect_pState::Reflect(Wr,  Grid.nfaceS(i, j+1, k));
                }
                if ( Grid.BCtypeS[i][k] == BC_NO_SLIP) {
                   Wl =  NavierStokes3D_ThermallyPerfect_pState::No_Slip(Wr, WoS[i][k], Grid.nfaceS(i, j+1, k),
                                                                         IPs.Pressure_Gradient,
                                                                         FIXED_TEMPERATURE_WALL);
                }
                if ( Grid.BCtypeS[i][k] == BC_MOVING_WALL) {
                   Wl =  NavierStokes3D_ThermallyPerfect_pState::Moving_Wall(Wr, WoS[i][k], Grid.nfaceS(i, j+1, k),
                                                                             IPs.Pressure_Gradient,
                                                                             IPs.Moving_Wall_Velocity,
                                                                             FIXED_TEMPERATURE_WALL);
                }

             } else if (j ==  JCu && 
                        ( Grid.BCtypeN[i][k] == BC_REFLECTION ||
                          Grid.BCtypeN[i][k] == BC_NO_SLIP||
                          Grid.BCtypeN[i][k] == BC_MOVING_WALL)) {
                
                dX =  Grid.xfaceN(i, j, k)- Grid.Cell[i][j][k].Xc;
                Wl =  W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x+
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y+
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                if ( Grid.BCtypeN[i][k] == BC_REFLECTION) {
                   Wr =  NavierStokes3D_ThermallyPerfect_pState::Reflect(Wl,  Grid.nfaceN(i, j, k));
                }
                if ( Grid.BCtypeN[i][k] == BC_NO_SLIP) {
                   Wr =  NavierStokes3D_ThermallyPerfect_pState::No_Slip(Wl, WoN[i][k], Grid.nfaceN(i, j, k),
                                                                         IPs.Pressure_Gradient,
                                                                         FIXED_TEMPERATURE_WALL );
                }
                if ( Grid.BCtypeN[i][k] == BC_MOVING_WALL) {
                   Wr =  NavierStokes3D_ThermallyPerfect_pState::Moving_Wall(Wl, WoN[i][k], Grid.nfaceN(i, j, k),
                                                                             IPs.Moving_Wall_Velocity,
                                                                             IPs.Pressure_Gradient,
                                                                             FIXED_TEMPERATURE_WALL );
                }
             } else {
                dX =  Grid.xfaceN(i, j, k)- Grid.Cell[i][j][k].Xc;
                Wl =  W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                dX =  Grid.xfaceS(i, j+1, k)- Grid.Cell[i][j+1][k].Xc;
                Wr =  W[i][j+1][k] +
                   ( phi[i][j+1][k]^ dWdx[i][j+1][k])*dX.x +
                   ( phi[i][j+1][k]^ dWdy[i][j+1][k])*dX.y +
                   ( phi[i][j+1][k]^ dWdz[i][j+1][k])*dX.z;
          } /* endif */

             switch(IPs.i_Flux_Function) {
             case FLUX_FUNCTION_HLLE :
                Flux =  NavierStokes3D_ThermallyPerfect_pState::FluxHLLE_n(Wl, Wr,  Grid.nfaceN(i, j, k));
                break;
             case FLUX_FUNCTION_ROE :
                            
                Flux =  NavierStokes3D_ThermallyPerfect_pState::FluxRoe_n(Wl, Wr,  Grid.nfaceN(i, j, k));

                break;
                
             } /* endswitch */
    
     
             // add viscous flux in j direction
             Flux -=  NavierStokes3D_ThermallyPerfect_pState::FluxViscous_n(Wl, Wr, W[i][j][k], W[i][j+1][k], 
                                   dWdx[i][j][k], dWdy[i][j][k],dWdz[i][j][k],
                                   dWdx[i][j+1][k], dWdy[i][j+1][k],dWdz[i][j+1][k],
                                   Grid.nfaceN(i, j, k),  Grid.Von(i, j, k),
                                   Grid.delta_on(i, j, k), Grid.volume(i, j, k), 
                                   Grid.volume(i, j+1, k));
             
             /* Evaluate cell-averaged solution changes. */
             dUdt[i][j][k][k_residual] -= 
                (IPs.CFL_Number* dt[i][j][k])*
                Flux* Grid.AfaceN(i, j, k)/
                Grid.volume(i, j, k);
              dUdt[i][j+1][k][k_residual] += 
                (IPs.CFL_Number* dt[i][j+1][k])*
                Flux* Grid.AfaceS(i, j+1, k)/
                 Grid.volume(i, j+1, k);

          /* Save south and north face boundary flux. */

          // if (j ==  JCl-1) {
//               FluxS[i] = -Flux* Grid.lfaceS(i, j+1);
//           } else if (j ==  JCu) {
//               FluxN[i] = Flux* Grid.lfaceN(i, j);
//           } /* endif */

          } /* endfor */
           dUdt[i][ JCl-1][k][k_residual] = U_VACUUM;
           dUdt[i][ JCu+1][k][k_residual] = U_VACUUM;
      } /* endfor */
   } /* endfor */
   
   // Add k-direction (gamma-direction) fluxes.
   for ( i =  ICl-1; i <=  ICu+1 ; ++i ){
      for ( j  =  JCl ; j <=  JCu ; ++j ){
          for ( k =  KCl-1 ; k <=  KCu ; ++k )  {
             
             /* Evaluate the cell interface j-direction fluxes. */
             if (k ==  KCl-1 && 
                 ( Grid.BCtypeB[i][j] == BC_REFLECTION  || 
                   Grid.BCtypeB[i][j] == BC_NO_SLIP ||
                   Grid.BCtypeB[i][j] == BC_MOVING_WALL)) {
                
                dX =  Grid.xfaceBot(i, j, k+1)- Grid.Cell[i][j][k+1].Xc;
                Wr =  W[i][j][k+1] +
                   ( phi[i][j][k+1]^ dWdx[i][j][k+1])*dX.x +
                   ( phi[i][j][k+1]^ dWdy[i][j][k+1])*dX.y+
                   ( phi[i][j][k+1]^ dWdz[i][j][k+1])*dX.z;
                if ( Grid.BCtypeB[i][j] == BC_REFLECTION) {
                   Wl =  NavierStokes3D_ThermallyPerfect_pState::Reflect(Wr,  Grid.nfaceBot(i, j, k+1));
                }
                if ( Grid.BCtypeB[i][j] == BC_NO_SLIP) {
                   Wl =  NavierStokes3D_ThermallyPerfect_pState::No_Slip(Wr, WoB[i][j], Grid.nfaceBot(i, j, k+1), 
                                                                         IPs.Pressure_Gradient,
                                                                         FIXED_TEMPERATURE_WALL);
                }
                if ( Grid.BCtypeB[i][j] == BC_MOVING_WALL) {
                   Wl =  NavierStokes3D_ThermallyPerfect_pState::Moving_Wall(Wr, WoB[i][j], Grid.nfaceBot(i, j, k+1),
                                                                             IPs.Moving_Wall_Velocity,
                                                                             IPs.Pressure_Gradient,
                                                                             FIXED_TEMPERATURE_WALL);
                }


             } else if (k ==  KCu && 
                        ( Grid.BCtypeT[i][j] == BC_REFLECTION ||
                          Grid.BCtypeT[i][j] == BC_NO_SLIP ||
                          Grid.BCtypeT[i][j] == BC_MOVING_WALL)) {
                
                dX =  Grid.xfaceTop(i, j, k)- Grid.Cell[i][j][k].Xc;
                Wl =  W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x+
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y+
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                if ( Grid.BCtypeT[i][j] == BC_REFLECTION) {
                   Wr =  NavierStokes3D_ThermallyPerfect_pState::Reflect(Wl,  Grid.nfaceTop(i, j, k));
                }
                if ( Grid.BCtypeT[i][j] == BC_NO_SLIP) {
                   Wr =  NavierStokes3D_ThermallyPerfect_pState::No_Slip(Wl, WoT[i][j], Grid.nfaceTop(i, j, k), 
                                
                                                                         IPs.Pressure_Gradient,FIXED_TEMPERATURE_WALL );
                }
                if ( Grid.BCtypeT[i][j] == BC_MOVING_WALL) {
                   Wr =  NavierStokes3D_ThermallyPerfect_pState::Moving_Wall(Wl, WoT[i][j], Grid.nfaceTop(i, j, k),
                                                                             IPs.Moving_Wall_Velocity,  
                                                                             IPs.Pressure_Gradient,
                                                                             FIXED_TEMPERATURE_WALL );
                }
                
             } else {
                dX =  Grid.xfaceTop(i, j, k)- Grid.Cell[i][j][k].Xc;
                Wl =  W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                dX =  Grid.xfaceBot(i, j, k+1)- Grid.Cell[i][j][k+1].Xc;
                Wr =  W[i][j][k+1] +
                   ( phi[i][j][k+1]^ dWdx[i][j][k+1])*dX.x +
                   ( phi[i][j][k+1]^ dWdy[i][j][k+1])*dX.y +
                   ( phi[i][j][k+1]^ dWdz[i][j][k+1])*dX.z;
                          
                
             } /* endif */
             
             switch(IPs.i_Flux_Function) {
             case FLUX_FUNCTION_HLLE :
                Flux =  NavierStokes3D_ThermallyPerfect_pState::FluxHLLE_n(Wl, Wr,  Grid.nfaceTop(i, j, k));
                break;
                
             case FLUX_FUNCTION_ROE :
                Flux =  NavierStokes3D_ThermallyPerfect_pState::FluxRoe_n(Wl, Wr,  Grid.nfaceTop(i, j, k));
                             
                break;
                
             } /* endswitch */
    
             // add viscous flux in k direction
             Flux -=  NavierStokes3D_ThermallyPerfect_pState::FluxViscous_n(Wl, Wr,W[i][j][k], W[i][j][k+1], 
                                   dWdx[i][j][k], dWdy[i][j][k],dWdz[i][j][k],
                                   dWdx[i][j][k+1], dWdy[i][j][k+1],dWdz[i][j][k+1],
                                   Grid.nfaceTop(i, j, k), Grid.Vot(i, j, k), 
                                   Grid.delta_ot(i, j, k), Grid.volume(i, j, k), 
                                   Grid.volume(i, j, k+1));
       
             /* Evaluate cell-averaged solution changes. */
             
             dUdt[i][j][k][k_residual] -= 
                (IPs.CFL_Number* dt[i][j][k])*
                Flux* Grid.AfaceTop(i, j, k)/
                Grid.volume(i, j, k);
             dUdt[i][j][k+1][k_residual] += 
                (IPs.CFL_Number* dt[i][j][k+1])*
                Flux* Grid.AfaceBot(i, j, k+1)/
                Grid.volume(i, j, k+1);
             
             /* Save top and bottom face boundary flux. */

             // if (j ==  JCl-1) {
//               FluxS[i] = -Flux* Grid.lfaceS(i, j+1);
//           } else if (j ==  JCu) {
//               FluxN[i] = Flux* Grid.lfaceN(i, j);
//           } /* endif */
             
          } /* endfor */
          dUdt[i][j][ KCl-1][k_residual] = U_VACUUM;
          dUdt[i][j][ KCu+1][k_residual] = U_VACUUM;
       } /* endfor */
   } /* endfor */
         
   /* Residual for the stage successfully calculated. */

    
   return (0);
    
}
