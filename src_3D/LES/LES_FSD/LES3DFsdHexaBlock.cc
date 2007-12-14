#ifndef  _LES3DFSD_HEXA_BLOCK_INCLUDED
#include "LES3DFsdHexaBlock.h"
#endif //  _LES3DFSD_HEXA_BLOCK_INCLUDED

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution values at the nodes of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template< >
void Hexa_Block<LES3DFsd_pState,LES3DFsd_cState>::
Output_Tecplot(Input_Parameters<LES3DFsd_pState,LES3DFsd_cState> &IPs,
               const int Number_of_Time_Steps,
               const double &Time,
               const int Block_Number,
               const int Output_Title,
               ostream &Out_File) {
   int i, j, k;
 
   LES3DFsd_pState W_node;
   
//     for (int k  = KCl; k <= KCu; ++k ) {
//      for (int j  = JCl; j <= JCu; ++j ) {
//       for (int i = ICl; i <= ICu; ++i ) {
//          Reconstruction_Second_Derivatives();	 
//       }
//      }
//     }

   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
   
   BCs(IPs);
   
   /* Output node solution data. */
   
   Out_File << setprecision(14);
   if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 3D Solution, "
               << "Time Step/Iteration Level = " << Number_of_Time_Steps
               << ", Time = " << Time
               << "\"" << "\n"
               << "VARIABLES = \"x\" \\ \n"
               << "\"y\" \\ \n"
               << "\"z\" \\ \n"
               << "\"rho\" \\ \n"
               << "\"u\" \\ \n"
               << "\"v\" \\ \n"
               << "\"w\" \\ \n"
               << "\"p\" \\ \n"
               << "\"C\" \\ \n"
               << "\"FSD/rho\" \\ \n"
               << "\"k\" \\ \n";
      //n species mass fractions names
      for(int i =0 ;i<W[0][0][0].ns ;i++){
         Out_File <<"\"c_"<<W[0][0][0].specdata[i].Speciesname()<<"\" \\ \n";
      }
      
      Out_File <<"\"T\" \\ \n"
               <<"\"FSD\" \\ \n"
               <<"\"Vorticity\" \\ \n"
               <<"\"Mx\" \\ \n"
               <<"\"My\" \\ \n"  
               <<"\"Mz\" \\ \n"  
               <<"\"Reaction_Rate_Fsd\" \\ \n"  
               <<"\"Resolved_Strain\" \\ \n"  
               <<"\"Resolved_Propagation_Curvature\" \\ \n"  
               <<"\"SFS_Strain\" \\ \n"  
               <<"\"SFS_Curvature\" \\ \n";
//                <<"\"Resolved_Curvature\" \\ \n"  
//                <<"\"Resolved_Propagation\" \\ \n"  
//                <<"\"Resolved_Convection_Progvar\" \\ \n"  
//   	          <<"\"Resolved_Convection_Fsd\" \\ \n"  
//                <<"\"NGT_Progvar\" \\ \n"  
//                <<"\"NGT_Fsd\" \\ \n"  
//                <<"\"SFS_Diffusion_Progvar\" \\ \n"  
//                <<"\"SFS_Diffusion_Fsd\" \\ \n"  
//                <<"\"Heat_Release_Strain\" \\ \n"  
//                <<"\"Net_Rate_Change_Progvar\" \\ \n"  
//                <<"\"Net_Rate_Change_Fsd\" \\ \n";

      
      Out_File<< "ZONE T =  \"Block Number = " << Block_Number
              << "\" \\ \n"
              << "I = " << Grid.INu -  Grid.INl + 1 << " \\ \n"
              << "J = " << Grid.JNu -  Grid.JNl + 1 << " \\ \n"
              << "K = " << Grid.KNu -  Grid.KNl + 1 << " \\ \n"
              << "DATAPACKING = POINT \n";
      
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
               << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
               << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";              
   } /* endif */
   
   
   for ( k  =  Grid.KNl ; k <=  Grid.KNu ; ++k ) {
      for ( j  =  Grid.JNl ; j <=  Grid.JNu ; ++j ) {
         for ( i =  Grid.INl ; i <=  Grid.INu ; ++i ) {
            W_node = Wn(i, j, k);
            Out_File << " "  << Grid.Node[i][j][k].X << W_node;
            Out_File.setf(ios::scientific);
            Out_File << " " << W_node.T() 
                     << " " << W_node.Fsd*W_node.rho
                     << " " << W_node.vorticity(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k])
                     << " " << W_node.M_x(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k])
                     << " " << W_node.M_y(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
                     << " " << W_node.M_z(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
                     << " " << W_node.Reaction_Rate_Fsd(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
                     << " " << W_node.Resolved_Strain(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
                     << " " << W_node.Resolved_Propagation_Curvature(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
                     << " " << W_node.SFS_Strain(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],Flow_Type,Grid.volume(i,j,k)) 
                     << " " << W_node.SFS_Curvature(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) << "\n";
//                      << " " << W_node.Resolved_Curvature(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                          d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                          d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k]) 
//                      << " " << W_node.Resolved_Propagation(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                          d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                          d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k])
//                      << " " << W_node.Resolved_Convection_Progvar(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
//                      << " " << W_node.Resolved_Convection_Fsd(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
//                      << " " << W_node.NGT_Progvar(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
//                      << " " << W_node.NGT_Fsd(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                               d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                               d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k]) 
//                      << " " << W_node.SFS_Diffusion_Progvar(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                             d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                             d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k],
//                                                             Flow_Type,Grid.volume(i,j,k)) 
//                      << " " << W_node.SFS_Diffusion_Fsd(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                         d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                         d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k],
//                                                         Flow_Type,Grid.volume(i,j,k)) 
//                      << " " << W_node.Heat_Release_Strain(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                           d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                           d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k]) 
//                      << " " << W_node.Net_Rate_Change_Progvar(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                               d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                               d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k],
//                                                               Flow_Type,Grid.volume(i,j,k)) 
//                      << " " << W_node.Net_Rate_Change_Fsd(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                           d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                           d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k],
//                                                           Flow_Type,Grid.volume(i,j,k)) << "\n; 

            Out_File.unsetf(ios::scientific);
         } /* endfor */
      } /* endfor */
   } /* endfor */

   Out_File << setprecision(6);
  
}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template< >
void Hexa_Block<LES3DFsd_pState,LES3DFsd_cState>::
Output_Cells_Tecplot(Input_Parameters<LES3DFsd_pState, 
                                      LES3DFsd_cState> &IPs,
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File) {

   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
   
   BCs(IPs);

   /* Output cell centred solution data. */

   int i, j, k;
   Out_File << setprecision(14);

//     for (int k  = KCl; k <= KCu; ++k ) {
//      for (int j  = JCl; j <= JCu; ++j ) {
//       for (int i = ICl; i <= ICu; ++i ) {
//          Reconstruction_Second_Derivatives();	 
//       }
//      }
//     }

   if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 3D Solution, "
               << "Time Step/Iteration Level = " << Number_of_Time_Steps
               << ", Time = " << Time
               << "\"" << "\n"
               << "VARIABLES = \"x\" \\ \n"
               << "\"y\" \\ \n"
               << "\"z\" \\ \n"
               << "\"rho\" \\ \n"
               << "\"u\" \\ \n"
               << "\"v\" \\ \n"
               << "\"w\" \\ \n"
               << "\"p\" \\ \n"
               << "\"C\" \\ \n"
               << "\"FSD/rho\" \\ \n"
               << "\"k\" \\ \n";
      //n species mass fractions names
      for( i =0; i<W[0][0][0].ns; i++){
         Out_File <<"\"c"<<W[0][0][0].specdata[i].Speciesname()<<"\" \\ \n";
      }
     
      Out_File <<"\"T\" \\ \n"
               <<"\"R\" \\ \n";
      
      Out_File  << "\"Mu_t\" \\ \n" 
                << "\"FSD\" \\ \n"
                << "\"Reaction_Rate_Fsd\" \\ \n"
                << "\"M_x\" \\ \n"
                << "\"M_y\" \\ \n"
                << "\"M_z\" \\ \n"
                << "\"Resolved_Strain\" \\ \n"
                << "\"Resolved_Propagation_Curvature\" \\ \n"
                << "\"SFS_Strain\" \\ \n"
                << "\"SFS_Curvature\" \\ \n";
//                 << "\"Resolved_Curvature\" \\ \n"
//                 << "\"Resolved_Propagation\" \\ \n"
//                 << "\"Resolved_Convection_Progvar\" \\ \n"
//                 << "\"Resolved_Convection_Fsd\" \\ \n"
//                 << "\"NGT_Progvar\" \\ \n"
//                 << "\"NGT_Fsd\" \\ \n"
//                 << "\"SFS_Diffusion_Progvar\" \\ \n"
//                 << "\"SFS_Diffusion_Fsd\" \\ \n"
//                 << "\"Heat_Release_Strain\" \\ \n"
//                 << "\"Net_Rate_Change_Progvar\" \\ \n"
//                 << "\"Net_Rate_Change_Fsd\" \\ \n";

      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " <<  ICu -  ICl + 2*Nghost +1 << " \\ \n"
               << "J = " <<  JCu -  JCl + 2*Nghost +1 << " \\ \n"
               << "K = " <<  KCu -  KCl + 2*Nghost +1 << " \\ \n"
               << "DATAPACKING = POINT \n";
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " <<  Grid.ICu -  Grid.ICl + 2* Nghost + 1 << " \\ \n"
               << "J = " <<  Grid.JCu -  Grid.JCl + 2* Nghost + 1 << " \\ \n"
               << "K = " <<  Grid.KCu -  Grid.KCl + 2* Nghost + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";
      
   } /* endif */

   for ( k =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k) {
      for ( j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) {
         for ( i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
            
            Out_File << " "  <<  Grid.Cell[i][j][k].Xc
                     <<  W[i][j][k];
            Out_File.setf(ios::scientific);
            Out_File << " " <<W[i][j][k].T() 
                     << " " <<W[i][j][k].Rtot()<<"\n ";

            Out_File << " " <<W[i][j][k].mu_t(dWdx[i][j][k],
                                              dWdy[i][j][k],
                                              dWdz[i][j][k],
                                              Flow_Type,
                                              Grid.volume(i,j,k)) 
                     << " " <<W[i][j][k].Fsd*W[i][j][k].rho
                     << " " <<W[i][j][k].Reaction_Rate_Fsd(dWdx[i][j][k],
                                                           dWdy[i][j][k],
                                                           dWdz[i][j][k])
                     << " " <<W[i][j][k].M_x(dWdx[i][j][k],
                                             dWdy[i][j][k],
                                             dWdz[i][j][k])
                     << " " <<W[i][j][k].M_y(dWdx[i][j][k],
                                             dWdy[i][j][k],
                                             dWdz[i][j][k])
                     << " " <<W[i][j][k].M_z(dWdx[i][j][k],
                                             dWdy[i][j][k],
                                             dWdz[i][j][k])
                     << " " <<W[i][j][k].Resolved_Strain(dWdx[i][j][k],
                                                         dWdy[i][j][k],
                                                         dWdz[i][j][k])
                     << " " <<W[i][j][k].Resolved_Propagation_Curvature(dWdx[i][j][k],
                                                                        dWdy[i][j][k],
                                                                        dWdz[i][j][k])
                     << " " <<W[i][j][k].SFS_Strain(dWdx[i][j][k],
                                                    dWdy[i][j][k],
                                                    dWdz[i][j][k],
                                                    Flow_Type,
                                                    Grid.volume(i,j,k))
                     << " " <<W[i][j][k].SFS_Curvature(dWdx[i][j][k],
                                                       dWdy[i][j][k],
                                                       dWdz[i][j][k]) << "\n";
//                      << " " <<W[i][j][k].Resolved_Curvature(dWdx[i][j][k],
//                                                             dWdy[i][j][k],
//                                                             dWdz[i][j][k],
//                                                             d_dWdx_dx[i][j][k],
//                                                             d_dWdy_dy[i][j][k],
//                                                             d_dWdz_dz[i][j][k],
//                                                             d_dWdx_dy[i][j][k],
//                                                             d_dWdx_dz[i][j][k],
//                                                             d_dWdy_dz[i][j][k])
//                      << " " <<W[i][j][k].Resolved_Propagation(dWdx[i][j][k],
//                                                               dWdy[i][j][k],
//                                                               dWdz[i][j][k],
//                                                               d_dWdx_dx[i][j][k],
//                                                               d_dWdy_dy[i][j][k],
//                                                               d_dWdz_dz[i][j][k],
//                                                               d_dWdx_dy[i][j][k],
//                                                               d_dWdx_dz[i][j][k],
//                                                               d_dWdy_dz[i][j][k])
//                      << " " <<W[i][j][k].Resolved_Convection_Progvar(dWdx[i][j][k],
//                                                                      dWdy[i][j][k],
//                                                                      dWdz[i][j][k])
//                      << " " <<W[i][j][k].Resolved_Convection_Fsd(dWdx[i][j][k],
//                                                                  dWdy[i][j][k],
//                                                                  dWdz[i][j][k])
//                      << " " <<W[i][j][k].NGT_Progvar(dWdx[i][j][k],
//                                                      dWdy[i][j][k],
//                                                      dWdz[i][j][k])
//                      << " " <<W[i][j][k].NGT_Fsd(dWdx[i][j][k],
//                                                  dWdy[i][j][k],
//                                                  dWdz[i][j][k],
//                                                  d_dWdx_dx[i][j][k],
//                                                  d_dWdy_dy[i][j][k],
//                                                  d_dWdz_dz[i][j][k],
//                                                  d_dWdx_dy[i][j][k],
//                                                  d_dWdx_dz[i][j][k],
//                                                  d_dWdy_dz[i][j][k])
//                      << " " <<W[i][j][k].SFS_Diffusion_Progvar(dWdx[i][j][k],
//                                                                dWdy[i][j][k],
//                                                                dWdz[i][j][k],
//                                                                d_dWdx_dx[i][j][k],
//                                                                d_dWdy_dy[i][j][k],
//                                                                d_dWdz_dz[i][j][k],
//                                                                d_dWdx_dy[i][j][k],
//                                                                d_dWdx_dz[i][j][k],
//                                                                d_dWdy_dz[i][j][k],
//                                                                Flow_Type,
//                                                                Grid.volume(i,j,k))
//                      << " " <<W[i][j][k].SFS_Diffusion_Fsd(dWdx[i][j][k],
//                                                            dWdy[i][j][k],
//                                                            dWdz[i][j][k],
//                                                            d_dWdx_dx[i][j][k],
//                                                            d_dWdy_dy[i][j][k],
//                                                            d_dWdz_dz[i][j][k],
//                                                            d_dWdx_dy[i][j][k],
//                                                            d_dWdx_dz[i][j][k],
//                                                            d_dWdy_dz[i][j][k],
//                                                            Flow_Type,
//                                                            Grid.volume(i,j,k))
//                      << " " <<W[i][j][k].Heat_Release_Strain(dWdx[i][j][k],
//                                                              dWdy[i][j][k],
//                                                              dWdz[i][j][k],
//                                                              d_dWdx_dx[i][j][k],
//                                                              d_dWdy_dy[i][j][k],
//                                                              d_dWdz_dz[i][j][k],
//                                                              d_dWdx_dy[i][j][k],
//                                                              d_dWdx_dz[i][j][k],
//                                                              d_dWdy_dz[i][j][k])
//                      << " " <<W[i][j][k].Net_Rate_Change_Progvar(dWdx[i][j][k],
//                                                                  dWdy[i][j][k],
//                                                                  dWdz[i][j][k],
//                                                                  d_dWdx_dx[i][j][k],
//                                                                  d_dWdy_dy[i][j][k],
//                                                                  d_dWdz_dz[i][j][k],
//                                                                  d_dWdx_dy[i][j][k],
//                                                                  d_dWdx_dz[i][j][k],
//                                                                  d_dWdy_dz[i][j][k],
//                                                                  Flow_Type,
//                                                                  Grid.volume(i,j,k))
//                      << " " <<W[i][j][k].Net_Rate_Change_Fsd(dWdx[i][j][k],
//                                                              dWdy[i][j][k],
//                                                              dWdz[i][j][k],
//                                                              d_dWdx_dx[i][j][k],
//                                                              d_dWdy_dy[i][j][k],
//                                                              d_dWdz_dz[i][j][k],
//                                                              d_dWdx_dy[i][j][k],
//                                                              d_dWdx_dz[i][j][k],
//                                                              d_dWdy_dz[i][j][k],
//                                                              Flow_Type,
//                                                              Grid.volume(i,j,k)) << "\n";
            Out_File.unsetf(ios::scientific);
         } /* endfor */
      } /* endfor */
   } /* endfor */
   
   Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution values at the nodes of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template< >
void Hexa_Block<LES3DFsd_pState,LES3DFsd_cState>::
Output_Nodes_Tecplot(Input_Parameters<LES3DFsd_pState,LES3DFsd_cState> &IPs,
               const int Number_of_Time_Steps,
               const double &Time,
               const int Block_Number,
               const int Output_Title,
               ostream &Out_File) {
   int i, j, k;
 
   LES3DFsd_pState W_node;
   
//     for (int k  = KCl; k <= KCu; ++k ) {
//      for (int j  = JCl; j <= JCu; ++j ) {
//       for (int i = ICl; i <= ICu; ++i ) {
//          Reconstruction_Second_Derivatives();	 
//       }
//      }
//     }

   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
   
   BCs(IPs);
   
   /* Output node solution data. */
   
   Out_File << setprecision(14);
   if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 3D Solution, "
               << "Time Step/Iteration Level = " << Number_of_Time_Steps
               << ", Time = " << Time
               << "\"" << "\n"
               << "VARIABLES = \"x\" \\ \n"
               << "\"y\" \\ \n"
               << "\"z\" \\ \n"
               << "\"rho\" \\ \n"
               << "\"u\" \\ \n"
               << "\"v\" \\ \n"
               << "\"w\" \\ \n"
               << "\"p\" \\ \n"
               << "\"C\" \\ \n"
               << "\"FSD/rho\" \\ \n"
               << "\"k\" \\ \n";
      //n species mass fractions names
      for(int i =0 ;i<W[0][0][0].ns ;i++){
         Out_File <<"\"c_"<<W[0][0][0].specdata[i].Speciesname()<<"\" \\ \n";
      }
      
      Out_File <<"\"T\" \\ \n"
               <<"\"FSD\" \\ \n"
               <<"\"Mx\" \\ \n"
               <<"\"My\" \\ \n"  
               <<"\"Mz\" \\ \n"  
               <<"\"Reaction_Rate_Fsd\" \\ \n"  
               <<"\"Resolved_Strain\" \\ \n"  
               <<"\"Resolved_Propagation_Curvature\" \\ \n"  
               <<"\"SFS_Strain\" \\ \n"  
               <<"\"SFS_Curvature\" \\ \n";
//                <<"\"Resolved_Curvature\" \\ \n"  
//                <<"\"Resolved_Propagation\" \\ \n"  
//                <<"\"Resolved_Convection_Progvar\" \\ \n"  
//                <<"\"Resolved_Convection_Fsd\" \\ \n"  
//                <<"\"NGT_Progvar\" \\ \n"  
//                <<"\"NGT_Fsd\" \\ \n"  
//                <<"\"SFS_Diffusion_Progvar\" \\ \n"  
//                <<"\"SFS_Diffusion_Fsd\" \\ \n"  
//                <<"\"Heat_Release_Strain\" \\ \n"  
//                <<"\"Net_Rate_Change_Progvar\" \\ \n"  
//                <<"\"Net_Rate_Change_Fsd\" \\ \n";

      
      Out_File<< "ZONE T =  \"Block Number = " << Block_Number
              << "\" \\ \n"
              << "I = " << Grid.INu -  Grid.INl + 1 << " \\ \n"
              << "J = " << Grid.JNu -  Grid.JNl + 1 << " \\ \n"
              << "K = " << Grid.KNu -  Grid.KNl + 1 << " \\ \n"
              << "DATAPACKING = POINT \n";
      
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
               << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
               << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";              
   } /* endif */
   
   
   for ( k  =  Grid.KNl ; k <=  Grid.KNu ; ++k ) {
      for ( j  =  Grid.JNl ; j <=  Grid.JNu ; ++j ) {
         for ( i =  Grid.INl ; i <=  Grid.INu ; ++i ) {
            W_node = Wn(i, j, k);
            Out_File << " "  << Grid.Node[i][j][k].X <<W_node;
            Out_File.setf(ios::scientific);
            Out_File << " " << W_node.T()
                     << " " << W_node.Fsd*W_node.rho
                     << " " << W_node.M_x(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k])
                     << " " << W_node.M_y(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
                     << " " << W_node.M_z(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
                     << " " << W_node.Reaction_Rate_Fsd(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
                     << " " << W_node.Resolved_Strain(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
                     << " " << W_node.Resolved_Propagation_Curvature(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
                     << " " << W_node.SFS_Strain(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],Flow_Type,Grid.volume(i,j,k)) 
                     << " " << W_node.SFS_Curvature(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) << "\n";
//                      << " " << W_node.Resolved_Curvature(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                          d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                          d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k]) 
//                      << " " << W_node.Resolved_Propagation(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                          d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                          d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k])
//                      << " " << W_node.Resolved_Convection_Progvar(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
//                      << " " << W_node.Resolved_Convection_Fsd(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
//                      << " " << W_node.NGT_Progvar(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) 
//                      << " " << W_node.NGT_Fsd(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                               d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                               d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k]) 
//                      << " " << W_node.SFS_Diffusion_Progvar(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                             d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                             d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k],
//                                                             Flow_Type,Grid.volume(i,j,k)) 
//                      << " " << W_node.SFS_Diffusion_Fsd(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                         d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                         d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k],
//                                                         Flow_Type,Grid.volume(i,j,k)) 
//                      << " " << W_node.Heat_Release_Strain(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                           d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                           d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k]) 
//                      << " " << W_node.Net_Rate_Change_Progvar(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                               d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                               d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k],
//                                                               Flow_Type,Grid.volume(i,j,k)) 
//                      << " " << W_node.Net_Rate_Change_Fsd(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
//                                                           d_dWdx_dx[i][j][k],d_dWdy_dy[i][j][k],d_dWdz_dz[i][j][k],
//                                                           d_dWdx_dy[i][j][k],d_dWdx_dz[i][j][k],d_dWdy_dz[i][j][k],
//                                                           Flow_Type,Grid.volume(i,j,k)) << "\n"; 

            Out_File.unsetf(ios::scientific);
         } /* endfor */
      } /* endfor */
   } /* endfor */

   Out_File << setprecision(6);
  
}

/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Apply initial conditions for the specified hexa      *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
template< >
int Hexa_Block<LES3DFsd_pState,LES3DFsd_cState>::
ICs(const int i_ICtype,
    Input_Parameters<LES3DFsd_pState,LES3DFsd_cState> &IPs){
   
   int i, j, k;
   
   double dpdx, dpdy, dpdz, delta_pres_x, delta_pres_y, delta_pres_z ;
   double zd, zz, di, Um;
   
   Flow_Type = IPs.i_Flow_Type;
   
   LES3DFsd_pState Wl, Wr;
   
   switch(i_ICtype) {
      
   case IC_UNIFORM :
      // Set the solution state to the initial state Wo[0].
      for (k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) 
         for ( j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) 
            for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
               W[i][j][k] = IPs.Wo;
               U[i][j][k] = W[i][j][k].U( );
               
            }
      break;
      
   case IC_TURBULENT_PREMIXED_FLAME :

    ifstream InFile;

      for (k  =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k ) 
         for (j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) 
            for ( i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
	      double xx = -Grid.Cell[i][j][k].Xc.x;//+(1.0/2.0)*IPs.Grid_IP.Box_Length;
	      double tau_fsd = 2218/298-1.0;//W[i][j][k].HeatRelease_Parameter();
       	      W[i][j][k].C = (erf(xx*4000.0)+1.0)/2.0;
              W[i][j][k].p = 101325.0;
	      W[i][j][k].premixed_mfrac();
       	      W[i][j][k].rho = 1.13*W[ICu][j][k].Rtot()/W[i][j][k].Rtot()/(1.0+tau_fsd*W[i][j][k].C);
              W[i][j][k].v.x = 1.13*0.3837/W[i][j][k].rho;
       	      W[i][j][k].Fsd = 3000.0*exp(-sqr(xx*4000.0))/sqrt(3.1415926)/W[i][j][k].rho;
	      W[i][j][k].k = 0.0;
              U[i][j][k] = W[i][j][k].U();

            } /* endfor */

      // Open data file for reading
      InFile.open("Initial_Turbulence_Fluctuations.dat", ios::in); 
      // Check to see if successful
      if(InFile.fail()){ 
        cerr<<"\nError opening file: Initial_Turbulence.dat to read" <<endl;
        exit(1); 
      } 

      bool interpolate_flag;
      double xx, yy, zz, dx, dy, dz, dd, uprime_x, uprime_y, uprime_z;

      InFile.setf(ios::skipws);
 
      for ( i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) 
         for (j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) 
	   for (k  =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k ) {

	 if ((i>=ICl && i<=ICu) && (j>=JCl && j<=JCu) && (k>=KCl && k<=KCu)) {

 	    interpolate_flag = true;

 	    do{
	      InFile >> xx >> yy >> zz >> uprime_x >> uprime_y >> uprime_z;

 	      if (fabs((xx - Grid.Cell[i][j][k].Xc.x)/Grid.Cell[i][j][k].Xc.x) < 1.0E-2  &&
		  fabs((yy - Grid.Cell[i][j][k].Xc.y)/Grid.Cell[i][j][k].Xc.y) < 1.0E-2  &&
		  fabs((zz - Grid.Cell[i][j][k].Xc.z)/Grid.Cell[i][j][k].Xc.z) < 1.0E-2) {

 		W[i][j][k].v.x += uprime_x;
 		W[i][j][k].v.y += uprime_y;
 		W[i][j][k].v.z += uprime_z;
 		interpolate_flag = false;
 		break;

  	      } else {
 		dx = Grid.Cell[i][j][k].Xc.x - xx;
 		dy = Grid.Cell[i][j][k].Xc.y - yy;
 		dz = Grid.Cell[i][j][k].Xc.z - zz;
 		dd = sqrt(dx*dx + dy*dy + dz*dz);

 		if ( dd == ZERO ||  dd > sqrt(sqr(IPs.Grid_IP.Box_Width) + sqr(IPs.Grid_IP.Box_Height) + sqr(IPs.Grid_IP.Box_Length) )) {
 		  cout << "\nYOU ARE NOT SUPPOSED TO REACH THIS LINE!!!";

  		} //else {
	     } // end if
	 }  while ( !InFile.eof() ); // end while
	  
	    // reset EOF flag
//  	    if ( InFile.eof() ) {
//  	      InFile.clear(); 
//  	      // Go to the beginning of the file 
//   	      InFile.seekg(0, ios::beg);
//  	    } 
	   
	 }
	    U[i][j][k] = W[i][j][k].U();
    }

      InFile.unsetf(ios::skipws);
      InFile.close();

   if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K){ 
     Linear_Reconstruction_LeastSquares(IPs.i_Limiter);
     for (k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) 
         for ( j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) 
            for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	      //       	      W[i][j][k].Fsd = 200.0*W[i][j][k].Fsd;
 		W[i][j][k].k = 0.005*sqr(W[i][j][k].filter_width(Grid.volume(i,j,k))*W[i][j][k].abs_strain_rate(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]));
  	        U[i][j][k] = W[i][j][k].U();
 	      }
   }

     break;  
      
   } //end of switch
   
   /* Set default values for the boundary conditions
      reference states. */

   for ( k = KCl-Nghost ; k<= KCu+Nghost; ++k )
      for ( j = JCl-Nghost ; j<= JCu+Nghost; ++j ){
         if ((k >= KCl && k <= KCu) && (j >= JCl && j <= JCu)) {
            WoW[j][k] = W[ICl][j][k];
            WoE[j][k] = W[ICu][j][k];
         } else if (j < JCl && k < KCl ) {
            WoW[j][k] = W[ICl][JCl][KCl];
            WoE[j][k] = W[ICu][JCl][KCl];
         } else if (j > JCu && k> KCu) {
            WoW[j][k] = W[ICl][JCu][KCu];
            WoE[j][k] = W[ICu][JCu][KCu];
         } else if(j < JCl &&(k >= KCl && k <= KCu)){
            WoW[j][k] = W[ICl][JCl][k];
            WoE[j][k] = W[ICu][JCl][k];
         } else if(j > JCu && (k >= KCl && k <= KCu)){
            WoW[j][k] = W[ICl][JCu][k];
            WoE[j][k] = W[ICu][JCu][k];
         } else if(k < KCl &&(j >= JCl && j <= JCu)){
            WoW[j][k] = W[ICl][j][KCl];
            WoE[j][k] = W[ICu][j][KCl];
         } else if(k > KCu && (j >= JCl && j <= JCu)){
            WoW[j][k] = W[ICl][j][KCu];
            WoE[j][k] = W[ICu][j][KCu];
         } else if(k > KCu && j < JCl ){
            WoW[j][k] = W[ICl][JCl][KCu];
            WoE[j][k] = W[ICu][JCl][KCu];
         } else if(k < KCl && j > JCu){
            WoW[j][k] = W[ICl][JCu][KCl];
            WoE[j][k] = W[ICu][JCu][KCl];
         }
      } /* endfor */ 
    
   for (  k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) 
      for (  i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
         if ((k >= KCl && k <= KCu) && (i >= ICl && i <= ICu)) {
            WoS[i][k] = W[i][JCl][k];
            WoN[i][k] = W[i][JCu][k];
         } else if (i < ICl && k< KCl) {
            WoS[i][k] = W[ICl][JCl][KCl];
            WoN[i][k] = W[ICl][JCu][KCl];
         } else if (i > ICu && k > KCu) {
            WoS[i][k] = W[ICu][JCl][KCu];
            WoN[i][k] = W[ICu][JCu][KCu];
         } else if (i<ICl && (k >= KCl && k <= KCu)){
            WoS[i][k] = W[ICl][JCl][k];
            WoN[i][k] = W[ICl][JCu][k];
         } else if (i>ICu && (k >= KCl && k <= KCu)){
            WoS[i][k] = W[ICu][JCl][k];
            WoN[i][k] = W[ICu][JCu][k];
         } else if ((i >= ICl && i <= ICu) && k< KCl) {
            WoS[i][k] = W[i][JCl][KCl];
            WoN[i][k] = W[i][JCu][KCl];
         } else if ((i >= ICl && i <= ICu) && k > KCu) {
            WoS[i][k] = W[i][JCl][KCu];
            WoN[i][k] = W[i][JCu][KCu];
         } else if (i < ICl  && k > KCu) {
            WoS[i][k] = W[ICl][JCl][KCu];
            WoN[i][k] = W[ICl][JCu][KCu];
         } else if (i >ICu  && k < KCl) {
            WoS[i][k] = W[ICu][JCl][KCl];
            WoN[i][k] = W[ICu][JCu][KCl];
         }
      } /* endfor */

   for (  j = JCl-Nghost ; j <= JCu+Nghost ; ++j ) 
      for (  i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
          if ((j >= JCl && j <= JCu) && (i >= ICl && i <= ICu)) {
             WoT[i][j] = W[i][j][KCu];
             WoB[i][j] = W[i][j][KCl];
          }  else if (i < ICl &&  j< JCl) {
             WoT[i][j] = W[ICl][JCl][KCu];
             WoB[i][j] = W[ICl][JCl][KCl];
          } else if(i > ICu &&  j > JCu) {
             WoT[i][j] = W[ICu][JCu][KCu];
             WoB[i][j] = W[ICu][JCu][KCl];
          }else if (i < ICl && (j >= JCl && j <= JCu)) {
             WoT[i][j] = W[ICl][j][KCu];
             WoB[i][j] = W[ICl][j][KCl];
          }else if (i > ICu && (j >= JCl && j <= JCu)) {
             WoT[i][j] = W[ICu][j][KCu];
             WoB[i][j] = W[ICu][j][KCl];
          } else if ((i >= ICl && i <= ICu) &&  j< JCl) {
             WoT[i][j] = W[i][JCl][KCu];
             WoB[i][j] = W[i][JCl][KCl];
          } else if ((i >= ICl && i <= ICu) &&  j> JCu) {
             WoT[i][j] = W[i][JCu][KCu];
             WoB[i][j] = W[i][JCu][KCl];
          } else if (i > ICu && j < JCl) {
             WoT[i][j] = W[ICu][JCl][KCu];
             WoB[i][j] = W[ICu][JCl][KCl];
          } else if (i < ICl && j > JCu) {
             WoT[i][j] = W[ICl][JCu][KCu];
             WoB[i][j] = W[ICl][JCu][KCl];
          }
       } /* endfor */
      
    return (0);
    
}

template< >
void Hexa_Block<LES3DFsd_pState,LES3DFsd_cState>::
BCs(Input_Parameters<LES3DFsd_pState,LES3DFsd_cState> &IPs){
   
   int i, j, k;
   double dpdx, dpdy, dpdz;
   Vector3D dX;
   Vector3D MOVING_WALL_VELOCITY = IPs.Moving_Wall_Velocity;

   for ( k =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k) 
      for ( j =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) {
         
         // Prescribe West boundary conditions.
         switch( Grid.BCtypeW[j][k]) {
            
         case BC_NONE :
            break;
            
         case BC_REFLECTION :
            W[ICl-1][j][k] = LES3DFsd_pState::Reflect(W[ICl][j][k],Grid.nfaceW(ICl,j,k));
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = LES3DFsd_pState::Reflect(W[ICl+1][j][k],Grid.nfaceW(ICl,j,k));
            U[ICl-2][j][k] = W[ICl-2][j][k].U();
            break;

         case BC_NO_SLIP :
            W[ICl-1][j][k] = LES3DFsd_pState::NoSlip(W[ICl][j][k],
                                                     WoW[j][k],
                                                     Grid.nfaceW(ICl,j,k),
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = LES3DFsd_pState::NoSlip(W[ICl+1][j][k],
                                                     WoW[j][k],
                                                     Grid.nfaceW(ICl,j,k),
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[ICl-2][j][k] = W[ICl-2][j][k].U();
            break;

         case BC_MOVING_WALL :
            W[ICl-1][j][k] = LES3DFsd_pState::MovingWall(W[ICl][j][k],
                                                         WoW[j][k],
                                                         Grid.nfaceW(ICl,j,k),
                                                         MOVING_WALL_VELOCITY,
                                                         IPs.Pressure_Gradient,
                                                         FIXED_TEMPERATURE_WALL);
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = LES3DFsd_pState::MovingWall(W[ICl+1][j][k], 
                                                         WoW[j][k], 
                                                         Grid.nfaceW(ICl,j,k),
                                                         MOVING_WALL_VELOCITY, 
                                                         IPs.Pressure_Gradient,
                                                         FIXED_TEMPERATURE_WALL);
            U[ICl-2][j][k] = W[ICl-2][j][k].U();
            break;

         case BC_FIXED_PRESSURE :
            W[ICl-1][j][k] = W[ICl][j][k];
            W[ICl-1][j][k].p = WoW[j][k].p;
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = W[ICl][j][k] ;
            W[ICl-2][j][k].p = WoW[j][k].p;
            U[ICl-2][j][k] = W[ICl-2][j][k].U();
            break;

         case BC_CONSTANT_EXTRAPOLATION :
            W[ICl-1][j][k] = W[ICl][j][k];
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = W[ICl][j][k] ;
            U[ICl-2][j][k] = W[ICl-2][j][k].U();
            break;

         case BC_CHANNEL_INFLOW:
            dpdx = IPs.Pressure_Gradient.x; 
            //for turbulent channel flow
            // p linearly varys based on constant pressure gradient 
            dX = Grid.Cell[ICl-1][j][k].Xc - Grid.Cell[ICl][j][k].Xc; 
            W[ICl-1][j][k] = WoW[j][k];
            W[ICl-1][j][k].v.x = W[ICl][j][k].v.x;
            W[ICl-1][j][k].p = WoW[j][k].p - dpdx*dX.x ;
            U[ICl-1][j][k] = W[ICl-1][j][k].U( );
            
            dX = Grid.Cell[ICl-2][j][k].Xc - Grid.Cell[ICl][j][k].Xc;
            W[ICl-2][j][k] = WoW[j][k];
            W[ICl-2][j][k].v.x = W[ICl][j][k].v.x;
            W[ICl-2][j][k].p = WoW[j][k].p - dpdx*dX.x;
            //W[ICl-2][j][k].p = (WoW[j][k].p - 2.0/3.0*W[ICl-2][j][k].rho*W[ICl-2][j][k].k) - dpdx*dX.x;
            U[ICl-2][j][k] = W[ICl-2][j][k].U( );
            break;
            
         case BC_PERIODIC :  
            W[ICl-1][j][k] = W[ICu-1][j][k];
            U[ICl-1][j][k] = U[ICu-1][j][k];
            W[ICl-2][j][k] = W[ICu-2][j][k];
            U[ICl-2][j][k] = U[ICu-2][j][k];
            break;

         case BC_INFLOW_SUBSONIC :
         // all fixed except pressure which is constant extrapolation
	    W[ICl-1][j][k] = WoW[j][k];
            W[ICl-1][j][k].p = W[ICl][j][k].p;
	    U[ICl-1][j][k] = W[ICl-1][j][k].U();
		
	    W[ICl-2][j][k] = WoW[j][k];
            W[ICl-2][j][k].p = W[ICl][j][k].p;
	    U[ICl-2][j][k] = W[ICl-2][j][k].U();
	    break;

         case BC_OUTFLOW_SUBSONIC :
	// all constant extrapolation except pressure which is fixed.
	    W[ICl-1][j][k] = W[ICl][j][k]; 
	    W[ICl-1][j][k].p = WoW[j][k].p;
	    U[ICl-1][j][k] = W[ICl-1][j][k].U();
	
	    W[ICl-2][j][k] = W[ICl][j][k]; 
	    W[ICl-2][j][k].p = WoW[j][k].p;
	    U[ICl-2][j][k] = W[ICl-2][j][k].U();
	    break;

         } /* endswitch */
         
         // Prescribe East boundary conditions.
         switch( Grid.BCtypeE[j][k]) {

         case BC_NONE :
            break;

         case BC_REFLECTION :
            W[ICu+1][j][k] = LES3DFsd_pState::Reflect(W[ICu][j][k],Grid.nfaceE(ICu,j,k));
            U[ICu+1][j][k] = W[ICu+1][j][k].U();
            W[ICu+2][j][k] = LES3DFsd_pState::Reflect(W[ICu-1][j][k],Grid.nfaceE(ICu,j,k));
            U[ICu+2][j][k] = W[ICu+2][j][k].U();
            break;

         case BC_NO_SLIP :
            W[ICu+1][j][k] = LES3DFsd_pState::NoSlip(W[ICu][j][k],
                                                     WoE[j][k],
                                                     Grid.nfaceE(ICu,j,k),
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[ICu+1][j][k] = W[ICu+1][j][k].U();
            W[ICu+2][j][k] = LES3DFsd_pState::NoSlip(W[ICu-1][j][k],
                                                     WoE[j][k],
                                                     Grid.nfaceE(ICu,j,k),
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[ICu+2][j][k] = W[ICu+2][j][k].U();
            break;

         case BC_MOVING_WALL :
            W[ICu+1][j][k] = LES3DFsd_pState::MovingWall(W[ICu][j][k],
                                                         WoE[j][k],
                                                         Grid.nfaceE(ICu,j,k),
                                                         MOVING_WALL_VELOCITY,
                                                         IPs.Pressure_Gradient,
                                                         FIXED_TEMPERATURE_WALL);
            U[ICu+1][j][k] = W[ICu+1][j][k].U();
            W[ICu+2][j][k] = LES3DFsd_pState::MovingWall(W[ICu-1][j][k], 
                                                         WoW[j][k], 
                                                         Grid.nfaceW(ICu,j,k),
                                                         MOVING_WALL_VELOCITY, 
                                                         IPs.Pressure_Gradient,
                                                         FIXED_TEMPERATURE_WALL);
            U[ICu+2][j][k] = W[ICu+2][j][k].U();
            break;

         case BC_FIXED_PRESSURE :
            W[ICu+1][j][k] = W[ICu][j][k];
            W[ICu+1][j][k].p = WoE[j][k].p;
            U[ICu+1][j][k] = W[ICu+1][j][k].U();
            
            W[ICu+2][j][k] = W[ICu-1][j][k];
            W[ICu+2][j][k].p = WoE[j][k].p; 
            U[ICu+2][j][k] = W[ICu+2][j][k].U();
            break;

         case BC_CONSTANT_EXTRAPOLATION :
            W[ICu+1][j][k] = W[ICu][j][k];
            U[ICu+1][j][k] = W[ICu+1][j][k].U();
            W[ICu+2][j][k] = W[ICu][j][k];
            U[ICu+2][j][k] = W[ICu+2][j][k].U();
            break;

         case BC_CHANNEL_OUTFLOW:
            dpdx = IPs.Pressure_Gradient.x; 
            // all constant extrapolation except pressure specified which linearly varys if there is pressure gradient
            dX = Grid.Cell[ICu+1][j][k].Xc - Grid.Cell[ICu][j][k].Xc; 
            W[ICu+1][j][k] = W[ICu][j][k]; 
            W[ICu+1][j][k].p = WoE[j][k].p-dpdx*dX.x;
            //      W[ICu+1][j][k].p = (WoE[j][k].p-2.0/3.0*W[ICu+1][j][k].rho*W[ICu+1][j][k].k)-dpdx*dX.x;
            U[ICu+1][j][k] = W[ICu+1][j][k].U( );
            
            dX = Grid.Cell[ICu+2][j][k].Xc - Grid.Cell[ICu][j][k].Xc; 
            W[ICu+2][j][k] = W[ICu][j][k];
            W[ICu+2][j][k].p = WoE[j][k].p -dpdx*dX.x; 	
            //W[ICu+2][j][k].p = (WoE[j][k].p -2.0/3.0*W[ICu+2][j][k].rho*W[ICu+2][j][k].k)-dpdx*dX.x; 	
            U[ICu+2][j][k] = W[ICu+2][j][k].U( );
            break;
            
         case BC_PERIODIC :
            W[ICu+1][j][k] = W[ICl+1][j][k];
            U[ICu+1][j][k] = U[ICl+1][j][k];
            W[ICu+2][j][k] = W[ICl+2][j][k];
            U[ICu+2][j][k] = U[ICl+2][j][k];
            break;

 	 case BC_INFLOW_SUBSONIC :
          // all fixed except pressure which is constant extrapolation
	    W[ICu+1][j][k] = WoE[j][k];
            W[ICu+1][j][k].p = W[ICu][j][k].p;
	    U[ICu+1][j][k] = W[ICu+1][j][k].U();

	    W[ICu+2][j][k] = WoE[j][k];
            W[ICu+2][j][k].p = W[ICu][j][k].p;
	    U[ICu+2][j][k] = W[ICu+2 ][j][k].U();
	    break;

	 case BC_OUTFLOW_SUBSONIC :
	  // all constant extrapolation except pressure which is fixed.
	    W[ICu+1][j][k] = W[ICu][j][k]; 
	    W[ICu+1][j][k].p = WoE[j][k].p;
	    U[ICu+1][j][k] = W[ICu+1][j][k].U();

	    W[ICu+2][j][k] = W[ICu][j][k]; 
	    W[ICu+2][j][k].p = WoE[j][k].p;
	    U[ICu+2][j][k] = W[ICu+2][j][k].U();
	    break;
            
         }//endofeastface
         
      } /* endfor */
       
 for ( k =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k )
      for ( i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
              
  // Prescribe North boundary conditions.
         switch( Grid.BCtypeN[i][k]) {

         case BC_NONE :
            break;
            
         case BC_REFLECTION :
            W[i][JCu+1][k] = LES3DFsd_pState::Reflect(W[i][JCu][k],Grid.nfaceN(i,JCu,k));
            U[i][JCu+1][k] = W[i][ JCu+1][k].U();
            W[i][JCu+2][k] = LES3DFsd_pState::Reflect(W[i][JCu-1][k],Grid.nfaceN(i,JCu,k));
            U[i][JCu+2][k] = W[i][JCu+2][k].U();
            break;

         case BC_NO_SLIP :
            W[i][JCu+1][k] = LES3DFsd_pState::NoSlip(W[i][JCu][k],
                                                     WoN[i][k],
                                                     Grid.nfaceN(i,JCu,k),
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[i][JCu+1][k] = W[i][JCu+1][k].U();
            W[i][JCu+2][k] = LES3DFsd_pState::NoSlip(W[i][JCu-1][k], 
                                                     WoN[i][k],
                                                     Grid.nfaceN(i,JCu,k),
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[i][JCu+2][k] = W[i][JCu+2][k].U();
            break;
            
         case BC_MOVING_WALL :
            W[i][JCu+1][k] = LES3DFsd_pState::MovingWall(W[i][JCu][k], 
                                                         WoN[i][k],
                                                         Grid.nfaceN(i,JCu,k),
                                                         MOVING_WALL_VELOCITY,
                                                         IPs.Pressure_Gradient,
                                                         FIXED_TEMPERATURE_WALL);
            U[i][JCu+1][k] = W[i][JCu+1][k].U();
            W[i][JCu+2][k] = LES3DFsd_pState::MovingWall(W[i][JCu-1][k], 
                                                         WoN[i][k],
                                                         Grid.nfaceN(i,JCu,k),
                                                         MOVING_WALL_VELOCITY,
                                                         IPs.Pressure_Gradient,
                                                         FIXED_TEMPERATURE_WALL);
            U[i][JCu+2][k] = W[i][JCu+2][k].U();
            break;

         case BC_FIXED_PRESSURE :
            W[i][ JCu+1][k] = W[i][JCu][k];
            W[i][ JCu+1][k].p = WoN[i][k].p;
            U[i][ JCu+1][k] = W[i][JCu+1][k].U();
          
            W[i][ JCu+2][k] = W[i][JCu-1][k];
            W[i][ JCu+2][k].p = WoN[i][k].p;
            U[i][ JCu+2][k] = W[i][JCu+2][k].U();
            break;
            
         case BC_CHANNEL_INFLOW:
            dpdy = IPs.Pressure_Gradient.y; 
            // all constant extrapolation except pressure specified which linearly varys if there is pressure gradient
            dX = Grid.Cell[i][JCu+1][k].Xc - Grid.Cell[i][JCu][k].Xc; 
            W[i][JCu+1][k].v.y = W[i][JCu][k].v.y; 
            W[i][JCu+1][k].p = WoN[i][k].p;
            //W[i][JCu+1][k].p = (WoN[i][k].p-2.0/3.0*W[i][JCu+1][k].rho*W[i][JCu+1][k].k)-dpdy*dX.y;
            
            U[i][JCu+1][k] = W[i][JCu+1][k].U( );
            
            dX = Grid.Cell[i][JCu+2][k].Xc - Grid.Cell[i][JCu][k].Xc; 
            W[i][JCu+2][k] = W[i][JCu][k];
            W[i][JCu+2][k].p = WoN[i][k].p;
            //W[i][JCu+2][k].p = (WoN[i][k].p -2.0/3.0*W[i][JCu+2][k].rho*W[i][JCu+2][k].k)-dpdy*dX.y; 	
            U[i][JCu+2][k] = W[i][JCu+2][k].U( );
            break;
            
         case BC_CONSTANT_EXTRAPOLATION :
            W[i][JCu+1][k] = W[i][JCu][k];
            U[i][JCu+1][k] = W[i][JCu+1][k].U();
            W[i][JCu+2][k] = W[i][JCu][k];
            U[i][JCu+2][k] = W[i][JCu+2][k].U();
            break;

         case BC_PERIODIC :
            W[i][JCu+1][k] = W[i][JCl+1][k];
            U[i][JCu+1][k] = U[i][JCl+1][k];
            W[i][JCu+2][k] = W[i][JCl+2][k];
            U[i][JCu+2][k] = U[i][JCl+2][k];
            break;

         case BC_INFLOW_SUBSONIC :
	   // all fixed except v.x (u) which is constant extrapolation
	    W[i][JCl-1][k] = WoS[i][k];
	    W[i][JCl-1][k].v.y = W[i][JCl][k].v.y;
	    U[i][JCl-1][k] = W[i][JCl-1][k].U();
	
	    W[i][JCl-2][k] = WoS[i][k];
	    W[i][JCl-2][k].v.y = W[i][JCl][k].v.y;
	    U[i][JCl-2][k] = W[i][JCl-2][k].U();
	    break;

         case BC_OUTFLOW_SUBSONIC :
	  // all constant extrapolation except pressure which is fixed.
	    W[i][JCl-1][k] = W[i][JCl][k];
	    W[i][JCl-1][k].p = WoS[i][k].p;
	    U[i][JCl-1][k] = W[i][JCl-1][k].U();
	  
	    W[i][JCl-2][k] = W[i][JCl][k];
	    W[i][JCl-2][k].p = WoS[i][k].p;
	    U[i][JCl-2][k] = W[i][JCl-2][k].U();
	    break;

         } /* endswitch */
    
         // Prescribe South boundary conditions.
         switch( Grid.BCtypeS[i][k]) {

         case BC_NONE :
            break;

          case BC_REFLECTION :
            W[i][JCl-1][k] = LES3DFsd_pState::Reflect(W[i][JCl][k],Grid.nfaceS(i,JCl,k));
            U[i][JCl-1][k] = W[i][JCl-1][k].U();
            W[i][JCl-2][k] = LES3DFsd_pState::Reflect(W[i][JCl+1][k],Grid.nfaceS(i,JCl,k));
            U[i][JCl-2][k] = W[i][JCl-2][k].U();
            break;

         case BC_NO_SLIP :
            W[i][JCl-1][k] = LES3DFsd_pState::NoSlip(W[i][JCl][k], 
                                                     WoS[i][k], 
                                                     Grid.nfaceS(i,JCl,k),
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[i][JCl-1][k] = W[i][JCl-1][k].U();
            W[i][JCl-2][k] = LES3DFsd_pState::NoSlip(W[i][JCl+1][k], 
                                                     WoS[i][k],
                                                     Grid.nfaceS(i,JCl,k),   
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[i][JCl-2][k] = W[i][JCl-2][k].U();
            break;
 
        case BC_MOVING_WALL :
            W[i][JCl-1][k] = LES3DFsd_pState::MovingWall(W[i][JCl][k], 
                                                         WoS[i][k], 
                                                         Grid.nfaceS(i,JCl,k),
                                                         MOVING_WALL_VELOCITY,
                                                         IPs.Pressure_Gradient, 
                                                         FIXED_TEMPERATURE_WALL);
            U[i][JCl-1][k] = W[i][JCl-1][k].U();
            W[i][JCl-2][k] = LES3DFsd_pState::MovingWall(W[i][JCl+1][k],
                                                         WoS[i][k],
                                                         Grid.nfaceS(i,JCl,k),
                                                         MOVING_WALL_VELOCITY, 
                                                         IPs.Pressure_Gradient,
                                                         FIXED_TEMPERATURE_WALL);
            U[i][JCl-2][k] =  W[i][JCl-2][k].U();
            break;

         case BC_FIXED_PRESSURE :
            W[i][JCl-1][k] = W[i][JCl][k];
            W[i][JCl-1][k].p = WoS[i][k].p;
            U[i][JCl-1][k] =  W[i][JCl-1][k].U();
            W[i][JCl-2][k] = W[i][JCl+1][k];
            W[i][JCl-2][k].p = WoS[i][k].p;
            U[i][JCl-2][k] =  W[i][JCl-2][k].U();
	    break;

         case BC_CONSTANT_EXTRAPOLATION :
            W[i][JCl-1][k] = W[i][JCl][k];
            U[i][JCl-1][k] = W[i][JCl-1][k].U();
            W[i][JCl-2][k] = W[i][JCl][k];
            U[i][JCl-2][k] = W[i][JCl-2][k].U();
	    break;

          case BC_CHANNEL_OUTFLOW:
            dpdy = IPs.Pressure_Gradient.y; 
            //for turbulent channel flow
            // k and omega are constant extrapolation, p linearly varys based on constant pressure gradient 
            dX = Grid.Cell[i][JCl-1][k].Xc - Grid.Cell[i][JCl][k].Xc; 
            W[i][JCl-1][k] = WoS[i][k];
            W[i][JCl-1][k].p = WoS[i][k].p -dpdy*dX.y;
            //W[i][JCl-1][k].p = (WoS[i][k].p )-dpdy*dX.y;
            U[i][JCl-1][k] = W[i][JCl-1][k].U( );
            
            dX = Grid.Cell[i][JCl-2][k].Xc - Grid.Cell[i][JCl][k].Xc;
            W[i][JCl-2][k].p = WoS[i][k].p - dpdy*dX.y;
            //W[i][JCl-2][k].p = (WoS[i][k].p - 2.0/3.0*W[i][JCl-2][k].rho*W[i][JCl-2][k].k) - dpdy*dX.y;
            U[i][JCl-2][k] = W[i][JCl-2][k].U( );
            break;

         case BC_PERIODIC :
            W[i][JCl-1][k] = W[i][JCu-1][k];
            U[i][JCl-1][k] = U[i][JCu-1][k];
            W[i][JCl-2][k] = W[i][JCu-2][k];
            U[i][JCl-2][k] = U[i][JCu-2][k];
            break;

         case BC_INFLOW_SUBSONIC :
	// all fixed except v.x (u) which is constant extrapolation
	    W[i][JCl-1][k] = WoS[i][k];
	    W[i][JCl-1][k].v.y = W[i][JCl][k].v.y;
	    U[i][JCl-1][k] = W[i][JCl-1][k].U();
	
	    W[i][JCl-2][k] = WoS[i][k];
	    W[i][JCl-2][k].v.y = W[i][JCl][k].v.y;
	    U[i][JCl-2][k] = W[i][JCl-2][k].U();
	    break;

         case BC_OUTFLOW_SUBSONIC :
	// all constant extrapolation except pressure which is fixed.
	    W[i][JCl-1][k] = W[i][JCl][k];
	    W[i][JCl-1][k].p = WoS[i][k].p;
	    U[i][JCl-1][k] = W[i][JCl-1][k].U();
	  
	    W[i][JCl-2][k] = W[i][JCl][k];
	    W[i][JCl-2][k].p = WoS[i][k].p;
	    U[i][JCl-2][k] = W[i][JCl-2][k].U();
	    break;
            
         } /* endswitch */
  } /* endfor */
 
 for ( j =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j )
      for ( i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {

         // Prescribe Bottom boundary conditions.
         switch( Grid.BCtypeB[i][j]) {
       
         case BC_NONE :
            break;

         case BC_REFLECTION :
            W[i][j][KCl-1] = LES3DFsd_pState::Reflect(W[i][j][KCl],Grid.nfaceBot(i,j,KCl));
            U[i][j][KCl-1] = W[i][j][KCl-1].U();
            W[i][j][KCl-2] = LES3DFsd_pState::Reflect(W[i][j][KCl+1],Grid.nfaceBot(i,j,KCl));
            U[i][j][KCl-2] = W[i][j][KCl-2].U();
            break;

         case BC_NO_SLIP :
            W[i][j][KCl-1] = LES3DFsd_pState::NoSlip(W[i][j][KCl], 
                                                     WoB[i][j],
                                                     Grid.nfaceBot(i,j,KCl),
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-1] = W[i][j][KCl-1].U(); 
            W[i][j][KCl-2] = LES3DFsd_pState::NoSlip(W[i][j][KCl+1], 
                                                     WoB[i][j],
                                                     Grid.nfaceBot(i,j,KCl),
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-2] = W[i][j][KCl-2].U();
            break;
            
         case BC_MOVING_WALL :
            W[i][j][KCl-1] = LES3DFsd_pState::MovingWall(W[i][j][KCl], 
                                                         WoB[i][j],
                                                         Grid.nfaceBot(i,j,KCl),
                                                         MOVING_WALL_VELOCITY,
                                                         IPs.Pressure_Gradient,
                                                         FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-1] = W[i][j][KCl-1].U();
            W[i][j][KCl-2] = LES3DFsd_pState::MovingWall(W[i][j][KCl+1], 
                                                         WoB[i][j],
                                                         Grid.nfaceBot(i,j ,KCl),
                                                         MOVING_WALL_VELOCITY,
                                                         IPs.Pressure_Gradient,
                                                         FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-2] =  W[i][j][KCl-2].U();
            break;

         case BC_FIXED_PRESSURE :
            W[i][j][KCl-1] = W[i][j][KCl];
            W[i][j][KCl-1].p = WoB[i][j].p;
            U[i][j][KCl-1] = W[i][j][KCl-1].U();
            W[i][j][KCl-2] = W[i][j][KCl+1];
            W[i][j][KCl-2].p = WoB[i][j].p;
            U[i][j][KCl-2] = W[i][j][KCl-2].U();
            break;

         case BC_CONSTANT_EXTRAPOLATION :
            W[i][j][KCl-1] = W[i][j][KCl];
            U[i][j][KCl-1] = W[i][j][KCl-1].U();
            W[i][j][KCl-2] = W[i][j][KCl];
            U[i][j][KCl-2] = W[i][j][KCl-2].U();
            break;

         case BC_CHANNEL_INFLOW:
            dpdz = IPs.Pressure_Gradient.z; 
            //for turbulent channel flow
            // k and omega are constant extrapolation, p linearly varys based on constant pressure gradient 
            dX = Grid.Cell[i][j][KCl-1].Xc - Grid.Cell[i][j][KCl].Xc; 
    
            if((j==JCl-1 || j==JCl-2 || j== JCu+1 || j== JCu+2 ) && 
               ((Grid.BCtypeN[i][KCl-1] == BC_NO_SLIP) ||
                (Grid.BCtypeN[i][KCl-2] == BC_NO_SLIP) ||
                (Grid.BCtypeS[i][KCl-1] == BC_NO_SLIP) ||
                (Grid.BCtypeS[i][KCl-2] == BC_NO_SLIP))){
               //do not overwrite the solid no-slip wall boundary condition 
            }else{
               
               W[i][j][KCl-1] = WoB[i][j];
               W[i][j][KCl-1].v.z = W[i][j][KCl].v.z;
               W[i][j][KCl-1].p = WoB[i][j].p - dpdz*dX.z;
               // W[i][j][KCl-1].p = (WoB[i][j].p  - 2.0/3.0*W[i][j][KCl-1].rho* W[i][j][KCl-1].k)-dpdz*dX.z;
               U[i][j][KCl-1] = W[i][j][KCl-1].U( );
               
               dX = Grid.Cell[i][j][KCl-2].Xc - Grid.Cell[i][j][KCl].Xc;
               W[i][j][KCl-2] = WoB[i][j];
               W[i][j][KCl-2].v.z = W[i][j][KCl].v.z;
               W[i][j][KCl-2].p = WoB[i][j].p - dpdz*dX.z;
               // W[i][j][KCl-2].p = (WoB[i][j].p - 2.0/3.0*W[i][j][KCl-2].rho* W[i][j][KCl-2].k ) - dpdz*dX.z;
               U[i][j][KCl-2] = W[i][j][KCl-2].U( );
               
            }
            
            break;

         case BC_PERIODIC :
            W[i][j][KCl-1] = W[i][j][KCu-1];
            U[i][j][KCl-1] = U[i][j][KCu-1];
            W[i][j][KCl-2] = W[i][j][KCu-2];
            U[i][j][KCl-2] = U[i][j][KCu-2];
            break;

         case BC_INFLOW_SUBSONIC :
	// all fixed except v.x (u) which is constant extrapolation
            W[i][j][KCl-1] = W[i][j][KCu-1];
            W[i][j][KCl-1].v.z = WoB[i][j].v.z;
            U[i][j][KCl-1] = W[i][j][KCu-1].U();
            W[i][j][KCl-2] = W[i][j][KCu-2];
            W[i][j][KCl-2].v.z = WoB[i][j].v.z;
            U[i][j][KCl-2] = W[i][j][KCu-2].U();
            break;

         case BC_OUTFLOW_SUBSONIC :
	// all constant extrapolation except pressure which is fixed.
            W[i][j][KCl-1] = W[i][j][KCu-1];
            W[i][j][KCl-1].p = WoB[i][j].p;
            U[i][j][KCl-1] = W[i][j][KCu-1].U();
            W[i][j][KCl-2] = W[i][j][KCu-2];
            W[i][j][KCl-2].p = WoB[i][j].p;
            U[i][j][KCl-2] = W[i][j][KCu-2].U();
            break;

         } /* endswitch */
       
  // Prescribe Top boundary conditions.
         switch( Grid.BCtypeT[i][j]) {
       
         case BC_NONE :
            break;
            
         case BC_REFLECTION :
            W[i][j][KCu+1] = LES3DFsd_pState::Reflect(W[i][j][KCu],Grid.nfaceTop(i,j,KCu));
            U[i][j][KCu+1] = W[i][j][KCu+1].U();
            W[i][j][KCu+2] = LES3DFsd_pState::Reflect(W[i][j][KCu-1],Grid.nfaceTop(i,j,KCu));
            U[i][j][KCu+2] = W[i][j][KCu+2].U();
            break;

         case BC_NO_SLIP :
            W[i][j][KCu+1] = LES3DFsd_pState::NoSlip(W[i][j][KCu], 
                                                     WoT[i][j],
                                                     Grid.nfaceTop(i,j,KCu),
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[i][j][KCu+1] = W[i][j][KCu+1].U();
            W[i][j][KCu+2] = LES3DFsd_pState::NoSlip(W[i][j][KCu-1], 
                                                     WoT[i][j],
                                                     Grid.nfaceTop(i,j,KCu),
                                                     IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
            U[i][j][KCu+2] = W[i][j][KCu+2].U();
            break;
            
         case BC_MOVING_WALL :
            W[i][j][KCu+1] = LES3DFsd_pState::MovingWall(W[i][j][KCu], 
                                                         WoT[i][j],
                                                         Grid.nfaceTop(i,j,KCu),
                                                         MOVING_WALL_VELOCITY,
                                                         IPs.Pressure_Gradient,
                                                         FIXED_TEMPERATURE_WALL);
            U[i][j][KCu+1] = W[i][j][KCu+1].U();
            W[i][j][KCu+2] = LES3DFsd_pState::MovingWall(W[i][j][KCu-1], 
                                                         WoT[i][j],
                                                         Grid.nfaceTop(i,j,KCu),
                                                         MOVING_WALL_VELOCITY,
                                                         IPs.Pressure_Gradient,
                                                         FIXED_TEMPERATURE_WALL);
            U[i][j][KCu+2] = W[i][j][KCu+2].U();
            break;

         case BC_FIXED_PRESSURE :
            W[i][j][KCu+1] = W[i][j][KCu];
            W[i][j][KCu+1].p = WoT[i][j].p;
            U[i][j][KCu+1] =  W[i][j][KCu+1].U();
            W[i][j][KCu+2] =  W[i][j][KCu-1];
            W[i][j][KCu+2].p = WoT[i][j].p;
            U[i][j][KCu+2] =  W[i][j][KCu+2].U();
            break;

         case BC_CONSTANT_EXTRAPOLATION :
            W[i][j][KCu+1] = W[i][j][KCu];
            U[i][j][KCu+1] = W[i][j][KCu+1].U();
            W[i][j][KCu+2] = W[i][j][KCu];
            U[i][j][KCu+2] = W[i][j][KCu+2].U();
            break;

         case BC_CHANNEL_OUTFLOW:
            dpdz = IPs.Pressure_Gradient.z; 
            // all constant extrapolation except pressure specified which linearly varys if there is pressure gradient
            dX = Grid.Cell[i][j][KCu+1].Xc - Grid.Cell[i][j][KCu].Xc; 
            W[i][j][KCu+1] = W[i][j][KCu]; 
            W[i][j][KCu+1].p = (WoT[i][j].p)-dpdz*dX.z;
            U[i][j][KCu+1] = W[i][j][KCu+1].U( );
            
            dX = Grid.Cell[i][j][KCu+2].Xc - Grid.Cell[i][j][KCu].Xc; 
            W[i][j][KCu+2] = W[i][j][KCu];
            W[i][j][KCu+2].p = (WoT[i][j].p)-dpdz*dX.z; 	
            U[i][j][KCu+2] = W[i][j][KCu+2].U( );
            break;

         case BC_PERIODIC :
            W[i][j][KCu+1] = W[i][j][KCl+1];
            U[i][j][KCu+1] = U[i][j][KCl+1];
            W[i][j][KCu+2] = W[i][j][KCl+2];
            U[i][j][KCu+2] = U[i][j][KCl+2];
            break;

         case BC_INFLOW_SUBSONIC :
	// all fixed except v.x (u) which is constant extrapolation
            W[i][j][KCl+1] = W[i][j][KCu+1];
            W[i][j][KCl+1].v.z = WoT[i][j].v.z;
            U[i][j][KCl+1] = W[i][j][KCu+1].U();
            W[i][j][KCl+2] = W[i][j][KCu+2];
            W[i][j][KCl+2].v.z = WoT[i][j].v.z;
            U[i][j][KCl+2] = W[i][j][KCu+2].U();
            break;

         case BC_OUTFLOW_SUBSONIC :
	// all constant extrapolation except pressure which is fixed.
            W[i][j][KCl+1] = W[i][j][KCu+1];
            W[i][j][KCl+1].p = WoT[i][j].p;
            U[i][j][KCl+1] = W[i][j][KCu+1].U();
            W[i][j][KCl+2] = W[i][j][KCu+2];
            W[i][j][KCl+2].p = WoT[i][j].p;
            U[i][j][KCl+2] = W[i][j][KCu+2].U();
            break;

         } /* endswitch */
      } /* endfor */
}

/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Euler time stepping scheme) for the    *
 * specified hexa solution block according to           *
 * the Courant-Friedrichs-Lewy condition.               *
 *                                                      *
 ********************************************************/
template<>
double Hexa_Block<LES3DFsd_pState,LES3DFsd_cState>::
CFL(Input_Parameters<LES3DFsd_pState,LES3DFsd_cState> &IPs){
   
   int i, j, k;
   double dtMin, d_i, d_j, d_k, v_i, v_j, v_k, a, dt_vis, nv, nv_t;
   double mr, aa_i, aa_j, aa_k;
   double length_n, delta_n, dTime;
   
   dtMin = MILLION;
   
   for ( k  =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k ) 
      for ( j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) 
         for ( i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
            if (i <  ICl || i >  ICu ||
                j <  JCl || j >  JCu || 
                k <  KCl || k >  KCu) {
               dt[i][j][k] = ZERO;
            } else {
               d_i = TWO*(Grid.volume(i,j,k)/
                          (Grid.AfaceE(i, j, k)+ Grid.AfaceW(i, j, k)));
               d_j = TWO*( Grid.volume(i, j, k)/
                           (Grid.AfaceN(i, j, k)+ Grid.AfaceS(i, j, k)));
               d_k = TWO*( Grid.volume(i, j, k)/
                           (Grid.AfaceTop(i, j, k)+ Grid.AfaceBot(i, j, k)));
               v_i = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                           (Grid.nfaceE(i, j, k)- Grid.nfaceW(i, j, k)));
               v_j = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                           ( Grid.nfaceN(i, j, k)- Grid.nfaceS(i, j, k)));
               v_k = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                            (Grid.nfaceTop(i, j, k)- Grid.nfaceBot(i, j, k)));
               
                  
               length_n = max(max(IPs.Grid_IP.Box_Length, IPs.Grid_IP.Box_Width), IPs.Grid_IP.Box_Height);  
	       delta_n = min(min(fabs(d_i),fabs(d_j)),fabs(d_k));

	     //no preconditioning
             if(IPs.Preconditioning == 0){
      	    
               a =  W[i][j][k].a();
               dt[i][j][k] = min(min(d_i/(a+fabs(v_i)), d_j/(a+fabs(v_j))),
                                 d_k/(a+fabs(v_k)));
               
	       //Low Mach Number Preconditioning
 	     } else if(IPs.Preconditioning == 1) { 
	       dt[i][j][k] = min(min(d_i/W[i][j][k].u_plus_aprecon(fabs(v_i),
     							           delta_n,
                                                                   length_n,
                                                                   dTime),
	       			     d_j/W[i][j][k].u_plus_aprecon(fabs(v_j),
					  		           delta_n,
                                                                   length_n,
								   dTime)),
	       			     d_k/W[i][j][k].u_plus_aprecon(fabs(v_k),
					  		           delta_n,
                                                                   length_n,
                                                                   dTime));
	     } /* endif */
                     
               if (IPs.i_Flow_Type != FLOWTYPE_INVISCID) {  
                  nv = W[i][j][k].mu()/W[i][j][k].rho;

               if (IPs.i_Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K) {  
		 nv_t = W[i][j][k].mu_t(dWdx[i][j][k],
                                        dWdy[i][j][k],
                                        dWdz[i][j][k],
                                        Flow_Type,
                                        Grid.volume(i, j, k)) / W[i][j][k].rho; 
                     nv = max(nv, nv_t);
	       }
                  
                  dt_vis = min(min((d_i*d_i)/(3.0*nv), (d_j*d_j)/(3.0*nv)), (d_k*d_k)/(3.0*nv)); 
                  dt[i][j][k]  = min(dt_vis, dt[i][j][k]);
		  /*
             int nn = NumVar();
             DenseMatrix dSwdU(nn,nn);
             dSwdU.zero();
             double max_diagonal = ONE;

             W[i][j][k].SemiImplicitSourceJacobi(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
                                                 d_dWdx_dW[i][j][k][0],d_dWdy_dW[i][j][k][0],
                                                 d_dWdz_dW[i][j][k][0],dSwdU,Grid.volume(i,j,k)); 

             if(IPs.Preconditioning == 1){
 	      DenseMatrix Pinv(nn,nn);
 	      Pinv.zero();
 	      W[i][j][k].Low_Mach_Number_Preconditioner_Inverse(Pinv, delta_n, length_n, dTime);	 
 	      dSwdU = Pinv*dSwdU;
             }
             for(int ii=0; ii<nn; ii++){
              max_diagonal = max(max_diagonal,fabs(dSwdU(ii,ii)));
             }
 	       dt[i][j][k] = min(HALF/max_diagonal, dt[i][j][k]);
		  */
            }       
             dtMin = min(dtMin,  dt[i][j][k]);
               
	    } /* endif */
            
	 } /* endfor */
   
   for ( k  =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k ) 
      for ( j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) 
         for ( i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
            if (i <  ICl || i >  ICu ||
                j <  JCl || j >  JCu || 
                k <  KCl || k >  KCu) {
               dt[i][j][k] = dtMin;
            } /* endif */
         } /* endfor */
   
   /* Return the global time step. */
   
   return (dtMin);
}

template<>
int Hexa_Block<LES3DFsd_pState,LES3DFsd_cState>::
dUdt_Multistage_Explicit(const int i_stage,
                         Input_Parameters<LES3DFsd_pState,LES3DFsd_cState> &IPs) {

   int i, j, k,  k_residual;
   double omega;
   Vector3D dX;
   
   LES3DFsd_pState Wl, Wr;
   LES3DFsd_cState Flux, Temp;

   LES3DFsd_cState U_VACUUM;
   U_VACUUM.Vacuum();
   LES3DFsd_pState W_VACUUM;
   W_VACUUM.Vacuum();

   double delta_n;

   Wl = W_VACUUM;
   Wr = W_VACUUM;
   Flux = U_VACUUM;
   Temp = U_VACUUM;


   // Evaluate the time step fraction and residual storage location for the stage. 
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


   // compute y+ etc.;
   //   Wall_Shear( );
   
    /* Compute viscous stresses and heat conduction vector if using cell-centered methods. */
//     if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
//         Input_Parameters.i_Viscous_Flux_Evaluation == VISCOUS_RECONSTRUCTION_ARITHMETIC) {
//        Viscous_Calculations(SolnBlk);
//     }

   /* Evaluate the time rate of change of the solution
      (i.e., the solution residuals) using a second-order
      limited upwind scheme with a variety of flux functions. */
   
   // Add i-direction (zeta-direction) fluxes.
   for ( k = KCl-1; k <= KCu+1; ++k ) {
      for ( j = JCl-1; j <= JCu+1; ++j ) {
         if ( i_stage == 1 ) {
            Uo[ICl-1][j][k] = U[ICl-1][j][k];
            dUdt[ICl-1][j][k][k_residual] = U_VACUUM;
         } else {
            dUdt[ICl-1][j][k][k_residual] = U_VACUUM;
         } /* endif */
         
         for ( i = ICl-1; i <= ICu; ++i ) {
            if ( i_stage == 1 ) {
               Uo[i+1][j][k] = U[i+1][j][k];
               dUdt[i+1][j][k][k_residual] = U_VACUUM;
            } else if (( j > JCl-1 && j < JCu+1 ) 
                       && (k > KCl-1 && k < KCu+1)) {
               switch(IPs.i_Time_Integration) {
               case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR:
                  //  dUdt[i+1][j][k][k_residual] = 
                  //  dUdt[i+1][j][k][k_residual];
                  break;
               case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA:
		  if (IPs.N_Stage == 2) {
                     // dUdt[i+1][j][k_residual] = 
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

            if (( j > JCl-1 && j < JCu+1 ) && 
                ( k > KCl-1 && k < KCu+1 )) {
               /* Evaluate the cell interface i-direction fluxes. */
               if (i ==  ICl-1 && (Grid.BCtypeW[j][k] == BC_REFLECTION ||
                                   Grid.BCtypeW[j][k] == BC_NO_SLIP||
                                   Grid.BCtypeW[j][k] == BC_MOVING_WALL) ) {
                  dX = Grid.xfaceW(i+1,j,k)- Grid.Cell[i+1][j][k].Xc;
                  Wr = W[i+1][j][k] + 
                     ( phi[i+1][j][k]^dWdx[i+1][j][k])*dX.x +
                     ( phi[i+1][j][k]^dWdy[i+1][j][k])*dX.y +
                     ( phi[i+1][j][k]^dWdz[i+1][j][k])*dX.z;
                                 
                  if ( Grid.BCtypeW[j][k] == BC_REFLECTION) {
                     Wl = LES3DFsd_pState::Reflect(Wr,Grid.nfaceW(i+1,j,k));
                  }
                  if ( Grid.BCtypeW[j][k] == BC_NO_SLIP) {
                     Wl = LES3DFsd_pState::NoSlip(Wr,WoW[j][k],
                                                  Grid.nfaceW(i+1,j,k), 
                                                  IPs.Pressure_Gradient,
                                                  FIXED_TEMPERATURE_WALL);
                                                         
                  }
                  if ( Grid.BCtypeW[j][k] == BC_MOVING_WALL) {
                     Wl = LES3DFsd_pState::MovingWall(Wr,WoW[j][k],
                                                      Grid.nfaceW(i+1,j,k),
                                                      IPs.Moving_Wall_Velocity,
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
                  }
               } else if (i ==  ICu && 
                         (Grid.BCtypeE[j][k] == BC_REFLECTION||
                          Grid.BCtypeE[j][k] == BC_NO_SLIP||
                          Grid.BCtypeE[j][k] == BC_MOVING_WALL)) {
                  dX = Grid.xfaceE(i,j,k)- Grid.Cell[i][j][k].Xc;
                  Wl = W[i][j][k] + 
                     ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                     ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                     ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                  if ( Grid.BCtypeE[j][k] == BC_REFLECTION) {
                     Wr = LES3DFsd_pState::Reflect(Wl,Grid.nfaceE(i,j,k));
                  } 
                  if ( Grid.BCtypeE[j][k] == BC_NO_SLIP) {
                     Wr = LES3DFsd_pState::NoSlip(Wl,WoE[j][k],
                                                  Grid.nfaceE(i,j,k), 
                                                  IPs.Pressure_Gradient, 
                                                  FIXED_TEMPERATURE_WALL);
                  } 
                  if ( Grid.BCtypeE[j][k] == BC_MOVING_WALL) {
                     Wr = LES3DFsd_pState::MovingWall(Wl,WoE[j][k],
                                                      Grid.nfaceE(i,j,k),  
                                                      IPs.Moving_Wall_Velocity,
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
                  }
               } else { 
                  dX = Grid.xfaceE(i,j,k)- Grid.Cell[i][j][k].Xc;
                  Wl = W[i][j][k] + 
                     ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                     ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                     ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                  dX = Grid.xfaceW(i+1,j,k)- Grid.Cell[i+1][j][k].Xc;
                  Wr = W[i+1][j][k] + 
                     ( phi[i+1][j][k]^ dWdx[i+1][j][k])*dX.x +
                     ( phi[i+1][j][k]^ dWdy[i+1][j][k])*dX.y +
                     ( phi[i+1][j][k]^ dWdz[i+1][j][k])*dX.z;
               } /* endif */
               
	  // Spacing for Preconditioner 
	  if(IPs.Flow_Type != FLOWTYPE_INVISCID && IPs.Preconditioning){
 	     delta_n = min(min(TWO*(Grid.volume(i,j,k)/(Grid.AfaceE(i,j,k)+Grid.AfaceW(i,j,k))),
  			       TWO*(Grid.volume(i,j,k)/(Grid.AfaceN(i,j,k)+Grid.AfaceS(i,j,k)))),
                               TWO*(Grid.volume(i,j,k)/(Grid.AfaceTop(i,j,k)+Grid.AfaceBot(i,j,k))));
	  }
               switch(IPs.i_Flux_Function) {
               case FLUX_FUNCTION_HLLE :
                  Flux = LES3DFsd_pState::FluxHLLE_n(Wl, Wr, Grid.nfaceE(i,j,k));
                  break;
               case FLUX_FUNCTION_ROE :
                  Flux = LES3DFsd_pState::FluxRoe_n(Wl, Wr, Grid.nfaceE(i,j,k));
                  break;
  	       case FLUX_FUNCTION_AUSM_PLUS_UP :
 		  Flux = LES3DFsd_pState::FluxAUSMplus_up_n(Wl, Wr, Grid.nfaceE(i,j,k));
 	          break;
               } /* endswitch */
               // add viscous flux in i direction
               Flux -=  LES3DFsd_pState::FluxViscous_n(Wl, Wr,
                                                       W[i][j][k],W[i+1][j][k], 
                                                       dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
                                                       dWdx[i+1][j][k],dWdy[i+1][j][k],dWdz[i+1][j][k],
                                                       Grid.nfaceE(i,j,k),Grid.Voe(i,j,k),  
                                                       Grid.delta_oe(i,j,k),Grid.volume(i,j,k), 
	                                               Grid.volume(i+1,j,k),Flow_Type);
               /* Evaluate cell-averaged solution changes. */
               dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*Flux*Grid.AfaceE(i,j,k)/Grid.volume(i,j,k);
               dUdt[i+1][j][k][k_residual] += (IPs.CFL_Number*dt[i+1][j][k])*Flux*Grid.AfaceW(i+1,j,k)/Grid.volume(i+1,j,k);
               /* Include source terms associated with turbulence model */
                dUdt[i][j][k][k_residual] += (IPs.CFL_Number*dt[i][j][k])*LES3DFsd_pState::Sturbulence(W[i][j][k],
                                                                                                       dWdx[i][j][k],
                                                                                                       dWdy[i][j][k], 
                                                                                                       dWdz[i][j][k],
                                                                                                       Flow_Type,
                                                                                                       Grid.volume(i,j,k));
          /* Include physical time derivative for dual time stepping */
//           if (IPs.Dual_Time_Stepping) {
//             //cout << "Source[k] dual time" << endl;
//             dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*
//                                               W[i][j][k].S_dual_time_stepping(U[i][j][k], Ut[i][j][k], Uold[i][j][k], 
//                                                                                    dTime, IPs.first_step);
//           }

	   /* Save west and east face boundary flux. */
               
               //    if (i ==  ICl-1) {
//                     FluxW[j] = -Flux* Grid.lfaceW(i+1, j);
//                 } else if (i ==  ICu) {
//                     FluxE[j] = Flux* Grid.lfaceE(i, j);
//                 } /* endif */ 
               

            } /* endif */
         } /* endfor */
          
         if (( j > JCl-1 && j < JCu+1 ) && 
             ( k > KCl-1 && k < KCu+1 )) {
            dUdt[ICl-1][j][k][k_residual] = U_VACUUM;
            dUdt[ICu+1][j][k][k_residual] = U_VACUUM;
         } /* endif */
      } /* endfor */
   } /* endfor */
 
   // Add j-direction (eta-direction) fluxes.
   for ( k = KCl; k <= KCu; ++k ) {
      for ( i = ICl; i <= ICu; ++i ) {
         for ( j = JCl-1; j <= JCu; ++j ) {
            /* Evaluate the cell interface j-direction fluxes. */
            if (j == JCl-1 && 
                ( Grid.BCtypeS[i][k] == BC_REFLECTION ||
                  Grid.BCtypeS[i][k] == BC_NO_SLIP||
                  Grid.BCtypeS[i][k] == BC_MOVING_WALL)) {
               dX = Grid.xfaceS(i,j+1,k)- Grid.Cell[i][j+1][k].Xc;
               Wr = W[i][j+1][k] +
                  ( phi[i][j+1][k]^ dWdx[i][j+1][k])*dX.x +
                  ( phi[i][j+1][k]^ dWdy[i][j+1][k])*dX.y+
                  ( phi[i][j+1][k]^ dWdz[i][j+1][k])*dX.z;
               if ( Grid.BCtypeS[i][k] == BC_REFLECTION) {
                  Wl = LES3DFsd_pState::Reflect(Wr,Grid.nfaceS(i,j+1,k));
                }
                if ( Grid.BCtypeS[i][k] == BC_NO_SLIP) {
                  Wl = LES3DFsd_pState::NoSlip(Wr,WoS[i][k],
                                               Grid.nfaceS(i,j+1,k),
                                               IPs.Pressure_Gradient,
                                               FIXED_TEMPERATURE_WALL);
                }
                if ( Grid.BCtypeS[i][k] == BC_MOVING_WALL) {
                  Wl = LES3DFsd_pState::MovingWall(Wr,WoS[i][k],
                                                   Grid.nfaceS(i,j+1,k),
                                                   IPs.Moving_Wall_Velocity,
                                                   IPs.Pressure_Gradient,
                                                   FIXED_TEMPERATURE_WALL);
                }
             } else if ( j ==  JCu && 
                       ( Grid.BCtypeN[i][k] == BC_REFLECTION ||
                         Grid.BCtypeN[i][k] == BC_NO_SLIP||
                         Grid.BCtypeN[i][k] == BC_MOVING_WALL)) {
                dX = Grid.xfaceN(i,j,k)- Grid.Cell[i][j][k].Xc;
                Wl = W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x+
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y+
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                if ( Grid.BCtypeN[i][k] == BC_REFLECTION) {
                   Wr = LES3DFsd_pState::Reflect(Wl,Grid.nfaceN(i,j,k));
                }
                if ( Grid.BCtypeN[i][k] == BC_NO_SLIP) {
                   Wr = LES3DFsd_pState::NoSlip(Wl,WoN[i][k],
                                                Grid.nfaceN(i,j,k),
                                                IPs.Pressure_Gradient,
                                                FIXED_TEMPERATURE_WALL);
                }
                if ( Grid.BCtypeN[i][k] == BC_MOVING_WALL) {
                   Wr = LES3DFsd_pState::MovingWall(Wl,WoN[i][k],
                                                    Grid.nfaceN(i,j,k),
                                                    IPs.Moving_Wall_Velocity,
                                                    IPs.Pressure_Gradient,
                                                    FIXED_TEMPERATURE_WALL);
                }
             } else {
                dX = Grid.xfaceN(i,j,k)- Grid.Cell[i][j][k].Xc;
                Wl = W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                dX =  Grid.xfaceS(i,j+1,k)- Grid.Cell[i][j+1][k].Xc;
                Wr =  W[i][j+1][k] +
                   ( phi[i][j+1][k]^ dWdx[i][j+1][k])*dX.x +
                   ( phi[i][j+1][k]^ dWdy[i][j+1][k])*dX.y +
                   ( phi[i][j+1][k]^ dWdz[i][j+1][k])*dX.z;
          } /* endif */

	  // Spacing for Preconditioner 
	  if(IPs.Flow_Type != FLOWTYPE_INVISCID && IPs.Preconditioning){
	    delta_n = min(min(TWO*(Grid.volume(i,j,k)/(Grid.AfaceE(i,j,k)+Grid.AfaceW(i,j,k))),
			      TWO*(Grid.volume(i,j,k)/(Grid.AfaceN(i,j,k)+Grid.AfaceS(i,j,k)))),
                              TWO*( Grid.volume(i,j,k)/(Grid.AfaceTop(i,j,k)+Grid.AfaceBot(i,j,k))));
	  }
             switch(IPs.i_Flux_Function) {
             case FLUX_FUNCTION_HLLE :
                Flux =  LES3DFsd_pState::FluxHLLE_n(Wl, Wr, Grid.nfaceN(i,j,k));
                break;
             case FLUX_FUNCTION_ROE :
                Flux =  LES3DFsd_pState::FluxRoe_n(Wl, Wr, Grid.nfaceN(i,j,k));
                break;
             case FLUX_FUNCTION_AUSM_PLUS_UP :
		Flux =  LES3DFsd_pState::FluxAUSMplus_up_n(Wl, Wr, Grid.nfaceN(i,j,k));
	        break;
             } /* endswitch */
             // add viscous flux in j direction
             Flux -= LES3DFsd_pState::FluxViscous_n(Wl,Wr, 
                                                    W[i][j][k],W[i][j+1][k],
                                                    dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
                                                    dWdx[i][j+1][k],dWdy[i][j+1][k],dWdz[i][j+1][k],
                                                    Grid.nfaceN(i,j,k),Grid.Von(i,j,k),
                                                    Grid.delta_on(i,j,k),Grid.volume(i,j,k), 
                                                    Grid.volume(i,j+1,k),Flow_Type);
             /* Evaluate cell-averaged solution changes. */
             dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*Flux*Grid.AfaceN(i,j,k)/Grid.volume(i,j,k);
             dUdt[i][j+1][k][k_residual] += (IPs.CFL_Number*dt[i][j+1][k])*Flux*Grid.AfaceS(i,j+1,k)/Grid.volume(i,j+1,k);
          /* Save south and north face boundary flux. */
          // if (j ==  JCl-1) {
//               FluxS[i] = -Flux* Grid.lfaceS(i, j+1);
//           } else if (j ==  JCu) {
//               FluxN[i] = Flux* Grid.lfaceN(i, j);
//           } /* endif */

	 } /* endfor */
           dUdt[i][JCl-1][k][k_residual] = U_VACUUM;
           dUdt[i][JCu+1][k][k_residual] = U_VACUUM;
      } /* endfor */
   } /* endfor */

   // Add k-direction (gamma-direction) fluxes.
   for ( i = ICl; i <= ICu; ++i ) {
      for ( j = JCl; j <= JCu; ++j ){
         for ( k = KCl-1; k <= KCu; ++k )  {
            /* Evaluate the cell interface j-direction fluxes. */
            if ( k == KCl-1 && 
               ( Grid.BCtypeB[i][j] == BC_REFLECTION  || 
                 Grid.BCtypeB[i][j] == BC_NO_SLIP ||
                 Grid.BCtypeB[i][j] == BC_MOVING_WALL)) {
                
               dX = Grid.xfaceBot(i,j,k+1)- Grid.Cell[i][j][k+1].Xc;
               Wr = W[i][j][k+1] +
                  ( phi[i][j][k+1]^ dWdx[i][j][k+1])*dX.x +
                  ( phi[i][j][k+1]^ dWdy[i][j][k+1])*dX.y+
                  ( phi[i][j][k+1]^ dWdz[i][j][k+1])*dX.z;

               if ( Grid.BCtypeB[i][j] == BC_REFLECTION) {
                  Wl = LES3DFsd_pState::Reflect(Wr,Grid.nfaceBot(i,j,k+1));
               }
               if ( Grid.BCtypeB[i][j] == BC_NO_SLIP) {
                  Wl = LES3DFsd_pState::NoSlip(Wr,WoB[i][j],
                                               Grid.nfaceBot(i,j,k+1), 
                                               IPs.Pressure_Gradient,
                                               FIXED_TEMPERATURE_WALL);
               }
               if ( Grid.BCtypeB[i][j] == BC_MOVING_WALL) {
                  Wl = LES3DFsd_pState::MovingWall(Wr,WoB[i][j],
                                                   Grid.nfaceBot(i,j,k+1),
                                                   IPs.Moving_Wall_Velocity,
                                                   IPs.Pressure_Gradient,
                                                   FIXED_TEMPERATURE_WALL);
               }
            } else if ( k == KCu && 
                      ( Grid.BCtypeT[i][j] == BC_REFLECTION ||
                        Grid.BCtypeT[i][j] == BC_NO_SLIP ||
                        Grid.BCtypeT[i][j] == BC_MOVING_WALL)) {
                
                dX = Grid.xfaceTop(i,j,k)- Grid.Cell[i][j][k].Xc;
                Wl = W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x+
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y+
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;

                if ( Grid.BCtypeT[i][j] == BC_REFLECTION) {
                   Wr = LES3DFsd_pState::Reflect(Wl,Grid.nfaceTop(i,j,k));
                }
                if ( Grid.BCtypeT[i][j] == BC_NO_SLIP) {
                   Wr = LES3DFsd_pState::NoSlip(Wl,WoT[i][j],
                                                Grid.nfaceTop(i,j,k), 
                                                IPs.Pressure_Gradient,
                                                FIXED_TEMPERATURE_WALL );
                }
                if ( Grid.BCtypeT[i][j] == BC_MOVING_WALL) {
                   Wr = LES3DFsd_pState::MovingWall(Wl,WoT[i][j],
                                                    Grid.nfaceTop(i,j,k),
                                                    IPs.Moving_Wall_Velocity,
                                                    IPs.Pressure_Gradient, 
                                                    FIXED_TEMPERATURE_WALL );
                }
             } else {
                dX = Grid.xfaceTop(i,j,k)- Grid.Cell[i][j][k].Xc;
                Wl = W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                dX = Grid.xfaceBot(i,j,k+1)- Grid.Cell[i][j][k+1].Xc;
                Wr = W[i][j][k+1] +
                   ( phi[i][j][k+1]^ dWdx[i][j][k+1])*dX.x +
                   ( phi[i][j][k+1]^ dWdy[i][j][k+1])*dX.y +
                   ( phi[i][j][k+1]^ dWdz[i][j][k+1])*dX.z;
             } /* endif */
             
	  // Spacing for Preconditioner 
	  if(IPs.Flow_Type != FLOWTYPE_INVISCID && IPs.Preconditioning){
 	     delta_n = min(min(TWO*(Grid.volume(i,j,k)/(Grid.AfaceE(i,j,k)+Grid.AfaceW(i,j,k))),
			       TWO*(Grid.volume(i,j,k)/(Grid.AfaceN(i,j,k)+Grid.AfaceS(i,j,k)))),
                               TWO*(Grid.volume(i,j,k)/(Grid.AfaceTop(i,j,k)+Grid.AfaceBot(i,j,k))));
	  }
             switch(IPs.i_Flux_Function) {
             case FLUX_FUNCTION_HLLE :
                Flux =  LES3DFsd_pState::FluxHLLE_n(Wl, Wr, Grid.nfaceTop(i,j,k));
                break;
             case FLUX_FUNCTION_ROE :
                Flux =  LES3DFsd_pState::FluxRoe_n(Wl, Wr, Grid.nfaceTop(i,j,k));
                break;
             case FLUX_FUNCTION_AUSM_PLUS_UP :
		Flux =  LES3DFsd_pState::FluxAUSMplus_up_n(Wl, Wr, Grid.nfaceTop(i,j,k));
	        break;
             } /* endswitch */
             // add viscous flux in k direction
             Flux -=  LES3DFsd_pState::FluxViscous_n(Wl, Wr,
                                                     W[i][j][k],W[i][j][k+1], 
                                                     dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
                                                     dWdx[i][j][k+1],dWdy[i][j][k+1],dWdz[i][j][k+1],
                                                     Grid.nfaceTop(i,j,k),Grid.Vot(i,j,k), 
                                                     Grid.delta_ot(i,j,k),Grid.volume(i,j,k), 
 						     Grid.volume(i,j,k+1),Flow_Type);
             /* Evaluate cell-averaged solution changes. */
             dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*Flux*Grid.AfaceTop(i,j,k)/Grid.volume(i,j,k);
             dUdt[i][j][k+1][k_residual] += (IPs.CFL_Number*dt[i][j][k+1])*Flux*Grid.AfaceBot(i,j,k+1)/Grid.volume(i,j,k+1);
             /* Save top and bottom face boundary flux. */
             // if (j ==  JCl-1) {
//               FluxS[i] = -Flux* Grid.lfaceS(i, j+1);
//           } else if (j ==  JCu) {
//               FluxN[i] = Flux* Grid.lfaceN(i, j);
//           } /* endif */
             
          } /* endfor */
          dUdt[i][j][KCl-1][k_residual] = U_VACUUM;
          dUdt[i][j][KCu+1][k_residual] = U_VACUUM;
      } /* endfor */
   } /* endfor */

   return (0);
}

template<>
int Hexa_Block<LES3DFsd_pState,LES3DFsd_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<LES3DFsd_pState,LES3DFsd_cState> &IPs){
   
   int i, j, k,  k_residual;
   double omega;
   
   // Memory for linear system solver.
   //   LES3DFsd_cState dU_precon;

  /* Additional variables for dual time stepping. */
  double dTime = ZERO;          // Physical time step
  double residual_denominator;  // Improves convergence for inner iterations
  double length_n = max(max(IPs.Grid_IP.Box_Length, IPs.Grid_IP.Box_Width), IPs.Grid_IP.Box_Height);    

//   if (IPs.Dual_Time_Stepping) {
//     dTime = IPs.dTime;
//   }

   /* Evaluate the time step fraction and residual storage 
      location for the stage. */
  
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
   /* Update solution variables for this stage. */
   for ( k  =  KCl ; k <=  KCu ; ++k ) {
      for ( j  =  JCl ; j <=  JCu ; ++j ) {
         for ( i =  ICl ; i <=  ICu ; ++i ) {
            
            if (IPs.Local_Time_Stepping == 
                GLOBAL_TIME_STEPPING || 
                IPs.Local_Time_Stepping == 
                SCALAR_LOCAL_TIME_STEPPING) {
               U[i][j][k] = Uo[i][j][k] + omega* dUdt[i][j][k][k_residual];
	       if (IPs.i_Flow_Type==FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY){
		 U[i][j][k].rhok=0.0;
	       }
	       U[i][j][k].premixed_mfrac();
            }

            // Check physical validity of update solution state
            if (IPs.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
	      if (!U[i][j][k].Realizable_Solution_Check()) {
                cout << "\n " << CFFC_Name() 
                     << " ERROR: Negative Density, Progress Variable, FSD, Kinetic Energy, and/or Sensible Energy: \n"
                     << " cell = (" << i << ", " << j <<", "<< k << ") " 
                     << " X = " <<  Grid.Cell[i][j][k].Xc 
                     << "\n U = " <<  U[i][j][k] 
                     << "\n dUdt = " << dUdt[i][j][k][k_residual] << "\n";
		return (1);
              } /* endif */

            } else {
	      if (!U[i][j][k].Realizable_Solution_Check()) {
                cout << "\n " << CFFC_Name() 
                     << " ERROR: Negative Density, Progress Variable, FSD, Kinetic Energy, and/or Sensible Energy: \n"
                     << " cell = (" << i << ", " << j <<", "<< k << ") " 
                     << " X = " <<  Grid.Cell[i][j][k].Xc 
                     << "\n U = " <<  U[i][j][k] 
                     << "\n dUdt = " << dUdt[i][j][k][k_residual] << "\n";
		return (1);
              } /* endif */
	    } /* endif */

        /************ FORM LHS FOR DUAL TIME STEPPING SIMI-IMPLICIT WITH PRECONDITIONING ***************/
//         if (Input_Parameters.Local_Time_Stepping == DUAL_SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
//           double dual_coef_I;
// 	  dual_coef_I = (Input_Parameters.first_step  ?  ONE:THREE/TWO);
          
//           // dSdU
//           dSdU = SemiImplicitBlockJacobi(SolnBlk,i, j,Input_Parameters.Simple_Chemistry);
//           // I
// 	  LinSys.A.identity();
//           // P
//           SolnBlk.Uo[i][j].Low_Mach_Number_Preconditioner(Precon,
// 							  SolnBlk.Flow_Type,
// 							  delta_n,
//                                                           length_n,
//                                                           dTime,
//                                                           Input_Parameters.Simple_Chemistry);
   
//           // SolnBlk.Uo[i][j].Low_Mach_Number_Preconditioner_Inverse(Precon_Inv,
// // 	         						 SolnBlk.Flow_Type,
// // 								 delta_n,
// //                                                                  Input_Parameters.Simple_Chemistry);

//           //                             dt
//           //  omega*CFL* dual_coef_I * ------- * I
//           //                            dTime
//           LinSys.A = dual_coef_I*(/*omega*/Input_Parameters.CFL_Number*SolnBlk.dt[i][j]/dTime)*LinSys.A;
//           //                            dt
//           //  omega*CFL*dual_coef_I *  ----- * I - omega*CFL*dt*dSdU 
//           //                           dTime
             
//           LinSys.A -= (/*omega*/Input_Parameters.CFL_Number*SolnBlk.dt[i][j]*dSdU);
    
//           //                              dt
//           //  P + omega*CFL*dual_coef_I  ----- * I - omega*CFL*dt*dSdU 
//           //                             dTime
//           LinSys.A += Precon;   
// 	}
            W[i][j][k] = U[i][j][k].W();
         } /* endfor */    	 
      } /* endfor */    
   } /* endfor */

   /* Solution successfully updated. */
    
   return (0);   
}
 

