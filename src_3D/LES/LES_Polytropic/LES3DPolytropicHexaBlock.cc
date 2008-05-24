
#ifndef  _LES3D_POLYTROPIC_HEXA_BLOCK_INCLUDED
#include "LES3DPolytropicHexaBlock.h"
#endif // _LES3D_POLYTROPIC_HEXA_BLOCK_INCLUDED

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution values at the nodes of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template<>
void Hexa_Block<LES3D_Polytropic_pState,LES3D_Polytropic_cState>::
Output_Tecplot(Input_Parameters<LES3D_Polytropic_pState,LES3D_Polytropic_cState> &IPs,
               const int Number_of_Time_Steps,
               const double &Time,
               const int Block_Number,
               const int Output_Title,
               ostream &Out_File) {
    
    LES3D_Polytropic_pState W_node;
    
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
        << "\"T\" \\ \n"
        << "\"M\" \\ \n"
        << "\"k_SFS\" \\ \n"
        << "\"Q_criterion\" \\ \n";
        
        Out_File<< "ZONE T =  \"Block Number = " << Block_Number
        << "\" \\ \n"
        << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
        << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
        << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
        << "DATAPACKING = POINT \n";
        
    } else {
        Out_File << "ZONE T =  \"Block Number = " << Block_Number
        << "\" \\ \n"
        << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
        << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
        << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
        << "DATAPACKING = POINT \n";              
    } /* endif */
    
    for (int k = Grid.KNl ; k <= Grid.KNu ; ++k) {
        for (int j = Grid.JNl ; j <= Grid.JNu ; ++j) {
            for (int i = Grid.INl ; i <= Grid.INu ; ++i) {
                W_node = Wn(i, j, k);
                Out_File << " "  << Grid.Node[i][j][k].X << W_node;
                Out_File.setf(ios::scientific);
                Out_File << " " << W_node.T()
                << " " << W_node.M()
                << " " << SFS_Kinetic_Energy_n(*this, i ,j ,k)
                << " " << Q_criterion_n(*this, i, j, k) << "\n"; 
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
template<>
void Hexa_Block<LES3D_Polytropic_pState,LES3D_Polytropic_cState>::
Output_Cells_Tecplot(Input_Parameters<LES3D_Polytropic_pState, 
                     LES3D_Polytropic_cState> &IPs,
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File) {
    
    /* Ensure boundary conditions are updated before
     evaluating solution at the nodes. */
    
    BCs(IPs);
    
    /* Output cell centred solution data. */
    
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
        << "\"T\" \\ \n"
        << "\"M\" \\ \n"
        << "\"k_SFS\" \\ \n"
        << "\"Q_criterion\" \\ \n";
        
        Out_File << "ZONE T =  \"Block Number = " << Block_Number
        << "\" \\ \n"
        << "I = " << Grid.ICu - Grid.ICl + 2*Grid.Nghost + 1 << " \\ \n"
        << "J = " << Grid.JCu - Grid.JCl + 2*Grid.Nghost + 1 << " \\ \n"
        << "K = " << Grid.KCu - Grid.KCl + 2*Grid.Nghost + 1 << " \\ \n"
        << "DATAPACKING = POINT \n";
    } else {
        Out_File << "ZONE T =  \"Block Number = " << Block_Number
        << "\" \\ \n"
        << "I = " << Grid.ICu - Grid.ICl + 2*Grid.Nghost + 1 << " \\ \n"
        << "J = " << Grid.JCu - Grid.JCl + 2*Grid.Nghost + 1 << " \\ \n"
        << "K = " << Grid.KCu - Grid.KCl + 2*Grid.Nghost + 1 << " \\ \n"
        << "DATAPACKING = POINT \n";
        
    } /* endif */
    
    for (int k = Grid.KCl-Grid.Nghost; k <= Grid.KCu+Grid.Nghost; ++k) {
        for (int j  = Grid.JCl-Grid.Nghost; j <= Grid.JCu+Grid.Nghost; ++j ) {
            for (int i = Grid.ICl-Grid.Nghost; i <= Grid.ICu+Grid.Nghost; ++i ) {
                
                Out_File << " "  <<  Grid.Cell[i][j][k].Xc
                <<  W[i][j][k];
                Out_File.setf(ios::scientific);
                Out_File << " " << W[i][j][k].T() 
                << " " << W[i][j][k].M()
                << " " << W[i][j][k].SFS_Kinetic_Energy(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],Grid.Cell[i][j][k].V)
                << " " << W[i][j][k].Q_criterion(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k]) <<"\n ";
                Out_File.unsetf(ios::scientific);
            } /* endfor */
        } /* endfor */
    } /* endfor */
    
    Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the solution values at the nodes of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template<>
void Hexa_Block<LES3D_Polytropic_pState,LES3D_Polytropic_cState>::
Output_Nodes_Tecplot(Input_Parameters<LES3D_Polytropic_pState,LES3D_Polytropic_cState> &IPs,
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File) {
    
    
    LES3D_Polytropic_pState W_node;
    
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
        << "\"T\" \\ \n"
        << "\"M\" \\ \n"
        << "\"k_SFS\" \\ \n"
        << "\"Q_criterion\" \\ \n";
        
        Out_File<< "ZONE T =  \"Block Number = " << Block_Number
        << "\" \\ \n"
        << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
        << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
        << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
        << "DATAPACKING = POINT \n";
        
    } else {
        Out_File << "ZONE T =  \"Block Number = " << Block_Number
        << "\" \\ \n"
        << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
        << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
        << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
        << "DATAPACKING = POINT \n";              
    } /* endif */
    
    for (int k = Grid.KNl ; k <= Grid.KNu ; ++k) {
        for (int j = Grid.JNl ; j <= Grid.JNu ; ++j) {
            for (int i = Grid.INl ; i <= Grid.INu ; ++i) {
                W_node = Wn(i, j, k);
                Out_File << " "  << Grid.Node[i][j][k].X << W_node;
                Out_File.setf(ios::scientific);
                Out_File << " " << W_node.T()
                << " " << W_node.M()
                << " " << SFS_Kinetic_Energy_n(*this, i ,j ,k)
                << " " << Q_criterion_n(*this, i, j, k) << "\n"; 
                Out_File.unsetf(ios::scientific);
            } /* endfor */
        } /* endfor */
    } /* endfor */
    
    Out_File << setprecision(6);
}




/********************************************************
 * Routine: ICs_Specializations                         *
 *                                                      *
 * Apply initial conditions for the k-equation model    *
 *                                                      *
 ********************************************************/
template<>
int Hexa_Block<LES3D_Polytropic_pState,LES3D_Polytropic_cState>::
ICs_Specializations(Input_Parameters<LES3D_Polytropic_pState,LES3D_Polytropic_cState> &IPs){
    
//    if (false) { //(IPs.Turbulence_IP.i_filter_type != FILTER_TYPE_IMPLICIT)*/ {
//        if (CFFC_Primary_MPI_Processor()) {
//            cout << endl;
//            cout << " ------------------------------------------------" << endl;
//            cout << "    Explicitly filtering the initial condition   " << endl;
//            cout << " ------------------------------------------------" << endl;        
//        }
//        LES3D_Polytropic_cState *** (Hexa_Block<LES3D_Polytropic_pState,LES3D_Polytropic_cState>::*U_ptr) = &Hexa_Block<LES3D_Polytropic_pState,LES3D_Polytropic_cState>::U;
//        double (LES3D_Polytropic_pState::*p_ptr) = p_ptr = &LES3D_Polytropic_pState::p; 
//        LES_Filter<LES3D_Polytropic_pState,LES3D_Polytropic_cState> Explicit_Filter(*this,IPs);
//        Explicit_Filter.filter(U_ptr);
        for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for ( int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
                for ( int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                    W[i][j][k] = U[i][j][k].W();
                }	  
            }
        }
//        if (CFFC_Primary_MPI_Processor()) {
//            Explicit_Filter.transfer_function();
//        }
//        Explicit_Filter.filter(p_ptr);
//        for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
//            for ( int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
//                for ( int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
//                    U[i][j][k] = W[i][j][k].U();
//                }	  
//            }
//        }
//        if (CFFC_Primary_MPI_Processor()) {
//            cout << "    Finished explicit filtering " << endl;
//        }
//    }
    
    Linear_Reconstruction_LeastSquares(IPs.i_Limiter);
    for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
        for ( int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
            for ( int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                U[i][j][k] = W[i][j][k].U();
            } /* endfor */ 	  
        } /* endfor */
    } /* endfor */
    
    /* Set default values for the boundary conditions
     reference states. */
    
    for (int k = KCl-Nghost ; k<= KCu+Nghost; ++k) {
        for (int j = JCl-Nghost ; j<= JCu+Nghost; ++j){
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
    } /* endfor */
    
    for (int  k = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
        for (int  i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
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
            } /* endif */
        } /* endfor */
    } /* endfor */
    
    for (int  j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
        for (int  i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
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
            } /* endif */
        } /* endfor */
    } /* endfor */
    return (0);
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
double Hexa_Block<LES3D_Polytropic_pState,LES3D_Polytropic_cState>::
CFL(Input_Parameters<LES3D_Polytropic_pState,LES3D_Polytropic_cState> &IPs){
    
    double dt_acoustic(MILLION), dt_viscous(MILLION);
    
    double dtMin, d_i, d_j, d_k, v_i, v_j, v_k, a, dt_vis, nu;
    double mr, aa_i, aa_j, aa_k;
    double length_n, delta_n, dTime;
    
    dtMin = MILLION;
    
    for (int k = KCl- Nghost; k <= KCu + Nghost; ++k) {
        for (int j =  JCl- Nghost; j <= JCu + Nghost; ++j) {
            for (int i =  ICl- Nghost; i <= ICu+ Nghost; ++i) {
                if (i <  ICl || i >  ICu ||
                    j <  JCl || j >  JCu || 
                    k <  KCl || k >  KCu) {
                    dt[i][j][k] = ZERO;
                } else {
                    d_i = TWO*(Grid.Cell[i][j][k].V/
                               (Grid.AfaceE(i, j, k)+ Grid.AfaceW(i, j, k)));
                    d_j = TWO*( Grid.Cell[i][j][k].V/
                               (Grid.AfaceN(i, j, k)+ Grid.AfaceS(i, j, k)));
                    d_k = TWO*( Grid.Cell[i][j][k].V/
                               (Grid.AfaceTop(i, j, k)+ Grid.AfaceBot(i, j, k)));
                    v_i = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                                (Grid.nfaceE(i, j, k)- Grid.nfaceW(i, j, k)));
                    v_j = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                                ( Grid.nfaceN(i, j, k)- Grid.nfaceS(i, j, k)));
                    v_k = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                                (Grid.nfaceTop(i, j, k)- Grid.nfaceBot(i, j, k)));
                    
                    
                    length_n = max(max(IPs.Grid_IP.Box_Length, 
                                       IPs.Grid_IP.Box_Width), 
                                   IPs.Grid_IP.Box_Height);  
                    delta_n = min(min(fabs(d_i),fabs(d_j)),fabs(d_k));
                    
                    /* ---------- Acoustic limitation ----------- */
                    if(IPs.Preconditioning == 0){
                        a =  W[i][j][k].a();
                        dt[i][j][k] = min(min(d_i/(a+fabs(v_i)), d_j/(a+fabs(v_j))),
                                          d_k/(a+fabs(v_k)));
                        
                    }  /* endif */
                    
                    /* ---------- Viscous limitation ------------ */
                    nu  = W[i][j][k].nu();
                    nu += W[i][j][k].nu_t(dWdx[i][j][k], dWdy[i][j][k], dWdz[i][j][k], Grid.Cell[i][j][k].V);
                    dt_vis = min(min((d_i*d_i)/nu, (d_j*d_j)/nu), (d_k*d_k)/nu)/THREE; 
                    dt[i][j][k]  = min(dt_vis, dt[i][j][k]);                        
                    dtMin = min(dtMin,  dt[i][j][k]);
                    
                } /* endif */
                
            } /* endfor */
        } /* endfor */
    } /* endfor */
    
    
    for (int k  =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k ) {
        for (int j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) {
            for (int i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
                if (i <  ICl || i >  ICu ||
                    j <  JCl || j >  JCu || 
                    k <  KCl || k >  KCu) {
                    dt[i][j][k] = dtMin;
                } /* endif */
            } /* endfor */
        } /* endfor */
    } /* endfor */
    
    /* Return the global time step. */
    
    return (dtMin);
    
}

/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine determines the solution residuals for a *
 * given stage of a variety of multi-stage explicit     *
 * time integration schemes for a given solution block. *
 *                                                      *
 ********************************************************/
template<>
int Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState>::
dUdt_Multistage_Explicit(const int i_stage,
                         Input_Parameters<LES3D_Polytropic_pState, LES3D_Polytropic_cState> &IPs) {
    
    int i, j, k,  k_residual;
    double omega;
    Vector3D dX;
    
    LES3D_Polytropic_pState Wl, Wr;
    LES3D_Polytropic_cState Flux, Temp;
    
    LES3D_Polytropic_cState U_VACUUM;
    U_VACUUM.Vacuum();
    LES3D_Polytropic_pState W_VACUUM;
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
                } 
            }
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
                                dUdt[i+1][j][k][k_residual] = dUdt[i+1][j][k][0] + 
                                TWO*dUdt[i+1][j][k][1] +
                                TWO*dUdt[i+1][j][k][2];
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
                            Wl = LES3D_Polytropic_pState::Reflect(Wr,Grid.nfaceW(i+1,j,k));
                        } /* endif */
                        
                        if ( Grid.BCtypeW[j][k] == BC_NO_SLIP) {
                            Wl = LES3D_Polytropic_pState::NoSlip(Wr,WoW[j][k],
                                                                 Grid.nfaceW(i+1,j,k), 
                                                                 IPs.Pressure_Gradient,
                                                                 FIXED_TEMPERATURE_WALL);
                            
                        } /* endif */
                        
                        if ( Grid.BCtypeW[j][k] == BC_MOVING_WALL) {
                            Wl = LES3D_Polytropic_pState::MovingWall(Wr,WoW[j][k],
                                                                     Grid.nfaceW(i+1,j,k),
                                                                     IPs.Moving_Wall_Velocity,
                                                                     IPs.Pressure_Gradient,
                                                                     FIXED_TEMPERATURE_WALL);
                        } /* endif */
                        
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
                            Wr = LES3D_Polytropic_pState::Reflect(Wl,Grid.nfaceE(i,j,k));
                        } /* endif */
                        
                        if ( Grid.BCtypeE[j][k] == BC_NO_SLIP) {
                            Wr = LES3D_Polytropic_pState::NoSlip(Wl,WoE[j][k],
                                                                 Grid.nfaceE(i,j,k), 
                                                                 IPs.Pressure_Gradient, 
                                                                 FIXED_TEMPERATURE_WALL);
                        } /* endif */
                        
                        if ( Grid.BCtypeE[j][k] == BC_MOVING_WALL) {
                            Wr = LES3D_Polytropic_pState::MovingWall(Wl,WoE[j][k],
                                                                     Grid.nfaceE(i,j,k),  
                                                                     IPs.Moving_Wall_Velocity,
                                                                     IPs.Pressure_Gradient,
                                                                     FIXED_TEMPERATURE_WALL);
                        } /* endif */
                        
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
                        delta_n = min(min(TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceE(i,j,k)+Grid.AfaceW(i,j,k))),
                                          TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceN(i,j,k)+Grid.AfaceS(i,j,k)))),
                                      TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceTop(i,j,k)+Grid.AfaceBot(i,j,k))));
                    } /* endif */
                    
                    // Evaluate inviscid numerical flux
                    switch(IPs.i_Flux_Function) {
                        case FLUX_FUNCTION_HLLE :
                            Flux = LES3D_Polytropic_pState::FluxHLLE_n(Wl, Wr, Grid.nfaceE(i,j,k));
                            break;
                        case FLUX_FUNCTION_ROE :
                            Flux = LES3D_Polytropic_pState::FluxRoe_n(Wl, Wr, Grid.nfaceE(i,j,k));
                            break;
                        case FLUX_FUNCTION_AUSM_PLUS_UP :
                            Flux = LES3D_Polytropic_pState::FluxAUSMplus_up_n(Wl, Wr, Grid.nfaceE(i,j,k));
                            break;
                    } /* endswitch */
                    
                    // Add viscous flux in i direction
                    Flux -=  LES3D_Polytropic_pState::FluxViscous_n(Wl, Wr,
                                                                    W[i][j][k], W[i+1][j][k], 
                                                                    dWdx[i][j][k], dWdy[i][j][k], dWdz[i][j][k],
                                                                    dWdx[i+1][j][k], dWdy[i+1][j][k], dWdz[i+1][j][k],
                                                                    Grid.nfaceE(i,j,k), Grid.Voe(i,j,k),  
                                                                    Grid.delta_oe(i,j,k), Grid.Cell[i][j][k].V, 
                                                                    Grid.Cell[i+1][j][k].V);
                    
                    /* Evaluate cell-averaged solution changes. */
                    dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*
                    Flux*Grid.AfaceE(i,j,k)/Grid.Cell[i][j][k].V;
                    dUdt[i+1][j][k][k_residual] += (IPs.CFL_Number*dt[i+1][j][k])*
                    Flux*Grid.AfaceW(i+1,j,k)/Grid.Cell[i+1][j][k].V;
                    
                    /* Include physical time derivative for dual time stepping. */
                    //                if (IPs.Dual_Time_Stepping) {
                    //                   dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*
                    //                                                W[i][j][k].S_dual_time_stepping(U[i][j][k], Ut[i][j][k], Uold[i][j][k], 
                    //                                                                                dTime, IPs.first_step);
                    //                } /* endif */
                    
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
                        Wl = LES3D_Polytropic_pState::Reflect(Wr,Grid.nfaceS(i,j+1,k));
                    } /* endif */
                    
                    if ( Grid.BCtypeS[i][k] == BC_NO_SLIP) {
                        Wl = LES3D_Polytropic_pState::NoSlip(Wr,WoS[i][k],
                                                             Grid.nfaceS(i,j+1,k),
                                                             IPs.Pressure_Gradient,
                                                             FIXED_TEMPERATURE_WALL);
                    } /* endif */
                    
                    if ( Grid.BCtypeS[i][k] == BC_MOVING_WALL) {
                        Wl = LES3D_Polytropic_pState::MovingWall(Wr,WoS[i][k],
                                                                 Grid.nfaceS(i,j+1,k),
                                                                 IPs.Moving_Wall_Velocity,
                                                                 IPs.Pressure_Gradient,
                                                                 FIXED_TEMPERATURE_WALL);
                    } /* endif */
                    
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
                        Wr = LES3D_Polytropic_pState::Reflect(Wl,Grid.nfaceN(i,j,k));
                    } /* endif */
                    
                    if ( Grid.BCtypeN[i][k] == BC_NO_SLIP) {
                        Wr = LES3D_Polytropic_pState::NoSlip(Wl,WoN[i][k],
                                                             Grid.nfaceN(i,j,k),
                                                             IPs.Pressure_Gradient,
                                                             FIXED_TEMPERATURE_WALL);
                    } /* endif */
                    
                    if ( Grid.BCtypeN[i][k] == BC_MOVING_WALL) {
                        Wr = LES3D_Polytropic_pState::MovingWall(Wl,WoN[i][k],
                                                                 Grid.nfaceN(i,j,k),
                                                                 IPs.Moving_Wall_Velocity,
                                                                 IPs.Pressure_Gradient,
                                                                 FIXED_TEMPERATURE_WALL);
                    } /* endif */
                    
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
                if (IPs.Flow_Type != FLOWTYPE_INVISCID && IPs.Preconditioning){
                    delta_n = min(min(TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceE(i,j,k)+Grid.AfaceW(i,j,k))),
                                      TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceN(i,j,k)+Grid.AfaceS(i,j,k)))),
                                  TWO*( Grid.Cell[i][j][k].V/(Grid.AfaceTop(i,j,k)+Grid.AfaceBot(i,j,k))));
                } /* endif */
                
                // Evaluate inviscid numerical flux
                switch(IPs.i_Flux_Function) {
                    case FLUX_FUNCTION_HLLE :
                        Flux = LES3D_Polytropic_pState::FluxHLLE_n(Wl, Wr, Grid.nfaceN(i,j,k));
                        break;
                    case FLUX_FUNCTION_ROE :
                        Flux = LES3D_Polytropic_pState::FluxRoe_n(Wl, Wr, Grid.nfaceN(i,j,k));
                        break;
                    case FLUX_FUNCTION_AUSM_PLUS_UP :
                        Flux = LES3D_Polytropic_pState::FluxAUSMplus_up_n(Wl, Wr, Grid.nfaceN(i,j,k));
                        break;
                } /* endswitch */
                
                // Add viscous flux in j direction
                Flux -= LES3D_Polytropic_pState::FluxViscous_n(Wl, Wr, 
                                                               W[i][j][k], W[i][j+1][k],
                                                               dWdx[i][j][k], dWdy[i][j][k], dWdz[i][j][k],
                                                               dWdx[i][j+1][k], dWdy[i][j+1][k], dWdz[i][j+1][k],
                                                               Grid.nfaceN(i,j,k), Grid.Von(i,j,k),
                                                               Grid.delta_on(i,j,k), Grid.Cell[i][j][k].V, 
                                                               Grid.Cell[i][j+1][k].V);
                
                /* Evaluate cell-averaged solution changes. */
                dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*
                Flux*Grid.AfaceN(i,j,k)/Grid.Cell[i][j][k].V;
                dUdt[i][j+1][k][k_residual] += (IPs.CFL_Number*dt[i][j+1][k])*
                Flux*Grid.AfaceS(i,j+1,k)/Grid.Cell[i][j+1][k].V;
                
                /* Save south and north face boundary flux. */
                
                //           if (j ==  JCl-1) {
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
                        Wl = LES3D_Polytropic_pState::Reflect(Wr,Grid.nfaceBot(i,j,k+1));
                    } /* endif */
                    
                    if ( Grid.BCtypeB[i][j] == BC_NO_SLIP) {
                        Wl = LES3D_Polytropic_pState::NoSlip(Wr,WoB[i][j],
                                                             Grid.nfaceBot(i,j,k+1), 
                                                             IPs.Pressure_Gradient,
                                                             FIXED_TEMPERATURE_WALL);
                    } /* endif */
                    
                    if ( Grid.BCtypeB[i][j] == BC_MOVING_WALL) {
                        Wl = LES3D_Polytropic_pState::MovingWall(Wr,WoB[i][j],
                                                                 Grid.nfaceBot(i,j,k+1),
                                                                 IPs.Moving_Wall_Velocity,
                                                                 IPs.Pressure_Gradient,
                                                                 FIXED_TEMPERATURE_WALL);
                    } /* endif */
                    
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
                        Wr = LES3D_Polytropic_pState::Reflect(Wl,Grid.nfaceTop(i,j,k));
                    } /* endif */
                    
                    if ( Grid.BCtypeT[i][j] == BC_NO_SLIP) {
                        Wr = LES3D_Polytropic_pState::NoSlip(Wl,WoT[i][j],
                                                             Grid.nfaceTop(i,j,k), 
                                                             IPs.Pressure_Gradient,
                                                             FIXED_TEMPERATURE_WALL );
                    } /* endif */
                    
                    if ( Grid.BCtypeT[i][j] == BC_MOVING_WALL) {
                        Wr = LES3D_Polytropic_pState::MovingWall(Wl,WoT[i][j],
                                                                 Grid.nfaceTop(i,j,k),
                                                                 IPs.Moving_Wall_Velocity,
                                                                 IPs.Pressure_Gradient, 
                                                                 FIXED_TEMPERATURE_WALL );
                    } /* endif */
                    
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
                if (IPs.Flow_Type != FLOWTYPE_INVISCID && IPs.Preconditioning) {
                    delta_n = min(min(TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceE(i,j,k)+Grid.AfaceW(i,j,k))),
                                      TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceN(i,j,k)+Grid.AfaceS(i,j,k)))),
                                  TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceTop(i,j,k)+Grid.AfaceBot(i,j,k))));
                } /* endif */
                
                // Evaluate inviscid numerical flux
                switch(IPs.i_Flux_Function) {
                    case FLUX_FUNCTION_HLLE :
                        Flux = LES3D_Polytropic_pState::FluxHLLE_n(Wl, Wr, Grid.nfaceTop(i,j,k));
                        break;
                    case FLUX_FUNCTION_ROE :
                        Flux = LES3D_Polytropic_pState::FluxRoe_n(Wl, Wr, Grid.nfaceTop(i,j,k));
                        break;
                    case FLUX_FUNCTION_AUSM_PLUS_UP :
                        Flux = LES3D_Polytropic_pState::FluxAUSMplus_up_n(Wl, Wr, Grid.nfaceTop(i,j,k));
                        break;
                } /* endswitch */
                
                // Add viscous flux in k direction
                Flux -=  LES3D_Polytropic_pState::FluxViscous_n(Wl, Wr,
                                                                W[i][j][k], W[i][j][k+1], 
                                                                dWdx[i][j][k], dWdy[i][j][k], dWdz[i][j][k],
                                                                dWdx[i][j][k+1], dWdy[i][j][k+1], dWdz[i][j][k+1],
                                                                Grid.nfaceTop(i,j,k), Grid.Vot(i,j,k), 
                                                                Grid.delta_ot(i,j,k), Grid.Cell[i][j][k].V, 
                                                                Grid.Cell[i][j][k+1].V);
                
                /* Evaluate cell-averaged solution changes. */
                dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*
                Flux*Grid.AfaceTop(i,j,k)/Grid.Cell[i][j][k].V;
                dUdt[i][j][k+1][k_residual] += (IPs.CFL_Number*dt[i][j][k+1])*
                Flux*Grid.AfaceBot(i,j,k+1)/Grid.Cell[i][j][k+1].V;
                
                /* Save top and bottom face boundary flux. */
                //           if (j ==  JCl-1) {
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

/********************************************************
 * Routine: Update_Solution_Multistage_Explicit         *
 *                                                      *
 * This routine updates solution states of the given    *
 * solution block for a variety of multi-stage explicit *
 * time integration schemes.                            *
 *                                                      *
 ********************************************************/
template<>
int Hexa_Block<LES3D_Polytropic_pState,LES3D_Polytropic_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<LES3D_Polytropic_pState,LES3D_Polytropic_cState> &IPs){
    
    int i, j, k,  k_residual;
    double omega;
    
    // Memory for linear system solver.
    //   LES3D_Polytropic_cState dU_precon;
    
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
                
                // Update conserved solution state
                if (IPs.Local_Time_Stepping == GLOBAL_TIME_STEPPING || 
                    IPs.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
                    U[i][j][k] = Uo[i][j][k] + omega* dUdt[i][j][k][k_residual];
                } /* endif */
                
                // Check physical validity of update solution state
                if(U[i][j][k].Unphysical_Properties()) {
                    cout << "\n " << CFFC_Name() 
                    << " ERROR: Negative Density, Total Energy, and/or Specific Internal Energy: \n"
                    << " cell = (" << i << ", " << j <<", "<< k << ") " 
                    << " X = " <<  Grid.Cell[i][j][k].Xc 
                    << "\n U = " <<  U[i][j][k] 
                    << "\n W = " <<  U[i][j][k].W()
                    << "\n dUdt = " << dUdt[i][j][k][k_residual] << "\n";
                    return (1);
                }
                
                // Update primitive solution state
                W[i][j][k] = U[i][j][k].W();
                
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
            } /* endfor */    	 
        } /* endfor */    
    } /* endfor */
    
    /* Solution successfully updated. */
    
    return (0);   
    
}
