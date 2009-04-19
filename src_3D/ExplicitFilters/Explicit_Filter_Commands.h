/*
 *  Explicit_Filter_Commands.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/03/08.
 *
 */

#ifndef _EXPLICIT_FILTER_COMMANDS_INCLUDED
#define _EXPLICIT_FILTER_COMMANDS_INCLUDED

#ifndef _HEXA_SOLVER_CLASSES_INCLUDED
#include  "../HexaBlock/HexaSolverClasses.h"
#endif  //_HEXA_SOLVER_CLASSES_INCLUDED

#ifndef _EXPLICIT_FILTER_CONSTANTS_INCLUDED
#include "Explicit_Filter_Constants.h"
#endif

/* ------------------------------------------------------------------------------------------------------------------------------ */
/**
 * CLASS : Explicit_Filter_Commands
 */

class Explicit_Filter_Commands {
public:
    
    template <typename Soln_pState, typename Soln_cState>
    static int Initialize_Filters(HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data);

    template <typename Soln_pState, typename Soln_cState>
    static int Filter_Solution(HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data, Explicit_Filter_Constants::Filter_Number filter);

    template <typename Soln_pState, typename Soln_cState>
    static int Filter_Residual(int i_stage, HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data);

    template <typename Soln_pState, typename Soln_cState>
    static int Transfer_Function(HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data);
    
    template <typename Filter_Variable, typename Soln_pState, typename Soln_cState>
    static int Filter(Filter_Variable filter_variable, HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data);

};

template <typename Soln_pState, typename Soln_cState>
int Explicit_Filter_Commands::Initialize_Filters(HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data)
{
    int error_flag(0);

    if (Solution_Data.Input.Turbulence_IP.i_filter_type != Explicit_Filter_Constants::IMPLICIT_FILTER) {
        
        Hexa_Block<Soln_pState,Soln_cState> *Soln_Blks = Solution_Data.Local_Solution_Blocks.Soln_Blks;
        for (int nBlk = 0; nBlk < Solution_Data.Local_Solution_Blocks.Number_of_Soln_Blks; nBlk++ ) {
            if (Solution_Data.Local_Solution_Blocks.Block_Used[nBlk]) {        
                Soln_Blks[nBlk].Explicit_Filter.Initialize(Data.batch_flag,Solution_Data.Input);
                Soln_Blks[nBlk].Explicit_Secondary_Filter.Initialize_Secondary(Data.batch_flag,Solution_Data.Input);
            }
        }
    }
    return error_flag;
}

template <typename Soln_pState, typename Soln_cState>
int Explicit_Filter_Commands::Filter_Solution(HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data, Explicit_Filter_Constants::Filter_Number filter)
{
    
    int error_flag(0);
    if (filter != Explicit_Filter_Constants::PRIMARY_FILTER && filter != Explicit_Filter_Constants::SECONDARY_FILTER){
        cout << endl << "Invalid filter number passed. " << endl;
    }

        
    if (Solution_Data.Input.Turbulence_IP.i_filter_type != Explicit_Filter_Constants::IMPLICIT_FILTER) {

        
        Soln_cState *** (Hexa_Block<Soln_pState,Soln_cState>::*U_ptr) = &Hexa_Block<Soln_pState,Soln_cState>::U;
        Hexa_Block<Soln_pState,Soln_cState> *Soln_Blks = Solution_Data.Local_Solution_Blocks.Soln_Blks;
        
        for (int nBlk = 0; nBlk < Solution_Data.Local_Solution_Blocks.Number_of_Soln_Blks; nBlk++ ) {
            if (Solution_Data.Local_Solution_Blocks.Block_Used[nBlk]) {
                
                if (filter==Explicit_Filter_Constants::PRIMARY_FILTER) {
                    Soln_Blks[nBlk].Explicit_Filter.filter(U_ptr);   
                } else if (filter==Explicit_Filter_Constants::SECONDARY_FILTER) {
                    Soln_Blks[nBlk].Explicit_Secondary_Filter.filter(U_ptr);
                }
                
                // Calculate primitive variables from conservative
                for (int k  = Soln_Blks[nBlk].KCl-Soln_Blks[nBlk].Nghost ; k <= Soln_Blks[nBlk].KCu+Soln_Blks[nBlk].Nghost ; ++k ) {
                    for ( int j  = Soln_Blks[nBlk].JCl-Soln_Blks[nBlk].Nghost ; j <= Soln_Blks[nBlk].JCu+Soln_Blks[nBlk].Nghost ; ++j ) {
                        for ( int i = Soln_Blks[nBlk].ICl-Soln_Blks[nBlk].Nghost ; i <= Soln_Blks[nBlk].ICu+Soln_Blks[nBlk].Nghost ; ++i ) {
                            Soln_Blks[nBlk].W[i][j][k] =  Soln_Blks[nBlk].U[i][j][k].W();
                        }	  
                    }
                }
            }
        }
        
        //    double Soln_pState::*p_ptr = p_ptr = &Soln_pState::p; 
        //    Explicit_Filter.filter(p_ptr);    
        //    for (int nBlk = 0; nBlk < Number_of_Soln_Blks; ++nBlk) {
        //        if (Block_Used[nBlk]) {
        //            for (int k  = Soln_Blks[nBlk].KCl-Soln_Blks[nBlk].Nghost ; k <= Soln_Blks[nBlk].KCu+Soln_Blks[nBlk].Nghost ; ++k ) {
        //                for ( int j  = Soln_Blks[nBlk].JCl-Soln_Blks[nBlk].Nghost ; j <= Soln_Blks[nBlk].JCu+Soln_Blks[nBlk].Nghost ; ++j ) {
        //                    for ( int i = Soln_Blks[nBlk].ICl-Soln_Blks[nBlk].Nghost ; i <= Soln_Blks[nBlk].ICu+Soln_Blks[nBlk].Nghost ; ++i ) {
        //                         Soln_Blks[nBlk].U[i][j][k] =  Soln_Blks[nBlk].W[i][j][k].U();
        //                    }	  
        //                }
        //            }
        //        }
        //    }            
    }
    return (error_flag);     
}

template <typename Soln_pState, typename Soln_cState>
int Explicit_Filter_Commands::Filter_Residual(int i_stage, HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data)
{
    
    int error_flag(0);
    
    if (Solution_Data.Input.Turbulence_IP.i_filter_type != Explicit_Filter_Constants::IMPLICIT_FILTER) {
        if (Solution_Data.Input.Turbulence_IP.Filter_Method == Explicit_Filter_Constants::FILTER_RESIDUALS) {
            
            int residual_index = 0;
            
            switch(Solution_Data.Input.i_Time_Integration) {
                case TIME_STEPPING_EXPLICIT_EULER :
                    residual_index = 0;
                    break;
                    
                case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
                    residual_index = 0;
                    break;
                    
                case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
                    residual_index = 0;
                    if (Solution_Data.Input.N_Stage == 4) {
                        if (i_stage == 4) {
                            residual_index = 0;
                        } else {
                            residual_index = i_stage - 1;
                        } /* endif */
                    } /* endif */
                    break;
                    
                case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
                    residual_index = 0;
                    break;
                    
                default:
                    residual_index = 0;
                    break;
                    
            } /* endswitch */
            
            
            
            Solution_Data.Local_Solution_Blocks.BCs_dUdt(Solution_Data.Input,residual_index);
            error_flag = Send_Messages_Residual<Hexa_Block<Soln_pState,Soln_cState> >
            (Solution_Data.Local_Solution_Blocks.Soln_Blks,
             Data.Local_Adaptive_Block_List,
             residual_index);
            
            Soln_cState **** (Hexa_Block<Soln_pState,Soln_cState>::*dUdt_ptr) = &Hexa_Block<Soln_pState,Soln_cState>::dUdt;
            Hexa_Block<Soln_pState,Soln_cState> *Soln_Blks = Solution_Data.Local_Solution_Blocks.Soln_Blks;
            
            for (int nBlk = 0; nBlk < Solution_Data.Local_Solution_Blocks.Number_of_Soln_Blks; nBlk++ ) {
                if (Solution_Data.Local_Solution_Blocks.Block_Used[nBlk]) {
                    Soln_Blks[nBlk].Explicit_Filter.filter(dUdt_ptr,residual_index);   
                }
            }
            
        }       
    }
    return (error_flag);     
}


template <typename Soln_pState, typename Soln_cState>
int Explicit_Filter_Commands::Transfer_Function(HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data)
{
    
    int error_flag(0);
    
    if (CFFC_Primary_MPI_Processor()) {
        
        Hexa_Block<Soln_pState,Soln_cState> *Soln_Blks = Solution_Data.Local_Solution_Blocks.Soln_Blks;
        
        bool first_block = true;
        
        for (int nBlk = 0; nBlk < Solution_Data.Local_Solution_Blocks.Number_of_Soln_Blks; nBlk++ ) {
            if (Solution_Data.Local_Solution_Blocks.Block_Used[nBlk]) {
                if (first_block)
                    Soln_Blks[nBlk].Explicit_Filter.transfer_function(Explicit_Filter_Constants::MIDDLE_CELL);   
                first_block = false;
            }
        }
    }

    return (error_flag);     
}

template <typename Filter_Variable, typename Soln_pState, typename Soln_cState>
int Explicit_Filter_Commands::Filter(Filter_Variable filter_variable, HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data)
{
    int error_flag(0);
    
    if (Solution_Data.Input.Turbulence_IP.i_filter_type != Explicit_Filter_Constants::IMPLICIT_FILTER) {
            
            Hexa_Block<Soln_pState,Soln_cState> *Soln_Blks = Solution_Data.Local_Solution_Blocks.Soln_Blks;
            
            for (int nBlk = 0; nBlk < Solution_Data.Local_Solution_Blocks.Number_of_Soln_Blks; nBlk++ ) {
                if (Solution_Data.Local_Solution_Blocks.Block_Used[nBlk]) {
                    
                    Soln_Blks[nBlk].Explicit_Filter.filter(filter_variable);   
                }
            }
    }
    return (error_flag); 
}



#endif
