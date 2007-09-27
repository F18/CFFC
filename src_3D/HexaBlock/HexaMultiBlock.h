/* HexaMultiBlock.h: Header file for list of multiple solution block. */

#ifndef _HEXA_MULTIBLOCK_INCLUDED
#define _HEXA_MULTIBLOCK_INCLUDED

/* Include various CFFC header files. */

#ifndef _HEXA_BLOCK_INCLUDED
#include "HexaBlock.h"
#endif // _HEXA_BLOCK_INCLUDED

#ifndef _GRID3D_HEXA_MULTIBLOCK_INCLUDED
#include "../Grid/Grid3DHexaMultiBlock.h"
#endif // _GRID3D_HEXA_MULTIBLOCK_INCLUDED

#ifndef _ADAPTIVEBLOCK3D_INCLUDED
#include "../AMR/AdaptiveBlock3D.h"
#endif //_ADAPTIVEBLOCK3D_INCLUDED

// a list of solution blocks on a processor
template<class HEXA_BLOCK> class Hexa_Multi_Block{
   
  private:
  public:
     
   HEXA_BLOCK *Soln_Blks;          // Array of hexahedral solution blocks.
   int Number_of_Soln_Blks;        // Number or size of array of hexahedral solution blocks. 
   int *Block_Used;                // Solution block usage indicator.
   int Allocated;                  // Indicates if the solution blocks have been allocated or not.

   /* Creation constructors. */
   Hexa_Multi_Block(void) {
      Number_of_Soln_Blks = 0; Allocated = 0;
      Soln_Blks = NULL; Block_Used = NULL;
   }

   Hexa_Multi_Block(const int Nblk) {
      Allocate(Nblk);
   }

   Hexa_Multi_Block(Input_Parameters<typename HEXA_BLOCK::Soln_pState,
 		                     typename HEXA_BLOCK::Soln_cState> &IPs) {
      Allocate(IPs.Number_of_Blocks_Per_Processor);
   }

   Hexa_Multi_Block(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                     typename HEXA_BLOCK::Soln_cState> &IPs,
                    Grid3D_Hexa_Multi_Block Multi_Block_Mesh) {
      Allocate(IPs.Number_of_Blocks_Per_Processor);
      int nblk = 0;
      for (int  k = 0; k < Multi_Block_Mesh.NBlk_Kdir; ++k) {
         for (int  j = 0; j < Multi_Block_Mesh.NBlk_Jdir; ++j) {
            for (int i = 0 ; i < Multi_Block_Mesh.NBlk_Jdir; ++i) {
               if (Multi_Block_Mesh.Grid_Blks[i][j][k].Allocated) {
	          Soln_Blks[nblk].Create_Block(Multi_Block_Mesh.Grid_Blks[i][j][k]);
                  Soln_Blks[nblk].Flow_Type = IPs.i_Flow_Type;
                  Block_Used[nblk] = HEXA_BLOCK_USED;
                  ++nblk;
               } /* endif */
            } /* endfor */
         } /* endfor */
      } /* endfor */
   }

   /* Destructor. */
   ~Hexa_Multi_Block() {
      Deallocate();
   }

   /* Define various member functions. */

   void Allocate(const int Nblk);

   void Deallocate(void);

   void Copy(Hexa_Multi_Block &Solution2);

   void Broadcast(void);

   int Read_Restart_Solution(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                              typename HEXA_BLOCK::Soln_cState> &Input,
                             AdaptiveBlock3D_List &Soln_Block_List,
                             int &Number_of_Time_Steps,
                             double &Time,
                             CPUTime &CPU_Time);

   int Write_Restart_Solution(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                               typename HEXA_BLOCK::Soln_cState> &Input,
                              AdaptiveBlock3D_List &Local_Adaptive_Block_List,
                              const int Number_of_Time_Steps,
                              const double &Time,
                              const CPUTime &CPU_Time);

   int Output_Tecplot(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                       typename HEXA_BLOCK::Soln_cState> &Input,
                      AdaptiveBlock3D_List &Local_Adaptive_Block_List,
                      const int Number_of_Time_Steps,
                      const double &Time);

   int Output_Cells_Tecplot(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                             typename HEXA_BLOCK::Soln_cState> &Input,
                            AdaptiveBlock3D_List &Local_Adaptive_Block_List,
                            const int Number_of_Time_Steps,
                            const double &Time);

   int Output_Nodes_Tecplot(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                              typename HEXA_BLOCK::Soln_cState> &Input,
                            AdaptiveBlock3D_List &Local_Adaptive_Block_List,
                            const int Number_of_Time_Steps,
                            const double &Time);

   double L1_Norm_Residual(void);

   double L2_Norm_Residual(void);

   double Max_Norm_Residual(void);

   void Set_Global_TimeStep(const double &Dt_min);

   void ICs(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                             typename HEXA_BLOCK::Soln_cState> &Input);

   void BCs(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                             typename HEXA_BLOCK::Soln_cState> &Input);

   int WtoU(void);

   double CFL(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                               typename HEXA_BLOCK::Soln_cState> &Input);

   int dUdt_Multistage_Explicit(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                typename HEXA_BLOCK::Soln_cState> &Input,
                                const int I_Stage);

   int Update_Solution_Multistage_Explicit(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                                            typename HEXA_BLOCK::Soln_cState> &Input, 
                                           const int I_Stage);

};

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for a list of 3D hexahedral          *
 * solution blocks.                                     *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
void Hexa_Multi_Block<HEXA_BLOCK>::Allocate(const int Nblk) {

   if (Nblk >= 1 && !Allocated) {
      Number_of_Soln_Blks = Nblk;

      Soln_Blks = new HEXA_BLOCK[Number_of_Soln_Blks];
      Block_Used = new int[Number_of_Soln_Blks];
      for (int usage = 0; usage<Number_of_Soln_Blks; ++usage) {
         Block_Used[usage] = HEXA_BLOCK_NOT_USED;
      } /* endfor */

      Allocated = 1;
   } /* endif */

}

/********************************************************
 * Routine: Deallocate                                  *
 *                                                      *
 * Deallocate memory for a list of 3D hexahedral        *
 * solution blocks.                                     *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
void Hexa_Multi_Block<HEXA_BLOCK>::Deallocate(void) {

   if (Number_of_Soln_Blks >= 1 && Allocated) {
      delete []Soln_Blks; delete []Block_Used;
      Soln_Blks = NULL; Block_Used = NULL;

      Number_of_Soln_Blks = 0;

      Allocated = 0;
   } /* endif */

}

/********************************************************
 * Routine: Copy                                        *
 *                                                      *
 * Make a copy of list of hexahedral solution blocks.   *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
void Hexa_Multi_Block<HEXA_BLOCK>::Copy(Hexa_Multi_Block<HEXA_BLOCK> & Solution2) {

  if (Solution2.Allocated) {

    /* Ensure the solution block arrays have same dimensions. */

    if (Allocated && (Number_of_Soln_Blks != Solution2.Number_of_Soln_Blks)) {
      Deallocate();
      Allocate(Solution2.Number_of_Soln_Blks);
    } else if (!Allocated) {
      Allocate(Solution2.Number_of_Soln_Blks);
    } /* endif */

    /* Copy each of the used solution blocks in Solution2. */

    for (int  i = 0 ; i < Number_of_Soln_Blks ; ++i ) {
       if (Solution2.Block_Used[i]) Soln_Blks[i].Copy(Solution2.Soln_Blks[i]);
    } /* endfor */

  } /* endif */

}

/********************************************************
 * Routine: Broadcast                                   *
 *                                                      *
 * Broadcast list of hexahedral solution blocks.        *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
void Hexa_Multi_Block<HEXA_BLOCK>::Broadcast(void) {

  if (Allocated) {

    /* Broadcast each of the hexahedral solution blocks. */

    for (int  i = 0 ; i < Number_of_Soln_Blks ; ++i ) {
      if (Block_Used[i]) Soln_Blks[i].Broadcast();
    } /* endfor */

  } /* endif */

}

/********************************************************
 * Routine: L1_Norm_Residual                            *
 *                                                      *
 * Determines the L1-norm of the solution residual for  *
 * a 1D array of 3D hexahedrial multi-block solution    *
 * blocks.  Useful for monitoring convergence of the    *
 * solution for steady state problems.                  *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
double Hexa_Multi_Block<HEXA_BLOCK>::L1_Norm_Residual(void) {
   
   double l1_norm;
   l1_norm = ZERO;
   
   /* Calculate the L1-norm.
      Sum the L1-norm for each solution block. */
   
   for (int nblk = 0; nblk < Number_of_Soln_Blks; ++nblk) {
      if (Block_Used[nblk]) {
         l1_norm += Soln_Blks[nblk].L1_Norm_Residual();
      } /* endif */
   }  /* endfor */

   /* Return the L1-norm. */

   return (l1_norm);

}

/********************************************************
 * Routine: L2_Norm_Residual                            *
 *                                                      *
 * Determines the L2-norm of the solution residual for  *
 * a 1D array of 3D hexahedrial multi-block solution  *
 * blocks.  Useful for monitoring convergence of the    *
 * solution for steady state problems.                  *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
double Hexa_Multi_Block<HEXA_BLOCK>::L2_Norm_Residual(void) {

   
   double l2_norm;

   l2_norm = ZERO;
   
   /* Sum the square of the L2-norm for each solution block. */
      
   for (int nblk = 0; nblk < Number_of_Soln_Blks; ++nblk) {
      if (Block_Used[nblk]) {
         l2_norm += Soln_Blks[nblk].L2_Norm_Residual();
      } /* endif */
   }  /* endfor */
   
   /* Calculate the L2-norm for all blocks. */
   
   l2_norm = sqrt(l2_norm);
   
   /* Return the L2-norm. */
   
   return (l2_norm);
   
}

/********************************************************
 * Routine: Max_Norm_Residual                           *
 *                                                      *
 * Determines the maximum norm of the solution residual *
 * for a 1D array of 3D hexahedrial multi-block       *
 * solution blocks.  Useful for monitoring convergence  *
 * of the solution for steady state problems.           *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
double Hexa_Multi_Block<HEXA_BLOCK>::Max_Norm_Residual(void) {
   
   double max_norm;
   
   max_norm = ZERO;
   
   /* Find the maximum norm for all solution blocks. */
   
   for (int nblk = 0; nblk < Number_of_Soln_Blks; ++nblk) {
      if(Block_Used[nblk]){
         max_norm = max(max_norm, (Soln_Blks[nblk].Max_Norm_Residual()));
      } /* endif */
   }  /* endfor */
   
   /* Return the maximum norm. */
   
   return (max_norm);
   
}

/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables of a 1D array of 3D hexahedrial *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
void Hexa_Multi_Block<HEXA_BLOCK>::ICs(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                                        typename HEXA_BLOCK::Soln_cState> &Input) {
   
   /* Assign initial data for each solution block. */

   for (int nblk = 0; nblk < Number_of_Soln_Blks; ++nblk) {
      if (Block_Used[nblk]) {
         Soln_Blks[nblk].ICs(Input.i_ICs, Input);
      } /* endif */
   }  /* endfor */
   
}

/********************************************************
 * Routine: BCs                                         *
 *                                                      *
 * Apply boundary conditions at boundaries of a 1D      *
 * array of 3D hexahedral multi-block solution          *
 * blocks.                                              *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
void Hexa_Multi_Block<HEXA_BLOCK>::BCs(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                                        typename HEXA_BLOCK::Soln_cState> &Input) {
   
   /* Prescribe boundary data for each solution block. */

   for (int nblk = 0; nblk < Number_of_Soln_Blks; ++nblk) {
      if (Block_Used[nblk]) {
         Soln_Blks[nblk].BCs(Input);
      } /* endif */
   }  /* endfor */
   
}

/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Euler time stepping scheme) for a 1D   *
 * array of 3D hexahedral multi-block solution          *
 * blocks according to the Courant-Friedrichs-Lewy      *
 * condition.                                           *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
double Hexa_Multi_Block<HEXA_BLOCK>::CFL(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                                          typename HEXA_BLOCK::Soln_cState> &Input) {
  
   double dtMin;

   dtMin = MILLION;
  
   /* Determine the allowable time step for each solution block. */
   /* Prescribe boundary data for each solution block. */

   for (int nblk = 0; nblk < Number_of_Soln_Blks; ++nblk) { 
      if (Block_Used[nblk]) {
         dtMin = min(dtMin, Soln_Blks[nblk].CFL(Input));
      } /* endif */
   }  /* endfor */
   
   /* Return the global time step. */
   
   return (dtMin);
    
}
/********************************************************
 * Routine: Set_Global_TimeStep                         *
 *                                                      *
 * Assigns global time step to a 1D array of 2D         *
 * hexadrial multi-block solution blocks for            *
 * time-accurate calculations.                          *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
void Hexa_Multi_Block<HEXA_BLOCK>::Set_Global_TimeStep(const double &Dt_min){

   for (int nblk = 0; nblk < Number_of_Soln_Blks; ++nblk) { 
      if (Block_Used[nblk]) {
         Soln_Blks[nblk].Set_Global_TimeStep(Dt_min);
      } /* endif */
   }  /* endfor */

}

/********************************************************
 * Routine: WtoU                                        *
 *                                                      *
 * Convert primitive solution vector to conservative    *
 * solution vector.                                     *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
int Hexa_Multi_Block<HEXA_BLOCK>::WtoU(void){
   
   int i, error_flag;
   error_flag = 0;
   
   /* Convert U to W for each solution block. */

   for(int nblk = 0; nblk < Number_of_Soln_Blks; ++nblk){
      if (Block_Used[nblk]) {
         error_flag = Soln_Blks[nblk].WtoU();
         if (error_flag) return (error_flag);
      } /* endif */
   }  /* endfor */
   
   return(error_flag);
   
}

/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine evaluates the stage solution residual   *
 * for a 1D array of 3D hexahedrial multi-block       *
 * solution blocks.  A variety of multistage explicit   *
 * time integration and upwind finite-volume spatial    *
 * discretization procedures can be used depending on   *
 * the specified input values.                          *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
int Hexa_Multi_Block<HEXA_BLOCK>::dUdt_Multistage_Explicit(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                                                            typename HEXA_BLOCK::Soln_cState> &Input,
                                                           const int I_Stage) {
   
   int i, error_flag;
   error_flag = 0;
   
   /* Evaluate the solution residual for each solution block. */

   for (int nblk = 0; nblk<Number_of_Soln_Blks; ++nblk) {
      if (Block_Used[nblk]) {
         error_flag =  Soln_Blks[nblk].dUdt_Multistage_Explicit(I_Stage, Input);
         if (error_flag) return (error_flag);
      } /* endif */
   }  /* endfor */
   
   return(error_flag);
   
}

/********************************************************
 * Routine: Update_Solution_Multistage_Explicit         *
 *                                                      *
 * This routine updates the solution for a 1D array of  *
 * 3D hexahedrial multi-block solution blocks.  A       *
 * variety of multistage explicit time integration      *
 * and upwind finite-volume spatial discretization      *
 * procedures can be used depending on the specified    *
 * input values.                                        *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
int Hexa_Multi_Block<HEXA_BLOCK>::Update_Solution_Multistage_Explicit(Input_Parameters<typename HEXA_BLOCK::Soln_pState,
                                                                                       typename HEXA_BLOCK::Soln_cState> &Input,
                                                                      const int I_Stage) {
   
   int i, error_flag;
   error_flag = 0;
   
   /* Update the solution for each solution block. */

   for (int nblk = 0; nblk < Number_of_Soln_Blks; ++nblk) {
      if (Block_Used[nblk]) {
         error_flag =  Soln_Blks[nblk].Update_Solution_Multistage_Explicit(I_Stage, Input);
         if (error_flag) return (error_flag);
      } /* endif */
   }  /* endfor */
   
   return(error_flag);
   
}

/********************************************************
 * Routine: Read_Restart_Solution                       *
 *                                                      *
 * Reads restart solution file(s) and assigns values to *
 * the solution variables of a 1D array of 3D           *
 * quadrilateral multi-block solution blocks.           *
 * Returns a non-zero value if cannot read any of the   *
 * restart solution files.                              *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
int Hexa_Multi_Block<HEXA_BLOCK>::Read_Restart_Solution(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                                                         typename HEXA_BLOCK::Soln_cState> &Input,
                                                        AdaptiveBlock3D_List &Local_Adaptive_Block_List,
                                                        int &Number_of_Time_Steps,
                                                        double &Time,
                                                        CPUTime &CPU_Time) {
   
   int i, i_new_time_set, nsteps, nblk;
   char prefix[256], extension[256], restart_file_name[256], line[256];
   char *restart_file_name_ptr;
   ifstream restart_file;
   double time0;
   CPUTime cpu_time0;
   
   i = 0;
   while (1) {
      if (Input.Restart_File_Name[i] == ' ' ||
          Input.Restart_File_Name[i] == '.') break;
      prefix[i]=Input.Restart_File_Name[i];
      i = i + 1;
      if (i > strlen(Input.Restart_File_Name) ) break;
   } /* endwhile */
   prefix[i] = '\0';
   strcat(prefix, "_blk");
   
   /* Read the initial data for each solution block. */

   i_new_time_set = 0;

   for (nblk = 0; nblk < Local_Adaptive_Block_List.Nblk; ++nblk) { 
      if (Local_Adaptive_Block_List.Block[nblk].used == ADAPTIVEBLOCK3D_USED){
         sprintf(extension, "%.6d", Local_Adaptive_Block_List.Block[nblk].gblknum);
         strcat(extension, ".soln");
         strcpy(restart_file_name, prefix);
         strcat(restart_file_name, extension);
         restart_file_name_ptr = restart_file_name;
          
         // Open restart file.
         restart_file.open(restart_file_name_ptr, ios::in);
         if (restart_file.bad()) return (1);
          
         // Read iteration/time data.
         restart_file.setf(ios::skipws);
         restart_file >> nsteps >> time0 >> cpu_time0;
         restart_file.unsetf(ios::skipws);

         if (!i_new_time_set) {
            Number_of_Time_Steps = nsteps;
            Input.Maximum_Number_of_Time_Steps += Number_of_Time_Steps;
            Time = time0;
            CPU_Time.cput = cpu_time0.cput;
            i_new_time_set = 1;
         } /* endif */

         // Read reference solution states.
         Input.Read_Reference_Solution_States(restart_file);

         // Read solution block data.
         restart_file >> Soln_Blks[nblk];
              
         // Close restart file.
         restart_file.close();

         Soln_Blks[nblk].Wall_Shear();
       
       }  /* endif */
    } /* endfor */
    
    /* Reading of restart files complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Write_Restart_Solution                      *
 *                                                      *
 * Writes restart solution file(s) for a 1D array of 3D *
 * hexahedrial multi-block solution blocks.           *
 * Returns a non-zero value if cannot write any of the  *
 * restart solution files.                              *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
int Hexa_Multi_Block<HEXA_BLOCK>::Write_Restart_Solution(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                                                          typename HEXA_BLOCK::Soln_cState> &Input,
                                                         AdaptiveBlock3D_List &Local_Adaptive_Block_List,
                                                         const int Number_of_Time_Steps,
                                                         const double &Time,
                                                         const CPUTime &CPU_Time) {
   
   int i, nblk;
   char prefix[256], extension[256], restart_file_name[256];
   char *restart_file_name_ptr;
   ofstream restart_file;
   
   /* Determine prefix of restart file names. */
   nblk = 0;
     
   i = 0;
   while (1) {
      if (Input.Restart_File_Name[i] == ' ' ||
          Input.Restart_File_Name[i] == '.') break;
      prefix[i]=Input.Restart_File_Name[i];
      i = i + 1;
      if (i > strlen(Input.Restart_File_Name) ) break;
   } /* endwhile */
   prefix[i] = '\0';
   strcat(prefix, "_blk");
    
   for (nblk = 0; nblk < Local_Adaptive_Block_List.Nblk; ++nblk){
      if (Local_Adaptive_Block_List.Block[nblk].used == ADAPTIVEBLOCK3D_USED){    
         sprintf(extension, "%.6d", Local_Adaptive_Block_List.Block[nblk].gblknum);
         strcat(extension, ".soln");
         strcpy(restart_file_name, prefix);
         strcat(restart_file_name, extension);
         restart_file_name_ptr = restart_file_name;
         
         // Open restart file.
         restart_file.open(restart_file_name_ptr, ios::out);
         if (restart_file.bad()) return (1);
      
         // Write iteration/time data.
         restart_file.setf(ios::scientific);
         restart_file << setprecision(14) << Number_of_Time_Steps 
                      << " " << Time << " " << CPU_Time << "\n";
         restart_file.unsetf(ios::scientific);
         
         // Write reference solution states.
         Input.Write_Reference_Solution_States(restart_file);
         
         // Write solution block data.
         restart_file << setprecision(14) << Soln_Blks[nblk];
         
         // Close restart file.
         restart_file.close();

      }  /* endif */

   } /* endfor */
   
   /* Writing of restart files complete.  Return zero value. */

   return(0);
      
}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the nodal solution values for a 1D            *
 * array of 3D hexahedral multi-block solution          *
 * blocks to the specified output data file(s) in a     *
 * format suitable for plotting with TECPLOT.           *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
int Hexa_Multi_Block<HEXA_BLOCK>::Output_Tecplot(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                                                  typename HEXA_BLOCK::Soln_cState> &Input,
                                                 AdaptiveBlock3D_List &Local_Adaptive_Block_List,
                                                 const int Number_of_Time_Steps,
                                                 const double &Time) {
   
    int i, i_output_title, nblk;
    char prefix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    
   
    /* Determine prefix of output data file names. */

    i = 0;
    while (1) {
       if (Input.Output_File_Name[i] == ' ' ||
           Input.Output_File_Name[i] == '.') break;
       prefix[i]=Input.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(Input.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_cpu");
   
    /* Determine output data file name for this processor. */

    sprintf(extension, "%.6d", Local_Adaptive_Block_List.ThisCPU);
    strcat(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data for each solution block. */

    i_output_title = 1;
    for (nblk = 0; nblk<Local_Adaptive_Block_List.Nblk; ++nblk) {
       if (Local_Adaptive_Block_List.Block[nblk].used == ADAPTIVEBLOCK3D_USED) {    
          Soln_Blks[nblk].Output_Tecplot(Input,
                                         Number_of_Time_Steps, 
                                         Time,
                                         Local_Adaptive_Block_List.Block[nblk].gblknum,
                                         i_output_title,
                                         output_file);
          if (i_output_title) i_output_title = 0;
       } /* endif */
    }  /* endfor */

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cell centred solution values for a 1D     *
 * array of 3D hexahedral multi-block solution          *
 * blocks to the specified output data file(s) in a     *
 * format suitable for plotting with TECPLOT.           *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
int Hexa_Multi_Block<HEXA_BLOCK>::Output_Cells_Tecplot(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                                                        typename HEXA_BLOCK::Soln_cState> &Input,
                                                        AdaptiveBlock3D_List &Local_Adaptive_Block_List,
                                                        const int Number_of_Time_Steps,
                                                        const double &Time) {
   
    int i, i_output_title, nblk;
    char prefix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    

    /* Determine prefix of output data file names. */
 
    i = 0;
    while (1) {
       if (Input.Output_File_Name[i] == ' ' ||
           Input.Output_File_Name[i] == '.') break;
       prefix[i]=Input.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(Input.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_cells_cpu");
       /* Determine output data file name for this processor. */

    sprintf(extension, "%.6d", Local_Adaptive_Block_List.ThisCPU);
    strcat(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data for each solution block. */

    i_output_title = 1;
    for (nblk = 0; nblk < Local_Adaptive_Block_List.Nblk; ++nblk) {
       if (Local_Adaptive_Block_List.Block[nblk].used == ADAPTIVEBLOCK3D_USED) {    
         Soln_Blks[nblk].Output_Cells_Tecplot(Input,
                                              Number_of_Time_Steps, 
                                              Time,
                                              Local_Adaptive_Block_List.Block[nblk].gblknum,
                                              i_output_title,
                                              output_file);
         if (i_output_title) i_output_title = 0;
       } /* endif */
    }  /* endfor */

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the nodal solution values for a 1D            *
 * array of 3D hexahedral multi-block solution          *
 * blocks to the specified output data file(s) in a     *
 * format suitable for plotting with TECPLOT.           *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
int Hexa_Multi_Block<HEXA_BLOCK>::Output_Nodes_Tecplot(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                                                        typename HEXA_BLOCK::Soln_cState> &Input,
                                                       AdaptiveBlock3D_List &Local_Adaptive_Block_List,
                                                       const int Number_of_Time_Steps,
                                                       const double &Time) {
   
    int i, i_output_title, nblk;
    char prefix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    

    /* Determine prefix of output data file names. */
 
    i = 0;
    while (1) {
       if (Input.Output_File_Name[i] == ' ' ||
           Input.Output_File_Name[i] == '.') break;
       prefix[i]=Input.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(Input.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_nodes_cpu");
   
    /* Determine output data file name for this processor. */

    sprintf(extension, "%.6d", Local_Adaptive_Block_List.ThisCPU);
    strcat(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data for each solution block. */

    i_output_title = 1;
    for (nblk = 0; nblk < Local_Adaptive_Block_List.Nblk; ++nblk) {
       if (Local_Adaptive_Block_List.Block[nblk].used == ADAPTIVEBLOCK3D_USED) {    
         Soln_Blks[nblk].Output_Nodes_Tecplot(Input,
                                              Number_of_Time_Steps, 
                                              Time,
                                              Local_Adaptive_Block_List.Block[nblk].gblknum,
                                              i_output_title,
                                              output_file);
         if (i_output_title) i_output_title = 0;
       } /* endif */
    }  /* endfor */

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);

}

#endif // _HEXA_MULTIBLOCK_INCLUDED
