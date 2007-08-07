/* HexaMultiBlock.h:  Header file creating multiblock list. */

#ifndef _HEXA_MULTIBLOCK_INCLUDED
#define _HEXA_MULTIBLOCK_INCLUDED

/* Include required C++ libraries. */

using namespace std;

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
template<class HEXA_BLOCK>
class Hexa_MultiBlock{
   
  private:
  public:
     
   HEXA_BLOCK **Hexa_Block_List; 
   int Size_of_Block_List;
   int *Block_Used;
   
   // constructor ...
   Hexa_MultiBlock(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                   typename HEXA_BLOCK::Soln_cState> &IPs,
                   Grid3D_Hexa_Block ****Grid_ptr);

   Hexa_MultiBlock(Input_Parameters<typename HEXA_BLOCK::Soln_pState,
                    typename HEXA_BLOCK::Soln_cState> &IPs);
 
   //destructor ...
   ~Hexa_MultiBlock(){
      
      delete []Hexa_Block_List;
      Hexa_Block_List = NULL;
      
      
   }// endofdestructor...
   

   int Read_Restart_Solution(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                             typename HEXA_BLOCK::Soln_cState> &IPs,
                             AdaptiveBlock3D_List &Soln_Block_List,
                             int &Number_of_Time_Steps,
                             double &Time,
                             CPUTime &CPU_Time);
   int Write_Restart_Solution(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                              typename HEXA_BLOCK::Soln_cState> &IPs,
                              AdaptiveBlock3D_List &Soln_Block_List,
                              const int Number_of_Time_Steps,
                              const double &Time,
                              const CPUTime &CPU_Time);
   int Output_Cells_Tecplot(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                            typename HEXA_BLOCK::Soln_cState> &IPs,
                            AdaptiveBlock3D_List &Soln_Block_List,
                            const int Number_of_Time_Steps,
                            const double &Time);
   int Output_Tecplot(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                      typename HEXA_BLOCK::Soln_cState> &IPs,
                      AdaptiveBlock3D_List &Soln_Block_List,
                      const int Number_of_Time_Steps,
                      const double &Time);
   double L1_Norm_Residual(void);
   double L2_Norm_Residual(void);
   double Max_Norm_Residual(void);
   void Set_Global_TimeStep(const double &Dt_min);
   void ICs(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
            typename HEXA_BLOCK::Soln_cState> &IPs);
   void BCs(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
            typename HEXA_BLOCK::Soln_cState> &IPs);
   int WtoU(void);
   double CFL(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
              typename HEXA_BLOCK::Soln_cState> &IPs);
   int dUdt_Multistage_Explicit(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                                typename HEXA_BLOCK::Soln_cState> &IPs,
                                const int I_Stage);
   int Update_Solution_Multistage_Explicit
      (Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
       typename HEXA_BLOCK::Soln_cState> &IPs, 
       const int I_Stage);
   // flow type dependent memeber functions 
   int Create_Wall_Data(void); // calculate wall data for serial code (all blocks on one processor)
   // the other Create wall data template function is in Turbulence/*.h for parallel version
/* Copy the boundaries of different blocks on on processor */
/* This was written for testing multiblock serial code, i.e.
   running the code with multiblock but on only one processor */
   int Copy_MultiBlk_Boundary_Info(AdaptiveBlock3D_List &Soln_Block_List,
                                   Grid3D_Input_Parameters &IPs);
   
/*    int Send_All_Boundary_Info(MultiBlk_Connectivity &Multiblock_conn, */
/*                               Input_Parameters &Input_Parameters); */
/*    void Message_Passing_Datatype(MultiBlk_Connectivity &MultiBlock_Connectivity, */
/*                                  Input_Parameters &Input_Parameters); */
      
};

  
// create a multiblock list on one processor
template<class HEXA_BLOCK>
Hexa_MultiBlock<HEXA_BLOCK>::Hexa_MultiBlock(
   Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
   typename HEXA_BLOCK::Soln_cState> &IPs, 
   Grid3D_Hexa_Block ****Grid_ptr){
   
   HEXA_BLOCK *SBptr ; // pointer
   
   int nblk = 0;
   
   // assign the number of blocks on this processor
   Size_of_Block_List = IPs.Number_of_Blocks_Per_Processor;
   // dynamically allocate the memory for the array of pointers 
   // (the solution block pointers) based on the size
   Hexa_Block_List = new HEXA_BLOCK *[Size_of_Block_List];
   Block_Used = new int [Size_of_Block_List];

   for(int usage = 0; usage<Size_of_Block_List; ++usage){
      Block_Used[usage] = 0;
      
   }//0 not assigning a real solution block 1 assigning real blocks
   
   
   for (int  k = 0 ; k < IPs.IP_Grid.KBlk ; ++k ) {
      for (int  j = 0; j <IPs.IP_Grid.JBlk ; ++j ) {
         for (int i = 0 ; i <IPs.IP_Grid.IBlk; ++i ) {
            if (Grid_ptr[i][j][k]->Node != NULL) {
               
               SBptr = new HEXA_BLOCK
                  (IPs.IP_Grid.ICells,
                   IPs.IP_Grid.JCells,
                   IPs.IP_Grid.KCells,
                   IPs.IP_Grid.Nghost,
                   i, j, k, 
                   IPs.i_Flow_Type, 
                   Grid_ptr[i][j][k]);
               // indexing the solution
               Hexa_Block_List[nblk] = SBptr;
               Block_Used[nblk] = 1;
               
               ++nblk;
               
               
               // solution block address in the list
            }
               
         }//endofi
         
      }//endofj
      
   }//endofk
   
}
template<class HEXA_BLOCK>
Hexa_MultiBlock<HEXA_BLOCK>::Hexa_MultiBlock(
		   Input_Parameters<typename HEXA_BLOCK::Soln_pState,
		      typename HEXA_BLOCK::Soln_cState> &IPs){

	   HEXA_BLOCK *SBptr ; // pointer

	      int nblk = 0;

	         // assign the number of blocks on this processor
	   Size_of_Block_List = IPs.Number_of_Blocks_Per_Processor;
	  // dynamically allocate the memory for the array of pointers 
	  // (the solution block pointers) based on the size
	   Hexa_Block_List = new HEXA_BLOCK *[Size_of_Block_List];
	   Block_Used = new int [Size_of_Block_List];
	   for(int usage = 0; usage<Size_of_Block_List; ++usage){
	    Block_Used[usage] = 0;
	   }//0 not assigning a real solution block 1 assigning real blocks
	      //
	                                                                                                                                                                                                                                                                                                
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
   double Hexa_MultiBlock<HEXA_BLOCK>::L1_Norm_Residual(void) {
   
   
   double l1_norm;
   l1_norm = ZERO;
   
   /* Calculate the L1-norm.
      Sum the L1-norm for each solution block. */
   
   for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){
      if(Block_Used[nblk]){
         l1_norm += Hexa_Block_List[nblk]->L1_Norm_Residual();
      }
      
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
   double Hexa_MultiBlock<HEXA_BLOCK>::L2_Norm_Residual(void) {

   
   double l2_norm;

   l2_norm = ZERO;
   
   /* Sum the square of the L2-norm for each solution block. */
      
   for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){
      if(Block_Used[nblk]){
         
         l2_norm += Hexa_Block_List[nblk]->L2_Norm_Residual();
      }
      
      
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
   double Hexa_MultiBlock<HEXA_BLOCK>::Max_Norm_Residual(void) {
   
   double max_norm;
   
   max_norm = ZERO;
   
   /* Find the maximum norm for all solution blocks. */
   
   
   for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){
      
      if(Block_Used[nblk]){
         
         max_norm = max(max_norm, (Hexa_Block_List[nblk]->Max_Norm_Residual()));
      }
      
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
   void Hexa_MultiBlock<HEXA_BLOCK>:: 
ICs(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                   typename HEXA_BLOCK::Soln_cState> &IPs) {
   
   
   /* Assign initial data for each solution block. */
   for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){
      if(Block_Used[nblk]){
         // Set initial data.
         Hexa_Block_List[nblk]->ICs(IPs.i_ICs, IPs);
      }
      
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
void Hexa_MultiBlock<HEXA_BLOCK>:: BCs(
Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                   typename HEXA_BLOCK::Soln_cState> &IPs) {
   
   /* Prescribe boundary data for each solution block. */
   for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){
      if(Block_Used[nblk]){
         Hexa_Block_List[nblk]->BCs(IPs);
      }
      
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
   double Hexa_MultiBlock<HEXA_BLOCK>::
CFL(Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                   typename HEXA_BLOCK::Soln_cState> &IPs) {

  
   double dtMin;

   dtMin = MILLION;
  
   /* Determine the allowable time step for each solution block. */
   /* Prescribe boundary data for each solution block. */
   for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){ 
      if(Block_Used[nblk]){
         dtMin = min(dtMin,  Hexa_Block_List[nblk]->CFL(IPs));
      }
      
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
   void Hexa_MultiBlock<HEXA_BLOCK>::
Set_Global_TimeStep(const double &Dt_min){

   for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){ 
      if(Block_Used[nblk]){
         Hexa_Block_List[nblk]->Set_Global_TimeStep(Dt_min);
      }
      
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
   int Hexa_MultiBlock<HEXA_BLOCK>::WtoU(void){
   
   
   int i, error_flag;
   
   error_flag = 0;
   
   /* Convert U to W for each solution block. */
   for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){
      if(Block_Used[nblk]){
         error_flag = Hexa_Block_List[nblk]->WtoU();
         if (error_flag) return (error_flag);
      }
      
   }  /* endfor */
   
   /* Hexahedrial multi-block solution blocks
      successfully updated.  Return. */
   
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
   int Hexa_MultiBlock<HEXA_BLOCK>::
dUdt_Multistage_Explicit(
   Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                   typename HEXA_BLOCK::Soln_cState> &IPs,
   const int I_Stage) {
   
   int i, error_flag;

   error_flag = 0;
   
   /* Evaluate the stage solution residual for each solution block. */
   /* Prescribe boundary data for each solution block. */
   for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){
      if(Block_Used[nblk]){
         error_flag =  Hexa_Block_List[nblk]->dUdt_Multistage_Explicit(
            I_Stage, IPs);
         if (error_flag) return (error_flag);
      }
      
   }  /* endfor */
   
   /* Residuals for each hexahedrial multi-block solution block
      successfully calcualted.  Return. */
   
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
int Hexa_MultiBlock<HEXA_BLOCK>::
Update_Solution_Multistage_Explicit(
   Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
   typename HEXA_BLOCK::Soln_cState> &IPs,
   const int I_Stage) {
   
   int i, error_flag;
   
   error_flag = 0;
   
   /* Update the solution for each solution block. */
   for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){
      if(Block_Used[nblk]){
         error_flag =  Hexa_Block_List[nblk]->Update_Solution_Multistage_Explicit(
            I_Stage, IPs);
      
         if (error_flag) return (error_flag);
      }
      
   }  /* endfor */
   
   /* Hexahedrial multi-block solution blocks
      successfully updated.  Return. */
   
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
int Hexa_MultiBlock<HEXA_BLOCK>::
Read_Restart_Solution(
   Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
   typename HEXA_BLOCK::Soln_cState> &IPs,
   AdaptiveBlock3D_List &Soln_Block_List,
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
      if (IPs.Restart_File_Name[i] == ' ' ||
          IPs.Restart_File_Name[i] == '.') break;
      prefix[i]=IPs.Restart_File_Name[i];
      i = i + 1;
      if (i > strlen(IPs.Restart_File_Name) ) break;
   } /* endwhile */
   prefix[i] = '\0';
   strcat(prefix, "_blk");
   
   /* Read the initial data for each solution block. */

    i_new_time_set = 0;

    // index for the solution block number
    nblk = 0;
    
  /*   for(nblk = 0; nblk<Size_of_Block_List; ++nblk){ */
/*        if(Block_Used[nblk]){ */
/*        sprintf(extension, "%.6d", nblk); */

    for(nblk = 0; nblk<Soln_Block_List.Nblk; ++nblk){ 
       if(Soln_Block_List.Block[nblk].used == ADAPTIVEBLOCK3D_USED){
          sprintf(extension, "%.6d", Soln_Block_List.Block[nblk].gblknum);
          strcat(extension, ".soln");
          strcpy(restart_file_name, prefix);
          strcat(restart_file_name, extension);
          restart_file_name_ptr = restart_file_name;
          
          // Open restart file.
          
          restart_file.open(restart_file_name_ptr, ios::in);
          if (restart_file.bad()) return (1);
          
          // Read solution block data.
          restart_file.setf(ios::skipws);
          restart_file >> nsteps >> time0 >> cpu_time0;
          restart_file.unsetf(ios::skipws);

       //----------------
       /*Multispecies*/
       restart_file.getline(line,sizeof(line)); 
       // get reaction set name
       restart_file >>IPs.react_name;
       IPs.Wo.React.set_reactions(IPs.react_name);
       
       // multispecies but no reactions
       if( IPs.react_name == "NO_REACTIONS"){
          restart_file.setf(ios::skipws);
          int num_species;
          //get number of species
          restart_file >> num_species; 
          string *species = new string[num_species];
          //get species names 
          for(int k=0; k<num_species; k++){
             restart_file >> species[k];
          } 
          restart_file.unsetf(ios::skipws);  
          IPs.Wo.React.set_species(species,num_species);
          delete[] species; 
       }
	  

       if (!i_new_time_set) {
          Number_of_Time_Steps = nsteps;
          IPs.Maximum_Number_of_Time_Steps += Number_of_Time_Steps;
          Time = time0;
          CPU_Time.cput = cpu_time0.cput;
        
          i_new_time_set = 1;
       } /* endif */
      
       restart_file >>(*(Hexa_Block_List[nblk]));
              
       // Close restart file.
       restart_file.close();

       Hexa_Block_List[nblk]->Wall_Shear();
       
       }  /* enditerating */
    
          /* Reading of restart files complete.  Return zero value. */
    
    }
    
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
int Hexa_MultiBlock<HEXA_BLOCK>::
Write_Restart_Solution(
   Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
   typename HEXA_BLOCK::Soln_cState> &IPs,
   AdaptiveBlock3D_List &Soln_Block_List,
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
      if (IPs.Restart_File_Name[i] == ' ' ||
          IPs.Restart_File_Name[i] == '.') break;
      prefix[i]=IPs.Restart_File_Name[i];
      i = i + 1;
      if (i > strlen(IPs.Restart_File_Name) ) break;
   } /* endwhile */
   prefix[i] = '\0';
   strcat(prefix, "_blk");
    
/*    for(nblk = 0; nblk<Size_of_Block_List; ++nblk){ */
      
/*       if(Block_Used[nblk]){ */
         
/*          sprintf(extension, "%.6d", nblk); */


         for(nblk = 0; nblk<Soln_Block_List.Nblk; ++nblk){
            if(Soln_Block_List.Block[nblk].used == ADAPTIVEBLOCK3D_USED){    
               sprintf(extension, "%.6d", Soln_Block_List.Block[nblk].gblknum);
         strcat(extension, ".soln");
         strcpy(restart_file_name, prefix);
         strcat(restart_file_name, extension);
         restart_file_name_ptr = restart_file_name;
         
         // Open restart file.
         restart_file.open(restart_file_name_ptr, ios::out);
         if (restart_file.bad()) return (1);
         
      
         // Write solution block data.
         restart_file.setf(ios::scientific);
         restart_file << setprecision(14) << Number_of_Time_Steps 
                      << " " << Time << " " << CPU_Time << "\n";
         restart_file.unsetf(ios::scientific);
         
         restart_file << IPs.react_name << "\n";
         if(IPs.react_name == "NO_REACTIONS"){
            restart_file << IPs.Wo.ns <<" ";
            for(int k=0; k< IPs.Wo.ns; k++){ 
               restart_file << IPs.multispecies[k] <<" ";
            }
            restart_file<<endl;
         }
         
         restart_file << setprecision(14) << (*(Hexa_Block_List[nblk]));
         
         // Close restart file.
         restart_file.close();


         
      }  /* endfor */


      
      /* Writing of restart files complete.  Return zero value. */
   }
   
   return(0);
      
}



/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cell centred solution values for a 1D     *
 * array of 3D hexahedrial   multi-block solution       *
 * blocks to the specified output data file(s) in a     *
 * format suitable for plotting with TECPLOT.           *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
int Hexa_MultiBlock<HEXA_BLOCK>::
                      Output_Cells_Tecplot(
Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
                   typename HEXA_BLOCK::Soln_cState> &IPs,
AdaptiveBlock3D_List &Soln_Block_List,
   const int Number_of_Time_Steps,
   const double &Time) {
   
   int i, i_output_title, nblk;
   char prefix[256], extension[256], output_file_name[256];
   char *output_file_name_ptr;
   ofstream output_file;    

   /* Determine prefix of output data file names. */
 
   nblk = 0;
     
   i = 0;
   
   while (1) {
      if (IPs.Output_File_Name[i] == ' ' ||
          IPs.Output_File_Name[i] == '.') break;
      prefix[i]=IPs.Output_File_Name[i];
      i = i + 1;
      if (i > strlen(IPs.Output_File_Name) ) break;
   } /* endwhile */
   prefix[i] = '\0';
   strcat(prefix, "_cells_cpu");
   
   /* Determine output data file name for this processor. */
   //only generate one file for all the solution blocks on a processor.
   
   /* Write the solution data for each solution block. */
   i_output_title = 1;
 /*   for(nblk = 0; nblk<Size_of_Block_List; ++nblk){ */
         
/*       if(Block_Used[nblk]){ */
              
/*          // generating a file for each block */
/*          sprintf(extension, "%.6d", nblk); */

  for(nblk = 0; nblk<Soln_Block_List.Nblk; ++nblk){
            if(Soln_Block_List.Block[nblk].used == ADAPTIVEBLOCK3D_USED){    
               sprintf(extension, "%.6d", Soln_Block_List.Block[nblk].gblknum);

         strcat(extension, ".dat");
         strcpy(output_file_name, prefix);
         strcat(output_file_name, extension);
         output_file_name_ptr = output_file_name;
         
         // Open the output data file.
         output_file.open(output_file_name_ptr, ios::out);
         if (output_file.bad()) return (1);
         i_output_title = 1;
         
         /* Write the solution data for each solution block. */
         
         Hexa_Block_List[nblk]->Output_Cells_Tecplot(Number_of_Time_Steps, 
                                                     Time,
                                                     Soln_Block_List.Block[nblk].gblknum,
                                                     i_output_title,
                                                     output_file);
         
         if (i_output_title) i_output_title = 0;
         
         
         output_file.close();
         
      }
   }
   

   /* Close the output data file. */
   // output_file.close();
   /* Writing of output data files complete.  Return zero value. */
   
   return(0);

}
/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the nodal solution values for a 1D            *
 * array of 3D hexahedrial   multi-block solution       *
 * blocks to the specified output data file(s) in a     *
 * format suitable for plotting with TECPLOT.           *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
template<class HEXA_BLOCK>
int Hexa_MultiBlock<HEXA_BLOCK>::
Output_Tecplot(
   Input_Parameters<typename HEXA_BLOCK::Soln_pState, 
   typename HEXA_BLOCK::Soln_cState> &IPs,
   AdaptiveBlock3D_List &Soln_Block_List,
   const int Number_of_Time_Steps,
   const double &Time) {
   
   int i, i_output_title, nblk;
   char prefix[256], extension[256], output_file_name[256];
   char *output_file_name_ptr;
   ofstream output_file;    
   
   /* Determine prefix of output data file names. */
   
   nblk = 0;
   
   i = 0;
   
   while (1) {
      if (IPs.Output_File_Name[i] == ' ' ||
          IPs.Output_File_Name[i] == '.') break;
      prefix[i]=IPs.Output_File_Name[i];
      i = i + 1;
      if (i > strlen(IPs.Output_File_Name) ) break;
   } /* endwhile */
   prefix[i] = '\0';
   strcat(prefix, "_nodes_cpu");
   
   /* Determine output data file name for this processor. */
   //only generate one file for all the solution blocks on a processor.
   
   /* Write the solution data for each solution block. */
   i_output_title = 1;

   for(nblk = 0; nblk<Soln_Block_List.Nblk; ++nblk){
      if(Soln_Block_List.Block[nblk].used == ADAPTIVEBLOCK3D_USED){    
         sprintf(extension, "%.6d", Soln_Block_List.Block[nblk].gblknum);
         
         strcat(extension, ".dat");
         strcpy(output_file_name, prefix);
         strcat(output_file_name, extension);
         output_file_name_ptr = output_file_name;
         
         // Open the output data file.
         output_file.open(output_file_name_ptr, ios::out);
         if (output_file.bad()) return (1);
         i_output_title = 1;
         
         /* Write the solution data for each solution block. */
              
         Hexa_Block_List[nblk]->Output_Tecplot(IPs,
                                               Number_of_Time_Steps, 
                                               Time,
                                               Soln_Block_List.Block[nblk].gblknum,
                                               i_output_title,
                                               output_file);
         
         if (i_output_title) i_output_title = 0;
         
         output_file.close();
         
      }
   }
   

   /* Close the output data file. */
   // output_file.close();
   /* Writing of output data files complete.  Return zero value. */
   
   return(0);

}


template<class HEXA_BLOCK>
int Hexa_MultiBlock<HEXA_BLOCK>::Create_Wall_Data(void){
   
   if( Hexa_Block_List[0]->Flow_Type != FLOWTYPE_TURBULENT_RANS_K_OMEGA){
      
      return (0);
      
   }else{
      
      int error_flag;
      double y_wall_temp;
      Vector3D X_wall_temp, n_wall_temp;
      int BC_wall_temp;
      
      y_wall_temp = std::numeric_limits<double>::max();
      for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){
         if(Block_Used[nblk]){
         
            
            // compute wall data ... ... 
            for ( int k = Hexa_Block_List[nblk]->Grid->KCl- Hexa_Block_List[nblk]->Grid->Nghost ; k <=   Hexa_Block_List[nblk]->Grid->KCu+  Hexa_Block_List[nblk]->Grid->Nghost ; ++k )
               for ( int j = Hexa_Block_List[nblk]->Grid->JCl- Hexa_Block_List[nblk]->Grid->Nghost ; j <=   Hexa_Block_List[nblk]->Grid->JCu+  Hexa_Block_List[nblk]->Grid->Nghost ; ++j ) 
                  for ( int i = Hexa_Block_List[nblk]->Grid->ICl-  Hexa_Block_List[nblk]->Grid->Nghost ; i <=   Hexa_Block_List[nblk]->Grid->ICu+  Hexa_Block_List[nblk]->Grid->Nghost ; ++i ) {
                     
                 
                     Hexa_Block_List[nblk]->Wall[i][j][k].ywall =  1.0e70;
                           
                     for(int iblk = 0; iblk<Size_of_Block_List; ++iblk){
                        if(Block_Used[iblk]){
                           
                           error_flag = Wall_Distance(Hexa_Block_List[iblk],
                                                      Hexa_Block_List[nblk]->Grid->Cell[i][j][k].Xc,
                                                      X_wall_temp, n_wall_temp,
                                                      y_wall_temp, BC_wall_temp);

                           if (y_wall_temp < Hexa_Block_List[nblk]->Wall[i][j][k].ywall ) {
                              Hexa_Block_List[nblk]->Wall[i][j][k].ywall = y_wall_temp;
                              Hexa_Block_List[nblk]->Wall[i][j][k].Xwall = X_wall_temp;
                              Hexa_Block_List[nblk]->Wall[i][j][k].nwall = n_wall_temp;
                              Hexa_Block_List[nblk]->Wall[i][j][k].BCwall = BC_wall_temp;
                           }
                        }
                        
                        
                     }
                     
                         
                     
                  }
                    
            // compute y+ for each cell
                      
         }// compute for used blocks
      }//loop through all blocks
           
   }//turbulent case
   

   return 0;
}

/* Copy the boundaries of different blocks on on processor */
/* This was written for testing multiblock serial code, i.e.
   running the code with multiblock but on only one processor */
template<class HEXA_BLOCK>
int Hexa_MultiBlock<HEXA_BLOCK>::Copy_MultiBlk_Boundary_Info(  AdaptiveBlock3D_List &Soln_Block_List,
                                                               Grid3D_Input_Parameters &IPs){
   
   
  
    
   HEXA_BLOCK *SBptr; // pointer - the target solution block
   HEXA_BLOCK *SBptr_n; // neighbor blocks of this target block

   int iCellLocal, iCellDonor, iCellNum;
   int jCellLocal, jCellDonor, jCellNum;
   int kCellLocal, kCellDonor, kCellNum;
   int icloc, jcloc, kcloc, icdon, jcdon,kcdon;
   
   for (int  bi = 0 ; bi <IPs.IBlk ; ++bi ) {
      for (int  bj = 0; bj <IPs.JBlk ; ++bj ) {
         for (int bk = 0 ; bk <IPs.KBlk ; ++bk ) {  

            for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){ 
               if(Block_Used[nblk]){
                  if( Hexa_Block_List[nblk]->Iindex == bi && 
                      Hexa_Block_List[nblk]->Jindex == bj &&
                      Hexa_Block_List[nblk]->Kindex == bk)
                     SBptr  =   Hexa_Block_List[nblk];
               }
            }//get the target solution block info
           
            // loop through the neigbour blocks
           for (int  nbi = bi-1 ; nbi <=bi+1 ; ++nbi ) {
              for (int  nbj = bj-1; nbj <=bj+1 ; ++nbj ) {
                 for (int nbk = bk-1 ; nbk <=bk+1 ; ++nbk ) {
/*                       test the [nbi, nbj, nbk] exists and that  */
/*                       [nbi, nbj, nbk] ! = [bi, bj, bk] */
/*                       loop through the list to get the index for  */
/*                       the neighour solution blocks */
/*                       or check if this neighbour block exists... */
                    if( (nbi<0 || nbi>=IPs.IBlk)
                        ||(nbj<0 || nbj>=IPs.JBlk)
                        ||(nbk<0 || nbk>=IPs.KBlk)
                        ||(nbi == bi && nbj == bj && nbk == bk)){
                       
                       // do nothing
                    }else{
                       // checking the index for the neighboring blocks
                       for(int nblk = 0; nblk<Size_of_Block_List; ++nblk){ 
                          if(Block_Used[nblk]){
                             if( Hexa_Block_List[nblk]->Iindex == nbi && 
                                 Hexa_Block_List[nblk]->Jindex == nbj &&
                                 Hexa_Block_List[nblk]->Kindex == nbk)
                                SBptr_n  =   Hexa_Block_List[nblk];
                          }
                       }//get the target solution block info
           
                        
                       if(nbi > bi){
                          iCellLocal = SBptr->ICu+1;
                          iCellDonor = SBptr_n->ICl;
                          iCellNum =  SBptr->Nghost;
                           
                       }else if(nbi < bi){
                          iCellLocal = SBptr->ICl - SBptr->Nghost;
                          iCellDonor = SBptr_n->ICu - SBptr_n->Nghost +1 ;
                          iCellNum =  SBptr->Nghost;
                           
                       }else{
                          iCellLocal = SBptr->ICl;
                          iCellDonor = SBptr_n->ICl;
                          iCellNum =  IPs.ICells;
                           
                       }// i direction
                       if(nbj > bj){
                          jCellLocal = SBptr->JCu+1;
                          jCellDonor = SBptr_n->JCl;
                          jCellNum =  SBptr->Nghost;
                       }else if(nbj < bj){
                          jCellLocal = SBptr->JCl - SBptr->Nghost;
                          jCellDonor = SBptr_n->JCu - SBptr_n->Nghost+1 ;
                          jCellNum =  SBptr->Nghost;
                       }else{
                          jCellLocal = SBptr->JCl;
                          jCellDonor = SBptr_n->JCl;
                          jCellNum =  IPs.JCells;
                       }// j direction
                       if(nbk > bk){
                          kCellLocal = SBptr->KCu+1;
                          kCellDonor = SBptr_n->KCl;
                          kCellNum =  SBptr->Nghost;
                       }else if(nbk < bk){
                          kCellLocal = SBptr->KCl - SBptr->Nghost;
                          kCellDonor = SBptr_n->KCu - SBptr_n->Nghost +1 ;
                          kCellNum =  SBptr->Nghost;
                       }else{
                          kCellLocal = SBptr->KCl;
                          kCellDonor = SBptr_n->KCl;
                          kCellNum = IPs.KCells;
                       }// k direction
                       icloc = iCellLocal;
                       icdon = iCellDonor;
                       // copying the donor cell information to the local cells
                       for (int  i = 0 ; i < iCellNum ; ++i ) {
                          jcloc = jCellLocal;
                          jcdon = jCellDonor;
                          for (int  j = 0 ; j < jCellNum ; ++j ) {
                             kcloc = kCellLocal;
                             kcdon = kCellDonor;
                             for (int k = 0 ; k < kCellNum ; ++k ) {
                                SBptr->U[icloc][jcloc][kcloc] =
                                   SBptr_n->U[icdon][jcdon][kcdon];
                                SBptr->W[icloc][jcloc][kcloc] =
                                   SBptr_n->W[icdon][jcdon][kcdon];

                                ++kcloc;
                                ++kcdon;
                                 
                             }
                             ++jcloc;
                             ++jcdon;
                          }
                          ++icloc;
                          ++icdon;
                       }//endofcopying information from donor cells to local cells
                    }//executable situation (exists) 
                 }
              }
            }//endofchecking the neighbour
        }
     }
   } //endofloopingthroughlocalblocks
   
   return 0;
   
}

#endif /* _HEXA_MULTIBLOCK_INCLUDED  */
