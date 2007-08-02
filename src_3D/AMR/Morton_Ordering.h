/* Morton_Ordering.h:  Header file defining subroutines used for
                       applying a Morton re-ordering to the octree
                       solution blocks. */

/************************************************************************
 *                                                                      *
 *                         Written by Tim Blair                         *
 *                      email: blair@oddjob.utias.utoronto.ca           *
 *         University of Toronto Institute for Aerospace Studies        *
 *                           Summer 2004                                *
 *                                                                      * 
 ************************************************************************/

#ifndef _MORTON_ORDERING_INCLUDED
#define _MORTON_ORDERING_INCLUDED

#ifndef _HEXA_MULTIBLOCK_INCLUDED
#include "../HexaBlock/HexaMultiBlock.h"
#endif //_HEXA_MULTIBLOCK_INCLUDED

/*******************************************************************
 * "int* OcTreeBlock_DataStructure_Morton(const OcTreeBlock_DataStructure &QTD, int *array ) "
 *
 * array = a pointer to an empty array 
 * (of length equal to the number of blocks) 
 * QTD = a pointer variable a OcTreeBlock_DataStructure.
 * 
 * The function  returns a pointer to the last element in the array position in the array 
 * as well as filling the array with the global block numbers in a morton order
 * For use with Morton_Order() below
 * Writen by Tim Blair - May 31 2004
 *******************************************************************/
class my_place
{
public:
  int x,y,z;
  unsigned long long mornum;
};


inline unsigned long long morton_number(int depth,unsigned int x, unsigned int y, unsigned int z)
{
  unsigned long long mask = 1 << (depth - 1);
  unsigned long long result = 0;
  int b;
  for ( b = depth; b--; ) {
	result |= x & mask; result <<= 1;
	result |= y & mask; result <<= 1;
	result |= z & mask; mask >>= 1;
  }
  return result;
}

inline int* OcTreeBlock_Morton(OcTreeBlock &QTB, int *array_ptr)
{
  if (QTB.block.used){               //no children then assign value
    array_ptr[0] = QTB.block.gblknum;       //assign Global adaptive block number into the correct position
    return (array_ptr +1);           //Shift the array pointer to point to an empty spot (array[1])
  }
  else{  //Go through all 8 children and assign positions in the array ptr
    array_ptr = OcTreeBlock_Morton(*QTB.childBSW_ptr,array_ptr);
    array_ptr = OcTreeBlock_Morton(*QTB.childBSE_ptr,array_ptr);
    array_ptr = OcTreeBlock_Morton(*QTB.childBNW_ptr,array_ptr);
    array_ptr = OcTreeBlock_Morton(*QTB.childBNE_ptr,array_ptr);
    array_ptr = OcTreeBlock_Morton(*QTB.childTSW_ptr,array_ptr);
    array_ptr = OcTreeBlock_Morton(*QTB.childTSE_ptr,array_ptr);
    array_ptr = OcTreeBlock_Morton(*QTB.childTNW_ptr,array_ptr);
    array_ptr = OcTreeBlock_Morton(*QTB.childTNE_ptr,array_ptr);
  }
  return array_ptr;
}


// Sort Function arranges the array into order of INCREASING morton number 
inline void sort(my_place *block_array, int length)
{
  int flag = 0;
  my_place temp;
  while(flag == 0){
    flag = 1;
    for(int i=0; i<length-1; i++)
      if(block_array[i].mornum > block_array[i+1].mornum) //Then switch and change flag
	{
	  temp = block_array[i];
	  block_array[i]=block_array[i+1];
	  block_array[i+1]= temp;
	  flag = 0;
	}
  }
}


inline int* OcTreeBlock_DataStructure_Morton(OcTreeBlock_DataStructure &QTD, int *array_ptr ) 
{
 int num_of_blocks = QTD.countUsedBlocks(); 


 my_place *block_array;
 block_array = new my_place[num_of_blocks];
 int roots_count = 0;
 int maxn;
  for(unsigned int x = 0; x<QTD.NRi; x++)
    for (unsigned int y = 0; y<QTD.NRj; y++)
      for (unsigned int z = 0; z<QTD.NRk; z++){
      if(QTD.Roots[x][y][z].childTNW_ptr != NULL || QTD.Roots[x][y][z].block.used == 1 ) {
	//If it has no children and is not being used then it must be a 'hole'
	maxn = max(QTD.NRi,max(QTD.NRj,QTD.NRk));
	block_array[roots_count].x = x;
	block_array[roots_count].y = y;                
	block_array[roots_count].z = z;                
	block_array[roots_count].mornum = morton_number(maxn,x,y,z);
	roots_count++;
	}
      
    }
  num_of_blocks = roots_count; 

  sort(&block_array[0], roots_count);

  for(roots_count = 0; roots_count<num_of_blocks; roots_count++){
    array_ptr = OcTreeBlock_Morton(QTD.Roots[block_array[roots_count].x][block_array[roots_count].y][block_array[roots_count].z], array_ptr);

  }

  return (array_ptr-1); //returns a pointer pointing to the last element in the array.

}


/*******************************************************************
 * This is the main function that rearanges the block structure on the parallel processors. 
 * 
 * Tim Blair June 9 2004
 ********************************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
   int Morton_ReOrdering_of_Solution_Blocks(OcTreeBlock_DataStructure &OcTree,
                                            AdaptiveBlock3D_List        &List_of_Local_Solution_Blocks,
					    Hexa_MultiBlock<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Hexa_MultiBlock_List,
                                            Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs,
                                            int                         &number_of_time_steps, 
                                            double                      &Time, 
                                            CPUTime                     &processor_cpu_time) {
   
  //if(OcTree.Ncpu==1)  return (0); //No reason to do this if there is only one processor

  int BLK = 0;
  int CPU = 0;
  int error_flag = 0;
    
    /********************************************************  
     * Write the restart files                             *
     ********************************************************/
   
    // Write restart files.
    error_flag = Write_OcTree(OcTree,
			      IPs);

    if (error_flag) {
      cout << "\n  ERROR: Unable to open  octree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif 
//    error_flag = CFDkit_OR_MPI(error_flag);
//    if (error_flag) return (error_flag);

    error_flag = Hexa_MultiBlock_List.Write_Restart_Solution(IPs,
                                                             List_of_Local_Solution_Blocks,
							     number_of_time_steps,
   					                     Time,
 					                     processor_cpu_time);
    if (error_flag) {
      cout << "\n  ERROR: Unable to open  restart output data file(s) "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */
//    error_flag = CFDkit_OR_MPI(error_flag);
//    if (error_flag) return (error_flag);
    
    
    /********************************************************  
     * MORTON ORDERING FUNCTIONS                            *
     ********************************************************/
    
    int num_of_blocks;    
    num_of_blocks = OcTree.countUsedBlocks();
    int *morton_array;
    morton_array = new int[num_of_blocks];
     int index;

    /********************************************************
     *REORDERING THE OC TREE Info                           *
     ********************************************************/
    // cout<<"\n     *REORDERING THE OC TREE Info*";cout.flush();
    OcTreeBlock_DataStructure_Morton(OcTree, morton_array );

    //Array to determine how many blocks go to each CPU
    int *cpu_array;
    cpu_array = new int[OcTree.Ncpu];
    for (int i=0;i<OcTree.Ncpu;i++)cpu_array[i]=0;
    index=0;
    for (int b=0; b<OcTree.Nblk;b++)
      for (int c = 0; c<OcTree.Ncpu; c++){
	if(index < num_of_blocks){
	  cpu_array[c] = cpu_array[c]+1;
	  index++;
	}
      }
     
    CPU = 0;
    BLK = 0;
    for (index=0; index <num_of_blocks; index++)
      for ( int oldCPU = 0; oldCPU < OcTree.Ncpu; oldCPU++ )
	for ( int oldBLK = 0; oldBLK < OcTree.Nblk; oldBLK++){
	  if (OcTree.Blocks[oldCPU][oldBLK]!=NULL && morton_array[index] == OcTree.Blocks[oldCPU][oldBLK]->block.gblknum){
	    OcTree.Blocks[oldCPU][oldBLK]->block.info.cpu = CPU;  
	    OcTree.Blocks[oldCPU][oldBLK]->block.info.blknum = BLK;
	    BLK++;
	    if(BLK == cpu_array[CPU]){
	      BLK = 0;
	      CPU++;
	    }
	  }
	}
 
    
    // cout<<"\n** OcTree.assign_block_pointers()    ";cout.flush();
   
    //This reassigns the .Blocks OcTreeBlock pointers according to the cpu# and local blk # recorded in .info of each block.
    OcTree.assign_block_pointers();
    
    //manipulate Octree & List_of_Local_Solution_Blocks, 
    
    //********** Need to recall these functions even though they were called within Read_OcTree *********
   
    // Find the neighbours of the root blocks.
    //cout<<"\n ** Find_Neighbours_of_Root_Solution_Blocks(   ";cout.flush();
    OcTreeBlock_DataStructure::Find_Neighbours_of_Root_Solution_Blocks(OcTree);
    
    // Modify block neighbours for grid geometries with 
    //periodic boundaries, etc. 
    //cout<<"\n ** Modify_Neighbours_of_Root_Solution_Blocks   ";cout.flush();
    OcTreeBlock_DataStructure::Modify_Neighbours_of_Root_Solution_Blocks(OcTree, IPs.IP_Grid.i_Grid);
  
    // Determine the neighbouring blocks of all used (active)
    //solution blocks in the octree data structure. This will
    //also copy block information to local processor solution block	list. 
   //cout<<"\n  **Find_Neighbours   ";cout.flush();
    OcTreeBlock_DataStructure::Find_Neighbours(OcTree, List_of_Local_Solution_Blocks);

    //******************************************************
    
    
    
//    error_flag = CFDkit_OR_MPI(error_flag);
//    if (error_flag) return (error_flag);
  
    /*Allocates memory for all message passing buffers used  *
     * to send solution information between neighbouring      *
     * adaptive blocks.  */                                  
    
    //cout<<"\n  ** Allocate_Message_Buffers   ";cout.flush();
     AdaptiveBlock3D_List::Allocate_Message_Buffers(List_of_Local_Solution_Blocks,
                             Hexa_MultiBlock_List.Hexa_Block_List[0]->NumVar()
                             +NUM_COMP_VECTOR3D);

//    error_flag = CFDkit_OR_MPI(error_flag);
//    if (error_flag) return (error_flag);


    /* Reads restart solution file(s) and assigns values to *
     * the solution variables of a 1D array of 3D           *
     * hexalateral multi-block solution blocks.             *
     * Returns a non-zero value if cannot read any of the   *
     * restart solution files.                              */
 

    error_flag = Hexa_MultiBlock_List.Read_Restart_Solution(IPs,
                                                            List_of_Local_Solution_Blocks,
				                            number_of_time_steps,
				                            Time,
				                            processor_cpu_time);

    if (error_flag) {
      cout << "\n  ERROR: Unable to open  restart input data file(s) "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif 

/*     error_flag = CFDkit_OR_MPI(error_flag); */
/*     if (error_flag) return (error_flag); */
/*     // Ensure each processor has the correct time and time!!! */
/*     number_of_time_steps = CFDkit_Maximum_MPI(number_of_time_steps); */
/*     Time = CFDkit_Maximum_MPI(Time); */
/*     processor_cpu_time.cput = CFDkit_Maximum_MPI(processor_cpu_time.cput); */
    
    
    /* Send solution information between neighbouring blocks to complete
       prescription of initial data. */
    
/*     CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization. */
   
   
/*    cout <<" \n MO- Send_All_Messages() ON ";cout.flush();
      error_flag = Send_All_Messages(Local_SolnBlk, 
				   List_of_Local_Solution_Blocks,
				   NUM_COMP_VECTOR3D,
				   ON);
     cout<<"\n MO- Send_All_Messages() OFF ";cout.flush();
   if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk, 
						    List_of_Local_Solution_Blocks,
						    Local_SolnBlk[0].NumVar(),
						    OFF);
   if (error_flag) {
      cout << "\n  ERROR: Message passing error during  solution intialization "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif 
  
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
*/    

    /* Prescribe boundary data consistent with initial data. */
//     Hexa_MultiBlock_List.BCs(IPs);

    return (0);
}

/******************************************************************************************
   *This next bit outputs a curve which connects all the blocks in the cpu/blknum order     *
   ****************************************************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void Morton_SFC_Output_Tecplot3D( Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>  &IPs,
                                  Hexa_MultiBlock<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Hexa_MultiBlock_List,
				   AdaptiveBlock3D_List        &Soln_Block_List) {

  int i, BLK, i_centre, j_centre, k_centre;

  int num_of_blocks = 0;
  for ( BLK = 0; BLK <= Soln_Block_List.Nblk-1; BLK++){
    if (Soln_Block_List.Block[BLK].used == ADAPTIVEBLOCK3D_USED) num_of_blocks++;
  } /* endfor */
 
  ofstream Output_File;

  if (num_of_blocks > 0) {
     char prefix[256], extension[256], output_file_name[256];
     //prefix
     i = 0;
     while (1) {
        if (IPs.Output_File_Name[i] == ' ' ||
            IPs.Output_File_Name[i] == '.') break;
        prefix[i]=IPs.Output_File_Name[i];
        i = i + 1;
        if (i > strlen(IPs.Output_File_Name) ) break;
     } /* endwhile */
     prefix[i] = '\0';
     strcat(prefix, "_Morton_cpu");

     //extension
     sprintf(extension, "%.4d", Soln_Block_List.ThisCPU);
     strcat(extension, ".dat");

     //complete file name
     strcpy(output_file_name, prefix);
     strcat(output_file_name, extension);
     char *Output_File_Name = output_file_name;
  
     //open output file
     Output_File.open(Output_File_Name,ios::out | ios::trunc);
     Output_File << setprecision(14);
     Output_File << "TITLE = \"Morton SFC \"" << "\n"<< "VARIABLES = \"X\", \"Y\", \"Z\", \n";
     
     //Zone Floor is the dummy zone used to view the morton spacing curve individually
     //used for Sizing the 'Floor' for tecplot
     Vector3D V; //, minvec;

/*   //minvec =  Local_SolnBlk[0].Grid.Node[Local_SolnBlk[0].Grid.INl ][Local_SolnBlk[0].Grid.JNl].X; 

     Output_File << "ZONE T =\"Floor\", I=2, J=2, F=POINT"<<"\n"
 		 <<minvec.x<<" "<<minvec.y<<" "<<minvec.z<<" "<<minvec.x<<" "<<minvec.y<<" "<<minvec.z<<" "
	         <<minvec.x<<" "<<minvec.y<<" "<<minvec.z<<" "<<minvec.x<<" "<<minvec.y<<" "<<minvec.z<<"\n";

		 <<0<<" "<<0<<" "<<0<<" "
		 <<0<<" "<<1<<" "<<0<<" "
		 <<1<<" "<<1<<" "<<0<<" "
		 <<1<<" "<<0<<" "<<0<<" ";
*/

//num_of_blocks = IPs.Number_of_Blocks_Per_Processor;
     Output_File << "GEOMETRY X=0, Y=0, Z=0, CS=GRID, C=BLUE, T=LINE3D, F=POINT"<<"\n"
	         <<"1"<<"\n"<< num_of_blocks<<"\n"; 
   
//     for ( BLK = 0; BLK <= num_of_blocks-1; BLK++){
     for ( BLK = 0; BLK <= Soln_Block_List.Nblk-1; BLK++){
       if (Soln_Block_List.Block[BLK].used == ADAPTIVEBLOCK3D_USED) {

          i_centre = Hexa_MultiBlock_List.Hexa_Block_List[BLK]-> Grid -> INl
                     +(Hexa_MultiBlock_List.Hexa_Block_List[BLK]->Grid -> INu
                     - Hexa_MultiBlock_List.Hexa_Block_List[BLK]->Grid -> INl)/2;

	  j_centre = Hexa_MultiBlock_List.Hexa_Block_List[BLK]->Grid -> JNl
                     +(Hexa_MultiBlock_List.Hexa_Block_List[BLK]->Grid -> JNu
                     - Hexa_MultiBlock_List.Hexa_Block_List[BLK]->Grid -> JNl)/2;
          k_centre = Hexa_MultiBlock_List.Hexa_Block_List[BLK]->Grid -> KNl
                     +(Hexa_MultiBlock_List.Hexa_Block_List[BLK]->Grid -> KNu
                     - Hexa_MultiBlock_List.Hexa_Block_List[BLK]->Grid -> KNl)/2;
          V = Hexa_MultiBlock_List.Hexa_Block_List[BLK]->Grid -> Node[i_centre][j_centre][k_centre].X;
          Output_File << V.x <<" "<< V.y <<" "<<V.z<<"\n"; //Output the centre node location for each used block

       } /* endif */
     } /* endfor */
     Output_File.close();
  } /* endif */
}

#endif /*_MORTON_ORDERING_INCLUDED*/
 
