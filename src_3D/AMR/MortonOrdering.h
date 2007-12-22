/* MortonOrdering.h:  Header file defining subroutines used for
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
 * "int* OctreeBlock_DataStructure_Morton(const OctreeBlock_DataStructure &QTD, int *array ) "
 *
 * array = a pointer to an empty array 
 * (of length equal to the number of blocks) 
 * QTD = a pointer variable a OctreeBlock_DataStructure.
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

inline unsigned long long morton_number(int depth,unsigned int x, unsigned int y, unsigned int z) {
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

inline int* OctreeBlock_Morton(OctreeBlock &QTB, int *array_ptr) {
  if (QTB.block.used){               //no children then assign value
    array_ptr[0] = QTB.block.gblknum;       //assign Global adaptive block number into the correct position
    return (array_ptr +1);           //Shift the array pointer to point to an empty spot (array[1])
  }
  else{  //Go through all 8 children and assign positions in the array ptr
    array_ptr = OctreeBlock_Morton(*QTB.childBSW_ptr,array_ptr);
    array_ptr = OctreeBlock_Morton(*QTB.childBSE_ptr,array_ptr);
    array_ptr = OctreeBlock_Morton(*QTB.childBNW_ptr,array_ptr);
    array_ptr = OctreeBlock_Morton(*QTB.childBNE_ptr,array_ptr);
    array_ptr = OctreeBlock_Morton(*QTB.childTSW_ptr,array_ptr);
    array_ptr = OctreeBlock_Morton(*QTB.childTSE_ptr,array_ptr);
    array_ptr = OctreeBlock_Morton(*QTB.childTNW_ptr,array_ptr);
    array_ptr = OctreeBlock_Morton(*QTB.childTNE_ptr,array_ptr);
  }
  return array_ptr;
}


// Sort Function arranges the array into order of INCREASING morton number 
inline void sort(my_place *block_array, int length) {
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


inline int* Octree_DataStructure_Morton(Octree_DataStructure &QTD, int *array_ptr ) {
 int num_of_blocks = QTD.countUsedBlocks(); 
 my_place *block_array;
 block_array = new my_place[num_of_blocks];
 int roots_count = 0;
 int maxn;
  for(unsigned int x = 0; x<QTD.NRi; x++)
    for (unsigned int y = 0; y<QTD.NRj; y++)
      for (unsigned int z = 0; z<QTD.NRk; z++){
      if(QTD.Roots[x*y*z].childTNW_ptr != NULL || QTD.Roots[x*y*z].block.used == 1 ) {
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
    array_ptr = OctreeBlock_Morton(QTD.Roots[roots_count], 
                                   array_ptr);

  }

  return (array_ptr-1); //returns a pointer pointing to the last element in the array.
}

/*******************************************************************
 * This is the main function that rearanges the block structure on 
 * the parallel processors. 
 * 
 * Tim Blair June 9 2004
 ********************************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Morton_ReOrdering_of_Solution_Blocks(Octree_DataStructure                                    &Octree,
                                         AdaptiveBlock3D_List                                    &Local_Adaptive_Block_List,
					 Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
                                         Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>              &IPs,
                                         int                                                     &number_of_time_steps, 
                                         double                                                  &Time, 
                                         CPUTime                                                 &processor_cpu_time) {
   
     int BLK = 0;
     int CPU = 0;
     int error_flag = 0;

    /********************************************************  
     * Write the restart files                             *
     ********************************************************/
   
    // Write restart files.
    error_flag = Write_Octree(Octree,
			      IPs);

    if (error_flag) {
      cout << "\n  ERROR: Unable to open  octree data file "
	   << "on processor "
	   << Local_Adaptive_Block_List.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif 
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    error_flag = Local_Solution_Blocks.Write_Restart_Solution(IPs,
                                                              Local_Adaptive_Block_List,
							      number_of_time_steps,
   					                      Time,
 					                      processor_cpu_time);
    if (error_flag) {
      cout << "\n  ERROR: Unable to open  restart output data file(s) "
	   << "on processor "
	   << Local_Adaptive_Block_List.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    /********************************************************  
     * MORTON ORDERING FUNCTIONS                            *
     ********************************************************/
    
    int num_of_blocks;    
    num_of_blocks = Octree.countUsedBlocks();
    int *morton_array;
    morton_array = new int[num_of_blocks];
     int index;

    /********************************************************
     *REORDERING THE OCTREE Info                           *
     ********************************************************/
    Octree_DataStructure_Morton(Octree, morton_array);

    //Array to determine how many blocks go to each CPU
    int *cpu_array;
    cpu_array = new int[Octree.Ncpu];
    for (int i=0;i<Octree.Ncpu;i++)cpu_array[i]=0;
    index=0;
    for (int b=0; b<Octree.Nblk;b++)
      for (int c = 0; c<Octree.Ncpu; c++){
	if(index < num_of_blocks){
	  cpu_array[c] = cpu_array[c]+1;
	  index++;
	}
      }
     
    CPU = 0;
    BLK = 0;
    for (index=0; index <num_of_blocks; index++)
      for ( int oldCPU = 0; oldCPU < Octree.Ncpu; oldCPU++ )
	for ( int oldBLK = 0; oldBLK < Octree.Nblk; oldBLK++){
	  if (Octree.Blocks[oldCPU][oldBLK]!=NULL && morton_array[index] == Octree.Blocks[oldCPU][oldBLK]->block.gblknum){
	    Octree.Blocks[oldCPU][oldBLK]->block.info.cpu = CPU;  
	    Octree.Blocks[oldCPU][oldBLK]->block.info.blknum = BLK;
	    BLK++;
	    if(BLK == cpu_array[CPU]){
	      BLK = 0;
	      CPU++;
	    }
	  }
	}
 
    //This reassigns the .Blocks OctreeBlock pointers according to the cpu# and local blk # recorded in .info of each block.
    Octree.assign_block_pointers();
    
    //manipulate Octree & Local_Adaptive_Block_List, 
    
    //********** Need to recall these functions even though they were called within Read_Octree *********

   // For unstructured root blocks, the information of neighbours of root solution blocks
    // has been assigned in the create intial solution blocks. Xinfeng Gao Oct.18 2007.

/*     // Find the neighbours of the root blocks. */
/*     Octree_DataStructure::Find_Neighbours_of_Root_Solution_Blocks(Octree); */
    
    // Modify block neighbours for grid geometries with 
    //periodic boundaries, etc. 
    Octree_DataStructure::Modify_Neighbours_of_Root_Solution_Blocks(Octree, IPs.Grid_IP.i_Grid);
    
    // Determine the neighbouring blocks of all used (active)
    //solution blocks in the octree data structure. This will
    //also copy block information to local processor solution block list. 
    // Octree_DataStructure::Find_Neighbours(Octree, Local_Adaptive_Block_List);

    //******************************************************
    
    /*Allocates memory for all message passing buffers used  *
     * to send solution information between neighbouring      *
     * adaptive blocks.  */                                  
    
//     AdaptiveBlock3D_List::Allocate_Message_Buffers(Local_Adaptive_Block_List,
//                                                    Local_Solution_Blocks.Soln_Blks[0].NumVar()+NUM_COMP_VECTOR3D);

    /* Reads restart solution file(s) and assigns values to *
     * the solution variables of a 1D array of 3D           *
     * hexalateral multi-block solution blocks.             *
     * Returns a non-zero value if cannot read any of the   *
     * restart solution files.                              */
 

    error_flag = Local_Solution_Blocks.Read_Restart_Solution(IPs,
                                                             Local_Adaptive_Block_List,
                                                             number_of_time_steps,
				                             Time,
				                             processor_cpu_time);

    if (error_flag) {
      cout << "\n  ERROR: Unable to open  restart input data file(s) "
	   << "on processor "
	   << Local_Adaptive_Block_List.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif 
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    // Ensure each processor has the correct time and time!!!
    number_of_time_steps = CFFC_Maximum_MPI(number_of_time_steps);
    Time = CFFC_Maximum_MPI(Time);
    processor_cpu_time.cput = CFFC_Maximum_MPI(processor_cpu_time.cput);
    
    /* Send solution information between neighbouring blocks to complete
       prescription of initial data. */
    
/*     CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.  */
/*     // First send mesh and geometry information. */
/*     error_flag = Send_Messages_Mesh_Geometry_Only<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > */
/*                     (Local_Solution_Blocks,  */
/*                      Local_Adaptive_Block_List); */
/*     if (error_flag) { */
/*         cout << "\n ERROR: Message passing error during geometry intialization " */
/*              << "on processor " */
/*              << Local_Adaptive_Block_List.ThisCPU */
/*              << ".\n"; */
/*         cout.flush(); */
/*     } /\* endif *\/ */
/*     error_flag = CFFC_OR_MPI(error_flag); */
/*     if (error_flag) return (error_flag); */
/*     // Correct exterior nodes to match with message passed geometry information. */
/*     Solution_Data.Local_Solution_Blocks.Correct_Grid_Exterior_Nodes(Data.Local_Adaptive_Block_List);\ */
/*     // Now send solution information and data. */
/*     error_flag = Send_Messages<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > */
/*                     (Local_Solution_Blocks,  */
/*                      Local_Adaptive_Block_List); */
/*     if (error_flag) { */
/*        cout << "\n ERROR: Message passing error during solution intialization " */
/*             << "on processor " */
/*             << Local_Adaptive_Block_List.ThisCPU */
/*             << ".\n"; */
/*        cout.flush(); */
/*     } /\* endif *\/ */
/*     error_flag = CFFC_OR_MPI(error_flag); */
/*     if (error_flag) return (error_flag); */

    /* Prescribe boundary data consistent with initial data. */
    Local_Solution_Blocks.BCs(IPs);

    return (0);

}

/******************************************************************************************
   *This next bit outputs a curve which connects all the blocks in the cpu/blknum order     *
   ****************************************************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void Morton_SFC_Output_Tecplot3D(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>               &IPs,
                                 Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >  &Local_Solution_Blocks,
				 AdaptiveBlock3D_List                                     &Local_Adaptive_Block_List) {

  int i, BLK, i_centre, j_centre, k_centre;

  int num_of_blocks = 0;
  for ( BLK = 0; BLK <= Local_Adaptive_Block_List.Nblk-1; BLK++){
    if (Local_Adaptive_Block_List.Block[BLK].used == ADAPTIVEBLOCK3D_USED) num_of_blocks++;
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
     sprintf(extension, "%.4d", Local_Adaptive_Block_List.ThisCPU);
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

     Output_File << "GEOMETRY X=0, Y=0, Z=0, CS=GRID, C=BLUE, T=LINE3D, F=POINT"<<"\n"
	         <<"1"<<"\n"<< num_of_blocks <<"\n"; 
   
     for ( BLK = 0; BLK <= Local_Adaptive_Block_List.Nblk-1; BLK++){
       if (Local_Adaptive_Block_List.Block[BLK].used == ADAPTIVEBLOCK3D_USED) {

          i_centre = Local_Solution_Blocks.Soln_Blks[BLK].Grid.INl
                     +(Local_Solution_Blocks.Soln_Blks[BLK].Grid.INu
                     - Local_Solution_Blocks.Soln_Blks[BLK].Grid.INl)/2;

	  j_centre = Local_Solution_Blocks.Soln_Blks[BLK].Grid.JNl
                     +(Local_Solution_Blocks.Soln_Blks[BLK].Grid.JNu
                     - Local_Solution_Blocks.Soln_Blks[BLK].Grid.JNl)/2;
          k_centre = Local_Solution_Blocks.Soln_Blks[BLK].Grid.KNl
                     +(Local_Solution_Blocks.Soln_Blks[BLK].Grid.KNu
                     - Local_Solution_Blocks.Soln_Blks[BLK].Grid.KNl)/2;
          V = Local_Solution_Blocks.Soln_Blks[BLK].Grid. Node[i_centre][j_centre][k_centre].X;
          Output_File << V.x <<" "<< V.y <<" "<<V.z<<"\n"; //Output the centre node location for each used block

       } /* endif */
     } /* endfor */
     Output_File.close();
  } /* endif */
}

#endif /*_MORTON_ORDERING_INCLUDED*/
 
