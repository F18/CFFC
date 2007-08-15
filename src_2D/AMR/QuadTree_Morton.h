/* QuadTree_Morton.h:  Header file defining subroutines used for
                       applying a Morton re-ordering to the quadtree
                       solution blocks. */

/************************************************************************
 *                                                                      *
 *                         Written by Tim Blair                         *
 *                      email: blair@oddjob.utias.utoronto.ca           *
 *         University of Toronto Institute for Aerospace Studies        *
 *                           Summer 2004                                *
 *                                                                      * 
 ************************************************************************/

#ifndef _QUADTREE_MORTON_INCLUDED
#define _QUADTREE_MORTON_INCLUDED

/*******************************************************************
 * "int* QuadTreeBlock_DataStructure_Morton(const QuadTreeBlock_DataStructure &QTD, int *array ) "
 *
 * array = a pointer to an empty array 
 * (of length equal to the number of blocks) 
 * QTD = a pointer variable a QuadTreeBlock_DataStructure.
 * 
 * The function  returns a pointer to the last element in the array position in the array 
 * as well as filling the array with the global block numbers in a morton order
 * For use with Morton_Order() below
 * Writen by Tim Blair - May 31 2004
 *******************************************************************/
class my_place
{
public:
  int x,y;
  unsigned long long mornum;
};


inline unsigned long long morton_number(int depth,unsigned int x, unsigned int y)
{
  unsigned long long mask = 1 << (depth - 1);
  unsigned long long result = 0;
  int b;
  for ( b = depth; b--; ) {
    result |= y & mask; result <<= 1;
    result |= x & mask; mask >>= 1;
  }
  return result;
}

inline int* QuadTreeBlock_Morton(QuadTreeBlock &QTB, int *array_ptr)
{
  if (QTB.block.used){               //no children then assign value
    array_ptr[0] = QTB.block.gblknum;       //assign Global adaptive block number into the correct position
    return (array_ptr +1);           //Shit the array pointer to point to an empty spot (array[1])
  }
  else{  //Go through all 4 children and assign positions in the array ptr
    array_ptr = QuadTreeBlock_Morton(*QTB.childSW_ptr,array_ptr);
    array_ptr = QuadTreeBlock_Morton(*QTB.childSE_ptr,array_ptr);
    array_ptr = QuadTreeBlock_Morton(*QTB.childNW_ptr,array_ptr);
    array_ptr = QuadTreeBlock_Morton(*QTB.childNE_ptr,array_ptr);
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


inline int* QuadTreeBlock_DataStructure_Morton(QuadTreeBlock_DataStructure &QTD, int *array_ptr ) 
{
 int num_of_blocks = QTD.countUsedBlocks(); 


 my_place *block_array;
 block_array = new my_place[num_of_blocks];
 int roots_count = 0;

  for(unsigned int x = 0; x<QTD.NRi; x++)
    for (unsigned int y = 0; y<QTD.NRj; y++){
      if(QTD.Roots[x][y].childNW_ptr != NULL || QTD.Roots[x][y].block.used == 1 ) {
	//If it has no children and is not begin used then it must be a 'hole'
	block_array[roots_count].x = x;
	block_array[roots_count].y = y;                
        block_array[roots_count].mornum = morton_number(max(QTD.NRj,QTD.NRi),x,y);
	roots_count++;
	}
      
    }
  num_of_blocks = roots_count; 

  sort(&block_array[0], roots_count);

  for(roots_count = 0; roots_count<num_of_blocks; roots_count++){
    array_ptr = QuadTreeBlock_Morton(QTD.Roots[block_array[roots_count].x][block_array[roots_count].y], array_ptr);

  }

  return (array_ptr-1); //returns a pointer pointing to the last element in the array.

}


/*******************************************************************
 * This is the main function that rearanges the block structure on the parallel processors. 
 * 
 * Tim Blair June 9 2004
 ********************************************************************/

template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
  int Morton_ReOrdering_of_Solution_Blocks(QuadTreeBlock_DataStructure &QuadTree,
		                           AdaptiveBlockResourceList   &List_of_Global_Solution_Blocks,
		                           AdaptiveBlock2D_List        &List_of_Local_Solution_Blocks,
		                           Quad_Soln_Block             *Local_SolnBlk, 
                                           Quad_Soln_Input_Parameters  &Input_Parameters,
		                           int                         &number_of_time_steps, 
                                           double                      &Time, 
                                           CPUTime                     &processor_cpu_time) {

  //if(QuadTree.Ncpu==1)  return (0); //No reason to do this if there is only one processor

  int BLK = 0;
  int CPU = 0;
  int error_flag = 0;

    
    /********************************************************  
     * Write the restart files                             *
     ********************************************************/
    
    
    // Write restart files.
    error_flag = Write_QuadTree(QuadTree,
				Input_Parameters);
    if (error_flag) {
      cout << "\n Euler2D ERROR: Unable to open Euler2D quadtree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    error_flag = Write_Restart_Solution(Local_SolnBlk, 
					List_of_Local_Solution_Blocks, 
					Input_Parameters,
					number_of_time_steps,
					Time,
					processor_cpu_time);
    if (error_flag) {
      cout << "\n Euler2D ERROR: Unable to open Euler2D restart output data file(s) "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    
    /********************************************************  
     * TIM'S MORTON ORDERING FUNCTIONS                            *
     ********************************************************/
    
    int num_of_blocks;    
    num_of_blocks = QuadTree.countUsedBlocks();
    int *morton_array;
    morton_array = new int[num_of_blocks];
     int index;
    
    
    /********************************************************  
     *REORDERING THE QUAD TREE Info  -Tim                   *
     ********************************************************/
      QuadTreeBlock_DataStructure_Morton(QuadTree, morton_array );

    //Array to determine how many blocks go to each CPU
    int *cpu_array;
    cpu_array = new int[QuadTree.Ncpu];
    for (int i=0;i<QuadTree.Ncpu;i++)cpu_array[i]=0;
    index=0;
    for (int b=0; b<QuadTree.Nblk;b++)
      for (int c = 0; c<QuadTree.Ncpu; c++){
	if(index < num_of_blocks){
	  cpu_array[c] = cpu_array[c]+1;
	  index++;
	}}
     
    CPU = 0;
    BLK = 0;
    for (index=0; index <num_of_blocks; index++)
      for ( int oldCPU = 0; oldCPU < QuadTree.Ncpu; oldCPU++ )
	for ( int oldBLK = 0; oldBLK < QuadTree.Nblk; oldBLK++){
	  if (QuadTree.Blocks[oldCPU][oldBLK]!=NULL && morton_array[index] == QuadTree.Blocks[oldCPU][oldBLK]->block.gblknum){
	    QuadTree.Blocks[oldCPU][oldBLK]->block.info.cpu = CPU;  
	    QuadTree.Blocks[oldCPU][oldBLK]->block.info.blknum = BLK;
	    BLK++;
	    if(BLK == cpu_array[CPU]){
	      BLK = 0;
	      CPU++;
	    }
	  }
	}
    
    
    
    //This reassigns the .Blocks QuadTreeBlock pointers according to the cpu# and local blk # recorded in .info of each block.
    QuadTree.assign_block_pointers();
    
    //manipulate Quadtree & List_of_Local_Solution_Blocks, 
    
    //********** Need to recall these functions even though they were called within Read_QuadTree *********
    
    // Find the neighbours of the root blocks.
    Find_Neighbours_of_Root_Solution_Blocks(QuadTree);
    
    // Modify block neighbours for grid geometries with 
    //periodic boundaries, etc. 
    Modify_Neighbours_of_Root_Solution_Blocks(QuadTree, Input_Parameters.i_Grid);
    
    // Determine the neighbouring blocks of all used (active)
    //solution blocks in the quadtree data structure. This will
    //also copy block information to local processor solution block	list. 
    Find_Neighbours(QuadTree, List_of_Local_Solution_Blocks);
    //******************************************************
    
    
    
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    /*Allocates memory for all message passing buffers used  *
     * to send solution information between neighbouring      *
     * adaptive blocks.  */                                  
    
    Allocate_Message_Buffers(List_of_Local_Solution_Blocks, Local_SolnBlk[0].NumVar()+NUM_COMP_VECTOR2D);
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    
    /* Reads restart solution file(s) and assigns values to *
     * the solution variables of a 1D array of 2D           *
     * quadrilateral multi-block solution blocks.           *
     * Returns a non-zero value if cannot read any of the   *
     * restart solution files.     */
    error_flag = Read_Restart_Solution(Local_SolnBlk, 
				       List_of_Local_Solution_Blocks, 
				       Input_Parameters,
				       number_of_time_steps,
				       Time,
				       processor_cpu_time);
    if (error_flag) {
      cout << "\n Euler2D ERROR: Unable to open Euler2D restart input data file(s) "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
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
    
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    
    error_flag = Send_All_Messages(Local_SolnBlk, 
				   List_of_Local_Solution_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);
    if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk, 
						    List_of_Local_Solution_Blocks,
						    Local_SolnBlk[0].NumVar(),
						    OFF);
    if (error_flag) {
      cout << "\n Euler2D ERROR: Message passing error during Euler2D solution intialization "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */
  
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    /* Prescribe boundary data consistent with initial data. */
    
    BCs(Local_SolnBlk, 
	List_of_Local_Solution_Blocks,
	Input_Parameters);
    
    /***********************************
     *Change Global Block Numbers    *
     **********************************/
    /* 
    if(1){//Enter flag for changing global numbers so that they may be reread correctly
      
    Euler2D_Quad_Block *temp_SolnBlk;
    num_of_blocks = QuadTree.countUsedBlocks();
    temp_SolnBlk = new Euler2D_Quad_Block[num_of_blocks];
    
    
    cout << "Order of CPU's and global block num's: ";
    for ( CPU = 0; CPU < QuadTree.Ncpu; CPU++ ){
      cout << "\nCPU "<<CPU<<" Has Blocks ";
      for ( BLK = 0; BLK < QuadTree.Nblk; BLK++){
      if(QuadTree.Blocks[CPU][BLK] != NULL)
      cout<<QuadTree.Blocks[CPU][BLK]->block.gblknum << " ";
      }
    }
    cout<<endl;
    

    for (int i=0;i<num_of_blocks;i++)
      temp_SolnBlk[i]=Local_SolnBlk[i];

     // Renumber all solution blocks, assigning a unique global block number. 
    Renumber_Solution_Blocks(QuadTree,List_of_Local_Solution_Blocks);
  }*/

    return (0);
}

/******************************************************************************************
   *This next bit outputs a curve which connects all the blocks in the cpu/blknum order     *
   ****************************************************************************************/

template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
  void Morton_SFC_Output_Tecplot(Quad_Soln_Block             *Local_SolnBlk,
                                 Quad_Soln_Input_Parameters  &Input_Parameters,
                                 AdaptiveBlock2D_List        &Soln_Block_List) {
 
  int i, BLK, i_centre, j_centre;

  int num_of_blocks = 0;
  for ( BLK = 0; BLK <= Soln_Block_List.Nblk-1; BLK++){
    if (Soln_Block_List.Block[BLK].used == ADAPTIVEBLOCK2D_USED) num_of_blocks++;
  } /* endfor */

  ofstream Output_File;

  if (num_of_blocks > 0) { 
     char prefix[256], extension[256], output_file_name[256];

     //prefix
     i = 0;
     while (1) {
        if (Input_Parameters.Output_File_Name[i] == ' ' ||
            Input_Parameters.Output_File_Name[i] == '.') break;
        prefix[i]=Input_Parameters.Output_File_Name[i];
        i = i + 1;
        if (i > strlen(Input_Parameters.Output_File_Name) ) break;
     } /* endwhile */
     prefix[i] = '\0';
     strcat(prefix, "_Morton_SFC_cpu");

     //extension
     sprintf(extension, "%.6d", Soln_Block_List.ThisCPU);
     strcat(extension, ".dat");

     //complete file name
     strcpy(output_file_name, prefix);
     strcat(output_file_name, extension);
     char *Output_File_Name = output_file_name;
  
     //open output file
     Output_File.open(Output_File_Name,ios::out | ios::trunc);
     Output_File << setprecision(14);
     Output_File << "TITLE = \"Morton SFC \"" << "\n"<< "VARIABLES = \"X\", \"Y\", \n";

     //used for Sizing the 'Floor' for tecplot
     Vector2D V, minvec;
     //minvec =  Local_SolnBlk[0].Grid.Node[Local_SolnBlk[0].Grid.INl ][Local_SolnBlk[0].Grid.JNl].X; 
     Output_File << "ZONE T =\"Floor\", I=2, J=2, F=POINT"<<"\n"
	         <<minvec.x<<" "<<minvec.y<<" "<<minvec.x<<" "<<minvec.y<<" "
	         <<minvec.x<<" "<<minvec.y<<" "<<minvec.x<<" "<<minvec.y<<"\n";
  
     Output_File << "GEOMETRY X=0, Y=0, CS=GRID, C=BLUE, T=LINE, F=POINT"<<"\n"
	         <<"1"<<"\n"<<num_of_blocks<<"\n"; 
   
     for ( BLK = 0; BLK <= Soln_Block_List.Nblk-1; BLK++){
       if (Soln_Block_List.Block[BLK].used == ADAPTIVEBLOCK2D_USED) {
          i_centre = Local_SolnBlk[BLK].Grid.INl+(Local_SolnBlk[BLK].Grid.INu-Local_SolnBlk[BLK].Grid.INl)/2;
          j_centre = Local_SolnBlk[BLK].Grid.JNl+(Local_SolnBlk[BLK].Grid.JNu-Local_SolnBlk[BLK].Grid.JNl)/2;
          V = Local_SolnBlk[BLK].Grid.Node[i_centre][j_centre].X;
          Output_File << V.x <<" "<<V.y<<"\n"; //Output the centre node location for each used block
       } /* endif */
     } /* endfor */
     Output_File.close();
  } /* endif */
  
}

#endif /*_QUADTREE_MORTON_INCLUDED*/
 
