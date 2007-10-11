/* Grid3DHexaMultiBlock.h:  Header file defining 
                            3D hexahedral multiblock grids. */

#ifndef _GRID3D_HEXA_MULTIBLOCK_INCLUDED
#define _GRID3D_HEXA_MULTIBLOCK_INCLUDED

/* Include 3D grid/mesh input header file. */

#ifndef _GRID3D_INPUT_INCLUDED
#include "Grid3DInput.h"
#endif // _GRID3D_INPUT_INCLUDED

/* Include 3D hexahedral block grid header file. */

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "Grid3DHexaBlock.h"
#endif // _GRID3D_HEXA_BLOCK_INCLUDED


#ifndef _BLOCK_ORIENTATION_INFO_INCLUDED
#include "../AMR/BlockOrientationInfo.h"
#endif // _BLOCK_ORIENTATION_INFO_INCLUDED


#define	GRID3D_NO_NEIGHBOUR             -1


/* Define the 3D hexahedral grid multiblock classes. */
class Grid3D_Hexa_Multi_Block_Connectivity{
   
  public:
   
   int      num_neighT,num_neighB,num_neighN; // Number of neighbouring grid blocks in each direction
   int      num_neighS,num_neighE,num_neighW; // for each of the block in the multi-block grid
   int      num_neighNW,num_neighNE,num_neighSE,num_neighSW;
   int      num_neighTN,num_neighTS,num_neighTE,num_neighTW;
   int      num_neighTNW,num_neighTSW,num_neighTNE,num_neighTSE;
   int      num_neighBN,num_neighBS,num_neighBE,num_neighBW;
   int      num_neighBNW,num_neighBSW,num_neighBNE,num_neighBSE;
   

   int      neighT,neighB,neighN; // Possible 26 neighbouring grid blocks in each direction
   int      neighS,neighE,neighW; // for each of the block in the multi-block grid
   int      neighNW[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighNE[BLOCK_ORIENTATION_MAX_NEIGHBOUR],
            neighSE[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighSW[BLOCK_ORIENTATION_MAX_NEIGHBOUR];
   int      neighTN[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighTS[BLOCK_ORIENTATION_MAX_NEIGHBOUR],
            neighTE[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighTW[BLOCK_ORIENTATION_MAX_NEIGHBOUR];
   int      neighTNW[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighTSW[BLOCK_ORIENTATION_MAX_NEIGHBOUR],
            neighTNE[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighTSE[BLOCK_ORIENTATION_MAX_NEIGHBOUR];
   int      neighBN[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighBS[BLOCK_ORIENTATION_MAX_NEIGHBOUR],
            neighBE[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighBW[BLOCK_ORIENTATION_MAX_NEIGHBOUR];
   int      neighBNW[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighBSW[BLOCK_ORIENTATION_MAX_NEIGHBOUR],
            neighBNE[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighBSE[BLOCK_ORIENTATION_MAX_NEIGHBOUR];
   

   Block_Orientation_Info  neighT_info,neighB_info,neighN_info; // Compact transformation matrix for each of the
   Block_Orientation_Info  neighS_info,neighE_info,neighW_info; // possible 26 neighbouring grid blocks in each direction
   Block_Orientation_Info  neighNW_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighNE_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],
                           neighSE_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighSW_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR]; // for each of the block in the multi-block grid
   Block_Orientation_Info  neighTN_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighTS_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],
                           neighTE_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighTW_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR];
   Block_Orientation_Info  neighTNW_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighTSW_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],
                           neighTNE_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighTSE_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR];
   Block_Orientation_Info  neighBN_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighBS_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],
                           neighBE_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighBW_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR];
   Block_Orientation_Info  neighBNW_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighBSW_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],
                           neighBNE_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR],neighBSE_info[BLOCK_ORIENTATION_MAX_NEIGHBOUR];


  //Creation constructors.
   Grid3D_Hexa_Multi_Block_Connectivity(void):
      neighT(GRID3D_NO_NEIGHBOUR), neighB(GRID3D_NO_NEIGHBOUR),
      neighN(GRID3D_NO_NEIGHBOUR), neighS(GRID3D_NO_NEIGHBOUR),
      neighE(GRID3D_NO_NEIGHBOUR), neighW(GRID3D_NO_NEIGHBOUR){
      
      for(int nN = 0; nN < BLOCK_ORIENTATION_MAX_NEIGHBOUR; nN++){
         neighNW[nN] = 0;  neighNE[nN] = 0; 
         neighSE[nN] = 0;  neighSW[nN] = 0; 
         neighTN[nN] = 0;  neighTS[nN] = 0; 
         neighTE[nN] = 0;  neighTW[nN] = 0; 
         neighTNW[nN] = 0;  neighTSW[nN] = 0; 
         neighTNE[nN] = 0;  neighTSE[nN] = 0; 
         neighBN[nN] = 0;  neighBS[nN] = 0; 
         neighBE[nN] = 0;  neighBW[nN] = 0; 
         neighBNW[nN] = 0;  neighBSW[nN] = 0; 
         neighBNE[nN] = 0;  neighBSE[nN] = 0; 
      }
         
      
   }

      
      
   void Deallocate(void);


};


class Grid3D_Hexa_Multi_Block_List{
  public:
    Grid3D_Hexa_Block     *Grid_Blks; // one dimensional array of grid block.
    int                    NBlk;
    int                    NBlk_Idir, 
                           NBlk_Jdir, 
                           NBlk_Kdir; // Number of blocks in i, j and k directions.
    int                    Allocated; // Indicates if the grid blocks have been allocated or not.
   
    Grid3D_Hexa_Multi_Block_Connectivity  *Connectivity; 
   
    /* Creation constructors. */
    Grid3D_Hexa_Multi_Block_List(void): 
       NBlk(0), NBlk_Idir(0), NBlk_Jdir(0), NBlk_Kdir(0), Grid_Blks(NULL), 
       Connectivity(NULL), Allocated(0){ }

    Grid3D_Hexa_Multi_Block_List(const int N) {
       Allocate(N);
    }

    Grid3D_Hexa_Multi_Block_List(const int Ni, 
                                 const int Nj, 
                                 const int Nk) {
       Allocate(Ni, Nj, Nk);
    }

    Grid3D_Hexa_Multi_Block_List(Grid3D_Input_Parameters &Input){
      // create various multiblock multiblock grid depending on input parameters
      switch(Input.i_Grid) {
        case GRID_CUBE :
          Create_Grid_Cube(Input);
          break;
        case GRID_CHANNEL_XDIR :
        case GRID_CHANNEL_YDIR:
        case GRID_CHANNEL_ZDIR:
          Create_Grid_Channel(Input);
          break;
        case GRID_COUETTE_XDIR :
        case GRID_COUETTE_YDIR:
        case GRID_COUETTE_ZDIR:
          Create_Grid_Couette(Input);
          break;
        case GRID_PIPE :
          Create_Grid_Pipe(Input);
          break;
        case GRID_BLUFF_BODY_BURNER :
          Create_Grid_Bluff_Body_Burner(Input);
          break;
        case GRID_ICEMCFD :
          //Create_Grid_ICEMCFD(Input);
          break;
        default:
          Create_Grid_Cube(Input);
          break;
      } /* endswitch */
    } 
   
    /* Destructor. */
    ~Grid3D_Hexa_Multi_Block_List(void) {
        Deallocate();
    }

    /* Other member functions  */

    void Allocate(const int Ni, const int Nj, const int Nk);

    void Allocate(const int N);

    void Deallocate(void);

    void Copy(Grid3D_Hexa_Multi_Block_List &Grid2);

    void Broadcast(void);

    void Output(ostream &Out_File);

    void Output_Tecplot(ostream &Out_File);

    void Output_Nodes_Tecplot(ostream &Out_File);

    void Output_Cells_Tecplot(ostream &Out_File);

    void Output_Gnuplot(ostream &Out_File);

    void Create_Grid(Grid3D_Input_Parameters &Input);

    void Create_Grid_Cube(Grid3D_Input_Parameters &Input);

    void Create_Grid_Channel(Grid3D_Input_Parameters &Input);

    void Create_Grid_Couette(Grid3D_Input_Parameters &Input);

    void Create_Grid_Pipe(Grid3D_Input_Parameters &Input);

    void Create_Grid_Bluff_Body_Burner(Grid3D_Input_Parameters &Input);

    void Find_Neighbours(Grid3D_Input_Parameters &Input);

    //void Create_Grid_ICEMCFD(Grid3D_Input_Parameters &Input);

  private:
    //copy and assignment are not permitted
    Grid3D_Hexa_Multi_Block_List(const Grid3D_Hexa_Multi_Block_List &G);
    Grid3D_Hexa_Multi_Block_List &operator = (const Grid3D_Hexa_Multi_Block_List &G);
   
};

class Grid3D_Hexa_Multi_Block{
  public:
    Grid3D_Hexa_Block   ***Grid_Blks; // 3 dimensional array of grid block.
    int                    NBlk_Idir, 
                           NBlk_Jdir, 
                           NBlk_Kdir; // Number of blocks in i, j and k directions.
    int                    Allocated; // Indicates if the grid blocks have been allocated or not.
   
    /* Creation constructors. */
    Grid3D_Hexa_Multi_Block(void){
       NBlk_Idir = 0; NBlk_Jdir = 0; NBlk_Kdir = 0;
       Grid_Blks = NULL; Allocated = 0;
    }

    Grid3D_Hexa_Multi_Block(const int Ni, 
                            const int Nj, 
                            const int Nk) {
       Allocate(Ni, Nj, Nk);
    }

    Grid3D_Hexa_Multi_Block(Grid3D_Input_Parameters &Input){
      // create various multiblock multiblock grid depending on input parameters
      switch(Input.i_Grid) {
        case GRID_CUBE :
          Create_Grid_Cube(Input);
          break;
        case GRID_CHANNEL_XDIR :
        case GRID_CHANNEL_YDIR:
        case GRID_CHANNEL_ZDIR:
          Create_Grid_Channel(Input);
          break;
        case GRID_COUETTE_XDIR :
        case GRID_COUETTE_YDIR:
        case GRID_COUETTE_ZDIR:
          Create_Grid_Couette(Input);
          break;
        case GRID_PIPE :
          Create_Grid_Pipe(Input);
          break;
        case GRID_BLUFF_BODY_BURNER :
          Create_Grid_Bluff_Body_Burner(Input);
          break;
        case GRID_ICEMCFD :
          Create_Grid_ICEMCFD(Input);
          break;
        default:
          Create_Grid_Cube(Input);
          break;
      } /* endswitch */

    
    } 
   
    /* Destructor. */
    ~Grid3D_Hexa_Multi_Block(void) {
        Deallocate();
    }

    /* Other member functions  */

    void Allocate(const int Ni, const int Nj, const int Nk);

    void Deallocate(void);

    void Copy(Grid3D_Hexa_Multi_Block &Grid2);

    void Broadcast(void);

    void Output(ostream &Out_File);

    void Output_Tecplot(ostream &Out_File);

    void Output_Nodes_Tecplot(ostream &Out_File);

    void Output_Cells_Tecplot(ostream &Out_File);

    void Output_Gnuplot(ostream &Out_File);

    void Create_Grid(Grid3D_Input_Parameters &Input);

    void Create_Grid_Cube(Grid3D_Input_Parameters &Input);

    void Create_Grid_Channel(Grid3D_Input_Parameters &Input);

    void Create_Grid_Couette(Grid3D_Input_Parameters &Input);

    void Create_Grid_Pipe(Grid3D_Input_Parameters &Input);

    void Create_Grid_Bluff_Body_Burner(Grid3D_Input_Parameters &Input);

    void Create_Grid_ICEMCFD(Grid3D_Input_Parameters &Input);

  private:
    //copy and assignment are not permitted
    Grid3D_Hexa_Multi_Block(const Grid3D_Hexa_Multi_Block &G);
    Grid3D_Hexa_Multi_Block &operator = (const Grid3D_Hexa_Multi_Block &G);
   
};

#endif // _GRID3D_HEXA_MULTIBLOCK_INCLUDED
