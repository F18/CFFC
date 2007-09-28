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

#define Mesh_Orientation_Matrix int

/* Define the 3D hexahedral grid multiblock classes. */

class Grid3D_Hexa_Multi_Block_List{
  public:
    Grid3D_Hexa_Block     *Grid_Blks; // one dimensional array of grid block.
    int                    NBlk;
    int                    NBlk_Idir, 
                           NBlk_Jdir, 
                           NBlk_Kdir; // Number of blocks in i, j and k directions.
    int                    Allocated; // Indicates if the grid blocks have been allocated or not.
    int      *neighT,*neighB,*neighN; // Possible 26 neighbouring grid blocks in each direction
    int      *neighS,*neighE,*neighW; // for each of the block in the multi-block grid 
    int      *neighNW,*neighNE,*neighSE,*neighSW; 
    int      *neighTN,*neighTS,*neighTE,*neighTW;  
    int      *neighTNW,*neighTSW,*neighTNE,*neighTSE;  
    int      *neighBN,*neighBS,*neighBE,*neighBW;  
    int      *neighBNW,*neighBSW,*neighBNE,*neighBSE;
   
    Mesh_Orientation_Matrix  *neighT_ctm,*neighB_ctm,*neighN_ctm; // Compact transformation matrix for each of the
    Mesh_Orientation_Matrix  *neighS_ctm,*neighE_ctm,*neighW_ctm; // possible 26 neighbouring grid blocks in each direction
    Mesh_Orientation_Matrix  *neighNW_ctm,*neighNE_ctm,*neighSE_ctm,*neighSW_ctm; // for each of the block in the multi-block grid
    Mesh_Orientation_Matrix  *neighTN_ctm,*neighTS_ctm,*neighTE_ctm,*neighTW_ctm; 
    Mesh_Orientation_Matrix  *neighTNW_ctm,*neighTSW_ctm,*neighTNE_ctm,*neighTSE_ctm;  
    Mesh_Orientation_Matrix  *neighBN_ctm,*neighBS_ctm,*neighBE_ctm,*neighBW_ctm;  
    Mesh_Orientation_Matrix  *neighBNW_ctm,*neighBSW_ctm,*neighBNE_ctm,*neighBSE_ctm;

    /* Creation constructors. */
    Grid3D_Hexa_Multi_Block_List(void): 
      NBlk(0), NBlk_Idir(0), NBlk_Jdir(0), NBlk_Kdir(0), Grid_Blks(NULL), Allocated(0),
      neighT(NULL), neighB(NULL), neighN(NULL), neighS(NULL), neighE(NULL), neighW(NULL),
      neighNW(NULL), neighNE(NULL), neighSE(NULL), neighSW(NULL), 
      neighTN(NULL), neighTS(NULL), neighTE(NULL), neighTW(NULL),   
      neighTNW(NULL), neighTSW(NULL), neighTNE(NULL), neighTSE(NULL),   
      neighBN(NULL), neighBS(NULL), neighBE(NULL), neighBW(NULL),   
      neighBNW(NULL), neighBSW(NULL), neighBNE(NULL), neighBSE(NULL),
      neighT_ctm(NULL), neighB_ctm(NULL), neighN_ctm(NULL), neighS_ctm(NULL), neighE_ctm(NULL), neighW_ctm(NULL),
      neighNW_ctm(NULL), neighNE_ctm(NULL), neighSE_ctm(NULL), neighSW_ctm(NULL), 
      neighTN_ctm(NULL), neighTS_ctm(NULL), neighTE_ctm(NULL), neighTW_ctm(NULL),   
      neighTNW_ctm(NULL), neighTSW_ctm(NULL), neighTNE_ctm(NULL), neighTSE_ctm(NULL),   
      neighBN_ctm(NULL), neighBS_ctm(NULL), neighBE_ctm(NULL), neighBW_ctm(NULL),   
      neighBNW_ctm(NULL), neighBSW_ctm(NULL), neighBNE_ctm(NULL), neighBSE_ctm(NULL) { 
    }

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

    void Find_Neighbours(void);

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
