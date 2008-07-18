/* Grid3DHexaMultiBlock.h:  Header file defining 
                            3D hexahedral multiblock grids. */

#ifndef _GRID3D_HEXA_MULTIBLOCK_INCLUDED
#define _GRID3D_HEXA_MULTIBLOCK_INCLUDED

/* Include 3D grid/mesh input header file. */

#ifndef _GRID3D_INPUT_INCLUDED
#include "../../Reconstruction/Reconstruction3D/Reconstruct3DInput.h"
#endif // _GRID3D_INPUT_INCLUDED

/* Include 3D hexahedral block grid header file. */

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "Grid3DHexaBlock.h"
#endif // _GRID3D_HEXA_BLOCK_INCLUDED

/* Include adaptive block orientation information header file. */

//#ifndef _BLOCK_ORIENTATION_INFO_INCLUDED
//#include "../../../../src_3D/AMR/BlockOrientationInfo.h"
//#endif // _BLOCK_ORIENTATION_INFO_INCLUDED

//! 
#define GRID3D_HEXA_MULTI_BLOCK_NAX_NUMBER_OF_DIRECTIONS_TO_NEIGHBOURS   26

//! Maximum number of grid block neighbours in each direction
#define	GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS_PER_DIRECTION  4
#define	GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS  2

//! No grid block neighbour block  
#define	GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR   -1

/* Define the 3D hexahedral grid multiblock classes. */

//class Grid3D_Hexa_Multi_Block_Connectivity {
//  public:
//   //! Number of neigbouring grid blocks in each of the 26 directions
//   int      number_neighbours[GRID3D_HEXA_MULTI_BLOCK_NAX_NUMBER_OF_DIRECTIONS_TO_NEIGHBOURS];
//
//   //! Number of neighbouring grid blocks sharing each of the six (6) faces of the block 
//   int      num_neighT,num_neighB,num_neighN,num_neighS,num_neighE,num_neighW;
//
//   // ! Number of neighbouring grid blocks sharing each of the twelve (12) edges of the block
//   int      num_neighNW,num_neighNE,num_neighSE,num_neighSW,
//            num_neighTN,num_neighTS,num_neighTE,num_neighTW,
//            num_neighBN,num_neighBS,num_neighBE,num_neighBW;
//
//   // ! Number of neighbouring grid blocks sharing each of the eight (8) corners of the block
//   int      num_neighTNW,num_neighTSW,num_neighTNE,num_neighTSE,
//            num_neighBNW,num_neighBSW,num_neighBNE,num_neighBSE;
//
//   //! Grid block numbers for each of the neigbouring grid blocks in each of the 26 directions
//   int      neighbour_blknum[GRID3D_HEXA_MULTI_BLOCK_NAX_NUMBER_OF_DIRECTIONS_TO_NEIGHBOURS]
//                            [GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS_PER_DIRECTION];
//
//   //! Grid block numbers for neighbouring grid blocks sharing each of the six (6) faces of the block
//   int      neighT,neighB,neighN,neighS,neighE,neighW;
//
//   // ! Grid block numbers for neighbouring grid blocks sharing each of the twelve (12) edges of the block
//   int      neighNW[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighNE[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighSE[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighSW[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighTN[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighTS[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighTE[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighTW[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighBN[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighBS[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighBE[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighBW[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS];
//
//   // ! Grid block numbers for grid blocks sharing each of the eight (8) corners of the block
//   int      neighTNW[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighTSW[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighTNE[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighTSE[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighBNW[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighBSW[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighBNE[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//            neighBSE[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS];
//   
//   //! Compact transformation matrix for each of the neigbouring grid blocks in each of the 26 directions
//   int      neighbour_orient_info[GRID3D_HEXA_MULTI_BLOCK_NAX_NUMBER_OF_DIRECTIONS_TO_NEIGHBOURS]
//                                 [GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS_PER_DIRECTION];
//
//   //! Compact transformation matrix for neighbouring grid blocks sharing each of the six (6) faces of the block
//   Block_Orientation_Info  neighT_info,neighB_info,neighN_info,neighS_info,neighE_info,neighW_info; 
//
//   //! Compact transformation matrix for neighbouring grid blocks sharing each of the twelve (12) edges of the block
//   Block_Orientation_Info  neighNW_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighNE_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighSE_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighSW_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS], 
//                           neighTN_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighTS_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighTE_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighTW_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighBN_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighBS_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighBE_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighBW_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS];
//
//   //! Compact transformation matrix for grid blocks sharing each of the eight (8) corners of the block
//   Block_Orientation_Info  neighTNW_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighTSW_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighTNE_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighTSE_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighBNW_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighBSW_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighBNE_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS],
//                           neighBSE_info[GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS];
//
//   //! Compact transformation matrix for grid blocks sharing each of the eight (8) corners of the block
//   Block_Boundary_Elements_on_Domain_Extent be;
//   
//   //! Creation constructor
//   Grid3D_Hexa_Multi_Block_Connectivity(void) :
//      num_neighT(0), num_neighB(0), num_neighN(0), num_neighS(0), num_neighE(0), num_neighW(0),
//      num_neighNE(0), num_neighNW(0), num_neighSE(0), num_neighSW(0),
//      num_neighTN(0), num_neighTS(0), num_neighTE(0), num_neighTW(0),
//      num_neighBN(0), num_neighBS(0), num_neighBE(0), num_neighBW(0),
//      num_neighTNW(0), num_neighTSW(0), num_neighTNE(0), num_neighTSE(0),
//      num_neighBNW(0), num_neighBSW(0), num_neighBNE(0), num_neighBSE(0),
//      neighT(GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR), neighB(GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR),
//      neighN(GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR), neighS(GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR),
//      neighE(GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR), neighW(GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR) {
//      for (int i_neigh = 0; i_neigh < GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS; i_neigh++) {
//         neighNW[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR; 
//         neighNE[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR; 
//         neighSE[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR;
//         neighSW[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR; 
//         neighTN[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR;
//         neighTS[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR; 
//         neighTE[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR;
//         neighTW[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR; 
//         neighBN[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR;
//         neighBS[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR; 
//         neighBE[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR;
//         neighBW[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR; 
//         neighTNW[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR;
//         neighTSW[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR; 
//         neighTNE[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR;
//         neighTSE[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR;
//         neighBNW[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR;
//         neighBSW[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR; 
//         neighBNE[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR;
//         neighBSE[i_neigh] = GRID3D_HEXA_MULTI_BLOCK_NO_NEIGHBOUR; 
//      } /* endfor */
//   }
//
//   //! Broadcast
//   void Broadcast(void);
//};


//class Grid3D_Hexa_Multi_Block_List {
//  public:
//    Grid3D_Hexa_Block     *Grid_Blks; // one dimensional array of grid block.
//    int                    NBlk;
//    int                    NBlk_Idir, 
//                           NBlk_Jdir, 
//                           NBlk_Kdir; // Number of blocks in i, j and k directions.
//    int                    Allocated; // Indicates if the grid blocks have been allocated or not.
// 
//    // Grid block connectivity information  
//    Grid3D_Hexa_Multi_Block_Connectivity  *Connectivity; 
//   
//    /* Creation constructors. */
//    Grid3D_Hexa_Multi_Block_List(void) : 
//       NBlk(0), NBlk_Idir(0), NBlk_Jdir(0), NBlk_Kdir(0), Grid_Blks(NULL), 
//       Connectivity(NULL), Allocated(0) { }
//
//    Grid3D_Hexa_Multi_Block_List(const int N) {
//       Allocate(N);
//    }
//
//    Grid3D_Hexa_Multi_Block_List(const int Ni, 
//                                 const int Nj, 
//                                 const int Nk) {
//       Allocate(Ni, Nj, Nk);
//    }
//
//    /* Destructor. */
//    ~Grid3D_Hexa_Multi_Block_List(void) {
//        Deallocate();
//    }
//
//    /* Other member functions  */
//
//    void Allocate(const int Ni, const int Nj, const int Nk);
//
//    void Allocate(const int N);
//
//    void Deallocate(void);
//
//    void Copy(Grid3D_Hexa_Multi_Block_List &Grid2);
//
//    void Broadcast(void);
//
//    void Output(ostream &Out_File);
//
//    void Output_Tecplot(ostream &Out_File);
//
//    void Output_Nodes_Tecplot(ostream &Out_File);
//
//    void Output_Cells_Tecplot(ostream &Out_File);
//
//    void Output_Gnuplot(ostream &Out_File);
//
//    void Create_Grid(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Cube(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Periodic_Box(Reconstruct3D_Input_Parameters &Input);
//    
//    void Create_Grid_Flat_Plate(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Turbulence_Box(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Bunsen_Inflow(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Bunsen_Burner(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Bunsen_Box(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Channel(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Couette(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Pipe(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Bluff_Body_Burner(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Bump_Channel_Flow(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_ICEMCFD(Reconstruct3D_Input_Parameters &Input);
//
//    void Find_Neighbours(Reconstruct3D_Input_Parameters &Input);
//
//  private:
//    //copy and assignment are not permitted
//    Grid3D_Hexa_Multi_Block_List(const Grid3D_Hexa_Multi_Block_List &G);
//    Grid3D_Hexa_Multi_Block_List &operator = (const Grid3D_Hexa_Multi_Block_List &G);
//};

class Grid3D_Hexa_Multi_Block {
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

    void Create_Grid(Reconstruct3D_Input_Parameters &Input);

    void Create_Grid_Cube(Reconstruct3D_Input_Parameters &Input);

//    void Create_Grid_Channel(Reconstruct3D_Input_Parameters &Input);
//    
//    void Create_Grid_Flat_Plate(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Couette(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Pipe(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Bump_Channel_Flow(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Bluff_Body_Burner(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Turbulence_Box(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Bunsen_Inflow(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Bunsen_Burner(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_Bunsen_Box(Reconstruct3D_Input_Parameters &Input);
//
//    void Create_Grid_ICEMCFD(Reconstruct3D_Input_Parameters &Input);

  private:
    //copy and assignment are not permitted
    Grid3D_Hexa_Multi_Block(const Grid3D_Hexa_Multi_Block &G);
    Grid3D_Hexa_Multi_Block &operator = (const Grid3D_Hexa_Multi_Block &G);
};

#endif // _GRID3D_HEXA_MULTIBLOCK_INCLUDED
