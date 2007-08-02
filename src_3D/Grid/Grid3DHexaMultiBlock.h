/* Grid3DHexaMultiBlock.h:  Header file defining 
                   3D Grid Multiblock. */

#ifndef _GRID3D_HEXA_MULTIBLOCK_INCLUDED
#define _GRID3D_HEXA_MULTIBLOCK_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include math macro, CFD, 3D vector, 3D cell*/

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif //_VECTOR3D_INCLUDED

#ifndef _CELL3D_INCLUDED
#include "../Grid/Cell3D.h"
#endif // _CELL3D_INCLUDED

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "Grid3DHexaBlock.h"
#endif // _GRID3D_HEXA_BLOCK_INCLUDED

#define GRID_INPUT_PARAMETER_LENGTH 128
// input parameters for grid generation
struct Grid3D_Input_Parameters{
  public:
   // number of blocks and cells in each direction
   int IBlk, JBlk, KBlk;
   int ICells, JCells, KCells;
   int Nghost;
   int i_Grid;
   char Grid_Type[GRID_INPUT_PARAMETER_LENGTH];
   // demensions for flow geometry
   double Box_Length, Box_Width, Box_Height;
   double I_stretching_factor, J_stretching_factor, 
      K_stretching_factor;
   int geometry_index;
   
   
  //constructor : set some default values
   Grid3D_Input_Parameters(void):
      IBlk(1), JBlk(1), KBlk(1),
      ICells(6), JCells(6), KCells(6),
      Nghost(2), Box_Length(ONE), Box_Width(ONE),
      Box_Height(ONE), i_Grid(GRID_CUBE), I_stretching_factor(0),
      J_stretching_factor(0), K_stretching_factor(0), 
      geometry_index(1){ }
   
   //destructor
   ~Grid3D_Input_Parameters(void){
      IBlk = 0; JBlk = 0; KBlk = 0;
      ICells = 0; JCells = 0; KCells = 0;
      Nghost = 0; geometry_index = 0;}
};

/* Define the HEXA 3D grid multiblock class. */

class Grid3D_Hexa_Multi_Block{
   
  public:
   Grid3D_Hexa_Block   ****Grid_ptr; // 3 dimensional array of pointers.
   int  NBI, NBJ, NBK; // Number of blocks in i, j and k directions;
   
   /* Creation constructors. */
   Grid3D_Hexa_Multi_Block(const Grid3D_Input_Parameters  &IP_Grid):
      NBI(IP_Grid.IBlk), 
      NBJ(IP_Grid.JBlk), NBK(IP_Grid.KBlk) {

      // allocate memory;
      Grid_ptr = new Grid3D_Hexa_Block***[NBI];
      for (int i = 0 ; i < NBI ; ++i ) {
         Grid_ptr[i] = new Grid3D_Hexa_Block**[NBJ];
         for (int j = 0 ; j <NBJ ; ++j ) {
            Grid_ptr[i][j] = new Grid3D_Hexa_Block*[NBK];
         }  /* endfor */
      }/* endfor */ 
      
      // create the pointers (pointing to the single grid block)
    
      int Nghost = IP_Grid.Nghost;
      int NNI = IP_Grid.ICells;
      int NNJ = IP_Grid.JCells;
      int NNK = IP_Grid.KCells;
   
      for (int i = 0 ; i <NBI ; ++i ) 
         for (int j = 0 ; j <NBJ ; ++j ) 
            for (int k = 0 ; k <NBK ; ++k ){
               
               Grid_ptr[i][j][k] = new Grid3D_Hexa_Block(
                  NNI,  NNJ , NNK, Nghost);
               
            }
     
      // create the multiblock grid
      switch(IP_Grid.i_Grid) {
         
      case GRID_CUBE :
         
         Grid_Cube(IP_Grid);
	 
         break;
         
      case GRID_CHANNEL :
         
         Grid_Channel(IP_Grid);
         break;

      case GRID_COUETTE :
         
         Grid_Couette(IP_Grid);
         
         break;
         
         
      } /* endswitch */

      
   } //end of creation constructor
   
   ~Grid3D_Hexa_Multi_Block(void) {
      
      Deallocate_Unused_Grid_Blocks();
      
      /* Deallocate memory. */
      for (int i = 0 ; i < NBI ; ++i ) {
         for ( int j=0 ; j <NBJ ; ++j ) {
            delete []Grid_ptr[i][j];
            Grid_ptr[i][j]=NULL;
         }  /* endfor */
         
         delete []Grid_ptr[i];
         Grid_ptr[i] = NULL;
      }  /* endfor */
      delete []Grid_ptr;
      Grid_ptr = NULL;


      NBI = 0;  NBJ=0;  NBK= 0;
         
  
   }


/* member functions  */
   void Grid_Cube(const Grid3D_Input_Parameters  &IP_Grid);
   void Grid_Channel(const Grid3D_Input_Parameters  &IP_Grid);
   void Grid_Couette(const Grid3D_Input_Parameters  &IP_Grid);
   void Deallocate_Unused_Grid_Blocks(void);
   void Write_Multi_Block_Grid(ostream &Out_File);
   void Output_Tecplot(ostream &Out_File);
   void Output_Nodes_Tecplot(ostream &Out_File);
   void Output_Cells_Tecplot(ostream &Out_File);
   void Output_Gnuplot(ostream &Out_File);
 
   
};

#endif // _GRID3D_HEXA_MULTIBLOCK_INCLUDED
