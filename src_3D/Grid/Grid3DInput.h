/* Grid3DInput.h:  Header file defining 
                   3D grid/mesh input data class. */

#ifndef _GRID3D_INPUT_INCLUDED
#define _GRID3D_INPUT_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include required CFFC header files. */

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#define GRID_INPUT_PARAMETER_LENGTH 128

/* Define the 3D hexahedral grid input class for controlling
   mesh selection and grid generation. */

struct Grid3D_Input_Parameters{
  public:
    // number of blocks and cells in each direction
    int NBlk_Idir, NBlk_Jdir, NBlk_Kdir;
    int NCells_Idir, NCells_Jdir, NCells_Kdir;
    int Nghost;

    // Grid type
    int i_Grid;
    char Grid_Type[GRID_INPUT_PARAMETER_LENGTH];

    // demensions for flow geometry
    double Box_Length, Box_Width, Box_Height;
    int Stretching_Type_Idir, Stretching_Type_Jdir, 
        Stretching_Type_Kdir; 
    double Stretching_Factor_Idir, Stretching_Factor_Jdir, 
           Stretching_Factor_Kdir;
   
   //constructor : set some default values
   Grid3D_Input_Parameters(void){
       NBlk_Idir = 1; NBlk_Jdir = 1; NBlk_Kdir = 1;
       NCells_Idir = 10; NCells_Jdir = 10; NCells_Kdir = 10; 
       Nghost = 2;
       Box_Length = ONE; Box_Width = ONE; Box_Height = ONE;
       i_Grid = GRID_CUBE; 
       Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
       Stretching_Type_Jdir = STRETCHING_FCN_LINEAR; 
       Stretching_Type_Kdir = STRETCHING_FCN_LINEAR;
       Stretching_Factor_Idir = 1.01;
       Stretching_Factor_Jdir = 1.01;
       Stretching_Factor_Kdir = 1.01;
   }
   
   //destructor
   ~Grid3D_Input_Parameters(void){
      NBlk_Idir = 0; NBlk_Jdir = 0; NBlk_Kdir = 0;
      NCells_Idir = 0; NCells_Jdir = 0; NCells_Kdir = 0;
      Nghost = 0;
   }

};

#endif // _GRID3D_INPUT_INCLUDED
