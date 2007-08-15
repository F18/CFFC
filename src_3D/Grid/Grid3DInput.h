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

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

#define GRID_INPUT_PARAMETER_LENGTH 128

/* Define the 3D hexahedral grid input class for controlling
   mesh selection and grid generation. */

struct Grid3D_Input_Parameters{
  public:
    // Number of blocks in multiblock mesh in each direction
    int NBlk_Idir, NBlk_Jdir, NBlk_Kdir;

    // Default number of cells in multiblock mesh in each direction
    int NCells_Idir, NCells_Jdir, NCells_Kdir;
    int Nghost;

    // Grid type indicator
    int i_Grid;
    char Grid_Type[GRID_INPUT_PARAMETER_LENGTH];

    // ICEMCFD Filenames
    char **ICEMCFD_FileNames;

    // Dimensions for basic flow geometry
    double Box_Length, Box_Width, Box_Height;

    // Mesh stretching parameters
    int Stretching_Type_Idir, Stretching_Type_Jdir, 
        Stretching_Type_Kdir; 
    double Stretching_Factor_Idir, Stretching_Factor_Jdir, 
           Stretching_Factor_Kdir;

   // Pipe mesh parameters
   double Pipe_Length, Pipe_Radius;

   // Bluff body burner mesh parameters
   double Radius_Fuel_Line, Radius_Bluff_Body, Radius_Coflow_Inlet_Pipe,
          Length_Coflow_Inlet_Pipe, Length_Combustor_Tube; 

   // Constructor: set some default values
   Grid3D_Input_Parameters(void){
       NBlk_Idir = 1; NBlk_Jdir = 1; NBlk_Kdir = 1;
       NCells_Idir = 10; NCells_Jdir = 10; NCells_Kdir = 10; 
       Nghost = 2;
       Box_Length = ONE; Box_Width = ONE; Box_Height = ONE;
       i_Grid = GRID_CUBE; 
       ICEMCFD_FileNames = ICEMCFD_get_filenames();
       Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
       Stretching_Type_Jdir = STRETCHING_FCN_LINEAR; 
       Stretching_Type_Kdir = STRETCHING_FCN_LINEAR;
       Stretching_Factor_Idir = 1.10;
       Stretching_Factor_Jdir = 1.10;
       Stretching_Factor_Kdir = 1.10;
       Pipe_Length = ONE; Pipe_Radius = 0.1;
       Radius_Fuel_Line = 1.80e-03;
       Radius_Bluff_Body = 25.04e-03;
       Radius_Coflow_Inlet_Pipe = 0.1;
       Length_Coflow_Inlet_Pipe = 0.127; 
       Length_Combustor_Tube = 0.508;
   }
   
   // Use default destructor
   //~Grid3D_Input_Parameters(void){
   //}

};

#endif // _GRID3D_INPUT_INCLUDED
