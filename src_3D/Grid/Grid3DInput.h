/* Grid3DInput.h:  Header file defining 3D grid/mesh input data class. */

#ifndef _GRID3D_INPUT_INCLUDED
#define _GRID3D_INPUT_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include required CFFC header files. */

#ifndef _CFD_INCLUDE
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif //_VECTOR3D_INCLUDED

#ifndef _ICEMCFD_INCLUDED
#include "../ICEMCFD/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

//#ifndef _GRID3D_HO_EXECUTIONMODE_INCLUDED
//#include "Grid3DHighOrderExecutionMode.h"
//#endif //_GRID3D_HO_EXECUTIONMODE_INCLUDED

#define GRID_INPUT_PARAMETER_LENGTH 256

/*!
 * Class: Grid3D_Input_Parameters
 *
 * @brief Input Parameters for 3D hexahedral mesh generation.
 *
 * This class defines and handles the input variables related to the
 * the 3D hexahedral mesh selection and grid generation.
 *
 */
class Grid3D_Input_Parameters{
  public:
    //@{ @name Number of blocks in multiblock mesh in each direction:
    int NBlk_Idir, NBlk_Jdir, NBlk_Kdir;
    //@}

    //@{ @name Default number of cells in multiblock mesh in each direction:
    int NCells_Idir, NCells_Jdir, NCells_Kdir;
    int Nghost;
    //@}

    //@{ @name Default number of cells in turbulence multiblock mesh in each direction:
    int NCells_Turbulence_Idir, NCells_Turbulence_Jdir, NCells_Turbulence_Kdir;
    //@}

    //@{ @name Grid type indicator:
    int i_Grid;
    char Grid_Type[GRID_INPUT_PARAMETER_LENGTH];
    //@}

    //@{ @name Multi-block mesh input file:
    char Grid_File_Name[GRID_INPUT_PARAMETER_LENGTH];
    //@}

    //@{ @name ICEMCFD Filenames:
    char **ICEMCFD_FileNames;
    //@}

    //@{ @name Dimensions for basic flow geometry:
    double Box_Length, Box_Width, Box_Height;
    //@}

    //@{ @name Mesh stretching parameters:
    int Mesh_Stretching;
    int Stretching_Type_Idir, Stretching_Type_Jdir, 
        Stretching_Type_Kdir; 
    double Stretching_Factor_Idir, Stretching_Factor_Jdir, 
           Stretching_Factor_Kdir;
    int Mesh_Smoothing;
    //@}

    //@{ @name Mesh shifting, scaling and rotation parameters:
    Vector3D X_Shift;
    double X_Scale, X_Rotate;
    //@}

    //@{ @name Mesh distortion:
    int Disturb_Interior_Nodes;
    //@}
    
    //@{ @name Flat Plate mesh parameters:
    double Plate_Length;
    //@}
    
    //@{ @name Pipe mesh parameters:
    double Pipe_Length, Pipe_Radius;
    //@}

    //@{ @name Bluff body burner mesh parameters:
    double Radius_Fuel_Line, Radius_Bluff_Body, Radius_Coflow_Inlet_Pipe,
           Length_Coflow_Inlet_Pipe, Length_Combustor_Tube; 
    //@}

    //@{ @name Bunsen burner mesh parameters:
    double Radius_Bunsen_Burner_Fuel_Line, Radius_Bunsen_Burner, Height_Bunsen_Burner;
    //@}

    //@{ @name Dimensions for slot burner mesh parameters:
    double Slot_Width;
    //@}

    //@{ @name Dimensions for turbulence box:
    double Turbulence_Box_Length, Turbulence_Box_Width, Turbulence_Box_Height;
    //@}

    //@{ @name Reconstruction type indicator:
    char Reconstruction_Type_For_Grid_Info[GRID_INPUT_PARAMETER_LENGTH];
    //@}

    //@{ @name Constructors and desctructors:
    //! Constructor (assign default values)
    Grid3D_Input_Parameters(void){
       // Basic grid parameters:
       NBlk_Idir = 1; NBlk_Jdir = 1; NBlk_Kdir = 1;
       NCells_Idir = 10; NCells_Jdir = 10; NCells_Kdir = 10; 
       Nghost = 2;
       NCells_Turbulence_Idir = NCells_Idir;  
       NCells_Turbulence_Jdir = NCells_Jdir;
       NCells_Turbulence_Kdir = NCells_Kdir;
       Box_Length = ONE; Box_Width = ONE; Box_Height = ONE;
       strcpy(Grid_Type,"Cube");
       i_Grid = GRID_CUBE; 
       strcpy(Grid_File_Name,"gridfile.grid");
       Mesh_Stretching = OFF;
       Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
       Stretching_Type_Jdir = STRETCHING_FCN_LINEAR; 
       Stretching_Type_Kdir = STRETCHING_FCN_LINEAR;
       Stretching_Factor_Idir = 1.10;
       Stretching_Factor_Jdir = 1.10;
       Stretching_Factor_Kdir = 1.10;
       Mesh_Smoothing = OFF;
       X_Shift = Vector3D_ZERO;
       X_Scale = ONE;
       X_Rotate = ZERO;
       Disturb_Interior_Nodes = OFF;
       // Pipe parameters:
       Pipe_Length = ONE; Pipe_Radius = 0.1234;
       // Flat plate parameters:
       Plate_Length = Box_Width/2;
       // Bluff body burner parameters:
       Radius_Fuel_Line = 0.0018;
       Radius_Bluff_Body = 0.025;
       Radius_Coflow_Inlet_Pipe = 0.07;
       Length_Coflow_Inlet_Pipe = 0.1; 
       Length_Combustor_Tube = 0.3;
       // Bunsen burner parameters:
       Radius_Bunsen_Burner_Fuel_Line = 0.0056;
       Radius_Bunsen_Burner = 0.025;
       Height_Bunsen_Burner = 0.075;
       // Slot burner parameters:
       Slot_Width = 0.025;
       // Turbulence box parameters:
       Turbulence_Box_Length = Box_Length; 
       Turbulence_Box_Width = Box_Width; 
       Turbulence_Box_Height = Box_Height;
       //ICEM Filenames:
       ICEMCFD_FileNames = ICEMCFD_get_filenames();
       // Set default values in class Grid3D_HO_Execution_Mode
       // Grid3D_HO_Execution_Mode::SetDefaults();
    }
   
    //! Destructor
    ~Grid3D_Input_Parameters(void){}
    //@}
 
    //@{ @name Other Member functions:
    //! Broadcast input parameters to all processors
    void Broadcast(void);
    //! Parse next input line
    int Parse_Next_Input_Control_Parameter(char *code, stringstream &value);
    //! Check validity of specified input parameters
    int Check_Inputs(void);
    //@}


    //@{ @name Input-output operators:
    friend ostream &operator << (ostream &out_file,
                                 const Grid3D_Input_Parameters &IP);
    friend istream &operator >> (istream &in_file,
                                 Grid3D_Input_Parameters &IP);
    void Output(ostream &out_file) const;
    //@}

};

#endif // _GRID3D_INPUT_INCLUDED
