/* Grid3DInput.cc: Definition of Grid3D_Input_Parameters class member functions. */

/* Include the Grid3DInput header file. */

#ifndef _GRID3D_INPUT_INCLUDED
#include "Grid3DInput.h"
#endif // _GRID3D_INPUT_INCLUDED

/* Define member functions. */

/***************************************************************************
 * Grid3D_Input_Parameters::Broadcast -- Broadcast to all processors.      *
 ***************************************************************************/
void Grid3D_Input_Parameters::Broadcast(void) {

#ifdef _MPI_VERSION
   // Basic Grid parameters:
   MPI::COMM_WORLD.Bcast(Grid_Type,
                         GRID_INPUT_PARAMETER_LENGTH,
                         MPI::CHAR, 0);
   MPI::COMM_WORLD.Bcast(&(i_Grid),
                         1,
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(Grid_File_Name,
                         GRID_INPUT_PARAMETER_LENGTH,
                         MPI::CHAR, 0);
   MPI::COMM_WORLD.Bcast(&(NBlk_Idir),
                         1,
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(NBlk_Jdir),
                         1,
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(NBlk_Kdir),
                         1,
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(NCells_Idir),
                         1,
                         MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(NCells_Jdir),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(NCells_Kdir),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Nghost),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(NCells_Turbulence_Idir),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(NCells_Turbulence_Jdir),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(NCells_Turbulence_Kdir),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Stretching_Type_Idir),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Mesh_Stretching),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Stretching_Type_Jdir),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Stretching_Type_Kdir),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Stretching_Factor_Idir),
                          1,
                          MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(Stretching_Factor_Jdir),
                          1,
                          MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(Stretching_Factor_Kdir),
                          1,
                          MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(Mesh_Smoothing),
                          1,
                          MPI::INT,0);
   MPI::COMM_WORLD.Bcast(&(X_Shift.x), 
                         1, 
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(X_Shift.y), 
                         1, 
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(X_Shift.z), 
                         1, 
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(X_Scale),
                         1,
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(X_Rotate),
                         1,
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Disturb_Interior_Nodes),
                         1,
                         MPI::INT,0);

   // Grid cube and box dimensions:
   MPI::COMM_WORLD.Bcast(&(Box_Length),
                         1,
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Box_Width),
                         1,
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Box_Height),
                         1,
                         MPI::DOUBLE, 0);

   // Pipe parameters:
   MPI::COMM_WORLD.Bcast(&(Pipe_Radius),
                         1,
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Pipe_Length),
                         1,
                         MPI::DOUBLE, 0);

   // Bluff body burner parameters:
   MPI::COMM_WORLD.Bcast(&(Radius_Fuel_Line),
                         1,
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Radius_Bluff_Body),
                          1,
                          MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Radius_Coflow_Inlet_Pipe),
                         1,
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Length_Coflow_Inlet_Pipe),
                         1,
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Length_Combustor_Tube),
                         1,
                         MPI::DOUBLE, 0);

   // Bunsen burner parameters:
   MPI::COMM_WORLD.Bcast(&(Radius_Bunsen_Burner_Fuel_Line),
                          1,
                          MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Radius_Bunsen_Burner),
                          1,
                          MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Height_Bunsen_Burner),
                          1,
                          MPI::DOUBLE, 0);

   // Slot burner parameters:
   MPI::COMM_WORLD.Bcast(&(Slot_Width),
                          1,
                          MPI::DOUBLE, 0);

   // Turbulence box parameters:
   MPI::COMM_WORLD.Bcast(&(Turbulence_Box_Length),
                         1,
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Turbulence_Box_Width),
                         1,
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(Turbulence_Box_Height),
                         1,
                         MPI::DOUBLE, 0);

   //ICEM Filenames:
   MPI::COMM_WORLD.Bcast(ICEMCFD_FileNames[0],
                         20,
                         MPI::CHAR, 0);
   MPI::COMM_WORLD.Bcast(ICEMCFD_FileNames[1],
                         20,
                         MPI::CHAR, 0);
   MPI::COMM_WORLD.Bcast(ICEMCFD_FileNames[2],
                         20,
                         MPI::CHAR, 0);
#endif

}


/*********************************************************************************
 * Grid3D_Input_Parameters::Parse_Next_Input_Control_Parameter - Parse input.    *
 *********************************************************************************/
int Grid3D_Input_Parameters::Parse_Next_Input_Control_Parameter(char *code, 
                                                                stringstream &value) {

// Returns:
//  - INVALID_INPUT_VALUE if code is valid but value is invalid
//  - INVALID_INPUT_CODE  if unknown code

  int i_command = INVALID_INPUT_CODE;
  string value_string;

  /* Parse next control parameter with Grid3D_HO_Execution_Mode parser */
  //Grid3D_HO_Execution_Mode::Parse_Next_Input_Control_Parameter(*this, i_command, code, value);

  if (strcmp(code, "Grid_Type") == 0) {
     i_command = 3001;
     value >> value_string;
     strcpy(Grid_Type, value_string.c_str());

     if (strcmp(Grid_Type, "Cartesian") == 0) {
        i_Grid = GRID_CARTESIAN_UNIFORM;
        Box_Length = ONE;
        Box_Width = ONE;
        Box_Height = ONE;

     } else if (strcmp(Grid_Type, "Cube") == 0) {
        i_Grid = GRID_CUBE;
        Box_Length = ONE;
        Box_Width = ONE;
        Box_Height = ONE;
         
     } else if (strcmp(Grid_Type, "Flat_Plate") == 0) {
         i_Grid = GRID_FLAT_PLATE;
         Box_Length = HALF*Plate_Length;
         Box_Width = TWO*Plate_Length;
         Box_Height = Plate_Length;

     } else if (strcmp(Grid_Type, "Periodic_Box") == 0) {
        i_Grid = GRID_PERIODIC_BOX;
        Box_Length = ONE;
        Box_Width = ONE;
        Box_Height = ONE;

     } else if (strcmp(Grid_Type, "Periodic_Box_With_Inflow") == 0) {
        i_Grid = GRID_PERIODIC_BOX_WITH_INFLOW;
        Box_Length = ONE;
        Box_Width = ONE;
        Box_Height = ONE;

     } else if (strcmp(Grid_Type, "Bunsen_Burner") == 0) {
        i_Grid = GRID_BUNSEN_BURNER;
        Radius_Bunsen_Burner_Fuel_Line = 0.0056;
        Radius_Bunsen_Burner = 0.025;
        Height_Bunsen_Burner = 0.075;
	Turbulence_Box_Length = 0.2;
        Turbulence_Box_Width = 0.14;
        Turbulence_Box_Height = 0.14;

     } else if (strcmp(Grid_Type, "Bunsen_Box") == 0) {
        i_Grid = GRID_BUNSEN_BOX;
        Box_Length = 0.2;
        Box_Width = 0.2;
        Box_Height = 0.2;
	Slot_Width = 0.025;
	Turbulence_Box_Length = 0.2;
        Turbulence_Box_Width = 0.14;
        Turbulence_Box_Height = 0.14;

     } else if (strcmp(Grid_Type, "Bunsen_Inflow") == 0) {
        i_Grid = GRID_BUNSEN_INFLOW;
        Radius_Bunsen_Burner_Fuel_Line = 0.0056;
        Radius_Bunsen_Burner = 0.07;
        Height_Bunsen_Burner = 0.2;
	Turbulence_Box_Length = 0.2;
        Turbulence_Box_Width = 0.14;
        Turbulence_Box_Height = 0.14;

     } else if (strcmp(Grid_Type, "Turbulence_Box") == 0) {
        i_Grid = GRID_TURBULENCE_BOX;
        Turbulence_Box_Length = ONE;
        Turbulence_Box_Width = ONE;
        Turbulence_Box_Height = ONE;

     } else if (strcmp(Grid_Type, "Channel") == 0) {
        i_Grid = GRID_CHANNEL_ZDIR;
        Pipe_Length = 0.1;
        Pipe_Radius = 0.2;

     } else if (strcmp(Grid_Type, "Channel_X") == 0) {
        i_Grid = GRID_CHANNEL_XDIR;
        Box_Length = 0.001;
        Box_Width = 0.2;
        Box_Height = 0.001;

     } else if (strcmp(Grid_Type, "Channel_Y") == 0) {
        i_Grid = GRID_CHANNEL_YDIR;
        Box_Length = 0.001;
        Box_Width = 0.001;
        Box_Height = 0.2;

     } else if (strcmp(Grid_Type, "Channel_Z") == 0) {
        i_Grid = GRID_CHANNEL_ZDIR;
        Box_Length = 0.2;
        Box_Width = 0.001;
        Box_Height = 0.001;

     } else if (strcmp(Grid_Type, "Couette") == 0) {
        i_Grid = GRID_COUETTE_ZDIR;
        Box_Length = 0.2;
        Box_Width  = 0.001;
        Box_Height = 0.001;

     } else if (strcmp(Grid_Type, "Couette_X") == 0) {
        i_Grid = GRID_COUETTE_XDIR;
        Box_Length = 0.001;
        Box_Width = 0.2;
        Box_Height = 0.001;

     } else if (strcmp(Grid_Type, "Couette_Y") == 0) {
        i_Grid = GRID_COUETTE_YDIR;
        Box_Length = 0.001;
        Box_Width = 0.001;
        Box_Height = 0.2;

     } else if (strcmp(Grid_Type, "Couette_Z") == 0) {
        i_Grid = GRID_COUETTE_ZDIR;
        Box_Length = 0.2;
        Box_Width = 0.001;
        Box_Height = 0.001;

     } else if (strcmp(Grid_Type, "Turbulent_Channel") == 0) {
        i_Grid = GRID_CHANNEL_ZDIR;
        Box_Length = 1.524;
        Box_Width  = 0.127;
        Box_Height = 0.127;

     } else if (strcmp(Grid_Type, "Pipe") == 0) {
        i_Grid = GRID_PIPE;

     } else if (strcmp(Grid_Type, "Bump_Channel_Flow") == 0) {
       i_Grid = GRID_BUMP_CHANNEL_FLOW;

     } else if (strcmp(Grid_Type, "Bluff_Body_Burner") == 0) {
        i_Grid = GRID_BLUFF_BODY_BURNER;

     } else if (strcmp(Grid_Type, "ICEMCFD") == 0) {
        i_Grid = GRID_ICEMCFD;
        ICEMCFD_FileNames = ICEMCFD_get_filenames();

     } else {
        i_command = INVALID_INPUT_VALUE;
     } /* endif */

  } else if (strcmp(code, "Grid_File_Name") == 0) {
     i_command = 3002;
     value >> value_string;
     strcpy(Grid_File_Name, value_string.c_str());
     strcat(Grid_File_Name, ".grid");

  } else if (strcmp(code, "Number_of_Cells_Idir") == 0) {
     i_command = 3003;
     value >> NCells_Idir;
     if (NCells_Idir < 1) i_command = INVALID_INPUT_VALUE;
     if (NCells_Turbulence_Idir < NCells_Idir) NCells_Turbulence_Idir = NCells_Idir;

  } else if (strcmp(code, "Number_of_Cells_Jdir") == 0) {
     i_command = 3004;
     value >> NCells_Jdir;
     if (NCells_Jdir <1) i_command = INVALID_INPUT_VALUE;
     if (NCells_Turbulence_Jdir < NCells_Jdir) NCells_Turbulence_Jdir = NCells_Jdir;

  } else if (strcmp(code, "Number_of_Cells_Kdir") == 0) {
     i_command = 3005;
     value >> NCells_Kdir;
     if (NCells_Kdir < 1) i_command = INVALID_INPUT_VALUE;
     if (NCells_Turbulence_Kdir < NCells_Kdir) NCells_Turbulence_Kdir = NCells_Kdir;

  } else if (strcmp(code, "Number_of_Ghost_Cells") == 0) {
     i_command = 3006;
     value >> Nghost;

  } else if (strcmp(code, "Number_of_Blocks_Idir") == 0) {
     i_command = 3007;
     value >> NBlk_Idir;
     if (NBlk_Idir < 1) i_command = INVALID_INPUT_VALUE;
 
  } else if (strcmp(code, "Number_of_Blocks_Jdir") == 0) {
     i_command = 3008;
     value >> NBlk_Jdir;
     if (NBlk_Jdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Number_of_Blocks_Kdir") == 0) {
     i_command = 3009;
     value >> NBlk_Kdir;
     if (NBlk_Kdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Mesh_Stretching") == 0) {
      i_command = 3009;
      value >> value_string;
      if (strcmp("ON", value_string.c_str()) == 0) {
          Mesh_Stretching = ON;
      } else if (strcmp("OFF", value_string.c_str()) == 0) {
          Mesh_Stretching = OFF;
      } else {
          i_command = INVALID_INPUT_VALUE;
      }
            
  } else if (strcmp(code, "Stretching_Factor_Idir") == 0) {
     i_command = 3010;
     value >> Stretching_Factor_Idir;
     if (Stretching_Factor_Idir < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Stretching_Factor_Jdir") == 0) {
     i_command = 3011;
     value >> Stretching_Factor_Jdir;
     if (Stretching_Factor_Jdir < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Stretching_Factor_Kdir") == 0) {
     i_command = 3012;
     value >> Stretching_Factor_Kdir;
     if (Stretching_Factor_Kdir < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Stretching_Type_Idir") == 0) {
     i_command = 3013;
     value >> Stretching_Type_Idir;
     if (Stretching_Type_Idir < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Stretching_Type_Jdir") == 0) {
     i_command = 3014;
     value >> Stretching_Type_Jdir;
     if (Stretching_Type_Jdir < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Stretching_Type_Kdir") == 0) {
     i_command = 3015;
     value >> Stretching_Type_Kdir;
     if (Stretching_Type_Kdir < 0) i_command = INVALID_INPUT_VALUE;
      
  } else if (strcmp(code, "Mesh_Smoothing") == 0) {
      i_command = 3018;
      value >> Mesh_Smoothing;
      
  } else if (strcmp(code, "X_Shift") == 0) {
     i_command = 3016;
     value >> X_Shift;

  } else if (strcmp(code, "X_Scale") == 0) {
     i_command = 3017;
     value >> X_Scale;
       
  } else if (strcmp(code, "X_Rotate") == 0) {
     i_command = 3018;
     value >> X_Rotate;
      
  } else if (strcmp(code, "Disturb_Interior_Nodes") == 0) {
      i_command = 3018;
      value >> Disturb_Interior_Nodes;
      
  } else if (strcmp(code, "Plate_Length") == 0) {
     i_command = 3019;
     value >> Plate_Length;
     if (Plate_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Pipe_Length") == 0) {
     i_command = 3020;
     value >> Pipe_Length;
     if (Pipe_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Pipe_Radius") == 0) {
     i_command = 3021;
     value >> Pipe_Radius;
     if (Pipe_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Radius_Bluff_Body") == 0) {
     i_command = 3022;
     value >> Radius_Bluff_Body;
     if (Radius_Bluff_Body <ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Radius_Fuel_Line") == 0) {
     i_command = 3023;
     value >> Radius_Fuel_Line;
     if (Radius_Fuel_Line <ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Radius_Coflow_Inlet_Pipe") == 0) {
     i_command = 3024;
     value >> Radius_Coflow_Inlet_Pipe;
     if (Radius_Coflow_Inlet_Pipe <ZERO) i_command = INVALID_INPUT_VALUE;   

  } else if (strcmp(code, "Length_Coflow_Inlet_Pipe") == 0) {
     i_command = 3025;
     value >> Length_Coflow_Inlet_Pipe;
     if (Length_Coflow_Inlet_Pipe <ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Length_Combustor_Tube") == 0) {
     i_command = 3026;
     value >> Length_Combustor_Tube;
     if (Length_Combustor_Tube <ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "ICEMCFD_Topology_File") == 0) {
     i_command = 3027;
     value >> value_string;
     strcpy(ICEMCFD_FileNames[0], value_string.c_str());

  } else if (strcmp(code, "ICEMCFD_Family_Boco_File") == 0) {
     i_command = 3028;
     value >> value_string;
     strcpy(ICEMCFD_FileNames[1], value_string.c_str());

  } else if (strcmp(code, "ICEMCFD_Family_Topo_File") == 0) {
     i_command = 3029;
     value >> value_string;
     strcpy(ICEMCFD_FileNames[2], value_string.c_str());

  } else if (strcmp(code, "Box_Length") == 0) {
     i_command = 3030;
     value >> Box_Length;
     if (Box_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Box_Width") == 0) {
     i_command = 3031;
     value >> Box_Width;
     if (Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Box_Height") == 0) {
     i_command = 3032;
     value >> Box_Height;
     if (Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Radius_Bunsen_Burner_Fuel_Line") == 0) {
     i_command = 3033;
     value >> Radius_Bunsen_Burner_Fuel_Line;
     if (Radius_Bunsen_Burner <ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Radius_Bunsen_Burner") == 0) {
     i_command = 3034;
     value >> Radius_Bunsen_Burner;
     if (Radius_Bunsen_Burner <ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Height_Bunsen_Burner") == 0) {
     i_command = 3035;
     value >> Height_Bunsen_Burner;
     if (Height_Bunsen_Burner <ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Number_of_Cells_Turbulence_Idir") == 0) {
     i_command = 3036;
     value >> NCells_Turbulence_Idir;
     if (NCells_Turbulence_Idir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Number_of_Cells_Turbulence_Jdir") == 0) {
     i_command = 3037;
     value >> NCells_Turbulence_Jdir;
     if (NCells_Turbulence_Jdir <1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Number_of_Cells_Turbulence_Kdir") == 0) {
     i_command = 3038;
     value >> NCells_Turbulence_Kdir;
     if (NCells_Turbulence_Kdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Turbulence_Box_Length") == 0) {
     i_command = 3039;
     value >> Turbulence_Box_Length;
     if (Turbulence_Box_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Turbulence_Box_Width") == 0) {
     i_command = 3040;
     value >> Turbulence_Box_Width;
     if (Turbulence_Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Turbulence_Box_Height") == 0) {
     i_command = 3041;
     value >> Turbulence_Box_Height;
     if (Turbulence_Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;
    
  } else if (strcmp(code, "Slot_Width") == 0) {
     i_command = 3042;
     value >> Slot_Width;
     if (Slot_Width <= ZERO) i_command = INVALID_INPUT_VALUE;
  
  } else if (strcmp(code, "Number_of_Ghost_Cells") == 0) {
      i_command = 3035;
      value >> Nghost;
      if (Nghost <1) i_command = INVALID_INPUT_VALUE;

  } else {
     i_command = INVALID_INPUT_CODE;
  } /* endif */

  return i_command;
  
}

/***************************************************************************
 * Grid3D_Input_Parameters::Check_Inputs -- Check input values.            *
 ***************************************************************************/
int Grid3D_Input_Parameters::Check_Inputs(void) {

  // Input parameters are consistent.  Exit successfully.
  return 0;

}

/***************************************************************************
 * Grid3D_Input_Parameters -- Input-output operators.                      *
 ***************************************************************************/
ostream &operator << (ostream &out_file,
                      const Grid3D_Input_Parameters &IP) {

  IP.Output(out_file);
  return (out_file);

}

istream &operator >> (istream &in_file,
                      Grid3D_Input_Parameters &IP) {
 
   return in_file;

}

void Grid3D_Input_Parameters::Output(ostream &out_file) const {

   out_file << setprecision(6);

   out_file << "\n  -> Grid: "
            << Grid_Type;

   switch(i_Grid) {
     case GRID_ICEMCFD :
        out_file << 
        out_file << "\n  -> topology file : " << ICEMCFD_FileNames[0];
        out_file << "\n  -> family_boco file : " << ICEMCFD_FileNames[1];
        out_file << "\n  -> family_topo file : " << ICEMCFD_FileNames[2];
        break;
     case GRID_CUBE : 
     case GRID_PERIODIC_BOX : 
     case GRID_PERIODIC_BOX_WITH_INFLOW : 
        out_file << "\n  -> Length of Solution Domain (m): "
                 << Box_Length;
        out_file << "\n  -> Width of Solution Domain (m): "
                 << Box_Width;
        out_file << "\n  -> Height of Solution Domain (m): "
                 << Box_Height;
        break;
     case GRID_BUNSEN_INFLOW : 
        out_file << "\n  -> Radius of Bunsen Burner Fuel Domain (m): "
                 << Radius_Bunsen_Burner_Fuel_Line;
        out_file << "\n  -> Height of Bunsen Burner Solution Domain (m): "
                 << Height_Bunsen_Burner;
        out_file << "\n  -> Length of Bunsen Box Turbulence Domain (m): "
                 << Turbulence_Box_Length;
        out_file << "\n  -> Width of Bunsen Box Turbulence Domain (m): "
                 << Turbulence_Box_Width;
        out_file << "\n  -> Height of Bunsen Box Turbulence Domain (m): "
                 << Turbulence_Box_Height;
        break;
     case GRID_BUNSEN_BURNER : 
        out_file << "\n  -> Radius of Bunsen Burner Fuel Domain (m): "
                 << Radius_Bunsen_Burner_Fuel_Line;
        out_file << "\n  -> Radius of Bunsen Burner Solution Domain (m): "
                 << Radius_Bunsen_Burner;
        out_file << "\n  -> Height of Bunsen Burner Solution Domain (m): "
                 << Height_Bunsen_Burner;
        out_file << "\n  -> Length of Bunsen Box Turbulence Domain (m): "
                 << Turbulence_Box_Length;
        out_file << "\n  -> Width of Bunsen Box Turbulence Domain (m): "
                 << Turbulence_Box_Width;
        out_file << "\n  -> Height of Bunsen Box Turbulence Domain (m): "
                 << Turbulence_Box_Height;
        break;
     case GRID_BUNSEN_BOX : 
        out_file << "\n  -> Length of Solution Domain (m): "
                 << Box_Length;
        out_file << "\n  -> Width of Solution Domain (m): "
                 << Box_Width;
        out_file << "\n  -> Height of Solution Domain (m): "
                 << Box_Height;
	out_file << "\n  -> Length of Bunsen Box Turbulence Domain (m): "
                 << Turbulence_Box_Length;
        out_file << "\n  -> Width of Bunsen Box Turbulence Domain (m): "
                 << Turbulence_Box_Width;
        out_file << "\n  -> Height of Bunsen Box Turbulence Domain (m): "
                 << Turbulence_Box_Height;
        out_file << "\n  -> Slot Width of Slot Burner (m): "
                 << Slot_Width;
        break;
     case GRID_COUETTE :
        out_file << "\n  -> Length of Solution Domain (m): "
                 << Box_Length;
        out_file << "\n  -> Width of Solution Domain (m): "
                 << Box_Width;
        out_file << "\n  -> Height of Solution Domain (m): "
                 << Box_Height;
        break;
     default:
        out_file << "\n  -> Length of Solution Domain (m): "
                 << Box_Length;
        out_file << "\n  -> Width of Solution Domain (m): "
                 << Box_Width;
        out_file << "\n  -> Height of Solution Domain (m): "
                 << Box_Height;
        break;
   } /* endswitch */

   out_file << "\n  -> Number of Blocks i-direction: "
            << NBlk_Idir;
   out_file << "\n  -> Number of Blocks j-direction: " 
            << NBlk_Jdir;
   out_file << "\n  -> Number of Blocks k-direction: "
            << NBlk_Kdir;
   out_file << "\n  -> Number of Cells i-direction: "
            << NCells_Idir;
   out_file << "\n  -> Number of Cells j-direction: " 
            << NCells_Jdir;
   out_file << "\n  -> Number of Cells k-direction: " 
            << NCells_Kdir;
   out_file << "\n  -> Number of Ghost Cells on each side: "
            << Nghost;

   out_file << "\n  -> Mesh shift, scale, and rotate: " 
            << X_Shift << " " << X_Scale << " " << X_Rotate;

}
