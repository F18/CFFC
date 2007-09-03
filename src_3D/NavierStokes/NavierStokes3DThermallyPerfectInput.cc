/****************** NavierStokes3DInput.cc ************************************
  Constructors for the Input class.

  associated header file:  Input.h
***********************************************************************/

// Include required CFFC header files

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // INPUT_INCLUDED

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED
#include "NavierStokes3DThermallyPerfectState.h"
#endif // NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED

/*************************************************************
 * NavierStokes3D_Input_Parameters -- Input-output operators.       *
 *************************************************************/
template< >
ostream &operator << (ostream &out_file,
                      const Input_Parameters<NavierStokes3D_ThermallyPerfect_pState,
                                             NavierStokes3D_ThermallyPerfect_cState> &IP) {
   
   out_file << setprecision(6);
   
   out_file << "\n\n Solving 3D MulitSpecies";
   if (IP.i_Flow_Type ==  FLOWTYPE_INVISCID){
      out_file<<" Euler (Inviscid) ";
   } else {
      out_file<<" Navier-Stokes (Viscous) ";
   }
   out_file<<"equations (IBVP/BVP) "; 
   if (IP.i_Flow_Type ==  FLOWTYPE_INVISCID) {
   } else if (IP.i_Flow_Type == FLOWTYPE_LAMINAR) {
      out_file << "\n  -> Laminar flow";
   } else if (IP.i_Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      out_file << "\n  -> Turbulent flow: RANS with k-epsilon turbulence model";
      if (IP.Wall_Boundary_Treatments==0){
         out_file <<" (Automatic wall boundary treatments)";
      } else {
         if(IP.Wall_Boundary_Treatments==2){
            out_file <<" (Low-Reynolds-number formulation)";
         }else{
            out_file <<" (Wall-function formulation)";
         }
      }
   } else if (IP.i_Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      out_file << "\n  -> Turbulent flow: RANS with k-oemga turbulence model";
      if (IP.Wall_Boundary_Treatments==0){
         out_file <<" (Automatic wall boundary treatments)";
      } else {
         if(IP.Wall_Boundary_Treatments==2){
            out_file <<" (Low-Reynolds-number formulation)";
         }else{
            out_file <<" (Wall-function formulation)";
         }
      }
      
   } else if (IP.i_Flow_Type == FLOWTYPE_TURBULENT_LES) {
      out_file << "\n  -> Turbulent flow: LES ";
   } else if (IP.i_Flow_Type == FLOWTYPE_TURBULENT_DES_K_OMEGA) {
      out_file << "\n  -> Turbulent flow: DES with k-oemga SGS turbulence model ";
   } else if (IP.i_Flow_Type == FLOWTYPE_TURBULENT_DNS) {
      out_file << "\n  -> Turbulent flow: DNS ";
   }

   out_file << "\n  -> Input File Name: " 
            << IP.Input_File_Name;
   if (IP.Time_Accurate) { 
      out_file << "\n  -> Time Accurate (Unsteady) Solution";
   } else {
      out_file << "\n  -> Time Invariant (Steady-State) Solution";
   }
   
   if(IP.Gravity) {
      out_file << "\n  -> With Gravity (-z)";
   }

   if(IP.debug_level){
      out_file << "\n  -> Debug level: "
	       << IP.debug_level;
   }
   out_file << "\n  -> Time Integration: " 
            << IP.Time_Integration_Type;
   out_file << "\n  -> Number of Stages in Multi-Stage Scheme: " 
            << IP.N_Stage;
   if (IP.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
      out_file << "\n  -> Global Time Stepping";
   } else if (IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Scalar Local Time Stepping";
   } else if (IP.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Matrix Local Time Stepping";
   } else if (IP.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
      out_file << "\n  -> Low-Mach-Number Local Preconditioning (Weiss-Smith)";
   } else if (IP.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Semi-Implicit Local Time Stepping";
   } else if (IP.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
      out_file << "\n  -> Semi-Implicit Low-Mach-Number Local Preconditioned Time Stepping";
   }
   if (IP.Preconditioning == 1){
      out_file <<"\n  -> Mach Number Reference "<<IP.Mach_Number_Reference;
   }
   
   if (IP.Residual_Smoothing) {
      out_file << "\n  -> Residual Smoothing";
      out_file << "\n  -> Epsilon: " 
               << IP.Residual_Smoothing_Epsilon;
      out_file << "\n  -> Gauss_Seidel_Iterations: " 
               << IP.Residual_Smoothing_Gauss_Seidel_Iterations;
   }

   out_file << "\n  -> Reconstruction: " 
            << IP.Reconstruction_Type;
   out_file << "\n  -> Limiter: " 
            << IP.Limiter_Type;   
   if (IP.Limiter_Type != LIMITER_ZERO && IP.Freeze_Limiter) {
      out_file << "\n  -> Freeze Limiter when L2-norm of residual is < "
	       << IP.Freeze_Limiter_Residual_Level;
   } 
   out_file << "\n  -> Flux Function: " 
            << IP.Flux_Function_Type;
   
   out_file << "\n  -> Reaction Mechanism: " 
            << IP.react_name;
   out_file << "\n  -> Species: "<<IP.Wo.ns
            << "\n  -> Initial mass fractions: ";
   for(int i=0; i<IP.Wo.ns; i++){
      out_file  <<"c["<<IP.multispecies[i]<<"]= ";
      out_file  << IP.Wo.spec[i].c<<", ";
   } 
//    if(IP.i_Flow_Type != FLOWTYPE_INVISCID){
//       out_file << "\n  -> Schmidt Numbers for Viscous flow: ";
//       for(int i=0; i <IP.Wo.ns; i++){
//          out_file  <<"Sc["<<IP.multispecies[i]<<"]= ";
//          out_file << IP.Schmidt[i]<<", ";
//       }
//    }
   out_file << "\n  -> Initial Conditions: " 
            << IP.ICs_Type;
   switch(IP.i_ICs) {
   case IC_CONSTANT :
      out_file << "\n  -> Pressure (kPa): " 
               << IP.Pressure/THOUSAND;
      out_file << "\n  -> Temperature (K): " 
               << IP.Temperature;
      out_file << "\n  -> Mach Number: " 
               << IP.Mach_Number;
      out_file << "\n  -> Flow Angle: " 
               << IP.Flow_Angle;
      break;
   case IC_UNIFORM :
      out_file << "\n  -> Pressure (kPa): " 
               << IP.Pressure/THOUSAND;
      out_file << "\n  -> Temperature (K): " 
               << IP.Temperature;
      out_file << "\n  -> Mach Number: " 
               << IP.Mach_Number;
      out_file << "\n  -> Flow Angle (degrees): " 
               << IP.Flow_Angle;
      break;
   case IC_SOD_XDIR :
      break;
   case IC_SOD_YDIR :
      break;
   case IC_GROTH_XDIR :
      break;
   case IC_GROTH_YDIR :
      break;
   case IC_EINFELDT_XDIR :
      break;
   case IC_EINFELDT_YDIR :
      break;
   case IC_SHOCK_BOX :
      break;
   case IC_HIGH_PRESSURE_RESERVOIR :
      break;
   case IC_LOW_PRESSURE_RESERVOIR :
      break;
   case IC_RIEMANN :
      break;
   case IC_RIEMANN_XDIR :
      break;
   case IC_RIEMANN_YDIR :
      out_file << "\n  -> Pressure (kPa): " 
               << IP.Pressure/THOUSAND;
      out_file << "\n  -> Temperature (K): " 
               << IP.Temperature;
      out_file << "\n  -> Mach Number: " 
               << IP.Mach_Number;
      out_file << "\n  -> Flow Angle (degrees): " 
               << IP.Flow_Angle;
      break;
   case IC_GAS_MIX :
      break;
   case IC_CHEM_CORE_FLAME:
      break;
   case IC_CHEM_INVERSE_FLAME:
      break;
   case IC_CHEM_1DFLAME:
      break;
   case IC_PRESSURE_GRADIENT_X:
      break;
   case IC_PRESSURE_GRADIENT_Y:
      break;
   case IC_VISCOUS_COUETTE:
      break;
   case IC_VISCOUS_COUETTE_PRESSURE_GRADIENT: 
      break;
      
   default:
      break;
   } /* endswitch */
   out_file << "\n  -> Grid: "
            << IP.IP_Grid.Grid_Type;
   switch(IP.IP_Grid.i_Grid) {

   case GRID_CUBE :
      out_file << "\n  -> Length of Solution Domain (m): "
               << IP.IP_Grid.Box_Length;
      out_file << "\n  -> Width of Solution Domain (m): "
               << IP.IP_Grid.Box_Width;
      out_file << "\n  -> Height of Solution Domain (m): "
               << IP.IP_Grid.Box_Height;
      break;
   case GRID_COUETTE :
      out_file << "\n  -> Length of Solution Domain (m): "
               << IP.IP_Grid.Box_Length;
      out_file << "\n  -> Width of Solution Domain (m): "
               << IP.IP_Grid.Box_Width;
      out_file << "\n  -> Height of Solution Domain (m): "
               << IP.IP_Grid.Box_Height;
      out_file << "\n  -> Moving Wall Velocity in x direction (m/s): "
               << IP.Moving_wall_velocity.x;
      out_file << "\n  -> Moving Wall Velocity in y direction (m/s): "
               << IP.Moving_wall_velocity.y;
      out_file << "\n  -> Moving Wall Velocity in z direction (m/s): "
               << IP.Moving_wall_velocity.z;
      break;

  /*  case GRID_RECTANGULAR_BOX : */
/*       out_file << "\n  -> Width of Solution Domain (m): " */
/*                << IP.Box_Width; */
/*       out_file << "\n  -> Height of Solution Domain (m): " */
/*                << IP.Box_Height; */
/*       break; */
/*    case GRID_COUETTE : */
/*       out_file << "\n  -> Width of Solution Domain (m): " */
/*                << IP.Box_Width; */
/*       out_file << "\n  -> Height of Solution Domain (m): " */
/*                << IP.Box_Height; */
/*       out_file << "\n  -> Moving Wall Velocity (m/s): " */
/*                << IP.Moving_wall_velocity; */
/*       break; */
/*    case GRID_1DFLAME : */
/*       out_file << "\n  -> Width of Solution Domain (m): " */
/*                << IP.Box_Width; */
/*       out_file << "\n  -> Height of Solution Domain (m): " */
/*                << IP.Box_Height; */
/*       break; */
/*    case GRID_LAMINAR_FLAME : */
/*       out_file << "\n  -> Width of Solution Domain (m): " */
/*                << IP.Pipe_Length; */
/*       out_file << "\n  -> Height of Solution Domain (m): " */
/*                << IP.Pipe_Radius; */
/*       break; */
/*    case GRID_FLAT_PLATE : */
/*       out_file << "\n  -> Plate Length (m): " */
/*                << IP.Plate_Length; */
/*       break; */
/*    case GRID_PIPE : */
/*       out_file << "\n  -> Pipe Length (m): " */
/*                << IP.Pipe_Length; */
/*       out_file << "\n  -> Pipe Radius (m): " */
/*                << IP.Pipe_Radius; */
/*       break; */
/*    case GRID_CAVITY_FLOW : */
/*       out_file << "\n  -> Cavity Depth (m): " */
/*                << IP.Box_Height; */
/*       break; */
      
/*    case GRID_BLUNT_BODY : */
/*       out_file << "\n  -> Cylinder Radius (m): " */
/*                << IP.Blunt_Body_Radius; */
/*       break; */
/*    case GRID_BLUFF_BODY : */
/*       out_file << "\n  -> Shroud Length (m): " */
/*                << IP.Length_Shroud; */
/*       out_file << "\n  -> Shroud Radius (m): " */
/*                << IP.Radius_Shroud; */
/*       out_file << "\n  -> Bluff Body Length (m): " */
/*                << IP.Length_BluffBody; */
/*       out_file << "\n  -> Bluff Body Radius (m): " */
/*                << IP.Radius_BluffBody; */
/*       out_file << "\n  -> Fuel Orifice Radius (m): " */
/*                << IP.Radius_Orifice; */
/*       break; */
/*    case GRID_DUMP_COMBUSTOR : */
/*       out_file << "\n  -> Inlet Pipe Length (m): " */
/*                << IP.Length_Inlet_Pipe; */
/*       out_file << "\n  -> Inlet Pipe Radius (m): " */
/*                << IP.Radius_Inlet_Pipe; */
/*       out_file << "\n  -> Combustor Tube Length (m): " */
/*                << IP.Length_Combustor_Tube; */
/*       out_file << "\n  -> Combustor Tube Radius(m): " */
/*                << IP.Radius_Combustor_Tube; */
/*       break; */
/*    case GRID_ROCKET_MOTOR : */
/*       out_file << "\n  -> Length of Grain (m): " */
/*                << IP.Grain_Length; */
/*       out_file << "\n  -> Radius of Grain (m): " */
/*                << IP.Grain_Radius; */
/*       out_file << "\n  -> Distance from Grain to Nozzle Throat (m): " */
/*                << IP.Grain_To_Throat_Length; */
/*       out_file << "\n  -> Length of the Nozzle (m): " */
/*                << IP.Nozzle_Length; */
/*       out_file << "\n  -> Radius of the Nozzle at Throat (m): " */
/*                << IP.Nozzle_Radius_Throat; */
/*       out_file << "\n  -> Radius of the Nozzle at Exit (m): " */
/*                << IP.Nozzle_Radius_Exit; */
/*       break; */
/*    case GRID_CIRCULAR_CYLINDER : */
/*       out_file << "\n  -> Cylinder Radius (m): " */
/*                << IP.Cylinder_Radius; */
/*       break; */
/*    case GRID_ELLIPSE : */
/*       out_file << "\n  -> Width of Ellipse along x-axis (m): " */
/*                << IP.Ellipse_Length_X_Axis; */
/*       out_file << "\n  -> Height of Ellipse along y-axis (m): " */
/*                << IP.Ellipse_Length_Y_Axis; */
/*       break; */
/*    case GRID_NACA_AEROFOIL : */
/*       out_file << "\n  -> NACA " */
/*                << IP.NACA_Aerofoil_Type; */
/*       out_file << "\n  -> Chord Length (m): " */
/*                << IP.Chord_Length; */
/*       break; */
/*    case GRID_FREE_JET : */
/*       out_file << "\n  -> Orifice Radius (m): " */
/*                << IP.Orifice_Radius; */
/*       break; */
/*    case GRID_WEDGE : */
/*       out_file << "\n  -> Wedge Angle (degrees): " << IP.Wedge_Angle; */
/*       out_file << "\n  -> Wedge Length (m): " << IP.Wedge_Length; */
/*       break; */
/*    case GRID_ICEMCFD : */
/*       break; */
/*    case GRID_READ_FROM_DEFINITION_FILE : */
/*       break; */
/*    case GRID_READ_FROM_GRID_DATA_FILE : */
/*       break; */
   default:
      out_file << "\n  -> Length of Solution Domain (m): "
               << IP.IP_Grid.Box_Length;
      out_file << "\n  -> Width of Solution Domain (m): "
               << IP.IP_Grid.Box_Width;
      out_file << "\n  -> Height of Solution Domain (m): "
               << IP.IP_Grid.Box_Height;
      break;
   } /* endswitch */
 
   if (IP.Number_of_Initial_Mesh_Refinements > 0)
      out_file << "\n  -> Number of Initial Mesh Refinements : " 
               << IP.Number_of_Initial_Mesh_Refinements;
   if (IP.Number_of_Uniform_Mesh_Refinements > 0)
      out_file << "\n  -> Number of Uniform Mesh Refinements : " 
               << IP.Number_of_Uniform_Mesh_Refinements;
   if (IP.Number_of_Boundary_Mesh_Refinements > 0)
      out_file << "\n  -> Number of Boundary Mesh Refinements : " 
               << IP.Number_of_Boundary_Mesh_Refinements;
   out_file << "\n  -> Number of Blocks i-direction: "
            << IP.IP_Grid.NBlk_Idir;
   out_file << "\n  -> Number of Blocks j-direction: " 
            << IP.IP_Grid.NBlk_Jdir;
   out_file << "\n  -> Number of Blocks k-direction: "
            << IP.IP_Grid.NBlk_Kdir;
   out_file << "\n  -> Number of Cells i-direction: "
            << IP.IP_Grid.NCells_Idir;
   out_file << "\n  -> Number of Cells j-direction: " 
            << IP.IP_Grid.NCells_Jdir;
   out_file << "\n  -> Number of Cells k-direction: " 
            << IP.IP_Grid.NCells_Kdir;
   out_file << "\n  -> CFL Number: " 
            << IP.CFL_Number;
   out_file << "\n  -> Maximum Time (ms): " 
            << IP.Time_Max*THOUSAND;
   out_file << "\n  -> Maximum Number of Time Steps (Iterations): " 
            << IP.Maximum_Number_of_Time_Steps;
   out_file << "\n  -> Number of Processors: " 
            << IP.Number_of_Processors;
   out_file << "\n  -> Number of Blocks Per Processor: " 
            << IP.Number_of_Blocks_Per_Processor;
   out_file << "\n  -> Output File Name: " 
            << IP.Output_File_Name;
   out_file << "\n  -> Output Format: " 
            << IP.Output_Format_Type;
   out_file << "\n  -> Restart Solution Save Frequency: "
            << IP.Restart_Solution_Save_Frequency
            << " steps (iterations)"; 
   if (IP.Time_Accurate_Plot_Freq !=0 && IP.Time_Accurate){
      out_file << "\n  -> Time Accurate Solution Plot Frequency: "
	       << IP.Time_Accurate_Plot_Freq
	       << " steps (iterations)"; 
   }
   return (out_file);
}

istream  &operator >> (istream &in_file,
                       Input_Parameters<NavierStokes3D_ThermallyPerfect_pState, 
                                        NavierStokes3D_ThermallyPerfect_cState> &IP) {
   return (in_file);
}



