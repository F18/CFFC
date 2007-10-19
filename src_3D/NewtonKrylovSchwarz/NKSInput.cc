/* NKSInput.cc: Definition of NKS_Input_Parameters class member functions. */

/* Include the NKSInput header file. */

#ifndef _NKSINPUT_INCLUDED
#include "NKSInput.h"
#endif // _NKSINPUT_INCLUDED

/* Define member functions. */

/***************************************************************************
 * NKS_Input_Parameters::Broadcast -- Broadcast to all processors.         *
 ***************************************************************************/
void NKS_Input_Parameters::Broadcast(void) {

#ifdef _MPI_VERSION
    //Newton 
    MPI::COMM_WORLD.Bcast(&(Maximum_Number_of_NKS_Iterations), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Overall_Tolerance), 
			  1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Relaxation_multiplier), 
                          1, 
                          MPI::DOUBLE, 0);
   
    // Dual Time Stepping
    MPI::COMM_WORLD.Bcast(&(Dual_Time_Stepping), 
			  1, 
			  MPI::INT, 0);  //BOOL
    MPI::COMM_WORLD.Bcast(&(Physical_Time_Integration), 
			  1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Physical_Time_CFL_Number),
			  1, 
			  MPI::DOUBLE, 0);  
    MPI::COMM_WORLD.Bcast(&(Physical_Time_Step),
			  1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Maximum_Number_of_DTS_Steps),
			  1, 
			  MPI::INT, 0);

    MPI::COMM_WORLD.Bcast(&(Freeze_Limiter_Immediately), 
                          1, 
                          MPI::INT, 0);
 
    // Finite Time Step
    MPI::COMM_WORLD.Bcast(&(Finite_Time_Step), 
			  1, 
			  MPI::INT, 0);  //BOOL
    MPI::COMM_WORLD.Bcast(&(Finite_Time_Step_Initial_CFL), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Finite_Time_Step_Final_CFL), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Finite_Time_Step_Max_CFL), 
			  1, 
			  MPI::DOUBLE, 0);

    // GMRES
    MPI::COMM_WORLD.Bcast(&(Maximum_Number_of_GMRES_Iterations), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(GMRES_Restart), 
                          1, 
                          MPI::INT, 0); 
    MPI::COMM_WORLD.Bcast(&(GMRES_Overlap), 
			  1, 
			  MPI::INT, 0); 
    MPI::COMM_WORLD.Bcast(&(GMRES_Initial_Tolerance), 
                          1, 
                          MPI::DOUBLE, 0); 
    MPI::COMM_WORLD.Bcast(&(GMRES_Final_Tolerance), 
                          1, 
                          MPI::DOUBLE, 0); 
    MPI::COMM_WORLD.Bcast(&(Normalization), 
			  1, 
			  MPI::INT, 0); //BOOL
    MPI::COMM_WORLD.Bcast(&(GMRES_CHECK),
			  1, 
			  MPI::INT, 0); //BOOL
    MPI::COMM_WORLD.Bcast(&(GMRES_Frechet_Derivative_Order),
			  1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Epsilon_Naught), 1, MPI::DOUBLE, 0);

    // GMRES Preconditioner Type
    MPI::COMM_WORLD.Bcast(&(GMRES_Block_Preconditioner), 
			  1, 
			  MPI::INT, 0);
    //Preconditoner Jacobian Order
    MPI::COMM_WORLD.Bcast(&(Jacobian_Order), 
                          1, 
                          MPI::INT, 0);    
    // GMRES ILUK_Level_of_Fill
    MPI::COMM_WORLD.Bcast(&(GMRES_ILUK_Level_of_Fill), 
			  1, 
			  MPI::INT, 0);

    // Output control
    MPI::COMM_WORLD.Bcast(&(Detect_Convergence_Stall), 
                          1, 
                          MPI::INT, 0); // bool
    MPI::COMM_WORLD.Bcast(&(DCS_Window), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Output_Format), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Output_Precision), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Output_Width), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(NKS_Write_Output_Cells_Freq), 
                          1, 
                          MPI::INT, 0); 

    MPI::COMM_WORLD.Bcast(&(Min_Number_of_Newton_Steps_With_Zero_Limiter), 
			  1, 
			  MPI::INT, 0);    
    MPI::COMM_WORLD.Bcast(&(Min_Number_of_Newton_Steps_Requiring_Jacobian_Update), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&( Min_L2_Norm_Requiring_Jacobian_Update), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&( Min_Finite_Time_Step_Norm_Ratio), 
			  1, 
			  MPI::DOUBLE, 0);
#endif

}

/***************************************************************************
 * NKS_Input_Parameters::Parse_Next_Input_Control_Parameter - Parse input. *
 ***************************************************************************/
int NKS_Input_Parameters::Parse_Next_Input_Control_Parameter(char *code, 
                                                             stringstream &value) {

// Returns:
//  - INVALID_INPUT_VALUE if code is valid but value is invalid
//  - INVALID_INPUT_CODE  if unknown code

  int i_command = INVALID_INPUT_CODE;
  string value_string;

  // NEWTON 
  if (strcmp(code, "Maximum_Number_of_NKS_Iterations") == 0) {
    i_command = 1001;
    value >> Maximum_Number_of_NKS_Iterations;
  
  } else if (strcmp(code, "NKS_Overall_Tolerance") == 0) {
    i_command = 1002;
    value >> Overall_Tolerance;
  
  } else if (strcmp(code, "Relaxation_multiplier") == 0) {
    i_command = 1003;
    value >> Relaxation_multiplier;

    //DUAL TIME STEPPING 
  } else if (strcmp(code, "NKS_Dual_Time_Stepping") == 0) {
    i_command = 1004;
    value >> value_string;
    if (value_string == "ON" || 
        value_string == "TRUE" || 
        value_string == "1") {
        Dual_Time_Stepping = true;
    } else if (value_string == "OFF" || 
               value_string == "FALSE" || 
               value_string == "0") {
        Dual_Time_Stepping = false;
    } else {
       i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code, "NKS_Physical_Time_Integration") == 0) {
    i_command = 65;  
    value >> value_string;
    if ( value_string == "Implicit_Euler" ) {
      Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
    } else if( value_string == "Second_Order_Backwards") {
      Physical_Time_Integration = TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD;
    } else {
      Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
    }
    
  } else if (strcmp(code, "NKS_Physical_Time_CFL") == 0) {
    i_command = 65;
    value >> Physical_Time_CFL_Number;
   
  } else if (strcmp(code, "NKS_Physical_Time_Step") == 0) {
    i_command = 65;
    value >> Physical_Time_Step;

  } else if (strcmp(code, "Maximum_Number_of_NKS_DTS_Steps" ) == 0) {
    i_command = 64;
    value >> Maximum_Number_of_DTS_Steps;

  // FINITE TIME
  } else if (strcmp(code, "NKS_Finite_Time_Step") == 0) {
    i_command = 1007; 
    value >> value_string;
    if (value_string == "ON" || 
        value_string == "TRUE" || 
        value_string == "1") {
       Finite_Time_Step = true;
    } else if (value_string == "OFF" || 
               value_string == "FALSE" || 
               value_string == "0") {
       Finite_Time_Step = true;
    } else {
       i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code, "NKS_Finite_Time_Step_Initial_CFL") == 0) {
    i_command = 1008;
    value >> Finite_Time_Step_Initial_CFL;

  } else if (strcmp(code, "NKS_Finite_Time_Step_Final_CFL") == 0) {
    i_command = 1009;
    value >> Finite_Time_Step_Final_CFL;

  } else if (strcmp(code, "NKS_Finite_Time_Step_Max_CFL") == 0) {
    i_command = 1010;
    value >> Finite_Time_Step_Max_CFL;

  // GMRES 
  } else if (strcmp(code, "Maximum_Number_of_GMRES_Iterations") == 0) {
    i_command = 1011;
    value >> Maximum_Number_of_GMRES_Iterations;

  } else if (strcmp(code, "GMRES_Restart") == 0) {
    i_command = 1012;
    value >> GMRES_Restart;

  } else if (strcmp(code, "GMRES_Overlap") == 0) {
    i_command = 1013;
    value >> GMRES_Overlap;

  } else if (strcmp(code, "GMRES_Normalization") == 0) {
    i_command = 1014;
    value >> value_string;
    if (value_string == "ON" || 
        value_string == "TRUE" || 
        value_string == "1") {
       Normalization = true;
    } else if (value_string == "OFF" || 
               value_string == "FALSE" || 
               value_string == "0") {
       Normalization = false;
    } else {
       i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code, "GMRES_Tolerance") == 0) { // backward compatible
    i_command = 1015;
    value >> GMRES_Initial_Tolerance;
    GMRES_Final_Tolerance = GMRES_Initial_Tolerance; 

  } else if (strcmp(code, "GMRES_Initial_Tolerance") == 0) {
    i_command = 1016;
    value >> GMRES_Initial_Tolerance;

  } else if (strcmp(code, "GMRES_Final_Tolerance") == 0) {
    i_command = 1017;
    value >> GMRES_Final_Tolerance;

  } else if (strcmp(code, "GMRES_Check") == 0) {
    i_command = 1018;
    value >> value_string;
    if (value_string == "ON" || 
        value_string == "TRUE" || 
        value_string == "1") {
       GMRES_CHECK = true;
    } else if (value_string == "OFF" || 
               value_string == "FALSE" || 
               value_string == "0") {
       GMRES_CHECK = false;
    } else {
       i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code, "GMRES_Frechet_Derivative_Order") == 0) {
    i_command = 1019;
    value >> value_string;
    if (value_string == "First_Order") {
       GMRES_Frechet_Derivative_Order = FIRST_ORDER;
    } else if(value_string == "Second_Order") {
       GMRES_Frechet_Derivative_Order = SECOND_ORDER;
    } else {
       i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code, "GMRES_Epsilon_Naught") == 0) {
    i_command = 1020;
    value >> Epsilon_Naught;

  } else if (strcmp(code, "GMRES_Block_Preconditioner") == 0) {
    i_command = 1021;
    value >> value_string;
    if (value_string == "Diagonal") {
       GMRES_Block_Preconditioner = Block_Jacobi;
    } else if(value_string == "ILUK") {
       GMRES_Block_Preconditioner = Block_ILUK;
    } else {
       i_command = INVALID_INPUT_VALUE;  
    }

  } else if (strcmp(code, "GMRES_ILUK_Level_of_Fill") == 0) {
    i_command = 1022;
    value >> GMRES_ILUK_Level_of_Fill;

  } else if (strcmp(code, "Jacobian_Order") == 0) {
    i_command = 1023;
    value >> value_string;
    if (value_string == "Source_Terms_Only") {
       Jacobian_Order = SOURCE_TERMS_ONLY;
    } else if (value_string == "First_Order_Inviscid_HLLE") {
       Jacobian_Order = FIRST_ORDER_INVISCID_HLLE;
    } else if (value_string == "First_Order_Inviscid_Roe") {
       Jacobian_Order = FIRST_ORDER_INVISCID_ROE; 
    } else if (value_string == "First_Order_Inviscid_AUSM_plus_up") {
       Jacobian_Order = FIRST_ORDER_INVISCID_AUSMPLUSUP;   
    } else if (value_string == "Second_Order_Diamond_Path_with_HLLE") {
       Jacobian_Order = SECOND_ORDER_DIAMOND_WITH_HLLE;
    } else if (value_string == "Second_Order_Diamond_Path_with_Roe") {
       Jacobian_Order = SECOND_ORDER_DIAMOND_WITH_ROE;	
    } else if (value_string == "Second_Order_Diamond_Path_with_AUSM_plus_up") {
       Jacobian_Order = SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP;
    } else {
       i_command = INVALID_INPUT_VALUE; 
    }

  } else if (strcmp(code, "NKS_Detect_Convergence_Stall") == 0) {
    i_command = 1024;
    value >> value_string;
    if (value_string == "ON" || 
        value_string == "TRUE" || 
        value_string == "1") {
       Detect_Convergence_Stall = true;
    } else if (value_string == "OFF" || 
               value_string == "FALSE" || 
               value_string == "0") {
       Detect_Convergence_Stall = false;
    } else {
       i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code, "NKS_DCS_Window") == 0) {
    i_command = 1025;
    value >> DCS_Window;

  } else if (strcmp(code, "NKS_Output_Format") == 0) {
    i_command = 1026;
    value >> value_string;
    Output_Format = 0;

  } else if (strcmp(code, "NKS_Output_Precision") == 0) {
    i_command = 1027;
    value >> Output_Precision;

  } else if (strcmp(code, "NKS_Output_Width") == 0) {
    i_command = 1028;
    value >> Output_Width;

  } else if (strcmp(code, "NKS_Freeze_Limiter_Immediately") == 0) {
    value >> i_command;
    switch (i_command) {
       case 1: 
         Freeze_Limiter_Immediately = FLI_NO;  
         break;
       case 2: 
         Freeze_Limiter_Immediately = FLI_YES; 
         break;
       default: 
         Freeze_Limiter_Immediately = FLI_NOT_USED; 
         break;
    }
    i_command = 1029;

  } else if (strcmp(code, "NKS_Write_Output_Cells_Freq") == 0) {
    i_command = 1030;
    value >> NKS_Write_Output_Cells_Freq;

  } else if (strcmp(code, "NKS_Min_Number_of_Newton_Steps_With_Zero_Limiter") == 0) {
    i_command = 1034;
    value >> Min_Number_of_Newton_Steps_With_Zero_Limiter;  

  } else if (strcmp(code, "NKS_Min_Number_of_Newton_Steps_Requiring_Jacobian_Update") == 0) {
    i_command = 1035;
    value >> Min_Number_of_Newton_Steps_Requiring_Jacobian_Update;

  } else if (strcmp(code, "NKS_Min_L2_Norm_Requiring_Jacobian_Update") == 0) {
    i_command = 1036;
    value >> Min_L2_Norm_Requiring_Jacobian_Update;
    
  } else if (strcmp(code, "NKS_Min_Finite_Time_Step_Norm_Ratio") == 0) {
    i_command = 1037;
    value >> Min_Finite_Time_Step_Norm_Ratio;
 
  } else {
    i_command = INVALID_INPUT_CODE;

  } /* endif */

  return i_command;
  
}

/***************************************************************************
 * NKS_Input_Parameters::Check_Inputs -- Check input values.               *
 ***************************************************************************/
int NKS_Input_Parameters::Check_Inputs(void) {

  // Input parameters are consistent.  Exit successfully.
  return 0;

}

/***************************************************************************
 * NKS_Input_Parameters::Memory Estimates - Estimate memory usage.         *
 ***************************************************************************/
void NKS_Input_Parameters::Memory_Estimates(const int &blocksize, 
                                            const int &block_mat_size,
                                            const int &used_blocks){

  cerr <<" \n NKS_Input_Parameters::Memory_Estimates not working yet ";

//   int INT = sizeof(int);  //4 bytes
//   int DOUBLE = sizeof(double); //8bytes
//   double MB = 1024.0*1024.0;

//   double GMRES = DOUBLE*( (blocksize+1) + 
// 			  (GMRES_Restart +1) +
// 			  (GMRES_Restart*2) +
// 			  (block_mat_size*blocksize) +
// 			  (GMRES_Restart +1)*(GMRES_Restart) +
// 			  (GMRES_Restart)*(block_mat_size*blocksize) +
// 			  (GMRES_Restart+1)*(block_mat_size*blocksize) );
  
//   int nnz = 5*block_mat_size;  //ONLY FOR 1st order, 13 for 2nd
  
//   int JACOBIAN = DOUBLE*(nnz*blocksize*blocksize) +
//                  INT * ( (block_mat_size +1)*3 + nnz );

				    				    
//   //ONLY FOR ILU
//   int upper = nnz*blocksize*blocksize / 2 - block_mat_size;
//   int lower = upper;

//   int PRECON = DOUBLE*( (block_mat_size*blocksize*blocksize) +
// 			( GMRES_ILUK_Level_of_Fill + 2)*(upper) +
// 			( GMRES_ILUK_Level_of_Fill + 2)*(upper) + block_mat_size ) +
//                INT * (  (block_mat_size +1)*2 +
// 			( GMRES_ILUK_Level_of_Fill + 2)*(upper) +
// 			( GMRES_ILUK_Level_of_Fill + 2)*(upper) + block_mat_size );

  
//   //Output
//   cout<<" NKS Memory Requirement Estimate (MB)";
//   cout<<"\n Preconditioner (ILU) = "<<PRECON/MB; 
//   cout<<"\n GMRES                = "<<GMRES/MB; 
//   cout<<"\n Jacobian             = "<<JACOBIAN/MB; 
//   cout<<"\n Total                = "<<(PRECON+JACOBIAN+GMRES)/MB; 
//   cout<<endl;
//   for (int star=0;star<75;star++) {cout<<"*";}
//   cout <<endl;

} 

/***************************************************************************
 * NKS_Input_Parameters -- Input-output operators.                         *
 ***************************************************************************/
ostream &operator << (ostream &out_file,
                      const NKS_Input_Parameters &IP) {
  
  IP.Output(out_file);
  return (out_file);

}

istream &operator >> (istream &in_file,
                      NKS_Input_Parameters &IP) {

  return in_file;

}

void NKS_Input_Parameters::Output(ostream &fout) const {

  //fout.setf(ios::scientific);
  fout.unsetf(ios::scientific);
  fout << " " << endl;
  for (int star=0;star<75;star++){fout <<"*";}
  fout << "\n********                   Newton-Krylov-Schwarz                 **********" << endl;   
  for (int star=0;star<75;star++){fout<<"*";}
  
  fout <<"\n Overall Tolerance     ====> " << Overall_Tolerance << endl;

  fout << " Relaxation Multiplier ====> " << Relaxation_multiplier << endl;

  //DTS
  fout <<" Time Accurate (DTS)   ====> ";
  if (Dual_Time_Stepping) { 
     fout << "ON\n"; 
     if ( Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {
       fout<<" DTS Time Integration  ====> Implicit Euler \n";
     } else if ( Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {
       fout<<" DTS Time Integration  ====> Second Order Backwards \n";
     } 
     if( Physical_Time_Step > ZERO){
       fout <<" DTS Fixed Time Step   ====> " << Physical_Time_Step << endl; 
     } else { 
       fout <<" DTS CFL Number        ====> " << Physical_Time_CFL_Number << endl;   
     }
     fout <<" DTS Max Steps         ====> " << Maximum_Number_of_DTS_Steps  << endl;        
  } else { 
    fout << "OFF\n"; 
  }
  
  //Finite Time Step
  if (Finite_Time_Step == ON) {     
    fout <<" Finite Time Step      ====> ON" << endl;
    fout <<" Initial_CFL           ====> " << Finite_Time_Step_Initial_CFL << endl;
		// But not everyone uses all of these ...
    fout <<" Final_CFL             ====> " << Finite_Time_Step_Final_CFL << endl;
    fout <<" Max_CFL               ====> " << Finite_Time_Step_Max_CFL << endl;
  } else {
    fout <<" Finite Time Step      ====> OFF" << endl; 
  } /* endif */ 
  
  // GMRES 
  fout <<" Maximum GMRES Its.    ====> " << Maximum_Number_of_GMRES_Iterations<< endl;
  fout <<" GMRES Restart Its.    ====> " << GMRES_Restart << endl;
  fout <<" Level of Overlap      ====> " << GMRES_Overlap << endl;
	if (fabs(GMRES_Initial_Tolerance-GMRES_Final_Tolerance) < 1e-10) {
  fout <<" GMRES Tolerance       ====> " << GMRES_Initial_Tolerance << endl;
	} else {
	 fout <<" GMRES Initial Tol     ====> " << GMRES_Initial_Tolerance << endl;
	 fout <<" GMRES Final Tol       ====> " << GMRES_Final_Tolerance << endl;
	}
  
  if (Normalization == ON) {
    fout <<" Normalization         ====> ON" << endl;
  } else {
    fout <<" Normalization         ====> OFF" << endl; 
  } 
  
  if (GMRES_CHECK) {
    fout <<" GMRES Check           ====> ON" << endl;
  }
  if( GMRES_Frechet_Derivative_Order == FIRST_ORDER) {
    fout<<   " Frechet Derivative    ====> First Order "<< endl;
  } else if ( GMRES_Frechet_Derivative_Order == SECOND_ORDER) {
    fout<<   " Frechet Derivative    ====> Second Order "<< endl;
  }
  fout <<" Matrix Free Epsilon0  ====> " << Epsilon_Naught << endl;
  
  // Precondtioner
  fout<< " Approximate Jacobian  ====>";
  if(Jacobian_Order == SOURCE_TERMS_ONLY) {
    fout<<" Source Terms Only "<<endl;
  } else if(Jacobian_Order == FIRST_ORDER_INVISCID_HLLE){
   fout<<" First Order Inviscid HLLE"<<endl;
  } else if (Jacobian_Order == FIRST_ORDER_INVISCID_ROE){
    fout<<" First Order Inviscid Roe"<<endl;
  } else if (Jacobian_Order == FIRST_ORDER_INVISCID_AUSMPLUSUP){
    fout<<" First Order Inviscid AUSM plus up"<<endl;   
  } else if (Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_HLLE){
    fout<<" Second Order Diamond Path with HLLE"<<endl;
  } else if (Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_ROE){
    fout<<" Second Order Diamond Path with Roe"<<endl;   
  } else if (Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP){
    fout<<" Second Order Diamond Path with AUSM plus up"<<endl;   
  }
  
  if (GMRES_Block_Preconditioner == Block_ILUK) {        
    fout << " Local Preconditioner  ====> ILU("<< GMRES_ILUK_Level_of_Fill <<")" << endl;
  } else if (GMRES_Block_Preconditioner == Block_Jacobi){   // Diagonal
    fout << " Local Preconditioner  ====> Diagonal" << endl; 
  } 

  if (Detect_Convergence_Stall) {
    fout << " Detect Conv. Stall    ====> ON";
    fout << " (Window: " << DCS_Window << ")" << endl;
  }
  switch (Freeze_Limiter_Immediately) {
    case FLI_NO:  fout <<" Freeze Lim Immediately ===> OFF" << endl; break;
    case FLI_YES: fout <<" Freeze Lim Immediately ===> ON"  << endl; break;
    default: break;
  }

  if (NKS_Write_Output_Cells_Freq > 0) {
    fout << " Write Output Freq     ====> " << NKS_Write_Output_Cells_Freq << endl; 
  }
  
  //End 
  for (int star=0;star<75;star++){fout<<"*";}
  fout << endl;

}
