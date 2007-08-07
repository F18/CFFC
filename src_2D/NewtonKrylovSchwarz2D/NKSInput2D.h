/* NKS2DInput.h:  Header file declaring related Inputs and functions */

#ifndef _NKS2DINPUT_INCLUDED
#define _NKS2DINPUT_INCLUDED
   
enum Block_Preconditioners { Block_ILUK, 
			     Block_Jacobi, 
			     Something_else };

enum Jacobian_Orders { SOURCE_TERMS_ONLY,
                       FIRST_ORDER_INVISCID_HLLE,
		       FIRST_ORDER_INVISCID_ROE,
		       FIRST_ORDER_INVISCID_AUSMPLUSUP,
		       SECOND_ORDER_DIAMOND_WITH_HLLE,
		       SECOND_ORDER_DIAMOND_WITH_ROE,
		       SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP,
                       SECOND_ORDER_OTHER };

enum Frechet_Derivatives { FIRST_ORDER = 1,
			   SECOND_ORDER = 2 };


/************************************************
 * Class: NKS_Input_Parameters                  *
 ***********************************************/ 
class NKS_Input_Parameters{
 private:
 public:
  // Newton Parameters
  int    Maximum_Number_of_NKS_Iterations;  //Outer
  double Overall_Tolerance;

  // Implicit Euler Parameters
  bool   Finite_Time_Step;   
  double Finite_Time_Step_Initial_CFL;
  double Finite_Time_Step_Max_CFL;

  // GMRES parameters 
  int    Maximum_Number_of_GMRES_Iterations; //Inner
  int    GMRES_Restart;
  int    GMRES_Overlap;
  double GMRES_Tolerance; 
  bool   Normalization;   

  // Matrix Free
  bool   GMRES_CHECK;
  int    GMRES_Frechet_Derivative_Order;

  // Precondtioner 
  int    GMRES_Block_Preconditioner; 
  int    Jacobian_Order;  
  int    GMRES_ILUK_Level_of_Fill;  
  
  // Default Constructor 
  NKS_Input_Parameters() {
    Maximum_Number_of_NKS_Iterations = 0;
    Overall_Tolerance = 1e-5;      
 
    Finite_Time_Step = false;
    Finite_Time_Step_Initial_CFL = 1.0; 
    Finite_Time_Step_Max_CFL = 1.0e12;
     
    Maximum_Number_of_GMRES_Iterations = 0;
    GMRES_Restart = 30;
    GMRES_Overlap = 0;  
    GMRES_Tolerance = 1e-5;
    Normalization = true;
 
    GMRES_CHECK = false;
    GMRES_Frechet_Derivative_Order = FIRST_ORDER; 
    
    GMRES_Block_Preconditioner = Block_Jacobi;
    Jacobian_Order = FIRST_ORDER_INVISCID_HLLE;
    GMRES_ILUK_Level_of_Fill = 0;   
    
  };

  void Output();
  int Parse_Next_Input_Control_Parameter(char *code, char *value);
  void Memory_Estimates(const int &, const int &, const int &);
  
  // Using Default Destructor 
  //~NKS_Input_Parameters() {};

  // Broadcast 
  void Broadcast_Input_Parameters() {

#ifdef _MPI_VERSION
    //Newton 
    MPI::COMM_WORLD.Bcast(&(Maximum_Number_of_NKS_Iterations), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Overall_Tolerance), 
			  1, 
                          MPI::DOUBLE, 0);
 
    // Finite Time Step
    MPI::COMM_WORLD.Bcast(&(Finite_Time_Step), 
			  1, 
			  MPI::INT, 0);  //BOOL
    MPI::COMM_WORLD.Bcast(&(Finite_Time_Step_Initial_CFL), 
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
    MPI::COMM_WORLD.Bcast(&(GMRES_Tolerance), 
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
#endif
  };


#ifdef _MPI_VERSION
  void Broadcast_Input_Parameters(MPI::Intracomm &Communicator,
					    const int Source_Rank){
    //Newton 0
    Communicator.Bcast(&(Maximum_Number_of_NKS_Iterations), 
                          1, 
                          MPI::INT, Source_Rank);
    Communicator.Bcast(&(Overall_Tolerance), 
			  1, 
                          MPI::DOUBLE, Source_Rank);
 
    // Finite Time Step
    Communicator.Bcast(&(Finite_Time_Step), 
			  1, 
			  MPI::INT, Source_Rank);  //BOOL
    Communicator.Bcast(&(Finite_Time_Step_Initial_CFL), 
			  1, 
			  MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Finite_Time_Step_Max_CFL), 
			  1, 
			  MPI::DOUBLE, Source_Rank);

    // GMRES
    Communicator.Bcast(&(Maximum_Number_of_GMRES_Iterations), 
                          1, 
                          MPI::INT, Source_Rank);
    Communicator.Bcast(&(GMRES_Restart), 
                          1, 
                          MPI::INT, Source_Rank); 
    Communicator.Bcast(&(GMRES_Overlap), 
			  1, 
			  MPI::INT, Source_Rank); 
    Communicator.Bcast(&(GMRES_Tolerance), 
			  1, 
                          MPI::DOUBLE, Source_Rank); 
    Communicator.Bcast(&(Normalization), 
			  1, 
			  MPI::INT, Source_Rank); //BOOL
    Communicator.Bcast(&(GMRES_CHECK),
			  1, 
			  MPI::INT, Source_Rank); //BOOL
    Communicator.Bcast(&(GMRES_Frechet_Derivative_Order),
			  1, 
			  MPI::INT, Source_Rank);

    // GMRES Preconditioner Type
    Communicator.Bcast(&(GMRES_Block_Preconditioner), 
			  1, 
			  MPI::INT, Source_Rank);
    //Preconditoner Jacobian Order
    Communicator.Bcast(&(Jacobian_Order), 
                          1, 
                          MPI::INT, Source_Rank);    
    // GMRES ILUK_Level_of_Fill
    Communicator.Bcast(&(GMRES_ILUK_Level_of_Fill), 
			  1, 
			  MPI::INT, Source_Rank);
  };
#endif

  // IO operators  
  friend ostream &operator << (ostream &out_file,
		               const NKS_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       NKS_Input_Parameters &IP);

};


// Returns:
//  - INVALID_INPUT_VALUE if code is valid but value is invalid
//  - INVALID_INPUT_CODE  if unknown code
inline int NKS_Input_Parameters::
Parse_Next_Input_Control_Parameter(char *code, char *value)
{
  int i_command = INVALID_INPUT_CODE;
  char *ptr = NULL;

  // NEWTON 
  if (strcmp(code, "Maximum_Number_of_NKS_Iterations") == 0) {
    i_command = 64;
    Maximum_Number_of_NKS_Iterations = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }
  
  } else if (strcmp(code, "NKS_Overall_Tolerance") == 0) {
    i_command = 61;
    Overall_Tolerance = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }
 

    // FINITE TIME
  } else if (strcmp(code, "NKS_Finite_Time_Step") == 0) {
    i_command = 66; 
//     if (strlen(value) > 1) {
//       for (unsigned ii = 0; ii < strlen(value); ii++) {
// 	value[ii] = tolower(value[ii]);
//       }
//     }
    if (strcmp(value, "OFF") == 0 || strcmp(value, "0") == 0) {
      Finite_Time_Step = false;
    } else {
      Finite_Time_Step = true;
    }

  } else if (strcmp(code, "NKS_Finite_Time_Step_Initial_CFL") == 0) {
    i_command = 70;
    Finite_Time_Step_Initial_CFL = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "NKS_Finite_Time_Step_Max_CFL") == 0) {
    i_command = 71;
    Finite_Time_Step_Max_CFL = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }


    // GMRES 
  } else if (strcmp(code, "Maximum_Number_of_GMRES_Iterations") == 0) {
    i_command = 65;
    Maximum_Number_of_GMRES_Iterations = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "GMRES_Restart") == 0) {
    i_command = 58;
    GMRES_Restart = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "GMRES_Overlap") == 0) {
    i_command = 59;
    GMRES_Overlap = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "GMRES_Tolerance") == 0) {
    i_command = 60;
    GMRES_Tolerance = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "GMRES_Normalization") == 0) {
    i_command = 67;
//     if (strlen(value) > 1) {
//       for (unsigned ii = 0; ii < strlen(value); ii++) {
// 	value[ii] = tolower(value[ii]);
//       }
//     }
    if (strcmp(value, "OFF") == 0 || strcmp(value, "0") == 0) {
      Normalization = false;
    } else {
      Normalization = true;
    }

  } else if (strcmp(code, "GMRES_Check") == 0) {
    i_command = 67;
  //   if (strlen(value) > 1) {
//       for (unsigned ii = 0; ii < strlen(value); ii++) {
// 	value[ii] = tolower(value[ii]);
//       }
//     }
    if (strcmp(value, "ON") == 0 || strcmp(value, "1") == 0) {
      GMRES_CHECK = true;
    } else {
      GMRES_CHECK = false;
    }

  } else if (strcmp(code, "GMRES_Frechet_Derivative_Order") == 0) {
    i_command = 62;
    if (strcmp(value, "First_Order") == 0) {
      GMRES_Frechet_Derivative_Order = FIRST_ORDER;
    } else if(strcmp(value, "Second_Order") == 0) {
      GMRES_Frechet_Derivative_Order = SECOND_ORDER;
    } else {
      GMRES_Frechet_Derivative_Order = FIRST_ORDER;
    }

  } else if (strcmp(code, "GMRES_Block_Preconditioner") == 0) {
    i_command = 62;
    if (strcmp(value, "Diagonal") == 0) {
      GMRES_Block_Preconditioner = Block_Jacobi;
    } else if(strcmp(value, "ILUK") == 0) {
      GMRES_Block_Preconditioner = Block_ILUK;
    } else {
      GMRES_Block_Preconditioner = Block_ILUK;
    }

  } else if (strcmp(code, "GMRES_ILUK_Level_of_Fill") == 0) {
    i_command = 63;
    GMRES_ILUK_Level_of_Fill = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "Jacobian_Order") == 0) {
    i_command = 631;
    if (strcmp(value, "Source_Terms_Only") == 0) {
      Jacobian_Order = SOURCE_TERMS_ONLY;
    } else if (strcmp(value, "First_Order_Inviscid_HLLE") == 0) {
      Jacobian_Order = FIRST_ORDER_INVISCID_HLLE;
    } else if (strcmp(value, "First_Order_Inviscid_Roe") == 0) {
      Jacobian_Order = FIRST_ORDER_INVISCID_ROE; 
    } else if (strcmp(value, "First_Order_Inviscid_AUSM_plus_up") == 0) {
      Jacobian_Order = FIRST_ORDER_INVISCID_AUSMPLUSUP;   
    } else if (strcmp(value, "Second_Order_Diamond_Path_with_HLLE") == 0) {
      Jacobian_Order = SECOND_ORDER_DIAMOND_WITH_HLLE;
    } else if (strcmp(value, "Second_Order_Diamond_Path_with_Roe") == 0) {
      Jacobian_Order = SECOND_ORDER_DIAMOND_WITH_ROE;	
    } else if (strcmp(value, "Second_Order_Diamond_Path_with_AUSM_plus_up") == 0) {
      Jacobian_Order = SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP;
    } else {
      Jacobian_Order = FIRST_ORDER_INVISCID_HLLE;
    }

  } else {
    i_command = INVALID_INPUT_CODE;
  }
  return i_command;
  
}

/***************************************************************
 * NKS_Input_Parameters -- Display Output Operator             *
 ***************************************************************/
inline void NKS_Input_Parameters::Output(){
  //cout.setf(ios::scientific);
  cout.unsetf(ios::scientific);
  cout << " " << endl;
  for (int star=0;star<75;star++){cout <<"*";}
  cout << "\n********                   Newton-Krylov-Schwarz                 **********" << endl;   
  for (int star=0;star<75;star++){cout<<"*";}
  
  cout <<"\n Overall Tolerance     ====> " << Overall_Tolerance << endl;
  
  if (Finite_Time_Step == ON) {     
    cout <<" Finite Time Step      ====> ON" << endl;
    cout <<" Initial_CFL           ====> " << Finite_Time_Step_Initial_CFL << endl;
    cout <<" Max_CFL               ====> " << Finite_Time_Step_Max_CFL << endl;
  } else {
    cout <<" Finite Time Step      ====> OFF" << endl; 
  } /* endif */ 
  
  // GMRES 
  cout <<" Maximum GMRES Its.    ====> " << Maximum_Number_of_GMRES_Iterations<< endl;
  cout <<" GMRES Restart Its.    ====> " << GMRES_Restart << endl;
  cout <<" Level of Overlap      ====> " << GMRES_Overlap << endl;
  cout <<" GMRES Tolerance       ====> " << GMRES_Tolerance << endl;
  
  if (Normalization == ON) {
    cout <<" Normalization         ====> ON" << endl;
  } else {
    cout <<" Normalization         ====> OFF" << endl; 
  } 
  
  if (GMRES_CHECK) {
    cout <<" GMRES Check           ====> ON" << endl;
  }
  if( GMRES_Frechet_Derivative_Order == FIRST_ORDER) {
    cout<<   " Frechet Derivative    ====> First Order "<< endl;
  } else if ( GMRES_Frechet_Derivative_Order == SECOND_ORDER) {
    cout<<   " Frechet Derivative    ====> Second Order "<< endl;
  }
  
  // Precondtioner
  cout<< " Approximate Jacobian  ====>";
  if(Jacobian_Order == SOURCE_TERMS_ONLY) {
    cout<<" Source Terms Only "<<endl;
  } else if(Jacobian_Order == FIRST_ORDER_INVISCID_HLLE){
   cout<<" First Order Inviscid HLLE"<<endl;
  } else if (Jacobian_Order == FIRST_ORDER_INVISCID_ROE){
    cout<<" First Order Inviscid Roe"<<endl;
  } else if (Jacobian_Order == FIRST_ORDER_INVISCID_AUSMPLUSUP){
    cout<<" First Order Inviscid AUSM plus up"<<endl;   
  } else if (Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_HLLE){
    cout<<" Second Order Diamond Path with HLLE"<<endl;
  } else if (Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_ROE){
    cout<<" Second Order Diamond Path with Roe"<<endl;   
  } else if (Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP){
    cout<<" Second Order Diamond Path with AUSM plus up"<<endl;   
  }
  
  if (GMRES_Block_Preconditioner == Block_ILUK) {        
    cout << " Local Preconditioner  ====> ILU("<< GMRES_ILUK_Level_of_Fill <<")" << endl;
  } else if (GMRES_Block_Preconditioner == Block_Jacobi){   // Diagonal
    cout << " Local Preconditioner  ====> Diagonal" << endl; 
  } 
  
  //End 
  for (int star=0;star<75;star++){cout<<"*";}
  cout << endl;
}

/***************************************************************
 * NKS_Input_Parameters -- Memory Estimates                    *
 * These are NOT!! right yet, especially the PRECON            *
 ***************************************************************/
inline void NKS_Input_Parameters::Memory_Estimates(const int &blocksize, 
						   const int &block_mat_size,
						   const int &used_blocks){

  cerr<<" \n NKS_Input_Parameters::Memory_Estimates not working yet ";
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


/***************************************************************
 * NKS_Input_Parameters -- Input-output operators.       *
 ***************************************************************/
inline ostream &operator << (ostream &out_file,
			     const NKS_Input_Parameters &IP) {    
    return (out_file);
}


inline istream &operator >> (istream &in_file,
			     NKS_Input_Parameters &IP) {
    return (in_file);
}

#endif // _NKS2DINPUT_INCLUDED
