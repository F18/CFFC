/* CompDomain2D.cc: Subroutines for CompDomain2D Classes*/

/* Include the header files */

#ifndef _COMPDOMAIN2D_INCLUDED
#include "CompDomain2D.h"
#endif // _COMPDOMAIN2D_INCLUDED

#ifndef _TESTFUNCTIONS_2D_INCLUDED
#include "TestFunctions_2D.h"
#endif // _TESTFUNCTIONS_2D_INCLUDED

#ifndef _CFD_INCLUDED
#include "CFD.h"
#endif // _CFD_INCLUDED

#include <cstring>
#include <cstdlib>
using namespace std;

/*************************************************************
 * _COMPDOMAIN2D_INCLUDED  -- External subroutines. *
 *************************************************************/
//*****
/*************************************************************
 ******  CompDomain2D_Cartesian -- External subrutines  ********
 ************************************************************/

/********************************************************
 * Routine: Generate_Mesh for Cartesian_Grid            *
 *                                                      *
 * Generates the mesh.                                  *
 *                                                      *
 ********************************************************/

CompDomain2D_Cartesian* Generate_Mesh(CompDomain2D_Cartesian *Grid,
				      Reconstruct2D_Input_Parameters &IP){

  int NGC; // the number of gostcells
  // it takes into account the reconstruction order
  int order; //reconstruction order
  int N; // the minimum number of cells necessary in the reconstruction 
  int error_flag;

  Grid = new CompDomain2D_Cartesian;

  order = IP.Reconstruction_Order;

  N = (order+1)*(order+2)/2;
  if (N <= 1)
    NGC = 0; // !!!! check for N<=3
  else if (N <= 12)
    NGC = 2; 
  else if (N <= 49)
    NGC = 3;
  else if (N <= 81)
    NGC = 4;
  else {
    cout << "\nThe program doesn't perform reconstruction for order "
	 << order;
    return (NULL);
  }

  // allocate memory for the grid
  Grid->allocate(IP.Number_of_Cells_Idir, IP.Number_of_Cells_Jdir, NGC);

  // generate the geometry of the mesh
  error_flag = SetGeometry(Grid, IP);
  if (error_flag) {
    Grid->deallocate();
    return NULL;
  }

  // initialize the cell average quantity for the specific function
  Compute_Cell_Average(Grid, IP);

  return Grid;

}

/********************************************************
 * Routine: SetGeometry for Cartesian_Grid              *
 *                                                      *
 * Sets the location of the cell centers, the values of *
 * the cell lengths and the locations of the subgrid    * 
 * points.                                              * 
 *                                                      *
 ********************************************************/

int SetGeometry(CompDomain2D_Cartesian *Grid,
		 Reconstruct2D_Input_Parameters &IP){

  //local variable
  Vector2D loc, delta;  // loc -- location of the first cell center
                        // delta -- delta_X and delta_Y of the cell
  Vector2D new_loc;    // new_loc -- the next cell center location
  int NbSubCell_x;      // Number of SubCell Points in Idir
  int NbSubCell_y;      // Number of SubCell Points in Jdir
  int ro;             // Reconstruction Order

  // initialize the parameters for each CompCell

  switch(IP.i_Grid){

  case GRID_READ_FROM_GRID_DATA_FILE:
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cout << "\n Grid_Read_From_Data_File subrutine wasn't implemented yet\a"
	 << endl;
    cout.flush();
    return 1;
    break;
  case GRID_READ_FROM_DEFINITION_FILE:
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cout << "\n Grid_Read_From_Definition_File subrutine wasn't implemented yet\a"
	 << endl;
    cout.flush();
    return 1;
    break;
  default:
    Vector2D temp(IP.Delta_X_Cell, IP.Delta_Y_Cell);
    delta = temp;
    Vector2D P_min(IP.X_min, IP.Y_min); // minim point
    Vector2D Translate(0.5-Grid->ICl, 0.5-Grid->JCl); //translation vector
    Translate.x = delta.x * Translate.x;
    Translate.y = delta.y * Translate.y;
    loc = P_min + Translate;
    NbSubCell_x = IP.Number_of_SubGrid_Points_Idir;
    NbSubCell_y = IP.Number_of_SubGrid_Points_Jdir;
    ro = IP.Reconstruction_Order;

    // initialize the geometry
    for (int i = 0; i<=Grid->NCi-1; i++){
      Translate.x = (double) i * delta.x;
      for (int j = 0; j<=Grid->NCj-1; j++){
	Translate.y = (double) j * delta.y;
	new_loc = loc + Translate;
	Grid->Mesh[i][j].SetCellParameters(new_loc, delta, NbSubCell_x,
					   NbSubCell_y, ro);
      }
    }
  } //endswitch

  return 0;
}

/********************************************************
 * Routine: Compute_Cell_Average for Cartesian_Grid     *
 *                                                      *
 * Computes the cell average quantity for the specified *
 * function.                                            * 
 ********************************************************/

void Compute_Cell_Average(CompDomain2D_Cartesian *Grid,
			  Reconstruct2D_Input_Parameters &IP){

  double F_Integral;
  Node2D SW, SE, NW, NE;
  double area;
  int index_i, index_j;
  int NumericalIntegration = ON;

  if (NumericalIntegration)
    cout << endl
	 << "********* Numerical integration *********" << endl;
  else
    cout << endl
	 << "********* Analitical integration *********" << endl;

  for (int i = 0; i<=Grid->NCi-1; i++)
    for (int j = 0; j<=Grid->NCj-1; j++){

      index_i = Grid->Mesh[i][j].Nx-1;
      index_j = Grid->Mesh[i][j].Ny-1;
      SW.setloc(Grid->Mesh[i][j].SubGrid[0][0]);
      SE.setloc(Grid->Mesh[i][j].SubGrid[index_i][0]);
      NW.setloc(Grid->Mesh[i][j].SubGrid[0][index_j]);
      NE.setloc(Grid->Mesh[i][j].SubGrid[index_i][index_j]);
  
      area = Grid->Mesh[i][j].geom.dx.x * Grid->Mesh[i][j].geom.dx.y;
      if (NumericalIntegration){
	F_Integral = quad2dAdaptiveGaussianQuadrature(IP.TestF,
						      SW.X.x, SE.X.x,
						      SW.X.y, NW.X.y,
						      11);
      }
      else{
	F_Integral = IP.IntTestF(SW.X.x, SE.X.x, SW.X.y, NW.X.y);
      }
      Grid->Mesh[i][j].U_cell = F_Integral/area;
    }
}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution values at the cell center for 2D *
 * reconstruction to the specified output data file(s)  *
 * in a format suitable for plotting with TECPLOT.      *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
int Output_Tecplot(CompDomain2D_Cartesian *Grid,
                   Reconstruct2D_Input_Parameters &IP) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256], output_directory[30];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Determine prefix of output data file names.
  
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  
  // Determine output data file name.
  strcpy(output_directory, "../Plots/");
  strcpy(output_file_name, output_directory);
  strcat(output_file_name, prefix);
  strcat(output_file_name, ".dat");
  output_file_name_ptr = output_file_name;
  
  // Open the output data file.
  
  output_file.open(output_file_name_ptr, ios::out);
  if (output_file.bad()) return (1);
  
  // Write the solution data.
  
  output_file << setprecision(14);
  output_file << "TITLE = \"" << IP.Reconstruction_Order;
  switch (IP.Reconstruction_Order){
    case 1:
      output_file << "st";
      break;
    case 2:
      output_file << "nd";
      break;
     case 3:
      output_file << "rd";
      break;
    default:
      output_file << "th";
    }
    output_file <<" Order 2D Reconstruction Solution for function "
                << IP.Function_Type << "\"" << "\n"
	        << "VARIABLES =" 
		<< "\n\"x\""
		<< "\n \"y\""
		<< "\n \"U\" \n";

//      for ( i = 0; i<= Grid->NCi-1 ; i++){
//         for (int j = 0; j<= Grid->NCj-1; j++){
                
    for ( i = Grid->ICl; i<= Grid->ICu ; i++){
      for (int j = Grid->JCl; j<= Grid->JCu; j++){
	Grid->Mesh[i][j].Output_Tecplot_Cell_Solution(output_file, i, j);
	output_file << endl;
      } // endfor
    } // endfor
    output_file << setprecision(6);
    
    // Close the output data file.
    
    output_file.close();
    
    // Writing of output data files complete.  Return zero value.
    
    return(0);
}


/********************************************************
 * Routine: Output_Derivatives_Tecplot                  *
 *                                                      *
 * Writes the derivative values at the subnodes of each *
 * cell for 2D reconstruction to the specified output   *
 * data file(s) in a format suitable for plotting with  *
 * TECPLOT.                                             *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
int Output_Derivatives_Tecplot(CompDomain2D_Cartesian *Grid,
			       Reconstruct2D_Input_Parameters &IP) {
  return 0;
}

/********************************************************
 * Routine: Output_Cell_Center_Values                   *
 *                                                      *
 * Writes the theoretical value of the function at the  *
 * cell center and the numerical one for                *
 * 1D reconstruction to the specified output            *
 * data file(s) in ASCII format.                        *
 * This file can be also imported in TECPLOT            *
 * Returns a non-zero value if cannot write any of the  *
 * solution files.                                      *
 *                                                      *
 ********************************************************/

int Output_Cell_Center_Values(CompDomain2D_Cartesian *Grid,
			      Reconstruct2D_Input_Parameters &IP){

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256], output_directory[30];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Determine prefix of output data file names.
  
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  
  // Determine output data file name.
  strcpy(output_directory, "../Plots/");
  strcpy(output_file_name, output_directory);
  strcat(output_file_name, prefix);
  strcat(output_file_name, "_CellCenter.dat");
  output_file_name_ptr = output_file_name;
  
  // Open the output data file.
  
  output_file.open(output_file_name_ptr, ios::out);
  if (output_file.bad()) return (1);
  
  // Write the solution data.
  
  output_file << setprecision(14);
  output_file << "TITLE = \"" << IP.Reconstruction_Order;
  switch (IP.Reconstruction_Order){
    case 1:
      output_file << "st";
      break;
    case 2:
      output_file << "nd";
      break;
     case 3:
      output_file << "rd";
      break;
    default:
      output_file << "th";
    }

    output_file <<" Order 2D Reconstruction Solution for function "
                << IP.Function_Type << "\"" << "\n"
	        << "VARIABLES =" 
		<< "\n\"x\""
		<< "\n \"y\""
		<< "\n \"U\" \n";

    output_file << "ZONE I=" 
		<< IP.Number_of_Cells_Idir
		<< ",  J="
		<< IP.Number_of_Cells_Jdir
		<< ",   F=Point\n"
		<< "T= CellCenter\n";
    
    for (int j = Grid->JCl; j<= Grid->JCu; j++){
      for (  i = Grid->ICl; i<= Grid->ICu; i++){
	output_file << Grid->Mesh[i][j].geom.xc.x << "\t" 
		    << Grid->Mesh[i][j].geom.xc.y << "\t"
		    << Grid->Mesh[i][j].U_cell << endl;
      } // endfor
    } // endfor
    output_file << setprecision(6);
    
    // Close the output data file.
    
    output_file.close();
    
    // Writing of output data files complete.  Return zero value.
    
    return(0);
  
}  

/********************************************************
 * Routine: Output_Function_Graph_Tecplot               *
 *                                                      *
 * Writes the representation of the function using its  *
 * mathematical definition to the specified output data *
 * file(s) in a format suitable for plotting            *
 * with TECPLOT.                                        *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
int Output_Function_Graph_Tecplot(CompDomain2D_Cartesian *Grid,
				  Reconstruct2D_Input_Parameters &IP) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256], output_directory[30];
  char *output_file_name_ptr;
  ofstream output_file;    
  
  // Determine prefix of output data file names.
  
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  
  // Determine output data file name.
  strcpy(output_directory, "../Plots/");
  strcpy(extension,"_Function_");
  strcat(extension,IP.Function_Type);
  strcpy(output_file_name, output_directory);
  strcat(output_file_name, prefix);
  strcat(output_file_name, extension);
  strcat(output_file_name, ".dat");
  output_file_name_ptr = output_file_name;
  
  // Open the output data file.
  
  output_file.open(output_file_name_ptr, ios::out);
  if (output_file.bad()) return (1);
  
  // Write the solution data.
  
  output_file << setprecision(14);
  output_file << "TITLE = \"" << IP.Function_Type
	      << " representation on the domain: ";
  output_file << IP.X_min << " :  " << IP.X_min + IP.X_Length
	      << " -- " << IP.Y_min << " :  " << IP.Y_min + IP.Y_Length 
	      <<"\""<< endl
	      << "VARIABLES = \"x\"  \"y\"  \"U\" \n";

  
  int N_standard = 100;
  // the number of points used for representing a domain length equal to 1
  int N_points_x = (int) (IP.Delta_X_Cell*N_standard);
  int N_points_y = (int) (IP.Delta_Y_Cell*N_standard);
  
  if (N_points_x == 0)
    N_points_x = 5;
  if (N_points_y == 0)
    N_points_y = 5;

  Vector2D FunctionGrid[N_points_x][N_points_y];
  double FunctionValues[N_points_x][N_points_y];
  
  double delta_x, delta_y;
  
  delta_x = IP.Delta_X_Cell/(N_points_x-1);
  delta_y = IP.Delta_Y_Cell/(N_points_y-1);

  for (int i_cell = Grid->ICl; i_cell<=Grid->ICu; i_cell++)
    for (int j_cell = Grid->JCl; j_cell<= Grid->JCu; j_cell++){  // for each cell

      // create grid for exact function representation
      FunctionGrid[0][0].x = Grid->Mesh[i_cell][j_cell].SubGrid[0][0].x;
      FunctionGrid[0][0].y = Grid->Mesh[i_cell][j_cell].SubGrid[0][0].y;

      
      for (i=0; i<=N_points_x-1; i++)
	for (int j=0; j<=N_points_y-1; j++){
	  FunctionGrid[i][j].x = FunctionGrid[0][0].x + i*delta_x;
	  FunctionGrid[i][j].y = FunctionGrid[0][0].y + j*delta_y;
	}
  
      // compute function values in the grid points
      for (i=0; i<=N_points_x-1; i++)
	for (int j=0; j<=N_points_y-1; j++){
	  FunctionValues[i][j]= IP.TestF(FunctionGrid[i][j].x,
					 FunctionGrid[i][j].y);
	  
	}//endfor
  
      // write the data to the output data file 
      output_file << "ZONE I=" 
		  << N_points_x
		  << ",  J="
		  << N_points_y
		  << ",   F=Point\n"
		  << "T=\"Cell [" << i_cell <<"," << j_cell << "]\"\n";
      
      output_file.setf(ios::scientific);
      
      for (int j=0; j<= N_points_y-1; j++)      
	for ( i=0; i<= N_points_x-1; i++ ){ 
	  output_file << FunctionGrid[i][j]<< "\t" 
		      << FunctionValues[i][j] << endl;
	}
      output_file.unsetf(ios::scientific);
  
    } // endfor

  output_file << setprecision(6);
  
  // Close the output data file.
  
  output_file.close();
  
  // Writing of output data files complete.  Return zero value.
  
  return(0);
}

/********************************************************
 * Routine: Output_Matlab                               *
 *                                                      *
 * Writes the solution values at the cell centers for 2D*
 * reconstruction to the specified output data file(s)  *
 * in a format suitable for plotting with MATLAB.       *
 * Returns a non-zero value if cannot write any of the  *
 * MATLAB solution files.                               *
 *                                                      *
 ********************************************************/
int Output_Matlab(CompDomain2D_Cartesian *Grid,
                  Reconstruct2D_Input_Parameters &IP) {

  return 0;
}

/********************************************************
 * Routine: Output_Derivatives_Matlab                   *
 *                                                      *
 * Writes the derivative values at the nodes of each    *
 * cell for 2D reconstruction to the specified output   *
 * data file(s) in a format suitable for plotting with  *
 * MATLAB.                                              *
 * Returns a non-zero value if cannot write any of the  *
 * MATLAB solution files.                               *
 *                                                      *
 ********************************************************/

int Output_Derivatives_Matlab(CompDomain2D_Cartesian *Grid,
			      Reconstruct2D_Input_Parameters &IP) {
  return 0;  
}

/********************************************************
 * Routine: Output_Function_Graph_Matlab                *
 *                                                      *
 * Writes the representation of the function using its  *
 * mathematical definition to the specified output data *
 * file(s) in a format suitable for plotting            *
 * with MATLAB.                                         *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/

int Output_Function_Graph_Matlab(Reconstruct2D_Input_Parameters &IP) {
 
  return 0;
}

double L1_Norm(CompDomain2D_Cartesian *Grid,
	       Reconstruct2D_Input_Parameters &IP){
/* Compute the L1 norm for the whole domain*/

  double L1_norm=0.0;

  for (int i=Grid->ICl; i<=Grid->ICu; i++){
    for (int j=Grid->JCl; j<=Grid->JCu; j++)
      L1_norm += Grid->Mesh[i][j].ErrorCellL1Norm(IP.TestF);
  }
  return L1_norm/(IP.Number_of_Cells_Idir*IP.Number_of_Cells_Jdir*
                  Grid->Mesh[1][1].geom.dx.x*Grid->Mesh[1][1].geom.dx.y);
}


double L2_Norm(CompDomain2D_Cartesian *Grid,
	       Reconstruct2D_Input_Parameters &IP){
/* Compute the L2 norm for the whole domain*/

  double L2_norm=0.0;

  for (int i=Grid->ICl; i<=Grid->ICu; i++){
    for (int j=Grid->JCl; j<=Grid->JCu; j++)
      L2_norm += Grid->Mesh[i][j].ErrorCellL2Norm(IP.TestF);
  }
  L2_norm = sqrt(L2_norm/(IP.Number_of_Cells_Idir*IP.Number_of_Cells_Jdir));

  return L2_norm;
}

/********************************************************
 * Routine: Output_L1_Norm_Tecplot                      *
 *                                                      *
 * Writes the L1_Norm of the reconstruction             *
 * to the specified output data                         *
 * file(s) in a format suitable for plotting            *
 * with TECPLOT.                                        *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 * The first call of the function writes the header file*
 * and the next calls only append to the end of the file*
 *                                                      *
 ********************************************************/

static int WrittenHeader = OFF;

int Output_L1_Norm_Tecplot(CompDomain2D_Cartesian *Grid,
			   Reconstruct2D_Input_Parameters &IP){

     int i, i_output_title;
    char prefix[256], extension[256], output_file_name[256], output_directory[30];
    char order[5];
    char *output_file_name_ptr;
    ofstream output_file;    
    double L1_norm;

    /* Determine prefix of output data file names. */

    i = 0;
    while (1) {
       if (IP.Output_File_Name[i] == ' ' ||
           IP.Output_File_Name[i] == '.') break;
       prefix[i]=IP.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(IP.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';

    /* Determine output data file name. */
    switch(IP.Reconstruction_Order){
    case 1:
      strcpy(order,"1");
      break;
    case 2:
      strcpy(order,"2");
      break;
    case 3:
      strcpy(order,"3");
      break;
    case 4:
      strcpy(order,"4");
      break;
    case 5:
      strcpy(order,"5");
      break;
    case 6:
      strcpy(order,"6");
      break;
    case 7:
      strcpy(order,"7");
      break;
    default:
      strcpy(order,"0");
    }
    strcpy(output_directory, "../Plots/");
    strcpy(extension,"_Function_");
    strcat(extension,IP.Function_Type);
    strcat(extension,"_L1Norm_Order");
    strcat(extension,order);
    strcat(extension,"_");
    strcat(extension,IP.Data_Dependent_Weighting);
    strcpy(output_file_name, output_directory);
    strcat(output_file_name, prefix);
    strcat(output_file_name, extension);
    strcat(output_file_name, ".dat");
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    if (!WrittenHeader){
      output_file.open(output_file_name_ptr, ios::trunc);
      if (output_file.bad()) return (1);
      // write the header
      output_file << "TITLE = \"" << IP.Function_Type
		  << " L1Norm representation on the domain: ";
      output_file << IP.X_min << " :  " << IP.X_min + IP.X_Length
		  << " -- " << IP.Y_min << " :  " << IP.Y_min + IP.Y_Length 
		  <<"\""<< endl
		  << "VARIABLES = \"#Points\" \"L1Norm\" \n"
		  << "ZONE \n";
      WrittenHeader = ON;
    }
    else{
      output_file.open(output_file_name_ptr, ios::app);
      if (output_file.bad()) return (1);
    }
    
    /* Write the solution data. */
    L1_norm = L1_Norm(Grid,IP);

    output_file << setprecision(14);
    output_file << " " 
		<< IP.Number_of_Cells_Idir*IP.Number_of_Cells_Jdir
		<< "  " << L1_norm << "\n";
    output_file.unsetf(ios::scientific);
    output_file << setprecision(6);

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);

}

static int WrittenHeaderL2 = OFF;

int Output_L2_Norm_Tecplot(CompDomain2D_Cartesian *Grid,
			   Reconstruct2D_Input_Parameters &IP){

     int i, i_output_title;
    char prefix[256], extension[256], output_file_name[256], output_directory[30];
    char order[5];
    char *output_file_name_ptr;
    ofstream output_file;    
    double L2_norm;

    /* Determine prefix of output data file names. */

    i = 0;
    while (1) {
       if (IP.Output_File_Name[i] == ' ' ||
           IP.Output_File_Name[i] == '.') break;
       prefix[i]=IP.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(IP.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';

    /* Determine output data file name. */
    switch(IP.Reconstruction_Order){
    case 1:
      strcpy(order,"1");
      break;
    case 2:
      strcpy(order,"2");
      break;
    case 3:
      strcpy(order,"3");
      break;
    case 4:
      strcpy(order,"4");
      break;
    case 5:
      strcpy(order,"5");
      break;
    case 6:
      strcpy(order,"6");
      break;
    case 7:
      strcpy(order,"7");
      break;
    default:
      strcpy(order,"0");
    }

    strcpy(output_directory, "../Plots/");
    strcpy(extension,"_Function_");
    strcat(extension,IP.Function_Type);
    strcat(extension,"_L2Norm_Order");
    strcat(extension,order);
    strcat(extension,"_");
    strcat(extension,IP.Data_Dependent_Weighting);
    strcpy(output_file_name, output_directory);
    strcat(output_file_name, prefix);
    strcat(output_file_name, extension);
    strcat(output_file_name, ".dat");
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    if (!WrittenHeaderL2){
      output_file.open(output_file_name_ptr, ios::trunc);
      if (output_file.bad()) return (1);
      // write the header
      output_file << "TITLE = \"" << IP.Function_Type
		  << " L2Norm representation on the domain: ";
      output_file << IP.X_min << " :  " << IP.X_min + IP.X_Length
		  << " -- " << IP.Y_min << " :  " << IP.Y_min + IP.Y_Length 
		  <<"\""<< endl
		  << "VARIABLES = \"#Points\" \"L2Norm\" \n"
		  << "ZONE \n";
      WrittenHeaderL2 = ON;
    }
    else{
      output_file.open(output_file_name_ptr, ios::app);
      if (output_file.bad()) return (1);
    }
    
    /* Write the solution data. */
    L2_norm = L2_Norm(Grid,IP);

    output_file << setprecision(14);
    output_file << " " 
		<< IP.Number_of_Cells_Idir*IP.Number_of_Cells_Jdir
		<< "  " << L2_norm << "\n";
    output_file.unsetf(ios::scientific);
    output_file << setprecision(6);

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);

}
