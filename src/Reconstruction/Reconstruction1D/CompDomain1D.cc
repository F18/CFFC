/* CompDomain1D.cc: Subroutines for CompDomain1D Classes*/

/* Include the header files */

#ifndef _COMPDOMAIN1D_INCLUDED
#include "CompDomain1D.h"
#endif // _COMPDOMAIN1D_INCLUDED

/*************************************************************
 * _COMPDOMAIN1D_INCLUDED  -- External subroutines. *
 *************************************************************/

/*************************************************************
 ******  CompDomain1D_Uniform -- External subrutines  ********
 ************************************************************/

/********************************************************
 * Routine: Generate_Mesh for Uniform_Grid              *
 *                                                      *
 * Generates the mesh.                                  *
 *                                                      *
 ********************************************************/

void Generate_Mesh(CompDomain1D_Uniform &Grid,
		   const Reconstruct1D_Input_Parameters &IP){

  int NC; // the number of cells for initializing the computational domain
  // it takes into account the reconstruction order

  if (IP.Reconstruction_Order%2 == 0)
    Grid.NCS = (IP.Reconstruction_Order)/2;
  else
    Grid.NCS = (IP.Reconstruction_Order + 1)/2;

  NC = IP.Number_of_Cells_Idir + 2*Grid.NCS;

  // allocate memory for the grid
  Grid.allocate(NC);

  // generate the geometry of the mesh
  SetGeometry(Grid, IP);

  // initialize the cell average quantity for the specific function
  Compute_Cell_Average(Grid, IP);
}

/********************************************************
 * Routine: SetGeometry for Uniform_Grid                *
 *                                                      *
 * Sets the location of the cell centers, the values of *
 * the cell lengths and the locations of the subgrid    * 
 * points.                                              * 
 *                                                      *
 ********************************************************/

void SetGeometry(CompDomain1D_Uniform &Grid,
		 const Reconstruct1D_Input_Parameters &IP){

  //local variable
  double loc, delta;  // loc -- location of the first cell center
                      // delta -- delta_X of the cell
  double new_loc;     // new_loc -- the next cell center location
  int NbSubCell;      // Number of SubCell Points
  int ro;             // Reconstruction Order

  // initialize the parameters for each CompCell
  delta = IP.Delta_Cell;
  loc = new_loc = IP.X_min+delta*(0.5-Grid.NCS);
  NbSubCell = IP.Number_of_SubGrid_Points;
  ro = IP.Reconstruction_Order;

  for (int i = 0; i<=Grid.NCi-1; i++){

    // initialize the geometry
    Grid.Mesh[i].SetCellParameters(new_loc, delta, NbSubCell, ro);
    new_loc = loc + (i+1)*delta;
  }
}

/********************************************************
 * Routine: Compute_Cell_Average for Uniform_Grid       *
 *                                                      *
 * Computes the cell average quantity for the specified *
 * function.                                            * 
 ********************************************************/

void Compute_Cell_Average(CompDomain1D_Uniform &Grid,
			  const Reconstruct1D_Input_Parameters &IP){

  double F_Integral;
  double x_inf, x_sup;
  double length;

  for (int i = 0; i<=Grid.NCi-1; i++){

    x_inf = Grid.Mesh[i].SubGrid[0];
    x_sup = Grid.Mesh[i].SubGrid[Grid.Mesh[i].N-1];
    length = Grid.Mesh[i].geom.dx;
    F_Integral = IP.IntTestF(x_sup, x_inf);
    Grid.Mesh[i].U_cell = 1.0/length*F_Integral;
  }
}

/****************************************************************
 ******  CompDomain1D_NonUniform -- External subrutines  ********
 ***************************************************************/
// Copy constructor 1
CompDomain1D_NonUniform::CompDomain1D_NonUniform(const CompDomain1D_NonUniform &G):
  NCi(G.NCi), NCS(G.NCS){
  Mesh = new CompCell1D_NonUniform[G.NCi];
  if (Mesh == NULL){
    cout << "Not enough memory!\n";
    exit(1);
  }
  for (int i=0; i <= NCi-1; i++){
    Mesh[i] = G.Mesh[i];
  }
}

double CompDomain1D_NonUniform::ComputeValueInPoint(double x){
  double value;
  int flag = OFF;
  for (int i=NCS; i<=NCi-NCS-1; i++){
    if ((Mesh[i].geom.x >= x)&&(Mesh[i].geom.x - Mesh[i].geom.dx/2 <= x)){
      value = Mesh[i].SolutionAtCoordinates(x);
      flag = ON;
      break;
    }
  }
  if (flag==OFF)
    cout << "The x coordinate is out of the computational domain!\n";
  return value;
}

/********************************************************
 * Routine: Generate_Mesh for NonUniform_Grid           *
 *                                                      *
 * Generates the mesh.                                  *
 *                                                      *
 ********************************************************/
CompDomain1D_NonUniform* Generate_Mesh(CompDomain1D_NonUniform *Grid,
				       const Reconstruct1D_Input_Parameters &IP){

  int NC; // the number of cells for initializing the computational domain
  // it takes into account the reconstruction order
  int error_flag;

  Grid = new CompDomain1D_NonUniform;

  if (IP.Reconstruction_Order == 0)
    Grid->NCS = 0;
  else if (IP.Reconstruction_Order%2 == 0){
    Grid->NCS = IP.Reconstruction_Order/2 + 1;
  }
  else
    Grid->NCS = (IP.Reconstruction_Order + 1)/2;

//   if (IP.Reconstruction_Order == 3){
//     Grid->NCS = 3;
//     Print(Grid->NCS)
//   }

  NC = IP.Number_of_Cells_Idir + 2*Grid->NCS;
  // allocate memory for the grid
  Grid->allocate(NC);

  // generate the geometry of the mesh
  error_flag = SetGeometry(Grid, IP);

  if (error_flag) {
    Grid->deallocate();
    return NULL;
  }

  // initialize the cell average quantity for the specific function
  Compute_Cell_Average(Grid, IP);

  // initialize D's
  for (int cell=0; cell<= Grid->NCi-1; cell++){
    Grid->Mesh[cell].D[0]= Grid->Mesh[cell].U_cell;
    Grid->Mesh[cell].SetSubDomainSolution();
  }

  return Grid;
}

/********************************************************
 * Routine: SetGeometry for NonUniform_Grid             *
 *                                                      *
 * Sets the location of the cell centers, the values of *
 * the cell lengths and the locations of the subgrid    * 
 * points.                                              * 
 *                                                      *
 ********************************************************/

int SetGeometry(CompDomain1D_NonUniform *Grid,
		 const Reconstruct1D_Input_Parameters &IP){

  //local variable
  double loc, delta;  // loc -- location of the first cell center
                      // delta -- delta_X of the cell
  double new_loc;     // new_loc -- the next cell center location
  int NbSubCell;      // Number of SubCell Points
  int ro;             // Reconstruction Order

  // initialize the parameters for each CompCell

  switch(IP.i_Grid){

  case GRID_READ_FROM_GRID_DATA_FILE:
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cout << "\n Grid_Read_From_Data_File subrutine wasn't implemented yet\a" << endl;
    cout.flush();
    return 1; break;
  case GRID_READ_FROM_DEFINITION_FILE:
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cout << "\n Grid_Read_From_Definition_File subrutine wasn't implemented yet\a" << endl;
    cout.flush();
    return 1; break;
  case GRID_NONUNIFORM:
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cout << "\n Grid_Nonuniform subrutine wasn't implemented yet\a" << endl;
    cout.flush();
    return 1; break;
  default:
    delta = IP.LengthX()/IP.iCell();
    loc = IP.X_min+delta*(0.5-Grid->NCS);
    NbSubCell = IP.Number_of_SubGrid_Points;
    ro = IP.Reconstruction_Order;
    // initialize the geometry
    for (int i = 0; i<=Grid->NCi-1; i++){
      new_loc = loc + (i)*delta;
      Grid->Mesh[i].SetCellParameters(new_loc, delta, NbSubCell, ro);
    }

  } //endswitch

  return 0;
}

/********************************************************
 * Routine: Compute_Cell_Average for NonUniform_Grid    *
 *                                                      *
 * Computes the cell average quantity for the specified *
 * function.                                            * 
 ********************************************************/
void Compute_Cell_Average(CompDomain1D_NonUniform *Grid,
			  const Reconstruct1D_Input_Parameters &IP){

  double F_Integral;
  double x_inf, x_sup;
  double length;

  for (int i = 0; i<=Grid->NCi-1; ++i){
    x_inf = Grid->Mesh[i].SubGrid[0];
    x_sup = Grid->Mesh[i].SubGrid[Grid->Mesh[i].N-1];
    length = Grid->Mesh[i].geom.dx;

    if (strcmp(IP.Integration_Type,"Numeric")==0)
      F_Integral = AdaptiveGaussianQuadrature(IP.TestF,x_inf, x_sup,14,F_Integral);
    else{
      F_Integral = IP.IntTestF(x_inf, x_sup);
    }

    Grid->Mesh[i].U_cell = F_Integral/length;
  }//endif

  // Find the maximum and minimum solution
  Grid->MaxSolution = Grid->Mesh[0].U_cell;
  Grid->MinSolution = Grid->Mesh[0].U_cell;
  for (int i = 1; i<=Grid->NCi-1; ++i){
    Grid->MaxSolution = max(Grid->MaxSolution, Grid->Mesh[i].U_cell);
    Grid->MinSolution = min(Grid->MinSolution, Grid->Mesh[i].U_cell);
  }

  // Set MaxDeltaSolutionOverDomain
  Grid->MaxDeltaSolutionOverDomain = Grid->MaxSolution - Grid->MinSolution;

}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution values at the nodes for 1D       *
 * reconstruction to the specified output data file(s)  *
 * in a format suitable for plotting with TECPLOT.      *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
int Output_Tecplot(const CompDomain1D_NonUniform *Grid,
                   const Reconstruct1D_Input_Parameters &IP) {

    int i;
    char prefix[256], output_file_name[256], output_directory[30];
    char *output_file_name_ptr;
    ofstream output_file;    

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
    strcpy(output_directory, "../Results/");
    strcpy(output_file_name, output_directory);
    strcat(output_file_name, prefix);
    strcat(output_file_name, ".dat");
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data. */

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
    output_file <<" Order 1D Reconstruction Solution for function "
                << IP.Function_Type << "\"" << "\n"
	        << "VARIABLES = \"x\" \"U\" \n"
                << "ZONE \n";

    int I_min  = Grid->NCS;
    int I_max = Grid->NCi - Grid->NCS - 1;

    for ( i = I_min; i<= I_max; i++ ) {
       for (int j = 0 ; j <= Grid->Mesh[i].N-1 ; j++ ) {

           output_file << " " << Grid->Mesh[i].SubGrid[j] 
		       << "  " << Grid->Mesh[i].U_points[j] << "\n";
           output_file.unsetf(ios::scientific);
       } /* endfor */
    } /* endfor */
    output_file << setprecision(6);

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);
}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution values at the nodes for 1D       *
 * IC data to the specified output data file(s)  *
 * in a format suitable for plotting with TECPLOT.      *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
int Output_Tecplot_Domain(const CompDomain1D_NonUniform *Grid,
			  const Reconstruct1D_Input_Parameters &IP) {


    int i;
    char prefix[256], output_file_name[256], output_directory[30];
    char *output_file_name_ptr;
    ofstream output_file;    

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
    strcpy(output_directory, "../Results/");
    strcpy(output_file_name, output_directory);
    strcat(output_file_name, prefix);
    strcat(output_file_name, "_grid.dat");
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data. */

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
    output_file <<" Order 1D Reconstruction Solution for function "
                << IP.Function_Type << "\"" << "\n"
	        << "VARIABLES = \"x\" \"U\" \n"
                << "ZONE \n";

    int I_min  = 0;
    int I_max = Grid->NCi - 1;

    for ( i = I_min; i<= I_max; i++ ) {
       for (int j = 0 ; j <= Grid->Mesh[i].N-1 ; j++ ) {

           output_file << " " << Grid->Mesh[i].SubGrid[j] 
		       << "  " << Grid->Mesh[i].U_points[j] << "\n";
           output_file.unsetf(ios::scientific);
       } /* endfor */
    } /* endfor */
    output_file << setprecision(6);

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);
}

/********************************************************
 * Routine: Output_Matlab                               *
 *                                                      *
 * Writes the solution values at the nodes for 1D       *
 * reconstruction to the specified output data file(s)  *
 * in a format suitable for plotting with MATLAB.       *
 * Returns a non-zero value if cannot write any of the  *
 * MATLAB solution files.                               *
 *                                                      *
 ********************************************************/
int Output_Matlab(const CompDomain1D_NonUniform *Grid,
                  const Reconstruct1D_Input_Parameters &IP) {

    int i;
    char prefix[256], output_file_name_data[256], output_directory[30];
    char output_file_name_script[256];
    char *output_file_name_ptr;
    ofstream output_file;    

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

    strcpy(output_directory, "../Results/");
    strcpy(output_file_name_data, output_directory);
    strcat(output_file_name_data, prefix);
    strcat(output_file_name_data, "_mtb.dat");
    output_file_name_ptr = output_file_name_data;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data. */

    output_file << setprecision(14);
    output_file << "%TITLE = \"" << IP.Reconstruction_Order;
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
    } //endswitch
    output_file <<" Order 1D Reconstruction Solution for function "
                << IP.Function_Type << "\"" << "\n"
	        << "%VARIABLES \n" 
		<< " x \t U \n";

    int I_min  = Grid->NCS;
    int I_max = Grid->NCi - Grid->NCS - 1;

    for ( i=I_min; i <=I_max; i++ ) {
       for (int j=0 ; j<= Grid->Mesh[i].N-1 ; j++ ) {

           output_file << " " << Grid->Mesh[i].SubGrid[j] 
		       << "  " << Grid->Mesh[i].U_points[j] << "\n";
           output_file.unsetf(ios::scientific);
       } /* endfor */
    } /* endfor */
    output_file << setprecision(6);

    /* Close the output data file. */

    output_file.close();

    /* Write the Matlab script file. */
    strcpy(output_file_name_script, output_directory);
    strcat(output_file_name_script, prefix);
    strcat(output_file_name_script, ".m");
    output_file_name_ptr = output_file_name_script;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    output_file << "%Matlab file for visualizing "<< output_file_name_data
		<< " file." <<endl << endl;

    output_file << "[x U] = textread('"<< output_file_name_data
		<<"','%f %f',...\n" 
		<< "'headerlines',3);" << endl;

    output_file << "figure(3) \n"
		<< "h = plot (x, U,'*','MarkerEdgeColor','r',"
		<< "'MarkerFaceColor','y', 'MarkerSize',5);" <<endl 
		<< "set (h, 'LineWidth', .5, {'LineStyle'}, {'-';})" <<endl
		<< "set (h, {'Color'}, {'g'})" << endl 
                //<< "axis equal tight" << endl !!! Just in case you need !!!
		<< "grid on" << endl
		<< "xlabel ('Domain', 'FontSize',12)" << endl
		<< "ylabel ('Function Distribution', 'FontSize',12,"
		<< "'Rotation',90)" << endl
		<< "legend (h, 'Function reprezentation' ) "<< endl
		<< "title ('" << IP.Reconstruction_Order;
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
    } //endswitch
    output_file <<" Order 1D Reconstruction Solution for function "
                << IP.Function_Type << "',..." << endl
		<< "'FontSize', 14)";

    /* Close the output data file. */

    output_file.close();

    /* Writing of Matlab commands file complete.  Return zero value. */

    return(0);
}

/********************************************************
 * Routine: Output_Derivatives_Tecplot                  *
 *                                                      *
 * Writes the derivative values at the subnodes of each *
 * cell for 1D reconstruction to the specified output   *
 * data file(s) in a format suitable for plotting with  *
 * TECPLOT.                                             *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
int Output_Derivatives_Tecplot(const CompDomain1D_NonUniform *Grid,
			       const Reconstruct1D_Input_Parameters &IP) {

    int i;
    char prefix[256], extension[256], output_file_name[256], output_directory[30];
    char *output_file_name_ptr;
    ofstream output_file;    
    double *Deriv_T;
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
    strcat(prefix,"_Deriv");

    if (IP.Reconstruction_Order == 0){
      cout << "\n MESSAGE "<< IP.Message_Number
	   << ": The derivatives are printed ONLY for Reconstruction_Order higher than 0.\a";
      //      IP.Message_Number++;
      return(0);
    }

    for (int order=1; order<=IP.Reconstruction_Order; order++){ // for each order

      // Determine output data file name for this reconstruction order.
      sprintf(extension, "%.1u", order);
      strcat(extension,".dat");
      strcpy(output_directory, "../Results/");
      strcpy(output_file_name, output_directory);
      strcat(output_file_name, prefix);
      strcat(output_file_name, extension);
      output_file_name_ptr = output_file_name;
      
      // Open the output data file.
      output_file.open(output_file_name_ptr, ios::out);
      if (output_file.bad()) return (1);

      // Write the solution data.
      output_file << setprecision(14);
      output_file << "TITLE = \"" << order;
      switch (order){
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
      output_file <<" Derivative Solution for 1D Reconstruction of function "
		  << IP.Function_Type << "\"" << "\n"
		  << "VARIABLES = \"x\" \"U_numeric\" \"U_theoretical\" \n"
		  << "ZONE\n";

      int I_min  = Grid->NCS;
      int I_max = Grid->NCi - Grid->NCS - 1;

      for ( i=I_min; i<=I_max; i++){ //for each cell in the Domain

	// allocate memory for theoretical derivative vector
	Deriv_T = new double[Grid->Mesh[i].N];
	
	// Compute the numeric solution for the specific derivative
	Grid->Mesh[i].SetDerivative_Representation(order);
	
	// Compute the theoretical solution for the specific derivative
       	Cell_Theoretic_Derives(Grid->Mesh[i], IP, order, Deriv_T);
	
	// print out the cell solution
	for (int j = 0 ; j <= Grid->Mesh[i].N-1 ; j++ ) {
	  output_file << " " << Grid->Mesh[i].SubGrid[j] 
		      << " " << Grid->Mesh[i].U_deriv_num[j]
		      << " " << Deriv_T[j] << "\n";
	  output_file.unsetf(ios::scientific);
	} // endfor 
      } // endfor
      output_file << setprecision(6);
      
      // Close the output data file.
      
      output_file.close();
      
      // deallocate memory for theoretical derivative vector

      delete [] Deriv_T; Deriv_T = NULL;
      
    } // endfor order 
    
    // Writing of output data files complete.  Return zero value.
    
    return(0);
}

/********************************************************
 * Routine: Output_Derivatives_Matlab                   *
 *                                                      *
 * Writes the derivative values at the nodes of each    *
 * cell for 1D reconstruction to the specified output   *
 * data file(s) in a format suitable for plotting with  *
 * MATLAB.                                              *
 * Returns a non-zero value if cannot write any of the  *
 * MATLAB solution files.                               *
 *                                                      *
 ********************************************************/

int Output_Derivatives_Matlab(const CompDomain1D_NonUniform *Grid,
			      const Reconstruct1D_Input_Parameters &IP) {

    int i;
    char prefix[256], extension[256], output_file_name_data[256],  output_directory[30];
    char output_file_name_script[256];
    char *output_file_name_ptr;
    ofstream output_file;    
    double *Deriv_T;
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
    strcat(prefix,"_Deriv");

    // Determine output data file name.

    if (IP.Reconstruction_Order == 0){
      cout << "\n MESSAGE "<< IP.Message_Number
	   << ": The derivatives are printed ONLY for Reconstruction_Order higher than 0.\a";
      //      IP.Message_Number++;
      return(0);
    }
    for (int order=1; order<=IP.Reconstruction_Order; order++){

      // Determine output data file names for this reconstruction order.
      sprintf(extension, "%.1u", order);
      strcpy(output_directory, "../Results/");
      strcpy(output_file_name_data, output_directory);
      strcat(output_file_name_data, prefix);
      strcat(output_file_name_data, extension);
      strcat(output_file_name_data, ".dat");
      output_file_name_ptr = output_file_name_data;

      // Open the output data file.

      output_file.open(output_file_name_ptr, ios::out);
      if (output_file.bad()) return (1);
      
      // Write the solution data.
      
      output_file << setprecision(14);
      output_file << "%TITLE = \"" << order;
      switch (order){
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
      } //endswitch
      output_file <<" Derivative Solution for 1D Reconstruction of function "
		  << IP.Function_Type << "\"" << "\n"
		  << "%VARIABLES \n" 
		  << " x \t U_numerical \t U_theoretical \n";

      int I_min  = Grid->NCS;
      int I_max = Grid->NCi - Grid->NCS - 1;
      
      for ( i=I_min; i<=I_max; i++){ //for each cell in the Domain
      
	// allocate memory for theoretical derivative vector
	Deriv_T = new double[1];
	
	// Compute the numeric solution for the specific derivative
	Grid->Mesh[i].SetDerivative_Representation(order);
	
	// Compute the theoretical solution for the specific derivative
       	Cell_Theoretic_Derives(Grid->Mesh[i], IP, order, Deriv_T);
	
	// print out the cell solution
	for (int j = 0 ; j <= Grid->Mesh[i].N-1 ; j++ ) {
	    output_file << " " << Grid->Mesh[i].SubGrid[j] 
			<< "  " << Grid->Mesh[i].U_deriv_num[j]
			<< " " << Deriv_T[j] << "\n";
	    output_file.unsetf(ios::scientific);
	} // endfor
      } // endfor
      output_file << setprecision(6);
      
      // Close the output data file.
      
      output_file.close();
      
      // Write the Matlab script file.
      strcpy(output_file_name_script, output_directory);
      strcat(output_file_name_script, prefix);
      strcat(output_file_name_script, extension);
      strcat(output_file_name_script, ".m");
      output_file_name_ptr = output_file_name_script;
      
      // Open the output data file.
      
      output_file.open(output_file_name_ptr, ios::out);
      if (output_file.bad()) return (1);
      
      output_file << "%Matlab file for visualizing "<< output_file_name_data
		  << " file." <<endl << endl;
      
      output_file << "[x U_numerical U_theoretical] = textread('"<< output_file_name_data
		  <<"','%f %f %f',...\n" 
		  << "'headerlines',3);" << endl;
      
      output_file << "figure(3) \n"
		  << "h = plot (x, U_numerical,'*','MarkerEdgeColor','b',"
		  << "'MarkerFaceColor','y', 'MarkerSize',5);" <<endl 
		  << "set (h, 'LineWidth', .5, {'LineStyle'}, {'-';})" <<endl
		  << "set (h, {'Color'}, {'g'})" << endl 
		  << "hold on" << endl
		  << "h1 = plot (x, U_theoretical,'+','MarkerEdgeColor','r',"
		  << "'MarkerFaceColor','y', 'MarkerSize',5);" <<endl 
		  << "set (h, 'LineWidth', .5, {'LineStyle'}, {'--';})" <<endl
		  << "set (h, {'Color'}, {'g'})" << endl 
		  << "set (h1, 'LineWidth', .5, {'LineStyle'}, {'-';})" <<endl
		  << "set (h1, {'Color'}, {'r'})" << endl 
	//<< "axis equal tight" << endl !!! Just in case you need !!!
		  << "grid on" << endl
		  << "xlabel ('Domain', 'FontSize',12)" << endl
		  << "ylabel ('Derivative Variation', 'FontSize',12,"
		  << "'Rotation',90)" << endl
		  << "legend ('Numerical Derivative', 'Analytical Derivative')"
		  << endl
		  << "title ('" << order;
      switch (order){
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
      } //endswitch
      output_file << " Derivative Solution for 1D Reconstruction of function "
		  << IP.Function_Type << "',..." << endl
		  << "'FontSize', 14)";
      
      // Close the output data file.
      
      output_file.close();
      
    }
    // Writing of Matlab commands file complete.  Return zero value.
    
    return(0);
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

int Output_Cell_Center_Values(const CompDomain1D_NonUniform *Grid,
			      const Reconstruct1D_Input_Parameters &IP){


  int i;
  char prefix[256], output_file_name[256];
  char output_directory[30];
  char *output_file_name_ptr;
  ofstream output_file;    
  double F_val; // the value of the function at the cell center  

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
  strcat(prefix,"_CellCenterValues.dat");
  strcpy(output_directory, "../Output_Data/");
  strcpy(output_file_name, output_directory);
  strcat(output_file_name, prefix);
  output_file_name_ptr = output_file_name;
    
  // Open the output data file.
    
  output_file.open(output_file_name_ptr, ios::out);
  if (output_file.bad()) return (1);

  // Write the solution data.
  output_file << "TITLE = \"The cell center values for function"
	      << IP.Function_Type << "\"" << endl;
  output_file << "VARIABLES = "
	      << setw(8) << "\"x\""
	      << setw(20) << "\"U numeric\""
	      << setw(20) << "\"U theoretical\"\n"
	      << "ZONE \n";
  output_file << setprecision(14);

  int I_min  = Grid->NCS;
  int I_max = Grid->NCi - Grid->NCS - 1;

  for (i=I_min; i<= I_max; i++) {
    F_val = IP.TestF(Grid->Mesh[i].geom.x ); // pointer to the asked functin during runtime
    output_file << " " << setw(20) << Grid->Mesh[i].geom.x 
		<< "  " << setw(20) << Grid->Mesh[i].U_cell
		<< "  " << setw(20) << F_val << "\n";
    output_file.unsetf(ios::scientific);
  } // endfor
  output_file << setprecision(6);

  // Close the output data file.
      
  output_file.close();
  
  // Writing of solution file complete.  Return zero value.

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
int Output_Function_Graph_Tecplot(const Reconstruct1D_Input_Parameters &IP) {

    int i;
    char prefix[256], extension[256], output_file_name[256], output_directory[30];
    char *output_file_name_ptr;
    ofstream output_file;    

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
    strcpy(output_directory, "../Results/");
    strcpy(extension,"_Function_");
    strcat(extension,IP.Function_Type);
    strcpy(output_file_name, output_directory);
    strcat(output_file_name, prefix);
    strcat(output_file_name, extension);
    strcat(output_file_name, ".dat");
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data. */

    output_file << setprecision(14);
    output_file << "TITLE = \"" << IP.Function_Type << " representation on the domain: ";
    output_file << IP.X_min << " :  " << IP.X_min + IP.LengthX() <<"\""<< endl
	        << "VARIABLES = \"x\" \"U\" \n"
                << "ZONE \n";

    int N_standard = 4000; // the number of points used for representing a domain length
                          // equal to 1
    int N_points = (int) IP.LengthX()*N_standard;

    double* FunctionGrid;
    double* FunctionValues;
    FunctionGrid = new double [N_points];
    FunctionValues = new double [N_points];

    double delta;

    delta = IP.LengthX()/(N_points-1);
    FunctionGrid[0] = IP.X_min;
    for (i=1; i<=N_points-1; i++)
      FunctionGrid[i] = FunctionGrid[0] + i*delta;

    for (int i=0; i<=N_points-1; i++){
      FunctionValues[i]= IP.TestF(FunctionGrid[i]); // pointer to the function run at runtime
    }

    for ( i = 0; i<= N_points-1; i++ ) {
           output_file << " " << FunctionGrid[i] 
		       << "  " << FunctionValues[i] << "\n";
           output_file.unsetf(ios::scientific);
       } /* endfor */
    output_file << setprecision(6);

    delete [] FunctionGrid; FunctionGrid = NULL;
    delete [] FunctionValues; FunctionValues = NULL;

    /* close the output data file. */

    output_file.close();       

    /* Writing of output data files complete.  Return zero value. */

    return(0);

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

int Output_Function_Graph_Matlab(const Reconstruct1D_Input_Parameters &IP) {
 
    int i;
    char prefix[256], extension[256], output_file_name_data[256], output_directory[30];
    char output_file_name_script[256];
    char *output_file_name_ptr;
    ofstream output_file;    

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

    strcpy(output_directory, "../Results/");
    strcpy(extension,"_Function_");
    strcat(extension,IP.Function_Type);
    strcpy(output_file_name_data, output_directory);
    strcat(output_file_name_data, prefix);
    strcat(output_file_name_data, extension);
    strcat(output_file_name_data, "_mtb.dat");
    output_file_name_ptr = output_file_name_data;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data. */

    output_file << setprecision(14);
    output_file << "%TITLE = \"" << IP.Function_Type << " representation on the domain: ";
    output_file << IP.X_min << " :  " << IP.X_min + IP.LengthX() <<"\""<< endl
 	        << "%VARIABLES \n" 
 		<< " x \t U \n";

    int N_standard = 100; // the number of points used for representing a domain length equal to 1
    int N_points = (int) IP.LengthX()*N_standard;

    double* FunctionGrid;
    double* FunctionValues;
    FunctionGrid = new double [N_points];
    FunctionValues = new double [N_points];

    double delta;

    delta = IP.LengthX()/(N_points-1);
    FunctionGrid[0] = IP.X_min;
    for (i=1; i<=N_points-1; i++)
      FunctionGrid[i] = FunctionGrid[0] + i*delta;

    for (int i=0; i<=N_points-1; i++){
      FunctionValues[i]= IP.TestF(FunctionGrid[i]);
    }

    for ( i = 0; i<= N_points-1; i++ ) {
           output_file << " " << FunctionGrid[i] 
		       << "  " << FunctionValues[i] << "\n";
           output_file.unsetf(ios::scientific);
       } /* endfor */
    output_file << setprecision(6);

    delete [] FunctionGrid; FunctionGrid = NULL;
    delete [] FunctionValues; FunctionValues = NULL;

    /* Close the output data file. */

    output_file.close();

    /* Write the Matlab script file. */
    strcpy(output_file_name_script, output_directory);
    strcat(output_file_name_script, prefix);
    strcat(output_file_name_script, extension);
    strcat(output_file_name_script, ".m");
    output_file_name_ptr = output_file_name_script;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    output_file << "%Matlab file for visualizing "<< output_file_name_data
 		<< " file." <<endl << endl;

    output_file << "[x U] = textread('"<< output_file_name_data
 		<<"','%f %f',...\n" 
 		<< "'headerlines',3);" << endl;

    output_file << "figure(3) \n"
		<< "h = plot (x, U,'*','MarkerEdgeColor','r',"
 		<< "'MarkerFaceColor','y', 'MarkerSize',5);" <<endl 
 		<< "set (h, 'LineWidth', .5, {'LineStyle'}, {'-';})" <<endl
 		<< "set (h, {'Color'}, {'g'})" << endl 
                //<< "axis equal tight" << endl !!! Just in case you need !!!
 		<< "grid on" << endl
 		<< "xlabel ('Domain', 'FontSize',12)" << endl
 		<< "ylabel ('Function Distribution', 'FontSize',12,"
 		<< "'Rotation',90)" << endl
 		<< "legend (h, 'Function reprezentation' ) "<< endl
 		<< "title ('" << IP.Function_Type
		<< " representation on the domain: ";
    output_file << IP.X_min << " :  " << IP.X_min + IP.LengthX() 
		<<"',..."<< endl
 		<< "'FontSize', 14)";

       /* Close the output data file. */

       output_file.close();

       /* Writing of Matlab commands file complete.  Return zero value. */

    return(0);
}

/********************************************************
 * Routine: Output_Derivatives_Tecplot                  *
 *                                                      *
 * Writes the derivative values at the subnodes of each *
 * cell for 1D reconstruction to the specified output   *
 * data file(s) in a format suitable for plotting with  *
 * TECPLOT.                                             *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/

void Cell_Theoretic_Derives(const CompCell1D_NonUniform &Cell,
			    const Reconstruct1D_Input_Parameters &IP,
			    int Order, // The order of the derivative
			    double* D  // The derivatives' vector
			    ){
 
  for (int sg=0; sg<=Cell.N-1; sg++){
    switch(Order){
    case 1:
      switch(IP.i_Function){
      case FUNCTION_DEFAULT:
	D[sg] = Test_Default1D_FirstDeriv(Cell.SubGrid[sg]);
	break;
      case FUNCTION_EXAMPLE_1:
	D[sg] = Test_Example1_FirstDeriv(Cell.SubGrid[sg]);
	break;
      case FUNCTION_EXAMPLE_2:
	D[sg] = Test_Example2_FirstDeriv(Cell.SubGrid[sg]);
	break;
      default: // for an undefined derivative function
	D[sg] = 0.0;
      }
      break;
    case 2:
      switch(IP.i_Function){
      case FUNCTION_DEFAULT:
	D[sg] = Test_Default1D_SecondDeriv(Cell.SubGrid[sg]);
	break;
      case FUNCTION_EXAMPLE_1:
	D[sg] = Test_Example1_SecondDeriv(Cell.SubGrid[sg]);
	break;
      case FUNCTION_EXAMPLE_2:
	D[sg] = Test_Example2_SecondDeriv(Cell.SubGrid[sg]);
	break;
      default: // for an undefined derivative function
	D[sg] = 0.0;
      }
      break;
    case 3:
      switch(IP.i_Function){
      case FUNCTION_DEFAULT:
	D[sg] = Test_Default1D_ThirdDeriv(Cell.SubGrid[sg]);
	break;
      case FUNCTION_EXAMPLE_1:
	D[sg] = Test_Example1_ThirdDeriv(Cell.SubGrid[sg]);
	break;
      case FUNCTION_EXAMPLE_2:
	D[sg] = Test_Example2_ThirdDeriv(Cell.SubGrid[sg]);
	break;
      default: // for an undefined derivative function
	D[sg] = 0.0;
      }
      break;
    case 4:
      switch(IP.i_Function){
      case FUNCTION_DEFAULT:
	D[sg] = Test_Default1D_FourthDeriv(Cell.SubGrid[sg]);
	break;
      case FUNCTION_EXAMPLE_1:
	D[sg] = Test_Example1_FourthDeriv(Cell.SubGrid[sg]);
	break;
      case FUNCTION_EXAMPLE_2:
	D[sg] = Test_Example2_FourthDeriv(Cell.SubGrid[sg]);
	break;
      default: // for an undefined derivative function
	D[sg] = 0.0;
      }
      break;
    default:  // for an undefined derivative order
      D[sg] = 0.0;
    }//endswitch
  }//endfor
}

double L1_Norm(const CompDomain1D_NonUniform *Grid,
	       const Reconstruct1D_Input_Parameters &IP){
/* Compute the L1 norm for the whole domain*/

  int ICl = Grid->NCS;
  int ICu = Grid->NCi-Grid->NCS-1;
  double L1_norm=0.0;

  for (int i=ICl; i<=ICu; i++){
    L1_norm += Grid->Mesh[i].ErrorRec;
  }
  return L1_norm/(IP.Number_of_Cells_Idir);
}


double L2_Norm(const CompDomain1D_NonUniform *Grid,
	       const Reconstruct1D_Input_Parameters &IP){
/* Compute the L2 norm for the whole domain*/

  int ICl = Grid->NCS;
  int ICu = Grid->NCi-Grid->NCS-1;
  double L2_norm=0.0;

  for (int i=ICl; i<=ICu; i++){
    L2_norm += pow(Grid->Mesh[i].ErrorRec,2);
  }
  L2_norm = sqrt(L2_norm/IP.Number_of_Cells_Idir);

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

int Output_L1_Norm_Tecplot(const CompDomain1D_NonUniform *Grid,
			   const Reconstruct1D_Input_Parameters &IP){

    static int WrittenHeader = OFF;
    static char OutputFile[256] = "START";
    int i;
    char prefix[256], extension[256], output_file_name[256], output_directory[30];
    char order[30];
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
      strcpy(order,"1"); break;
    case 2:
      strcpy(order,"2"); break;
    case 3:
      strcpy(order,"3"); break;
    case 4:
      strcpy(order,"4"); break;
    case 5:
      strcpy(order,"5"); break;
    case 6:
      strcpy(order,"6"); break;
    case 7:
      strcpy(order,"7"); break;
    default:
      strcpy(order,"0");
    }

    strcpy(output_directory, "../Results/");
    strcpy(extension,"_Function_");
    strcat(extension,IP.Function_Type);
    strcat(extension,"_L1Norm_Order");
    strcat(extension,order);
    strcat(extension,"_");
    strcpy(output_file_name, output_directory);
    strcat(output_file_name, prefix);
    strcat(output_file_name, extension);
    strcat(output_file_name, ".dat");
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    if (strcmp(OutputFile,output_file_name)){
      WrittenHeader = OFF;
      strcpy(OutputFile,output_file_name);
    }

    if (!WrittenHeader){
      output_file.open(output_file_name_ptr, ios::trunc);
      if (output_file.bad()) return (1);
      // write the header
      output_file << "TITLE = \"" << IP.Function_Type
		  << " L1Norm representation on the domain: ";
      output_file << IP.X_min << " :  " << IP.X_min + IP.LengthX()
		  <<"\""<< endl
		  << "VARIABLES = \"#Points\" \"L1Norm\" \n"
		  << "ZONE \n";
      WrittenHeader = ON;
    }
    else{
      output_file.open(output_file_name_ptr, ios::app);
      if (output_file.bad()) return (1);
    }

    cout << "Output L1 Norm in file: \n" << "\t" << output_file_name_ptr << endl;

    /* Write the solution data. */
    L1_norm = L1_Norm(Grid,IP);

    output_file.setf(ios::showpoint,ios::right);
    output_file << " " << setw(8) <<  IP.Number_of_Cells_Idir << setprecision(14)
		<< "   " << L1_norm << "\n";
    output_file.unsetf(ios::scientific);
    output_file << setprecision(6);

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Output_L2_Norm_Tecplot                      *
 *                                                      *
 * Writes the L2_Norm of the reconstruction             *
 * to the specified output data                         *
 * file(s) in a format suitable for plotting            *
 * with TECPLOT.                                        *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 * The first call of the function writes the header file*
 * and the next calls only append to the end of the file*
 *                                                      *
 ********************************************************/

int Output_L2_Norm_Tecplot(const CompDomain1D_NonUniform *Grid,
			   const Reconstruct1D_Input_Parameters &IP){

    static int WrittenHeaderL2 = OFF;
    static char OutputFileL2[256] = "START"; 
    int i;
    char prefix[256], extension[256], output_file_name[256], output_directory[30];
    char order[30];
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

    strcpy(output_directory, "../Results/");
    strcpy(extension,"_Function_");
    strcat(extension,IP.Function_Type);
    strcat(extension,"_L2Norm_Order");
    strcat(extension,order);
    strcat(extension,"_");
    strcpy(output_file_name, output_directory);
    strcat(output_file_name, prefix);
    strcat(output_file_name, extension);
    strcat(output_file_name, ".dat");
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    if (strcmp(OutputFileL2,output_file_name)){
      WrittenHeaderL2 = OFF;
      strcpy(OutputFileL2,output_file_name);
    }


    if (!WrittenHeaderL2){
      output_file.open(output_file_name_ptr, ios::trunc);
      if (output_file.bad()) return (1);
      // write the header
      output_file << "TITLE = \"" << IP.Function_Type
		  << " L2Norm representation on the domain: ";
      output_file << IP.X_min << " :  " << IP.X_min + IP.LengthX()
		  <<"\""<< endl
		  << "VARIABLES = \"#Points\" \"L2Norm\" \n"
		  << "ZONE \n";
      WrittenHeaderL2 = ON;
    }
    else{
      output_file.open(output_file_name_ptr, ios::app);
      if (output_file.bad()) return (1);
    }
    
    cout << "Output L2 Norm in file: \n" << "\t" << output_file_name_ptr << endl;

    /* Compute the solution data */
    L2_norm = L2_Norm(Grid,IP);

    /* Write the solution data. */
    output_file.setf(ios::showpoint,ios::right);
    output_file << " " << setw(8) <<  IP.Number_of_Cells_Idir << setprecision(14)
		<< "   " << L2_norm << "\n";
    output_file.unsetf(ios::scientific);
    output_file << setprecision(6);

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);
}
