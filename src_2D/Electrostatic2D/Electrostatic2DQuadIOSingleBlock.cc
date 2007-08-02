/**********************************************************************
 * Electrostatic2DQuadIOSingleBlock.cc                                *
 *                                                                    *
 * Single-block versions of input and output subroutines for 2D       *
 * Electrostatic multi-block quadrilateral mesh solution classes.     *
 *                                                                    *
 **********************************************************************/

// Include 2D Electrostatic quadrilateral mesh solution header file.

#ifndef _ELECTROSTATIC2D_QUAD_INCLUDED
#include "Electrostatic2DQuad.h"
#endif // _ELECTROSTATIC2D_QUAD_INCLUDED

/**********************************************************************
 * Electrostatic2D_Quad_Block -- IO Single Block External             *
 *                               Subroutines.                         *
 **********************************************************************/

/**********************************************************************
 * Routine: Write_Solution_Block                                      *
 *                                                                    *
 * Writes the cell centred solution values of the specified           *
 * quadrilateral solution block to the specified output stream for    *
 * restart purposes.                                                  *
 *                                                                    *
 **********************************************************************/
void Write_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk,
	                  ostream &Out_File) {

  Out_File << setprecision(14) << SolnBlk << setprecision(6);

}

/**********************************************************************
 * Routine: Read_Solution_Block                                       *
 *                                                                    *
 * Reads the cell centred solution values for the specified           *
 * quadrilateral solution block from the specified input stream as    *
 * required for restart purposes.                                     *
 *                                                                    *
 **********************************************************************/
void Read_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk,
	                 istream &In_File) {

  In_File >> SolnBlk;

}

/**********************************************************************
 * Routine: Output_Tecplot                                            *
 *                                                                    *
 * Writes the solution values at the nodes of the specified           *
 * quadrilateral solution block to the specified output stream        *
 * suitable for plotting with TECPLOT.                                *
 *                                                                    *
 **********************************************************************/
void Output_Tecplot(Electrostatic2D_Quad_Block &SolnBlk,
		    Electrostatic2D_Input_Parameters &IP,
		    const int Number_of_Time_Steps,
		    const double &Time,
		    const int Block_Number,
		    const int Output_Title,
		    ostream &Out_File) {

  Electrostatic2DState U_node;

  // Ensure boundary conditions are updated before evaluating
  // solution at the nodes.
  BCs(SolnBlk,IP);

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFDkit_Name() << ": 2D Electrostatic Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"Ex\" \\ \n"
	     << "\"Ey\" \\ \n"
	     << "\"V\" \\ \n";
  }

  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	   << "F = POINT \\ \n";

  for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
    for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
      U_node = SolnBlk.Un(i,j);
      Out_File.setf(ios::scientific);
      Out_File << " " << SolnBlk.Grid.Node[i][j].X
	       << " " << U_node.E.x
	       << " " << U_node.E.y
	       << " " << U_node.V
               << endl;
    }
  }

  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Output_Cells_Tecplot                                      *
 *                                                                    *
 * Writes the cell centred solution values of the specified           *
 * quadrilateral solution block to the specified output stream        *
 * suitable for plotting with TECPLOT.                                *
 *                                                                    *
 **********************************************************************/
void Output_Cells_Tecplot(Electrostatic2D_Quad_Block &SolnBlk,
			  Electrostatic2D_Input_Parameters &IP,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int CPU_Number,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFDkit_Name() << ": 2D Electrostatic Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"Ex\" \\ \n"
	     << "\"Ey\" \\ \n"
	     << "\"V\" \\ \n";
  }
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 2*SolnBlk.Nghost + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 2*SolnBlk.Nghost + 1 << " \\ \n"
	   << "F = POINT \n";

  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      Out_File.setf(ios::scientific);
      Out_File << " " << SolnBlk.Grid.Cell[i][j].Xc
	       << " " << SolnBlk.U[i][j].E.x
	       << " " << SolnBlk.U[i][j].E.y
	       << " " << SolnBlk.U[i][j].V
	       << endl;
    }
  }
  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Output_Quasi3D_Tecplot                                    *
 *                                                                    *
 * Writes the solution values at the nodes of the specified active    *
 * quadrilateral solution block to the specified output stream        *
 * suitable for plotting with TECPLOT in a quasi-3D format.           *
 *                                                                    *
 **********************************************************************/
void Output_Quasi3D_Tecplot(Electrostatic2D_Quad_Block &SolnBlk,
			    Electrostatic2D_Input_Parameters &IP,
			    const int Number_of_Time_Steps,
			    const double &Time,
			    const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File) {

  Electrostatic2DState U_node;
  int nrr, numberofrotations = 360/15;
  int numberofnodes = ((SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1)*
		       (SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1));
  int numberofcells = ((SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1)*
		       (SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1));

  // Ensure boundary conditions are updated before evaluating
  // solution at the nodes.
  BCs(SolnBlk,IP);

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFDkit_Name() << ": 2D Electrostatic Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"z\" \\ \n"
	     << "\"Ex\" \\ \n"
	     << "\"Ey\" \\ \n"
	     << "\"Ez\" \\ \n"
	     << "\"V\" \\ \n";
  }

  for (int nr = 0; nr < numberofrotations; nr++) {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number << nr
	     << "\" \\ \n"
	     << "N = " << 2*numberofnodes << " \\ \n"
	     << "E = " << numberofcells << " \\ \n"
	     << "F = FEPOINT \\ \n"
	     << "ET = BRICK \n";
    for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
      for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
	U_node = SolnBlk.Un(i,j);
	Out_File.setf(ios::scientific);
	Out_File << " " << SolnBlk.Grid.Node[i][j].X.x
		 << " " << SolnBlk.Grid.Node[i][j].X.y*sin(TWO*PI*double(nr)/double(numberofrotations))
		 << " " << SolnBlk.Grid.Node[i][j].X.y*cos(TWO*PI*double(nr)/double(numberofrotations))
		 << " " << U_node.E.x
		 << " " << U_node.E.y*sin(TWO*PI*double(nr)/double(numberofrotations))
		 << " " << U_node.E.y*cos(TWO*PI*double(nr)/double(numberofrotations))
		 << " " << U_node.V
		 << endl;
      }
    }
    if (nr < numberofrotations-1) nrr = nr + 1;
    else nrr = 0;
    for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
      for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
	U_node = SolnBlk.Un(i,j);
	Out_File.setf(ios::scientific);
	Out_File << " " << SolnBlk.Grid.Node[i][j].X.x
		 << " " << SolnBlk.Grid.Node[i][j].X.y*sin(TWO*PI*double(nrr)/double(numberofrotations))
		 << " " << SolnBlk.Grid.Node[i][j].X.y*cos(TWO*PI*double(nrr)/double(numberofrotations))
		 << " " << U_node.E.x
		 << " " << U_node.E.y*sin(TWO*PI*double(nrr)/double(numberofrotations))
		 << " " << U_node.E.y*cos(TWO*PI*double(nrr)/double(numberofrotations))
		 << " " << U_node.V
		 << endl;
      }
    }

    // Connectivity table.
    for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
      for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
	Out_File << (j-2)*(SolnBlk.Grid.INu-1) + i - 1                         << " "
		 << (j-2)*(SolnBlk.Grid.INu-1) + i                             << " "
		 << (j-1)*(SolnBlk.Grid.INu-2) + i + 1 + (j-2)                 << " "
		 << (j-1)*(SolnBlk.Grid.INu-2) + i     + (j-2)                 << " "
		 << (j-2)*(SolnBlk.Grid.INu-1) + i - 1         + numberofnodes << " "
		 << (j-2)*(SolnBlk.Grid.INu-1) + i             + numberofnodes << " "
		 << (j-1)*(SolnBlk.Grid.INu-2) + i + 1 + (j-2) + numberofnodes << " "
		 << (j-1)*(SolnBlk.Grid.INu-2) + i     + (j-2) + numberofnodes << endl;
      }
    }

  }

  Out_File << setprecision(6);

}
