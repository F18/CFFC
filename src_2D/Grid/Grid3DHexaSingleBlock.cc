/* Grid3DHexaSingleBlock.cc:  Single-block subroutines for 
                              3D hexahedral block grid class. */

/* Include 3D hexahedral block grid type header file. */

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "Grid3DHexa.h"
#endif // _GRID3D_HEXA_BLOCK_INCLUDED

/*************************************************************************
 * Grid3D_Hexa_Block -- External subroutines for single grid block.      *
 *************************************************************************/

/********************************************************
 * Routine: Write_Hexa_Block                            *
 *                                                      *
 * Writes 3D hexahedral    grid block to the            *
 * specified output stream for retrieval and re-use     *
 * purposes.                                            *
 *                                                      *
 ********************************************************/
void Write_Hexa_Block(Grid3D_Hexa_Block &Grid,
	              ostream &Out_File) {

    Out_File << setprecision(14) << Grid << setprecision(6);

}

/********************************************************
 * Routine: Read_Hexa_Block                             *
 *                                                      *
 * Reads 3D hexahedral    grid block from the           *
 * specified input stream.                              *
 *                                                      *
 ********************************************************/
void Read_Hexa_Block(Grid3D_Hexa_Block &Grid,
	             istream &In_File) {

    In_File >> Grid;

}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the nodes of the hexahedral    mesh to the    *
 * specified output stream in a format suitable for     *
 * plotting the grid with TECPLOT.                      *
 *                                                      *
 ********************************************************/
void Output_Tecplot(Grid3D_Hexa_Block &Grid,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFDkit_Name()
                << ": 3D Structured Curvilinear Grid Block (Node Locations)"
                << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
	        << "\"z\" \n"  
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
                << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
		<< "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
                << "F = POINT \n";

    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
                << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
		<< "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    for ( int k = Grid.KNl ; k <= Grid.KNu ; ++k){
      for ( int j = Grid.JNl ; j <= Grid.JNu ; ++j ) {
	for ( int i = Grid.INl ; i <= Grid.INu ; ++i ) {
	  Out_File << " " << Grid.Node[i][j][k].X << "\n";
	} /* endfor */
      } /* endfor */
    }/*endfor*/
    Out_File << setprecision(6);

}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the nodes of the quadrilateral mesh to the    *
 * specified output stream in a format suitable for     *
 * plotting the grid with TECPLOT.  Includes boundary   *
 * nodes.                                               *
 *                                                      *
 ********************************************************/
void Output_Nodes_Tecplot(Grid3D_Hexa_Block &Grid,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFDkit_Name()
                << ": 3D Structured Curvilinear Grid Block (Node Locations)"
                << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
		<< "\"z\" \n"
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.INu - Grid.INl + 1 + 4 << " \\ \n"
                << "J = " << Grid.JNu - Grid.JNl + 1 + 4 << " \\ \n"
		<< "K = " << Grid.KNu - Grid.KNl + 1 + 4 << " \\ \n"
                << "F = POINT \n";

    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.INu - Grid.INl + 1 + 4 << " \\ \n"
                << "J = " << Grid.JNu - Grid.JNl + 1 + 4 << " \\ \n"
		<< "K = " << Grid.KNu - Grid.KNl + 1 + 4 << " \\ \n" 
                << "F = POINT \n";
    } /* endif */
    
    for(int k = Grid.KNl-2 ; k <= Grid.KNu+2 ; ++k) {
      for (int j  = Grid.JNl-2 ; j <= Grid.JNu+2 ; ++j ) {
	for (int i = Grid.INl-2 ; i <= Grid.INu+2 ; ++i ) {
	  Out_File << " " << Grid.Node[i][j][k].X << "\n";
	} /* endfor */
      } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cells of the quadrilateral mesh to the    *
 * specified output stream in a format suitable for     *
 * plotting the grid with TECPLOT.                      *
 *                                                      *
 ********************************************************/
void Output_Cells_Tecplot(Grid3D_Hexa_Block &Grid,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {


    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFDkit_Name()
                << ": 3D Structured Curvilinear Grid Block (Cell Locations)"
                << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
		<< "\"z\" \n"
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.ICu - Grid.ICl + 1 << " \\ \n"
                << "J = " << Grid.JCu - Grid.JCl + 1 << " \\ \n"
		<< "K = " << Grid.KCu - Grid.KCl + 1 << " \\ \n"
                << "F = POINT \n";

    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.ICu - Grid.ICl + 1 << " \\ \n"
                << "J = " << Grid.JCu - Grid.JCl + 1 << " \\ \n"
		<< "K = " << Grid.KCu - Grid.KCl + 1 << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    for (int k = Grid.KCl ; k <= Grid.KCu ; ++k) {
      for (int j  = Grid.JCl ; j <= Grid.JCu ; ++j ) {
	for (int i = Grid.ICl ; i <= Grid.ICu ; ++i ) {
	  Out_File << " " << Grid.Cell[i][j][k].Xc << "\n";
	} /* endfor */
      } /* endfor */
    }/*endfor*/
    Out_File << setprecision(6);

}

/********************************************************
 * Routine: Output_Gnuplot                              *
 *                                                      *
 * Writes the nodes of the hexahedral    mesh to the    *
 * specified output stream in a format suitable for     *
 * plotting the grid with GNUPLOT.                      *
 *                                                      *
 ********************************************************/
void Output_Gnuplot(Grid3D_Hexa_Block &Grid,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {


    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "# " << CFDkit_Name()
                << ": 3D Structured Curvilinear Grid Block (Node Locations)"
                << "\n"
	        << "# x(m), y(m), z(m)\n";
    } /* endif */
    
    for( int k = Grid.KNl ; k <= Grid.KNu ; ++k) {
      for (int j  = Grid.JNl ; j <= Grid.JNu ; ++j ) {
	for (int i = Grid.INl ; i <= Grid.INu ; ++i ) {
	  Out_File << " " << Grid.Node[i][j][k].X << "\n";
	} /* endfor */
	Out_File << "\n";
      } /* endfor */
      Out_File << "\n";
    }/*endfor*/
    for (int i  = Grid.INl ; i <= Grid.INu ; ++i ) {
      for (int j = Grid.JNl ; j <= Grid.JNu ; ++j ) {
	for (int k = Grid.KNl ; k <= Grid.KNu ; ++k ){
	  Out_File << " " << Grid.Node[i][j][k].X << "\n";
	} /* endfor */
	Out_File << "\n";
      } /* endfor */
      Out_File << "\n";
    }/*endfor*/
    Out_File << setprecision(6);
    
}

/*****************************************************************
 * Routine: Update Cells
 * 
 * Updates the cell information for the quadrilateal mesh block
 *                                                           
 *****************************************************************/

void Update_Cells(Grid3D_Hexa_Block &Grid) {
  for(int k = Grid.KCl-Grid.Nghost; k <=Grid.KCu +Grid.Nghost ; ++k){ 
    for( int j = Grid.JCl-Grid.Nghost ; j <= Grid.JCu +Grid.Nghost; ++j) {
      for ( int i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
	Grid.Cell[i][j][k].I = i ;
	Grid.Cell[i][j][k].J = j ;
	Grid.Cell[i][j][k].K = k ;
	Grid.Cell[i][j][k].Xc = Grid.centroid(i, j, k);
       	Grid.Cell[i][j][k].V = Grid.volume(Grid.Cell[i][j][k]);
      } /* endfor */
    } /* endfor */
  } /* endfor */
}

/**************************************************************
Update exterior nodes
**************************************************************/
void Update_Exterior_Nodes(Grid3D_Hexa_Block &Grid){
  Vector3D norm_dir, X_norm, X_tan;
  Vector3D norm_dir1,X_norm1,X_tan1;
  for(int k = Grid.KNl; k<=Grid.KNu; ++k){
    for(int j = Grid.JNl; j<=Grid.JNu; ++j){
      for(int z=1; z<=Grid.Nghost; z++){
	if(j>Grid.JNl && j < Grid.JNu){
	  //create ghost cells on west side
	  norm_dir = - HALF*(Grid.nfaceW(Grid.ICl, j, k) +
			     Grid.nfaceW(Grid.ICl, j-1, k));
	  X_norm = ((Grid.Node[Grid.INl+z][j][k].X-
		     Grid.Node[Grid.INl][j][k].X)*norm_dir)*norm_dir;
	  X_tan  = (Grid.Node[Grid.INl+z][j][k].X - 
		    Grid.Node[Grid.INl][j][k].X)-X_norm;
	  
	  Grid.Node[Grid.INl-z][j][k].X = Grid.Node[Grid.INl][j][k].X - X_norm + X_tan; 
	  //create ghost cells on east side
	  norm_dir1 = - HALF*(Grid.nfaceE(Grid.ICl, j, k) +
			     Grid.nfaceE(Grid.ICl, j-1, k));
	  X_norm1 = ((Grid.Node[Grid.INu][j][k].X-
		     Grid.Node[Grid.INu-z][j][k].X)*norm_dir1)*norm_dir1;
	  X_tan1  = (Grid.Node[Grid.INu][j][k].X - 
		    Grid.Node[Grid.INu-z][j][k].X)-X_norm1;
	  
	  Grid.Node[Grid.INu+z][j][k].X = Grid.Node[Grid.INu][j][k].X + X_norm1 - X_tan1; 
	}
      }
    }
  }
  
  for(int j = Grid.JNl; j<=Grid.JNu; ++j){
    for(int i = Grid.INl; i<=Grid.INu; ++i){
      for(int z=1; z<=Grid.Nghost; z++){
	if(i>Grid.INl && i < Grid.INu){
	  //create ghost cells on bottom side
	  norm_dir = - HALF*(Grid.nfaceBot(i, j, Grid.KCl) +
			     Grid.nfaceBot(i-1, j, Grid.KCl));
	  X_norm = ((Grid.Node[i][j][Grid.KNl+z].X-
		     Grid.Node[i][j][Grid.KNl].X)*norm_dir)*norm_dir;
	  X_tan  = (Grid.Node[i][j][Grid.KNl+z].X - 
		    Grid.Node[i][j][Grid.KNl].X)-X_norm;
	  
	  Grid.Node[i][j][Grid.KNl-z].X = Grid.Node[i][j][Grid.KNl].X - X_norm + X_tan; 
	  //create ghost cells on top side
	  norm_dir1 = - HALF*(Grid.nfaceTop(i, j, Grid.KCl) +
			     Grid.nfaceTop(i-1, j, Grid.KCl));
	  X_norm1 = ((Grid.Node[i][j][Grid.KNu].X-
		     Grid.Node[i][j][Grid.KNu-z].X)*norm_dir1)*norm_dir1;
	  X_tan1  = (Grid.Node[i][j][Grid.KNu].X - 
		    Grid.Node[i][j][Grid.KNu-z].X)-X_norm1;
	  
	  Grid.Node[i][j][Grid.KNu+z].X = Grid.Node[i][j][Grid.KNu].X + X_norm1 - X_tan1; 
	}	
      }
    }
  }
  for(int i = Grid.INl; i<=Grid.INu; ++i){
    for(int k = Grid.KNl; k<=Grid.KNu; ++k){
      for(int z=1; z<=Grid.Nghost; z++){
	if(k >Grid.KNl && k < Grid.KNu){
	//create ghost cells on south side
	  norm_dir = - HALF*(Grid.nfaceS(i, Grid.JCl, k) +
			     Grid.nfaceS(i, Grid.JCl, k-1));
	  X_norm = ((Grid.Node[i][Grid.JNl+z][k].X-
		     Grid.Node[i][Grid.JNl][k].X)*norm_dir)*norm_dir;
	  X_tan  = (Grid.Node[i][Grid.JNl+z][k].X - 
		    Grid.Node[i][Grid.JNl][k].X)-X_norm;
	  
	  Grid.Node[i][Grid.JNl-z][k].X = Grid.Node[i][Grid.JNl][k].X - X_norm + X_tan; 
	  //create ghost cells on north side
	  norm_dir1 = - HALF*(Grid.nfaceN(i, Grid.JCl, k) +
			     Grid.nfaceN(i, Grid.JCl, k-1));
	  X_norm1 = ((Grid.Node[i][Grid.JNu][k].X-
		     Grid.Node[i][Grid.JNu-z][k].X)*norm_dir1)*norm_dir1;
	  X_tan1  = (Grid.Node[i][Grid.JNu][k].X - 
		    Grid.Node[i][Grid.JNu-z][k].X)-X_norm1;
	  
	  Grid.Node[i][Grid.JNu+z][k].X = Grid.Node[i][Grid.JNu][k].X + X_norm1 - X_tan1; 
	}
      }
    }
  }
  
  //the following code deals with placing ghost cells at the edges and corners
  for(int k = Grid.KNl-Grid.Nghost; k<=Grid.KNu+Grid.Nghost; ++k){
    for(int j = Grid.JNu+1; j<=Grid.JNu+Grid.Nghost; ++j){
      for(int z = 1; z<=Grid.Nghost;++z){
	Grid.Node[Grid.INu+z][j][k]= (Grid.Node[Grid.INu][j][k].X +
				      (Grid.Node[Grid.INu][j][k].X - Grid.Node[Grid.INu-z][j][k].X));
	Grid.Node[Grid.INl-z][j][k]= (Grid.Node[Grid.INl][j][k].X +
				      (Grid.Node[Grid.INl][j][k].X - Grid.Node[Grid.INl+z][j][k].X));
      }
    }
  }
  
  for(int k = Grid.KNl-Grid.Nghost; k<=Grid.KNu+Grid.Nghost; ++k){
    for(int j = Grid.JNl-Grid.Nghost; j<=Grid.JNu-1; ++j){
      for(int z = 1; z<=Grid.Nghost;++z){
	Grid.Node[Grid.INu+z][j][k]= (Grid.Node[Grid.INu][j][k].X +
				      (Grid.Node[Grid.INu][j][k].X - Grid.Node[Grid.INu-z][j][k].X));
	Grid.Node[Grid.INl-z][j][k]= (Grid.Node[Grid.INl][j][k].X +
				      (Grid.Node[Grid.INl][j][k].X - Grid.Node[Grid.INl+z][j][k].X));
      }
    }
  }

  for(int i = Grid.INl; i <=Grid.INu; ++i){
    for(int k = Grid.KNl-Grid.Nghost; k<=Grid.KNu-1; ++k){
      for(int z=1; z<=Grid.Nghost; ++z){
	Grid.Node[i][Grid.JNu+z][k]= (Grid.Node[i][Grid.JNu][k].X +
				      (Grid.Node[i][Grid.JNu][k].X - Grid.Node[i][Grid.JNu-z][k].X));
	Grid.Node[i][Grid.KNu-z][k]= (Grid.Node[i][Grid.JNl][k].X +
				      (Grid.Node[i][Grid.JNl][k].X - Grid.Node[i][Grid.JNl+z][k].X));
      }
    }
  }
  
  for(int i = Grid.INl; i <=Grid.INu; ++i){
    for(int k = Grid.KNu+1; k<=Grid.KNu+Grid.Nghost; ++k){
      for(int z=1; z<=Grid.Nghost; ++z){
	Grid.Node[i][Grid.JNu+z][k]= (Grid.Node[i][Grid.JNu][k].X +
				      (Grid.Node[i][Grid.JNu][k].X - Grid.Node[i][Grid.JNu-z][k].X));
	Grid.Node[i][Grid.KNu-z][k]= (Grid.Node[i][Grid.JNl][k].X +
				      (Grid.Node[i][Grid.JNl][k].X - Grid.Node[i][Grid.JNl+z][k].X));
      }
    }
  }
  
  for(int j = Grid.JNl; j <=Grid.JNu; ++j){
    for(int k = Grid.KNl-Grid.Nghost; k<=Grid.KNu-1; ++k){
      for(int z=1; z<=Grid.Nghost; ++z){
	Grid.Node[Grid.INu+z][j][k]= (Grid.Node[Grid.INu][j][k].X +
				      (Grid.Node[Grid.INu][j][k].X - Grid.Node[Grid.INu-z][j][k].X));
	Grid.Node[Grid.INu-z][Grid.KNu-z][k]= (Grid.Node[j][Grid.JNl][k].X +
					       (Grid.Node[Grid.INl][j][k].X - Grid.Node[Grid.INl+z][j][k].X));
      }
    }
  }
  
  for(int j = Grid.JNl; j <=Grid.JNu; ++j){
    for(int k = Grid.KNu+1; k<=Grid.KNu+Grid.Nghost; ++k){
      for(int z=1; z<=Grid.Nghost; ++z){
	Grid.Node[Grid.INu+z][j][k]= (Grid.Node[Grid.INu][j][k].X +
				      (Grid.Node[Grid.INu][j][k].X - Grid.Node[Grid.INu-z][j][k].X));
	Grid.Node[Grid.INu-z][Grid.KNu-z][k]= (Grid.Node[j][Grid.JNl][k].X +
					       (Grid.Node[Grid.INl][j][k].X - Grid.Node[Grid.INl+z][j][k].X));
      }
    }
  }
  
}

/********************************************************
 * Routine: Rotate_Hexa_Block                           *
 *                                                      *
 * Rotates the quadrilateral grid block.                *
 *                                                      *
 ********************************************************/
void Rotate_Hexa_Block(Grid3D_Hexa_Block &Grid,
                       const double &Angle, 
		       const double &Angle1, 
		       const double &Angle2) {
  double cos_angle, sin_angle; Vector3D X;
  cos_angle = cos(Angle); sin_angle = sin(Angle);
  //rotation about the Z-axis
  for(int k = Grid.KNl - Grid.Nghost; k <= Grid.KNu + Grid.Nghost; ++k){
    for (int j = Grid.JNl-Grid.Nghost ; j <= Grid.JNu+Grid.Nghost; ++j ) {
      for (int i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost; ++i ) {
	X.x = Grid.Node[i][j][k].X.x*cos_angle -
	  Grid.Node[i][j][k].X.y*sin_angle +
	  0*Grid.Node[i][j][k].X.z;
	X.y = Grid.Node[i][j][k].X.x*sin_angle +
	  Grid.Node[i][j][k].X.y*cos_angle +
	  0*Grid.Node[i][j][k].X.z;
	X.z = 0*Grid.Node[i][j][k].X.x*cos_angle +
	  0*Grid.Node[i][j][k].X.y*sin_angle +
	  Grid.Node[i][j][k].X.z; 
	Grid.Node[i][j][k].X = X;
      } /* endfor */
    } /* endfor */
  } /*endfor*/ 

  cos_angle=cos(Angle1);sin_angle= sin(Angle1);
  //rotation about the X-axis
  for(int k = Grid.KNl - Grid.Nghost; k <= Grid.KNu + Grid.Nghost; ++k){
    for (int j = Grid.JNl-Grid.Nghost ; j <= Grid.JNu+Grid.Nghost; ++j ) {
      for (int i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost; ++i ) {
	X.x = Grid.Node[i][j][k].X.x +
	  0*Grid.Node[i][j][k].X.y +
	  0*Grid.Node[i][j][k].X.z;
	X.y = 0*Grid.Node[i][j][k].X.x +
	  Grid.Node[i][j][k].X.y*cos_angle -
	  Grid.Node[i][j][k].X.z*sin_angle;
	X.z = 0*Grid.Node[i][j][k].X.x +
	  Grid.Node[i][j][k].X.y*sin_angle +
	  Grid.Node[i][j][k].X.z*cos_angle; 
	Grid.Node[i][j][k].X = X;
      } /* endfor */
    } /* endfor */
  } /*endfor*/

  cos_angle=cos(Angle2); sin_angle= sin(Angle2);
  //rotation about the Y-axis
  for(int k = Grid.KNl - Grid.Nghost; k <= Grid.KNu + Grid.Nghost; ++k){
    for (int j = Grid.JNl-Grid.Nghost ; j <= Grid.JNu+Grid.Nghost; ++j ) {
      for (int i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost; ++i ) {
	X.x = Grid.Node[i][j][k].X.x*cos_angle +
	  0*Grid.Node[i][j][k].X.y +
	  Grid.Node[i][j][k].X.z*sin_angle;
	X.y = 0*Grid.Node[i][j][k].X.x +
	  Grid.Node[i][j][k].X.y +
	  0*Grid.Node[i][j][k].X.z;
	X.z = -Grid.Node[i][j][k].X.x*sin_angle +
	  0*Grid.Node[i][j][k].X.y +
	  Grid.Node[i][j][k].X.z*cos_angle; 
	Grid.Node[i][j][k].X = X;
      } /* endfor */
    } /* endfor */
  } /*endfor*/
  
  for( int k = Grid.KCl-Grid.Nghost; k<= Grid.KCu + Grid.Nghost; ++k){
    for (int j = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j) {
      for (int i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
	Grid.Cell[i][j][k].Xc = Grid.centroid(i, j, k);
      } /* endfor */
    } /* endfor */
  } /*endfor */ 
 
}
