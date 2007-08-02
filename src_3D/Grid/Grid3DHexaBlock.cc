/* Grid3DHexaSingleBlock.cc:  Single-block subroutines for 
                              3D hexahedral block grid class. */

/* Include 3D hexahedral block grid type header file. */

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "Grid3DHexaBlock.h"
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
void Grid3D_Hexa_Block::Write_Hexa_Block(Grid3D_Hexa_Block &Grid, 
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
void Grid3D_Hexa_Block::Read_Hexa_Block(Grid3D_Hexa_Block &Grid,
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
void  Grid3D_Hexa_Block::Output_Tecplot(const int Block_Number,
                                        const int Output_Title,
                                        ostream &Out_File) {

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name()
                << ": 3D Structured Curvilinear Grid Block (Node Locations)"
                << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
	        << "\"z\" \n"  
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << INu - INl + 1 << " \\ \n"
                << "J = " << JNu - JNl + 1 << " \\ \n"
		<< "K = " << KNu - KNl + 1 << " \\ \n"
                << "F = POINT \n";

    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << INu - INl + 1 << " \\ \n"
                << "J = " << JNu - JNl + 1 << " \\ \n"
		<< "K = " << KNu - KNl + 1 << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    for ( int k = KNl ; k <= KNu ; ++k){
      for ( int j = JNl ; j <= JNu ; ++j ) {
	for ( int i = INl ; i <= INu ; ++i ) {
	  Out_File << " " << Node[i][j][k].X << "\n";
	} /* endfor */
      } /* endfor */
    }/*endfor*/
    Out_File << setprecision(6);

}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the nodes of the hexahedral mesh to the       *
 * specified output stream in a format suitable for     *
 * plotting the grid with TECPLOT.  Includes boundary   *
 * nodes.                                               *
 *                                                      *
 ********************************************************/
void  Grid3D_Hexa_Block::Output_Nodes_Tecplot(const int Block_Number,
                                              const int Output_Title,
                                              ostream &Out_File) {

   Out_File << setprecision(14);
   if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name()
               << ": 3D Structured Curvilinear Grid Block (Node Locations)"
               << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
		<< "\"z\" \n"
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << INu - INl + 1 + 4 << " \\ \n"
                << "J = " << JNu - JNl + 1 + 4 << " \\ \n"
		<< "K = " << KNu - KNl + 1 + 4 << " \\ \n"
                << "F = POINT \n";

    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << INu - INl + 1 + 4 << " \\ \n"
                << "J = " << JNu - JNl + 1 + 4 << " \\ \n"
		<< "K = " << KNu - KNl + 1 + 4 << " \\ \n" 
                << "F = POINT \n";
    } /* endif */
    
    for(int k = KNl-2 ; k <= KNu+2 ; ++k) {
      for (int j  = JNl-2 ; j <= JNu+2 ; ++j ) {
	for (int i = INl-2 ; i <= INu+2 ; ++i ) {
	  Out_File << " " << Node[i][j][k].X << "\n";
	} /* endfor */
      } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cells of the hexahedral mesh to the       *
 * specified output stream in a format suitable for     *
 * plotting the grid with TECPLOT.                      *
 *                                                      *
 ********************************************************/
void  Grid3D_Hexa_Block::Output_Cells_Tecplot(const int Block_Number,
                                              const int Output_Title,
                                              ostream &Out_File) {


    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name()
                << ": 3D Structured Curvilinear Grid Block (Cell Locations)"
                << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
		<< "\"z\" \n"
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << ICu - ICl + 1 << " \\ \n"
                << "J = " << JCu - JCl + 1 << " \\ \n"
		<< "K = " << KCu - KCl + 1 << " \\ \n"
                << "F = POINT \n";

    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << ICu - ICl + 1 << " \\ \n"
                << "J = " << JCu - JCl + 1 << " \\ \n"
		<< "K = " << KCu - KCl + 1 << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    for (int k = KCl ; k <= KCu ; ++k) {
      for (int j  = JCl ; j <= JCu ; ++j ) {
	for (int i = ICl ; i <= ICu ; ++i ) {
	  Out_File << " " << Cell[i][j][k].Xc << "\n";
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
void  Grid3D_Hexa_Block::Output_Gnuplot(const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {


    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "# " << CFFC_Name()
                << ": 3D Structured Curvilinear Grid Block (Node Locations)"
                << "\n"
	        << "# x(m), y(m), z(m)\n";
    } /* endif */
    
    for( int k = KNl ; k <= KNu ; ++k) {
      for (int j  = JNl ; j <= JNu ; ++j ) {
	for (int i = INl ; i <= INu ; ++i ) {
	  Out_File << " " << Node[i][j][k].X << "\n";
	} /* endfor */
	Out_File << "\n";
      } /* endfor */
      Out_File << "\n";
    }/*endfor*/
    for (int i  = INl ; i <= INu ; ++i ) {
      for (int j = JNl ; j <= JNu ; ++j ) {
	for (int k = KNl ; k <= KNu ; ++k ){
	  Out_File << " " << Node[i][j][k].X << "\n";
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
 * Updates the cell information for the hexahedral mesh block
 *                                                           
 *****************************************************************/

void  Grid3D_Hexa_Block::Update_Cells(void) {
  for(int k = KCl-Nghost; k <=KCu +Nghost ; ++k){ 
    for( int j = JCl-Nghost ; j <= JCu +Nghost; ++j) {
      for ( int i = ICl-Nghost ; i <= ICu+Nghost; ++i) {
	Cell[i][j][k].I = i ;
	Cell[i][j][k].J = j ;
	Cell[i][j][k].K = k ;
	Cell[i][j][k].Xc = centroid(i, j, k);
        Cell[i][j][k].V = volume(Cell[i][j][k]);

        
      } /* endfor */
    } /* endfor */
  } /* endfor */
}

/**************************************************************
Update exterior nodes
**************************************************************/
void Grid3D_Hexa_Block::Update_Exterior_Nodes(){
  Vector3D norm_dir, X_norm, X_tan;
  Vector3D norm_dir1,X_norm1,X_tan1;
  for(int k = KNl; k<=KNu; ++k){
    for(int j = JNl; j<=JNu; ++j){
      for(int z=1; z<=Nghost; z++){
	if(j>JNl && j < JNu){
	  //create ghost cells on west side
	  norm_dir = - HALF*(nfaceW(ICl, j, k) +
			     nfaceW(ICl, j-1, k));
	  X_norm = ((Node[INl+z][j][k].X-
		     Node[INl][j][k].X)*norm_dir)*norm_dir;
	  X_tan  = (Node[INl+z][j][k].X - 
		    Node[INl][j][k].X)-X_norm;
	  
	  Node[INl-z][j][k].X = Node[INl][j][k].X - X_norm + X_tan; 
	  //create ghost cells on east side
	  norm_dir1 = - HALF*(nfaceE(ICl, j, k) +
			     nfaceE(ICl, j-1, k));
	  X_norm1 = ((Node[INu][j][k].X-
		     Node[INu-z][j][k].X)*norm_dir1)*norm_dir1;
	  X_tan1  = (Node[INu][j][k].X - 
		    Node[INu-z][j][k].X)-X_norm1;
	  
	  Node[INu+z][j][k].X = Node[INu][j][k].X + X_norm1 - X_tan1; 
	}
      }
    }
  }
  
  for(int j = JNl; j<=JNu; ++j){
    for(int i = INl; i<=INu; ++i){
      for(int z=1; z<=Nghost; z++){
	if(i>INl && i < INu){
	  //create ghost cells on bottom side
	  norm_dir = - HALF*(nfaceBot(i, j, KCl) +
			     nfaceBot(i-1, j, KCl));
	  X_norm = ((Node[i][j][KNl+z].X-
		     Node[i][j][KNl].X)*norm_dir)*norm_dir;
	  X_tan  = (Node[i][j][KNl+z].X - 
		    Node[i][j][KNl].X)-X_norm;
	  
	  Node[i][j][KNl-z].X = Node[i][j][KNl].X - X_norm + X_tan; 
	  //create ghost cells on top side
	  norm_dir1 = - HALF*(nfaceTop(i, j, KCl) +
			     nfaceTop(i-1, j, KCl));
	  X_norm1 = ((Node[i][j][KNu].X-
		     Node[i][j][KNu-z].X)*norm_dir1)*norm_dir1;
	  X_tan1  = (Node[i][j][KNu].X - 
		    Node[i][j][KNu-z].X)-X_norm1;
	  
	  Node[i][j][KNu+z].X = Node[i][j][KNu].X + X_norm1 - X_tan1; 
	}	
      }
    }
  }
  for(int i = INl; i<=INu; ++i){
    for(int k = KNl; k<=KNu; ++k){
      for(int z=1; z<=Nghost; z++){
	if(k >KNl && k < KNu){
	//create ghost cells on south side
	  norm_dir = - HALF*(nfaceS(i, JCl, k) +
			     nfaceS(i, JCl, k-1));
	  X_norm = ((Node[i][JNl+z][k].X-
		     Node[i][JNl][k].X)*norm_dir)*norm_dir;
	  X_tan  = (Node[i][JNl+z][k].X - 
		    Node[i][JNl][k].X)-X_norm;
	  
	  Node[i][JNl-z][k].X = Node[i][JNl][k].X - X_norm + X_tan; 
	  //create ghost cells on north side
	  norm_dir1 = - HALF*(nfaceN(i, JCl, k) +
			     nfaceN(i, JCl, k-1));
	  X_norm1 = ((Node[i][JNu][k].X-
		     Node[i][JNu-z][k].X)*norm_dir1)*norm_dir1;
	  X_tan1  = (Node[i][JNu][k].X - 
		    Node[i][JNu-z][k].X)-X_norm1;
	  
	  Node[i][JNu+z][k].X = Node[i][JNu][k].X + X_norm1 - X_tan1; 
	}
      }
    }
  }
  
  //the following code deals with placing ghost cells at the edges and corners
  for(int k = KNl-Nghost; k<=KNu+Nghost; ++k){
    for(int j = JNu+1; j<=JNu+Nghost; ++j){
      for(int z = 1; z<=Nghost;++z){
	Node[INu+z][j][k]= (Node[INu][j][k].X +
				      (Node[INu][j][k].X - Node[INu-z][j][k].X));
	Node[INl-z][j][k]= (Node[INl][j][k].X +
				      (Node[INl][j][k].X - Node[INl+z][j][k].X));
      }
    }
  }
  
  for(int k = KNl-Nghost; k<=KNu+Nghost; ++k){
    for(int j = JNl-Nghost; j<=JNu-1; ++j){
      for(int z = 1; z<=Nghost;++z){
	Node[INu+z][j][k]= (Node[INu][j][k].X +
				      (Node[INu][j][k].X - Node[INu-z][j][k].X));
	Node[INl-z][j][k]= (Node[INl][j][k].X +
				      (Node[INl][j][k].X - Node[INl+z][j][k].X));
      }
    }
  }

  for(int i = INl; i <=INu; ++i){
    for(int k = KNl-Nghost; k<=KNu-1; ++k){
      for(int z=1; z<=Nghost; ++z){
	Node[i][JNu+z][k]= (Node[i][JNu][k].X +
				      (Node[i][JNu][k].X - Node[i][JNu-z][k].X));
	Node[i][KNu-z][k]= (Node[i][JNl][k].X +
				      (Node[i][JNl][k].X - Node[i][JNl+z][k].X));
      }
    }
  }
  
  for(int i = INl; i <=INu; ++i){
    for(int k = KNu+1; k<=KNu+Nghost; ++k){
      for(int z=1; z<=Nghost; ++z){
	Node[i][JNu+z][k]= (Node[i][JNu][k].X +
				      (Node[i][JNu][k].X - Node[i][JNu-z][k].X));
	Node[i][KNu-z][k]= (Node[i][JNl][k].X +
				      (Node[i][JNl][k].X - Node[i][JNl+z][k].X));
      }
    }
  }
  
  for(int j = JNl; j <=JNu; ++j){
    for(int k = KNl-Nghost; k<=KNu-1; ++k){
      for(int z=1; z<=Nghost; ++z){
	Node[INu+z][j][k]= (Node[INu][j][k].X +
				      (Node[INu][j][k].X - Node[INu-z][j][k].X));
	Node[INu-z][KNu-z][k]= (Node[j][JNl][k].X +
					       (Node[INl][j][k].X - Node[INl+z][j][k].X));
      }
    }
  }
  
  for(int j = JNl; j <=JNu; ++j){
    for(int k = KNu+1; k<=KNu+Nghost; ++k){
      for(int z=1; z<=Nghost; ++z){
	Node[INu+z][j][k]= (Node[INu][j][k].X +
				      (Node[INu][j][k].X - Node[INu-z][j][k].X));
	Node[INu-z][KNu-z][k]= (Node[j][JNl][k].X +
					       (Node[INl][j][k].X - Node[INl+z][j][k].X));
      }
    }
  }
  
}

/********************************************************
 * Routine: Rotate_Hexa_Block                           *
 *                                                      *
 * Rotates the hexahedral grid block.                   *
 *                                                      *
 ********************************************************/
void  Grid3D_Hexa_Block::Rotate_Hexa_Block(const double &Angle, 
                                           const double &Angle1, 
                                           const double &Angle2) {
   double cos_angle, sin_angle; Vector3D X;
  cos_angle = cos(Angle); sin_angle = sin(Angle);
  //rotation about the Z-axis
  for(int k = KNl - Nghost; k <= KNu + Nghost; ++k){
    for (int j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
      for (int i = INl-Nghost ; i <= INu+Nghost; ++i ) {
	X.x = Node[i][j][k].X.x*cos_angle -
	  Node[i][j][k].X.y*sin_angle +
	  0*Node[i][j][k].X.z;
	X.y = Node[i][j][k].X.x*sin_angle +
	  Node[i][j][k].X.y*cos_angle +
	  0*Node[i][j][k].X.z;
	X.z = 0*Node[i][j][k].X.x*cos_angle +
	  0*Node[i][j][k].X.y*sin_angle +
	  Node[i][j][k].X.z; 
	Node[i][j][k].X = X;
      } /* endfor */
    } /* endfor */
  } /*endfor*/ 

  cos_angle=cos(Angle1);sin_angle= sin(Angle1);
  //rotation about the X-axis
  for(int k = KNl - Nghost; k <= KNu + Nghost; ++k){
    for (int j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
      for (int i = INl-Nghost ; i <= INu+Nghost; ++i ) {
	X.x = Node[i][j][k].X.x +
	  0*Node[i][j][k].X.y +
	  0*Node[i][j][k].X.z;
	X.y = 0*Node[i][j][k].X.x +
	  Node[i][j][k].X.y*cos_angle -
	  Node[i][j][k].X.z*sin_angle;
	X.z = 0*Node[i][j][k].X.x +
	  Node[i][j][k].X.y*sin_angle +
	  Node[i][j][k].X.z*cos_angle; 
	Node[i][j][k].X = X;
      } /* endfor */
    } /* endfor */
  } /*endfor*/

  cos_angle=cos(Angle2); sin_angle= sin(Angle2);
  //rotation about the Y-axis
  for(int k = KNl - Nghost; k <= KNu + Nghost; ++k){
    for (int j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
      for (int i = INl-Nghost ; i <= INu+Nghost; ++i ) {
	X.x = Node[i][j][k].X.x*cos_angle +
	  0*Node[i][j][k].X.y +
	  Node[i][j][k].X.z*sin_angle;
	X.y = 0*Node[i][j][k].X.x +
	  Node[i][j][k].X.y +
	  0*Node[i][j][k].X.z;
	X.z = -Node[i][j][k].X.x*sin_angle +
	  0*Node[i][j][k].X.y +
	  Node[i][j][k].X.z*cos_angle; 
	Node[i][j][k].X = X;
      } /* endfor */
    } /* endfor */
  } /*endfor*/
  
  for( int k = KCl-Nghost; k<= KCu + Nghost; ++k){
    for (int j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
      for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
	Cell[i][j][k].Xc = centroid(i, j, k);
      } /* endfor */
    } /* endfor */
  } /*endfor */ 
 
}



/*************************************************************************
 * Grid3D_Hexa_Block -- External subroutines for single grid block.      *
 *************************************************************************/
/********************************************************
 * Routine: Create_Hexa_Block                           *
 *                                                      *
 * Create hexa- grid block for a Cartesian              *
 * mesh defined on a cube with given length, width,     *
 *  and height (origion 0,0,0)                          *
 ********************************************************/
void  Grid3D_Hexa_Block::Create_Hexa_Block(const double Length,
                                           const double Width,
                                           const double Height,
                                           const double x_orig,
                                           const double y_orig,
                                           const double z_orig,
                                           const double alpha,
                                           const double beta,
                                           const double gamma,
                                           const int  BCtype_top,
                                           const int  BCtype_bottom,
                                           const int  BCtype_north,
                                           const int  BCtype_south,
                                           const int  BCtype_west,
                                           const int  BCtype_east,
                                           const int Number_of_Cells_Idir,
                                           const int Number_of_Cells_Jdir,
                                           const int Number_of_Cells_Kdir,
                                           const int Number_of_Ghost_Cells) {

   int i, j, k;
   double dx, dy, dz;
   double xx = 0.0;
   
   dx = Length/double(Number_of_Cells_Idir);
   dy = Width/double(Number_of_Cells_Jdir);
   dz = Height/double(Number_of_Cells_Kdir);
  
   // Generate the mesh
   for ( k = KNl-Nghost ; k <= KNu+Nghost ; ++k) 
      for ( j = JNl-Nghost ; j <= JNu+Nghost ; ++j) 
         for ( i = INl-Nghost ; i <= INu+Nghost ; ++i) {
            
            // stretching to both x ends
            if(alpha > 0.0){
               // stretching to both walls in z direction. overwrite  Node[i][j][k].X.z.
               xx = double(i-INl)/double(INu-INl);
               Node[i][j][k].X.x = x_orig + Length*StretchingFcn(xx, alpha, ZERO, STRETCHING_FCN_MINMAX_CLUSTERING);
	    }else{
               Node[i][j][k].X.x = x_orig + (i-Nghost)*dx;
                            
	    }
            // stretching to both y ends
            if(beta > 0.0){
               
               // stretching to both walls in z direction. overwrite  Node[i][j][k].X.z.
               xx = double(j-JNl)/double(JNu-JNl);
               Node[i][j][k].X.y = y_orig + Width*StretchingFcn(xx, beta, ZERO, STRETCHING_FCN_MINMAX_CLUSTERING);
               
            }else{
               Node[i][j][k].X.y = y_orig + (j-Nghost)*dy;
	    }
            //stretching to both z ends
            if(gamma > 0.0){
               // stretching to both walls in z direction. overwrite  Node[i][j][k].X.z.
               xx = double(k-KNl)/double(KNu-KNl);
               Node[i][j][k].X.z = z_orig + Height*StretchingFcn(xx, gamma, ZERO, STRETCHING_FCN_MINMAX_CLUSTERING);
	    }else{
               Node[i][j][k].X.z = z_orig + (k-Nghost)*dz;
               
	    }
         }
   
   // Update the ghost nodes coordinates according to the stretching factors;
   if(alpha > 0.0){
      for ( k = KNl-Nghost ; k <= KNu+Nghost ; ++k) 
         for ( j = JNl-Nghost ; j <= JNu+Nghost ; ++j) 
            for ( i = INl-Nghost ; i <= INu+Nghost ; ++i) {
               // update the ghost nodes coordinates;
               if(i==INu+1){
                  Node[i][j][k].X.x =  Node[INu][j][k].X.x +( Node[INu][j][k].X.x - Node[INu-1][j][k].X.x);
               }
               if(i==INu+2){
                  Node[i][j][k].X.x =  Node[INu][j][k].X.x +( Node[INu][j][k].X.x - Node[INu-2][j][k].X.x);
               }
               
               if(i==INl - 1){
                  Node[i][j][k].X.x =  Node[INl][j][k].X.x +( Node[INl][j][k].X.x - Node[INl+1][j][k].X.x);
               }
               if(i==INl - 2){
                  Node[i][j][k].X.x =  Node[INl][j][k].X.x +( Node[INl][j][k].X.x - Node[INl+2][j][k].X.x);
               }// 
            }
   }
   if(beta > 0.0){
      for ( k = KNl-Nghost ; k <= KNu+Nghost ; ++k) 
         for ( j = JNl-Nghost ; j <= JNu+Nghost ; ++j) 
            for ( i = INl-Nghost ; i <= INu+Nghost ; ++i) {
               // update the ghost nodes coordinates;
               if(j==JNu+1){
                  Node[i][j][k].X.y =  Node[i][JNu][k].X.y +( Node[i][JNu][k].X.y - Node[i][JNu-1][k].X.y);
               }
               if(j==JNu+2){
                  Node[i][j][k].X.y =  Node[i][JNu][k].X.y +( Node[i][JNu][k].X.y - Node[i][JNu-2][k].X.y);
               }
               
               if(j==JNl - 1){
                  Node[i][j][k].X.y =  Node[i][JNl][k].X.y +( Node[i][JNl][k].X.y - Node[i][JNl+1][k].X.y);
               }
               if(j==JNl - 2){
                  Node[i][j][k].X.y =  Node[i][JNl][k].X.y +( Node[i][JNl][k].X.y - Node[i][JNl+2][k].X.y);
               }// 
            }
   }
              
   if(gamma > 0.0){
      for ( k = KNl-Nghost ; k <= KNu+Nghost ; ++k) 
         for ( j = JNl-Nghost ; j <= JNu+Nghost ; ++j) 
            for ( i = INl-Nghost ; i <= INu+Nghost ; ++i) {
               // update the ghost nodes coordinates;
               if(k==KNu+1){
                  Node[i][j][k].X.z =  Node[i][j][KNu].X.z +( Node[i][j][KNu].X.z - Node[i][j][KNu-1].X.z);
               }
               if(k==KNu+2){
                  Node[i][j][k].X.z =  Node[i][j][KNu].X.z +( Node[i][j][KNu].X.z - Node[i][j][KNu-2].X.z);
               }
               
               if(k==KNl - 1){
                  Node[i][j][k].X.z =  Node[i][j][KNl].X.z +( Node[i][j][KNl].X.z - Node[i][j][KNl+1].X.z);
               }
               if(k==KNl - 2){
                  Node[i][j][k].X.z =  Node[i][j][KNl].X.z +( Node[i][j][KNl].X.z - Node[i][j][KNl+2].X.z);
               }// 
            }
   }
   

   // update cell coordinates
   Update_Cells();
   
   //store the boundary conditions
   for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) 
      for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
         
         BCtypeT[i][j] =  BCtype_top;
         BCtypeB[i][j] =  BCtype_bottom;
         
      }
   for ( k = KCl-Nghost ; k <= KCu+Nghost ; ++k) 
      for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j){
         BCtypeW[j][k] =  BCtype_west;
         BCtypeE[j][k] =  BCtype_east;
               
      }
   for ( k = KCl-Nghost ; k <= KCu+Nghost ; ++k) 
      for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i){
         BCtypeN[i][k] =  BCtype_north;
         BCtypeS[i][k] =  BCtype_south;
         
      }
   
  
   
}

/********************************************************
 * Routine: Copy_Quad_Block                             *
 *                                                      *
 * Copies the node and cell locations of Hexahedral     *
 * mesh block Grid2 to Grid1.                           *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Copy_Hexa_Block(Grid3D_Hexa_Block *Grid1,
                     Grid3D_Hexa_Block *Grid2) {
   

       int i, j,k, ni, nj, nk;
 
    /* Allocate (re-allocate) memory for the cells and nodes 
       of the hexarilateral mesh block Grid1 as necessary. */
 
      
    if (Grid1->NNi != Grid2->NNi || Grid1->NNj != Grid2->NNj ||
        Grid1->NNk != Grid2->NNk || Grid1->NCk != Grid2->NCk ||
        Grid1->NCi != Grid2->NCi || Grid1->NCj != Grid2->NCj ) { 
 
       if (Grid1->Node != NULL && Grid1->Cell != NULL) { 
          Grid1->deallocate();
       } else if (Grid1->Node != NULL) {
          Grid1->deallocateNodes();
       } else if (Grid1->Cell != NULL) {
          Grid1->deallocateCells();
       } /* endif */

 
       if (Grid2->Node != NULL) Grid1->allocate(Grid2->NCi- 2*Grid2->Nghost, 
                                              Grid2->NCj-2*Grid2->Nghost, 
                                              Grid2->NCk-2*Grid2->Nghost,
                                              Grid2->Nghost);
 

    } /* endif */

    /* Copy the node locations of grid block Grid2 to Grid1-> */

    if (Grid2->Node != NULL) {
      for (k  = Grid2->KNl- Grid2->Nghost ; k <= Grid2->KNu+ Grid2->Nghost ; ++k ) 
	for (j  = Grid2->JNl- Grid2->Nghost ; j <= Grid2->JNu+ Grid2->Nghost ; ++j ) {
	  for ( i = Grid2->INl- Grid2->Nghost ; i <= Grid2->INu+ Grid2->Nghost ; ++i ) {
	    Grid1->Node[i][j][k].X = Grid2->Node[i][j][k].X;
	  } /* endfor */
	} /* endfor */
    } /* endif */
    
    /* Copy the cell values of grid block Grid2 to Grid1-> */

    if (Grid2->Node != NULL) {
      for ( k = Grid2->KCl- Grid2->Nghost; k <= Grid2->KCu+ Grid2->Nghost ; ++k) 
        for ( j = Grid2->JCl- Grid2->Nghost; j <= Grid2->JCu+ Grid2->Nghost ; ++j) {
	  for ( i = Grid2->ICl- Grid2->Nghost ; i <= Grid2->ICu+ Grid2->Nghost ; ++i) {
	    Grid1->Cell[i][j][k].I  = Grid2->Cell[i][j][k].I;
	    Grid1->Cell[i][j][k].J  = Grid2->Cell[i][j][k].J;
	    Grid1->Cell[i][j][k].Xc = Grid2->Cell[i][j][k].Xc;
	    Grid1->Cell[i][j][k].V  = Grid2->Cell[i][j][k].V;
	  } /* endfor */
        } /* endfor */
    } /* endif */

    /* Copy the boundary condition type info of grid block 
       Grid2 to Grid1-> */

    if (Grid2->Node != NULL) {
      for ( i = Grid2->ICl- Grid2->Nghost ; i <= Grid2->ICu+ Grid2->Nghost ; ++i) 
	for ( k = Grid2->KCl- Grid2->Nghost ; k <= Grid2->KCu+ Grid2->Nghost ; ++k) {
	  Grid1->BCtypeN[i][k] = Grid2->BCtypeN[i][k];
	  Grid1->BCtypeS[i][k] = Grid2->BCtypeS[i][k];
        } /* endfor */
      for ( j = Grid2->JCl- Grid2->Nghost ; j <= Grid2->JCu+ Grid2->Nghost ; ++j) {
        for ( k = Grid2->KCl- Grid2->Nghost ; k <= Grid2->KCu+ Grid2->Nghost ; ++k) {
	  Grid1->BCtypeE[j][k] = Grid2->BCtypeE[j][k];
	  Grid1->BCtypeW[j][k] = Grid2->BCtypeW[j][k];
        } /* endfor */
        for ( i = Grid2->ICl- Grid2->Nghost ; i <= Grid2->ICu+ Grid2->Nghost ; ++i) 
	  for ( j = Grid2->JCl- Grid2->Nghost ; j <= Grid2->JCu+ Grid2->Nghost ; ++j) {
	    Grid1->BCtypeT[i][j] = Grid2->BCtypeT[i][j];
	    Grid1->BCtypeB[i][j] = Grid2->BCtypeB[i][j];
	  } /* endfor */
// 	Grid1->BCs_N = Grid2->BCs_N;
// 	Grid1->BCs_S = Grid2->BCs_S;
// 	Grid1->BCs_E = Grid2->BCs_E;
// 	Grid1->BCs_W = Grid2->BCs_W;
// 	Grid1->BCs_T = Grid2->BCs_T;
// 	Grid1->BCs_B = Grid2->BCs_B;
	
      } /* endif */
      
    }

}
