/*! \file    Grid3DHexaBlock.cc
 *  \brief   Defines member functions for 3D hexahedral grid block class. */


/* Include 3D hexahedral block grid type header file. */

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "Grid3DHexaBlock.h"
#endif // _GRID3D_HEXA_BLOCK_INCLUDED

/*************************************************************************
 * Grid3D_Hexa_Block -- Member functions for single grid block.          *
 *************************************************************************/

/********************************************************
 * Routine: Create_Block                                *
 *                                                      *
 * Creates a hexahedral grid block for a Cartesian      *
 * mesh defined on a cube with given length, width,     *
 * and height (origion 0,0,0)                           *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Create_Block(const double &Length,
                                     const double &Width,
                                     const double &Height,
                                     const double &x_orig,
                                     const double &y_orig,
                                     const double &z_orig,
                                     const double &alpha,
                                     const double &beta,
                                     const double &gamma,
                                     const int BCtype_top,
                                     const int BCtype_bottom,
                                     const int BCtype_north,
                                     const int BCtype_south,
                                     const int BCtype_west,
                                     const int BCtype_east,
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
   for (k = KNl-Nghost ; k <= KNu+Nghost ; ++k) {
      for (j = JNl-Nghost ; j <= JNu+Nghost ; ++j) {
         for (i = INl-Nghost ; i <= INu+Nghost ; ++i) {
            // stretching to both x ends
            if (alpha > 0.0){
               // stretching to both walls in z direction. overwrite  Node[i][j][k].X.z.
               xx = double(i-INl)/double(INu-INl);
               Node[i][j][k].X.x = x_orig + Length*
                                   StretchingFcn(xx, alpha, ZERO, STRETCHING_FCN_MINMAX_CLUSTERING);
	    } else {
               Node[i][j][k].X.x = x_orig + (i-Nghost)*dx;
                            
	    }

            // stretching to both y ends
            if (beta > 0.0) {
               // stretching to both walls in z direction. overwrite  Node[i][j][k].X.z.
               xx = double(j-JNl)/double(JNu-JNl);
               Node[i][j][k].X.y = y_orig + Width*
                                   StretchingFcn(xx, beta, ZERO, STRETCHING_FCN_MINMAX_CLUSTERING);
            } else {
               Node[i][j][k].X.y = y_orig + (j-Nghost)*dy;
	    }

            //stretching to both z ends
            if (gamma > 0.0) {
               // stretching to both walls in z direction. overwrite  Node[i][j][k].X.z.
               xx = double(k-KNl)/double(KNu-KNl);
               Node[i][j][k].X.z = z_orig + Height*
                                   StretchingFcn(xx, gamma, ZERO, STRETCHING_FCN_MINMAX_CLUSTERING);
	    } else {
               Node[i][j][k].X.z = z_orig + (k-Nghost)*dz;
	    }
         } /* endfor */
      } /* endfor */
   } /* endfor */
   
   // Update the ghost nodes coordinates according to the stretching factors;
   if (alpha > 0.0) {
      for (k = KNl-Nghost ; k <= KNu+Nghost ; ++k) {
         for (j = JNl-Nghost ; j <= JNu+Nghost ; ++j) {
            for (i = INl-Nghost ; i <= INu+Nghost ; ++i) {
               if(i==INu+1){
                  Node[i][j][k].X.x = Node[INu][j][k].X.x +
                                      (Node[INu][j][k].X.x - Node[INu-1][j][k].X.x);
               }
               if(i==INu+2){
                  Node[i][j][k].X.x = Node[INu][j][k].X.x +
                                      (Node[INu][j][k].X.x - Node[INu-2][j][k].X.x);
               }
               
               if(i==INl - 1){
                  Node[i][j][k].X.x = Node[INl][j][k].X.x +
                                      (Node[INl][j][k].X.x - Node[INl+1][j][k].X.x);
               }
               if(i==INl - 2){
                  Node[i][j][k].X.x = Node[INl][j][k].X.x +
                                      (Node[INl][j][k].X.x - Node[INl+2][j][k].X.x);
               }
            } /* endfor */
	 } /* endfor */
      } /* endfor */
   } /* endif */

   if (beta > 0.0) {
      for (k = KNl-Nghost ; k <= KNu+Nghost ; ++k) {
         for (j = JNl-Nghost ; j <= JNu+Nghost ; ++j) {
            for (i = INl-Nghost ; i <= INu+Nghost ; ++i) {
               if(j==JNu+1){
                  Node[i][j][k].X.y = Node[i][JNu][k].X.y +
                                      (Node[i][JNu][k].X.y - Node[i][JNu-1][k].X.y);
               }
               if(j==JNu+2){
                  Node[i][j][k].X.y = Node[i][JNu][k].X.y +
                                      (Node[i][JNu][k].X.y - Node[i][JNu-2][k].X.y);
               }
               
               if(j==JNl - 1){
                  Node[i][j][k].X.y = Node[i][JNl][k].X.y +
                                      (Node[i][JNl][k].X.y - Node[i][JNl+1][k].X.y);
               }
               if(j==JNl - 2){
                  Node[i][j][k].X.y = Node[i][JNl][k].X.y +
                                      (Node[i][JNl][k].X.y - Node[i][JNl+2][k].X.y);
               } 
            } /* endfor */
	 } /* endfor */
      } /* endfor */
   } /* endif */
              
   if (gamma > 0.0) {
      for (k = KNl-Nghost ; k <= KNu+Nghost ; ++k) {
	 for (j = JNl-Nghost ; j <= JNu+Nghost ; ++j) {
            for (i = INl-Nghost ; i <= INu+Nghost ; ++i) {
               if(k==KNu+1){
                  Node[i][j][k].X.z = Node[i][j][KNu].X.z +
                                      (Node[i][j][KNu].X.z - Node[i][j][KNu-1].X.z);
               }
               if(k==KNu+2){
                  Node[i][j][k].X.z = Node[i][j][KNu].X.z +
                                      (Node[i][j][KNu].X.z - Node[i][j][KNu-2].X.z);
               }
               
               if(k==KNl - 1){
                  Node[i][j][k].X.z = Node[i][j][KNl].X.z +
                                      (Node[i][j][KNl].X.z - Node[i][j][KNl+1].X.z);
               }
               if(k==KNl - 2){
                  Node[i][j][k].X.z = Node[i][j][KNl].X.z +
                                      (Node[i][j][KNl].X.z - Node[i][j][KNl+2].X.z);
               } 
            } /* endfor */
	 } /* endfor */
      } /* endfor */
   } /* endif */
   

   // update cell coordinates
   Update_Cells();
   
   //store the boundary conditions
   for (j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
      for (i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
         BCtypeT[i][j] =  BCtype_top;
         BCtypeB[i][j] =  BCtype_bottom;
      } /* endfor */
   } /* endfor */

   for (k = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
      for (j = JCl-Nghost ; j <= JCu+Nghost ; ++j){
         BCtypeW[j][k] =  BCtype_west;
         BCtypeE[j][k] =  BCtype_east;
      } /* endfor */
   } /* endfor */

   for (k = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
      for (i = ICl-Nghost ; i <= ICu+Nghost ; ++i){
         BCtypeN[i][k] =  BCtype_north;
         BCtypeS[i][k] =  BCtype_south;
      } /* endfor */
   } /* endfor */

}

/********************************************************
 * Routine: Copy                                        *
 *                                                      *
 * Copy the node and cell locations of hexahedral       *
 * mesh block Grid2 to current grid block.              *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Copy(Grid3D_Hexa_Block &Grid2) {

    if (Grid2.Allocated) {

       /* Allocate (re-allocate) memory for the cells and nodes 
          of the hexahedral mesh block as necessary. */
      
       if (NNi != Grid2.NNi || NNj != Grid2.NNj ||
           NNk != Grid2.NNk || NCi != Grid2.NCi ||
           NCj != Grid2.NCj || NCk != Grid2.NCk ||
           Nghost != Grid2.Nghost) { 
          if (Allocated) { 
             deallocate();
          } /* endif */

          allocate(Grid2.NCi-2*Grid2.Nghost, 
                   Grid2.NCj-2*Grid2.Nghost, 
                   Grid2.NCk-2*Grid2.Nghost,
                   Grid2.Nghost);
       } /* endif */

       /* Copy the node locations of grid block Grid2. */

       for (int k = Grid2.KNl-Grid2.Nghost; k <= Grid2.KNu+Grid2.Nghost; ++k) {
	  for (int j = Grid2.JNl-Grid2.Nghost; j <= Grid2.JNu+Grid2.Nghost; ++j) {
	     for (int i = Grid2.INl-Grid2.Nghost; i <= Grid2.INu+Grid2.Nghost; ++i) {
	        Node[i][j][k].X = Grid2.Node[i][j][k].X;
	     } /* endfor */
	  } /* endfor */
       } /* endfor */
    
       /* Copy the cell values of grid block Grid2. */

       for (int k = Grid2.KCl-Grid2.Nghost; k <= Grid2.KCu+Grid2.Nghost; ++k) {
          for (int j = Grid2.JCl-Grid2.Nghost; j <= Grid2.JCu+Grid2.Nghost; ++j) {
	     for (int i = Grid2.ICl-Grid2.Nghost; i <= Grid2.ICu+Grid2.Nghost; ++i) {
	        Cell[i][j][k].I  = Grid2.Cell[i][j][k].I;
	        Cell[i][j][k].J  = Grid2.Cell[i][j][k].J;
	        Cell[i][j][k].K  = Grid2.Cell[i][j][k].K;
	        Cell[i][j][k].Xc = Grid2.Cell[i][j][k].Xc;
	        Cell[i][j][k].V  = Grid2.Cell[i][j][k].V;
             Cell[i][j][k].Jacobian  = Grid2.Cell[i][j][k].Jacobian;

	     } /* endfor */
          } /* endfor */
       } /* endfor */

       /* Copy the boundary condition type info of grid block Grid2. */

       for (int k = Grid2.KCl-Grid2.Nghost; k <= Grid2.KCu+Grid2.Nghost; ++k) {
          for (int i = Grid2.ICl-Grid2.Nghost; i <= Grid2.ICu+Grid2.Nghost; ++i) {
	     BCtypeN[i][k] = Grid2.BCtypeN[i][k];
	     BCtypeS[i][k] = Grid2.BCtypeS[i][k];
          } /* endfor */
       } /* endfor */

       for (int k = Grid2.KCl-Grid2.Nghost; k <= Grid2.KCu+Grid2.Nghost; ++k) {
          for (int j = Grid2.JCl-Grid2.Nghost; j <= Grid2.JCu+Grid2.Nghost; ++j) {
	     BCtypeE[j][k] = Grid2.BCtypeE[j][k];
	     BCtypeW[j][k] = Grid2.BCtypeW[j][k];
          } /* endfor */
       } /* endfor */


       for (int j = Grid2.JCl-Grid2.Nghost; j <= Grid2.JCu+Grid2.Nghost; ++j) {
          for (int i = Grid2.ICl-Grid2.Nghost; i <= Grid2.ICu+Grid2.Nghost; ++i) {
	     BCtypeT[i][j] = Grid2.BCtypeT[i][j];
	     BCtypeB[i][j] = Grid2.BCtypeB[i][j];
	   } /* endfor */
       } /* endfor */
    } /* endif */

    // Copy the number of Gauss quadrature points
    CopyNumberOfGaussQuadraturePoints(Grid2.getNumGQP());

}

/********************************************************
 * Routine: Broadcast                                   *
 *                                                      *
 * Broadcast the node and cell locations of hexahedral  *
 * mesh block.                                          *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Broadcast(void) {

#ifdef _MPI_VERSION

    int ni, nj, nk, ng, mesh_allocated, buffer_size;
    int *i_buffer;
    double *buffer;

    /* Broadcast the number of regular and ghost cells 
       for the grid block. */

    if (CFFC_Primary_MPI_Processor()) {
      ni = NCi;
      nj = NCj;
      nk = NCk;
      ng = Nghost;
      mesh_allocated = Allocated;
    } /* endif */

    MPI::COMM_WORLD.Bcast(&ni, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nj, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nk, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&ng, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&mesh_allocated, 1, MPI::INT, 0);

    /* On non-primary MPI processors, allocate (re-allocate) 
       memory for the cells and nodes of the hexihedral
       mesh block as necessary. */

    if (!CFFC_Primary_MPI_Processor()) {
       if (mesh_allocated && (NCi != ni || NCj != nj || NCk != nk || Nghost != ng)) { 
          if (Allocated) { 
             deallocate();
          } /* endif */

          allocate(ni-2*ng, nj-2*ng, nk-2*ng, ng);

       } /* endif */
    } /* endif */

    /* Broadcast the node locations for grid block. */

    if (mesh_allocated) {
       ni = (INu+Nghost) - (INl-Nghost) + 1;
       nj = (JNu+Nghost) - (JNl-Nghost) + 1;
       nk = (KNu+Nghost) - (KNl-Nghost) + 1;
       buffer = new double[3*ni*nj*nk];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (int k  = KNl-Nghost; k <= KNu+Nghost; ++k) {
             for (int j  = JNl-Nghost; j <= JNu+Nghost; ++j) {
                for (int i = INl-Nghost; i <= INu+Nghost; ++i) {
 	           buffer[buffer_size] = Node[i][j][k].X.x;
 	           buffer[buffer_size+1] = Node[i][j][k].X.y;
 	           buffer[buffer_size+2] = Node[i][j][k].X.z;
                   buffer_size = buffer_size + 3;
                } /* endfor */
             } /* endfor */
	  } /* endfor */
       } /* endif */

       buffer_size = 3*ni*nj*nk;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (int k  = KNl-Nghost; k <= KNu+Nghost; ++k) {
             for (int j  = JNl-Nghost; j <= JNu+Nghost; ++j) {
                for (int i = INl-Nghost; i <= INu+Nghost; ++i) {
 	           Node[i][j][k].X.x = buffer[buffer_size];
 	           Node[i][j][k].X.y = buffer[buffer_size+1];
 	           Node[i][j][k].X.z = buffer[buffer_size+2];
                   buffer_size = buffer_size + 3;
                } /* endfor */
             } /* endfor */
	  } /* endfor */
       } /* endif */

       delete []buffer;
       buffer = NULL;
    } /* endif */

    /*  On non-primary MPI processors, update the cells
        of the hexahedral mesh block. */

    if (mesh_allocated && !CFFC_Primary_MPI_Processor()) {
       Update_Cells();
    } /* endif */

    /* Broadcast boundary condition type info. */

    if (mesh_allocated) {
       ni = (INu+Nghost) - (INl-Nghost) + 1;
       nj = (JNu+Nghost) - (JNl-Nghost) + 1;
       nk = (KNu+Nghost) - (KNl-Nghost) + 1;
       i_buffer = new int[2*ni*nk+2*nj*nk+2*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          
          for (int k = KCl-Nghost; k <= KCu+Nghost; ++k) {
             for (int i = ICl-Nghost; i <= ICu+Nghost; ++i) {
                i_buffer[buffer_size] = BCtypeN[i][k];
	        i_buffer[buffer_size+1] = BCtypeS[i][k];
                buffer_size = buffer_size + 2;
             } /* endfor */
          } /* endfor */

          for (int k = KCl-Nghost; k <= KCu+Nghost; ++k) {
             for (int j = JCl-Nghost; j <= JCu+Nghost; ++j) {
	        i_buffer[buffer_size] = BCtypeE[j][k];
	        i_buffer[buffer_size+1] = BCtypeW[j][k];
                buffer_size = buffer_size + 2;
             } /* endfor */
          } /* endfor */

	  for (int j = JCl-Nghost; j <= JCu+Nghost; ++j) {
             for (int i = ICl-Nghost; i <= ICu+Nghost; ++i) {
	        i_buffer[buffer_size] = BCtypeT[i][j];
	        i_buffer[buffer_size+1] = BCtypeB[i][j];
                buffer_size = buffer_size + 2;
	     } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = 2*ni*nk+2*nj*nk+2*ni*nj;
       MPI::COMM_WORLD.Bcast(i_buffer, buffer_size, MPI::INT, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;

          for (int k = KCl-Nghost; k <= KCu+Nghost; ++k) {
             for (int i = ICl-Nghost; i <= ICu+Nghost; ++i) {
                BCtypeN[i][k] = i_buffer[buffer_size];
	        BCtypeS[i][k] = i_buffer[buffer_size+1];
                buffer_size = buffer_size + 2;
             } /* endfor */
          } /* endfor */

          for (int k = KCl-Nghost; k <= KCu+Nghost; ++k) {
             for (int j = JCl-Nghost; j <= JCu+Nghost; ++j) {
	        BCtypeE[j][k] = i_buffer[buffer_size];
	        BCtypeW[j][k] = i_buffer[buffer_size+1];
                buffer_size = buffer_size + 2;
             } /* endfor */
          } /* endfor */

	  for (int j = JCl-Nghost; j <= JCu+Nghost; ++j) {
             for (int i = ICl-Nghost; i <= ICu+Nghost; ++i) {
	        BCtypeT[i][j] = i_buffer[buffer_size];
	        BCtypeB[i][j] = i_buffer[buffer_size+1];
                buffer_size = buffer_size + 2;
	     } /* endfor */
          } /* endfor */
       } /* endif */

       delete []i_buffer;
       i_buffer = NULL;
    } /* endif */

#endif

}

#ifdef _MPI_VERSION
/********************************************************
 * Routine: Broadcast                                   *
 *                                                      *
 * Broadcast the node and cell locations of hexahedral  *
 * mesh block.                                          *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Broadcast(MPI::Intracomm &Communicator) {

}
#endif

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the nodes of the hexahedral mesh to the       *
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
                << "I = " << INu - INl + 1 + 2*Nghost << " \\ \n"
                << "J = " << JNu - JNl + 1 + 2*Nghost << " \\ \n"
		<< "K = " << KNu - KNl + 1 + 2*Nghost << " \\ \n"
                << "F = POINT \n";

    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << INu - INl + 1 + 2*Nghost << " \\ \n"
                << "J = " << JNu - JNl + 1 + 2*Nghost << " \\ \n"
		<< "K = " << KNu - KNl + 1 + 2*Nghost << " \\ \n" 
                << "F = POINT \n";
    } /* endif */
    
    for(int k = KNl-Nghost ; k <= KNu+Nghost ; ++k) {
      for (int j  = JNl-Nghost ; j <= JNu+Nghost ; ++j ) {
	for (int i = INl-Nghost ; i <= INu+Nghost ; ++i ) {
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
                << "I = " << ICu - ICl + 1 + 2*Nghost << " \\ \n"
                << "J = " << JCu - JCl + 1 + 2*Nghost << " \\ \n"
		<< "K = " << KCu - KCl + 1 + 2*Nghost << " \\ \n"
                << "F = POINT \n";

    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << ICu - ICl + 1 + 2*Nghost << " \\ \n"
                << "J = " << JCu - JCl + 1 + 2*Nghost << " \\ \n"
		<< "K = " << KCu - KCl + 1 + 2*Nghost << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
      for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
	for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	  Out_File << " " << Cell[i][j][k].Xc << "\n";
	} /* endfor */
      } /* endfor */
    }/*endfor*/
    Out_File << setprecision(6);

}

/********************************************************
 * Routine: Output_Gnuplot                              *
 *                                                      *
 * Writes the nodes of the hexahedral mesh to the       *
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

/**************************************************************
 * Routine: Update_Exterior_Nodes                             *
 *                                                            *
 * Update the exterior nodes for the hexahedral mesh block.   *
 *                                                            *
 **************************************************************/
void Grid3D_Hexa_Block::Update_Exterior_Nodes(void) {

  Vector3D norm_dir, X_norm, X_tan;
  Vector3D norm_dir1,X_norm1,X_tan1;

  for (int k = KNl; k <= KNu; ++k){
    for (int j = JNl; j <= JNu; ++j){
      for (int z = 1; z <= Nghost; z++) {
	 //create ghost cells on west side
	 if (j > JNl && j < JNu && k > KNl && k < KNu) {
	    norm_dir = - QUARTER*(nfaceW(ICl, j, k) +
		  	          nfaceW(ICl, j-1, k)+
                                  nfaceW(ICl, j, k-1)+
                                  nfaceW(ICl, j-1, k-1));
	 } else if (j > JNl && j < JNu) {
	    norm_dir = - HALF*(nfaceW(ICl, j, k) +
                               nfaceW(ICl, j-1, k));
	 } else if (k > KNl && k < KNu) {
	    norm_dir = - HALF*(nfaceW(ICl, j, k) +
                               nfaceW(ICl, j, k-1));
         } else {
            norm_dir = - nfaceW(ICl, j, k);
	 } /* endif */
	 X_norm = ((Node[INl+z][j][k].X-
		    Node[INl][j][k].X)*norm_dir)*norm_dir;
	 X_tan  = (Node[INl+z][j][k].X - 
		   Node[INl][j][k].X)-X_norm;
	 Node[INl-z][j][k].X = Node[INl][j][k].X - X_norm + X_tan; 

	 //create ghost cells on east side
	 if (j > JNl && j < JNu && k > KNl && k < KNu) {
	    norm_dir1 = QUARTER*(nfaceE(ICu, j, k) +
		                 nfaceE(ICu, j-1, k)+
                                 nfaceE(ICu, j, k-1)+
                                 nfaceE(ICu, j-1, k-1));
	 } else if (j > JNl && j < JNu) {
	    norm_dir1 = HALF*(nfaceE(ICu, j, k) +
                              nfaceE(ICu, j-1, k));
	 } else if (k > KNl && k < KNu) {
	    norm_dir1 = HALF*(nfaceE(ICu, j, k) +
                              nfaceE(ICu, j, k-1));
         } else {
            norm_dir1 = nfaceE(ICu, j, k);
	 } /* endif */
	 X_norm1 = ((Node[INu][j][k].X-
		     Node[INu-z][j][k].X)*norm_dir1)*norm_dir1;
	 X_tan1  = (Node[INu][j][k].X - 
		    Node[INu-z][j][k].X)-X_norm1;
	 Node[INu+z][j][k].X = Node[INu][j][k].X + X_norm1 - X_tan1; 
      }
    }
  }
  
  for (int k = KNl; k <= KNu; ++k) {
    for (int i = INl; i <= INu; ++i) {
      for (int z = 1; z <= Nghost; z++) {
 	 //create ghost cells on south side
	 if (i > INl && i < INu && k > KNl && k < KNu) {
	    norm_dir = - QUARTER*(nfaceS(i, JCl, k) +
		  	          nfaceS(i-1, JCl, k)+
                                  nfaceS(i, JCl, k-1)+
                                  nfaceS(i-1, JCl, k-1));
	 } else if (i > INl && i < INu) {
	    norm_dir = - HALF*(nfaceS(i, JCl, k) +
		 	       nfaceS(i-1, JCl, k));
	 } else if (k > KNl && k < KNu) {
	    norm_dir = - HALF*(nfaceS(i, JCl, k) +
			       nfaceS(i, JCl, k-1));
         } else {
            norm_dir = - nfaceS(i, JCl, k);
	 } /* endif */
	 X_norm = ((Node[i][JNl+z][k].X-
		    Node[i][JNl][k].X)*norm_dir)*norm_dir;
	 X_tan  = (Node[i][JNl+z][k].X - 
	           Node[i][JNl][k].X)-X_norm;
	 Node[i][JNl-z][k].X = Node[i][JNl][k].X - X_norm + X_tan; 

   	 //create ghost cells on north side
	 if (i > INl && i < INu && k > KNl && k < KNu) {
	    norm_dir1 = QUARTER*(nfaceN(i, JCu, k) +
		  	         nfaceN(i-1, JCu, k)+
                                 nfaceN(i, JCu, k-1)+
                                 nfaceN(i-1, JCu, k-1));
	 } else if (i > INl && i < INu) {
	    norm_dir1 = HALF*(nfaceN(i, JCu, k) +
		 	      nfaceN(i-1, JCu, k));
	 } else if (k > KNl && k < KNu) {
	    norm_dir1 = HALF*(nfaceN(i, JCu, k) +
			      nfaceN(i, JCu, k-1));
         } else {
            norm_dir1 = nfaceN(i, JCu, k);
	 } /* endif */
	 X_norm1 = ((Node[i][JNu][k].X-
		     Node[i][JNu-z][k].X)*norm_dir1)*norm_dir1;
	 X_tan1  = (Node[i][JNu][k].X - 
	            Node[i][JNu-z][k].X)-X_norm1;
	 Node[i][JNu+z][k].X = Node[i][JNu][k].X + X_norm1 - X_tan1; 
      }
    }
  }

  for (int j = JNl; j <= JNu; ++j) {
    for (int i = INl; i <= INu; ++i) {
      for (int z = 1; z <= Nghost; z++) {
	 //create ghost cells on bottom side
	 if (i > INl && i < INu && j > JNl && j < JNu) {
	    norm_dir = - QUARTER*(nfaceBot(i, j, KCl) +
		  	          nfaceBot(i-1, j, KCl)+
                                  nfaceBot(i, j-1, KCl)+
                                  nfaceBot(i-1, j-1, KCl));
	 } else if (i > INl && i < INu) {
	    norm_dir = - HALF*(nfaceBot(i, j, KCl) +
		 	       nfaceBot(i-1, j, KCl));
	 } else if (j > JNl && j < JNu) {
	    norm_dir = - HALF*(nfaceBot(i, j, KCl) +
			       nfaceBot(i, j-1, KCl));
         } else {
            norm_dir = - nfaceBot(i, j, KCl);
	 } /* endif */
	 X_norm = ((Node[i][j][KNl+z].X-
		    Node[i][j][KNl].X)*norm_dir)*norm_dir;
	 X_tan  = (Node[i][j][KNl+z].X - 
		   Node[i][j][KNl].X)-X_norm;
	 Node[i][j][KNl-z].X = Node[i][j][KNl].X - X_norm + X_tan; 

	 //create ghost cells on top side
	 if (i > INl && i < INu && j > JNl && j < JNu) {
	    norm_dir1 = QUARTER*(nfaceTop(i, j, KCu) +
		  	         nfaceTop(i-1, j, KCu)+
                                 nfaceTop(i, j-1, KCu)+
                                 nfaceTop(i-1, j-1, KCu));
	 } else if (i > INl && i < INu) {
	    norm_dir1 = HALF*(nfaceTop(i, j, KCu) +
		 	      nfaceTop(i-1, j, KCu));
	 } else if (j > JNl && j < JNu) {
	    norm_dir1 = HALF*(nfaceTop(i, j, KCu) +
			      nfaceTop(i, j-1, KCu));
         } else {
            norm_dir1 = nfaceTop(i, j, KCu);
	 } /* endif */
	 X_norm1 = ((Node[i][j][KNu].X-
		     Node[i][j][KNu-z].X)*norm_dir1)*norm_dir1;
	 X_tan1  = (Node[i][j][KNu].X - 
		    Node[i][j][KNu-z].X)-X_norm1;
	 Node[i][j][KNu+z].X = Node[i][j][KNu].X + X_norm1 - X_tan1; 
      }
    }
  }
  
  //the following code deals with placing ghost cells at the edges and corners
  for (int k = KNl; k <= KNu; ++k){
    for (int j = JNu+1; j <= JNu+Nghost; ++j){
       for (int z = 1; z <= Nghost; ++z){
	  Node[INu+z][j][k] = (Node[INu][j][k].X +
		  	      (Node[INu][j][k].X - Node[INu-z][j][k].X));
	  Node[INl-z][j][k] = (Node[INl][j][k].X +
			      (Node[INl][j][k].X - Node[INl+z][j][k].X));
      }
    }
  }
  
  for (int k = KNl; k <= KNu; ++k){
    for (int j = JNl-Nghost; j <= JNl-1; ++j){
       for (int z = 1; z <= Nghost; ++z){
	  Node[INu+z][j][k] = (Node[INu][j][k].X +
			      (Node[INu][j][k].X - Node[INu-z][j][k].X));
	  Node[INl-z][j][k] = (Node[INl][j][k].X +
			      (Node[INl][j][k].X - Node[INl+z][j][k].X));
      }
    }
  }

  for (int i = INl; i <= INu; ++i){
     for (int j = JNu+1; j <= JNu+Nghost; ++j){
        for (int z = 1; z <= Nghost; ++z){
	   Node[i][j][KNu+z] = (Node[i][j][KNu].X +
			       (Node[i][j][KNu].X - Node[i][j][KNu-z].X));
	   Node[i][j][KNl-z] = (Node[i][j][KNl].X +
			       (Node[i][j][KNl].X - Node[i][j][KNl+z].X));
      }
    }
  }
  
  for (int i = INl; i <= INu; ++i){
     for (int j = JNl-Nghost; j <= JNl-1; ++j){
        for (int z = 1; z <= Nghost; ++z){
	   Node[i][j][KNu+z] = (Node[i][j][KNu].X +
			       (Node[i][j][KNu].X - Node[i][j][KNu-z].X));
	   Node[i][j][KNl-z] = (Node[i][j][KNl].X +
			       (Node[i][j][KNl].X - Node[i][j][KNl+z].X));
      }
    }
  }
  
  for (int j = JNl-Nghost; j <= JNl-1; ++j){
    for (int k = KNl-Nghost; k <= KNl-1; ++k){
      for (int z = 1; z <= Nghost; ++z){
	  Node[INu+z][j][k] = (Node[INu][j][k].X +
		 	      (Node[INu][j][k].X - Node[INu-z][j][k].X));
	  Node[INl-z][j][k] = (Node[INl][j][k].X +
		  	      (Node[INl][j][k].X - Node[INl+z][j][k].X));
      }
    }
  }
  
  for (int j = JNu+1; j <= JNu+Nghost; ++j){
    for (int k = KNl-Nghost; k <= KNl-1; ++k){
      for (int z = 1; z <= Nghost; ++z){
	  Node[INu+z][j][k] = (Node[INu][j][k].X +
		 	      (Node[INu][j][k].X - Node[INu-z][j][k].X));
	  Node[INl-z][j][k] = (Node[INl][j][k].X +
		  	      (Node[INl][j][k].X - Node[INl+z][j][k].X));
      }
    }
  }

  for (int j = JNl-Nghost; j <= JNl-1; ++j){
     for (int k = KNu+1; k <= KNu+Nghost; ++k){
        for (int z = 1; z <= Nghost; ++z){
	   Node[INu+z][j][k] = (Node[INu][j][k].X +
			       (Node[INu][j][k].X - Node[INu-z][j][k].X));
	   Node[INl-z][j][k] = (Node[INl][j][k].X +
			       (Node[INl][j][k].X - Node[INl+z][j][k].X));
      }
    }
  }

  for (int j = JNu+1; j <= JNu+Nghost; ++j){
     for (int k = KNu+1; k <= KNu+Nghost; ++k){
        for (int z = 1; z <= Nghost; ++z){
	   Node[INu+z][j][k] = (Node[INu][j][k].X +
			       (Node[INu][j][k].X - Node[INu-z][j][k].X));
	   Node[INl-z][j][k] = (Node[INl][j][k].X +
			       (Node[INl][j][k].X - Node[INl+z][j][k].X));
      }
    }
  }

//   for (int i = INu+1; i <= INu+Nghost; ++i){
//     for (int k = KNl-Nghost; k <= KNl-1; ++k){
//       for (int z = 1; z <= Nghost; ++z){
// 	  Node[i][JNu+z][k] = (Node[i][JNu][k].X +
// 			      (Node[i][JNu][k].X - Node[i][JNu-z][k].X));
// 	  Node[i][JNl-z][k] = (Node[i][JNl][k].X +
// 		              (Node[i][JNl][k].X - Node[i][JNl+z][k].X));
//       }
//     }
//   }
  
//   for (int i = INl-Nghost; i <= INl-1; ++i){
//      for (int k = KNu+1; k <= KNu+Nghost; ++k){
//        for (int z = 1; z <= Nghost; ++z){
// 	  Node[i][JNu+z][k] = (Node[i][JNu][k].X +
// 			      (Node[i][JNu][k].X - Node[i][JNu-z][k].X));
// 	  Node[i][JNl-z][k] = (Node[i][JNl][k].X +
// 			      (Node[i][JNl][k].X - Node[i][JNl+z][k].X));
//       }
//     }
//   }
  
}

/**************************************************************
 * Routine: Update_Exterior_Nodes_Zdir                        *
 *                                                            *
 * Update the exterior nodes for the hexahedral mesh block.   *
 *                                                            *
 **************************************************************/
void Grid3D_Hexa_Block::Update_Exterior_Nodes_Zdir(void) {

  Vector3D norm_dir, X_norm, X_tan;

  for (int j = JNl-Nghost; j <= JNu+Nghost; ++j) {
    for (int i = INl-Nghost; i <= INu+Nghost; ++i) {
      for (int z = 1; z <= Nghost; z++) {
	 //create ghost cells on bottom side
	 if (i == INu+Nghost && j == JNu+Nghost) {
            norm_dir = - nfaceBot(i-1, j-1, KCl);
	 } else if (i == INu+Nghost && j == JNl-Nghost) {
            norm_dir = - nfaceBot(i-1, j, KCl);
         } else if (i == INl-Nghost && j == JNu+Nghost) {
            norm_dir = - nfaceBot(i, j-1, KCl);
         } else if (i == INl-Nghost && j == JNl-Nghost) {
            norm_dir = - nfaceBot(i, j, KCl);
         } else if (i == INu+Nghost) {
	    norm_dir = - HALF*(nfaceBot(i-1, j, KCl) +
		 	       nfaceBot(i-1, j-1, KCl));
         } else if (j == JNu+Nghost) {
	    norm_dir = - HALF*(nfaceBot(i, j-1, KCl) +
		 	       nfaceBot(i-1, j-1, KCl));
         } else if (i == INl-Nghost) {
	    norm_dir = - HALF*(nfaceBot(i, j, KCl) +
		 	       nfaceBot(i, j-1, KCl));
         } else if (j == JNl-Nghost) {
	    norm_dir = - HALF*(nfaceBot(i, j, KCl) +
		 	       nfaceBot(i-1, j, KCl));
         } else {
	    norm_dir = - QUARTER*(nfaceBot(i, j, KCl) +
		  	          nfaceBot(i-1, j, KCl)+
                                  nfaceBot(i, j-1, KCl)+
                                  nfaceBot(i-1, j-1, KCl));
         } /* endif */
	 X_norm = ((Node[i][j][KNl+z].X-
		    Node[i][j][KNl].X)*norm_dir)*norm_dir;
	 X_tan  = (Node[i][j][KNl+z].X - 
		   Node[i][j][KNl].X)-X_norm;
	 Node[i][j][KNl-z].X = Node[i][j][KNl].X - X_norm + X_tan; 

	 //create ghost cells on top side
         if (i == INu+Nghost && j == JNu+Nghost) {
            norm_dir = nfaceTop(i-1, j-1, KCu);
	 } else if (i == INu+Nghost && j == JNl-Nghost) {
            norm_dir = nfaceTop(i-1, j, KCu);
         } else if (i == INl-Nghost && j == JNu+Nghost) {
            norm_dir = nfaceTop(i, j-1, KCu);
         } else if (i == INl-Nghost && j == JNl-Nghost) {
            norm_dir = nfaceTop(i, j, KCu);
         } else if (i == INu+Nghost) {
	    norm_dir = HALF*(nfaceTop(i-1, j, KCu) +
		 	     nfaceTop(i-1, j-1, KCu));
         } else if (j == JNu+Nghost) {
	    norm_dir = HALF*(nfaceTop(i, j-1, KCu) +
		 	     nfaceTop(i-1, j-1, KCu));
         } else if (i == INl-Nghost) {
	    norm_dir = HALF*(nfaceTop(i, j, KCu) +
		 	     nfaceTop(i, j-1, KCu));
         } else if (j == JNl-Nghost) {
	    norm_dir = HALF*(nfaceTop(i, j, KCu) +
		             nfaceTop(i-1, j, KCu));
         } else {
	    norm_dir = QUARTER*(nfaceTop(i, j, KCu) +
		  	        nfaceTop(i-1, j, KCu)+
                                nfaceTop(i, j-1, KCu)+
                                nfaceTop(i-1, j-1, KCu));
	 } /* endif */
	 X_norm = ((Node[i][j][KNu].X-
		    Node[i][j][KNu-z].X)*norm_dir)*norm_dir;
	 X_tan  = (Node[i][j][KNu].X - 
		   Node[i][j][KNu-z].X)-X_norm;
	 Node[i][j][KNu+z].X = Node[i][j][KNu].X + X_norm - X_tan; 
      }
    }
  }  
  
}

/**************************************************************
 * Routine: Correct_Exterior_Nodes                            *
 *                                                            *
 * Correct the exterior nodes of the hexahedral mesh block.   *
 *                                                            *
 **************************************************************/
void Grid3D_Hexa_Block::Correct_Exterior_Nodes(const int ii, 
                                               const int jj, 
                                               const int kk, 
                                               const int *be) {

   int igs, ige, igd, jgs, jge, jgd, kgs, kge, kgd;
   int im, jm, km;
   int face_not_boundary = 0;
   int i_be, j_be, k_be;
   
   im = 0;
   jm = 0;
   km = 0;
   
   if (ii) {
      if (ii>0) {
         igs = INu +1;
         ige = INu +Nghost;
         igd = 1;
      } else {
         igs = INl-1;
         ige = INl-Nghost;
         igd = -1;
      }
      i_be = 9*(ii+1) + 3 +1;
      if (be[i_be]) {
         im = -ii;
      } else {
         ++face_not_boundary;
      }
   } else {
      igs = INl;
      ige = INu;
      igd = 1;
   } /* endif */

   if (jj) {
      if (jj>0){
         jgs = JNu +1;
         jge = JNu +Nghost;
         jgd = 1;
      } else {
         jgs = JNl-1;
         jge = JNl-Nghost;
         jgd = -1;
      }
      j_be = 9 + 3*(jj+1) +1;
      if (be[j_be]) {
         jm = -jj;
      } else {
         ++face_not_boundary;
      }
   } else {
      jgs = JNl;
      jge = JNu;
      jgd = 1;
   } /* endif */
   
   if (kk) {
      if (kk>0) {
         kgs = KNu +1;
         kge = KNu +Nghost;
         kgd = 1;
      } else {
         kgs = KNl - 1;
         kge = KNl - Nghost;
         kgd = -1;
      }
      k_be = 9 + 3 +(kk+1);
      if (be[k_be]) {
         km =   -kk;
      } else {
         ++face_not_boundary;
      }
   } else {
      kgs = KNl;
      kge = KNu;
      kgd = 1;
   } /* endif */

   // A special case of those edges whose two faces are interior, but the edge is 
   // at the domain extent. 
   if (face_not_boundary == abs(ii) + abs(jj) + abs(kk)) {
      im = -ii;
      jm = -jj;
      km = -kk;
   } /* endif */
   
   for (int i = igs; (i-igs)*(i-ige)<=0; i+= igd) {
      for (int j = jgs; (j-jgs)*(j-jge)<=0; j+= jgd) {
         for (int k = kgs; (k-kgs)*(k-kge)<=0; k+= kgd){
            Node[i][j][k].X = 2*Node[i+im*(abs(i-igs)+1)][j+jm*(abs(j-jgs)+1)][k+km*(abs(k-kgs)+1)].X -
                              Node[i+2*im*((abs(i-igs))+1)][j+2*jm*((abs(j-jgs))+1)][k+2*km*((abs(k-kgs))+1)].X;
                      
         } /* endfor */
      } /* endfor */
   } /* endfor */
   
}

/**************************************************************
 * Routine: Fix_Corner_Cells_for_3_Blks_Abutting              *
 *                                                            *
 * For those three blocks abutting each other, each block     *
 * has no corner nodes. The corner nodes geometry and         *
 * solutons don't have real physical meaning. This situation  *
 * will corrupte the gradient reconstruction. Also the        *
 * output soluion will have these unphysical regions, which   *
 * might confuse the analysis.  The most convenient way       *
 * to fix those nodes are that just make them coincide with   *
 * the nearest phyiscal ones, and all the reconstructions     *
 * and outputs remain the general format.                     *
 **************************************************************/
int Grid3D_Hexa_Block::Fix_Corner_Cells_for_3_Blks_Abutting(const int i_elem, 
                                                             const int j_elem, 
                                                             const int k_elem, 
                                                            const int numNeigh,
                                                            const int be) {

   int execute_this_prog = 0;
   
   if( ((abs(i_elem) && abs(j_elem) && !(k_elem)) ||
        (abs(i_elem) && abs(k_elem) && !(j_elem)) ||
        (abs(j_elem) && abs(k_elem) && !(i_elem))) && (!numNeigh && !be) ) {
      
      execute_this_prog = 1;
   }
   
   if(!execute_this_prog) return 0;
   
   // execute this program for the true situation.
   int i_nearest, j_nearest, k_nearest;
   int i_inc, j_inc, k_inc;
   
   // default values
   i_nearest = Nghost;
   j_nearest = Nghost;
   k_nearest = Nghost;
   
   i_inc = 0;
   j_inc = 0;
   k_inc = 0;
   
   // set the corner's nearest nodes based on where the element is located.
   if(i_elem <0){
      i_nearest = Nghost;
      i_inc = -1;
   }else{
      i_nearest = ICu;
      i_inc = 1;
   }
   if(j_elem <0){
      j_nearest = Nghost;
      j_inc = -1;
   }else{
      j_nearest = JCu;
      j_inc = 1;
   }
   if(k_elem <0){
      k_nearest = Nghost;
      k_inc = -1;
   }else{
      k_nearest = KCu;
      k_inc = 1;
   }
   
   // coincide the ghost corners with the nearest physical ones.
   if(abs(i_elem) && abs(j_elem) && !(k_elem)){
      for (int kDir = KCl-Nghost; kDir<= KCu+Nghost-1; ++kDir){
         for(int jDir = 1; jDir <= Nghost; ++jDir){
            for(int iDir = 1; iDir <= Nghost; ++iDir){
               Cell[i_nearest + i_inc*iDir][j_nearest + j_inc*jDir][kDir].Xc = 
                  Cell[i_nearest][j_nearest][kDir].Xc;
               
            }
         }
      }
   }
   if(abs(i_elem) && abs(k_elem) && !(j_elem)){
      for (int jDir = JCl-Nghost; jDir<= JCu+Nghost-1; ++jDir){
         for(int iDir = 1; iDir <= Nghost; ++iDir){
            for(int kDir = 1; kDir <= Nghost; ++kDir){
               Cell[i_nearest + i_inc*iDir][jDir][k_nearest + k_inc*kDir].Xc = 
                  Cell[i_nearest][jDir][k_nearest].Xc;
            }
         }
      }
   }

   if (abs(j_elem) && abs(k_elem) && !(i_elem)){
      for (int iDir = ICl-Nghost; iDir<= ICu+Nghost-1; ++iDir){
         for(int jDir = 1; jDir <= Nghost; ++jDir){
            for(int kDir = 1; kDir <= Nghost; ++kDir){
               Cell[iDir][j_nearest + j_inc*jDir][k_nearest + k_inc*kDir].Xc = 
                  Cell[iDir][j_nearest][k_nearest].Xc;
               
            }
         }
      }
   }
   
   return 0;
   
}


/*****************************************************************
 * Routine: Update Cells                                         *
 *                                                               *
 * Updates the cell information for the hexahedral mesh block.   *
 *                                                               *
 *****************************************************************/
void Grid3D_Hexa_Block::Update_Cells(void) {

  for(int k = KCl-Nghost; k <=KCu+Nghost ; ++k){ 
    for(int j = JCl-Nghost ; j <= JCu+Nghost; ++j) {
      for (int i = ICl-Nghost ; i <= ICu+Nghost; ++i) {
	Cell[i][j][k].I = i ;
	Cell[i][j][k].J = j ;
	Cell[i][j][k].K = k ;
	Cell[i][j][k].Xc = centroid(i, j, k);
    Cell[i][j][k].V = volume(Cell[i][j][k]);

      } /* endfor */
    } /* endfor */
  } /* endfor */
    
    for(int k = KCl-Nghost; k <=KCu+Nghost ; ++k){ 
        for(int j = JCl-Nghost ; j <= JCu+Nghost; ++j) {
            for (int i = ICl-Nghost ; i <= ICu+Nghost; ++i) {
                /* calculate jacobian to 4th order */
                Cell[i][j][k].Jacobian = jacobian(i,j,k,4);
            } /* endfor */
        } /* endfor */
    } /* endfor */    
}

/*****************************************************************
 * Routine: Update Cells HighOrder                               *
 *                                                               *
 * Updates the cell information for the hexahedral mesh block.   *
 *                                                               *
 *****************************************************************/
void Grid3D_Hexa_Block::Update_Cells_HighOrder(void) {
    
    Update_Cells();
    
    // If using CENO Reconstruction, then we need 
    // to compute and store the geometric coefficients
    // -----------------------------------------------
    if(Grid3D_HO_Execution_Mode::USE_HO_CENO_GRID){

        for(int k = KCl-Nghost; k <=KCu+Nghost ; ++k){ 
            for(int j = JCl-Nghost ; j <= JCu+Nghost; ++j) {
                for (int i = ICl-Nghost ; i <= ICu+Nghost; ++i) {
                    ComputeGeometricCoefficients(i,j,k);
                } /* endfor */
            } /* endfor */
        } /* endfor */
    }
}

/*****************************************************************
 * Routine: Update Ghost Cells                                   *
 *                                                               *
 * Updates the ghost cells information of the hexahedral mesh    *
 * block.                                                        *
 *                                                               *
 *****************************************************************/
void Grid3D_Hexa_Block::Update_Ghost_Cells(void) {

  for(int k = KCl-Nghost; k <=KCu+Nghost ; ++k){ 
    for(int j = JCl-Nghost ; j <= JCu+Nghost; ++j) {
      for (int i = ICl-Nghost ; i <= ICu+Nghost; ++i) {
         if (k < KCl || k > KCu || j < JCl || j > JCu || i < ICl || i > ICu) {
            Cell[i][j][k].I = i ;
            Cell[i][j][k].J = j ;
            Cell[i][j][k].K = k ;
            Cell[i][j][k].Xc = centroid(i, j, k);
            Cell[i][j][k].V = volume(Cell[i][j][k]);
         } /* endif */

      } /* endfor */
    } /* endfor */
  } /* endfor */
    
    for(int k = KCl-Nghost; k <=KCu+Nghost ; ++k){ 
        for(int j = JCl-Nghost ; j <= JCu+Nghost; ++j) {
            for (int i = ICl-Nghost ; i <= ICu+Nghost; ++i) {
                if (k < KCl || k > KCu || j < JCl || j > JCu || i < ICl || i > ICu) {
                    /* calculate jacobian to 4th order */
                    Cell[i][j][k].Jacobian = jacobian(i,j,k,4);
                } /* endif */
            } /* endfor */
        } /* endfor */
    } /* endfor */

}


/*****************************************************************
 * Routine: Update Ghost Cells HighOrder                         *
 *                                                               *
 * Updates the ghost cells information of the hexahedral mesh    *
 * block.                                                        *
 *                                                               *
 *****************************************************************/
void Grid3D_Hexa_Block::Update_Ghost_Cells_HighOrder(void) {
    
    Update_Ghost_Cells();
    // If using CENO Reconstruction, then we need 
    // to compute and store the geometric coefficients
    // -----------------------------------------------
    if(Grid3D_HO_Execution_Mode::USE_HO_CENO_GRID){
        for(int k = KCl-Nghost; k <=KCu+Nghost ; ++k){ 
            for(int j = JCl-Nghost ; j <= JCu+Nghost; ++j) {
                for (int i = ICl-Nghost ; i <= ICu+Nghost; ++i) {
                    if (k < KCl || k > KCu || j < JCl || j > JCu || i < ICl || i > ICu) {
                        ComputeGeometricCoefficients(i,j,k);
                    } /* endif */
                } /* endfor */
            } /* endfor */
        } /* endfor */
    } /* endif */
}

/********************************************************
 * Routine: Set_BCs_Xdir                                *
 *                                                      *
 * Set boundary conditions for east & west boundaries   *
 * of hexahedral grid block.                            *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Set_BCs_Xdir(const int BCtype_east_boundary,
                                     const int BCtype_west_boundary) {

    for (int k = KCl-Nghost; k <= KCu+Nghost; ++k) {
       for (int j = JCl-Nghost; j <= JCu+Nghost; ++j) {
          BCtypeE[j][k] = BCtype_east_boundary;
          BCtypeW[j][k] = BCtype_west_boundary;
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: Set_BCs_Ydir                                *
 *                                                      *
 * Set boundary conditions for north & south boundaries *
 * of hexahedral grid block.                            *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Set_BCs_Ydir(const int BCtype_north_boundary,
                                     const int BCtype_south_boundary) {

    for (int k = KCl-Nghost; k <= KCu+Nghost; ++k) {
       for (int i = ICl-Nghost; i <= ICu+Nghost; ++i) {
 	  BCtypeN[i][k] = BCtype_north_boundary;
          BCtypeS[i][k] = BCtype_south_boundary;
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: Set_BCs_Zdir                                *
 *                                                      *
 * Set boundary conditions for top & bottom boundaries  *
 * of hexahedral grid block.                            *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Set_BCs_Zdir(const int BCtype_top_boundary,
                                     const int BCtype_bottom_boundary) {

    for (int j = JCl-Nghost; j <= JCu+Nghost; ++j) {
       for (int i = ICl-Nghost; i <= ICu+Nghost; ++i) {
 	  BCtypeT[i][j] = BCtype_top_boundary;
          BCtypeB[i][j] = BCtype_bottom_boundary;
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: Set_BCs                                     *
 *                                                      *
 * Set boundary conditions for all six boundaries       *
 * of hexahedral grid block.                            *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Set_BCs(const int BCtype_east_boundary,
                                const int BCtype_west_boundary,
                                const int BCtype_north_boundary,
                                const int BCtype_south_boundary,
                                const int BCtype_top_boundary,
                                const int BCtype_bottom_boundary) {

    Set_BCs_Xdir(BCtype_east_boundary,
                 BCtype_west_boundary);
    Set_BCs_Ydir(BCtype_north_boundary,
                 BCtype_south_boundary);
    Set_BCs_Zdir(BCtype_top_boundary,
                 BCtype_bottom_boundary);

}

/********************************************************
 * Routine: Rotate                                      *
 *                                                      *
 * Rotates the hexahedral grid block.                   *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Rotate(const double &Angle, 
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

/********************************************************
 * Routine: Extrude                                     *
 *                                                      *
 * Extrudes a 3D hexahedral grid block from 2D          *
 * quadrilateral grid block.                            *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Extrude(Grid2D_Quad_Block &Grid2D_XYplane, 
                                const int Nk,
                                const int Stretching_Flag,
                                const int i_Stretching_Kdir,
                                const double &Stretching_Kdir,
                                const double &Z_min,
                                const double &Z_max) {
    int Stretching_type_Kdir = i_Stretching_Kdir;
    if (Stretching_Flag == OFF) {
        Stretching_type_Kdir = STRETCHING_FCN_LINEAR;
    }
    
    Extrude(Grid2D_XYplane,
            Nk,
            Stretching_type_Kdir,
            Stretching_Kdir,
            Z_min,
            Z_max);
            
}

void Grid3D_Hexa_Block::Extrude(Grid2D_Quad_Block &Grid2D_XYplane, 
           			const int Nk,
                                const int i_Stretching_Kdir,
			        const double &Stretching_Kdir,
                                const double &Z_min,
                                const double &Z_max) {

  double S_k, Z;

  /* Allocate memory for the grid. */
 
  if (Allocated) deallocate();
  allocate(Grid2D_XYplane.NNi-2*Grid2D_XYplane.Nghost-1, 
           Grid2D_XYplane.NNj-2*Grid2D_XYplane.Nghost-1, 
           Nk,
           Grid2D_XYplane.Nghost);

  /* Extrude the grid points from the 2D quadrilateral grid block. */

  for (int k = KNl; k <= KNu; k++) {
    for (int j =  JNl-Nghost; j <= JNu+Nghost; j++) {
      for (int i = INl-Nghost; i <= INu+Nghost; i++) {
        S_k = double(k-KNl)/double(KNu-KNl);
        S_k = StretchingFcn(S_k, Stretching_Kdir, ZERO, i_Stretching_Kdir);
        Z = S_k * (Z_max-Z_min) + Z_min;
        Node[i][j][k].X.x = Grid2D_XYplane.Node[i][j].X.x;
        Node[i][j][k].X.y = Grid2D_XYplane.Node[i][j].X.y;
        Node[i][j][k].X.z = Z;
      } /* endfor */
    } /* endfor */
  }/* endfor */

  /* Process the extruded 3D grid. */

  //Update_Exterior_Nodes();
  Update_Exterior_Nodes_Zdir();
  Update_Cells();

  /* Assign East, West, North, and South boundary conditions. */

  for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
     for (int j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
        BCtypeW[j][k] = Grid2D_XYplane.BCtypeW[j];
        BCtypeE[j][k] = Grid2D_XYplane.BCtypeE[j];
     } /* endfor */
  } /* endfor */

  for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
     for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
        BCtypeN[i][k] = Grid2D_XYplane.BCtypeN[i];
        BCtypeS[i][k] = Grid2D_XYplane.BCtypeS[i];
     } /* endfor */ 
  } /* endfor */

}

Vector3D Grid3D_Hexa_Block::Delta_minimum(void) {
    double dx(1e10), dy(1e10), dz(1e10);
    double D;
    for (int k = KNl-Nghost; k <= KNu+Nghost-1; k++) {
        for (int j =  JNl-Nghost; j <= JNu+Nghost-1; j++) {
            for (int i = INl-Nghost; i <= INu+Nghost-1; i++) {
                D = Node[i+1][j][k].X.x-Node[i][j][k].X.x;
                if (D < dx)     dx = D;
                
                D = Node[i][j+1][k].X.y-Node[i][j][k].X.y;
                if (D < dy)     dy = D;      
                
                D = Node[i][j][k+1].X.z-Node[i][j][k].X.z;
                if (D < dz)     dz = D;
            }
        }
    }
    return Vector3D(dx,dy,dz);
}
//! Routine: ComputeGeometricCoefficients
//  -----------------------------------------------------------------------------
/*! Purpose: Computes the geometric moments with respect to the cell
 *           centroid (xCC,yCC,zCC) of cell (ii,jj,kk).
 *
 *  \note    The geometric moment is the integral over the cell domain
 *           (per unit volume) of a polynomial function with the form:
 *           (x-xCC)^l*(y-yCC)^m*(z-zCC)^n
 *
 *  \warning This subroutine is for cartesian grids (no curved boundaries).
 *
 *///----------------------------------------------------------------------------
void Grid3D_Hexa_Block::ComputeGeometricCoefficients(const int &ii, const int &jj, const int &kk){

  int p1,p2,p3;
  double DummyParam;

  // ----- Define the polynomial function for the cell ----- //

  GeneralizedPolynomialFunctionOfThreeVariables Polynom(0,0,0,
							centroid(ii,jj,kk).x,
							centroid(ii,jj,kk).y,
							centroid(ii,jj,kk).z);

  // ----- Compute The Geometric Coefficients ----- //
  
  // Note: When the method of integration is used the moments associated with the
  // ----  first order are not quite accurate

  for (int i =  Cell[ii][jj][kk].GeomCoeff().FirstElem(); i<=Cell[ii][jj][kk].GeomCoeff().LastElem(); ++i){

    p1 = Cell[ii][jj][kk].GeomCoeff(i).P1();
    p2 = Cell[ii][jj][kk].GeomCoeff(i).P2();
    p3 = Cell[ii][jj][kk].GeomCoeff(i).P3();

    Polynom.ChangePowersTo(p1,p2,p3);

    Cell[ii][jj][kk].GeomCoeffValue(i) = IntegratePolynomialOverTheCell(Polynom,14,DummyParam,ii,jj,kk) / volume(ii,jj,kk);

  } /* endfor */
}

//! Routine: IntegratePolynomialOverTheCell
//  -----------------------------------------------------------------------------
/*! Purpose: Performs the integration of a polynomial over the given
 *           cell (ii,jj,kk) using Adaptive Gaussian Quadrature
 *  
 *  \param [in]  FuncObj       the polynomial to be integrated
 *  \param [in]  digits        the number of digits required for convergence
 *  \param [in]  ii            the i-index of the cell
 *  \param [in]  jj            the j-index of the cell
 *  \param [in]  kk            the k-index of the cell
 *  \param [in] _dummy_param   used only to provide the return type
 *  
 *  \return ReturnType  the result of the integration
 *  
 *///----------------------------------------------------------------------------
template<typename FO, class ReturnType>
ReturnType Grid3D_Hexa_Block::IntegratePolynomialOverTheCell(FO FuncObj,
							     const int & digits,
							     const ReturnType & _dummy_param, 
							     const int &ii, 
							     const int &jj, 
							     const int &kk){

  return AdaptiveGaussianQuadrature(FuncObj, 
				    nodeSWBot(ii,jj,kk).X.x, nodeSEBot(ii,jj,kk).X.x, 
				    nodeSWBot(ii,jj,kk).X.y, nodeNWBot(ii,jj,kk).X.y,
				    nodeSWBot(ii,jj,kk).X.z, nodeSWTop(ii,jj,kk).X.z, 
				    digits,
				    _dummy_param);
}

#define DXDI 1
#define DYDI 2
#define DZDI 3
#define DXDJ 4
#define DYDJ 5
#define DZDJ 6
#define DXDK 7
#define DYDK 8
#define DZDK 9

/* Finite difference for equally spaced samples */
double Grid3D_Hexa_Block::Central_Finite_Difference(const int i, const int j, const int k, const int derivative, const double &dt, int order) {
    
    int n=int(ceil(order/2.0));
    int N=2*n;
    RowVector coefficients(N);
    ColumnVector samples(N);
    
    for(int p=1; p<=n; p++) {
        switch (derivative) {
            case DXDI:
                samples(n-p)   = Cell[i-p][j][k].Xc.x;
                samples(n+p-1) = Cell[i+p][j][k].Xc.x;
                break;
            case DYDI:
                samples(n-p)   = Cell[i-p][j][k].Xc.y;
                samples(n+p-1) = Cell[i+p][j][k].Xc.y;
                break;
            case DZDI:
                samples(n-p)   = Cell[i-p][j][k].Xc.z;
                samples(n+p-1) = Cell[i+p][j][k].Xc.z;
                break;
            case DXDJ:
                samples(n-p)   = Cell[i][j-p][k].Xc.x;
                samples(n+p-1) = Cell[i][j+p][k].Xc.x;
                break;
            case DYDJ:
                samples(n-p)   = Cell[i][j-p][k].Xc.y;
                samples(n+p-1) = Cell[i][j+p][k].Xc.y;
                break;
            case DZDJ:
                samples(n-p)   = Cell[i][j-p][k].Xc.z;
                samples(n+p-1) = Cell[i][j+p][k].Xc.z;
                break;
            case DXDK:
                samples(n-p)   = Cell[i][j][k-p].Xc.x;
                samples(n+p-1) = Cell[i][j][k+p].Xc.x;
                break;
            case DYDK:
                samples(n-p)   = Cell[i][j][k-p].Xc.y;
                samples(n+p-1) = Cell[i][j][k+p].Xc.y;
                break;
            case DZDK:
                samples(n-p)   = Cell[i][j][k-p].Xc.z;
                samples(n+p-1) = Cell[i][j][k+p].Xc.z;
                break;
        }
    }
    
    switch (N) {
        case 2:
            /* 2nd order */
            coefficients(0) = -1.0;
            coefficients(1) = 1.0;
            coefficients /= 2.0;
            break;
        case 4:
            /* 4th order */
            coefficients(0) =  1.0;
            coefficients(1) = -8.0;
            coefficients(2) =  8.0;
            coefficients(3) = -1.0;
            coefficients /= 12.0;
            break;
        case 6:
            /* 6th order */
            coefficients(0) = -1.0;
            coefficients(1) =  9.0;
            coefficients(2) = -45.0;
            coefficients(3) =  45.0;
            coefficients(4) = -9.0;
            coefficients(5) =  1.0;
            coefficients /= 60.0;
            break;
        case 8:
            /* 8th order */
            coefficients(0) =  3.0;
            coefficients(1) = -32.0;
            coefficients(2) =  168.0;
            coefficients(3) = -672.0;
            coefficients(4) =  672.0;
            coefficients(5) = -168.0;
            coefficients(6) =  32.0;
            coefficients(7) = -3.0;
            coefficients /= 840.0;
            break;
        case 10:
            /* 10th order */
            coefficients(0) = -2.0;
            coefficients(1) =  25.0;
            coefficients(2) = -150.0;
            coefficients(3) =  600.0;
            coefficients(4) = -2100.0;
            coefficients(5) =  2100.0;
            coefficients(6) = -600.0;
            coefficients(7) =  150.0;
            coefficients(8) = -25.0;
            coefficients(9) =  2.0;
            coefficients /= 2520.0;
            break;
    }
    
    return coefficients*samples/dt;
}

/* Finite difference for equally spaced samples */
double Grid3D_Hexa_Block::Forward_Finite_Difference(const int i, const int j, const int k, const int derivative, const double &dt, int order) {
    
    
    int n = int(ceil(order/2.0));
    order = 2*n;
    int N=2*n+1;
    RowVector coefficients(N);
    ColumnVector samples(N);
    //double sample;
    for(int p=1; p<=N; p++) {
        switch (derivative) {
            case DXDI:
                samples(N-p) = Cell[i-p][j][k].Xc.x;
                //sample = samples(N-p);
                break;
            case DYDI:
                samples(N-p) = Cell[i-p][j][k].Xc.y;
                //sample = samples(N-p);
                break;
            case DZDI:       
                samples(N-p) = Cell[i-p][j][k].Xc.z;
                //sample = samples(N-p);
                break;
            case DXDJ:
                samples(N-p) = Cell[i][j-p][k].Xc.x;
                //sample = samples(N-p);
                break;
            case DYDJ:
                samples(N-p) = Cell[i][j-p][k].Xc.y;
                //sample = samples(N-p);
                break;
            case DZDJ:
                samples(N-p) = Cell[i][j-p][k].Xc.z;
                //sample = samples(N-p);
                break;
            case DXDK:
                samples(N-p) = Cell[i][j][k-p].Xc.x;
                //sample = samples(N-p);
                break;
            case DYDK:
                samples(N-p) = Cell[i][j][k-p].Xc.y;
                //sample = samples(N-p);
                break;
            case DZDK:
                samples(N-p) = Cell[i][j][k-p].Xc.z;
                //sample = samples(N-p);
                break;
        }
    }
    
    switch (order) {
        case 2:
            /* 2nd order */
            coefficients(0) =  1.0;
            coefficients(1) = -4.0;
            coefficients(2) =  3.0;
            coefficients /= 2.0;
            break;
        case 4:
            /* 4th order */
            coefficients(0) =  3.0;
            coefficients(1) = -16.0;
            coefficients(2) =  36.0;
            coefficients(3) = -48.0;
            coefficients(4) =  25.0;
            coefficients /= 12.0;
            break;
        case 6:
            /* 6th order */
            coefficients(0) =  10.0;
            coefficients(1) = -72.0;
            coefficients(2) =  225.0;
            coefficients(3) = -400.0;
            coefficients(4) =  450.0;
            coefficients(5) = -360.0;
            coefficients(6) =  147.0;
            coefficients /= 60.0;
            break;
        case 8:
            /* 8th order */
            coefficients(0) =  105.0;
            coefficients(1) = -960.0;
            coefficients(2) =  3920.0;
            coefficients(3) = -9408.0;
            coefficients(4) =  14700.0;
            coefficients(5) = -15680.0;
            coefficients(6) =  11760.0;
            coefficients(7) = -6720.0;
            coefficients(8) =  2283.0;
            coefficients /= 840.0;
            break;
        case 10:
            /* 10th order */
            coefficients(0)  =  252.0;
            coefficients(1)  = -2800.0;
            coefficients(2)  =  14175.0;
            coefficients(3)  = -43200.0;
            coefficients(4)  =  88200.0;
            coefficients(5)  = -127008.0;
            coefficients(6)  =  132300.0;
            coefficients(7)  = -100800.0;
            coefficients(8)  =  56700.0;
            coefficients(9)  = -25200.0;
            coefficients(10) =  7381.0;
            coefficients /= 2520.0;
            break;
    }
    
    return coefficients*samples/dt;
}


/* Finite difference for equally spaced samples */
double Grid3D_Hexa_Block::Backward_Finite_Difference(const int i, const int j, const int k, const int derivative, const double &dt, int order) {
    
    int n=int(ceil(order/2.0));
    order=2*n;
    int N=2*n+1;
    RowVector coefficients(N);
    ColumnVector samples(N);
    //double sample ;
    for(int p=0; p<N; p++) {
        switch (derivative) {
            case DXDI:
                samples(p) = Cell[i+p][j][k].Xc.x;
                //cout << "Cell("<<i+p<<","<<j<<","<<k<<") = " << Cell[i+p][j][k].Xc << endl;
                //sample = samples(p);
                break;
            case DYDI:
                samples(p) = Cell[i+p][j][k].Xc.y;
                //sample = samples(p);
                break;
            case DZDI:       
                samples(p) = Cell[i+p][j][k].Xc.z;
                //sample = samples(p);
                break;
            case DXDJ:
                samples(p) = Cell[i][j+p][k].Xc.x;
                //cout << "Cell("<<i<<","<<j<<","<<k<<") = " << Cell[i+p][j][k].Xc << endl;
                //sample = samples(p);
                break;
            case DYDJ:
                samples(p) = Cell[i][j+p][k].Xc.y;
                //sample = samples(p);
                break;
            case DZDJ:
                samples(p) = Cell[i][j+p][k].Xc.z;
                //sample = samples(p);
                break;
            case DXDK:
                samples(p) = Cell[i][j][k+p].Xc.x;
                //sample = samples(p);
                break;
            case DYDK:
                samples(p) = Cell[i][j][k+p].Xc.y;
                //sample = samples(p);
                break;
            case DZDK:
                samples(p) = Cell[i][j][k+p].Xc.z;
                //sample = samples(p);
                break;
        }
    }
    
    switch (order) {
        case 2:
            /* 2nd order */
            coefficients(0) = -3.0;
            coefficients(1) =  4.0;
            coefficients(2) = -1.0;
            coefficients /= 2.0;
            break;
        case 4:
            /* 4th order */
            coefficients(0) = -25.0;
            coefficients(1) =  48.0;
            coefficients(2) = -36.0;
            coefficients(3) =  16.0;
            coefficients(4) = -3.0;
            coefficients /= 12.0;
            break;
        case 6:
            /* 6th order */
            coefficients(0) = -147.0;
            coefficients(1) =  360.0;
            coefficients(2) = -450.0;
            coefficients(3) =  400.0;
            coefficients(4) = -225.0;
            coefficients(5) =  72.0;
            coefficients(6) = -10.0;
            coefficients /= 60.0;
            break;
        case 8:
            /* 8th order */
            coefficients(0) =  2283.0;
            coefficients(1) = -6720.0;
            coefficients(2) =  11760.0;
            coefficients(3) = -15680.0;
            coefficients(4) =  14700.0;
            coefficients(5) = -9408.0;
            coefficients(6) =  3920.0;
            coefficients(7) = -960.0;
            coefficients(8) =  105.0;
            coefficients /= 840.0;
            break;
        case 10:
            /* 10th order */
            coefficients(0)  = -7381.0;
            coefficients(1)  =  25200.0;
            coefficients(2)  = -56700.0;
            coefficients(3)  =  100800.0;
            coefficients(4)  = -132300.0;
            coefficients(5)  =  127008.0;
            coefficients(6)  = -88200.0;
            coefficients(7)  =  43200.0;
            coefficients(8)  = -14175.0;
            coefficients(9)  =  2800.0;
            coefficients(10) = -252.0;
            coefficients /= 2520.0;
            break;
    }
    
    return coefficients*samples/dt;
}

double Grid3D_Hexa_Block::Finite_Difference(const int i, const int j, const int k, const int derivative, const double &dt, int order) {
    
    int n=int(ceil(order/2.0));
    
    int index, last_index;
    switch(derivative) {
        case DXDI:
        case DYDI:
        case DZDI:
            index = i;
            last_index = NCi-1;
            break;
        case DXDJ:
        case DYDJ:
        case DZDJ:
            index = j;
            last_index = NCj-1;
            break;
        case DXDK:
        case DYDK:
        case DZDK:
            index = k;
            last_index = NCk-1;
            break;
    }
    
    if (n > index) {
        return Backward_Finite_Difference(i,j,k, derivative, dt, order);
    } else if (n > last_index - index) {
        return  Forward_Finite_Difference(i,j,k, derivative, dt, order);
    } else {
        return  Central_Finite_Difference(i,j,k, derivative, dt, order);
    }
}


double Grid3D_Hexa_Block::jacobian(const Cell3D &theCell, int order) {
    return jacobian(theCell.I,theCell.J,theCell.K,order);
}

double Grid3D_Hexa_Block::jacobian(const int i, const int j, const int k, int order) {
    double dt = 1.0;
    
    double dxdi, dxdj, dxdk, dydi, dydj, dydk, dzdi, dzdj, dzdk, detJm1, detJ;
    
    /* calculate  dxdi */
    dxdi = Finite_Difference(i,j,k,DXDI, dt, order);
    dxdj = Finite_Difference(i,j,k,DXDJ, dt, order);
    dxdk = Finite_Difference(i,j,k,DXDK, dt, order);
    dydi = Finite_Difference(i,j,k,DYDI, dt, order);
    dydj = Finite_Difference(i,j,k,DYDJ, dt, order);
    dydk = Finite_Difference(i,j,k,DYDK, dt, order);
    dzdi = Finite_Difference(i,j,k,DZDI, dt, order);
    dzdj = Finite_Difference(i,j,k,DZDJ, dt, order);
    dzdk = Finite_Difference(i,j,k,DZDK, dt, order);
    
    detJm1 = -dxdk*dydj*dzdi + dxdj*dydk*dzdi + dxdk*dydi*dzdj 
        - dxdi*dydk*dzdj - dxdj*dydi*dzdk + dxdi*dydj*dzdk;
    detJ = ONE/detJm1;
      
    Cell[i][j][k].dXc.x = dxdi + dxdj + dxdk;
    Cell[i][j][k].dXc.y = dydi + dydj + dydk;
    Cell[i][j][k].dXc.z = dzdi + dzdj + dzdk;
    
    return detJ;
}

void Grid3D_Hexa_Block::Jacobian_Matrix(DenseMatrix &J, const int i, const int j, const int k, int order) {
    double dt = 1.0;
        
    DenseMatrix Jm1(3,3);
    
    double dxdi, dxdj, dxdk, dydi, dydj, dydk, dzdi, dzdj, dzdk, detJm1, detJ;
    
    /* calculate  dxdi */
    dxdi = Finite_Difference(i,j,k,DXDI, dt, order);
    dxdj = Finite_Difference(i,j,k,DXDJ, dt, order);
    dxdk = Finite_Difference(i,j,k,DXDK, dt, order);
    dydi = Finite_Difference(i,j,k,DYDI, dt, order);
    dydj = Finite_Difference(i,j,k,DYDJ, dt, order);
    dydk = Finite_Difference(i,j,k,DYDK, dt, order);
    dzdi = Finite_Difference(i,j,k,DZDI, dt, order);
    dzdj = Finite_Difference(i,j,k,DZDJ, dt, order);
    dzdk = Finite_Difference(i,j,k,DZDK, dt, order);
    
    Jm1(0,0) = dxdi;    Jm1(0,1) = dxdj;    Jm1(0,2) = dxdk;
    Jm1(1,0) = dydi;    Jm1(1,1) = dydj;    Jm1(1,2) = dydk;
    Jm1(2,0) = dzdi;    Jm1(2,1) = dzdj;    Jm1(2,2) = dzdk;
    
    detJm1 = -dxdk*dydj*dzdi + dxdj*dydk*dzdi + dxdk*dydi*dzdj 
    - dxdi*dydk*dzdj - dxdj*dydi*dzdk + dxdi*dydj*dzdk;
    detJ = ONE/detJm1;
    
    /* calculate inverse */
    
    J(0,0) = detJ*(-(dydk*dzdj) + dydj*dzdk);
    J(0,1) = detJ*(dydk*dzdi - dydi*dzdk);
    J(0,2) = detJ*(-(dydj*dzdi) + dydi*dzdj);
    J(1,0) = detJ*(dxdk*dzdj - dxdj*dzdk);
    J(1,1) = detJ*(-(dxdk*dzdi) + dxdi*dzdk);
    J(1,2) = detJ*(dxdj*dzdi - dxdi*dzdj);
    J(2,0) = detJ*(-(dxdk*dydj) + dxdj*dydk);
    J(2,1) = detJ*(dxdk*dydi - dxdi*dydk);
    J(2,2) = detJ*(-(dxdj*dydi) + dxdi*dydj);
    
}


void Grid3D_Hexa_Block::Disturb_Interior_Nodes(const int Number_of_Iterations) {
    double MinDistance, phi, theta , Displacement;
	
	/* Displace the interior nodes of the quadrilateral mesh block without affecting the boundary nodes */
	for (int num_iter=1; num_iter<=Number_of_Iterations; ++num_iter){
		
		//srand48(num_iter);        // set the seed for the pseudo-random number generator
		
		for (int i=INl+1; i<=INu-1 ; i++) {
			for (int j=JNl+1; j<=JNu-1; j++) {
                for (int k=KNl+1; k<=KNu-1; k++) {
                    
                    // Determine the minimum distance between the Node[i][j] and all the neighbour edges
                    MinDistance = MinimumNodeFaceDistance(i,j,k);
                    
                    // Generate a random angle for the displacement direction
                    phi   = 2.0*PI*drand48();
                    theta =     PI*drand48();
                    
                    // Calculate the displacement -> 5% of the MinDistance
                    Displacement = 0.05 * MinDistance;
                    
                    // Calculate the new node location
                    Node[i][j][k].X.x += Displacement*sin(theta)*cos(phi);
                    Node[i][j][k].X.y += Displacement*sin(theta)*sin(phi);
                    Node[i][j][k].X.z += Displacement*cos(theta);
                    
                    
                    
                }    // endfor (i)
            } // endfor (j)
        } // endfor (k)
        
    } //endfor (num_iter)
}

/********************************************************\
 * Routine: MinimumNodeFaceDistance                       *
 *                                                        *
 * Determines the minimum distance between the Node[i][j] *
 * and all the neighbour edges, for a quadrilateral mesh. *
 \********************************************************/
double Grid3D_Hexa_Block::MinimumNodeFaceDistance(const int i, const int j, const int k) 
{
	
	// Obs. For the quadrilateral mesh there are 8 edges and 4 diagonals that must be checked for determining the min distance.
	/***************************************\
     *          6      5                     *
     *       o-----o------o                  *
     *       |            |                  *
     *     7 |   (i,j)    | 4                *
     *       o     o      o                  *
     *       |            |                  *
     *     8 |            | 3                *
     *       o-----o------o                  *
     *          1     2                      *
     \***************************************/
	
	double MinDistance;
	
    // Face (i-1/2,j,k)
	
	MinDistance = DistanceFromPointToFace(Node[i][j][k], Node[i-1][j-1][k-1], Node[i-1][j-1][k], Node[i-1][j][k], Node[i-1][j][k-1]);
	MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i-1][j-1][k], Node[i-1][j-1][k+1], Node[i-1][j][k+1], Node[i-1][j][k]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i-1][j][k], Node[i-1][j][k+1], Node[i-1][j+1][k+1], Node[i-1][j+1][k]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i-1][j][k-1], Node[i-1][j][k], Node[i-1][j+1][k], Node[i-1][j+1][k-1]));
    
    // Face (i+1/2,j,k)
	
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i+1][j-1][k-1], Node[i+1][j-1][k], Node[i+1][j][k], Node[i+1][j][k-1]));
	MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i+1][j-1][k], Node[i+1][j-1][k+1], Node[i+1][j][k+1], Node[i+1][j][k]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i+1][j][k], Node[i+1][j][k+1], Node[i+1][j+1][k+1], Node[i+1][j+1][k]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i+1][j][k-1], Node[i+1][j][k], Node[i+1][j+1][k], Node[i+1][j+1][k-1]));
    
    
    // Face (i,j-1/2,k)
	
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i-1][j-1][k-1], Node[i-1][j-1][k], Node[i][j-1][k], Node[i][j-1][k-1]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i-1][j-1][k], Node[i-1][j-1][k+1], Node[i][j-1][k+1], Node[i][j-1][k]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i][j-1][k], Node[i][j-1][k+1], Node[i+1][j-1][k+1], Node[i+1][j-1][k]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i][j-1][k-1], Node[i][j-1][k], Node[i+1][j-1][k], Node[i+1][j-1][k-1]));
    
    // Face (i,j+1/2,k)
	
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i-1][j+1][k-1], Node[i-1][j+1][k], Node[i][j+1][k], Node[i][j+1][k-1]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i-1][j+1][k], Node[i-1][j+1][k+1], Node[i][j+1][k+1], Node[i][j+1][k]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i][j+1][k], Node[i][j+1][k+1], Node[i+1][j+1][k+1], Node[i+1][j+1][k]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i][j+1][k-1], Node[i][j+1][k], Node[i+1][j+1][k], Node[i+1][j+1][k-1]));
    
    
    
    // Face (i,j,k-1/2)
	
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i-1][j-1][k-1], Node[i-1][j][k-1], Node[i][j][k-1], Node[i][j-1][k-1]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i-1][j][k-1], Node[i-1][j+1][k-1], Node[i][j+1][k-1], Node[i][j][k-1]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i][j][k-1], Node[i][j+1][k-1], Node[i+1][j+1][k-1], Node[i+1][j][k-1]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i][j-1][k-1], Node[i][j][k-1], Node[i+1][j][k-1], Node[i+1][j-1][k-1]));
    
    // Face (i,j,k+1/2)
	
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i-1][j-1][k+1], Node[i-1][j][k+1], Node[i][j][k+1], Node[i][j-1][k+1]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i-1][j][k+1], Node[i-1][j+1][k+1], Node[i][j+1][k+1], Node[i][j][k+1]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i][j][k+1], Node[i][j+1][k+1], Node[i+1][j+1][k+1], Node[i+1][j][k+1]));
    MinDistance = min(MinDistance, 
                      DistanceFromPointToFace(Node[i][j][k], Node[i][j-1][k+1], Node[i][j][k+1], Node[i+1][j][k+1], Node[i+1][j-1][k+1]));
    
	return MinDistance;
}

/**********************************************************
 * Routine: DistanceFromPointToLine                       *
 *                                                        *
 * Determines the distance between the Point and the face *
 *  by using 4 points in the face and the given point     *
 **********************************************************/
double Grid3D_Hexa_Block::DistanceFromPointToFace(const Node3D &Point, const Node3D &node1, const Node3D &node2, const Node3D &node3, const Node3D &node4) {
    Vector3D Xp = (node1.X + node2.X + node3.X + node4.X)/FOUR;
    Vector3D n1,n2,n3,n4;
    
	/***************************************\
     *     2                3              *
     *       o------------o                *
     *       |            |                *
     *       |     Xp     |                *
     *       |     o      |                *
     *       |            |                *
     *       |            |                *
     *       o------------o                *
     *     1                4              *
     \***************************************/
    
    n1 = (Xp-node1.X)^(node2.X-node1.X);
    n2 = (Xp-node2.X)^(node3.X-node2.X);
    n3 = (Xp-node3.X)^(node4.X-node3.X);
    n4 = (Xp-node4.X)^(node1.X-node4.X);
    if(abs(n1)!=0)
        n1=n1/abs(n1);
    if(abs(n2)!=0)
        n2=n2/abs(n2);
    if(abs(n3)!=0)
        n3=n3/abs(n3);
    if(abs(n4)!=0)
        n4=n4/abs(n4);
    
    double A1, A2, A3, A4;
    A1 = abs(((Xp-node1.X)^(node1.X-node2.X)))/2;
    A2 = abs(((Xp-node2.X)^(node2.X-node3.X)))/2;
    A3 = abs(((Xp-node3.X)^(node3.X-node4.X)))/2;
    A4 = abs(((Xp-node4.X)^(node4.X-node1.X)))/2;
    
    
    Vector3D n =((A1*n1+A2*n2+A3*n3+A4*n4));
    if(abs(n)!=0)
        n=n/abs(n); 
    
    return DistanceFromPointToFace(Point, Xp, n);
}

/**********************************************************
 * Routine: DistanceFromPointToLine                       *
 *                                                        *
 * Determines the distance between the Point and the face *
 * defined by a point on the face and the normal.         *
 *********************************************************/
double Grid3D_Hexa_Block::DistanceFromPointToFace(const Node3D &Point, const Vector3D &Xp, const Vector3D &n) {
    Vector3D dX = Point.X - Xp;
    return fabs(n*dX)/n.abs();
}
