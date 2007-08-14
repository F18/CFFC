/* Grid3DHexaBlock.cc:  Single-block subroutines for 
                             3D hexahedral block grid class. */

/* Include 3D hexahedral block grid type header file. */

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "Grid3DHexaBlock.h"
#endif // _GRID3D_HEXA_BLOCK_INCLUDED

/*************************************************************************
 * Grid3D_Hexa_Block -- Member functions for single grid block.          *
 *************************************************************************/

/********************************************************
 * Routine: Copy                                        *
 *                                                      *
 * Copies the node and cell locations of hexahedral     *
 * mesh block Grid2 to current grid block.              *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Copy(Grid3D_Hexa_Block &Grid2) {

    int i, j,k, ni, nj, nk;
 
    if (Grid2.Used) {

       /* Allocate (re-allocate) memory for the cells and nodes 
          of the hexahedral mesh block as necessary. */
      
       if (NNi != Grid2.NNi || NNj != Grid2.NNj ||
           NNk != Grid2.NNk || NCk != Grid2.NCk ||
           NCi != Grid2.NCi || NCj != Grid2.NCj ) { 
          if (Used) { 
             deallocate();
          } /* endif */

          allocate(Grid2.NCi-2*Grid2.Nghost, 
                   Grid2.NCj-2*Grid2.Nghost, 
                   Grid2.NCk-2*Grid2.Nghost,
                   Grid2.Nghost);
       } /* endif */

       /* Copy the node locations of grid block Grid2. */

       for (k = Grid2.KNl-Grid2.Nghost; k <= Grid2.KNu+Grid2.Nghost; ++k) {
	  for (j = Grid2.JNl-Grid2.Nghost; j <= Grid2.JNu+Grid2.Nghost; ++j) {
	     for (i = Grid2.INl-Grid2.Nghost; i <= Grid2.INu+Grid2.Nghost; ++i) {
	        Node[i][j][k].X = Grid2.Node[i][j][k].X;
	     } /* endfor */
	  } /* endfor */
       } /* endfor */
    
       /* Copy the cell values of grid block Grid2. */

       for (k = Grid2.KCl-Grid2.Nghost; k <= Grid2.KCu+Grid2.Nghost; ++k) {
          for (j = Grid2.JCl-Grid2.Nghost; j <= Grid2.JCu+Grid2.Nghost; ++j) {
	     for (i = Grid2.ICl-Grid2.Nghost; i <= Grid2.ICu+Grid2.Nghost; ++i) {
	        Cell[i][j][k].I  = Grid2.Cell[i][j][k].I;
	        Cell[i][j][k].J  = Grid2.Cell[i][j][k].J;
	        Cell[i][j][k].Xc = Grid2.Cell[i][j][k].Xc;
	        Cell[i][j][k].V  = Grid2.Cell[i][j][k].V;
	     } /* endfor */
          } /* endfor */
       } /* endfor */

       /* Copy the boundary condition type info of grid block Grid2. */

       for (i = Grid2.ICl-Grid2.Nghost; i <= Grid2.ICu+Grid2.Nghost; ++i) {
	  for (k = Grid2.KCl-Grid2.Nghost; k <= Grid2.KCu+Grid2.Nghost; ++k) {
	     BCtypeN[i][k] = Grid2.BCtypeN[i][k];
	     BCtypeS[i][k] = Grid2.BCtypeS[i][k];
          } /* endfor */
       } /* endfor */

       for (j = Grid2.JCl-Grid2.Nghost; j <= Grid2.JCu+Grid2.Nghost; ++j) {
          for (k = Grid2.KCl-Grid2.Nghost; k <= Grid2.KCu+Grid2.Nghost; ++k) {
	      BCtypeE[j][k] = Grid2.BCtypeE[j][k];
	      BCtypeW[j][k] = Grid2.BCtypeW[j][k];
          } /* endfor */
       } /* endfor */

       for (i = Grid2.ICl-Grid2.Nghost; i <= Grid2.ICu+Grid2.Nghost; ++i) {
	   for (j = Grid2.JCl-Grid2.Nghost; j <= Grid2.JCu+Grid2.Nghost; ++j) {
	     BCtypeT[i][j] = Grid2.BCtypeT[i][j];
	     BCtypeB[i][j] = Grid2.BCtypeB[i][j];
	   } /* endfor */
       } /* endfor */

    } /* endif */

}

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
void Grid3D_Hexa_Block::Update_Exterior_Nodes(){

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
void Grid3D_Hexa_Block::Update_Exterior_Nodes_Zdir(){

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

/*****************************************************************
 * Routine: Update Cells                                         *
 *                                                               *
 * Updates the cell information for the hexahedral mesh block.   *
 *                                                               *
 *****************************************************************/
void Grid3D_Hexa_Block::Update_Cells(void) {

  for(int k = KCl-Nghost; k <=KCu+Nghost ; ++k){ 
    for( int j = JCl-Nghost ; j <= JCu+Nghost; ++j) {
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
                                const int i_Stretching_Kdir,
			        const double &Stretching_Kdir,
                                const double &Z_min,
                                const double &Z_max) {

  double S_k, Z;

  /* Allocate memory for the grid. */
 
  if (Used) deallocate();
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

/********************************************************
 * Routine: Create_Grid_Cube                            *
 *                                                      *
 * Creates a hexahedral grid block for a Cartesian      *
 * mesh defined on a cube with given length, width,     *
 * and height (origion 0,0,0)                           *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Block::Create_Grid_Cube(const double &Length,
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
   for ( k = KNl-Nghost ; k <= KNu+Nghost ; ++k) 
      for ( j = JNl-Nghost ; j <= JNu+Nghost ; ++j) 
         for ( i = INl-Nghost ; i <= INu+Nghost ; ++i) {
            
            // stretching to both x ends
            if(alpha > 0.0){
               // stretching to both walls in z direction. overwrite  Node[i][j][k].X.z.
               xx = double(i-INl)/double(INu-INl);
               Node[i][j][k].X.x = x_orig + Length*
                                   StretchingFcn(xx, alpha, ZERO, STRETCHING_FCN_MINMAX_CLUSTERING);
	    }else{
               Node[i][j][k].X.x = x_orig + (i-Nghost)*dx;
                            
	    }

            // stretching to both y ends
            if(beta > 0.0){
               
               // stretching to both walls in z direction. overwrite  Node[i][j][k].X.z.
               xx = double(j-JNl)/double(JNu-JNl);
               Node[i][j][k].X.y = y_orig + Width*
                                   StretchingFcn(xx, beta, ZERO, STRETCHING_FCN_MINMAX_CLUSTERING);
               
            }else{
               Node[i][j][k].X.y = y_orig + (j-Nghost)*dy;
	    }

            //stretching to both z ends
            if(gamma > 0.0){
               // stretching to both walls in z direction. overwrite  Node[i][j][k].X.z.
               xx = double(k-KNl)/double(KNu-KNl);
               Node[i][j][k].X.z = z_orig + Height*
                                   StretchingFcn(xx, gamma, ZERO, STRETCHING_FCN_MINMAX_CLUSTERING);
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
   for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
      for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
         BCtypeT[i][j] =  BCtype_top;
         BCtypeB[i][j] =  BCtype_bottom;
      }
   }

   for ( k = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
      for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j){
         BCtypeW[j][k] =  BCtype_west;
         BCtypeE[j][k] =  BCtype_east;
      }
   }

   for ( k = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
      for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i){
         BCtypeN[i][k] =  BCtype_north;
         BCtypeS[i][k] =  BCtype_south;
      }
   }

}
