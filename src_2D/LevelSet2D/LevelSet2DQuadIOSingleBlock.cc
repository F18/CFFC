/******************************************************************//**
 * \file LevelSet2DQuadIOSingleBlock.cc                               
 *                                                                    
 * Single-block versions of input and output subroutines for 2D Level 
 * Set multi-block quadrilateral mesh solution class.                 
 *                                                                    
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Input and Output Single Block External    *
 *                          Subroutines.                              *
 **********************************************************************/

/******************************************************************//**
 * Routine: Write_Solution_Block                                      
 *                                                                    
 * Writes the cell centred solution values of the specified           
 * quadrilateral solution block to the specified output stream for    
 * restart purposes.                                                  
 *                                                                    
 **********************************************************************/
void Write_Solution_Block(LevelSet2D_Quad_Block &SolnBlk,
	                  ostream &Out_File) {
  
  Out_File << setprecision(14) << SolnBlk << setprecision(6);
  
}

/******************************************************************//**
 * Routine: Read_Solution_Block                                       
 *                                                                    
 * Reads the cell centred solution values for the specified           
 * quadrilateral solution block from the specified input stream as    
 * required for restart purposes.                                     
 *                                                                    
 **********************************************************************/
void Read_Solution_Block(LevelSet2D_Quad_Block &SolnBlk,
	                 istream &In_File) {
  
  In_File >> SolnBlk;
  
}

/******************************************************************//**
 * Routine: Output_Tecplot                                            
 *                                                                    
 * Writes the solution values at the nodes of the specified           
 * quadrilateral solution block to the specified output stream        
 * suitable for plotting with TECPLOT.                                
 *                                                                    
 **********************************************************************/
void Output_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
                    const int Number_of_Time_Steps,
                    const double &Time,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {

  // Ensure boundary conditions are updated before evaluating
  // solution at the nodes.
  BCs(SolnBlk);

  // Output node solution data.
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D LevelSet Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"psi\" \\ \n"
	     << "\"F\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n";
  }
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	   << "F = POINT \n";
  for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
    for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
      Out_File << SolnBlk.Grid.Node[i][j].X
	       << SolnBlk.Un(i,j);
      Out_File << endl;
    }
  }
  Out_File << setprecision(6);
 
}

/******************************************************************//**
 * Routine: Output_Cells_Tecplot                                      
 *                                                                    
 * Writes the cell centred solution values of the specified           
 * quadrilateral solution block to the specified output stream        
 * suitable for plotting with TECPLOT.                                
 *                                                                    
 **********************************************************************/
void Output_Cells_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int CPU_Number,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

  int ng = ON;  // Turn ON to print ghost cells.
  double ddU;

  // Ensure boundary conditions are updated before evaluating
  // solution at the nodes.
  BCs(SolnBlk);

  Linear_Reconstruction_LeastSquares(SolnBlk,LIMITER_VENKATAKRISHNAN);

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D LevelSet Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"psi\" \\ \n"
	     << "\"F\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"sign\" \\ \n"
	     << "\"dpsidx\" \\ \n"
	     << "\"dpsidy\" \\ \n"
	     << "\"kappa\" \\ \n"
	     << "\"gradMag\" \\ \n"
	     << "\"dFdx\" \\ \n"
	     << "\"dFdy\" \\ \n"
	     << "\"ddpsi\" \\ \n"
	     << "\"Eikonal\" \\ \n";
  }
  Out_File << "ZONE T =  \"Block Number = " << CPU_Number << "::" << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.ICu - SolnBlk.ICl + 2*SolnBlk.Nghost*ng + 1 << " \\ \n"
	   << "J = " << SolnBlk.JCu - SolnBlk.JCl + 2*SolnBlk.Nghost*ng + 1 << " \\ \n"
	   << "F = POINT \n";

  for (int j = SolnBlk.JCl-SolnBlk.Nghost*ng; j <= SolnBlk.JCu+SolnBlk.Nghost*ng; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost*ng; i <= SolnBlk.ICu+SolnBlk.Nghost*ng; i++) {
      if (i >= SolnBlk.ICl && i <= SolnBlk.ICu && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	Laplacian_Reconstruction_GreenGauss(SolnBlk,i,j,1,ddU);
      } else {
	ddU = ZERO;
      }
      Out_File << " " << SolnBlk.Grid.Cell[i][j].Xc 
	       << " " << SolnBlk.U[i][j]
	       << " " << SolnBlk.sign[i][j]
	       << " " << SolnBlk.dUdx[i][j].psi
	       << " " << SolnBlk.dUdy[i][j].psi
 	       << " " << SolnBlk.kappa[i][j].psi
	       << " " << SolnBlk.gradMag[i][j].psi
	       << " " << SolnBlk.dUdx[i][j].F
	       << " " << SolnBlk.dUdy[i][j].F
// 	       << " " << sqrt(sqr(SolnBlk.dUdx[i][j].psi) + sqr(SolnBlk.dUdy[i][j].psi)) - ONE
	       << " " << abs( sqrt(max(sqr(max(SolnBlk.dUdxm[i][j].psi,ZERO)),
				       sqr(min(SolnBlk.dUdxp[i][j].psi,ZERO))) +
				   max(sqr(max(SolnBlk.dUdym[i][j].psi,ZERO)),
				       sqr(min(SolnBlk.dUdyp[i][j].psi,ZERO)))) - ONE )
	       << " " << ddU
	       << endl;
    }
  }
  Out_File << setprecision(6);

}

/******************************************************************//**
 * Routine: Output_Nodes_Tecplot                                      
 *                                                                    
 * Writes the node locations of the specified quadrilateral solution  
 * block to the specified output stream suitable for plotting with    
 * TECPLOT.                                                           
 *                                                                    
 **********************************************************************/
void Output_Nodes_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name()
	     << ": 2D LevelSet Nodes\"\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n";
  }
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 2*SolnBlk.Nghost + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 2*SolnBlk.Nghost + 1 << " \\ \n"
	   << "F = POINT \n";

  for (int j = SolnBlk.Grid.JNl-SolnBlk.Nghost; j <= SolnBlk.Grid.JNu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.Grid.INl-SolnBlk.Nghost; i <= SolnBlk.Grid.INu+SolnBlk.Nghost; i++) {
      Out_File << " " << SolnBlk.Grid.Node[i][j].X << endl;
    }
  }

  Out_File << setprecision(6);

}

/******************************************************************//**
 * Routine: Output_Interface_Tecplot                                  
 *                                                                    
 * Writes the interface spline nodes to the specified output stream   
 * suitable for plotting with TECPLOT.                                
 *                                                                    
 **********************************************************************/
void Output_Interface_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
			      const int Number_of_Time_Steps,
			      const double &Time,
			      const int Block_Number,
			      const int Output_Title,
			      ostream &Out_File) {

  // Output node solution data.
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Level Set Interface Nodes "
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"F\" \\ \n"
	     << "\"Fx\" \\ \n"
	     << "\"Fy\" \\ \n"
	     << "\"nx\" \\ \n"
	     << "\"ny\" \\ \n";

  for (int ni = 1; ni <= SolnBlk.Interface_List.Ni; ni++) {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number << ni
	     << "\" \\ \n"
	     << "I = " << SolnBlk.Interface_List[ni].Spline.np << " \\ \n"
	     << "J = " << 1 << " \\ \n"
	     << "F = POINT \n";
    Out_File << setprecision(14);
    for (int np = 0; np < SolnBlk.Interface_List[ni].Spline.np; np++) {
      Out_File.setf(ios::scientific);
      Out_File << " "  << SolnBlk.Interface_List[ni].Spline.Xp[np];
      Out_File << " "  << SolnBlk.Interface_List[ni].Fn(np)
	       << " "  << SolnBlk.Interface_List[ni].F[np].x
	       << " "  << SolnBlk.Interface_List[ni].F[np].y
	       << SolnBlk.Interface_List[ni].normal(np)
	       << endl;
      Out_File.unsetf(ios::scientific);
    }
    Out_File << setprecision(6);
  }

}

/******************************************************************//**
 * Routine: Output_Circle_Tecplot                                     
 *                                                                    
 * Writes the comparison of the exact and computed solutions for a    
 * circle interface to the specified output stream suitable for       
 * plotting with TECPLOT.                                             
 *                                                                    
 **********************************************************************/
void Output_Circle_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
			   LevelSet2D_Input_Parameters &IP,
			   const double &Time,
			   const int Block_Number,
			   const int Output_Title,
			   ostream &Out_File,
			   double &l1_norm,
			   double &l2_norm,
			   double &max_norm) {

  double psi_exact, a, b, r;

  // Determine the coefficients of the exact signed-distance field.
  a = SolnBlk.Interface_List[1].Xref.x;
  b = SolnBlk.Interface_List[1].Xref.y;
  r = SolnBlk.Interface_List[1].Length1;

  if (SolnBlk.Interface_List[1].Motion == INTERFACE_MOTION_LEVELSET_EXPAND) {
    r += Time*SolnBlk.Interface_List[1].Speed.x;
  }

  if (IP.i_BulkFlowField_Type == INTERFACE_BULKFLOWFIELD_UNIFORM) {
    a += Time*IP.V.x;
    b += Time*IP.V.y;
  } else if (IP.i_BulkFlowField_Type == INTERFACE_BULKFLOWFIELD_SWIRL) {
    //a = ;
    //b = ;
  }

  // Determine the error norms and output the solution data.
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name()
	     << ": 2D Level Set solution comparison for a circle interface"
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"psi\" \\ \n"
	     << "\"psi_exact\" \\ \n"
	     << "\"psi_error\" \\ \n"
	     << "\"eikonal_test\" \\ \n";
  Out_File << "ZONE T =  \"Block Number = " << Block_Number << "\" \\ \n"
	   << "I = " << SolnBlk.ICu - SolnBlk.ICl + 1 << " \\ \n"
	   << "J = " << SolnBlk.JCu - SolnBlk.JCl + 1 << " \\ \n"
	   << "F = POINT \\ \n";
  Out_File.setf(ios::scientific);
  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      Linear_Reconstruction_LeastSquares(SolnBlk,i,j,LIMITER_ZERO);
      psi_exact = sqrt(sqr(SolnBlk.Grid.Cell[i][j].Xc.x-a) +
		       sqr(SolnBlk.Grid.Cell[i][j].Xc.y-b)) - r;
      l1_norm += fabs(psi_exact - SolnBlk.U[i][j].psi);
      l2_norm += sqr(psi_exact - SolnBlk.U[i][j].psi);
      max_norm = max(max_norm,fabs(psi_exact - SolnBlk.U[i][j].psi));
      Out_File << SolnBlk.Grid.Cell[i][j].Xc << " " 
	       << SolnBlk.U[i][j].psi << " " 
	       << psi_exact << " "
	       << fabs(psi_exact - SolnBlk.U[i][j].psi) << " "
	       << sqrt(sqr(SolnBlk.dUdx[i][j].psi) + sqr(SolnBlk.dUdy[i][j].psi))
	       << endl;
    }
  }
  Out_File.unsetf(ios::scientific);

}

/******************************************************************//**
 * Routine: Output_Ellipse_Tecplot                                    
 *                                                                    
 * Writes the comparison of the exact and computed solutions for an   
 * ellipse interface to the specified output stream suitable for      
 * plotting with TECPLOT.                                             
 *                                                                    
 **********************************************************************/
void Output_Ellipse_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
			    LevelSet2D_Input_Parameters &IP,
			    const double &Time,
			    const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File,
			    double &l1_norm,
			    double &l2_norm,
			    double &max_norm) {

}

/******************************************************************//**
 * Routine: Output_Zalesaks_Disk_Tecplot                              
 *                                                                    
 * Writes the comparison of the exact and computed solutions for      
 * Zalesak's disk to the specified output stream suitable for         
 * plotting with TECPLOT.                                             
 *                                                                    
 **********************************************************************/
void Output_Zalesaks_Disk_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
				  LevelSet2D_Input_Parameters &IP,
				  const double &Time,
				  const int Block_Number,
				  const int Output_Title,
				  ostream &Out_File,
				  double &l1_norm,
				  double &l2_norm,
				  double &max_norm) {

}
