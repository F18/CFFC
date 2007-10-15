/******************************************************************//**
 * \file LevelSet2DQuadMultiBlock.cc                                  
 *                                                                    
 * Multi-block versions of subroutines for 2D Level Set multi-block   
 * quadrilateral mesh solution classes.                               
 *                                                                    
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Multiple Block External Subroutines.      *
 **********************************************************************/

/******************************************************************//**
 * Routine: Allocate                                                  
 *                                                                    
 * Allocate memory for 1D array of 2D quadrilateral multi-block       
 * solution blocks.                                                   
 *                                                                    
 **********************************************************************/
LevelSet2D_Quad_Block* Allocate(LevelSet2D_Quad_Block *Soln_ptr,
				LevelSet2D_Input_Parameters &Input_Parameters) {
  
  // Allocate memory.
  Soln_ptr = new LevelSet2D_Quad_Block[Input_Parameters.Number_of_Blocks_Per_Processor];
  
  // Return memory location.
  return Soln_ptr;

}

/******************************************************************//**
 * Routine: Deallocate                                                
 *                                                                    
 * Deallocate memory for 1D array of 2D quadrilateral multi-block     
 * solution blocks.                                                   
 *                                                                    
 **********************************************************************/
LevelSet2D_Quad_Block* Deallocate(LevelSet2D_Quad_Block *Soln_ptr,
				  LevelSet2D_Input_Parameters &Input_Parameters) {
  
  // Deallocate memory.
  for (int nb = 0; nb < Input_Parameters.Number_of_Blocks_Per_Processor; nb++)
    if (Soln_ptr[nb].U != NULL) Soln_ptr[nb].deallocate();
  
  delete []Soln_ptr;
  Soln_ptr = NULL;
  
  // Return memory location.
  return Soln_ptr;
  
}

/******************************************************************//**
 * Routine: ICs                                                       
 *                                                                    
 * Initializes the interface(s) on each of the solution blocks in the 
 * 1D array of 2D quadrilateral multi-block solution blocks.          
 *                                                                    
 **********************************************************************/
void ICs(LevelSet2D_Quad_Block *Soln_ptr,
	 AdaptiveBlock2D_List &Soln_Block_List,
	 LevelSet2D_Input_Parameters &Input_Parameters) {
  
  int error_flag;

  error_flag = Initialize_Interfaces(Soln_ptr,Soln_Block_List,Input_Parameters);

}

/******************************************************************//**
 * Routine: Construct_Bulk_Flow_Field                                 
 *                                                                    
 * Constructs the bulk flow field data to the solution variables of a 
 * 1D array of 2D quadrilateral multi-block solution blocks.          
 *                                                                    
 **********************************************************************/
int Construct_Bulk_Flow_Field(LevelSet2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      LevelSet2D_Input_Parameters &Input_Parameters) {

  int error_flag = 0;

  // Assign initial data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Construct_Bulk_Flow_Field(Soln_ptr[nb],Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  // Return error_flag;
  return error_flag;

}

/******************************************************************//**
 * Routine: BCs                                                       
 *                                                                    
 * Apply boundary conditions at boundaries of a 1D array of 2D        
 * quadrilateral multi-block solution blocks.                         
 *                                                                    
 **********************************************************************/
void BCs(LevelSet2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
	 LevelSet2D_Input_Parameters &Input_Parameters) {

  // Prescribe boundary data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      BCs(Soln_ptr[nb]);
    }
  }

}

/******************************************************************//**
 * Routine: Linear_Reconstruction                                     
 *                                                                    
 * This routine determines the linear reconstruction of the solution  
 * gradients for a 1D array of 2D quadrilateral multi-block Level Set 
 * solution blocks.                                                   
 *                                                                    
 **********************************************************************/
// int Linear_Reconstruction(LevelSet2D_Quad_Block *Soln_ptr,
// 			  AdaptiveBlock2D_List &Soln_Block_List,
// 			  LevelSet2D_Input_Parameters &Input_Parameters) {

//   int error_flag = 0;

//   // Perform the linear reconstruction within each cell of the 
//   // computational grid for this stage to determine the gradient of 
//   // the level set function.    
//   for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
//     if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
//       switch(Input_Parameters.i_Reconstruction) {
//       case RECONSTRUCTION_GREEN_GAUSS :
// 	Linear_Reconstruction_GreenGauss(Soln_ptr[nb],Input_Parameters.i_Limiter);    
// 	break;
//       case RECONSTRUCTION_LEAST_SQUARES :
// 	Linear_Reconstruction_LeastSquares(Soln_ptr[nb],Input_Parameters.i_Limiter);
// 	break;
//       default:
// 	Linear_Reconstruction_LeastSquares(Soln_ptr[nb],Input_Parameters.i_Limiter);
// 	break;
//       }
//     }
//   }
  
//   // Quadrilateral multi-block solution blocks successfully updated.
//   return error_flag;
 
// }

/******************************************************************//**
 * Routine: Reconstruction_EssentiallyNonOscillatory                  
 *                                                                    
 * This routine determines the linear reconstruction of the solution  
 * gradients of a given variable, n, for a 1D array of 2D Cartesian   
 * solution blocks using an essentially non-oscillatory scheme.       
 *                                                                    
 **********************************************************************/
int Reconstruction_EssentiallyNonOscillatory(LevelSet2D_Quad_Block *Soln_ptr,
					     AdaptiveBlock2D_List &Soln_Block_List,
					     LevelSet2D_Input_Parameters &Input_Parameters,
					     const int n) {

  int error_flag = 0;
 
  // Perform the linear reconstruction of the specified state variable,
  // n, within each cell of the computational grid.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Reconstruction_EssentiallyNonOscillatory(Soln_ptr[nb],
					       n,
					       Input_Parameters.i_Reconstruction);
    }
  }

  // Cartesian multi-block solution blocks successfully updated.
  return error_flag;

}

/******************************************************************//**
 * Routine: Reconstruction_WeightedEssentiallyNonOscillatory          
 *                                                                    
 * This routine determines the linear reconstruction of the solution  
 * gradients of a given variable, n, for a 1D array of 2D Cartesian   
 * solution blocks using an essentially non-oscillatory scheme.       
 *                                                                    
 **********************************************************************/
int Reconstruction_WeightedEssentiallyNonOscillatory(LevelSet2D_Quad_Block *Soln_ptr,
						     AdaptiveBlock2D_List &Soln_Block_List,
						     const int n) {

  int error_flag = 0;
 
  // Perform the linear reconstruction of the specified state variable,
  // n, within each cell of the computational grid.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Reconstruction_WeightedEssentiallyNonOscillatory(Soln_ptr[nb],
						       n);
    }
  }

  // Cartesian multi-block solution blocks successfully updated.
  return error_flag;

}

/******************************************************************//**
 * Routine: Reconstruction_Curvature                                  
 *                                                                    
 * This routine determines the curvature of a given variable, n, for  
 * a 1D array of 2D Cartesian solution blocks using one of three      
 * methods.                                                           
 *                                                                    
 **********************************************************************/
int Reconstruction_Curvature(LevelSet2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     LevelSet2D_Input_Parameters &Input_Parameters,
			     const int n) {
  int error_flag = 0;

  /* Perform the reconstruction of the curvature of the specified
     state variable, n, within each block of the block list. */
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Reconstruction_Curvature(Soln_ptr[nb],
			       Input_Parameters,
			       n);
    }
  }

  // Cartesian multi-block solution blocks successfully updated.
  return error_flag;

}

/******************************************************************//**
 * Routine: Set_Global_TimeStep                                       
 *                                                                    
 * Assigns global time step to a 1D array of 2D quadrilateral         
 * multi-block solution blocks for time-accurate calculations.        
 *                                                                    
 **********************************************************************/
void Set_Global_TimeStep(LevelSet2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         const double &Dt_min) {
  
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Set_Global_TimeStep(Soln_ptr[nb],Dt_min);
    }
  }
  
}

/******************************************************************//**
 * Routine: L1_Norm_Residual                                          
 *                                                                    
 * Determines the L1-norm of the solution residual for a 1D array of  
 * 2D quadrilateral multi-block solution blocks.  Useful for          
 * monitoring convergence of the solution for steady state problems.  
 *                                                                    
 **********************************************************************/
double L1_Norm_Residual(LevelSet2D_Quad_Block *Soln_ptr,
			AdaptiveBlock2D_List &Soln_Block_List,
			const int &n) {
  
  double l1_norm = ZERO;
  
  // Calculate the L1-norm.  Sum the L1-norm for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_norm +=  L1_Norm_Residual(Soln_ptr[nb],n);
    }
  }

  // Return the L1-norm.
  return l1_norm;
  
}

/******************************************************************//**
 * Routine: L2_Norm_Residual                                          
 *                                                                    
 * Determines the L2-norm of the solution residual for a 1D array of  
 * 2D quadrilateral multi-block solution blocks.  Useful for          
 * monitoring convergence of the solution for steady state problems.  
 *                                                                    
 **********************************************************************/
double L2_Norm_Residual(LevelSet2D_Quad_Block *Soln_ptr,
			AdaptiveBlock2D_List &Soln_Block_List,
			const int &n) {

  double l2_norm = ZERO, l2_norm_get = ZERO;
  
  // Sum the square of the L2-norm for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l2_norm_get = L2_Norm_Residual(Soln_ptr[nb],n);
      l2_norm += sqr(l2_norm_get);
    }
  }

  // Calculate the L2-norm for all blocks.
  l2_norm = sqrt(l2_norm);
  
  // Return the L2-norm.
  return l2_norm;
  
}

/******************************************************************//**
 * Routine: Max_Norm_Residual                                         
 *                                                                    
 * Determines the maximum norm of the solution residual for a 1D      
 * array of 2D quadrilateral multi-block solution blocks.  Useful for 
 * monitoring convergence of the solution for steady state problems.  
 *                                                                    
 **********************************************************************/
double Max_Norm_Residual(LevelSet2D_Quad_Block *Soln_ptr,
			 AdaptiveBlock2D_List &Soln_Block_List,
			 const int &n) {

  double max_norm = ZERO, max_norm_get = ZERO;
  
  // Find the maximum norm for all solution blocks.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      max_norm_get = Max_Norm_Residual(Soln_ptr[nb],n);
      max_norm = max(max_norm,max_norm_get);
    }
  }

  // Return the maximum norm.
  return max_norm;
  
}

/******************************************************************//**
 * Routine: Store_Initial_Eikonal_Solution                            
 *                                                                    
 * This routine makes a copy of the solution prior to solving the     
 * Eikonal equation. This initial solution is used in solving the     
 * equation.                                                          
 *                                                                    
 **********************************************************************/
int Store_Initial_Eikonal_Solution(LevelSet2D_Quad_Block *Soln_ptr,
				   AdaptiveBlock2D_List &Soln_Block_List) {

  int error_flag;

  // Copy the solution for all blocks.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Store_Initial_Eikonal_Solution(Soln_ptr[nb]);
      if (error_flag) return error_flag;
    }
  }

  // Initial Eikonal solution copied successfully.
  return 0;

}

/******************************************************************//**
 * Routine: Calculate_Sign_Function                                   
 *                                                                    
 * This routine determines the sign function of the signed distance   
 * function for a 1D array of 2D quadrilateral multi-block solution   
 * blocks.  Required for the solution of the Eikonal equation for     
 * redistancing the distance function and the scalar extension        
 * equation.                                                          
 *                                                                    
 **********************************************************************/
int Calculate_Sign_Function(LevelSet2D_Quad_Block *Soln_ptr,
			    AdaptiveBlock2D_List &Soln_Block_List,
			    LevelSet2D_Input_Parameters &Input_Parameters) {

  int error_flag;

  // Calculate the sign function for all solution blocks.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Calculate_Sign_Function(Soln_ptr[nb],
					   Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  // Sign function calculated successfully.
  return 0;

}

/******************************************************************//**
 * Routine: Initial_Adaptive_Mesh_Refinement                          
 *                                                                    
 * This routine calls the AMR routine for initial refinement.  This   
 * routine is required to facilitate reapplication of the initial     
 * conditions.                                                         
 *                                                                    
 **********************************************************************/
int Initial_Adaptive_Mesh_Refinement(LevelSet2D_Quad_Block *Soln_ptr,
				     LevelSet2D_Input_Parameters &Input_Parameters,
				     QuadTreeBlock_DataStructure &QuadTree,
				     AdaptiveBlockResourceList &GlobalSolnBlockList,
				     AdaptiveBlock2D_List &Soln_Block_List) {

  if (!Input_Parameters.Number_of_Initial_Mesh_Refinements) return 0;

  int error_flag;

  for (int number_of_initial_mesh_refinements = 1; 
       number_of_initial_mesh_refinements <= Input_Parameters.Number_of_Initial_Mesh_Refinements; 
       number_of_initial_mesh_refinements++) {

    // Call AMR.
    error_flag = AMR(Soln_ptr,
		     Input_Parameters,
		     QuadTree,
		     GlobalSolnBlockList,
		     Soln_Block_List,
		     ON,ON);
    if (error_flag) return error_flag;

    // Reinitialize the interfaces.
    error_flag = Initialize_Interfaces(Soln_ptr,
				       Soln_Block_List,
				       Input_Parameters);
    if (error_flag) return error_flag;

    // Resolve geometric extension problem.
    if (Input_Parameters.i_Initial_Distance_Type == LEVELSET_INITIAL_EXTENSION_GEOMETRIC) {
      error_flag = Geometric_Extension_Problem(Soln_ptr,
					       Soln_Block_List,
					       Input_Parameters);
    } else if (Input_Parameters.i_Initial_Distance_Type == LEVELSET_INITIAL_EXTENSION_EXACT) {
      error_flag = Exact_Initial_Extension(Soln_ptr,
					   Soln_Block_List,
					   Input_Parameters);
    }
    if (error_flag) return error_flag;

    // Redistance level set function to be a signed distance function.
    error_flag = Explicit_Eikonal_Equation(Soln_ptr,
  					   Input_Parameters,
  					   QuadTree,
  					   GlobalSolnBlockList,
  					   Soln_Block_List,
  					   OFF);
    if (error_flag) return error_flag;

    // Resolve scalar (front speed) extension equation.
    error_flag = Explicit_Scalar_Extension_Equation(Soln_ptr,
						    Input_Parameters,
						    QuadTree,
						    GlobalSolnBlockList,
						    Soln_Block_List);
    if (error_flag) return error_flag;

    // Update ghostcell information and prescribe boundary conditions 
    // to ensure that the solution is consistent on each block.
    error_flag = Send_All_Messages(Soln_ptr,
				   Soln_Block_List,
				   Soln_ptr[0].NumVar(),
				   OFF);
    if (error_flag) return error_flag;

    // Print statistics.
    if (CFFC_Primary_MPI_Processor())
      cout << "\n Refinement Level #"      << number_of_initial_mesh_refinements
	   << " : Number of Blocks = "     << QuadTree.countUsedBlocks()
	   << ", Number of Cells = "       << QuadTree.countUsedCells()
	   << ", Refinement Efficiency = " << QuadTree.efficiencyRefinement();

  }

  // Initial interface refined successfully.
  return 0;

}

int Initial_Adaptive_Mesh_Refinement(LevelSet2D_Quad_Block *Soln_ptr,
				     LevelSet2D_Input_Parameters &Input_Parameters,
				     QuadTreeBlock_DataStructure &QuadTree,
				     AdaptiveBlockResourceList &GlobalSolnBlockList,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     const Interface2D_List &Interface_List) {

  if (!Input_Parameters.Number_of_Initial_Mesh_Refinements) return 0;

  int error_flag;

  for (int number_of_initial_mesh_refinements = 1; 
       number_of_initial_mesh_refinements <= Input_Parameters.Number_of_Initial_Mesh_Refinements; 
       number_of_initial_mesh_refinements++) {

    // Call AMR.
    error_flag = AMR(Soln_ptr,
		     Input_Parameters,
		     QuadTree,
		     GlobalSolnBlockList,
		     Soln_Block_List,
		     ON,ON);
    if (error_flag) return error_flag;

    // Reinitialize the interfaces.
    error_flag = Initialize_Interfaces(Soln_ptr,
				       Soln_Block_List,
				       Input_Parameters,
				       Interface_List);
    if (error_flag) return error_flag;

    // Resolve geometric extension problem.
    if (Input_Parameters.i_Initial_Distance_Type == LEVELSET_INITIAL_EXTENSION_GEOMETRIC) {
      error_flag = Geometric_Extension_Problem(Soln_ptr,
					       Soln_Block_List,
					       Input_Parameters);
    } else if (Input_Parameters.i_Initial_Distance_Type == LEVELSET_INITIAL_EXTENSION_EXACT) {
      error_flag = Exact_Initial_Extension(Soln_ptr,
					   Soln_Block_List,
					   Input_Parameters);
    }
    if (error_flag) return error_flag;

    // Redistance level set function to be a signed distance function.
    error_flag = Explicit_Eikonal_Equation(Soln_ptr,
  					   Input_Parameters,
  					   QuadTree,
  					   GlobalSolnBlockList,
  					   Soln_Block_List,
  					   OFF);
    if (error_flag) return error_flag;

    // Resolve scalar (front speed) extension equation.
    error_flag = Explicit_Scalar_Extension_Equation(Soln_ptr,
						    Input_Parameters,
						    QuadTree,
						    GlobalSolnBlockList,
						    Soln_Block_List);
    if (error_flag) return error_flag;

    // Update ghostcell information and prescribe boundary conditions 
    // to ensure that the solution is consistent on each block.
    error_flag = Send_All_Messages(Soln_ptr,
				   Soln_Block_List,
				   Soln_ptr[0].NumVar(),
				   OFF);
    if (error_flag) return error_flag;

    // Print statistics.
    if (CFFC_Primary_MPI_Processor())
      cout << "\n Refinement Level #"      << number_of_initial_mesh_refinements
	   << " : Number of Blocks = "     << QuadTree.countUsedBlocks()
	   << ", Number of Cells = "       << QuadTree.countUsedCells()
	   << ", Refinement Efficiency = " << QuadTree.efficiencyRefinement();

  }

  // Initial interface refined successfully.
  return 0;

}

/******************************************************************//**
 * Routine: Uniform_Adaptive_Mesh_Refinement                          
 *                                                                    
 * This routine calls the AMR routine for uniform refinement.  This   
 * routine is required to facilitate reapplication of the uniform     
 * conditions.                                                         
 *                                                                    
 **********************************************************************/
int Uniform_Adaptive_Mesh_Refinement(LevelSet2D_Quad_Block *Soln_ptr,
				     LevelSet2D_Input_Parameters &Input_Parameters,
				     QuadTreeBlock_DataStructure &QuadTree,
				     AdaptiveBlockResourceList &GlobalSolnBlockList,
				     AdaptiveBlock2D_List &Soln_Block_List) {

  if (!Input_Parameters.Number_of_Uniform_Mesh_Refinements) return 0;

  int error_flag;

  // Call Uniform AMR.
  error_flag = Uniform_AMR(Soln_ptr,
			   Input_Parameters,
			   QuadTree,
			   GlobalSolnBlockList,
			   Soln_Block_List);
  if (error_flag) return error_flag;

  // Reinitialize the interfaces.
  error_flag = Initialize_Interfaces(Soln_ptr,
				     Soln_Block_List,
				     Input_Parameters);
  if (error_flag) return error_flag;

  // Resolve geometric extension problem.
  error_flag = Geometric_Extension_Problem(Soln_ptr,
					   Soln_Block_List,
					   Input_Parameters);
  if (error_flag) return error_flag;

  // Redistance level set function to be a signed distance function.
  error_flag = Explicit_Eikonal_Equation(Soln_ptr,
					 Input_Parameters,
					 QuadTree,
					 GlobalSolnBlockList,
					 Soln_Block_List,
					 OFF);
  if (error_flag) return error_flag;

  // Resolve scalar (front speed) extension equation.
  error_flag = Explicit_Scalar_Extension_Equation(Soln_ptr,
						  Input_Parameters,
						  QuadTree,
						  GlobalSolnBlockList,
						  Soln_Block_List);
  if (error_flag) return error_flag;

  // Update ghostcell information and prescribe boundary conditions 
  // to ensure that the solution is consistent on each block.
  error_flag = Send_All_Messages(Soln_ptr,
				 Soln_Block_List,
				 Soln_ptr[0].NumVar(),
				 OFF);
  if (error_flag) return error_flag;

  // Uniform mesh refinement was successfull.
  return 0;

}

int Uniform_Adaptive_Mesh_Refinement(LevelSet2D_Quad_Block *Soln_ptr,
				     LevelSet2D_Input_Parameters &Input_Parameters,
				     QuadTreeBlock_DataStructure &QuadTree,
				     AdaptiveBlockResourceList &GlobalSolnBlockList,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     const Interface2D_List &Interface_List) {

  if (!Input_Parameters.Number_of_Uniform_Mesh_Refinements) return 0;

  int error_flag;

  // Call Uniform AMR.
  error_flag = Uniform_AMR(Soln_ptr,
			   Input_Parameters,
			   QuadTree,
			   GlobalSolnBlockList,
			   Soln_Block_List);
  if (error_flag) return error_flag;

  // Reinitialize the interfaces.
  error_flag = Initialize_Interfaces(Soln_ptr,
				     Soln_Block_List,
				     Input_Parameters,
				     Interface_List);
  if (error_flag) return error_flag;

  // Resolve geometric extension problem.
  error_flag = Geometric_Extension_Problem(Soln_ptr,
					   Soln_Block_List,
					   Input_Parameters);
  if (error_flag) return error_flag;

  // Redistance level set function to be a signed distance function.
  error_flag = Explicit_Eikonal_Equation(Soln_ptr,
					 Input_Parameters,
					 QuadTree,
					 GlobalSolnBlockList,
					 Soln_Block_List,
					 OFF);
  if (error_flag) return error_flag;

  // Resolve scalar (front speed) extension equation.
  error_flag = Explicit_Scalar_Extension_Equation(Soln_ptr,
						  Input_Parameters,
						  QuadTree,
						  GlobalSolnBlockList,
						  Soln_Block_List);
  if (error_flag) return error_flag;

  // Update ghostcell information and prescribe boundary conditions 
  // to ensure that the solution is consistent on each block.
  error_flag = Send_All_Messages(Soln_ptr,
				 Soln_Block_List,
				 Soln_ptr[0].NumVar(),
				 OFF);
  if (error_flag) return error_flag;

  // Uniform mesh refinement was successful.
  return 0;

}
