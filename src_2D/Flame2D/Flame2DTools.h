/////////////////////////////////////////////////////////////////////
///
/// \file Flame2DTools.cc
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the function definition for 
///        the Flame2DTools.  Basically functions that do not fit
///        anywhere else.
///
/////////////////////////////////////////////////////////////////////
#ifndef _FLAME2D_TOOLS_INCLUDED
#define _FLAME2D_TOOLS_INCLUDED

// inludes
#include "Flame2DQuad.h"

/////////////////////////////////////////////////////////////////////
/// Time Accurate Ouput
/////////////////////////////////////////////////////////////////////
int Open_Time_Accurate_File(ofstream &Time_Accurate_File,
			    char *File_Name,
			    const int Append_to_Fileconst,
			    const Flame2D_pState &Soln);
  
int Close_Time_Accurate_File(ofstream &Time_Accurate_File);
  
void Output_to_Time_Accurate_File(ostream &Time_Accurate_File,
				  const double &Time,
				  Flame2D_pState &Soln);


/////////////////////////////////////////////////////////////////////
/// Exact Solutions
/////////////////////////////////////////////////////////////////////
int Output_Viscous_Channel(Flame2D_Quad_Block *Soln_ptr,
			   AdaptiveBlock2D_List &Soln_Block_List,
			   Flame2D_Input_Parameters &IP);
void Output_Viscous_Channel(Flame2D_Quad_Block &SolnBlk,
			    const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File,
			    double &l1_norm,
			    double &l2_norm,
			    double &max_norm,
			    double &Vwall,
			    const double dp);
  
int Output_Flat_Plate(Flame2D_Quad_Block *Soln_ptr,
		      AdaptiveBlock2D_List &Soln_Block_List,
		      Flame2D_Input_Parameters &IP);
void Output_Flat_Plate(Flame2D_Quad_Block &SolnBlk,
		       const int Block_Number,
		       const int Output_Title_Soln,
		       ostream &Out_File_Soln,
		       const int Output_Title_Skin,
		       ostream &Out_File_Skin,
		       Flame2D_pState &Winf,
		       double &l1_norm,
		       double &l2_norm,
		       double &max_norm);
  
int Output_Driven_Cavity_Flow(Flame2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      Flame2D_Input_Parameters &IP);
  
void Output_Driven_Cavity_Flow(Flame2D_Quad_Block &SolnBlk,
			       const int Block_Number,
			       const int Output_Title,
			       ostream &Out_File_u,
			       ostream &Out_File_v,
			       const double &Re,
			       const double &Vwall,
			       const double &length);
  
  
#endif //_FLAME2D_DRDU_INCLUDED
