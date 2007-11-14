/*! \file AdvectDiffuse2DState.cc
  @brief Subroutines for 2D Advection Diffusion Equation Solution State Class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "New_AdvectDiffuse2DState.h" /* Include 2D advection diffusion equation solution state header file. */

// ===  Static member variables ===
AdvectDiffuse2D_State_New::AdvectionVelocityType AdvectDiffuse2D_State_New::V = VelocityFields::Quiescent_Flow;
AdvectDiffuse2D_State_New::DiffusionFieldType AdvectDiffuse2D_State_New::k = DiffusionFields::Zero_Diffusion;
SourceTermFields* AdvectDiffuse2D_State_New::SourceTerm = NULL;

