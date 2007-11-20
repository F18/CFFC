/*!\file New_AdvectDiffuse2DQuad.h
  \brief Temporary header file defining 2D Advection Diffusion Equation Quadrilateral Mesh Solution Classes. */

#ifndef _NEW_ADVECTDIFFUSE2D_QUAD_INCLUDED
#define _NEW_ADVECTDIFFUSE2D_QUAD_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "New_AdvectDiffuse2DState.h"  /* Include 2D advection diffusion equation solution state header file */
#include "../Grid/Cell2D.h"        /* Include 2D cell header file */
#include "../Grid/Grid2DQuad.h"    /* Include 2D quadrilateral multiblock grid header file */
#include "../AMR/QuadTree.h"       /* Include quadtree header file */
#include "../AMR/AMR.h"            /* Include AMR header file */
#include "AdvectDiffuse2DInput.h"  /* Include 2D advection diffusion equation header file */
#include "../ICEM/ICEMCFD.h"       /* Include ICEMCFD header file. */

/* Define the structures and classes. */

#define	NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D    3

/*!
 * Class: AdvectDiffuse2D_Quad_Block_New
 *
 * @brief Class definition of the 2D advection-diffusion solution blocks.
 *
 * \verbatim
 * Member functions
 *       U      -- Return solution variable for the block (cell average).
 *    Grid      -- Return the solution block quadrilateral grid.
 *      dt      -- Return local time step for the solution block.
 *    dudt      -- Return the solution residuals.
 *    dudx      -- Return the unlimited solution gradients (x-direction)
 *                 for block.
 *    dudy      -- Return the unlimited solution gradients (y-direction)
 *                 for block.
 *     phi      -- Return the solution slope limiter.
 *    FluxN     -- Return array of north boundary solution fluxes.
 *    FluxS     -- Return array of south boundary solution fluxes.
 *    FluxE     -- Return array of east boundary solution fluxes.
 *    FluxW     -- Return array of west boundary solution fluxes.
 *     NCi      -- Return number of cells in
 *                 the i-direction (zeta-direction).
 *     ICl      -- Return lower index for cells in 
 *                 the i-direction (zeta-direction).
 *     ICu      -- Return upper index for cells in 
 *                 the i-direction (zeta-direction).
 *     NCj      -- Return number of cells in 
 *                 the j-direction (eta-direction).
 *     JCl      -- Return lower index for cells in 
 *                 the j-direction (eta-direction).
 *     JCu      -- Return upper index for cells in 
 *                 the j-direction (eta-direction).
 *  Nghost      -- Return number of ghost (halo or 
 *                 overlap) cells.
 * Axisymmetric -- Return axisymmetric geometry
 *                 indicator (=1 for axisymmetric,
 *                            =0 for planar geometry).
 * Freeze_Limiter -- Return limiter freezing indicator
 *                 (=1 for limiter freezing on,
 *                 (=0 for limiter freezing off).
 *     UoN      -- Return array of reference states for the application
 *                 of boundary conditions at the north boundary of the
 *                 solution block.
 *     UoS      -- Return array of reference states for the application
 *                 of boundary conditions at the south boundary of the
 *                 solution block.
 *     UoE      -- Return array of reference states for the application
 *                 of boundary conditions at the east boundary of the
 *                 solution block.
 *     UoW      -- Return array of reference states for the application
 *                 of boundary conditions at the west boundary of the
 *                 solution block.
 *   allocate   -- Allocate memory for structured quadrilateral
 *                 solution block.
 *   deallocate -- Deallocate memory for structured quadrilateral
 *                 solution block.
 *      Un      -- Return solution state at the specified node.
 *    UnNW      -- Return solution state at the north-west node.
 *    UnNE      -- Return solution state at the north-east node.
 *    UnSW      -- Return solution state at the south-west node.
 *    UnSE      -- Return solution state at the south-east node.
 *      un      -- Return solution value at the specified node.
 *    unNW      -- Return solution value at the north-west node.
 *    unNE      -- Return solution value at the north-east node.
 *    unSW      -- Return solution value at the south-west node.
 *    unSE      -- Return solution value at the south-east node.
 * evaluate_limiters -- Set flag to evaluate limiters.
 * freeze_limiters -- Set flag to freeze limiters.
 * Member functions required for message passing.
 *    NumVar    -- Returns number of solution variables in primitive
 *                 and conserved solution state vectors.
 * LoadSendBuffer -- Loads send buffer.
 * LoadSendBuffer_F2C -- Loads send buffer for fine to coarse block
 *                       messages.
 * UnloadReceiveBuffer -- Unloads receive buffer.
 * UnloadReceiveBuffer_F2C -- Unloads receive buffer for fine to
 *                            coarse block messages.
 * SubcellReconstruction -- Performs subcell solution reconstruction
 *                          used in adaptive mesh refinement.
 * LoadSendBuffer_Flux_F2C -- Loads send buffer for sending conservative
 *                            flux corrections from fine to coarse
 *                            solution blocks.
 * UnLoadSendBuffer_Flux_F2C -- Loads send buffer for sending 
 *                              conservative flux corrections from fine 
 *                              to coarse solution blocks.
 *
 * Member operators
 *      S -- a 2D advection diffusion equation solution
 *
 * S = S;
 * cout << S; (output function)
 * cin  >> S; (input function)
 * \endverbatim
 */
class AdvectDiffuse2D_Quad_Block_New{
private:
  AdvectDiffuse2D_Quad_Block_New(const AdvectDiffuse2D_Quad_Block_New &Soln); //!< Private copy constructor
  AdvectDiffuse2D_Quad_Block_New operator = (const AdvectDiffuse2D_Quad_Block_New &Soln);   //!< Private assignment operator

public:
  //! @name Defined datatypes
  //@{
  typedef Vector2D (* Exact_Gradient_Function) (const double, const double);
  //@}

  //! @name Solution state arrays:
  //@{
  AdvectDiffuse2D_State_New    **U; //!< Solution state.
  //@}

  //! @name Grid block information:
  //@{
  int                      NCi, //!< Total number of i-direction cells.
    ICl, //!< First i-direction non-ghost cell counter.
    ICu; //!< Final i-direction non-ghost cell counter.
  int                      NCj, //!< Total number of j-direction cells.
    JCl, //!< First j-direction non-ghost cell counter.
    JCu; //!< Final j-direction non-ghost cell counter.
  int                   Nghost; //!< Number of ghost cells.
  Grid2D_Quad_Block       Grid; //!< 2D quadrilateral grid geometry.
  //@}

  //! @name Residual and time-stepping arrays:
  //@{
  double                         **dt; //!< Local time step.
  AdvectDiffuse2D_State_New   ***dUdt; //!< Solution residual.
  AdvectDiffuse2D_State_New      **Uo; //!< Initial solution.
  static int residual_variable; //!< Static integer that indicates which variable is used for residual calculations.
  //@}

  //! @name Solution gradient arrays:
  //@{
  AdvectDiffuse2D_State_New    **dUdx; //!< Unlimited solution gradient (x-direction).
  AdvectDiffuse2D_State_New    **dUdy; //!< Unlimited solution gradient (y-direction).
  AdvectDiffuse2D_State_New     **phi; //!< Solution slope limiter.
  //@}

  //! @name Boundary solution flux arrays:
  //@{
  AdvectDiffuse2D_State_New   *FluxN, //!< North boundary solution flux.
                              *FluxS, //!< South boundary solution flux.
                              *FluxE, //!< East boundary solution flux.
                              *FluxW; //!< West boundary solution flux.
  //@}

  //! @name Problem indicator flags:
  //@{
  int             Axisymmetric; //!< Axisymmetric geometry indicator.
  int           Freeze_Limiter; //!< Limiter freezing indicator.
  static char*   solutionTitle; //!< Solution title info
  //@}

  //! @name Boundary condtion reference states:
  //@{
  AdvectDiffuse2D_State_New   *UoN, //!< Boundary condition reference states for north boundary.
                              *UoS, //!< Boundary condition reference states for south boundary.
                              *UoE, //!< Boundary condition reference states for east boundary.
                              *UoW; //!< Boundary condition reference states for west boundary.
  //@}

  //! @name Pointers to exact solutions
  //@{
  static Exact_Gradient_Function ExactGrad; //!< Exact gradient for equations with analytic solution
  static FunctionType2D ExactSoln;          //!< Exact solution for equations with analytic solution
  //@}
	      
  //! @name Creation, copy, and assignment constructors.
  //@{
  //! Default constructor.
  AdvectDiffuse2D_Quad_Block_New(void);

  /* Destructor. */
  ~AdvectDiffuse2D_Quad_Block_New(void){ deallocate(); }
  //@}

  //! @name Allocate and deallocate functions.
  //@{
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);

  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //@}

  //! @name Bilinear interplation (Zingg & Yarrow).
  //@{
  //! Return solution state at specified node.
  AdvectDiffuse2D_State_New Un(const int ii, const int jj);

  AdvectDiffuse2D_State_New UnNW(const int ii, const int jj); //!< Return solution state at cell NW node.
  AdvectDiffuse2D_State_New UnNE(const int ii, const int jj); //!< Return solution state at cell NE node.
  AdvectDiffuse2D_State_New UnSE(const int ii, const int jj); //!< Return solution state at cell SE node.
  AdvectDiffuse2D_State_New UnSW(const int ii, const int jj); //!< Return solution state at cell SW node.

  //! Return solution state at specified node.
  double un(const int ii, const int jj);

  double unNW(const int ii, const int jj); //!< Return solution state at cell NW node.
  double unNE(const int ii, const int jj); //!< Return solution state at cell NE node.
  double unSE(const int ii, const int jj); //!< Return solution state at cell SE node.
  double unSW(const int ii, const int jj); //!< Return solution state at cell SW node.
  //@}

  //! @name Evaluate diffusive flux for the cell.
  //@{
  //@}

  //! @name Evaluate velocity and diffusion coefficient at the centroid of a specified cell.
  //@{
  Vector2D VelocityAtCellCentroid(const int & ii, const int & jj) const;     //!< Calculate velocity at the cell centroid
  double DiffusionCoeffAtCellCentroid(const int & ii, const int & jj) const; //!< Calculate diffusion coeff. at the cell centroid
  //@}

  //! Determine the stability limit imposed by the source term
  double SourceTermStabilityLimit(const int & ii, const int & jj) const;

  //! @name Member functions for limiter freezing.
  //@{
  void evaluate_limiters(void);
  void freeze_limiters(void);
  //@}

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file,
			       const AdvectDiffuse2D_Quad_Block_New &Soln);
  friend istream &operator >> (istream &in_file,
			       AdvectDiffuse2D_Quad_Block_New &Soln);
  //@}

  //! @name Member functions required for message passing.
  //@{
  //! Number of solution state variables.
  int NumVar(void);
  //! Load send message passing buffer.
  int LoadSendBuffer(double *buffer,
		     int &buffer_count,
		     const int buffer_size,
		     const int i_min, 
		     const int i_max,
		     const int i_inc,
		     const int j_min, 
		     const int j_max,
		     const int j_inc);
  //! Load F2C send message passing buffer.
  int LoadSendBuffer_F2C(double *buffer,
			 int &buffer_count,
			 const int buffer_size,
			 const int i_min, 
			 const int i_max,
			 const int i_inc,
			 const int j_min, 
			 const int j_max,
			 const int j_inc);
  //! Load C2F send message passing buffer.
  int LoadSendBuffer_C2F(double *buffer,
			 int &buffer_count,
			 const int buffer_size,
			 const int i_min, 
			 const int i_max,
			 const int i_inc,
			 const int j_min, 
			 const int j_max,
			 const int j_inc,
			 const int face,
			 const int sector);
  //! Unload receive message passing buffer.
  int UnloadReceiveBuffer(double *buffer,
			  int &buffer_count,
			  const int buffer_size,
			  const int i_min, 
			  const int i_max,
			  const int i_inc,
			  const int j_min, 
			  const int j_max,
			  const int j_inc);
  //! Unload F2C receive message passing buffer.
  int UnloadReceiveBuffer_F2C(double *buffer,
			      int &buffer_count,
			      const int buffer_size,
			      const int i_min, 
			      const int i_max,
			      const int i_inc,
			      const int j_min, 
			      const int j_max,
			      const int j_inc);
  //! Unload C2F receive message passing buffer.
  int UnloadReceiveBuffer_C2F(double *buffer,
			      int &buffer_count,
			      const int buffer_size,
			      const int i_min, 
			      const int i_max,
			      const int i_inc,
			      const int j_min, 
			      const int j_max,
			      const int j_inc);
  //! Subcell solution reconstruction within given computational cell.
  void SubcellReconstruction(const int i,
			     const int j,
			     const int Limiter);
  //! Load F2C conservative flux message passing buffer.
  int LoadSendBuffer_Flux_F2C(double *buffer,
			      int &buffer_count,
			      const int buffer_size,
			      const int i_min, 
			      const int i_max,
			      const int i_inc,
			      const int j_min, 
			      const int j_max,
			      const int j_inc);
  //! Unload F2C conservative flux message passing buffer.
  int UnloadReceiveBuffer_Flux_F2C(double *buffer,
				   int &buffer_count,
				   const int buffer_size,
				   const int i_min, 
				   const int i_max,
				   const int i_inc,
				   const int j_min, 
				   const int j_max,
				   const int j_inc);
  //@}

};

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::Un?? -- Get cell node solution states.     *
 **************************************************************************/
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_Quad_Block_New::UnNW(const int ii, const int jj) {
  return (Un(ii, jj+1));
}

inline AdvectDiffuse2D_State_New AdvectDiffuse2D_Quad_Block_New::UnNE(const int ii, const int jj) {
  return (Un(ii+1, jj+1));
}

inline AdvectDiffuse2D_State_New AdvectDiffuse2D_Quad_Block_New::UnSE(const int ii, const int jj) {
  return (Un(ii+1, jj));
}

inline AdvectDiffuse2D_State_New AdvectDiffuse2D_Quad_Block_New::UnSW(const int ii, const int jj) {
  return (Un(ii, jj));
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::un?? -- Get cell node solution values.     *
 **************************************************************************/
inline double AdvectDiffuse2D_Quad_Block_New::unNW(const int ii, const int jj) {
  return (un(ii, jj+1));
}

inline double AdvectDiffuse2D_Quad_Block_New::unNE(const int ii, const int jj) {
  return (un(ii+1, jj+1));
}

inline double AdvectDiffuse2D_Quad_Block_New::unSE(const int ii, const int jj) {
  return (un(ii+1, jj));
}

inline double AdvectDiffuse2D_Quad_Block_New::unSW(const int ii, const int jj) {
  return (un(ii, jj));
}

/**********************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::evaluate_limiters -- Set flag to evaluate limiters.*
 **********************************************************************************/
inline void AdvectDiffuse2D_Quad_Block_New::evaluate_limiters(void) {
  Freeze_Limiter = OFF; 
}

/**********************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::freeze_limiters -- Set flag to freeze limiters.    *
 **********************************************************************************/
inline void AdvectDiffuse2D_Quad_Block_New::freeze_limiters(void) {
  Freeze_Limiter = ON; 
}

/*******************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::NumVar -- Returns number of state variables.    *
 *******************************************************************************/
inline int AdvectDiffuse2D_Quad_Block_New::NumVar(void) {
  return (int(NUM_VAR_ADVECTDIFFUSE2D));
}


/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New -- Single Block External Subroutines.       *
 **************************************************************************/
extern void Write_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
	                        istream &In_File);

extern void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                     MPI::Intracomm &Communicator, 
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk1,
		                AdvectDiffuse2D_Quad_Block_New &SolnBlk2);

extern int Prolong_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk_Fine,
				  AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk_Coarse,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_SW,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_SE,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_NW,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_NE);

extern void BCs(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern void Output_Tecplot(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
			   AdvectDiffuse2D_Input_Parameters &IP,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
				 AdvectDiffuse2D_Input_Parameters &IP,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nodes_Tecplot(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void ICs(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
 	        const int i_ICtype,
                AdvectDiffuse2D_State_New *Wo);

extern double CFL(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(AdvectDiffuse2D_Quad_Block_New &SolnBlk, 
                                const double &Dt_min);

extern double L1_Norm_Residual(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern double L2_Norm_Residual(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern double Max_Norm_Residual(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern void Linear_Reconstruction_GreenGauss(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
					       const int Limiter);

extern void Residual_Smoothing(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                               const int k_residual,
			       double &epsilon, 
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  AdvectDiffuse2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern void Fix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block_New SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block_New SolnBlk);

extern void Apply_Boundary_Flux_Corrections(AdvectDiffuse2D_Quad_Block_New,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New SolnBlk,
                                                                const int i_stage,
                                                                const int n_stage,
                                                                const double &CFL_Number,
                                                                const int Time_Integration_Type,
                                                                const int Reconstruction_Type,
                                                                const int Limiter_Type,
                                                                const int Number_Neighbours_North_Boundary,
                                                                const int Number_Neighbours_South_Boundary,
                                                                const int Number_Neighbours_East_Boundary,
                                                                const int Number_Neighbours_West_Boundary);

extern int dUdt_Residual_Evaluation(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
				    AdvectDiffuse2D_Input_Parameters &Input_Parameters);

#if 0

extern int dUdt_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
   	                            const int i_stage,
                                    AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Update_Solution_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
   	                                       const int i_stage,
                                               AdvectDiffuse2D_Input_Parameters &Input_Parameters);

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New -- Multiple Block External Subroutines.     *
 **************************************************************************/

extern AdvectDiffuse2D_Quad_Block_New* Allocate(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
						AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern AdvectDiffuse2D_Quad_Block_New* Deallocate(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
						  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void ICs(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
		AdaptiveBlock2D_List &Soln_Block_List,
		AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		QuadTreeBlock_DataStructure &QuadTree);

extern void Set_Advection_Velocity_Field(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                         AdaptiveBlock2D_List &Soln_Block_List,
                                         AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Analytical_Solution(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
				    const AdaptiveBlock2D_List &Soln_Block_List,
				    const AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Read_Restart_Solution(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Nodes_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Mesh_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern void BCs(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern double CFL(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
                  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, 
                                const double &Dt_min);

extern double L1_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void Evaluate_Limiters(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Residual_Smoothing(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Local_Soln_Block_List,
                                    AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

extern int Initial_AMR2(AdvectDiffuse2D_Quad_Block_New       *Soln_ptr,
                        AdvectDiffuse2D_Input_Parameters &InputParameters,
                        QuadTreeBlock_DataStructure      &QuadTree,
                        AdaptiveBlockResourceList        &GlobalSolnBlockList,
                        AdaptiveBlock2D_List             &LocalSolnBlockList);

extern int Uniform_AMR2(AdvectDiffuse2D_Quad_Block_New       *Soln_ptr,
                        AdvectDiffuse2D_Input_Parameters &InputParameters,
                        QuadTreeBlock_DataStructure      &QuadTree,
                        AdaptiveBlockResourceList        &GlobalSolnBlockList,
                        AdaptiveBlock2D_List             &LocalSolnBlockList);

extern int Boundary_AMR2(AdvectDiffuse2D_Quad_Block_New       *Soln_ptr,
                         AdvectDiffuse2D_Input_Parameters &InputParameters,
                         QuadTreeBlock_DataStructure      &QuadTree,
                         AdaptiveBlockResourceList        &GlobalSolnBlockList,
                         AdaptiveBlock2D_List             &LocalSolnBlockList);

/********************************************************************************
 * AdvectDiffuse2D_Quad_Block_New -- Multiple Block External Subroutines for Mesh.  *
 ********************************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                            AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                 AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters);

#endif


#endif /* _ADVECTDIFFUSE2D_QUAD_INCLUDED  */
