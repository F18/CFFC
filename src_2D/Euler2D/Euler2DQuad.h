/* Euler2DQuad.h:  Header file defining 
                   2D Euler Quadrilateral Mesh Solution Classes. */

#ifndef _EULER2D_QUAD_INCLUDED
#define _EULER2D_QUAD_INCLUDED

/* Include 2D Euler state, 2D cell, 2D quadrilateral multiblock 
   grid, quadtree, AMR, NASA rotor, and 2D Euler input header files. */

#include "Euler2DState.h"
#include "../Grid/Cell2D.h"
#include "../Grid/Grid2DQuad.h"
#include "../Grid/HO_Grid2DQuad.h"     /* Include 2D quadrilateral block grid header file */
#include "../AMR/QuadTree.h"
#include "../AMR/AMR.h"
#include "../Grid/NASARotor37.h"
#include "../Grid/NASARotor67.h"
#include "Euler2DInput.h"
#include "../Math/LinearSystems.h"  /* Include the linear systems header file. */
#include "../ICEM/ICEMCFD.h"        /* Include ICEMCFD input header file. */
#include "../System/System_Linux.h"    /* Include System Linux header file. */
#include "../HighOrderReconstruction/AccuracyAssessment2D.h" /* Include 2D accuracy assessment header file. */
#include "../HighOrderReconstruction/HighOrder2D.h" /* Include 2D high-order template class header file. */
#include "../HighOrderReconstruction/Cauchy_BoundaryConditions.h" /* Include 2D high-order boundary conditions header file. */
#include "Euler2DExactSolutions.h"  /* Include 2D Euler exact solutions header file */

/* Define the structures and classes. */

#define	NUMBER_OF_RESIDUAL_VECTORS_EULER2D    3

/*!
 * Class: Euler2D_Quad_Block
 *
 * @brief Class definition of the 2D Euler solution blocks.
 *
 * \verbatim
 * Member functions
 *       W      -- Return primitive variable solution for the block
 *                 (cell average).
 *       U      -- Return conserved variable solution for the block
 *                 (cell average).
 *    Grid      -- Return the solution block quadrilateral grid or mesh.
 *      dt      -- Return local time step for the solution block.
 *    dUdt      -- Return the solution block residuals.
 *    dWdx      -- Return the unlimited primitive variable solution
 *                 gradients (x-direction) for the block.
 *    dWdy      -- Return the unlimited primitive variable solution
 *                 gradients (y-direction) for the block.
 *     phi      -- Return the solution slope limiters.
 *      Uo      -- Return initial solution state.
 *    FluxN     -- Return array of north boundary solution fluxes.
 *    FluxS     -- Return array of south boundary solution fluxes.
 *    FluxE     -- Return array of east boundary solution fluxes.
 *    FluxW     -- Return array of west boundary solution fluxes.
 *     NCi      -- Return number of cells in the i-direction 
 *                 (zeta-direction).
 *     ICl      -- Return lower index for cells in the i-direction 
 *                 (zeta-direction).
 *     ICu      -- Return upper index for cells in the i-direction 
 *                 (zeta-direction).
 *     NCj      -- Return number of cells in the j-direction 
 *                 (eta-direction).
 *     JCl      -- Return lower index for cells in the j-direction 
 *                 (eta-direction).
 *     JCu      -- Return upper index for cells in the j-direction 
 *                 (eta-direction).
 *  Nghost      -- Number of ghost cells.
 * Axisymmetric -- Return axisymmetric flow indicator (= 1 for 
 *                 axisymmetric flow, = 0 for planar flow).
 * Freeze_Limiter -- Return limiter freezing indicator
 *                   (=1 for limiter freezing on,
 *                   (=0 for limiter freezing off).
 *     WoN      -- Return array of reference states for the
 *                 application of boundary conditions at the north
 *                 boundary of the solution block.
 *     WoS      -- Return array of reference states for the
 *                 application of boundary conditions at the south
 *                 boundary of the solution block.
 *     WoE      -- Return array of reference states for the
 *                 application of boundary conditions at the east
 *                 boundary of the solution block.
 *     WoW      -- Return array of reference states for the
 *                 application of boundary conditions at the west
 *                 boundary of the solution block.
 *   allocate   -- Allocate memory for structured quadrilateral
 *                 solution block.
 *   deallocate -- Deallocate memory for structured quadrilateral
 *                 solution block.
 *      Wn      -- Return primitive variable solution at the
 *                 specified node.
 *      Un      -- Return conserved variable solution at the
 *                 specified node.
 *    WnNW      -- Return primitive variable solution at the
 *                 north-west node.
 *    WnNE      -- Return primitive variable solution at the
 *                 north-east node.
 *    WnSW      -- Return primitive variable solution at the
 *                 south-west node.
 *    WnSE      -- Return primitive variable solution at the
 *                 south-east node.
 *    UnNW      -- Return conserved variable solution at the
 *                 north-west node.
 *    UnNE      -- Return conserved variable solution at the
 *                 north-east node.
 *    UnSW      -- Return conserved variable solution at the
 *                 south-west node.
 *    UnSE      -- Return conserved variable solution at the
 *                 south-east node.
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
 * QuadraticSubcellReconstruction -- Performs quadratic subcell solution
 *                                   reconstruction used in adaptive
 *                                   mesh refinement.
 * LoadSendBuffer_Flux_F2C -- Loads send buffer for sending
 *                            conservative flux corrections from fine
 *                            to coarse solution blocks.
 * UnLoadSendBuffer_Flux_F2C -- Loads send buffer for sending
 *                              conservative flux corrections from
 *                              fine to coarse solution blocks.
 *
 * Member operators
 *      S -- a 2D NavierStokes solution
 *
 * S = S;
 * cout << S; (output function)
 * cin  >> S; (input function)
 * \endverbatim
 */
class Euler2D_Quad_Block{
public:

  //! @name Defined public types:
  //@{
  typedef Euler2D_ExactSolutions Exact_Solution_Type;
  typedef AccuracyAssessment2D<Euler2D_Quad_Block> Accuracy_Assessment_Type;
  typedef Euler2D_pState Soln_State;
  typedef HighOrder2D<Soln_State> HighOrderType; //!< high-order variable data type. Use primitive variables for reconstruction.
  typedef Cauchy_BCs<Soln_State> BC_Type;
  typedef std::vector<double> DoubleArrayType;
  //@}

  //@{ @name Solution state arrays:
  Euler2D_pState           **W; //!< Primitive solution state.
  Euler2D_cState           **U; //!< Conserved solution state.
  //@}

  //@{ @name Grid block information:
  int                      NCi, //!< Total number of i-direction cells.
                           ICl, //!< First i-direction non-ghost cell counter.
                           ICu; //!< Final i-direction non-ghost cell counter.
  int                      NCj, //!< Total number of j-direction cells.
                           JCl, //!< First j-direction non-ghost cell counter.
                           JCu; //!< Final j-direction non-ghost cell counter.
  int                   Nghost; //!< Number of ghost cells.
  Grid2D_Quad_Block_HO    Grid; //!< 2D quadrilateral grid geometry.
  //@}

  //@{ @name Residual and time-stepping arrays:
  double                  **dt; //!< Local time step.
  Euler2D_cState       ***dUdt; //!< Solution residual.
  Euler2D_cState          **Uo; //!< Initial solution state.
  static int residual_variable; //!< Static integer that indicates which variable is used for residual calculations.  
  static int Number_of_Residual_Norms; //!< How many Residual norms to plot?
  //@}

  //@{ @name Solution gradient arrays:
  Euler2D_pState        **dWdx; //!< Unlimited solution gradient (x-direction).
  Euler2D_pState        **dWdy; //!< Unlimited solution gradient (y-direction).
  Euler2D_pState     **ddWdxdx; //!< Unlimited solution gradient (x-direction).
  Euler2D_pState     **ddWdxdy; //!< Unlimited solution gradient (x-direction).
  Euler2D_pState     **ddWdydy; //!< Unlimited solution gradient (x-direction).
  Euler2D_pState         **phi; //!< Solution slope limiter.
  //@}

  //@{ @name Boundary solution flux arrays:
  Euler2D_cState        *FluxN, //!< North boundary solution flux.
                        *FluxS, //!< South boundary solution flux.
                        *FluxE, //!< East boundary solution flux.
                        *FluxW; //!< West boundary solution flux.
  //@}

  //@{ @name Problem indicator flags:
  int             Axisymmetric; //!< Axisymmetric flow indicator.
  int           Freeze_Limiter; //!< Limiter freezing indicator.
  static int         Flow_Type; //!< Flow type flag (always inviscid).
  static char*   solutionTitle; //!< Solution title info
  //@}

  //@{ @name Boundary condtion reference states:
  Euler2D_pState          *WoN, //!< Boundary condition reference states for north boundary.
                          *WoS, //!< Boundary condition reference states for south boundary.
                          *WoE, //!< Boundary condition reference states for east boundary.
                          *WoW; //!< Boundary condition reference states for west boundary.

  //! Reference values for north and south boundary conditon reference states
  Euler2D_pState Ref_State_BC_North, Ref_State_BC_South; 
  //! Reference values for east and west boundary conditon reference states
  Euler2D_pState Ref_State_BC_East, Ref_State_BC_West;
  //@}

  //! @name Accuracy assessment data:
  //@{
  static Exact_Solution_Type *ExactSoln;          //!< Pointer to the exact solution
  static Exact_Solution_Type *ExactSolution(void){ return ExactSoln;} //!< Return the exact solution pointer
  Accuracy_Assessment_Type AssessAccuracy;   //!< Variable to provide access to accuracy assessment subroutines
  //@}


  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Euler2D_Quad_Block(void);

  /* Destructor. */
  ~Euler2D_Quad_Block(void){ deallocate(); }
  //@}

  //@{ @name Allocate and deallocate functions.
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);

  //! Allocate memory for high-order variables
  void allocate_HighOrder(const int & NumberOfReconstructions,
			  const vector<int> & ReconstructionOrders);

  //! Allocate memory for high-order boundary conditions
  void allocate_HighOrder_BoundaryConditions(void);

  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);

  //! Deallocate high-order variable memory for structured quadrilateral solution block.
  void deallocate_HighOrder(void);

  //! Deallocate high-order boundary conditions memory.
  void deallocate_HighOrder_BoundaryConditions(void);
  //@}

  //@{ @name Bilinear interplation (Zingg & Yarrow).
  //! Return primitive solution state at specified node.
  Euler2D_pState Wn(const int ii, const int jj);

  //! Return conserverd solution state at specified node.
  Euler2D_cState Un(const int ii, const int jj);

  Euler2D_pState WnNW(const int ii, const int jj); //!< Return primitive solution state at cell NW node.
  Euler2D_pState WnNE(const int ii, const int jj); //!< Return primitive solution state at cell NE node.
  Euler2D_pState WnSE(const int ii, const int jj); //!< Return primitive solution state at cell SE node.
  Euler2D_pState WnSW(const int ii, const int jj); //!< Return primitive solution state at cell SW node.

  Euler2D_cState UnNW(const int ii, const int jj); //!< Return conserved solution state at cell NW node.
  Euler2D_cState UnNE(const int ii, const int jj); //!< Return conserved solution state at cell NE node.
  Euler2D_cState UnSE(const int ii, const int jj); //!< Return conserved solution state at cell SE node.
  Euler2D_cState UnSW(const int ii, const int jj); //!< Return conserved solution state at cell SW node.
  //@}

  //! @name Field access
  //@{
  const Euler2D_pState& CellSolution(const int &ii, const int &jj) const { return W[ii][jj]; }

  //! @name High-order variables
  //@{
  //! Return the array of high-order variable of the current block
  const HighOrderType * HighOrderVariables(void) const { return HO_Ptr; }
  //! Return the high-order variable in the "Pos" position of the current block
  HighOrderType & HighOrderVariable(const unsigned short int & Pos) { return HO_Ptr[Pos]; }
  const HighOrderType & HighOrderVariable(const unsigned short int & Pos) const { return HO_Ptr[Pos]; }
  //! Return the number of high-order variables
  const unsigned short int & NumberOfHighOrderObjects(void) const { return NumberOfHighOrderVariables; }
  //@}

  //! @name High-order boundary conditions (used mostly for constrained reconstruction)
  //@{
  const BC_Type * BC_NorthCell(void) { return HO_WoN;}
  BC_Type & BC_NorthCell(const int &iCell){ return HO_WoN[iCell]; }
  const BC_Type & BC_NorthCell(const int &iCell) const { return HO_WoN[iCell]; }
  const BC_Type * BC_SouthCell(void) { return HO_WoS;}
  BC_Type & BC_SouthCell(const int &iCell){ return HO_WoS[iCell]; }
  const BC_Type & BC_SouthCell(const int &iCell) const { return HO_WoS[iCell]; }
  const BC_Type * BC_EastCell(void) { return HO_WoE;}
  BC_Type & BC_EastCell(const int &jCell){ return HO_WoE[jCell]; }
  const BC_Type & BC_EastCell(const int &jCell) const { return HO_WoE[jCell]; }
  const BC_Type * BC_WestCell(void) { return HO_WoW;}
  BC_Type & BC_WestCell(const int &jCell){ return HO_WoW[jCell]; }
  const BC_Type & BC_WestCell(const int &jCell) const { return HO_WoW[jCell]; }

  void BCs_HighOrder(void);
  //@}

  //@} //end(Field access)

  //! @name Normalization related functions:
  //@{
  //! Get normalization state
  const Euler2D_pState getNormalizationState(const int &ii, const int &jj) const { return RefU; }
  //! Set normalization state
  static void Set_Normalization_Reference_State(const Euler2D_pState & State){ RefU = State; }
  //@}

  //@{ @name Member functions for limiter freezing.
  void evaluate_limiters(void); //!< Set flags for limiter evaluation.
  void freeze_limiters(void);   //!< Set flags for limiter freezing.
  //@}

  //! @name Member functions to compute the piecewise linear solution at a particular location
  //@{
  void SetPiecewiseLinearReconstructionStencil(const int &i, const int &j,
					       int i_index[], int j_index[],
					       int & n_ptr);
  Euler2D_pState PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
						 const double &DeltaXToCentroid,
						 const double &DeltaYToCentroid) const;
  double PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
					 const double &DeltaXToCentroid,
					 const double &DeltaYToCentroid,
					 const unsigned int &parameter) const;
  Euler2D_pState PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
						   const Vector2D &CalculationPoint) const;
  double PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
					   const Vector2D &CalculationPoint,
					   const unsigned int &parameter) const;

  Euler2D_pState UnlimitedPiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
							  const double &DeltaXToCentroid,
							  const double &DeltaYToCentroid) const;
  Euler2D_pState UnlimitedPiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
							    const Vector2D &CalculationPoint) const;
  //@}

  //! @name Member functions to set boundary states
  //@{
  //! Set reference values for boundary reference states
  void Set_Reference_Values_For_Boundary_States(const Euler2D_pState & Ref_North,
						const Euler2D_pState & Ref_South,
						const Euler2D_pState & Ref_East,
						const Euler2D_pState & Ref_West);
  //! Set boundary reference states based on user's input data
  void Set_Boundary_Reference_States_Based_On_Input(const Euler2D_Input_Parameters &IP);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Euler2D_Quad_Block &Soln);
  friend istream &operator >> (istream &in_file, Euler2D_Quad_Block &Soln);
  //@}

  //@{ @name Member functions required for message passing.
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

  //! @name Member functions used for plotting.
  //@{
  void Output_Nodes_Tecplot_HighOrder(const int &Number_of_Time_Steps,
				      const double &Time,
				      const int &Block_Number,
				      const int &Output_Title,
				      ostream &Out_File,
				      const int & StartI_CellIndex,
				      const int & EndI_CellIndex,
				      const int & StartJ_CellIndex,
				      const int & EndJ_CellIndex,
				      const int &IndexHO = 0);

  void Output_Tecplot_HighOrder(const int &Number_of_Time_Steps,
				const double &Time,
				const int &Block_Number,
				const int &Output_Title,
				ostream &Out_File,
				const int &IndexHO = 0);

  void Output_Tecplot_HighOrder_Debug_Mode(AdaptiveBlock2D_List &Soln_Block_List,
					   const Euler2D_Input_Parameters &P,
					   const int &Block_Number,
					   const int &IndexHO = 0);

  void Output_Nodes_Tecplot_HighOrder(const int &Number_of_Time_Steps,
				      const double &Time,
				      const int &Block_Number,
				      const int &Output_Title,
				      ostream &Out_File,
				      const int &IndexHO = 0);

  void Output_Nodes_Tecplot_HighOrder_Debug_Mode(AdaptiveBlock2D_List &Soln_Block_List,
						 const Euler2D_Input_Parameters &P,
						 const int &Block_Number,
						 const int &IndexHO = 0);

  void Output_Cells_Tecplot_HighOrder(const int &Number_of_Time_Steps,
				      const double &Time,
				      const int &Block_Number,
				      const int &Output_Title,
				      ostream &Out_File,
				      const int &IndexHO = 0);

  void Output_Cells_Tecplot_HighOrder_Debug_Mode(AdaptiveBlock2D_List &Soln_Block_List,
						 const Euler2D_Input_Parameters &P,
						 const int &Block_Number,
						 const int &IndexHO = 0);
  //@}


private:
  Euler2D_Quad_Block(const Euler2D_Quad_Block &Soln);  //! Private copy constructor.
  Euler2D_Quad_Block & operator = (const Euler2D_Quad_Block &Soln);   //!< Private assignment operator

  static Euler2D_pState RefU;	//!< reference state for normalizing the solution

  //! @name High-order variables and member functions
  //@{
  HighOrderType* HO_Ptr;  //!< pointer to an array of high-order variables
  unsigned short int NumberOfHighOrderVariables; //!< counter for the total number of high-order variables

  //! Allocate high-order variables array.
  void allocate_HighOrder_Array(const int & NumberOfReconstructions);
  //@}

  //! @name High-order boundary conditions (used mostly for constrained reconstruction)
  //@{
  BC_Type  *HO_WoN, 		//!< High-order boundary condition reference states for North boundary
    *HO_WoS,            	//!< High-order boundary condition reference states for South boundary
    *HO_WoE,			//!< High-order boundary condition reference states for East boundary
    *HO_WoW;            	//!< High-order boundary condition reference states for West boundary
  //@}

};

/****************************************//**
 * Default constructor.
 *****************************************/
inline Euler2D_Quad_Block::Euler2D_Quad_Block(void):
  AssessAccuracy(this),
  Ref_State_BC_North(0.0), Ref_State_BC_South(0.0),
  Ref_State_BC_East(0.0), Ref_State_BC_West(0.0),
  HO_Ptr(NULL), NumberOfHighOrderVariables(0)
{

  // Grid size and variables:
  NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
  // Solution variables:
  W = NULL; U = NULL; dt = NULL; dUdt = NULL; 
  dWdx = NULL; dWdy = NULL; phi = NULL; Uo = NULL;
  FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
  WoN = NULL; WoS = NULL; WoE = NULL; WoW = NULL;
  HO_WoN = NULL; HO_WoS = NULL; HO_WoE = NULL; HO_WoW = NULL;
  Axisymmetric = 0; Freeze_Limiter = OFF;
  
  // Get access to the Euler2D_ExactSolutions object
  ExactSoln = &Euler2D_ExactSolutions::getInstance();

}

/****************************************//**
 * Private copy constructor. (shallow copy)
 *****************************************/
inline Euler2D_Quad_Block::Euler2D_Quad_Block(const Euler2D_Quad_Block &Soln):
  AssessAccuracy(this), HO_Ptr(NULL)
{
  NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu; 
  NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu; Nghost = Soln.Nghost;
  Grid = Soln.Grid; W = Soln.W; U = Soln.U; dt = Soln.dt; dUdt = Soln.dUdt; 
  dWdx = Soln.dWdx; dWdy = Soln.dWdy; phi = Soln.phi; Uo = Soln.Uo;
  FluxN = Soln.FluxN; FluxS = Soln.FluxS; FluxE = Soln.FluxE; FluxW = Soln.FluxW;
  WoN = Soln.WoN; WoS = Soln.WoS; WoE = Soln.WoE; WoW = Soln.WoW;
  HO_WoN = Soln.HO_WoN; HO_WoS = Soln.HO_WoS; HO_WoE = Soln.HO_WoE; HO_WoW = Soln.HO_WoW;
  Axisymmetric = 0; Freeze_Limiter = Soln.Freeze_Limiter;

  Ref_State_BC_North = Soln.Ref_State_BC_North;
  Ref_State_BC_South = Soln.Ref_State_BC_South;
  Ref_State_BC_East = Soln.Ref_State_BC_East;
  Ref_State_BC_West = Soln.Ref_State_BC_West;

  HO_Ptr = Soln.HO_Ptr;
  NumberOfHighOrderVariables = Soln.NumberOfHighOrderVariables;
}

/***********************//**
 * Allocate memory.                       
 ***************************/
inline void Euler2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
   int i, j, k; 
   assert(Ni > 1 && Nj > 1 && Ng > 1);

   // Check to see if the current block dimensions differ from the required ones.
   if ( (Nghost != Ng) || (NCi != Ni+2*Ng) || (NCj != Nj+2*Ng) ){ 

     // free the memory if there is memory allocated
     deallocate();

     // allocate new memory
     Grid.allocate(Ni, Nj, Ng);
     NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1; 
     NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
     W = new Euler2D_pState*[NCi]; U = new Euler2D_cState*[NCi];
     dt = new double*[NCi]; dUdt = new Euler2D_cState**[NCi];
     dWdx = new Euler2D_pState*[NCi]; dWdy = new Euler2D_pState*[NCi];
     phi = new Euler2D_pState*[NCi]; Uo = new Euler2D_cState*[NCi];
     for ( i = 0; i <= NCi-1 ; ++i ) {
       W[i] = new Euler2D_pState[NCj]; U[i] = new Euler2D_cState[NCj];
       dt[i] = new double[NCj]; dUdt[i] = new Euler2D_cState*[NCj];
       for ( j = 0; j <= NCj-1 ; ++j ) 
	 { dUdt[i][j] = new Euler2D_cState[NUMBER_OF_RESIDUAL_VECTORS_EULER2D]; }
       dWdx[i] = new Euler2D_pState[NCj]; dWdy[i] = new Euler2D_pState[NCj];
       phi[i] = new Euler2D_pState[NCj]; Uo[i] = new Euler2D_cState[NCj];
     } /* endfor */
     FluxN = new Euler2D_cState[NCi]; FluxS = new Euler2D_cState[NCi];
     FluxE = new Euler2D_cState[NCj]; FluxW = new Euler2D_cState[NCj];
     WoN = new Euler2D_pState[NCi]; WoS = new Euler2D_pState[NCi];
     WoE = new Euler2D_pState[NCj]; WoW = new Euler2D_pState[NCj];
     // Set the solution residuals, gradients, limiters, and other values to zero.
     for (j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
       for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	 for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_EULER2D-1 ; ++k ) {
	   dUdt[i][j][k] = Euler2D_U_VACUUM;
	 } /* endfor */
	 dWdx[i][j].Vacuum(); dWdy[i][j].Vacuum();
	 phi[i][j].Vacuum(); Uo[i][j].Vacuum();
	 dt[i][j] = ZERO;
       } /* endfor */
     } /* endfor */

   }/* endif */
}

/*****************************************//**
 * Allocate memory for high-order variables.
 ********************************************/
inline void Euler2D_Quad_Block::allocate_HighOrder(const int & NumberOfReconstructions,
						   const vector<int> & ReconstructionOrders){

  bool _pseudo_inverse_allocation_(false);
  int i;

  // Decide whether to allocate the pseudo-inverse
  if (CENO_Execution_Mode::CENO_SPEED_EFFICIENT){
    _pseudo_inverse_allocation_ = true;
  }

  // Re-allocate new memory if necessary
  if (NumberOfReconstructions != NumberOfHighOrderVariables){

    // allocate the high-order array
    allocate_HighOrder_Array(NumberOfReconstructions);
    
    // set the reconstruction order of each high-order object
    for (i=0; i<NumberOfHighOrderVariables; ++i){
      HO_Ptr[i].InitializeVariable(ReconstructionOrders[i],
				   Grid,
				   _pseudo_inverse_allocation_);
    }

  } else {
    // check the reconstruction orders
    for (i=0; i<ReconstructionOrders.size(); ++i){
      if (HighOrderVariable(i).RecOrder() != ReconstructionOrders[i]){
	// change the reconstruction order of the high-order object
	HO_Ptr[i].SetReconstructionOrder(ReconstructionOrders[i]);
      }
    } // endfor
  }// endif
}

/***************************************************//**
 * Allocate memory for high-order boundary conditions
 * if the corresponding boundary reconstruction is 
 * constrained.
 * Assume that Grid and block indexes have been setup!
 *******************************************************/
inline void Euler2D_Quad_Block::allocate_HighOrder_BoundaryConditions(void){

  int i,j;

  // allocate North BCs
  if ( Grid.IsNorthBoundaryReconstructionConstrained() ){

    if (HO_WoN != NULL){
      // deallocate memory
      delete [] HO_WoN; HO_WoN = NULL;
    }

    // allocate new memory
    HO_WoN = new BC_Type[NCi];

    // allocate BC memory for each flux calculation point
    for (i=ICl; i<=ICu; ++i){
      BC_NorthCell(i).allocate(Grid.NumOfConstrainedGaussQuadPoints_North(i,JCu));
    }

  } else if ( HO_WoN != NULL){
    // deallocate memory
    delete [] HO_WoN; HO_WoN = NULL;
  }

  // allocate South BCs
  if ( Grid.IsSouthBoundaryReconstructionConstrained() ){

    if (HO_WoS != NULL){
      // deallocate memory
      delete [] HO_WoS; HO_WoS = NULL;
    }

    // allocate new memory    
    HO_WoS = new BC_Type[NCi];
    
    // allocate BC memory for each flux calculation point
    for (i=ICl; i<=ICu; ++i){
      BC_SouthCell(i).allocate(Grid.NumOfConstrainedGaussQuadPoints_South(i,JCl));
    }    
  } else if (HO_WoS != NULL){
    // deallocate memory
    delete [] HO_WoS; HO_WoS = NULL;
  }

  // allocate East BCs
  if ( Grid.IsEastBoundaryReconstructionConstrained() ){

    if (HO_WoE != NULL){
      // deallocate memory
      delete [] HO_WoE; HO_WoE = NULL;
    }

    // allocate new memory    
    HO_WoE = new BC_Type[NCj];

    // allocate BC memory for each flux calculation point
    for (j=JCl; j<=JCu; ++j){
      BC_EastCell(j).allocate(Grid.NumOfConstrainedGaussQuadPoints_East(ICu,j));
    }
  } else if (HO_WoE != NULL){
    // deallocate memory
    delete [] HO_WoE; HO_WoE = NULL;
  }

  // allocate West BCs
  if ( Grid.IsWestBoundaryReconstructionConstrained() ){

    if (HO_WoW != NULL){
      // deallocate memory
      delete [] HO_WoW; HO_WoW = NULL;
    }

    // allocate new memory    
    HO_WoW = new BC_Type[NCj];

    // allocate BC memory for each flux calculation point
    for (j=JCl; j<=JCu; ++j){
      BC_WestCell(j).allocate(Grid.NumOfConstrainedGaussQuadPoints_West(ICl,j));
    }
  } else if (HO_WoW != NULL){
    // deallocate memory
    delete [] HO_WoW; HO_WoW = NULL;
  }

}

/**************************************************************************
 * Euler2D_Quad_Block::deallocate -- Deallocate memory.                   *
 **************************************************************************/
inline void Euler2D_Quad_Block::deallocate(void) {
  if (U != NULL){
    /* free the memory if there is memory allocated */
    int i, j; Grid.deallocate(); 
    for ( i = 0; i <= NCi-1 ; ++i ) {
      delete []W[i]; W[i] = NULL; delete []U[i]; U[i] = NULL;
      delete []dt[i]; dt[i] = NULL; 
      for ( j = 0; j <= NCj-1 ; ++j ) { delete []dUdt[i][j]; dUdt[i][j] = NULL; }
      delete []dUdt[i]; dUdt[i] = NULL;
      delete []dWdx[i]; dWdx[i] = NULL; delete []dWdy[i]; dWdy[i] = NULL;
      delete []phi[i]; phi[i] = NULL; delete []Uo[i]; Uo[i] = NULL;
    } /* endfor */
    delete []W; W = NULL; delete []U; U = NULL;
    delete []dt; dt = NULL; delete []dUdt; dUdt = NULL;
    delete []dWdx; dWdx = NULL; delete []dWdy; dWdy = NULL; 
    delete []phi; phi = NULL; delete []Uo; Uo = NULL;
    delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
    delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
    delete []WoN; WoN = NULL; delete []WoS; WoS = NULL;
    delete []WoE; WoE = NULL; delete []WoW; WoW = NULL;
    NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;

    deallocate_HighOrder();
    deallocate_HighOrder_BoundaryConditions();   
  }/*endif*/
}

/******************************************//**
 * Deallocate memory for high-order variables
 *********************************************/
inline void Euler2D_Quad_Block::deallocate_HighOrder(void) {
  delete []HO_Ptr; HO_Ptr = NULL;
  NumberOfHighOrderVariables = 0;
}

/****************************************************//**
 * Deallocate memory for high-order boundary conditions
 *******************************************************/
inline void Euler2D_Quad_Block::deallocate_HighOrder_BoundaryConditions(void) {
  if (HO_WoN != NULL){
    delete [] HO_WoN; HO_WoN = NULL;
  }
  if (HO_WoS != NULL){
    delete [] HO_WoS; HO_WoS = NULL;
  }
  if (HO_WoE != NULL){
    delete [] HO_WoE; HO_WoE = NULL;
  }
  if (HO_WoW != NULL){
    delete [] HO_WoW; HO_WoW = NULL;
  }
}

/**************************************************************************
 * Euler2D_Quad_Block::Wn -- Node primitive variable solution state.      *
 **************************************************************************/
inline Euler2D_pState Euler2D_Quad_Block::Wn(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    Euler2D_pState A, B, C, D;
    x=Grid.Node[ii][jj].X.x; y=Grid.Node[ii][jj].X.y;
    ax=Grid.Cell[ii-1][jj-1].Xc.x;
    bx=Grid.Cell[ii-1][jj].Xc.x-Grid.Cell[ii-1][jj-1].Xc.x; 
    cx=Grid.Cell[ii][jj-1].Xc.x-Grid.Cell[ii-1][jj-1].Xc.x; 
    dx=Grid.Cell[ii][jj].Xc.x+Grid.Cell[ii-1][jj-1].Xc.x-
       Grid.Cell[ii-1][jj].Xc.x-Grid.Cell[ii][jj-1].Xc.x;
    ay=Grid.Cell[ii-1][jj-1].Xc.y; 
    by=Grid.Cell[ii-1][jj].Xc.y-Grid.Cell[ii-1][jj-1].Xc.y; 
    cy=Grid.Cell[ii][jj-1].Xc.y-Grid.Cell[ii-1][jj-1].Xc.y; 
    dy=Grid.Cell[ii][jj].Xc.y+Grid.Cell[ii-1][jj-1].Xc.y-
       Grid.Cell[ii-1][jj].Xc.y-Grid.Cell[ii][jj-1].Xc.y;
    aa=bx*dy-dx*by; bb=dy*(ax-x)+bx*cy-cx*by+dx*(y-ay); cc=cy*(ax-x)+cx*(y-ay);
    if (fabs(aa) < TOLER*TOLER) {
      if (fabs(bb) >= TOLER*TOLER) { zeta1=-cc/bb; }
      else { zeta1 = -cc/sgn(bb)*(TOLER*TOLER); } 
      if (fabs(cy+dy*zeta1) >= TOLER*TOLER) { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
      else { eta1 = HALF; } zeta2=zeta1; eta2=eta1;
    } else {
      if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { zeta1=HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; }
      else { zeta1 = -HALF*bb/aa; } 
      if (fabs(cy+dy*zeta1) < TOLER*TOLER) { eta1=-ONE; } 
      else { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
      if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { zeta2=HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa; }
      else { zeta2 = -HALF*bb/aa; }
      if (fabs(cy+dy*zeta2) < TOLER*TOLER) { eta2=-ONE; } 
      else { eta2=(y-ay-by*zeta2)/(cy+dy*zeta2); }
    } /* end if */
    if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
        eta1  > -TOLER && eta1  < ONE + TOLER) {
      zeta=zeta1; eta=eta1;
    } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
               eta2  > -TOLER && eta2  < ONE + TOLER) {
      zeta=zeta2; eta=eta2;
    } else {
      zeta=HALF; eta=HALF;
    } /* endif */
    A=W[ii-1][jj-1]; B=W[ii-1][jj]-W[ii-1][jj-1]; C=W[ii][jj-1]-W[ii-1][jj-1];
    D=W[ii][jj]+W[ii-1][jj-1]-W[ii-1][jj]-W[ii][jj-1];
    return (A+B*zeta+C*eta+D*zeta*eta);
}

/**************************************************************************
 * Euler2D_Quad_Block::Un -- Node conserved variable solution state.      *
 **************************************************************************/
inline Euler2D_cState Euler2D_Quad_Block::Un(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    Euler2D_cState A, B, C, D;
    x=Grid.Node[ii][jj].X.x; y=Grid.Node[ii][jj].X.y;
    ax=Grid.Cell[ii-1][jj-1].Xc.x;
    bx=Grid.Cell[ii-1][jj].Xc.x-Grid.Cell[ii-1][jj-1].Xc.x; 
    cx=Grid.Cell[ii][jj-1].Xc.x-Grid.Cell[ii-1][jj-1].Xc.x; 
    dx=Grid.Cell[ii][jj].Xc.x+Grid.Cell[ii-1][jj-1].Xc.x-
       Grid.Cell[ii-1][jj].Xc.x-Grid.Cell[ii][jj-1].Xc.x;
    ay=Grid.Cell[ii-1][jj-1].Xc.y; 
    by=Grid.Cell[ii-1][jj].Xc.y-Grid.Cell[ii-1][jj-1].Xc.y; 
    cy=Grid.Cell[ii][jj-1].Xc.y-Grid.Cell[ii-1][jj-1].Xc.y; 
    dy=Grid.Cell[ii][jj].Xc.y+Grid.Cell[ii-1][jj-1].Xc.y-
       Grid.Cell[ii-1][jj].Xc.y-Grid.Cell[ii][jj-1].Xc.y;
    aa=bx*dy-dx*by; bb=dy*(ax-x)+bx*cy-cx*by+dx*(y-ay); cc=cy*(ax-x)+cx*(y-ay);
    if (fabs(aa) < TOLER*TOLER) {
      if (fabs(bb) >= TOLER*TOLER) { zeta1=-cc/bb; }
      else { zeta1 = -cc/sgn(bb)*(TOLER*TOLER); } 
      if (fabs(cy+dy*zeta1) >= TOLER*TOLER) { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
      else { eta1 = HALF; } zeta2=zeta1; eta2=eta1;
    } else {
      if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { zeta1=HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; }
      else { zeta1 = -HALF*bb/aa; } 
      if (fabs(cy+dy*zeta1) < TOLER*TOLER) { eta1=-ONE; } 
      else { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
      if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { zeta2=HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa; }
      else { zeta2 = -HALF*bb/aa; }
      if (fabs(cy+dy*zeta2) < TOLER*TOLER) { eta2=-ONE; } 
      else { eta2=(y-ay-by*zeta2)/(cy+dy*zeta2); }
    } /* end if */
    if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
        eta1  > -TOLER && eta1  < ONE + TOLER) {
      zeta=zeta1; eta=eta1;
    } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
               eta2  > -TOLER && eta2  < ONE + TOLER) {
      zeta=zeta2; eta=eta2;
    } else {
      zeta=HALF; eta=HALF;
    } /* endif */
    A=U[ii-1][jj-1]; B=U[ii-1][jj]-U[ii-1][jj-1]; C=U[ii][jj-1]-U[ii-1][jj-1];
    D=U[ii][jj]+U[ii-1][jj-1]-U[ii-1][jj]-U[ii][jj-1];
    return (A+B*zeta+C*eta+D*zeta*eta);
}

/**************************************************************************
 * Euler2D_Quad_Block::Wn?? -- Get cell node primitive solution states.   *
 **************************************************************************/
inline Euler2D_pState Euler2D_Quad_Block::WnNW(const int ii, const int jj) {
  return (Wn(ii, jj+1));
}

inline Euler2D_pState Euler2D_Quad_Block::WnNE(const int ii, const int jj) {
  return (Wn(ii+1, jj+1));
}

inline Euler2D_pState Euler2D_Quad_Block::WnSE(const int ii, const int jj) {
  return (Wn(ii+1, jj));
}

inline Euler2D_pState Euler2D_Quad_Block::WnSW(const int ii, const int jj) {
  return (Wn(ii, jj));
}

/**************************************************************************
 * Euler2D_Quad_Block::Un?? -- Get cell node conserved solution states.   *
 **************************************************************************/
inline Euler2D_cState Euler2D_Quad_Block::UnNW(const int ii, const int jj) {
  return (Un(ii, jj+1));
}

inline Euler2D_cState Euler2D_Quad_Block::UnNE(const int ii, const int jj) {
  return (Un(ii+1, jj+1));
}

inline Euler2D_cState Euler2D_Quad_Block::UnSE(const int ii, const int jj) {
  return (Un(ii+1, jj));
}

inline Euler2D_cState Euler2D_Quad_Block::UnSW(const int ii, const int jj) {
  return (Un(ii, jj));
}

/**************************************************************************
 * Euler2D_Quad_Block::evaluate_limiters -- Set flag to evaluate limiters.*
 **************************************************************************/
inline void Euler2D_Quad_Block::evaluate_limiters(void) {
  Freeze_Limiter = OFF; 
}

/**************************************************************************
 * Euler2D_Quad_Block::freeze_limiters -- Set flag to freeze limiters.    *
 **************************************************************************/
inline void Euler2D_Quad_Block::freeze_limiters(void) {
  Freeze_Limiter = ON; 
}

/*!
 * Determine the number of neighbouring cells to
 * be used in the piecewise linear reconstruction
 * procedure and set the cell indexes in the 
 * provided containers.
 *
 * \param [in] i the i-index of the reconstructed cell
 * \param [in] j the j-index of the reconstructed cell
 * \param [out] i_index the "i" indexes of the neighbouring cells that form the stencil
 * \param [out] j_index the "j" indexes of the neighbouring cells that form the stencil
 * \param [out] n_ptr the number of neighbouring cells that form the stencil.
 * This number is typically 8, but it can vary for different boundary conditions.
 */
inline void Euler2D_Quad_Block::SetPiecewiseLinearReconstructionStencil(const int &i, const int &j,
									int i_index[], int j_index[],
									int &n_pts){

  if (i <= ICl-Nghost || i >= ICu+Nghost || j <= JCl-Nghost || j >= JCu+Nghost) {
    n_pts = 0;
  } else if ((i == ICl-1) && (Grid.BCtypeW[j] != BC_NONE)) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeW[j] == BC_PERIODIC ||
	       Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	       Grid.BCtypeW[j] == BC_NEUMANN ||
	       Grid.BCtypeW[j] == BC_ROBIN) {
      if (j == JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (j == JCl) {
	n_pts = 3;
	i_index[0] = i+1; j_index[0] = j  ;
	i_index[1] = i  ; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j+1;
      } else if (j == JCu) {
	n_pts = 3;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } /* endif */
    } /* endif */           
  } else if ((i == ICu+1) && (Grid.BCtypeE[j] != BC_NONE)) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeE[j] == BC_PERIODIC ||
	       Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeE[j] == BC_CHARACTERISTIC ||	       
	       Grid.BCtypeE[j] == BC_NEUMANN ||
	       Grid.BCtypeE[j] == BC_ROBIN) {
      if (j == JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (j == JCl) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i  ; j_index[2] = j+1;
      } else if (j == JCu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } /* endif */
    } /* endif */
  } else if ((j == JCl-1) && (Grid.BCtypeS[i] != BC_NONE)) {
    if (i == ICl-1 || i == ICu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeS[i] == BC_PERIODIC ||
	       Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeS[i] == BC_CHARACTERISTIC ||	       
	       Grid.BCtypeS[i] == BC_NEUMANN ||
	       Grid.BCtypeS[i] == BC_ROBIN) {
      if (i == ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (i == ICl) {
	n_pts = 3;
	i_index[0] = i+1; j_index[0] = j  ;
	i_index[1] = i  ; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j+1;
      } else if (i == ICu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i  ; j_index[2] = j+1;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } /* endif */
    } /* endif */
  } else if ((j == JCu+1) && (Grid.BCtypeN[i] != BC_NONE)) {
    if (i == ICl-1 || i == ICu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeN[i] == BC_PERIODIC ||
	       Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeN[i] == BC_CHARACTERISTIC ||	       
	       Grid.BCtypeN[i] == BC_NEUMANN ||
	       Grid.BCtypeN[i] == BC_ROBIN) {
      if (i == ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (i == ICl) {
	n_pts = 3;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
      } else if (i == ICu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } /* endif */
    } /* endif */
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  } /* endif */
}

/**************************************************************************
 * Euler2D_Quad_Block -- Input-output operators.                          *
 **************************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Euler2D_Quad_Block &SolnBlk) {
  int i, j; 
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " 
	   << SolnBlk.ICu << " " << SolnBlk.Nghost << "\n";
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << "\n";
  out_file << SolnBlk.Axisymmetric << "\n";
  if (SolnBlk.NCi == 0 || SolnBlk.NCj == 0) return(out_file);
  for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
     for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
         out_file << SolnBlk.U[i][j] << "\n";
     } /* endfor */
  } /* endfor */
  for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
     out_file << SolnBlk.WoW[j] << "\n";
     out_file << SolnBlk.WoE[j] << "\n";
  } /* endfor */
  for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
     out_file << SolnBlk.WoS[i] << "\n";
     out_file << SolnBlk.WoN[i] << "\n";
  } /* endfor */

  // Output the high-order variables
  out_file << SolnBlk.NumberOfHighOrderVariables << "\n";
  for (int k = 1; k <= SolnBlk.NumberOfHighOrderVariables; ++k){
    out_file << SolnBlk.HighOrderVariable(k-1);
  }/* endfor */

  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Euler2D_Quad_Block &SolnBlk) {
  int i, j, k, ni, il, iu, ng, nj, jl, ju, n_HO;
  Grid2D_Quad_Block_HO New_Grid;
  in_file >> New_Grid; 
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu >> ng;
  in_file >> nj >> jl >> ju;
  in_file >> SolnBlk.Axisymmetric;
  in_file.unsetf(ios::skipws);
  if (ni == 0 || nj == 0) {
      SolnBlk.deallocate(); return(in_file);
  } /* endif */
  if (SolnBlk.U == NULL || SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) {
      if (SolnBlk.U != NULL) SolnBlk.deallocate(); 
      SolnBlk.allocate(ni - 2*ng, nj - 2*ng, ng);
  } /* endif */

  // Copy the temporary mesh into the grid of the current solution block
  SolnBlk.Grid = New_Grid;

  // Read the solution & Initialize some data structures
  for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
     for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
         in_file >> SolnBlk.U[i][j];
         SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
         for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_EULER2D-1 ; ++k ) {
	   SolnBlk.dUdt[i][j][k].Vacuum();
         } /* endfor */
	 SolnBlk.dWdx[i][j].Vacuum();
	 SolnBlk.dWdy[i][j].Vacuum();
	 SolnBlk.phi[i][j].Vacuum();
	 SolnBlk.Uo[i][j].Vacuum();
	 SolnBlk.dt[i][j] = ZERO;
     } /* endfor */
  } /* endfor */
  for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
     in_file >> SolnBlk.WoW[j];
     in_file >> SolnBlk.WoE[j];
  } /* endfor */
  for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
     in_file >> SolnBlk.WoS[i];
     in_file >> SolnBlk.WoN[i];
  } /* endfor */

  // Allocate memory for the high-order boundary conditions (the values don't need to be read)
  SolnBlk.allocate_HighOrder_BoundaryConditions();

  in_file.setf(ios::skipws);

  // Read the high-order variables
  in_file >> n_HO;
  SolnBlk.allocate_HighOrder_Array(n_HO);

  for (int k = 1; k <= SolnBlk.NumberOfHighOrderVariables; ++k){
    // Read the current variable
    in_file >> SolnBlk.HighOrderVariable(k-1);
    // Associate this variable to the current grid
    SolnBlk.HighOrderVariable(k-1).AssociateGeometry(SolnBlk.Grid);
  }/* endfor */

  in_file.setf(ios::skipws);

  return (in_file);
}

/*******************************************************************************
 * Euler2D_Quad_Block::NumVar -- Returns number of state variables.            *
 *******************************************************************************/
inline int Euler2D_Quad_Block::NumVar(void) {
  return (int(NUM_VAR_EULER2D));
}

/*!
 * Compute the solution provided by the piecewise linear reconstruction
 * at a location relative to the cell centroid.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param DeltaXToCentroid the X-delta between the point of interest and the centroid
 * \param DeltaYToCentroid the Y-delta between the point of interest and the centroid
 */
inline Euler2D_pState Euler2D_Quad_Block::PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
									  const double &DeltaXToCentroid,
									  const double &DeltaYToCentroid) const{
  return W[ii][jj] + (phi[ii][jj]^(dWdx[ii][jj]*DeltaXToCentroid + dWdy[ii][jj]*DeltaYToCentroid) );
}

/*!
 * Compute the solution provided by the piecewise linear reconstruction
 * for a particular state parameter at a location relative to the cell centroid.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param DeltaXToCentroid the X-delta between the point of interest and the centroid
 * \param DeltaYToCentroid the Y-delta between the point of interest and the centroid
 * \param parameter the parameter for which the solution is computed
 */
inline double Euler2D_Quad_Block::PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
								  const double &DeltaXToCentroid,
								  const double &DeltaYToCentroid,
								  const unsigned int &parameter) const{
  return W[ii][jj][parameter] + (phi[ii][jj][parameter]*(dWdx[ii][jj][parameter]*DeltaXToCentroid + 
							 dWdy[ii][jj][parameter]*DeltaYToCentroid) );
}

/*!
 * Compute the solution provided by the piecewise linear reconstruction
 * at a particular location.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param CalculationPoint the Cartesian coordinates of the point of interest
 */
inline Euler2D_pState Euler2D_Quad_Block::PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
									    const Vector2D &CalculationPoint) const{
  // calculate the distance between the point of interest and the centroid of the cell
  Vector2D dX(CalculationPoint - Grid.Cell[ii][jj].Xc);	
  return PiecewiseLinearSolutionForDelta(ii,jj,dX.x,dX.y);
}

/*!
 * Compute the solution provided by the piecewise linear reconstruction
 * for a particular state parameter at a particular location.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param CalculationPoint the Cartesian coordinates of the point of interest
 * \param parameter the parameter for which the solution is computed
 */
inline double Euler2D_Quad_Block::PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
								    const Vector2D &CalculationPoint,
								    const unsigned int &parameter) const{
  // calculate the distance between the point of interest and the centroid of the cell
  Vector2D dX(CalculationPoint - Grid.Cell[ii][jj].Xc);	  
  return PiecewiseLinearSolutionForDelta(ii,jj,dX.x,dX.y,parameter);
}

/*!
 * Compute the solution provided by the UNLIMITED (i.e. no limiter!) piecewise 
 * linear reconstruction at a location relative to the cell centroid.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param DeltaXToCentroid the X-delta between the point of interest and the centroid
 * \param DeltaYToCentroid the Y-delta between the point of interest and the centroid
 */
inline Euler2D_pState Euler2D_Quad_Block::UnlimitedPiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
										   const double &DeltaXToCentroid,
										   const double &DeltaYToCentroid) const{

  return W[ii][jj] + dWdx[ii][jj]*DeltaXToCentroid + dWdy[ii][jj]*DeltaYToCentroid;
}

/*!
 * Compute the solution provided by the UNLIMITED (i.e. no limiter!) piecewise 
 * linear reconstruction at a particular location.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param CalculationPoint the Cartesian coordinates of the point of interest
 */
inline Euler2D_pState Euler2D_Quad_Block::UnlimitedPiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
										     const Vector2D &CalculationPoint) const{

  // calculate the distance between the point of interest and the centroid of the cell
  Vector2D dX(CalculationPoint - Grid.Cell[ii][jj].Xc);	
  return UnlimitedPiecewiseLinearSolutionForDelta(ii,jj,dX.x,dX.y);
}

/*!
 * Allocate the high-order array.
 */
inline void Euler2D_Quad_Block::allocate_HighOrder_Array(const int & NumberOfReconstructions){

  if (NumberOfReconstructions != NumberOfHighOrderVariables){
    
    // deallocate memory
    deallocate_HighOrder();
    
    // allocate new memory for high-order objects
    HO_Ptr = new HighOrderType [NumberOfReconstructions];
    
    // check for successful memory allocation
    if (HO_Ptr != NULL){
      NumberOfHighOrderVariables = NumberOfReconstructions;
    }
  }
}

/*******************************************************************************
 *                                                                             *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                              *
 *                                                                             *
 *******************************************************************************/

/*******************************************************************************
 * Euler2D_Quad_Block::LoadSendBuffer -- Loads send message buffer.            *
 *******************************************************************************/
inline int Euler2D_Quad_Block::LoadSendBuffer(double *buffer,
                                              int &buffer_count,
                                              const int buffer_size,
                                              const int i_min, 
                                              const int i_max,
                                              const int i_inc,
                                              const int j_min, 
                                              const int j_max,
                                              const int j_inc) {
  int i, j, k;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = U[i][j][k];
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Euler2D_Quad_Block::LoadSendBuffer_F2C -- Loads send message buffer for     *
 *                                           fine to coarse block message      *
 *                                           passing.                          *
 *******************************************************************************/
inline int Euler2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
                                                  int &buffer_count,
                                                  const int buffer_size,
                                                  const int i_min, 
                                                  const int i_max,
                                                  const int i_inc,
                                                  const int j_min, 
                                                  const int j_max,
                                                  const int j_inc) {
  int i, j, k;
  for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*U[i  ][j  ][k]+
                                   Grid.Cell[i+1][j  ].A*U[i+1][j  ][k]+
                                   Grid.Cell[i  ][j+1].A*U[i  ][j+1][k]+
                                   Grid.Cell[i+1][j+1].A*U[i+1][j+1][k])/
                                  (Grid.Cell[i  ][j  ].A+
                                   Grid.Cell[i+1][j  ].A+
                                   Grid.Cell[i  ][j+1].A+
                                   Grid.Cell[i+1][j+1].A);
/*            buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*U[i  ][j  ][k]+ */
/*                                    Grid.Cell[i+1][j  ].A*U[i+1][j  ][k]+ */
/*                                    Grid.Cell[i  ][j+1].A*U[i  ][j+1][k]+ */
/*                                    Grid.Cell[i+1][j+1].A*U[i+1][j+1][k]); */
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Euler2D_Quad_Block::LoadSendBuffer_C2F -- Loads send message buffer for     *
 *                                           coarse to fine block message      *
 *                                           passing.                          *
 *******************************************************************************/
inline int Euler2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
                                                  int &buffer_count,
                                                  const int buffer_size,
                                                  const int i_min, 
                                                  const int i_max,
                                                  const int i_inc,
                                                  const int j_min, 
                                                  const int j_max,
                                                  const int j_inc,
						  const int face,
						  const int sector) {
  int i, j, k;
  Vector2D dX;
  Euler2D_pState Wfine;
  Euler2D_cState Ufine;

  if (j_inc > 0) {
    if (i_inc > 0) {
      for ( j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j, LIMITER_VENKATAKRISHNAN);
	  // Evaluate SW sub (fine) cell values if required.
	  if (!(face == NORTH && sector == WEST && Nghost%2 && j == j_min) &&
	      !(face == NORTH && sector == EAST && Nghost%2 && (i == i_min || j == j_min)) &&
	      !(face == SOUTH && sector == EAST && Nghost%2 && i == i_min) &&
	      !(face == EAST && sector == NORTH && Nghost%2 && (i == i_min || j == j_min)) &&
	      !(face == EAST && sector == SOUTH && Nghost%2 && i == i_min) &&
	      !(face == WEST && sector == NORTH && Nghost%2 && j == j_min) &&
	      !(face == NORTH_EAST && Nghost%2 && (i == i_min || j == j_min)) &&
	      !(face == NORTH_WEST && Nghost%2 && j == j_min) &&
	      !(face == SOUTH_EAST && Nghost%2 && i == i_min)) {
	    dX = Grid.centroidSW(i,j) - Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
	  } /* endif */
	  // Evaluate SE sub (fine) cell values if required.
	  if (!(face == NORTH && sector == WEST && Nghost%2 && (i == i_max || j == j_min)) &&
	      !(face == NORTH && sector == EAST && Nghost%2 && j == j_min) &&
	      !(face == SOUTH && sector == WEST && Nghost%2 && i == i_max) &&
	      !(face == EAST && sector == NORTH && Nghost%2 && j == j_min) &&
	      !(face == WEST && sector == NORTH && Nghost%2 && (i == i_max || j == j_min)) &&
	      !(face == WEST && sector == SOUTH && Nghost%2 && i == i_max) &&
	      !(face == NORTH_EAST && Nghost%2 && j == j_min) &&
	      !(face == NORTH_WEST && Nghost%2 && (i == i_max || j == j_min)) &&
	      !(face == SOUTH_WEST && Nghost%2 && i == i_max)) {
	    dX = Grid.centroidSE(i,j) - Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
	  } /* endif */
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate NW sub (fine) cell values if required.
	  if (!(face == NORTH && sector == EAST && Nghost%2 && i == i_min) &&
	      !(face == SOUTH && sector == EAST && Nghost%2 && (i == i_min || j == j_max)) &&
	      !(face == SOUTH && sector == WEST && Nghost%2 && j == j_max) &&
	      !(face == EAST && sector == NORTH && Nghost%2 && i == i_min) &&
	      !(face == EAST && sector == SOUTH && Nghost%2 && (i == i_min || j == j_max)) &&
	      !(face == WEST && sector == SOUTH && Nghost%2 && j == j_max) &&
	      !(face == NORTH_EAST && Nghost%2 && i == i_min) &&
	      !(face == SOUTH_EAST && Nghost%2 && (i == i_min || j == j_max)) &&
	      !(face == SOUTH_WEST && Nghost%2 && j == j_max)) {
	    dX = Grid.centroidNW(i,j) - Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
	  } /* endif */
	  // Evaluate NE sub (fine) cell values if required.
	  if (!(face == NORTH && sector == WEST && Nghost%2 && i == i_max) &&
	      !(face == SOUTH && sector == EAST && Nghost%2 && j == j_max) &&
	      !(face == SOUTH && sector == WEST && Nghost%2 && (i == i_max || j == j_max)) &&
	      !(face == EAST && sector == SOUTH && Nghost%2 && j == j_max) &&
	      !(face == WEST && sector == NORTH && Nghost%2 && i == i_max) &&
	      !(face == WEST && sector == SOUTH && Nghost%2 && (i == i_max || j == j_max)) &&
	      !(face == NORTH_WEST && Nghost%2 && i == i_max) &&
	      !(face == SOUTH_EAST && Nghost%2 && j == j_max) &&
	      !(face == SOUTH_WEST && Nghost%2 && (i == i_max || j == j_max))) {
	    dX = Grid.centroidNE(i,j) - Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
	  } /* endif */
	} /* endfor */
      } /* endfor */

      return 0;

    } /* endif */
  } /* endif */

  // Load send message buffer for the coarse-to-fine grid for cases in
  // which one (or both) of the increments is negative.  Only for two
  // ghost cells.

  if (j_min == j_max) { // North or south boundary.
     // Four different orderings to consider depending on the value of i_inc & j_inc.
     if (j_inc > 0) {
        if (i_inc > 0) {
 	   return 1;
        } else {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    Grid.Node[i+1][j_min].X+
                    Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i][j_min].X+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
                    Grid.Node[i+1][j_min+1].X)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc+
                    Grid.Node[i][j_min+1].X+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } else {
        if (i_inc > 0) {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc+
                    Grid.Node[i][j_min+1].X+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
                    Grid.Node[i+1][j_min+1].X)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i][j_min].X+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    Grid.Node[i+1][j_min].X+
                    Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } else {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
                    Grid.Node[i+1][j_min+1].X)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc+
                    Grid.Node[i][j_min+1].X+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    Grid.Node[i+1][j_min].X+
                    Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i][j_min].X+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */
  } else { // East or west boundary.
     // Four different orderings to consider depending on the value of i_inc & j_inc.
     if (j_inc > 0) {
        if (i_inc > 0) {
 	   return 1;
        } else {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    Grid.Node[i_min+1][j].X+
                    Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i_min][j].X+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
             dX = (Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
                    Grid.Node[i_min+1][j+1].X)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc+
                    Grid.Node[i_min][j+1].X+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } else {
        if (i_inc > 0) {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc+
                    Grid.Node[i_min][j+1].X+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
                    Grid.Node[i_min+1][j+1].X)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i_min][j].X+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    Grid.Node[i_min+1][j].X+
                    Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } else {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
                    Grid.Node[i_min+1][j+1].X)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc+
                    Grid.Node[i_min][j+1].X+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    Grid.Node[i_min+1][j].X+
                    Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i_min][j].X+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */
  } /* endif */
  return(0);
}

/*******************************************************************************
 * Euler2D_Quad_Block::UnloadReceiveBuffer -- Unloads receive message buffer.  *
 *******************************************************************************/
inline int Euler2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
                                                   int &buffer_count,
                                                   const int buffer_size,
                                                   const int i_min, 
                                                   const int i_max,
                                                   const int i_inc,
                                                   const int j_min, 
                                                   const int j_max,
                                                   const int j_inc) {
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
  	buffer_count = buffer_count + NUM_VAR_EULER2D;
        if (buffer_count >= buffer_size) return(1);
        U[i][j] = Euler2D_cState(buffer[buffer_count-3],
                                 buffer[buffer_count-2],
                                 buffer[buffer_count-1],
                                 buffer[buffer_count]);
        W[i][j] = U[i][j].W();
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Euler2D_Quad_Block::UnloadReceiveBuffer_F2C -- Unloads receive message      *
 *                                                buffer for fine to coarse    *
 *                                                block message passing.       *
 *******************************************************************************/
inline int Euler2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) {
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
  	buffer_count = buffer_count + NUM_VAR_EULER2D;
        if (buffer_count >= buffer_size) return(1);
        U[i][j] = Euler2D_cState(buffer[buffer_count-3],
                                 buffer[buffer_count-2],
                                 buffer[buffer_count-1],
                                 buffer[buffer_count]);
/*         U[i][j] = Euler2D_cState(buffer[buffer_count-3], */
/*                                  buffer[buffer_count-2], */
/*                                  buffer[buffer_count-1], */
/*                                  buffer[buffer_count])/Grid.Cell[i][j].A; */
        W[i][j] = U[i][j].W();
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Euler2D_Quad_Block::UnloadReceiveBuffer_C2F -- Unloads receive message      *
 *                                                buffer for coarse to fine    *
 *                                                block message passing.       *
 *******************************************************************************/
inline int Euler2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) {
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
  	buffer_count = buffer_count + NUM_VAR_EULER2D;
        if (buffer_count >= buffer_size) return(1);
        U[i][j] = Euler2D_cState(buffer[buffer_count-3],
                                 buffer[buffer_count-2],
                                 buffer[buffer_count-1],
                                 buffer[buffer_count]);
        W[i][j] = U[i][j].W();
     } /* endfor */
  } /* endfor */
  return(0);
}

/**************************************************************************
 * Euler2D_Quad_Block::SubcellReconstruction --                           *
 *               Performs the subcell reconstruction of solution state    *
 *               within a given cell (i,j) of the computational mesh for  *
 *               the specified quadrilateral solution block.              *
 **************************************************************************/
inline void Euler2D_Quad_Block::SubcellReconstruction(const int i, 
                                                      const int j,
                                                      const int Limiter) {

  int n, n2, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Euler2D_pState DU, DUDx_ave, DUDy_ave;

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  if (i == ICl-2 || i == ICu+2 ||
      j == JCl-2 || j == JCu+2) {
    n_pts = 0;
  } else if ((i == ICl-1) && 
             (Grid.BCtypeW[j] != BC_NONE)) {
    if (j == JCl-1 || j == JCu+1) {
       n_pts = 0;
    } else if (Grid.BCtypeW[j] == BC_PERIODIC ||
               Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
       if (j == JCl) {
          n_pts = 5;
          i_index[0] = i-1; j_index[0] = j  ;
          i_index[1] = i+1; j_index[1] = j  ;
          i_index[2] = i-1; j_index[2] = j+1;
          i_index[3] = i  ; j_index[3] = j+1;
          i_index[4] = i+1; j_index[4] = j+1;
       } else if (j == JCu) {
          n_pts = 5;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j-1;
          i_index[3] = i-1; j_index[3] = j  ;
          i_index[4] = i+1; j_index[4] = j  ;
       } else {
          n_pts = 8;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j-1;
          i_index[3] = i-1; j_index[3] = j  ;
          i_index[4] = i+1; j_index[4] = j  ;
          i_index[5] = i-1; j_index[5] = j+1;
          i_index[6] = i  ; j_index[6] = j+1;
          i_index[7] = i+1; j_index[7] = j+1;
       } /* endif */
    } else {
       if (j == JCl) {
          n_pts = 3;
          i_index[0] = i+1; j_index[0] = j  ;
          i_index[1] = i  ; j_index[1] = j+1;
          i_index[2] = i+1; j_index[2] = j+1;
       } else if (j == JCu) {
          n_pts = 3;
          i_index[0] = i  ; j_index[0] = j-1;
          i_index[1] = i+1; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j  ;
       } else {
          n_pts = 5;
          i_index[0] = i  ; j_index[0] = j-1;
          i_index[1] = i+1; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j  ;
          i_index[3] = i  ; j_index[3] = j+1;
          i_index[4] = i+1; j_index[4] = j+1;
       } /* endif */
    } /* endif */           
  } else if ((i == ICu+1) && 
             (Grid.BCtypeE[j] != BC_NONE)) {
    if (j == JCl-1 || j == JCu+1) {
       n_pts = 0;
    } else if (Grid.BCtypeE[j] == BC_PERIODIC ||
               Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
       if (j == JCl) {
          n_pts = 5;
          i_index[0] = i-1; j_index[0] = j  ;
          i_index[1] = i+1; j_index[1] = j  ;
          i_index[2] = i-1; j_index[2] = j+1;
          i_index[3] = i  ; j_index[3] = j+1;
          i_index[4] = i+1; j_index[4] = j+1;
       } else if (j == JCu) {
          n_pts = 5;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j-1;
          i_index[3] = i-1; j_index[3] = j  ;
          i_index[4] = i+1; j_index[4] = j  ;
       } else {
          n_pts = 8;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j-1;
          i_index[3] = i-1; j_index[3] = j  ;
          i_index[4] = i+1; j_index[4] = j  ;
          i_index[5] = i-1; j_index[5] = j+1;
          i_index[6] = i  ; j_index[6] = j+1;
          i_index[7] = i+1; j_index[7] = j+1;
       } /* endif */
    } else {
       if (j == JCl) {
          n_pts = 3;
          i_index[0] = i-1; j_index[0] = j  ;
          i_index[1] = i-1; j_index[1] = j+1;
          i_index[2] = i  ; j_index[2] = j+1;
       } else if (j == JCu) {
          n_pts = 3;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i-1; j_index[2] = j  ;
       } else {
          n_pts = 5;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i-1; j_index[2] = j  ;
          i_index[3] = i-1; j_index[3] = j+1;
          i_index[4] = i  ; j_index[4] = j+1;
       } /* endif */
    } /* endif */
  } else if ((j == JCl-1) && 
             (Grid.BCtypeS[i] != BC_NONE)) {
    if (i == ICl-1 || i == ICu+1) {
       n_pts = 0;
    } else if (Grid.BCtypeS[i] == BC_PERIODIC ||
               Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
       if (i == ICl) {
          n_pts = 5;
          i_index[0] = i  ; j_index[0] = j-1;
          i_index[1] = i+1; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j  ;
          i_index[3] = i  ; j_index[3] = j+1;
          i_index[4] = i+1; j_index[4] = j+1;
       } else if (i == ICu) {
          n_pts = 5;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i-1; j_index[2] = j  ;
          i_index[3] = i-1; j_index[3] = j+1;
          i_index[4] = i  ; j_index[4] = j+1;
       } else {
          n_pts = 8;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j-1;
          i_index[3] = i-1; j_index[3] = j  ;
          i_index[4] = i+1; j_index[4] = j  ;
          i_index[5] = i-1; j_index[5] = j+1;
          i_index[6] = i  ; j_index[6] = j+1;
          i_index[7] = i+1; j_index[7] = j+1;
       } /* endif */
    } else {
       if (i == ICl) {
          n_pts = 3;
          i_index[0] = i+1; j_index[0] = j  ;
          i_index[1] = i  ; j_index[1] = j+1;
          i_index[2] = i+1; j_index[2] = j+1;
       } else if (i == ICu) {
          n_pts = 3;
          i_index[0] = i-1; j_index[0] = j  ;
          i_index[1] = i-1; j_index[1] = j+1;
          i_index[2] = i  ; j_index[2] = j+1;
       } else {
          n_pts = 5;
          i_index[0] = i-1; j_index[0] = j  ;
          i_index[1] = i+1; j_index[1] = j  ;
          i_index[2] = i-1; j_index[2] = j+1;
          i_index[3] = i  ; j_index[3] = j+1;
          i_index[4] = i+1; j_index[4] = j+1;
       } /* endif */
    } /* endif */
  } else if ((j == JCu+1) && 
             (Grid.BCtypeN[i] != BC_NONE)) {
    if (i == ICl-1 || i == ICu+1) {
       n_pts = 0;
    } else if (Grid.BCtypeN[i] == BC_PERIODIC ||
               Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
       if (i == ICl) {
          n_pts = 5;
          i_index[0] = i  ; j_index[0] = j-1;
          i_index[1] = i+1; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j  ;
          i_index[3] = i  ; j_index[3] = j+1;
          i_index[4] = i+1; j_index[4] = j+1;
       } else if (i == ICu) {
          n_pts = 5;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i-1; j_index[2] = j  ;
          i_index[3] = i-1; j_index[3] = j+1;
          i_index[4] = i  ; j_index[4] = j+1;
       } else {
          n_pts = 8;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j-1;
          i_index[3] = i-1; j_index[3] = j  ;
          i_index[4] = i+1; j_index[4] = j  ;
          i_index[5] = i-1; j_index[5] = j+1;
          i_index[6] = i  ; j_index[6] = j+1;
          i_index[7] = i+1; j_index[7] = j+1;
       } /* endif */
    } else {
       if (i == ICl) {
          n_pts = 3;
          i_index[0] = i  ; j_index[0] = j-1;
          i_index[1] = i+1; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j  ;
       } else if (i == ICu) {
          n_pts = 3;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i-1; j_index[2] = j  ;
       } else {
          n_pts = 5;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j-1;
          i_index[3] = i-1; j_index[3] = j  ;
          i_index[4] = i+1; j_index[4] = j  ;
       } /* endif */
    } /* endif */
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  } /* endif */
  
  if (n_pts > 0) {
      DUDx_ave = Euler2D_W_VACUUM;
      DUDy_ave = Euler2D_W_VACUUM;
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;
  
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
          dX = Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
               Grid.Cell[i][j].Xc;
          DU = W[ i_index[n2] ][ j_index[n2] ] - 
               W[i][j];
          DUDx_ave += DU*dX.x;
          DUDy_ave += DU*dX.y;
          DxDx_ave += dX.x*dX.x;
          DxDy_ave += dX.x*dX.y;
          DyDy_ave += dX.y*dX.y;
      } /* endfor */
  					    
      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);
      dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
      // Calculate slope limiters. 
      if (!Freeze_Limiter) {
         for ( n = 1 ; n <= NUM_VAR_EULER2D ; ++n ) {
            u0Min = W[i][j][n];
            u0Max = u0Min;
            for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
               u0Min = min(u0Min, W[ i_index[n2] ][ j_index[n2] ][n]);
               u0Max = max(u0Max, W[ i_index[n2] ][ j_index[n2] ][n]);
            } /* endfor */
    
            dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
            uQuad[0] = W[i][j][n] + 
                       dWdx[i][j][n]*dX.x +
                       dWdy[i][j][n]*dX.y ;
            dX = Grid.xfaceW(i, j)-Grid.Cell[i][j].Xc;
            uQuad[1] = W[i][j][n] + 
                       dWdx[i][j][n]*dX.x +
                       dWdy[i][j][n]*dX.y ;
            dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
            uQuad[2] = W[i][j][n] + 
                       dWdx[i][j][n]*dX.x +
                       dWdy[i][j][n]*dX.y ;
            dX = Grid.xfaceS(i, j)-Grid.Cell[i][j].Xc;
            uQuad[3] = W[i][j][n] + 
                       dWdx[i][j][n]*dX.x +
                       dWdy[i][j][n]*dX.y ;
    
            switch(Limiter) {
              case LIMITER_ONE :
                phi_n = ONE;
                break;
              case LIMITER_ZERO :
                phi_n = ZERO;
                break;
              case LIMITER_BARTH_JESPERSEN :
                phi_n = Limiter_BarthJespersen(uQuad, W[i][j][n], 
                                               u0Min, u0Max, 4);
                break;
              case LIMITER_VENKATAKRISHNAN :
                phi_n = Limiter_Venkatakrishnan(uQuad, W[i][j][n], 
                                                u0Min, u0Max, 4);
                break;
              case LIMITER_VANLEER :
                phi_n = Limiter_VanLeer(uQuad, W[i][j][n], 
                                        u0Min, u0Max, 4);
                break;
              case LIMITER_VANALBADA :
                phi_n = Limiter_VanAlbada(uQuad, W[i][j][n], 
                                          u0Min, u0Max, 4);
                break;
              default:
                phi_n = Limiter_BarthJespersen(uQuad, W[i][j][n], 
                                               u0Min, u0Max, 4);
                break;
            } /* endswitch */

	    phi[i][j][n] = phi_n;

         } /* endfor */
      } /* endif */
  } else {
      dWdx[i][j] = Euler2D_W_VACUUM;
      dWdy[i][j] = Euler2D_W_VACUUM; 
      phi[i][j]  = Euler2D_W_VACUUM;
  } /* endif */
    
}

/*******************************************************************************
 * Euler2D_Quad_Block::LoadSendBuffer_Flux_F2C -- Loads send message buffer for*
 *                                                fine to coarse block message *
 *                                                passing of conservative      *
 *                                                solution fluxes.             *
 *******************************************************************************/
inline int Euler2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) {
  int i, j, k;
  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxS[i  ][k]+
                                   FluxS[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxN[i  ][k]+
                                   FluxN[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxW[j][k]+
                                   FluxW[j+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxE[j][k]+
                                   FluxE[j+1][k]);
        } /* endfor */
     } /* endfor */
  } /* endif */
  return(0);
}

/*******************************************************************************
 * Euler2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message *
 *                                                buffer for fine to coarse    *
 *                                                block message passing of     *
 *                                                conservative solution fluxes.*
 *******************************************************************************/
inline int Euler2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
                                                            int &buffer_count,
                                                            const int buffer_size,
                                                            const int i_min, 
                                                            const int i_max,
                                                            const int i_inc,
                                                            const int j_min, 
                                                            const int j_max,
                                                            const int j_inc) {
  int i, j;
  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
  	buffer_count = buffer_count + NUM_VAR_EULER2D;
        if (buffer_count >= buffer_size) return(1);
        FluxS[i] = -Euler2D_cState(buffer[buffer_count-3],
                                   buffer[buffer_count-2],
                                   buffer[buffer_count-1],
                                   buffer[buffer_count])
                   -FluxS[i];
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
  	buffer_count = buffer_count + NUM_VAR_EULER2D;
        if (buffer_count >= buffer_size) return(1);
        FluxN[i] = -Euler2D_cState(buffer[buffer_count-3],
                                   buffer[buffer_count-2],
                                   buffer[buffer_count-1],
                                   buffer[buffer_count])
                   -FluxN[i];
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
  	buffer_count = buffer_count + NUM_VAR_EULER2D;
        if (buffer_count >= buffer_size) return(1);
        FluxW[j] = -Euler2D_cState(buffer[buffer_count-3],
                                   buffer[buffer_count-2],
                                   buffer[buffer_count-1],
                                   buffer[buffer_count])
                   -FluxW[j];
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
  	buffer_count = buffer_count + NUM_VAR_EULER2D;
        if (buffer_count >= buffer_size) return(1);
        FluxE[j] = -Euler2D_cState(buffer[buffer_count-3],
                                   buffer[buffer_count-2],
                                   buffer[buffer_count-1],
                                   buffer[buffer_count])
                   -FluxE[j];
     } /* endfor */
  } /* endif */
  return(0);
}

/**************************************************************************
 * Euler2D_Quad_Block -- Single Block External Subroutines.               *
 **************************************************************************/

extern void Write_Solution_Block(Euler2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(Euler2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Broadcast_Solution_Block(Euler2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(Euler2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator, 
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(Euler2D_Quad_Block &SolnBlk1,
		                Euler2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(Euler2D_Quad_Block &SolnBlk_Fine,
				  Euler2D_Quad_Block &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(Euler2D_Quad_Block &SolnBlk_Coarse,
				   Euler2D_Quad_Block &SolnBlk_Original_SW,
				   Euler2D_Quad_Block &SolnBlk_Original_SE,
				   Euler2D_Quad_Block &SolnBlk_Original_NW,
				   Euler2D_Quad_Block &SolnBlk_Original_NE);

extern void Output_Tecplot(Euler2D_Quad_Block &SolnBlk,
			   Euler2D_Input_Parameters &IP,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(Euler2D_Quad_Block &SolnBlk,
		                 Euler2D_Input_Parameters &IP,
                                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nodes_Tecplot(Euler2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Gradients_Tecplot(Euler2D_Quad_Block &SolnBlk,
				     const int Number_of_Time_Steps,
				     const double &Time,
				     const int Block_Number,
				     const int Output_Title,
				     ostream &Out_File);

extern void Output_Tecplot_Quasi3D(Euler2D_Quad_Block &SolnBlk,
				   Euler2D_Input_Parameters &IP,
				   const int Number_of_Time_Steps,
				   const double &Time,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File);

extern void Output_Ringleb_Flow_Solution(Euler2D_Quad_Block &SolnBlk,
					 const int Block_Number,
					 const int Output_Title,
					 ostream &Out_File);

extern void Output_Ringleb_Flow_Error(Euler2D_Quad_Block &SolnBlk,
				      double &l1_norm,
				      double &l2_norm,
				      double &max_norm);

extern void ICs(Euler2D_Quad_Block &SolnBlk,
		Euler2D_Input_Parameters &IP,
                Euler2D_pState *Wo);

extern void BCs(Euler2D_Quad_Block &SolnBlk,
		Euler2D_Input_Parameters &IP);

extern double CFL(Euler2D_Quad_Block &SolnBlk,
                  Euler2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Euler2D_Quad_Block &SolnBlk, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Euler2D_Quad_Block &SolnBlk, const int &norm);

extern double L2_Norm_Residual(Euler2D_Quad_Block &SolnBlk, const int &norm);

extern double Max_Norm_Residual(Euler2D_Quad_Block &SolnBlk, const int &norm);

extern void Linear_Reconstruction_GreenGauss(Euler2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Euler2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Euler2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_2(Euler2D_Quad_Block &SolnBlk,
                                                 const int i,
                                                 const int j,
					         const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Euler2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Linear_Reconstruction(Euler2D_Quad_Block &SolnBlk,
				  const int & Reconstruction_Type,
				  const int & Limiter);

extern void Residual_Smoothing(Euler2D_Quad_Block &SolnBlk,
                               const int k_residual,
			       double &epsilon, 
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  Euler2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          Euler2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(Euler2D_Quad_Block &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(Euler2D_Quad_Block &SolnBlk);

extern void Apply_Boundary_Flux_Corrections(Euler2D_Quad_Block &SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Euler2D_Quad_Block &SolnBlk,
                                                                const int i_stage,
                                                                const int n_stage,
                                                                const double &CFL_Number,
                                                                const int Time_Integration_Type,
                                                                const int Reconstruction_Type,
                                                                const int Limiter_Type,
                                                                const int Flux_Function_Type,
                                                                const int Number_Neighbours_North_Boundary,
                                                                const int Number_Neighbours_South_Boundary,
                                                                const int Number_Neighbours_East_Boundary,
                                                                const int Number_Neighbours_West_Boundary);

extern int dUdt_Residual_Evaluation(Euler2D_Quad_Block &SolnBlk,
				    Euler2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(Euler2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    Euler2D_Input_Parameters &Input_Parameters);

extern int Update_Solution_Multistage_Explicit(Euler2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
                                               Euler2D_Input_Parameters &Input_Parameters);

extern int Determine_Maximum_Surface_Pressure(Euler2D_Quad_Block &SolnBlk,
					      double &maximum_pressure);

extern int Determine_Mach_Stem_Height(Euler2D_Quad_Block &SolnBlk,
				      Euler2D_Input_Parameters &Input_Parameters,
				      const double &ms_pressure,
				      const int Number_of_Time_Steps,
				      const double &Time,
				      const int Block_Number,
				      const int Output_Title,
				      ostream &Out_File);

extern void Output_Wedge_Solution_Distribution_Tecplot(Euler2D_Quad_Block &SolnBlk,
						       const double &wedge_angle,
						       const int Number_of_Time_Steps,
						       const double &Time,
						       const int Block_Number,
						       const int Output_Title,
						       ostream &Out_File);

/**************************************************************************
 * Euler2D_Quad_Block -- Multiple Block External Subroutines.             *
 **************************************************************************/

extern Euler2D_Quad_Block* Allocate(Euler2D_Quad_Block *Soln_ptr,
                                    Euler2D_Input_Parameters &Input_Parameters);

extern Euler2D_Quad_Block* Deallocate(Euler2D_Quad_Block *Soln_ptr,
                                      Euler2D_Input_Parameters &Input_Parameters);

extern Euler2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block_HO **InitMeshBlk,
						       Euler2D_Quad_Block *Soln_ptr,
						       Euler2D_Input_Parameters &Input_Parameters,
						       QuadTreeBlock_DataStructure &QuadTree,
						       AdaptiveBlockResourceList &GlobalSolnBlockList,
						       AdaptiveBlock2D_List &LocalSolnBlockList);

extern void ICs(Euler2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                Euler2D_Input_Parameters &Input_Parameters);

extern int Read_Restart_Solution(Euler2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 Euler2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(Euler2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  Euler2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(Euler2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          Euler2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(Euler2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                Euler2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Nodes_Tecplot(Euler2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                Euler2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Gradients_Tecplot(Euler2D_Quad_Block *Soln_ptr,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    Euler2D_Input_Parameters &Input_Parameters,
				    const int Number_of_Time_Steps,
				    const double &Time);

extern int Output_Tecplot_Quasi3D(Euler2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  Euler2D_Input_Parameters &Input_Parameters,
				  const int Number_of_Time_Steps,
				  const double &Time);

extern int Output_Mesh_Tecplot(Euler2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Euler2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(Euler2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Euler2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Ringleb_Flow(Euler2D_Quad_Block *Soln_ptr,
			       AdaptiveBlock2D_List &Soln_Block_List,
			       Euler2D_Input_Parameters &Input_Parameters);

extern void Linear_Reconstruction(Euler2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  Euler2D_Input_Parameters &Input_Parameters);

extern void BCs(Euler2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		Euler2D_Input_Parameters &Input_Parameters);

extern double CFL(Euler2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
                  Euler2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Euler2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Euler2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(Euler2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(Euler2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void L1_Norm_Residual(Euler2D_Quad_Block *Soln_ptr, 
			     AdaptiveBlock2D_List &Soln_Block_List,
			     double *l1_norm);

extern void L2_Norm_Residual(Euler2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     double *l2_norm);

extern void Max_Norm_Residual(Euler2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      double *max_norm);

extern void Evaluate_Limiters(Euler2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(Euler2D_Quad_Block *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Residual_Smoothing(Euler2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Euler2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(Euler2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Euler2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                Euler2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Multistage_Explicit(Euler2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Soln_Block_List,
                                    Euler2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(Euler2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               Euler2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

extern int Determine_Maximum_Surface_Pressure(Euler2D_Quad_Block *Soln_ptr,
					      AdaptiveBlock2D_List &Soln_Block_List,
					      Euler2D_Input_Parameters &Input_Parameters,
					      double &maximum_pressure);

extern int Determine_Mach_Stem_Height(Euler2D_Quad_Block *Soln_ptr,
				      AdaptiveBlock2D_List &Soln_Block_List,
				      Euler2D_Input_Parameters &Input_Parameters,
				      const double &ms_pressure,
				      const int Number_of_Time_Steps,
				      const double &Time);

extern int Output_Wedge_Solution_Distribution_Tecplot(Euler2D_Quad_Block *Soln_ptr,
						      AdaptiveBlock2D_List &Soln_Block_List,
						      Euler2D_Input_Parameters &Input_Parameters,
						      const int Number_of_Time_Steps,
						      const double &Time);

/**************************************************************************
 * Euler2D_Quad_Block -- Multiple Block External Subroutines for Mesh.    *
 **************************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            Euler2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      Euler2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             Euler2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                            Euler2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  Euler2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                 Euler2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          Euler2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Euler2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Euler2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid3D_Hexa_Block ***Grid_ptr,
                          Euler2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Blocks_Kdir);

/**************************************************************************
 * Euler2D_Quad_Block -- Solvers.                                         *
 **************************************************************************/

extern int Euler2DQuadSolver(char *Input_File_Name_ptr,
                             int batch_flag);


/*
 * Include specializations CFFC header files.
 * Must be included at the end of the file!!!
 */
#include "Euler2DAccuracyAssessment.h" /* Include 2D accuracy assessment header 
					  file with specializations for Euler solution. */

#endif /* _EULER2D_QUAD_INCLUDED  */
