/*!\file BGK1D.h
  \brief Header file defining 1D BGK Solution Classes. */

#ifndef _BGK1D_INCLUDED
#define _BGK1D_INCLUDED

/* Include 1D BGK state and cell header files. */

#ifndef _BGK1D_VECTOR_INCLUDED
#include "BGK1DVector.h"
#endif // _BGK1D_VECTOR_INCLUDED

#ifndef _CELL1D_INCLUDED
#include "../Grid/Cell1D.h"
#endif // _CELL1D_INCLUDED

#ifndef _CFD1DINPUT_INCLUDED
#include "../CFD/CFD1DInput.h"
#endif // _CFD1DINPUT_INCLUDED

using namespace std;

/* Define the classes. */

/******************************************************//**
 * Class: BGK1D_UniformMesh
 *
 * @brief Class definition of the 1D BGK computational cell (i.e. solution + geometry)
 *
 * Member functions
 * -     V       -- Return vector of distribution function discretization.
 * -     X       -- Return cell geometry.
 * -     dt      -- Return local time step.
 * -   dVdt      -- Return the solution residual.
 * -   dVdx      -- Return the unlimited solution
 * -                gradient.
 * -    phi      -- Return the solution slope limiters.
 * -     Vo      -- Return initial solution state.
 * - Nghost      -- Number of ghost cells (!= 1 for
 * -                high-order)
 *
 * - cout << S; (output function)
 * - cin  >> S; (input function)
 *
 ********************************************************/
class BGK1D_UniformMesh{
 public:

  BGK1D_Vector    V;   //!< discretization of distribution function
  Cell1D_Uniform  X;   //!< Cell geometry.
  BGK1D_Vector dVdt;   //!< Solution residual.
  BGK1D_Vector dVdx;   //!< Unlimited solution gradient.
  BGK1D_Vector  phi;   //!< Solution slope limiter.
  BGK1D_Vector   Vo;   //!< Initial solution state.
  // Made public so can access them.

  int Nghost;            //!< Number of ghost cells(!= 1 for high-order)
  int ICl, ICu;	         //!< Indexes (Start & End)

  /* Creation, copy, and assignment constructors. */
  BGK1D_UniformMesh(void);
  BGK1D_UniformMesh(const BGK1D_UniformMesh &Soln);

  /* Field access */
  const double & CellCenter(void) const {return X.x;} //!< return cell center
  const double & CellDelta (void) {return X.dx;} //!< return cell delta

  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file,
			       const BGK1D_UniformMesh &Soln);
  friend istream &operator >> (istream &in_file,
			       BGK1D_UniformMesh &Soln);

 private:

};

/*******************************************************
 *  BGK1D_UniformMesh::MemberFunctions()         *
 ******************************************************/
/* Constructor */
inline BGK1D_UniformMesh::BGK1D_UniformMesh(void){
  V.zero(); Vo.zero();
  X = Cell1D_Uniform_ONE;
  dVdt.zero(); dVdx.zero(); phi.zero();
}

inline BGK1D_UniformMesh::BGK1D_UniformMesh(const BGK1D_UniformMesh &Soln) {
  V = Soln.V;
  dVdt = Soln.dVdt; dVdx = Soln.dVdx;
  phi = Soln.phi; Vo = Soln.Vo;
  Nghost = Soln.Nghost; ICl = Soln.ICl; ICu = Soln.ICu;
}

/********************************************************
 * BGK1D_UniformMesh -- Input-output operators.       *
 ********************************************************/
inline ostream &operator << (ostream &out_file,
			     const BGK1D_UniformMesh &Soln) {
  out_file << Soln.X << Soln.V;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     BGK1D_UniformMesh &Soln) {
  in_file >> Soln.X >> Soln.V;
  return (in_file);
}

/********************************************************
 * BGK1D_UniformMesh -- External subroutines.         *
 ********************************************************/

extern BGK1D_UniformMesh* Allocate(BGK1D_UniformMesh *Soln_ptr,
				   const CFD1D_Input_Parameters &IP);

extern BGK1D_UniformMesh* Deallocate(BGK1D_UniformMesh *Soln_ptr);

extern void Output_Tecplot(BGK1D_UniformMesh *Soln,
  		           const int Number_of_Cells,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

void Output_Tecplot(BGK1D_UniformMesh *Soln,
		    const CFD1D_Input_Parameters &IP,
		    const int Number_of_Time_Steps,
		    const double &Time,
		    ostream &out_file);

extern void Grid(BGK1D_UniformMesh *Soln,
                 const double &xmin,
		 const double &xmax,
		 const int Number_of_Cells);

extern void ICs(BGK1D_UniformMesh *Soln,
                char *gas_ptr,
		const int i_ICtype,
		const int Number_of_Cells,
		CFD1D_Input_Parameters &IP);

extern void Linear_Reconstruction_MUSCL(BGK1D_UniformMesh *Soln,
                                        const int Number_of_Cells,
					const int Limiter);

extern void Linear_Reconstruction_GreenGauss(BGK1D_UniformMesh *Soln,
                                             const int Number_of_Cells,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(BGK1D_UniformMesh *Soln,
                                               const int Number_of_Cells,
					       const int Limiter);

extern int dVdt_explicitEuler_upwind(BGK1D_UniformMesh *Soln,
				     const int Number_of_Cells,
				     double &dtMin,
				     const double &CFL_Number,
				     const int Flux_Function_Type,
				     const int Local_Time_Stepping);

extern int dVdt_2stage_2ndOrder_upwind(BGK1D_UniformMesh *Soln,
				       const CFD1D_Input_Parameters &IP,
			               double &dtMin,
				       const double &CFL_Number,
                                       const int Reconstruction_Type,
                                       const int Limiter_Type,
				       const int Flux_Function_Type,
			               const int Local_Time_Stepping);

void LimitedLinearReconstructionOverDomain(BGK1D_UniformMesh *Soln, const CFD1D_Input_Parameters &IP);

/********************************************************
 * BGK1D_UniformMesh -- Solvers.                      *
 ********************************************************/

extern int BGK1DSolver(char *Input_File_Name_ptr,
		       int batch_flag);

#endif /* _BGK1D_INCLUDED  */
