/* HyperHeat1D.h:  Header file defining 1D Hyperbolic Heat Equations
                   (Maxwell-Cattaneo Equation) Solution Classes. */

#ifndef _HYPERHEAT1D_INCLUDED
#define _HYPERHEAT1D_INCLUDED

/* Include 1D hyperbolic heat equations solution state and
   cell header files. */

#ifndef _HYPERHEAT1D_STATE_INCLUDED
#include "HyperHeat1DState.h"
#endif // _HYPERHEAT1D_STATE_INCLUDED

#ifndef _CELL1D_INCLUDED
#include "../Grid/Cell1D.h"
#endif // _CELL1D_INCLUDED

/* Define the classes. */

/********************************************************
 * Class: HyperHeat1D_UniformMesh                       *
 *                                                      *
 * Member functions                                     *
 *      U       -- Return solution state.               *
 *      X       -- Return cell geometry.                *
 *      dt      -- Return local time step.              *
 *    dUdt      -- Return the solution residual.        *
 *    dUdx      -- Return the unlimited solution        *
 *                 gradient.                            *
 *     phi      -- Return the solution slope limiters.  *
 *      Uo      -- Return initial solution state.       *
 *                                                      *
 * Member operators                                     *
 *      S -- a 1D hyperbolic heat equation solution     *
 *                                                      *
 * S = S;                                               *
 * S = S + S;                                           *
 * S = S - S;                                           *
 * S = +S;                                              *
 * S = -S;                                              *
 * S += S;                                              *
 * S -= S;                                              *
 * S == S;                                              *
 * S != S;                                              *
 * cout << S; (output function)                         *
 * cin  >> S; (input function)                          *
 *                                                      *
 ********************************************************/
class HyperHeat1D_UniformMesh{
  private:
  public:
    HyperHeat1D_State     U;   // Conserved solution state.
    Cell1D_Uniform        X;   // Cell geometry.
    double               dt;   // Local time step.
    HyperHeat1D_State  dUdt;   // Solution residual.
    HyperHeat1D_State  dUdx;   // Unlimited solution gradient.
    HyperHeat1D_State   phi;   // Solution slope limiter.
    HyperHeat1D_State    Uo;   // Initial solution state.
 	                       // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    HyperHeat1D_UniformMesh(void) {
       U = HyperHeat1D_State(ZERO,ZERO); X = Cell1D_Uniform_ONE; dt = ZERO;
       dUdt = HyperHeat1D_State(ZERO,ZERO); dUdx = HyperHeat1D_State(ZERO,ZERO);
       phi = HyperHeat1D_State(ZERO,ZERO); Uo = HyperHeat1D_State(ZERO,ZERO);
    }

    HyperHeat1D_UniformMesh(const HyperHeat1D_UniformMesh &Soln) {
       U = Soln.U; X = Soln.X; dt = Soln.dt;
       dUdt = Soln.dUdt; dUdx = Soln.dUdx; phi = Soln.phi; Uo = Soln.Uo;
    }

    HyperHeat1D_UniformMesh(const HyperHeat1D_State &U0,
			    const Cell1D_Uniform &X0) {
       U = U0; X = X0; dt = ZERO;
       dUdt = HyperHeat1D_State(ZERO,ZERO); dUdx = HyperHeat1D_State(ZERO,ZERO);
       phi = HyperHeat1D_State(ZERO,ZERO); Uo = HyperHeat1D_State(ZERO,ZERO);
    }

    /* Destructor. */
    // ~HyperHeat1D_UniformMesh(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // HyperHeat1D_UniformMesh operator = (const HyperHeat1D_UniformMesh &Soln);
    // Use automatically generated assignment operator.

    /* Binary arithmetic operators. */
//    friend HyperHeat1D_UniformMesh operator +(const HyperHeat1D_UniformMesh &Soln1,
//					        const HyperHeat1D_UniformMesh &Soln2);
//    friend HyperHeat1D_UniformMesh operator -(const HyperHeat1D_UniformMesh &Soln1,
//					        const HyperHeat1D_UniformMesh &Soln2);

    /* Unary arithmetic operators. */
//    friend HyperHeat1D_UniformMesh operator +(const HyperHeat1D_UniformMesh &Soln);
//    friend HyperHeat1D_UniformMesh operator -(const HyperHeat1D_UniformMesh &Soln);

    /* Shortcut arithmetic operators. */
//    friend HyperHeat1D_UniformMesh &operator +=(HyperHeat1D_UniformMesh &Soln1,
//					          const HyperHeat1D_UniformMesh &Soln2);
//    friend HyperHeat1D_UniformMesh &operator -=(HyperHeat1D_UniformMesh &Soln1,
//					          const HyperHeat1D_UniformMesh &Soln2);
    
    /* Relational operators. */
    friend int operator ==(const HyperHeat1D_UniformMesh &Soln1,
			   const HyperHeat1D_UniformMesh &Soln2);
    friend int operator !=(const HyperHeat1D_UniformMesh &Soln1,
			   const HyperHeat1D_UniformMesh &Soln2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const HyperHeat1D_UniformMesh &Soln);
    friend istream &operator >> (istream &in_file,
				 HyperHeat1D_UniformMesh &Soln);
    
};

/********************************************************
 * HyperHeat1D_UniformMesh -- Relational operators.     *
 ********************************************************/
inline int operator ==(const HyperHeat1D_UniformMesh &Soln1,
		       const HyperHeat1D_UniformMesh &Soln2) {
    return (Soln1.U == Soln2.U && Soln1.X == Soln2.X && Soln1.dt == Soln2.dt &&
	    Soln1.dUdt == Soln2.dUdt && Soln1.dUdx == Soln2.dUdx &&
	    Soln1.phi == Soln2.phi && Soln1.Uo == Soln2.Uo);
}

inline int operator !=(const HyperHeat1D_UniformMesh &Soln1,
		       const HyperHeat1D_UniformMesh &Soln2) {
    return (Soln1.U != Soln2.U || Soln1.X != Soln2.X || Soln1.dt != Soln2.dt ||
	    Soln1.dUdt != Soln2.dUdt || Soln1.dUdx != Soln2.dUdx ||
	    Soln1.phi != Soln2.phi || Soln1.Uo == Soln2.Uo);
}

/********************************************************
 * HyperHeat1D_UniformMesh -- Input-output operators.   *
 ********************************************************/
inline ostream &operator << (ostream &out_file,
			     const HyperHeat1D_UniformMesh &Soln) {
    out_file << Soln.X << Soln.U;
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     HyperHeat1D_UniformMesh &Soln) {
    in_file >> Soln.X >> Soln.U;
    return (in_file);
}

/********************************************************
 * HyperHeat1D_UniformMesh -- External subroutines.     *
 ********************************************************/

extern HyperHeat1D_UniformMesh* Allocate(HyperHeat1D_UniformMesh *Soln_ptr,
				         const int Number_of_Cells);

extern HyperHeat1D_UniformMesh* Deallocate(HyperHeat1D_UniformMesh *Soln_ptr,
				           const int Number_of_Cells);

extern void Output_Gnuplot(HyperHeat1D_UniformMesh *Soln,
  		           const int Number_of_Cells,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Output_Tecplot(HyperHeat1D_UniformMesh *Soln,
  		           const int Number_of_Cells,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Grid(HyperHeat1D_UniformMesh *Soln,
                 const double &xmin,
		 const double &xmax,
		 const int Number_of_Cells);

extern void ICs(HyperHeat1D_UniformMesh *Soln,
		const int i_ICtype,
                const double &Kappa,
                const double &Tau,
		const int Number_of_Cells);

extern double CFL(HyperHeat1D_UniformMesh *Soln,
  		  const int Number_of_Cells);

extern void Linear_Reconstruction_MUSCL(HyperHeat1D_UniformMesh *Soln,
                                        const int Number_of_Cells,
                                        const int Limiter);

extern void Linear_Reconstruction_GreenGauss(HyperHeat1D_UniformMesh *Soln,
                                             const int Number_of_Cells,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(HyperHeat1D_UniformMesh *Soln,
                                               const int Number_of_Cells,
					       const int Limiter);

extern int dUdt_explicitEuler_upwind(HyperHeat1D_UniformMesh *Soln,
    	                             const int Number_of_Cells,
                                     const int BC_type_left,
				     const int BC_type_right,
                                     const HyperHeat1D_State U_left,
                                     const HyperHeat1D_State U_right,
			             double &dtMin,
				     const double &CFL_Number,
                                     const int Flux_Function_Type,
			             const int Local_Time_Stepping);

extern int dUdt_semiimplicitEuler_upwind(HyperHeat1D_UniformMesh *Soln,
    	                                 const int Number_of_Cells,
                                         const int BC_type_left,
					 const int BC_type_right,
                                         const HyperHeat1D_State U_left,
                                         const HyperHeat1D_State U_right,
			                 double &dtMin,
					 const double &CFL_Number,
                                         const int Flux_Function_Type,
			                 const int Local_Time_Stepping);

extern int dUdt_2stage_2ndOrder_upwind(HyperHeat1D_UniformMesh *Soln,
    	                               const int Number_of_Cells,
                                       const int BC_type_left,
				       const int BC_type_right,
                                       const HyperHeat1D_State U_left,
                                       const HyperHeat1D_State U_right,
			               double &dtMin,
				       const double &CFL_Number,
                                       const int Reconstruction_Type,
                                       const int Limiter_Type,
                                       const int Flux_Function_Type,
			               const int Local_Time_Stepping);

/********************************************************
 * HyperHeat1D_UniformMesh -- Solvers.                  *
 ********************************************************/

extern int HyperHeat1DSolver(char *Input_File_Name_ptr,
                             int batch_flag);

#endif /* _HYPERHEAT1D_INCLUDED  */
