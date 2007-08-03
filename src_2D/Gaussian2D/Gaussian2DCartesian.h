/* Gaussian2DCartesian.h:  Header file defining 
                           2D Gaussian Cartesian Mesh Solution Classes. */

#ifndef _GAUSSIAN2D_CARTESIAN_INCLUDED
#define _GAUSSIAN2D_CARTESIAN_INCLUDED

/* Include 2D Gaussian state, 2D cell, and 2D Gaussian input header files. */

#ifndef _GAUSSIAN2D_STATE_INCLUDED
#include "Gaussian2DState.h"
#endif // _GAUSSIAN2D_STATE_INCLUDED

#ifndef _CELL2D_INCLUDED
#include "../Grid/Cell2D.h"
#endif // _CELL2D_INCLUDED

#ifndef _GAUSSIAN2D_INPUT_INCLUDED
#include "Gaussian2DInput.h"
#endif // _GAUSSIAN2D_INPUT_INCLUDED

/* Define the structures and classes. */

/********************************************************
 * Class: Gaussian2D_Cartesian_Cell                     *
 *                                                      *
 * Member functions                                     *
 *      W       -- Return primitive solution state.     *
 *      U       -- Return conserved solution state.     *
 *      X       -- Return 2D cell geometry.             *
 *      dt      -- Return local time step.              *
 *    dUdt      -- Return the solution residual.        *
 *    dWdx      -- Return the unlimited solution        *
 *                 gradient (x-direction)               *
 *    dWdy      -- Return the unlimited solution        *
 *                 gradient (y-direction)               *
 *     phi      -- Return the solution slope limiters.  *
 *      Uo      -- Return initial solution state.       *
 *                                                      *
 * Member operators                                     *
 *      S -- a 2D Gaussian solution                     *
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
class Gaussian2D_Cartesian_Cell{
  private:
  public:
    Gaussian2D_pState         W;  // Primitive solution state.
    Gaussian2D_cState         U;  // Conserved solution state.
    Cell2D_Cartesian          X;  // 2D Cell geometry.
    double                   dt;  // Local time step.
    Gaussian2D_cState      dUdt;  // Solution residual.
    Gaussian2D_pState      dWdx;  // Unlimited solution gradient
                                  // (x-direction).
    Gaussian2D_pState      dWdy;  // Unlimited solution gradient
                                  // (y-direction).
    Gaussian2D_pState       phi;  // Solution slope limiter.
    Gaussian2D_cState        Uo;  // Initial solution state.
                                  // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Gaussian2D_Cartesian_Cell(void) {
       W = Gaussian2D_W_STDATM; U = Gaussian2D_U_STDATM;
       X = Cell2D_Cartesian_ONE; dt = ZERO;
       dUdt = Gaussian2D_U_VACUUM; dWdx = Gaussian2D_W_VACUUM;
       dWdy = Gaussian2D_W_VACUUM; phi = Gaussian2D_W_VACUUM;
       Uo = Gaussian2D_U_VACUUM;
    }

    Gaussian2D_Cartesian_Cell(const Gaussian2D_Cartesian_Cell &Soln) {
       W = Soln.W; U = Soln.U; X = Soln.X;
       dt = Soln.dt; dUdt = Soln.dUdt; dWdx = Soln.dWdx;
       dWdy = Soln.dWdy; phi = Soln.phi; Uo = Soln.Uo;
    }

    Gaussian2D_Cartesian_Cell(const Gaussian2D_pState &W0,
		   	      const Gaussian2D_cState &U0,
			      const Cell2D_Cartesian &X0) {
       W = W0; U = U0; X = X0;
       dt = ZERO; dUdt = Gaussian2D_U_VACUUM;
       dWdx = Gaussian2D_W_VACUUM; dWdy = Gaussian2D_W_VACUUM;
       phi = Gaussian2D_W_VACUUM; Uo = Gaussian2D_U_VACUUM;
    }

    Gaussian2D_Cartesian_Cell(const Gaussian2D_pState &W0,
			      const Cell2D_Cartesian &X0) {
       W = W0; U = W0.U(); X = X0;
       dt = ZERO; dUdt = Gaussian2D_U_VACUUM;
       dWdx = Gaussian2D_W_VACUUM; dWdy = Gaussian2D_W_VACUUM;
       phi = Gaussian2D_W_VACUUM; Uo = Gaussian2D_U_VACUUM;
    }

    Gaussian2D_Cartesian_Cell(const Gaussian2D_cState &U0,
			      const Cell2D_Cartesian &X0) {
       W = U0.W(); U = U0; X = X0;
       dt = ZERO; dUdt = Gaussian2D_U_VACUUM;
       dWdx = Gaussian2D_W_VACUUM; dWdy = Gaussian2D_W_VACUUM;
       phi = Gaussian2D_W_VACUUM; Uo = Gaussian2D_U_VACUUM;
    }
    
    /* Destructor. */
    // ~Gaussian2D_Cartesian_Cell(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // Gaussian2D_Cartesian_Cell operator = (const Gaussian2D_Cartesian_Cell &Soln);
    // Use automatically generated assignment operator.

    /* Relational operators. */
    friend int operator ==(const Gaussian2D_Cartesian_Cell &Soln1,
			   const Gaussian2D_Cartesian_Cell &Soln2);
    friend int operator !=(const Gaussian2D_Cartesian_Cell &Soln1,
			   const Gaussian2D_Cartesian_Cell &Soln2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Gaussian2D_Cartesian_Cell &Soln);
    friend istream &operator >> (istream &in_file,
				 Gaussian2D_Cartesian_Cell &Soln);
    
};

/*************************************************************
 * Gaussian2D_Cartesian_Cell -- Relational operators.           *
 *************************************************************/
inline int operator ==(const Gaussian2D_Cartesian_Cell &Soln1,
		       const Gaussian2D_Cartesian_Cell &Soln2) {
    return (Soln1.W == Soln2.W && Soln1.U == Soln2.U &&
	    Soln1.X == Soln2.X && Soln1.dt == Soln2.dt &&
	    Soln1.dUdt == Soln2.dUdt && Soln1.dWdx == Soln2.dWdx &&
	    Soln1.dWdy == Soln2.dWdy && Soln1.phi == Soln2.phi &&
	    Soln1.Uo == Soln2.Uo);
}

inline int operator !=(const Gaussian2D_Cartesian_Cell &Soln1,
		       const Gaussian2D_Cartesian_Cell &Soln2) {
    return (Soln1.W != Soln2.W || Soln1.U != Soln2.U ||
	    Soln1.X != Soln2.X || Soln1.dt != Soln2.dt ||
	    Soln1.dUdt != Soln2.dUdt || Soln1.dWdx != Soln2.dWdx ||
	    Soln1.dWdy != Soln2.dWdy || Soln1.phi != Soln2.phi ||
	    Soln1.Uo == Soln2.Uo);
}

/*************************************************************
 * Gaussian2D_Cartesian_Cell -- Input-output operators.         *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Gaussian2D_Cartesian_Cell &Soln) {
    out_file << Soln.X << Soln.W;
    out_file.setf(ios::scientific);
    out_file << " " << Soln.W.T();
    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Gaussian2D_Cartesian_Cell &Soln) {
    in_file >> Soln.X >> Soln.W;
    Soln.U = Soln.W.U();
    return (in_file);
}

/*************************************************************
 * Gaussian2D_Cartesian_Cell -- External subroutines.           *
 *************************************************************/

extern Gaussian2D_Cartesian_Cell** Allocate(Gaussian2D_Cartesian_Cell **Soln_ptr,
			  		    const int Number_of_Cells_Xdir,
		                            const int Number_of_Cells_Ydir);

extern Gaussian2D_Cartesian_Cell** Deallocate(Gaussian2D_Cartesian_Cell **Soln_ptr,
					      const int Number_of_Cells_Xdir,
		                              const int Number_of_Cells_Ydir);

extern void Output_Gnuplot(Gaussian2D_Cartesian_Cell **Soln,
  		           const int Number_of_Cells_Xdir,
  		           const int Number_of_Cells_Ydir,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Output_Tecplot(Gaussian2D_Cartesian_Cell **Soln,
  		           const int Number_of_Cells_Xdir,
  		           const int Number_of_Cells_Ydir,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Grid(Gaussian2D_Cartesian_Cell **Soln,
                 const double &xmin,
		 const double &xmax,
                 const double &ymin,
		 const double &ymax,
		 const int Number_of_Cells_Xdir,
		 const int Number_of_Cells_Ydir);

extern void ICs(Gaussian2D_Cartesian_Cell **Soln,
                char *gas_ptr,
		const int i_ICtype,
		const int Number_of_Cells_Xdir,
		const int Number_of_Cells_Ydir);

extern double CFL(Gaussian2D_Cartesian_Cell **Soln,
   		  const int Number_of_Cells_Xdir,
                  const int Number_of_Cells_Ydir);

extern void Linear_Reconstruction_GreenGauss(Gaussian2D_Cartesian_Cell **Soln,
   		                             const int Number_of_Cells_Xdir,
                                             const int Number_of_Cells_Ydir,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Gaussian2D_Cartesian_Cell **Soln,
   		                               const int Number_of_Cells_Xdir,
                                               const int Number_of_Cells_Ydir,
					       const int Limiter);

extern int dUdt_explicitEuler_upwind(Gaussian2D_Cartesian_Cell **Soln,
   		                     const int Number_of_Cells_Xdir,
                                     const int Number_of_Cells_Ydir,
				     const int BC_type_left,
	      			     const int BC_type_right,
				     const int BC_type_bottom,
	      			     const int BC_type_top,
                                     const Gaussian2D_pState W_left,
                                     const Gaussian2D_pState W_right,
                                     const Gaussian2D_pState W_bottom,
                                     const Gaussian2D_pState W_top,
			             double &dtMin,
				     const double &CFL_Number,
				     const int Flux_Function_Type,
			             const int Local_Time_Stepping);

extern int dUdt_2stage_2ndOrder_upwind(Gaussian2D_Cartesian_Cell **Soln,
   		                       const int Number_of_Cells_Xdir,
                                       const int Number_of_Cells_Ydir,
				       const int BC_type_left,
	      			       const int BC_type_right,
				       const int BC_type_bottom,
	      			       const int BC_type_top,
                                       const Gaussian2D_pState W_left,
                                       const Gaussian2D_pState W_right,
                                       const Gaussian2D_pState W_bottom,
                                       const Gaussian2D_pState W_top,
			               double &dtMin,
				       const double &CFL_Number,
                                       const int Reconstruction_Type,
                                       const int Limiter_Type,
				       const int Flux_Function_Type,
			               const int Local_Time_Stepping);


/*************************************************************
 * Gaussian2D_Cartesian_Cell -- Solvers.                     *
 *************************************************************/

extern int Gaussian2DCartesianSolver(char *Input_File_Name_ptr,
                                     int batch_flag);

#endif /* _GAUSSIAN2D_CARTESIAN_INCLUDED  */











