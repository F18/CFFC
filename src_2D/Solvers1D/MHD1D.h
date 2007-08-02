/* MHD1D.h:  Header file defining 1D MHD Solution Classes. */

#ifndef _MHD1D_INCLUDED
#define _MHD1D_INCLUDED

/* Include 1D MHD state and cell header files. */

#ifndef _MHD1D_STATE_INCLUDED
#include "MHD1DState.h"
#endif // _MHD1D_STATE_INCLUDED

#ifndef _CELL1D_INCLUDED
#include "../Grid/Cell1D.h"
#endif // _CELL1D_INCLUDED

/* Define the classes. */

/********************************************************
 * Class: MHD1D_UniformMesh                             *
 *                                                      *
 * Member functions                                     *
 *      W       -- Return primitive solution state.     *
 *      U       -- Return conserved solution state.     *
 *      X       -- Return cell geometry.                *
 *      dt      -- Return local time step.              *
 *    dUdt      -- Return the solution residual.        *
 *                                                      *
 * Member operators                                     *
 *      S -- a 1D MHD solution                          *
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
class MHD1D_UniformMesh{
  private:
  public:
    MHD1D_pState      W;   // Primitive solution state.
    MHD1D_cState      U;   // Conserved solution state.
    Cell1D_Uniform    X;   // Cell geometry.
    double           dt;   // Local time step.
    MHD1D_cState   dUdt;   // Solution residual.
 	                   // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    MHD1D_UniformMesh(void) {
       W = MHD1D_W_REF; U = MHD1D_U_REF;
       X = Cell1D_Uniform_ONE; dt = ZERO;
       dUdt = MHD1D_U_ZERO;
    }

    MHD1D_UniformMesh(const MHD1D_UniformMesh &Soln) {
       W = Soln.W; U = Soln.U; X = Soln.X;
       dt = Soln.dt; dUdt = Soln.dUdt;
    }

    MHD1D_UniformMesh(const MHD1D_pState &W0,
		      const MHD1D_cState &U0,
		      const Cell1D_Uniform &X0) {
       W = W0; U = U0; X = X0;
       dt = ZERO; dUdt = MHD1D_U_ZERO;
    }

    MHD1D_UniformMesh(const MHD1D_pState &W0,
		      const Cell1D_Uniform &X0) {
       W = W0; U = W0.U(); X = X0;
       dt = ZERO; dUdt = MHD1D_U_ZERO;
    }

    MHD1D_UniformMesh(const MHD1D_cState &U0,
		      const Cell1D_Uniform &X0) {
       W = U0.W(); U = U0; X = X0;
       dt = ZERO; dUdt = MHD1D_U_ZERO;
    }
    
    /* Destructor. */
    // ~MHD1D_UniformMesh(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // MHD1D_UniformMesh operator = (const MHD1D_UniformMesh &Soln);
    // Use automatically generated assignment operator.

    /* Relational operators. */
    friend int operator ==(const MHD1D_UniformMesh &Soln1,
			   const MHD1D_UniformMesh &Soln2);
    friend int operator !=(const MHD1D_UniformMesh &Soln1,
			   const MHD1D_UniformMesh &Soln2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const MHD1D_UniformMesh &Soln);
    friend istream &operator >> (istream &in_file,
				 MHD1D_UniformMesh &Soln);
    
};

/********************************************************
 * MHD1D_UniformMesh -- Relational operators.           *
 ********************************************************/
inline int operator ==(const MHD1D_UniformMesh &Soln1,
		       const MHD1D_UniformMesh &Soln2) {
    return (Soln1.W == Soln2.W && Soln1.U == Soln2.U &&
	    Soln1.X == Soln2.X && Soln1.dt == Soln2.dt &&
	    Soln1.dUdt == Soln2.dUdt);
}

inline int operator !=(const MHD1D_UniformMesh &Soln1,
		       const MHD1D_UniformMesh &Soln2) {
    return (Soln1.W != Soln2.W || Soln1.U != Soln2.U ||
	    Soln1.X != Soln2.X || Soln1.dt != Soln2.dt ||
	    Soln1.dUdt != Soln2.dUdt);
}

/********************************************************
 * MHD1D_UniformMesh -- Input-output operators.         *
 ********************************************************/
inline ostream &operator << (ostream &out_file,
			     const MHD1D_UniformMesh &Soln) {
    out_file << Soln.X << Soln.W;
    out_file.setf(ios::scientific);
    out_file << " " << Soln.W.T();
    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     MHD1D_UniformMesh &Soln) {
    in_file >> Soln.X >> Soln.W;
    Soln.U = Soln.W.U();
    return (in_file);
}

/********************************************************
 * MHD1D_UniformMesh -- External subroutines.           *
 ********************************************************/

extern MHD1D_UniformMesh* Allocate(MHD1D_UniformMesh *Soln_ptr,
				   const int Number_of_Cells);

extern MHD1D_UniformMesh* Deallocate(MHD1D_UniformMesh *Soln_ptr,
				     const int Number_of_Cells);

extern void Output_Gnuplot(MHD1D_UniformMesh *Soln,
  		           const int Number_of_Cells,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Output_Tecplot(MHD1D_UniformMesh *Soln,
  		           const int Number_of_Cells,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Grid(MHD1D_UniformMesh *Soln,
                 const double &xmin,
		 const double &xmax,
		 const int Number_of_Cells);

extern void ICs(MHD1D_UniformMesh *Soln,
                char *gas_ptr,
		const int i_ICtype,
		const int Number_of_Cells);

extern double CFL(MHD1D_UniformMesh *Soln,
  		  const int Number_of_Cells);

extern int dUdt_explicitEuler_upwind(MHD1D_UniformMesh *Soln,
    	                             const int Number_of_Cells,
			             double &dtMin,
				     const double &CFL_Number,
				     const int Flux_Function_Type,
			             const int Local_Time_Stepping);

extern int dUdt_LaxFriedrichs(MHD1D_UniformMesh *Soln,
    	                      const int Number_of_Cells,
			      double &dtMin,
		              const double &CFL_Number,
			      const int Local_Time_Stepping);

/********************************************************
 * MHD1D_UniformMesh -- Solvers.                        *
 ********************************************************/

extern int MHD1DSolver(char *Input_File_Name_ptr,
                       int batch_flag);

#endif /* _MHD1D_INCLUDED  */
