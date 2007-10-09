/* Euler1D.h:  Header file defining 1D Euler Solution Classes. */

#ifndef _EULER1D_INCLUDED
#define _EULER1D_INCLUDED

/* Include 1D Euler state and cell header files. */

#ifndef _EULER1D_STATE_INCLUDED
#include "Euler1DState.h"
#endif // _EULER1D_STATE_INCLUDED

#ifndef _CELL1D_INCLUDED
#include "../Grid/Cell1D.h"
#endif // _CELL1D_INCLUDED

#include "../CFD/CFD1DInput.h"

/* Define the classes. */

/********************************************************
 * Class: Euler1D_UniformMesh                           *
 *                                                      *
 * Member functions                                     *
 *      W       -- Return primitive solution state.     *
 *      U       -- Return conserved solution state.     *
 *      X       -- Return cell geometry.                *
 *      dt      -- Return local time step.              *
 *    dUdt      -- Return the solution residual.        *
 *    dUdx      -- Return the unlimited solution        *
 *                 gradient.                            *
 *     phi      -- Return the solution slope limiters.  *
 *      Uo      -- Return initial solution state.       *
 *                                                      *
 * Member operators                                     *
 *      S -- a 1D Euler solution                        *
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
class Euler1D_UniformMesh{
  private:
  public:
    Euler1D_pState    W;   // Primitive solution state.
    Euler1D_cState    U;   // Conserved solution state.
    Euler1D_pState CharactVar;    // Characteristic solution state.
    Cell1D_Uniform    X;   // Cell geometry.
    double           dt;   // Local time step.
    Euler1D_cState dUdt;   // Solution residual.
    Euler1D_pState dWdx;   // Unlimited solution gradient.
    Euler1D_pState  phi;   // Solution slope limiter.
    Euler1D_cState   Uo;   // Initial solution state.
                           // Made public so can access them.
    int Nghost;     // Number of ghost cells(!= 1 for high-order)

    /* Field access */
    const double & CellCenter(void) const {return X.x;}
    const double & CellDelta (void) {return X.dx;}

    /* Primitive variables */
    const Euler1D_pState & CellSolutionPrimVar(void) const {return W;}
    Euler1D_pState & CellSolutionPrimVar(void) {return W;}
    const double & CellSolutionPrimVar(int & VarPosition) const { return W[VarPosition];}
    double & CellSolutionPrimVar(int & VarPosition){ return W[VarPosition];}
    /* Conservative variables */
    const Euler1D_cState & CellSolutionConsVar(void) const {return U;}
    Euler1D_cState & CellSolutionConsVar(void) {return U;}
    const double & CellSolutionConsVar(int & VarPosition) const { return U[VarPosition];}
    double & CellSolutionConsVar(int & VarPosition) { return U[VarPosition];}
    /* Characteristic variables */
    const Euler1D_pState & CellSolutionCharactVar(void) const {return CharactVar;}
    Euler1D_pState & CellSolutionCharactVar(void) {return CharactVar;}
    const double & CellSolutionCharactVar(int & VarPosition) const { return CharactVar[VarPosition];}
    double & CellSolutionCharactVar(int & VarPosition) { return CharactVar[VarPosition];}
  
	      
    /* Creation, copy, and assignment constructors. */
    Euler1D_UniformMesh(void) {
       W = Euler1D_W_STDATM; U = Euler1D_U_STDATM;
       X = Cell1D_Uniform_ONE; dt = ZERO;
       dUdt = Euler1D_U_VACUUM;
       dWdx = Euler1D_W_VACUUM; phi = Euler1D_W_VACUUM;
       Uo = Euler1D_U_VACUUM;
    }

    Euler1D_UniformMesh(const Euler1D_UniformMesh &Soln) {
       W = Soln.W; U = Soln.U; X = Soln.X;
       dt = Soln.dt; dUdt = Soln.dUdt; dWdx = Soln.dWdx;
       phi = Soln.phi; Uo = Soln.Uo;
    }

    Euler1D_UniformMesh(const Euler1D_pState &W0,
			const Euler1D_cState &U0,
			const Cell1D_Uniform &X0) {
       W = W0; U = U0; X = X0;
       dt = ZERO; dUdt = Euler1D_U_VACUUM;
       dWdx = Euler1D_W_VACUUM; phi = Euler1D_W_VACUUM;
       Uo = Euler1D_U_VACUUM;
    }

    Euler1D_UniformMesh(const Euler1D_pState &W0,
			const Cell1D_Uniform &X0) {
       W = W0; U = W0.U(); X = X0;
       dt = ZERO; dUdt = Euler1D_U_VACUUM;
       dWdx = Euler1D_W_VACUUM; phi = Euler1D_W_VACUUM;
       Uo = Euler1D_U_VACUUM;
    }

    Euler1D_UniformMesh(const Euler1D_cState &U0,
			const Cell1D_Uniform &X0) {
       W = U0.W(); U = U0; X = X0;
       dt = ZERO; dUdt = Euler1D_U_VACUUM;
       dWdx = Euler1D_W_VACUUM; phi = Euler1D_W_VACUUM;
       Uo = Euler1D_U_VACUUM;
    }
    
    /* Destructor. */
    // ~Euler1D_UniformMesh(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // Euler1D_UniformMesh operator = (const Euler1D_UniformMesh &Soln);
    // Use automatically generated assignment operator.

    /* Binary arithmetic operators. */
//    friend Euler1D_UniformMesh operator +(const Euler1D_UniformMesh &Soln1,
//					  const Euler1D_UniformMesh &Soln2);
//    friend Euler1D_UniformMesh operator -(const Euler1D_UniformMesh &Soln1,
//					  const Euler1D_UniformMesh &Soln2);

    /* Unary arithmetic operators. */
//    friend Euler1D_UniformMesh operator +(const Euler1D_UniformMesh &Soln);
//    friend Euler1D_UniformMesh operator -(const Euler1D_UniformMesh &Soln);

    /* Shortcut arithmetic operators. */
//    friend Euler1D_UniformMesh &operator +=(Euler1D_UniformMesh &Soln1,
//					    const Euler1D_UniformMesh &Soln2);
//    friend Euler1D_UniformMesh &operator -=(Euler1D_UniformMesh &Soln1,
//					    const Euler1D_UniformMesh &Soln2);
    
    /* Relational operators. */
    friend int operator ==(const Euler1D_UniformMesh &Soln1,
			   const Euler1D_UniformMesh &Soln2);
    friend int operator !=(const Euler1D_UniformMesh &Soln1,
			   const Euler1D_UniformMesh &Soln2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Euler1D_UniformMesh &Soln);
    friend istream &operator >> (istream &in_file,
				 Euler1D_UniformMesh &Soln);
    
};

/********************************************************
 * Euler1D_UniformMesh -- Relational operators.         *
 ********************************************************/
inline int operator ==(const Euler1D_UniformMesh &Soln1,
		       const Euler1D_UniformMesh &Soln2) {
    return (Soln1.W == Soln2.W && Soln1.U == Soln2.U &&
	    Soln1.X == Soln2.X && Soln1.dt == Soln2.dt &&
	    Soln1.dUdt == Soln2.dUdt && Soln1.dWdx == Soln2.dWdx &&
	    Soln1.phi == Soln2.phi && Soln1.Uo == Soln2.Uo);
}

inline int operator !=(const Euler1D_UniformMesh &Soln1,
		       const Euler1D_UniformMesh &Soln2) {
    return (Soln1.W != Soln2.W || Soln1.U != Soln2.U ||
	    Soln1.X != Soln2.X || Soln1.dt != Soln2.dt ||
	    Soln1.dUdt != Soln2.dUdt || Soln1.dWdx != Soln2.dWdx ||
	    Soln1.phi != Soln2.phi || Soln1.Uo == Soln2.Uo);
}

/********************************************************
 * Euler1D_UniformMesh -- Input-output operators.       *
 ********************************************************/
inline ostream &operator << (ostream &out_file,
			     const Euler1D_UniformMesh &Soln) {
    out_file << Soln.X << Soln.W;
    out_file.setf(ios::scientific);
    out_file << " " << Soln.W.T();
    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Euler1D_UniformMesh &Soln) {
    in_file >> Soln.X >> Soln.W;
    Soln.U = Soln.W.U();
    return (in_file);
}

/********************************************************
 * Euler1D_UniformMesh -- External subroutines.         *
 ********************************************************/

extern Euler1D_UniformMesh* Allocate(Euler1D_UniformMesh *Soln_ptr,
				     const int Number_of_Cells);

extern Euler1D_UniformMesh* Deallocate(Euler1D_UniformMesh *Soln_ptr,
				       const int Number_of_Cells);

extern void Output_Gnuplot(Euler1D_UniformMesh *Soln,
  		           const int Number_of_Cells,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Output_Tecplot(Euler1D_UniformMesh *Soln,
  		           const int Number_of_Cells,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Grid(Euler1D_UniformMesh *Soln,
                 const double &xmin,
		 const double &xmax,
		 const int Number_of_Cells);

extern void ICs(Euler1D_UniformMesh *Soln,
                char *gas_ptr,
		const int i_ICtype,
		const int Number_of_Cells);

extern double CFL(Euler1D_UniformMesh *Soln,
  		  const int Number_of_Cells);

extern void Linear_Reconstruction_MUSCL(Euler1D_UniformMesh *Soln,
                                        const int Number_of_Cells,
					const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Euler1D_UniformMesh *Soln,
                                             const int Number_of_Cells,
					     const int Limiter);

extern void Linear_Reconstruction_Characteristic(Euler1D_UniformMesh *Soln,
                                                 const int Number_of_Cells,
					         const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Euler1D_UniformMesh *Soln,
                                               const int Number_of_Cells,
					       const int Limiter);

extern int dUdt_explicitEuler_upwind(Euler1D_UniformMesh *Soln,
    	                             const int Number_of_Cells,
			             double &dtMin,
				     const double &CFL_Number,
				     const int Flux_Function_Type,
			             const int Local_Time_Stepping);

extern int dUdt_LaxFriedrichs(Euler1D_UniformMesh *Soln,
    	                      const int Number_of_Cells,
			      double &dtMin,
			      const double &CFL_Number,
			      const int Local_Time_Stepping);

extern int dUdt_LaxWendroff(Euler1D_UniformMesh *Soln,
    	                    const int Number_of_Cells,
			    double &dtMin,
			    const double &CFL_Number,
			    const int Local_Time_Stepping);

extern int dUdt_MacCormack(Euler1D_UniformMesh *Soln,
    	                   const int Number_of_Cells,
			   double &dtMin,
			   const double &CFL_Number,
			   const int Local_Time_Stepping);

extern int dUdt_Hancock(Euler1D_UniformMesh *Soln,
    	                const int Number_of_Cells,
			double &dtMin,
			const double &CFL_Number,
                        const int Limiter_Type,
                        const int Flux_Function_Type,
			const int Local_Time_Stepping);

extern int dUdt_2stage_2ndOrder_upwind(Euler1D_UniformMesh *Soln,
    	                               const int Number_of_Cells,
			               double &dtMin,
				       const double &CFL_Number,
                                       const int Reconstruction_Type,
                                       const int Limiter_Type,
				       const int Flux_Function_Type,
			               const int Local_Time_Stepping);

/********************************************************
 * Euler1D_UniformMesh -- Solvers.                      *
 ********************************************************/

extern int Euler1DSolver(char *Input_File_Name_ptr,
                         int batch_flag);

#endif /* _EULER1D_INCLUDED  */
