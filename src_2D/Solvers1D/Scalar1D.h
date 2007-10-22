/* Scalar1D.h:  Header file defining 1D Scalar Advection Solution Classes. */

#ifndef _SCALAR1D_INCLUDED
#define _SCALAR1D_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>

using namespace std;

/* Include CFD and 1D cell header files. */

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _CELL1D_INCLUDED
#include "../Grid/Cell1D.h"
#endif // _CELL1D_INCLUDED

#include "../CFD/CFD1DInput.h"

/* Define the classes. */

/********************************************************
 * Class: Scalar1D_UniformMesh                          *
 *                                                      *
 * Member functions                                     *
 *      u       -- Return solution.                     *
 *      a       -- Return advection velocity.           *
 *      X       -- Return cell geometry.                *
 *      dt      -- Return local time step.              *
 *    dudt      -- Return the solution residual.        *
 *    dudx      -- Return the unlimited solution        *
 *                 gradient.                            *
 *     phi      -- Return the solution slope limiters.  *
 *      uo      -- Return initial solution state.       *
 *                                                      *
 * Member operators                                     *
 *      S -- a 1D Scalar solution                       *
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
class Scalar1D_UniformMesh{
  private:
  public:
    double            u;   // Solution variable.
    double            a;   // Advection velocity.
    Cell1D_Uniform    X;   // Cell geometry.
    double           dt;   // Local time step.
    double         dudt;   // Solution residual.
    double         dudx;   // Unlimited solution gradient.
    double          phi;   // Solution slope limiter.
    double           uo;   // Initial solution state.
	                   // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Scalar1D_UniformMesh(void) {
       u = ZERO; a = ONE; X = Cell1D_Uniform_ONE; dt = ZERO;
       dudt = ZERO; dudx = ZERO; phi = ZERO; uo = ZERO;
    }

    Scalar1D_UniformMesh(const Scalar1D_UniformMesh &Soln) {
       u = Soln.u; a = Soln.a; X = Soln.X;
       dt = Soln.dt; dudt = Soln.dudt; dudx = Soln.dudx;
       phi = Soln.phi; uo = Soln.uo;
    }

    Scalar1D_UniformMesh(const double &u0,
			 const double &a0,
			 const Cell1D_Uniform &X0) {
       u = u0; a = a0; X = X0; dt = ZERO; dudt = ZERO; 
       dudx = ZERO; phi = ZERO; uo = ZERO;
    }

    /* Destructor. */
    // ~Scalar1D_UniformMesh(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // Scalar1D_UniformMesh operator = (const Scalar1D_UniformMesh &Soln);
    // Use automatically generated assignment operator.

    /* Relational operators. */
    friend int operator ==(const Scalar1D_UniformMesh &Soln1,
			   const Scalar1D_UniformMesh &Soln2);
    friend int operator !=(const Scalar1D_UniformMesh &Soln1,
			   const Scalar1D_UniformMesh &Soln2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Scalar1D_UniformMesh &Soln);
    friend istream &operator >> (istream &in_file,
				 Scalar1D_UniformMesh &Soln);
    
};

/********************************************************
 * Scalar1D_UniformMesh -- Relational operators.        *
 ********************************************************/
inline int operator ==(const Scalar1D_UniformMesh &Soln1,
		       const Scalar1D_UniformMesh &Soln2) {
    return (Soln1.u == Soln2.u && Soln1.a == Soln2.a &&
	    Soln1.X == Soln2.X && Soln1.dt == Soln2.dt &&
	    Soln1.dudt == Soln2.dudt && Soln1.dudx == Soln2.dudx &&
	    Soln1.phi == Soln2.phi && Soln1.uo == Soln2.uo);
}

inline int operator !=(const Scalar1D_UniformMesh &Soln1,
		       const Scalar1D_UniformMesh &Soln2) {
    return (Soln1.u != Soln2.u || Soln1.a != Soln2.a ||
	    Soln1.X != Soln2.X || Soln1.dt != Soln2.dt ||
	    Soln1.dudt != Soln2.dudt || Soln1.dudx != Soln2.dudx ||
	    Soln1.phi != Soln2.phi || Soln1.uo == Soln2.uo);
}

/********************************************************
 * Scalar1D_UniformMesh -- Input-output operators.      *
 ********************************************************/
inline ostream &operator << (ostream &out_file,
			     const Scalar1D_UniformMesh &Soln) {
    out_file << Soln.X;
    out_file.setf(ios::scientific);
    out_file << " " << Soln.u << " " << Soln.a;
    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Scalar1D_UniformMesh &Soln) {
    in_file >> Soln.X;
    in_file.setf(ios::skipws);
    in_file >> Soln.u >> Soln.a;
    in_file.unsetf(ios::skipws);
    return (in_file);
}

/********************************************************
 * Scalar1D_UniformMesh -- External subroutines.        *
 ********************************************************/

extern Scalar1D_UniformMesh* Allocate(Scalar1D_UniformMesh *Soln_ptr,
				      const int Number_of_Cells);

extern Scalar1D_UniformMesh* Deallocate(Scalar1D_UniformMesh *Soln_ptr,
					const int Number_of_Cells);

extern void Output_Gnuplot(Scalar1D_UniformMesh *Soln,
  		           const int Number_of_Cells,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Output_Tecplot(Scalar1D_UniformMesh *Soln,
  		           const int Number_of_Cells,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Grid(Scalar1D_UniformMesh *Soln,
                 const double &xmin,
		 const double &xmax,
		 const int Number_of_Cells);

extern void ICs(Scalar1D_UniformMesh *Soln,
		const int i_ICtype,
		const int Number_of_Cells);

extern double CFL(Scalar1D_UniformMesh *Soln,
  		  const int Number_of_Cells);

extern void Linear_Reconstruction_MUSCL(Scalar1D_UniformMesh *Soln,
                                        const int Number_of_Cells,
					const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Scalar1D_UniformMesh *Soln,
                                             const int Number_of_Cells,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Scalar1D_UniformMesh *Soln,
                                               const int Number_of_Cells,
					       const int Limiter);

extern int dudt_explicitEuler_upwind(Scalar1D_UniformMesh *Soln,
    	                             const int Number_of_Cells,
			             double &dtMin,
				     const double &CFL_Number,
			             const int Local_Time_Stepping);

extern int dudt_LaxFriedrichs(Scalar1D_UniformMesh *Soln,
    	                      const int Number_of_Cells,
			      double &dtMin,
			      const double &CFL_Number,
			      const int Local_Time_Stepping);

extern int dudt_LaxWendroff(Scalar1D_UniformMesh *Soln,
    	                    const int Number_of_Cells,
			    double &dtMin,
			    const double &CFL_Number,
			    const int Local_Time_Stepping);

extern int dudt_MacCormack(Scalar1D_UniformMesh *Soln,
    	                   const int Number_of_Cells,
			   double &dtMin,
			   const double &CFL_Number,
			   const int Local_Time_Stepping);

extern int dUdt_2stage_2ndOrder_upwind(Scalar1D_UniformMesh *Soln,
	                               const int Number_of_Cells,
			               double &dtMin,
				       const double &CFL_Number,
                                       const int Reconstruction_Type,
                                       const int Limiter_Type,
			               const int Local_Time_Stepping);

/********************************************************
 * Scalar1D_UniformMesh -- Solvers.                     *
 ********************************************************/

extern int Scalar1DSolver(char *Input_File_Name_ptr,
                          int batch_flag);

#endif /* _SCALAR1D_INCLUDED  */
