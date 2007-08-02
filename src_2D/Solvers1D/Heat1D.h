/* Heat1D.h:  Header file defining 1D Heat Equation Solution Classes. */

#ifndef _HEAT1D_INCLUDED
#define _HEAT1D_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>

using namespace std;

/* Include linear systems, CFD, and 1D cell header files. */

#ifndef _LINEARSYSTEMS_INCLUDED
#include "../Math/LinearSystems.h"
#endif // _LINEARSYSTEMS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _CELL1D_INCLUDED
#include "../Grid/Cell1D.h"
#endif // _CELL1D_INCLUDED

/* Define the classes. */

/********************************************************
 * Class: Heat1D_UniformMesh                            *
 *                                                      *
 * Member functions                                     *
 *      T       -- Return temperature.                  *
 *      qx      -- Return heat flux.                    *
 *      kappa   -- Return thermal conductivity.         *
 *      X       -- Return node geometry.                *
 *      dt      -- Return local time step.              *
 *    dTdt      -- Return the solution residual.        *
 *                                                      *
 * Member operators                                     *
 *      S -- a 1D heat equation solution                *
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
class Heat1D_UniformMesh{
  private:
  public:
    double             T;   // Temperature.
    double            qx;   // Heat flux.
    double         kappa;   // Thermal conductivity.
    Cell1D_Uniform     X;   // Node geometry.
    double            dt;   // Local time step.
    double          dTdt;   // Solution residual.
	                    // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Heat1D_UniformMesh(void) {
       T = ZERO; qx = ZERO; kappa = ONE; X = Cell1D_Uniform_ONE; dt = ZERO;
       dTdt = ZERO;
    }

    Heat1D_UniformMesh(const Heat1D_UniformMesh &Soln) {
       T = Soln.T; qx = Soln.qx; kappa = Soln.kappa; X = Soln.X;
       dt = Soln.dt; dTdt = Soln.dTdt;
    }

    Heat1D_UniformMesh(const double &T0,
                       const double &qx0,
		       const double &kappa0,
		       const Cell1D_Uniform &X0) {
       T = T0; qx = qx0; kappa = kappa0; X = X0;
       dt = ZERO; dTdt = ZERO;
    }

    /* Destructor. */
    // ~Heat1D_UniformMesh(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // Heat1D_UniformMesh operator = (const Heat1D_UniformMesh &Soln);
    // Use automatically generated assignment operator.

    /* Relational operators. */
    friend int operator ==(const Heat1D_UniformMesh &Soln1,
			   const Heat1D_UniformMesh &Soln2);
    friend int operator !=(const Heat1D_UniformMesh &Soln1,
			   const Heat1D_UniformMesh &Soln2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Heat1D_UniformMesh &Soln);
    friend istream &operator >> (istream &in_file,
				 Heat1D_UniformMesh &Soln);
    
};

/********************************************************
 * Heat1D_UniformMesh -- Relational operators.          *
 ********************************************************/
inline int operator ==(const Heat1D_UniformMesh &Soln1,
		       const Heat1D_UniformMesh &Soln2) {
    return (Soln1.T == Soln2.T && Soln1.qx == Soln2.qx &&
	    Soln1.kappa == Soln2.kappa &&
	    Soln1.X == Soln2.X && Soln1.dt == Soln2.dt &&
	    Soln1.dTdt == Soln2.dTdt);
}

inline int operator !=(const Heat1D_UniformMesh &Soln1,
		       const Heat1D_UniformMesh &Soln2) {
    return (Soln1.T != Soln2.T || Soln1.qx != Soln2.qx ||
	    Soln1.kappa != Soln2.kappa ||
	    Soln1.X != Soln2.X || Soln1.dt != Soln2.dt ||
	    Soln1.dTdt != Soln2.dTdt);
}

/********************************************************
 * Heat1D_UniformMesh -- Input-output operators.        *
 ********************************************************/
inline ostream &operator << (ostream &out_file,
			     const Heat1D_UniformMesh &Soln) {
    out_file << Soln.X;
    out_file.setf(ios::scientific);
    out_file << " " << Soln.T << " " << Soln.qx << " " << Soln.kappa;
    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Heat1D_UniformMesh &Soln) {
    in_file >> Soln.X;
    in_file.setf(ios::skipws);
    in_file >> Soln.T >> Soln.qx >> Soln.kappa;
    in_file.unsetf(ios::skipws);
    return (in_file);
}

/********************************************************
 * Heat1D_UniformMesh -- External subroutines.          *
 ********************************************************/

extern Heat1D_UniformMesh* Allocate(Heat1D_UniformMesh *Soln_ptr,
				    const int Number_of_Cells);

extern Heat1D_UniformMesh* Deallocate(Heat1D_UniformMesh *Soln_ptr,
				      const int Number_of_Cells);

extern void Output_Gnuplot(Heat1D_UniformMesh *Soln,
  		           const int Number_of_Nodes,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Output_Tecplot(Heat1D_UniformMesh *Soln,
  		           const int Number_of_Nodes,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Grid(Heat1D_UniformMesh *Soln,
                 const double &xmin,
		 const double &xmax,
		 const int Number_of_Nodes);

extern void ICs(Heat1D_UniformMesh *Soln,
		const int i_ICtype,
                const double &Kappa,
		const int Number_of_Nodes);

extern void AnalyticSoln(Heat1D_UniformMesh *Soln,
                         const int i_ICtype,
                         const double &Kappa,
                         const double &Time,
		         const int Number_of_Nodes);

extern void HeatFlux(Heat1D_UniformMesh *Soln,
                     const int Number_of_Nodes);

extern double Rstability(Heat1D_UniformMesh *Soln,
  		         const int Number_of_Nodes);

extern int dTdt_SimpleExplicit(Heat1D_UniformMesh *Soln,
    	                       const int Number_of_Nodes,
			       double &dtMin,
			       const double &CFL_Number,
			       const int Local_Time_Stepping);

extern int dTdt_SimpleImplicit(Heat1D_UniformMesh *Soln,
    	                       const int Number_of_Nodes,
			       double &dtMin,
			       const double &CFL_Number,
			       const int Local_Time_Stepping,
                               TriDiagonalSystemLinEqs &Ax_equal_b);

extern int dTdt_CrankNicolson(Heat1D_UniformMesh *Soln,
    	                      const int Number_of_Nodes,
			      double &dtMin,
			      const double &CFL_Number,
                              const int Local_Time_Stepping,
                              TriDiagonalSystemLinEqs &Ax_equal_b);

extern int dTdt_ADE(Heat1D_UniformMesh *Soln,
    	            const int Number_of_Nodes,
		    double &dtMin,
		    const double &CFL_Number,
		    const int Local_Time_Stepping);

/********************************************************
 * Heat1D_UniformMesh -- Solvers.                       *
 ********************************************************/

extern int Heat1DSolver(char *Input_File_Name_ptr,
                        int batch_flag);

#endif /* _HEAT1D_INCLUDED  */
