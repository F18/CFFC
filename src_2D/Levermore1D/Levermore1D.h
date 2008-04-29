/*!\file Levermore1D.h
  \brief Header file defining 1D Levermore Solution Classes. */

#ifndef _LEVERMORE1D_INCLUDED
#define _LEVERMORE1D_INCLUDED

/* Include 1D Levermore state and cell header files. */

#ifndef _LEVERMORE1D_STATE_INCLUDED
#include "Levermore1DState.h"
#endif // _LEVERMORE1D_STATE_INCLUDED

#ifndef _CELL1D_INCLUDED
#include "../Grid/Cell1D.h"
#endif // _CELL1D_INCLUDED

#ifndef _CFD1DINPUT_INCLUDED
#include "../CFD/CFD1DInput.h"
#endif // _CFD1DINPUT_INCLUDED

#define   PREDICTED_MOMENT_NUMBER (Levermore1D_Vector::get_length())

/* Define the classes. */

/******************************************************//**
 * Class: Levermore1D_UniformMesh
 *
 * @brief Class definition of the 1D Levermore computational cell (i.e. solution + geometry)
 *
 * Member functions
 * -     W       -- Return primitive solution state.
 * -     U       -- Return conserved solution state.
 * -     A       -- Return closure weights.
 * -   dUdA      -- Return Hessian of density potential
 * - dUdA_inv    -- Return inverse Hessian of density potential
 * -   dFdA      -- Return Hessian of flux potential
 * -   dFdU      -- Return Flux Jacobian
 * -     X       -- Return cell geometry.
 * -     dt      -- Return local time step.
 * -   dUdt      -- Return the solution residual.
 * -   dUdx      -- Return the unlimited solution
 * -                gradient.
 * -    phi      -- Return the solution slope limiters.
 * -     Uo      -- Return initial solution state.
 * -     Uo      -- Return initial weights.
 * - Nghost      -- Number of ghost cells (!= 1 for
 * -                high-order)
 *
 * - cout << S; (output function)
 * - cin  >> S; (input function)
 *
 ********************************************************/
class Levermore1D_UniformMesh{
public:

  Levermore1D_pState    W;   //!< Primitive solution state.
  Levermore1D_cState    U;   //!< Conserved solution state.
  Levermore1D_weights   A;   //!< Closure weights
  DenseMatrix        dUdA;   //!< Hessian of density potential.
  DenseMatrix    dUdA_inv;   //!< Inverse Hessian of density potential.
  DenseMatrix        dFdA;   //!< Hessian of flux potential.
  DenseMatrix        dFdU;   //!< Flux Jacobian
  ColumnVector     Lambda;   //!< Vector containing Eigenvalues.
  double       lambda_max;   //!< Max eigenvalue.
  double       lambda_min;   //!< Min eigenvalue.
  Cell1D_Uniform        X;   //!< Cell geometry.
  double               dt;   //!< Local time step.
  Levermore1D_cState dUdt;   //!< Solution residual.
  //Levermore1D_cState TotaldUdt; //!< Solution residual. Used for RK4
  Levermore1D_pState dWdx;   //!< Unlimited solution gradient.
  Levermore1D_pState  phi;   //!< Solution slope limiter.
  Levermore1D_cState   Uo;   //!< Initial solution state.
  Levermore1D_weights  Ao;   //!< Initial weights.
  double         detector;   //!< The detector to determine syncronization of U and A
  double       detector_l;   //!< The detector to determine syncronization of U and A at left boundary
  double       detector_r;   //!< The detector to determine syncronization of U and A at right boundary
  double predicted_moment;   //!< Predicted moment used in detector
  double predicted_moment_l;   //!< Predicted moment used in detector at left boundary
  double predicted_moment_r;   //!< Predicted moment used in detector at right boundary
  int   number_of_resyncs;   //!< The number of times a resync was done
  Levermore1D_cState Ul_old, Ur_old;  //!< Remember left and right reconstructed states
  Levermore1D_weights Al_old, Ar_old; //!< Remember left and right reconstructed weights
  // Made public so can access them.

  int Nghost;            //!< Number of ghost cells(!= 1 for high-order)
  int ICl, ICu;	         //!< Indexes (Start & End)

  /* Creation, copy, and assignment constructors. */
  Levermore1D_UniformMesh(void);
  Levermore1D_UniformMesh(const Levermore1D_UniformMesh &Soln);
  Levermore1D_UniformMesh(const Levermore1D_pState &W0, const Levermore1D_cState &U0,
			  const Levermore1D_weights &A0, const Cell1D_Uniform &X0);
  Levermore1D_UniformMesh(const Levermore1D_pState &W0, const Cell1D_Uniform &X0);
  Levermore1D_UniformMesh(const Levermore1D_cState &U0, const Cell1D_Uniform &X0);

  /* Field access */
  const double & CellCenter(void) const {return X.x;} //!< return cell center
  const double & CellDelta (void) {return X.dx;} //!< return cell delta

  /* Primitive variables */
  const Levermore1D_pState & CellSolutionPrimVar(void) const {return W;}
  Levermore1D_pState & CellSolutionPrimVar(void) {return W;} //!< return the primitive solution state
  const double & CellSolutionPrimVar(int & VarPosition) const { return W[VarPosition];}
  double & CellSolutionPrimVar(int & VarPosition){ return W[VarPosition];} //!< return component of the primitive solution state
  /* Conservative variables */
  const Levermore1D_cState & CellSolutionConsVar(void) const {return U;}
  Levermore1D_cState & CellSolutionConsVar(void) {return U;} //!< return the conserved solution state
  const double & CellSolutionConsVar(int & VarPosition) const { return U[VarPosition];}
  double & CellSolutionConsVar(int & VarPosition) { return U[VarPosition];} //!< return component of the conserved solution state
  /* Closure Weights */
  const Levermore1D_weights & CellSolutionWeights(void) const {return A;}
  Levermore1D_weights & CellSolutionWeights(void) {return A;} //!< return the characteristic solution state
  const double & CellSolutionWeights(int & VarPosition) const { return A[VarPosition];}
  double & CellSolutionWeights(int & VarPosition) { return A[VarPosition];} //!< return component of the characteristic solution state
//  /* High-order variables */
//  HighOrderType & CellHighOrder(void) { return HO_dWdx; } //!< return the high-order variable for the primitive solution state
//  const HighOrderType & CellHighOrder(void) const { return HO_dWdx; }

  /* Operating functions */
//  double SolutionAtCoordinates_PWL (const double & X_Coord, const unsigned parameter) ;

//  /* Relational operators. */
//  friend int operator ==(const Levermore1D_UniformMesh &Soln1,
//			 const Levermore1D_UniformMesh &Soln2);
//  friend int operator !=(const Levermore1D_UniformMesh &Soln1,
//			 const Levermore1D_UniformMesh &Soln2);

  /* Set gas state */
  void set_state(const Levermore1D_pState &W0) {
    W = W0; U = Levermore1D_cState(W); A = Levermore1D_weights(U);
    detector = 0.0; detector_l = 0.0; detector_r = 0.0;
    Ul_old = U; Ur_old = U; Al_old = A; Ar_old = A;
    reset_predicted_moments();
  }
  void set_state(const Levermore1D_cState &U0) {
    U = U0; W = Levermore1D_pState(U); A = Levermore1D_weights(U);
    detector = 0.0; detector_l = 0.0; detector_r = 0.0;
    Ul_old = U; Ur_old = U; Al_old = A; Ar_old = A;
    reset_predicted_moments();
  }
//  void set_state(const Levermore1D_weights &A0) {
//    A = A0; U = Levermore1D_cState(A); W = Levermore1D_pState(U);
//  }
  void set_state(const Levermore1D_pState &W0,
		 const Levermore1D_cState &U0,
		 const Levermore1D_weights &A0) {
    W = W0; U = U0; A = A0;
    detector = 0.0; detector_l = 0.0; detector_r = 0.0;
    Ul_old = U; Ur_old = U; Al_old = A; Ar_old = A;
    reset_predicted_moments();
  }

  void calculate_Hessians() {
    double us(U[2]/U[1]);
    dUdA = U.d2hda2(A,us);
    dFdA = U.d2jda2(A,us);
    dUdA_inv = dUdA.inverse();
    dFdU = dFdA*dUdA_inv;
    Lambda = dFdU.eigenvalues();
    lambda_max = Lambda(0);
    lambda_min = Lambda(0);
    for(int i=1; i<Levermore1D_Vector::get_length();++i) {
      lambda_max = max(lambda_max,Lambda(i));
      lambda_min = min(lambda_min,Lambda(i));
    }
  }

  void reset_predicted_moment() {
    double us = U[2]/U[1];
    predicted_moment = U.moment(PREDICTED_MOMENT_NUMBER,A,us);
  }

  void reset_predicted_moment_l() {
    double us = Ul_old[2]/Ur_old[1];
    predicted_moment_l = Ul_old.moment(PREDICTED_MOMENT_NUMBER,Al_old,us);
  }

  void reset_predicted_moment_r() {
    double us = Ur_old[2]/Ur_old[1];
    predicted_moment_r = Ur_old.moment(PREDICTED_MOMENT_NUMBER,Ar_old,us);
  }
  void reset_predicted_moments() {
    reset_predicted_moment();
    reset_predicted_moment_l();
    reset_predicted_moment_r();
  }

  void update_predicted_moment(const ColumnVector &Update) {
    static Levermore1D_Vector temp;
    double us = U[2]/U[1];
    for(int i=1; i<=Levermore1D_Vector::get_length(); ++i) {
      temp[i] = U.moment(i+PREDICTED_MOMENT_NUMBER-1,A,us);
    }
    predicted_moment += temp*Update;
  }

  void update_predicted_moment_l(const ColumnVector &Update) {
    static Levermore1D_Vector temp;
    double us = Ul_old[2]/Ul_old[1];
    for(int i=1; i<=Levermore1D_Vector::get_length(); ++i) {
      temp[i] = Ul_old.moment(i+PREDICTED_MOMENT_NUMBER-1,Al_old,us);
    }
    predicted_moment_l += temp*Update;
  }

  void update_predicted_moment_r(const ColumnVector &Update) {
    static Levermore1D_Vector temp;
    double us = Ur_old[2]/Ur_old[1];
    for(int i=1; i<=Levermore1D_Vector::get_length(); ++i) {
      temp[i] = Ur_old.moment(i+PREDICTED_MOMENT_NUMBER-1,Ar_old,us);
    }
    predicted_moment_r += temp*Update;
  }

  void calculate_detector() {
    detector = U.detector_value(A,predicted_moment);
  }

  void calculate_detector_l() {
    detector_l = Ul_old.detector_value(Al_old,predicted_moment_l);
  }

  void calculate_detector_r() {
    detector_r = Ur_old.detector_value(Ar_old,predicted_moment_r);
  }

  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file,
			       const Levermore1D_UniformMesh &Soln);
  friend istream &operator >> (istream &in_file,
			       Levermore1D_UniformMesh &Soln);

private:
//  HighOrderType   HO_dWdx; //!< High-order derivatives container for the primitive solution state

};

/*******************************************************
 *  Levermore1D_UniformMesh::MemberFunctions()         *
 ******************************************************/
/* Constructor */
inline Levermore1D_UniformMesh::Levermore1D_UniformMesh(void){
  W.zero(); U.zero(); A.zero();
  X = Cell1D_Uniform_ONE; dt = ZERO;
  dUdt.zero(); dWdx.zero(); phi.zero();
  Uo.zero(); Ao.zero();
  Ul_old = U; Ur_old = U; Al_old = A; Ar_old = A;
  detector = 0.0; detector_l = 0.0; detector_r = 0.0;
  number_of_resyncs = 0;
  reset_predicted_moments();
}

inline Levermore1D_UniformMesh::Levermore1D_UniformMesh(const Levermore1D_UniformMesh &Soln) {
  W = Soln.W; U = Soln.U; A = Soln.A; X = Soln.X;
  dt = Soln.dt; dUdt = Soln.dUdt; dWdx = Soln.dWdx;
  phi = Soln.phi; Uo = Soln.Uo; Ao = Soln.Ao;
  Nghost = Soln.Nghost; ICl = Soln.ICl; ICu = Soln.ICu;
  lambda_max = Soln.lambda_max;
  lambda_min = Soln.lambda_min;
  detector = Soln.detector;
  detector_l = Soln.detector_l;
  detector_r = Soln.detector_r;
  number_of_resyncs = Soln.number_of_resyncs;
  predicted_moment = Soln.predicted_moment;
  predicted_moment_l = Soln.predicted_moment_l;
  predicted_moment_r = Soln.predicted_moment_r;
  Ul_old = U; Ur_old = U; Al_old = A; Ar_old = A;
}

inline Levermore1D_UniformMesh::Levermore1D_UniformMesh(const Levermore1D_pState &W0,
							const Levermore1D_cState &U0,
							const Levermore1D_weights &A0,
							const Cell1D_Uniform &X0) {
  W = W0; U = U0; A = A0; X = X0;
  dt = ZERO; dUdt.zero();
  dWdx.zero(); phi.zero(); Uo.zero(); Ao.zero();
  Ul_old = U; Ur_old = U; Al_old = A; Ar_old = A;
  detector = 0.0; detector_l = 0.0; detector_r = 0.0;
  number_of_resyncs = 0;
  reset_predicted_moments();
}

inline Levermore1D_UniformMesh::Levermore1D_UniformMesh(const Levermore1D_pState &W0,
							const Cell1D_Uniform &X0) {
  W = W0; U = Levermore1D_cState(W0); A = Levermore1D_weights(W0);
  X = X0; dt = ZERO; dUdt.zero();
  dWdx.zero(); phi.zero(); Uo.zero(); Ao.zero();
  Ul_old = U; Ur_old = U; Al_old = A; Ar_old = A;
  detector = 0.0; number_of_resyncs = 0;
  number_of_resyncs = 0;
  reset_predicted_moments();
}

inline Levermore1D_UniformMesh::Levermore1D_UniformMesh(const Levermore1D_cState &U0,
							const Cell1D_Uniform &X0) {
  W = Levermore1D_pState(U0); U = U0; A = Levermore1D_weights(U0);
  X = X0; dt = ZERO; dUdt.zero();
  dWdx.zero(); phi.zero(); Uo.zero(); Ao.zero();
  Ul_old = U; Ur_old = U; Al_old = A; Ar_old = A;
  detector = 0.0; detector_l = 0.0; detector_r = 0.0;
  number_of_resyncs = 0;
  reset_predicted_moments();
}

////! Return the solution of the piecewise limited linear reconstruction at the coordinate X_Coord,
////  for the required parameter.
//inline double Levermore1D_UniformMesh::SolutionAtCoordinates_PWL (const double & X_Coord, const unsigned parameter) {
//  SOLN_pSTATE W_Temp;
//  double Distance(X_Coord - X.x);
//  W_Temp = W + (phi^dWdx)*Distance;
//  return W_Temp[parameter];
//}

///********************************************************
// * Levermore1D_UniformMesh -- Relational operators.         *
// ********************************************************/
//inline int operator ==(const Levermore1D_UniformMesh &Soln1,
//		       const Levermore1D_UniformMesh &Soln2) {
//  return (Soln1.W == Soln2.W && Soln1.U == Soln2.U &&
//	  Soln1.X == Soln2.X && Soln1.dt == Soln2.dt &&
//	  Soln1.dUdt == Soln2.dUdt && Soln1.dWdx == Soln2.dWdx &&
//	  Soln1.phi == Soln2.phi && Soln1.Uo == Soln2.Uo &&
//	  Soln1.CellSolutionCharactVar() == Soln2.CellSolutionCharactVar() &&
//	  Soln1.CellHighOrder() == Soln2.CellHighOrder() );
//}
//
//inline int operator !=(const Levermore1D_UniformMesh &Soln1,
//		       const Levermore1D_UniformMesh &Soln2) {
//  return (Soln1.W != Soln2.W || Soln1.U != Soln2.U ||
//	  Soln1.X != Soln2.X || Soln1.dt != Soln2.dt ||
//	  Soln1.dUdt != Soln2.dUdt || Soln1.dWdx != Soln2.dWdx ||
//	  Soln1.phi != Soln2.phi || Soln1.Uo != Soln2.Uo ||
//	  Soln1.CellSolutionCharactVar() != Soln2.CellSolutionCharactVar() ||
//	  Soln1.CellHighOrder() != Soln2.CellHighOrder());
//}

/********************************************************
 * Levermore1D_UniformMesh -- Input-output operators.       *
 ********************************************************/
inline ostream &operator << (ostream &out_file,
			     const Levermore1D_UniformMesh &Soln) {
  out_file << Soln.X << Soln.W;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Levermore1D_UniformMesh &Soln) {
  in_file >> Soln.X >> Soln.W;
  Soln.U = Levermore1D_cState(Soln.W);
  Soln.A = Levermore1D_weights(Soln.W);
  return (in_file);
}

/********************************************************
 * Levermore1D_UniformMesh -- External subroutines.         *
 ********************************************************/

extern Levermore1D_UniformMesh* Allocate(Levermore1D_UniformMesh *Soln_ptr,
					 const CFD1D_Input_Parameters &IP);

extern Levermore1D_UniformMesh* Deallocate(Levermore1D_UniformMesh *Soln_ptr);

extern void Output_Gnuplot(Levermore1D_UniformMesh *Soln,
  		           const int Number_of_Cells,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

extern void Output_Tecplot(Levermore1D_UniformMesh *Soln,
  		           const int Number_of_Cells,
		           const int Number_of_Time_Steps,
                           const double &Time,
	                   ostream &out_file);

void Output_Tecplot(Levermore1D_UniformMesh *Soln,
		    const CFD1D_Input_Parameters &IP,
		    const int Number_of_Time_Steps,
		    const double &Time,
		    ostream &out_file);

//void Output_Tecplot_HighOrder(Levermore1D_UniformMesh *Soln,
//			      const int Number_of_Cells,
//			      const int Number_of_Time_Steps,
//			      const double &Time,
//			      ostream &out_file);

extern void Grid(Levermore1D_UniformMesh *Soln,
                 const double &xmin,
		 const double &xmax,
		 const int Number_of_Cells);

extern void ICs(Levermore1D_UniformMesh *Soln,
                char *gas_ptr,
		const int i_ICtype,
		const int Number_of_Cells,
		CFD1D_Input_Parameters &IP);

extern double CFL(Levermore1D_UniformMesh *Soln,
  		  const int Number_of_Cells);

extern void Linear_Reconstruction_MUSCL(Levermore1D_UniformMesh *Soln,
                                        const int Number_of_Cells,
					const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Levermore1D_UniformMesh *Soln,
                                             const int Number_of_Cells,
					     const int Limiter);

//extern void Linear_Reconstruction_Characteristic(Levermore1D_UniformMesh *Soln,
//                                                 const int Number_of_Cells,
//					         const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Levermore1D_UniformMesh *Soln,
                                               const int Number_of_Cells,
					       const int Limiter);

extern int dUdt_explicitEuler_upwind(Levermore1D_UniformMesh *Soln,
					 const int Number_of_Cells,
					 double &dtMin,
					 const double &CFL_Number,
					 const int Flux_Function_Type,
					 const int Local_Time_Stepping);

extern int dUdt_LaxFriedrichs(Levermore1D_UniformMesh *Soln,
    	                      const int Number_of_Cells,
			      double &dtMin,
			      const double &CFL_Number,
			      const int Local_Time_Stepping);

extern int dUdt_LaxWendroff(Levermore1D_UniformMesh *Soln,
    	                    const int Number_of_Cells,
			    double &dtMin,
			    const double &CFL_Number,
			    const int Local_Time_Stepping);

extern int dUdt_MacCormack(Levermore1D_UniformMesh *Soln,
    	                   const int Number_of_Cells,
			   double &dtMin,
			   const double &CFL_Number,
			   const int Local_Time_Stepping);

extern int dUdt_Hancock(Levermore1D_UniformMesh *Soln,
    	                const int Number_of_Cells,
			double &dtMin,
			const double &CFL_Number,
                        const int Limiter_Type,
                        const int Flux_Function_Type,
			const int Local_Time_Stepping);

extern int dUdt_2stage_2ndOrder_upwind(Levermore1D_UniformMesh *Soln,
				       const CFD1D_Input_Parameters &IP,
			               double &dtMin,
				       const double &CFL_Number,
                                       const int Reconstruction_Type,
                                       const int Limiter_Type,
				       const int Flux_Function_Type,
			               const int Local_Time_Stepping);

void LimitedLinearReconstructionOverDomain(Levermore1D_UniformMesh *Soln, const CFD1D_Input_Parameters &IP);

/********************************************************
 * Levermore1D_UniformMesh -- Solvers.                      *
 ********************************************************/

extern int Levermore1DSolver(char *Input_File_Name_ptr,
			     int batch_flag);

#endif /* _LEVERMORE1D_INCLUDED  */
