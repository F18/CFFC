/************************************************************************
  15Step.h:  Header file defining a variety of fortran functions
             required to compute the reaction rates of a 15 step
             CH4 mechanism.
************************************************************************/


extern "C"  {

  // Reaction rates for CH4-15step mechanism based on GRI2.11
  void ckwyp15step211_( double& P, double& T, double* Y, double* Wdot);
			    
  // Complexified version of the subroutine above
  void cplx15step211_( cplx& P, cplx& T, cplx* Y, cplx* Wdot);
			    
  // Reaction rates for CH4-15step mechanism based on GRI3
  void ckwyp15step30_( double& P, double& T, double* Y, double* Wdot);

  // Complexified version of the subroutine above
  void cplx15step30_( cplx& P, cplx& T, cplx* Y, cplx* Wdot);

}
			    
