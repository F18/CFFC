#ifndef _NKS_DTS_INCLUDED
#define _NKS_DTS_INCLUDED


/*! *****************************************************************************************
 * class: DTS_NKS_Quad_Block
 *
 * @brief Dual-time-stepping NKS solution block data.
 *
 * This class creates and destroys the data arrays required for the
 * dual-time-stepping calculations for nonlinear partial
 * differential equations.
 *
 * \begin{verbatim}
 *      Un   -- Conserved solution state at time-step n.
 *      Uo   -- Conserved solution state at time-step n-1.
 * \end{verbatim}
 ********************************************************************************************/
template <typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
class DTS_NKS_Quad_Block {
 public:
  int            NCi,        //!< Number of i-direction cells (including ghost cells).
                 NCj,        //!< Number of j-direction cells (including ghost cells).
                 blocksize;  //!< Number of variables (equations being solved)
  double      ***Un,         //!< Conserved solution state at time-step n.
              ***Unminus1;   //!< Conserved solution state at time-step n-1.
  double         DTS_dTime;      //!< Physical Time step.
  INPUT_TYPE    *Input_Parameters; //!< Input Parameters Pointer

 
  //! Default constructor.
  DTS_NKS_Quad_Block(void): NCi(0), NCj(0), blocksize(0), DTS_dTime(ZERO),
			    Un(NULL), Unminus1(NULL), Input_Parameters(NULL) {}

  //! Constructor.
  DTS_NKS_Quad_Block(const int &Ni, const int &Nj, const int &blocksize, INPUT_TYPE &IPs) {
    allocate(Ni,Nj,blocksize,IPs);
  }

  //! Destructor.
  ~DTS_NKS_Quad_Block(void) {
    deallocate();
  }

  //! Memory allocation.
  void allocate(const int &Ni, const int &Nj, const int &_blocksize, INPUT_TYPE &IPs) {
    Input_Parameters = &(IPs);
    NCi = Ni; NCj = Nj; blocksize=_blocksize;
    Un = new double**[NCi];
    Unminus1 = new double**[NCi];
    for( int i=0; i<NCi; i++){
      Un[i] = new double*[NCj];
      Unminus1[i] = new double*[NCj];
      for (int j = 0; j < NCj ; j++ ){
	Un[i][j] = new double[blocksize];
	Unminus1[i][j] = new double[blocksize];
      }
    }
  }

  //! Memory deallocation.
  void deallocate(void) { 
    for (int i=0; i < NCi ; i++ ) {
      for (int j = 0; j < NCj ; j++ ){
	delete[] Un[i][j];       Un[i][j]=NULL; 
	delete[] Unminus1[i][j]; Unminus1[i][j]=NULL; 
      }
      delete[] Un[i];            Un[i]=NULL;  
      delete[] Unminus1[i];      Unminus1[i]=NULL;  
    }
    if (Un != NULL) { delete []Un; Un = NULL; }
    if (Unminus1 != NULL) { delete []Unminus1; Unminus1 = NULL; }
  }
 
  void Store_Previous(SOLN_BLOCK_TYPE &, const double &);
};

/*! *******************************************************************
 * DTS_NKS_Quad_Block::Store_Previous_Solution                        *
 *                                                                    *
 * Store the previous solution.                                       *
 *                                                                    *
 **********************************************************************/
template <typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void DTS_NKS_Quad_Block<SOLN_BLOCK_TYPE,INPUT_TYPE>::
Store_Previous(SOLN_BLOCK_TYPE &SolnBlk, const double &_dTime){

  for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      // Implicit Euler only uses Un so could save some time here ???
      if (Input_Parameters->NKS_IP.Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {
	for (int k = 0; k < blocksize; k++) {
	  Un[i][j][k] = SolnBlk.U[i][j][k+1];
	}
	// BDF2
      } else if (Input_Parameters->NKS_IP.Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {
	for (int k = 0; k < blocksize; k++) {
	  Unminus1[i][j][k] =Un[i][j][k];        // t(n-1) 
	  Un[i][j][k] = SolnBlk.U[i][j][k+1];    // t(n)	    
	}
      } 	      
    }
  }
  // Stdore Physical time step
  DTS_dTime = _dTime;

}


/******************** TEMPLATED FUNCTION **************************************************************/

/*! *******************************************************************
 * LHS_Time                                                           *
 *                                                                    *
 * Left Hand Side time componenet                                     *
 *        [ I*LHS_Time + J ] = -R(Un)                                 *
 *                                                                    *
 **********************************************************************/
template <typename INPUT_TYPE>
inline double LHS_Time(INPUT_TYPE &Input_Parameters, double& d_tau, const double &DTS_dTime){
  
  // Dual Time Stepping 
  if (Input_Parameters.NKS_IP.Dual_Time_Stepping){
    // Implicit Euler
    if (Input_Parameters.NKS_IP.Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {
      return (ONE/d_tau + ONE/DTS_dTime);
      //BDF2
    } else if (Input_Parameters.NKS_IP.Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {
       return (ONE/d_tau + THREE/(TWO*DTS_dTime));
    }
  }
   
  //Standard
  return ONE/(d_tau);
  
}

/*! *************************************************************************
 *  Overloading dUdt_Residual_Evaluation function to include Dual Time 
 *  Stepping source term.
 *
 ****************************************************************************/
template <typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
int dUdt_Residual_Evaluation_NKS(SOLN_BLOCK_TYPE &SolnBlk,
				 DTS_NKS_Quad_Block<SOLN_BLOCK_TYPE,INPUT_TYPE> &DTS_SolnBlk,
				 INPUT_TYPE &Input_Parameters){

  int error_flag = dUdt_Residual_Evaluation(SolnBlk,Input_Parameters);

  // Add dual time stepping Source Term to Residual ie. dUdt[i][j][0]
  if (Input_Parameters.NKS_IP.Dual_Time_Stepping) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
	for (int k = 0; k <  DTS_SolnBlk.blocksize; k++) {
	  //Implicit Euler  R(U_n)* = R(U_n) + (U - Un)/dt
	  if (Input_Parameters.NKS_IP.Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER) {		  
	    SolnBlk.dUdt[i][j][0][k+1] -= (SolnBlk.U[i][j][k+1] - DTS_SolnBlk.Un[i][j][k])/(DTS_SolnBlk.DTS_dTime); 
	    // BDF2         R(U_n)* = R(U_n) + (3U - 4Un +  Un-1)/2dt  
	  } else if (Input_Parameters.NKS_IP.Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {     
	    SolnBlk.dUdt[i][j][0][k+1] -= (THREE*SolnBlk.U[i][j][k+1] - FOUR*DTS_SolnBlk.Un[i][j][k]
					   + DTS_SolnBlk.Unminus1[i][j][k])/(TWO*DTS_SolnBlk.DTS_dTime);
	  }
	} 	      
      }
    }
  }  
     
  return error_flag;
}




#endif //_NKS_DTS_INCLUDED
