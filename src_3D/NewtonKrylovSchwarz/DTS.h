#ifndef _DTS_INCLUDED
#define _DTS_INCLUDED


/*! *****************************************************************************************
 * class: DTS_Quad_Block
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
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
class DTS_Hexa_Block {
 public:
  int            NCi,        //!< Number of i-direction cells (including ghost cells). (WHY AM I STORING GHOST CELLS???) 
                 NCj,        //!< Number of j-direction cells (including ghost cells).
                 NCk;        //!< Number of k-direction cells (including ghost cells).
  SOLN_cSTATE ***Un,         //!< Conserved solution state at time-step n.
              ***Unminus1;   //!< Conserved solution state at time-step n-1.
  double         DTS_dTime;  //!< Physical Time step.
  
 
  //! Default constructor.
  DTS_Hexa_Block(void): NCi(0), NCj(0), NCk(0), DTS_dTime(ZERO),
			Un(NULL), Unminus1(NULL)  {}

  //! Constructor.
  DTS_Hexa_Block(const int i, const int j, const int k){
    allocates(i,j,k);
  }

  //! Destructor.
  ~DTS_Hexa_Block(void) {
    deallocate();
  }

  //! Memory allocation.
  void allocate(const int i, const int j, const int k){
    NCi = i; NCj = j; NCk = k;      

    Un = new SOLN_cSTATE**[NCi];
    Unminus1 = new SOLN_cSTATE**[NCi];
    for( int i=0; i<NCi; i++){
      Un[i] = new SOLN_cSTATE*[NCj];
      Unminus1[i] = new SOLN_cSTATE*[NCj];
      for (int j = 0; j < NCj ; j++ ){
	Un[i][j] = new SOLN_cSTATE[NCk];
	Unminus1[i][j] = new SOLN_cSTATE[NCk];
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
 
  void Store_Previous();
};

/*! *******************************************************************
 * DTS_Hexa_Block::Store_Previous_Solution                        *
 *                                                                    *
 * Store the previous solution.                                       *
 *                                                                    *
 **********************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
void DTS_Hexa_Block<SOLN_pSTATE,SOLN_cSTATE>::
Store_Previous(){

//   for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
//     for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
//       for (int k = 0; k < blocksize; k++) {
// 	Unminus1[i][j][k] =Un[i][j][k];        // t(n-1) 
// 	Un[i][j][k] = SolnBlk.U[i][j][k+1];    // t(n)	          
//       }
//     } 	      
//   }
 
}


#endif //_DTS_INCLUDED
