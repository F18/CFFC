#ifndef _GMRES3D_INCLUDED
#define _GMRES3D_INCLUDED

#include <iostream>
#include <cstdio>
#include <cmath> 

#ifndef _BLAS1_H_     
#include "blas1.h"   //Wrappers for BLAS Fortran libraries, so need to link with -lblas
#endif // _BLAS1_H_

#include "Block_Preconditioner.h"
#include "NKS_DTS.h"

using namespace std;

/***********************************************************
 * Class: GMRES_Block                                      *
 *                                                         *
 * Member functions                                        *
 *       s      -- Return residual                         *
 *      cs      -- Return cos vector                       *
 *      sn      -- Return sin vector                       *
 *       W      -- Return Az vector                        *
 *       z      -- Return inversion of preconditioner      *
 *                  times v vector                         *
 *       b      -- Return RHS vector                       *
 *       x      -- Return solution vector, delta u         *
 *       V      -- Return Krylov search direction vector   *
 *       H      -- Return Hessenberg matrix                *
 *   restart    -- Return restart                          *
 * overlap      -- Return level of overlap                 *
 * blocksize    -- Return number of variables              *
 * scalar_dim   -- Return NCi * NCj * NCk * blocksize      *
 *                   iterations                            *
 *                                                         *
 * ** MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING **     *
 * vector_swtich -- Return vector switch                   *
 *                  Note: 0 for x vector (solution)        *
 *                        1 for z vector (z = Minv * v)    *
 * search_directions -- Return number of search            *
 *                      directions                         *
 *     NCi      -- Return number of cells in               *
 *                 the i-direction (zeta-direction)        *
 *     ICl      -- Return lower index for cells in         *
 *                 the i-direction (zeta-direction)        *
 *     ICu      -- Return upper index for cells in         *
 *                 the i-direction (zeta-direction)        *
 *     NCj      -- Return number of cells in               *
 *                 the j-direction (eta-direction)         *
 *     JCl      -- Return lower index for cells in         *
 *                 the j-direction (eta-direction)         *
 *     JCu      -- Return upper index for cells in         *
 *                 the j-direction (eta-direction)         *
 *    Grid      -- Return a dummy object.                  *
 *                 (Required by message passing routines)  *   
 *  Nghost      -- Return number of ghost (halo or         *
 *                 overlap) cells.                         *
 *  Hexa_Block_ptr -- Pointer to solution block.           *
 *                                                         *
 *                                                         *
 ***********************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE>
class GMRES_Block{
private:

  /*********************************************/
  double *normalize_valuesR;         //vector(s) of normalizing values corresponding to SOLN_cSTATE
  double *normalize_valuesU;
  void set_normalize_values(void);   //set normalization values specific to SOLN_cSTATE
  double normalizeR(const double value, const int variable)     //normalizes "R" ie dUdt ie. norm(U)*norm(dt)
  {return (value/normalize_valuesR[variable]);}  
  double normalizeU(const double value, const int variable)     
  {return (value/normalize_valuesU[variable]);} 
  double normalizeUtoR(const double value, const int variable)     
  {return (value*normalize_valuesU[variable]/normalize_valuesR[variable]); }
  double denormalizeR(const double value, const int variable)           
  {return (value*normalize_valuesR[variable]);}  
  double denormalizeU(const double value, const int variable)   // Used for adding Uo + denormalize(detaUdt)
  {return (value*normalize_valuesU[variable]);} 
  /**************************************************/

  //3D Solution Block i,j,k,[var] indexing conversion to 1D Vector 
  int index(const int i, const int j, const int k,const int var) {return ((k*NCi*NCj+j*NCi+i)*blocksize+var);}
  int index(const int i, const int j, const int k) {return ((k*NCi*NCj+j*NCi+i)*blocksize);}

public:
  /* Solutition INFORMATION THAT DOESN'T CHANGE , could be static ???*/
  int                     restart; // number of gmres iterations before restart
  int                     overlap; // level of overlap  
  int                   blocksize; // number of variables            //Soln_ptr->NumVar();
  int                  scalar_dim; // xpts * ypts * blocksize

  /* USED INTERNAL TO GMRES ROUTINE */
  int               vector_switch; // to select the specified vector for message passing
  int           search_directions;    // number of search directions

  /* GMRES SOLUTION VECTORS ie DATA */
  double *                      s; // residual vector -> I think
  double *                     cs; // cos vector
  double *                     sn; // sin vector
  double *                      V; // Krylov search direction vector
  double *                      W; // A*z -> M^(-1)*x
  double *                      H; // Hessenberg matrix
  double *                      b; // RHS vector R(U)
  double *                      x; // initial guess of delta u
    
  /* MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING. */
  int                 NCi,ICl,ICu; // i-direction cell counters.
  int                 NCj,JCl,JCu; // j-direction cell counters.
  int                 NCk,KCl,KCu; // k-direction cell counters.
  int                      Nghost; // Number of ghost cells.

//   Grid3D_Quad_Block          Grid; // dummy pointer.
  Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>            *Hexa_Block_ptr;                                                    
  Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>      *IP;

  /****************************************************************/
  /************** Creation and copy constructors. *****************/
  /****************************************************************/  
  GMRES_Block(void):  restart(0), overlap(0), blocksize(0),
		      scalar_dim(0), search_directions(0),
		      normalize_valuesR(NULL), normalize_valuesU(NULL),
		      s(NULL), cs(NULL), sn(NULL), b(NULL), x(NULL),
		      H(NULL), W(NULL), V(NULL), Hexa_Block_ptr(NULL), IP(NULL),
		      vector_switch(1), NCi(0), ICl(0), ICu(0), 
		      NCj(0), JCl(0), JCu(0), NCk(0), KCl(0), KCu(0), Nghost(0) {}
  
  // GMRES_Block(const GMRES_Block &G);                //FIX so that it actually copies, not just passes pointers!!!!!!!!!!!!
  // GMRES_Block operator = (const GMRES_Block &G);    //Setup proper assignment operator
  
  /* Allocate and deallocate memory for GMRES_Block variables. */
  void allocate( Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs, 
		 Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>  *_Hexa_Block_ptr, 
		 const int &_blocksize);  

  void deallocate();
  ~GMRES_Block() { deallocate(); }
  void Initialize(void);

  //Member functions to manipulate data between GMRES vectors and Solution Block datastructure
  void calculate_perturbed_residual(const double &epsilon); 
  void calculate_perturbed_residual_2nd(const double &epsilon); 
  void calculate_perturbed_residual_Restart(const double &epsilon);
  void calculate_perturbed_residual_2nd_Restart(const double &epsilon);
  void calculate_Matrix_Free(const double &epsilon);
  void calculate_Matrix_Free_Restart(const double &epsilon);


  //Norm Calculations and Vector Operations
  double L2_Norm(const double* v);
  double Dotproduct(const double *v1, const double *v2);

  //Return denormalized deltaU
  double deltaU(const int i, const int j, const int k, const int var){ 
    return ( denormalizeU( x[index(i,j,k) ], var) ); }  

//   //Return info for testing/degugging.
//   double deltaU_test(const int i, const int j, const int k){ return x[index(i,j,k)]; }   
//   double b_test(const int i, const int j, const int k){ return b[index(i,j,k)]; }   

  /**************************************************************
   * Message Passing Member Functions for Message Passing       *
   * based on Solution_Block Data Structures                    *
   **************************************************************/
  /* Number of solution state variables. */
  int NumVar(void);

  /* Load send message passing buffer. */
  int LoadSendBuffer(double *buffer,
                     int &buffer_count,
                     const int buffer_size,
                     const int i_min, 
                     const int i_max,
                     const int i_inc,
                     const int j_min, 
                     const int j_max,
                     const int j_inc);
  int LoadSendBuffer_F2C(double *buffer,
                         int &buffer_count,
                         const int buffer_size,
                         const int i_min, 
                         const int i_max,
                         const int i_inc,
                         const int j_min, 
                         const int j_max,
                         const int j_inc);
  int LoadSendBuffer_C2F(double *buffer,
                         int &buffer_count,
                         const int buffer_size,
                         const int i_min, 
                         const int i_max,
                         const int i_inc,
                         const int j_min, 
                         const int j_max,
                         const int j_inc);
  /* Unload receive message passing buffer. */
  int UnloadReceiveBuffer(double *buffer,
                          int &buffer_count,
                          const int buffer_size,
                          const int i_min, 
                          const int i_max,
                          const int i_inc,
                          const int j_min, 
                          const int j_max,
                          const int j_inc);
  int UnloadReceiveBuffer_F2C(double *buffer,
                              int &buffer_count,
                              const int buffer_size,
                              const int i_min, 
                              const int i_max,
                              const int i_inc,
                              const int j_min, 
                              const int j_max,
                              const int j_inc);
  int UnloadReceiveBuffer_C2F(double *buffer,
                              int &buffer_count,
                              const int buffer_size,
                              const int i_min, 
                              const int i_max,
                              const int i_inc,
                              const int j_min, 
                              const int j_max,
                              const int j_inc);
  /* Subcell solution reconstruction within given computational cell. */
  void SubcellReconstruction(const int i,
                             const int j,
                             const int Limiter);

  // Load and unload conservative flux message passing buffer. 
  // NOT USED just added for compatibility with Message Passing Templates.
  int LoadSendBuffer_Flux_F2C(double *buffer,
                              int &buffer_count,
                              const int buffer_size,
                              const int i_min, 
                              const int i_max,
                              const int i_inc,
                              const int j_min, 
                              const int j_max,
                              const int j_inc){ return(0); }

  int UnloadReceiveBuffer_Flux_F2C(double *buffer,
				   int &buffer_count,
				   const int buffer_size,
				   const int i_min, 
				   const int i_max,
				   const int i_inc,
				   const int j_min, 
				   const int j_max,
				   const int j_inc) { return(0);} 
};


/*!*********************************************************************************
 *  Normalizing values and BCS needs to be Provided By _Quad_Block Specialization  *
 ***********************************************************************************/ 
template <typename SOLN_pSTATE,typename SOLN_cSTATE>
void GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
set_normalize_values(void) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF set_normalize_values for GMRES3D.h requried \n";
  exit(1);
}


/**************************************************************************
 * GMRES_Block::allocate -- Allocate memory.                              *
 **************************************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
void GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
allocate(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs, 
	 Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> *_Hexa_Block_ptr, 
	 const int &_blocksize)
{  
  //set pointer to solution block and set all necessary GMRES information
  Hexa_Block_ptr = _Hexa_Block_ptr;
  IP = &(IPs);

  NCi = Hexa_Block_ptr->NCi;
  ICl = Hexa_Block_ptr->ICl;
  ICu = Hexa_Block_ptr->ICu;
  NCj = Hexa_Block_ptr->NCj;
  JCl = Hexa_Block_ptr->JCl;
  JCu = Hexa_Block_ptr->JCu;
  NCk = Hexa_Block_ptr->NCk;
  KCl = Hexa_Block_ptr->KCl;
  KCu = Hexa_Block_ptr->KCu;
  Nghost = Hexa_Block_ptr->Nghost;

  overlap = IPs.NKS_IP.GMRES_Overlap;
  restart = IPs.NKS_IP.GMRES_Restart;

  //Check if a valid solution block and GMRES parameters
  assert(restart > 1);  assert(NCi > 1);  assert(NCj > 1); assert(NCk >1);
  blocksize = _blocksize; 
  scalar_dim = NCi * NCj * NCk * blocksize;

  // Allocate Memory
  normalize_valuesU = new double[blocksize];
  normalize_valuesR = new double[blocksize];
   s = new double[restart+1];
  cs = new double[restart];
  sn = new double[restart];
   b = new double[scalar_dim];
   x = new double[scalar_dim]; 
   H = new double[restart*(restart+1)];
   W = new double[restart * scalar_dim];
   V = new double[(restart + 1) * scalar_dim];

   //Setup Normalizing Values Based on Solution Block Type
   if(IPs.NKS_IP.Normalization){
     set_normalize_values();
   } else { //Set to 1
     for (int i=0; i< blocksize; i++){
       normalize_valuesU[i] = ONE;
       normalize_valuesR[i] = ONE;
     }
   }
}

/**************************************************************************
 * GMRES_Block::deallocate -- Deallocate memory.                          *
 **************************************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline void GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
deallocate() 
{
  if(normalize_valuesR !=NULL) delete[] normalize_valuesR; normalize_valuesR = NULL;
  if(normalize_valuesU !=NULL) delete[] normalize_valuesU; normalize_valuesU = NULL;
  if(s != NULL)      delete [] s;      s = NULL;
  if(cs != NULL)     delete [] cs;    cs = NULL;
  if(sn != NULL)     delete [] sn;    sn = NULL; 
  if(b != NULL)      delete [] b;      b = NULL; 
  if(x != NULL)      delete [] x;      x = NULL; 
  if(H != NULL)      delete [] H;      H = NULL;     
  if(W != NULL)      delete [] W;      W = NULL; 
  if(V != NULL)      delete [] V;      V = NULL; 
}

/**************************************************************************
 * GMRES_Block::Initialize()  Setup/Reset GMRES variables                 *
 **************************************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline void GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
Initialize(void)
{
  /* Initialize all GMRES variables except b, which is copied from dUdt */  
  for (int i=0;i<restart;++i) {  cs[i] = ZERO; sn[i] = ZERO; s[i] = ZERO;  }
  for (int i=0;i<scalar_dim;++i){  x[i] = ZERO; b[i] = ZERO; }
  for (int i=0;i<(restart*(restart+1));++i)         H[i] = ZERO;
  for (int i=0;i< restart*scalar_dim;++i)           W[i] = ZERO;
  for (int i=0;i<((restart + 1) * scalar_dim);++i)  V[i] = ZERO;
  
  // Load b ie RHS ie dUodt and normalize 
  for (int k = KCl;  k <= KCu; k++){
    for (int j = JCl;  j <= JCu; j++){
      for (int i = ICl;  i <= ICu; i++){
	for(int varindex =0; varindex < blocksize; varindex++){		
	  b[index(i,j,k,varindex)] = normalizeR(Hexa_Block_ptr->dUdt[i][j][k][0][varindex+1],varindex);
	}    	      
      } 
    }     
  } 
}

// 3D TO HERE !!

/**************************************************************************
 * Routine: calculate_pertubed_residual                                   *
 **************************************************************************/
// Calculate Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * W(i) )
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline void GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
calculate_perturbed_residual(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	Hexa_Block_ptr->U[i][j][varindex+1] = Hexa_Block_ptr->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*W[search_directions*scalar_dim+index(i,j,varindex)], varindex); 
      }   
      /* Update primitive variables. */
      Hexa_Block_ptr->W[i][j] = Hexa_Block_ptr->U[i][j].W();      
    }
  }  
}

// Calculate Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * x(i) )
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline void GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
calculate_perturbed_residual_Restart(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	Hexa_Block_ptr->U[i][j][varindex+1] = Hexa_Block_ptr->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  
      /* Update primitive variables. */
      Hexa_Block_ptr->W[i][j] = Hexa_Block_ptr->U[i][j].W();      
    }
  }  
}

/********************************************************
 * Routine: calculate_Matrix_Free                       *
 ********************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline void GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
calculate_Matrix_Free(const double &epsilon)
{
  //Taking into acount NKS overlap
  int JCl_overlap = 0; int JCu_overlap = 0;
  int ICu_overlap = 0; int ICl_overlap = 0;		  
  if(overlap){	
    if ( Hexa_Block_ptr->Grid.BCtypeS[ICl] == BC_NONE)  JCl_overlap = overlap; 
    if ( Hexa_Block_ptr->Grid.BCtypeN[ICu] == BC_NONE)  JCu_overlap = overlap;
    if ( Hexa_Block_ptr->Grid.BCtypeE[JCu] == BC_NONE)  ICu_overlap = overlap;
    if ( Hexa_Block_ptr->Grid.BCtypeW[JCl] == BC_NONE)  ICl_overlap = overlap;
  }
	   
  // Non-Overlap Ghost Cells R(U) already set to zero by dUdt calculation  
  /* V(i+1) = ( R(U+epsilon*W) - b) / epsilon - z / h */
  for (int j = JCl - JCl_overlap; j <= JCu + JCu_overlap; j++) {
    for (int i = ICl - ICl_overlap; i <= ICu + ICu_overlap; i++) {
      for(int k =0; k < blocksize; k++){	
	int iter = index(i,j,k);		
	//Matrix Free V(i+1) 
	V[(search_directions+1)*scalar_dim+iter] = 
	  ( normalizeR(Hexa_Block_ptr->dUdt[i][j][0][k+1],k) - b[iter]) / epsilon 
	  -  normalizeUtoR( W[(search_directions)*scalar_dim + iter] / ( Hexa_Block_ptr->dt[i][j]) , k );

// #ifdef _NKS_VERBOSE_NAN_CHECK
// 	// nan check most commonly caused by nans in dUdt !!!!
// 	if (V[(search_directions+1)*scalar_dim+iter] != V[(search_directions+1)*scalar_dim+iter] ){
// 	  cout<<"\n nan in V[ "<<(search_directions+1)*scalar_dim+iter<<"] at "<<i<<" "<<j<<" "<<k
// 	      <<" dUdt "<<  normalizeR(Hexa_Block_ptr->dUdt[i][j][0][k+1],k) <<" b "<< b[iter]
// 	      <<" z "<<W[(search_directions)*scalar_dim + iter]<< " h "<<( Hexa_Block_ptr->dt[i][j]*ao);
// 	}
// #endif 
      }      
    } 
  } 
}

template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline void GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
calculate_Matrix_Free_Restart(const double &epsilon)
{
  //Taking into acount NKS overlap
  int JCl_overlap = 0; int JCu_overlap = 0;
  int ICu_overlap = 0; int ICl_overlap = 0;		  
  if(overlap){	
    if ( Hexa_Block_ptr->Grid.BCtypeS[ICl] == BC_NONE)  JCl_overlap = overlap; 
    if ( Hexa_Block_ptr->Grid.BCtypeN[ICu] == BC_NONE)  JCu_overlap = overlap;
    if ( Hexa_Block_ptr->Grid.BCtypeE[JCu] == BC_NONE)  ICu_overlap = overlap;
    if ( Hexa_Block_ptr->Grid.BCtypeW[JCl] == BC_NONE)  ICl_overlap = overlap;
  }
	  
  // Non-Overlap Ghost Cells R(U) already set to zero by dUdt 
  
  /* V(0) = ( R(U + epsilon*W) - b) / epsilon - x / h */
  for (int j = JCl - JCl_overlap; j <= JCu + JCu_overlap; j++) {
    for (int i = ICl - ICl_overlap; i <= ICu + ICu_overlap; i++) {
      for(int k =0; k < blocksize; k++){	
	int iter = index(i,j,k);		
	//Matrix Free V(0) 
	V[iter] =  ( normalizeR(Hexa_Block_ptr->dUdt[i][j][0][k+1],k) - b[iter]) / epsilon 
	  - normalizeUtoR( x[iter] / ( Hexa_Block_ptr->dt[i][j]) , k);
      }      
    } 
  } 
}

/********************************************************
 * Routine: L2_norms                                    *
 ********************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline double GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
L2_Norm(const double* v) {
  integer inc = 1;
  return ( F77NAME(dnrm2)( &scalar_dim, v, &inc) );
}

/********************************************************
 * Routine: Dot_Product                                 *
 ********************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline double GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
Dotproduct(const double *v1, const double *v2) {
  integer inc = 1;
  return (F77NAME(ddot)( &scalar_dim, v1, &inc, v2, &inc));	   
}


/*******************************************************************************
 *                                                                             *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                              *
 *                                                                             *
 *******************************************************************************/

/*******************************************************************************
 * GMRES_Block::NumVar -- Returns number of state variables.                   *
 *******************************************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline int GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
NumVar(void) {
  return (blocksize);             
}

/*******************************************************************************
 * GMRES_Block::LoadSendBuffer -- Loads send message buffer.                   *
 *******************************************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline int GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
LoadSendBuffer(double *buffer,
	       int &buffer_count,
	       const int buffer_size,
	       const int i_min, 
	       const int i_max,
	       const int i_inc,
	       const int j_min, 
	       const int j_max,
	       const int j_inc) {
  
  for (int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    for (int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for (int  k = 0 ; k < blocksize; ++ k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return(1);
	if (vector_switch) {
	  buffer[buffer_count] = W[(search_directions)*scalar_dim + index(i,j,k)];
	} else {
	  buffer[buffer_count] = x[index(i,j,k)];
	}
      } 
    } 
  } 
  return(0);
} 

/*******************************************************************************
 * GMRES_Block::LoadSendBuffer_F2C -- Loads send message buffer for            *
 *                                    fine to coarse block message             *
 *                                    passing.                                 *
 *******************************************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline int GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
LoadSendBuffer_F2C(double *buffer,
		   int &buffer_count,
		   const int buffer_size,
		   const int i_min, 
		   const int i_max,
		   const int i_inc,
		   const int j_min, 
		   const int j_max,
		   const int j_inc) {
  
  for (int j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
    for (int i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
      for (int k = 0 ; k < blocksize; ++ k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return(1);
	if (vector_switch) {
	  buffer[buffer_count] = (Hexa_Block_ptr->Grid.Cell[i  ][j  ].A*W[(search_directions)*scalar_dim + index(i  ,j  ,k)]+
				  Hexa_Block_ptr->Grid.Cell[i+1][j  ].A*W[(search_directions)*scalar_dim + index(i+1,j  ,k)]+
				  Hexa_Block_ptr->Grid.Cell[i  ][j+1].A*W[(search_directions)*scalar_dim + index(i  ,j+1,k)]+
				  Hexa_Block_ptr->Grid.Cell[i+1][j+1].A*W[(search_directions)*scalar_dim + index(i+1,j+1,k)])/
	                         (Hexa_Block_ptr->Grid.Cell[i  ][j  ].A+
				  Hexa_Block_ptr->Grid.Cell[i+1][j  ].A+
				  Hexa_Block_ptr->Grid.Cell[i  ][j+1].A+
				  Hexa_Block_ptr->Grid.Cell[i+1][j+1].A);
	} else {
	  buffer[buffer_count] = (Hexa_Block_ptr->Grid.Cell[i  ][j  ].A*x[index(i  ,j  ,k)]+
				  Hexa_Block_ptr->Grid.Cell[i+1][j  ].A*x[index(i+1,j  ,k)]+
				  Hexa_Block_ptr->Grid.Cell[i  ][j+1].A*x[index(i  ,j+1,k)]+
				  Hexa_Block_ptr->Grid.Cell[i+1][j+1].A*x[index(i+1,j+1,k)])/
	                         (Hexa_Block_ptr->Grid.Cell[i  ][j  ].A+
				  Hexa_Block_ptr->Grid.Cell[i+1][j  ].A+
				  Hexa_Block_ptr->Grid.Cell[i  ][j+1].A+
				  Hexa_Block_ptr->Grid.Cell[i+1][j+1].A);
	} 
      } 
    } 
  } 
  return(0);
}

/*******************************************************************************
 * GMRES_Block::LoadSendBuffer_C2F -- Loads send message buffer for            *
 *                                    coarse to fine block message             *
 *                                    passing.                                 *
 *******************************************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline int GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
LoadSendBuffer_C2F(double *buffer,
		   int &buffer_count,
		   const int buffer_size,
		   const int i_min, 
		   const int i_max,
		   const int i_inc,
		   const int j_min, 
		   const int j_max,
		   const int j_inc) {
  int i, j, k;
  Vector2D dX;
  SOLN_pSTATE Wcoarse, Wfine;
  int LIMITER = LIMITER_ZERO; //LIMITER_VENKATAKRISHNAN 

  if (j_min == j_max) { // North or south boundary.
     // Four different orderings to consider depending on the value of i_inc & j_inc.
    if (j_inc > 0) {             
      if (i_inc > 0) {

	/******************************* CASE #1 ***************************************/
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j_min, LIMITER);                    
	  
	  // Evaluate SW sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i,j_min,k)];
	    } 
	  } 
	  dX = (Hexa_Block_ptr->Grid.Node[i][j_min].X+
		HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min].X)+
		HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i][j_min+1].X)+
		Hexa_Block_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
	        Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
	    (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;

	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } 
	  
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min].X)+
		Hexa_Block_ptr->Grid.Node[i+1][j_min].X + Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
		HALF*(Hexa_Block_ptr->Grid.Node[i+1][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
	        Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
	    (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
            
	  for ( k = 0 ; k < blocksize; ++ k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  }
	} 
	
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate NW sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i,j_min,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i][j_min+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
                    Hexa_Block_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min+1].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
          
	      for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i+1][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min+1].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)+
                    Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
	
	/******************************* CASE #2 ***************************************/
      } else {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER);
              // Evaluate SE sub (fine) cell values.
              for (k = 0 ; k < blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	         } else {
	            Wcoarse[k+1] = x[index(i,j_min,k)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min].X)+
                    Hexa_Block_ptr->Grid.Node[i+1][j_min].X+
                    Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i+1][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Hexa_Block_ptr->Grid.Node[i][j_min].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i][j_min+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for ( k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate NE sub (fine) cell values.
              for (k = 0 ; k < blocksize; ++ k) {
   	         if (vector_switch) {
	           Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	         } else {
	            Wcoarse[k+1] = x[index(i,j_min,k)];
	         } /* endif */
	      } /* endfor */
              dX = (Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i+1][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min+1].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)+
                    Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i][j_min+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
                    Hexa_Block_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min+1].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
        } /* endif */

      /******************************* CASE #3 ***************************************/
    } else {
        if (i_inc > 0) {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER);
              // Evaluate NW sub (fine) cell values.
              for (k = 0 ; k < blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	         } else {
	            Wcoarse[k+1] = x[index(i,j_min,k)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i][j_min+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
                    Hexa_Block_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min+1].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i+1][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min+1].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)+
                    Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SW sub (fine) cell values.
             for (k = 0 ; k < blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	         } else {
	            Wcoarse[k+1] = x[index(i,j_min,k)];
	         } /* endif */
	      } /* endfor */
              dX = (Hexa_Block_ptr->Grid.Node[i][j_min].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i][j_min+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for ( k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min].X)+
                    Hexa_Block_ptr->Grid.Node[i+1][j_min].X+
                    Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i+1][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
	   /******************************* CASE #4 ***************************************/
        } else {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER);
              // Evaluate NE sub (fine) cell values.
              for (k = 0 ; k < blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	         } else {
	            Wcoarse[k+1] = x[index(i,j_min,k)];
	         } /* endif */
	      } /* endfor */
              dX = (Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i+1][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min+1].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)+
                    Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;

	      // Wfine = Wfine.U();  ??MISTAKE MADE IN KALVINS ORIGINAL ?????

             for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i][j_min+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
                    Hexa_Block_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min+1].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SE sub (fine) cell values.
              for (k = 0 ; k < blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	         } else {
	            Wcoarse[k+1] = x[index(i,j_min,k)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min].X)+
                    Hexa_Block_ptr->Grid.Node[i+1][j_min].X+
                    Hexa_Block_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i+1][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for ( k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Hexa_Block_ptr->Grid.Node[i][j_min].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i][j_min].X+Hexa_Block_ptr->Grid.Node[i][j_min+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdx[i][j_min])*dX.x +
                      (Hexa_Block_ptr->phi[i][j_min]^Hexa_Block_ptr->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */


    /******************************* CASE #5 ***************************************/
  } else { // East or west boundary.
     // Four different orderings to consider depending on the value of i_inc & j_inc.
     if (j_inc > 0) {
        if (i_inc > 0) {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER);
              // Evaluate SW sub (fine) cell values.
              for (k = 0 ; k < blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i_min,j,k)];
	         } else {
	            Wcoarse[k+1] = x[index(i_min,j,k)];
	         } /* endif */
	      } /* endfor */
              dX = (Hexa_Block_ptr->Grid.Node[i_min][j].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min][j+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j].X)+
                    Hexa_Block_ptr->Grid.Node[i_min+1][j].X+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min+1][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min][j+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    Hexa_Block_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j+1].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min+1][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j+1].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)+
                    Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
	   /******************************* CASE #6 ***************************************/
        } else {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER);
              // Evaluate SE sub (fine) cell values.
              for (k = 0 ; k < blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i_min,j,k)];
	         } else {
	            Wcoarse[k+1] = x[index(i_min,j,k)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j].X)+
                    Hexa_Block_ptr->Grid.Node[i_min+1][j].X+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min+1][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Hexa_Block_ptr->Grid.Node[i_min][j].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min][j+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min+1][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j+1].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)+
                    Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min][j+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    Hexa_Block_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j+1].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
        } /* endif */	
	/******************************* CASE #7 ***************************************/
     } else {
        if (i_inc > 0) {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER);
              // Evaluate NW sub (fine) cell values.
              for (k = 0 ; k < blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i_min,j,k)];
	         } else {
	            Wcoarse[k+1] = x[index(i_min,j,k)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min][j+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    Hexa_Block_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j+1].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min+1][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j+1].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)+
                    Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Hexa_Block_ptr->Grid.Node[i_min][j].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min][j+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j].X)+
                    Hexa_Block_ptr->Grid.Node[i_min+1][j].X+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min+1][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */\

	   /******************************* CASE #8 ***************************************/
        } else {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER);
              // Evaluate NE sub (fine) cell values.
              for (k = 0 ; k < blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i_min,j,k)];
	         } else {
	            Wcoarse[k+1] = x[index(i_min,j,k)];
	         } /* endif */
	      } /* endfor */
              dX = (Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min+1][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j+1].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)+
                    Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min][j+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    Hexa_Block_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j+1].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j].X)+
                    Hexa_Block_ptr->Grid.Node[i_min+1][j].X+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min+1][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Hexa_Block_ptr->Grid.Node[i_min][j].X+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(Hexa_Block_ptr->Grid.Node[i_min][j].X+Hexa_Block_ptr->Grid.Node[i_min][j+1].X)+
                    Hexa_Block_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   Hexa_Block_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdx[i_min][j])*dX.x +
                      (Hexa_Block_ptr->phi[i_min][j]^Hexa_Block_ptr->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */
  } /* endif */
  return(0);

}

/*******************************************************************************
 * GMRES_Block::UnloadReceiveBuffer -- Unloads receive message buffer.         *
 *******************************************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline int GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
UnloadReceiveBuffer(double *buffer,
		    int &buffer_count,
		    const int buffer_size,
		    const int i_min, 
		    const int i_max,
		    const int i_inc,
		    const int j_min, 
		    const int j_max,
		    const int j_inc) {
  int i, j, k;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for ( k = 0 ; k < blocksize; ++ k) {
	  buffer_count++;
           if (buffer_count >= buffer_size) return(1);
	   if (vector_switch) {
	     W[(search_directions)*scalar_dim + index(i,j,k)] = buffer[buffer_count];  
	   } else {
	     x[index(i,j,k)] = buffer[buffer_count];
	   }
        } 
     } 
  } 
  return(0);
}

/*******************************************************************************
 * GMRES_Block::UnloadReceiveBuffer_F2C -- Unloads receive message             *
 *                                         buffer for fine to coarse           *
 *                                         block message passing.              *
 *******************************************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline int GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
UnloadReceiveBuffer_F2C(double *buffer,
			int &buffer_count,
			const int buffer_size,
			const int i_min, 
			const int i_max,
			const int i_inc,
			const int j_min, 
			const int j_max,
			const int j_inc) {

  for (int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for (int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for (int k = 0 ; k < blocksize; ++ k) {
	  buffer_count++;
           if (buffer_count >= buffer_size) return(1);
	   if (vector_switch) {
	      W[(search_directions)*scalar_dim + index(i,j,k)] = buffer[buffer_count];
	   } else {
	      x[index(i,j,k)] = buffer[buffer_count];
	   } 
        } 
     } 
  } 

  return(0);
}

/*******************************************************************************
 * GMRES_Block::UnloadReceiveBuffer_C2F -- Unloads receive message             *
 *                                         buffer for coarse to fine           *
 *                                         block message passing.              *
 *******************************************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline int GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
UnloadReceiveBuffer_C2F(double *buffer,
			int &buffer_count,
			const int buffer_size,
			const int i_min, 
			const int i_max,
			const int i_inc,
			const int j_min, 
			const int j_max,
			const int j_inc) {

  for (int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for (int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for (int k = 0 ; k < blocksize; ++ k) {
	  buffer_count++;
	  if (buffer_count >= buffer_size) return(1);
	  if (vector_switch) {
	    W[(search_directions)*scalar_dim + index(i,j,k)] = buffer[buffer_count];
	  } else {
	    x[index(i,j,k)] = buffer[buffer_count];
	  }	  
        } 
     } 
  }

  return(0);
}

/**************************************************************************
 * GMRES_Block::SubcellReconstruction --                                  *
 *              Performs the subcell reconstruction of solution state     *
 *              within a given cell (i,j) of the computational mesh for   *
 *              the specified quadrilateral solution block.               *
 **************************************************************************/
template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
inline void GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
SubcellReconstruction(const int i, 
		      const int j,
		      const int Limiter) {
  
  int n, n2, n_pts, i_index[8], j_index[8], k;
  double u0, u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  SOLN_pSTATE U0, DU, DUDx_ave, DUDy_ave, W_VACUUM;
  W_VACUUM.Vacuum();
  
  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

//   if (i == ICl-Nghost || i == ICu+Nghost ||
//       j == JCl-Nghost || j == JCu+Nghost) {
//     n_pts = 0;
//   } else if ((i == ICl-Nghost+1) && 
// 	     (Grid.BCtypeW[j] != BC_NONE)) {
//     if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
//       n_pts = 0;
//     } else if (Hexa_Block_ptr->Grid.BCtypeW[j] == BC_PERIODIC ||
//                  Hexa_Block_ptr->Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
//                  Hexa_Block_ptr->Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//                  Hexa_Block_ptr->Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
//          if (j == JCl) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i+1; j_index[1] = j  ;
//             i_index[2] = i-1; j_index[2] = j+1;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (j == JCu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (j == JCl) {
//             n_pts = 3;
//             i_index[0] = i+1; j_index[0] = j  ;
//             i_index[1] = i  ; j_index[1] = j+1;
//             i_index[2] = i+1; j_index[2] = j+1;
//          } else if (j == JCu) {
//             n_pts = 3;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//          } else {
//             n_pts = 5;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } /* endif */
//       } /* endif */           
//     } else if ((i == ICu+Nghost-1) && 
//                (Hexa_Block_ptr->Grid.BCtypeE[j] != BC_NONE)) {
//       if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
//          n_pts = 0;
//       } else if (Hexa_Block_ptr->Grid.BCtypeE[j] == BC_PERIODIC ||
//                  Hexa_Block_ptr->Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//                  Hexa_Block_ptr->Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//                  Hexa_Block_ptr->Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
//          if (j == JCl) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i+1; j_index[1] = j  ;
//             i_index[2] = i-1; j_index[2] = j+1;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (j == JCu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (j == JCl) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i-1; j_index[1] = j+1;
//             i_index[2] = i  ; j_index[2] = j+1;
//          } else if (j == JCu) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//          } else {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//             i_index[3] = i-1; j_index[3] = j+1;
//             i_index[4] = i  ; j_index[4] = j+1;
//          } /* endif */
//       } /* endif */
//     } else if ((j == JCl-Nghost+1) && 
//                (Hexa_Block_ptr->Grid.BCtypeS[i] != BC_NONE)) {
//       if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
//          n_pts = 0;
//       } else if (Hexa_Block_ptr->Grid.BCtypeS[i] == BC_PERIODIC ||
//                  Hexa_Block_ptr->Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//                  Hexa_Block_ptr->Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//                  Hexa_Block_ptr->Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
//          if (i == ICl) {
//             n_pts = 5;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (i == ICu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//             i_index[3] = i-1; j_index[3] = j+1;
//             i_index[4] = i  ; j_index[4] = j+1;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (i == ICl) {
//             n_pts = 3;
//             i_index[0] = i+1; j_index[0] = j  ;
//             i_index[1] = i  ; j_index[1] = j+1;
//             i_index[2] = i+1; j_index[2] = j+1;
//          } else if (i == ICu) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i-1; j_index[1] = j+1;
//             i_index[2] = i  ; j_index[2] = j+1;
//          } else {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i+1; j_index[1] = j  ;
//             i_index[2] = i-1; j_index[2] = j+1;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } /* endif */
//       } /* endif */
//     } else if ((j == JCu+Nghost-1) && 
//                (Hexa_Block_ptr->Grid.BCtypeN[i] != BC_NONE)) {
//       if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
//          n_pts = 0;
//       } else if (Hexa_Block_ptr->Grid.BCtypeN[i] == BC_PERIODIC ||
//                  Hexa_Block_ptr->Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//                  Hexa_Block_ptr->Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//                  Hexa_Block_ptr->Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
//          if (i == ICl) {
//             n_pts = 5;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (i == ICu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//             i_index[3] = i-1; j_index[3] = j+1;
//             i_index[4] = i  ; j_index[4] = j+1;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (i == ICl) {
//             n_pts = 3;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//          } else if (i == ICu) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//          } else {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//          } /* endif */
//       } /* endif */

  //FOR VISCOUS -> CHANGED TO USE ALL 8
  if (i == ICl-Nghost || i == ICu+Nghost ||
      j == JCl-Nghost || j == JCu+Nghost) {
    n_pts = 0;
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  } /* endif */

  if (n_pts > 0) {
      DUDx_ave = W_VACUUM;
      DUDy_ave = W_VACUUM;
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;

      for ( k = 0 ; k < blocksize; ++ k) {
	if (vector_switch) {
	  U0[k+1] = W[(search_directions)*scalar_dim + index(i,j,k)];
	} else {
	  U0[k+1] = x[index(i,j,k)];
	}
      } 

      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	dX = Hexa_Block_ptr->Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - Hexa_Block_ptr->Grid.Cell[i][j].Xc;
	for ( k = 0 ; k < blocksize; ++ k) {
	  if (vector_switch) {
	    DU[k+1] = W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] , k)] -  U0[k+1];
	  } else {
	    DU[k+1] = x[index( i_index[n2] , j_index[n2] , k)] -  U0[k+1];
	  } /* endif */
	} /* endfor */
	
	DUDx_ave += DU*dX.x;
	DUDy_ave += DU*dX.y;
	DxDx_ave += dX.x*dX.x;
	DxDy_ave += dX.x*dX.y;
	DyDy_ave += dX.y*dX.y;
      } /* endfor */
  					    
      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);

      Hexa_Block_ptr->dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                            (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      Hexa_Block_ptr->dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                            (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
      // Calculate slope limiters. 
      if (!Hexa_Block_ptr->Freeze_Limiter) {
         for ( n = 1 ; n <= blocksize ; ++n ) {
	    u0 = U0[n];
            u0Min = U0[n];
            u0Max = u0Min;
            for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
               if (vector_switch) {
                  u0Min = min(u0Min, W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] ,n-1)]);
                  u0Max = max(u0Max, W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] ,n-1)]);
               } else {
                  u0Min = min(u0Min, x[index(i_index[n2] , j_index[n2] , n-1)]);
                  u0Max = max(u0Max, x[index(i_index[n2] , j_index[n2] , n-1)]);
               } /* endif */
            } /* endfor */
    
            dX = Hexa_Block_ptr->Grid.xfaceE(i, j)-Hexa_Block_ptr->Grid.Cell[i][j].Xc;
            uQuad[0] = u0 + 
                       Hexa_Block_ptr->dWdx[i][j][n]*dX.x +
                       Hexa_Block_ptr->dWdy[i][j][n]*dX.y ;
            dX = Hexa_Block_ptr->Grid.xfaceW(i, j)-Hexa_Block_ptr->Grid.Cell[i][j].Xc;
            uQuad[1] = u0 + 
                       Hexa_Block_ptr->dWdx[i][j][n]*dX.x +
                       Hexa_Block_ptr->dWdy[i][j][n]*dX.y ;
            dX = Hexa_Block_ptr->Grid.xfaceN(i, j)-Hexa_Block_ptr->Grid.Cell[i][j].Xc;
            uQuad[2] = u0 + 
                       Hexa_Block_ptr->dWdx[i][j][n]*dX.x +
                       Hexa_Block_ptr->dWdy[i][j][n]*dX.y ;
            dX = Hexa_Block_ptr->Grid.xfaceS(i, j)-Hexa_Block_ptr->Grid.Cell[i][j].Xc;
            uQuad[3] = u0 + 
                       Hexa_Block_ptr->dWdx[i][j][n]*dX.x +
                       Hexa_Block_ptr->dWdy[i][j][n]*dX.y ;
    
            switch(Limiter) {
              case LIMITER_ONE :
                phi_n = ONE;
                break;
              case LIMITER_ZERO :
                phi_n = ZERO;
                break;
              case LIMITER_BARTH_JESPERSEN :
                phi_n = Limiter_BarthJespersen(uQuad, u0, 
                                               u0Min, u0Max, 4);
                break;
              case LIMITER_VENKATAKRISHNAN :
                phi_n = Limiter_Venkatakrishnan(uQuad, u0, 
                                                u0Min, u0Max, 4);
                break;
              case LIMITER_VANLEER :
                phi_n = Limiter_VanLeer(uQuad, u0, 
                                        u0Min, u0Max, 4);
                break;
              case LIMITER_VANALBADA :
                phi_n = Limiter_VanAlbada(uQuad, u0, 
                                          u0Min, u0Max, 4);
                break;
              default:
                phi_n = Limiter_BarthJespersen(uQuad, u0, 
                                               u0Min, u0Max, 4);
                break;
            } /* endswitch */

	    Hexa_Block_ptr->phi[i][j][n] = phi_n;

         } /* endfor */
      } /* endif */
  } else {
      Hexa_Block_ptr->dWdx[i][j] = W_VACUUM;
      Hexa_Block_ptr->dWdy[i][j] = W_VACUUM; 
      Hexa_Block_ptr->phi[i][j]  = W_VACUUM;
  } /* endif */

}

// /*******************************************************************************
//  * GMRES_Block::LoadSendBuffer_Flux_F2C -- Loads send message buffer for       *
//  *                                         fine to coarse block message        *
//  *                                         passing of conservative             *
//  *                                         solution fluxes.                    *
//  *******************************************************************************/
// template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
// inline int GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
// LoadSendBuffer_Flux_F2C(double *buffer,
// 			int &buffer_count,
// 			const int buffer_size,
// 			const int i_min, 
// 			const int i_max,
// 			const int i_inc,
// 			const int j_min, 
// 			const int j_max,
// 			const int j_inc) {
//   return(0); 
// }

// /*******************************************************************************
//  * GMRES_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message        *
//  *                                              buffer for fine to coarse      *
//  *                                              block message passing of       *
//  *                                              conservative solution fluxes.  *
//  *******************************************************************************/
// template <typename SOLN_pSTATE,typename SOLN_cSTATE> 
// inline int GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>::
// UnloadReceiveBuffer_Flux_F2C(double *buffer,
// 			     int &buffer_count,
// 			     const int buffer_size,
// 			     const int i_min, 
// 			     const int i_max,
// 			     const int i_inc,
// 			     const int j_min, 
// 			     const int j_max,
// 			     const int j_inc) {
//   return(0);
// }





///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////


/*!*********************************************************
 * Class: GMRES_RightPrecon_MatrixFree                     *
 *                                                         *
 * Performs GMRES using each of the blocks ....            *
 *                                                         *
 *                                                         *
 ***********************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
class GMRES_RightPrecon_MatrixFree {

private:
  
  //Pointers to Solution Info and DATA
  Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >  *Hexa_Blocks;                                                    
  Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>              *IP;

  //GMRES LOCAL DYNAMIC DATA
  int Number_of_GMRES_Iterations;
  double relative_residual;
  double global_time_step_size;
  GMRES_Block<SOLN_pSTATE,SOLN_cSTATE> *G;                  

  //private GMRES helper functions
  void GeneratePlaneRotation(const double &dx,const double &dy, double &cs, double &sn);
  void ApplyPlaneRotation(double &dx, double &dy,const double &cs,const double &sn);

public:
  
  //default constructor's 
  GMRES_RightPrecon_MatrixFree(void):
    Hexa_Blocks(NULL), IP(NULL), G(NULL), 
    Number_of_GMRES_Iterations(0), global_time_step_size(ZERO), relative_residual(ZERO) {}
  
  GMRES_RightPrecon_MatrixFree(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs, 
			       Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Hexa_Multi_Block_List, 
			       const int &blocksize);

  // Proper functions not built yet.
  // GMRES_RightPrecon_MatrixFree(const &GMRES_RightPrecon_MatrixFree) {}  
  // GMRES_RightPrecon_MatrixFree operator= (const &GMRES_RightPrecon_MatrixFree) {}

  //Memory Handlers
  void allocate() { G = new GMRES_Block<SOLN_pSTATE,SOLN_cSTATE>[Hexa_Blocks->Size_of_Block_List]; } 
  void deallocate() {  if(G != NULL) delete[] G;  G = NULL; }
  ~GMRES_RightPrecon_MatrixFree() { deallocate(); }

  // Constructors
  int solve(Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE> *Block_precon, 
	    const int Number_of_Newton_Steps);
  
  // CHANGE FOR 3D !!
  double deltaU(const int Bcount, const int i, const int j, const int k,const int var){
    return (G[Bcount].deltaU(i,j,k,var)); }
  double deltaU_test(const int Bcount, const int i, const int j, const int k){
    return (G[Bcount].deltaU_test(i,j,k)); }
  double b_test(const int Bcount, const int i, const int j, const int k){
    return (G[Bcount].b_test(i,j,k)); }  

};

/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline GMRES_RightPrecon_MatrixFree<SOLN_pSTATE,SOLN_cSTATE>::
GMRES_RightPrecon_MatrixFree(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs, 
			     Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Hexa_Multi_Block_List, 			     
			     const int &blocksize):
  Hexa_Blocks(&Hexa_Multi_Block_List),
  IP(&IPs),
  Number_of_GMRES_Iterations(0),
  relative_residual(0)  
{
  //Setup Memory for GMRES_Block's 
  allocate(); 

  // Setup GMRES for each Block that is used on this processor
  for (int  Bcount = 0 ; Bcount < Hexa_Multi_Block_List.Size_of_Block_List; ++Bcount ) {
    if ( Hexa_Multi_Block_List.Block_Used[Bcount] ) {
      G[Bcount].allocate(IPs,Hexa_Multi_Block_List.Hexa_Block_List[Bcount],blocksize);
    } 
  } 
}


/********************************************************
 * Routine: GMRES Plane Rotation Functions              *
 ********************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void GMRES_RightPrecon_MatrixFree<SOLN_pSTATE,SOLN_cSTATE>::
GeneratePlaneRotation(const double &dx, const double &dy, double &cs, double &sn) {

  if (dy == ZERO) {  //EQUALITIES WITH FP!!! -> NOT RELIABLE
    cs = ONE;
    sn = ZERO;
  } else if (fabs(dy) > fabs(dx)) {
    double temp = dx / dy;
    sn = ONE / sqrt( ONE + temp*temp );
    cs = temp * sn;
  } else {
    double temp = dy / dx;
    cs = ONE / sqrt( ONE + temp*temp );
    sn = temp * cs;
  }
}
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void GMRES_RightPrecon_MatrixFree<SOLN_pSTATE,SOLN_cSTATE>::
ApplyPlaneRotation(double &dx, double &dy, const double &cs,const double &sn) {
  double temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx =  temp; //cs * dx + sn * dy;
}


 /****************************************************************************
  * GMRES_RightPrecon_MatrixFree::solve -- apply right-preconditioned        *
  *                                              matrix-free GMRES.          *
  ****************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
int GMRES_RightPrecon_MatrixFree<SOLN_pSTATE,SOLN_cSTATE>::
solve(Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE> *Block_precon, 
      const int Number_of_Newton_Steps) {   
  
  int error_flag;

//   int m1 = Input_Parameters->NKS_IP.GMRES_Restart+1; //restart+1 used in H calculation
//   int error_flag(0);  
//   double resid0(ZERO);
//   double beta(ZERO);
//   double EPSILON_PARAM(1.0e-6);   //  1e0-8 original, 1e-10 fixed minimum, Viscous stable @ 1e-6 
//   double epsilon(EPSILON_PARAM);

//   //FORTRAN NAMES
//   integer inc(1);  // vector stride is always 1
//   doublereal temp;
//   Number_of_GMRES_Iterations = 0;

//   /**************************************************************************/
//   //Setup/Reset GMRES variables
//   for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {            
//     if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
//       G[Bcount].Initialize();
//     }
//   }
//   /**************************************************************************/
  
//   /**************************************************************************/
//   /******************* BEGIN GMRES CALCULATION ******************************/
//   /**************************************************************************/
//   bool i_first_time_through = true;
//   bool i_restart_flag = false;

//   do {

//     /**************************************************************************/
//     /************* FIRST TIME THROUGH - NO RESTART NEEDED *********************/
//     /**************************************************************************/    
//     // V(0) already set to ZERO when initialized, so nothing to do if not restart. 
      
//     /**************************************************************************/
//     /********** NOT FIRST TIME THROUGH - RESTART APPLIED **********************/
//     /**************************************************************************/
//     if(!i_first_time_through) {                                                     

//       if (CFDkit_Primary_MPI_Processor() && i_restart_flag){   
// 	cout << "\n GMRES Restarted at -> GMRES (Inner Iterations) = " << Number_of_GMRES_Iterations; 	 
//       } 
      
//       /************************************************************************/
//       /* Set vector switch for message passing or x (solution) instead of W. */
//       for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){  
// 	  G[Bcount].vector_switch = 0;
// 	}
//       }
//       /************************************************************************/

//       /************************************************************************/
//       /* MPI barrier to ensure processor synchronization. */
//       CFDkit_Barrier_MPI();  
      
//       /* Send "x" solution information between neighbouring blocks.*/          
//       error_flag = Send_All_Messages(G, 
// 				     *List_of_Local_Solution_Blocks,
// 				     G[0].NumVar(), 
// 				     OFF);
//       if (error_flag) {
// 	cout << "\n GMRES2D ERROR: GMRES message passing error on processor "
// 	     << List_of_Local_Solution_Blocks->ThisCPU
// 	     << ".\n";
// 	cout.flush();
//       } 
//       error_flag = CFDkit_OR_MPI(error_flag);
//       /************************************************************************/

//       /**************************************************************************/
//       /* Calculate global epsilon based on 2-norm of x. */  
//       double total_norm_x = ZERO;
//       for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
// 	  total_norm_x += sqr( G[Bcount].L2_Norm(G[Bcount].x) );	  
// 	}
//       } 
//       total_norm_x = sqrt(CFDkit_Summation_MPI(total_norm_x));      
//       epsilon = EPSILON_PARAM/sqrt(total_norm_x);       
//       /**************************************************************************/



//       /**************************************************************************/
//       /******* BEGIN MATRIX-FREE FOR RESTART ************************************/
//       /**************************************************************************/
//       for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  
// 	  //Calculate R(U+epsilon*(Minv*x(i))) -> Soln_ptr.U =  Soln_ptr.Uo + epsilon * W(i)
// 	  G[Bcount].calculate_perturbed_residual_Restart(epsilon);	  
// 	  BCs(Soln_ptr[Bcount],*Input_Parameters);
// 	  dUdt_Residual_Evaluation(Soln_ptr[Bcount],*Input_Parameters);	  
// 	} 
//       } 
//       /**************************************************************************/

//       /**************************************************************************/
//       // Send boundary flux corrections at block interfaces with resolution changes.
//       error_flag = Send_Conservative_Flux_Corrections(Soln_ptr,
// 						      *List_of_Local_Solution_Blocks,
// 						      Soln_ptr[0].NumVar());
//       if (error_flag) {
// 	cout << "\n GMRES ERROR: Solution flux correction message passing error on processor "
// 	     << List_of_Local_Solution_Blocks->ThisCPU << ".\n";
// 	cout.flush();
//       } 
//       error_flag = CFDkit_OR_MPI(error_flag);
      
//       // Apply boundary flux corrections to residual to ensure that method is conservative.
//       Apply_Boundary_Flux_Corrections(Soln_ptr, *List_of_Local_Solution_Blocks);
//       /**************************************************************************/
      
//       /////////////////////////// 2ND /////////////////////////////////// 
//       if(Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == SECOND_ORDER) {
// 	for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	  if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {	    
// 	    //Calculate R(U+epsilon*(Minv*x(i))) -> Soln_ptr.U =  Soln_ptr.Uo + epsilon * W(i)
// 	    G[Bcount].calculate_perturbed_residual_2nd_Restart(epsilon);	  
// 	    BCs(Soln_ptr[Bcount],*Input_Parameters);
// 	    dUdt_Residual_Evaluation(Soln_ptr[Bcount],*Input_Parameters);	  
// 	  } 
// 	} 
// 	/**************************************************************************/
	
// 	/**************************************************************************/
// 	// Send boundary flux corrections at block interfaces with resolution changes.
// 	error_flag = Send_Conservative_Flux_Corrections(Soln_ptr,
// 							*List_of_Local_Solution_Blocks,
// 							Soln_ptr[0].NumVar());
// 	if (error_flag) {
// 	  cout << "\n GMRES ERROR: Solution flux correction message passing error on processor "
// 	       << List_of_Local_Solution_Blocks->ThisCPU << ".\n";
// 	  cout.flush();
// 	} 
// 	error_flag = CFDkit_OR_MPI(error_flag);
	
// 	// Apply boundary flux corrections to residual to ensure that method is conservative.
// 	Apply_Boundary_Flux_Corrections(Soln_ptr, *List_of_Local_Solution_Blocks);
// 	/**************************************************************************/
//       }
//       //////////////////////////////////////////////////////////////////////////
            
//       /**************************************************************************/
//       //Calculate Matrix Free V(0) = ( R(U+epsilon*x) - b) / epsilon - x / h */
//       for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {	  
// 	  G[Bcount].calculate_Matrix_Free_Restart(epsilon);	  	   	  
// 	}
//       }   
//       /**************************************************************************/
//       /****************** END OF MATRIX-FREE FOR RESTART ************************/
//       /**************************************************************************/


//     }  
//     /**************************************************************************/
//     /***********  END OF RESTART APPLIED  *************************************/
//     /**************************************************************************/


//     /**************************************************************************/
//     // CALCULATE NORM OF FIRST SEARCH VECTOR, V(0).
//     beta= ZERO;
//     for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
//       if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
// 	/* V(0) = ONE*(-b) + V(0) = Ax - b */
// 	temp = ONE;
// 	F77NAME(daxpy)(&G[Bcount].scalar_dim, &temp, G[Bcount].b, &inc, 
// 		       &(G[Bcount].V[(0)*G[Bcount].scalar_dim]), &inc);

// 	/* beta = norm(V(0)) */	
// 	beta += sqr(G[Bcount].L2_Norm(&(G[Bcount].V[(0)*G[Bcount].scalar_dim])));
				       
//       }
//     }

//     beta = sqrt(CFDkit_Summation_MPI(beta));
//     /**************************************************************************/


//     /**************************************************************************/
//     //GMRES check compares computed maxtrix-vector products to real Ax by
//     //using the "restart" Ax.
//     if(Input_Parameters->NKS_IP.GMRES_CHECK && !i_first_time_through && !i_restart_flag) {    
//       // is ||beta||/||b|| =  to final relative_residual from GMRES
//       break;
//     }
//     /**************************************************************************/
    

//     /**************************************************************************/
//     // RESSCALE V(0) USING NORM.
//     for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
//       if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	
// 	/* V(0) = -V(0)/beta */
// 	temp = -ONE/beta;
// 	F77NAME(dscal)(&G[Bcount].scalar_dim, &temp, &(G[Bcount].V[(0)*G[Bcount].scalar_dim]), &inc);
	
// 	//reset residual's to zero, restart
// 	for (int i = 1; i <  G[Bcount].restart+1; i++)
// 	  G[Bcount].s[i] = ZERO;  
	
// 	G[Bcount].s[0] = beta;
	
// 	//Set vector switch to pass W 
// 	G[Bcount].vector_switch  = 1;  

//       } 
//     } 
//     /**************************************************************************/
    
//     /* save the very first residual norm */
//     if (Number_of_GMRES_Iterations == 0) {
//       resid0 = beta;
//     }

//     /*************************************************************************/
//     /**************** Begin Primary GMRES Loop *******************************/
//     /*************************************************************************/
//     int search_direction_counter = -1;

//     do {

//       search_direction_counter++;
//       Number_of_GMRES_Iterations++;

//       /**************************************************************************/
//       // Calculate z vector using preconditioner.
//       for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  
// 	  /* Set search direction counter. */
// 	  G[Bcount].search_directions = search_direction_counter;  
	  
// 	  //Apply Block Precondtioner  z = Minv * V(i) -> stored in W(i). 
// 	  Block_precon[Bcount].Apply_Preconditioner(G[Bcount].scalar_dim, 1, &(G[Bcount].V[(search_direction_counter)*G[Bcount].scalar_dim]), 
// 						    G[Bcount].scalar_dim, &(G[Bcount].W[(search_direction_counter)*G[Bcount].scalar_dim]), 
// 						    G[Bcount].scalar_dim);	  
// 	} 
//       } 
//       /**************************************************************************/

    
//       /**************************************************************************/
//       /* MPI barrier to ensure processor synchronization. */
//       CFDkit_Barrier_MPI();  
      
//       /* Send solution information between neighbouring blocks.*/
//       /* Passes "W = Minv * V(i) = z"  information in ghost cells */       
//       error_flag = Send_All_Messages(G, 
// 				     *List_of_Local_Solution_Blocks,
// 				     G[0].NumVar(), 
// 				     OFF); if (error_flag) {
// 	cout << "\n GMRES2D ERROR: GMRES message passing error on processor "
// 	     << List_of_Local_Solution_Blocks->ThisCPU
// 	     << ".\n";
// 	cout.flush();
//       } /* endif */
//       error_flag = CFDkit_OR_MPI(error_flag);
//       /**************************************************************************/
  

//       /**************************************************************************/
//       /* Calculate global epsilon base on 2-norm of z. */
//       double total_norm_z = ZERO; 
//       for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
// 	  total_norm_z += sqr(G[Bcount].L2_Norm(&(G[Bcount].W[(search_direction_counter)*G[Bcount].scalar_dim])));
//  	} 
//       }       
//       total_norm_z = sqrt(CFDkit_Summation_MPI(total_norm_z));     	                      
//       epsilon = EPSILON_PARAM/sqrt(total_norm_z);            
//       /**************************************************************************/



//       /**************************************************************************/
//       /***************** BEGIN MATRIX-FREE FOR PRIMARY GMRES LOOP ***************/
//       /**************************************************************************/ 
           
//       // Calculate perturbed Residual R(U+epsilon*(Minv*V(i)))
//       for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) { 
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) { 
//  	  //Calculate R(U+epsilon*(Minv*V(i))) -> Soln_ptr.U =  Soln_ptr.Uo + epsilon * W(i)
//  	  G[Bcount].calculate_perturbed_residual(epsilon);
// 	  //Apply Regular Soln_ptr BC'S 
// 	  BCs(Soln_ptr[Bcount],*Input_Parameters);
// 	  //modified to calculate in "overlap" cells as well
// 	  dUdt_Residual_Evaluation(Soln_ptr[Bcount],*Input_Parameters);	  
// 	 } 
//       }
//       /**************************************************************************/
    
//       /**************************************************************************/
//       // Send boundary flux corrections at block interfaces with resolution changes. (changes to dUdt)
//       error_flag = Send_Conservative_Flux_Corrections(Soln_ptr,
// 	              			              *List_of_Local_Solution_Blocks,
// 						      Soln_ptr[0].NumVar());
//       if (error_flag) {
// 	cout << "\n GMRES ERROR: Solution flux correction message passing error on processor "
// 	     << List_of_Local_Solution_Blocks->ThisCPU << ".\n"; cout.flush();
//       } 
//       error_flag = CFDkit_OR_MPI(error_flag);
	  
//       // Apply boundary flux corrections to residual to ensure that method is conservative.
//       Apply_Boundary_Flux_Corrections(Soln_ptr,*List_of_Local_Solution_Blocks);
      
//       /**************************************************************************/


//       ////////////////// 2nd ORDER ////////////////////////////////////////
//       if(Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == SECOND_ORDER) {
	 
// 	// Calculate perturbed Residual R(U+epsilon*(Minv*V(i)))
// 	for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) { 
// 	  if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) { 
// 	    //Calculate R(U+epsilon*(Minv*V(i))) -> Soln_ptr.U =  Soln_ptr.Uo - epsilon * W(i)
// 	    G[Bcount].calculate_perturbed_residual_2nd(epsilon);
// 	    //Apply Regular Soln_ptr BC'S 
// 	    BCs(Soln_ptr[Bcount],*Input_Parameters);
// 	    //modified to calculate in "overlap" cells as well
// 	    dUdt_Residual_Evaluation(Soln_ptr[Bcount],*Input_Parameters);	  
// 	  } 
// 	}
// 	/**************************************************************************/
    
// 	/**************************************************************************/
// 	// Send boundary flux corrections at block interfaces with resolution changes.
// 	error_flag = Send_Conservative_Flux_Corrections(Soln_ptr,
// 							*List_of_Local_Solution_Blocks,
// 							Soln_ptr[0].NumVar());
// 	if (error_flag) {
// 	  cout << "\n GMRES ERROR: Solution flux correction message passing error on processor "
// 	       << List_of_Local_Solution_Blocks->ThisCPU << ".\n"; cout.flush();
// 	} 
// 	error_flag = CFDkit_OR_MPI(error_flag);
	
// 	// Apply boundary flux corrections to residual to ensure that method is conservative.
// 	Apply_Boundary_Flux_Corrections(Soln_ptr,*List_of_Local_Solution_Blocks);
	
//       }
//       ////////////////////////////////////////////////////////////////////////////


//       /**************************************************************************/
//       // Calculate Matrix Free V(i+1) = ( R(U + epsilon*W) - R(U) ) / epsilon - z / h
//       for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) { 
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) { 
// 	  G[Bcount].calculate_Matrix_Free(epsilon);	  	  
// 	} 
//       }

//       /**************************************************************************/
//       /******************* END OF MATRIX-FREE ***********************************/
//       /**************************************************************************/


//       /**************************************************************************/
//       // H norm calculation
//       double total_norm_H = ZERO;

//       for (int k = 0; k <= search_direction_counter; k++) {
// 	total_norm_H = ZERO;
// 	for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	  if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
// 	    /* total_H = dot(W,V(k)) */
// 	    total_norm_H += G[Bcount].Dotproduct( &(G[Bcount].V[(search_direction_counter+1)*G[Bcount].scalar_dim]),
// 						  &(G[Bcount].V[(k)*G[Bcount].scalar_dim]) );
// 	  }
// 	}	
// 	total_norm_H = CFDkit_Summation_MPI(total_norm_H);

// 	for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	  if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	    
// 	    /* V(i+1) = V(i+1) -  H(k, i) * V(k) */
// 	    G[Bcount].H[(search_direction_counter)*m1+(k)] = total_norm_H;

// 	    temp = -ONE * G[Bcount].H[(search_direction_counter)*m1+(k)]; //store in temp for passing by address

// 	    F77NAME(daxpy)(&G[Bcount].scalar_dim, &temp, &(G[Bcount].V[(k)*G[Bcount].scalar_dim]), &inc,
// 			   &(G[Bcount].V[(search_direction_counter+1)*G[Bcount].scalar_dim]), &inc);
// 	  } 
// 	} 
//       } 
//       /**************************************************************************/

// #ifdef _NKS_VERBOSE
//       if (CFDkit_Primary_MPI_Processor()) { 
// 	cout << "\n  GMRES (Inner Iterations) = " << Number_of_GMRES_Iterations << " total_norm_H1 " <<   total_norm_H;	 
//       } 
// #endif
 
//       /**************************************************************************/
//       total_norm_H = ZERO;
//       /* Calculate 2-norm of V(i+1) -> H(i+1,1) */
//       for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
// 	  total_norm_H += sqr(G[Bcount].L2_Norm(&(G[Bcount].V[(search_direction_counter+1)*G[Bcount].scalar_dim]))); 	  
// 	} 
//       } 
      	
//       total_norm_H = sqrt(CFDkit_Summation_MPI(total_norm_H)); 
//       /**************************************************************************/

// #ifdef _NKS_VERBOSE   
//       if (CFDkit_Primary_MPI_Processor()) { 
// 	cout << "\n  GMRES (Inner Iterations) = " << Number_of_GMRES_Iterations << " total_norm_H2 " <<   total_norm_H;	 
//       } 
// #endif

//       /**************************************************************************/
//       for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	  
// 	  G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter+1)] = total_norm_H;

// 	  temp = ONE / G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter+1)]; // total_norm_H;
	  
// 	  /* V(i+1) = V(i+1) / H(i+1, i) */
// 	  F77NAME(dscal)(&G[Bcount].scalar_dim, &temp, &(G[Bcount].V[(search_direction_counter+1)*G[Bcount].scalar_dim]), &inc);
	  
// 	  } 
//       } 
//       /**************************************************************************/

      
//       /**************************************************************************/
//       for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	  
// 	  //POSSIBLY USE BLAS1 drtog, drotmg, drot, etc to do these faster ???
	  
// 	  for (int k = 0; k < search_direction_counter; k++) {
// 	    ApplyPlaneRotation(G[Bcount].H[(search_direction_counter)*m1+(k)], 
// 			       G[Bcount].H[(search_direction_counter)*m1+(k+1)], 
// 			       G[Bcount].cs[k], G[Bcount].sn[k]);
// 	  } 
	  
// 	  GeneratePlaneRotation(G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter)],
// 				G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter+1)], 
// 				G[Bcount].cs[search_direction_counter], G[Bcount].sn[search_direction_counter]);
	  
	  
// 	  ApplyPlaneRotation(G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter)],
// 			     G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter+1)], 
// 			     G[Bcount].cs[search_direction_counter], G[Bcount].sn[search_direction_counter]);
	  
	  
// 	  ApplyPlaneRotation(G[Bcount].s[search_direction_counter], G[Bcount].s[search_direction_counter+1], 
// 			     G[Bcount].cs[search_direction_counter], G[Bcount].sn[search_direction_counter]);
	  
// 	} 
//       } 
//       /**************************************************************************/


//       /**************************************************************************/
//       relative_residual = ZERO;
//       for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
// 	  relative_residual = max( fabs(G[Bcount].s[search_direction_counter+1]) , relative_residual); 
// 	}
//       }      
      
//       relative_residual = CFDkit_Maximum_MPI(relative_residual)/resid0;    

//       /**************************************************************************/
      
      
//       /**************************************************************************/
//       // Output progress

//       //Verbose Output for CHECKING
//       if (CFDkit_Primary_MPI_Processor() && Number_of_GMRES_Iterations%5 == 0 && Input_Parameters->NKS_IP.GMRES_CHECK) { 
// 	if(Number_of_GMRES_Iterations == 5){	
// 	  cout << "\n  GMRES Iter.  \t   resid0 \t   resid \t  rel_resid  \t   L2||z||   \t  epsilon ";
// 	} 
// 	cout << "\n \t" << Number_of_GMRES_Iterations << "\t" << resid0 << "\t" 
// 	     << relative_residual*resid0 << "\t" << relative_residual
// 	     <<"\t"<< total_norm_z << "\t  "<<epsilon;	
//       } 

//       //Standard Output
//       if (relative_residual <= Input_Parameters->NKS_IP.GMRES_Tolerance || 
// 	  Number_of_GMRES_Iterations+1 > Input_Parameters->NKS_IP.Maximum_Number_of_GMRES_Iterations ) {

// 	if (CFDkit_Primary_MPI_Processor()) { 
// 	  cout << "\n Finished GMRES with (Inner Iterations) = " << Number_of_GMRES_Iterations << " resid0 = " 
// 	       << resid0 << " resid = " << relative_residual*resid0 << " relative_residual = " << relative_residual<<endl;
// 	} 
// 	break;
//       }       
//       /**************************************************************************/

//     } while (search_direction_counter+1 < Input_Parameters->NKS_IP.GMRES_Restart && 
// 	     Number_of_GMRES_Iterations+1 <= Input_Parameters->NKS_IP.Maximum_Number_of_GMRES_Iterations);
//     /**************************************************************************/
//     /******************* END PRIMARY GMRES LOOP *******************************/
//     /**************************************************************************/


//     /************************************************************************/
//     // UPDATE SOLUTION. 
//     for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
//       if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){  
		
// 	/* solve upper triangular system in place */	  
// 	for (int j = search_direction_counter; j >= 0; j--) { 
// 	  G[Bcount].s[j] /= G[Bcount].H[(j)*m1+(j)];
// 	  for (int k = j-1; k >= 0; k--)
// 	    G[Bcount].s[k] -= G[Bcount].H[(j)*m1+(k)] * G[Bcount].s[j];
// 	}
	
// 	/* update solution, x */ 
// 	for (int j = 0; j <= search_direction_counter; j++) {
// 	  /* x = x + s[j] * W(j) */
// 	  F77NAME(daxpy)(&G[Bcount].scalar_dim, &G[Bcount].s[j], &(G[Bcount].W[(j)*G[Bcount].scalar_dim]), 
// 			 &inc, G[Bcount].x, &inc); 
// 	}   
//       } 
//     } 
//     /************************************************************************/

//     i_first_time_through = false;

//     // flag for restart
//     if( relative_residual > Input_Parameters->NKS_IP.GMRES_Tolerance && 
//         Number_of_GMRES_Iterations+1 <= Input_Parameters->NKS_IP.Maximum_Number_of_GMRES_Iterations) {
//       i_restart_flag = true;
//     } else {
//       i_restart_flag = false;
//     }

//   } while ( i_restart_flag || Input_Parameters->NKS_IP.GMRES_CHECK);
 
//   /**************************************************************************/
//   /********************* END GMRES CALCULATION ******************************/
//   /**************************************************************************/


//   /**************************************************************************/
//   // Check GMRES for accurate matrix vector products 
//   // ie. does ||(Ax-b)||/||b|| ~= relative_residual
  
//   if(Input_Parameters->NKS_IP.GMRES_CHECK){

//     int blocksize = G[0].blocksize;
//     double total_norm_b(ZERO);
//     double *total_norm_eqn_b = new double[blocksize];
//     double *total_norm_eqn_r = new double[blocksize];
//     for(int i=0; i< blocksize; i++){  total_norm_eqn_r[i] = ZERO; total_norm_eqn_b[i] =ZERO; }

//     //L2 norms 
//     for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
//       if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
// 	total_norm_b += sqr(G[Bcount].L2_Norm( G[Bcount].b));	
// 	for(int i=0; i< blocksize; i++){	
// 	  total_norm_eqn_b[i] += sqr(G[Bcount].L2_Norm(i, G[Bcount].b));	
// 	  total_norm_eqn_r[i] += sqr(G[Bcount].L2_Norm(i, G[Bcount].V));		  
// 	}
//       }
//     } 
     
//     total_norm_b = sqrt(CFDkit_Summation_MPI(total_norm_b));
//     for(int i=0; i< blocksize; i++){
//       total_norm_eqn_r[i] = sqrt(CFDkit_Summation_MPI(total_norm_eqn_r[i]));
//       total_norm_eqn_b[i] = sqrt(CFDkit_Summation_MPI(total_norm_eqn_b[i]));
//     }

//     //Output Reduction for systems & per equation
//     if (CFDkit_Primary_MPI_Processor()) {

//       cout<<" GMRES Check total   ||(Ax-b)||/||b|| = "<<beta<<" / "<<total_norm_b<<" = "<<beta/total_norm_b
// 	  <<" /relative_residual (ideal 1) -> " <<((beta/total_norm_b)/relative_residual)<<endl;

//        for(int i=0; i< blocksize; i++){	
// 	 cout<<" Reduction for Eqn "<<i<<" ||(Ax-b)||/||b|| = "<<total_norm_eqn_r[i]
// 	     <<" / "<<total_norm_eqn_b[i]<<" = "<<total_norm_eqn_r[i]/total_norm_eqn_b[i]<<endl;
//        }
//     }   

//     delete[] total_norm_eqn_b; delete[] total_norm_eqn_r; 
//   }
//   /**************************************************************************/


//   /**************************************************************************/
//   if ( Number_of_GMRES_Iterations >= Input_Parameters->NKS_IP.Maximum_Number_of_GMRES_Iterations ) {
//     if (CFDkit_Primary_MPI_Processor()) {  
//       cout << "\n GMRES2D - GMRES ERROR: Unable to reach the specified convergence tolerance.";
//       cout << "\n Final gmrestol = " << relative_residual << endl;
//       cout.flush();
//     }
//   }
//   /**************************************************************************/


  /**************************************************************************/
  return error_flag;

} /* End of GMRES_RightPrecon_MatrixFree::solve. */ 


#endif // _GMRES3D_INCLUDED


// #ifdef _NKS_VERBOSE 
//       G[0].Output_W(); 
//       G[0].Output_V();

//       double total_norm_b =ZERO;
//       double total_norm_Vp1 =ZERO;      
     
//       for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
// 	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
// 	  total_norm_b += sqr(G[Bcount].L2_Norm( G[Bcount].b));
// 	  total_norm_Vp1 += sqr(G[Bcount].L2_Norm(&(G[Bcount].V[(search_direction_counter+1)*G[Bcount].scalar_dim])));	
// 	} 
//       } 
      
//       total_norm_b = sqrt(CFDkit_Summation_MPI(total_norm_b));
//       total_norm_Vp1 = sqrt(CFDkit_Summation_MPI(total_norm_Vp1));     	    
      
//       if (CFDkit_Primary_MPI_Processor()) { 
//  	cout << "\n  GMRES (Inner Iterations) = " << Number_of_GMRES_Iterations << " total_norm_b " <<  total_norm_b;	
// 	cout << "\n  GMRES (Inner Iterations) = " << Number_of_GMRES_Iterations << " total_norm_Vp1 " <<  total_norm_Vp1;	 
//       } 
// #endif

