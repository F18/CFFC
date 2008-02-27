#ifndef _GMRES2D_INCLUDED
#define _GMRES2D_INCLUDED

#include <iostream>
#include <cstdio>
#include <cmath> 

#ifndef _BLAS1_H_     
#include "blas1.h"   //Wrappers for BLAS Fortran libraries, so need to link with -lblas
#endif // _BLAS1_H_

#ifndef _BLOCK_PRECONDITONER_INCLUDED 
#include "Block_Preconditioner2D.h"
#endif 

#ifndef _NKS_DTS_INCLUDED
#include "DTS_NKS.h"
#endif

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
 * scalar_dim   -- Return NCi * NCj * blocksize            *
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
 *  SolnBlk  -- Return pointer to solution block.          *
 *                                                         *
 *                                                         *
 ***********************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
class GMRES_Block{
private:

  /*********************************************/
  double *normalize_valuesR;      //vector(s) of normalizing values corresponding to SOLN_BLOCK_TYPE
  double *normalize_valuesU;
  void set_normalize_values(void);   //set normalization values specific to SOLN_BLOCK_TYPE
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

  //Solution Block i,j,k indexing conversion to 1D Vector 
  int index(int i, int j, int k) {return ((j*NCi+i)*blocksize+k);}
  int index(int i, int j) {return ((j*NCi+i)*blocksize);}

public:
  /* Solutition INFORMATION THAT DOESN'T CHANGE , could be static ???*/
  int                     restart; // number of gmres iterations before restart
  int                     overlap; // level of overlap  
  int                   blocksize; // number of variables            //Soln_ptr->NumVar();
  int                  scalar_dim; // xpts * ypts * blocksize

  /* USED INTERNAL TO GMRES ROUTINE */
  int               vector_switch; // to select the specified vector for message passing
  int           search_directions; // number of search directions

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
  int                      Nghost; // Number of ghost cells.
  Grid2D_Quad_Block          Grid; // dummy pointer.
  SOLN_BLOCK_TYPE  *SolnBlk;      // Pointer to solution block. 
  INPUT_TYPE *Input_Parameters;

  // Dual Time Stepping
  DTS_NKS_Quad_Block<SOLN_BLOCK_TYPE,INPUT_TYPE> *DTS_ptr;


  /****************************************************************/
  /************** Creation and copy constructors. *****************/
  /****************************************************************/  
  GMRES_Block(void)  {
    restart = 0;   overlap = 0;  blocksize = 0;
    scalar_dim = 0; search_directions = 0; 
    normalize_valuesR = NULL; normalize_valuesU = NULL; 
    s = NULL;  cs = NULL;  sn = NULL;  
    b = NULL;  x = NULL; H = NULL;  W = NULL;  V = NULL; 
    vector_switch = 100; NCi = 0;  ICl = 0; ICu = 0; 
    NCj = 0;  JCl = 0; JCu = 0; Nghost = 0;
    SolnBlk = NULL;      
    DTS_ptr = NULL;
  }
  
  // GMRES_Block(const GMRES_Block &G); //FIX so that it actually copies, not just passes pointers!!!!!!!!!!!!
  // GMRES_Block operator = (const GMRES_Block &G); //Setup proper assignment operator
  
  /* Allocate and deallocate memory for GMRES_Block variables. */
  void allocate( const int m, const int overlap_cells, 
		 const bool normalize, SOLN_BLOCK_TYPE &Soln_Block_ptr, 
		 INPUT_TYPE &IP, const int &_blocksize, 
		 DTS_NKS_Quad_Block<SOLN_BLOCK_TYPE,INPUT_TYPE> &DTS);  
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

  //TESTING ie _NKS_VERBOSE functions
  void Output_W(void) { for (int i=0;i< restart*scalar_dim;++i){ cout<<"\n W(i) "<<W[i]; }}
  void Output_V(void) { for (int i=0;i<((restart + 1) * scalar_dim);++i) { cout<<"\n V(i) "<<V[i]; }}
  void Output_U(int what);
  void Output_GMRES_vars_Tecplot(const int Number_of_Time_Steps,
			         const int Block_Number,
			         bool print_title,
			         const double &l2_norm,
			         const double &l2_norm_rel,
			         ofstream &fout);

  //Norm Calculations  and Vector Operations
  double L2_Norm(const double* v);
  double L2_Norm(int N, const double* v, int inc);
  double L2_Norm(int k, const double* v);
  double Dotproduct(const double *v1, const double *v2);

  // norm of normalized solution vector
  double L1_Norm_Unorm(void);

  //Return denormalized deltaU
  double deltaU(const int i, const int j, const int k){
    return ( denormalizeU(x[index(i,j,k)],k) ); }  

  //Return info for testing/degugging.
  double deltaU_test(const int i, const int j, const int k){ return x[index(i,j,k)]; }   
  double b_test(const int i, const int j, const int k){ return b[index(i,j,k)]; }   

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
                         const int j_inc,
			 const int face,
			 const int sector);
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
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
set_normalize_values(void) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF set_normalize_values for GMRES2D.h requried \n";
  exit(1);
}

/**************************************************************************
 * Routine: TESTING FUNCITON                                              *
 **************************************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Output_U(int what) {    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      if(what == 1)cout<<"\n U ("<<i<<","<<j<<") = "<<SolnBlk->U[i][j];      
      if(what == 2)cout<<"\n W ("<<i<<","<<j<<") = "<<SolnBlk->W[i][j];
      if(what == 3)cout<<"\n dUdt ("<<i<<","<<j<<") = "<<SolnBlk->dUdt[i][j][0];
    }
  }  
}

template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Output_GMRES_vars_Tecplot (const int Number_of_Time_Steps,
		           const int Block_Number,
		           bool print_title,
		           const double &l2_norm,
		           const double &l2_norm_rel,
		           ofstream &fout) {

  // whether or not to print the ghost cells. Use 1 or 0.
  int ng = 1;
  int kvar = Input_Parameters->i_Residual_Variable - 1;

  if (print_title) {
    fout << "TITLE = \"";
    fout << "GMRES deltaU: Iteration = " << Number_of_Time_Steps
	 << scientific << setprecision(10) << ", L2-Norm = " << l2_norm
	 << scientific << setprecision(3) << ", L2-Norm Rel = " << l2_norm_rel << fixed
	 << setprecision(14);
    fout << "\"" << "\n"
	 << "VARIABLES = \"x\" \\ \n"
	 << "\"y\" \\ \n";
    fout << "\"b_k" << kvar	<< "\" \\ \n";
    for (int q = 0; q < blocksize; q++) {
	fout << "\"deltaU_k" << q << "\" \\ \n";
    }
    for (int q = 0; q < blocksize; q++) {
	fout << "\"V0_" << q << "\" \\ \n";
    }
  }

  fout << setprecision(14);

  fout << "ZONE T =  \"Block " << Block_Number << " CPU" << CFFC_MPI::This_Processor_Number
       << "\" \\ \n"
       << "I = " << SolnBlk->Grid.ICu - SolnBlk->Grid.ICl + 2*SolnBlk->Nghost*ng + 1 << " \\ \n"
       << "J = " << SolnBlk->Grid.JCu - SolnBlk->Grid.JCl + 2*SolnBlk->Nghost*ng + 1 << " \\ \n"
       << "F = POINT \n";

  fout.setf(ios::scientific);
  for (int j = JCl - Nghost*ng; j <= JCu + Nghost*ng; j++) {
      for (int i = ICl - Nghost*ng; i <= ICu + Nghost*ng; i++) {
          fout << " " << SolnBlk->Grid.Cell[i][j].Xc;
	  fout << " " << fabs(b[index(i,j,kvar)]); 
	  for (int q = 0; q < blocksize; q++) {
	      fout << " " << fabs(deltaU(i,j,q)); 
	  }
	  for (int q = 0; q < blocksize; q++) {
	      fout << " " << fabs(V[0 + index(i, j, q)]);
	  }
      fout << endl;
      }
  }

  fout << setprecision(6);
}

/**************************************************************************
 * GMRES_Block::allocate -- Allocate memory.                              *
 **************************************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
allocate( const int m, const int overlap_cells, 
	  const bool normalize, SOLN_BLOCK_TYPE &Soln_Block_ptr, 
	  INPUT_TYPE &IP, const int &_blocksize, 
	  DTS_NKS_Quad_Block<SOLN_BLOCK_TYPE,INPUT_TYPE> &DTS)
{  
  //set pointer to solution block and set all necessary GMRES information
  SolnBlk = &(Soln_Block_ptr);
  Input_Parameters = &(IP);
  DTS_ptr = &(DTS);
  NCi = SolnBlk->NCi;
  ICl = SolnBlk->ICl;
  ICu = SolnBlk->ICu;
  NCj = SolnBlk->NCj;
  JCl = SolnBlk->JCl;
  JCu = SolnBlk->JCu;
  Nghost = SolnBlk->Nghost;
  overlap = overlap_cells;  
  restart = m;


  //Check if a valid solution block and GMRES parameters
  assert(restart > 1);  assert(NCi > 1);  assert(NCj > 1);
  blocksize = _blocksize; 
	// scalar_dim could be a static variable, yes...?
  scalar_dim = NCi * NCj * blocksize;

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
   if(normalize){
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
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
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
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Initialize(void)
{
  /* Initialize all GMRES variables except b, which is copied from dUdt */  
  for (int i=0;i<restart;++i) {  cs[i] = ZERO; sn[i] = ZERO; s[i] = ZERO;  }
  for (int i=0;i<scalar_dim;++i){  x[i] = ZERO; b[i] = ZERO; }
  for (int i=0;i<(restart*(restart+1));++i)         H[i] = ZERO;
  for (int i=0;i< restart*scalar_dim;++i)           W[i] = ZERO;
  for (int i=0;i<((restart + 1) * scalar_dim);++i)  V[i] = ZERO;
  
  // Load b ie RHS ie dUodt and normalize
  for (int j = JCl;  j <= JCu; j++){
    for (int i = ICl;  i <= ICu; i++){
      for(int varindex =0; varindex < blocksize; varindex++){		
	b[index(i,j,varindex)] = normalizeR(SolnBlk->dUdt[i][j][0][varindex+1],varindex);
      }    	      
    } 
  }      
}

/**************************************************************************
 * Routine: calculate_pertubed_residual                                   *
 **************************************************************************/
// Calculate Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * W(i) )
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
calculate_perturbed_residual(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*W[search_directions*scalar_dim+index(i,j,varindex)], varindex); 
      }   
      /* Update primitive variables. */
      SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
    }
  }  
}

// Copy forward difference & calculate backwards for 2nd order derivative 
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
calculate_perturbed_residual_2nd(const double &epsilon) {    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      //copy back R + epsilon * W(i)
      SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];

      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] - 
	  denormalizeU( epsilon*W[search_directions*scalar_dim+index(i,j,varindex)], varindex); 
      }   
      /* Update primitive variables. */
      SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
    }
  }  
}

// Calculate Restart Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * x(i) )
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
calculate_perturbed_residual_Restart(const double &epsilon) {    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  
      /* Update primitive variables. */
      SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
    }
  }  
}

// Copy forward difference & calculate backwards for 2nd order derivative  
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
calculate_perturbed_residual_2nd_Restart(const double &epsilon) {    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) { 
      //copy back R + epsilon * W(i)
      SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];
      
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] -
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  
      /* Update primitive variables. */
      SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
    }
  }  
}

/********************************************************
 * Routine: calculate_Matrix_Free                       *
 ********************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
calculate_Matrix_Free(const double &epsilon) {
  //Taking into acount NKS overlap
  int JCl_overlap = 0; int JCu_overlap = 0;
  int ICu_overlap = 0; int ICl_overlap = 0;		  
  if(overlap){	
    if ( SolnBlk->Grid.BCtypeS[ICl] == BC_NONE)  JCl_overlap = overlap; 
    if ( SolnBlk->Grid.BCtypeN[ICu] == BC_NONE)  JCu_overlap = overlap;
    if ( SolnBlk->Grid.BCtypeE[JCu] == BC_NONE)  ICu_overlap = overlap;
    if ( SolnBlk->Grid.BCtypeW[JCl] == BC_NONE)  ICl_overlap = overlap;
  }
	   
  // Non-Overlap Ghost Cells R(U) already set to zero by dUdt calculation  
  /* V(i+1) = ( R(U+epsilon*W) - b) / epsilon - z / h */
  for (int j = JCl - JCl_overlap; j <= JCu + JCu_overlap; j++) {
    for (int i = ICl - ICl_overlap; i <= ICu + ICu_overlap; i++) {
      for(int k =0; k < blocksize; k++){	
	int iter = index(i,j,k);		
	//Matrix Free V(i+1) 
	if( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == FIRST_ORDER ){
	  //forwards differenceing R(U+epsilon) - R(U) / epsilon
	  V[(search_directions+1)*scalar_dim+iter] = 
	    ( normalizeR(SolnBlk->dUdt[i][j][0][k+1],k) - b[iter]) / epsilon ;
	} else if ( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == SECOND_ORDER ){
	  //2nd order R(U+epsilon) - R(U-epsilon) / 2*epsilon
	  V[(search_directions+1)*scalar_dim+iter] = 
	    normalizeR( SolnBlk->dUdt[i][j][1][k+1] - SolnBlk->dUdt[i][j][0][k+1],k )/(TWO*epsilon);
	}
	V[(search_directions+1)*scalar_dim+iter] -=  
	  normalizeUtoR( W[(search_directions)*scalar_dim + iter] * 
			 LHS_Time<INPUT_TYPE>(*Input_Parameters,SolnBlk->dt[i][j],DTS_ptr->DTS_dTime), k);

// #ifdef _NKS_VERBOSE_NAN_CHECK
// 	// nan check most commonly caused by nans in dUdt !!!!
// 	if (V[(search_directions+1)*scalar_dim+iter] != V[(search_directions+1)*scalar_dim+iter] ){
// 	  cout<<"\n nan in V[ "<<(search_directions+1)*scalar_dim+iter<<"] at "<<i<<" "<<j<<" "<<k
// 	      <<" dUdt "<<  normalizeR(SolnBlk->dUdt[i][j][0][k+1],k) <<" b "<< b[iter]
// 	      <<" z "<<W[(search_directions)*scalar_dim + iter]<< " h "<<( SolnBlk->dt[i][j]*ao);
// 	}
// #endif 
      }      
    } 
  } 
}

template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
calculate_Matrix_Free_Restart(const double &epsilon) {
  //Taking into acount NKS overlap
  int JCl_overlap = 0; int JCu_overlap = 0;
  int ICu_overlap = 0; int ICl_overlap = 0;		  
  if(overlap){	
    if ( SolnBlk->Grid.BCtypeS[ICl] == BC_NONE)  JCl_overlap = overlap; 
    if ( SolnBlk->Grid.BCtypeN[ICu] == BC_NONE)  JCu_overlap = overlap;
    if ( SolnBlk->Grid.BCtypeE[JCu] == BC_NONE)  ICu_overlap = overlap;
    if ( SolnBlk->Grid.BCtypeW[JCl] == BC_NONE)  ICl_overlap = overlap;
  }
	  
  // Non-Overlap Ghost Cells R(U) already set to zero by dUdt   
  /* V(0) = ( R(U + epsilon*W) - b) / epsilon - x / h */
  for (int j = JCl - JCl_overlap; j <= JCu + JCu_overlap; j++) {
    for (int i = ICl - ICl_overlap; i <= ICu + ICu_overlap; i++) {
      for(int k =0; k < blocksize; k++){	
	int iter = index(i,j,k);		
	//Matrix Free V(0) 
	if( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == FIRST_ORDER ){
	  //forwards differenceing R(U+epsilon) - R(U) / epsilon
	  V[iter] = (normalizeR(SolnBlk->dUdt[i][j][0][k+1],k) - b[iter]) / epsilon ;
	} else if ( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == SECOND_ORDER ){
	  //2nd order R(U+epsilon) - R(U-epsilon) / 2*epsilon
	  V[iter] = normalizeR( SolnBlk->dUdt[i][j][1][k+1] - SolnBlk->dUdt[i][j][0][k+1],k)/(TWO*epsilon);
	}
	V[iter] -= normalizeUtoR( x[iter] * LHS_Time<INPUT_TYPE>(*Input_Parameters,
								 SolnBlk->dt[i][j],
								 DTS_ptr->DTS_dTime), k); 
      }      
    } 
  } 
}

/********************************************************
 * Routine: L2_norms                                    *
 ********************************************************/
// template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
// inline double GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
// L2_Norm_noghost(const double* v) {   
//   int ind;
//   double l2_norm =ZERO;
//   for (int j = JCl ; j <= JCu ; j++) {
//     for (int i = ICl ; i <= ICu ; i++) {
//       for(int k = 1; k <= blocksize; k++){	
// 	ind = index(i,j,k);
// 	l2_norm += v[ind]*v[ind]; 
//       }        
//     }
//   }  
//   return(sqrt(l2_norm));
// }

template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline double GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
L2_Norm(const double* v) {
  integer inc = 1;
  return ( F77NAME(dnrm2)( &scalar_dim, v, &inc) );
}

template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline double GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
L2_Norm(int N, const double* v, int inc) {
  return ( F77NAME(dnrm2)( &N, v, &inc) );
}

template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline double GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
L2_Norm(int k, const double* v){

  int ind;
  double l2_norm(ZERO);
  for (int j = JCl ; j <= JCu ; j++) {
    for (int i = ICl ; i <= ICu ; i++) {
      ind = index(i,j,k);
      l2_norm += v[ind]*v[ind];             
    }
  }  
  return(sqrt(l2_norm));

}

/********************************************************
 * Routine: Dot_Product                                 *
 ********************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline double GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Dotproduct(const double *v1, const double *v2) {
  integer inc = 1;
  return (F77NAME(ddot)( &scalar_dim, v1, &inc, v2, &inc));	   
}


/********************************************************
 * L1 Norm of normalized solution variable.             *
 ********************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline double GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
L1_Norm_Unorm(void) {
  double l1_norm_u(ZERO); 
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++)
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++)
      for(int varindex = 0; varindex < blocksize; varindex++)
	l1_norm_u += fabs(normalizeU( SolnBlk->U[i][j][varindex+1], varindex));
  return l1_norm_u;
}

/*******************************************************************************
 *                                                                             *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                              *
 *                                                                             *
 *******************************************************************************/

/*******************************************************************************
 * GMRES_Block::NumVar -- Returns number of state variables.                   *
 *******************************************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline int GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
NumVar(void) {
  return (blocksize);             
}

/*******************************************************************************
 * GMRES_Block::LoadSendBuffer -- Loads send message buffer.                   *
 *******************************************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline int GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
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
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline int GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
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
	  buffer[buffer_count] = (SolnBlk->Grid.Cell[i  ][j  ].A*W[(search_directions)*scalar_dim + index(i  ,j  ,k)]+
				  SolnBlk->Grid.Cell[i+1][j  ].A*W[(search_directions)*scalar_dim + index(i+1,j  ,k)]+
				  SolnBlk->Grid.Cell[i  ][j+1].A*W[(search_directions)*scalar_dim + index(i  ,j+1,k)]+
				  SolnBlk->Grid.Cell[i+1][j+1].A*W[(search_directions)*scalar_dim + index(i+1,j+1,k)])/
	                         (SolnBlk->Grid.Cell[i  ][j  ].A+
				  SolnBlk->Grid.Cell[i+1][j  ].A+
				  SolnBlk->Grid.Cell[i  ][j+1].A+
				  SolnBlk->Grid.Cell[i+1][j+1].A);
	} else {
	  buffer[buffer_count] = (SolnBlk->Grid.Cell[i  ][j  ].A*x[index(i  ,j  ,k)]+
				  SolnBlk->Grid.Cell[i+1][j  ].A*x[index(i+1,j  ,k)]+
				  SolnBlk->Grid.Cell[i  ][j+1].A*x[index(i  ,j+1,k)]+
				  SolnBlk->Grid.Cell[i+1][j+1].A*x[index(i+1,j+1,k)])/
	                         (SolnBlk->Grid.Cell[i  ][j  ].A+
				  SolnBlk->Grid.Cell[i+1][j  ].A+
				  SolnBlk->Grid.Cell[i  ][j+1].A+
				  SolnBlk->Grid.Cell[i+1][j+1].A);
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
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline int GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
LoadSendBuffer_C2F(double *buffer,
		   int &buffer_count,
		   const int buffer_size,
		   const int i_min, 
		   const int i_max,
		   const int i_inc,
		   const int j_min, 
		   const int j_max,
		   const int j_inc,
		   const int face,
		   const int sector) {
  int i, j, k;
  Vector2D dX;
  SOLN_VAR_TYPE Wcoarse, Wfine;
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
	  dX = (SolnBlk->Grid.Node[i][j_min].X+
		HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
		HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
		SolnBlk->Grid.Cell[i][j_min].Xc)/FOUR -
	        SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;

	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } 
	  
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
		SolnBlk->Grid.Node[i+1][j_min].X + SolnBlk->Grid.Cell[i][j_min].Xc+
		HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
	        SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
            
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
	  dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
                    SolnBlk->Grid.Cell[i][j_min].Xc+
                    SolnBlk->Grid.Node[i][j_min+1].X+
                    HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
          
	      for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolnBlk->Grid.Cell[i][j_min].Xc+
                    HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
                    SolnBlk->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
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
              dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
                    SolnBlk->Grid.Node[i+1][j_min].X+
                    SolnBlk->Grid.Cell[i][j_min].Xc+
                    HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolnBlk->Grid.Node[i][j_min].X+
                    HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
                    HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
                    SolnBlk->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
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
              dX = (SolnBlk->Grid.Cell[i][j_min].Xc+
                    HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
                    SolnBlk->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
                    SolnBlk->Grid.Cell[i][j_min].Xc+
                    SolnBlk->Grid.Node[i][j_min+1].X+
                    HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
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
              dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
                    SolnBlk->Grid.Cell[i][j_min].Xc+
                    SolnBlk->Grid.Node[i][j_min+1].X+
                    HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolnBlk->Grid.Cell[i][j_min].Xc+
                    HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
                    SolnBlk->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
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
              dX = (SolnBlk->Grid.Node[i][j_min].X+
                    HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
                    HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
                    SolnBlk->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
              for ( k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
                    SolnBlk->Grid.Node[i+1][j_min].X+
                    SolnBlk->Grid.Cell[i][j_min].Xc+
                    HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
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
              dX = (SolnBlk->Grid.Cell[i][j_min].Xc+
                    HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
                    SolnBlk->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;

	      // Wfine = Wfine.U();  ??MISTAKE MADE IN KALVINS ORIGINAL ?????

             for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
                    SolnBlk->Grid.Cell[i][j_min].Xc+
                    SolnBlk->Grid.Node[i][j_min+1].X+
                    HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
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
              dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
                    SolnBlk->Grid.Node[i+1][j_min].X+
                    SolnBlk->Grid.Cell[i][j_min].Xc+
                    HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
              for ( k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolnBlk->Grid.Node[i][j_min].X+
                    HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
                    HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
                    SolnBlk->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolnBlk->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
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
              dX = (SolnBlk->Grid.Node[i_min][j].X+
                    HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
                    HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
                    SolnBlk->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
                    SolnBlk->Grid.Node[i_min+1][j].X+
                    SolnBlk->Grid.Cell[i_min][j].Xc+
                    HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
                    SolnBlk->Grid.Cell[i_min][j].Xc+
                    SolnBlk->Grid.Node[i_min][j+1].X+
                    HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolnBlk->Grid.Cell[i_min][j].Xc+
                    HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
                    SolnBlk->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
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
              dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
                    SolnBlk->Grid.Node[i_min+1][j].X+
                    SolnBlk->Grid.Cell[i_min][j].Xc+
                    HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolnBlk->Grid.Node[i_min][j].X+
                    HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
                    HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
                    SolnBlk->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolnBlk->Grid.Cell[i_min][j].Xc+
                    HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
                    SolnBlk->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
                    SolnBlk->Grid.Cell[i_min][j].Xc+
                    SolnBlk->Grid.Node[i_min][j+1].X+
                    HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
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
              dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
                    SolnBlk->Grid.Cell[i_min][j].Xc+
                    SolnBlk->Grid.Node[i_min][j+1].X+
                    HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolnBlk->Grid.Cell[i_min][j].Xc+
                    HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
                    SolnBlk->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolnBlk->Grid.Node[i_min][j].X+
                    HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
                    HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
                    SolnBlk->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
                    SolnBlk->Grid.Node[i_min+1][j].X+
                    SolnBlk->Grid.Cell[i_min][j].Xc+
                    HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
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
              dX = (SolnBlk->Grid.Cell[i_min][j].Xc+
                    HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
                    SolnBlk->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
                    SolnBlk->Grid.Cell[i_min][j].Xc+
                    SolnBlk->Grid.Node[i_min][j+1].X+
                    HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
                    SolnBlk->Grid.Node[i_min+1][j].X+
                    SolnBlk->Grid.Cell[i_min][j].Xc+
                    HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolnBlk->Grid.Node[i_min][j].X+
                    HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
                    HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
                    SolnBlk->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolnBlk->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
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
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline int GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
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
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline int GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
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
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline int GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
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
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
inline void GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
SubcellReconstruction(const int i, 
		      const int j,
		      const int Limiter) {
  
  int n, n2, n_pts, i_index[8], j_index[8], k;
  double u0, u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  SOLN_VAR_TYPE U0, DU, DUDx_ave, DUDy_ave, W_VACUUM;
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
//     } else if (SolnBlk->Grid.BCtypeW[j] == BC_PERIODIC ||
//                  SolnBlk->Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
//                  SolnBlk->Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk->Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
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
//                (SolnBlk->Grid.BCtypeE[j] != BC_NONE)) {
//       if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
//          n_pts = 0;
//       } else if (SolnBlk->Grid.BCtypeE[j] == BC_PERIODIC ||
//                  SolnBlk->Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//                  SolnBlk->Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk->Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
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
//                (SolnBlk->Grid.BCtypeS[i] != BC_NONE)) {
//       if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
//          n_pts = 0;
//       } else if (SolnBlk->Grid.BCtypeS[i] == BC_PERIODIC ||
//                  SolnBlk->Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//                  SolnBlk->Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk->Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
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
//                (SolnBlk->Grid.BCtypeN[i] != BC_NONE)) {
//       if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
//          n_pts = 0;
//       } else if (SolnBlk->Grid.BCtypeN[i] == BC_PERIODIC ||
//                  SolnBlk->Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//                  SolnBlk->Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk->Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
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
	dX = SolnBlk->Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - SolnBlk->Grid.Cell[i][j].Xc;
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

      SolnBlk->dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                            (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      SolnBlk->dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                            (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
      // Calculate slope limiters. 
      if (!SolnBlk->Freeze_Limiter) {
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
    
            dX = SolnBlk->Grid.xfaceE(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[0] = u0 + 
                       SolnBlk->dWdx[i][j][n]*dX.x +
                       SolnBlk->dWdy[i][j][n]*dX.y ;
            dX = SolnBlk->Grid.xfaceW(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[1] = u0 + 
                       SolnBlk->dWdx[i][j][n]*dX.x +
                       SolnBlk->dWdy[i][j][n]*dX.y ;
            dX = SolnBlk->Grid.xfaceN(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[2] = u0 + 
                       SolnBlk->dWdx[i][j][n]*dX.x +
                       SolnBlk->dWdy[i][j][n]*dX.y ;
            dX = SolnBlk->Grid.xfaceS(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[3] = u0 + 
                       SolnBlk->dWdx[i][j][n]*dX.x +
                       SolnBlk->dWdy[i][j][n]*dX.y ;
    
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

	    SolnBlk->phi[i][j][n] = phi_n;

         } /* endfor */
      } /* endif */
  } else {
      SolnBlk->dWdx[i][j] = W_VACUUM;
      SolnBlk->dWdy[i][j] = W_VACUUM; 
      SolnBlk->phi[i][j]  = W_VACUUM;
  } /* endif */

}

// /*******************************************************************************
//  * GMRES_Block::LoadSendBuffer_Flux_F2C -- Loads send message buffer for       *
//  *                                         fine to coarse block message        *
//  *                                         passing of conservative             *
//  *                                         solution fluxes.                    *
//  *******************************************************************************/
// template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
// inline int GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
// LoadSendBuffer_Flux_F2C(double *buffer,
// 			int &buffer_count,
// 			const int buffer_size,
// 			const int i_min, 
// 			const int i_max,
// 			const int i_inc,
// 			const int j_min, 
// 			const int j_max,
// 			const int j_inc) {
//   return(0); //SOMETHING MISSING ?????
// }

// /*******************************************************************************
//  * GMRES_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message        *
//  *                                              buffer for fine to coarse      *
//  *                                              block message passing of       *
//  *                                              conservative solution fluxes.  *
//  *******************************************************************************/
// template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE> 
// inline int GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
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
 *                                                         *
 *                                                         *
 *                                                         *
 ***********************************************************/
template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
class GMRES_RightPrecon_MatrixFree {

private:
  
  //Pointers to Solution Info and DATA
  SOLN_BLOCK_TYPE  *Soln_ptr;
  AdaptiveBlock2D_List *List_of_Local_Solution_Blocks;
  INPUT_TYPE *Input_Parameters;

  DTS_NKS_Quad_Block<SOLN_BLOCK_TYPE,INPUT_TYPE> *DTS_ptr;

  //GMRES LOCAL DYNAMIC DATA
  int Number_of_GMRES_Iterations;
  double relative_residual;
  double global_time_step_size;
  GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> *G;                  

  //private GMRES helper functions
  void GeneratePlaneRotation(const double &dx,const double &dy, double &cs, double &sn);
  void ApplyPlaneRotation(double &dx, double &dy,const double &cs,const double &sn);

public:
  
  //default constructor's 
  GMRES_RightPrecon_MatrixFree(void):
    Soln_ptr(NULL), List_of_Local_Solution_Blocks(NULL), Input_Parameters(NULL),
    G(NULL), Number_of_GMRES_Iterations(0), global_time_step_size(ZERO), 
    relative_residual(ZERO),DTS_ptr(NULL) {}
  
  GMRES_RightPrecon_MatrixFree(SOLN_BLOCK_TYPE *Soln_ptr, 
			       AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
			       INPUT_TYPE &Input_Parameters, const int &blocksize, 
			       DTS_NKS_Quad_Block<SOLN_BLOCK_TYPE,INPUT_TYPE> *DTS);

  // Proper functions not built yet.
  // GMRES_RightPrecon_MatrixFree(const &GMRES_RightPrecon_MatrixFree) {}  
  // GMRES_RightPrecon_MatrixFree operator= (const &GMRES_RightPrecon_MatrixFree) {}

  //Memory Handlers
  void allocate() { G = new GMRES_Block<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>[List_of_Local_Solution_Blocks->Nblk]; } 
  void deallocate() {  if(G != NULL) delete[] G;  G = NULL; }
  ~GMRES_RightPrecon_MatrixFree() { deallocate(); }

  // Constructors
  int solve(Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> *Block_precon, double tol,
			bool *GMRES_restarted_at_least_once, bool *GMRES_failed, int *GMRES_iters,
			double *res_cputime, int *res_nevals);

  double deltaU(const int Bcount, const int i, const int j, const int k){
    return (G[Bcount].deltaU(i,j,k)); }
  double deltaU_test(const int Bcount, const int i, const int j, const int k){
    return (G[Bcount].deltaU_test(i,j,k)); }
  double b_test(const int Bcount, const int i, const int j, const int k){
    return (G[Bcount].b_test(i,j,k)); }  

  int Output_GMRES_vars_Tecplot(const int Number_of_Time_Steps,
			        const double &l2_norm,
			        const double &l2_norm_rel);

  // calulate perturbation parameter
  void calculate_epsilon_restart(double &epsilon); 
  void calculate_epsilon(double &epsilon, double &l2_norm_z, 
			 const int &search_direction_counter); 

  // norm calculations
  double L1_Norm_z(const int &search_direction_counter);
  double L2_Norm_z(const int &search_direction_counter);
  double L1_Norm_x(void);
  double L2_Norm_x(void);
};

/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/
template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
GMRES_RightPrecon_MatrixFree(SOLN_BLOCK_TYPE *SolnBlk, 
			     AdaptiveBlock2D_List &List_of_Local_Blocks,
			     INPUT_TYPE &IP, const int &blocksize, 
			     DTS_NKS_Quad_Block<SOLN_BLOCK_TYPE,INPUT_TYPE> *DTS) :
  
  Soln_ptr(SolnBlk), 
  List_of_Local_Solution_Blocks(&List_of_Local_Blocks),   
  Input_Parameters(&IP),
  DTS_ptr(DTS),
  Number_of_GMRES_Iterations(0),
  relative_residual(0)  
{
  //Setup Memory for GMRES_Block's 
  allocate(); 
  // Setup GMRES for each Block 
  for (int  Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk; ++Bcount ) {
    if (List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      G[Bcount].allocate(Input_Parameters->NKS_IP.GMRES_Restart,
			 Input_Parameters->NKS_IP.GMRES_Overlap,
			 Input_Parameters->NKS_IP.Normalization,
			 Soln_ptr[Bcount],IP,blocksize, DTS_ptr[Bcount]);
    } 
  } 
}

/********************************************************
 * Routine: GMRES Plane Rotation Functions              *
 ********************************************************/
template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
GeneratePlaneRotation(const double &dx, const double &dy, double &cs, double &sn) {

  if (dy == ZERO) {  //NOT A GOOD IDEA .....
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

template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
ApplyPlaneRotation(double &dx, double &dy, const double &cs,const double &sn) {
  double temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx =  temp; //cs * dx + sn * dy;
}

/****************************************************************************
 * GMRES_RightPrecon_MatrixFree::solve -- apply right-preconditioned        *
 *                                              matrix-free GMRES.          *
 ****************************************************************************/
template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
int GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
solve(Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> *Block_precon, 
      double tol,
      bool *GMRES_restarted_at_least_once, 
      bool *GMRES_failed, 
      int *GMRES_iters,
      double *res_cputime, 
      int *res_nevals) {

  int m1 = Input_Parameters->NKS_IP.GMRES_Restart+1; //restart+1 used in H calculation
  int error_flag(0);  
  double resid0(ZERO);
  double beta(ZERO);
  double epsilon(ZERO);
  double total_norm_z(ZERO);

  int met_tol(0);
  bool do_one_more_iter_for_check(false);

  //FORTRAN NAMES
  integer inc(1);  // vector stride is always 1
  doublereal temp;
  Number_of_GMRES_Iterations = 0;

  clock_t t0;

  // Once GMRES_restarted_at_least_once is set to true, it remains true.
  // That is, it is for information only (for the caller) and so should 
  // not be tested anywhere in a logical statement.
  *GMRES_restarted_at_least_once = false;
  *GMRES_failed = false;

  /**************************************************************************/
  //Setup/Reset GMRES variables
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {            
    if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      G[Bcount].Initialize();
    }
  }
  /**************************************************************************/
  
  /**************************************************************************/
  /******************* BEGIN GMRES CALCULATION ******************************/
  /**************************************************************************/
  bool i_first_time_through = true;

  do {

    /**************************************************************************/
    /************* FIRST TIME THROUGH - NO RESTART NEEDED *********************/
    /**************************************************************************/    
    // V(0) already set to ZERO when initialized, so nothing to do if not restart. 
      
    /**************************************************************************/
    /********** NOT FIRST TIME THROUGH - RESTART APPLIED **********************/
    /**************************************************************************/
    if(!i_first_time_through) {

      if (CFFC_Primary_MPI_Processor() && !do_one_more_iter_for_check){   
	 switch (Input_Parameters->NKS_IP.output_format) {
	   case OF_SCOTT:
	     cout << "\n GMRES Restarted at -> GMRES (Inner Iterations) = " ;
	     cout << Number_of_GMRES_Iterations; 	 
	     break;
	   case OF_ALISTAIR: break; // picked up below
	     default: break;
	 }
      } 
      
      /************************************************************************/
      /* Set vector switch for message passing or x (solution) instead of W. */
      for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){  
	  G[Bcount].vector_switch = 0;
	}
      }
      /************************************************************************/

      /************************************************************************/
      /* MPI barrier to ensure processor synchronization. */
      CFFC_Barrier_MPI();  
      
      /* Send "x" solution information between neighbouring blocks.*/          
      error_flag = Send_All_Messages(G, 
				     *List_of_Local_Solution_Blocks,
				     G[0].NumVar(), 
				     OFF);
      if (error_flag) {
	cout << "\n GMRES2D ERROR: GMRES message passing error on processor "
	     << List_of_Local_Solution_Blocks->ThisCPU
	     << ".\n";
	cout.flush();
      } 
      error_flag = CFFC_OR_MPI(error_flag);
      /************************************************************************/

      /**************************************************************************/
      /* Calculate global epsilon based on 2-norm of x. */  
      calculate_epsilon_restart(epsilon);

      /**************************************************************************/
      /******* BEGIN MATRIX-FREE FOR RESTART ************************************/
      /**************************************************************************/
      for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  
	  //Calculate R(U+epsilon*(Minv*x(i))) -> Soln_ptr.U =  Soln_ptr.Uo + epsilon * W(i)
	  G[Bcount].calculate_perturbed_residual_Restart(epsilon);	  
	  BCs(Soln_ptr[Bcount],*Input_Parameters);
		
	  t0 = clock();
	  
	  error_flag = dUdt_Residual_Evaluation_NKS<SOLN_BLOCK_TYPE,INPUT_TYPE>
	               (Soln_ptr[Bcount],DTS_ptr[Bcount],*Input_Parameters);
	  
	  *res_cputime += double(clock() - t0) / double(CLOCKS_PER_SEC); 
	} 
      } 
      /**************************************************************************/

      t0 = clock();

      /**************************************************************************/
      // Send boundary flux corrections at block interfaces with resolution changes.
      error_flag = Send_Conservative_Flux_Corrections(Soln_ptr,
						      *List_of_Local_Solution_Blocks,
						      Soln_ptr[0].NumVar());
      if (error_flag) {
	cout << "\n GMRES ERROR: Solution flux correction message passing error on processor "
	     << List_of_Local_Solution_Blocks->ThisCPU << ".\n";
	cout.flush();
      } 
      error_flag = CFFC_OR_MPI(error_flag);
      
      // Apply boundary flux corrections to residual to ensure that method is conservative.
      Apply_Boundary_Flux_Corrections(Soln_ptr, *List_of_Local_Solution_Blocks);
      /**************************************************************************/

			*res_cputime += double(clock() - t0) / double(CLOCKS_PER_SEC); 
      (*res_nevals)++;
      
      /////////////////////////// 2ND /////////////////////////////////// 
      if(Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == SECOND_ORDER) {
	for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
	  if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {	    
	    //Calculate R(U+epsilon*(Minv*x(i))) -> Soln_ptr.U =  Soln_ptr.Uo + epsilon * W(i)
	    G[Bcount].calculate_perturbed_residual_2nd_Restart(epsilon);	  
	    BCs(Soln_ptr[Bcount],*Input_Parameters);
	   
	    t0 = clock();
	    
	    error_flag = dUdt_Residual_Evaluation_NKS<SOLN_BLOCK_TYPE,INPUT_TYPE>
	      (Soln_ptr[Bcount],DTS_ptr[Bcount],*Input_Parameters);

	    *res_cputime += double(clock() - t0) / double(CLOCKS_PER_SEC); 
	  } 
	} 
	/**************************************************************************/
	
	/**************************************************************************/
	// Send boundary flux corrections at block interfaces with resolution changes.
	t0 = clock();

	error_flag = Send_Conservative_Flux_Corrections(Soln_ptr,
							*List_of_Local_Solution_Blocks,
							Soln_ptr[0].NumVar());
	if (error_flag) {
	  cout << "\n GMRES ERROR: Solution flux correction message passing error on processor "
	       << List_of_Local_Solution_Blocks->ThisCPU << ".\n";
	  cout.flush();
	} 
	error_flag = CFFC_OR_MPI(error_flag);
	
	// Apply boundary flux corrections to residual to ensure that method is conservative.
	Apply_Boundary_Flux_Corrections(Soln_ptr, *List_of_Local_Solution_Blocks);

			*res_cputime += double(clock() - t0) / double(CLOCKS_PER_SEC); 
      (*res_nevals)++;
	/**************************************************************************/
      }
      //////////////////////////////////////////////////////////////////////////
            
      /**************************************************************************/
      //Calculate Matrix Free V(0) = ( R(U+epsilon*x) - b) / epsilon - x / h */
      for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {	  
	  G[Bcount].calculate_Matrix_Free_Restart(epsilon);	  	   	  
	}
      }   
      /**************************************************************************/
      /****************** END OF MATRIX-FREE FOR RESTART ************************/
      /**************************************************************************/


    }  
    /**************************************************************************/
    /***********  END OF RESTART APPLIED  *************************************/
    /**************************************************************************/


    /**************************************************************************/
    // CALCULATE NORM OF FIRST SEARCH VECTOR, V(0).
    beta= ZERO;
    for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
      if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	/* V(0) = ONE*(-b) + V(0) = Ax - b */
	temp = ONE;
	F77NAME(daxpy)(&G[Bcount].scalar_dim, &temp, G[Bcount].b, &inc, 
		       &(G[Bcount].V[(0)*G[Bcount].scalar_dim]), &inc);

	/* beta = norm(V(0)) */	
	beta += sqr(G[Bcount].L2_Norm(&(G[Bcount].V[(0)*G[Bcount].scalar_dim])));
      }
    }

    beta = sqrt(CFFC_Summation_MPI(beta));
    /**************************************************************************/


    /**************************************************************************/
    //GMRES check compares computed maxtrix-vector products to real Ax by
    //using the "restart" Ax.
    if(Input_Parameters->NKS_IP.GMRES_CHECK && do_one_more_iter_for_check) {
      // is ||beta||/||b|| =  to final relative_residual from GMRES
      break;
    }
    /**************************************************************************/
    

    /**************************************************************************/
    // RESSCALE V(0) USING NORM.
    for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
      if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	/* V(0) = -V(0)/beta */
	temp = -ONE/beta;
	F77NAME(dscal)(&G[Bcount].scalar_dim, &temp, &(G[Bcount].V[(0)*G[Bcount].scalar_dim]), &inc);
	
	//reset residual's to zero, restart
	for (int i = 1; i <  G[Bcount].restart+1; i++)
	  G[Bcount].s[i] = ZERO;  
	
	G[Bcount].s[0] = beta;
	
	//Set vector switch to pass W 
	G[Bcount].vector_switch  = 1;  
      } 
    } 
    /**************************************************************************/
    
    /* save the very first residual norm */
    if (Number_of_GMRES_Iterations == 0) {
      resid0 = beta;
    }

    /*************************************************************************/
    /**************** Begin Primary GMRES Loop *******************************/
    /*************************************************************************/
    int search_direction_counter = -1;

    do {

      search_direction_counter++;
      Number_of_GMRES_Iterations++;

      /**************************************************************************/
      // Calculate z vector using preconditioner.
      for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  
	  /* Set search direction counter. */
	  G[Bcount].search_directions = search_direction_counter;  
	  
	  //Apply Block Precondtioner  z = Minv * V(i) -> stored in W(i). 
	  Block_precon[Bcount].Apply_Preconditioner(G[Bcount].scalar_dim, 1, 
                                                    &(G[Bcount].V[(search_direction_counter)*G[Bcount].scalar_dim]), 
						    G[Bcount].scalar_dim, 
                                                    &(G[Bcount].W[(search_direction_counter)*G[Bcount].scalar_dim]), 
						    G[Bcount].scalar_dim);	  
	} 
      } 
      /**************************************************************************/

    
      /**************************************************************************/
      /* MPI barrier to ensure processor synchronization. */
      CFFC_Barrier_MPI();  
      
      /* Send solution information between neighbouring blocks.*/
      /* Passes "W = Minv * V(i) = z"  information in ghost cells */       
      error_flag = Send_All_Messages(G, 
				     *List_of_Local_Solution_Blocks,
				     G[0].NumVar(), 
				     OFF); if (error_flag) {
	cout << "\n GMRES2D ERROR: GMRES message passing error on processor "
	     << List_of_Local_Solution_Blocks->ThisCPU
	     << ".\n";
	cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      /**************************************************************************/
  

      /**************************************************************************/
      /* Calculate global epsilon base on 2-norm of z. */
      calculate_epsilon(epsilon, total_norm_z, search_direction_counter);

     /**************************************************************************/

      /**************************************************************************/
      /***************** BEGIN MATRIX-FREE FOR PRIMARY GMRES LOOP ***************/
      /**************************************************************************/ 
           
      // Calculate perturbed Residual R(U+epsilon*(Minv*V(i)))
      for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) { 
	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) { 
 	  //Calculate R(U+epsilon*(Minv*V(i))) -> Soln_ptr.U =  Soln_ptr.Uo + epsilon * W(i)
 	  G[Bcount].calculate_perturbed_residual(epsilon);
	  //Apply Regular Soln_ptr BC'S 
	  BCs(Soln_ptr[Bcount],*Input_Parameters);

	  t0 = clock();

	  //modified to calculate in "overlap" cells as well
	  error_flag = dUdt_Residual_Evaluation_NKS<SOLN_BLOCK_TYPE,INPUT_TYPE>
	               (Soln_ptr[Bcount],DTS_ptr[Bcount],*Input_Parameters);

	  *res_cputime += double(clock() - t0) / double(CLOCKS_PER_SEC); 
	 } 
      }
      /**************************************************************************/

      t0 = clock();
    
      /**************************************************************************/
      // Send boundary flux corrections at block interfaces with resolution changes. (changes to dUdt)
      error_flag = Send_Conservative_Flux_Corrections(Soln_ptr,
	              			              *List_of_Local_Solution_Blocks,
						      Soln_ptr[0].NumVar());
      if (error_flag) {
	cout << "\n GMRES ERROR: Solution flux correction message passing error on processor "
	     << List_of_Local_Solution_Blocks->ThisCPU << ".\n"; cout.flush();
      } 
      error_flag = CFFC_OR_MPI(error_flag);
	  
      // Apply boundary flux corrections to residual to ensure that method is conservative.
      Apply_Boundary_Flux_Corrections(Soln_ptr,*List_of_Local_Solution_Blocks);

			*res_cputime += double(clock() - t0) / double(CLOCKS_PER_SEC); 
			(*res_nevals)++;
      
      /**************************************************************************/

      ////////////////// 2nd ORDER ////////////////////////////////////////
      if(Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == SECOND_ORDER) {
	 
	// Calculate perturbed Residual R(U+epsilon*(Minv*V(i)))
	for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) { 
	  if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) { 
	    //Calculate R(U+epsilon*(Minv*V(i))) -> Soln_ptr.U =  Soln_ptr.Uo - epsilon * W(i)
	    G[Bcount].calculate_perturbed_residual_2nd(epsilon);
	    //Apply Regular Soln_ptr BC'S 
	    BCs(Soln_ptr[Bcount],*Input_Parameters);
	    
	    t0 = clock();	 
	    
	    error_flag = dUdt_Residual_Evaluation_NKS<SOLN_BLOCK_TYPE,INPUT_TYPE>
	      (Soln_ptr[Bcount],DTS_ptr[Bcount],*Input_Parameters);
	    
	    *res_cputime += double(clock() - t0) / double(CLOCKS_PER_SEC); 
	  } 
	}
	/**************************************************************************/
    
	/**************************************************************************/
	// Send boundary flux corrections at block interfaces with resolution changes.
			t0 = clock();
	error_flag = Send_Conservative_Flux_Corrections(Soln_ptr,
							*List_of_Local_Solution_Blocks,
							Soln_ptr[0].NumVar());
	if (error_flag) {
	  cout << "\n GMRES ERROR: Solution flux correction message passing error on processor "
	       << List_of_Local_Solution_Blocks->ThisCPU << ".\n"; cout.flush();
	} 
	error_flag = CFFC_OR_MPI(error_flag);
	
	// Apply boundary flux corrections to residual to ensure that method is conservative.
	Apply_Boundary_Flux_Corrections(Soln_ptr,*List_of_Local_Solution_Blocks);
			*res_cputime += double(clock() - t0) / double(CLOCKS_PER_SEC); 
			(*res_nevals)++;
	
      }
      ////////////////////////////////////////////////////////////////////////////


      /**************************************************************************/
      // Calculate Matrix Free V(i+1) = ( R(U + epsilon*W) - R(U) ) / epsilon - z / h
      for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) { 
	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) { 
	  G[Bcount].calculate_Matrix_Free(epsilon);	  	  
	} 
      }

      /**************************************************************************/
      /******************* END OF MATRIX-FREE ***********************************/
      /**************************************************************************/

      /**************************************************************************/
      // H norm calculation
      double total_norm_H = ZERO;

      for (int k = 0; k <= search_direction_counter; k++) {
	total_norm_H = ZERO;
	for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
	  if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	    /* total_H = dot(W,V(k)) */
	    total_norm_H += G[Bcount].Dotproduct( &(G[Bcount].V[(search_direction_counter+1)*G[Bcount].scalar_dim]),
						  &(G[Bcount].V[(k)*G[Bcount].scalar_dim]) );
	  }
	}	
	total_norm_H = CFFC_Summation_MPI(total_norm_H);

	for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
	  if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	    
	    /* V(i+1) = V(i+1) -  H(k, i) * V(k) */
	    G[Bcount].H[(search_direction_counter)*m1+(k)] = total_norm_H;

	    temp = -ONE * G[Bcount].H[(search_direction_counter)*m1+(k)]; //store in temp for passing by address

	    F77NAME(daxpy)(&G[Bcount].scalar_dim, &temp, &(G[Bcount].V[(k)*G[Bcount].scalar_dim]), &inc,
			   &(G[Bcount].V[(search_direction_counter+1)*G[Bcount].scalar_dim]), &inc);
	  } 
	} 
      } 
      /**************************************************************************/

#ifdef _NKS_VERBOSE
      if (CFFC_Primary_MPI_Processor()) { 
	cout << "\n  GMRES (Inner Iterations) = " << Number_of_GMRES_Iterations << " total_norm_H1 " <<   total_norm_H;	 
      } 
#endif
 
      /**************************************************************************/
      total_norm_H = ZERO;
      /* Calculate 2-norm of V(i+1) -> H(i+1,1) */
      for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	  total_norm_H += sqr(G[Bcount].L2_Norm(&(G[Bcount].V[(search_direction_counter+1)*G[Bcount].scalar_dim]))); 	  
	} 
      } 
      	
      total_norm_H = sqrt(CFFC_Summation_MPI(total_norm_H)); 
      /**************************************************************************/

#ifdef _NKS_VERBOSE   
      if (CFFC_Primary_MPI_Processor()) { 
	cout << "\n  GMRES (Inner Iterations) = " << Number_of_GMRES_Iterations << " total_norm_H2 " <<   total_norm_H;	 
      } 
#endif

      /**************************************************************************/
      for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	  
	  G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter+1)] = total_norm_H;

	  temp = ONE / G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter+1)]; // total_norm_H;
	  
	  /* V(i+1) = V(i+1) / H(i+1, i) */
	  F77NAME(dscal)(&G[Bcount].scalar_dim, &temp, &(G[Bcount].V[(search_direction_counter+1)*G[Bcount].scalar_dim]), &inc);
	  
	  } 
      } 
      /**************************************************************************/

      
      /**************************************************************************/
      for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	  
	  //POSSIBLY USE BLAS1 drtog, drotmg, drot, etc to do these faster ???
	  
	  for (int k = 0; k < search_direction_counter; k++) {
	    ApplyPlaneRotation(G[Bcount].H[(search_direction_counter)*m1+(k)], 
			       G[Bcount].H[(search_direction_counter)*m1+(k+1)], 
			       G[Bcount].cs[k], G[Bcount].sn[k]);
	  } 
	  
	  GeneratePlaneRotation(G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter)],
				G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter+1)], 
				G[Bcount].cs[search_direction_counter], G[Bcount].sn[search_direction_counter]);
	  
	  
	  ApplyPlaneRotation(G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter)],
			     G[Bcount].H[(search_direction_counter)*m1+(search_direction_counter+1)], 
			     G[Bcount].cs[search_direction_counter], G[Bcount].sn[search_direction_counter]);
	  
	  
	  ApplyPlaneRotation(G[Bcount].s[search_direction_counter], G[Bcount].s[search_direction_counter+1], 
			     G[Bcount].cs[search_direction_counter], G[Bcount].sn[search_direction_counter]);
	  
	} 
      } 
      /**************************************************************************/


      /**************************************************************************/
      relative_residual = ZERO;
      for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
	if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	  relative_residual = max( fabs(G[Bcount].s[search_direction_counter+1]) , relative_residual); 
	}
      }      
      
      relative_residual = CFFC_Maximum_MPI(relative_residual)/resid0;    

      /**************************************************************************/
      
      // Output progress

      //Verbose Output for CHECKING
      if (CFFC_Primary_MPI_Processor() && Number_of_GMRES_Iterations%5 == 0 && 
          Input_Parameters->NKS_IP.GMRES_CHECK) { 
        if(Number_of_GMRES_Iterations == 5){    
          cout << "\n  GMRES Iter.  \t   resid0 \t   resid \t  rel_resid  \t   L2||z||   \t  epsilon ";
        } 
        cout << "\n \t" << Number_of_GMRES_Iterations << "\t" << resid0 << "\t" 
             << relative_residual*resid0 << "\t" << relative_residual
             <<"\t"<< total_norm_z << "\t  "<<epsilon;   
      } 

      if (relative_residual <= tol) { 
	 met_tol = 1; 
      }
      // Now that the tolerance is calculated using math
      // functions (see the call to GMRES::solve() in NKS2D.h), 
      // round-off error on different computers
      // could mean that one CPU meets the tolerance while
      // another does not. This MPI-OR statement avoids any
      // such problems.
      //   -- Alistair Wood Wed Feb 28 2007 
      met_tol = CFFC_OR_MPI(met_tol);

      if (met_tol) {
	break;
      } else if (Number_of_GMRES_Iterations >= 
	         Input_Parameters->NKS_IP.Maximum_Number_of_GMRES_Iterations) {
	*GMRES_failed = true;
	break;
      } else if (search_direction_counter >= Input_Parameters->NKS_IP.GMRES_Restart-1) {
	*GMRES_restarted_at_least_once = true;
	break;
      }

    } while (1);

    /**************************************************************************/
    /******************* END PRIMARY GMRES LOOP *******************************/
    /**************************************************************************/

    /************************************************************************/
    // UPDATE SOLUTION. 
    for (int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
      if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED){  
		
	/* solve upper triangular system in place */	  
	for (int j = search_direction_counter; j >= 0; j--) { 
	  G[Bcount].s[j] /= G[Bcount].H[(j)*m1+(j)];
	  for (int k = j-1; k >= 0; k--)
	    G[Bcount].s[k] -= G[Bcount].H[(j)*m1+(k)] * G[Bcount].s[j];
	}
	
	/* update solution, x */ 
	for (int j = 0; j <= search_direction_counter; j++) {
	  /* x = x + s[j] * W(j) */
	  F77NAME(daxpy)(&G[Bcount].scalar_dim, &G[Bcount].s[j], &(G[Bcount].W[(j)*G[Bcount].scalar_dim]), 
			 &inc, G[Bcount].x, &inc); 
	}   
      } 
    } 

    if (met_tol || *GMRES_failed) { // We have met the tolerance or failed so we will stop here ...
	if (Input_Parameters->NKS_IP.GMRES_CHECK) { // ... but if GMRES_CHECK ...
 	   do_one_more_iter_for_check = true; // ... then do one more iteration to check GMRES residual. 
	} else {
	   break;
	}
    }
    // and if we did not meet the tolerance and we have not failed then set:
    i_first_time_through = false;
    // and go back to the start for a restart.
  } while (1);
 
  /**************************************************************************/
  /********************* END GMRES CALCULATION ******************************/
  /**************************************************************************/

  /**************************************************************************/
  // Check GMRES for accurate matrix vector products 
  // ie. does ||(Ax-b)||/||b|| ~= relative_residual
  
  if(Input_Parameters->NKS_IP.GMRES_CHECK){

    int blocksize = G[0].blocksize;
    double total_norm_b(ZERO);
    double *total_norm_eqn_b = new double[blocksize];
    double *total_norm_eqn_r = new double[blocksize];
    for(int i=0; i< blocksize; i++){  total_norm_eqn_r[i] = ZERO; total_norm_eqn_b[i] =ZERO; }

    //L2 norms 
    for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
      if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	total_norm_b += sqr(G[Bcount].L2_Norm( G[Bcount].b));	
	for(int i=0; i< blocksize; i++){	
	  total_norm_eqn_b[i] += sqr(G[Bcount].L2_Norm(i, G[Bcount].b));	
	  total_norm_eqn_r[i] += sqr(G[Bcount].L2_Norm(i, G[Bcount].V));		  
	}
      }
    } 
     
    total_norm_b = sqrt(CFFC_Summation_MPI(total_norm_b));
    for(int i=0; i< blocksize; i++){
      total_norm_eqn_r[i] = sqrt(CFFC_Summation_MPI(total_norm_eqn_r[i]));
      total_norm_eqn_b[i] = sqrt(CFFC_Summation_MPI(total_norm_eqn_b[i]));
    }

    //Output Reduction for systems & per equation
    if (CFFC_Primary_MPI_Processor()) {

      cout<<"\n GMRES Check total   ||(Ax-b)||/||b|| = "<<beta<<" / "<<total_norm_b<<" = "<<beta/total_norm_b
	  <<" /relative_residual (ideal 1) -> " <<((beta/total_norm_b)/relative_residual)<<endl;

       for(int i=0; i< blocksize; i++){	
	 if (total_norm_eqn_b[i]) // -> watch divide by zero
	   cout<<" Reduction for Eqn "<<i<<" ||(Ax-b)||/||b|| = "<<total_norm_eqn_r[i]
	       <<" / "<<total_norm_eqn_b[i]<<" = "<<total_norm_eqn_r[i]/total_norm_eqn_b[i]<<endl;
	 else
	   cout<<" Reduction for Eqn "<<i<<" ||(Ax-b)||/||b|| = "<<total_norm_eqn_r[i]
	       <<" / "<<total_norm_eqn_b[i]<<" = "<<total_norm_eqn_r[i]/PICO<<endl;	   
       }
    }   

    delete[] total_norm_eqn_b; delete[] total_norm_eqn_r; 
  }

  if (CFFC_Primary_MPI_Processor() && !Input_Parameters->NKS_IP.Dual_Time_Stepping ) { //
    switch (Input_Parameters->NKS_IP.output_format) {
    case OF_SCOTT:
      cout << "\n Finished GMRES with (Inner Iterations) = " << Number_of_GMRES_Iterations;
      cout << " resid0 = " << resid0;
      cout << " resid = " << relative_residual*resid0;
      cout << " relative_residual = " << relative_residual << endl;
      if (*GMRES_failed) {
	cout << "\n GMRES2D - GMRES ERROR: Unable to reach the specified convergence tolerance.";
	cout << "\n Final gmrestol = " << relative_residual << endl;
      }
      break;
    case OF_ALISTAIR: {
      int output_width = Input_Parameters->NKS_IP.output_width; 
      cout << setw(3) << Number_of_GMRES_Iterations << "  ";
      cout.unsetf(ios::scientific); cout.setf(ios::fixed);
      cout << setw(output_width) << relative_residual * 1000.0;
      cout.unsetf(ios::fixed); cout.setf(ios::scientific);
      cout << setw(output_width) << resid0;
      if (*GMRES_restarted_at_least_once) { 
	cout << "  R"; 
      } else { 
             cout << "   "; 
      }
      if (*GMRES_failed) { 
	cout << "  F"; 
      } else { 
	cout << "   "; 
      }
    }
      break;
    default:
      break;
    } 
  } 
  
  *GMRES_iters = Number_of_GMRES_Iterations;
  return error_flag;

} /* End of GMRES_RightPrecon_MatrixFree::solve. */ 

template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
int GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Output_GMRES_vars_Tecplot(const int Number_of_Time_Steps,
 		          const double &l2_norm,
 		          const double &l2_norm_rel) {

  if (List_of_Local_Solution_Blocks->Nused() == 0) return 0;
  int i = 0;
  bool print_title = true;
  char prefix[256], output_file_name[256];
  ofstream output_file;    

  for (i = 0; i < strlen(Input_Parameters->Output_File_Name); i++) {
     if (Input_Parameters->Output_File_Name[i] == ' ' ||
	Input_Parameters->Output_File_Name[i] == '.') {
	break;
     }
  }
  strncpy(prefix, Input_Parameters->Output_File_Name, i);
  prefix[i] = '\0';
	
  sprintf(output_file_name, "%s_n1%.4d_gmres_cpu%.3d.dat",
          prefix, Number_of_Time_Steps, List_of_Local_Solution_Blocks->ThisCPU);
  
  output_file.open(output_file_name, ios::out);
  if (!output_file.good()) { return 1; }

  for (int blk = 0; blk < List_of_Local_Solution_Blocks->Nblk; blk++) {
     if (List_of_Local_Solution_Blocks->Block[blk].used == ADAPTIVEBLOCK2D_USED) {
        G[blk].Output_GMRES_vars_Tecplot(Number_of_Time_Steps,
		 	                 blk,
					 print_title,
					 l2_norm,
					 l2_norm_rel,
					 output_file);
        if (print_title) { print_title = false; }
     }
  }

  output_file.close();
  return 0;

}

/**************************************************************************
 * Routine: calculate_epsilon                                             *
 **************************************************************************/
//
// Restart vesion
//
template <typename SOLN_VAR_TYPE,
	  typename SOLN_BLOCK_TYPE, 
	  typename INPUT_TYPE>
inline void GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,
					 SOLN_BLOCK_TYPE,
					 INPUT_TYPE>::
calculate_epsilon_restart(double &epsilon)
{    
  // Calculate global epsilon based on 2-norm of x.
  // FIXED - extra square root of L2_Norm_x here, should be eps = sqrt(eps_mach)/||x||_2
  epsilon = Input_Parameters->NKS_IP.Epsilon_Naught/L2_Norm_x();
}


//
// Non-Restart version
//
template <typename SOLN_VAR_TYPE,
	  typename SOLN_BLOCK_TYPE, 
	  typename INPUT_TYPE>
inline void GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,
					 SOLN_BLOCK_TYPE,
					 INPUT_TYPE>::
calculate_epsilon(double &epsilon, double &l2_norm_z, const int &search_direction_counter)
{    
  // compute l2 norm of z
  l2_norm_z = L2_Norm_z(search_direction_counter);

  // Calculate global epsilon base on 2-norm of z.
  // FIXED - extra square root of L2_Norm_z here, should be eps = sqrt(eps_mach)/||z||_2
  epsilon = Input_Parameters->NKS_IP.Epsilon_Naught/l2_norm_z;
}


/**************************************************************************
 * Some necessary norm calculations.                                      *
 **************************************************************************/
//
// 2-norm of z.
//
template <typename SOLN_VAR_TYPE,
	  typename SOLN_BLOCK_TYPE, 
	  typename INPUT_TYPE>
inline double GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,
					   SOLN_BLOCK_TYPE,
					   INPUT_TYPE>::
L2_Norm_z(const int &search_direction_counter)
{    
  double total_norm_z(ZERO); 
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
    if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      total_norm_z += sqr(G[Bcount].L2_Norm(&(G[Bcount].W[(search_direction_counter)*
							  G[Bcount].scalar_dim])));
    } 
  }       
  return sqrt(CFFC_Summation_MPI(total_norm_z));
}

//
// 1-norm of z.
//
template <typename SOLN_VAR_TYPE,
	  typename SOLN_BLOCK_TYPE, 
	  typename INPUT_TYPE>
inline double GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,
					   SOLN_BLOCK_TYPE,
					   INPUT_TYPE>::
L1_Norm_z(const int &search_direction_counter)
{    
  double total_norm_z(ZERO); 
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
    if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      total_norm_z += fabs(G[Bcount].L2_Norm(&(G[Bcount].W[(search_direction_counter)*
							   G[Bcount].scalar_dim])));
    } 
  }       
  return CFFC_Summation_MPI(total_norm_z);
}


//
// 2-norm of x.
//
template <typename SOLN_VAR_TYPE,
	  typename SOLN_BLOCK_TYPE, 
	  typename INPUT_TYPE>
inline double GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,
					   SOLN_BLOCK_TYPE,
					   INPUT_TYPE>::
L2_Norm_x(void)
{
  double total_norm_x(ZERO);
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
    if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      total_norm_x += sqr( G[Bcount].L2_Norm(G[Bcount].x) );	  
    }
  } 
  return sqrt(CFFC_Summation_MPI(total_norm_x));      
}

//
// 1-norm of x.
//
template <typename SOLN_VAR_TYPE,
	  typename SOLN_BLOCK_TYPE, 
	  typename INPUT_TYPE>
inline double GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,
					   SOLN_BLOCK_TYPE,
					   INPUT_TYPE>::
L1_Norm_x(void)
{
  double total_norm_x(ZERO);
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
    if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      total_norm_x += fabs( G[Bcount].L2_Norm(G[Bcount].x) );	  
    }
  } 
  return CFFC_Summation_MPI(total_norm_x);      
}


#endif // _GMRES2D_INCLUDED

