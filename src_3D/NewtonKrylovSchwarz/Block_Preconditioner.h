#ifndef _BLOCK_PRECONDITONER_INCLUDED 
#define _BLOCK_PRECONDITONER_INCLUDED 

/* BPKIT data structures and precondtioners along with BLAS 
   FORTRAN wrappers */ 

#ifndef _BLOCKMAT_H_
#include "BlockMat.h"
#endif //_BLOCKMAT_H_

#ifndef _DENSEMAT_H_
#include "DenseMat.h"
#endif //_DENSEMAT_H_

#ifndef _BILUK_H_
#include "BILUK.h"
#endif //_BILUK_H_

#ifndef _BRELAX_H_  //For BJacobi
#include "BRelax.h"
#endif //_BRELAX_H_

#include "../Math/Matrix.h"
#include "../Math/Vector2D.h"

#include "NKS_DTS.h"

enum ThreeDStencils { FIRST_ORDER_STENCIL = 7,
		  SECOND_ORDER_STENCIL = 27 }; 

enum Locations {  STENCIL_CENTER = 0,		  
                  STENCIL_NORTH = 1,
		  STENCIL_SOUTH = 2,
		  STENCIL_EAST  = 3,
		  STENCIL_WEST  = 4,
		  STENCIL_BOTTOM = 5,
		  STENCIL_TOP   = 6 };

/********************************************************
 * Class: Block_Preconditioner                          *
 *                                                      *
 * Member functions                                     * 
 *                                                      *
 ********************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
class Block_Preconditioner {
  private:
  int blocksize;
  int Jacobian_stencil_size; 

  void Get_Block_Index(const int &i, int *);
  void Get_Block_Index(const int &,const int &,const int &, int *, int *);
  void Setup_Jacobian_approximation();          //called from Create_Preconditioner
  void Setup_Preconditioner();                  //called from Update_Jacobian
  void Implicit_Euler(const int&,const int&,const int&, DenseMatrix*, const double&);      //called from Update_Jacobian       
  void First_Order_Inviscid_Jacobian_HLLE(const int&,const int&,const int&, DenseMatrix*); //called from Update_Jacobian
  void First_Order_Inviscid_Jacobian_Roe(const int&,const int&,const int&, DenseMatrix*);  //called from Update_Jacobian
  void First_Order_Inviscid_Jacobian_AUSM_plus_up(const int&,const int&,const int&, DenseMatrix*);  //called from Update_Jacobian
  void Second_Order_Viscous_Jacobian(const int&,const int&,const int&, DenseMatrix*);      //called from Update_Jacobian

  DenseMatrix Rotation_Matrix_3D(const Vector3D &nface, const int &A_matrix); //Used in Jacobian_LocalBlock
  DenseMatrix Rotation_Matrix_2D(const Vector2D &nface, const int &A_matrix); //Used in Jacobian_LocalBlock
  DenseMat DenseMatrix_to_DenseMat(const DenseMatrix &B); 

  protected:
  public:
  //Address of corresponding Solution block data
  Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>           *SolnBlk;                                                    
  Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>     *Input;

  // BPKIT BlockMat storage of Approximate Jacobian for use with Preconditioner  
  BlockMat Block_Jacobian_approx;

  // BPKIT BpPrecon Preconditioner Types
  BILUK   *ILUK_Precon;
  BJacobi *Jacobi_Precon;

  /*******************************************************/
  //default constructors
  Block_Preconditioner(void):
    ILUK_Precon(NULL), Jacobi_Precon(NULL), SolnBlk(NULL), Input(NULL) {}
  Block_Preconditioner( Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> *SolnBlk_ptr, 
			Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &Input_ptr, 			
			const int &_blocksize )
  { Create_Preconditioner(Input_ptr,SolnBlk_ptr, _blocksize); }
  
  //Constructors
  //These should be generalized for any precondtioner, not just "BPkit" types....
  void Create_Preconditioner( Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> &SolnBlk_ptr,
			      Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &Input_ptr, 			      
			      const int &_blocksize);
  void Update_Jacobian_and_Preconditioner(const double &DTS_dTime);  
  void Apply_Preconditioner(int A, int B, const double *C, int D, double *E, int F); 
  
  //Equation Specific Specializations required for Inviscid, Viscous, and Source Term Jacobians
  //Inviscid
  void Preconditioner_dFIdU(DenseMatrix &dFdU,SOLN_pSTATE W); 
  void Preconditioner_dFIdU_Roe(DenseMatrix &, const int, const int, const int,const int);  
  void Preconditioner_dFIdU_AUSM_plus_up(DenseMatrix &, const int, const int, const int);  
  //Viscous
  void Preconditioner_dFVdU(DenseMatrix &, const int, const int, const int, const int, const int ,const int);
  //Source 
  void Preconditioner_dSdU(int,int,int,DenseMatrix&);
  void normalize_Preconditioner_dFdU(DenseMatrix &dFdU);

  //Destructor / memory cleanup
  ~Block_Preconditioner(){ deallocate(); deallocate_Precon(); }
  void deallocate_Precon(){
    if (ILUK_Precon != NULL) delete ILUK_Precon; ILUK_Precon=NULL;
    if (Jacobi_Precon != NULL) delete Jacobi_Precon; Jacobi_Precon=NULL;
  }
  void deallocate(){}

};

/*!***********************************************************
 *  dFdU needs to be Provided By _Quad_Block Specialization  *
 *  otherwise the following errors will be shown.            *
 *************************************************************/ 
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Preconditioner_dFIdU(DenseMatrix &dFdU, SOLN_pSTATE W) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF Preconditioner_dFIdU for Block_Preconditioner2D.h requried \n";
  exit(1);
}

template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Preconditioner_dFIdU_Roe(DenseMatrix &dFdU, const int, const int, const int, const int) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF Preconditioner_dFIdU_Roe for Block_Preconditioner2D.h requried \n";
  exit(1);
}

template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Preconditioner_dFIdU_AUSM_plus_up(DenseMatrix &dFdU, const int, const int, const int) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF Preconditioner_dFIdU_AUSM_plus_up for Block_Preconditioner2D.h requried \n";
  exit(1);
}

template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Preconditioner_dFVdU(DenseMatrix &dFvdU, const int Rii, const int Rjj, 
		     const int Wii, const int Wjj, const int Orient_face, const int Orient_cell) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF Preconditioner_dFVdU for Block_Preconditioner2D.h requried \n";
  exit(1);
}

// BLANK UNLESS OTHERWISE SPECIALIZED
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Preconditioner_dSdU( int cell_index_i, int cell_index_j, int cell_index_k, DenseMatrix &Jacobian) {
}

template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF normalize_Preconditioner_dFdU for Block_Preconditioner2D.h requried \n";
  exit(1);
}


/*!**************************************************************
 * Generate block matrix indexing for a given matrix row "i"    *
 ****************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Get_Block_Index(const int &i, int *Block_index_j)
{
  //Determine Block_Indicies for Block Matrix formation
  if( Jacobian_stencil_size == FIRST_ORDER_STENCIL){        

    Block_index_j[STENCIL_CENTER] = i;   
    if ( i > SolnBlk->NCi) {  
      Block_index_j[STENCIL_NORTH] = (i%(SolnBlk->NCi*SolnBlk->NCj) >= SolnBlk->NCi) ? i - SolnBlk->NCi : -1;  
      Block_index_j[STENCIL_SOUTH] = ( (i+SolnBlk->NCi)%(SolnBlk->NCi*SolnBlk->NCj) >= SolnBlk->NCi) ? i + SolnBlk->NCi : -1;   
    } else {
      Block_index_j[STENCIL_NORTH] = i - SolnBlk->NCi;
      Block_index_j[STENCIL_SOUTH] = i + SolnBlk->NCi;
    }
    Block_index_j[STENCIL_EAST] = ( i%SolnBlk->NCi != 0) ? i-1 : -1;
    Block_index_j[STENCIL_WEST] = ( (i+1)%SolnBlk->NCi != 0) ? i+1 : -1;
    Block_index_j[STENCIL_BOTTOM] = i - SolnBlk->NCi* SolnBlk->NCj;   
    Block_index_j[STENCIL_TOP] =  SolnBlk->NCi* SolnBlk->NCj + i ;
        
  }  else if ( Jacobian_stencil_size == SECOND_ORDER_STENCIL) {
//     Block_index_j[NE] = ( (i - SolnBlk->NCi)%SolnBlk->NCi != 0) ? i - SolnBlk->NCi - 1 : -1;
//     Block_index_j[NORTH] = i - SolnBlk->NCi;
//     Block_index_j[NW] = ( (i - SolnBlk->NCi+1)%SolnBlk->NCi != 0) ? i - SolnBlk->NCi + 1 : -1;
//     Block_index_j[EAST] = ( i%SolnBlk->NCi != 0) ? i-1 : -1;
//     Block_index_j[CENTER] = i;
//     Block_index_j[WEST] = ( (i+1)%SolnBlk->NCi != 0) ? i+1 : -1;
//     Block_index_j[SE] = ( (i + SolnBlk->NCi)%SolnBlk->NCi != 0) ? i + SolnBlk->NCi - 1 : -1;
//     Block_index_j[SOUTH] = i + SolnBlk->NCi;
//     Block_index_j[SW] = ( (i + SolnBlk->NCi+1)%SolnBlk->NCi != 0) ? i + SolnBlk->NCi + 1 : -1;
  }  
}

/****************************************************************
 * Generate block matrix indexing for a given cell (i,j)        *
 ****************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Get_Block_Index(const int &cell_index_i, const int &cell_index_j, const int &cell_index_k,
		int *Block_index_i, int *Block_index_j)
{ 
  //J index 
  for( int i=0; i<Jacobian_stencil_size; i++)  
    Block_index_j[i] = cell_index_k*(SolnBlk->NCi*SolnBlk->NCj) + 
                       cell_index_j*SolnBlk->NCi + cell_index_i; 

  //I index 
  if( Jacobian_stencil_size == FIRST_ORDER_STENCIL){        
    /*! Determine 1st order Block_Indicies for Cell (i,j)   
     *                    
     *            ---    5
     *           | 1 | / 
     *        --- --- ---
     *       | 4 | 0 | 3 |                
     *        --- --- ---
     *         / | 2 |
     *       6    ---
     */  
    Block_index_i[STENCIL_CENTER] = cell_index_k*(SolnBlk->NCi*SolnBlk->NCj) +
                                    cell_index_j*SolnBlk->NCi+cell_index_i;      
    Block_index_i[STENCIL_NORTH]  = Block_index_i[STENCIL_CENTER] - SolnBlk->NCi;
    Block_index_i[STENCIL_SOUTH]  = Block_index_i[STENCIL_CENTER] + SolnBlk->NCi;   
    Block_index_i[STENCIL_EAST]   = Block_index_i[STENCIL_CENTER] - 1;
    Block_index_i[STENCIL_WEST]   = Block_index_i[STENCIL_CENTER] + 1;   
    Block_index_i[STENCIL_BOTTOM] = Block_index_i[STENCIL_CENTER] - (SolnBlk->NCi*SolnBlk->NCj);
    Block_index_i[STENCIL_TOP]    = Block_index_i[STENCIL_CENTER] + (SolnBlk->NCi*SolnBlk->NCj);
    
  }  else if ( Jacobian_stencil_size == SECOND_ORDER_STENCIL) {
//     /*! Determine 2nd order Block_Indicies for Cell (i,j)   
//      *
//      *        --- --- ---
//      *       | 5 | 1 | 8 |
//      *        --- --- --- 
//      *       | 2 | 0 | 4 |               //2D -> Needs to be replace by 3D (27)
//      *       --- --- --- -
//      *       | 6 | 3 | 7 |
//      *        --- --- ---
//      */ 
//     Block_index_i[NE] = (cell_index_j - 1)*SolnBlk->NCi+cell_index_i-1;    //NE  
//     Block_index_i[NORTH] = (cell_index_j - 1)*SolnBlk->NCi+cell_index_i;   //NORTH
//     Block_index_i[NW] = (cell_index_j - 1)*SolnBlk->NCi+cell_index_i+1;    //NW   
//     Block_index_i[EAST] = cell_index_j*SolnBlk->NCi+cell_index_i-1;        //EAST      
//     Block_index_i[CENTER] = cell_index_j*SolnBlk->NCi+cell_index_i;        //CENTER
//     Block_index_i[WEST] = cell_index_j*SolnBlk->NCi+cell_index_i+1;        //WEST
//     Block_index_i[SE] = (cell_index_j + 1)*SolnBlk->NCi+cell_index_i-1;    //SE     
//     Block_index_i[SOUTH] = (cell_index_j + 1)*SolnBlk->NCi+cell_index_i;   //SOUTH    
//     Block_index_i[SW] = (cell_index_j + 1)*SolnBlk->NCi+cell_index_i+1;    //SW  
  }  

}


/*********************************************************
 *  Create the Preconditoner as required by BPKIT        *
 *********************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Create_Preconditioner(Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> &SolnBlk_ptr,
		      Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &Input_ptr, 
		      const int &_blocksize)
{
  SolnBlk = &(SolnBlk_ptr);
  Input = &(Input_ptr);
  blocksize = _blocksize;       //number of equations

  if (Input->NKS_IP.Jacobian_Order == SOURCE_TERMS_ONLY ||
      Input->NKS_IP.Jacobian_Order == FIRST_ORDER_INVISCID_HLLE ||
      Input->NKS_IP.Jacobian_Order == FIRST_ORDER_INVISCID_ROE || 
      Input->NKS_IP.Jacobian_Order == FIRST_ORDER_INVISCID_AUSMPLUSUP) {
    Jacobian_stencil_size = FIRST_ORDER_STENCIL; 
  } else if (Input->NKS_IP.Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_HLLE ||
             Input->NKS_IP.Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_ROE ||
             Input->NKS_IP.Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP) {
    Jacobian_stencil_size = SECOND_ORDER_STENCIL;  
    cerr<<"\n 2nd order Preconditioner Jacobian not finished yet! ";
    exit(1);
  } else {
    cerr<<"\n Invalid Jacobian Preconditioner Order "<<Input->NKS_IP.Jacobian_Order;
    exit(1);
  }
    
  //Sets up memory for Approximate Jacobian (Block_Jacobian_approx)
  Setup_Jacobian_approximation(); 
}

/********************************************************************************
 *  Jacobian (dR/dU) Approximation Stencil Allocation - only called on STARTUP   *
 *********************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Setup_Jacobian_approximation(){

  int block_mat_size = SolnBlk->NCi * SolnBlk->NCj * SolnBlk->NCk;   //block matrix size based on icells*jcells
  int nnz;

  //! 1st order 3D Jacobian, Sparse Block Matrix Memory Allocation of non-zero block entries
  if( Jacobian_stencil_size == FIRST_ORDER_STENCIL){        
    nnz = 2*( 2*( 8 + 5*(SolnBlk->NCi-2)) + (10 + 6*(SolnBlk->NCi-2))*(SolnBlk->NCj-2) ) 
       + (SolnBlk->NCk-2)*(  2*(10 + 6*(SolnBlk->NCi-2)) + (12 + 7*(SolnBlk->NCi-2))*(SolnBlk->NCj-2) );

    //! 2nd order 3D Jacobian
  }  else if ( Jacobian_stencil_size == SECOND_ORDER_STENCIL) {
    //nnz = 9*block_mat_size - 6*(SolnBlk->NCi + SolnBlk->NCj) + 4;    
  }      

  //cout<<"\n\n nnz = "<<nnz<<" for NCi= "<<SolnBlk->NCi<<" NCj= "<<SolnBlk->NCj<<" NCk= "<<SolnBlk->NCk;

  //!Temporary Storage Arrays
  int *i_index = new int[nnz];      //location of dense local blocks, in global sparse block Matrix   
  int *j_index = new int[nnz]; 
  double *Data = new double [nnz*blocksize*blocksize];
  int *block_j = new int[Jacobian_stencil_size];
  
  //!Setup Stencil for approximate Jacobian
  int nnz_count(0);   
  //! Loop through each row of Block Matrix
  for (int i=0; i< block_mat_size; i++){       
    Get_Block_Index(i,block_j);
    int stencil = 0; int local=0;
    //! Determine which entries correspond to 1st, 2nd, etc. order stencil    
    // cout<<"\n i "<<i<<" j(s) "; cout.flush();
    while( stencil < Jacobian_stencil_size) { 
      int j = block_j[stencil]; 
      //cout<<" "<<j; cout.flush();
      if( j >= 0 && j < block_mat_size){	
	i_index[nnz_count] = i+1;  // bpkit assumes Fortran indexing ie. starting at 1 
	j_index[nnz_count] = j+1;
	for(int k=0; k<blocksize*blocksize; k++){ //initialize Block Matrix as identity matrix
	  if( i == j && k%(blocksize+1) == 0){
	    Data[nnz_count*blocksize*blocksize + k] = ONE;
	  } else {
	    Data[nnz_count*blocksize*blocksize + k] = ZERO;
	  }
	}	 
	local++;
	nnz_count++;
      }
      stencil++; 
    } //cout<<"  tot = "<<local;     
  }

  if( nnz != nnz_count){ cerr<<"\n Number of nonzero blocks mismatch error in approximate Jacobian formation "
			     <<nnz<<" != "<<nnz_count<<"\n"; exit(1); }

  //! Create Sparse bpkit "BlockMat" Block Matrix 
  Block_Jacobian_approx.setup(block_mat_size,nnz,i_index,j_index,Data,blocksize);                  
  
  //! Clean up local memory
  delete[] i_index; delete[] j_index;  delete[] Data;  delete[] block_j;

}


/**********************************************************
 *  Update BlockMat with approximation to the Jacobian    *
 **********************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Update_Jacobian_and_Preconditioner(const double &DTS_dTime)
{
  
  //!Local Variables and Temporary Storage
  int block_mat_size = SolnBlk->NCi*SolnBlk->NCj*SolnBlk->NCk; 
  DenseMatrix *Jacobian_Data = new DenseMatrix[Jacobian_stencil_size];
  for(int i=0; i<Jacobian_stencil_size; i++) { Jacobian_Data[i] = DenseMatrix(blocksize,blocksize,ZERO); }
  int *block_i = new int[Jacobian_stencil_size]; 
  int *block_j = new int[Jacobian_stencil_size]; 

  //! Initially assume no overlap
  int JCl_overlap = 0, JCu_overlap = 0, ICu_overlap =0, ICl_overlap = 0, KCl_overlap = 0, KCu_overlap = 0;
  
//   //! If overlap determine which block boundaries are internal, ie. BC_NONE
//   if(Input->NKS_IP.GMRES_Overlap){	
//     if (SolnBlk->Grid.BCtypeS[SolnBlk->ICl] == BC_NONE)  JCl_overlap = Input->NKS_IP.GMRES_Overlap;
//     if (SolnBlk->Grid.BCtypeN[SolnBlk->ICu] == BC_NONE)  JCu_overlap = Input->NKS_IP.GMRES_Overlap;
//     if (SolnBlk->Grid.BCtypeE[SolnBlk->JCu] == BC_NONE)  ICu_overlap = Input->NKS_IP.GMRES_Overlap;
//     if (SolnBlk->Grid.BCtypeW[SolnBlk->JCl] == BC_NONE)  ICl_overlap = Input->NKS_IP.GMRES_Overlap;
//   }

  //*********************************************************************************//
  /*! Calculate Jacobians for each cell and Update Global Jacobian Block Matrix
   * loop through all non-ghost cells, including overlap cells.  Ghost cells already
   * set to Zero or Identity by initialization.
   **********************************************************************************/
  for(int i= SolnBlk->ICl - ICl_overlap; i<= SolnBlk->ICu + ICu_overlap; i++){    
    for(int j= SolnBlk->JCl - JCl_overlap; j<= SolnBlk->JCu + ICu_overlap; j++){  
      for(int k= SolnBlk->KCl - KCl_overlap; k<= SolnBlk->KCu + KCu_overlap; k++){  

	//--------------------------------------------------------------------------//
	//! Calculate Local Approximate Jacobian                        
	switch(Input->NKS_IP.Jacobian_Order){
	case SOURCE_TERMS_ONLY :
	  Implicit_Euler(i,j,k, Jacobian_Data,DTS_dTime);
	  Preconditioner_dSdU(i,j,k,Jacobian_Data[STENCIL_CENTER]);
	  break;
	case FIRST_ORDER_INVISCID_HLLE : 
	  Implicit_Euler(i,j,k, Jacobian_Data,DTS_dTime);
 	  First_Order_Inviscid_Jacobian_HLLE(i,j,k, Jacobian_Data);
	  Preconditioner_dSdU(i,j,k,Jacobian_Data[STENCIL_CENTER]);                       
	  break;
	case FIRST_ORDER_INVISCID_ROE : 
	  Implicit_Euler(i,j,k, Jacobian_Data,DTS_dTime);
	  First_Order_Inviscid_Jacobian_Roe(i,j,k, Jacobian_Data);
	  Preconditioner_dSdU(i,j,k,Jacobian_Data[STENCIL_CENTER]);
	  break;
	case FIRST_ORDER_INVISCID_AUSMPLUSUP : 
	  Implicit_Euler(i,j,k, Jacobian_Data,DTS_dTime);
	  First_Order_Inviscid_Jacobian_AUSM_plus_up(i,j,k, Jacobian_Data);
	  Preconditioner_dSdU(i,j,k,Jacobian_Data[STENCIL_CENTER]);
	  break;
	case SECOND_ORDER_DIAMOND_WITH_HLLE:
	  Implicit_Euler(i,j,k, Jacobian_Data,DTS_dTime);
	  First_Order_Inviscid_Jacobian_HLLE(i,j,k, Jacobian_Data);   
	  Second_Order_Viscous_Jacobian(i,j,k, Jacobian_Data);    
	  Preconditioner_dSdU(i,j,k,Jacobian_Data[STENCIL_CENTER]);   
	  break;
	case SECOND_ORDER_DIAMOND_WITH_ROE :	
	  Implicit_Euler(i,j,k, Jacobian_Data,DTS_dTime);
	  First_Order_Inviscid_Jacobian_Roe(i,j,k, Jacobian_Data);
	  Second_Order_Viscous_Jacobian(i,j,k, Jacobian_Data);    
	  Preconditioner_dSdU(i,j,k,Jacobian_Data[STENCIL_CENTER]); 
	  break;
	case SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP :	
	  Implicit_Euler(i,j,k, Jacobian_Data,DTS_dTime);
	  First_Order_Inviscid_Jacobian_AUSM_plus_up(i,j,k, Jacobian_Data);
	  Second_Order_Viscous_Jacobian(i,j,k, Jacobian_Data);    
	  Preconditioner_dSdU(i,j,k,Jacobian_Data[STENCIL_CENTER]);  
	  break;
	}
	
	//--------------------------------------------------------------------------//
	//! Get Block Matrix locations that have components from a given Cell(i,j)
	
	Get_Block_Index(i,j,k, block_i, block_j);
// 	cout<<"\n i "<<i<<" j "<<j<<" k "<<k<<" j "<<block_j[0]<<" i ";
// 	for( int block = 0; block < Jacobian_stencil_size; block++) cout<<" "<<block_i[block];	
	
	//HARD CODED FOR 2 GHOST CELLS & ZERO OVERLAP ???

	//fudge for iGhost Cell reset to zero   //jGhost cell already zero                
	if(block_i[STENCIL_NORTH] < k*(SolnBlk->NCi*SolnBlk->NCj) + TWO*SolnBlk->NCi)    { Jacobian_Data[STENCIL_NORTH].zero(); /* cout<<"\n NZERO";*/}
	if(block_i[STENCIL_SOUTH] > (k+1)*(SolnBlk->NCi*SolnBlk->NCj) - TWO*SolnBlk->NCi){ Jacobian_Data[STENCIL_SOUTH].zero();  /*cout<<"\n SZERO";*/}
	if(block_i[STENCIL_BOTTOM] < TWO*(SolnBlk->NCi*SolnBlk->NCj))                    { Jacobian_Data[STENCIL_BOTTOM].zero(); /*cout<<"\n BZERO";*/}
	if(block_i[STENCIL_TOP] > block_mat_size - TWO*(SolnBlk->NCi*SolnBlk->NCj))      { Jacobian_Data[STENCIL_TOP].zero();    /*cout<<"\n TZERO";*/}

//       if(Jacobian_stencil_size == SECOND_ORDER_STENCIL){		
// 	if(block_i[NE] < TWO*SolnBlk->NCi)    Jacobian_Data[NE].zero();
// 	if(block_i[NW] < TWO*SolnBlk->NCi)    Jacobian_Data[NW].zero(); 
// 	if(block_i[SE] > block_mat_size - TWO*SolnBlk->NCi)    Jacobian_Data[SE].zero();
// 	if(block_i[SW] > block_mat_size - TWO*SolnBlk->NCi)    Jacobian_Data[SW].zero();
//       }
// 	//--------------------------------------------------------------------------//
	//! Update BlockMat with Local Approximate Jacobians 
	for( int block = 0; block < Jacobian_stencil_size; block++){
	  // Normalize

// 	  cout<<endl<<Jacobian_Data[block];

	  normalize_Preconditioner_dFdU(Jacobian_Data[block]);
	  
	  //can be sped up by more intelligent logic in bkpkit (BlockMat.cc  "setblock")
	  Block_Jacobian_approx.setblock( block_i[block], block_j[block], DenseMatrix_to_DenseMat(Jacobian_Data[block]));
	  
	  Jacobian_Data[block].zero(); //Just in case to avoid +=/-= issues
	}     

      }
    }
  }
  
  // Local Memory cleanup
  delete[] Jacobian_Data; delete[] block_i; delete[] block_j;

  //Setup appropriate Preconditioner after Jacobian has been formed/Updated
  Setup_Preconditioner();
}

/*****************************************************************************
 *  Add finite time step to dia gonal.                                        *
 *****************************************************************************/ 
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Implicit_Euler(const int &cell_index_i,const int &cell_index_j,const int &cell_index_k,
	       DenseMatrix* Jacobian,const double& DTS_dTime){   

  DenseMatrix Diag(blocksize,blocksize);      
  Diag.identity();    
  //Cacluate LHS depeneding on Steady State of Dual Time Stepping
  Diag *= LHS_Time<SOLN_pSTATE,SOLN_cSTATE>(*Input, SolnBlk->dt[cell_index_i][cell_index_j][cell_index_k],DTS_dTime);  
  Jacobian[STENCIL_CENTER] -= Diag;

}

/*****************************************************************************
 *  Calculate First Order Local Jacobian Block(s) Coresponding to Cell(i,j)  *
 *  using HLLE                                                               *
 *****************************************************************************/ 
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
First_Order_Inviscid_Jacobian_HLLE(const int &cell_index_i,const int &cell_index_j,const int &cell_index_k, 
				   DenseMatrix* Jacobian){              
  
//   DenseMatrix test(blocksize,blocksize,cell_index_i+cell_index_j+cell_index_k); 
//   for( int block = 0; block < Jacobian_stencil_size; block++) Jacobian[block] += test;   

  //! Caculate normal vectors -> in Vector2D format. 
  Vector3D nface_N = SolnBlk->Grid.nfaceN(cell_index_i,cell_index_j-1,cell_index_k);
  Vector3D nface_S = SolnBlk->Grid.nfaceS(cell_index_i,cell_index_j+1,cell_index_k);      
  Vector3D nface_E = SolnBlk->Grid.nfaceE(cell_index_i-1,cell_index_j,cell_index_k);
  Vector3D nface_W = SolnBlk->Grid.nfaceW(cell_index_i+1,cell_index_j,cell_index_k);
  Vector3D nface_Bot = SolnBlk->Grid.nfaceBot(cell_index_i,cell_index_j,cell_index_k-1);
  Vector3D nface_Top = SolnBlk->Grid.nfaceTop(cell_index_i,cell_index_j,cell_index_k+1);

  //! Calculate wavespeeds using solutions in the rotated frame -> in Vector2D format.  
  Vector2D lambdas_N = SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].HLLE_wavespeeds(SolnBlk->W[cell_index_i][cell_index_j-1][cell_index_k], 
				       SolnBlk->W[cell_index_i][cell_index_j][cell_index_k], nface_N);
  Vector2D lambdas_S = SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].HLLE_wavespeeds(SolnBlk->W[cell_index_i][cell_index_j+1][cell_index_k], 
				       SolnBlk->W[cell_index_i][cell_index_j][cell_index_k], nface_S);  
  Vector2D lambdas_E = SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].HLLE_wavespeeds(SolnBlk->W[cell_index_i-1][cell_index_j][cell_index_k], 
				       SolnBlk->W[cell_index_i][cell_index_j][cell_index_k], nface_E);
  Vector2D lambdas_W = SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].HLLE_wavespeeds(SolnBlk->W[cell_index_i+1][cell_index_j][cell_index_k], 
				       SolnBlk->W[cell_index_i][cell_index_j][cell_index_k], nface_W);
  Vector2D lambdas_Bot = SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].HLLE_wavespeeds(SolnBlk->W[cell_index_i][cell_index_j][cell_index_k-1], 
					 SolnBlk->W[cell_index_i][cell_index_j][cell_index_k], nface_Bot);
  Vector2D lambdas_Top = SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].HLLE_wavespeeds(SolnBlk->W[cell_index_i][cell_index_j][cell_index_k+1], 
					 SolnBlk->W[cell_index_i][cell_index_j][cell_index_k], nface_Top);
 
  //! Calculate constants gamma and beta -> scalar values. 
  double gamma_N = (lambdas_N.x*lambdas_N.y)/(lambdas_N.y-lambdas_N.x);
  double beta_N  = - lambdas_N.x/(lambdas_N.y-lambdas_N.x);
  double gamma_S = (lambdas_S.x*lambdas_S.y)/(lambdas_S.y-lambdas_S.x);
  double beta_S  = - lambdas_S.x/(lambdas_S.y-lambdas_S.x);
  double gamma_E = (lambdas_E.x*lambdas_E.y)/(lambdas_E.y-lambdas_E.x);
  double beta_E  = - lambdas_E.x/(lambdas_E.y-lambdas_E.x);
  double gamma_W = (lambdas_W.x*lambdas_W.y)/(lambdas_W.y-lambdas_W.x);
  double beta_W  = - lambdas_W.x/(lambdas_W.y-lambdas_W.x);  
  double gamma_Bot = (lambdas_Bot.x*lambdas_Bot.y)/(lambdas_Bot.y-lambdas_Bot.x);
  double beta_Bot  = - lambdas_Bot.x/(lambdas_Bot.y-lambdas_Bot.x);
  double gamma_Top = (lambdas_Top.x*lambdas_Top.y)/(lambdas_Top.y-lambdas_Top.x);
  double beta_Top  = - lambdas_Top.x/(lambdas_Top.y-lambdas_Top.x);

  //! Obtain rotation matrices with normal vector -> matrices in DenseMatrix format. 
  DenseMatrix A_N( Rotation_Matrix_3D(nface_N, 1) );
  DenseMatrix AI_N( Rotation_Matrix_3D(nface_N, 0));
  DenseMatrix A_S( Rotation_Matrix_3D(nface_S, 1) );
  DenseMatrix AI_S( Rotation_Matrix_3D(nface_S, 0) );
  DenseMatrix A_E( Rotation_Matrix_3D(nface_E, 1) );
  DenseMatrix AI_E( Rotation_Matrix_3D(nface_E, 0) );
  DenseMatrix A_W( Rotation_Matrix_3D(nface_W, 1) );
  DenseMatrix AI_W( Rotation_Matrix_3D(nface_W, 0) );
  DenseMatrix A_Bot( Rotation_Matrix_3D(nface_Bot, 1) );
  DenseMatrix AI_Bot( Rotation_Matrix_3D(nface_Bot, 0) );
  DenseMatrix A_Top( Rotation_Matrix_3D(nface_Top, 1) );
  DenseMatrix AI_Top( Rotation_Matrix_3D(nface_Top, 0) );

  //! Calculate dFdU using solutions in the rotated frame -> matrix in DenseMatrix format. 
  DenseMatrix dFdU_N(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_S(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_E(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_W(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_Top(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_Bot(blocksize,blocksize,ZERO); 

  //Solution Rotate provided in pState 
  Preconditioner_dFIdU( dFdU_N, SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].Rotate(nface_N));   
  Preconditioner_dFIdU( dFdU_S, SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].Rotate(nface_S));
  Preconditioner_dFIdU( dFdU_E, SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].Rotate(nface_E));
  Preconditioner_dFIdU( dFdU_W, SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].Rotate(nface_W));
  Preconditioner_dFIdU( dFdU_Bot, SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].Rotate(nface_Bot));
  Preconditioner_dFIdU( dFdU_Top, SolnBlk->W[cell_index_i][cell_index_j][cell_index_k].Rotate(nface_Top));
  
  DenseMatrix II(blocksize,blocksize);  II.identity();     

  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  //North
  Jacobian[STENCIL_NORTH] = (SolnBlk->Grid.AfaceN(cell_index_i,cell_index_j-1,cell_index_k) 
		 * AI_N * (beta_N * dFdU_N + gamma_N * II) * A_N); 

  //South
  Jacobian[STENCIL_SOUTH] = (SolnBlk->Grid.AfaceS(cell_index_i,cell_index_j+1,cell_index_k) 
		 * AI_S * (beta_S * dFdU_S + gamma_S * II) * A_S);

  //East
  Jacobian[STENCIL_EAST] = (SolnBlk->Grid.AfaceE(cell_index_i-1,cell_index_j,cell_index_k) 
		 * AI_E * (beta_E * dFdU_E + gamma_E * II) * A_E);

  //West
  Jacobian[STENCIL_WEST] = (SolnBlk->Grid.AfaceW(cell_index_i+1,cell_index_j,cell_index_k) 
		 * AI_W * (beta_W * dFdU_W + gamma_W * II) * A_W);

  //Bottom
  Jacobian[STENCIL_BOTTOM] = (SolnBlk->Grid.AfaceBot(cell_index_i+1,cell_index_j,cell_index_k-1) 
			      * AI_Bot * (beta_Bot * dFdU_Bot + gamma_Bot * II) * A_Bot);

  //Top
  Jacobian[STENCIL_TOP] = (SolnBlk->Grid.AfaceTop(cell_index_i+1,cell_index_j,cell_index_k+1) 
			   * AI_Top * (beta_Top * dFdU_Top + gamma_Top * II) * A_Top);

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left) 
  Jacobian[STENCIL_CENTER] += (Jacobian[STENCIL_NORTH] + Jacobian[STENCIL_SOUTH] + Jacobian[STENCIL_EAST]  + Jacobian[STENCIL_WEST] +
			       Jacobian[STENCIL_BOTTOM] + Jacobian[STENCIL_TOP])/SolnBlk->Grid.volume(cell_index_i,cell_index_j,cell_index_k);

  Jacobian[STENCIL_NORTH] = -Jacobian[STENCIL_NORTH]/SolnBlk->Grid.volume(cell_index_i,cell_index_j-1,cell_index_k);
  Jacobian[STENCIL_SOUTH] = -Jacobian[STENCIL_SOUTH]/SolnBlk->Grid.volume(cell_index_i,cell_index_j+1,cell_index_k);
  Jacobian[STENCIL_EAST] = -Jacobian[STENCIL_EAST]/SolnBlk->Grid.volume(cell_index_i-1,cell_index_j,cell_index_k);
  Jacobian[STENCIL_WEST] = -Jacobian[STENCIL_WEST]/SolnBlk->Grid.volume(cell_index_i+1,cell_index_j,cell_index_k);
  Jacobian[STENCIL_BOTTOM] = -Jacobian[STENCIL_BOTTOM]/SolnBlk->Grid.volume(cell_index_i,cell_index_j,cell_index_k-1);
  Jacobian[STENCIL_TOP] = -Jacobian[STENCIL_TOP]/SolnBlk->Grid.volume(cell_index_i,cell_index_j,cell_index_k+1);

}

/*****************************************************************************
 *  Calculate First Order Local Jacobian Block(s) Coresponding to Cell(i,j)  *
 *  using Roe                                                                *
 *****************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
First_Order_Inviscid_Jacobian_Roe(const int &cell_index_i,const int &cell_index_j,const int &cell_index_k,  
				   DenseMatrix* Jacobian){              
    
//   //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
//   Preconditioner_dFIdU_Roe(Jacobian[STENCIL_NORTH],cell_index_i,cell_index_j,STENCIL_NORTH);
//   Preconditioner_dFIdU_Roe(Jacobian[STENCIL_SOUTH],cell_index_i,cell_index_j,STENCIL_SOUTH); 
//   Preconditioner_dFIdU_Roe(Jacobian[STENCIL_EAST],cell_index_i,cell_index_j,STENCIL_EAST);        
//   Preconditioner_dFIdU_Roe(Jacobian[STENCIL_WEST],cell_index_i,cell_index_j,STENCIL_WEST); 

//   //Center calculated from neighbours
//   //! Using the fact that dF/dU(right) = - dF/dU(left) 
//   Jacobian[STENCIL_CENTER] += (Jacobian[STENCIL_NORTH] + Jacobian[STENCIL_SOUTH] + Jacobian[STENCIL_EAST]  + Jacobian[STENCIL_WEST])
//                      /SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A;

//   Jacobian[STENCIL_NORTH] = -Jacobian[STENCIL_NORTH]/SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A;
//   Jacobian[STENCIL_SOUTH] = -Jacobian[STENCIL_SOUTH]/SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A;
//   Jacobian[STENCIL_EAST] = -Jacobian[STENCIL_EAST]/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A;
//   Jacobian[STENCIL_WEST] = -Jacobian[STENCIL_WEST]/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A;

}

/*****************************************************************************
 *  Calculate First Order Local Jacobian Block(s) Coresponding to Cell(i,j)  *
 *  using AUSM_plus_up                                                       *
 ****************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
First_Order_Inviscid_Jacobian_AUSM_plus_up(const int &cell_index_i,const int &cell_index_j,const int &cell_index_k,  
					   DenseMatrix* Jacobian){              

//   //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
//   Preconditioner_dFIdU_AUSM_plus_up(Jacobian[STENCIL_NORTH],cell_index_i,cell_index_j,STENCIL_NORTH);
//   Preconditioner_dFIdU_AUSM_plus_up(Jacobian[STENCIL_SOUTH],cell_index_i,cell_index_j,STENCIL_SOUTH); 
//   Preconditioner_dFIdU_AUSM_plus_up(Jacobian[STENCIL_EAST],cell_index_i,cell_index_j,STENCIL_EAST);        
//   Preconditioner_dFIdU_AUSM_plus_up(Jacobian[STENCIL_WEST],cell_index_i,cell_index_j,STENCIL_WEST); 

//   //! Using the fact that dF/dU(right) = - dF/dU(left) 
//   Jacobian[STENCIL_CENTER] += (Jacobian[STENCIL_NORTH] + Jacobian[STENCIL_SOUTH] + Jacobian[STENCIL_EAST]  + Jacobian[STENCIL_WEST])
//     /SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A;

//   Jacobian[STENCIL_NORTH] = -Jacobian[STENCIL_NORTH]/SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A;
//   Jacobian[STENCIL_SOUTH] = -Jacobian[STENCIL_SOUTH]/SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A;
//   Jacobian[STENCIL_EAST] = -Jacobian[STENCIL_EAST]/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A;
//   Jacobian[STENCIL_WEST] = -Jacobian[STENCIL_WEST]/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A;
  
}

/****************************************************************************
 *  Calculate Second Order Local Jacobian Block(s) Coresponding to Cell(i,j) *                      
 ****************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Second_Order_Viscous_Jacobian(const int &cell_index_i,const int &cell_index_j,const int &cell_index_k, DenseMatrix* Jacobian){

//   // A real cludge with all the DenseMatrices and recalculations, 
//   //  but just to test, need to change for performance....
//   DenseMatrix JacobianN(blocksize,blocksize,ZERO); 
//   DenseMatrix JacobianS(blocksize,blocksize,ZERO); 
//   DenseMatrix JacobianE(blocksize,blocksize,ZERO); 
//   DenseMatrix JacobianW(blocksize,blocksize,ZERO); 

//   //Also should rewrite to minimize dR/dU calls and just 
//   //call dWdU, but this needs to be done in Preconditioner_dFVdU
  
//   //***************** dR(i,j)/dU(i,j) *********************************************/
//   //STENCIL_CENTER
//   Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j,
// 		       cell_index_i,cell_index_j,STENCIL_NORTH,STENCIL_CENTER);
//   Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j,
// 		       cell_index_i,cell_index_j,STENCIL_SOUTH,STENCIL_CENTER);
//   Preconditioner_dFVdU(JacobianE,cell_index_i,cell_index_j,
// 		       cell_index_i,cell_index_j,STENCIL_EAST,STENCIL_CENTER);
//   Preconditioner_dFVdU(JacobianW,cell_index_i,cell_index_j,
// 		       cell_index_i,cell_index_j,STENCIL_WEST,STENCIL_CENTER);

//   Jacobian[STENCIL_CENTER] += (JacobianN + JacobianS + JacobianE  + JacobianW) 
//     /SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A;


//   /***************** dR(i,j-1)/dU(i,j) *********************************************/
//   //STENCIL_NORTH                           
//   JacobianN.zero();
//   Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j-1,
// 		       cell_index_i,cell_index_j,STENCIL_EAST,STENCIL_NORTH); 
//   Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j-1,
// 		       cell_index_i,cell_index_j,STENCIL_NORTH,STENCIL_NORTH); 
//   Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j-1,
// 		       cell_index_i,cell_index_j,STENCIL_WEST,STENCIL_NORTH); 
 
//   /***************** dR(i,j+1)/dU(i,j) *********************************************/
//   //STENCIL_SOUTH 
//   JacobianS.zero();
//   Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j+1,
// 		       cell_index_i,cell_index_j,STENCIL_EAST,STENCIL_SOUTH);
//   Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j+1,
// 		       cell_index_i,cell_index_j,STENCIL_SOUTH,STENCIL_SOUTH);
//   Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j+1,
// 		       cell_index_i,cell_index_j,STENCIL_WEST,STENCIL_SOUTH);
  
//   /***************** dR(i-1,j)/dU(i,j) *********************************************/
//   //STENCIL_EAST 
//   JacobianE.zero();
//   Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j,
// 		       cell_index_i,cell_index_j,STENCIL_NORTH,STENCIL_EAST);
//   Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j,
// 		       cell_index_i,cell_index_j,STENCIL_EAST,STENCIL_EAST);
//   Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j,
// 		       cell_index_i,cell_index_j,STENCIL_SOUTH,STENCIL_EAST);

//   /***************** dR(i+1,j)/dU(i,j) *********************************************/
//   //STENCIL_WEST
//   JacobianW.zero(); 
//   Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j,
// 		       cell_index_i,cell_index_j, STENCIL_NORTH,STENCIL_WEST);
//   Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j,
// 		       cell_index_i,cell_index_j, STENCIL_WEST,STENCIL_WEST);
//   Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j,
// 		       cell_index_i,cell_index_j, STENCIL_SOUTH,STENCIL_WEST);

//   Jacobian[STENCIL_NORTH] += JacobianN/SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A;
//   Jacobian[STENCIL_SOUTH] += JacobianS/SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A;
//   Jacobian[STENCIL_EAST] += JacobianE/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A;
//   Jacobian[STENCIL_WEST] += JacobianW/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A; 

//   /********************************************************************************/
//   //CORNERS reuse matrices
//   JacobianN.zero();  JacobianS.zero();
//   JacobianE.zero();  JacobianW.zero();

//   /***************** dR(i+1,j+1)/dU(i,j) *********************************************/
//   //STENCIL_SOUTHSTENCIL_WEST
//   Preconditioner_dFVdU(JacobianS,cell_index_i+1,cell_index_j+1,
// 		       cell_index_i,cell_index_j,STENCIL_SOUTH,SW);
//   Preconditioner_dFVdU(JacobianS,cell_index_i+1,cell_index_j+1,
// 		       cell_index_i,cell_index_j, STENCIL_WEST,SW);  
//   Jacobian[SW] += JacobianS/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j+1].A;
  
//   /***************** dR(i-1,j+1)/dU(i,j) *********************************************/
//   //STENCIL_SOUTHSTENCIL_EAST
//   Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j+1,
// 		       cell_index_i,cell_index_j,STENCIL_SOUTH,SE);
//   Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j+1,
// 		       cell_index_i,cell_index_j,STENCIL_EAST,SE);  
//   Jacobian[SE] += JacobianE/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j+1].A;

//   /***************** dR(i+1,j-1)/dU(i,j) *********************************************/
//   //STENCIL_NORTHSTENCIL_WEST
//   Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j-1,
// 		       cell_index_i,cell_index_j, STENCIL_NORTH,NW);
//   Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j-1,
// 		       cell_index_i,cell_index_j, STENCIL_WEST,NW);  
//   Jacobian[NW] += JacobianW/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j-1].A;

//   /***************** dR(i-1,j-1)/dU(i,j) *********************************************/
//   //STENCIL_NORTHSTENCIL_EAST
//   Preconditioner_dFVdU(JacobianN,cell_index_i-1,cell_index_j-1,
// 		       cell_index_i,cell_index_j, STENCIL_NORTH,NE);
//   Preconditioner_dFVdU(JacobianN,cell_index_i-1,cell_index_j-1,
// 		       cell_index_i,cell_index_j, STENCIL_EAST,NE);  
//   Jacobian[NE] += JacobianN/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j-1].A;


}

/*********************************************************
 *  Setup preconditioner from Jacobian Approx            *
 *********************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Setup_Preconditioner()
{
  deallocate_Precon(); // cludge of a workaround to fix up how BpPrecon uses memory in .setup
  
  ILUK_Precon = new BILUK;
  Jacobi_Precon = new BJacobi;
  
  if (Input->NKS_IP.GMRES_Block_Preconditioner == Block_ILUK ) {      
    /* ILU(k) */   
    ILUK_Precon->growth = Input->NKS_IP.GMRES_ILUK_Level_of_Fill+1; // BPKIT HACK -> increases storage needed for high level of fill
    ILUK_Precon->localprecon(LP_INVERSE);  // currently using "exact" local inverse, many other options which may be cheaper !!!!!
    ILUK_Precon->setup(Block_Jacobian_approx, Input->NKS_IP.GMRES_ILUK_Level_of_Fill);
  } else if (Input->NKS_IP.GMRES_Block_Preconditioner == Block_Jacobi) {
    /* Diagonal */    
    Jacobi_Precon->localprecon(LP_INVERSE);
    Jacobi_Precon->setup(Block_Jacobian_approx);
  }

}

/*********************************************************
 * Apply preconditioner - in GMRES                       *
 *********************************************************/  //FYI Issues with passing by address &
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline void Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Apply_Preconditioner(int A, int B, const double *C, int D, double *E, int F)
{
  /* z = Minv * V(i) -> stored in W(i). */
  if (Input->NKS_IP.GMRES_Block_Preconditioner == Block_ILUK ) {
    ILUK_Precon->apply(A, B, C, D , E, F);
  } else if (Input->NKS_IP.GMRES_Block_Preconditioner == Block_Jacobi ) {
    Jacobi_Precon->apply(A, B, C, D , E, F);
  }
}

/********************************************************
 * Block_Preconditioner::DenseMat_to_DenseMatrix.       *
 ********************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline DenseMat Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
DenseMatrix_to_DenseMat(const DenseMatrix &B) {
  DenseMat A(B.dim(0),B.dim(1));
  for (int i=0; i<B.dim(0); i++) {
    for( int j=0; j<B.dim(1); j++) {
      A(i,j) = B(i,j); 
    }
  }
  return A;
} 


/*********************************************************
 * Routine: Rotation_Matrix_3D                           *
 *                                                       *
 * This function returns either the rotation matrix, A,  *
 * or the inverse of A.                                  *
 *                                                       *
 * Note: A_matrix = 1 for returning A.                   *
 *                = 0 for returning inverse of A.        *
 *                                                       *
 *********************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline DenseMatrix Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Rotation_Matrix_3D(const Vector3D &nface, const int &A_matrix) 
{

  // for a 3D unit normal rotated to align with the x-axis
  double cos_angle = nface.x;
  double sin_angle = sqrt( nface.y*nface.y + nface.z*nface.z);

  Vector3D rot_axis(0,nface.z,-nface.y);

  DenseMatrix mat(blocksize,blocksize); //TEMP
  mat.identity();

  if (A_matrix) {             
    // Rotation Matrix, A                 
    mat(1,1) = cos_angle;
    mat(1,2) = -rot_axis.z*sin_angle;
    mat(1,3) = rot_axis.y*sin_angle;

    mat(2,1) = rot_axis.z*sin_angle;
    mat(2,2) = rot_axis.y*rot_axis.y*(ONE-cos_angle)+cos_angle;
    mat(2,3) = rot_axis.y*rot_axis.z*(ONE-cos_angle);    

    mat(3,1) = -rot_axis.y*sin_angle;
    mat(3,2) = rot_axis.y*rot_axis.z*(ONE-cos_angle);
    mat(3,3) = rot_axis.z*rot_axis.z*(ONE-cos_angle)+cos_angle;    

    //Inverse
  }else {
    mat(1,1) = cos_angle;
    mat(1,2) = rot_axis.z*sin_angle;
    mat(1,3) = -rot_axis.y*sin_angle;

    mat(2,1) = -rot_axis.z*sin_angle;
    mat(2,2) = rot_axis.y*rot_axis.y*(ONE-cos_angle)+cos_angle;
    mat(2,3) = rot_axis.y*rot_axis.z*(ONE-cos_angle);    

    mat(3,1) = rot_axis.y*sin_angle;
    mat(3,2) = rot_axis.y*rot_axis.z*(ONE-cos_angle);
    mat(3,3) = rot_axis.z*rot_axis.z*(ONE-cos_angle)+cos_angle;     

    //mat.pseudo_inverse_override();
  } 

  return mat;

} /* End of Rotation_Matrix. */

/*********************************************************
 * Routine: Rotation_Matrix                              *
 *                                                       *
 * This function returns either the rotation matrix, A,  *
 * or the inverse of A.                                  *
 *                                                       *
 * Note: A_matrix = 1 for returning A.                   *
 *                = 0 for returning inverse of A.        *
 *                                                       *
 *********************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
inline DenseMatrix Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>::
Rotation_Matrix_2D(const Vector2D &nface, const int &A_matrix) 
{
  
  double cos_angle = nface.x; 
  double sin_angle = nface.y;

  DenseMatrix mat(blocksize,blocksize); //TEMP
  mat.identity();

  if (A_matrix) {             
    // Rotation Matrix, A                  
    mat(1,1) = cos_angle;
    mat(1,2) = sin_angle;
    mat(2,1) = -sin_angle;
    mat(2,2) = cos_angle;
  } else {
    // Rotation Matrix, Inv A 
    mat(1,1) = cos_angle;
    mat(1,2) = -sin_angle;
    mat(2,1) = sin_angle;
    mat(2,2) = cos_angle;
  } /* endif */

  return mat;

} /* End of Rotation_Matrix. */

#endif  // _BLOCK_PRECONDITONER_INCLUDED 





  /************* DEBUGGING FOR FIRST ORDER *********************
  cout<<"\n CELL "<<cell_index_i<<" "<<cell_index_j;  
  cout<<"\n faces \n N"<< nface_N <<"\n S "<<nface_S<<"\n E "<<nface_E<<"\n W "<<nface_W;
  cout<<"\n lambdas \n N"<< lambdas_N <<"\n S "<<lambdas_S<<"\n E "<<lambdas_E<<"\n W "<<lambdas_W;
  cout<<"\n dFdUs \n N \n"<< dFdU_N <<"\n S \n"<<dFdU_S<<"\n E \n"<<dFdU_E<<"\n W \n"<<dFdU_W;
  cout<<"\n Jacobians \n N \n"<<Jacobian[1]<<"\n S \n"<<Jacobian[3]
      <<"\n E \n"<<Jacobian[4]<<"\n W \n"<<Jacobian[2]<<"\n C \n"<<Jacobian[0];
  *********************************************/


//   /************************************************************************/
//   /*********************** TEST OF BPKIT DATA STRUCTURE *******************/
//   blocksize = 2;
//   int nrow = 4;
//   int nnz = 6;
//   int *row= new int[nnz];   //row index of block
//   int *col = new int[nnz];  //column index of block
//   double *A = new double[nnz*blocksize*blocksize]; //1D array of values corresponding to (row,col) indices
  
//   for(int i=0; i<blocksize*blocksize; i++){
//     A[i] = ONE;
//     A[i+blocksize*blocksize] = 2.0;
//     A[i+2*blocksize*blocksize] = 3.0;
//     A[i+3*blocksize*blocksize] = 4.0;
//     A[i+4*blocksize*blocksize] = 5.0;
//     A[i+5*blocksize*blocksize] = 6.0;
//   }  
//   col[0] = 1;
//   col[1] = 2;
//   col[2] = 3;
//   col[3] = 4;
//   col[4] = 3;
//   col[5] = 1;
//   row[0] = 1;
//   row[1] = 2;
//   row[2] = 3;
//   row[3] = 4;
//   row[4] = 1;
//   row[5] = 3;
  
//   Block_Jacobian_approx.setup(nrow,nnz,row,col,A,blocksize);                

//   cout<<"\n Test of BlockMat\n "<<Block_Jacobian_approx;

//   delete[] row; delete[] col; delete[] A;
 
//   DenseMat Q(2,2);
//   Q(0,0) = 1.5;
//   Q(1,0) = 2.5;
//   Q(0,1) = 3.5;
//   Q(1,1) = 4.5;
//   cout<<"\n  TEST replacement block at 2,0 \n"<<Q;

//   Block_Jacobian_approx.setblock(2,0,Q);
  
//   cout<<"\n Test of BlockMat after modification\n "<<Block_Jacobian_approx;
//   cout<<"\n End of Test \n\n"; exit(1); 

//   /*********************** TEST OF BPKIT DATA STRUCTURE *******************/
//   /************************************************************************/

