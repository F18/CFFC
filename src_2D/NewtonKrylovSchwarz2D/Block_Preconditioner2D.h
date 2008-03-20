#ifndef _BLOCK_PRECONDITONER_INCLUDED 
#define _BLOCK_PRECONDITONER_INCLUDED 

/* BPKIT data structures and precondtioners along with BLAS 
   FORTRAN wrappers */ 

#include "BlockMat.h"
#include "DenseMat.h"
#include "BILUK.h"
#include "BRelax.h" //For BJacobi
#include "../Math/Matrix.h" //For DenseMatrix

#ifndef _NKS_DTS_INCLUDED
#include "DTS_NKS.h"
#endif

/********************************************************
 * Class: Block_Preconditioner                          *
 *                                                      *
 * Member functions                                     * 
 *                                                      *
 ********************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
class Block_Preconditioner {
  private:
  int blocksize;
  int Jacobian_stencil_size; 

  void Get_Block_Index(const int &i, int *);
  void Get_Block_Index(const int &,const int &, int *, int *);
  void Pre_Precon_SolnBlk_Init(void);
  void Setup_Jacobian_approximation();          //called from Create_Preconditioner
  void Setup_Preconditioner();                  //called from Update_Jacobian
  void Implicit_Euler(const int&,const int&, DenseMatrix*, const double& );     //called from Update_Jacobian       
  void First_Order_Inviscid_Jacobian_HLLE(const int&,const int&, DenseMatrix*); //called from Update_Jacobian
  void First_Order_Inviscid_Jacobian_Roe(const int&,const int&, DenseMatrix*);  //called from Update_Jacobian
  void First_Order_Inviscid_Jacobian_AUSM_plus_up(const int&,const int&, DenseMatrix*);  //called from Update_Jacobian
  void Second_Order_Viscous_Jacobian(const int&,const int&, DenseMatrix*);      //called from Update_Jacobian

  DenseMatrix Rotation_Matrix(const Vector2D &nface, const int &A_matrix); //Used in Jacobian_LocalBlock
  DenseMat DenseMatrix_to_DenseMat(const DenseMatrix &B); 

  protected:
  public:
  //Address of corresponding Solution block data
  SOLN_BLOCK_TYPE *SolnBlk;
  INPUT_TYPE *Input_Parameters;  

  // BPKIT BlockMat storage of Approximate Jacobian for use with Preconditioner  
  BlockMat Block_Jacobian_approx;

  // BPKIT BpPrecon Preconditioner Types
  BILUK   *ILUK_Precon;
  BJacobi *Jacobi_Precon;

  /*******************************************************/
  //default constructors
  Block_Preconditioner(void):Jacobian_stencil_size(5), blocksize(1),
    SolnBlk(NULL),Input_Parameters(NULL), ILUK_Precon(NULL), Jacobi_Precon(NULL) {}
  Block_Preconditioner(SOLN_BLOCK_TYPE  &Soln_ptr,INPUT_TYPE &IP, const int &_blocksize)		      
  { Create_Preconditioner(Soln_ptr,IP,_blocksize); }
  
  //Constructors
  //These should be generalized for any precondtioner, not just "BPkit" types....
  void Create_Preconditioner(SOLN_BLOCK_TYPE  &Soln_ptr, INPUT_TYPE &IP, const int &_blocksize);
  void Update_Jacobian_and_Preconditioner(const double &DTS_dTime);  
  void Apply_Preconditioner(int A, int B, const double *C, int D, double *E, int F); 
  
  //Equation Specific Specializations required for Inviscid, Viscous, and Source Term Jacobians
  //Inviscid
  void Preconditioner_dFIdU(DenseMatrix &dFdU, SOLN_VAR_TYPE W); 
  void Preconditioner_dFIdU_Roe(DenseMatrix &, const int, const int, const int);  
  void Preconditioner_dFIdU_AUSM_plus_up(DenseMatrix &, const int, const int, const int);  
  //Viscous
  void Preconditioner_dFVdU(DenseMatrix &, const int, const int, const int, const int, const int ,const int);
  //Source 
  void Preconditioner_dSdU(int cell_index_i, int cell_index_j, DenseMatrix &Jacobian);
  void normalize_Preconditioner_dFdU(DenseMatrix &dFdU);

  //Destructor / memory cleanup
  ~Block_Preconditioner(){ deallocate();  deallocate_Precon(); }   //??ISSUES WITH !=NULL DTS WITH AMR???
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
template<typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Pre_Precon_SolnBlk_Init(void) {
  // Do Nothing
}

template<typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Preconditioner_dFIdU(DenseMatrix &dFdU, SOLN_VAR_TYPE W) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF Preconditioner_dFIdU for Block_Preconditioner2D.h requried \n";
  exit(1);
}

template<typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Preconditioner_dFIdU_Roe(DenseMatrix &dFdU, const int, const int, const int) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF Preconditioner_dFIdU_Roe for Block_Preconditioner2D.h requried \n";
  exit(1);
}

template<typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Preconditioner_dFIdU_AUSM_plus_up(DenseMatrix &dFdU, const int, const int, const int) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF Preconditioner_dFIdU_AUSM_plus_up for Block_Preconditioner2D.h requried \n";
  exit(1);
}

//  Determines the Jacobian of the part of the equation written for
//  the "Orient_face" face of cell (ri, rj) with respect to the
//  variables of cell (wi, wj). The cell at (wi, wj) is to the
//  "Orient_cell" of the cell at (ri, rj). Remember that j increases
//  to the north while i increases to the east.
//  
//  For example, suppose that 
//  (ri, rj) = (ci-1, cj)
//  (wi, wj) = (ci  , cj)
//  Orient_face = North
//  Orient_cell = East
//  
//  Here, the equation written for the north face of the cell to the
//  west (ci-1, cj) will contain terms that depend on (ci, cj). We
//  will determine the derivate of those terms with respect to the
//  variables in cell (ci, cj). Notice that (ci, cj) is to the east
//  of (ci-1, cj).
//  
//  For another example, suppose that 
//  (ri, rj) = (ci, cj-1)
//  (wi, wj) = (ci, cj  )
//  Orient_face = North
//  Orient_cell = North
//  
//  Here, the equation written for the north face of the cell to the
//  south (ci, cj-1) will contain terms that depend on (ci, cj).
//  Notice that (ci, cj) is to the north of (ci, cj-1).  
//  
//  -- Code by Scott Northrup, this text by Alistair Wood.
template<typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Preconditioner_dFVdU(DenseMatrix &dFvdU, const int Rii, const int Rjj, 
		     const int Wii, const int Wjj, const int Orient_face, const int Orient_cell) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF Preconditioner_dFVdU for Block_Preconditioner2D.h requried \n";
  exit(1);
}

// BLANK UNLESS OTHERWISE SPECIALIZED
template<typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Preconditioner_dSdU( int cell_index_i, int cell_index_j, DenseMatrix &Jacobian) {
//   cerr<<"\n EXPLICIT SPECIALIZATION OF Preconditioner_dSdU for Block_Preconditioner2D.h requried \n";
//   exit(1);
}

template<typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) {
  cerr<<"\n EXPLICIT SPECIALIZATION OF normalize_Preconditioner_dFdU for Block_Preconditioner2D.h requried \n";
  exit(1);
}


/*!**************************************************************
 * Generate block matrix indexing for a given matrix row "i"    *
 ****************************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Get_Block_Index(const int &i, int *Block_index_j)
{
  //Determine Block_Indicies for Block Matrix formation
  if( Jacobian_stencil_size == 5){        
    Block_index_j[NORTH] = i - SolnBlk->NCi ; 
    Block_index_j[EAST] = ( i%SolnBlk->NCi != 0) ? i-1 : -1;
    Block_index_j[CENTER] = i;     
    Block_index_j[WEST] = ( (i+1)%SolnBlk->NCi != 0) ? i+1 : -1;
    Block_index_j[SOUTH] = i + SolnBlk->NCi ;   
  }  else if ( Jacobian_stencil_size == 9) {
    Block_index_j[NORTH_EAST] = ( (i - SolnBlk->NCi)%SolnBlk->NCi != 0) ? i - SolnBlk->NCi - 1 : -1;
    Block_index_j[NORTH] = i - SolnBlk->NCi;
    Block_index_j[NORTH_WEST] = ( (i - SolnBlk->NCi+1)%SolnBlk->NCi != 0) ? i - SolnBlk->NCi + 1 : -1;
    Block_index_j[EAST] = ( i%SolnBlk->NCi != 0) ? i-1 : -1;
    Block_index_j[CENTER] = i;
    Block_index_j[WEST] = ( (i+1)%SolnBlk->NCi != 0) ? i+1 : -1;
    Block_index_j[SOUTH_EAST] = ( (i + SolnBlk->NCi)%SolnBlk->NCi != 0) ? i + SolnBlk->NCi - 1 : -1;
    Block_index_j[SOUTH] = i + SolnBlk->NCi;
    Block_index_j[SOUTH_WEST] = ( (i + SolnBlk->NCi+1)%SolnBlk->NCi != 0) ? i + SolnBlk->NCi + 1 : -1;
  }  
}

/****************************************************************
 * Generate block matrix indexing for a given cell (i,j)        *
 ****************************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Get_Block_Index(const int &cell_index_i,const int &cell_index_j, int *Block_index_i, int *Block_index_j)
{ 
  //J index 
  for( int i=0; i<Jacobian_stencil_size; i++)  
    Block_index_j[i] = cell_index_j*SolnBlk->NCi+cell_index_i; 

  //I index 
  if( Jacobian_stencil_size == 5){        
    /*! Determine 1st order Block_Indicies for Cell (i,j)   
     *
     *            ---  
     *           | 1 |
     *        --- --- ---
     *       | 2 | 0 | 4 |                
     *        --- --- ---
     *           | 3 |
     *            ---
     */
    Block_index_i[NORTH] = (cell_index_j - 1)*SolnBlk->NCi+cell_index_i;  //NORTH  
    Block_index_i[EAST] = cell_index_j*SolnBlk->NCi+cell_index_i-1;       //EAST     
    Block_index_i[CENTER] = cell_index_j*SolnBlk->NCi+cell_index_i;       //CENTER   
    Block_index_i[WEST] = cell_index_j*SolnBlk->NCi+cell_index_i+1;       //WEST         
    Block_index_i[SOUTH] = (cell_index_j + 1)*SolnBlk->NCi+cell_index_i;  //SOUTH    
    
  }  else if ( Jacobian_stencil_size == 9) {
    /*! Determine 2nd order Block_Indicies for Cell (i,j)   
     *
     *        --- --- ---
     *       | 6 | 1 | 5 |
     *        --- --- --- 
     *       | 2 | 0 | 4 | 
     *       --- --- --- -
     *       | 8 | 3 | 7 |
     *        --- --- ---
     */ 
    Block_index_i[NORTH_EAST] = (cell_index_j - 1)*SolnBlk->NCi+cell_index_i-1;    //NORTH_EAST  
    Block_index_i[NORTH] = (cell_index_j - 1)*SolnBlk->NCi+cell_index_i;   //NORTH
    Block_index_i[NORTH_WEST] = (cell_index_j - 1)*SolnBlk->NCi+cell_index_i+1;    //NORTH_WEST   
    Block_index_i[EAST] = cell_index_j*SolnBlk->NCi+cell_index_i-1;        //EAST      
    Block_index_i[CENTER] = cell_index_j*SolnBlk->NCi+cell_index_i;        //CENTER
    Block_index_i[WEST] = cell_index_j*SolnBlk->NCi+cell_index_i+1;        //WEST
    Block_index_i[SOUTH_EAST] = (cell_index_j + 1)*SolnBlk->NCi+cell_index_i-1;    //SOUTH_EAST     
    Block_index_i[SOUTH] = (cell_index_j + 1)*SolnBlk->NCi+cell_index_i;   //SOUTH    
    Block_index_i[SOUTH_WEST] = (cell_index_j + 1)*SolnBlk->NCi+cell_index_i+1;    //SOUTH_WEST  
  }  

}


/*********************************************************
 *  Create the Preconditoner as required by BPKIT        *
 *********************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Create_Preconditioner( SOLN_BLOCK_TYPE  &Soln_ptr, INPUT_TYPE &IP, const int &_blocksize)
{
  SolnBlk = &(Soln_ptr);
  Input_Parameters= &(IP);
  blocksize = _blocksize;       //number of equations

  if (IP.NKS_IP.Jacobian_Order == SOURCE_TERMS_ONLY ||
      IP.NKS_IP.Jacobian_Order == FIRST_ORDER_INVISCID_HLLE ||
      IP.NKS_IP.Jacobian_Order == FIRST_ORDER_INVISCID_ROE || 
      IP.NKS_IP.Jacobian_Order == FIRST_ORDER_INVISCID_AUSMPLUSUP) {
    Jacobian_stencil_size = 5;
  } else if (IP.NKS_IP.Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_HLLE ||
             IP.NKS_IP.Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_ROE ||
             IP.NKS_IP.Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP) {
    Jacobian_stencil_size = 9;
  } else {
    cerr<<"\n Invalid Jacobian Preconditioner Order "<<IP.NKS_IP.Jacobian_Order;
    exit(1);
  }
    
  //Sets up memory for Approximate Jacobian (Block_Jacobian_approx)
  Setup_Jacobian_approximation(); 
}

/********************************************************************************
 *  Jacobian (dR/dU) Approximation Stencil Allocation - only called on STARTUP   *
 *********************************************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Setup_Jacobian_approximation(){

  int block_mat_size = SolnBlk->NCi*SolnBlk->NCj;       //block matrix size based on icells*jcells
  int nnz;

  //! 1st order 2D Jacobian, Sparse Block Matrix Memory Allocation of non-zero block entries
  if( Jacobian_stencil_size == 5){        
    nnz = 5*block_mat_size - 2*(SolnBlk->NCi + SolnBlk->NCj); 
    //! 2nd order 2D Jacobian
  }  else if ( Jacobian_stencil_size == 9) {
    nnz = 9*block_mat_size - 6*(SolnBlk->NCi + SolnBlk->NCj) + 4;    
  }      

  //!Temporary Storage Arrays
  int *i_index = new int[nnz];      //location of dense local blocks, in global sparse block Matrix   
  int *j_index = new int[nnz];                            //TEMP VAR      
  double *Data = new double [nnz*blocksize*blocksize];    //TEMP VAR
  int *block_i = new int[Jacobian_stencil_size];          //TEMP VAR
  
  //!Setup Stencil for approximate Jacobian
  int nnz_count = 0;   
  //! Loop through each row of Block Matrix
  for (int i=0; i< block_mat_size; i++){       
    Get_Block_Index(i,block_i);
    int stencil = 0;
    //! Determine which entries correspond to 1st, 2nd, etc. order stencil    
    //cout<<"\n i "<<i;
    while( stencil < Jacobian_stencil_size) { 
      int j = block_i[stencil];
      //cout<<" "<<j;
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
	nnz_count++;
      }
      stencil++;	
    }      
  }

  if( nnz != nnz_count){ cerr<<"\n Number of nonzero blocks mismatch error in approximate Jacobian formation "
			     <<nnz<<" != "<<nnz_count<<"\n"; exit(1); }

  //! Create Sparse bpkit "BlockMat" Block Matrix 
  Block_Jacobian_approx.setup(block_mat_size,nnz,i_index,j_index,Data,blocksize);                  
  
  //! Clean up local memory
  delete[] i_index; delete[] j_index; delete[] Data;  delete[] block_i;

}


/**********************************************************
 *  Update BlockMat with approximation to the Jacobian    *
 **********************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Update_Jacobian_and_Preconditioner(const double &DTS_dTime)
{

  //!Local Variables and Temporary Storage
  int block_mat_size = SolnBlk->NCi*SolnBlk->NCj; 
  DenseMatrix *Jacobian_Data = new DenseMatrix[Jacobian_stencil_size];        //TEMP VAR  
  for(int i=0; i<Jacobian_stencil_size; i++) { Jacobian_Data[i] = DenseMatrix(blocksize,blocksize,ZERO); }
  int *block_i = new int[Jacobian_stencil_size];                                //TEMP VAR
  int *block_j = new int[Jacobian_stencil_size];                                  //TEMP VAR
  
  //! Initially assume no overlap
  int JCl_overlap = 0, JCu_overlap = 0, ICu_overlap =0, ICl_overlap = 0;
  
  //! If overlap determine which block boundaries are internal, ie. BC_NONE
  if(Input_Parameters->NKS_IP.GMRES_Overlap){	
    if (SolnBlk->Grid.BCtypeS[SolnBlk->ICl] == BC_NONE)  JCl_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
    if (SolnBlk->Grid.BCtypeN[SolnBlk->ICu] == BC_NONE)  JCu_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
    if (SolnBlk->Grid.BCtypeE[SolnBlk->JCu] == BC_NONE)  ICu_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
    if (SolnBlk->Grid.BCtypeW[SolnBlk->JCl] == BC_NONE)  ICl_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
  }

  //! Perform any specialized initialization of the solution block
  //! before building the preconditioner
  Pre_Precon_SolnBlk_Init();

  //*********************************************************************************//
  /*! Calculate Jacobians for each cell and Update Global Jacobian Block Matrix
   * loop through all non-ghost cells, including overlap cells.  Ghost cells already
   * set to Zero or Identity by initialization.
   **********************************************************************************/
  for(int i= SolnBlk->ICl - ICl_overlap; i<= SolnBlk->ICu + ICu_overlap; i++){    
    for(int j= SolnBlk->JCl - JCl_overlap; j<= SolnBlk->JCu + ICu_overlap; j++){  
         
      //--------------------------------------------------------------------------//
      //! Calculate Local Approximate Jacobian                        
      switch(Input_Parameters->NKS_IP.Jacobian_Order){
      case SOURCE_TERMS_ONLY :
 	Implicit_Euler(i,j, Jacobian_Data,DTS_dTime);
	Preconditioner_dSdU(i,j,Jacobian_Data[CENTER]);
	break;
      case FIRST_ORDER_INVISCID_HLLE : 
	Implicit_Euler(i,j, Jacobian_Data,DTS_dTime);
	First_Order_Inviscid_Jacobian_HLLE(i,j, Jacobian_Data);
	Preconditioner_dSdU(i,j,Jacobian_Data[CENTER]);                       
	break;
      case FIRST_ORDER_INVISCID_ROE :  
	Implicit_Euler(i,j, Jacobian_Data,DTS_dTime); 
	First_Order_Inviscid_Jacobian_Roe(i,j, Jacobian_Data);
	Preconditioner_dSdU(i,j,Jacobian_Data[CENTER]);
	break;
      case FIRST_ORDER_INVISCID_AUSMPLUSUP : 
	Implicit_Euler(i,j, Jacobian_Data,DTS_dTime);
	First_Order_Inviscid_Jacobian_AUSM_plus_up(i,j, Jacobian_Data);
	Preconditioner_dSdU(i,j,Jacobian_Data[CENTER]);
	break;
      case SECOND_ORDER_DIAMOND_WITH_HLLE:
	Implicit_Euler(i,j, Jacobian_Data,DTS_dTime);
	First_Order_Inviscid_Jacobian_HLLE(i,j, Jacobian_Data);   
	Second_Order_Viscous_Jacobian(i,j, Jacobian_Data);    
	Preconditioner_dSdU(i,j,Jacobian_Data[CENTER]);   
	break;
      case SECOND_ORDER_DIAMOND_WITH_ROE :	
	Implicit_Euler(i,j, Jacobian_Data,DTS_dTime);
	First_Order_Inviscid_Jacobian_Roe(i,j, Jacobian_Data);
	Second_Order_Viscous_Jacobian(i,j, Jacobian_Data);    
	Preconditioner_dSdU(i,j,Jacobian_Data[CENTER]); 
	break;
      case SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP :	
	Implicit_Euler(i,j, Jacobian_Data,DTS_dTime);
	First_Order_Inviscid_Jacobian_AUSM_plus_up(i,j, Jacobian_Data);
	Second_Order_Viscous_Jacobian(i,j, Jacobian_Data);    
	Preconditioner_dSdU(i,j,Jacobian_Data[CENTER]);  
	break;
      }

      //--------------------------------------------------------------------------//
      //! Get Block Matrix locations that have components from a given Cell(i,j)
      Get_Block_Index(i,j, block_i, block_j);
      
      //fudge for iGhost Cell reset to zero   //jGhost cell already zero                
      if(block_i[NORTH] < TWO*SolnBlk->NCi) Jacobian_Data[NORTH].zero();
      if(block_i[SOUTH] > block_mat_size - TWO*SolnBlk->NCi) Jacobian_Data[SOUTH].zero();

      if(Jacobian_stencil_size == 9){		
	if(block_i[NORTH_EAST] < TWO*SolnBlk->NCi)    Jacobian_Data[NORTH_EAST].zero();
	if(block_i[NORTH_WEST] < TWO*SolnBlk->NCi)    Jacobian_Data[NORTH_WEST].zero(); 
	if(block_i[SOUTH_EAST] > block_mat_size - TWO*SolnBlk->NCi)    Jacobian_Data[SOUTH_EAST].zero();
	if(block_i[SOUTH_WEST] > block_mat_size - TWO*SolnBlk->NCi)    Jacobian_Data[SOUTH_WEST].zero();
      }
      //--------------------------------------------------------------------------//
      //! Update BlockMat with Local Approximate Jacobians 
      for( int block = 0; block < Jacobian_stencil_size; block++){
	// Normalize
	normalize_Preconditioner_dFdU(Jacobian_Data[block]);

	//can be sped up by more intelligent logic in bkpkit (BlockMat.cc  "setblock")
	Block_Jacobian_approx.setblock( block_i[block], block_j[block], DenseMatrix_to_DenseMat(Jacobian_Data[block]));

	Jacobian_Data[block].zero(); //Just in case to avoid +=/-= issues
      }     
    }
  }

  //Local Memory cleanup
  delete[] Jacobian_Data; delete[] block_i; delete[] block_j;

  //Setup appropriate Preconditioner after Jacobian has been formed/Updated
  Setup_Preconditioner();
}

/*****************************************************************************
 *  Add finite time step to diagonal.                                        *
 *****************************************************************************/ 
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Implicit_Euler(const int &cell_index_i,const int &cell_index_j, DenseMatrix* Jacobian,const double& DTS_dTime)
{   
  DenseMatrix Diag(blocksize,blocksize);      
  Diag.identity();    
  //Cacluate LHS depeneding on Steady State of Dual Time Stepping
  Diag *= LHS_Time<INPUT_TYPE>(*Input_Parameters, SolnBlk->dt[cell_index_i][cell_index_j],DTS_dTime);  
  Jacobian[CENTER] -= Diag;
}

/*****************************************************************************
 *  Calculate First Order Local Jacobian Block(s) Coresponding to Cell(i,j)  *
 *  using HLLE                                                               *
 *****************************************************************************/ 
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
First_Order_Inviscid_Jacobian_HLLE(const int &cell_index_i,const int &cell_index_j, 
				   DenseMatrix* Jacobian){              
  
  //! Caculate normal vectors -> in Vector2D format. 
  Vector2D nface_N = SolnBlk->Grid.nfaceN(cell_index_i,cell_index_j-1);
  Vector2D nface_S = SolnBlk->Grid.nfaceS(cell_index_i,cell_index_j+1);
  Vector2D nface_E = SolnBlk->Grid.nfaceE(cell_index_i-1,cell_index_j);
  Vector2D nface_W = SolnBlk->Grid.nfaceW(cell_index_i+1,cell_index_j);

  //! Calculate wavespeeds using solutions in the rotated frame -> in Vector2D format.
  Vector2D lambdas_N = HLLE_wavespeeds(SolnBlk->W[cell_index_i][cell_index_j-1], 
				       SolnBlk->W[cell_index_i][cell_index_j], nface_N);
  Vector2D lambdas_S = HLLE_wavespeeds(SolnBlk->W[cell_index_i][cell_index_j+1], 
				       SolnBlk->W[cell_index_i][cell_index_j], nface_S);  
  Vector2D lambdas_E = HLLE_wavespeeds(SolnBlk->W[cell_index_i-1][cell_index_j], 
				       SolnBlk->W[cell_index_i][cell_index_j], nface_E);
  Vector2D lambdas_W = HLLE_wavespeeds(SolnBlk->W[cell_index_i+1][cell_index_j], 
				       SolnBlk->W[cell_index_i][cell_index_j], nface_W);

  //checks necessary ????
  //if ((lambdas_W.y-lambdas_W.x) == ZERO) cout << " WEST : HLLE_wavespeeds " << endl;
 
  //! Calculate constants gamma and beta -> scalar values. 
  double gamma_N = (lambdas_N.x*lambdas_N.y)/(lambdas_N.y-lambdas_N.x);
  double beta_N  = - lambdas_N.x/(lambdas_N.y-lambdas_N.x);
  double gamma_S = (lambdas_S.x*lambdas_S.y)/(lambdas_S.y-lambdas_S.x);
  double beta_S  = - lambdas_S.x/(lambdas_S.y-lambdas_S.x);
  double gamma_E = (lambdas_E.x*lambdas_E.y)/(lambdas_E.y-lambdas_E.x);
  double beta_E  = - lambdas_E.x/(lambdas_E.y-lambdas_E.x);
  double gamma_W = (lambdas_W.x*lambdas_W.y)/(lambdas_W.y-lambdas_W.x);
  double beta_W  = - lambdas_W.x/(lambdas_W.y-lambdas_W.x);

  //! Obtain rotation matrices with normal vector -> matrices in DenseMatrix format. 
  DenseMatrix A_N( Rotation_Matrix(nface_N, 1) );           //TEMP VAR
  DenseMatrix AI_N( Rotation_Matrix(nface_N, 0) );            //TEMP VAR
  DenseMatrix A_S( Rotation_Matrix(nface_S, 1) );            //TEMP VAR
  DenseMatrix AI_S( Rotation_Matrix(nface_S, 0) );            //TEMP VAR
  DenseMatrix A_E( Rotation_Matrix(nface_E, 1) );             //TEMP VAR 
  DenseMatrix AI_E( Rotation_Matrix(nface_E, 0) );            //TEMP VAR
  DenseMatrix A_W( Rotation_Matrix(nface_W, 1) );             //TEMP VAR  
  DenseMatrix AI_W( Rotation_Matrix(nface_W, 0) );            //TEMP VAR      

  //! Calculate dFdU using solutions in the rotated frame -> matrix in DenseMatrix format. 
  DenseMatrix dFdU_N(blocksize,blocksize,ZERO);            //TEMP VAR  
  DenseMatrix dFdU_S(blocksize,blocksize,ZERO);            //TEMP VAR
  DenseMatrix dFdU_E(blocksize,blocksize,ZERO);             //TEMP VAR
  DenseMatrix dFdU_W(blocksize,blocksize,ZERO);             //TEMP VAR  

  //Solution Rotate provided in pState 
  Preconditioner_dFIdU(dFdU_N, Rotate(SolnBlk->W[cell_index_i][cell_index_j], nface_N)); 
  Preconditioner_dFIdU(dFdU_S, Rotate(SolnBlk->W[cell_index_i][cell_index_j], nface_S));
  Preconditioner_dFIdU(dFdU_E, Rotate(SolnBlk->W[cell_index_i][cell_index_j], nface_E));
  Preconditioner_dFIdU(dFdU_W, Rotate(SolnBlk->W[cell_index_i][cell_index_j], nface_W));
  
  DenseMatrix II(blocksize,blocksize);  II.identity();    

  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  //North
  Jacobian[NORTH] = (SolnBlk->Grid.lfaceN(cell_index_i,cell_index_j-1) 
		 * AI_N * (beta_N * dFdU_N + gamma_N * II) * A_N); 

  //South
  Jacobian[SOUTH] = (SolnBlk->Grid.lfaceS(cell_index_i,cell_index_j+1) 
		 * AI_S * (beta_S * dFdU_S + gamma_S * II) * A_S);

  //East
  Jacobian[EAST] = (SolnBlk->Grid.lfaceE(cell_index_i-1,cell_index_j) 
		 * AI_E * (beta_E * dFdU_E + gamma_E * II) * A_E);

  //West
  Jacobian[WEST] = (SolnBlk->Grid.lfaceW(cell_index_i+1,cell_index_j) 
		 * AI_W * (beta_W * dFdU_W + gamma_W * II) * A_W);

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left) 
  Jacobian[CENTER] += (Jacobian[NORTH] + Jacobian[SOUTH] + Jacobian[EAST]  + Jacobian[WEST])
    /SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A;

  Jacobian[NORTH] = -Jacobian[NORTH]/SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A;
  Jacobian[SOUTH] = -Jacobian[SOUTH]/SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A;
  Jacobian[EAST] = -Jacobian[EAST]/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A;
  Jacobian[WEST] = -Jacobian[WEST]/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A;

}

/*****************************************************************************
 *  Calculate First Order Local Jacobian Block(s) Coresponding to Cell(i,j)  *
 *  using Roe                                                                *
 *****************************************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
First_Order_Inviscid_Jacobian_Roe(const int &cell_index_i,const int &cell_index_j, 
				   DenseMatrix* Jacobian){              
    
  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  Preconditioner_dFIdU_Roe(Jacobian[NORTH],cell_index_i,cell_index_j,NORTH);
  Preconditioner_dFIdU_Roe(Jacobian[SOUTH],cell_index_i,cell_index_j,SOUTH); 
  Preconditioner_dFIdU_Roe(Jacobian[EAST],cell_index_i,cell_index_j,EAST);        
  Preconditioner_dFIdU_Roe(Jacobian[WEST],cell_index_i,cell_index_j,WEST); 

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left) 
  Jacobian[CENTER] += (Jacobian[NORTH] + Jacobian[SOUTH] + Jacobian[EAST]  + Jacobian[WEST])
                     /SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A;

  Jacobian[NORTH] = -Jacobian[NORTH]/SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A;
  Jacobian[SOUTH] = -Jacobian[SOUTH]/SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A;
  Jacobian[EAST] = -Jacobian[EAST]/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A;
  Jacobian[WEST] = -Jacobian[WEST]/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A;

}

/*****************************************************************************
 *  Calculate First Order Local Jacobian Block(s) Coresponding to Cell(i,j)  *
 *  using AUSM_plus_up                                                       *
 ****************************************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
First_Order_Inviscid_Jacobian_AUSM_plus_up(const int &cell_index_i,const int &cell_index_j, 
					   DenseMatrix* Jacobian){              

  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  Preconditioner_dFIdU_AUSM_plus_up(Jacobian[NORTH],cell_index_i,cell_index_j,NORTH);
  Preconditioner_dFIdU_AUSM_plus_up(Jacobian[SOUTH],cell_index_i,cell_index_j,SOUTH); 
  Preconditioner_dFIdU_AUSM_plus_up(Jacobian[EAST],cell_index_i,cell_index_j,EAST);        
  Preconditioner_dFIdU_AUSM_plus_up(Jacobian[WEST],cell_index_i,cell_index_j,WEST); 

  //! Using the fact that dF/dU(right) = - dF/dU(left) 
  Jacobian[CENTER] += (Jacobian[NORTH] + Jacobian[SOUTH] + Jacobian[EAST]  + Jacobian[WEST])
    /SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A;

  Jacobian[NORTH] = -Jacobian[NORTH]/SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A;
  Jacobian[SOUTH] = -Jacobian[SOUTH]/SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A;
  Jacobian[EAST] = -Jacobian[EAST]/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A;
  Jacobian[WEST] = -Jacobian[WEST]/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A;
  
}

/****************************************************************************
 *  Calculate Second Order Local Jacobian Block(s) Coresponding to Cell(i,j) *                      
 ****************************************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Second_Order_Viscous_Jacobian(const int &cell_index_i,const int &cell_index_j, DenseMatrix* Jacobian){

  // A real cludge with all the DenseMatrices and recalculations, 
  //  but just to test, need to change for performance....
  DenseMatrix JacobianN(blocksize,blocksize,ZERO);   //TEMP VAR
  DenseMatrix JacobianS(blocksize,blocksize,ZERO);    //TEMP VAR
  DenseMatrix JacobianE(blocksize,blocksize,ZERO);     //TEMP VAR 
  DenseMatrix JacobianW(blocksize,blocksize,ZERO);     //TEMP VAR 

  //Also should rewrite to minimize dR/dU calls and just 
  //call dWdU, but this needs to be done in Preconditioner_dFVdU
  
  //***************** dR(i,j)/dU(i,j) *********************************************/
  //CENTER
  Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j,
		       cell_index_i,cell_index_j,NORTH,CENTER);
  Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j,
		       cell_index_i,cell_index_j,SOUTH,CENTER);
  Preconditioner_dFVdU(JacobianE,cell_index_i,cell_index_j,
		       cell_index_i,cell_index_j,EAST,CENTER);
  Preconditioner_dFVdU(JacobianW,cell_index_i,cell_index_j,
		       cell_index_i,cell_index_j,WEST,CENTER);

  Jacobian[CENTER] += (JacobianN + JacobianS + JacobianE  + JacobianW) 
    /SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A;


  /***************** dR(i,j-1)/dU(i,j) *********************************************/
  //NORTH                           
  JacobianN.zero();
  Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j-1,
		       cell_index_i,cell_index_j,EAST,NORTH); 
  Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j-1,
		       cell_index_i,cell_index_j,NORTH,NORTH); 
  Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j-1,
		       cell_index_i,cell_index_j,WEST,NORTH); 
 
  /***************** dR(i,j+1)/dU(i,j) *********************************************/
  //SOUTH 
  JacobianS.zero();
  Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j+1,
		       cell_index_i,cell_index_j,EAST,SOUTH);
  Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j+1,
		       cell_index_i,cell_index_j,SOUTH,SOUTH);
  Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j+1,
		       cell_index_i,cell_index_j,WEST,SOUTH);
  
  /***************** dR(i-1,j)/dU(i,j) *********************************************/
  //EAST 
  JacobianE.zero();
  Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j,
		       cell_index_i,cell_index_j,NORTH,EAST);
  Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j,
		       cell_index_i,cell_index_j,EAST,EAST);
  Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j,
		       cell_index_i,cell_index_j,SOUTH,EAST);

  /***************** dR(i+1,j)/dU(i,j) *********************************************/
  //WEST
  JacobianW.zero(); 
  Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j,
		       cell_index_i,cell_index_j, NORTH,WEST);
  Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j,
		       cell_index_i,cell_index_j, WEST,WEST);
  Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j,
		       cell_index_i,cell_index_j, SOUTH,WEST);

  Jacobian[NORTH] += JacobianN/SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A;
  Jacobian[SOUTH] += JacobianS/SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A;
  Jacobian[EAST] += JacobianE/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A;
  Jacobian[WEST] += JacobianW/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A; 

  /********************************************************************************/
  //CORNERS reuse matrices
  JacobianN.zero();  JacobianS.zero();
  JacobianE.zero();  JacobianW.zero();

  /***************** dR(i+1,j+1)/dU(i,j) *********************************************/
  //SOUTHWEST
  Preconditioner_dFVdU(JacobianS,cell_index_i+1,cell_index_j+1,
		       cell_index_i,cell_index_j,SOUTH,SOUTH_WEST);
  Preconditioner_dFVdU(JacobianS,cell_index_i+1,cell_index_j+1,
		       cell_index_i,cell_index_j, WEST,SOUTH_WEST);  
  Jacobian[SOUTH_WEST] += JacobianS/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j+1].A;
  
  /***************** dR(i-1,j+1)/dU(i,j) *********************************************/
  //SOUTHEAST
  Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j+1,
		       cell_index_i,cell_index_j,SOUTH,SOUTH_EAST);
  Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j+1,
		       cell_index_i,cell_index_j,EAST,SOUTH_EAST);  
  Jacobian[SOUTH_EAST] += JacobianE/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j+1].A;

  /***************** dR(i+1,j-1)/dU(i,j) *********************************************/
  //NORTHWEST
  Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j-1,
		       cell_index_i,cell_index_j, NORTH,NORTH_WEST);
  Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j-1,
		       cell_index_i,cell_index_j, WEST,NORTH_WEST);  
  Jacobian[NORTH_WEST] += JacobianW/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j-1].A;

  /***************** dR(i-1,j-1)/dU(i,j) *********************************************/
  //NORTHEAST
  Preconditioner_dFVdU(JacobianN,cell_index_i-1,cell_index_j-1,
		       cell_index_i,cell_index_j, NORTH,NORTH_EAST);
  Preconditioner_dFVdU(JacobianN,cell_index_i-1,cell_index_j-1,
		       cell_index_i,cell_index_j, EAST,NORTH_EAST);  
  Jacobian[NORTH_EAST] += JacobianN/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j-1].A;


}

/*********************************************************
 *  Setup preconditioner from Jacobian Approx            *
 *********************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Setup_Preconditioner()
{
  deallocate_Precon(); // cludge of a workaround to fix up how BpPrecon uses memory in .setup
 
  ILUK_Precon = new BILUK;         //TEMP VAR          
  Jacobi_Precon = new BJacobi;     //TEMP VAR
  
  if (Input_Parameters->NKS_IP.GMRES_Block_Preconditioner == Block_ILUK ) {      
    /* ILU(k) */   
    ILUK_Precon->localprecon(LP_INVERSE);  // currently using "exact" local inverse, many other options which may be cheaper !!!!!
    ILUK_Precon->setup(Block_Jacobian_approx, Input_Parameters->NKS_IP.GMRES_ILUK_Level_of_Fill);
  } else if (Input_Parameters->NKS_IP.GMRES_Block_Preconditioner == Block_Jacobi) {
    /* Diagonal */    
    Jacobi_Precon->localprecon(LP_INVERSE);
    Jacobi_Precon->setup(Block_Jacobian_approx);
  }

}

/*********************************************************
 * Apply preconditioner - in GMRES                       *
 *********************************************************/  //FYI Issues with passing by address &
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline void Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Apply_Preconditioner(int A, int B, const double *C, int D, double *E, int F)
{
  /* z = Minv * V(i) -> stored in W(i). */
  if (Input_Parameters->NKS_IP.GMRES_Block_Preconditioner == Block_ILUK ) {
    ILUK_Precon->apply(A, B, C, D , E, F);
  } else if (Input_Parameters->NKS_IP.GMRES_Block_Preconditioner == Block_Jacobi ) {
    Jacobi_Precon->apply(A, B, C, D , E, F);
  }
}

/********************************************************
 * Block_Preconditioner::DenseMat_to_DenseMatrix.       *
 ********************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline DenseMat Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
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
 * Routine: Rotation_Matrix                              *
 *                                                       *
 * This function returns either the rotation matrix, A,  *
 * or the inverse of A.                                  *
 *                                                       *
 * Note: A_matrix = 1 for returning A.                   *
 *                = 0 for returning inverse of A.        *
 *                                                       *
 *********************************************************/
template <typename SOLN_VAR_TYPE, typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
inline DenseMatrix Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>::
Rotation_Matrix(const Vector2D &nface, const int &A_matrix) 
{
  double cos_angle = nface.x; 
  double sin_angle = nface.y;
    
  DenseMatrix mat(blocksize,blocksize);
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

