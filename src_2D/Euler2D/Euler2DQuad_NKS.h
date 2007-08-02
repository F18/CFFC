#ifndef _EULER2D_QUAD_NKS_INCLUDED
#define _EULER2D_QUAD_NKS_INCLUDED

#ifndef _EULER2D_QUAD_INCLUDED
#include "Euler2DQuad.h"
#endif // _EULER2D_QUAD_INCLUDED

#ifndef _EULER2D_QUAD_GMRES_INCLUDED
#include "Euler2DQuad_GMRES.h"
#endif // _EULER2D_QUAD_GMRES_INCLUDED

/********************************************************
 * Class: Index                                         *
 *                                                      *
 * Member functions                                     * 
 *       i      -- Return i index in the block matrix.  *
 *       j      -- Return j index in the block matrix.  *
 *                                                      *
 ********************************************************/
class Index {
 public:
  int x;
  int y;

  Index() {int x = -1; int y = -1;}

  /* Destructor. */
  // ~index(void);
  // Use automatically generated destructor.

}; /* End of Index class. */

/********************************************************
 * Class: HBTMat_Class                                  *
 *                                                      *
 * Member functions                                     * 
 *       deallocate_mem -- Deallocate memory for        *
 *                         all member variables.        *
 *                                                      *
 ********************************************************/
class HBTMat_Class {
 public:
  int *ncolv;  
  int *nrowv;  
  int *rowind;
  
  double *valv;
  double *rhs;  
  double *guess;
  double *exact;

  int blocksize;
  int fill_pattern;

  HBTMat_Class(int BLOCKSIZE, int FILL_PATTERN) {
    ncolv = NULL; nrowv = NULL; rowind = NULL;
    valv = NULL; rhs = NULL; guess = NULL; exact = NULL;
    blocksize = BLOCKSIZE; fill_pattern = FILL_PATTERN; }

  /* Destructor. */
  // ~HBTMat_Class(void);
  // Use automatically generated destructor.

  void deallocate_mem() {
    delete [] ncolv; 
    delete [] nrowv; 
    delete [] rowind;
    delete [] valv;  
    delete [] rhs;   
    delete [] guess;  
    delete [] exact;
  }
};

/********************************************************
 * Class: BlockMatrix                                   *
 *                                                      *
 * Member functions                                     * 
 *      create  -- creates block matrix object          *
 *    getvalue  -- Return a value in a specified block. *
 *    setvalue  -- Re-set a value in a specified block. *
 *    getblock  -- Return a specified block in Dense    *
 *                 format.                              *
 *    setblock  -- Re-set a specified block.            *
 *    get_nx    -- Return the number of points in x-dir.*
 *    get_ny    -- Return the number of points in y-dir.*
 *    get_nz    -- Return the number of points in z-dir.* 
 *    get_blocksize    -- Return the number of          *
 *                        variables.                    *
 *    get_fill_pattern -- Return fill pattern number.   * 
 *    Mindex_2D -- Return i and j indices of the block  * 
 *                 matrix given cell location, i and j. * 
 *                 Note : returning_vector[0] = i index *
 *                        returning_vector[1] = j index *
 *    Mindex_3D -- Return i and j indices of the block  * 
 *                 matrix given cell location, i, j     *
 *                 and k.                               * 
 *    Gindex_2D -- Return i and j indices of the grid   * 
 *                 given cell location, i and j.        * 
 *                 Note : returning_vector[0] = i index *
 *                        returning_vector[1] = j index *
 *    DenseMat_to_DenseMatrix -- Return dense matrix in *
 *                               DenseMatrix format.    *    
 *    DenseMatrix_to_DenseMat -- Return dense matrix in *
 *                               DenseMat format.       *  
 *                                                      *
 ********************************************************/
class BlockMatrix : public BlockMat
{
 private: 
 public:
  int nx, ny, nz;
  int blocksize;
  int fill_pattern;

  BlockMatrix() { nx = ny = nz = 0; blocksize = 0; fill_pattern = 0; }
  BlockMatrix(int nx, int ny, int nz, HBTMat &model, 
	      int blocksize, int fill_pattern);
  ~BlockMatrix() {}
  
  void create(int nx, int ny, int nz, HBTMat &model, int blocksize, int fill_pattern);

  double getvalue(int BLK_i, int BLK_j, int LOCAL_i, int LOCAL_j);
  void   setvalue(int BLK_i, int BLK_j, int LOCAL_i, int LOCAL_j, double value);

  DenseMatrix getblock(int BLK_i, int BLK_j, int check = 1);
  void setblock(int BLK_i, int BLK_j, DenseMatrix &A);

  int get_nx() {return nx;}
  int get_ny() {return ny;}
  int get_nz() {return nz;}

  int get_blocksize()    {return blocksize;}
  int get_fill_pattern() {return fill_pattern;}

  Index Mindex_2D(int base_cell_i, int base_cell_j, int location);
  Index Mindex_3D(int base_cell_i, int base_cell_j, int base_cell_k, int location);

  friend Index Gindex_2D(int base_cell_i, int base_cell_j, int location);

  friend DenseMatrix DenseMat_to_DenseMatrix(const DenseMat    &A);
  friend DenseMat    DenseMatrix_to_DenseMat(const DenseMatrix &A);

}; /* End of BlockMatrix class. */

/********************************************************
 * BlockMatrix::BlockMatrix -- Constructor.             *
 ********************************************************/
inline BlockMatrix::BlockMatrix(int nx, 
                                int ny, 
                                int nz, 
                                HBTMat &model, 
                                int blocksize, 
                                int fill_pattern) : BlockMat(model, blocksize, DENSE)
{
  this->nx = nx;
  this->ny = ny;
  this->nz = nz;
 
  this->blocksize    = blocksize;
  this->fill_pattern = fill_pattern;
  
} /* End of BlockMatrix::BlockMatrix. */

/********************************************************
 * BlockMatrix::create -- creates block matrix object   *
 ********************************************************/
inline void BlockMatrix::create(int nx, 
                                int ny, 
                                int nz, 
                                HBTMat &model,
			        int blocksize, 
                                int fill_pattern) {
  
  setup(model, blocksize, DENSE);
  this->nx = nx;
  this->ny = ny;
  this->nz = nz;
  
  this->blocksize    = blocksize;
  this->fill_pattern = fill_pattern;
  
} /* End of setup. */

/********************************************************
 * BlockMatrix::getvalue -- a value in the specified    *
 *                          block.                      *
 ********************************************************/
inline double BlockMatrix::getvalue(int BLK_i, 
                                    int BLK_j, 
                                    int LOCAL_i, 
                                    int LOCAL_j) {

  DenseMatrix densematrix(blocksize,blocksize,0.0);
  densematrix = getblock(BLK_i, BLK_j, 0);

  if (OFF) {
     cout << "BlockMatrix::getvalue at -> Local("<<LOCAL_i<<","
                                                 <<LOCAL_j<<") = "
                                                 <<densematrix(LOCAL_i, LOCAL_j)
                                                 <<endl;
  }

  return densematrix(LOCAL_i, LOCAL_j);

} /* End of BlockMatrix::getvalue. */ 

/********************************************************
 * BlockMatrix::setvalue -- a value in the specified    *
 *                          block.                      *
 ********************************************************/
inline void  BlockMatrix::setvalue(int BLK_i, 
                                   int BLK_j, 
                                   int LOCAL_i, 
                                   int LOCAL_j, 
                                   double value) {

  DenseMat *dmat;
  int index = -99;
  
  dmat = (DenseMat *) &(val(0));
  for (int t=row_ptr(BLK_i); t<row_ptr(BLK_i+1); t++) {
    if (col_ind(t) == BLK_j) {
      index = t;
    }
  } 
  
  if (index < ZERO) { 
    cout << "BlockMatrix::index_2D - unable to set the corresponding value in the specified block of the Jacobian matrix!! " << endl;
  } else {
    dmat[index](LOCAL_i, LOCAL_j) = value;
    if (OFF) {
    cout << "\nBlockMatrix::setvalue at -> Block("<<BLK_i<<","<<BLK_j<<")"
	 << ", Non-Zero Block #" <<index+1<<" of "<< numnz() << endl;
    cout << "                            Local("<<LOCAL_i<<","
	                                        <<LOCAL_j<<")"<< endl;
    }
  }
  
} /* End of BlockMatrix::setvalue. */

/********************************************************
 * BlockMatrix::getblock -- get the specified block.    *
 ********************************************************/
inline DenseMatrix BlockMatrix::getblock(int BLK_i, 
                                         int BLK_j, 
                                         int check) {
    
  double *value = NULL;
  DenseMatrix densematrix(blocksize, blocksize, 0.0);
  
  DenseMat *dmat;
  dmat = (DenseMat *) &(val(0));
  
  for (int t=row_ptr(BLK_i); t<row_ptr(BLK_i+1); t++) {
    if (col_ind(t) == BLK_j) {
      value = dmat[t].Data();
      for (int i=0;i<blocksize*blocksize;i++) {
	densematrix(i%blocksize,int(i/blocksize)) = value[i];
      }
      cout << "\nBlockMatrix::getblock at -> Block("<<BLK_i<<","<<BLK_j<<")"
	   << ", Non-Zero Block #" <<t+1<<" of "<< numnz() << endl;      
      if (check) { 
	cout << densematrix << endl;
      }
      return densematrix;

    }
  } 
  cout << "BlockMatrix::index_2D - unable to return the corresponding block from the Jacobian matrix!! " << endl;
  
  return densematrix;

} /* End of BlockMatrix::getblock. */

/********************************************************
 * BlockMatrix::setblock -- set the specified block.    *
 ********************************************************/
inline void BlockMatrix::setblock(int BLK_i, 
                                  int BLK_j, 
                                  DenseMatrix &A) {

  DenseMat *dmat;
  int index = -99;

  dmat = (DenseMat *) &(val(0));
  for (int t=row_ptr(BLK_i); t<row_ptr(BLK_i+1); t++) {
    if (col_ind(t) == BLK_j) {
      index = t;
    }
  } 

  if (index < ZERO) { 
    cout << "BlockMatrix::index_2D - unable to set the corresponding block in the Jacobian matrix!! " << endl;
  } else {
    for (int i=0; i<blocksize; i++) {
      for( int j=0; j<blocksize; j++) {
	dmat[index](i,j) = A(i,j);
      }
    }
    if (OFF) {
       cout << "\nBlockMatrix::setblock at -> Block("<<BLK_i<<","
	    <<BLK_j<<")" << ", non-zero block #" <<index+1<<" of "
	    << numnz() << endl;
    }
  }

} /* End of BlockMatrix::setblock. */

/********************************************************
 * BlockMatrix::Mindex_2D -- find indices in the block  * 
 *                           matrix given cell location * 
 *                           (i, j) and location.       *
 ********************************************************/
inline Index BlockMatrix::Mindex_2D(int base_cell_i, 
                                    int base_cell_j, 
                                    int location)
{
  Index index_ij;
  index_ij.x = 999;
  index_ij.y = 999;
  int index_i = 0;
  int index_j = 0;

  index_ij.x = base_cell_j * nx + base_cell_i;

  if (fill_pattern == 1) {
    switch (location) {
    case CENTER:
      index_i = base_cell_i;
      index_j = base_cell_j;
      index_ij.y = index_j * nx + index_i;
      break;
    case NORTH:
      index_i = base_cell_i;
      index_j = base_cell_j + 1;
      index_ij.y = index_j * nx + index_i;
      break;
    case SOUTH:
      index_i = base_cell_i;
      index_j = base_cell_j - 1;
      index_ij.y = index_j * nx + index_i;
      break;
    case EAST:
      index_i = base_cell_i + 1;
      index_j = base_cell_j;
      index_ij.y = index_j * nx + index_i;
      break;
    case WEST:
      index_i = base_cell_i - 1;
      index_j = base_cell_j;
      index_ij.y = index_j * nx + index_i;
      break;
    default:
      cout << "BlockMatrix::index_2D - unable to return the corresponding i and j indices in the Jacobian matrix!! " << endl;
    }
  } else {
    switch (location) {
    case CENTER:
      index_i = base_cell_i;
      index_j = base_cell_j;
      index_ij.y = index_j * nx + index_i;
      break;
    case NORTH:
      index_i = base_cell_i;
      index_j = base_cell_j + 1;
      index_ij.y = index_j * nx + index_i;
      break;
    case SOUTH:
      index_i = base_cell_i;
      index_j = base_cell_j - 1;
      index_ij.y = index_j * nx + index_i;
      break;
    case EAST:
      index_i = base_cell_i + 1;
      index_j = base_cell_j;
      index_ij.y = index_j * nx + index_i;
      break;
    case WEST:
      index_i = base_cell_i - 1;
      index_j = base_cell_j;
      index_ij.y = index_j * nx + index_i;
      break;
    case SNORTH:
      index_i = base_cell_i;
      index_j = base_cell_j + 2;
      index_ij.y = index_j * nx + index_i;
      break;
    case SSOUTH:
      index_i = base_cell_i;
      index_j = base_cell_j - 2;
      index_ij.y = index_j * nx + index_i;
      break;
    case SEAST:
      index_i = base_cell_i + 2;
      index_j = base_cell_j;
      index_ij.y = index_j * nx + index_i;
      break;
    case SWEST:
      index_i = base_cell_i - 2;
      index_j = base_cell_j;
      index_ij.y = index_j * nx + index_i;
      break;
    default:
      cout << "BlockMatrix::index_2D - unable to return the corresponding i and j indices in the Jacobian matrix!! " << endl;
    }
  }

  for (int t=row_ptr(int(index_ij.x)); t<row_ptr(int(index_ij.x+1)); t++) {
    if (col_ind(t) == index_ij.y) {
      if (OFF) {
	cout << " Mindex_2D.x = " << index_ij.x 
	     << " , Mindex_2D.y = " << index_ij.y << endl;   
      }
    return index_ij;
    }
  }
  
  index_ij.x = -1;
  index_ij.y = -1;
  
  return index_ij;

}; /* End of BlockMatrix::index_2D. */

/********************************************************
 * BlockMatrix::Mindex_3D -- find indices in the block  * 
 *                           matrix given cell location * 
 *                           (i, j, k) and location.    *
 ********************************************************/
inline Index BlockMatrix::Mindex_3D(int base_cell_i, 
                                    int base_cell_j, 
                                    int base_cell_k, 
                                    int location) {
  Index index_ijk;

  index_ijk.x = base_cell_i;
  index_ijk.y = base_cell_j+base_cell_k;

  return index_ijk;

} /* End of BlockMatrix::Mindex_3D. */

/********************************************************
 * BlockMatrix::Gindex_2D -- find indices in the grid   * 
 *                           given cell location (i, j) * 
 *                           and location.              *
 ********************************************************/
inline Index Gindex_2D(int base_cell_i, 
                       int base_cell_j, 
                       int location) {
  Index index_ij;
  index_ij.x = 99999;
  index_ij.y = 99999;

  switch(location) { 
  case CENTER:
    index_ij.x = base_cell_i;
    index_ij.y = base_cell_j;
    break;
  case NORTH:
    index_ij.x = base_cell_i;
    index_ij.y = base_cell_j+1;
    break;
  case SOUTH:
    index_ij.x = base_cell_i;
    index_ij.y = base_cell_j-1;
    break;
  case EAST:
    index_ij.x = base_cell_i+1;
    index_ij.y = base_cell_j;
    break;
  case WEST:
    index_ij.x = base_cell_i-1;
    index_ij.y = base_cell_j;
    break;
  default:
    cout << "Problem :  BlockMatrix::Gindex_2D......" << endl;
    break;
  }

  return index_ij;

} /* End of BlockMatrix::Gindex_2D. */

/********************************************************
 * BlockMatrix::DenseMat_to_DenseMatrix.                *
 ********************************************************/
inline DenseMatrix DenseMat_to_DenseMatrix(const DenseMat &A) {
  DenseMatrix dm(A.numrow(), A.numcol(), 0.0);

  for (int i=0; i<A.numrow(); i++) {
    for( int j=0; j<A.numcol(); j++) {
      dm(i,j) = A(i,j);
    }
  }

  return(dm);

} /* End of BlockMatrix::DenseMat_to_DenseMatrix. */

/********************************************************
 * BlockMatrix::DenseMatrix_to_DenseMat.                *
 ********************************************************/
inline DenseMat DenseMatrix_to_DenseMat(const DenseMatrix &A) {
  DenseMat dm(A.dim(0), A.dim(1));

  for (int i=0; i<A.dim(0); i++) {
    for( int j=0; j<A.dim(1); j++) {
      dm(i,j) = A(i,j);
    }
  }
 
  return(dm);
  
} /* End of BlockMatrix::DenseMatrix_to_DenseMat. */

/**************************************************************************
 * Euler2D_Quad_NKS -- External Subroutines.                              *
 **************************************************************************/

extern int Newton_Krylov_Schwarz_Solver(CPUTime &NKS_time,
					ostream &Progress_File,
                                        int &Number_of_Startup_Interations,
					Euler2D_Quad_Block *Soln_ptr,
					AdaptiveBlock2D_List &Soln_Block_List,
					Euler2D_Input_Parameters &Input_Parameters);

extern DenseMatrix Jacobian_LocalBlock(int x_pt, 
                                       int y_pt, 
                                       Euler2D_Quad_Block &SolnBlk, 
                                       int location, 
                                       Index M_ind, 
                                       int blocksize, 
                                       int normalization, 
                                       double dTime = 0.0);

extern HBTMat Create_Block_Matrix(int nx, int ny, int nz, HBTMat_Class &HC);

extern HBTMat Create_BlockMat_Object_5 (int xpts, int ypts, HBTMat_Class &HC);

extern HBTMat Create_BlockMat_Object_13(int xpts, int ypts, HBTMat_Class &HC);

extern void normalize_Jacobian(DenseMatrix &dFdU);

#endif  /* _EULER2D_QUAD_NKS_INCLUDED */
