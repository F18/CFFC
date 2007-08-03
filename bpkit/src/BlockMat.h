//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _BLOCKMAT_H_
#define _BLOCKMAT_H_

#include <iostream>

using namespace std;

#include "BpMatrix.h"
#include "DenseMat.h"

class LocalMat;
class HBTMat;

enum BlockType
{
    DENSE = 0,
    CSR   = 1
};

class BlockMat : public BpMatrix
{
 private:
  LocalMat **a;
  int *ja;
  int *ia;
  
  int nrow;
  int ncol;
  int nnz;
  
  int *kvstr;
  int *kvstc;
  
  int contig_memory;
  
  // helping functions
  void construct_Dense (const HBTMat& A);
  void construct_Sparse(const HBTMat& A);

 public: 
  /* Added BlockMat() by K. Tsang on March 27, 2003. */
  BlockMat() {a = NULL; ja = NULL; ia = NULL;
              nrow = 0; ncol = 0; nnz  = 0;
              kvstr = NULL; kvstc = NULL; contig_memory = 0;
	      //cout << " End : BlockMat constructor......" << endl;
             }
  BlockMat(const HBTMat&, int, BlockType);
  BlockMat(const HBTMat&, int, const int *, BlockType);
  BlockMat(int nrow, int nnz, int *row, int *col, double *A, int blocksize);
  ~BlockMat();
  
  void setup(const HBTMat&, int, BlockType);

  const LocalMat& val(unsigned int i) const {return *a[i];}
  LocalMat& val(unsigned int i) {return *a[i];}
  const int& row_ptr(unsigned int i) const {return ia[i];}
  const int& col_ind(unsigned int i) const {return ja[i];}

  //added by SN, probably temporary !
  void setblock(const int, const int, const DenseMat&);
  void setup(int nrow, int nnz, int *row, int *col, double *A, int blocksize);


  int numrow() const {return nrow;}
  int numcol() const {return ncol;}
  int numnz()  const {return nnz;}
  
  int dimrow() const {return kvst_row(nrow);} // scalar dimension
  int dimcol() const {return kvst_col(ncol);}
  
  const int& kvst_row(unsigned int i) const {return kvstr[i];}
  const int& kvst_col(unsigned int i) const {return kvstc[i];}
  
  void mult(int, int, const double *, int, double *, int) const;
  void trans_mult(int, int, const double *, int, double *, int) const;
};

ostream& operator << (ostream& os, const BlockMat& mat);

#endif // _BLOCKMAT_H_
