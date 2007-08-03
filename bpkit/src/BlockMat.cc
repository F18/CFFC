//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include <new>

using namespace std;

#include "BlockMat.h"
#include "BlockVec.h"
//#include "DenseMat.h"
#include "CSRMat.h"
#include "BPKIT.h"

#include "HBTMat.h"

BlockMat::~BlockMat()
{
    if (a != NULL)
    {
        if (contig_memory)
        {
            delete [] a[0]->Data();       // values
            for (int i=0; i<nnz; i++)
                a[i]->Data() = NULL;      // so DenseMat destructor safe
            delete [] (DenseMat *) a[0];  // must be array of DenseMat
        }
        else
        {
            delete [] (CSRMat *) a[0];
#if 0
            for (int i=0; i<nnz; i++)
                delete a[i];          // delete each matrix
#endif
        }
    }
    delete [] a;
    delete [] ja;
    delete [] ia;
    delete [] kvstr;
    delete [] kvstc;

}

void BlockMat::mult(int nr, int nc, const double *B, int ldu, 
    double *C, int ldv) const
{
    BlockVec Bb(nr, nc, B, ldu, kvstc);
    BlockVec Cb(dimrow(), nc, C, ldv, kvstr);
    Cb.VecSetToZero();

    for (int i=0; i<nrow; i++)
    {
        for (int j=ia[i]; j<ia[i+1]; j++)
	{
             a[j]->Mat_Vec_Mult(Bb(ja[j]), Cb(i), 1.0, 1.0);
	}
    }
}

void BlockMat::trans_mult(int nr, int nc, const double *B, int ldu,
    double *C, int ldv) const
{
    BlockVec Bb(nr, nc, B, ldu, kvstc);
    BlockVec Cb(dimcol(), nc, C, ldv, kvstr);
    Cb.VecSetToZero();

    for (int i=0; i<nrow; i++)
    {
        for (int j=ia[i]; j<ia[i+1]; j++)
	{
             a[j]->Mat_Trans_Vec_Mult(Bb(i), Cb(ja[j]), 1.0, 1.0);
	}
    }
}

// helping function that assumes partitionings already set
// assumes already set: nrow, ncol, kvstr, kvstc
// Storage for dense block matrices will be contiguous across the blocks.

void BlockMat::construct_Dense(const HBTMat& A)
{
    int i, j, k, ii, jj;

    // inverse kvst_col
    int *ikc = new int[A.dimcol()];
    for (i=0; i<ncol; i++)
        for (j=kvst_col(i); j<kvst_col(i+1); j++)
            ikc[j] = i;

    // count number of block nonzeros and number of scalar stored entries

    int *iwk = new int[ncol];
    for (i=0; i<ncol; i++)
        iwk[i] = 0;
    nnz = 0;
    int nnzs = 0; // number of scalar nonzeros

    for (i=0; i<nrow; i++)
    {
        int neqr = kvst_row(i+1)-kvst_row(i);

        // loop on all the elements in the block row to determine block sparsity
        for (jj=A.row_ptr(kvst_row(i)); jj<A.row_ptr(kvst_row(i+1)); jj++)
        {
            iwk[ikc[A.col_ind(jj)]] = 1;
        }

        // count and set iwk back to zero
        for (j=0; j<ncol; j++)
        {
            if (iwk[j])
            {
                nnz++;
                nnzs = nnzs + neqr * (kvst_col(j+1)-kvst_col(j));
                iwk[j] = 0;
            }
        }
    }
    // cerr << nnz << " block nonzeros.\n";
    // cerr << nnzs << " stored scalar entries.\n";

    // allocate memory
    contig_memory = TRUE;
    ia = new int[nrow+1];
    ja = new int[nnz];
    a  = new LocalMatp[nnz];   // array of pointers to blocks

    DenseMat *p = new DenseMat[nnz]; // array of blocks

    for (k=0; k<nnz; k++)
        a[k] = &p[k];

    double *pp = new double[nnzs];  // array of values

    // fill entries and ja -----------------------------------------------------

    int a0 = 0;
    k = 0;
    ia[0] = 0;

    // loop on block rows
    for (i=0; i<nrow; i++)
    {
        int neqr = kvst_row(i+1)-kvst_row(i);

        // loop on all the elements in the block row to determine block sparsity
        for (jj=A.row_ptr(kvst_row(i)); jj<A.row_ptr(kvst_row(i+1)); jj++)
        {
            iwk[ikc[A.col_ind(jj)]] = 1;
        }

        // count and set iwk back to zero
        for (j=0; j<ncol; j++)
        {
            if (iwk[j])
            {
                ja[k] = j;
                p[k].set(neqr,kvst_col(j+1)-kvst_col(j),pp);
                pp = pp + neqr * (kvst_col(j+1)-kvst_col(j));
                k++;
                iwk[j] = 0;
            }
        }
        ia[i+1] = k;

        // loop on scalar rows in block row
        for (ii=0; ii<neqr; ii++)
        {
            double *b0 = a[ia[i]]->Data() + ii;  
                                             // loc of first nonz in this
                                             // scalar row in the block format
            // loop on block columns
            for (j=ia[i]; j<ia[i+1]; j++)
            {
                // loop on scalar columns within block column
                for (jj=kvst_col(ja[j]); jj<kvst_col(ja[j]+1); jj++)
                {
                    if (a0 >= A.row_ptr(kvst_row(i)+ii+1))
                        *b0 = 0.0;
                    else
                        if (jj == A.col_ind(a0))
                            *b0 = A.val(a0++);
                        else
                            *b0 = 0.0;
                    b0 = b0 + neqr;
                }
            }
        }
    }

    delete [] ikc;
    delete [] iwk;
}


// conversion to blockmat with sparse blocks
// Storage for sparse block matrices will not be contiguous across the blocks.

void BlockMat::construct_Sparse(const HBTMat& A)
{
    int i, j, k, ii, jj;

    // inverse kvst_col
    int *ikc = new int[A.dimcol()];
    for (i=0; i<ncol; i++)
        for (j=kvst_col(i); j<kvst_col(i+1); j++)
            ikc[j] = i;

    // count number of block nonzeros and number of scalar stored entries

    int *iwk = new int[ncol];
    for (i=0; i<ncol; i++)
        iwk[i] = 0;
    nnz = 0;

    for (i=0; i<nrow; i++)
    {
        // loop on all the elements in the block row to determine block sparsity
        for (jj=A.row_ptr(kvst_row(i)); jj<A.row_ptr(kvst_row(i+1)); jj++)
        {
            iwk[ikc[A.col_ind(jj)]] = 1;
        }

        // count and set iwk back to zero
        for (j=0; j<ncol; j++)
        {
            if (iwk[j])
            {
                nnz++;
                iwk[j] = 0;
            }
        }
    }
    // cerr << nnz << " block nonzeros.\n";

    // allocate memory
    contig_memory = FALSE;
    ia = new int[nrow+1];
    ja = new int[nnz];
    a  = new LocalMatp[nnz];   // array of pointers to blocks

    CSRMat *p = new CSRMat[nnz]; // array of blocks

    for (k=0; k<nnz; k++)
        a[k] = &p[k];

    // fill entries and ja -----------------------------------------------------

    k = 0;
    ia[0] = 0;

    // loop on block rows
    for (i=0; i<nrow; i++)
    {
        int neqr = kvst_row(i+1)-kvst_row(i);

        // loop on all the elements in the block row to determine block sparsity
        for (jj=A.row_ptr(kvst_row(i)); jj<A.row_ptr(kvst_row(i+1)); jj++)
        {
            iwk[ikc[A.col_ind(jj)]]++;
        }

        // count and set iwk back to zero
        for (j=0; j<ncol; j++)
        {
            if (iwk[j])
            {
                ja[k] = j;

		double *a_ = new double[iwk[j]];
		int *ja_ = new int[iwk[j]];
		int *ia_ = new int[neqr+1];

                p[k].set(a_, ja_, ia_, neqr, kvst_col(j+1)-kvst_col(j), iwk[j]);

		int kk = 0;
		ia_[0] = 0;
		for (ii=0; ii<neqr; ii++) // scalar rows
		{
                    for (jj=A.row_ptr(kvst_row(i)+ii); 
			 jj<A.row_ptr(kvst_row(i)+ii+1); jj++)
		    {
			if (A.col_ind(jj) >= kvst_col(j) &&
			    A.col_ind(jj) <  kvst_col(j+1))
			{
			    ja_[kk] = A.col_ind(jj) - kvst_col(j);
			     a_[kk++] = A.val(jj);
			}
		    }
		    ia_[ii+1] = kk;
		}

                // cerr << (CSRMat) p[k] << endl << endl;

                k++;
                iwk[j] = 0;
            }
        }
        ia[i+1] = k;
    }

    delete [] ikc;
    delete [] iwk;
}

void BlockMat::setup(const HBTMat& A, int blocksize, BlockType blocktype)
{
    a = (LocalMat **) NULL;
    ja = (int *) NULL;
    ia = (int *) NULL;
    kvstr = NULL;
    kvstc = NULL;
    contig_memory = FALSE;

    int i;

    nrow = A.dimrow() / blocksize;
    if (A.dimrow() % blocksize)
        nrow++;

    ncol = A.dimcol() / blocksize;
    if (A.dimcol() % blocksize)
        ncol++;

    kvstr = new int[nrow+1];
    kvstc = new int[ncol+1];

    for (i=0; i<nrow; i++)
        kvstr[i] = i*blocksize;
    for (i=0; i<ncol; i++)
        kvstc[i] = i*blocksize;

    // set the last one explicitly, in case it is uneven
    kvstr[nrow] = A.dimrow();
    kvstc[ncol] = A.dimcol();

    if (blocktype == CSR)
        construct_Sparse(A);
    else
        construct_Dense(A);
}

BlockMat::BlockMat(const HBTMat& A, int blocksize, BlockType blocktype)
{
    a = (LocalMat **) NULL;
    ja = (int *) NULL;
    ia = (int *) NULL;
    kvstr = NULL;
    kvstc = NULL;
    contig_memory = FALSE;

    int i;

    nrow = A.dimrow() / blocksize;
    if (A.dimrow() % blocksize)
        nrow++;

    ncol = A.dimcol() / blocksize;
    if (A.dimcol() % blocksize)
        ncol++;

    kvstr = new int[nrow+1];
    kvstc = new int[ncol+1];

    for (i=0; i<nrow; i++)
        kvstr[i] = i*blocksize;
    for (i=0; i<ncol; i++)
        kvstc[i] = i*blocksize;

    // set the last one explicitly, in case it is uneven
    kvstr[nrow] = A.dimrow();
    kvstc[ncol] = A.dimcol();

    if (blocktype == CSR)
        construct_Sparse(A);
    else
        construct_Dense(A);
}

BlockMat::BlockMat(const HBTMat& A, int numblocks,
    const int *partit, BlockType blocktype)
{
    a = (LocalMat **) NULL;
    ja = (int *) NULL;
    ia = (int *) NULL;
    kvstr = NULL;
    kvstc = NULL;
    contig_memory = FALSE;

    int scalar_dim = A.dimrow();
    nrow = numblocks;
    ncol = numblocks;

    kvstr = new int[nrow+1];
    kvstc = new int[nrow+1];

    for (int i=0; i<nrow+1; i++)
    {
        kvstr[i] = partit[i];
        kvstc[i] = partit[i];
    }

    if (blocktype == CSR)
        construct_Sparse(A);
    else
        construct_Dense(A);
}

// Non-member functions

ostream& operator << (ostream& os, const BlockMat& mat)
{
        int M = mat.numrow();
        int N = mat.numcol();
        int rowp1, colp1;
        int flag = 0;

#ifdef _GNU_GCC_3
        ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
        ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
#else
#ifdef _GNU_GCC_296
        long olda = os.setf(ios::right,ios::adjustfield);
        long oldf = os.setf(ios::scientific,ios::floatfield);
#else
        long olda = os.setf(ios::right,ios::adjustfield);
        long oldf = os.setf(ios::scientific,ios::floatfield);
#endif
#endif

        int oldp = os.precision(12);

        for (int i = 0; i < M ; i++)
           for (int j=mat.row_ptr(i);j<mat.row_ptr(i+1);j++)
           {   
              rowp1 =  i ;
              colp1 =  mat.col_ind(j) ;
	      os <<" Block ("<< rowp1<<","<<colp1<<")"<<endl;
#if 1
              //os.width(20); 
              os << (DenseMat &) mat.val(j); // *(mat.val(j).Data());
#endif
              os << endl;
           }

#if 0
        if (flag == 0)
        {
           os.width(14);
           os <<  M ; os << "    " ;
           os.width(14);
           os <<  N ; os << "    " ;
           os.width(20);
           os <<  mat(M-1,N-1) << "\n";
        }
#endif

        os.setf(olda,ios::adjustfield);
        os.setf(oldf,ios::floatfield);
        os.precision(oldp);

        return os;
}

// Convert a matrix in block coordinate format to the BlockMat data structure.
// nrow = number of block rows
// nnz = number of block nonzeros
// note that the rows in the resultant data structure
// may not be sorted; they are in the order in which they were given

BlockMat::BlockMat(int nrow_, int nnz_, int *row, int *col, double *A, 
  int blocksize)
{
    a = (LocalMat **) NULL;
    ja = (int *) NULL;
    ia = (int *) NULL;
    kvstr = NULL;
    kvstc = NULL;
    contig_memory = FALSE;

    int i;

    nrow = nrow_;
    ncol = nrow_;
    nnz  = nnz_;

    kvstr = new int[nrow+1];
    kvstc = new int[ncol+1];

    for (i=0; i<=nrow; i++)
        kvstr[i] = i*blocksize;
    for (i=0; i<=ncol; i++)
        kvstc[i] = i*blocksize;

    // allocate memory
    contig_memory = TRUE;
    ia = new int[nrow+1];
    ja = new int[nnz];
    a  = new LocalMatp[nnz];   // array of pointers to blocks

    DenseMat *p = new DenseMat[nnz]; // array of blocks

    for (i=0; i<nnz; i++)
        a[i] = &p[i];

    double *pp = new double[nnz*blocksize*blocksize];  // array of values

    // construct ia ja structure

    // count the nonzeros in each row and store in ia array
    for (i=0; i<=nrow; i++)
        ia[i] = 0;

    for (i=0; i<nnz; i++)
        ia[row[i]-1]++;

    // put starting position of each row in ia array
    int prev, temp;
    prev = 0;
    for (i=0; i<=nrow; i++)
    {
        temp = ia[i];
        ia[i] = prev;
        prev = prev + temp;
    }

    // fill output matrix
    double *b;
    int j;
    for (i=0; i<nnz; i++)
    {
        int irow = row[i]-1;
        ja[ia[irow]] = col[i] - 1; // 0-based

        p[ia[irow]].set(blocksize,blocksize,pp);
	pp += blocksize*blocksize;

        b = a[ia[irow]]->Data();
        for (j=0; j<blocksize*blocksize; j++)
            *b++ = *A++;

        ia[irow]++;
    }

    // shift back ia
    for (i=nrow-1; i>=0; i--)
        ia[i+1] = ia[i];
    ia[0] = 0;

#if 0
for (i=0; i<=nrow; i++)
   printf("ia(%d) = %d\n", i, ia[i]);
for (i=0; i<ia[nrow]; i++)
   printf("ja(%d) = %d\n", i, ja[i]);
for (i=0; i<nnz; i++)
   printf("a(%d) = %f\n", i, *(a[i]->Data()));
fflush(NULL);
#endif
}

// Convert a matrix in block coordinate format to the BlockMat data structure.
// nrow = number of block rows
// nnz = number of block nonzeros
// note that the rows in the resultant data structure
// may not be sorted; they are in the order in which they were given

void BlockMat::setup(int nrow_, int nnz_, int *row, int *col, double *A, 
  int blocksize)
{
    a = (LocalMat **) NULL;
    ja = (int *) NULL;
    ia = (int *) NULL;
    kvstr = NULL;
    kvstc = NULL;
    contig_memory = FALSE;

    int i;

    nrow = nrow_;
    ncol = nrow_;
    nnz  = nnz_;

    kvstr = new int[nrow+1];
    kvstc = new int[ncol+1];

    for (i=0; i<=nrow; i++)
        kvstr[i] = i*blocksize;
    for (i=0; i<=ncol; i++)
        kvstc[i] = i*blocksize;

    // allocate memory
    contig_memory = TRUE;
    ia = new int[nrow+1];
    ja = new int[nnz];
    a  = new LocalMatp[nnz];   // array of pointers to blocks

    DenseMat *p = new DenseMat[nnz]; // array of DENSE blocks  

    for (i=0; i<nnz; i++)
        a[i] = &p[i];

    double *pp = new double[nnz*blocksize*blocksize];  // array of values

    // construct ia ja structure

    // count the nonzeros in each row and store in ia array
    for (i=0; i<=nrow; i++)
        ia[i] = 0;

    for (i=0; i<nnz; i++)
        ia[row[i]-1]++;

    // put starting position of each row in ia array
    int prev, temp;
    prev = 0;
    for (i=0; i<=nrow; i++)
    {
        temp = ia[i];
        ia[i] = prev;
        prev = prev + temp;
    }

    // fill output matrix
    double *b;
    int j;
    for (i=0; i<nnz; i++)
    {
        int irow = row[i]-1;
        ja[ia[irow]] = col[i] - 1; // 0-based

        p[ia[irow]].set(blocksize,blocksize,pp);
	pp += blocksize*blocksize;

        b = a[ia[irow]]->Data();
        for (j=0; j<blocksize*blocksize; j++)
            *b++ = *A++;

        ia[irow]++;
    }

    // shift back ia
    for (i=nrow-1; i>=0; i--)
        ia[i+1] = ia[i];
    ia[0] = 0;

#if 0
for (i=0; i<=nrow; i++)
   printf("ia(%d) = %d\n", i, ia[i]);
for (i=0; i<ia[nrow]; i++)
   printf("ja(%d) = %d\n", i, ja[i]);
for (i=0; i<nnz; i++)
   printf("a(%d) = %f\n", i, *(a[i]->Data()));
fflush(NULL);
#endif
}



/********************************************************
 * BlockMat::setblock -- set the specified block.    *
 ********************************************************/
void BlockMat::setblock(const int BLK_i,const int BLK_j,const DenseMat &A) {

   int index = -99;
   //find equivalent Localmat location in block CSR form (SIMPLER/BETTER WAY ???)
   for (int t = ia[BLK_i]; t < ia[BLK_i+1]; t++) {
     if (ja[t] == BLK_j) {
       index = t;
     }
   } 
   //Copy A into BlockMat data structure at (I,J)
   a[index][0].MatCopy(A);

} 
