//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>

using namespace std;

#include "HBTMat.h"
#include "BPKIT.h"
#include "SparseUtil.h"
#include "BpResource.h"

#ifndef RAND_MAX
#define RAND_MAX        32767
#endif

// load resources
int HBTMat::rhs_resource = BpGetInt(HBTMAT_RHS, HBTMAT_RHS_RAND);
int HBTMat::transpose_resource = BpGetInt(HBTMAT_TRANSPOSE, 0);
int HBTMat::scale_resource = BpGetInt(HBTMAT_SCALE, 0);

HBTMat::~HBTMat()
{
    if (owndata == TRUE)
    {
        delete [] a;
        delete [] ja;
        delete [] ia;
        delete [] rhs;
        delete [] guess;
        delete [] exact;
    }
}

void HBTMat::FreeMat()
{
    if (owndata == TRUE)
    {
        delete [] a;
        delete [] ja;
        delete [] ia;
        a = NULL;
        ja = NULL;
        ia = NULL;
    }
}

HBTMat::HBTMat(const char *filename)
{
    owndata = TRUE;

    FILE *infile = fopen(filename, "r");
    if (infile == NULL)
    {
       cerr << filename << endl;
       bperror("HBTMat: cannot open file", 0);
    }

    char line[82];
    char Type[4];
    int dummy, Rhscrd, Nrhs;
    double *rhsdata;

    fgets(line, 82, infile);
    sscanf(line, "%72c %8c", title, key);
    fgets(line, 82, infile);
    sscanf(line, "%i %i %i %i %i", &dummy, &dummy, &dummy, &dummy, &Rhscrd);

    fgets(line, 82, infile);
    sscanf(line, "%3c %i %i %i %i", Type, &nrow, &ncol, &nnz, &dummy);

    // check this is an RUA file
    if (strncmp(Type, "RUA", 3) != 0 && strncmp(Type, "rua", 3) != 0)
    {
	cerr << filename << endl;
        bperror("HBTMat: HB file is not in RUA format", 0);
    }

    Nrhs = 0;
    rhsdata = NULL;
    if (Rhscrd != 0)
    {
        fgets(line, 82, infile);
        fgets(line, 82, infile);
        sscanf(line, "%3c %i %i", Type, &Nrhs, &dummy);
        rhsdata = new double[3*nrow];
    }
    fclose(infile);

    a     = new double[nnz];
    ja    = new int[nnz];
    ia    = new int[ncol+1];
    rhs   = new double[nrow];
    guess = new double[nrow];
    exact = new double[nrow];

    int flen = strlen(filename);
    F77NAME(readhb)(filename, &flen, &nrow, &nnz, &Nrhs, 
        a, ja, ia, rhs, guess, exact, rhsdata);

    //
    // process right-hand sides
    //
    int i;
    if (Nrhs != 0)
    {
	// cerr << "rhs provided" << endl;
	for (i=0; i<nrow; i++)
	    rhs[i] = rhsdata[i];
    }
    else
    {
        if (rhs_resource == HBTMAT_RHS_ONES)
	    for (i=0; i<nrow; i++)
	        rhs[i] = 1.0;
	else
	    for (i=0; i<nrow; i++)
	        ;//rhs[i] = (double) (2*rand()) / (double) RAND_MAX - 1.0;
    }

    if (Type[1] == 'G' || Type[1] == 'g')
    {
	// cerr << "Guess provided" << endl;
	for (i=0; i<nrow; i++)
	    guess[i] = rhsdata[i+nrow];
    }
    else
    {
	for (i=0; i<nrow; i++)
	    guess[i] = 0.0;
    }

    if (Type[2] == 'X' || Type[2] == 'x')
    {
	// cerr << "Exact solution provided" << endl;
	for (i=0; i<nrow; i++)
	    exact[i] = rhsdata[i+2*nrow];
    }
    else
    {
	for (i=0; i<nrow; i++)
	    exact[i] = 0.0;
    }

    // convert to 0-based indexing
    for (i=0; i<=nrow; i++)
	ia[i]--;
    for (i=0; i<ia[nrow]; i++)
	ja[i]--;

    // scale columns and rows in CSC format
    if (scale_resource)
    {
        double *factors = new double[nrow];
        double norm;
	int j;

        // scale columns and scale back exact solution
        for (i=0; i<nrow; i++)
        {
            norm = 0.0;
            for (j=ia[i]; j<ia[i+1]; j++)
                norm += a[j]*a[j];

            norm = 1.0/sqrt(norm);
            for (j=ia[i]; j<ia[i+1]; j++)
               a[j] *= norm;

            exact[i] /= norm;
        }
      
        // scale rows and rhs
        for (i=0; i<nrow; i++)
            factors[i] = 0.0;
      
        // first pass:  compute 2-norm of each column
        for (i=0; i<ia[nrow]; i++)
            factors[ja[i]] += a[i]*a[i];

        for (i=0; i<nrow; i++)
	{
            factors[i] = 1.0 / sqrt(factors[i]);
            rhs[i] *= factors[i];
	}
      
        // second pass: scale
        for (i=0; i<ia[nrow]; i++)
	       a[i] *= factors[ja[i]];
	
        delete [] factors;
    }
    
    // transpose
    if (transpose_resource)
    {
	assert (nrow == ncol);
        csrtrans(nrow, a, ja, ia);
    }

    delete [] rhsdata;
}

// This creates an object with data provided by the user.
// Assumes 0-based indexing.
// Will transpose, if necessary.
// Do not use FreeMat when this constructor is used.

HBTMat::HBTMat(int n, double *a_, int *ja_, int *ia_, 
  double *rhs_, double *guess_, double *exact_)
{
    owndata = FALSE;
    strcpy(title, "None");
    strcpy(key, "None");
    nrow = n;
    ncol = n;
    nnz = ia_[nrow];
    a = a_;
    ja = ja_;
    ia = ia_;
    rhs = rhs_;
    guess = guess_;
    exact = exact_;

    // Commented by K.Tsang on April 2, 2003.
    // Reason: Transpose is not required.
//      if (transpose_resource)
//          csrtrans(n, a, ja, ia);
}
