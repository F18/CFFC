//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _HBTMAT_H_
#define _HBTMAT_H_

#include <cstring>
#include <iostream>

using namespace std;

// Harwell-Boeing Transposed format

class HBTMat
{
private:
    char title[73];
    char key[9];

    double *a;
    int *ja;
    int *ia;

    int nrow;
    int ncol;
    int nnz;

    double *rhs;
    double *guess;
    double *exact;

    int owndata;

public:
    /* Added HBTMat() by K. Tsang on March 27, 2003. */
    HBTMat() {//cout << " HBTMat constructor......" << endl;
              a = NULL; ja = NULL; ia = NULL; 
              nrow = ncol = nnz = 0; owndata = 0;
	      rhs = NULL; guess = NULL; exact = NULL;
	      strcpy(title, "None"); strcpy(key, "None");}

    HBTMat(const char *);
    HBTMat(int, double *, int *, int *, double *, double *, double *);
   ~HBTMat();
    void FreeMat();

    const double&  val(int i) const { return a[i]; }
    const int& row_ptr(int i) const { return ia[i]; }
    const int& col_ind(int i) const { return ja[i]; }

    int dimrow() const {return nrow;}
    int dimcol() const {return ncol;}
    int numnz() const {return nnz;}

    double *get_rhs() const {return rhs;}
    double *get_guess() const {return guess;}
    double *get_exact() const {return exact;}

    static int rhs_resource;
    static int transpose_resource;
    static int scale_resource;
};

#define HBTMAT_RHS               "HBTMat.rhs"
#define HBTMAT_RHS_RAND          0
#define HBTMAT_RHS_ONES          1

#define HBTMAT_TRANSPOSE         "HBTMat.transpose"

#define HBTMAT_SCALE             "HBTMat.scale"

#endif // _HBTMAT_H_
