//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _CSRMAT_H_
#define _CSRMAT_H_

#include <cmath>
#include <cassert>

using namespace std;

#include "LocalMat.h"
#include "BlockVec.h"
#include "BpResource.h"
#include "SparseUtil.h"
#include "BPKIT.h"

class CSRMat : public LocalMat
{
private:

public:
    double *a;
    int    *ja;
    int    *ia;

    int nrow;
    int ncol;
    int nnz;

    CSRMat() {a = NULL; ja = ia = NULL; nrow = ncol = nnz = 0;}
   ~CSRMat() {delete [] a; delete [] ja; delete [] ia;}

    CSRMat(const int r, const int c);

    const double& val(unsigned int i) const {return a[i];}
    const int& row_ptr(unsigned int i) const {return ia[i];}
    const int& col_ind(unsigned int i) const {return ja[i];}

    int numrow() const {return nrow;}
    int numcol() const {return ncol;}
    int numnz()  const {return nnz;}

    void set(double *a_, int *ja_, int *ia_, int nrow_, int ncol_, int nnz_)
        {a=a_; ja=ja_; ia=ia_; nrow=nrow_; ncol=ncol_; nnz=nnz_;}

    // virtual functions

    double *& Data() {return a;}
    const double *Data() const {return a;}
    LocalMat *CreateEmpty() const {return new CSRMat();}
    LocalMat *CreateInv(LocalPrecon&) const;
    void SetToZero(int, int);
    void MatCopy(const LocalMat&);
    void Print(ostream&) const;

    // virtual mathematical functions

    void Mat_Trans(LocalMat *B) const;
    void Mat_Mat_Add(const LocalMat *B, LocalMat *C,
        double alpha) const;
    void Mat_Mat_Mult(const LocalMat *B, LocalMat *C, 
        double alpha, double beta) const;
    void Mat_Vec_Mult(const BlockVec& B, BlockVec& C,
        double alpha, double beta) const;
    void Mat_Trans_Vec_Mult(const BlockVec& B, BlockVec& C,
        double alpha, double beta) const;
    void Mat_Vec_Solve(const BlockVec& B, BlockVec& C) const;
    void Mat_Trans_Vec_Solve(const BlockVec& B, BlockVec& C)const;
};

/*
 *  Non-member function declarations
 */

ostream& operator << (ostream& os, const CSRMat& mat);

// resource names

#define APINV0_TRANSPOSE "APINV0.right" // 0=no, 1=yes
#define APINV_BANDED_TRANSPOSE "APINV_BANDED.right" // 0=no, 1=yes

#define APINVS_METHOD    "APINVS.method"    // 0=orig, 1=new
#define APINVS_GUESS     "APINVS.guess"     // 0=ident, 1=transp
#define APINVS_NOUTER    "APINVS.nouter"    //
#define APINVS_NINNER    "APINVS.ninner"    //
#define APINVS_TRANSPOSE "APINVS.right"     // 0=no, 1=yes
#define APINVS_EPSILON   "APINVS.epsilon"   //
#define APINVS_PRINT     "APINVS.print"     //

#define RILUK_GROWTH     "RILUK.growth"     //
#define RILUK_THRESH     "RILUK.thresh"     //
#define RILUK_REL        "RILUK.rel"        //
#define RILUK_ABS        "RILUK.abs"        //
#define RILUK_STAT       "RILUK.stat"       //

#define ILUT_PERMTOL     "ILUT.permtol"     //
#define ILUT_PBLOCK      "ILUT.pblock"      //
#define ILUT_STAT        "ILUT.stat"        //

/*
 *  Derivations
 */

class CSRMat_INVERSE : public CSRMat
{
private:
public:
    CSRMat_INVERSE(const CSRMat&);
    void Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
        {solve_is_mult(B, X);}
};
    
class CSRMat_LU : public CSRMat
{
private:
    double *al;    // lower triangular factor
    int    *jal;
    int    *ial;
    double *au;    // upper triangular factor
    int    *jau;
    int    *iau;

public:
    CSRMat_LU(const CSRMat&);
   ~CSRMat_LU();
    void Mat_Vec_Solve(const BlockVec&, BlockVec&) const;
};
    
class CSRMat_RILUK : public CSRMat
{
private:
    double *al;    // lower triangular factor
    int    *jal;
    int    *ial;
    double *au;    // upper triangular factor
    int    *jau;
    int    *iau;

public:
    CSRMat_RILUK(const CSRMat&, int, double);
   ~CSRMat_RILUK();
    void Mat_Vec_Solve(const BlockVec&, BlockVec&) const;

    static int growth;
    static double stab_thresh;
    static double stab_rel;
    static double stab_abs;
    static int stat;
};
    
class CSRMat_ILUT : public CSRMat
{
private:
    int    numrow;
    double *alu;
    int    *jlu;
    int    *ju;

public:
    CSRMat_ILUT(const CSRMat&, int, double);
   ~CSRMat_ILUT();
    void Mat_Vec_Solve(const BlockVec&, BlockVec&) const;

    static double permtol;
    static int    pblock;
    static int    stat;
};

class CSRMat_APINV0 : public CSRMat
{
public:
    CSRMat_APINV0(const CSRMat&);
    void Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
        {solve_is_mult(B, X);}

    static int transpose;
};

class CSRMat_APINVS : public CSRMat
{
public:
    CSRMat_APINVS(const CSRMat&, int, int);
    void Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
        {solve_is_mult(B, X);}

    // resources
    static int method;
    static int guess;
    static int nouter;
    static int ninner;
    static int transpose;
    static double epsilon;
    static int print;
};

class CSRMat_DIAG : public CSRMat
{
private:
public:
    CSRMat_DIAG(const CSRMat&);
    void Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
        {solve_is_mult(B, X);}
};

class CSRMat_TRIDIAG : public CSRMat
{
private:
    double *val;  // ncol*3 matrix in sparse diagonal format
public:
    CSRMat_TRIDIAG(const CSRMat&);
   ~CSRMat_TRIDIAG() {delete [] val;}
    void Mat_Vec_Solve(const BlockVec&, BlockVec&) const;
};

class CSRMat_SOR : public CSRMat
{
private:
    const CSRMat *Ap;
    int iterations_;
    double omega_;
    int *idiag;      // pointer to diagonal element
public:
    CSRMat_SOR(const CSRMat&, const double, const int);
   ~CSRMat_SOR() {delete [] idiag;}
    void Mat_Vec_Solve(const BlockVec&, BlockVec&) const;
};

class CSRMat_SSOR : public CSRMat
{
private:
    const CSRMat *Ap;
    int iterations_;
    double omega_;
    int *idiag;      // pointer to diagonal element
public:
    CSRMat_SSOR(const CSRMat&, const double, const int);
   ~CSRMat_SSOR() {delete [] idiag;}
    void Mat_Vec_Solve(const BlockVec&, BlockVec&) const;
};

class CSRMat_GMRES : public CSRMat
{
private:
    const CSRMat *Ap;
    int iterations_;
    double threshold_;
public:
    CSRMat_GMRES(const CSRMat&, const int, const double);
    void Mat_Vec_Solve(const BlockVec&, BlockVec&) const;
};

class CSRMat_APINV_BANDED : public CSRMat
{
public:
    CSRMat_APINV_BANDED(const CSRMat&, int);
    void Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
        {solve_is_mult(B, X);}

    static int transpose;
};

class CSRMat_APINV_TRUNC : public CSRMat
{
private:
public:
    CSRMat_APINV_TRUNC(const CSRMat&, int);
    void Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
        {solve_is_mult(B, X);}
};

#endif // _CSRMAT_H_
