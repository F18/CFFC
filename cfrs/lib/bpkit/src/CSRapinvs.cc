//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include <cmath>

using namespace std;

#include "BPKIT.h"
#include "CSRMat.h"

// prototypes for this file

static void apinvb0(int n, int mcol, 
  const int *ia, const int *ja, const double *a, 
  int *im, int *jm, double *m, 
  const int *ib, const int *jb, const double *b,
  int lfil, int self_precon, int nouter, int ninner);

static void apinvb(int n, int mcol, 
  const int *ia, const int *ja, const double *a, 
  int *im, int *jm, double *m, 
  const int *ib, const int *jb, const double *b,
  int lfil, int self_precon, int nouter, double epsilon);

static void matscale3(int n, int *ia, int *ja, double *a, double& alpha);

static void matscale3t(int n, int *ia, int *ja, double *a, double& alpha);

static void init_trans(int n, double *ao, int *jao, int *iao,
  double *a, int *ja, int *ia, int lfil, double alpha);

static void spminres(const int nrhs, const int *irhs, const double *rhs,
  int& ns, int *is, double *s, const double *a, const int *ja, const int *ia,
  const double *m, const int *jm, const int *im, int lfil, int *iwk,
  int self_precon, double *r, int *ir, double *t, int *it, double *q, int *iq);

static double spminres2(const int nrhs, const int *irhs, const double *rhs,
  int& ns, int *is, double *s, const double *a, const int *ja, const int *ia,
  const double *m, const int *jm, const int *im, int lfil, int *iwk,
  int self_precon, double *r, int *ir, double *t, int *it, double *q, int *iq);

static void spamuxc(const int *ia, const int *ja, const double *a,
  const int nv, const int *iv, const double *v,
  int& nb, int *ib, double *b, int *iwk);

static void spamuxc2(int lfil, 
  const int *ia, const int *ja, const double *a,
  const int nv, const int *iv, const double *v,
  int& nb, int *ib, double *b, int *iwk);

static void spdot(const int na, const int *ia, const double *a,
  const int nb, const int *ib, const double *b, double &prod, int *iwk);

static void spaxpy(int& na, int *ia, double *a,
  const int nb, const int *ib, const double *b, double w, int *iwk);

static void dropt(double *t, int *it, int& nt, 
  const int *ix, const int nx, int lfil, int *iwk);

static void drop(double *t, int *it, int& nt, int lfil);

// resources

int CSRMat_APINVS::method     = BpGetInt(APINVS_METHOD,     1);
int CSRMat_APINVS::guess      = BpGetInt(APINVS_GUESS,      0);
int CSRMat_APINVS::nouter     = BpGetInt(APINVS_NOUTER,     1);
int CSRMat_APINVS::ninner     = BpGetInt(APINVS_NINNER,     1);
int CSRMat_APINVS::transpose  = BpGetInt(APINVS_TRANSPOSE,  1);
double CSRMat_APINVS::epsilon = BpGetDouble(APINVS_EPSILON, 0.0);
int CSRMat_APINVS::print      = BpGetInt(APINVS_PRINT, 0);

CSRMat_APINVS::CSRMat_APINVS(const CSRMat& A, int lfil, int self_precon)
{
    int i;
    double alpha;

    int n = nrow = ncol = A.nrow;
    nnz = n * lfil;

    // pointers to transpose
    double *ao;
    int *jao;
    int *iao;

    // transpose matrix if want right apinv, or transpose init guess
    if (transpose || guess)
    {
        ao = new double[A.nnz];
        jao = new int[A.nnz];
        iao = new int[n+1];
	csrtrans(n, A.a, A.ja, A.ia, ao, jao, iao);
    }

    // allocate space for approximate inverse
    a = new double[nnz];
    ja = new int[nnz];
    ia = new int[n+1];

    if (!guess)
    {
	// scaled identity initial guess
        matscale3(n, A.ia, A.ja, A.a, alpha);
        for (i=0; i<n; i++)
        {
             a[i*lfil] = alpha;
            ja[i*lfil] = i;
            ia[i] = 1;
        }
    }
    else
    {
	// scaled transpose initial guess
        matscale3t(n, A.ia, A.ja, A.a, alpha);

	init_trans(n, ao, jao, iao, a, ja, ia, lfil, alpha);
    }

    // form matrix of right-hand sides
    double *b = new double[n];
    int *jb = new int[n];
    int *ib = new int[n+1];

    for (i=0; i<n; i++)
    {
        b[i] = 1.0;
        jb[i] = i;
        ib[i] = i;
    }
    ib[n] = n;

    if (method == 0)
    {
	if (transpose)
            apinvb0(n, n, iao, jao, ao, ia, ja, a, ib, jb, b,
                lfil, self_precon, nouter, ninner);
	else
            apinvb0(n, n, A.ia, A.ja, A.a, ia, ja, a, ib, jb, b,
                lfil, self_precon, nouter, ninner);
    }
    else
    {
        // call with nouter = lfil
	if (transpose)
            apinvb(n, n, iao, jao, ao, ia, ja, a, ib, jb, b,
                lfil, self_precon, lfil, epsilon);
	else
            apinvb(n, n, A.ia, A.ja, A.a, ia, ja, a, ib, jb, b,
                lfil, self_precon, lfil, epsilon);
    }

    delete [] b;
    delete [] jb;
    delete [] ib;

    if (transpose || guess)
    {
        delete [] ao;
        delete [] jao;
        delete [] iao;
    }

    // convert to regular indexing, squeeze out the zeros

    int j, k = 0;
    for (i=0; i<n; i++)
    {
	for (j=i*lfil; j<i*lfil+ia[i]; j++)
	{
	    a[k] =  a[j];
	   ja[k++] = ja[j];
	}
    }
    for (i=n-1; i>=0; i--)
        ia[i+1] = ia[i];
    ia[0] = 0;
    for (i=0; i<n; i++)
	ia[i+1] = ia[i] + ia[i+1];

    // transpose the solution in place
    if (transpose)
        csrtrans(n, a, ja, ia);

    if (print)
       cout << *this << endl;
}

// Return alpha so that rho(I - alpha A) is minimized.
// The scale factor is tr(A)/tr(A'A) = tr(A)/(norm(A,fro)**2)

static void matscale3(
  int n, 
  int *ia, 
  int *ja, 
  double *a, 
  double& alpha)
{
      int i, j;
      double tr, fnorm2;

      // first compute trace
      tr = 0.0;
      for (i=0; i<n; i++)
	 for (j=ia[i]; j<ia[i+1]; j++)
	    if (ja[j] == i)
	       tr += a[j];

      // Compute Frobenius norm squared.
      fnorm2 = 0.0;
      for (i=0; i<ia[n]; i++)
	 fnorm2 += a[i]*a[i];

      alpha = tr/fnorm2;
}

// Return alpha so that rho(I - alpha C) is minimized.
// The scale factor is tr(C)/tr(C'C) = tr(C)/(norm(C,fro)**2)
// where C = A A'

static void matscale3t(
  int n, 
  int *ia, 
  int *ja, 
  double *a, 
  double& alpha)
{
      int i, j;
      double tr, fnorm2, prod;

      int *iwk = new int[n];
      for (i=0; i<n; i++)
          iwk[i] = 0;

      // Multiply A*A', producing rows one at a time.

      tr = 0.0;
      fnorm2 = 0.0;
      for (i=0; i<n; i++)    // rows of A
      {
         for (j=0; j<n; j++) // cols of A'
	 {
            spdot(ia[i+1]-ia[i], &ja[ia[i]], &a[ia[i]], 
                ia[j+1]-ia[j], &ja[ia[j]], &a[ia[j]], prod, iwk);

            if (i == j) 
                tr += prod;
            fnorm2 += prod*prod;
         }
      }
      alpha = tr/fnorm2;
      delete [] iwk;
}

// the transpose is available in ao, jao, iao

static void init_trans(int n, double *ao, int *jao, int *iao, 
  double *a, int *ja, int *ia, int lfil, double alpha)
{
    int i, j, ip;
    double *ax = new double[n];
    int *jax = new int[n];

    for (i=0; i<n; i++)
    {
        ia[i] = iao[i+1] - iao[i];

        if (ia[i] > lfil)
        {
            ip = 0;
            for (j=iao[i]; j<iao[i+1]; j++, ip++)
            {
                 ax[ip] = alpha * ao[j];
                jax[ip] = jao[j];
            }
            drop(ax, jax, ia[i], lfil);
            ip = i*lfil;
            for (j=0; j<ia[i]; j++, ip++)
            {
                 a[ip] =  ax[j];
                ja[ip] = jax[j];
            }
        }
        else
        {
            int ip = i*lfil;
            for (j=iao[i]; j<iao[i+1]; j++, ip++)
            {
                 a[ip] = alpha * ao[j];
                ja[ip] = jao[j];
            }
        }
    }
    delete [] jax;
    delete [] ax;
}

// calls original version
static void apinvb0(
  int n, 
  int mcol, 
  const int *ia, 
  const int *ja, 
  const double *a, 
  int *im, 
  int *jm, 
  double *m, 
  const int *ib, 
  const int *jb, 
  const double *b,
  int lfil, 
  int self_precon,
  int nouter,
  int ninner)
{
    int iouter, iinner, i;

    // initialize subprogram workspace
    int *iwk = new int[n];
    for (i=0; i<n; i++)
        iwk[i] = -1;

    // allocate r, t, q
    double *r = new double[n];
    double *t = new double[n];
    double *q = new double[n];
    int *ir = new int[n];
    int *it = new int[n];
    int *iq = new int[n];

    // outer iterations
    for (iouter=0; iouter<nouter; iouter++)
    {
        // loop on the columns in rhs
        for (i=0; i<mcol; i++)
        {
	    // inner iterations
            for (iinner=0; iinner<ninner; iinner++)
            {
                spminres(ib[i+1]-ib[i], &jb[ib[i]], &b[ib[i]], 
                    im[i], &jm[i*lfil], &m[i*lfil],
                    a, ja, ia, m, jm, im, lfil, iwk, self_precon,
                    r, ir, t, it, q, iq);
            }
        }
    }

    delete [] r;
    delete [] t;
    delete [] q;
    delete [] ir;
    delete [] it;
    delete [] iq;
    delete [] iwk;
}

static void apinvb(
  int n, 
  int mcol, 
  const int *ia, 
  const int *ja, 
  const double *a, 
  int *im, 
  int *jm, 
  double *m, 
  const int *ib, 
  const int *jb, 
  const double *b,
  int lfil, 
  int self_precon,
  int nouter,
  double epsilon)
{
    int iouter, i;

    int *done = new int[n];
    for (i=0; i<n; i++)
        done[i] = FALSE;

    // initialize subprogram workspace
    int *iwk = new int[n];
    for (i=0; i<n; i++)
        iwk[i] = -1;

    // allocate r, t, q
    double *r = new double[n];
    double *t = new double[n];
    double *q = new double[n];
    int *ir = new int[n];
    int *it = new int[n];
    int *iq = new int[n];

    // outer iterations
    for (iouter=0; iouter<nouter; iouter++)
    {
        // loop on the columns in rhs
        for (i=0; i<mcol; i++)
        {
	    if (!done[i])
	    {
                double res;
                res = spminres2(ib[i+1]-ib[i], &jb[ib[i]], &b[ib[i]], 
                    im[i], &jm[i*lfil], &m[i*lfil],
                    a, ja, ia, m, jm, im, lfil, iwk, self_precon,
                    r, ir, t, it, q, iq);
	        if (res <= epsilon)
	   	    done[i] = TRUE;
            }
        }
    }

    delete [] r;
    delete [] t;
    delete [] q;
    delete [] ir;
    delete [] it;
    delete [] iq;
    delete [] iwk;
    delete [] done;
}

// original version
static void spminres(
  const int nrhs,
  const int *irhs,
  const double *rhs,
  int& ns,                  // in and out
  int *is,
  double *s,
  const double *a,
  const int *ja,
  const int *ia,
  const double *b,
  const int *jb,
  const int *ib,
  int lfil,
  int *iwk,
  int self_precon,
  double *r, int *ir,
  double *t, int *it,
  double *q, int *iq)
{
    int j;
    double dot1, dot2, wstep;

    int nr, nt, nq;

    // compute residual
    spamuxc(ia, ja, a, ns, is, s, nr, ir, r, iwk);
    spaxpy(nr, ir, r, nrhs, irhs, rhs, -1.0, iwk);
    for (j=0; j<nr; j++)
        r[j] = -r[j];

    // call self-preconditioner, if necessary
    if (self_precon)
    {
        spamuxc2(lfil, ib, jb, b, nr, ir, r, nt, it, t, iwk);
    }
    else
    {
        nt = nr;
        for (j=0; j<nr; j++)
        {
             t[j] =  r[j];
            it[j] = ir[j];
        }
    }

    // q = A*t
    spamuxc(ia, ja, a, nt, it, t, nq, iq, q, iwk);

    // dot1 = r dot q
    spdot(nr, ir, r, nq, iq, q, dot1, iwk);

    // dot2 = q dot q
    dot2 = 0.0;
    for (j=0; j<nq; j++)
        dot2 += q[j]*q[j];
    if (dot2 == 0.0)
        return;

    wstep = dot1/dot2;

    // first copy s into full length workspace, q
    nq = ns;
    for (j=0; j<ns; j++)
    {
         q[j] =  s[j];
        iq[j] = is[j];
    }

    // s = s + wstep*t
    spaxpy(nq, iq, q, nt, it, t, wstep, iwk);

    // drop
    drop(q, iq, nq, lfil);

    // copy from q back to s
    ns = nq;
    for (j=0; j<nq; j++)
    {
         s[j] =  q[j];
        is[j] = iq[j];
    }
}


// smart version
static double spminres2(
  const int nrhs,
  const int *irhs,
  const double *rhs,
  int& ns,                  // in and out
  int *is,
  double *s,
  const double *a,
  const int *ja,
  const int *ia,
  const double *b,
  const int *jb,
  const int *ib,
  int lfil,
  int *iwk,
  int self_precon,
  double *r, int *ir,
  double *t, int *it,
  double *q, int *iq)
{
    int j;
    double dot1, dot2, wstep;

    int nr, nt, nq;

    // compute residual
    spamuxc(ia, ja, a, ns, is, s, nr, ir, r, iwk);
    spaxpy(nr, ir, r, nrhs, irhs, rhs, -1.0, iwk);
    for (j=0; j<nr; j++)
        r[j] = -r[j];

    // call self-preconditioner, if necessary
    if (self_precon)
    {
        spamuxc2(lfil, ib, jb, b, nr, ir, r, nt, it, t, iwk);
    }
    else
    {
        nt = nr;
        for (j=0; j<nr; j++)
        {
             t[j] =  r[j];
            it[j] = ir[j];
        }
    }

    // drop in t
    dropt(t, it, nt, is, ns, lfil, iwk);

    // q = A*t
    spamuxc(ia, ja, a, nt, it, t, nq, iq, q, iwk);

    // dot1 = r dot q
    spdot(nr, ir, r, nq, iq, q, dot1, iwk);

    // dot2 = q dot q
    dot2 = 0.0;
    for (j=0; j<nq; j++)
        dot2 += q[j]*q[j];
    if (dot2 == 0.0)
        return 0.0;

    wstep = dot1/dot2;

    // this is the only place s is modified, in place
    // s = s + wstep*t
    spaxpy(ns, is, s, nt, it, t, wstep, iwk);

    // update the residual; compute and return its norm
    spaxpy(nr, ir, r, nq, iq, q, -wstep, iwk);

    dot2 = 0.0;
    for (j=0; j<nr; j++)
        dot2 += r[j]*r[j];

    return dot2;
}

static void spamuxc(
  const int *ia,
  const int *ja,
  const double *a,
  const int nv,
  const int *iv,
  const double *v,
  int& nb,
  int *ib,
  double *b,
  int *iwk)
{
    int i, j, k;
    double vmult;

    nb = 0;
    for (i=0; i<nv; i++)
    {
        vmult = v[i];
        for (j=ia[iv[i]]; j<ia[iv[i]+1]; j++)
        {
            k = iwk[ja[j]];
	    if (k != -1)
	    {
		b[k] = b[k] + a[j]*vmult;
	    }
	    else
	    {
		ib[nb] = ja[j];
		b[nb] = a[j]*vmult;
		iwk[ja[j]] = nb++;
	    }
	}
    }

    // reset full-length work vector back to -1
    for (i=0; i<nb; i++)
	iwk[ib[i]] = -1;
}

static void spamuxc2(
  int lfil,
  const int *ia,
  const int *ja,
  const double *a,
  const int nv,
  const int *iv,
  const double *v,
  int& nb,
  int *ib,
  double *b,
  int *iwk)
{
    int i, j, k;
    double vmult;

    nb = 0;
    for (i=0; i<nv; i++)
    {
	vmult = v[i];

	// next three lines differ from spamuxc only
	int jstrt = iv[i]*lfil;
	int jstop = jstrt + ia[iv[i]];
	for (j=jstrt; j<jstop; j++)
	{
	    k = iwk[ja[j]];
	    if (k != -1)
	    {
		b[k] = b[k] + a[j]*vmult;
	    }
	    else
	    {
		ib[nb] = ja[j];
		b[nb] = a[j]*vmult;
		iwk[ja[j]] = nb++;
	    }
	}
    }

    // reset full-length work vector back to -1
    for (i=0; i<nb; i++)
	iwk[ib[i]] = -1;
}

static void spdot(
  const int na,
  const int *ia,
  const double *a,
  const int nb,
  const int *ib,
  const double *b,
  double &prod,
  int *iwk)
{
    int i, k;

    // scatter index of b into full-length work vector
    for (i=0; i<nb; i++)
       iwk[ib[i]] = i;

    // accumulate dot product using sparsity of a
    prod = 0.0;
    for (i=0; i<na; i++)
    {
       k = iwk[ia[i]];
       if (k != -1) 
          prod += a[i]* b[k];
    }

    // reset full-length work vector back to -1
    for (i=0; i<nb; i++)
       iwk[ib[i]] = -1;
}

static void spaxpy(
  int& na,
  int *ia,
  double *a,
  const int nb,
  const int *ib,
  const double *b,
  double w,
  int *iwk)
{
    int i, j;

    // scatter index of b into full-length work vector
    for (i=0; i<nb; i++)
       iwk[ib[i]] = i;

    // do daxpy using sparsity of a
    for (i=0; i<na; i++)
    {
       j = ia[i];
       if (iwk[j] != -1)
       {
          a[i] += w*b[iwk[j]];
          iwk[j] = -1;
       }
    }

    // do daxpy using sparsity of b, not in a
    for (i=0; i<nb; i++)
    {
       j = ib[i];
       if (iwk[j] != -1)
       {
          a[na] = w*b[iwk[j]];
          ia[na++] = j;
          // reset element in full-length work vector back to -1
          iwk[j] = -1;
       }
    }
}

// n is number of elements
// ncut is number desired

static void qsplit(double *a, int *ind, int n, int ncut)
{
        double tmp, abskey;
        int itmp, first, last, mid;

	ncut--;
        first = 0;
        last = n-1;
        if (ncut < first || ncut > last) return;

//    outer loop -- while mid .ne. ncut do

      while (1)
      {
        mid = first;
        abskey = ABS(a[mid]);
        for (int j=first+1; j<=last; j++)
	{
           if (ABS(a[j]) > abskey)
           {
              mid = mid+1;
              // interchange
              tmp = a[mid];
              itmp = ind[mid];
              a[mid] = a[j];
              ind[mid] = ind[j];
              a[j]  = tmp;
              ind[j] = itmp;
	   }
	}

//      interchange

        tmp = a[mid];
        a[mid] = a[first];
        a[first]  = tmp;

        itmp = ind[mid];
        ind[mid] = ind[first];
        ind[first] = itmp;

//      test for while loop

        if (mid == ncut) return;
        if (mid > ncut)
           last = mid-1;
        else
           first = mid+1;

     } // endwhile
}

static void drop(double *t, int *it, int& nt, int lfil)
{
    if (nt > lfil)
    {
	qsplit(t, it, nt, lfil);
	nt = lfil;
    }
}

//    Drop in t, not in sparsity of x
//    plus max element if not full lfil
static void dropt(
  double *t,
  int *it,
  int& nt,   // in and out
  const int *ix,
  const int nx,
  int lfil,
  int *iwk)
{
      int i, k, imax; //, imax2;
      double  tmax; //, tmax2;

      // scatter indices of x
      for (i=0; i<nx; i++)
         iwk[ix[i]] = 1;

      // copy over original t
      k = 0;
      imax = 0;
      tmax = 0.0;
      // imax2 = 0;
      // tmax2 = 0.0;
      for (i=0; i<nt; i++)
      {
	 // in the pattern
         if (iwk[it[i]] != -1) 
	 {
             t[k] =  t[i];
            it[k] = it[i];
            k++;
	 }
         else
	 {
	    // keep track of two largest
            if (ABS(t[i]) > ABS(tmax)) 
	    {
#if 0
               tmax2 = tmax;
               imax2 = imax;
#endif
               tmax = t[i];
               imax = it[i];
	    }
#if 0
            else if (ABS(t[i]) > ABS(tmax2)) 
	    {
               tmax2 = t[i];
               imax2 = it[i];
            }
#endif
         }
      }
      nt = k;

      // assert(nt <= nx);
      // add extra element
      if (nx < lfil && tmax != 0.0)   // nt changed to nx on 09/23/95
      {
          t[nt] = tmax;
         it[nt] = imax;
         nt++;
      }

#if 0 // only needed for exchange stage
      // see above fix
      // add next extra element on end
      if (nt < lfil && tmax2 != 0.0) 
      {
          t[nt] = tmax2;
         it[nt] = imax2;
      }
#endif

      // reset workspace
      for (i=0; i<nx; i++)
         iwk[ix[i]] = -1;
}
