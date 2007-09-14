//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include "BPKIT.h"
#include "SparseUtil.h"
#include "stdio.h" // kludge

// Transpose square matrix in place.
void csrtrans(int n, double *a, int *ja, int *ia)
{
    int i, j, k, next;

    double *ao  = new double[ia[n]];
    int    *jao = new int[ia[n]];
    int    *iao = new int[n+1];

    // initialize
    for (i=0; i<=n; i++)
        iao[i] = 0;

    // compute lengths of rows of transpose of A
    for (i=0; i<n; i++)
        for (k=ia[i]; k<ia[i+1]; k++)
            iao[ja[k]+1]++;

    // compute pointers from lengths
    iao[0] = 0;
    for (i=0; i<n; i++)
        iao[i+1] = iao[i] + iao[i+1];

    // now do the actual copying
    for (i=0; i<n; i++)
    {
        for (k=ia[i]; k<ia[i+1]; k++)
	{
            j = ja[k];
            next = iao[j];
            ao[next] = a[k];
            jao[next] = i;
            iao[j] = next+1;
        }
    }

    // reshift iao into ia
    for (i=n-1; i>=0; i--)
        ia[i+1] = iao[i];
    ia[0] = 0;

    // copy values and indices
    for (i=0; i<ia[n]; i++)
    {
	a[i] = ao[i];
	ja[i] = jao[i];
    }

    delete [] ao;
    delete [] jao;
    delete [] iao;
}

// Transpose square matrix, not in place.
void csrtrans(int n, double *a, int *ja, int *ia,
                     double *ao, int *jao, int *iao)
{
    int i, j, k, next;

    // initialize
    for (i=0; i<=n; i++)
        iao[i] = 0;

    // compute lengths of rows of transpose of A
    for (i=0; i<n; i++)
        for (k=ia[i]; k<ia[i+1]; k++)
            iao[ja[k]+1]++;

    // compute pointers from lengths
    iao[0] = 0;
    for (i=0; i<n; i++)
        iao[i+1] = iao[i] + iao[i+1];

    // now do the actual copying
    for (i=0; i<n; i++)
    {
        for (k=ia[i]; k<ia[i+1]; k++)
	{
            j = ja[k];
            next = iao[j];
            ao[next] = a[k];
            jao[next] = i;
            iao[j] = next+1;
        }
    }

    // reshift iao into ia
    for (i=n-1; i>=0; i--)
        iao[i+1] = iao[i];
    iao[0] = 0;
}

// shell sort
// stable, so it is fast if already sorted
// sorts x[0:n-1] in place, ascending order.

void shell_sort(
  const int n,
  int x[])
{
    int m, max, j, k, itemp;
    
    m = n/2;

    while (m > 0) {
        max = n - m;
        for (j=0; j<max; j++)
        {
            for (k=j; k>=0; k-=m)
            {
                if (x[k+m] >= x[k])
                    break;
                itemp = x[k+m];
                x[k+m] = x[k];
                x[k] = itemp;
            }
        }   
        m = m/2;
    }
}

// allocate space for integer data for level ILU

void allocate_ilu(
  const int levfill,                 // level of fill
  const int n,                       // order of matrix
  int *nzl, int *nzu,                // space allocated
  const int ia[], const int ja[],    // input
  int *ial[], int *jal[],            // output lower factor structure
  int *iau[], int *jau[],            // output upper factor structure
  int growth)                        // storage parameter
{
    int i, j, nzla, nzua;
    int maxl, maxu;

    maxu = n*(n+1)/2;
    maxl = n*n - maxu;

    // count number of entries in lower and upper triangular parts

    nzla = 0;
    nzua = 0;
    for (i=0; i<n; i++)
    {
        for (j=ia[i]; j<ia[i+1]; j++)
        {
            if (ja[j] < i)
               nzla++;
            else
               nzua++;
        }
    }

    *ial = new int[n+1];
    *iau = new int[n+1];

    if (levfill == 0) // ILU(0)
    {
        *nzl = nzla;
        *nzu = nzua + n;
    }
    else if (levfill == n) // full factorization
    {
        *nzl = maxl;
        *nzu = maxu;
    }
    else
    {
        *nzl = MIN((levfill+growth)*nzla, maxl);
        *nzu = MIN((levfill+growth)*nzua + n, maxu);
    }

    *jal = new int[*nzl];
    *jau = new int[*nzu];

    if (levfill != 0)
    {
        // cerr << "nnz in unfactored mat: " << nzla + nzua << endl;
        // cerr << "nnz allocated for ILU: " << *nzl + *nzu << endl;
    }
}

// symbolic level ILU
// fortran-style data structures and algorithms
// factors into separate upper and lower parts
// assumes rows are sorted  *** symbolic ILU
// assumes no zero rows

void symbolic_ilu(
  const int levfill,                 // level of fill
  const int n,                       // order of matrix
  int *nzl,                          // input-output
  int *nzu,                          // input-output
  const int ia[], const int ja[],    // input
  int ial[], int jal[],              // output lower factor structure
  int iau[], int jau[])              // output upper factor structure
{
#if 0
    // copy if levfill is 0
    if (levfill == 0)
    {
        int kl = 0;
        int ku = 0;
        ial[0] = 0;
        iau[0] = 0;
        for (int i=0; i<n; i++)
        {
            for (int j=ia[i]; j<ia[i+1]; j++)
            {
                if (ja[j] < i)
                   jal[kl++] = ja[j];
                else
                   jau[ku++] = ja[j];
            }
            ial[i+1] = kl;
            iau[i+1] = ku;
            shell_sort(ial[i+1]-ial[i], &jal[ial[i]]);
            shell_sort(iau[i+1]-iau[i], &jau[iau[i]]);
        }
        assert (kl <= *nzl); 
        assert (ku <= *nzu);
        *nzl = kl;
        *nzu = ku;
    }
    else
#endif
    {

    int *lnklst = new int[n];
    int *curlev = new int[n];
    int *levels = new int[*nzu];
    int *iwork = new int[n];

    int knzl = 0;
    int knzu = 0;

    ial[0] = 0;
    iau[0] = 0;

    for (int i=0; i<n; i++)
    {
        int first, next, j;

        // copy column indices of row into workspace and sort them

        int len = ia[i+1] - ia[i];
        next = 0;
        for (j=ia[i]; j<ia[i+1]; j++)
            iwork[next++] = ja[j];
        shell_sort(len, iwork);

        // construct implied linked list for row

        first = iwork[0];
        curlev[first] = 0;

        for (j=0; j<=len-2; j++)
        {
            lnklst[iwork[j]] = iwork[j+1];
            curlev[iwork[j]] = 0;
        }

        lnklst[iwork[len-1]] = n;
        curlev[iwork[len-1]] = 0;

        // merge with rows in U

        next = first;
        while (next < i)
        {
            int oldlst = next;
            int nxtlst = lnklst[next];
            int row = next;
            int ii;

            // scan row

            for (ii=iau[row]+1; ii<iau[row+1]; /*nop*/)
            {
                if (jau[ii] < nxtlst)
                {
                    // new fill-in
                    int newlev = curlev[row] + levels[ii] + 1;
                    if (newlev <= levfill)
                    {
                        lnklst[oldlst]  = jau[ii];
                        lnklst[jau[ii]] = nxtlst;
                        oldlst = jau[ii];
                        curlev[jau[ii]] = newlev;
                    }
                    ii++;
                }
                else if (jau[ii] == nxtlst)
                {
                    oldlst = nxtlst;
                    nxtlst = lnklst[oldlst];
                    int newlev = curlev[row] + levels[ii] + 1;
                    curlev[jau[ii]] = MIN(curlev[jau[ii]], newlev);
                    ii++;
                }
                else // (jau[ii] > nxtlst)
                {
                    oldlst = nxtlst;
                    nxtlst = lnklst[oldlst];
                }
            }
            next = lnklst[next];
        }
        
        // gather the pattern into L and U

        next = first;
        while (next < i)
        {
           if (knzl >= *nzl)
	     bperror("Not enough space allocated for symbolic factor", 0);

            jal[knzl++] = next;
            next = lnklst[next];
        }
        ial[i+1] = knzl;

        if (next != i)
        {
            // cerr << i << "  U has zero on diag, forcing nonzero" << endl;
            if (knzu >= *nzu)
	       bperror("Not enough space allocated for symbolic factor", 0);

            levels[knzu] = 2*n; // infinity
            jau[knzu++] = i;
        }

        while (next < n)
        {
            if (knzu >= *nzu)
	       bperror("Not enough space allocated for symbolic factor", 0);

            levels[knzu] = curlev[next];
            jau[knzu++] = next;
            next = lnklst[next];
        }
        iau[i+1] = knzu;
    }

    delete [] lnklst;
    delete [] curlev;
    delete [] levels;
    delete [] iwork;

    *nzl = knzl;
    *nzu = knzu;

    // cerr << "Actual nnz for ILU: " << *nzl + *nzu << endl;
    }
}

// assumes diag is first element in U

void numeric_ilu(
  const double omega,
  const int n, 
  const int ia[], const int ja[], const double a[],
  int ial[], int jal[], double al[],
  int iau[], int jau[], double au[],
  double thresh_resource, double rel_resource, double abs_resource)
{
    int i, k, kk, id, idd;
    double mult, modif;

    double *row = new double[n];
    int *marker = new int[n];

    for (i=0; i<n; i++)
    {
        row[i] = 0.0;
        marker[i] = 0;
    }

    for (i=0; i<n; i++)
    {
        // scatter row of A
        for (k=ia[i]; k<ia[i+1]; k++)
            row[ja[k]] = a[k];

        // scatter data structure of L and U
        for (k=ial[i]; k<ial[i+1]; k++)
            marker[jal[k]] = 1;
        for (k=iau[i]; k<iau[i+1]; k++)
            marker[jau[k]] = 1;
        
        modif = 0.0;

        // eliminate the elements in L in order
        for (k=ial[i]; k<ial[i+1]; k++)
        {
            id = jal[k];
            mult = row[id] / au[iau[id]];
            row[id] = mult;

            for (kk=iau[id]+1; kk<iau[id+1]; kk++)
            {
                idd = jau[kk];
                if (marker[idd]) 
                    row[idd] -= mult*au[kk];
		else if (omega != 0.0)
		    modif -= mult*au[kk];
            }
        }

        // gather resulting rows in L and U
        for (k=ial[i]; k<ial[i+1]; k++)
        {
            al[k] = row[jal[k]];
            row[jal[k]] = 0.0;
            marker[jal[k]] = 0;
        }
        for (k=iau[i]; k<iau[i+1]; k++)
        {
            au[k] = row[jau[k]];
            row[jau[k]] = 0.0;
            marker[jau[k]] = 0;
        }

	// regular modification
        if (omega != 0.0)
            au[iau[i]] += omega*modif;

        // stabilization
        if (ABS(au[iau[i]]) < thresh_resource)
	{
	    double sign = SGN(au[iau[i]]);

            au[iau[i]] = rel_resource * au[iau[i]];
            au[iau[i]] = au[iau[i]] + sign*abs_resource;
	}

        // check if still zero
        if (au[iau[i]] == 0.0)
	    bperror("Zero pivot encountered in numeric ILU factorization", i);
    }

    delete [] marker;
    delete [] row;
}

// print some statistics about the factors
// pivot elements are ordered first in U
// diagonal element of L is not stored

void ilu_stat(
  const int n, 
  const int ia[], const int ja[], const double a[],
  int ial[], int jal[], double al[],
  int iau[], int jau[], double au[])
{
    int i, k, countdd = 0;
    double maxa = 0.0, maxlu = 0.0, minpiv = 1.e50;
    double sum;

    for (i=0; i<ia[n]; i++)
        maxa = MAX(maxa, ABS(a[i]));

    for (i=0; i<n; i++)
        minpiv = MIN(minpiv, ABS(au[iau[i]]));

    //for (i=0; i<n; i++)
    //   printf("pivot:  %.2e\n", ABS(au[iau[i]]));

    for (i=0; i<ial[n]; i++)
        maxlu = MAX(maxlu, ABS(al[i]));

    for (i=0; i<n; i++)
        for (k=iau[i]+1; k<iau[i+1]; k++)
            maxlu = MAX(maxlu, ABS(au[i]));

    // fraction diagonally dominant
    for (i=0; i<n; i++)
    {
        sum = 0.0;
        for (k=ial[i]; k<ial[i+1]; k++)
	    sum += ABS(al[k]);
	if (sum <= 1.0) countdd++;

        sum = 0.0;
        for (k=iau[i]+1; k<iau[i+1]; k++)
            sum += ABS(au[k]);
        if (sum/au[iau[i]] <= 1.0) countdd++;
    }

/*
    cout << "max abs A = " << maxa << endl;
    cout << "max abs L+U = " << maxlu << endl;
    cout << "min pivot = " << minpiv << endl;
    cout << "fraction diag dom = " << countdd/(2.0*n) << endl;
*/
    int *levnum = new int[n];
    for (i=0; i<n; i++)
        levnum[i] = 0;

    int levi, nlev = 0;
    for (i=0; i<n; i++)
    {
         levi = 0;
	 for (k=ial[i]; k<ial[i+1]; k++)
            levi = MAX(levi, levnum[jal[k]]);
         levi++;
         levnum[i] = levi;
         nlev = MAX(nlev,levi);
    }
    delete [] levnum;

    printf("storage for preconditioner: %6d\n", ial[n]+iau[n]-2+n);

    printf("A LU PIV     %.2e   %.2e   %.2e\n", maxa, maxlu, 1.0/minpiv);
}
