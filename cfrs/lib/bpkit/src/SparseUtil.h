//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _SPARSEUTIL_H_
#define _SPARSEUTIL_H_

void csrtrans(int n, double *a, int *ja, int *ia);

void csrtrans(int n, double *a, int *ja, int *ia,
    double *ao, int *jao, int *iao);

void shell_sort(
  const int n,
  int x[]);

void allocate_ilu(
  const int levfill,
  const int n,
  int *nzl, int *nzu,
  const int ia[], const int ja[],
  int *ial[], int *jal[],
  int *iau[], int *jau[],
  int growth);

void symbolic_ilu(
  const int levinc,
  const int n,
  int *nzl,
  int *nzu,
  const int ia[], const int ja[],
  int ial[], int jal[],
  int iau[], int jau[]);

void numeric_ilu(
  const double omega,
  const int n,
  const int ia[], const int ja[], const double a[],
  int ial[], int jal[], double al[],
  int iau[], int jau[], double au[],
  double thresh_resource, double rel_resource, double abs_resource);

void ilu_stat(
  const int n,
  const int ia[], const int ja[], const double a[],
  int ial[], int jal[], double al[],
  int iau[], int jau[], double au[]);

#endif // _SPARSEUTIL_H_
