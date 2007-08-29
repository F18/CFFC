
// miscellaneous functions

#ifndef _POLYFIT_H
#include "Polyfit.h"
#endif 


typedef std::vector<double>       vector_fp;


double polyfit(int n, double* x, double* y, double* w, 
	       int maxdeg, int& ndeg, double eps, double* r) {

  int nn = n;
  int mdeg = maxdeg;
  int ndg = ndeg;
  double epss = eps;
  int ierr;
  int worksize = 3*n + 3*maxdeg + 3;
  vector_fp awork(worksize,0.0);
  vector_fp coeffs(n+1, 0.0);
  double zer = 0.0;

  dpolft_(&nn, x, y, w, &mdeg, &ndg, &epss, &coeffs[0],
	  &ierr, &awork[0]);

  if (ierr != 1) std::cout << "\nDPOLFT returned error code IERR = " << ierr 
			   << " while attempting to fit " << n 
			   << " data points to a polynomial of degree " << maxdeg;

  ndeg = ndg;
  dpcoef_(&ndg, &zer, r, &awork[0]);

  return epss;
}   



