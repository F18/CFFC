
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


double polyfit_smoothing(int N, double *x, double *y, int order, double window, bool window_flag, double *y_smooth) {
    
    double k=window;
    int m = int(window/2.0)*2; // even integer number
    
    double *x_temp = new double [N];
    double *y_temp = new double [N];
    int *indexes = new int [N];
    double *weight = new double [N];
    double *coeff = new double [order+1];
    double eps;
    int index;
    int ndeg;
    double xmin, xmax;
    for (int i=0; i<N; i++) {
        
        if(window_flag) {    // use the elements in the physical space
            xmin = max(x[i]-k/2.0,x[0]);
            xmax = min(x[i]+k/2.0,x[N-1]);
            
            if(xmin==x[0])
                xmax = xmin+k;
            if(xmax==x[N-1])
                xmin = xmax-k;
            
            
            index = 0;
            for (int j=0; j<N; j++) {
                if (x[j]>=xmin && x[j]<=xmax){
                    x_temp[index] = x[j];
                    y_temp[index] = y[j];
                    weight[index] = -1;
                    indexes[index]=j;
                    index++;
                }
            }
            
            if (order >= index){
                int last_index = index-1;
                // add additional indexes
                if(xmin == x[0]) {
                    // add indexes only to right

                    for (int j=1; j<(order-last_index+1); j++){
                        x_temp[index] = x[indexes[last_index]+j];
                        y_temp[index] = y[indexes[last_index]+j];
                        weight[index] = -1;
                        indexes[index] = indexes[last_index]+j;
                        index++;
                    }
                } else if(xmax == x[N-1]) {
                    // add indexes only to left
                    for (int j=1; j<(order-last_index+1); j++){
                        x_temp[index] = x[indexes[0]-j];
                        y_temp[index] = y[indexes[0]-j];
                        weight[index] = -1;
                        indexes[index] = indexes[0]-j;
                        index++;
                    }
                } else {
                    // add indexes to left and right
                    int j = 1;
                    while (order >= index){

                        // add index to right
                        x_temp[index] = x[indexes[last_index]+j];
                        y_temp[index] = y[indexes[last_index]+j];
                        weight[index] = -1;
                        indexes[index] = indexes[last_index]+j;
                        index++;
                        
                        // add index to left
                        x_temp[index] = x[indexes[0]-j];
                        y_temp[index] = y[indexes[0]-j];
                        weight[index] = -1;
                        indexes[index] = indexes[0]-j;
                        index++;
                        
                        j++;
                    }
                }
                    
            }
            
        } else {            // use the elements in the index space
            index = 0;
            for (int j=max(i-m/2,0); j<=min(i+m/2,N-1); j++) {
                x_temp[index] = x[j];
                y_temp[index] = y[j];
                weight[index] = -1;
                index++;
            }
        }
        

        // double polyfit(int n, double* x, double* y, double* w, int maxdeg, int& ndeg, double eps, double* r);
        eps = max(eps, polyfit(index, x_temp, y_temp, weight, order, ndeg, 0, coeff) );
        
        y_smooth[i] = polyN(x[i],coeff,ndeg);
        
    }
    
    return eps;
    
}


double polyfit_smoothing_logscale(int N, double *x, double *y, int order, double window, bool window_flag, double *y_smooth) {
    
    double *logx = new double [N];
    double *logy = new double [N];
    double *logy_smooth = new double [N];
    double eps;

    if(x[0]==0){
        int Nlog = N-1;
        for (int i=1; i<N; i++) {
            logx[i-1] = log10(x[i]);
            logy[i-1] = log10(y[i]);
        }
        eps = polyfit_smoothing(Nlog,logx,logy,order,window,window_flag,logy_smooth);
        
        y_smooth[0]=y[0];
        for (int i=1; i<N; i++)
            y_smooth[i] = pow(10.0,logy_smooth[i-1]);
        
    } else {

        for (int i=0; i<N; i++) {
            logx[i] = log10(x[i]);
            logy[i] = log10(y[i]);
        }
        eps = polyfit_smoothing(N,logx,logy,order,window,window_flag,logy_smooth);
        
        for (int i=0; i<N; i++)
            y_smooth[i] = pow(10,logy_smooth[i]);
    }
    
    delete[] logx;
    delete[] logy;
    delete[] logy_smooth;
    
    return eps;
}


