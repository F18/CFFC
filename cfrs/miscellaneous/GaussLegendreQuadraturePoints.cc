#include <iostream>
#include <iomanip.h>

void gauleg(double x[], double w[], int n);

int main(){

  int n;
  double *x;
  double *w;
  cout << "This subroutine returns the abscissas and the weights" << endl
       << "of the Gauss-Legendre n-point quadrature formula." << endl
       << "There are both the negative and the positive values!" << endl;
  cout << "Introduce the number of points --> n= ";

  if (cin >> n) {

    x = new double [n+1];
    w = new double [n+1];
  
    gauleg(x,w,n);
  
    cout.precision(15);
    cout.setf(ios::left);
    for (int i=1; i<=n; i++){
      cout << "abscissa= " << setw(18) <<  x[i] << "\t"
	   << "weight= "<< setw(18) << w[i] << endl;
    }
  
    delete []x; x=NULL;
    delete []w; w=NULL; 
  }
  else
    cout << "Please introduce a number next time. Bye!" << endl;  
  return 0;
  
}
