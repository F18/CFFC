#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double linearregression(double *xx, double *yy, int ni, int ne);

int main (void) { 

  int num = 4;
  double N[num];
  double L1[num], L2[num], Max[num];

  // Number of nodes.
  N[0] = 1.0/sqrt(104.0);
  N[1] = 1.0/sqrt(410.0);
  N[2] = 1.0/sqrt(1641.0);
  N[3] = 1.0/sqrt(6556.0);
  //N[4] = 1.0/sqrt(23649.0);
 
  // L1 Norm.
  L1[0] = 0.0006089960;
  L1[1] = 0.0001359110;
  L1[2] = 0.0000353185;
  L1[3] = 0.0000178633;
  //L1[4] = 0;
  cout << "  L1-norm = " << linearregression(N,L1,0,num) << endl;
  cout << "  L1-norm (first 2) = " << linearregression(N,L1,0,2) << endl;
  cout << "  L1-norm (first 3) = " << linearregression(N,L1,0,3) << endl;
  cout << "  L1-norm (middle 2) = " << linearregression(N,L1,1,3) << endl;
  cout << "  L1-norm (last 3) = " << linearregression(N,L1,num-3,num) << endl;
  cout << "  L1-norm (last 2) = " << linearregression(N,L1,num-2,num) << endl;
  // L2 Norm.
  L2[0] = 0.0009675040;
  L2[1] = 0.0002379520;
  L2[2] = 0.0000816755;
  L2[3] = 0.0000644961;
  //L2[4] = 0;
  cout << "  L2-norm = " << linearregression(N,L2,0,num) << endl;
  cout << "  L2-norm (first 2) = " << linearregression(N,L2,0,2) << endl;
  cout << "  L2-norm (first 3) = " << linearregression(N,L2,0,3) << endl;
  cout << "  L2-norm (middle 2) = " << linearregression(N,L2,1,3) << endl;
  cout << "  L2-norm (last 3) = " << linearregression(N,L2,num-3,num) << endl;
  cout << "  L2-norm (last 2) = " << linearregression(N,L2,num-2,num) << endl;
  // Max Norm.
  Max[0] = 0.003866280;
  Max[1] = 0.001575550;
  Max[2] = 0.000959969;
  Max[3] = 0.005041300;
  //Max[4] = 0;
  cout << "  Max-norm = " << linearregression(N,Max,0,num) << endl;
  cout << "  Max-norm (first 2) = " << linearregression(N,Max,0,2) << endl;
  cout << "  Max-norm (first 3) = " << linearregression(N,Max,0,3) << endl;
  cout << "  Max-norm (middle 2) = " << linearregression(N,Max,1,3) << endl;
  cout << "  Max-norm (last 3) = " << linearregression(N,Max,num-3,num) << endl;
  cout << "  Max-norm (last 2) = " << linearregression(N,Max,num-2,num) << endl;
  cout << endl;

  ofstream output_file;    

  // Open the output data file.
  output_file.open("order.dat",ios::out);
  if (output_file.bad()) return 1;

  output_file.setf(ios::scientific);
  output_file << "TITLE = \": Ringleb's Flow Error Norms \"\n"
	      << "VARIABLES = \"N\" \\ \n"
	      << "\"L1\" \\ \n"
	      << "\"L2\" \\ \n"
	      << "\"Max\" \\ \n";
  output_file << "ZONE T =  \"Block Number = " << 0
	   << "\" \\ \n"
	   << "I = " << 5 << " \\ \n"
	   << "J = " << 1 << " \\ \n"
	   << "F = POINT \\ \n";

  output_file.setf(ios::scientific);
  for (int n = 0; n < num; n++) output_file << N[n] << " "
					    << L1[n] << " "
					    << L2[n] << " "
					    << Max[n] << endl;

  // Close the output data file.
  output_file.close();

  return 0;

}

double linearregression(double *xx, double *yy, int ni, int ne) {
  double xb, yb, xy, xb2, a, b, r, rx, ry, y1, yn;
  double x[ne-ni], y[ne-ni];
  for (int n = ni; n < ne; n++) {
    x[n-ni] = log(xx[n]);
    y[n-ni] = log(yy[n]);
  }
  xb = 0.0;
  yb = 0.0;
  for (int n = 0; n < ne-ni; n++) {
    xb += x[n]/double(ne-ni);
    yb += y[n]/double(ne-ni);
  }
  xy = 0.0;
  xb2 = 0.0;
  for (int n = 0; n < ne-ni; n++) {
    xy += x[n]*y[n];
    xb2 += x[n]*x[n];
  }
  xy = xy/double(ne-ni);
  xb2 = xb2/double(ne-ni);
  a = (yb*xb2 + xb*xy)/(xb2 - xb*xb);
  b = (xy - xb*yb)/(xb2 - xb*xb);
  r = 0.0;
  rx = 0.0;
  ry = 0.0;
  for (int n = 0; n < ne-ni; n++) {
    r += (x[n]-xb)*(y[n]-yb);
    rx += (x[n]-xb)*(x[n]-xb);
    ry += (y[n]-yb)*(y[n]-yb);
  }
  r /= (sqrt(rx)*sqrt(ry));

  y1 = b*x[0] + a;
  yn = b*x[ne-ni-1] + a;

 return (yn-y1)/(x[ne-ni-1]-x[0]);

}
