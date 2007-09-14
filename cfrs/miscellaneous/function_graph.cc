

#include <iostream>
#include <fstream.h>
#include "TestFunctions.h"

int main(){

  ofstream out;
  out.open("function.dat");

  int N;

  cout << "N= ";
  cin >> N;
  cout << endl;

  double *x;
  double x_inf, x_sup;

  x = new double [N];

  x_inf = -1.0;
  x_sup = 1.0;

  double delta;
  double *val;
  val = new double [N];
  delta = (x_sup-x_inf)/(N-1);
  for (int i=0; i<=N-1; i++)
    {
      x[i] = x_inf + i*delta;
      val[i] = Test_Default1D(x[i]);
    }

  out << "VARIABLES=\"x\" \"F\"" << endl << "ZONE" << endl;

  for (int j=0; j<= N-1; j++){
    out << x[j] <<"\t" << val[j] <<endl;
  }

  delete [] x;
  out.close();

  return 0;

}
