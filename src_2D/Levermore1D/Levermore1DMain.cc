#include<iostream>
#include"Levermore1DState.h"

using namespace std;

#define NUMBER 10

int main() {

  Levermore1D_pState<NUMBER> W;

  W[1] = 3.0;
  W[2] = 1.2;
  W[3] = 1000.0;
  W[4] = 1000000.0;
  W[5] = 1000000000.0;

  Levermore1D_cState<NUMBER> U(W);

  cout.precision(40);

  cout << W << endl << U << endl;
  cout << Levermore1D_pState<NUMBER>(U) << endl;

  return 0;
}
