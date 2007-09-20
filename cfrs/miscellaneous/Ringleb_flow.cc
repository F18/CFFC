/* This program computes the Ringleb flow representation */

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

struct Vector2D {
  double x;
  double y;
};

int main(){
  // Variables
//   int NK = 31;
//   int NQ = 50;
  int NK = 32;
  int NQ = 50;
  Vector2D point[NK][NQ];
  Vector2D KQ[NK][NQ];
  double delta_q;
  double delta_k;
  double k_init, q_init, q_final;
  double q,k, c, J;
  double rho[NK][NQ];

  delta_k = (1.5-0.75)/(double)(NK-1);
  k_init = 0.75;
  q_init = 0.5;
  for (int i=0; i<= NK-1; i++){
    k= k_init + i*delta_k;
    q_final = k; // condition y = 0
    delta_q = (q_final - q_init)/(NQ-1);
    for (int j=0; j<= NQ-1; j++){
      q = q_init + j*delta_q;
      KQ[i][j].x = k;
      KQ[i][j].y = q;
    }
  }

  // file open
  ofstream output;
  output.open("Ringleb.dat");

  for (int i=0; i<=NK-1; i++)
    for (int j=0; j<=NQ-1; j++){
      c = sqrt(1 - 0.2*KQ[i][j].y*KQ[i][j].y);
      rho[i][j] = pow(c, 5);
      J = 1.0/c + 1.0/(3.0*pow(c,3)) + 1.0/(5.0*pow(c,5)) -
	  0.5*log((1.0+c)/(1.0-c));
      point[i][j].x = (0.5/rho[i][j])*(2.0/(KQ[i][j].x*KQ[i][j].x) - 
				       1.0/(KQ[i][j].y*KQ[i][j].y)) - 0.5*J;
      point[i][j].y = (1.0/(KQ[i][j].x*rho[i][j]*KQ[i][j].y))*
	              sqrt(1- (KQ[i][j].y*KQ[i][j].y)/(KQ[i][j].x*KQ[i][j].x));
    }

  //output zone
  output << "VARIABLES=\"x\"  \"y\"  \"rho\"" << endl;
  output << "ZONE I=" << NK << ", \t J="<< NQ <<", \t F=POINT" << endl;
  for (int j=0; j<=NQ-1; j++){
    for (int i=0; i<=NK-1; i++)
      output << point[i][j].x << "\t" << point[i][j].y << "\t"<< rho[i][j]<< endl;
  }
}
