/* CompCell1D.cc: Subroutines and data for 1D computational cells.*/

/* Include 1D Comp cell header file. */
#ifndef _COMPCELL1D_INCLUDED
#include "CompCell1D.h"
#endif //_COMPCELL1D_INCLUDED 

/************************************************************
 *        CompCell1D_Uniform Member Functions               *
 ***********************************************************/

void CompCell1D_Uniform::SetCellParameters(double loc, double delta,int nr, 
					   int ro){

  geom.x = loc;
  geom.dx = delta;
  U_cell = 0.0;
  SetDynamicVariable(nr, ro);
}

void CompCell1D_Uniform::SetCellParameters(double loc, double delta, double delta_sg,
					   int ro){
  
  geom.x = loc;
  geom.dx = delta;
  SetDynamicVariable(delta_sg, ro);
}

void CompCell1D_Uniform::SetDynamicVariable(int nr, int ro){
  
  double x_start;
  int nr_local;
  
  nr_local = nr;
  if (nr_local < 3)
    nr_local = 3;
  if (nr_local%2 != 1)
    {
      nr_local ++; //only odd numbers
    }
  N = nr_local;
  // allocate memory
  allocate(N,ro);
  
  // initialize subgrid geometry and solution
  delta_sg = geom.dx/(N-1);
  x_start = geom.x - geom.dx/2;
  for (int i=0; i<N; i++)
    {
      SubGrid[i] = x_start + i*delta_sg;
    }
  
  // initialize reconstruction coefficients
  for (int i=0; i<=RO; i++)
    D[i] = 0.0;
  SetSubDomainSolution();
}

void CompCell1D_Uniform::SetDynamicVariable(double delta, int ro){
  
  double x_start, ratio;
  double diff;
  
  ratio = geom.dx/delta;
  diff = ratio-(int)ratio;
  if (diff < 0.5)
    N = (int)ratio;
  else 
    N = (int)ratio + 1;
  
  if (N%2 != 1){
    N ++; //only odd numbers 
  }
  
  if (N < 3){
    N = 3;
    delta_sg = geom.dx/(N-1);
    cout << "The maximum delta_x for subgrid must be: " << delta_sg
	 << endl;
  }
  
  // allocate memory
  allocate(N,ro);
  
  // initialize subgrid geometry and solution
  delta_sg = geom.dx/(N-1);
  if (delta_sg != delta)
    cout << "The new value for delta_x is: " << delta_sg << endl;
  
  x_start = geom.x - geom.dx/2;
  for (int i=0; i<=N-1; i++ )
    {
      SubGrid[i] = x_start + i*delta_sg;
    }
 
  // initialize reconstruction coefficients
  D[0] = U_cell;
  for (int i=1; i<=RO; i++)
    D[i] = 0.0;

  SetSubDomainSolution();
}

/****************************************************************
 * SetSubDomainSolution subrutine                               *
 *                                                              *
 * Truncated Taylor series expansion reconstruction for order k *
 *           k                                                  *
 *   U(x) = SUM alpha[p]*D[p]*(x - x_j)^p                       *
 *          p=0                                                 *
 **************************************************************/

void CompCell1D_Uniform::SetSubDomainSolution(){
  
  double diff;
  
  for (int i=0; i<=N-1; i++){
    diff = (SubGrid[i] - geom.x);
    U_points[i] = D[0];
    for (int j=1; j<=RO; j++)
      U_points[i] += pow(diff,j)*D[j];
  }
}

/************************************************************
 *********                                      *************
 *******  CompCell1D_NonUniform Member Functions  ***********
 *********                                      *************
 ***********************************************************/

void CompCell1D_NonUniform::SetCellParameters(double loc, double delta,int nr, 
					      int ro){

  geom.x = loc;
  geom.dx = delta;
  U_cell = 0.0;
  ErrorRec = 0.0;
  SetDynamicVariable(nr, ro);
}

void CompCell1D_NonUniform::SetCellParameters(double loc, double delta, double delta_sg,
					   int ro){
  
  geom.x = loc;
  geom.dx = delta;
  ErrorRec = 0.0;
  SetDynamicVariable(delta_sg, ro);
}

void CompCell1D_NonUniform::SetDynamicVariable(int nr, int ro){
  
  double x_start;
  int nr_local;

  nr_local = nr;
  if (nr_local < 3)
    nr_local = 3;
  if (nr_local%2 != 1)
    nr_local ++; //only odd numbers

  N = nr_local;
  // allocate memory
  allocate(N,ro);

  // initialize subgrid geometry and solution
  delta_sg = geom.dx/(N-1);
  x_start = geom.x - geom.dx/2;
  for (int i=0; i<N; i++){
    SubGrid[i] = x_start + i*delta_sg;
  }

  // initialize reconstruction coefficients
  D[0] = U_cell;
  for (int i=1; i<=RO; i++)
    D[i] = 0.0;

  SetSubDomainSolution();
}

void CompCell1D_NonUniform::SetDynamicVariable(double delta, int ro){
  
  double x_start, ratio;
  double diff;
  
  ratio = geom.dx/delta;
  diff = ratio-(int)ratio;
  if (diff < 0.5)
    N = (int)ratio;
  else 
    N = (int)ratio + 1;
  
  if (N%2 != 1){
    N ++; //only odd numbers 
  }
  
  if (N < 3){
    N = 3;
    delta_sg = geom.dx/(N-1);
    cout << "The minimum delta_x for subgrid must be: " << delta_sg
	 << endl;
  }
  
  // allocate memory
  allocate(N,ro);
  
  // initialize subgrid geometry and solution
  delta_sg = geom.dx/(N-1);
  if (delta_sg != delta)
    cout << "The new value for delta_x is: " << delta_sg << endl;
  
  x_start = geom.x - geom.dx/2;
  for (int i=0; i<=N-1; i++ )
    {
      SubGrid[i] = x_start + i*delta_sg;
    }
 
  // initialize reconstruction coefficients
  for (int i=0; i<=RO; i++)
    D[i] = 0.0;

  SetSubDomainSolution();
}

/****************************************************************
 * SetSubDomainSolution subrutine                               *
 *                                                              *
 * Truncated Taylor series expansion reconstruction for order k *
 *           k                                                  *
 *   U(x) = SUM alpha[p]*D[p]*(x - x_j)^p                       *
 *          p=0                                                 *
 **************************************************************/

void CompCell1D_NonUniform::SetSubDomainSolution(){
  
  double diff;

  for (int i=0; i<=N-1; i++){
    diff = SubGrid[i] - geom.x;
    U_points[i] = D[0];
    for (int j=1; j<=RO; j++){
      U_points[i] += pow(diff,j)*D[j];
    }
  }
}

void CompCell1D_NonUniform::SetSubDomainSolution(CompCell1D_NonUniform &Cell){
  
  double diff;
  
  for (int i=0; i<=N-1; i++){
    diff = SubGrid[i] - Cell.geom.x;
    U_points[i] = D[0];
    for (int j=1; j<=RO; j++)
      U_points[i] += pow(diff,j)*D[j];
  }
}

void CompCell1D_NonUniform::SetDerivative_Representation(int order){
  // order -- the order of the derivative

  double diff;
  double sum;
  int val;

  for (int i=0; i<=N-1; i++){
    if (order > RO)
      U_deriv_num[i]= 0;
    else {
      diff = SubGrid[i] - geom.x;  
      sum = 0.0;
      for (int p=order; p<=RO; p++){
	val = p-order;
	sum += D[p]*pow(diff,val);
      }
      U_deriv_num[i] = sum;
    }
  }
}

double CompCell1D_NonUniform::SolutionAtCoordinates(double x){

  // computes the value of the function in the point "x", using the
  // coefficients determined during the reconstruction 
  double diff = x - geom.x;
  double FuncValue = 0.0;
  for (int j=0; j<=RO; j++)
    FuncValue += pow(diff,j)*D[j];
  return FuncValue;
}

void CompCell1D_NonUniform::ComputeReconstructionError (FunctionType1D ExactFunction){

  _Error_<FunctionType1D,CompCell1D_NonUniform,double> ErrorFunction(ExactFunction,this);

  double xStart = geom.x-0.5*geom.dx;
  double xStop = geom.x+0.5*geom.dx;

  long double DummyParam = 0.0;

  ErrorRec = AdaptiveGaussianQuadrature(ErrorFunction,xStart,xStop,14,DummyParam)/geom.dx;

}
