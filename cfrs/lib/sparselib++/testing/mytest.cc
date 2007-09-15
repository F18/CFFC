#include <cstdlib>

using namespace std;

#include "compcol_double.h"
#include "comprow_double.h"
#include "ilupre_double.h"
#include "diagpre_double.h"
#include "icpre_double.h"
#include "iohb_double.h"
#include "iotext_double.h"

#include "mvv.h" 
#include "mvm.h"

int
main(int argc, char * argv[])
{

  int i, i_nonzero, i_row;

  double aa;
  MV_Vector_double GG;
  MV_Vector_double VV(10,1);
  MV_ColMat_double MM(10, 7, 0);
  VV[3]=0.00;
  VV(3) = 2.00;
  MM(5,1) = 10.00;
  for (i = 0; i < 7; ++i) {
    MM(i, i) = 1.00;
  }
  aa = dot(VV, VV);
  cout << "========================\n";
  cout << VV << "\n";
  cout << MM << "\n";
  cout << aa << "\n";
  cout << MM.size(0) << " " << MM.size(1) << " " << MM.lda() << "\n";
  cout << MM(2, MV_VecIndex(0, 6)) << "\n";
  cout << "========================\n";
  cout.flush();

  CompCol_Mat_double Ac;
  CompRow_Mat_double Ar;

  // Define nonzero entries of matrix.
  double val[64] = {-4.0,  1.0,  1.0,             // row  0:  0- 2
		     1.0, -4.0,  1.0,  1.0,       // row  1:  3- 6
		     1.0, -4.0,  1.0,  1.0,       // row  2:  7-10
		     1.0, -4.0,  1.0,             // row  3: 11-13
		     1.0, -4.0,  1.0,  1.0,       // row  4: 14-17
		     1.0,  1.0, -4.0,  1.0,  1.0, // row  5: 18-22
		     1.0,  1.0, -4.0,  1.0,  1.0, // row  6: 23-27
		     1.0,  1.0, -4.0,  1.0,       // row  7: 28-31
		     1.0, -4.0,  1.0,  1.0,       // row  8: 32-35
		     1.0,  1.0, -4.0,  1.0,  1.0, // row  9: 36-40
		     1.0,  1.0, -4.0,  1.0,  1.0, // row 10: 41-45
		     1.0,  1.0, -4.0,  1.0,       // row 11: 46-49
		     1.0, -4.0,  1.0,             // row 12: 50-52
	             1.0,  1.0, -4.0,  1.0,       // row 13: 53-56
	             1.0,  1.0, -4.0,  1.0,       // row 14: 57-60
	             1.0,  1.0, -4.0};            // row 15: 61-63
  int colind[64] = {  0,  1,  4,                  // row  0:  0- 2 
                      0,  1,  2,  5,              // row  1:  3- 6
                      1,  2,  3,  6,              // row  2:  7-10
		      2,  3,  7,                  // row  3: 11-13
                      0,  4,  5,  8,              // row  4: 14-17
                      1,  4,  5,  6,  9,          // row  5: 18-22
                      2,  5,  6,  7, 10,          // row  6: 23-27
                      3,  6,  7, 11,              // row  7: 28-31
                      4,  8,  9, 12,              // row  8: 32-35
                      5,  8,  9, 10, 13,          // row  9: 36-40
                      6,  9, 10, 11, 14,          // row 10: 41-45
                      7, 10, 11, 15,              // row 11: 46-49
                      8, 12, 13,                  // row 12: 50-52
                      9, 12, 13, 14,              // row 13: 53-56
                     10, 13, 14, 15,              // row 14: 57-60
                     11, 14, 15};                 // row 15: 61-63
  int rowptr[17] = { 0,  // row  0
                     3,  // row  1
                     7,  // row  2
                    11,  // row  3
                    14,  // row  4
                    18,  // row  5
                    23,  // row  6
                    28,  // row  7
                    32,  // row  8
                    36,  // row  9
                    41,  // row 10
                    46,  // row 11
                    50,  // row 12
                    53,  // row 13
                    57,  // row 14
		    61,  // row 15
                    64}; // row 16

  Ar = CompRow_Mat_double(16, 16, 64, val, rowptr, colind);
//    cout << Ar;
//    writetxtfile_mat("dummy", Ar);

  CompRow_ILUPreconditioner_double Mr(Ar);
  
//    int n = Ar.dim(0);
//    cout << n << "\n";
//    MV_Vector_double b(n), x(n);
//    b(MV_VecIndex(0, n-1)) = 0.00;
//    b(15) = 1.00;
//    x = Mr.solve(b);
//    cout << x(MV_VecIndex(0, n-1));  

  CompCol_Mat_double Bc;
  CompRow_Mat_double Br;

  // Define nonzero entries of matrix.

  double b_val[17] = {-2.0,  1.0,  2.0,             // row  0:  0- 2
		       1.0, -2.0,  1.0,  2.0,       // row  1:  3- 6
		       1.0, -2.0,  1.0,             // row  3:  7- 9
  		       1.0,  1.0, -2.0,  1.0,       // row  4: 10-13
		       1.0,  1.0, -2.0};            // row  5: 14-16
  int b_colind[17] = {  0,  1,  3,                  // row  0:  0- 2 
                        0,  1,  2,  4,              // row  1:  3- 6
                        1,  2,  3,                  // row  2:  7- 9
		        0,  2,  3,  4,              // row  3: 10-13
                        1,  3,  4};                 // row  4: 14-16
  int b_rowptr[6] = { 0,  // row  0
                      3,  // row  1
                      7,  // row  2
                     10,  // row  3
                     14,  // row  4
                     17}; // row  5

  Br = CompRow_Mat_double(5, 5, 17, b_val, b_rowptr, b_colind);
  cout << "B======\n";
  cout << Br;

  CompRow_ILUPreconditioner_double MBr(Br);
  
  cout << "L======\n";
  cout << MBr.l();
  cout << "U======\n";
  cout << MBr.u();

  int m = Br.dim(0);
  cout << m << "\n";
  MV_Vector_double bb(m), xb(m);
  bb(MV_VecIndex(0, m-1)) = 0.00;
  bb(1) = 1.00;
  xb = MBr.solve(bb);
  cout << xb(MV_VecIndex(0, m-1));

  DiagPreconditioner_double DBr(Br);
  cout << DBr.diag(0) << " " << DBr.diag(1) << "\n";

  return(0);
}
