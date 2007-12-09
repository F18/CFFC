/**********************************************************************
 * EllipticOperatorAnalysis2D.h: Header file defining 2D elliptic     *
 *                               operator analysis classes.           *
 **********************************************************************/

#ifndef _ELLIPTIC2D_INCLUDED
#define _ELLIPTIC2D_INCLUDED

// Include required C++ libraries.
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include math macro header file.

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

// Include the CFD header file.

#ifndef _CFD_INCLUDED
#include "CFD.h"
#endif // _CFD_INCLUDED

// Include the 2D vector header file.

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

// Include the 2D tensor header file.

#ifndef _TENSOR2D_INCLUDED
#include "../Math/Tensor2D.h"
#endif //_TENSOR2D_INCLUDED

/*!
 * Class: EllipticOperatorCoefficients2D
 *
 * Class storing the coefficients of the least-squares gradient
 * reconstruction.
 *
 */
class EllipticOperatorCoefficients2D {
private:
public:

  //! Stencil size.
  int nn;

  //! Stencil mid-point.
  int nc;

  //! Array of coefficients.
  double  **u;

  //! Creation constructor.
  EllipticOperatorCoefficients2D(void) { nn = 0; nc = 0; u = NULL; }

  //! Creation constructor.
  EllipticOperatorCoefficients2D(const int &N) { u = NULL; allocate(N); }

  //! Copy constructor.
  EllipticOperatorCoefficients2D(const EllipticOperatorCoefficients2D &U) { copy(U); }

  //! Destructor.
  ~EllipticOperatorCoefficients2D(void) { deallocate(); };

  //! Copy function.
  void copy(const EllipticOperatorCoefficients2D &U) {
    if (u == NULL) allocate(U.nn);
    for (int j = 0; j < nn; j++) for (int i = 0; i < nn; i++) u[i][j] = U.u[i][j];
  }

  //! Allocation function.
  void allocate(const int &N) {
    assert(N >= 3 && N%2 == 1);
    if (u != NULL) deallocate();
    nn = N; nc = (N-1)/2;
    u = new double*[nn];
    for (int i = 0; i < nn; i++) {
      u[i] = new double[nn];
      for (int j = 0; j < nn; j++) u[i][j] = ZERO;
    }
  }

  //! Deallocation function.
  void deallocate(void) {
    for (int i = 0; i < nn; i++) { delete []u[i]; u[i] = NULL; }
    delete []u; u = NULL;
    nn = 0; nc = 0;
  }

  //! Zero function.
  void Zero(void) {
    for (int j = 0; j < nn; j++) for (int i = 0; i < nn; i++) u[i][j] = ZERO;
  }

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const EllipticOperatorCoefficients2D &U) { return out_file; }
  friend istream &operator >> (istream &in_file, EllipticOperatorCoefficients2D &U) { return in_file; }
  //@}

};

/*****************************************************************************
 * Routine: Elliptic_Operator_Analysis_Average_Gradient_Linear_Least_Squares *
 *                                                                           *
 *****************************************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Elliptic_Operator_Analysis_Average_Gradient_Linear_Least_Squares(Quad_Soln_Block *Local_SolnBlk,
 								     AdaptiveBlock2D_List &Local_Solution_Block_List,
 								     Quad_Soln_Input_Parameters &IP,
 								     const int &correction_flag) {

  // Exit immediately if no solution blocks are being used.
  if (!Local_Solution_Block_List.Nused()) return 0;

  ///////////////////////////////////////////
  // Allocate reconstruction coefficients. //
  ///////////////////////////////////////////

  int n_gradient_stencil = 3, n_laplacian_stencil = 5;
  EllipticOperatorCoefficients2D **Gradient_dudx, **Gradient_dudy;
  EllipticOperatorCoefficients2D **Laplacian;

  Gradient_dudx = new EllipticOperatorCoefficients2D*[Local_SolnBlk[0].NCi];
  Gradient_dudy = new EllipticOperatorCoefficients2D*[Local_SolnBlk[0].NCi];
  Laplacian = new EllipticOperatorCoefficients2D*[Local_SolnBlk[0].NCi];
  for (int i = 0; i < Local_SolnBlk[0].NCi; i++) {
    Gradient_dudx[i] = new EllipticOperatorCoefficients2D[Local_SolnBlk[0].NCj];
    Gradient_dudy[i] = new EllipticOperatorCoefficients2D[Local_SolnBlk[0].NCj];
    Laplacian[i] = new EllipticOperatorCoefficients2D[Local_SolnBlk[0].NCj];
    for (int j = 0; j < Local_SolnBlk[0].NCj; j++) {
      Gradient_dudx[i][j].allocate(n_gradient_stencil);
      Gradient_dudy[i][j].allocate(n_gradient_stencil);
      Laplacian[i][j].allocate(n_laplacian_stencil);
    }
  }

  ///////////////////////////////////
  // Prepare the output data file. //
  ///////////////////////////////////

  int name;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream out;    

  // Determine prefix of output data file names.
  name = 0;
  while (1) {
    if (IP.Output_File_Name[name] == ' ' || IP.Output_File_Name[name] == '.') break;
    prefix[name] = IP.Output_File_Name[name];
    name++;
    if (name > strlen(IP.Output_File_Name)) break;
  }
  prefix[name] = '\0';
  if (!correction_flag) strcat(prefix,"_eoa_aglls_cpu");
  else strcat(prefix,"_eoa_caglls_cpu");

  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Local_Solution_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data file.
  out.open(output_file_name_ptr,ios::out);
  if (out.fail()) return 1;

  /////////////////////////////////////////////////////////////////
  // Compute the reconstruction coefficients for all of the used //
  // solution blocks and output the data to a file.              //
  /////////////////////////////////////////////////////////////////

  // Declare gradient and Laplacian reconstruction variables:
  double dl, det, xhat, yhat, xhat2, xhatyhat, yhat2;
  Vector2D dx, nhatN, nhatS, nhatE, nhatW, thatN, thatS, thatE, thatW;
  Vector2D nfaceN, nfaceS, nfaceE, nfaceW;

  // Perform the analysis on each solution block that is in use.
  for (int nb = 0; nb < Local_Solution_Block_List.Nblk; nb++) {
    if (Local_Solution_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Initialize the coefficients.
      for (int j = Local_SolnBlk[nb].JCl-1; j <= Local_SolnBlk[nb].JCu+1; j++) {
 	for (int i = Local_SolnBlk[nb].ICl-1; i <= Local_SolnBlk[nb].ICu+1; i++) {
 	  Gradient_dudx[i][j].Zero();
 	  Gradient_dudy[i][j].Zero();
 	  Laplacian[i][j].Zero();
 	}
      }

      // Compute the coefficients for the cell-centred least-squares
      // gradient reconstruction.
      for (int j = Local_SolnBlk[nb].JCl-1; j <= Local_SolnBlk[nb].JCu+1; j++) {
 	for (int i = Local_SolnBlk[nb].ICl-1; i <= Local_SolnBlk[nb].ICu+1; i++) {
 	  // Construct coefficients of the normal matrix for the current cell.
	  xhat2 = ZERO; xhatyhat = ZERO; yhat2 = ZERO;
 	  for (int jj = j-1; jj <= j+1; jj++) {
 	    for (int ii = i-1; ii <= i+1; ii++) {
	      //if (ii == i || jj == j) {
 	      xhat = Local_SolnBlk[nb].Grid.Cell[ii][jj].Xc.x - Local_SolnBlk[nb].Grid.Cell[i][j].Xc.x;
 	      yhat = Local_SolnBlk[nb].Grid.Cell[ii][jj].Xc.y - Local_SolnBlk[nb].Grid.Cell[i][j].Xc.y;
	      xhat2 += xhat*xhat; xhatyhat += xhat*yhat; yhat2 += yhat*yhat;
	      //}
 	    }
 	  }
	  // Determine the determinant of the coefficient matrix.
	  det = xhat2*yhat2 - xhatyhat*xhatyhat;
 	  // Determine the least squares coefficients.
 	  for (int jj = j-1; jj <= j+1; jj++) {
 	    for (int ii = i-1; ii <= i+1; ii++) {
	      if (!(ii == i && jj == j)) {
		xhat = Local_SolnBlk[nb].Grid.Cell[ii][jj].Xc.x - Local_SolnBlk[nb].Grid.Cell[i][j].Xc.x;
		yhat = Local_SolnBlk[nb].Grid.Cell[ii][jj].Xc.y - Local_SolnBlk[nb].Grid.Cell[i][j].Xc.y;
		//if (ii == i || jj == j) {
  		Gradient_dudx[i][j].u[ii-i+1][jj-j+1] = (xhat*yhat2 - yhat*xhatyhat)/det;
   		Gradient_dudy[i][j].u[ii-i+1][jj-j+1] = (yhat*xhat2 - xhat*xhatyhat)/det;
  		Gradient_dudx[i][j].u[1][1] -= (xhat*yhat2 - yhat*xhatyhat)/det;
  		Gradient_dudy[i][j].u[1][1] -= (yhat*xhat2 - xhat*xhatyhat)/det;
		//}
 	      }
 	    }
 	  }
	}
      }

      // Compute the coefficients for the Laplacian reconstruction.
      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
 	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	  // Summarize the cell face normals (not unit normals).
	  nfaceN = Local_SolnBlk[nb].Grid.nfaceN(i,j)*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	  nfaceS = Local_SolnBlk[nb].Grid.nfaceS(i,j)*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	  nfaceE = Local_SolnBlk[nb].Grid.nfaceE(i,j)*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	  nfaceW = Local_SolnBlk[nb].Grid.nfaceW(i,j)*Local_SolnBlk[nb].Grid.lfaceW(i,j);
 	  for (int jj = j-1; jj <= j+1; jj++) {
 	    for (int ii = i-1; ii <= i+1; ii++) {
	      // NORTH face:
	      if (Local_SolnBlk[nb].Grid.lfaceN(i,j) > ZERO) {
		Laplacian[i][j].u[ii-i+2][jj-j+3] += HALF*(Gradient_dudx[i][j+1].u[ii-i+1][jj-j+1]*nfaceN.x +
							   Gradient_dudy[i][j+1].u[ii-i+1][jj-j+1]*nfaceN.y);
		Laplacian[i][j].u[ii-i+2][jj-j+2] += HALF*(Gradient_dudx[i][j  ].u[ii-i+1][jj-j+1]*nfaceN.x +
							   Gradient_dudy[i][j  ].u[ii-i+1][jj-j+1]*nfaceN.y);
	      }
	      // SOUTH face:
	      if (Local_SolnBlk[nb].Grid.lfaceS(i,j) > ZERO) {
		Laplacian[i][j].u[ii-i+2][jj-j+1] += HALF*(Gradient_dudx[i][j-1].u[ii-i+1][jj-j+1]*nfaceS.x +
							   Gradient_dudy[i][j-1].u[ii-i+1][jj-j+1]*nfaceS.y);
		Laplacian[i][j].u[ii-i+2][jj-j+2] += HALF*(Gradient_dudx[i][j  ].u[ii-i+1][jj-j+1]*nfaceS.x +
							   Gradient_dudy[i][j  ].u[ii-i+1][jj-j+1]*nfaceS.y);
	      }
	      // EAST face:
	      if (Local_SolnBlk[nb].Grid.lfaceE(i,j) > ZERO) {
		Laplacian[i][j].u[ii-i+3][jj-j+2] += HALF*(Gradient_dudx[i+1][j].u[ii-i+1][jj-j+1]*nfaceE.x +
							   Gradient_dudy[i+1][j].u[ii-i+1][jj-j+1]*nfaceE.y);
		Laplacian[i][j].u[ii-i+2][jj-j+2] += HALF*(Gradient_dudx[i  ][j].u[ii-i+1][jj-j+1]*nfaceE.x +
							   Gradient_dudy[i  ][j].u[ii-i+1][jj-j+1]*nfaceE.y);
	      }
	      // WEST face:
	      if (Local_SolnBlk[nb].Grid.lfaceW(i,j) > ZERO) {
		Laplacian[i][j].u[ii-i+1][jj-j+2] += HALF*(Gradient_dudx[i-1][j].u[ii-i+1][jj-j+1]*nfaceW.x +
							   Gradient_dudy[i-1][j].u[ii-i+1][jj-j+1]*nfaceW.y);
		Laplacian[i][j].u[ii-i+2][jj-j+2] += HALF*(Gradient_dudx[i  ][j].u[ii-i+1][jj-j+1]*nfaceW.x +
							   Gradient_dudy[i  ][j].u[ii-i+1][jj-j+1]*nfaceW.y);
	      }
	    }
	  }
	  // Add the correction term if required.
	  if (correction_flag) {
	    // NORTH face:
	    if (Local_SolnBlk[nb].Grid.lfaceN(i,j) > ZERO) {
	      dx = Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      dl = dx.abs();
	      thatN = dx/dl;
	      nhatN = Local_SolnBlk[nb].Grid.nfaceN(i,j);
	      Laplacian[i][j].u[2][3] += (dot(nhatN,nfaceN)/dot(nhatN,thatN))/dl;
	      Laplacian[i][j].u[2][2] -= (dot(nhatN,nfaceN)/dot(nhatN,thatN))/dl;
	      for (int jj = j-1; jj <= j+1; jj++) {
		for (int ii = i-1; ii <= i+1; ii++) {
		  Laplacian[i][j].u[ii-i+2][jj-j+3] -= HALF*(Gradient_dudx[i][j+1].u[ii-i+1][jj-j+1]*thatN.x +
							     Gradient_dudy[i][j+1].u[ii-i+1][jj-j+1]*thatN.y)*dot(nhatN,nfaceN)/dot(nhatN,thatN);
		  Laplacian[i][j].u[ii-i+2][jj-j+2] -= HALF*(Gradient_dudx[i][j  ].u[ii-i+1][jj-j+1]*thatN.x + 
							     Gradient_dudy[i][j  ].u[ii-i+1][jj-j+1]*thatN.y)*dot(nhatN,nfaceN)/dot(nhatN,thatN);
		}
	      }
	    }
	    // SOUTH face:
	    if (Local_SolnBlk[nb].Grid.lfaceS(i,j) > ZERO) {
	      dx = Local_SolnBlk[nb].Grid.Cell[i][j-1].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      dl = dx.abs();
	      thatS = dx/dl;
	      nhatS = Local_SolnBlk[nb].Grid.nfaceS(i,j);
	      Laplacian[i][j].u[2][1] += (dot(nhatS,nfaceS)/dot(nhatS,thatS))/dl;
	      Laplacian[i][j].u[2][2] -= (dot(nhatS,nfaceS)/dot(nhatS,thatS))/dl;
	      for (int jj = j-1; jj <= j+1; jj++) {
		for (int ii = i-1; ii <= i+1; ii++) {
		  Laplacian[i][j].u[ii-i+2][jj-j+1] -= HALF*(Gradient_dudx[i][j-1].u[ii-i+1][jj-j+1]*thatS.x +
							     Gradient_dudy[i][j-1].u[ii-i+1][jj-j+1]*thatS.y)*dot(nhatS,nfaceS)/dot(nhatS,thatS);
		  Laplacian[i][j].u[ii-i+2][jj-j+2] -= HALF*(Gradient_dudx[i][j  ].u[ii-i+1][jj-j+1]*thatS.x +
							     Gradient_dudy[i][j  ].u[ii-i+1][jj-j+1]*thatS.y)*dot(nhatS,nfaceS)/dot(nhatS,thatS);
		}
	      }
	    }
	    // EAST face:
	    if (Local_SolnBlk[nb].Grid.lfaceE(i,j) > ZERO) {
	      dx = Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      dl = dx.abs();
	      thatE = dx/dl;
	      nhatE = Local_SolnBlk[nb].Grid.nfaceE(i,j);
	      Laplacian[i][j].u[3][2] += (dot(nhatE,nfaceE)/dot(nhatE,thatE))/dl;
	      Laplacian[i][j].u[2][2] -= (dot(nhatE,nfaceE)/dot(nhatE,thatE))/dl;
	      for (int jj = j-1; jj <= j+1; jj++) {
		for (int ii = i-1; ii <= i+1; ii++) {
		  Laplacian[i][j].u[ii-i+3][jj-j+2] -= HALF*(Gradient_dudx[i+1][j].u[ii-i+1][jj-j+1]*thatE.x +
							     Gradient_dudy[i+1][j].u[ii-i+1][jj-j+1]*thatE.y)*dot(nhatE,nfaceE)/dot(nhatE,thatE);
		  Laplacian[i][j].u[ii-i+2][jj-j+2] -= HALF*(Gradient_dudx[i  ][j].u[ii-i+1][jj-j+1]*thatE.x +
							     Gradient_dudy[i  ][j].u[ii-i+1][jj-j+1]*thatE.y)*dot(nhatE,nfaceE)/dot(nhatE,thatE);
		}
	      }
	    }
	    // WEST face:
	    if (Local_SolnBlk[nb].Grid.lfaceW(i,j) > ZERO) {
	      dx = Local_SolnBlk[nb].Grid.Cell[i-1][j].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      dl = dx.abs();
	      thatW = dx/dl;
	      nhatW = Local_SolnBlk[nb].Grid.nfaceW(i,j);
	      Laplacian[i][j].u[1][2] += (dot(nhatW,nfaceW)/dot(nhatW,thatW))/dl;
	      Laplacian[i][j].u[2][2] -= (dot(nhatW,nfaceW)/dot(nhatW,thatW))/dl;
	      for (int jj = j-1; jj <= j+1; jj++) {
		for (int ii = i-1; ii <= i+1; ii++) {
		  Laplacian[i][j].u[ii-i+1][jj-j+2] -= HALF*(Gradient_dudx[i-1][j].u[ii-i+1][jj-j+1]*thatW.x +
							     Gradient_dudy[i-1][j].u[ii-i+1][jj-j+1]*thatW.y)*dot(nhatW,nfaceW)/dot(nhatW,thatW);
		  Laplacian[i][j].u[ii-i+2][jj-j+2] -= HALF*(Gradient_dudx[i  ][j].u[ii-i+1][jj-j+1]*thatW.x +
							     Gradient_dudy[i  ][j].u[ii-i+1][jj-j+1]*thatW.y)*dot(nhatW,nfaceW)/dot(nhatW,thatW);
		}
	      }
	    }
	  }
	  // Divide by the cell area to complete stencil construction.
	  for (int m = 0; m < n_laplacian_stencil; m++) {
	    for (int n = 0; n < n_laplacian_stencil; n++) {
	      Laplacian[i][j].u[n][m] /= Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    }
	  }
 	}
      }

      // Output the coefficients of the Laplacian reconstruction to data files.
      out << setprecision(8);
      if (nb == 0) {
	out << " ELLIPTIC OPERATOR ANALYSIS::";
	if (!correction_flag) out << "AVERAGE CELL-CENTRED GRADIENT WITH LINEAR LEAST SQUARES RECONSTRUCTION";
	else out << "CORRECTED AVERAGE CELL-CENTRED GRADIENT WITH LINEAR LEAST SQUARES RECONSTRUCTION";
	out << endl;
      }
      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
 	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	  out << endl;
	  out << "========================================================================";
	  out << endl << endl << " Xc(" << i << "," << j << ") = ("
	      << Local_SolnBlk[nb].Grid.Cell[i][j].Xc.x << ","
	      << Local_SolnBlk[nb].Grid.Cell[i][j].Xc.y << ")";
	  out << endl << endl << " -> Reconstruction weights: ";
	  for (int m = n_laplacian_stencil-1; m >= 0; m--) {
	    if (m != n_laplacian_stencil-1) out << "                            ";
	    for (int n = 0; n < n_laplacian_stencil; n++) {
	      if (fabs(Laplacian[i][j].u[n][m]/Laplacian[i][j].u[2][2]) < 1.0e-12) Laplacian[i][j].u[n][m] = ZERO;
	      out << " " << -Laplacian[i][j].u[n][m]/Laplacian[i][j].u[2][2];
	    }
	    out << endl;
	  }
	  out << endl << " -> Taylor series coefficients:";
	  Output_Taylor_Series_Coefficients(Local_SolnBlk[nb],i,j,Laplacian[i][j],out);
	  out << endl << endl << " -> Raw coefficients: ";
	  for (int m = n_laplacian_stencil-1; m >= 0; m--) {
	    if (m != n_laplacian_stencil-1) out << "                      ";
	    for (int n = 0; n < n_laplacian_stencil; n++) {
	      out << " " << Laplacian[i][j].u[n][m];
	    }
	    out << endl;
	  }
	}
      }

    }
  }

  // Close the output data file.
  out.close();

  /////////////////////////////////////////////
  // Deallocate reconstruction coefficients. //
  /////////////////////////////////////////////

  for (int i = 0; i < Local_SolnBlk[0].NCi; i++) {
    for (int j = 0; j < Local_SolnBlk[0].NCj; j++) {
      Gradient_dudx[i][j].deallocate();
      Gradient_dudy[i][j].deallocate();
      Laplacian[i][j].deallocate();
    }
    delete []Gradient_dudx[i]; Gradient_dudx[i] = NULL;
    delete []Gradient_dudy[i]; Gradient_dudy[i] = NULL;
    delete []Laplacian[i]; Laplacian[i] = NULL;
  }
  delete []Gradient_dudx; Gradient_dudx = NULL;
  delete []Gradient_dudy; Gradient_dudy = NULL;
  delete []Laplacian; Laplacian = NULL;

  // Average gradient with least-squares coefficients calculated
  // successfully.
  return 0;

}


/**********************************************************************
 * Routine: Elliptic_Operator_Analysis_Directional_Derivative         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Elliptic_Operator_Analysis_Directional_Derivative(Quad_Soln_Block *Local_SolnBlk,
 						      AdaptiveBlock2D_List &Local_Solution_Block_List,
 						      Quad_Soln_Input_Parameters &IP) {

  // Exit immediately if no solution blocks are being used.
  if (!Local_Solution_Block_List.Nused()) return 0;

  ///////////////////////////////////////////
  // Allocate reconstruction coefficients. //
  ///////////////////////////////////////////

  int n_laplacian_stencil = 3;
  EllipticOperatorCoefficients2D **Laplacian;

  Laplacian = new EllipticOperatorCoefficients2D*[Local_SolnBlk[0].NCi];
  for (int i = 0; i < Local_SolnBlk[0].NCi; i++) {
    Laplacian[i] = new EllipticOperatorCoefficients2D[Local_SolnBlk[0].NCj];
    for (int j = 0; j < Local_SolnBlk[0].NCj; j++) {
      Laplacian[i][j].allocate(n_laplacian_stencil);
    }
  }

  ///////////////////////////////////
  // Prepare the output data file. //
  ///////////////////////////////////

  int name;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream out;    

  // Determine prefix of output data file names.
  name = 0;
  while (1) {
    if (IP.Output_File_Name[name] == ' ' || IP.Output_File_Name[name] == '.') break;
    prefix[name] = IP.Output_File_Name[name];
    name++;
    if (name > strlen(IP.Output_File_Name)) break;
  }
  prefix[name] = '\0';
  strcat(prefix,"_eoa_dd_cpu");

  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Local_Solution_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data file.
  out.open(output_file_name_ptr,ios::out);
  if (out.fail()) return 1;

  /////////////////////////////////////////////////////////////////
  // Compute the reconstruction coefficients for all of the used //
  // solution blocks and output the data to a file.              //
  /////////////////////////////////////////////////////////////////

  // Declare gradient and Laplacian reconstruction variables:
  double dl;
  Vector2D dx, nhatN, nhatS, nhatE, nhatW, thatN, thatS, thatE, thatW;
  Vector2D nfaceN, nfaceS, nfaceE, nfaceW;

  // Perform the analysis on each solution block that is in use.
  for (int nb = 0; nb < Local_Solution_Block_List.Nblk; nb++) {
    if (Local_Solution_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Compute the coefficients for the Laplacian reconstruction.
      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
 	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	  // Initialize the coefficients:
 	  Laplacian[i][j].Zero();
	  // NORTH face:
	  if (Local_SolnBlk[nb].Grid.lfaceN(i,j) > ZERO) {
	    dx = Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    dl = dx.abs();
	    thatN = dx/dl;
	    nhatN = Local_SolnBlk[nb].Grid.nfaceN(i,j);
	    nfaceN = nhatN*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	    Laplacian[i][j].u[1][2] = (dot(nhatN,nfaceN)/dot(nhatN,thatN))/dl;
	    Laplacian[i][j].u[1][1] -= (dot(nhatN,nfaceN)/dot(nhatN,thatN))/dl;
	  }
	  // SOUTH face:
	  if (Local_SolnBlk[nb].Grid.lfaceS(i,j) > ZERO) {
	    dx = Local_SolnBlk[nb].Grid.Cell[i][j-1].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    dl = dx.abs();
	    thatS = dx/dl;
	    nhatS = Local_SolnBlk[nb].Grid.nfaceS(i,j);
	    nfaceS = nhatS*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	    Laplacian[i][j].u[1][0] = (dot(nhatS,nfaceS)/dot(nhatS,thatS))/dl;
	    Laplacian[i][j].u[1][1] -= (dot(nhatS,nfaceS)/dot(nhatS,thatS))/dl;
	  }
	  // EAST face:
	  if (Local_SolnBlk[nb].Grid.lfaceE(i,j) > ZERO) {
	    dx = Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    dl = dx.abs();
	    thatE = dx/dl;
	    nhatE = Local_SolnBlk[nb].Grid.nfaceE(i,j);
	    nfaceE = nhatE*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    Laplacian[i][j].u[2][1] = (dot(nhatE,nfaceE)/dot(nhatE,thatE))/dl;
	    Laplacian[i][j].u[1][1] -= (dot(nhatE,nfaceE)/dot(nhatE,thatE))/dl;
	  }
	  // WEST face:
	  if (Local_SolnBlk[nb].Grid.lfaceW(i,j) > ZERO) {
	    dx = Local_SolnBlk[nb].Grid.Cell[i-1][j].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    dl = dx.abs();
	    thatW = dx/dl;
	    nhatW = Local_SolnBlk[nb].Grid.nfaceW(i,j);
	    nfaceW = nhatW*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	    Laplacian[i][j].u[0][1] = (dot(nhatW,nfaceW)/dot(nhatW,thatW))/dl;
	    Laplacian[i][j].u[1][1] -= (dot(nhatW,nfaceW)/dot(nhatW,thatW))/dl;
	  }
	  // Divide by the cell area to complete stencil construction.
	  for (int m = 0; m < n_laplacian_stencil; m++) {
	    for (int n = 0; n < n_laplacian_stencil; n++) {
	      Laplacian[i][j].u[n][m] /= Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    }
	  }
 	}
      }

      // Output the coefficients of the Laplacian reconstruction to data files.
      out << setprecision(8);
      if (nb == 0) {
	out << " ELLIPTIC OPERATOR ANALYSIS::DIRECTIONAL DERIVATIVE";
	out << endl;
      }
      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
 	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	  out << endl;
	  out << "========================================================================";
	  out << endl << endl << " Xc(" << i << "," << j << ") = ("
	      << Local_SolnBlk[nb].Grid.Cell[i][j].Xc.x << ","
	      << Local_SolnBlk[nb].Grid.Cell[i][j].Xc.y << ")";
	  out << endl << endl << " -> Reconstruction weights: ";
	  for (int m = n_laplacian_stencil-1; m >= 0; m--) {
	    if (m != n_laplacian_stencil-1) out << "                            ";
	    for (int n = 0; n < n_laplacian_stencil; n++) {
	      if (fabs(Laplacian[i][j].u[n][m]/Laplacian[i][j].u[1][1]) < 1.0e-12) Laplacian[i][j].u[n][m] = ZERO;
	      out << " " << -Laplacian[i][j].u[n][m]/Laplacian[i][j].u[1][1];
	    }
	    out << endl;
	  }
	  out << endl << " -> Taylor series coefficients:";
	  Output_Taylor_Series_Coefficients(Local_SolnBlk[nb],i,j,Laplacian[i][j],out);
	  out << endl << endl << " -> Raw coefficients: ";
	  for (int m = n_laplacian_stencil-1; m >= 0; m--) {
	    if (m != n_laplacian_stencil-1) out << "                      ";
	    for (int n = 0; n < n_laplacian_stencil; n++) {
	      out << " " << Laplacian[i][j].u[n][m];
	    }
	    out << endl;
	  }
	}
      }

    }
  }

  // Close the output data file.
  out.close();

  /////////////////////////////////////////////
  // Deallocate reconstruction coefficients. //
  /////////////////////////////////////////////

  for (int i = 0; i < Local_SolnBlk[0].NCi; i++) {
    for (int j = 0; j < Local_SolnBlk[0].NCj; j++) {
      Laplacian[i][j].deallocate();
    }
    delete []Laplacian[i]; Laplacian[i] = NULL;
  }
  delete []Laplacian; Laplacian = NULL;

  // Elliptic operator analysis with directional-derivative face
  // gradient reconstruction performed successfully.
  return 0;

}

/**********************************************************************
 * Routine: Elliptic_Operator_Analysis_Diamond_Path                   *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Elliptic_Operator_Analysis_Diamond_Path(Quad_Soln_Block *Local_SolnBlk,
 					    AdaptiveBlock2D_List &Local_Solution_Block_List,
 					    Quad_Soln_Input_Parameters &IP,
 					    const int &weighting_flag) {

  // Exit immediately if no solution blocks are being used.
  if (!Local_Solution_Block_List.Nused()) return 0;

  ///////////////////////////////////////////
  // Allocate reconstruction coefficients. //
  ///////////////////////////////////////////

  int n_laplacian_stencil = 3;
  EllipticOperatorCoefficients2D **Laplacian;
  double *wa, *wb;

  Laplacian = new EllipticOperatorCoefficients2D*[Local_SolnBlk[0].NCi];
  for (int i = 0; i < Local_SolnBlk[0].NCi; i++) {
    Laplacian[i] = new EllipticOperatorCoefficients2D[Local_SolnBlk[0].NCj];
    for (int j = 0; j < Local_SolnBlk[0].NCj; j++) {
      Laplacian[i][j].allocate(n_laplacian_stencil);
    }
  }
  wa = new double[4];
  wb = new double[4];
  for (int n = 0; n < 4; n++) { wa[n] = QUARTER; wb[n] = QUARTER; }

  ///////////////////////////////////
  // Prepare the output data file. //
  ///////////////////////////////////

  int name;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream out;    

  // Determine prefix of output data file names.
  name = 0;
  while (1) {
    if (IP.Output_File_Name[name] == ' ' || IP.Output_File_Name[name] == '.') break;
    prefix[name] = IP.Output_File_Name[name];
    name++;
    if (name > strlen(IP.Output_File_Name)) break;
  }
  prefix[name] = '\0';
  if (weighting_flag == 0) strcat(prefix,"_eoa_dpsw_cpu");
  else if (weighting_flag == 1) strcat(prefix,"_eoa_dplpwhc_cpu");
  else if (weighting_flag == 2) strcat(prefix,"_eoa_dplpwzy_cpu");
  else return 1;

  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Local_Solution_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data file.
  out.open(output_file_name_ptr,ios::out);
  if (out.fail()) return 1;

  /////////////////////////////////////////////////////////////////
  // Compute the reconstruction coefficients for all of the used //
  // solution blocks and output the data to a file.              //
  /////////////////////////////////////////////////////////////////

  // Declare gradient and Laplacian reconstruction variables:
  double ds, dl;
  Vector2D dx, nhatN, nhatS, nhatE, nhatW, thatN, thatS, thatE, thatW;
  Vector2D nfaceN, nfaceS, nfaceE, nfaceW;

  // Perform the analysis on each solution block that is in use.
  for (int nb = 0; nb < Local_Solution_Block_List.Nblk; nb++) {
    if (Local_Solution_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Compute the coefficients for the Laplacian reconstruction.
      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
 	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	  // Initialize the coefficients:
 	  Laplacian[i][j].Zero();
	  // NORTH face:
	  if (Local_SolnBlk[nb].Grid.lfaceN(i,j) > ZERO) {
	    if (weighting_flag == 1) {
	      Linear_Weighting_Coefficients_HC(Local_SolnBlk[nb],i-1,j,wa);
	      Linear_Weighting_Coefficients_HC(Local_SolnBlk[nb],i,j,wb);
	    } else if (weighting_flag == 2) {
	      Linear_Weighting_Coefficients_ZY(Local_SolnBlk[nb],i-1,j,wa);
	      Linear_Weighting_Coefficients_ZY(Local_SolnBlk[nb],i,j,wb);
	    }
	    dx = Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    ds = dx.abs();
	    thatN = dx/ds;
	    nhatN = Vector2D(-thatN.y,thatN.x);
	    nfaceN = Local_SolnBlk[nb].Grid.nfaceN(i,j);
	    dl = Local_SolnBlk[nb].Grid.lfaceN(i,j);
	    Laplacian[i][j].u[0][1] -= (dot(nhatN,nfaceN)/dot(nfaceN,thatN))*(wb[0]*dl/dl);
	    Laplacian[i][j].u[1][1] -= (dot(nfaceN/ds + nhatN*wb[2]/dl - nhatN*wa[3]/dl,nfaceN*dl)/dot(nfaceN,thatN));
	    Laplacian[i][j].u[2][1] += (dot(nhatN,nfaceN)/dot(nfaceN,thatN))*(wa[1]*dl/dl);
	    Laplacian[i][j].u[0][2] -= (dot(nhatN,nfaceN)/dot(nfaceN,thatN))*(wb[3]*dl/dl);
	    Laplacian[i][j].u[1][2] -= (dot(- nfaceN/ds + nhatN*wb[1]/dl - nhatN*wa[0]/dl,nfaceN*dl)/dot(nfaceN,thatN));
	    Laplacian[i][j].u[2][2] += (dot(nhatN,nfaceN)/dot(nfaceN,thatN))*(wa[2]*dl/dl);
	  }
	  // SOUTH face:
	  if (Local_SolnBlk[nb].Grid.lfaceS(i,j) > ZERO) {
	    if (weighting_flag == 1) {
	      Linear_Weighting_Coefficients_HC(Local_SolnBlk[nb],i-1,j-1,wa);
	      Linear_Weighting_Coefficients_HC(Local_SolnBlk[nb],i,j-1,wb);
	    } else if (weighting_flag == 2) {
	      Linear_Weighting_Coefficients_ZY(Local_SolnBlk[nb],i-1,j-1,wa);
	      Linear_Weighting_Coefficients_ZY(Local_SolnBlk[nb],i,j-1,wb);
	    }
	    dx = Local_SolnBlk[nb].Grid.Cell[i][j-1].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    ds = dx.abs();
	    thatS = dx/ds;
	    nhatS = Vector2D(-thatS.y,thatS.x);
	    nfaceS = Local_SolnBlk[nb].Grid.nfaceS(i,j);
	    dl = Local_SolnBlk[nb].Grid.lfaceS(i,j);
	    Laplacian[i][j].u[0][0] += (dot(nhatS,nfaceS)/dot(nfaceS,thatS))*(wa[0]*dl/dl);
	    Laplacian[i][j].u[1][0] -= (dot(- nfaceS/ds - nhatS*wa[1]/dl + nhatS*wb[0]/dl,nfaceS*dl)/dot(nfaceS,thatS));
	    Laplacian[i][j].u[2][0] -= (dot(nhatS,nfaceS)/dot(nfaceS,thatS))*(wb[1]*dl/dl);
	    Laplacian[i][j].u[0][1] += (dot(nhatS,nfaceS)/dot(nfaceS,thatS))*(wa[3]*dl/dl);
	    Laplacian[i][j].u[1][1] -= (dot(nfaceS/ds - nhatS*wa[2]/dl + nhatS*wa[3]/dl,nfaceS*dl)/dot(nfaceS,thatS));
	    Laplacian[i][j].u[2][1] -= (dot(nhatS,nfaceS)/dot(nfaceS,thatS))*(wb[2]*dl/dl);
	  }
 	  // EAST face:
	  if (Local_SolnBlk[nb].Grid.lfaceE(i,j) > ZERO) {
	    if (weighting_flag == 1) {
	      Linear_Weighting_Coefficients_HC(Local_SolnBlk[nb],i,j-1,wa);
	      Linear_Weighting_Coefficients_HC(Local_SolnBlk[nb],i,j,wb);
	    } else if (weighting_flag == 2) {
	      Linear_Weighting_Coefficients_ZY(Local_SolnBlk[nb],i,j-1,wa);
	      Linear_Weighting_Coefficients_ZY(Local_SolnBlk[nb],i,j,wb);
	    }
	    dx = Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    ds = dx.abs();
	    thatE = dx/ds;
	    nhatE = Vector2D(-thatE.y,thatE.x);
	    nfaceE = Local_SolnBlk[nb].Grid.nfaceE(i,j);
	    dl = Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    Laplacian[i][j].u[1][2] -= (dot(nhatE,nfaceE)/dot(nfaceE,thatE))*(wb[3]*dl/dl);
	    Laplacian[i][j].u[1][1] -= (dot(nfaceE/ds + nhatE*wb[0]/dl - nhatE*wa[3]/dl,nfaceE*dl)/dot(nfaceE,thatE));
	    Laplacian[i][j].u[1][0] += (dot(nhatE,nfaceE)/dot(nfaceE,thatE))*(wa[0]*dl/dl);
	    Laplacian[i][j].u[2][2] -= (dot(nhatE,nfaceE)/dot(nfaceE,thatE))*(wb[2]*dl/dl);
	    Laplacian[i][j].u[2][1] -= (dot(- nfaceE/ds + nhatE*wb[1]/dl - nhatE*wa[2]/dl,nfaceE*dl)/dot(nfaceE,thatE));
	    Laplacian[i][j].u[2][0] += (dot(nhatE,nfaceE)/dot(nfaceE,thatE))*(wa[1]*dl/dl);
	  }
 	  // WEST face:
	  if (Local_SolnBlk[nb].Grid.lfaceW(i,j) > ZERO) {
	    if (weighting_flag == 1) {
	      Linear_Weighting_Coefficients_HC(Local_SolnBlk[nb],i-1,j,wa);
	      Linear_Weighting_Coefficients_HC(Local_SolnBlk[nb],i-1,j-1,wb);
	    } else if (weighting_flag == 2) {
	      Linear_Weighting_Coefficients_ZY(Local_SolnBlk[nb],i-1,j,wa);
	      Linear_Weighting_Coefficients_ZY(Local_SolnBlk[nb],i-1,j-1,wb);
	    }
	    dx = Local_SolnBlk[nb].Grid.Cell[i-1][j].Xc - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    ds = dx.abs();
	    thatW = dx/ds;
	    nhatW = Vector2D(-thatW.y,thatW.x);
	    nfaceW = Local_SolnBlk[nb].Grid.nfaceW(i,j);
	    dl = Local_SolnBlk[nb].Grid.lfaceW(i,j);
	    Laplacian[i][j].u[1][0] -= (dot(nhatW,nfaceW)/dot(nfaceW,thatW))*(wb[1]*dl/dl);
	    Laplacian[i][j].u[1][1] -= (dot(nfaceW/ds + nhatW*wb[2]/dl - nhatW*wa[1]/dl,nfaceW*dl)/dot(nfaceW,thatW));
	    Laplacian[i][j].u[1][2] += (dot(nhatW,nfaceW)/dot(nfaceW,thatW))*(wa[2]*dl/dl);
	    Laplacian[i][j].u[0][0] -= (dot(nhatW,nfaceW)/dot(nfaceW,thatW))*(wb[0]*dl/dl);
	    Laplacian[i][j].u[0][1] -= (dot(- nfaceW/ds + nhatW*wb[3]/dl - nhatW*wa[0]/dl,nfaceW*dl)/dot(nfaceW,thatW));
	    Laplacian[i][j].u[0][2] += (dot(nhatW,nfaceW)/dot(nfaceW,thatW))*(wa[3]*dl/dl);
	  }
	  // Divide by the cell area to complete stencil construction.
	  for (int m = 0; m < n_laplacian_stencil; m++) {
	    for (int n = 0; n < n_laplacian_stencil; n++) {
	      Laplacian[i][j].u[n][m] /= Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    }
	  }
 	}
      }

      // Output the coefficients of the Laplacian reconstruction to data files.
      out << setprecision(8);
      if (nb == 0) {
	out << " ELLIPTIC OPERATOR ANALYSIS::DIAMOND PATH WITH";
	if (weighting_flag == 0) out << " SIMPLE WEIGHTING";
	else if (weighting_flag == 1) out << " LINEARITY PRESERVING WEIGHTING -- HOLMES AND CONNELL";
	else if (weighting_flag == 2) out << " LINEARITY PRESERVING WEIGHTING -- ZINGG AND YARROW";
	out << endl;
      }
      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
 	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	  out << endl;
	  out << "========================================================================";
	  out << endl << endl << " Xc(" << i << "," << j << ") = ("
	      << Local_SolnBlk[nb].Grid.Cell[i][j].Xc.x << ","
	      << Local_SolnBlk[nb].Grid.Cell[i][j].Xc.y << ")";
	  out << endl << endl << " -> Reconstruction weights: ";
	  for (int m = n_laplacian_stencil-1; m >= 0; m--) {
	    if (m != n_laplacian_stencil-1) out << "                            ";
	    for (int n = 0; n < n_laplacian_stencil; n++) {
	      if (fabs(Laplacian[i][j].u[n][m]/Laplacian[i][j].u[1][1]) < 1.0e-12) Laplacian[i][j].u[n][m] = ZERO;
	      out << " " << -Laplacian[i][j].u[n][m]/Laplacian[i][j].u[1][1];
	    }
	    out << endl;
	  }
	  out << endl << " -> Taylor series coefficients:";
	  Output_Taylor_Series_Coefficients(Local_SolnBlk[nb],i,j,Laplacian[i][j],out);
	  out << endl << endl << " -> Raw coefficients: ";
	  for (int m = n_laplacian_stencil-1; m >= 0; m--) {
	    if (m != n_laplacian_stencil-1) out << "                      ";
	    for (int n = 0; n < n_laplacian_stencil; n++) {
	      out << " " << Laplacian[i][j].u[n][m];
	    }
	    out << endl;
	  }
	}
      }

    }
  }

  // Close the output data file.
  out.close();

  /////////////////////////////////////////////
  // Deallocate reconstruction coefficients. //
  /////////////////////////////////////////////

  for (int i = 0; i < Local_SolnBlk[0].NCi; i++) {
    for (int j = 0; j < Local_SolnBlk[0].NCj; j++) {
      Laplacian[i][j].deallocate();
    }
    delete []Laplacian[i]; Laplacian[i] = NULL;
  }
  delete []Laplacian; Laplacian = NULL;

  delete []wa; wa = NULL;
  delete []wb; wb = NULL;

  // Elliptic operator analysis with diamond-path face gradient
  // reconstruction performed successfully.
  return 0;

}

/**********************************************************************
 * Routine: Output_Taylor_Series_Coefficients                         *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block>
void Output_Taylor_Series_Coefficients(Quad_Soln_Block SolnBlk,
				       const int &i, const int &j,
				       EllipticOperatorCoefficients2D &Laplacian,
				       ofstream &out) {

  int nc = Laplacian.nc;
  double xhat, yhat,
         alpha, alphadx, alphady, alphadx2, alphadxdy, alphady2,
         alphadx3, alphadx2dy, alphadxdy2, alphady3,
         alphadx4, alphadx3dy, alphadx2dy2, alphadxdy3, alphady4;

  alpha = ZERO; alphadx = ZERO; alphady = ZERO;
  alphadx2 = ZERO; alphadxdy = ZERO; alphady2 = ZERO;
  alphadx3 = ZERO; alphadx2dy = ZERO; alphadxdy2 = ZERO; alphady3 = ZERO;
  alphadx4 = ZERO; alphadx3dy = ZERO; alphadx2dy2 = ZERO; alphadxdy3 = ZERO; alphady4 = ZERO;

  for (int jj = j-nc; jj <= j+nc; jj++) {
    for (int ii = i-nc; ii <= i+nc; ii++) {
      xhat = SolnBlk.Grid.Cell[ii][jj].Xc.x - SolnBlk.Grid.Cell[i][j].Xc.x;
      yhat = SolnBlk.Grid.Cell[ii][jj].Xc.y - SolnBlk.Grid.Cell[i][j].Xc.y;
      alpha += Laplacian.u[ii-i+nc][jj-j+nc];
      alphadx += Laplacian.u[ii-i+nc][jj-j+nc]*xhat;
      alphady += Laplacian.u[ii-i+nc][jj-j+nc]*yhat;
      alphadx2 += Laplacian.u[ii-i+nc][jj-j+nc]*sqr(xhat);
      alphadxdy += Laplacian.u[ii-i+nc][jj-j+nc]*xhat*yhat;
      alphady2 += Laplacian.u[ii-i+nc][jj-j+nc]*sqr(yhat);
      alphadx3 += Laplacian.u[ii-i+nc][jj-j+nc]*cube(xhat);
      alphadx2dy += Laplacian.u[ii-i+nc][jj-j+nc]*sqr(xhat)*yhat;
      alphadxdy2 += Laplacian.u[ii-i+nc][jj-j+nc]*xhat*sqr(yhat);
      alphady3 += Laplacian.u[ii-i+nc][jj-j+nc]*cube(yhat);
      alphadx4 += Laplacian.u[ii-i+nc][jj-j+nc]*sqr(xhat)*sqr(xhat);
      alphadx3dy += Laplacian.u[ii-i+nc][jj-j+nc]*cube(xhat)*yhat;
      alphadx2dy2 += Laplacian.u[ii-i+nc][jj-j+nc]*sqr(xhat)*sqr(yhat);
      alphadxdy3 += Laplacian.u[ii-i+nc][jj-j+nc]*xhat*cube(yhat);
      alphady4 += Laplacian.u[ii-i+nc][jj-j+nc]*sqr(yhat)*sqr(yhat);
    }
  }
  if (fabs(alpha) < NANO) alpha = ZERO;
  if (fabs(alphadx) < NANO) alphadx = ZERO;
  if (fabs(alphady) < NANO) alphady = ZERO;
  if (fabs(alphadx2) < NANO) alphadx2 = ZERO;
  if (fabs(alphadxdy) < NANO) alphadxdy = ZERO;
  if (fabs(alphady2) < NANO) alphady2 = ZERO;
  if (fabs(alphadx3) < NANO) alphadx3 = ZERO;
  if (fabs(alphadx2dy) < NANO) alphadx2dy = ZERO;
  if (fabs(alphadxdy2) < NANO) alphadxdy2 = ZERO;
  if (fabs(alphady3) < NANO) alphady3 = ZERO;
  if (fabs(alphadx4) < NANO) alphadx4 = ZERO;
  if (fabs(alphadx3dy) < NANO) alphadx3dy = ZERO;
  if (fabs(alphadx2dy2) < NANO) alphadx2dy2 = ZERO;
  if (fabs(alphadxdy3) < NANO) alphadxdy3 = ZERO;
  if (fabs(alphady4) < NANO) alphady4 = ZERO;
  out << endl << "    L(u) = (" << alpha << ") + " << "(" << alphadx << ") dudx + (" << alphady << ") dudy + "
      << endl << "           (" << alphadx2 << ")/2 d2udx2 + (" << alphadxdy << ") d2udxdy + (" << alphady2 << ")/2 d2udy2 + "
      << endl << "           (" << alphadx3 << ")/6 d3udx3 + (" << alphadx2dy << ")/2 d3udx2dy + (" << alphadxdy2 << ")/2 d3udxdy2 + (" << alphady3 << ")/6 d3udy3 + ...";
//       << endl << "           (" << alphadx4 << ")/24 d4udx4 + (" << alphadx3dy << ")/6 d4udx3dy + (" << alphadx2dy2 << ")/4 d4udx2dy2 + (" << alphadxdy3 << ")/6 d4udxdy3 + (" << alphady4 << ")/24 d4udy4 + ... ";

}

/**********************************************************************
 * Routine: Linear_Weighting_Coefficients_HC                          *
 *                                                                    *
 * Holmes and Connell (AIAA Paper 1989-1932-CP)                       *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block>
void Linear_Weighting_Coefficients_HC(Quad_Soln_Block SolnBlk,
				      const int &ii, const int &jj,
				      double *w) {
  Vector2D lambda, R, X0, X1, X2, X3, X4;
  Tensor2D I;
  // Summarize cell-centres and states.
  X0 = SolnBlk.Grid.Node[ii][jj].X;
  X1 = SolnBlk.Grid.Cell[ii  ][jj  ].Xc;
  X2 = SolnBlk.Grid.Cell[ii+1][jj  ].Xc;
  X3 = SolnBlk.Grid.Cell[ii+1][jj+1].Xc;
  X4 = SolnBlk.Grid.Cell[ii  ][jj+1].Xc;
  // Determine weighting coefficients:
  R = X1 + X2 + X3 + X4 - FOUR*X0;
  I.xx = sqr(X1.x - X0.x) + sqr(X2.x - X0.x) + sqr(X3.x - X0.x) + sqr(X4.x - X0.x);
  I.xy = (X1.x - X0.x)*(X1.y - X0.y) + (X2.x - X0.x)*(X2.y - X0.y) +
         (X3.x - X0.x)*(X3.y - X0.y) + (X4.x - X0.x)*(X4.y - X0.y);
  I.yy = sqr(X1.y - X0.y) + sqr(X2.y - X0.y) + sqr(X3.y - X0.y) + sqr(X4.y - X0.y);
  lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
  lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
  // Determine the weights:
  w[0] = ONE + lambda*(X1 - X0);
  w[1] = ONE + lambda*(X2 - X0);
  w[2] = ONE + lambda*(X3 - X0);
  w[3] = ONE + lambda*(X4 - X0);
}

/**********************************************************************
 * Routine: Linear_Weighting_Coefficients_ZY                          *
 *                                                                    *
 * Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13 No. 3 1992)   *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block>
void Linear_Weighting_Coefficients_ZY(Quad_Soln_Block SolnBlk,
				      const int &ii, const int &jj,
				      double *w) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
         eta1, zeta1, eta2, zeta2, eta, zeta;
  x  = SolnBlk.Grid.Node[ii][jj].X.x;
  y  = SolnBlk.Grid.Node[ii][jj].X.y;
  ax = SolnBlk.Grid.Cell[ii-1][jj-1].Xc.x;
  bx = SolnBlk.Grid.Cell[ii-1][jj  ].Xc.x - SolnBlk.Grid.Cell[ii-1][jj-1].Xc.x; 
  cx = SolnBlk.Grid.Cell[ii  ][jj-1].Xc.x - SolnBlk.Grid.Cell[ii-1][jj-1].Xc.x; 
  dx = SolnBlk.Grid.Cell[ii  ][jj  ].Xc.x + SolnBlk.Grid.Cell[ii-1][jj-1].Xc.x -
       SolnBlk.Grid.Cell[ii-1][jj  ].Xc.x - SolnBlk.Grid.Cell[ii  ][jj-1].Xc.x;
  ay = SolnBlk.Grid.Cell[ii-1][jj-1].Xc.y;
  by = SolnBlk.Grid.Cell[ii-1][jj  ].Xc.y - SolnBlk.Grid.Cell[ii-1][jj-1].Xc.y; 
  cy = SolnBlk.Grid.Cell[ii  ][jj-1].Xc.y - SolnBlk.Grid.Cell[ii-1][jj-1].Xc.y; 
  dy = SolnBlk.Grid.Cell[ii  ][jj  ].Xc.y + SolnBlk.Grid.Cell[ii-1][jj-1].Xc.y -
       SolnBlk.Grid.Cell[ii-1][jj  ].Xc.y - SolnBlk.Grid.Cell[ii  ][jj-1].Xc.y;
  aa = bx*dy - dx*by;
  bb = dy*(ax-x) + bx*cy - cx*by+dx*(y-ay);
  cc = cy*(ax-x) + cx*(y-ay);
  if (fabs(aa) < TOLER*TOLER) {
    if (fabs(bb) >= TOLER*TOLER) zeta1 = -cc/bb;
    else zeta1 = -cc/sgn(bb)*(TOLER*TOLER);
    if (fabs(cy+dy*zeta1) >= TOLER*TOLER) eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1); 
    else eta1 = HALF;
    zeta2 = zeta1;
    eta2  = eta1;
  } else {
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) zeta1 = HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; 
    else zeta1 = -HALF*bb/aa;
    if (fabs(cy+dy*zeta1) < TOLER*TOLER) eta1 = -ONE;
    else eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1);
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) zeta2 = HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa; 
    else zeta2 = -HALF*bb/aa;
    if (fabs(cy+dy*zeta2) < TOLER*TOLER) eta2 = -ONE;
    else eta2 = (y-ay-by*zeta2)/(cy+dy*zeta2);
  }
  if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
      eta1  > -TOLER && eta1  < ONE + TOLER) {
    zeta = zeta1;
    eta  = eta1;
  } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
	     eta2  > -TOLER && eta2  < ONE + TOLER) {
    zeta = zeta2;
    eta  = eta2;
  } else {
    zeta = HALF;
    eta  = HALF;
  }
  // Determine the weights:
  w[0] = ONE - zeta - eta - zeta*eta;
  w[1] = eta - zeta*eta;
  w[2] = zeta*eta;
  w[3] = zeta - zeta*eta;
}

#endif // _ELLIPTIC2D_INCLUDED
