/* CompDomain1D.h: Header file defining 1D domanin types*/

#ifndef _COMPDOMAIN1D_INCLUDED
#define _COMPDOMAIN1D_INCLUDED

/* Include math macro header file. */

#ifndef _RECONSTRUCT1D_INPUT_INCLUDED
#include "Reconstruct1DInput.h"
#endif // _RECONSTRUCT1D_INPUT_INCLUDED

#ifndef _COMPCELL1D_INCLUDED
#include "CompCell1D.h"
#endif // _COMPCELL1D_INCLUDED

#ifndef _TESTFUNCTIONS_INCLUDED
#include "TestFunctions/TestFunctions.h"
#endif // _TESTFUNCTIONS_INCLUDED

#include <cstring>
using namespace std;

/* Define the classes */

/**********************************************************
 * Class: CompDomain1D_Uniform                            *
 *                                                        *
 * Member functions                                       *
 *     *Mesh     -- Pointer of type CompCell1D_Uniform    *
 *                  points to the vector of the cells     *
 *     NCi -- Number of cells in I direction              *
 *     NCS -- Half the number of cells of the stencil     *
 *            minus the central value                     *
 *       ex. NCS = 2 => the stencil is: j-2,j-1,j,j+1,j+2 *
 *********************************************************/

class CompDomain1D_Uniform{
 private:
 public:
  CompCell1D_Uniform* Mesh;
  int NCi;
  int NCS;

  /*Creation, copy, and assignment constructors. */
  CompDomain1D_Uniform(){
    NCi = 0;
    NCS = 0;
    Mesh = NULL;
  };
  CompDomain1D_Uniform(const CompDomain1D_Uniform &G){
    CompCell1D_Uniform temp;
    NCi = G.NCi;
    NCS = G.NCS;
    Mesh = new CompCell1D_Uniform[G.NCi];
    for (int i=0; i <= NCi-1; i++){
      temp = G.Mesh[i];
      Mesh[i] = temp;
    }
  };

  /*Destructor */
  //Use the default destructor
  //~CompDomain1D_Uniform();
  
  /* Allocate memory for Computational Domain. */
  void allocate(const int Ni);
  
  /* Deallocate memory for Computational Domain. */
  void deallocate(void);
};

/**************************************************************************
 * CompDomain1D_Uniform::allocate -- Allocate memory.                     *
 **************************************************************************/

inline void CompDomain1D_Uniform::allocate(const int Ni){

  NCi = Ni;
  Mesh = new CompCell1D_Uniform[NCi];
}

/**************************************************************************
 * CompDomain1D_Uniform::deallocate -- Deallocate memory.                 *
 **************************************************************************/

inline void CompDomain1D_Uniform::deallocate(void){

  delete [] Mesh; Mesh = NULL;
}

/********************************************************
 * CompDomain1D_Uniform -- External subrutines          *
 *******************************************************/

extern void Generate_Mesh(CompDomain1D_Uniform &Grid,
			  const Reconstruct1D_Input_Parameters &IP);

extern void SetGeometry(CompDomain1D_Uniform &Grid,
			const Reconstruct1D_Input_Parameters &IP);

extern void Compute_Cell_Average(CompDomain1D_Uniform &Grid,
				 const Reconstruct1D_Input_Parameters &IP);

/**********************************************************
 * Class: CompDomain1D_NonUniform                         *
 *                                                        *
 * Member functions                                       *
 *     *Mesh     -- Pointer of type CompCell1D_NonUniform *
 *                  points to the vector of the cells     *
 *     NCi -- Number of cells in I direction              *
 *     NCS -- Half the number of cells of the stencil     *
 *            minus the central value                     *
 *       ex. NCS = 2 => the stencil is: j-2,j-1,j,j+1,j+2 *
 *********************************************************/

class CompDomain1D_NonUniform{
 private:
 public:
  CompCell1D_NonUniform* Mesh;
  int NCi;
  int NCS;
  double MinSolution;
  double MaxSolution;
  double MaxDeltaSolutionOverDomain;

  /*Creation, copy, and assignment constructors. */

  CompDomain1D_NonUniform():
    NCi(0), NCS(0), Mesh(NULL){}

  // Copy constructor 1
  CompDomain1D_NonUniform(const CompDomain1D_NonUniform &);

  //Use the default destructor
  //~CompDomain1D_NonUniform();

  /* Allocate memory for Computational Domain. */
  void allocate(const int Ni);
  
  /* Deallocate memory for Computational Domain. */
  void deallocate(void);

  // Compute Value in Point
  double ComputeValueInPoint(double x);
};

/**************************************************************************
 * Compdomain1d_Nonuniform::allocate -- Allocate memory.                  *
 **************************************************************************/

inline void CompDomain1D_NonUniform::allocate(const int Ni){

  NCi = Ni;
  Mesh = new CompCell1D_NonUniform[NCi];
  if (Mesh == NULL){
    cout << "Not enough memory!\n";
    exit(1);
  }
}

/**************************************************************************
 * CompDomain1D_NonUniform::deallocate -- Deallocate memory.              *
 **************************************************************************/

inline void CompDomain1D_NonUniform::deallocate(void){

  for (int i=0; i<=NCi-1; i++)
    Mesh[i].deallocate();
  cout << "\n Deallocate Computational_Cells";
  delete [] Mesh; Mesh = NULL;
  cout << "\n Deallocate Computational_Domain1D -- Mesh";
}

/***********************************************************
 * CompDomain1D_NonUniform -- External subrutines          *
 **********************************************************/

extern CompDomain1D_NonUniform* Generate_Mesh(CompDomain1D_NonUniform *Grid,
					      const Reconstruct1D_Input_Parameters &IP);

extern int SetGeometry(CompDomain1D_NonUniform *Grid,
			const Reconstruct1D_Input_Parameters &IP);

extern void Compute_Cell_Average(CompDomain1D_NonUniform *Grid,
				 const Reconstruct1D_Input_Parameters &IP);

extern int Output_Tecplot(const CompDomain1D_NonUniform *Grid,
			  const Reconstruct1D_Input_Parameters &IP);

extern int Output_Tecplot_Domain(const CompDomain1D_NonUniform *Grid,
				 const Reconstruct1D_Input_Parameters &IP);

extern int Output_Matlab (const CompDomain1D_NonUniform *Grid,
			  const Reconstruct1D_Input_Parameters &IP);

extern int Output_Derivatives_Tecplot(const CompDomain1D_NonUniform *Grid,
				      const Reconstruct1D_Input_Parameters &IP);

extern int Output_Derivatives_Matlab(const CompDomain1D_NonUniform *Grid,
				     const Reconstruct1D_Input_Parameters &IP);

extern int Output_Cell_Center_Values(const CompDomain1D_NonUniform *Grid,
				     const Reconstruct1D_Input_Parameters &IP);

extern int Output_Function_Graph_Matlab(const Reconstruct1D_Input_Parameters &IP);

extern int Output_Function_Graph_Tecplot(const Reconstruct1D_Input_Parameters &IP);

extern void  Cell_Theoretic_Derives(const CompCell1D_NonUniform &Cell,
				    const Reconstruct1D_Input_Parameters &IP,
				    int Order, double* D);

/* Compute the L1 norm*/
extern double L1_Norm(const CompDomain1D_NonUniform *Grid,
		      const Reconstruct1D_Input_Parameters &IP);

/* Compute the L2 norm*/
extern double L2_Norm(const CompDomain1D_NonUniform *Grid,
		      const Reconstruct1D_Input_Parameters &IP);

extern int Output_L1_Norm_Tecplot(const CompDomain1D_NonUniform *Grid,
				  const Reconstruct1D_Input_Parameters &IP);

extern int Output_L2_Norm_Tecplot(const CompDomain1D_NonUniform *Grid,
				  const Reconstruct1D_Input_Parameters &IP);

#endif // _COMPDOMAIN1D_INCLUDED
