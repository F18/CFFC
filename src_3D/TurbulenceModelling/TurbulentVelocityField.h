#ifndef _TURBULENT_VELOCITY_FIELD_INCLUDED 
#define _TURBULENT_VELOCITY_FIELD_INCLUDED

/* Include required C++ header files. */

#include <cmath> 
#include <cassert>
#include <cstdlib>     // defines the drand48() function
#include <ctime>       // defines the time() function
#include <limits>
#include <complex>

/* Include required CFFC header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _NumericalLibrary_INCLUDED
#include "../Math/NumericalLibrary.h"
#endif // _NumericalLibrary_INCLUDED

#ifndef _SPLINEFIT_INCLUDED
#include "../Math/SplineFit.h"
#endif // _SPLINEFIT_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif // _VECTOR3D_INCLUDED

#ifndef _TENSOR3D_INCLUDED
#include "../Math/Tensor3D.h"
#endif // _TENSOR3D_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // _INPUT_INCLUDED

#ifndef _OCTREE_INCLUDED
#include "../AMR/Octree.h"
#endif // _OCTREE_INCLUDED

#ifndef _GRID3D_HEXA_MULTIBLOCK_INCLUDED
#include "../Grid/Grid3DHexaMultiBlock.h"
#endif // _GRID3D_HEXA_MULTIBLOCK_INCLUDED

#ifndef _ADAPTIVEBLOCK3D_INCLUDED
#include "../AMR/AdaptiveBlock3D.h"
#endif // _ADAPTIVEBLOCK3D_INCLUDED

#ifndef _FFTW_INCLUDED
#include "fftw3.h"
#endif // _FFTW_INCLUDED

#ifndef _GNUPLOT_INCLUDED
#include "../System/gnuplot.h"
#endif // _GNUPLOT_INCLUDED

#ifndef _BINARY_TREE_INCLUDED
#include "../Math/Binarytree.h"
#endif // _BINARY_TREE_INCLUDED


// Constants
const complex<double>  I(0.0, 1.0);      // sqrt(-1.0)

#define TURBULENT_VELOCITY_FIELD_DATA_USED        1

#define TURBULENT_VELOCITY_FIELD_DATA_NOT_USED    0




/*!
 * Class: Turbulent_Velocity_Field_Block
 *
 * \brief Class used to store turbulent velocity field for a single block.
 *
 */
class Turbulent_Velocity_Field_Block {
  public:
  int                   gblknum; // global turbulent velocity field number
    int             NCi,ICl,ICu; // i-direction turbulent velocity field cell counters
    int             NCj,JCl,JCu; // j-direction turbulent velocity field cell counters
    int             NCk,KCl,KCu; // k-direction turbulent velocity field cell counters
    int                  Nghost; // number of ghost cells
    Vector3D        ***Velocity; // array of turbulent velocity field vectors
    Vector3D        ***Position; // array of turbulent velocity field position vectors
    Vector3D   Node_INl_JNl_KNl, // diagonally opposite corners of the block 
               Node_INu_JNu_KNu;


    int Allocated; // Indicates whether or not the turbulent velocity field data has been allocated.
    
    /* Constructors. */
    Turbulent_Velocity_Field_Block(void) {
       NCi = 0; ICl = 0; ICu = 0; 
       NCj = 0; JCl = 0; JCu = 0;
       NCk = 0; KCl = 0; KCu = 0;
       Nghost = 0;
       gblknum = 0;
       Allocated = TURBULENT_VELOCITY_FIELD_DATA_NOT_USED;
       Velocity = NULL;
       Position = NULL;
    }

    Turbulent_Velocity_Field_Block(const int Ni, 
                                   const int Nj, 
                                   const int Nk,
                                   const int Ng) {
      allocate(Ni, Nj, Nk, Ng);
    }

    /* Destructor. */
    ~Turbulent_Velocity_Field_Block(void) {
       deallocate();
    }

    /* Allocate memory for velocity field data. */
    void allocate(const int Ni, 
                  const int Nj, 
                  const int Nk,
                  const int Ng);
    
    /* Deallocate memory for velocity field data. */
    void deallocate(void);

    /* Reconstruct velocity field data. */
    void LeastSquares_Reconstruction(const int &i,
				     const int &j,
				     const int &k,
				     Vector3D &dVdx,
				     Vector3D &dVdy,
				     Vector3D &dVdz);

    void Copy(Turbulent_Velocity_Field_Block &Block2);

    void Broadcast(void);

#ifdef _MPI_VERSION    
    void Broadcast(MPI::Intracomm &Communicator, 
                   const int Source_CPU);
#endif
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, 
                                 const Turbulent_Velocity_Field_Block &V);

    friend istream &operator >> (istream &in_file, 
                                 Turbulent_Velocity_Field_Block &V);
  
    /* Other useful member functions. */

  private:
    //copy and assignment are not permitted
    Turbulent_Velocity_Field_Block(const Turbulent_Velocity_Field_Block &V);
    Turbulent_Velocity_Field_Block &operator =(const Turbulent_Velocity_Field_Block &V);
};

/*************************************************************************
 * Turbulent_Velocity_Field_Block::allocate -- Allocate memory.          *
 *************************************************************************/
inline void Turbulent_Velocity_Field_Block::allocate(const int Ni, 
                                                     const int Nj, 
                                                     const int Nk,
                                                     const int Ng) {
   assert( Ni >= 1 && Nj >= 1 && Nk >= 1 && Ng >=1 && !Allocated);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
   NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1;
   NCk = Nk+2*Ng; KCl = Ng; KCu = Nk+Ng-1;
   Nghost = Ng;
   Allocated = TURBULENT_VELOCITY_FIELD_DATA_USED;

   Velocity = new Vector3D**[NCi];
   Position = new Vector3D**[NCi];
   for (int i = 0; i <= NCi-1; ++i ){
      Velocity[i] = new Vector3D*[NCj];
      Position[i] = new Vector3D*[NCj];
      for (int j = 0; j <= NCj-1; ++j ){
         Velocity[i][j] = new Vector3D[NCk];
	 Position[i][j] = new Vector3D[NCk];
      } /* endfor */
   } /* endfor */
}

/*************************************************************************
 * Turbulent_Velocity_Field_Block::deallocate -- Deallocate memory.      *
 *************************************************************************/
inline void Turbulent_Velocity_Field_Block::deallocate(void) {
   if (Allocated) {
      assert(NCi >= 1 && NCj >= 1 && NCk >= 1);
      for (int i = 0; i <= NCi-1 ; ++i ) {
         for ( int j = 0 ; j <= NCj-1 ; ++j) {
            delete []Velocity[i][j]; Velocity[i][j] = NULL;
	    delete []Position[i][j]; Position[i][j] = NULL;
         } /* endfor */
         delete []Velocity[i]; Velocity[i] = NULL;
	 delete []Position[i]; Position[i] = NULL;
      }/*endfor*/
      delete []Velocity; Velocity = NULL;
      delete []Position; Position = NULL;
  
      NCi = 0; ICl = 0; ICu = 0; 
      NCj = 0; JCl = 0; JCu = 0; 
      NCk = 0; KCl = 0; KCu = 0; 
      Nghost = 0;
      Allocated = TURBULENT_VELOCITY_FIELD_DATA_NOT_USED;
   } /* endif */
}

/*************************************************************************
 * Turbulent_Velocity_Field_Block::Copy -- Copy a block.                 *
 *************************************************************************/
inline void Turbulent_Velocity_Field_Block::Copy(Turbulent_Velocity_Field_Block &Block2){
  if (Block2.Allocated) {
    //  Allocate memory as required.
    if (NCi != Block2.NCi ||
	NCj != Block2.NCj ||
	NCk != Block2.NCk ||
	Nghost != Block2.Nghost) {
      if (Allocated) {
	deallocate();
      } /*endif */
            
      allocate(Block2.NCi-2*Block2.Nghost, 
	       Block2.NCj-2*Block2.Nghost, 
	       Block2.NCk-2*Block2.Nghost, 
	       Block2.Nghost);
    } /* endif */
        
    /* Assign initial turbulent velocity field to solution block. */
        
    for (int i = ICl ; i <= ICu ; i++) {
      for (int j = JCl ; j <= JCu ; j++) {
	for (int k = KCl ; k <= KCu ; k++) {
	  Velocity[i][j][k] = Block2.Velocity[i][j][k];
	  Position[i][j][k] = Block2.Position[i][j][k];
	} /* endfor */
      } /* endfor */
    } /* endfor */

    gblknum = Block2.gblknum;
    Node_INl_JNl_KNl = Block2.Node_INl_JNl_KNl;
    Node_INu_JNu_KNu = Block2.Node_INu_JNu_KNu;
        
  }
        
}
   



/*!
 * Class: Turbulent_Velocity_Field_Multi_Block
 *
 * \brief Class used to store turbulent velocity field for a
 *        1D array of turbulent velocity field solution blocks.
 *
 */
class Turbulent_Velocity_Field_Multi_Block_List {
  public:
    Turbulent_Velocity_Field_Block  *Vel_Blks; // one dimensional array of velocity block.
    int                                  NBlk;
    int                             NBlk_Idir, 
                                    NBlk_Jdir, 
                                    NBlk_Kdir; // Number of blocks in i, j and k directions.
    int                             Allocated; // Indicates if the velocity blocks have been allocated or not.
 
    /* Creation constructors. */
    Turbulent_Velocity_Field_Multi_Block_List(void) : 
       NBlk(0), NBlk_Idir(0), NBlk_Jdir(0), NBlk_Kdir(0), Vel_Blks(NULL), Allocated(0) { }

    Turbulent_Velocity_Field_Multi_Block_List(const int N) {
       Allocate(N);
    }

    Turbulent_Velocity_Field_Multi_Block_List(const int Ni, 
                                              const int Nj, 
                                              const int Nk) {
       Allocate(Ni, Nj, Nk);
    }

    /* Destructor. */
    ~Turbulent_Velocity_Field_Multi_Block_List(void) {
        Deallocate();
    }

    /* Other member functions  */

    void Allocate(const int Ni, const int Nj, const int Nk);

    void Allocate(const int N);

    void Deallocate(void);

    void Create(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                const Grid3D_Input_Parameters &Input);

    void Interpolate_Turbulent_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
				     Turbulent_Velocity_Field_Multi_Block_List &Interpolated_Velocity_Field);
    

  private:
    //copy and assignment are not permitted
    Turbulent_Velocity_Field_Multi_Block_List(const Turbulent_Velocity_Field_Multi_Block_List &V);
    Turbulent_Velocity_Field_Multi_Block_List &operator = (const Turbulent_Velocity_Field_Multi_Block_List &V);
};

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Allocate -- Allocate memory.   *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Allocate(const int Ni, 
                                                                const int Nj, 
                                                                const int Nk) {
   if (Ni >= 1 && Nj >= 1 && Nk >= 1 && !Allocated) {
      NBlk_Idir = Ni; 
      NBlk_Jdir = Nj; 
      NBlk_Kdir = Nk; 
      NBlk = Ni*Nj*Nk;
      Vel_Blks = new Turbulent_Velocity_Field_Block[NBlk];
      Allocated = 1;
   } /* endif */
}

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Allocate -- Allocate memory.   *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Allocate(const int N) {
   if (N >= 1 && !Allocated) {
      NBlk_Idir = N; 
      NBlk_Jdir = 1; 
      NBlk_Kdir = 1; 
      NBlk = N;
      Vel_Blks = new Turbulent_Velocity_Field_Block[NBlk];
      Allocated = 1;
   } /* endif */
}

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Allocate -- Deallocate memory. *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Deallocate(void) {
   if (NBlk >= 1 && Allocated) {
       delete []Vel_Blks;
       Vel_Blks = NULL;
       NBlk_Idir = 0; 
       NBlk_Jdir = 0; 
       NBlk_Kdir = 0;
       NBlk = 0;
       Allocated = 0;
   } /* endif */
}

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Create --                      *
 *       Create memory containers for velocity field data.                   *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Create(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                                              const Grid3D_Input_Parameters &Input) {
   if (Initial_Mesh.NBlk >= 1 && Initial_Mesh.Allocated) {
      Allocate(Initial_Mesh.NBlk_Idir, Initial_Mesh.NBlk_Jdir, Initial_Mesh.NBlk_Kdir);
      for (int nBlk = 0; nBlk <= NBlk-1; ++nBlk ) {
	if (Initial_Mesh.Grid_Blks[nBlk].Allocated) Vel_Blks[nBlk].allocate(Initial_Mesh.Grid_Blks[nBlk].NCi-
                                                                            2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                            Initial_Mesh.Grid_Blks[nBlk].NCj-
                                                                            2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                            Initial_Mesh.Grid_Blks[nBlk].NCk-
                                                                            2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                            Initial_Mesh.Grid_Blks[nBlk].Nghost);
      } /* endfor */
   } /* endif */
}


/********************************************************
 * Routine: Inflow_Turbulence_XY_Plane                  *
 *                                                      *
 * Feeds turbulent velocity fluctuations on a XY-plane  *
 * of a 3D grid. It applies a turbulent inflow BC.      *
 *                                                      *
 ********************************************************/
template<typename HEXA_BLOCK, typename SOLN_pSTATE, typename SOLN_cSTATE>
void Inflow_Turbulence_XY_Plane(HEXA_BLOCK &Solution_Block,
				const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field,
				const Input_Parameters<SOLN_pSTATE,SOLN_cSTATE> &IPs, 
				const double &Time) {

  int N(1);
  double Lz(IPs.Grid_IP.Turbulence_Box_Length);
  double L_convected = double(N)*Lz - Time*IPs.Mean_Velocity.z;

  if ( L_convected < ZERO ) {
    while ( L_convected < ZERO ) {
      N++;
      L_convected = double(N)*Lz - Time*IPs.Mean_Velocity.z;
    }
  }
  
  int nBlk, k, ii, jj, kk;
  double new_z, delta, dmin, dmin1, delta1, zmax, zmin, delta_z;  	      
  Vector3D Position_on_Slice, dX, dVdx, dVdy, dVdz;

  for (int n_ghost = 1; n_ghost <= Solution_Block.Nghost; ++n_ghost) {
    // loop over the cells on the boundary  of the solution block
    for (int i = Solution_Block.ICl; i <= Solution_Block.ICu; ++i) {
      for (int j = Solution_Block.JCl; j <= Solution_Block.JCu; ++j) {
	if (Solution_Block.Grid.BCtypeT[i][j] == BC_INFLOW_TURBULENCE ||
	    Solution_Block.Grid.BCtypeB[i][j] == BC_INFLOW_TURBULENCE ) {
	  
	  // Top boundary
	  if (Solution_Block.Grid.BCtypeT[i][j]) {
	    
	    k = Solution_Block.KCu + n_ghost;

	    delta_z = fabs(Solution_Block.Grid.Cell[i][j][k].Xc.z 
			   - Solution_Block.Grid.Cell[i][j][Solution_Block.KCu].Xc.z)
	              - Solution_Block.Grid.volume(i, j, Solution_Block.KCu)/
	                (Solution_Block.Grid.AfaceTop(i, j, Solution_Block.KCu) 
			 + Solution_Block.Grid.AfaceBot(i, j, Solution_Block.KCu));

	  // Bottom boundary
	  } else if (Solution_Block.Grid.BCtypeB[i][j]) {

	    k = Solution_Block.KCl - n_ghost;
	    
	    delta_z = fabs(Solution_Block.Grid.Cell[i][j][k].Xc.z 
			   - Solution_Block.Grid.Cell[i][j][Solution_Block.KCl].Xc.z)
	              - Solution_Block.Grid.volume(i, j, Solution_Block.KCl)/
	                (Solution_Block.Grid.AfaceTop(i, j, Solution_Block.KCl) 
			 + Solution_Block.Grid.AfaceBot(i, j, Solution_Block.KCl));
	  }

	  new_z = L_convected - delta_z;

	  Solution_Block.W[i][j][k].v = IPs.Mean_Velocity;


	  // find the appropriate turbulent field block to interpolate the fluctuations
	  for (nBlk = 0; nBlk < Velocity_Field.NBlk; ++nBlk) {
	    zmax = max(Velocity_Field.Vel_Blks[nBlk].Node_INl_JNl_KNl.z, 
		       Velocity_Field.Vel_Blks[nBlk].Node_INu_JNu_KNu.z);
	    zmin = min(Velocity_Field.Vel_Blks[nBlk].Node_INl_JNl_KNl.z, 
		       Velocity_Field.Vel_Blks[nBlk].Node_INu_JNu_KNu.z);

	    if ( new_z >= zmin  &&  new_z <= zmax ) break;

	    if ( (nBlk == Velocity_Field.NBlk-1) && (new_z < zmin  ||  new_z > zmax) ) {
	      cout << "\n --> n_ghost = " << n_ghost  << " i = " << i 
		   << " j = " << j << " k = " << k << flush << endl;
	      cerr << "\n new_z = " << new_z << " not found in turbulence blocks";
	      exit(1);
	    }   
	  }


	  //-----------------------------------------------------------------------------------
	  // determine ii, jj, kk of the velocity field block to interpolate the fluctuations
	  //-----------------------------------------------------------------------------------


	  // find index in k-direction for the closest point to new_z of the slice
	  // ---> This should be more general. 
	  dmin = 1E10;
	  for (int n = Velocity_Field.Vel_Blks[nBlk].KCl; n <= Velocity_Field.Vel_Blks[nBlk].KCu; ++n) {
	    delta = fabs(Velocity_Field.Vel_Blks[nBlk].Position[ Velocity_Field.Vel_Blks[nBlk].ICl ]
			                                       [ Velocity_Field.Vel_Blks[nBlk].JCl ]
			                                       [n].z - new_z);
	    if ( delta < dmin ) { 
	      dmin = delta;
	      kk = n;
	    }
	    if (dmin == 0.0) break;
	  }


	  dmin  = 1E10;
	  dmin1 = 1E10;
	  for (int p = Velocity_Field.Vel_Blks[nBlk].ICl; p <= Velocity_Field.Vel_Blks[nBlk].ICu; ++p) {
	    for (int q = Velocity_Field.Vel_Blks[nBlk].JCl; q <= Velocity_Field.Vel_Blks[nBlk].JCu; ++q) {

	      // find index in i-direction for the closest point to x of the slice
	      delta = fabs(Velocity_Field.Vel_Blks[nBlk].Position[p][q][kk].x 
			   - Solution_Block.Grid.Cell[i][j][k].Xc.x);

	      if ( delta < dmin )  { 
		dmin = delta;
		ii = p;
	      }
    
	      // find index in j-direction for the closest point to y of the slice
	      delta1 = fabs(Velocity_Field.Vel_Blks[nBlk].Position[p][q][kk].y 
			   - Solution_Block.Grid.Cell[i][j][k].Xc.y);

	      if ( delta1 < dmin1 )  { 
		dmin1 = delta1;
		jj = q;
	      }

	    }
	  }


	  //------------------------------------------------------------------
	  // use least squares to reconstruct the turbulent velocity field 
	  //
	  // ---> Gradients can be computed once and be stored!!!
	  //
	  //------------------------------------------------------------------
	
	  Velocity_Field.Vel_Blks[nBlk].LeastSquares_Reconstruction(ii, jj, kk, 
								    dVdx, dVdy, dVdz);
           
	  Position_on_Slice.x = Solution_Block.Grid.Cell[i][j][k].Xc.x;  
	  Position_on_Slice.y = Solution_Block.Grid.Cell[i][j][k].Xc.y;  
	  Position_on_Slice.z = new_z;

	  dX = Position_on_Slice - Velocity_Field.Vel_Blks[nBlk].Position[ii][jj][kk];
	
	  Solution_Block.W[i][j][k].v += Velocity_Field.Vel_Blks[nBlk].Velocity[ii][jj][kk]
	                                 + dVdx*dX.x + dVdy*dX.y + dVdz*dX.z;

	} /* end if */
      } /* end for */
    } /* end for */ 
  } /* end for */	  



  //---------------
  // ghost cells 
  //---------------

  // i < ICl
  for (int n_ghost = 1; n_ghost <= Solution_Block.Nghost; ++n_ghost) {
    for (int i = Solution_Block.ICl-Solution_Block.Nghost; i < Solution_Block.ICl; ++i) {
      for (int j = Solution_Block.JCl-Solution_Block.Nghost; j <= Solution_Block.JCu+Solution_Block.Nghost; ++j) {
	// Top boundary
	if (Solution_Block.Grid.BCtypeT[i][j]) {
	  k = Solution_Block.KCu + n_ghost;
	// Bottom boundary
	} else if (Solution_Block.Grid.BCtypeB[i][j]) {
	  k = Solution_Block.KCl - n_ghost;
	}
	
	if (j >= Solution_Block.JCl  &&  j <= Solution_Block.JCu) {
	  Solution_Block.W[i][j][k].v = Solution_Block.W[Solution_Block.ICl][j][k].v;
	// corners
	} else if (j < Solution_Block.JCl) {
	  Solution_Block.W[i][j][k].v = Solution_Block.W[Solution_Block.ICl][Solution_Block.JCl][k].v;
	} else if (j > Solution_Block.JCu) {
	  Solution_Block.W[i][j][k].v = Solution_Block.W[Solution_Block.ICl][Solution_Block.JCu][k].v;
	}

      }
    }
  }

  // i > ICu 
  for (int n_ghost = 1; n_ghost <= Solution_Block.Nghost; ++n_ghost) {
    for (int i = Solution_Block.ICu+1; i <= Solution_Block.ICu+Solution_Block.Nghost; ++i) {
      for (int j = Solution_Block.JCl-Solution_Block.Nghost; j <= Solution_Block.JCu+Solution_Block.Nghost; ++j) {
	// Top boundary
	if (Solution_Block.Grid.BCtypeT[i][j]) {
	  k = Solution_Block.KCu + n_ghost;
	// Bottom boundary
	} else if (Solution_Block.Grid.BCtypeB[i][j]) {
	  k = Solution_Block.KCl - n_ghost;
	}

	if (j >= Solution_Block.JCl  &&  j <= Solution_Block.JCu) {
	  Solution_Block.W[i][j][k].v = Solution_Block.W[Solution_Block.ICu][j][k].v;
	// corners
	} else if (j < Solution_Block.JCl) {
	  Solution_Block.W[i][j][k].v = Solution_Block.W[Solution_Block.ICu][Solution_Block.JCl][k].v;
	} else if (j > Solution_Block.JCu) {
	  Solution_Block.W[i][j][k].v = Solution_Block.W[Solution_Block.ICu][Solution_Block.JCu][k].v;
	}

      }
    }
  }

  // j < JCl
  for (int n_ghost = 1; n_ghost <= Solution_Block.Nghost; ++n_ghost) {
    for (int j = Solution_Block.JCl-Solution_Block.Nghost; j < Solution_Block.JCl; ++j) {
      for (int i = Solution_Block.ICl; i <= Solution_Block.ICu; ++i) {
	// Top boundary
	if (Solution_Block.Grid.BCtypeT[i][j]) {
	  k = Solution_Block.KCu + n_ghost;
	// Bottom boundary
	} else if (Solution_Block.Grid.BCtypeB[i][j]) {
	  k = Solution_Block.KCl - n_ghost;
	}

	Solution_Block.W[i][j][k].v = Solution_Block.W[i][Solution_Block.JCl][k].v;

      }
    }
  }
  
  // j > JCu 
  for (int n_ghost = 1; n_ghost <= Solution_Block.Nghost; ++n_ghost) {
    for (int j = Solution_Block.JCu+1; j <= Solution_Block.JCu+Solution_Block.Nghost; ++j) {
      for (int i = Solution_Block.ICl; i <= Solution_Block.ICu; ++i) {
	// Top boundary
	if (Solution_Block.Grid.BCtypeT[i][j]) {
	  k = Solution_Block.KCu + n_ghost;
	// Bottom boundary
	} else if (Solution_Block.Grid.BCtypeB[i][j]) {
	  k = Solution_Block.KCl - n_ghost;
	}

	Solution_Block.W[i][j][k].v = Solution_Block.W[i][Solution_Block.JCu][k].v;

      }
    }
  }
  

  // update SFS turbulence kinetic energy and conservative state variables
  for (int n_ghost = 1; n_ghost <= Solution_Block.Nghost; ++n_ghost) {
    for (int i = Solution_Block.ICl-Solution_Block.Nghost; i <= Solution_Block.ICu+Solution_Block.Nghost; ++i) {
      for (int j = Solution_Block.JCl-Solution_Block.Nghost; j <= Solution_Block.JCu+Solution_Block.Nghost; ++j) {
	if (Solution_Block.Grid.BCtypeT[i][j] == BC_INFLOW_TURBULENCE ||
	    Solution_Block.Grid.BCtypeB[i][j] == BC_INFLOW_TURBULENCE ) {
	  
	  // Top boundary
	  if (Solution_Block.Grid.BCtypeT[i][j]) {
	    k = Solution_Block.KCu + n_ghost;
	  // Bottom boundary
	  } else if (Solution_Block.Grid.BCtypeB[i][j]) {
	    k = Solution_Block.KCl - n_ghost;
	  }
	  
	  // reconstruct gradients
	  Solution_Block.Linear_Reconstruction_LeastSquares(i, j, k, IPs.i_Limiter);

	  // update k
	  Solution_Block.W[i][j][k].k = 0.005*sqr(Solution_Block.W[i][j][k].filter_width(Solution_Block.Grid.volume(i,j,k))*
						  Solution_Block.W[i][j][k].abs_strain_rate(Solution_Block.dWdx[i][j][k],
											    Solution_Block.dWdy[i][j][k],
											    Solution_Block.dWdz[i][j][k]));

	  // update the conservative state variables
	  Solution_Block.U[i][j][k] = Solution_Block.W[i][j][k].U();

	}
      }
    }
  }


}


template<typename HEXA_BLOCK, typename SOLN_pSTATE, typename SOLN_cSTATE>
int Inflow_Turbulence_XY_Plane(HEXA_BLOCK *Solution_Block,
			       const AdaptiveBlock3D_List &LocalSolnBlockList,
			       const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field,
			       const Input_Parameters<SOLN_pSTATE,SOLN_cSTATE> &IPs,
			       const double &Time) {
  int error_flag(0);

  for (int Nblk = 0 ; Nblk <= LocalSolnBlockList.Nblk-1 ; Nblk++ ) {
    if (LocalSolnBlockList.Block[Nblk].used == ADAPTIVEBLOCK3D_USED) {
      Inflow_Turbulence_XY_Plane(Solution_Block[Nblk], Velocity_Field, IPs, Time);
    }
  }

  return error_flag;

}


/********************************************************
 * Routine: IC_Assign_Turbulence_Fresh_Gas              *
 *                                                      *
 * Assigns turbulent velocity fluctuations to a volume  *
 * of a 3D domain. It is used with the slot and Bunsen  *  
 * burner configurations.                               *                         
 *                                                      *
 ********************************************************/
template<typename HEXA_BLOCK, typename SOLN_pSTATE, typename SOLN_cSTATE>
void IC_Assign_Turbulence_Fresh_Gas(HEXA_BLOCK &Solution_Block,
				    const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field,
				    const Input_Parameters<SOLN_pSTATE,SOLN_cSTATE> &IPs) {

  int nnBlk, ii, jj, kk, n;
  double delta, dmin;
  double xmax, xmin, ymax, ymin, zmax, zmin;
  Vector3D dVdx, dVdy, dVdz, dX, local_X;
  

  for (int i = Solution_Block.ICl; i <= Solution_Block.ICu; ++i) {
    for (int j = Solution_Block.JCl; j <= Solution_Block.JCu; ++j) {
      for (int k = Solution_Block.KCl; k <= Solution_Block.KCu; ++k) {

	local_X = Solution_Block.Grid.Cell[i][j][k].Xc;

	
	if (// Slot burner
	    (IPs.Grid_IP.i_Grid == GRID_BUNSEN_BOX  &&  
	     (local_X.z <= 0.02  &&  fabs(local_X.y) <= 0.0125)) ||
	    // Bunsen burner
	    (IPs.Grid_IP.i_Grid == GRID_BUNSEN_BURNER  &&  
	     (local_X.z <= 0.035  &&  (sqr(local_X.x) + sqr(local_X.y) <= sqr(0.0056))))
	    ) {
	 
	  // find  nnBlk, ii, jj and kk to perform the interpolation

	  for (nnBlk = 0; nnBlk < Velocity_Field.NBlk; ++nnBlk) {
	    xmax = max(Velocity_Field.Vel_Blks[nnBlk].Node_INl_JNl_KNl.x, 
		       Velocity_Field.Vel_Blks[nnBlk].Node_INu_JNu_KNu.x);
	    xmin = min(Velocity_Field.Vel_Blks[nnBlk].Node_INl_JNl_KNl.x, 
		       Velocity_Field.Vel_Blks[nnBlk].Node_INu_JNu_KNu.x);

	    ymax = max(Velocity_Field.Vel_Blks[nnBlk].Node_INl_JNl_KNl.y, 
		       Velocity_Field.Vel_Blks[nnBlk].Node_INu_JNu_KNu.y);
	    ymin = min(Velocity_Field.Vel_Blks[nnBlk].Node_INl_JNl_KNl.y, 
		       Velocity_Field.Vel_Blks[nnBlk].Node_INu_JNu_KNu.y);

	    zmax = max(Velocity_Field.Vel_Blks[nnBlk].Node_INl_JNl_KNl.z, 
		       Velocity_Field.Vel_Blks[nnBlk].Node_INu_JNu_KNu.z);
	    zmin = min(Velocity_Field.Vel_Blks[nnBlk].Node_INl_JNl_KNl.z, 
		       Velocity_Field.Vel_Blks[nnBlk].Node_INu_JNu_KNu.z);

	    if ( (local_X.x >= xmin  &&  local_X.x <= xmax) &&
		 (local_X.y >= ymin  &&  local_X.y <= ymax) &&
		 (local_X.z >= zmin  &&  local_X.z <= zmax) ) {

	      break;
	    }

	  } /* end for*/ 

	  
	  
	  // search in X-direction
	  dmin = 1E10;
	  for (n = Velocity_Field.Vel_Blks[nnBlk].ICl; n <= Velocity_Field.Vel_Blks[nnBlk].ICu; ++n) {
	    delta = fabs(Velocity_Field.Vel_Blks[nnBlk].Position[n]
			 [Velocity_Field.Vel_Blks[nnBlk].JCl]
			 [Velocity_Field.Vel_Blks[nnBlk].KCl].x -
			 local_X.x);
	    if ( delta < dmin )  { 
	      dmin = delta;
	      ii = n;
	    }
	    if (dmin == 0.0) break;
	  } /* end for */


	  // search in Y-direction
	  dmin = 1E10; 
	  for (n = Velocity_Field.Vel_Blks[nnBlk].JCl; n <= Velocity_Field.Vel_Blks[nnBlk].JCu; ++n) {
	    delta = fabs(Velocity_Field.Vel_Blks[nnBlk].Position[Velocity_Field.Vel_Blks[nnBlk].ICl]
			 [n]
			 [Velocity_Field.Vel_Blks[nnBlk].KCl].y -
			 local_X.y);
	    if ( delta < dmin )  { 
	      dmin = delta;
	      jj = n;
	    }
	    if (dmin == 0.0) break;
	  } /* end for */


	  // search in Z-direction
	  dmin = 1E10; 
	  for (n = Velocity_Field.Vel_Blks[nnBlk].KCl; n <= Velocity_Field.Vel_Blks[nnBlk].KCu; ++n) {
	    delta = fabs(Velocity_Field.Vel_Blks[nnBlk].Position[Velocity_Field.Vel_Blks[nnBlk].ICl]
			 [Velocity_Field.Vel_Blks[nnBlk].JCl]
			 [n].z - 
			 local_X.z);
	    if ( delta < dmin )  { 
	      dmin = delta;
	      kk = n;
	    }
	    if (dmin == 0.0) break;
	  } /* end for */


	  if (ii < Velocity_Field.Vel_Blks[nnBlk].ICl || ii > Velocity_Field.Vel_Blks[nnBlk].ICu || 
	      jj < Velocity_Field.Vel_Blks[nnBlk].JCl || jj > Velocity_Field.Vel_Blks[nnBlk].JCu || 
	      kk < Velocity_Field.Vel_Blks[nnBlk].KCl || kk > Velocity_Field.Vel_Blks[nnBlk].KCu) {
	    cout << "\n Index out of bound!!! -> ii = " << ii 
		 << "  jj = " << jj << "  kk = " << kk
		 << "\nICu = " << Velocity_Field.Vel_Blks[nnBlk].ICu 
		 << "  JCu = " << Velocity_Field.Vel_Blks[nnBlk].JCu 
		 << "  KCu = " << Velocity_Field.Vel_Blks[nnBlk].KCu; 
	  }


	  // use least squares to reconstruct the turbulent velocity field 
	  Velocity_Field.Vel_Blks[nnBlk].LeastSquares_Reconstruction(ii, jj, kk, 
								     dVdx, dVdy, dVdz);

	  dX = local_X - Velocity_Field.Vel_Blks[nnBlk].Position[ii][jj][kk];

	  Solution_Block.W[i][j][k].v +=  Velocity_Field.Vel_Blks[nnBlk].Velocity[ii][jj][kk] 
	                                  + dVdx*dX.x + dVdy*dX.y + dVdz*dX.z;

	  // update the conservative state variables
	  Solution_Block.U[i][j][k] = Solution_Block.W[i][j][k].U();

	} /* endif */

      } /* endfor */
    } /* endfor */
  } /* endfor */


}


template<typename HEXA_BLOCK, typename SOLN_pSTATE, typename SOLN_cSTATE>
int IC_Assign_Turbulence_Fresh_Gas(HEXA_BLOCK *Solution_Block,
				   const AdaptiveBlock3D_List &LocalSolnBlockList,
				   const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field,
				   const Input_Parameters<SOLN_pSTATE,SOLN_cSTATE> &IPs) {
  int error_flag(0);

  for (int Nblk = 0 ; Nblk <= LocalSolnBlockList.Nblk-1 ; Nblk++ ) {
    if (LocalSolnBlockList.Block[Nblk].used == ADAPTIVEBLOCK3D_USED) {
      IC_Assign_Turbulence_Fresh_Gas(Solution_Block[Nblk], Velocity_Field, IPs);
    }
  }

  return error_flag;

}




/*!
 * Class: K_container
 *
 * \brief Class defined to contain spectral information at one wave number
 *
 * A class to contain all the spectral information at one wave number.
 * Multiple combinations of (k1,k2,k3) lead to the same abs_k. Therefore needs
 * container will hold all the combinations of (k1,k2,k3) that give this abs_k as
 * an index. [ index = l + nz*( j + Ny*i ) ] \n
 * Relational operators and output operator are required to make use of a
 * binary tree sorting algorithm. (see class BinaryTree<itemtype>)
 */
class K_container {
public:
    K_container() : k(0), N(1), Ek(0), velocity_scaling_factor(1), indexcounter(0) {  }
    double k;
    int N;
    int indexcounter;
    int *indexes;
    int index;
    double Ek;
    double Emodel;
    double velocity_scaling_factor;
    // Relational operators
    //! W == W
    friend bool operator ==(const K_container &K1, const K_container &K2) {
        return ((fabs((K1.k - K2.k)/(K1.k+NANO)) <= TOLER));
    }
    //! W != W
    friend bool operator !=(const K_container &K1, const K_container &K2) {
        return (K1.k != K2.k);
    }
    //! W <= W
    friend bool operator <=(const K_container &K1, const K_container &K2) {
        return (K1.k <= K2.k);
    }
    //! W >= W
    friend bool operator >=(const K_container &K1, const K_container &K2) {
        return (K1.k >= K2.k);
    }
    //! W < W
    friend bool operator <(const K_container &K1, const K_container &K2) {
        return (K1.k < K2.k);
    }
    //! W > W
    friend bool operator >(const K_container &K1, const K_container &K2) {
        return (K1.k > K2.k);
    }
    
    // Input-output operators. 
    //! ostream << W
    friend ostream& operator << (ostream &out_file, const K_container &K) {
        out_file.setf(ios::scientific);
        out_file << " " << K.k  << "     " << K.N << "     " << K.Ek;
        out_file.unsetf(ios::scientific);
        return (out_file);
    }
    //! istream >> W
    //friend istream& operator >> (istream &in_file,  K_container &K); 
};


/*!
 * Class: K_BinaryTree
 *
 * \brief This class inherits from BinaryTree<itemtype> to adapt to perform on K_containers
 *
 * A class that inherits from BinaryTree<itemtype>. The functioncall InsertNode
 * is overloaded so that an item with the same value of k as an existing node 
 * won't generate a new node in the tree but will instead store the index
 * of the item in the existing node
 *
 * This class is needed in the class RandomFieldRogallo to efficiently sort and 
 * store all the K_containers.
 */
class K_BinaryTree : public BinaryTree<K_container> {
public:
    
    K_BinaryTree( ) : BinaryTree<K_container>( ) { }
    
    
    void InsertNode(K_container &newItem) {
        InsertNode(K_BinaryTree::root, newItem);
    }
    
    void InsertNode(treeNode<K_container> *&root, K_container &newItem) {
        int MAX_INDEXES = 500;
        if ( root == NULL ) {
            root = new treeNode<K_container>( newItem );
            root->item.indexes = new int[MAX_INDEXES];
            root->item.indexes[0] = newItem.index;
            root->item.indexcounter++;
            return;
        }
        else if ( newItem == root->item ) {
            int i = root->item.indexcounter;
            root->item.indexes[i] = newItem.index;
            root->item.indexcounter++;
            root->item.N++;
            if (root->item.indexcounter==MAX_INDEXES)
                cerr << "\n Maximum amount of indexes reached";
        }
        else if ( newItem < root->item ) {
            InsertNode( root->left, newItem );
        }
        else {
            InsertNode( root->right, newItem );
        }
    }
};


/*!
 * Class: RandomFieldRogallo
 *
 * \brief Class defined to generate random fluctuations using
 * Rogallo's procedure.
 *
 * A class to generate a turbulence field with Rogallo's procedure. \n
 * Parameters to control in inputfile are:
 *  - Spectrum_Model:  the model spectrum
 *  - Domain_Integral_Lengthscale_Ratio: to control the amount oflargest eddies in the domain.
 *    As a general rule of thumb for isotropic turbulence : L11 = 0.4 L where L is the specified lengthscale
 *  - Turbulent_Kinetic_Energy: It should be noted that this is a theoretical value
 *    on a continuous spectrum integrated from 0 to infinity
 *
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
class RandomFieldRogallo {   
public:
    string spectrum_name;
    int spectrum_flag;  //!< Turbulence kinetic energy spectrum flag.
    double L1, L2, L3;  //!< Grid dimensions
    int Nx, Ny, Nz;     //!< Number of grid nodes
    double Ls;          //!< Smallest domain length
    double Lp;          //!< Length scale with highest energy
    double k_L;         //!< wave number with highest energy 
    double TKE;         //!< Turbulent kinetic energy for continuous spectrum
    double nu;          //!< Kinematic viscosity
    double ReL;         //!< Turbulence Reynolds Number
    double eta;         //!< Kolmogorov scale
    double k_eta;       //!< Kolmogorov wave number
    double eps;         //!< Dissipation
    double up;          //!< RMS velocity
    double Rlambda;     //!< Taylor scale Reynolds Number
    double lambda_g;    //!< Taylor scale
    
    const char *File_Name;

    
    //! Constructor
    RandomFieldRogallo(Input_Parameters<SOLN_pSTATE,SOLN_cSTATE> &IPs) {
        File_Name = IPs.Output_File_Name;
        spectrum_flag = IPs.Turbulence_IP.i_spectrum;
        spectrum_name = IPs.Turbulence_IP.spectrum;
	if (IPs.Grid_IP.i_Grid == GRID_BUNSEN_BURNER || IPs.Grid_IP.i_Grid == GRID_BUNSEN_BOX) {
	  L1 = IPs.Grid_IP.Turbulence_Box_Width;
	  L2 = IPs.Grid_IP.Turbulence_Box_Height;
	  L3 = IPs.Grid_IP.Turbulence_Box_Length;
	} else {
	  L1 = IPs.Grid_IP.Box_Width;
	  L2 = IPs.Grid_IP.Box_Height;
	  L3 = IPs.Grid_IP.Box_Length;
	}
        Ls = min(min(L1,L2),L3);
        Lp = Ls/IPs.Turbulence_IP.LLR;
        TKE = IPs.Turbulence_IP.TKE;
        nu = IPs.Wo.nu();

        //Calculate :
        up  = sqrt(TKE*TWO/THREE);         
        eps = pow(TKE,THREE/TWO)/Lp;      
        ReL = sqrt(TKE)*Lp/nu;              
        eta = pow(ReL,-THREE/FOUR)*Lp;    
        k_eta = TWO*PI/eta;
        k_L = TWO*PI/Lp;
        Rlambda = sqrt(TWENTY/THREE * ReL); 
        lambda_g = nu*Rlambda/up;           
    }  

    //! k1 = n1 * 2*PI/L1
    double k_1(const int &n1) const {
        if (n1<=Nx/2)
            return  n1 * TWO*PI/L1; 
        else
            return  (Nx-n1) * TWO*PI/L1;
    }
    //! k2 = n2 * 2*PI/L2
    double k_2(const int &n2) const {
        if (n2<=Ny/2)
            return  n2 * TWO*PI/L2; 
        else
            return  (Ny-n2) * TWO*PI/L2;
    }
    //! k3 = n3 * 2*PI/L3
    double k_3(const int &n3) const {
        if (n3<=Nz/2)
            return  n3 * TWO*PI/L3; 
        else
            return  (Nz-n3) * TWO*PI/L3;
    }
    
    //! Generates a random double
    double random_double() const { return  drand48(); }

    //! The model for the energy spectrum
    double Energy_Spectrum_Value(const double &abs_wave_num) const;
    
/** @name Rogallo procedure functions */
/*        --------------------------- */    
//@{
    //! alpha as defined by Rogallo, 1981 for continuous spectra
    complex<double> alpha(const double &abs_wave_num, 
			  const double &theta1, 
			  const double &phi) const;
    
    //! beta as defined by Rogallo, 1981 for continuous spectra
    complex<double> beta(const double &abs_wave_num, 
			 const double &theta2, 
			 const double &phi) const;
    
    //! alpha as defined by Rogallo, 1981 for discrete spectra
    complex<double> alpha_discrete(const double &abs_wave_num, 
                                   const double &N,
                                   const double &theta1, 
                                   const double &phi) const;
    
    //! beta as defined by Rogallo, 1981 for discrete spectra
    complex<double> beta_discrete(const double &abs_wave_num, 
                                  const double &N,
                                  const double &theta2, 
                                  const double &phi) const;

    //! Create and assign a homogeneous turbulence velocity field
    int Create_Homogeneous_Turbulence_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                                     const Grid3D_Input_Parameters &IPs,
                                                     int &batch_flag,
                                                     Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field);
//@}
    
/** @name Write and read a Velocity Field */
/*        ------------------------------- */    
//@{
    void Write_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                        const Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
                                        const Grid3D_Input_Parameters &IPs) const;
       
    void Read_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                       Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
                                       Grid3D_Input_Parameters &IPs);
//@}
    
/** @name Output the spectrum to a gnuplot file */
/*        ------------------------------------- */    
//@{
    ofstream Turbulence_Spectrum_File;
    int Open_Turbulence_Spectrum_File( );
    int Output_Turbulence_Spectrum_to_File(const double &k,
                                           const double &model_Energy,
                                           const double &constructed_Energy);
    int Close_Turbulence_Spectrum_File( );
//@}    

};

//-----------------------------------------------------//
//         Members of RandomFieldRogallo class         //
//-----------------------------------------------------//


// Prescribed energy spectrum
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
double RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Energy_Spectrum_Value(const double &abs_wave_num) const {

    double k, kp, kd, u, EE;
    double C, A, alpha, a_s;
    int s;
    
    double epsVKP;
 
    // variable declaration for Pope:
    double p0, cL, beta, c_eta, shape_function_fL, shape_function_f_eta;
    
    k = abs_wave_num;
  
    switch (spectrum_flag) {
            
        /*****  Lee and Reynolds  *****/
        case SPECTRUM_LEE_REYNOLDS :
            C = 20.0;  kp = 80.0;   
            EE = (k <= kp)  ?  C*k*k : C*kp*kp*pow(k/kp, -5.0/3.0);
            break;

        /*****  Laval and Nazarenko paper  *****/
        case SPECTRUM_LAVAL_NAZARENKO :
            //double EE, kp, C;
            C = 1.0; // kp = 4.0;
            kp = TWO*PI/Lp;
            EE = C*k*exp(-pow(k/kp, 2.0));
            break;

        /*****   von Karman-Pao   *****/
        case SPECTRUM_VON_KARMAN_PAO :
            A = 1.5;    alpha = 1.5;
         
            kp = TWO*PI/Lp;
            kd = 362.0*kp 
            /* kp = 1, kd = 362.0*/;
            u = 0.8;    
            epsVKP = 2.73E-3; //2.73E-03;
            EE = (A/epsVKP)*pow(u, 5.0);
            EE = EE*pow(k/kp, 4.0) * exp(-3.0*alpha*pow(k/kd, 4.0/3.0)/2.0);
            EE = EE/pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
            break;

        /*****  Haworth and Poinsot paper  *****/
        case SPECTRUM_HAWORTH_POINSOT :
            //double EE, kp, 
            u = up;
            kp = TWO*PI/Lp;
            EE = (32.0/3.0) * sqrt(2.0/PI)* (u*u/kp) * pow(k/kp, 4.0) * exp(-2.0*(k/kp)*(k/kp));
            break;

        /*****  Chasnov paper 1996  *****/
        case SPECTRUM_CHASNOV :  
            //double a_s, EE, kp = 4.0 /*8.0 20.0 4.0*/, u = 0.095; /*0.1 28.3 21.21  0.001*/  
            //int s = 3;
            s = 3;
            kp = 4.0;
            u = 0.009; 
            // u = 0.001  ->  Re_lambda = 10 
            // u = 0.009  ->  Re_lambda = 98
            // u = 0.095  ->  Re_lambda = 950
   
            a_s = pow(2.0*double(s)+1.0, double(s)+1.0)/(factorial(s)*pow(2.0, double(s)));
            EE = (HALF*a_s*u*u/kp)*pow(k/kp, 2.0*double(s)+1.0);
            EE = EE*exp(-(double(s)+HALF)*pow(k/kp, 2.0));
            break;

        /*****   Bell & Day report   *****/
        case SPECTRUM_BELL_DAY :
            //  kd = 1/(2*dx)
            kp = 3.0;   kd = ONE/0.576E-3;
            EE = pow(k/kp, 4.0) * exp(-9.0*pow(k/kd, 4.0/3.0)/4.0);
            EE /= pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
            break;

        /*****   Pope   "Turbulent Flows" , Cambridge University Press *****/
        case SPECTRUM_POPE :
        default :
            // shape function for the large eddy wave number:
            p0 = 2.0;  // --- p0 = 4 gives same slope as "Von Karman Pao" in energy containing range
            cL = 6.78;
            shape_function_fL = pow( (k*Lp) / sqrt( sqr(k*Lp) + cL) , FIVE/THREE+p0 );
            
            // shape function for the kolmogorov wave number:
            beta = 5.2;
            c_eta = 0.40;
            shape_function_f_eta = exp(-beta * ( pow( pow(k*eta,FOUR) + pow(c_eta,FOUR) , ONE/FOUR) - c_eta) );
            
            // the resulting energy spectrum
            C = 1.5;
            EE = C * pow(eps,TWO/THREE) * pow(k,-FIVE/THREE) * shape_function_fL * shape_function_f_eta;
        
            break;    

    }  // end switch

    return (k == 0.0)  ?  0.0 : EE;
}


// alpha
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
alpha(const double &abs_wave_num, const double &theta1, const double &phi) const {
    
    double E, k;
    k = abs_wave_num;  E = Energy_Spectrum_Value(k);
    return (k == 0.0)  ?  0.0 : sqrt(TWO*E/(4.0*PI*k*k)) * exp(I*theta1) * cos(phi);
}

// beta
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
beta(const double &abs_wave_num, const double &theta2, const double &phi) const {
    
    double E, k;
    k = abs_wave_num;  E = Energy_Spectrum_Value(k);
    return (k == 0.0)  ?  0.0 : sqrt(TWO*E/(4.0*PI*k*k)) * exp(I*theta2) * sin(phi);    
}

// alpha
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
alpha_discrete(const double &abs_wave_num, const double &N, const double &theta1, const double &phi) const {
    
    double E, k;
    k = abs_wave_num;  E = Energy_Spectrum_Value(k);
    return (k == 0.0)  ?  0.0 : sqrt(E/N) * exp(I*theta1) * cos(phi);
}

// beta
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
beta_discrete(const double &abs_wave_num, const double &N, const double &theta2, const double &phi) const {
    
    double E, k;
    k = abs_wave_num;  E = Energy_Spectrum_Value(k);
    return (k == 0.0)  ?  0.0 : sqrt(E/N) * exp(I*theta2) * sin(phi);
}

 
/*!
 * Subroutine: Create_Homogeneous_Turbulence_Velocity_Field
 * 
 * Creates and assigns a turbulence velocity field
 *
 *
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Create_Homogeneous_Turbulence_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
					     const Grid3D_Input_Parameters &IPs,
					     int &batch_flag,
                                             Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field ) {

  int errorflag;
  errorflag = Open_Turbulence_Spectrum_File( );
  bool rescale = false;
  
  Nx = IPs.NCells_Turbulence_Idir*Initial_Mesh.NBlk_Idir;
  Ny = IPs.NCells_Turbulence_Jdir*Initial_Mesh.NBlk_Jdir;
  Nz = IPs.NCells_Turbulence_Kdir*Initial_Mesh.NBlk_Kdir; 

  double        *u, *v, *w;        // Arrays to store the velocity fluctuations in physical space
  fftw_complex  *uu, *vv, *ww;     // Arrays to store the velocity fluctuations in spectral space
  fftw_plan      physical;

  int index;
  int nz = Nz/2+1;

    // Allocation of arrays used in the transforms
    u = (double *) malloc(Nx*Ny*Nz * sizeof(double));
    v = (double *) malloc(Nx*Ny*Nz * sizeof(double));
    w = (double *) malloc(Nx*Ny*Nz * sizeof(double));
    uu = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
    vv = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
    ww = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));


    int seed = 1; // = time(NULL);   // assigns the current time to the seed
    srand48(seed);             // changes the seed for drand48()
    double k1, k2, k3, abs_k;  // Wave numbers
    
    double k0x, k0y, k0z, abs_k0, kcx, kcy, kcz, abs_kc;
    k0x = TWO*PI/L1;    // smallest wavenumber on grid in x-direction
    k0y = TWO*PI/L2;    // smallest wavenumber on grid in y-direction
    k0z = TWO*PI/L3;    // smallest wavenumber on grid in z-direction
    abs_k0 = min(min(k0x,k0y),k0z);
    kcx = Nx/2*k0x;     // largest wavenumber on grid in x-direction
    kcy = Ny/2*k0y;     // largest wavenumber on grid in y-direction
    kcz = Nz/2*k0z;     // largest wavenumber on grid in z-direction
    abs_kc = sqrt( sqr(kcx) + sqr(kcy) + sqr(kcz) );

    if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) {
        cout << "\n\n";
        cout << " ==========================================================================\n"; 
        cout << "         Generating Homogeneous Isotropic Turbulent Velocity Field \n";
        cout << " ==========================================================================\n";
        cout << endl;
        cout << "   Turbulence parameters in case of continuous spectrum: " << endl;
        cout << "    -->  TKE         = " << TKE << endl;
        cout << "    -->  u_RMS       = " << up << endl;
        cout << "    -->  L           = " << Lp << endl;
        cout << "    -->  L11 ~ 0.4 L = " << Lp*0.4 << endl;
        cout << "    -->  eps         = " << eps << endl;
        cout << "    -->  eta         = " << eta << endl;
        cout << "    -->  ReL         = " << ReL << endl;
        cout << "    -->  Rlambda     = " << Rlambda << endl;
        cout << "    -->  lambda_g    = " << lambda_g << endl;
    }
    
        
    K_container *K = new K_container [Nx*Ny*nz];  // Container of information for a wavenumber
    K_BinaryTree Ktree;                           // A binary tree that will order and store K_containers
    for (int i=0; i<Nx; ++i) {
        k1 = k_1(i);

        for (int j=0; j<Ny; ++j) {
            k2 = k_2(j);

            for (int l=0; l<nz; ++l) {
                k3 = k_3(l);

                abs_k = sqrt(k1*k1 + k2*k2 + k3*k3); // Wave number magnitude
                
                index = l + nz*(j+Ny*i);
                K[index].k = abs_k;
                K[index].index = index;
                Ktree.InsertNode(K[index]);
            }   
        }
    }
    delete[] K;

    int nK = Ktree.countNodes();    // The total number of unique values of abs_k
    
    K = Ktree.asArray();            // Make an array of the tree nodes in order.
    
    
    /* -------------------- gnuplot model -------------------- */
    #ifdef _GNUPLOT
    int npts = int((k_eta-abs_k0)/abs_k0);
        double kmin = abs_k0/10;
        double kmax = k_eta/4.0;
        dpoint *dp = new dpoint[npts];
        for (int i=0; i<npts; i++) {
            dp[i].x= kmin + i*(kmax-kmin)/(npts-1);
            dp[i].y= Energy_Spectrum_Value(dp[i].x);
        }
        Gnuplot_Control h1;
        h1.gnuplot_init();
        h1.gnuplot_setstyle("lines") ;
        h1.gnuplot_cmd("set logscale xy");
        h1.gnuplot_cmd("set grid");
        h1.gnuplot_set_xlabel("k");
        h1.gnuplot_set_ylabel("E(k)");
        h1.gnuplot_set_title("Turbulent kinetic energy spectrum");
        h1.gnuplot_plot1d_var2(dp,npts,"model");
        sleep(2);
    #endif

    
    double theta1, theta2, phi;         // random phases for Rogallo spectrum
    complex<double> aa, bb;             // complex energy components for Rogallo spectrum
    double deno;                        // denominator in Rogallo spectrum function
    
    double *EE;
    EE = new double [Nx*Ny*nz];   // Energy of one wavenumber vector
    
    
    for(int ii=0; ii<nK; ii++) {            // For every abs_k
        for(int jj=0; jj<K[ii].N; jj++) {       // For all grid points with this value of abs_k
            index = K[ii].indexes[jj];
            int i = index/(Ny*nz);              // convert index to (i,j,l) coordinates
            int j = (index - Ny*nz*i)/nz;
            int l = index - nz*(j+Ny*i);

            k1 = k_1(i);
            k2 = k_2(j);
            k3 = k_3(l);
            abs_k = K[ii].k;
            
            /* ----------- Calculate velocities with Rogallo's procedure ---------- */
            theta1 = TWO*PI*random_double();   // Random number (0, 2*PI)
            theta2 = TWO*PI*random_double();   // Random number (0, 2*PI)
            phi    = TWO*PI*random_double();   // Random number (0, 2*PI)
            
            if ( theta1 == theta2  && theta2 == phi ) {
                cerr << "\n theta1, theta2 and phi are all equal.";
            }
            
            //aa = alpha_discrete(abs_k, K[ii].N, theta1, phi);
            //bb =  beta_discrete(abs_k, K[ii].N, theta2, phi);
            
            aa = alpha(abs_k, theta1, phi);
            bb =  beta(abs_k, theta2, phi);

            deno = abs_k * sqrt(k1*k1 + k2*k2) + PICO;
            
            uu[index][0] = real( (aa*abs_k*k2 + bb*k1*k3)/deno ) ;
            uu[index][1] = imag( (aa*abs_k*k2 + bb*k1*k3)/deno ) ;
            vv[index][0] = real( (bb*k2*k3 - aa*abs_k*k1)/deno ) ;
            vv[index][1] = imag( (bb*k2*k3 - aa*abs_k*k1)/deno ) ;
            ww[index][0] = real( -( bb*(k1*k1 + k2*k2) )/deno  ) ;
            ww[index][1] = imag( -( bb*(k1*k1 + k2*k2) )/deno  ) ;
            
            /* ------------ Some exceptions to make FFT real ------------ */
            if ( l==0  ||  l==Nz/2) {
                if ( j>Ny/2  ||  ( i>Nx/2  &&  (j==0 || j==Ny/2) ) ) {
                    int iconj = (i==0  ?  0 : Nx-i);
                    int jconj = (j==0  ?  0 : Ny-j);
                    int index_conj = l + nz*(jconj+Ny*iconj);
                    
                    uu[index][0] =  uu[index_conj][0]; // complex conjugates
                    uu[index][1] = -uu[index_conj][1];
                    vv[index][0] =  vv[index_conj][0];
                    vv[index][1] = -vv[index_conj][1];
                    ww[index][0] =  ww[index_conj][0];
                    ww[index][1] = -ww[index_conj][1];
                    
                } else if ( (i==0 || i==Nx/2)  &&  (j==0 || j==Ny/2) ) { // real values at 8 corners
                    uu[index][0] = sqrt( sqr(uu[index][0]) + sqr(uu[index][1]) );
                    uu[index][1] = 0.0;    
                    vv[index][0] = sqrt( sqr(vv[index][0]) + sqr(vv[index][1]) );
                    vv[index][1] = 0.0;
                    ww[index][0] = sqrt( sqr(ww[index][0]) + sqr(ww[index][1]) );
                    ww[index][1] = 0.0; 
                }
            }

            /* ---------- Calculate energy corresponding with this index -------- */
            // Velocity correlation tensor for this wavenumber
            double Rxx = (sqr(uu[index][0]) + sqr(uu[index][1]));
            double Ryy = (sqr(vv[index][0]) + sqr(vv[index][1]));
            double Rzz = (sqr(ww[index][0]) + sqr(ww[index][1]));
            
            // Energy = HALF * trace(velocity correlation tensor Rij)
            EE[index] = HALF * (Rxx + Ryy + Rzz) ;
            
            /* ------- Calculate energy corresponding with this wavenumber ------ */
            K[ii].Ek += TWO * EE[index] * (TWO*PI*sqr(abs_k)) / K[ii].N;
            
        }
        K[ii].Emodel = Energy_Spectrum_Value(K[ii].k);

        if (rescale) {  // rescale is false now
            if (K[ii].Ek != ZERO)
                K[ii].velocity_scaling_factor = sqrt(K[ii].Emodel/(K[ii].Ek));
            else 
                K[ii].velocity_scaling_factor = ONE;
        }

    }
    

    /* ---------------- gnuplot Rogallo spectrum -------------- */
    double *k  = new double [nK];
    double *Ek = new double [nK];
    double *Emodel = new double [nK];

    for (int ii=1; ii <nK; ii++) {
        k[ii-1]  = K[ii].k;   // Don't wish to store K[0].k which is 0
        Ek[ii-1] = K[ii].Ek; 
        Emodel[ii-1] = K[ii].Emodel;
    }
    #ifdef _GNUPLOT
        h1.gnuplot_setstyle("points");
        h1.gnuplot_plot1d_var2(k,Ek,nK-1,"Rogallo");
    #endif
    
    double TKE_model   = CubicSplinesIntegration(k,Emodel,nK-1);
    double TKE_rogallo = CubicSplinesIntegration(k,Ek,nK-1);

    if (!batch_flag) {       
        cout << endl;
        cout << "   Total Kinetic Energy calculated in spectral space: " << endl;
        cout << "    -->  TKE of the model            = " << TKE_model << endl;
        cout << "    -->  TKE with Rogallo            = " << TKE_rogallo << endl;
    }
    
    
    /* ------- rescale the rogallo spectrum to match the model spectrum ------- */
    if (rescale) {   // rescale is false now
        for (int ii=0; ii<nK; ii++) {
            K[ii].Ek = 0;
            double abs_k = K[ii].k;
            for (int jj=0; jj<K[ii].N; jj++) {
                int index = K[ii].indexes[jj];
                double velocity_scaling_factor = K[ii].velocity_scaling_factor;
                uu[index][0] *= velocity_scaling_factor;
                uu[index][1] *= velocity_scaling_factor;
                vv[index][0] *= velocity_scaling_factor;
                vv[index][1] *= velocity_scaling_factor;
                ww[index][0] *= velocity_scaling_factor;
                ww[index][1] *= velocity_scaling_factor;
                
                double Rxx = (sqr(uu[index][0]) + sqr(uu[index][1]));
                double Ryy = (sqr(vv[index][0]) + sqr(vv[index][1]));
                double Rzz = (sqr(ww[index][0]) + sqr(ww[index][1]));
                EE[index] = HALF * (Rxx + Ryy + Rzz);
                
                K[ii].Ek += TWO * EE[index] * (TWO*PI*sqr(abs_k)) / K[ii].N;
            }
            if (ii>0) {
                Ek[ii-1] = K[ii].Ek; // Don't wish to store K[0].k which is 0
            } 
        } 
        
        TKE_rogallo = CubicSplinesIntegration(k,Ek,nK-1);

        if (!batch_flag)
            cout << "    -->  TKE with Rogallo rescaled   = " << TKE_rogallo << endl;
    }
    
    /* ------------------ write gnuplot spectrum files --------------- */
    for (int ii=0; ii <nK; ii++) {
        if (ii>0) {
            k[ii-1]  = K[ii].k;   // Don't wish to store K[0].k which is 0
            Ek[ii-1] = K[ii].Ek; 
            Emodel[ii-1] = K[ii].Emodel;
        }    
        Output_Turbulence_Spectrum_to_File(K[ii].k,K[ii].Emodel,K[ii].Ek);
    }
    
    /* ------------------ gnuplot the rescaled spectrum -------------- */
    #ifdef _GNUPLOT
        if (rescale)    // rescale is false now
            h1.gnuplot_plot1d_var2(k,Ek,nK,"Rogallo - rescaled");
    #endif

    /* ---------- Total kinetic energy for continuous spectrum ---------- */


 
    typedef double (RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::*Energy_Spectrum_Type) (const double &abs_wave_num) const;
    _Member_Function_Wrapper_<RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>,Energy_Spectrum_Type, double> Energy_function(this, &RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::Energy_Spectrum_Value);
    double dummy;
    double TKE_entire_range = AdaptiveGaussianQuadrature(Energy_function, 0.0, k_eta/2.0, dummy,5);
    
    if (!batch_flag) {       
        cout << "    -->  TKE of the entire range     = " << TKE_entire_range << endl;
        cout << "    -->  TKE percentage on the grid  = " << int(TKE_rogallo/TKE*HUNDRED) << " %" << endl;
	cout << "    -->  u_RMS with Rogallo          = " << sqrt(TWO/THREE*TKE_rogallo) << endl;
    }

    
    
    /* ----------------- calculate L11 ---------------------- */
    
    
    double *Ek_over_k = new double [nK];
    for (int ii=1; ii <nK; ii++) {
            Ek_over_k[ii-1] = K[ii].Ek/K[ii].k; 
    }
    double L11_spectral = THREE*PI/FOUR * CubicSplinesIntegration(k,Ek_over_k,nK-1) / TKE_rogallo;
    
    if (!batch_flag)      
        cout << "    -->  L11                         = " << L11_spectral << endl;

    /* ---------------------------------------------------------------------------- *
     *          Perform Fourier transform from spectral to physical space           *
     * ---------------------------------------------------------------------------- */

    physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, uu, u, FFTW_ESTIMATE);
    fftw_execute(physical); 
    fftw_destroy_plan(physical);
 
    physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, vv, v, FFTW_ESTIMATE);
    fftw_execute(physical); 
    fftw_destroy_plan(physical);
  
    physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, ww, w, FFTW_ESTIMATE);
    fftw_execute(physical); 
    fftw_destroy_plan(physical);
  
    /* --------------- Total kinetic energy in physical space --------------- */
    double TKE_physical = 0;
    double u2 = 0 , v2 = 0 , w2 = 0;
    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny; j++) {
            for (int l=0; l<Nz; l++) {
                index = l + Nz*(j+Ny*i);
                u[index] *= sqrt(pow(TWO*PI,3.0)/(L1*L2*L3));
                v[index] *= sqrt(pow(TWO*PI,3.0)/(L1*L2*L3));
                w[index] *= sqrt(pow(TWO*PI,3.0)/(L1*L2*L3));

                u2 += ONE/(Nx*Ny*Nz) * sqr(u[index]);   // average of u^2
                v2 += ONE/(Nx*Ny*Nz) * sqr(v[index]);   // average of v^2
                w2 += ONE/(Nx*Ny*Nz) * sqr(w[index]);   // average of w^2
            } /* endfor */
        } /* endfor */
    } /* endfor */
    TKE_physical = HALF*(u2+v2+w2);
    double uRMS = sqrt(TWO/THREE*TKE_physical);
    
    if (!batch_flag) {
        cout << endl;
        cout << "   Total Kinetic Energy calculated in physical space: " << endl;
        cout << "    -->  discrete TKE  = " << TKE_physical << endl;
        cout << "    -->  u_RMS         = " << uRMS << endl;
    }
    
    
    
    /* ---------- Longitudinal velocity correlation function fr ---------- */
    double *Rr = new double [Nx];
    double *fr = new double [Nx];
    double *r = new double [Nx];
    for (int ri=0; ri<Nx; ri++){ // loop for variation in r
        Rr[ri] = 0;
    }
    for (int l=0; l<Nz; l++) {
        for (int j=0; j<Ny; j++) {
            for (int i=0; i<Nx; i++) {
                for (int ii=i; ii<Nx; ii++) {
                    int indexr = l + Nz*( j + Ny*ii );
                    int index  = l + Nz*( j + Ny*i );
                    int ri = ii-i;
                    Rr[ri] += ONE/(Nx*Ny*Nz) * (u[indexr]*u[index]); 
                }
            }
        }
    }
    double dr = L1/(Nx-1.0);
    for (int ri = 0; ri<Nx; ri++) {
        r[ri] = ri*dr;
        fr[ri] = Rr[ri]/u2;
    }
    
    /* ---------------- L11 = Integral of fr ---------------- */
    double L11_physical = CubicSplinesIntegration(r,fr,Nx);
    if (!batch_flag)
        cout << "    -->  L11           = " << L11_physical << endl;

    
#ifdef _GNUPLOT
    Gnuplot_Control h2;
    h2.gnuplot_init();
    h2.gnuplot_setstyle("lines") ;
    h2.gnuplot_cmd("set grid");
    h2.gnuplot_set_xlabel("r");
    h2.gnuplot_set_ylabel("f(r)");
    h2.gnuplot_set_title("longitudinal velocity autocorrelation function");
    h2.gnuplot_plot1d_var2(r,fr,Nx,"");
//    sleep(2);    
//    h2.gnuplot_close();
#endif
    
    
//    /* ---------------------------------------------------------------------------- *
//     *          Perform Fourier transform from physical to spectral space           *
//     * ---------------------------------------------------------------------------- */
//    fftw_plan      spectral;
//
//    spectral = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, u, uu, FFTW_ESTIMATE);
//    fftw_execute(spectral); 
//    fftw_destroy_plan(spectral);
//    
//    spectral = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, v, vv, FFTW_ESTIMATE);
//    fftw_execute(physical); 
//    fftw_destroy_plan(spectral);
//    
//    spectral = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, w, ww, FFTW_ESTIMATE);
//    fftw_execute(spectral); 
//    fftw_destroy_plan(spectral);
//    
//    
//    
//    for (int ii=0; ii<nK; ii++) {
//        K[ii].Ek = 0;
//        double abs_k = K[ii].k;
//        for (int jj=0; jj<K[ii].N; jj++) {
//            int index = K[ii].indexes[jj];
//            
//            uu[index][0] *= ONE/(Nx*Ny*Nz)   / sqrt(pow(TWO*PI,3)/(L1*L2*L3));
//            uu[index][1] *= ONE/(Nx*Ny*Nz)   / sqrt(pow(TWO*PI,3)/(L1*L2*L3));
//            vv[index][0] *= ONE/(Nx*Ny*Nz)   / sqrt(pow(TWO*PI,3)/(L1*L2*L3));
//            vv[index][1] *= ONE/(Nx*Ny*Nz)   / sqrt(pow(TWO*PI,3)/(L1*L2*L3));
//            ww[index][0] *= ONE/(Nx*Ny*Nz)   / sqrt(pow(TWO*PI,3)/(L1*L2*L3));
//            ww[index][1] *= ONE/(Nx*Ny*Nz)   / sqrt(pow(TWO*PI,3)/(L1*L2*L3));
//
//            double Rxx = (sqr(uu[index][0]) + sqr(uu[index][1]));
//            double Ryy = (sqr(vv[index][0]) + sqr(vv[index][1]));
//            double Rzz = (sqr(ww[index][0]) + sqr(ww[index][1]));
//            EE[index] = HALF * (Rxx + Ryy + Rzz);
//            
//            K[ii].Ek += TWO * EE[index] * (TWO*PI*sqr(abs_k)) / K[ii].N;
//        }
//        if (ii>0) {
//            k[ii-1] = K[ii].k;
//            Ek[ii-1] = K[ii].Ek; // Don't wish to store K[0].k which is 0
//        } 
//    } 
//    
//    /* ------------------ gnuplot the rescaled spectrum -------------- */
//#ifdef _GNUPLOT
//        h1.gnuplot_plot1d_var2(k,Ek,nK-1,"from physical");
//#endif
//    
//    cout << " TKE from physical back to spectral = " << CubicSplinesIntegration(k,Ek,nK-1);
//    
  

  /* ------------------------------------ Assign turbulent velocity field ------------------------------------ */
  int nBlk, ix, iy, iz;
  int INl, INu, JNl, JNu, KNl, KNu;
  for (int kBlk = 0; kBlk <= Initial_Mesh.NBlk_Kdir-1; ++kBlk) {
     for (int jBlk = 0; jBlk <= Initial_Mesh.NBlk_Jdir-1; ++jBlk) {
        for (int iBlk = 0; iBlk <= Initial_Mesh.NBlk_Idir-1; ++iBlk) {
            nBlk = iBlk +
                   jBlk*Initial_Mesh.NBlk_Idir +
                   kBlk*Initial_Mesh.NBlk_Idir*Initial_Mesh.NBlk_Jdir;
            for (int i = Initial_Mesh.Grid_Blks[nBlk].ICl; i <= Initial_Mesh.Grid_Blks[nBlk].ICu; ++i) {
               for (int j = Initial_Mesh.Grid_Blks[nBlk].JCl; j <= Initial_Mesh.Grid_Blks[nBlk].JCu; ++j) {
                  for (int k = Initial_Mesh.Grid_Blks[nBlk].KCl; k <= Initial_Mesh.Grid_Blks[nBlk].KCu; ++k) {
		     ix = iBlk*IPs.NCells_Turbulence_Idir+(i-Initial_Mesh.Grid_Blks[nBlk].ICl);
		     iy = jBlk*IPs.NCells_Turbulence_Jdir+(j-Initial_Mesh.Grid_Blks[nBlk].JCl);
		     iz = kBlk*IPs.NCells_Turbulence_Kdir+(k-Initial_Mesh.Grid_Blks[nBlk].KCl);
	             index = iz + 
                             iy*Nz + 
                             ix*Ny*Nz;
                     Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k].x = u[index];
                     Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k].y = v[index];
                     Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k].z = w[index];
		     Initial_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k] = Initial_Mesh.Grid_Blks[nBlk].Cell[i][j][k].Xc;	     
                     
		  } /* endfor */
	       } /* endfor */
	    } /* endfor */  
	    INl = Initial_Mesh.Grid_Blks[nBlk].INl; 
	    INu = Initial_Mesh.Grid_Blks[nBlk].INu; 
	    JNl = Initial_Mesh.Grid_Blks[nBlk].JNl; 
	    JNu = Initial_Mesh.Grid_Blks[nBlk].JNu; 
	    KNl = Initial_Mesh.Grid_Blks[nBlk].KNl; 
	    KNu = Initial_Mesh.Grid_Blks[nBlk].KNu;
	    Initial_Velocity_Field.Vel_Blks[nBlk].Node_INl_JNl_KNl = Initial_Mesh.Grid_Blks[nBlk].Node[INl][JNl][KNl].X; 
	    Initial_Velocity_Field.Vel_Blks[nBlk].Node_INu_JNu_KNu = Initial_Mesh.Grid_Blks[nBlk].Node[INu][JNu][KNu].X;

	} /* endfor */
     } /* endfor */
  } /* endfor */
  
      /* --------------- Deallocations --------------- */
    for (int ii=0; ii <nK; ii++) {
        delete[] K[ii].indexes;
    }
    delete[] K;
    delete[] EE;
    delete[] k;
    delete[] Ek;

    delete[] Emodel;
    fftw_free(u);
    fftw_free(v);
    fftw_free(w);
    fftw_free(uu);
    fftw_free(vv);
    fftw_free(ww);
    
    #ifdef _GNUPLOT
        delete[] dp;
       // sleep(2);
       // h1.gnuplot_close();
    #endif
    
    errorflag = Close_Turbulence_Spectrum_File( );      
    return (0);
}


// Write_Turbulent_Velocity_Field
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Write_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
			       const Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
                               const Grid3D_Input_Parameters &IPs) const {

    ofstream out_file;
    out_file.open("Initial_Turbulence_Fluctuations.dat", ios::out);
    if (out_file.fail()) {
      cerr<<"\n Error opening file: Initial_Turbulence_Fluctuations.dat to write" << endl;
      exit(1);
    } /* endif */
  
    out_file.setf(ios::scientific);

    for (int nBlk = 0; nBlk <= Initial_Mesh.NBlk-1; ++nBlk ) {
       for (int i = Initial_Mesh.Grid_Blks[nBlk].ICl; i <= Initial_Mesh.Grid_Blks[nBlk].ICu; ++i) {
          for (int j = Initial_Mesh.Grid_Blks[nBlk].JCl; j <= Initial_Mesh.Grid_Blks[nBlk].JCu; ++j) {
             for (int k = Initial_Mesh.Grid_Blks[nBlk].KCl; k <= Initial_Mesh.Grid_Blks[nBlk].KCu; ++k) {
	      out_file << setprecision(10)
		       << Initial_Mesh.Grid_Blks[nBlk].Cell[i][j][k].Xc << " " 
                       << Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k];
	     } /* endfor */
	  } /* endfor */
       } /* endfor */
       out_file.setf(ios::scientific);
       out_file << endl;
    } /* endfor */

    out_file.close();

}

// Read_Turbulent_Velocity_Field
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Read_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
			      Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
                              Grid3D_Input_Parameters &IPs) {

   ifstream InFile;

   // Open data file for reading
   InFile.open("Initial_Turbulence_Fluctuations.dat", ios::in); 
   // Check to see if successful
   if(InFile.fail()){ 
     cerr<<"\n Error opening file: Initial_Turbulence.dat to read" <<endl;
     exit(1); 
   } 

   bool interpolate_flag;
   double xx, yy, zz, dx, dy, dz, dd, uprime_x, uprime_y, uprime_z;

   InFile.setf(ios::skipws);
 
   InFile.unsetf(ios::skipws);
   InFile.close();

}

/*!
 * Open_Turbulence_Spectrum_File
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Open_Turbulence_Spectrum_File( ) {
    
    int i;
    char prefix[256], extension[256], 
    turbulence_spectrum_file_name[256], gnuplot_file_name[256];
    char *turbulence_spectrum_file_name_ptr, *gnuplot_file_name_ptr;
    ofstream gnuplot_file;
    
    /* Determine the name of the turbulence spectrum file. */
    
    i = 0;
    while (1) {
        if (File_Name[i] == ' ' ||
            File_Name[i] == '.') break;
        prefix[i] = File_Name[i];
        i = i + 1;
        if (i > strlen(File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_energy_spectrum");
    
    strcpy(extension, ".dat");
    strcpy(turbulence_spectrum_file_name, prefix);
    strcat(turbulence_spectrum_file_name, extension);
    
    turbulence_spectrum_file_name_ptr = turbulence_spectrum_file_name;
    
    /* Open the turbulence spectrum file. */
    
    Turbulence_Spectrum_File.open(turbulence_spectrum_file_name_ptr, ios::out);
    if (Turbulence_Spectrum_File.bad()) return (1);
    
    /* Write the appropriate GNUPLOT command file for 
     plotting turbulence progress file information. */
    
    strcpy(extension, ".gplt");
    strcpy(gnuplot_file_name, prefix);
    strcat(gnuplot_file_name, extension);
    
    gnuplot_file_name_ptr = gnuplot_file_name;
    
    gnuplot_file.open(gnuplot_file_name_ptr, ios::out);
    if (gnuplot_file.fail()) return(1);
    
    gnuplot_file << "set title \"Turbulence energy spectrum "<< spectrum_name <<"\"\n"
    << "set xlabel \"wavenumber k \"\n"
    << "set ylabel \"EE\"\n" 
    << "set logscale xy\n"
    << "plot \"" << turbulence_spectrum_file_name_ptr << "\" using 1:2 \\\n"
    << "     title \"" << "model" << "\" with lines , \\\n"
    << "\"" << turbulence_spectrum_file_name_ptr << "\"  using 1:3 \\\n"
    << "     title \"" << "constructed" << "\" with lines\n"
    << "pause -1  \"Hit return to continue\"\n";
    
    gnuplot_file.close();
    
    /* Preparation of progress file complete.
     Return zero value. */    
    return(0);
}


/*!
 * Output_Turbulence_Spectrum_to_File
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Output_Turbulence_Spectrum_to_File(const double &k,
                                   const double &model_Energy,
                                   const double &constructed_Energy) {
    
    Turbulence_Spectrum_File << setprecision(6);
    Turbulence_Spectrum_File.setf(ios::scientific);
    
    Turbulence_Spectrum_File << k   << " " << model_Energy
    << " " << constructed_Energy 
    << "\n";
    
    Turbulence_Spectrum_File.unsetf(ios::scientific);
    Turbulence_Spectrum_File.flush();
    
    return(0);
}

/*!
 * Close_Turbulence_Spectrum_File
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Close_Turbulence_Spectrum_File( ) {
    Turbulence_Spectrum_File.close();
    return(0);
}


// Assign_Homogeneous_Turbulence_Velocity_Field
template<typename HEXA_BLOCK>
void Assign_Homogeneous_Turbulence_Velocity_Field(HEXA_BLOCK *Solution_Block,
                                                  const AdaptiveBlock3D_List &LocalSolnBlockList,
                                                  const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field) {

   /* Assign initial turbulent velocity field to each solution block. */

   for (int nBlk = 0 ; nBlk <= LocalSolnBlockList.Nblk-1 ; nBlk++) {
      if (LocalSolnBlockList.Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
	 if (Velocity_Field.Vel_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum].Allocated) {
  	    Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Block[nBlk],
                                                         Velocity_Field.Vel_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum]);
	 } /* endif */
      } /* endif */
   }  /* endfor */

}

template<typename HEXA_BLOCK>
void Assign_Homogeneous_Turbulence_Velocity_Field(HEXA_BLOCK &Solution_Block,
                                                  const Turbulent_Velocity_Field_Block &Velocity_Field) {

  /* Assign initial turbulent velocity field to solution block. */

  for (int i = Solution_Block.ICl ; i <= Solution_Block.ICu ; i++) {
     for (int j = Solution_Block.JCl ; j <= Solution_Block.JCu ; j++) {
        for (int k = Solution_Block.KCl ; k <= Solution_Block.KCu ; k++) {
           Solution_Block.W[i][j][k].v += Velocity_Field.Velocity[i][j][k];
 	   Solution_Block.U[i][j][k] = Solution_Block.W[i][j][k].U();
	} /* endfor */
     } /* endfor */
  } /* endfor */

}


/********************************************************
 *          Open_Turbulence_Progress_File               *
 ********************************************************/
inline int Open_Turbulence_Progress_File(ofstream &Turbulence_Progress_File,
		    		         char *File_Name,
				         const int Append_to_File) {

    int i;
    char prefix[256], extension[256], 
         turbulence_progress_file_name[256], gnuplot_file_name[256];
    char *turbulence_progress_file_name_ptr, *gnuplot_file_name_ptr;
    ofstream gnuplot_file;

    /* Determine the name of the turbulence progress file. */

    i = 0;
    while (1) {
       if (File_Name[i] == ' ' ||
           File_Name[i] == '.') break;
       prefix[i] = File_Name[i];
       i = i + 1;
       if (i > strlen(File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_turbulence_statistics");

    strcpy(extension, ".dat");
    strcpy(turbulence_progress_file_name, prefix);
    strcat(turbulence_progress_file_name, extension);

    turbulence_progress_file_name_ptr = turbulence_progress_file_name;

    /* Open the turbulence progress file. */

    if (Append_to_File) {
       Turbulence_Progress_File.open(turbulence_progress_file_name_ptr, ios::out|ios::app);
    } else {
       Turbulence_Progress_File.open(turbulence_progress_file_name_ptr, ios::out);
    } /* endif */
    if (Turbulence_Progress_File.bad()) return (1);

    /* Write the appropriate GNUPLOT command file for 
       plotting turbulence progress file information. */

    strcpy(extension, ".gplt");
    strcpy(gnuplot_file_name, prefix);
    strcat(gnuplot_file_name, extension);

    gnuplot_file_name_ptr = gnuplot_file_name;

    gnuplot_file.open(gnuplot_file_name_ptr, ios::out);
    if (gnuplot_file.fail()) return(1);

    gnuplot_file << "set title \"Turbulence parameters progress \"\n"
                 << "set xlabel \"Time \"\n"
                 << "set ylabel \"TKE/u_rms/enstrophy/Taylor/St\"\n" 
                 << "set logscale xy\n"
                 << "plot \"" << turbulence_progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%lf%*lf%*lf%*lf%*lf%*lf\" \\\n"
                 << "     title \"TKE\" with lines, \\\n"
                 << "\"" << turbulence_progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%lf%*lf%*lf%*lf%*lf\" \\\n"
                 << "     title \"u_rms\" with lines, \\\n"
                 << "\"" << turbulence_progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%*lf%lf%*lf%*lf%*lf\" \\\n"
                 << "     title \"enstrophy\" with lines, \\\n"
                 << "\"" << turbulence_progress_file_name_ptr << "\""
		 << " using 1:2 \"%lf%*lf%*lf%*lf%*lf%*lf%lf%*lf%*lf\" \\\n"
	         << "     title \"Taylor_scale\" with lines, \\\n"
                 << "\"" << turbulence_progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%lf\" \\\n"
                 << "     title \"St\" with lines\n"
                 << "pause -1  \"Hit return to continue\"\n";

    gnuplot_file.close();

    /* Preparation of progress file complete.
       Return zero value. */

    return(0);

}

/********************************************************
 *          Close_Turbulence_Progress_File              *
 ********************************************************/
inline int Close_Turbulence_Progress_File(ofstream &Turbulence_Progress_File) {

    Turbulence_Progress_File.close();

    return(0);

}

/********************************************************
 *          Output_Turbulence_Progress_to_File          *
 ********************************************************/
inline int Output_Turbulence_Progress_to_File(ostream &Turbulence_Progress_File,
		 			      const int Number_of_Time_Steps,
					      const double &Time,
					      const CPUTime &CPU_Time,
					      const double &Total_Energy,
                                              const double &u_rms,
                                              const double &Total_Enstrophy,
					      const double &Taylor_scale,
					      const double &viscosity,
                                              const double &turbulent_burning_rate) {

    Turbulence_Progress_File << setprecision(6);
    Turbulence_Progress_File << Time
                             << " " << Number_of_Time_Steps
                             << " " << CPU_Time.min();
    Turbulence_Progress_File.setf(ios::scientific);
    Turbulence_Progress_File << " " << Total_Energy
                             << " " << u_rms
                             << " " << Total_Enstrophy
                             << " " << Taylor_scale
			     << " " << viscosity
			     << " " << turbulent_burning_rate
                             << "\n";                       
    Turbulence_Progress_File.unsetf(ios::scientific);
    Turbulence_Progress_File.flush();

    return(0);

}



template<typename HEXA_BLOCK>
int Longitudinal_Correlation(const Octree_DataStructure &OcTree,
                             const AdaptiveBlock3D_ResourceList &Global_Soln_Block_List,
			     const AdaptiveBlock3D_List &LocalSolnBlockList,
			     HEXA_BLOCK *Solution_Block,
			     const Grid3D_Input_Parameters &IPs,
			     const double &u_ave, 
                             const double &sqr_u);



#endif // _TURBULENT_VELOCITY_FIELD_INCLUDED 
