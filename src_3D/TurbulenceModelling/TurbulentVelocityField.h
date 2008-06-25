
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
    int                 gblknum; // global turbulent velocity field number
    int             NCi,ICl,ICu; // i-direction turbulent velocity field cell counters
    int             NCj,JCl,JCu; // j-direction turbulent velocity field cell counters
    int             NCk,KCl,KCu; // k-direction turbulent velocity field cell counters
    int                  Nghost; // number of ghost cells
    Vector3D        ***Velocity; // array of turbulent velocity field vectors
    Vector3D        ***Position; // array of turbulent velocity field position vectors
    Vector3D   Node_INl_JNl_KNl, // diagonally opposite corners of the block 
               Node_INu_JNu_KNu;

    Vector3D            ***dVdx; // array of turbulent velocity field derivatives
    Vector3D            ***dVdy; // array of turbulent velocity field derivatives
    Vector3D            ***dVdz; // array of turbulent velocity field derivatives

    int  Allocated; // Indicates whether or not the turbulent velocity field data have been allocated.
    int _Allocated; // Indicates whether or not the turbulent velocity field gradients have been allocated.
    
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
       _Allocated = TURBULENT_VELOCITY_FIELD_DATA_NOT_USED;
       dVdx = NULL;
       dVdy = NULL;
       dVdz = NULL;      
    }

    Turbulent_Velocity_Field_Block(const int Ni, 
                                   const int Nj, 
                                   const int Nk,
                                   const int Ng) {
      allocate(Ni, Nj, Nk, Ng);
    }

    //! Destructor
    ~Turbulent_Velocity_Field_Block(void) {
       deallocate_gradients();
       deallocate();        
    }

    //! Allocate memory for velocity field data. 
    void allocate(const int Ni, 
                  const int Nj, 
                  const int Nk,
                  const int Ng);
    
    //! Deallocate memory for velocity field data.
    void deallocate(void);

    //! Allocate memory for velocity field gradients. 
    void allocate_gradients(void);

    //! Deallocate memory for velocity field gradients.
    void deallocate_gradients(void);

    //! Reconstruct velocity field data. 
    void LeastSquares_Reconstruction(const int &i,
				     const int &j,
				     const int &k,
				     Vector3D &dVdx,
				     Vector3D &dVdy,
				     Vector3D &dVdz);

    //! Least Squares Reconstruction
    void LeastSquares_Reconstruction(const int &i,
				     const int &j,
				     const int &k);

    //! Least Squares Reconstruction
    void LeastSquares_Reconstruction(void);

    //! Copy operator
    void Copy(Turbulent_Velocity_Field_Block &Block2);

    //! Assign indices
    void assign_indices(int &n_pts,
			const int &i,
			const int &j,
			const int &k,
			int i_index[],
			int j_index[],
			int k_index[]);

    //! Broadcast
    void Broadcast(void);

#ifdef _MPI_VERSION   
    //! Broadcast from a given processor
    void Broadcast(MPI::Intracomm &Communicator, 
                   const int Source_CPU);
#endif
    
    //! Output operator
    friend ostream &operator << (ostream &out_file, 
                                 const Turbulent_Velocity_Field_Block &V);
    //! Input operator
    friend istream &operator >> (istream &in_file, 
                                 Turbulent_Velocity_Field_Block &V);
  
    /* Other useful member functions. */

  private:
    //copy and assignment are not permitted
    Turbulent_Velocity_Field_Block(const Turbulent_Velocity_Field_Block &V);
    Turbulent_Velocity_Field_Block &operator =(const Turbulent_Velocity_Field_Block &V);
};

/*************************************************************************
 * Turbulent_Velocity_Field_Block -- Input-output operators.             *
 *************************************************************************/
inline ostream &operator << (ostream &out_file,
                             const Turbulent_Velocity_Field_Block &V){  
    
    out_file << V.NCi << " " << V.NCj << " " << V.NCk << "\n";
    out_file << V.Nghost << " " << V.gblknum << "\n";
    out_file << V.Node_INl_JNl_KNl << " " << V.Node_INu_JNu_KNu << "\n";
    
    if (V.NCi == 0 || V.NCj == 0 || V.NCk == 0 || V.Nghost == 0)  return(out_file);
    
    for (int k=V.KCl-V.Nghost; k<= V.KCu+V.Nghost; ++k) {
        for(int j= V.JCl-V.Nghost; j<= V.JCu+V.Nghost; ++j) {
            for(int i=V.ICl-V.Nghost; i<= V.ICu+V.Nghost; ++i) {
                out_file << V.Position[i][j][k] << " " << V.Velocity[i][j][k] << "\n";
            } 
        } 
    } /* endfor */ 
    
    return (out_file);
    
}

inline istream &operator >> (istream &in_file,
			     Turbulent_Velocity_Field_Block &V) {
   
  int ni, nj, nk, ng, gblknum;
  Vector3D Node_l, Node_u;
   
  in_file.setf(ios::skipws);
  in_file >> ni >> nj >> nk;
  in_file >> ng >> gblknum;
  in_file >> Node_l >> Node_u;
  in_file.unsetf(ios::skipws);

  if (ni == 0 || nj == 0 || nk == 0 || ng == 0) {
    if (V.Allocated) V.deallocate(); 
    return(in_file);
  } 

  if (!V.Allocated || V.NCi != ni || V.NCj != nj
      || V.NCk != nk || V.Nghost != ng) {
    if (V.Allocated) V.deallocate();
    V.allocate(ni-2*ng, nj-2*ng, nk-2*ng, ng);
  } 
   
  V.gblknum = gblknum;
  V.Node_INl_JNl_KNl = Node_l;
  V.Node_INu_JNu_KNu = Node_u;
   
  for (int k=V.KCl-V.Nghost; k<= V.KCu+V.Nghost ; ++k ) {
    for (int j=V.JCl-V.Nghost; j<= V.JCu+V.Nghost; ++j) {
      for (int i=V.ICl-V.Nghost; i<= V.ICu+V.Nghost; ++i) {
	in_file >> V.Position[i][j][k] >> V.Velocity[i][j][k];
      } 
    } 
  } /* endfor */
   
  return (in_file);

}



/**
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

    void Broadcast(void);

    void Create(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                const Grid3D_Input_Parameters &Input);

    void Interpolate_Turbulent_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
				     Turbulent_Velocity_Field_Multi_Block_List &Interpolated_Velocity_Field);
    /* Read and write turbulent velocity field */
    int Write_Turbulent_Velocity_Field(void) const;   
    int Read_Turbulent_Velocity_Field(void);
    
    void Collect_Blocks_from_all_processors(Octree_DataStructure &OcTree,
                                            AdaptiveBlock3D_List &LocalSolnBlockList,
                                            AdaptiveBlock3D_ResourceList &Global_Soln_Block_List) ;
    
    void Assign_Local_Velocity_Field_Blocks(Turbulent_Velocity_Field_Multi_Block_List &Local_List);
    

    
  private:
    //copy and assignment are not permitted
    Turbulent_Velocity_Field_Multi_Block_List(const Turbulent_Velocity_Field_Multi_Block_List &V);
    Turbulent_Velocity_Field_Multi_Block_List &operator = (const Turbulent_Velocity_Field_Multi_Block_List &V);
};


/**
 * Create_Local_Velocity_Field_Multi_Block_List
 *
 * Creates and allocates a Velocity_Field_Multi_Block_List
 * containing velocity field blocks corresponding to local blocks
 *
 * \param [out] Local_List          The list of local velocity field blocks
 * \param [in]  Solution_Block      The local solution blocks
 * \param [in]  LocalSolnBlockList  The local AdaptiveBlock3D_List
 */
template <typename HEXA_BLOCK>
inline void Create_Local_Velocity_Field_Multi_Block_List(Turbulent_Velocity_Field_Multi_Block_List &Local_List, 
                                                         HEXA_BLOCK *Solution_Block, 
                                                         AdaptiveBlock3D_List &LocalSolnBlockList) {
    if (LocalSolnBlockList.Nused() >= 1) {
        Local_List.Allocate(LocalSolnBlockList.Nused());
        for (int nBlk = 0; nBlk <= LocalSolnBlockList.Nused(); ++nBlk ) {
            if (LocalSolnBlockList.Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
                Local_List.Vel_Blks[nBlk].allocate(Solution_Block[nBlk].NCi - 2*Solution_Block[nBlk].Nghost,
                                                   Solution_Block[nBlk].NCj - 2*Solution_Block[nBlk].Nghost,
                                                   Solution_Block[nBlk].NCk - 2*Solution_Block[nBlk].Nghost,
                                                   Solution_Block[nBlk].Nghost);
                Local_List.Vel_Blks[nBlk].gblknum = LocalSolnBlockList.Block[nBlk].info.gblknum;
            }         
        } /* endfor */
    } /* endif */
}



/**
 * Get_Local_Homogeneous_Turbulence_Velocity_Field
 *
 * Fills all the local velocity field blocks with velocity fluctuations around
 * a domain averaged velocity with data from the local solution blocks
 *
 * \param [in]  Solution_Block              All the local Solution Hexa_Blocks
 * \param [in]  LocalSolnBlockList          The AdaptiveBlock3D_List
 * \param [in]  v_average                   The average velocity of the whole domain
 * \param [out] Local_Velocity_Field_List   The list of local Turbulent_Velocity_Field_Blocks
 */
template<typename HEXA_BLOCK>
void Get_Local_Homogeneous_Turbulence_Velocity_Field(HEXA_BLOCK *Solution_Block,
                                                     const AdaptiveBlock3D_List &LocalSolnBlockList,
                                                     const Vector3D &v_average,
                                                     Turbulent_Velocity_Field_Multi_Block_List &Local_Velocity_Field_List) {
    /* Assign initial turbulent velocity field to each solution block. */
    for (int nBlk = 0 ; nBlk <= LocalSolnBlockList.Nblk-1 ; nBlk++) {
        if (LocalSolnBlockList.Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
            if (Local_Velocity_Field_List.Vel_Blks[nBlk].Allocated) {
                Get_Local_Homogeneous_Turbulence_Velocity_Field(Solution_Block[nBlk],
                                                                v_average,
                                                                Local_Velocity_Field_List.Vel_Blks[nBlk]);
            } /* endif */
        } /* endif */
    }  /* endfor */
}


/**
 * Get_Local_Homogeneous_Turbulence_Velocity_Field
 *
 * Copies velocity fluctuations from a local solution block to a
 * local velocity field block.
 *
 * \param [in]  Solution_Block  A local Solution Hexa_Block
 * \param [in]  v_average       The average velocity of the whole domain
 * \param [out] Velocity_Field  A local Turbulent_Velocity_Field_Block
 */
template<typename HEXA_BLOCK>
void Get_Local_Homogeneous_Turbulence_Velocity_Field(HEXA_BLOCK &Solution_Block,
                                                     const Vector3D &v_average,
                                                     Turbulent_Velocity_Field_Block &Velocity_Field) {
    
    /* Assign initial turbulent velocity field to solution block. */
    
    for (int i = Velocity_Field.ICl ; i <= Velocity_Field.ICu ; i++) {
        for (int j = Velocity_Field.JCl ; j <= Velocity_Field.JCu ; j++) {
            for (int k = Velocity_Field.KCl ; k <= Velocity_Field.KCu ; k++) {
                Velocity_Field.Velocity[i][j][k] = Solution_Block.W[i][j][k].v - v_average;
                Velocity_Field.Position[i][j][k] = Solution_Block.Grid.Cell[i][j][k].Xc;
            } /* endfor */
        } /* endfor */
    } /* endfor */
}

/**
 * Assign_Homogeneous_Turbulence_Velocity_Field
 *
 * Superimposes turbulence velocity fluctuations on all the local solution blocks
 * with velocities from the global Turbulent_Velocity_Field_Multi_Block_List
 *
 * \param [in]  Solution_Block        The local Solution Hexa_Blocks
 * \param [in]  LocalSolnBlockList    The local AdaptiveBlock3D_List
 * \param [out] Velocity_Field        The list of all Turbulent_Velocity_Field_Blocks
 */
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

/**
 * Assign_Homogeneous_Turbulence_Velocity_Field
 *
 * Superimposes velocities from a Turbulent_Velocity_Field_Block to
 * a local solution block.
 *
 * \param [out]  Solution_Block   A local Solution Hexa_Block
 * \param [in]   Velocity_Field   A Turbulent_Velocity_Field_Block
 */
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

/**
 * Assign_Homogeneous_Turbulence_Velocity_Field
 *
 * Sets velocities from all local solution block to domain averaged velocity,
 * superimposed with turbulence velocity fluctuations from the global
 * Turbulent_Velocity_Field_Multi_Block_List.
 *
 * \param [in]  Solution_Block        The local Solution Hexa_Blocks
 * \param [in]  LocalSolnBlockList    The local AdaptiveBlock3D_List
 * \param [in]  v_average             The domain averaged velocity
 * \param [out] Velocity_Field        The list of all Turbulent_Velocity_Field_Blocks
 */
template<typename HEXA_BLOCK>
void Assign_Homogeneous_Turbulence_Velocity_Field(HEXA_BLOCK *Solution_Block,
                                                  const AdaptiveBlock3D_List &LocalSolnBlockList,
                                                  const Vector3D &v_average,
                                                  const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field) {
    
    /* Assign initial turbulent velocity field to each solution block. */
    
    for (int nBlk = 0 ; nBlk <= LocalSolnBlockList.Nblk-1 ; nBlk++) {
        if (LocalSolnBlockList.Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
            if (Velocity_Field.Vel_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum].Allocated) {
                Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Block[nBlk],
                                                             Velocity_Field.Vel_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum],
                                                             v_average);
            } /* endif */
        } /* endif */
    }  /* endfor */
    
}

/**
 * Assign_Homogeneous_Turbulence_Velocity_Field
 *
 * Set velocities from a local solution block to the velocities of 
 * a Turbulent_Velocity_Field_Block superimposed with a given average velocity.
 *
 * \param [out]  Solution_Block   A local Solution Hexa_Block
 * \param [in]   Velocity_Field   A Turbulent_Velocity_Field_Block
 * \param [in]   v_average        The domain averaged velocity
 */
template<typename HEXA_BLOCK>
void Assign_Homogeneous_Turbulence_Velocity_Field(HEXA_BLOCK &Solution_Block,
                                                  const Turbulent_Velocity_Field_Block &Velocity_Field,
                                                  const Vector3D &v_average) {
    
    /* Assign initial turbulent velocity field to solution block. */
    
    for (int i = Solution_Block.ICl ; i <= Solution_Block.ICu ; i++) {
        for (int j = Solution_Block.JCl ; j <= Solution_Block.JCu ; j++) {
            for (int k = Solution_Block.KCl ; k <= Solution_Block.KCu ; k++) {
                Solution_Block.W[i][j][k].v = Velocity_Field.Velocity[i][j][k] + v_average;
                Solution_Block.U[i][j][k] = Solution_Block.W[i][j][k].U();
            } /* endfor */
        } /* endfor */
    } /* endfor */
}


/********************************************************
 * Routine: Inflow_Turbulence_XY_Plane                  *
 *                                                      *
 * Feeds turbulent velocity fluctuations on a XY-plane  *
 * of a 3D grid. It applies a turbulent inflow BC.      *
 *                                                      *
 ********************************************************/
template<typename HEXA_BLOCK>
void Inflow_Turbulence_XY_Plane(HEXA_BLOCK &Solution_Block,
				const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field,
				const Input_Parameters<typename HEXA_BLOCK::Soln_pState,
				                       typename HEXA_BLOCK::Soln_cState> &IPs, 
				const double &Time) {

  int N(1);
  double Lz(IPs.Grid_IP.Turbulence_Box_Length);
  double L_convected = double(N)*Lz - Time*IPs.Mean_Velocity.z;

  if ( L_convected < ZERO ) {
    while (true) {
      N++;
      L_convected = double(N)*Lz - Time*IPs.Mean_Velocity.z;
      if (L_convected >= ZERO) break;
    }
  }

  int nBlk, k, ii, jj, kk;
  double new_z, delta, dmin, dmin1, delta1, zmax, zmin, delta_z;  	      
  Vector3D Position_on_Slice, dX;

  Vector3D dVdx, dVdy, dVdz;

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
	              - Solution_Block.Grid.Cell[i][j][Solution_Block.KCu].V/
	                (Solution_Block.Grid.AfaceTop(i, j, Solution_Block.KCu) 
			 + Solution_Block.Grid.AfaceBot(i, j, Solution_Block.KCu));
	  // Bottom boundary
	  } else if (Solution_Block.Grid.BCtypeB[i][j]) {

	    k = Solution_Block.KCl - n_ghost;
	    
	    delta_z = fabs(Solution_Block.Grid.Cell[i][j][k].Xc.z 
			   - Solution_Block.Grid.Cell[i][j][Solution_Block.KCl].Xc.z)
	              - Solution_Block.Grid.Cell[i][j][Solution_Block.KCl].V/
	                (Solution_Block.Grid.AfaceTop(i, j, Solution_Block.KCl) 
			 + Solution_Block.Grid.AfaceBot(i, j, Solution_Block.KCl));
	  }
	  new_z = L_convected - delta_z;
	  if (new_z < NANO) new_z = ZERO;

	  Solution_Block.W[i][j][k].v = IPs.Mean_Velocity;

	  // find the appropriate turbulent field block to interpolate the fluctuations
	  for (nBlk = 0; nBlk < Velocity_Field.NBlk; ++nBlk) {
	    zmax = max(Velocity_Field.Vel_Blks[nBlk].Node_INl_JNl_KNl.z, 
		       Velocity_Field.Vel_Blks[nBlk].Node_INu_JNu_KNu.z);
	    zmin = min(Velocity_Field.Vel_Blks[nBlk].Node_INl_JNl_KNl.z, 
		       Velocity_Field.Vel_Blks[nBlk].Node_INu_JNu_KNu.z);

	    if ( new_z >= zmin  &&  new_z <= zmax ) {

	      // Reconstruct gradients of the velocity field if necessary
/* 	      if ( !Velocity_Field.Vel_Blks[nBlk]._Allocated ) { */
/* 		Velocity_Field.Vel_Blks[nBlk].allocate_gradients(); */
/* 		Velocity_Field.Vel_Blks[nBlk].LeastSquares_Reconstruction(); */
/* 	      } */
	      break;
	    }
	
	    if ( (nBlk == Velocity_Field.NBlk-1) && (new_z < zmin  ||  new_z > zmax) ) {
	      cout << "\n --> n_ghost = " << n_ghost  << " i = " << i 
		   << " j = " << j << " k = " << k;
	      cout << "\n WARNING: new_z = " << new_z << " not found in turbulence blocks.";
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
	  /* Solution_Block.Linear_Reconstruction_LeastSquares(i, j, k, IPs.i_Limiter); */

	  // update k
	  /* Solution_Block.W[i][j][k].k = 0.005*sqr(Solution_Block.W[i][j][k].filter_width()* */
/* 						  Solution_Block.W[i][j][k].abs_strain_rate(Solution_Block.dWdx[i][j][k], */
/* 											    Solution_Block.dWdy[i][j][k], */
/* 											    Solution_Block.dWdz[i][j][k])); */

	  // update the conservative state variables
	  Solution_Block.U[i][j][k] = Solution_Block.W[i][j][k].U();
	}
      }
    }
  }


}


template<typename HEXA_BLOCK>
int Inflow_Turbulence_XY_Plane(HEXA_BLOCK *Solution_Block,
			       const AdaptiveBlock3D_List &LocalSolnBlockList,
			       const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field,
			       const Input_Parameters<typename HEXA_BLOCK::Soln_pState,
				                      typename HEXA_BLOCK::Soln_cState> &IPs,
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
template<typename HEXA_BLOCK>
void IC_Assign_Turbulence_Fresh_Gas(HEXA_BLOCK &Solution_Block,
				    const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field,
				    const Input_Parameters<typename HEXA_BLOCK::Soln_pState,
				                           typename HEXA_BLOCK::Soln_cState> &IPs) {

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
	     (local_X.z <= IPs.Fresh_Gas_Height  &&  fabs(local_X.y) <= 0.0125)) ||
	    // Bunsen burner
	    (IPs.Grid_IP.i_Grid == GRID_BUNSEN_BURNER  &&  
	     (local_X.z <= IPs.Fresh_Gas_Height  &&  (sqr(local_X.x) + sqr(local_X.y) <= sqr(0.0056))))
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


template<typename HEXA_BLOCK>
int IC_Assign_Turbulence_Fresh_Gas(HEXA_BLOCK *Solution_Block,
				   const AdaptiveBlock3D_List &LocalSolnBlockList,
				   const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field,
				   const Input_Parameters<typename HEXA_BLOCK::Soln_pState,
				                          typename HEXA_BLOCK::Soln_cState> &IPs) {
  int error_flag(0);

  for (int Nblk = 0 ; Nblk <= LocalSolnBlockList.Nblk-1 ; Nblk++ ) {
    if (LocalSolnBlockList.Block[Nblk].used == ADAPTIVEBLOCK3D_USED) {
      IC_Assign_Turbulence_Fresh_Gas(Solution_Block[Nblk], Velocity_Field, IPs);
    }
  }

  return error_flag;

}




/**
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
    K_container() : MAX_INDEXES(100), k(0), N(1), Ek(0), velocity_scaling_factor(1), indexcounter(0) {  }
    double k;
    int N;
    int indexcounter;
    int *indexes;
    int index;
    int MAX_INDEXES;
    double Ek;
    double Ek_smooth;
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


/**
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
    
    K_BinaryTree( ) : BinaryTree<K_container>( ) {MAX_INDEXES=0; }
    
    
    void InsertNode(K_container &newItem) {
        InsertNode(K_BinaryTree::root, newItem);
    }
    
private: 
    int MAX_INDEXES;
    void InsertNode(treeNode<K_container> *&root, K_container &newItem) {
        if ( root == NULL ) {
            root = new treeNode<K_container>( newItem );
            
            root->item.MAX_INDEXES =  (root->item.MAX_INDEXES > MAX_INDEXES ? root->item.MAX_INDEXES : MAX_INDEXES);
            root->item.indexes = new int[root->item.MAX_INDEXES];
            root->item.indexes[0] = newItem.index;
            root->item.indexcounter++;
            return;
        }
        else if ( newItem == root->item ) {
            int i = root->item.indexcounter;
            root->item.indexes[i] = newItem.index;
            root->item.indexcounter++;
            root->item.N++;
            if (root->item.indexcounter==root->item.MAX_INDEXES) {
                int *indexes_copy = new int[root->item.MAX_INDEXES];
                for (int j=0; j<root->item.indexcounter; j++){
                    indexes_copy[j] = root->item.indexes[j];
                }
                delete[] root->item.indexes;
                
                root->item.MAX_INDEXES *= 5;
                root->item.indexes = new int[root->item.MAX_INDEXES];
                
                for (int j=0; j<root->item.indexcounter; j++){
                    root->item.indexes[j] = indexes_copy[j];
                }
                
                MAX_INDEXES =  (root->item.MAX_INDEXES > MAX_INDEXES ? root->item.MAX_INDEXES : MAX_INDEXES);
                
                delete[] indexes_copy;
            }
        }
        else if ( newItem < root->item ) {
            InsertNode( root->left, newItem );
        }
        else {
            InsertNode( root->right, newItem );
        }
    }
};


/**
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
    //! Constructor
    RandomFieldRogallo(Input_Parameters<SOLN_pSTATE,SOLN_cSTATE> &IPs);
    
    //! Destructor
    ~RandomFieldRogallo(void);
    
    const char *File_Name;
    
/** @name spectrum imposed members */
/*        ------------------------ */    
//@{    
    string spectrum_name;
    int spectrum_flag;  //!< Turbulence kinetic energy spectrum flag.
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
//@}  
    
/** @name Grid imposed members */
/*        -------------------- */    
//@{
    //! Grid dimensions
    double L1, L2, L3;  
    //!< Number of grid nodes
    int Nx, Ny, Nz, nz;
    //!< Number of cells per block
    int NCells_Idir, NCells_Jdir, NCells_Kdir;
    //! wave numbers
    double k0x, k0y, k0z, abs_k0, kcx, kcy, kcz, abs_kc;
    //! array of all wave numbers visible on grid
    K_container *K;
    //! number of wave numbers visible on grid
    int nK;
    //! k1 = n1 * 2*PI/L1
    double k_1(const int &n1) const {
        if (n1<=Nx/2)   return       n1 * TWO*PI/L1; 
        else            return  (Nx-n1) * TWO*PI/L1;
    }
    //! k2 = n2 * 2*PI/L2
    double k_2(const int &n2) const {
        if (n2<=Ny/2)   return       n2 * TWO*PI/L2; 
        else            return  (Ny-n2) * TWO*PI/L2;
    }
    //! k3 = n3 * 2*PI/L3
    double k_3(const int &n3) const {
        if (n3<=Nz/2)   return       n3 * TWO*PI/L3; 
        else            return  (Nz-n3) * TWO*PI/L3;
    }
    
    //! Assign wave numbers and K_container array
    void Calculate_wave_numbers(void);
    
//@}

    
    
/** @name Rogallo procedure functions */
/*        --------------------------- */    
//@{
    
    //! Generates a random double
    double random_double() const { return  drand48(); }
    
    //! alpha as defined by Rogallo, 1981 for continuous spectra
    complex<double> alpha(const double &abs_wave_num, 
                          const double &theta1, 
                          const double &phi) const;
    
    //! beta as defined by Rogallo, 1981 for continuous spectra
    complex<double> beta(const double &abs_wave_num, 
                         const double &theta2, 
                         const double &phi) const;
    
    void Rogallo_Procedure(void);
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
    
    
/** @name General member functions */
/*        ------------------------ */    
//@{
    
    //! Velocity field
    double        *u, *v, *w;        // Arrays to store the velocity fluctuations in physical space
    fftw_complex  *uu, *vv, *ww;     // Arrays to store the velocity fluctuations in Fourier space
    
    double *dudy;                    // derivative to y of velocity fluctuations in x-direction
    
    //! The model for the energy spectrum
    double Energy_Spectrum_Value(const double &abs_wave_num) const;
    
    //! Calculate the energy spectrum of a given velocity field
    int Get_Energy_Spectrum(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                            int &batch_flag,
                            Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field);
    
    //! Create and assign a homogeneous turbulence velocity field
    int Create_Homogeneous_Turbulence_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                                     int &batch_flag,
                                                     Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field);

    void Import_from_Velocity_Field_Blocks(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                           Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field_List);

    void Export_to_Velocity_Field_Blocks(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                         Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field_List);
    
    double Calculate_Energy_Spectrum(void);
    
    double Rescale_Energy_Spectrum_to_Initial_Spectrum(void);
    
    void FFT_spectral_to_physical(void);
    
    void FFT_physical_to_spectral(void);
    
//@}    
    
/** @name Output the spectrum to a gnuplot file */
/*        ------------------------------------- */    
//@{
    ofstream Turbulence_Spectrum_File;
    int Open_Turbulence_Spectrum_File(void);
    int Output_Turbulence_Spectrum_to_File(void);
    int Close_Turbulence_Spectrum_File(void);
//@}    
    
    
/** @name Statistical Parameters */
/*        ---------------------- */
//@{
    double TKE_model;
    double TKE_grid;
    double TKE_physical;
    double L11_spectral;
    double L11_physical;
    double L22_physical;
    double L33_physical;
    double u2, v2, w2, uRMS;
    
    double Longitudinal_Correlation_spectral(void);
    void Spatial_Averages(void);
    double Longitudinal_Correlation(void);
    double Transversal_Correlation_L22(void);
    double Transversal_Correlation_L33(void);
    double Dissipation(void);
//@}    
    

};

//---------------------------------------------------------------------------//
//                    Members of RandomFieldRogallo class                    //
//---------------------------------------------------------------------------//


/**
 * Constructor
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
RandomFieldRogallo(Input_Parameters<SOLN_pSTATE,SOLN_cSTATE> &IPs) {
    File_Name = IPs.Output_File_Name;
    spectrum_flag = IPs.Turbulence_IP.i_spectrum;
    spectrum_name = IPs.Turbulence_IP.spectrum;
    
    if (IPs.Grid_IP.i_Grid == GRID_BUNSEN_BURNER || IPs.Grid_IP.i_Grid == GRID_BUNSEN_BOX) {
        NCells_Idir = IPs.Grid_IP.NCells_Turbulence_Idir;
        NCells_Jdir = IPs.Grid_IP.NCells_Turbulence_Jdir;
        NCells_Kdir = IPs.Grid_IP.NCells_Turbulence_Kdir;
        
        L1 = IPs.Grid_IP.Turbulence_Box_Width;
        L2 = IPs.Grid_IP.Turbulence_Box_Height;
        L3 = IPs.Grid_IP.Turbulence_Box_Length;
        
    } else {
        NCells_Idir = IPs.Grid_IP.NCells_Idir;
        NCells_Jdir = IPs.Grid_IP.NCells_Jdir;
        NCells_Kdir = IPs.Grid_IP.NCells_Kdir;
        
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
};  

/**
 * Destructor
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
~RandomFieldRogallo(void) {
    for (int ii=0; ii <nK; ii++) {
        delete[] K[ii].indexes;
    }
    delete[] K;
    delete[] dudy;
    fftw_free(u);
    fftw_free(v);
    fftw_free(w);
    fftw_free(uu);
    fftw_free(vv);
    fftw_free(ww);
    // clean up accumulated data by the fftw plan
    fftw_cleanup();  
}


/**
 * Prescribed energy spectrum
 */
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


/**
 * alpha as defined in Rogallo's procedure
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
alpha(const double &abs_wave_num, const double &theta1, const double &phi) const {
    
    double E, k;
    k = abs_wave_num;  E = Energy_Spectrum_Value(k);
    return (k == 0.0)  ?  0.0 : sqrt(TWO*E/(4.0*PI*k*k)) * exp(I*theta1) * cos(phi);
}

/**
 * beta as defined in Rogallo's procedure
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
beta(const double &abs_wave_num, const double &theta2, const double &phi) const {
    
    double E, k;
    k = abs_wave_num;  E = Energy_Spectrum_Value(k);
    return (k == 0.0)  ?  0.0 : sqrt(TWO*E/(4.0*PI*k*k)) * exp(I*theta2) * sin(phi);    
}


/**
 * RandomFieldRogallo::Calculate_wave_numbers()
 *
 * Calculates grid characteristic wave numbers
 * and fills the K_container with all possible
 * absolute wave numbers on the grid, saving the
 * index of the each combination.
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Calculate_wave_numbers(void) {
    
    
    k0x = TWO*PI/L1;    // smallest wavenumber on grid in x-direction
    k0y = TWO*PI/L2;    // smallest wavenumber on grid in y-direction
    k0z = TWO*PI/L3;    // smallest wavenumber on grid in z-direction
    abs_k0 = min(min(k0x,k0y),k0z);
    kcx = Nx/2*k0x;     // largest wavenumber on grid in x-direction
    kcy = Ny/2*k0y;     // largest wavenumber on grid in y-direction
    kcz = Nz/2*k0z;     // largest wavenumber on grid in z-direction
    abs_kc = sqrt( sqr(kcx) + sqr(kcy) + sqr(kcz) );
    
    K = new K_container [Nx*Ny*nz];  // Container of information for a wavenumber
    K_BinaryTree Ktree;              // A binary tree that will order and store K_containers
    int index;
    double k1,k2,k3,abs_k;
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
    
    nK = Ktree.countNodes();        // The total number of unique values of abs_k
    
    K = Ktree.asArray();            // Make an array of the tree nodes in order.
    
}

/**
 * RandomFieldRogallo::Import_from_Velocity_Field_Blocks
 *
 * Assigns the vectors u, v, w used in FFT from the
 * Turbulent_Velocity_Field_Multi_Block_List
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Import_from_Velocity_Field_Blocks(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                  Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field_List) {
    int index;
    int nBlk, ix, iy, iz;
    Vector3D dVdx, dVdy, dVdz;
    for (int kBlk = 0; kBlk <= Initial_Mesh.NBlk_Kdir-1; ++kBlk) {
        for (int jBlk = 0; jBlk <= Initial_Mesh.NBlk_Jdir-1; ++jBlk) {
            for (int iBlk = 0; iBlk <= Initial_Mesh.NBlk_Idir-1; ++iBlk) {
                nBlk = iBlk + 
                jBlk*Initial_Mesh.NBlk_Idir + 
                kBlk*Initial_Mesh.NBlk_Idir*Initial_Mesh.NBlk_Jdir;
                for (int i = Initial_Mesh.Grid_Blks[nBlk].ICl; i <= Initial_Mesh.Grid_Blks[nBlk].ICu; ++i) {
                    for (int j = Initial_Mesh.Grid_Blks[nBlk].JCl; j <= Initial_Mesh.Grid_Blks[nBlk].JCu; ++j) {
                        for (int k = Initial_Mesh.Grid_Blks[nBlk].KCl; k <= Initial_Mesh.Grid_Blks[nBlk].KCu; ++k) {
                            ix = iBlk*NCells_Idir+(i-Initial_Mesh.Grid_Blks[nBlk].ICl);
                            iy = jBlk*NCells_Jdir+(j-Initial_Mesh.Grid_Blks[nBlk].JCl);
                            iz = kBlk*NCells_Kdir+(k-Initial_Mesh.Grid_Blks[nBlk].KCl);
                            index = iz + Nz*(iy + ix*Ny);
                            u[index] = Velocity_Field_List.Vel_Blks[nBlk].Velocity[i][j][k].x;
                            v[index] = Velocity_Field_List.Vel_Blks[nBlk].Velocity[i][j][k].y;
                            w[index] = Velocity_Field_List.Vel_Blks[nBlk].Velocity[i][j][k].z;
                            Velocity_Field_List.Vel_Blks[nBlk].LeastSquares_Reconstruction(i,j,k,dVdx,dVdy,dVdz);
                            dudy[index] = dVdy.x;
                        }
                    }
                }
            }
        }
    } 
};


/**
 * RandomFieldRogallo::Import_from_Velocity_Field_Blocks
 *
 * Assigns the Turbulent_Velocity_Field_Multi_Block_List
 * from the vectors u, v, w used in FFT
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Export_to_Velocity_Field_Blocks(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field_List) {
    int index;
    int nBlk, ix, iy, iz;
    int INl,INu,JNl,JNu,KNl,KNu;
    for (int kBlk = 0; kBlk <= Initial_Mesh.NBlk_Kdir-1; ++kBlk) {
        for (int jBlk = 0; jBlk <= Initial_Mesh.NBlk_Jdir-1; ++jBlk) {
            for (int iBlk = 0; iBlk <= Initial_Mesh.NBlk_Idir-1; ++iBlk) {
                nBlk = iBlk + 
                jBlk*Initial_Mesh.NBlk_Idir + 
                kBlk*Initial_Mesh.NBlk_Idir*Initial_Mesh.NBlk_Jdir;
                for (int i = Initial_Mesh.Grid_Blks[nBlk].ICl; i <= Initial_Mesh.Grid_Blks[nBlk].ICu; ++i) {
                    for (int j = Initial_Mesh.Grid_Blks[nBlk].JCl; j <= Initial_Mesh.Grid_Blks[nBlk].JCu; ++j) {
                        for (int k = Initial_Mesh.Grid_Blks[nBlk].KCl; k <= Initial_Mesh.Grid_Blks[nBlk].KCu; ++k) {
                            ix = iBlk*NCells_Idir+(i-Initial_Mesh.Grid_Blks[nBlk].ICl);
                            iy = jBlk*NCells_Jdir+(j-Initial_Mesh.Grid_Blks[nBlk].JCl);
                            iz = kBlk*NCells_Kdir+(k-Initial_Mesh.Grid_Blks[nBlk].KCl);
                            index = iz + Nz*(iy + ix*Ny);
                            Velocity_Field_List.Vel_Blks[nBlk].Velocity[i][j][k].x = u[index];
                            Velocity_Field_List.Vel_Blks[nBlk].Velocity[i][j][k].y = v[index];
                            Velocity_Field_List.Vel_Blks[nBlk].Velocity[i][j][k].z = w[index];
                            Velocity_Field_List.Vel_Blks[nBlk].Position[i][j][k] = Initial_Mesh.Grid_Blks[nBlk].Cell[i][j][k].Xc;	                                 
                        }
                    }
                }
                INl = Initial_Mesh.Grid_Blks[nBlk].INl; 
                INu = Initial_Mesh.Grid_Blks[nBlk].INu; 
                JNl = Initial_Mesh.Grid_Blks[nBlk].JNl; 
                JNu = Initial_Mesh.Grid_Blks[nBlk].JNu; 
                KNl = Initial_Mesh.Grid_Blks[nBlk].KNl; 
                KNu = Initial_Mesh.Grid_Blks[nBlk].KNu;
                Velocity_Field_List.Vel_Blks[nBlk].Node_INl_JNl_KNl = Initial_Mesh.Grid_Blks[nBlk].Node[INl][JNl][KNl].X; 
                Velocity_Field_List.Vel_Blks[nBlk].Node_INu_JNu_KNu = Initial_Mesh.Grid_Blks[nBlk].Node[INu][JNu][KNu].X;
            }
        } 
    } 
};

/**
 * RandomFieldRogallo::Spatial_Averages
 *
 * Calculates the Total Kinetic Energy and the
 * RMS velocity, assuming the grid is uniform!
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::Spatial_Averages(void) {
    int index;
    double dudy2 = 0;
    u2 = 0 , v2 = 0 , w2 = 0;
    double invN = ONE/(Nx*Ny*Nz);
    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny; j++) {
            for (int l=0; l<Nz; l++) {
                index = l + Nz*(j+Ny*i);                
                u2 += invN * sqr(u[index]);   // average of u^2
                v2 += invN * sqr(v[index]);   // average of v^2
                w2 += invN * sqr(w[index]);   // average of w^2
                dudy2 += invN * sqr(dudy[index]);
            } 
        } 
    } 
    TKE_physical = HALF*(u2+v2+w2);
    uRMS = sqrt(TWO/THREE*TKE_physical);
    lambda_g = sqrt(u2/(dudy2+PICO));
}


/**
 * RandomFieldRogallo::Longitudinal_Correlation
 *
 * Calculates the longitudinal velocity autocorrelation 
 * lengthscale L11 in physical space,
 * assuming the grid is uniform!
 * 
 * \return L11: longitudinal velocity autocorrelation lengthscale
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
double RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::Longitudinal_Correlation(void) {
    
    double *Rr = new double [Nx];
    double *fr = new double [Nx];
    double *r  = new double [Nx];
    for (int i=0; i<Nx; i++){ // loop for variation in r
        Rr[i] = 0;
    }
    int ri, indexr, index;
    for (int l=0; l<Nz; l++) {
        for (int j=0; j<Ny; j++) {
            for (int i=0; i<Nx; i++) {
                for (int ii=i; ii<Nx; ii++) {
                    indexr = l + Nz*( j + Ny*ii );
                    index  = l + Nz*( j + Ny*i );
                    ri = ii-i;
                    Rr[ri] += ONE/(Nx*Ny*Nz) * (u[indexr]*u[index]); 
                }
            }
        }
    }
    double dr = L1/(Nx-1.0);
    for (int i = 0; i<Nx; i++) {
        r[i] = i*dr;
        fr[i] = Rr[i]/u2;
    }
    
    L11_physical = CubicSplinesIntegration(r,fr,Nx);

    
#ifdef _GNUPLOT
    Gnuplot_Control h2;
    h2.gnuplot_init();
    h2.gnuplot_setstyle("lines") ;
    h2.gnuplot_cmd("set grid");
    h2.gnuplot_set_xlabel("r");
    h2.gnuplot_set_ylabel("f(r)");
    h2.gnuplot_set_title("longitudinal velocity autocorrelation function");
    h2.gnuplot_plot1d_var2(r,fr,Nx,"");
#endif
    
    
    delete[] Rr;
    delete[] fr;
    delete[] r;
    
    return L11_physical;    // Note: this is a class member
}


/**
 * RandomFieldRogallo::Transversal_Correlation_L22
 *
 * Calculates the transversal velocity autocorrelation 
 * lengthscale L22 in physical space,
 * assuming the grid is uniform!
 * 
 * \return L22: transversal velocity autocorrelation lengthscale
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
double RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::Transversal_Correlation_L22(void) {
    
    double *Rr = new double [Ny];
    double *fr = new double [Ny];
    double *r  = new double [Ny];
    for (int j=0; j<Ny; j++){ // loop for variation in r
        Rr[j] = 0;
    }
    int indexr, index, rj;
    for (int l=0; l<Nz; l++) {
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny; j++) {
                for (int jj=j; jj<Ny; jj++) {
                    indexr = l + Nz*( jj + Ny*i );
                    index  = l + Nz*( j + Ny*i );
                    rj = jj-j;
                    Rr[rj] += ONE/(Nx*Ny*Nz) * (v[indexr]*v[index]); 
                }
            }
        }
    }
    double dr = L2/(Ny-1.0);
    for (int j = 0; j<Ny; j++) {
        r[j] = j*dr;
        fr[j] = Rr[j]/v2;
    }
    
    L22_physical = CubicSplinesIntegration(r,fr,Ny);
    
    delete[] Rr;
    delete[] fr;
    delete[] r;
    
    return L22_physical;    // Note: this is a class member
}


/**
 * RandomFieldRogallo::Transversal_Correlation_L33
 *
 * Calculates the transversal velocity autocorrelation 
 * lengthscale L33 in physical space,
 * assuming the grid is uniform!
 * 
 * \return L33: transversal velocity autocorrelation lengthscale
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
double RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::Transversal_Correlation_L33(void) {
    
    double *Rr = new double [Nz];
    double *fr = new double [Nz];
    double *r  = new double [Nz];
    for (int l=0; l<Nz; l++){ // loop for variation in r
        Rr[l] = 0;
    }
    int indexr, index, rl;
    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny; j++) {
            for (int l=0; l<Nz; l++) {
                for (int ll=l; ll<Nz; ll++) {
                    indexr = ll + Nz*( j + Ny*i );
                    index  = l + Nz*( j + Ny*i );
                    rl = ll-l;
                    Rr[rl] += ONE/(Nx*Ny*Nz) * (w[indexr]*w[index]); 
                }
            }
        }
    }
    double dr = L3/(Nz-1.0);
    for (int l = 0; l<Nz; l++) {
        r[l] = l*dr;
        fr[l] = Rr[l]/w2;
    }
    
    L33_physical = CubicSplinesIntegration(r,fr,Nz);
    
    delete[] Rr;
    delete[] fr;
    delete[] r;
    
    return L33_physical;    // Note: this is a class member
}


/**
 * RandomFieldRogallo::Longitudinal_Correlation_spectral
 *
 * Calculates longitudinal velocity autocorrelation 
 * lengthscale in spectral space,
 * assuming the turbulence is isotropic!
 *
 * \return L11: longitudinal velocity autocorrelation lengthscale
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
double RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::Longitudinal_Correlation_spectral(void) {
    
    double *Ek_over_k = new double [nK];
    double *k = new double [nK];
    for (int ii=1; ii <nK; ii++) {
        k[ii-1] = K[ii].k;
        Ek_over_k[ii-1] = K[ii].Ek/K[ii].k; 
    }
    L11_spectral = THREE*PI/FOUR * CubicSplinesIntegration(k,Ek_over_k,nK-1) / TKE_grid;
    
    delete[] Ek_over_k;
    delete[] k;
    
    return L11_spectral;  // Note: this is a class member
}

 

/**
 * RandomFieldRogallo::Dissipation
 *
 * Calculates the dissipation in spectral space
 * which is visible on the grid only
 *
 * \return Dissipation
 */ 
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
double RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::Dissipation(void) {
    
    double *Dk_smooth = new double [nK];
    double *k = new double [nK];
    for (int ii=0; ii <nK; ii++) {
        k[ii] = K[ii].k;
        Dk_smooth[ii] = TWO*nu*sqr(K[ii].k)*K[ii].Ek_smooth; 
        // For non smooth data: Dk[ii] = TWO*nu*sqr(K[ii].k)*K[ii].Ek; 
        
    }
    double dissipation = CubicSplinesIntegration(k,Dk_smooth,nK);
    
#ifdef _GNUPLOT
    int npts = int((k_eta-abs_k0)/abs_k0);
    double kmin = abs_k0/10;
    double kmax = k_eta/4.0;
    dpoint *dp = new dpoint[npts];
    dpoint *Ep = new dpoint[npts];
    for (int i=0; i<npts; i++) {
        dp[i].x = kmin + i*(kmax-kmin)/(npts-1);
        dp[i].y = TWO*nu*sqr(dp[i].x)*Energy_Spectrum_Value(dp[i].x);
    }
    Gnuplot_Control h1;
    h1.gnuplot_init();
    h1.gnuplot_setstyle("lines") ;
    h1.gnuplot_cmd("set grid");
    h1.gnuplot_cmd("set logscale xy");
    h1.gnuplot_set_xlabel("k");
    h1.gnuplot_set_ylabel("D(k)");
    h1.gnuplot_set_title("Dissipation spectrum");
    
    // the initial dissipation spectrum
    h1.gnuplot_plot1d_var2(dp,npts,"dissipation model");
    delete[] dp;
    
    // the current dissipation spectrum
    h1.gnuplot_plot1d_var2(k,Dk_smooth,nK,"dissipation on grid");
#endif
    
    delete[] Dk_smooth;
    delete[] k;
    
    return dissipation;
}



/**
 * RandomFieldRogallo::Rogallo_Procedure
 *
 * The spectral velocity components that are used in FFT are 
 * assigned using Rogallo's procedure described in paper:
 * "Numerical Experiments in Homogeneous Turbulence" (1981)
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::Rogallo_Procedure(void) {
    
    /* --------------- Create spectral velocity fluctuations ---------------- */
    int seed = 1; // = time(NULL);      // assigns the current time to the seed
    srand48(seed);                      // changes the seed for drand48()
    double k1, k2, k3, abs_k;           // Wave numbers
    double theta1, theta2, phi;         // random angles for Rogallo spectrum
    complex<double> aa, bb;             // complex energy components for Rogallo spectrum
    double deno;                        // denominator in Rogallo spectrum function
    int index, i, j, l;                 // indexes
    int index_conj, iconj, jconj;
    
    for(int ii=0; ii<nK; ii++) {                // For every abs_k
        for(int jj=0; jj<K[ii].N; jj++) {       // For all grid points with this value of abs_k
            index = K[ii].indexes[jj];      // index corresponding to (i,j,l)
            i = index/(Ny*nz);              // convert index to (i,j,l) coordinates
            j = (index - Ny*nz*i)/nz;
            l = index - nz*(j+Ny*i);
            
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
            
            aa = alpha(abs_k, theta1, phi);
            bb = beta (abs_k, theta2, phi);
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
                    iconj = (i==0  ?  0 : Nx-i);
                    jconj = (j==0  ?  0 : Ny-j);
                    index_conj = l + nz*(jconj+Ny*iconj);
                    
                    // complex conjugates
                    uu[index][0] =  uu[index_conj][0];
                    uu[index][1] = -uu[index_conj][1];
                    vv[index][0] =  vv[index_conj][0];
                    vv[index][1] = -vv[index_conj][1];
                    ww[index][0] =  ww[index_conj][0];
                    ww[index][1] = -ww[index_conj][1];
                    
                } else if ( (i==0 || i==Nx/2)  &&  (j==0 || j==Ny/2) ) {
                    // real values at 8 corners
                    uu[index][0] = sqrt( sqr(uu[index][0]) + sqr(uu[index][1]) );
                    uu[index][1] = 0.0;    
                    vv[index][0] = sqrt( sqr(vv[index][0]) + sqr(vv[index][1]) );
                    vv[index][1] = 0.0;
                    ww[index][0] = sqrt( sqr(ww[index][0]) + sqr(ww[index][1]) );
                    ww[index][1] = 0.0; 
                }
            }                
        }
    }
};

/**
 * RandomFieldRogallo::Calculate_Energy_Spectrum
 *
 * The Kinetic Energy spectrum is calculated from the
 * spectral velocity components. The Total Kinetic Energy (TKE)
 * is calculated by integrating the spectrum.
 *
 * \return Total Kinetic Energy
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
double RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::Calculate_Energy_Spectrum(void) {
    
    double *k  = new double [nK];
    double *Ek = new double [nK];
    double *Ek_smooth = new double [nK];
    double *Emodel = new double [nK];
    
    double *EE = new double [Nx*Ny*nz];   // Energy of one wavenumber vector
    double Rxx, Ryy, Rzz;
    
    int index, i, j, l;                 // indexes
    double abs_k;
    
    for(int ii=0; ii<nK; ii++) {                // For every abs_k
        K[ii].Ek = ZERO;
        for(int jj=0; jj<K[ii].N; jj++) {       // For all grid points with this value of abs_k
            index = K[ii].indexes[jj];      // index corresponding to (i,j,l)
            i = index/(Ny*nz);              // convert index to (i,j,l) coordinates
            j = (index - Ny*nz*i)/nz;
            l = index - nz*(j+Ny*i);
            
            abs_k = K[ii].k;
            
            /* ------ Calculate energy corresponding with this wave number ------ */
            // Velocity correlation tensor for this wavenumber
            Rxx = (sqr(uu[index][0]) + sqr(uu[index][1]));
            Ryy = (sqr(vv[index][0]) + sqr(vv[index][1]));
            Rzz = (sqr(ww[index][0]) + sqr(ww[index][1]));
            
            // Energy = HALF * trace(velocitycorrelation tensor Rij)
            EE[index] = HALF * (Rxx + Ryy + Rzz) ;
            
            /* ------- Calculate energy corresponding with this wavenumber ------ */
            K[ii].Ek += TWO * EE[index] * (TWO*PI*sqr(abs_k)) / K[ii].N;
            
        }
        K[ii].Emodel = Energy_Spectrum_Value(K[ii].k);
        
        if (K[ii].Ek != ZERO)
            K[ii].velocity_scaling_factor = sqrt(K[ii].Emodel/(K[ii].Ek));
        else 
            K[ii].velocity_scaling_factor = ONE;
        
    }
    
    /* ------------------- Total Kinetic Energy ------------------- */
    for (int ii=0; ii<nK; ii++) {
        k[ii]  = K[ii].k;
        Ek[ii] = K[ii].Ek; 
        Emodel[ii] = K[ii].Emodel;
    }
    
    TKE_model = CubicSplinesIntegration(k,Emodel,nK);
    TKE_grid  = CubicSplinesIntegration(k,Ek,nK);
    
    /* --------------- Calculate Smoothed spectrum ---------------- */
    polyfit_smoothing_logscale(nK,k,Ek,3,0.1,true,Ek_smooth);
    for (int ii=0; ii<nK; ii++) {
        K[ii].Ek_smooth = Ek_smooth[ii];
    }
        
    delete[] k;
    delete[] EE;
    delete[] Ek;
    delete[] Ek_smooth;
    delete[] Emodel;
    
    return TKE_grid;  // Note: This is a class member
};


/**
 * RandomFieldRogallo::Rescale_Energy_Spectrum_to_Initial_Spectrum
 *
 * The spectral velocity components are being rescaled to match the same
 * spectral content as the initial Kinetic Energy Spectrum
 *
 * \return Total Kinetic Energy
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
double RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Rescale_Energy_Spectrum_to_Initial_Spectrum(void) {
    
    double *k  = new double [nK];
    double *Ek = new double [nK];
    double *Ek_smooth = new double [nK];
    
    double *EE = new double [Nx*Ny*nz];
    int index;
    double velocity_scaling_factor;
    double abs_k;
    double Rxx, Ryy, Rzz;
    for (int ii=0; ii<nK; ii++) {
        K[ii].Ek = 0;
        abs_k = K[ii].k;
        for (int jj=0; jj<K[ii].N; jj++) {
            index = K[ii].indexes[jj];
            velocity_scaling_factor = K[ii].velocity_scaling_factor;
            uu[index][0] *= velocity_scaling_factor;
            uu[index][1] *= velocity_scaling_factor;
            vv[index][0] *= velocity_scaling_factor;
            vv[index][1] *= velocity_scaling_factor;
            ww[index][0] *= velocity_scaling_factor;
            ww[index][1] *= velocity_scaling_factor;
            
            Rxx = (sqr(uu[index][0]) + sqr(uu[index][1]));
            Ryy = (sqr(vv[index][0]) + sqr(vv[index][1]));
            Rzz = (sqr(ww[index][0]) + sqr(ww[index][1]));
            EE[index] = HALF * (Rxx + Ryy + Rzz);
            
            K[ii].Ek += TWO * EE[index] * (TWO*PI*sqr(abs_k)) / K[ii].N;
        }
        
    } 
    
    /* ------------------- Total Kinetic Energy ------------------- */
    for (int ii=0; ii <nK; ii++) {
        k[ii]  = K[ii].k;   // Don't wish to store K[0].k which is 0
        Ek[ii] = K[ii].Ek; 
    }
    
    TKE_grid = CubicSplinesIntegration(k,Ek,nK);
    
    /* --------------- Calculate Smoothed spectrum ---------------- */
    polyfit_smoothing_logscale(nK,k,Ek,3,0.1,true,Ek_smooth);
    for (int ii=0; ii<nK; ii++) {
        K[ii].Ek_smooth = Ek_smooth[ii];
    }
    
    delete[] k;
    delete[] Ek;
    delete[] Ek_smooth;
    delete[] EE;
    
    return TKE_grid;  // Note: this is a class member
}

/**
 * RandomFieldRogallo::FFT_spectral_to_physical
 *
 * The physical velocity components are calculated
 * through FFT of the spectral velocity components.
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::FFT_spectral_to_physical(void) {
    
    fftw_plan physical;
    physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, uu, u, FFTW_ESTIMATE);
    fftw_execute(physical); 
    fftw_destroy_plan(physical);
    
    physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, vv, v, FFTW_ESTIMATE);
    fftw_execute(physical); 
    fftw_destroy_plan(physical);
    
    physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, ww, w, FFTW_ESTIMATE);
    fftw_execute(physical); 
    fftw_destroy_plan(physical);
    
    /* ---------- FFT scaling factor ---------- */
    double FFT_scaling_factor = sqrt(cube(TWO*PI)/(L1*L2*L3));
    int index;
    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny; j++) {
            for (int l=0; l<Nz; l++) {
                index = l + Nz*(j+Ny*i);
                u[index] *= FFT_scaling_factor;
                v[index] *= FFT_scaling_factor;
                w[index] *= FFT_scaling_factor;
            } 
        } 
    }
};

/**
 * RandomFieldRogallo::FFT_physical_to_spectral
 *
 * The spectral velocity components are calculated
 * through FFT of the physical velocity components.
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::FFT_physical_to_spectral(void) {
    
    fftw_plan spectral;
    spectral = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, u, uu, FFTW_ESTIMATE);
    fftw_execute(spectral); 
    fftw_destroy_plan(spectral);
    
    spectral = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, v, vv, FFTW_ESTIMATE);
    fftw_execute(spectral); 
    fftw_destroy_plan(spectral);
    
    spectral = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, w, ww, FFTW_ESTIMATE);
    fftw_execute(spectral); 
    fftw_destroy_plan(spectral);
    
    /* --------------- FFT scaling factor ---------------- */
    
    double FFT_scaling_factor = ONE/(Nx*Ny*Nz)   / sqrt(pow(TWO*PI,3)/(L1*L2*L3));
    int index;
    for (int ii=0; ii<nK; ii++) {
        for (int jj=0; jj<K[ii].N; jj++) {
            index = K[ii].indexes[jj];
            
            uu[index][0] *= FFT_scaling_factor;
            uu[index][1] *= FFT_scaling_factor;
            vv[index][0] *= FFT_scaling_factor;
            vv[index][1] *= FFT_scaling_factor;
            ww[index][0] *= FFT_scaling_factor;
            ww[index][1] *= FFT_scaling_factor;
        }
    }
};



/**
 * Subroutine: Create_Homogeneous_Turbulence_Velocity_Field
 * 
 * Creates and assigns a turbulence velocity field
 *
 *
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Create_Homogeneous_Turbulence_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                             int &batch_flag,
                                             Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field) {
    
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) {
        cout << "\n\n";
        cout << " ==========================================================================\n"; 
        cout << "         Generating Homogeneous Isotropic Turbulent Velocity Field \n";
        cout << " ==========================================================================\n";
        cout << endl;
        cout << "   Turbulence statistics of the initial spectrum (k=0:infinity) : " << endl;
        cout << "    -->  spectrum     = " << spectrum_name << endl; 
        cout << "    -->  TKE          = " << TKE << endl;
        cout << "    -->  u_RMS        = " << up << endl;
        cout << "    -->  L            = " << Lp << endl;
        cout << "    -->  L11 ~ 0.43 L = " << Lp*0.43 << endl;
        cout << "    -->  eps          = " << eps << endl;
        cout << "    -->  eta          = " << eta << endl;
        cout << "    -->  ReL          = " << ReL << endl;
        cout << "    -->  Rlambda      = " << Rlambda << endl;
        cout << "    -->  lambda_g     = " << lambda_g << endl;
    }
    
    
    int errorflag;
    errorflag = Open_Turbulence_Spectrum_File( );
    bool rescale = true;
    
    
    //  Allocate_spectrum();
    Nx = NCells_Idir * Initial_Mesh.NBlk_Idir;
    Ny = NCells_Jdir * Initial_Mesh.NBlk_Idir;
    Nz = NCells_Kdir * Initial_Mesh.NBlk_Kdir;
    nz = Nz/2+1;
    
    // Allocation of arrays used in the transforms
    u = (double *) malloc(Nx*Ny*Nz * sizeof(double));
    v = (double *) malloc(Nx*Ny*Nz * sizeof(double));
    w = (double *) malloc(Nx*Ny*Nz * sizeof(double));
    uu = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
    vv = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
    ww = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
    fftw_plan      physical;
    
    
    // Allocation of dudy (needs to be allocated)
    dudy = new double [Nx*Ny*Nz];
    
    
    
    Calculate_wave_numbers();
    
    
    
    
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
    h1.gnuplot_set_title("Kinetic energy spectrum");
    h1.gnuplot_plot1d_var2(dp,npts,"model");
    sleep(2);
    delete[] dp;
#endif
    
    
    /* --------------- Create spectral velocity fluctuations ---------------- */
    Rogallo_Procedure();
    
    
    /* --------------------------- Energy spectrum -------------------------- */
    Calculate_Energy_Spectrum();
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) {       
        cout << endl;
        cout << "   Turbulence statistics in spectral space: " << endl;
        cout << "    -->  TKE of the model            = " << TKE_model << endl;
        cout << "    -->  TKE with Rogallo            = " << TKE_grid << endl;
    }
    
    /* ---------------- gnuplot rogallo spectrum -------------- */
#ifdef _GNUPLOT
    double *k = new double [nK];
    double *Ek = new double [nK];
    for (int ii=0; ii <nK; ii++) {
        k[ii]  = K[ii].k;   
        Ek[ii] = K[ii].Ek; 
    }
    h1.gnuplot_setstyle("points");
    h1.gnuplot_plot1d_var2(k,Ek,nK,"Rogallo");
#endif
    
    
    
    rescale = true;
    /* ------- rescale the rogallo spectrum to match the model spectrum ------- */
    if (rescale) { 
        Rescale_Energy_Spectrum_to_Initial_Spectrum();
        
        if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) 
            cout << "    -->  TKE with Rogallo rescaled   = " << TKE_grid << endl;
    }
    
    /* ------------------ write gnuplot spectrum files --------------- */
    Output_Turbulence_Spectrum_to_File();
    
    
    /* ------------------ gnuplot the rescaled spectrum -------------- */
#ifdef _GNUPLOT
    if (rescale) {
        for (int ii=0; ii <nK; ii++) {
            k[ii]  = K[ii].k; 
            Ek[ii] = K[ii].Ek; 
        }
        h1.gnuplot_plot1d_var2(k,Ek,nK,"Rogallo - rescaled"); 
    }
    delete[] k;
    delete[] Ek;
#endif
    
    /* ------------ Percentage of TKE on the grid ----------- */
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) {
        cout << "    -->  TKE percentage on the grid  = " << int(TKE_grid/TKE*HUNDRED) << " %" << endl;
    }
    
    
    
    /* ----------------- calculate L11 ---------------------- */
    Longitudinal_Correlation_spectral();
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) {
        cout.setf(ios::fixed, ios::floatfield);
        cout.precision(4);    
        cout << "    -->  L11                         = " << setw(7) << L11_spectral
             << "      Lx/L11  = " << setw(7) << L1 / L11_spectral
             << "      L11/L  = " << setw(7) << L11_spectral  / Lp << endl;
        cout.unsetf(ios::fixed);
        cout.unsetf(ios::floatfield);
    }
    
    /* --------------------- dissipation -------------------- */
    double diss = Dissipation();
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) )
        cout << "    -->  eps on grid                 = " << diss << endl;

    
    
    /* ---------------------------------------------------------------------------- *
     *               Fourier transform from spectral to physical space              *
     * ---------------------------------------------------------------------------- */
    
    FFT_spectral_to_physical();
    
    Spatial_Averages();
    
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) {
        cout << endl;
        cout << "   Turbulence statistics in physical space:" << endl;
        cout << "    -->  TKE    = " << TKE_physical << endl;
        cout << "    -->  u_RMS  = " << uRMS << endl;
        cout << "    -->  <u^2>  = " << u2 << endl;
        cout << "    -->  <v^2>  = " << v2 << endl;
        cout << "    -->  <w^2>  = " << w2 << endl;
    }
    
    
    /* --------------- Velocity auto-correlation functions -------------- */
    Longitudinal_Correlation();
    Transversal_Correlation_L22();
    Transversal_Correlation_L33();
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) {
        cout.setf(ios::fixed, ios::floatfield);
        cout.precision(4);    
        cout << "    -->  L11    = " << setw(7) << L11_physical 
             << "      Lx/L11  = " << setw(7) << L1 / L11_physical
             << "      L11/L  = " << setw(7) << L11_physical  / Lp << endl;
        cout << "    -->  L22    = " << setw(7) << L22_physical 
             << "      Ly/L22  = " << setw(7) << L2 / L22_physical
             << "      L22/L  = " << setw(7) << L22_physical  / Lp << endl;
        cout << "    -->  L33    = " << setw(7) << L33_physical 
             << "      Lz/L33  = " << setw(7) << L3 / L33_physical
             << "      L33/L  = " << setw(7) << L33_physical  / Lp << endl;
        cout.unsetf(ios::fixed);
        cout.unsetf(ios::floatfield);
    }
    
    /* ----------------- Assign turbulent velocity field ------------------ */
    Export_to_Velocity_Field_Blocks(Initial_Mesh, Initial_Velocity_Field);
    
    
    
    errorflag = Close_Turbulence_Spectrum_File( );      
    return (0);
}



template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Get_Energy_Spectrum(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                    int &batch_flag,
                    Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field_List) {
    
    
    /* ------- assign grid parameters --------- */
    Nx = NCells_Idir * Initial_Mesh.NBlk_Idir;
    Ny = NCells_Jdir * Initial_Mesh.NBlk_Idir;
    Nz = NCells_Kdir * Initial_Mesh.NBlk_Kdir;
    nz = Nz/2+1;
    
    /* ------- allocate physical and spectral velocity fluctuation arrays ---------- */    
    u = (double *) malloc(Nx*Ny*Nz * sizeof(double));
    v = (double *) malloc(Nx*Ny*Nz * sizeof(double));
    w = (double *) malloc(Nx*Ny*Nz * sizeof(double));
    uu = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
    vv = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
    ww = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
    
    
    /* ------------ allocate dudy --------------- */
    dudy = new double[Nx*Ny*Nz];

    /* ----- assign the wave number list and some key wave numbers ----- */
    Calculate_wave_numbers();
    
    
    /* ------------ assign the physical velocity fluctuations ------------- */
    Import_from_Velocity_Field_Blocks(Initial_Mesh,
                                      Velocity_Field_List); 
    
    
    
    /* --------------- Turbulent kinetic energy in physical space --------------- */
    Spatial_Averages();
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) {
        cout << endl;
        cout << "   Turbulence statistics in physical space:" << endl;
        cout << "    -->  TKE      = " << TKE_physical << endl;
        cout << "    -->  u_RMS    = " << uRMS << endl;
        cout << "    -->  <u^2>    = " << u2 << endl;
        cout << "    -->  <v^2>    = " << v2 << endl;
        cout << "    -->  <w^2>    = " << w2 << endl;
        cout << "    -->  lambda_g = " << lambda_g << endl;
        
    }
    
    /* --------------- Velocity auto-correlation functions -------------- */
    Longitudinal_Correlation();
    Transversal_Correlation_L22();
    Transversal_Correlation_L33();
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) {
        cout.setf(ios::fixed, ios::floatfield);
        cout.precision(4);        
        cout << "    -->  L11      = " << setw(7) << L11_physical 
        << "    Lx/L11  = " << setw(6) << L1 / L11_physical << endl;
        cout << "    -->  L22      = " << setw(7) << L22_physical 
        << "    Ly/L22  = " << setw(6) << L2 / L22_physical << endl;
        cout << "    -->  L33      = " << setw(7) << L33_physical 
        << "    Lz/L33  = " << setw(6) << L3 / L33_physical << endl;
        cout.unsetf(ios::fixed);
        cout.unsetf(ios::floatfield);
    }
    
    /* ---------------------------------------------------------------------------- *
     *          Perform Fourier transform from physical to spectral space           *
     * ---------------------------------------------------------------------------- */
    
    FFT_physical_to_spectral();
    
    Calculate_Energy_Spectrum();
    
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) {       
        cout << endl;
        cout << "   Turbulence statistics in spectral space: " << endl;
        cout << "    -->  TKE on the grid  = " << TKE_grid << endl;
    }
    
    /* ----------------- calculate L11 ---------------------- */
    Longitudinal_Correlation_spectral();
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) ) {
        cout.setf(ios::fixed, ios::floatfield);
        cout.precision(4); 
        cout << "    -->  L11              = " << setw(6) << L11_spectral
             << "      Lx/L11  = " << setw(6) << L1 / L11_spectral << endl;
        cout.unsetf(ios::fixed);
        cout.unsetf(ios::floatfield);
    }
    
    /* --------------------- dissipation -------------------- */
    double diss = Dissipation();    
    if (CFFC_Primary_MPI_Processor() && !(batch_flag) )
        cout << "    -->  eps on grid      = " << diss << endl;

    
    /* --------------- gnuplot initial spectrum -------------- */
#ifdef _GNUPLOT
    int npts = int((k_eta-abs_k0)/abs_k0);
    double kmin = abs_k0/10;
    double kmax = k_eta/4.0;
    dpoint *energy = new dpoint[npts];
    for (int i=0; i<npts; i++) {
        energy[i].x= kmin + i*(kmax-kmin)/(npts-1);
        energy[i].y= Energy_Spectrum_Value(energy[i].x);
    }
    Gnuplot_Control h1;
    h1.gnuplot_init();
    h1.gnuplot_setstyle("lines") ;
    h1.gnuplot_cmd("set logscale xy");
    h1.gnuplot_cmd("set grid");
    h1.gnuplot_set_xlabel("k");
    h1.gnuplot_set_ylabel("E(k)");
    h1.gnuplot_set_title("Kinetic energy spectrum");
    h1.gnuplot_plot1d_var2(energy,npts,"initial");
    
    delete[] energy;
#endif
    
    /* ------------------ gnuplot the numerical spectrum -------------- */
#ifdef _GNUPLOT
    double *k  = new double [nK];
    double *Ek = new double [nK];
    for (int ii=0; ii <nK; ii++) {
        k[ii] = K[ii].k;
        Ek[ii] = K[ii].Ek; 
    }
    h1.gnuplot_setstyle("points");
    h1.gnuplot_plot1d_var2(k,Ek,nK,"on grid");
    sleep(2);
    delete[] k;
    delete[] Ek;
#endif
    
    Open_Turbulence_Spectrum_File();
    Output_Turbulence_Spectrum_to_File();
    Close_Turbulence_Spectrum_File();
    
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
    
    gnuplot_file 
    << "set title \"Turbulence kinetic energy spectrum "<< spectrum_name <<"\"\n"
    << "set xlabel \"wavenumber k \"\n"
    << "set ylabel \"EE\"\n" 
    << "set logscale xy\n"
    << "plot \"" << turbulence_spectrum_file_name_ptr << "\" using 1:2 \\\n"
    << "     title \"" << "initial condition"    << "\" with lines , \\\n"
    
    << "\""      << turbulence_spectrum_file_name_ptr << "\" using 1:3 \\\n"
    << "     title \"" << "LES"                  << "\" with points , \\\n"
    
    << "\""      << turbulence_spectrum_file_name_ptr << "\" using 1:4 \\\n"
    << "     title \"" << "smoothed"             << "\" with lines \n"
    
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
Output_Turbulence_Spectrum_to_File() {
    
    Turbulence_Spectrum_File << setprecision(6);
    Turbulence_Spectrum_File.setf(ios::scientific);
    
    for (int ii=0; ii <nK; ii++) {
        Turbulence_Spectrum_File    << K[ii].k  << " " << K[ii].Emodel
        << " " << K[ii].Ek 
        << " " << K[ii].Ek_smooth
        << "\n";
    } 
    
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
