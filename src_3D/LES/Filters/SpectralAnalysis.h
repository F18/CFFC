/*
 *  SpectralAnalysis.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 13/04/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

//# include "../../TurbulenceModelling/TurbulentVelocityField.h"

#include "fftw3.h"
#include "LES_Filters.h"
//typedef complex<double> Complex

/*!
 * Class: Scalar_Field_Block
 *
 * \brief Class used to store turbulent velocity field for a single block.
 *
 */
class Scalar_Field_Block {
public:
    int                 gblknum; // global turbulent velocity field number
    int             NCi,ICl,ICu; // i-direction turbulent velocity field cell counters
    int             NCj,JCl,JCu; // j-direction turbulent velocity field cell counters
    int             NCk,KCl,KCu; // k-direction turbulent velocity field cell counters
    int                  Nghost; // number of ghost cells
    double        ***s; // array of turbulent velocity field vectors
    

    int  Allocated; // Indicates whether or not the turbulent velocity field data have been allocated.
    
    /* Constructors. */
    Scalar_Field_Block(void) {
        NCi = 0; ICl = 0; ICu = 0; 
        NCj = 0; JCl = 0; JCu = 0;
        NCk = 0; KCl = 0; KCu = 0;
        Nghost = 0;
        gblknum = 0;
        Allocated = TURBULENT_VELOCITY_FIELD_DATA_NOT_USED;
        s = NULL;
 
    }
    
    Scalar_Field_Block(const int Ni, 
                                   const int Nj, 
                                   const int Nk,
                                   const int Ng) {
        allocate(Ni, Nj, Nk, Ng);
    }
    
    //! Destructor
    ~Scalar_Field_Block(void) {
        deallocate();        
    }
    
    //! Allocate memory for velocity field data. 
    void allocate(const int Ni, 
                  const int Nj, 
                  const int Nk,
                  const int Ng);
    
    //! Deallocate memory for velocity field data.
    void deallocate(void);
    
    //! Copy operator
    void Copy(Scalar_Field_Block &Block2);
    
    //! Broadcast
    void Broadcast(void);
    
#ifdef _MPI_VERSION   
    //! Broadcast from a given processor
    void Broadcast(MPI::Intracomm &Communicator, 
                   const int Source_CPU);
#endif
    
    //! Output operator
    friend ostream &operator << (ostream &out_file, 
                                 const Scalar_Field_Block &V);
    //! Input operator
    friend istream &operator >> (istream &in_file, 
                                 Scalar_Field_Block &V);
    
    /* Other useful member functions. */
    
private:
    //copy and assignment are not permitted
    Scalar_Field_Block(const Scalar_Field_Block &V);
    Scalar_Field_Block &operator =(const Scalar_Field_Block &V);
};

/*************************************************************************
 * Scalar_Field_Block -- Input-output operators.             *
 *************************************************************************/
inline ostream &operator << (ostream &out_file,
                             const Scalar_Field_Block &S){  
    
    out_file << S.NCi << " " << S.NCj << " " << S.NCk << "\n";
    out_file << S.Nghost << " " << S.gblknum << "\n";
    
    if (S.NCi == 0 || S.NCj == 0 || S.NCk == 0 || S.Nghost == 0)  return(out_file);
    
    for (int k=S.KCl-S.Nghost; k<= S.KCu+S.Nghost; ++k) {
        for(int j= S.JCl-S.Nghost; j<= S.JCu+S.Nghost; ++j) {
            for(int i=S.ICl-S.Nghost; i<= S.ICu+S.Nghost; ++i) {
                out_file << S.s[i][j][k] << "\n";
            } 
        } 
    } /* endfor */ 
    
    return (out_file);
    
}

inline istream &operator >> (istream &in_file,
                             Scalar_Field_Block &V) {
    
    int ni, nj, nk, ng, gblknum;
    Vector3D Node_l, Node_u;
    
    in_file.setf(ios::skipws);
    in_file >> ni >> nj >> nk;
    in_file >> ng >> gblknum;
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
    
    for (int k=V.KCl-V.Nghost; k<= V.KCu+V.Nghost ; ++k ) {
        for (int j=V.JCl-V.Nghost; j<= V.JCu+V.Nghost; ++j) {
            for (int i=V.ICl-V.Nghost; i<= V.ICu+V.Nghost; ++i) {
                in_file >> V.s[i][j][k];
            } 
        } 
    } /* endfor */
    
    return (in_file);
    
}



/**
 * Class: Turbulent_Scalar_Field_Multi_Block
 *
 * \brief Class used to store turbulent velocity field for a
 *        1D array of turbulent velocity field solution blocks.
 *
 */
class Scalar_Field_Multi_Block_List {
public:
    Scalar_Field_Block  *Scalar_Field_Blks; // one dimensional array of velocity block.
    int                                  NBlk;
    int                             NBlk_Idir, 
    NBlk_Jdir, 
    NBlk_Kdir; // Number of blocks in i, j and k directions.
    int                             Allocated; // Indicates if the velocity blocks have been allocated or not.
    
    /* Creation constructors. */
    Scalar_Field_Multi_Block_List(void) : 
    NBlk(0), NBlk_Idir(0), NBlk_Jdir(0), NBlk_Kdir(0), Scalar_Field_Blks(NULL), Allocated(0) { }
    
    Scalar_Field_Multi_Block_List(const int N) {
        Allocate(N);
    }
    
    Scalar_Field_Multi_Block_List(const int Ni, 
                                  const int Nj, 
                                  const int Nk) {
        Allocate(Ni, Nj, Nk);
    }
    
    /* Destructor. */
    ~Scalar_Field_Multi_Block_List(void) {
        Deallocate();
    }
    
    /* Other member functions  */
    
    void Allocate(const int Ni, const int Nj, const int Nk);
    
    void Allocate(const int N);
    
    void Deallocate(void);
    
    void Broadcast(void);
    
    void Create_Global(Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                       Grid3D_Input_Parameters &Input);
    
    template <typename Soln_pState, typename Soln_cState>
    void Create_Local(Hexa_Block<Soln_pState,Soln_cState> *Solution_Block, 
                             AdaptiveBlock3D_List &LocalSolnBlockList);
    
    void Assign_Local_Scalar_Field_Blocks(Scalar_Field_Multi_Block_List &Local_List);
    
    
private:
    //copy and assignment are not permitted
    Scalar_Field_Multi_Block_List(const Scalar_Field_Multi_Block_List &V);
    Scalar_Field_Multi_Block_List &operator = (const Scalar_Field_Multi_Block_List &V);
};


/* -------------------------------------------------------------------------- */
/*                      Scalar_Field_Block member functions                   */
/* -------------------------------------------------------------------------- */



/*************************************************************************
 * Scalar_Field_Block::allocate -- Allocate memory.          *
 *************************************************************************/
void Scalar_Field_Block::allocate(const int Ni, 
                                              const int Nj, 
                                              const int Nk,
                                              const int Ng) {
    assert( Ni >= 1 && Nj >= 1 && Nk >= 1 && Ng >=1 && !Allocated);
    NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
    NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1;
    NCk = Nk+2*Ng; KCl = Ng; KCu = Nk+Ng-1;
    Nghost = Ng;
    Allocated = TURBULENT_VELOCITY_FIELD_DATA_USED;
    
    s = new double**[NCi];
    for (int i = 0; i <= NCi-1; ++i ){
        s[i] = new double*[NCj];
        for (int j = 0; j <= NCj-1; ++j ){
            s[i][j] = new double[NCk];
        } /* endfor */
    } /* endfor */
    
}


/*************************************************************************
 * Scalar_Field_Block::deallocate -- Deallocate memory.      *
 *************************************************************************/
void Scalar_Field_Block::deallocate(void) {
    if (Allocated) {
        assert(NCi >= 1 && NCj >= 1 && NCk >= 1);
        for (int i = 0; i <= NCi-1 ; ++i ) {
            for ( int j = 0 ; j <= NCj-1 ; ++j) {
                delete []s[i][j]; s[i][j] = NULL;
            } /* endfor */
            delete []s[i]; s[i] = NULL;
        }/*endfor*/
        delete []s; s = NULL;
        
        NCi = 0; ICl = 0; ICu = 0; 
        NCj = 0; JCl = 0; JCu = 0; 
        NCk = 0; KCl = 0; KCu = 0; 
        Nghost = 0;
        Allocated = TURBULENT_VELOCITY_FIELD_DATA_NOT_USED;
    } /* endif */
}



/*************************************************************************
 * Scalar_Field_Block::Copy -- Copy a block.                 *
 *************************************************************************/
void Scalar_Field_Block::Copy(Scalar_Field_Block &Block2){
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
                    s[i][j][k] = Block2.s[i][j][k];
                } /* endfor */
            } /* endfor */
        } /* endfor */
        
        gblknum = Block2.gblknum;        
    }
    
}
















/*********************************************************************
 * Scalar_Field_Block::Broadcast                         *
 *                     -- Broadcast turbulent velocity field block.  *
 *********************************************************************/
void Scalar_Field_Block::Broadcast(void) {
    
#ifdef _MPI_VERSION
    
    int ni, nj, nk, ng, block_allocated, buffer_size;
    double *buffer;
        
    /* Broadcast the number of cells in each direction. */
    
    if (CFFC_Primary_MPI_Processor()) {
        ni = NCi;
        nj = NCj;
        nk = NCk;
        ng = Nghost;
        
        if (Allocated)
            block_allocated = 1;
        else
            block_allocated = 0;
        
    } /* endif */
    
    MPI::COMM_WORLD.Bcast(&ni, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nj, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nk, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&ng, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&block_allocated, 1, MPI::INT, 0);
    
    /* On non-primary MPI processors, allocate (re-allocate) 
     memory for the block as necessary. */
    
    if (!CFFC_Primary_MPI_Processor()) {
        if (NCi != ni || NCj != nj ||  NCk != nk || Nghost != ng ) { 
            if (block_allocated) allocate(ni-2*ng, nj-2*ng, nk-2*ng, ng); 
        }
    } 
    
    /* Broadcast the global block number */
    
    MPI::COMM_WORLD.Bcast(&gblknum, 1, MPI::INT, 0);
    
    /* Broadcast the scalar quantity */
    
    if (block_allocated) {
        ni = (ICu+Nghost) - (ICl-Nghost) + 1;
        nj = (JCu+Nghost) - (JCl-Nghost) + 1;
        nk = (KCu+Nghost) - (KCl-Nghost) + 1;
        
        buffer_size = (ni*nj*nk);
        buffer = new double[buffer_size]; 
        
        if (CFFC_Primary_MPI_Processor()) {
            int i_buffer = 0;
            for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
                for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
                    for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                        buffer[i_buffer] = s[i][j][k];
                        i_buffer ++;
                    }
                }
            } /* endfor */ 
        } /* endif */
	    
        MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);
        
        if (!CFFC_Primary_MPI_Processor()) {
            int i_buffer = 0;
            for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
                for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
                    for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
                        s[i][j][k] = buffer[i_buffer];
                        i_buffer ++;
                    }	 
                }
            } /* endfor */
        } /* endif */
        
        delete []buffer; 
        buffer = NULL;
        
    } /* endif */
#endif
    
}


/*********************************************************************
 * Scalar_Field_Block::Broadcast                         *
 *                     -- Broadcast turbulent velocity field block.  *
 *********************************************************************/
#ifdef _MPI_VERSION
void Scalar_Field_Block::Broadcast(MPI::Intracomm &Communicator, 
                                               const int Source_CPU){
    
    int Source_Rank = 0;
    int ni, nj, nk, ng, block_allocated, buffer_size;
    double *buffer;
        
    /* Broadcast the number of cells in each direction. */
    
    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
        ni = NCi;
        nj = NCj;
        nk = NCk;
        ng = Nghost; 
        
        if (Allocated)
            block_allocated = 1;
        else
            block_allocated = 0;
    } /* endif */
    
    Communicator.Bcast(&ni, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&nj, 1, MPI::INT, Source_Rank); 
    Communicator.Bcast(&nk, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&ng, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&block_allocated, 1, MPI::INT, Source_Rank);
    
    
    /* On non-source MPI processors, allocate (re-allocate) 
     memory for the block as necessary. */
    
    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
        if (NCi != ni || NCj != nj ||  NCk != nk || Nghost != ng) { 
            if (block_allocated) {
                allocate(ni-2*ng, nj-2*ng, nk-2*ng, ng); 
            }
        } /* endif */
    } /* endif */
    
    
    /* Broadcast the global block number */
    
    Communicator.Bcast(&gblknum, 1, MPI::INT, Source_Rank);

    /* Broadcast the scalar quantity */
    
    if (block_allocated) {
        ni = (ICu+Nghost) - (ICl-Nghost) + 1;
        nj = (JCu+Nghost) - (JCl-Nghost) + 1;
        nk = (KCu+Nghost) - (KCl-Nghost) + 1;
        
        buffer_size = (ni*nj*nk);
        buffer = new double[buffer_size];
        
        if (CFFC_MPI::This_Processor_Number == Source_CPU) {
            int i_buffer = 0;
            for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
                for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
                    for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
                        buffer[i_buffer] = s[i][j][k];
                        i_buffer ++;
                    }	 
                }
            } /* endfor */
        } /* endif */
        
        Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);
        
        if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
            int i_buffer = 0;
            for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
                for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
                    for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
                        s[i][j][k] = buffer[i_buffer];
                        i_buffer ++;
                    }	 
                }
            } /* endfor */
        } /* endif */
        
        delete []buffer; 
        buffer = NULL;
        
    } /* endif */ 
    
}
#endif 



/* -------------------------------------------------------------------------- */
/*          Scalar_Field_Multi_Block_List member functions        */
/* -------------------------------------------------------------------------- */



/*****************************************************************************
 * Scalar_Field_Multi_Block_List::Allocate -- Allocate memory.   *
 *****************************************************************************/
void Scalar_Field_Multi_Block_List::Allocate(const int Ni, 
                                                         const int Nj, 
                                                         const int Nk) {
    if (Ni >= 1 && Nj >= 1 && Nk >= 1 && !Allocated) {
        NBlk_Idir = Ni; 
        NBlk_Jdir = Nj; 
        NBlk_Kdir = Nk; 
        NBlk = Ni*Nj*Nk;
        Scalar_Field_Blks = new Scalar_Field_Block[NBlk];
        Allocated = 1;
    } /* endif */
}

/*****************************************************************************
 * Scalar_Field_Multi_Block_List::Allocate -- Allocate memory.   *
 *****************************************************************************/
void Scalar_Field_Multi_Block_List::Allocate(const int N) {
    if (N >= 1 && !Allocated) {
        NBlk_Idir = N; 
        NBlk_Jdir = 1; 
        NBlk_Kdir = 1; 
        NBlk = N;
        Scalar_Field_Blks = new Scalar_Field_Block[NBlk];
        Allocated = 1;
    } /* endif */
}

/*****************************************************************************
 * Scalar_Field_Multi_Block_List::Allocate -- Deallocate memory. *
 *****************************************************************************/
void Scalar_Field_Multi_Block_List::Deallocate(void) {
    if (NBlk >= 1 && Allocated) {
        delete []Scalar_Field_Blks;
        Scalar_Field_Blks = NULL;
        NBlk_Idir = 0; 
        NBlk_Jdir = 0; 
        NBlk_Kdir = 0;
        NBlk = 0;
        Allocated = 0;
    } /* endif */
}




/**
 * Scalar_Field_Multi_Block_List::Create
 * 
 * Allocate a full scalar field multi block list using the initial mesh
 * for grid information
 */
void Scalar_Field_Multi_Block_List::Create_Global(Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                                  Grid3D_Input_Parameters &Input) {
    if (Initial_Mesh.NBlk >= 1 && Initial_Mesh.Allocated) {
        Allocate(Initial_Mesh.NBlk_Idir, Initial_Mesh.NBlk_Jdir, Initial_Mesh.NBlk_Kdir);
        for (int nBlk = 0; nBlk <= NBlk-1; ++nBlk ) {
            if (Initial_Mesh.Grid_Blks[nBlk].Allocated) Scalar_Field_Blks[nBlk].allocate(Initial_Mesh.Grid_Blks[nBlk].NCi - 2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                                         Initial_Mesh.Grid_Blks[nBlk].NCj - 2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                                         Initial_Mesh.Grid_Blks[nBlk].NCk - 2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                                         Initial_Mesh.Grid_Blks[nBlk].Nghost);
        } 
    } 
}



/**
 * Create_Local_Scalar_Field_Multi_Block_List
 *
 * Creates and allocates a Scalar_Field_Multi_Block_List
 * containing velocity field blocks corresponding to local blocks
 *
 * \param [out] Local_List          The list of local velocity field blocks
 * \param [in]  Solution_Block      The local solution blocks
 * \param [in]  LocalSolnBlockList  The local AdaptiveBlock3D_List
 */
template <typename Soln_pState, typename Soln_cState>
inline void Scalar_Field_Multi_Block_List::Create_Local(Hexa_Block<Soln_pState,Soln_cState> *Solution_Block, 
                                                        AdaptiveBlock3D_List &LocalSolnBlockList) {
    if (LocalSolnBlockList.Nused() >= 1) {
        Allocate(LocalSolnBlockList.Nused());
        for (int nBlk = 0; nBlk <= LocalSolnBlockList.Nused(); ++nBlk ) {
            if (LocalSolnBlockList.Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
                Scalar_Field_Blks[nBlk].allocate(Solution_Block[nBlk].NCi - 2*Solution_Block[nBlk].Nghost,
                                                 Solution_Block[nBlk].NCj - 2*Solution_Block[nBlk].Nghost,
                                                 Solution_Block[nBlk].NCk - 2*Solution_Block[nBlk].Nghost,
                                                 Solution_Block[nBlk].Nghost);
                Scalar_Field_Blks[nBlk].gblknum = LocalSolnBlockList.Block[nBlk].info.gblknum;
            }         
        }
    }
}







/*****************************************************************************
 * Scalar_Field_Multi_Block_List::Broadcast --                   *
 *     Broadcast the turbulent velocity field data from the main processor.  *
 *****************************************************************************/
void Scalar_Field_Multi_Block_List::Broadcast() {
    
#ifdef _MPI_VERSION
    
    int nblki, nblkj, nblkk, Block_allocated;
    
    if ( CFFC_Primary_MPI_Processor() ) {
        nblki = NBlk_Idir;
        nblkj = NBlk_Jdir;
        nblkk = NBlk_Kdir;
        if (Allocated)
            Block_allocated = 1;
        else
            Block_allocated = 0;  
    }
    
    MPI::COMM_WORLD.Bcast(&nblki, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nblkj, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nblkk, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&Block_allocated, 1, MPI::INT, 0);
    
    /* On non-primary MPI processors, allocate (re-allocate)
     memory as necessary. */
    
    if ( !CFFC_Primary_MPI_Processor() ) {
        if (NBlk_Idir != nblki || NBlk_Jdir != nblkj ||  NBlk_Kdir != nblkk) {
            if (Block_allocated) Allocate(nblki, nblkj, nblkk);
        }
    } 
    
    for (int nBlk = 0; nBlk < NBlk; ++nBlk) {
        Scalar_Field_Blks[nBlk].Broadcast();
    }
    
#endif
}
























/* --------------------------------------------------------------------------------------------------------------- */
/*                                                 STANDALONES                                                     */
/* --------------------------------------------------------------------------------------------------------------- */






///**
// * Assign_Homogeneous_Turbulence_Scalar_Field
// *
// * Superimposes turbulence velocity fluctuations on all the local solution blocks
// * with velocities from the global Scalar_Field_Multi_Block_List
// *
// * \param [in]  Solution_Block        The local Solution Hexa_Blocks
// * \param [in]  LocalSolnBlockList    The local AdaptiveBlock3D_List
// * \param [out] Scalar_Field        The list of all Scalar_Field_Blocks
// */
//template<typename HEXA_BLOCK>
//void Assign_Scalar_Field(HEXA_BLOCK *Solution_Block,
//                         const AdaptiveBlock3D_List &LocalSolnBlockList,
//                         const Scalar_Field_Multi_Block_List &Scalar_Field) {
//    
//    /* Assign initial turbulent velocity field to each solution block. */
//    
//    for (int nBlk = 0 ; nBlk <= LocalSolnBlockList.Nblk-1 ; nBlk++) {
//        if (LocalSolnBlockList.Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
//            if (Scalar_Field.Scalar_Field_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum].Allocated) {
//                Assign_Scalar_Field(Solution_Block[nBlk],
//                                    Scalar_Field.Scalar_Field_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum]);
//            } /* endif */
//        } /* endif */
//    }  /* endfor */
//    
//}
//
///**
// * Assign_Homogeneous_Turbulence_Scalar_Field
// *
// * Superimposes velocities from a Scalar_Field_Block to
// * a local solution block.
// *
// * \param [out]  Solution_Block   A local Solution Hexa_Block
// * \param [in]   Scalar_Field   A Scalar_Field_Block
// */
//template<typename HEXA_BLOCK>
//void Assign_Scalar_Field(HEXA_BLOCK &Solution_Block,
//                         const Scalar_Field_Block &Scalar_Field) {
//    
//    /* Assign initial turbulent velocity field to solution block. */
//    
//    for (int i = Solution_Block.ICl ; i <= Solution_Block.ICu ; i++) {
//        for (int j = Solution_Block.JCl ; j <= Solution_Block.JCu ; j++) {
//            for (int k = Solution_Block.KCl ; k <= Solution_Block.KCu ; k++) {
//                Solution_Block.W[i][j][k].rho = Scalar_Field.s[i][j][k];
//            } /* endfor */
//        } /* endfor */
//    } /* endfor */
//}
//
///**
// * Assign_Scalar_Field
// *
// * Sets velocities from all local solution block to domain averaged velocity,
// * superimposed with turbulence velocity fluctuations from the global
// * Scalar_Field_Multi_Block_List.
// *
// * \param [in]  Solution_Block        The local Solution Hexa_Blocks
// * \param [in]  LocalSolnBlockList    The local AdaptiveBlock3D_List
// * \param [in]  v_average             The domain averaged velocity
// * \param [out] Scalar_Field        The list of all Scalar_Field_Blocks
// */
//template<typename HEXA_BLOCK>
//void Assign_Scalar_Field(HEXA_BLOCK *Solution_Block,
//                                                  const AdaptiveBlock3D_List &LocalSolnBlockList,
//                                                  const Vector3D &v_average,
//                                                  const Scalar_Field_Multi_Block_List &Scalar_Field) {
//    
//    /* Assign initial turbulent velocity field to each solution block. */
//    
//    for (int nBlk = 0 ; nBlk <= LocalSolnBlockList.Nblk-1 ; nBlk++) {
//        if (LocalSolnBlockList.Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
//            if (Scalar_Field.Scalar_Field_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum].Allocated) {
//                Assign_Homogeneous_Turbulence_Scalar_Field(Solution_Block[nBlk],
//                                                             Scalar_Field.Scalar_Field_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum],
//                                                             v_average);
//            } /* endif */
//        } /* endif */
//    }  /* endfor */
//    
//}
//
///**
// * Assign_Scalar_Field
// *
// * Set velocities from a local solution block to the velocities of 
// * a Scalar_Field_Block superimposed with a given average velocity.
// *
// * \param [out]  Solution_Block   A local Solution Hexa_Block
// * \param [in]   Scalar_Field   A Scalar_Field_Block
// * \param [in]   v_average        The domain averaged velocity
// */
//template<typename HEXA_BLOCK>
//void Assign_Scalar_Field(HEXA_BLOCK &Solution_Block,
//                         const Scalar_Field_Block &Scalar_Field,
//                         const Vector3D &v_average) {
//    
//    /* Assign initial turbulent velocity field to solution block. */
//    
//    for (int i = Solution_Block.ICl ; i <= Solution_Block.ICu ; i++) {
//        for (int j = Solution_Block.JCl ; j <= Solution_Block.JCu ; j++) {
//            for (int k = Solution_Block.KCl ; k <= Solution_Block.KCu ; k++) {
//                Solution_Block.W[i][j][k].rho = Scalar_Field.s[i][j][k];
//                Solution_Block.U[i][j][k] = Solution_Block.W[i][j][k].U();
//            } /* endfor */
//        } /* endfor */
//    } /* endfor */
//}


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
class Spectral_container {
public:
    Spectral_container() : MAX_INDEXES(100), k(0), N(1), absS(0), realS(0), imagS(0) , indexcounter(0) {  }
    double k;
    double absS;
    double realS;
    double imagS;
    double rescaling_factor;
    int N;
    int indexcounter;
    int *indexes;
    int index;
    int MAX_INDEXES;
    // Relational operators
    //! W == W
    friend bool operator ==(const Spectral_container &K1, const Spectral_container &K2) {
        return ((fabs((K1.k - K2.k)/(K1.k+NANO)) <= TOLER));
    }
    //! W != W
    friend bool operator !=(const Spectral_container &K1, const Spectral_container &K2) {
        return (K1.k != K2.k);
    }
    //! W <= W
    friend bool operator <=(const Spectral_container &K1, const Spectral_container &K2) {
        return (K1.k <= K2.k);
    }
    //! W >= W
    friend bool operator >=(const Spectral_container &K1, const Spectral_container &K2) {
        return (K1.k >= K2.k);
    }
    //! W < W
    friend bool operator <(const Spectral_container &K1, const Spectral_container &K2) {
        return (K1.k < K2.k);
    }
    //! W > W
    friend bool operator >(const Spectral_container &K1, const Spectral_container &K2) {
        return (K1.k > K2.k);
    }
    
    // Input-output operators. 
    //! ostream << W
    friend ostream& operator << (ostream &out_file, const Spectral_container &K) {
        out_file.setf(ios::scientific);
        out_file << " " << K.k  << "     " << K.N << "      " << K.absS ;
        out_file.unsetf(ios::scientific);
        return (out_file);
    }
    //! istream >> W
    //friend istream& operator >> (istream &in_file,  K_container &K); 
};


/**
 * Class: Spectral_BinaryTree
 *
 * \brief This class inherits from BinaryTree<itemtype> to adapt to perform on Spectral_container
 *
 * A class that inherits from BinaryTree<itemtype>. The functioncall InsertNode
 * is overloaded so that an item with the same value of k as an existing node 
 * won't generate a new node in the tree but will instead store the index
 * of the item in the existing node
 *
 * This class is needed in the class RandomFieldRogallo to efficiently sort and 
 * store all the Spectral_container.
 */
class Spectral_BinaryTree : public BinaryTree<Spectral_container> {
    
public:
    
    Spectral_BinaryTree( ) : BinaryTree<Spectral_container>( ) {MAX_INDEXES=0; }
    
    
    void InsertNode(Spectral_container &newItem) {
        InsertNode(Spectral_BinaryTree::root, newItem);
    }
    
private: 
    int MAX_INDEXES;
    void InsertNode(treeNode<Spectral_container> *&root, Spectral_container &newItem) {
        if ( root == NULL ) {
            root = new treeNode<Spectral_container>( newItem );
            
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







template<typename Soln_pState, typename Soln_cState>
class SpectralAnalysis {
public:
    
    AdaptiveBlock3D_List                    *LocalSolnBlkList_ptr;  // List with properties of SolnBlks
    Hexa_Block<Soln_pState,Soln_cState>     *Solution_Blocks_ptr;   // array of SolnBlks
    HexaSolver_Data                         *Data_ptr;
    Scalar_Field_Multi_Block_List           Global_Scalar_Field;
    Scalar_Field_Multi_Block_List           Local_Scalar_Field;
    
    /* ----- constructor ----- */
    SpectralAnalysis(HexaSolver_Data &Data,
                     HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data);    

    
    
    

    
    
    
    /** @name Grid imposed members */
    /*        -------------------- */    
    //@{
    //! Grid dimensions
    double L1, L2, L3, Ls;  
    //!< Number of grid nodes
    int Nx, Ny, Nz, nz;
    //!< Number of cells per block
    int NCells_Idir, NCells_Jdir, NCells_Kdir;
    //! wave numbers
    double k0x, k0y, k0z, abs_k0, kcx, kcy, kcz, abs_kc;
    //! array of all wave numbers visible on grid
    Spectral_container *K;
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
    
    

    
    
    /** @name General member functions */
    /*        ------------------------ */    
    //@{
    
    //! Velocity field
    double        *s;        // Array to store the scalar quantity in physical space
    fftw_complex  *ss;       // Array to store the scalar quantity in Fourier space
    
    void FFT_spectral_to_physical(void);
    void FFT_physical_to_spectral(void);
    
    //@}    
    
    /** @name Output the spectrum to a gnuplot file */
    /*        ------------------------------------- */    
    //@{
    const char *File_Name;
    const char *scalar_name;
    ofstream Spectrum_File;
    int Open_Spectrum_File(char *file_suffix);
    int Output_Spectrum_to_File(void);
    int Close_Spectrum_File(void);
    void Average_Spectrum();
    
    void Output_Spectrum_1D(char *file_suffix){
        Average_Spectrum();
        Open_Spectrum_File(file_suffix);
        Output_Spectrum_to_File();
        Close_Spectrum_File();
    }
    
    //@}    
    
    
    
    /** @name Functions working with Scalar_Field_Blocks */
    /*        ------------------------------------------ */
    //@{
    
    template <typename Function_type>
    void Assign_Local_Scalar_Field(Function_type func);
    
    void Assign_Local_Scalar_Field_Block(Hexa_Block<Soln_pState,Soln_cState> &Solution_Block,
                                         Scalar_Field_Block &Scalar_Field,
                                         double (Soln_pState::*&member));
    
    void Assign_Local_Scalar_Field_Block(Hexa_Block<Soln_pState,Soln_cState> &Solution_Block,
                                         Scalar_Field_Block &Scalar_Field,
                                         double (Soln_pState::*&function)(void));
    
    void Assign_Local_Scalar_Field_Block(Hexa_Block<Soln_pState,Soln_cState> &Solution_Block,
                                         Scalar_Field_Block &Scalar_Field,
                                         double (Soln_pState::*&function)(const Soln_pState&, const Soln_pState&, const Soln_pState&));
    
    void Assign_Local_Scalar_Field_Block(Hexa_Block<Soln_pState,Soln_cState> &Solution_Block,
                                         Scalar_Field_Block &Scalar_Field,
                                         double (Soln_pState::*&function)(const Soln_pState&, const Soln_pState&, const Soln_pState&, const double&));

    
    
    void Gather_from_all_processors(void);
    
    void Import_From_Global_Scalar_Field(void);
    void Export_To_Local_Scalar_Field(void);

    
    template <typename Function_type>
    void Assign_Scalar_Field(Function_type func) {
        Assign_Local_Scalar_Field(func);
        Gather_from_all_processors();
        Import_From_Global_Scalar_Field();
    }
    
    template <typename Function_type>
    void Get_Spectrum(Function_type func, char* scalar_name){
        Assign_Scalar_Field(func);
        Calculate_wave_numbers();
        FFT_physical_to_spectral();
        Output_Spectrum_1D(scalar_name);
        /* do it all */
    }
    
    template <typename Member_type>
    void Set_Spectrum(Member_type member) {
        Calculate_wave_numbers();
        Calculate_Spectrum();
        FFT_spectral_to_physical();
        Export_To_Local_Scalar_Field();
        Assign_Local_Solution_Blocks(member);


    }
    
    template <typename Member_type>
    void Assign_Local_Solution_Blocks(Member_type member);
    void Assign_Local_Solution_Block(Hexa_Block<Soln_pState,Soln_cState> &Solution_Block,
                                     Scalar_Field_Block &Scalar_Field,
                                     double (Soln_pState::*&member));
    
    
    
    void Calculate_Spectrum(void);
    Complex alpha_Rogallo(const double &theta);
    double random_double();
    
    //@}

    
};

template <typename Soln_pState, typename Soln_cState>
Complex SpectralAnalysis<Soln_pState,Soln_cState>::
alpha_Rogallo(const double &theta){
    return ONE*exp(I*theta);
}

// Random number generator
template <typename Soln_pState, typename Soln_cState>
double SpectralAnalysis<Soln_pState,Soln_cState>::
random_double(){
    return  drand48();
}


template <typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState,Soln_cState>::Calculate_Spectrum(void) {
    
    
    
    int index, i, j, l, iconj, jconj, index_conj;
    double theta;
    Complex aa;
 
    for (int ii=0; ii<nK; ii++) {
        
        K[ii].realS = 0;
        K[ii].imagS = 0;
        K[ii].absS  = 0;
        
        for (int jj=0; jj<K[ii].N; jj++) {       // For all grid points with this value of abs_k
            index = K[ii].indexes[jj];      // index corresponding to (i,j,l)
            i = index/(Ny*nz);              // convert index to (i,j,l) coordinates
            j = (index - Ny*nz*i)/nz;
            l = index - nz*(j+Ny*i);
            
            //k1 = k_1(i);
            //k2 = k_2(j);
            //k3 = k_3(l);
            //abs_k = K[ii].k;
            
            theta =  TWO*PI*random_double();   // Random number (0, 2*PI)
            aa = alpha_Rogallo(theta);
            ss[index][0] = real(aa);
            ss[index][1] = imag(aa);
                        
            
            /* ------------ Some exceptions to make FFT real ------------ */
            if ( l==0  ||  l==Nz/2) {
                if ( j>Ny/2  ||  ( i>Nx/2  &&  (j==0 || j==Ny/2) ) ) {
                    iconj = (i==0  ?  0 : Nx-i);
                    jconj = (j==0  ?  0 : Ny-j);
                    index_conj = l + nz*(jconj+Ny*iconj);
                    
                    // complex conjugates
                    ss[index][0] =  ss[index_conj][0];
                    ss[index][1] = -ss[index_conj][1];
                    
                } else if ( (i==0 || i==Nx/2)  &&  (j==0 || j==Ny/2) ) {
                    // real values at 8 corners
                    ss[index][0] = sqrt( sqr(ss[index][0]) + sqr(ss[index][1]) );
                    ss[index][1] = 0.0;
                }
            }     
            K[ii].realS += ss[index][0] / K[ii].N;
            K[ii].imagS += ss[index][1] / K[ii].N;
            K[ii].absS  += sqrt( sqr(ss[index][0]) + sqr(ss[index][1]) ) / K[ii].N;            
        }
        K[ii].rescaling_factor = ONE/K[ii].absS;
    }
    
    
    /* -- rescale -- */
    for (int ii=0; ii<nK; ii++) {
        for (int jj=0; jj<K[ii].N; jj++) {       // For all grid points with this value of abs_k
            index = K[ii].indexes[jj];      // index corresponding to (i,j,l)

            ss[index][0] *= K[ii].rescaling_factor;
            ss[index][1] *= K[ii].rescaling_factor;
                        
            K[ii].realS += ss[index][0] / K[ii].N;
            K[ii].imagS += ss[index][1] / K[ii].N;
            K[ii].absS  += sqrt( sqr(ss[index][0]) + sqr(ss[index][1]) ) / K[ii].N;            
        }
    }
    
    
    
    
    
    
    
}


template <typename Soln_pState, typename Soln_cState>
SpectralAnalysis<Soln_pState,Soln_cState>::
SpectralAnalysis(HexaSolver_Data &Data,
                 HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data) {
    
    Solution_Blocks_ptr  = Solution_Data.Local_Solution_Blocks.Soln_Blks;
    LocalSolnBlkList_ptr = &(Data.Local_Adaptive_Block_List);
    Data_ptr = &Data;
    
    // allocate scalar field
    Global_Scalar_Field.Create_Global(Data.Initial_Mesh,Solution_Data.Input.Grid_IP);
    Local_Scalar_Field.Create_Local(Solution_Blocks_ptr,*LocalSolnBlkList_ptr);

    File_Name = Solution_Data.Input.Output_File_Name;
    
    NCells_Idir = Solution_Data.Input.Grid_IP.NCells_Idir;
    NCells_Jdir = Solution_Data.Input.Grid_IP.NCells_Jdir;
    NCells_Kdir = Solution_Data.Input.Grid_IP.NCells_Kdir;
    
    L1 = Solution_Data.Input.Grid_IP.Box_Width;
    L2 = Solution_Data.Input.Grid_IP.Box_Height;
    L3 = Solution_Data.Input.Grid_IP.Box_Length;
    Ls = min(min(L1,L2),L3);
    
    
    
    //  Allocate_spectrum();
    Nx = NCells_Idir * Data.Initial_Mesh.NBlk_Idir;
    Ny = NCells_Jdir * Data.Initial_Mesh.NBlk_Idir;
    Nz = NCells_Kdir * Data.Initial_Mesh.NBlk_Kdir;
    nz = Nz/2+1;
    
    // Allocation of arrays used in the transforms
    s  = (double *) malloc(Nx*Ny*Nz * sizeof(double));
    ss = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));


}


/**
 * RandomFieldRogallo::Calculate_wave_numbers()
 *
 * Calculates grid characteristic wave numbers
 * and fills the K_container with all possible
 * absolute wave numbers on the grid, saving the
 * index of the each combination.
 */
template<typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState, Soln_cState>::
Calculate_wave_numbers(void) {
    
    
    k0x = TWO*PI/L1;    // smallest wavenumber on grid in x-direction
    k0y = TWO*PI/L2;    // smallest wavenumber on grid in y-direction
    k0z = TWO*PI/L3;    // smallest wavenumber on grid in z-direction
    abs_k0 = min(min(k0x,k0y),k0z);
    kcx = Nx/2*k0x;     // largest wavenumber on grid in x-direction
    kcy = Ny/2*k0y;     // largest wavenumber on grid in y-direction
    kcz = Nz/2*k0z;     // largest wavenumber on grid in z-direction
    abs_kc = sqrt( sqr(kcx) + sqr(kcy) + sqr(kcz) );
    
    K = new Spectral_container [Nx*Ny*nz];  // Container of information for a wavenumber
    Spectral_BinaryTree Ktree;              // A binary tree that will order and store K_containers
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


template<typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState, Soln_cState>::
Average_Spectrum(void) {

    int index = 0;
    for (int ii=0; ii <nK; ii++) {
        K[ii].realS = 0;
        K[ii].imagS = 0;
        K[ii].absS = 0;
        for(int jj=0; jj<K[ii].N; jj++) {
            index = K[ii].indexes[jj];
            K[ii].realS += ss[index][0];
            K[ii].imagS += ss[index][1];
            K[ii].absS += sqrt( sqr(ss[index][0]) + sqr(ss[index][1]) );

        }
        K[ii].realS /= K[ii].N;
        K[ii].imagS /= K[ii].N;
        K[ii].absS /= K[ii].N;

    } 

}



/**
 * RandomFieldRogallo::FFT_spectral_to_physical
 *
 * The physical velocity components are calculated
 * through FFT of the spectral velocity components.
 */
template<typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState, Soln_cState>::
FFT_spectral_to_physical(void) {
    
    fftw_plan physical;
    physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, ss, s, FFTW_ESTIMATE);
    fftw_execute(physical); 
    fftw_destroy_plan(physical);
    
    /* ---------- FFT scaling factor ---------- */
    double FFT_scaling_factor = sqrt(pow(TWO*PI,3.0)/(L1*L2*L3));
    int index;
    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny; j++) {
            for (int l=0; l<Nz; l++) {
                index = l + Nz*(j+Ny*i);
                s[index] *= FFT_scaling_factor;
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
template<typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState, Soln_cState>::
FFT_physical_to_spectral(void) {
    
    fftw_plan spectral;
    spectral = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, s, ss, FFTW_ESTIMATE);
    fftw_execute(spectral); 
    fftw_destroy_plan(spectral);
        
    /* --------------- FFT scaling factor ---------------- */
    
    double FFT_scaling_factor = ONE/(Nx*Ny*Nz)   / sqrt(pow(TWO*PI,3)/(L1*L2*L3));
    int index;
    for (int ii=0; ii<nK; ii++) {
        for (int jj=0; jj<K[ii].N; jj++) {
            index = K[ii].indexes[jj];
            ss[index][0] *= FFT_scaling_factor;
            ss[index][1] *= FFT_scaling_factor;
        }
    }
};




/**
 * Open_Spectrum_File
 */
template<typename Soln_pState, typename Soln_cState>
int SpectralAnalysis<Soln_pState, Soln_cState>::
Open_Spectrum_File(char *suffix) {
    
    scalar_name = suffix;
    
    int i;
    char prefix[256], extension[256], 
    spectrum_file_name[256], gnuplot_file_name[256];
    char *spectrum_file_name_ptr, *gnuplot_file_name_ptr;
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
    strcat(prefix, "_");
    strcat(prefix, suffix);
    strcat(prefix, "_spectrum");
    
    strcpy(extension, ".dat");
    strcpy(spectrum_file_name, prefix);
    strcat(spectrum_file_name, extension);
    
    spectrum_file_name_ptr = spectrum_file_name;
    
    /* Open the turbulence spectrum file. */
    
    Spectrum_File.open(spectrum_file_name_ptr, ios::out);
    if (Spectrum_File.bad()) return (1);
    
    /* Write the appropriate GNUPLOT command file for 
     plotting turbulence progress file information. */
    
    strcpy(extension, ".gplt");
    strcpy(gnuplot_file_name, prefix);
    strcat(gnuplot_file_name, extension);
    
    gnuplot_file_name_ptr = gnuplot_file_name;
    
    gnuplot_file.open(gnuplot_file_name_ptr, ios::out);
    if (gnuplot_file.fail()) return(1);
    
    gnuplot_file 
    << "set title \"Spectrum of "<< scalar_name <<"\"\n"
    << "set xlabel \"wavenumber k \"\n"
    << "set ylabel \"Fourier Transform\"\n" 
    << "set grid \n"
    //<< "set logscale xy\n"
    << "plot \"" << spectrum_file_name_ptr << "\" using 1:2 \\\n"
    << "     title \"" << "real part"    << "\" with lines , \\\n"
    << "\"" << spectrum_file_name_ptr << "\" using 1:3 \\\n"
    << "     title \"" << "imaginary part"    << "\" with lines , \\\n"
    << "\"" << spectrum_file_name_ptr << "\" using 1:4 \\\n"
    << "     title \"" << "absolute"    << "\" with lines \n"
    << "pause -1  \"Hit return to continue\"\n";
    
    gnuplot_file.close();
    
    /* Preparation of progress file complete.
     Return zero value. */    
    return(0);
}


/**
 * Output_Spectrum_to_File
 */
template<typename Soln_pState, typename Soln_cState>
int SpectralAnalysis<Soln_pState, Soln_cState>::
Output_Spectrum_to_File() {
    
    Spectrum_File << setprecision(6);
    Spectrum_File.setf(ios::scientific);
    
    for (int ii=0; ii <nK; ii++) {
        Spectrum_File    << K[ii].k  << " " << K[ii].realS 
                                     << " " << K[ii].imagS
                                     << " " << K[ii].absS
                         << "\n";
    } 
    
    Spectrum_File.unsetf(ios::scientific);
    Spectrum_File.flush();
    
    return(0);
}

/**
 * Close_Spectrum_File
 */
template<typename Soln_pState, typename Soln_cState>
int SpectralAnalysis<Soln_pState, Soln_cState>::
Close_Spectrum_File( ) {
    Spectrum_File.close();
    return(0);
}
















template <typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState,Soln_cState>::Import_From_Global_Scalar_Field(void) {
    
    Grid3D_Hexa_Multi_Block_List *Mesh_ptr = &Data_ptr->Initial_Mesh;
    
    int index;
    int nBlk, ix, iy, iz;
    for (int kBlk = 0; kBlk <= Mesh_ptr->NBlk_Kdir-1; ++kBlk) {
        for (int jBlk = 0; jBlk <= Mesh_ptr->NBlk_Jdir-1; ++jBlk) {
            for (int iBlk = 0; iBlk <= Mesh_ptr->NBlk_Idir-1; ++iBlk) {
                nBlk = iBlk + 
                jBlk*Mesh_ptr->NBlk_Idir + 
                kBlk*Mesh_ptr->NBlk_Idir*Mesh_ptr->NBlk_Jdir;
                for (int i = Mesh_ptr->Grid_Blks[nBlk].ICl; i <= Mesh_ptr->Grid_Blks[nBlk].ICu; ++i) {
                    for (int j = Mesh_ptr->Grid_Blks[nBlk].JCl; j <= Mesh_ptr->Grid_Blks[nBlk].JCu; ++j) {
                        for (int k = Mesh_ptr->Grid_Blks[nBlk].KCl; k <= Mesh_ptr->Grid_Blks[nBlk].KCu; ++k) {
                            ix = iBlk*NCells_Idir+(i-Mesh_ptr->Grid_Blks[nBlk].ICl);
                            iy = jBlk*NCells_Jdir+(j-Mesh_ptr->Grid_Blks[nBlk].JCl);
                            iz = kBlk*NCells_Kdir+(k-Mesh_ptr->Grid_Blks[nBlk].KCl);
                            index = iz + Nz*(iy + ix*Ny);
                            s[index] = Global_Scalar_Field.Scalar_Field_Blks[nBlk].s[i][j][k];
                        }
                    }
                }
            }
        }
    } 
}


template <typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState,Soln_cState>::Export_To_Local_Scalar_Field(void) {
    
    Grid3D_Hexa_Multi_Block_List *Mesh_ptr = &Data_ptr->Initial_Mesh;

    int index;
    int nBlk, ix, iy, iz;
    int INl,INu,JNl,JNu,KNl,KNu;
    for (int kBlk = 0; kBlk <= Mesh_ptr->NBlk_Kdir-1; ++kBlk) {
        for (int jBlk = 0; jBlk <= Mesh_ptr->NBlk_Jdir-1; ++jBlk) {
            for (int iBlk = 0; iBlk <= Mesh_ptr->NBlk_Idir-1; ++iBlk) {
                nBlk = iBlk + 
                jBlk*Mesh_ptr->NBlk_Idir + 
                kBlk*Mesh_ptr->NBlk_Idir*Mesh_ptr->NBlk_Jdir;
                for (int i = Mesh_ptr->Grid_Blks[nBlk].ICl; i <= Mesh_ptr->Grid_Blks[nBlk].ICu; ++i) {
                    for (int j = Mesh_ptr->Grid_Blks[nBlk].JCl; j <= Mesh_ptr->Grid_Blks[nBlk].JCu; ++j) {
                        for (int k = Mesh_ptr->Grid_Blks[nBlk].KCl; k <= Mesh_ptr->Grid_Blks[nBlk].KCu; ++k) {
                            ix = iBlk*NCells_Idir+(i-Mesh_ptr->Grid_Blks[nBlk].ICl);
                            iy = jBlk*NCells_Jdir+(j-Mesh_ptr->Grid_Blks[nBlk].JCl);
                            iz = kBlk*NCells_Kdir+(k-Mesh_ptr->Grid_Blks[nBlk].KCl);
                            index = iz + Nz*(iy + ix*Ny);
                            Local_Scalar_Field.Scalar_Field_Blks[nBlk].s[i][j][k] = s[index];                                 
                        }
                    }
                }
            }
        } 
    } 
}



/**
 * Scalar_Field_Multi_Block_List::Collect_Blocks_from_all_processors
 *
 * Collects each local Scalar_Field_Block from all the processors
 * and copies them in this (global) Scalar_Field_Multi_Block_List by
 * comparing the global block number.
 *
 * \param [in] OcTree                  The Octree_DataStructure
 * \param [in] LocalSolnBlockList      The local AdaptiveBlock3D_List
 * \param [in] Global_Soln_Block_List  The global AdaptiveBlock3D_ResourceList
 */
template <typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState,Soln_cState>::Gather_from_all_processors(void) {
    
    /* First copy the Local_Scalar_Field in the Global_Scalar_Field */
    for (int nBlk = 0; nBlk < Local_Scalar_Field.NBlk; nBlk++) {
        Global_Scalar_Field.Scalar_Field_Blks[Local_Scalar_Field.Scalar_Field_Blks[nBlk].gblknum].Copy(Local_Scalar_Field.Scalar_Field_Blks[nBlk]);
    }
    
    
    
    
    /* Then gather the rest of the Global_Scalar_Field */
    
    int *new_blocks_CPU, *CPUs_in_new_blocks;
    int  new_blocks_BLK;
    CPUs_in_new_blocks = new int[Data_ptr->Global_Adaptive_Block_List.Ncpu];
    new_blocks_CPU = new int[Data_ptr->Global_Adaptive_Block_List.Ncpu];
    
#ifdef _MPI_VERSION
    MPI::Intracomm new_comm;
    MPI::Group     big_group = MPI::COMM_WORLD.Get_group();
    MPI::Group     new_group;
#endif
    
    Scalar_Field_Block  Scalar_Field_Block_duplicated;
    OctreeBlock *octree_block_duplicated_ptr;
    
    for (int iCPU = 0; iCPU <= Data_ptr->Octree.Ncpu-1; ++iCPU ) { // Loop over available processors.
        for (int iBLK = 0; iBLK <= Data_ptr->Octree.Nblk-1; ++iBLK ) { // Loop over available blocks.
            if (Data_ptr->Octree.Blocks[iCPU][iBLK] != NULL) {
                if (Data_ptr->Octree.Blocks[iCPU][iBLK]->block.used) { // Check if the solution block is used.
                    
                    octree_block_duplicated_ptr = Data_ptr->Octree.Blocks[iCPU][iBLK];
                    new_blocks_CPU[0] = octree_block_duplicated_ptr->block.info.cpu;
                    new_blocks_BLK = octree_block_duplicated_ptr->block.info.gblknum;
                    
                    CPUs_in_new_blocks[0] = new_blocks_CPU[0];
                    
                    if (Data_ptr->Global_Adaptive_Block_List.Ncpu > 1) {
                        for (int iNEW = 1 ; iNEW < Data_ptr->Global_Adaptive_Block_List.Ncpu; ++iNEW ) {
                            if (iNEW != new_blocks_CPU[0]) {
                                new_blocks_CPU[iNEW] = iNEW;
                            } else {
                                new_blocks_CPU[iNEW] = 0;
                            }
                            CPUs_in_new_blocks[iNEW] = new_blocks_CPU[iNEW];
                        }
                    }
                    
                    if (Data_ptr->Local_Adaptive_Block_List.ThisCPU == new_blocks_CPU[0]) { 
                        Scalar_Field_Block_duplicated.Copy(Global_Scalar_Field.Scalar_Field_Blks[new_blocks_BLK]); 
                    } 
#ifdef _MPI_VERSION 
                    new_group = big_group.Incl(Data_ptr->Global_Adaptive_Block_List.Ncpu, CPUs_in_new_blocks);
                    new_comm = MPI::COMM_WORLD.Create(new_group);
                    
                    Scalar_Field_Block_duplicated.Broadcast(new_comm, new_blocks_CPU[0]); 
                    
                    if (new_comm != MPI::COMM_NULL) new_comm.Free();
                    new_group.Free();
#endif
                    Global_Scalar_Field.Scalar_Field_Blks[Scalar_Field_Block_duplicated.gblknum].Copy(Scalar_Field_Block_duplicated);
                    
                }
            }
        }
    }
    
    //delete octree_block_duplicated_ptr;
    delete[] CPUs_in_new_blocks;
    delete[] new_blocks_CPU;
    
    
}




/**
 * Get_Local_Homogeneous_Turbulence_Scalar_Field
 *
 * Fills all the local velocity field blocks with velocity fluctuations around
 * a domain averaged velocity with data from the local solution blocks
 *
 * \param [in]  Solution_Block              All the local Solution Hexa_Blocks
 * \param [in]  LocalSolnBlockList          The AdaptiveBlock3D_List
 * \param [in]  v_average                   The average velocity of the whole domain
 * \param [out] Local_Scalar_Field_List   The list of local Scalar_Field_Blocks
 */
template<typename Soln_pState, typename Soln_cState>
template<typename Function_type>
void SpectralAnalysis<Soln_pState,Soln_cState>::Assign_Local_Scalar_Field(Function_type func) {
    /* Assign initial turbulent velocity field to each solution block. */
    for (int nBlk = 0 ; nBlk < LocalSolnBlkList_ptr->Nblk ; nBlk++) {
        if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
            if (Local_Scalar_Field.Scalar_Field_Blks[nBlk].Allocated) {
                Assign_Local_Scalar_Field_Block(Solution_Blocks_ptr[nBlk],
                                                Local_Scalar_Field.Scalar_Field_Blks[nBlk],
                                                func);
            } 
        } 
    } 
}


///**
// * Get_Local_Homogeneous_Turbulence_Scalar_Field
// *
// * Copies velocity fluctuations from a local solution block to a
// * local velocity field block.
// *
// * \param [in]  Solution_Block  A local Solution Hexa_Block
// * \param [in]  v_average       The average velocity of the whole domain
// * \param [out] Scalar_Field  A local Scalar_Field_Block
// */
//template<typename Soln_pState, typename Soln_cState>
//template<typename Function_type>
//void SpectralAnalysis<Soln_pState,Soln_cState>::Assign_Local_Scalar_Field_Block(Hexa_Block<Soln_pState,Soln_cState> &Solution_Block,
//                                                                                Scalar_Field_Block &Scalar_Field,
//                                                                                Function_type func) {
//    
//    /* Assign initial turbulent velocity field to solution block. */
//    
//    for (int i = Scalar_Field.ICl ; i <= Scalar_Field.ICu ; i++) {
//        for (int j = Scalar_Field.JCl ; j <= Scalar_Field.JCu ; j++) {
//            for (int k = Scalar_Field.KCl ; k <= Scalar_Field.KCu ; k++) {
//                Scalar_Field.s[i][j][k] = (Solution_Block.W[i][j][k].*func)();
//            } 
//        }
//    }
//}


template<typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState,Soln_cState>::Assign_Local_Scalar_Field_Block(Hexa_Block<Soln_pState,Soln_cState> &Solution_Block,
                                                                                Scalar_Field_Block &Scalar_Field,
                                                                                double (Soln_pState::*&member)) {
    
    /* Assign initial turbulent velocity field to solution block. */
    
    for (int i = Scalar_Field.ICl ; i <= Scalar_Field.ICu ; i++) {
        for (int j = Scalar_Field.JCl ; j <= Scalar_Field.JCu ; j++) {
            for (int k = Scalar_Field.KCl ; k <= Scalar_Field.KCu ; k++) {
                Scalar_Field.s[i][j][k] = Solution_Block.W[i][j][k].*member;
            } 
        }
    }
}


template<typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState,Soln_cState>::Assign_Local_Scalar_Field_Block(Hexa_Block<Soln_pState,Soln_cState> &Solution_Block,
                                                                                Scalar_Field_Block &Scalar_Field,
                                                                                double (Soln_pState::*&func)(void)) {
    
    /* Assign initial turbulent velocity field to solution block. */
    
    for (int i = Scalar_Field.ICl ; i <= Scalar_Field.ICu ; i++) {
        for (int j = Scalar_Field.JCl ; j <= Scalar_Field.JCu ; j++) {
            for (int k = Scalar_Field.KCl ; k <= Scalar_Field.KCu ; k++) {
                Scalar_Field.s[i][j][k] = (Solution_Block.W[i][j][k].*func)();
            } 
        }
    }
}

template<typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState,Soln_cState>::Assign_Local_Scalar_Field_Block(Hexa_Block<Soln_pState,Soln_cState> &Solution_Block,
                                                                                Scalar_Field_Block &Scalar_Field,
                                                                                double (Soln_pState::*&func)(const Soln_pState&, const Soln_pState&, const Soln_pState&)) {
    
    /* Assign initial turbulent velocity field to solution block. */
    
    for (int i = Scalar_Field.ICl ; i <= Scalar_Field.ICu ; i++) {
        for (int j = Scalar_Field.JCl ; j <= Scalar_Field.JCu ; j++) {
            for (int k = Scalar_Field.KCl ; k <= Scalar_Field.KCu ; k++) {
                Scalar_Field.s[i][j][k] = (Solution_Block.W[i][j][k].*func)(Solution_Block.dWdx[i][j][k],Solution_Block.dWdy[i][j][k],Solution_Block.dWdz[i][j][k]);
            } 
        }
    }
}

template<typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState,Soln_cState>::Assign_Local_Scalar_Field_Block(Hexa_Block<Soln_pState,Soln_cState> &Solution_Block,
                                                                                Scalar_Field_Block &Scalar_Field,
                                                                                double (Soln_pState::*&func)(const Soln_pState&, const Soln_pState&, const Soln_pState&, const double&)) {
    
    /* Assign initial turbulent velocity field to solution block. */
    
    for (int i = Scalar_Field.ICl ; i <= Scalar_Field.ICu ; i++) {
        for (int j = Scalar_Field.JCl ; j <= Scalar_Field.JCu ; j++) {
            for (int k = Scalar_Field.KCl ; k <= Scalar_Field.KCu ; k++) {
                Scalar_Field.s[i][j][k] = (Solution_Block.W[i][j][k].*func)(Solution_Block.dWdx[i][j][k],Solution_Block.dWdy[i][j][k],Solution_Block.dWdz[i][j][k],Solution_Block.Grid.Cell[i][j][k].V);
            } 
        }
    }
}







template<typename Soln_pState, typename Soln_cState>
template<typename Member_type>
void SpectralAnalysis<Soln_pState,Soln_cState>::Assign_Local_Solution_Blocks(Member_type member) {
        
    for (int nBlk = 0 ; nBlk<LocalSolnBlkList_ptr->Nblk ; nBlk++) {
        if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
            if (Local_Scalar_Field.Scalar_Field_Blks[LocalSolnBlkList_ptr->Block[nBlk].info.gblknum].Allocated) {
                Assign_Local_Solution_Block(Solution_Blocks_ptr[nBlk],
                                            Local_Scalar_Field.Scalar_Field_Blks[LocalSolnBlkList_ptr->Block[nBlk].info.gblknum],
                                            member);
            } 
        } 
    } 
}


template<typename Soln_pState, typename Soln_cState>
void SpectralAnalysis<Soln_pState,Soln_cState>::Assign_Local_Solution_Block(Hexa_Block<Soln_pState,Soln_cState> &Solution_Block,
                                                                            Scalar_Field_Block &Scalar_Field,
                                                                            double (Soln_pState::*&member)) {
    
    for (int i = Scalar_Field.ICl ; i <= Scalar_Field.ICu ; i++) {
        for (int j = Scalar_Field.JCl ; j <= Scalar_Field.JCu ; j++) {
            for (int k = Scalar_Field.KCl ; k <= Scalar_Field.KCu ; k++) {
                Solution_Block.W[i][j][k].*member = Scalar_Field.s[i][j][k];
            } 
        }
    }
}







///**
// * Assign_Homogeneous_Turbulence_Scalar_Field
// *
// * Superimposes turbulence velocity fluctuations on all the local solution blocks
// * with velocities from the global Scalar_Field_Multi_Block_List
// *
// * \param [in]  Solution_Block        The local Solution Hexa_Blocks
// * \param [in]  LocalSolnBlockList    The local AdaptiveBlock3D_List
// * \param [out] Scalar_Field        The list of all Scalar_Field_Blocks
// */
//template<typename HEXA_BLOCK>
//void Assign_Scalar_Field(HEXA_BLOCK *Solution_Block,
//                         const AdaptiveBlock3D_List &LocalSolnBlockList,
//                         const Scalar_Field_Multi_Block_List &Scalar_Field) {
//    
//    /* Assign initial turbulent velocity field to each solution block. */
//    
//    for (int nBlk = 0 ; nBlk <= LocalSolnBlockList.Nblk-1 ; nBlk++) {
//        if (LocalSolnBlockList.Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
//            if (Scalar_Field.Scalar_Field_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum].Allocated) {
//                Assign_Scalar_Field(Solution_Block[nBlk],
//                                    Scalar_Field.Scalar_Field_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum]);
//            } /* endif */
//        } /* endif */
//    }  /* endfor */
//    
//}
//
///**
// * Assign_Homogeneous_Turbulence_Scalar_Field
// *
// * Superimposes velocities from a Scalar_Field_Block to
// * a local solution block.
// *
// * \param [out]  Solution_Block   A local Solution Hexa_Block
// * \param [in]   Scalar_Field   A Scalar_Field_Block
// */
//template<typename HEXA_BLOCK>
//void Assign_Scalar_Field(HEXA_BLOCK &Solution_Block,
//                         const Scalar_Field_Block &Scalar_Field) {
//    
//    /* Assign initial turbulent velocity field to solution block. */
//    
//    for (int i = Solution_Block.ICl ; i <= Solution_Block.ICu ; i++) {
//        for (int j = Solution_Block.JCl ; j <= Solution_Block.JCu ; j++) {
//            for (int k = Solution_Block.KCl ; k <= Solution_Block.KCu ; k++) {
//                Solution_Block.W[i][j][k].rho = Scalar_Field.s[i][j][k];
//            } /* endfor */
//        } /* endfor */
//    } /* endfor */
//}
//
///**
// * Assign_Scalar_Field
// *
// * Sets velocities from all local solution block to domain averaged velocity,
// * superimposed with turbulence velocity fluctuations from the global
// * Scalar_Field_Multi_Block_List.
// *
// * \param [in]  Solution_Block        The local Solution Hexa_Blocks
// * \param [in]  LocalSolnBlockList    The local AdaptiveBlock3D_List
// * \param [in]  v_average             The domain averaged velocity
// * \param [out] Scalar_Field        The list of all Scalar_Field_Blocks
// */
//template<typename HEXA_BLOCK>
//void Assign_Scalar_Field(HEXA_BLOCK *Solution_Block,
//                                                  const AdaptiveBlock3D_List &LocalSolnBlockList,
//                                                  const Vector3D &v_average,
//                                                  const Scalar_Field_Multi_Block_List &Scalar_Field) {
//    
//    /* Assign initial turbulent velocity field to each solution block. */
//    
//    for (int nBlk = 0 ; nBlk <= LocalSolnBlockList.Nblk-1 ; nBlk++) {
//        if (LocalSolnBlockList.Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
//            if (Scalar_Field.Scalar_Field_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum].Allocated) {
//                Assign_Homogeneous_Turbulence_Scalar_Field(Solution_Block[nBlk],
//                                                             Scalar_Field.Scalar_Field_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum],
//                                                             v_average);
//            } /* endif */
//        } /* endif */
//    }  /* endfor */
//    
//}
//
///**
// * Assign_Scalar_Field
// *
// * Set velocities from a local solution block to the velocities of 
// * a Scalar_Field_Block superimposed with a given average velocity.
// *
// * \param [out]  Solution_Block   A local Solution Hexa_Block
// * \param [in]   Scalar_Field   A Scalar_Field_Block
// * \param [in]   v_average        The domain averaged velocity
// */
//template<typename HEXA_BLOCK>
//void Assign_Scalar_Field(HEXA_BLOCK &Solution_Block,
//                         const Scalar_Field_Block &Scalar_Field,
//                         const Vector3D &v_average) {
//    
//    /* Assign initial turbulent velocity field to solution block. */
//    
//    for (int i = Solution_Block.ICl ; i <= Solution_Block.ICu ; i++) {
//        for (int j = Solution_Block.JCl ; j <= Solution_Block.JCu ; j++) {
//            for (int k = Solution_Block.KCl ; k <= Solution_Block.KCu ; k++) {
//                Solution_Block.W[i][j][k].rho = Scalar_Field.s[i][j][k];
//                Solution_Block.U[i][j][k] = Solution_Block.W[i][j][k].U();
//            } /* endfor */
//        } /* endfor */
//    } /* endfor */
//}

















//
///************************************************************************//**
// * \class _Member_Function_Wrapper_
// * \brief Adaptor for making a class member function look like an ordinary function
// *
// * This wrapper is useful for passing member functions to subroutines that 
// * take ordinary functions (e.g. numerical integration subroutines)
// ****************************************************************************/
//template <typename Soln_pState, typename member_type>
//class Member_of_Member{
//    
//public:
//    
//    /* constructor */
//    Member_of_Member(Vector3D (Soln_pState::*V_ptr),
//                     member_type V_mem_ptr): Vector_ptr(V_ptr), Vector_member_ptr(V_mem_ptr) {
//     
//        resulting_ptr = (Vector_ptr->*Vector_member_ptr);
//    }
//    
//    double (Soln_pState::*resulting_ptr);
//        
//private:
//    Vector3D (Soln_pState::*Vector_ptr);
//    //double (Vector3D::*Vector_member_ptr);
//    member_type Vector_member_ptr;
//
//};
//
