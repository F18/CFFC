/*
 *  SpectralAnalysis.cc
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/04/08.
 *
 */

#include "SpectralAnalysis.h"

/* -------------------------------------------------------------------------- */
/*                      Scalar_Field_Block member functions                   */
/* -------------------------------------------------------------------------- */



/**
 * Scalar_Field_Block::allocate -- Allocate memory.    
 */
void Scalar_Field_Block::allocate(const int Ni, 
                                  const int Nj, 
                                  const int Nk,
                                  const int Ng) {
    assert( Ni >= 1 && Nj >= 1 && Nk >= 1 && Ng >=1 && !Allocated);
    NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
    NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1;
    NCk = Nk+2*Ng; KCl = Ng; KCu = Nk+Ng-1;
    Nghost = Ng;
    Allocated = SCALAR_FIELD_DATA_USED;
    
    s = new double**[NCi];
    for (int i = 0; i <= NCi-1; ++i ){
        s[i] = new double*[NCj];
        for (int j = 0; j <= NCj-1; ++j ){
            s[i][j] = new double[NCk];
        } /* endfor */
    } /* endfor */
    
}


/**
 * Scalar_Field_Block::deallocate -- Deallocate memory.     
 */
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
        Allocated = SCALAR_FIELD_DATA_NOT_USED;
    } /* endif */
}



/**
 * Scalar_Field_Block::Copy -- Copy a block.  
 */
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
















/**
 * Scalar_Field_Block::Broadcast                       
 * Broadcast turbulent velocity field block. 
 */
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


/**
 * Scalar_Field_Block::Broadcast                        
 * Broadcast turbulent velocity field block.  
 */
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



/**
 * Scalar_Field_Multi_Block_List::Allocate
 * Allocate memory.  
 */
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

/**
 * Scalar_Field_Multi_Block_List::Allocate 
 * Allocate memory.  
 */
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

/**
 * Scalar_Field_Multi_Block_List::Deallocate 
 * Deallocate memory.
 */
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
 * Scalar_Field_Multi_Block_List::Create_Global
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
 * Scalar_Field_Multi_Block_List::Broadcast
 * Broadcast the turbulent velocity field data from the main processor. 
 */
inline void Scalar_Field_Multi_Block_List::Broadcast() {
    
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