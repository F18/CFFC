/* AdaptiveBlock.h:  Header file defining adaptive block classes. */

#ifndef _ADAPTIVEBLOCK_INCLUDED
#define _ADAPTIVEBLOCK_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include math macro, CFD, and MPI header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

/* Define the block in-use indicators. */

#define	ADAPTIVEBLOCK2D_USED                            1
#define	ADAPTIVEBLOCK2D_NOT_USED                        0

/* Define the maximum number of possible neighbouring
   adaptive blocks. */

#define	ADAPTIVEBLOCK2D_NUMBER_NEIGHBOURS_MAX           2
#define	ADAPTIVEBLOCK2D_NUMBER_CORNER_NEIGHBOURS_MAX    1

/* Define the block sector indicators for 
   block refinement and coarsening. */

#define	ADAPTIVEBLOCK2D_SECTOR_NONE                    -1
#define	ADAPTIVEBLOCK2D_SECTOR_SW                       0
#define	ADAPTIVEBLOCK2D_SECTOR_SE                       1
#define	ADAPTIVEBLOCK2D_SECTOR_NW                       2
#define	ADAPTIVEBLOCK2D_SECTOR_NE                       3

/* Define the block refinement indicators. */

#define	ADAPTIVEBLOCK2D_REFINE                          1
#define	ADAPTIVEBLOCK2D_NOCHANGE                        0
#define	ADAPTIVEBLOCK2D_COARSEN                        -1

/* Define the classes. */

/********************************************************
 * Class: AdaptiveBlockResourceList                     *
 *                                                      *
 * Member functions                                     *
 *      ThisCPU   -- Return the global MPI number or    *
 *                   rank for the CPU that is executing *
 *                   the task or process.               *
 *      Ncpu      -- Return number of available CPUs.   *
 *      Nblk      -- Return number of blocks per CPU.   *
 *      Ntotal    -- Return total number of available   *
 *                   adaptive blocks for computation.   *
 *      Nused     -- Return number of adaptive blocks   *
 *                   currently in use in computation.   *
 *      Nfree     -- Return number of free or open      *
 *                   blocks which may be added to       *
 *                   computation.                       *
 *      CPU       -- Return array of CPU numbers.       *
 *      Block     -- Return array of local block        *
 *                   numbers.                           *
 *      allocate  -- Allocate memory for resource list. *
 *      deallocate -- Deallocate memory for resource    *
 *                    list.                             *
 *      initialize -- Initialize resource list.         *
 *      nextCPU    -- Return CPU number of next         *
 *                    available block.                  *
 *      nextBlock  -- Return Block number of next       *
 *                    available block.                  *
 *      update_next -- Use next available block,        *
 *                     updating all counters.           *
 *      returnCPU -- Return CPU number to list of       *
 *                   available blocks.                  *
 *      returnBlock -- Return Block number to list of   *
 *                     available blocks.                *
 *      update_return -- Update the returned block, use *
 *                       as next available block,       *
 *                       updating all counters.         *
 *                                                      *
 * Member operators                                     *
 *      R -- adaptive block resource list               *
 *                                                      *
 * R = R;                                               *
 * cout << R; (output function)                         *
 * cin  >> R; (input function)                          *
 *                                                      *
 ********************************************************/
class AdaptiveBlockResourceList{
  private:
  public:
    int    ThisCPU;  // Number or rank of CPU executing the
                     // current process.
    int       Ncpu;  // Number of available processors.
    int       Nblk;  // Number of blocks per processor.
    int     Ntotal;  // Total number of available adaptive blocks.
    int      Nused;  // Number of blocks currently being used.
    int      Nfree;  // Number of blocks that are currently free.
    int       *CPU;  // Array of CPU numbers.
    int     *Block;  // Array of local block numbers.
                     // Made public so can access them.

    /* Creation, copy, and assignment constructors. */
    AdaptiveBlockResourceList(void) {
       ThisCPU = 0; Ncpu = 0; Nblk = 0; Ntotal = 0; Nused = 0; Nfree = 0;
       CPU = NULL; Block = NULL; 
    }

    AdaptiveBlockResourceList(const AdaptiveBlockResourceList &R) {
       ThisCPU = R.ThisCPU; Ncpu = R.Ncpu; Nblk = R.Nblk; Ntotal = R.Ntotal; 
       Nused = R.Nused; Nfree = R.Nfree; CPU = R.CPU; Block = R.Block;
    }

    /* Destructor. */
    // ~AdaptiveBlockResourceList(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // AdaptiveBlockResourceList operator = (const AdaptiveBlockResourceList &R);
    // Use automatically generated assignment operator.

    /* Allocate memory for the resource list. */
    void allocate(const int ncpu, const int nblk);

    /* Deallocate memory for the resource list. */
    void deallocate(void);

    /* Initialize resource list. */
    void initialize(void);

    /* Copy resource list. */
    void copy(const AdaptiveBlockResourceList &R);

    /* CPU number of next available block. */
    int nextCPU(void);

    /* Block number of next available block. */
    int nextBlock(void);

    /* Update next available block. */
    void update_next(void);

    /* Return CPU number to list of available blocks. */
    void returnCPU(const int cpu);

    /* Return Block number to list of available blocks. */
    void returnBlock(const int block);

    /* Update returned block. */
    void update_return(void);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, 
                                 const AdaptiveBlockResourceList &R);
    friend istream &operator >> (istream &in_file, 
                                 AdaptiveBlockResourceList &R);

};

/*******************************************************************
 * AdaptiveBlockResourceList::allocate -- Allocate memory.         *
 *******************************************************************/
inline void AdaptiveBlockResourceList::allocate(const int ncpu, 
                                                const int nblk) {
   assert( ncpu > 0 && nblk > 0); 
   Ncpu = ncpu; Nblk = nblk; Ntotal = ncpu*nblk; Nused = 0; Nfree = Ntotal;
   CPU = new int[Ntotal]; Block = new int[Ntotal];
}

/*******************************************************************
 * AdaptiveBlockResourceList::deallocate -- Deallocate memory.     *
 *******************************************************************/
inline void AdaptiveBlockResourceList::deallocate(void) {
   Ncpu = 0; Nblk = 0; Ntotal = 0; Nused = 0; Nfree = 0;
   delete []CPU; CPU = NULL; delete []Block; Block = NULL; 
}

/*******************************************************************
 * AdaptiveBlockResourceList::initialize -- Initialization.        *
 *******************************************************************/
inline void AdaptiveBlockResourceList::initialize(void) {
  int i, j;
  Nused = 0;
  for ( j = 0; j <= Nblk-1; ++j ) {
      for ( i = 0; i <= Ncpu-1; ++i ) {
	  CPU[Nused] = i;
          Block[Nused] = j;
          ++Nused;
      } /* endfor */
  } /* endfor */
  Nused = 0;
}

/*******************************************************************
 * AdaptiveBlockResourceList::initialize -- Copy.                  *
 *******************************************************************/
inline void AdaptiveBlockResourceList::copy(const AdaptiveBlockResourceList &R) {
  int nused;
  if (CPU==NULL && Block==NULL) deallocate();
  allocate(R.Ncpu,R.Nblk);
  ThisCPU = R.ThisCPU;
  Nused = R.Nused;
  Nfree = R.Nfree;
  nused = 0;
  for (int j = 0; j < Nblk; j++) {
    for (int i = 0; i < Ncpu; i++) {
      CPU[nused] = R.CPU[nused];
      Block[nused] = R.Block[nused];
      nused++;
    }
  }
}

/*******************************************************************
 * AdaptiveBlockResourceList::nextCPU -- Next CPU number.          *
 *******************************************************************/
inline int AdaptiveBlockResourceList::nextCPU(void) {
   return(CPU[Nused]);
}

/*******************************************************************
 * AdaptiveBlockResourceList::nextBlock -- Next Block number.      *
 *******************************************************************/
inline int AdaptiveBlockResourceList::nextBlock(void) {
   return(Block[Nused]);
}

/*******************************************************************
 * AdaptiveBlockResourceList::update_next -- Update next.          *
 *******************************************************************/
inline void AdaptiveBlockResourceList::update_next(void) {
   Nused += 1; Nfree -= 1;
}

/*******************************************************************
 * AdaptiveBlockResourceList::returnCPU -- Return CPU number.      *
 *******************************************************************/
inline void AdaptiveBlockResourceList::returnCPU(const int cpu) {
   CPU[Nused-1] = cpu;
}

/*******************************************************************
 * AdaptiveBlockResourceList::returnBlock -- Return Block number.  *
 *******************************************************************/
inline void AdaptiveBlockResourceList::returnBlock(const int block) {
   Block[Nused-1] = block;
}

/**********************************************************************
 * AdaptiveBlockResourceList::update_return -- Update returned block. *
 **********************************************************************/
inline void AdaptiveBlockResourceList::update_return(void) {
   Nused -= 1; Nfree += 1;
}

/*******************************************************************
 * AdaptiveBlockResourceList -- Input-output operators.            *
 *******************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const AdaptiveBlockResourceList &R) {
  int i;
  for ( i = 0; i <= R.Ntotal-1; ++i ) {
      out_file << " " << R.CPU[i] << " " << R.Block[i] << "\n";
  } /* endfor */
  return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             AdaptiveBlockResourceList &R) {
  int i;
  for ( i = 0; i <= R.Ntotal-1; ++i ) {
      in_file >> R.CPU[i] >> R.Block[i];
  } /* endfor */
  return (in_file);
}

/*************************************************************
 * AdaptiveBlockResourceList -- External subroutines.        *
 *************************************************************/

extern void Create_Block_Resource_List(AdaptiveBlockResourceList &List_of_Available_Blocks,
                                       const int Number_of_Processors,
                                       const int Number_of_Blocks_per_Processor);

/********************************************************
 * Class: AdaptiveBlock2D_Dimensions                    *
 *                                                      *
 * Member functions                                     *
 *         i    -- Return adaptive block dimension in   *
 *                 the i-direction.                     *
 *         j    -- Return adaptive block dimension in   *
 *                 the j-direction.                     *
 *     ghost    -- Return number of ghost (halo or      *
 *                 overlap) cells for adaptive block.   * 
 *                                                      *
 * Member operators                                     *
 *      D -- adaptive block dimensions                  *
 *                                                      *
 * D = D;                                               *
 * cout << D; (output function)                         *
 * cin  >> D; (input function)                          *
 *                                                      *
 ********************************************************/
class AdaptiveBlock2D_Dimensions{
  private:
  public:
    int           i,j;  // Adaptive block dimensions.
    int         ghost;  // Adaptive block number of ghost 
                        // (halo or overlap) cells.
	                // Made public so can access them.
    
    /* Creation, copy, and assignment constructors. */
    AdaptiveBlock2D_Dimensions(void) {
       i = 0; j = 0; ghost = 0;
    }

    AdaptiveBlock2D_Dimensions(const AdaptiveBlock2D_Dimensions &Blk_Dimen) {
       i = Blk_Dimen.i; j = Blk_Dimen.j; ghost = Blk_Dimen.ghost;
    }

    AdaptiveBlock2D_Dimensions(const int i_dimen,
	                       const int j_dimen,
                               const int ghost_dimen) {
       i = i_dimen; j = j_dimen; ghost = ghost_dimen;
    }

    /* Destructor. */
    // ~AdaptiveBlock2D_Dimensions(void);
    // Use automatically generated destructor.

    /* copy function */
    void copy(const AdaptiveBlock2D_Dimensions &Blk_Dimen) {
      i = Blk_Dimen.i; j = Blk_Dimen.j; ghost = Blk_Dimen.ghost;
    }


    /* Assignment operator. */
    // AdaptiveBlock2D_Dimensions operator = (const AdaptiveBlock2D_Dimensions &Blk_Dimen);
    // Use automatically generated assignment operator.

    //@{ @name Binary arithmetic operators.
    AdaptiveBlock2D_Dimensions operator *(const double &a) const;
    friend AdaptiveBlock2D_Dimensions operator *(const double &a, const AdaptiveBlock2D_Dimensions &Blk_Dimen);
    //@}

    //@{ Shortcut arithmetic operators.
    AdaptiveBlock2D_Dimensions &operator +=(const AdaptiveBlock2D_Dimensions &Blk_Dimen);
    //@}

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const AdaptiveBlock2D_Dimensions &Blk_Dimen);
    friend istream &operator >> (istream &in_file,
				 AdaptiveBlock2D_Dimensions &Blk_Dimen);

};

/**********************************************************************
 * AdaptiveBlock2D_Dimensions -- Binary arithmetic operators.         *
 **********************************************************************/
//inline AdaptiveBlock2D_Dimensions AdaptiveBlock2D_Dimensions::operator *(const double &a) const {
//  return AdaptiveBlock2D_Dimensions(i*a,j*a,ghost*a);
//}

//inline AdaptiveBlock2D_Dimensions operator *(const double &a, const AdaptiveBlock2D_Dimensions &Blk_Dimen) {
//  return AdaptiveBlock2D_Dimensions(Blk_Dimen.i*a,Blk_Dimen.j*a,Blk_Dimen.ghost*a);
//}

/**********************************************************************
 * AdaptiveBlock2D_Dimensions -- Shortcut arithmetic operators.       *
 **********************************************************************/
inline AdaptiveBlock2D_Dimensions& AdaptiveBlock2D_Dimensions::operator +=(const AdaptiveBlock2D_Dimensions &Blk_Dimen) {
  i += Blk_Dimen.i; j += Blk_Dimen.j; ghost += Blk_Dimen.ghost;
  return *this;
}

/*************************************************************
 * AdaptiveBlock2D_Dimensions -- Input-output operators.     *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const AdaptiveBlock2D_Dimensions &Blk_Dimen) {
  out_file << " " << Blk_Dimen.i << " " << Blk_Dimen.j << " " << Blk_Dimen.ghost;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     AdaptiveBlock2D_Dimensions &Blk_Dimen) {
  in_file.setf(ios::skipws);
  in_file >> Blk_Dimen.i >> Blk_Dimen.j >> Blk_Dimen.ghost;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Class: AdaptiveBlock2D_Info                          *
 *                                                      *
 * Member functions                                     *
 *       cpu    -- Return adaptive block processor      *
 *                 or processor element number.         *
 *    blknum    -- Return local adaptive block number   *
 *                 for the given processor.             *
 *     dimen    -- Return adaptive block dimensions.    *
 *    sector    -- Return adaptive block sector.        *
 *     level    -- Return adaptive block refinement     *
 *                 level.                               *
 *                                                      *
 * Member operators                                     *
 *      I -- adaptive block info                        *
 *                                                      *
 * I = I;                                               *
 * cout << I; (output function)                         *
 * cin  >> I; (input function)                          *
 *                                                      *
 ********************************************************/
class AdaptiveBlock2D_Info{
  private:
  public:
    int                            cpu; // Adaptive block processor number.
    int                         blknum; // Local adaptive block number.
    AdaptiveBlock2D_Dimensions   dimen; // Adaptive block dimensions. 
    int                         sector; // Adaptive block sector number.
    int                          level; // Adaptive block refinement level.
	                                // Made public so can access them.
    
    /* Creation, copy, and assignment constructors. */
    AdaptiveBlock2D_Info(void) {
       cpu = 0; blknum = 0; dimen.i = 0; dimen.j = 0; dimen.ghost = 0; 
       sector = ADAPTIVEBLOCK2D_SECTOR_NONE; level = 0;
    }

    AdaptiveBlock2D_Info(const AdaptiveBlock2D_Info &Blk_Info) {
       cpu = Blk_Info.cpu; blknum = Blk_Info.blknum;
       dimen = Blk_Info.dimen; sector = Blk_Info.sector; level = Blk_Info.level;
    }

    AdaptiveBlock2D_Info(const int i_processor,
	                 const int i_block,
                         const int i_dimen,
			 const int j_dimen,
                         const int ghost_dimen,
                         const int i_sector,
	                 const int i_level) {
       cpu = i_processor; blknum = i_block;
       dimen.i = i_dimen; dimen.j = j_dimen; dimen.ghost = ghost_dimen;
       sector = i_sector; level = i_level;
    }

    /* Destructor. */
    // ~AdaptiveBlock2D_Info(void);
    // Use automatically generated destructor.

    /* copy function */
    void copy(const AdaptiveBlock2D_Info &Blk_Info) {
       cpu = Blk_Info.cpu; blknum = Blk_Info.blknum;
       dimen.copy(Blk_Info.dimen); sector = Blk_Info.sector; level = Blk_Info.level;
    }

    /* Assignment operator. */
    // AdaptiveBlock2D_Info operator = (const AdaptiveBlock2D_Info &Blk_Info);
    // Use automatically generated assignment operator.

    /* Reset the block information. */
    void reset(void);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const AdaptiveBlock2D_Info &Blk_Info);
    friend istream &operator >> (istream &in_file,
				 AdaptiveBlock2D_Info &Blk_Info);

};

/*************************************************************
 * AdaptiveBlock2D_Info -- Reset the block information.      *
 *************************************************************/
inline void AdaptiveBlock2D_Info::reset(void) {
  cpu = 0; 
  blknum = 0; 
  dimen.i = 0; dimen.j = 0; dimen.ghost = 0; 
  sector = ADAPTIVEBLOCK2D_SECTOR_NONE;
  level = 0;
}

/*************************************************************
 * AdaptiveBlock2D_Info -- Input-output operators.           *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const AdaptiveBlock2D_Info &Blk_Info) {
  out_file << " " << Blk_Info.cpu 
           << " " << Blk_Info.blknum << Blk_Info.dimen 
           << " " << Blk_Info.sector
           << " " << Blk_Info.level;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     AdaptiveBlock2D_Info &Blk_Info) {
  in_file.setf(ios::skipws);
  in_file >> Blk_Info.cpu >> Blk_Info.blknum;
  in_file.unsetf(ios::skipws);
  in_file >> Blk_Info.dimen;
  in_file.setf(ios::skipws);
  in_file >> Blk_Info.sector >> Blk_Info.level;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/*************************************************************
 * AdaptiveBlock2D_Info -- External subroutines.             *
 *************************************************************/

extern void Broadcast_Adaptive_Block_Info(AdaptiveBlock2D_Info &Blk_Info);

/********************************************************
 * Class: AdaptiveBlock2D                               *
 *                                                      *
 * Member functions                                     *
 *      used     -- Return adaptive block usage.        *
 *   gblknum     -- Return global block number.         *
 *      info     -- Return adaptive block info.         *
 *        nN     -- Return number of north neighbour    *
 *                  blocks.                             *
 *        nS     -- Return number of south neighbour    *
 *                  blocks.                             *
 *        nE     -- Return number of east neighbour     *
 *                  blocks.                             *
 *        nW     -- Return number of west neighbour     *
 *                  blocks.                             *
 *       nNW     -- Return number of north-west corner  *
 *                  neighbour blocks.                   *
 *       nNE     -- Return number of north-east corner  *
 *                  neighbour blocks.                   *
 *       nSE     -- Return number of south-east corner  *
 *                  neighbour blocks.                   *
 *       nSW     -- Return number of south-west corner  *
 *                  neighbour blocks.                   *
 *     infoN     -- Return info for north neighbour     *
 *                  blocks.                             *
 *     infoS     -- Return info for south neighbour     *
 *                  blocks.                             *
 *     infoE     -- Return info for east neighbour      *
 *                  blocks.                             *
 *     infoW     -- Return info for west neighbour      *
 *                  blocks.                             *
 *    infoNW     -- Return info for north-west corner   *
 *                  neighbour blocks.                   *
 *    infoNE     -- Return info for north-east corner   *
 *                  neighbour blocks.                   *
 *    infoSE     -- Return info for south-east corner   *
 *                  neighbour blocks.                   *
 *    infoSW     -- Return info for south-west corner   *
 *                  neighbour blocks.                   * 
 *                                                      *
 * NOTE:                                                *
 * If a dimension of a neighbouring block is negative,  *
 * then the cell ordering is reversed for that          *
 * neighbouring block in the computational coordinate   *
 * direction corresponding to the negative dimension.   *
 *                                                      *
 * Member operators                                     *
 *      B -- adaptive block                             *
 *                                                      *
 * B = B;                                               *
 * cout << B; (output function)                         *
 * cin  >> B; (input function)                          *
 *                                                      *
 ********************************************************/
class AdaptiveBlock2D{
  private:
  public:
    int                    used;  // Adaptive block usage indicator.
    int                 gblknum;  // Global adaptive block number.
    AdaptiveBlock2D_Info   info;  // Adaptive block Info.
    int             nN,nS,nE,nW;  // Number of neighbouring blocks
                                  // to the north, south, east, and west.    
    AdaptiveBlock2D_Info infoN[ADAPTIVEBLOCK2D_NUMBER_NEIGHBOURS_MAX],
                         infoS[ADAPTIVEBLOCK2D_NUMBER_NEIGHBOURS_MAX],
                         infoE[ADAPTIVEBLOCK2D_NUMBER_NEIGHBOURS_MAX],
                         infoW[ADAPTIVEBLOCK2D_NUMBER_NEIGHBOURS_MAX];
                                  // Block info of neighbouring blocks
                                  // to the north, south, east, and west.
    int         nNW,nNE,nSE,nSW;  // Number of neighboring blocks
                                  // to the north-west, north-east,
                                  // south-east, and south-west corners.
    AdaptiveBlock2D_Info infoNW[ADAPTIVEBLOCK2D_NUMBER_CORNER_NEIGHBOURS_MAX],
                         infoNE[ADAPTIVEBLOCK2D_NUMBER_CORNER_NEIGHBOURS_MAX],
                         infoSE[ADAPTIVEBLOCK2D_NUMBER_CORNER_NEIGHBOURS_MAX],
                         infoSW[ADAPTIVEBLOCK2D_NUMBER_CORNER_NEIGHBOURS_MAX];
                                  // Block info of neighbouring blocks
                                  // to the north-west, north-east,
                                  // south-east, and south-west corners.
	                          // Made public so can access them.
    
    /* Creation, copy, and assignment constructors. */
    AdaptiveBlock2D(void) {
       used = ADAPTIVEBLOCK2D_NOT_USED; gblknum = 0;
       nN = 0; nS = 0; nE = 0; nW = 0;
       nNW = 0; nNE = 0; nSE = 0; nSW = 0;
    }

    AdaptiveBlock2D(const AdaptiveBlock2D &Blk) {
       int i;
       used = Blk.used; gblknum = Blk.gblknum; info = Blk.info;
       nN = Blk.nN; nS = Blk.nS; nE = Blk.nE; nW = Blk.nW;
       for (i = 0; i <= ADAPTIVEBLOCK2D_NUMBER_NEIGHBOURS_MAX-1; ++i) {
	 infoN[i] = Blk.infoN[i]; infoS[i] = Blk.infoS[i];
	 infoE[i] = Blk.infoE[i]; infoW[i] = Blk.infoW[i];
       } /* endfor */
       nNW = Blk.nNW; nNE = Blk.nNE; nSE = Blk.nSE; nSW = Blk.nSW;
       for (i = 0; i <= ADAPTIVEBLOCK2D_NUMBER_CORNER_NEIGHBOURS_MAX-1; ++i) {
	 infoNW[i] = Blk.infoNW[i]; infoNE[i] = Blk.infoNE[i];
	 infoSE[i] = Blk.infoSE[i]; infoSW[i] = Blk.infoSW[i];
       } /* endfor */
    }

    AdaptiveBlock2D(const AdaptiveBlock2D *Blk) {
       int i;
       used = Blk->used; gblknum = Blk->gblknum; info = Blk->info;
       nN = Blk->nN; nS = Blk->nS; nE = Blk->nE; nW = Blk->nW;
       for (i = 0; i <= ADAPTIVEBLOCK2D_NUMBER_NEIGHBOURS_MAX-1; ++i) {
	 infoN[i] = Blk->infoN[i]; infoS[i] = Blk->infoS[i];
	 infoE[i] = Blk->infoE[i]; infoW[i] = Blk->infoW[i];
       } /* endfor */
       nNW = Blk->nNW; nNE = Blk->nNE; nSE = Blk->nSE; nSW = Blk->nSW;
       for (i = 0; i <= ADAPTIVEBLOCK2D_NUMBER_CORNER_NEIGHBOURS_MAX-1; ++i) {
	 infoNW[i] = Blk->infoNW[i]; infoNE[i] = Blk->infoNE[i];
	 infoSE[i] = Blk->infoSE[i]; infoSW[i] = Blk->infoSW[i];
       } /* endfor */
    }

    /* Destructor. */
    // ~AdaptiveBlock2D(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // AdaptiveBlock2D operator =(const AdaptiveBlock2D &Blk);
    // Use automatically generated assignment operator.

    /* copy function */
    void copy(const AdaptiveBlock2D &Blk) {
      int i;
      used = Blk.used; gblknum = Blk.gblknum; info.copy(Blk.info);
       nN = Blk.nN; nS = Blk.nS; nE = Blk.nE; nW = Blk.nW;
       for (i = 0; i <= ADAPTIVEBLOCK2D_NUMBER_NEIGHBOURS_MAX-1; ++i) {
	 infoN[i].copy(Blk.infoN[i]); infoS[i].copy(Blk.infoS[i]);
	 infoE[i].copy(Blk.infoE[i]); infoW[i].copy(Blk.infoW[i]);
       } /* endfor */
       nNW = Blk.nNW; nNE = Blk.nNE; nSE = Blk.nSE; nSW = Blk.nSW;
       for (i = 0; i <= ADAPTIVEBLOCK2D_NUMBER_CORNER_NEIGHBOURS_MAX-1; ++i) {
	 infoNW[i].copy(Blk.infoNW[i]); infoNE[i].copy(Blk.infoNE[i]);
	 infoSE[i].copy(Blk.infoSE[i]); infoSW[i].copy(Blk.infoSW[i]);
       } /* endfor */
    }

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const AdaptiveBlock2D &Blk);
    friend istream &operator >> (istream &in_file,
				 AdaptiveBlock2D &Blk);

};

/*************************************************************
 * AdaptiveBlock2D -- Input-output operators.                *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const AdaptiveBlock2D &Blk) {
  int i;
  out_file << " " << Blk.used << " " << Blk.gblknum << Blk.info;
  out_file << "\n " << Blk.nN;
  for (i = 0; i <= Blk.nN-1; ++i) {
     out_file << " " << Blk.infoN[i];
  } /* endfor */
  out_file << "\n " << Blk.nS;
  for (i = 0; i <= Blk.nS-1; ++i) {
     out_file << " " << Blk.infoS[i];
  } /* endfor */
  out_file << "\n " << Blk.nE;
  for (i = 0; i <= Blk.nE-1; ++i) {
     out_file << " " << Blk.infoE[i];
  } /* endfor */
  out_file << "\n " << Blk.nW;
  for (i = 0; i <= Blk.nW-1; ++i) {
     out_file << " " << Blk.infoW[i];
  } /* endfor */
  out_file << "\n " << Blk.nNW;
  for (i = 0; i <= Blk.nNW-1; ++i) {
     out_file << " " << Blk.infoNW[i];
  } /* endfor */
  out_file << "\n " << Blk.nNE;
  for (i = 0; i <= Blk.nNE-1; ++i) {
     out_file << " " << Blk.infoNE[i];
  } /* endfor */
  out_file << "\n " << Blk.nSE;
  for (i = 0; i <= Blk.nSE-1; ++i) {
     out_file << " " << Blk.infoSE[i];
  } /* endfor */
  out_file << "\n " << Blk.nSW;
  for (i = 0; i <= Blk.nSW-1; ++i) {
     out_file << " " << Blk.infoSW[i];
  } /* endfor */
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     AdaptiveBlock2D &Blk) {
  int i;
  in_file.setf(ios::skipws); in_file >> Blk.used >> Blk.gblknum; 
  in_file.unsetf(ios::skipws);
  in_file >> Blk.info;
  in_file.setf(ios::skipws); in_file >> Blk.nN; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nN-1; ++i) {
     in_file >> Blk.infoN[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nS; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nS-1; ++i) {
     in_file >> Blk.infoS[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nE; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nE-1; ++i) {
     in_file >> Blk.infoE[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nW; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nW-1; ++i) {
     in_file >> Blk.infoW[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nNW; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nNW-1; ++i) {
     in_file >> Blk.infoNW[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nNE; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nNE-1; ++i) {
     in_file >> Blk.infoNE[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nSE; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nSE-1; ++i) {
     in_file >> Blk.infoSE[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nSW; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nSW-1; ++i) {
     in_file >> Blk.infoSW[i];
  } /* endfor */
  return (in_file);
}

/*************************************************************
 * AdaptiveBlock2D -- External subroutines.                  *
 *************************************************************/

extern void Broadcast_Adaptive_Block(AdaptiveBlock2D &Blk);

/********************************************************
 * Class: AdaptiveBlock2D_List                          *
 *                                                      *
 * Member functions                                     *
 *     ThisCPU -- Return the global MPI number or       *
 *                rank for the CPU that is executing    *
 *                the task or process.                  *
 *        Nblk -- Return number of local adaptive       *
 *                blocks in the list.                   *
 *       Block -- Return list of adaptive blocks.       *
 *       RefineFlag -- Return list of mesh refinement   *
 *                     flags for the adaptive blocks.   *
 *       message_noreschange_northface_sendbuf          *
 *       message_noreschange_northface_recbuf           *
 *       message_noreschange_southface_sendbuf          *
 *       message_noreschange_southface_recbuf           *
 *       message_noreschange_eastface_sendbuf           *
 *       message_noreschange_eastface_recbuf            *
 *       message_noreschange_westface_sendbuf           *
 *       message_noreschange_westface_recbuf            *
 *       message_noreschange_northwestcorner_sendbuf    *
 *       message_noreschange_northwestcorner_recbuf     *
 *       message_noreschange_northeastcorner_sendbuf    *
 *       message_noreschange_northeastcorner_recbuf     *
 *       message_noreschange_southeastcorner_sendbuf    *
 *       message_noreschange_southeastcorner_recbuf     *
 *       message_noreschange_southwestcorner_sendbuf    *
 *       message_noreschange_southwestcorner_recbuf     *
 *        -- Return message passing send and receive    *
 *           buffers for block faces and corners with   *
 *           no block resolution (refinement) change.   *
 *    allocate -- Allocate memory for block list.       *
 *  deallocate -- Deallocate memory for block list.     *
 *  nochangeAll-- Sets the mesh refinement flags to     *
 *                force no refinement or coarsening of  *
 *                adaptive blocks (default).            *
 *  refineAll  -- Sets the mesh refinement flags to     *
 *                force refinement of all adaptive      *
 *                blocks.                               *
 *  coarsenAll -- Sets the mesh refinement flags to     *
 *                force coasening of all adaptive       *
 *                blocks.                               *
 *  setRefineAll -- Sets the mesh refinement flags of   *
 *                  all blocks to specified value.      *
 *                                                      *
 * Member operators                                     *
 *      B -- adaptive block list                        *
 *                                                      *
 * B = B;                                               *
 * cout << B; (output function)                         *
 * cin  >> B; (input function)                          *
 *                                                      *
 ********************************************************/
class AdaptiveBlock2D_List{
  private:
  public:
    int              ThisCPU; // Number or rank of CPU executing the
                              // current process.
    int                 Nblk; // Number of local adaptive
                              // blocks in the list.
    AdaptiveBlock2D   *Block; // List of adaptive blocks.
    int          *RefineFlag; // Block refinement (coarsening) flag.
    double **message_noreschange_northface_sendbuf,
      **message_noreschange_northface_recbuf,
      **message_noreschange_southface_sendbuf,
      **message_noreschange_southface_recbuf,
      **message_noreschange_eastface_sendbuf,
      **message_noreschange_eastface_recbuf,
      **message_noreschange_westface_sendbuf,
      **message_noreschange_westface_recbuf,
      **message_noreschange_northwestcorner_sendbuf,
      **message_noreschange_northwestcorner_recbuf,
      **message_noreschange_northeastcorner_sendbuf,
      **message_noreschange_northeastcorner_recbuf,
      **message_noreschange_southeastcorner_sendbuf,
      **message_noreschange_southeastcorner_recbuf,
      **message_noreschange_southwestcorner_sendbuf,
      **message_noreschange_southwestcorner_recbuf;
                              // Message passing send and receive buffers
                              // for block faces and corners with no block
                              // resolution (refinement) change.
    double ***message_reschange_northface_sendbuf,
      ***message_reschange_northface_recbuf,
      ***message_reschange_southface_sendbuf,
      ***message_reschange_southface_recbuf,
      ***message_reschange_eastface_sendbuf,
      ***message_reschange_eastface_recbuf,
      ***message_reschange_westface_sendbuf,
      ***message_reschange_westface_recbuf,
      **message_reschange_northwestcorner_sendbuf,
      **message_reschange_northwestcorner_recbuf,
      **message_reschange_northeastcorner_sendbuf,
      **message_reschange_northeastcorner_recbuf,
      **message_reschange_southeastcorner_sendbuf,
      **message_reschange_southeastcorner_recbuf,
      **message_reschange_southwestcorner_sendbuf,
      **message_reschange_southwestcorner_recbuf;
                              // Message passing send and receive buffers
                              // for block faces and corners with block
                              // resolution (refinement) change.
                              // Made public so can access them.
    
    /* Creation, copy, and assignment constructors. */
    AdaptiveBlock2D_List(void) {
      ThisCPU = 0; Nblk = 0; Block = NULL; RefineFlag = NULL;
      // No resolution change message buffers.
      message_noreschange_northface_sendbuf = NULL;
      message_noreschange_northface_recbuf = NULL;
      message_noreschange_southface_sendbuf = NULL;
      message_noreschange_southface_recbuf = NULL;
      message_noreschange_eastface_sendbuf = NULL;
      message_noreschange_eastface_recbuf = NULL;
      message_noreschange_westface_sendbuf = NULL;
      message_noreschange_westface_recbuf = NULL;
      message_noreschange_northwestcorner_sendbuf = NULL;
      message_noreschange_northwestcorner_recbuf = NULL;
      message_noreschange_northeastcorner_sendbuf = NULL;
      message_noreschange_northeastcorner_recbuf = NULL;
      message_noreschange_southeastcorner_sendbuf = NULL;
      message_noreschange_southeastcorner_recbuf = NULL;
      message_noreschange_southwestcorner_sendbuf = NULL;
      message_noreschange_southwestcorner_recbuf = NULL;
      // Resolution change message buffers.
      message_reschange_northface_sendbuf = NULL;
      message_reschange_northface_recbuf = NULL;
      message_reschange_southface_sendbuf = NULL;
      message_reschange_southface_recbuf = NULL;
      message_reschange_eastface_sendbuf = NULL;
      message_reschange_eastface_recbuf = NULL;
      message_reschange_westface_sendbuf = NULL;
      message_reschange_westface_recbuf = NULL;
      message_reschange_northwestcorner_sendbuf = NULL;
      message_reschange_northwestcorner_recbuf = NULL;
      message_reschange_northeastcorner_sendbuf = NULL;
      message_reschange_northeastcorner_recbuf = NULL;
      message_reschange_southeastcorner_sendbuf = NULL;
      message_reschange_southeastcorner_recbuf = NULL;
      message_reschange_southwestcorner_sendbuf = NULL;
      message_reschange_southwestcorner_recbuf = NULL;
    }

    AdaptiveBlock2D_List(const AdaptiveBlock2D_List &Blk_List) {
      ThisCPU = Blk_List.ThisCPU; Nblk = Blk_List.Nblk; 
      Block = Blk_List.Block; RefineFlag = Blk_List.RefineFlag;
      // No resolution change message buffers.
      message_noreschange_northface_sendbuf = Blk_List.message_noreschange_northface_sendbuf;
      message_noreschange_northface_recbuf = Blk_List.message_noreschange_northface_recbuf;
      message_noreschange_southface_sendbuf = Blk_List.message_noreschange_southface_sendbuf;
      message_noreschange_southface_recbuf = Blk_List.message_noreschange_southface_recbuf;
      message_noreschange_eastface_sendbuf = Blk_List.message_noreschange_eastface_sendbuf;
      message_noreschange_eastface_recbuf = Blk_List.message_noreschange_eastface_recbuf;
      message_noreschange_westface_sendbuf = Blk_List.message_noreschange_westface_sendbuf;
      message_noreschange_westface_recbuf = Blk_List.message_noreschange_westface_recbuf;
      message_noreschange_northwestcorner_sendbuf = Blk_List.message_noreschange_northwestcorner_sendbuf;
      message_noreschange_northwestcorner_recbuf = Blk_List.message_noreschange_northwestcorner_recbuf;
      message_noreschange_northeastcorner_sendbuf = Blk_List.message_noreschange_northeastcorner_sendbuf;
      message_noreschange_northeastcorner_recbuf = Blk_List.message_noreschange_northeastcorner_recbuf;
      message_noreschange_southeastcorner_sendbuf = Blk_List.message_noreschange_southeastcorner_sendbuf;
      message_noreschange_southeastcorner_recbuf = Blk_List.message_noreschange_southeastcorner_recbuf;
      message_noreschange_southwestcorner_sendbuf = Blk_List.message_noreschange_southwestcorner_sendbuf;
      message_noreschange_southwestcorner_recbuf = Blk_List.message_noreschange_southwestcorner_recbuf;
      // Resolution change message buffers.
      message_reschange_northface_sendbuf = Blk_List.message_reschange_northface_sendbuf;
      message_reschange_northface_recbuf = Blk_List.message_reschange_northface_recbuf;
      message_reschange_southface_sendbuf = Blk_List.message_reschange_southface_sendbuf;
      message_reschange_southface_recbuf = Blk_List.message_reschange_southface_recbuf;
      message_reschange_eastface_sendbuf = Blk_List.message_reschange_eastface_sendbuf;
      message_reschange_eastface_recbuf = Blk_List.message_reschange_eastface_recbuf;
      message_reschange_westface_sendbuf = Blk_List.message_reschange_westface_sendbuf;
      message_reschange_westface_recbuf = Blk_List.message_reschange_westface_recbuf;
      message_reschange_northwestcorner_sendbuf = Blk_List.message_reschange_northwestcorner_sendbuf;
      message_reschange_northwestcorner_recbuf = Blk_List.message_reschange_northwestcorner_recbuf;
      message_reschange_northeastcorner_sendbuf = Blk_List.message_reschange_northeastcorner_sendbuf;
      message_reschange_northeastcorner_recbuf = Blk_List.message_reschange_northeastcorner_recbuf;
      message_reschange_southeastcorner_sendbuf = Blk_List.message_reschange_southeastcorner_sendbuf;
      message_reschange_southeastcorner_recbuf = Blk_List.message_reschange_southeastcorner_recbuf;
      message_reschange_southwestcorner_sendbuf = Blk_List.message_reschange_southwestcorner_sendbuf;
      message_reschange_southwestcorner_recbuf = Blk_List.message_reschange_southwestcorner_recbuf;
    }

    /* Destructor. */
    // ~AdaptiveBlock2D_List(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // AdaptiveBlock2D_List operator = (const AdaptiveBlock2D_List &Blk_List);
    // Use automatically generated assignment operator.

    /* Allocate memory for adaptive block list. */
    void allocate(const int N);

    /* Deallocate memory for adaptive block list. */
    void deallocate(void);

    /* Set mesh refinement and coarsening flags. */
    void nochangeAll(void);
    void refineAll(void);
    void coarsenAll(void);
    void setRefineAll(const int Flag);

    /* Returns the number of blocks used. */
    int Nused(void);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const AdaptiveBlock2D_List &Blk_List);
    friend istream &operator >> (istream &in_file,
				 AdaptiveBlock2D_List &Blk_List);

};

/*************************************************************
 * AdaptiveBlock2D_List::allocate -- Allocate memory.        *
 *************************************************************/
inline void AdaptiveBlock2D_List::allocate(const int N) {
   int i, j; assert( N > 0 ); Nblk = N;
   Block = new AdaptiveBlock2D[Nblk]; 
   RefineFlag = new int[Nblk]; nochangeAll();
   //
   // No resolution change message buffers.
   message_noreschange_northface_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      message_noreschange_northface_sendbuf[i] = new double[1];
      message_noreschange_northface_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_northface_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_northface_recbuf[i] = new double[1];
      message_noreschange_northface_recbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_southface_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_southface_sendbuf[i] = new double[1];
      message_noreschange_southface_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_southface_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      message_noreschange_southface_recbuf[i] = new double[1];
      message_noreschange_southface_recbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_eastface_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_eastface_sendbuf[i] = new double[1];
      message_noreschange_eastface_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_eastface_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_eastface_recbuf[i] = new double[1];
      message_noreschange_eastface_recbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_westface_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_westface_sendbuf[i] = new double[1];
      message_noreschange_westface_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_westface_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_westface_recbuf[i] = new double[1];
      message_noreschange_westface_recbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_northwestcorner_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_northwestcorner_sendbuf[i] = new double[1];
      message_noreschange_northwestcorner_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_northwestcorner_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_northwestcorner_recbuf[i] = new double[1];
      message_noreschange_northwestcorner_recbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_northeastcorner_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_northeastcorner_sendbuf[i] = new double[1];
      message_noreschange_northeastcorner_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_northeastcorner_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_northeastcorner_recbuf[i] = new double[1];
      message_noreschange_northeastcorner_recbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_southeastcorner_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_southeastcorner_sendbuf[i] = new double[1];
      message_noreschange_southeastcorner_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_southeastcorner_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_noreschange_southeastcorner_recbuf[i] = new double[1];
      message_noreschange_southeastcorner_recbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_southwestcorner_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      message_noreschange_southwestcorner_sendbuf[i] = new double[1];
      message_noreschange_southwestcorner_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_noreschange_southwestcorner_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      message_noreschange_southwestcorner_recbuf[i] = new double[1];
      message_noreschange_southwestcorner_recbuf[i][0] = ZERO;
   } /* endfor */
   //
   // Resolution change message buffers.
   message_reschange_northface_sendbuf = new double**[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      message_reschange_northface_sendbuf[i] = new double*[2];
      for (j = 0; j < 2 ; ++j) {
         message_reschange_northface_sendbuf[i][j] = new double[1];
         message_reschange_northface_sendbuf[i][j][0] = ZERO;
      } /* endfor */
   } /* endfor */
   message_reschange_northface_recbuf = new double**[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_reschange_northface_recbuf[i] = new double*[2];
      for (j = 0; j < 2 ; ++j) {
         message_reschange_northface_recbuf[i][j] = new double[1];
         message_reschange_northface_recbuf[i][j][0] = ZERO;
      } /* endfor */
   } /* endfor */
   message_reschange_southface_sendbuf = new double**[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      message_reschange_southface_sendbuf[i] = new double*[2];
      for (j = 0; j < 2 ; ++j) {
         message_reschange_southface_sendbuf[i][j] = new double[1];
         message_reschange_southface_sendbuf[i][j][0] = ZERO;
      } /* endfor */
   } /* endfor */
   message_reschange_southface_recbuf = new double**[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_reschange_southface_recbuf[i] = new double*[2];
      for (j = 0; j < 2 ; ++j) {
         message_reschange_southface_recbuf[i][j] = new double[1];
         message_reschange_southface_recbuf[i][j][0] = ZERO;
      } /* endfor */
   } /* endfor */
   message_reschange_eastface_sendbuf = new double**[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      message_reschange_eastface_sendbuf[i] = new double*[2];
      for (j = 0; j < 2 ; ++j) {
         message_reschange_eastface_sendbuf[i][j] = new double[1];
         message_reschange_eastface_sendbuf[i][j][0] = ZERO;
      } /* endfor */
   } /* endfor */
   message_reschange_eastface_recbuf = new double**[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_reschange_eastface_recbuf[i] = new double*[2];
      for (j = 0; j < 2 ; ++j) {
         message_reschange_eastface_recbuf[i][j] = new double[1];
         message_reschange_eastface_recbuf[i][j][0] = ZERO;
      } /* endfor */
   } /* endfor */
   message_reschange_westface_sendbuf = new double**[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      message_reschange_westface_sendbuf[i] = new double*[2];
      for (j = 0; j < 2 ; ++j) {
         message_reschange_westface_sendbuf[i][j] = new double[1];
         message_reschange_westface_sendbuf[i][j][0] = ZERO;
      } /* endfor */
   } /* endfor */
   message_reschange_westface_recbuf = new double**[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_reschange_westface_recbuf[i] = new double*[2];
      for (j = 0; j < 2 ; ++j) {
         message_reschange_westface_recbuf[i][j] = new double[1];
         message_reschange_westface_recbuf[i][j][0] = ZERO;
      } /* endfor */
   } /* endfor */
   message_reschange_northwestcorner_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_reschange_northwestcorner_sendbuf[i] = new double[1];
      message_reschange_northwestcorner_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_reschange_northwestcorner_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_reschange_northwestcorner_recbuf[i] = new double[1];
      message_reschange_northwestcorner_recbuf[i][0] = ZERO;
   } /* endfor */
   message_reschange_northeastcorner_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_reschange_northeastcorner_sendbuf[i] = new double[1];
      message_reschange_northeastcorner_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_reschange_northeastcorner_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_reschange_northeastcorner_recbuf[i] = new double[1];
      message_reschange_northeastcorner_recbuf[i][0] = ZERO;
   } /* endfor */
   message_reschange_southeastcorner_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_reschange_southeastcorner_sendbuf[i] = new double[1];
      message_reschange_southeastcorner_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_reschange_southeastcorner_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) { 
      message_reschange_southeastcorner_recbuf[i] = new double[1];
      message_reschange_southeastcorner_recbuf[i][0] = ZERO;
   } /* endfor */
   message_reschange_southwestcorner_sendbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      message_reschange_southwestcorner_sendbuf[i] = new double[1];
      message_reschange_southwestcorner_sendbuf[i][0] = ZERO;
   } /* endfor */
   message_reschange_southwestcorner_recbuf = new double*[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      message_reschange_southwestcorner_recbuf[i] = new double[1];
      message_reschange_southwestcorner_recbuf[i][0] = ZERO;
   } /* endfor */
}

/*************************************************************
 * AdaptiveBlock2D_List::deallocate -- Deallocate memory.    *
 *************************************************************/
inline void AdaptiveBlock2D_List::deallocate(void) {
   int i, j; delete []Block; Block = NULL;
   delete []RefineFlag; RefineFlag = NULL;
   //
   // No resolution change message buffers.
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_northface_sendbuf[i]; 
      message_noreschange_northface_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_northface_sendbuf;
   message_noreschange_northface_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_northface_recbuf[i]; 
      message_noreschange_northface_recbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_northface_recbuf;
   message_noreschange_northface_recbuf = NULL;
  for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_southface_sendbuf[i]; 
      message_noreschange_southface_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_southface_sendbuf;
   message_noreschange_southface_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_southface_recbuf[i]; 
      message_noreschange_southface_recbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_southface_recbuf;
   message_noreschange_southface_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_eastface_sendbuf[i]; 
      message_noreschange_eastface_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_eastface_sendbuf;
   message_noreschange_eastface_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_eastface_recbuf[i]; 
      message_noreschange_eastface_recbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_eastface_recbuf;
   message_noreschange_eastface_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_westface_sendbuf[i]; 
      message_noreschange_westface_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_westface_sendbuf;
   message_noreschange_westface_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_westface_recbuf[i]; 
      message_noreschange_westface_recbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_westface_recbuf;
   message_noreschange_westface_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_northwestcorner_sendbuf[i]; 
      message_noreschange_northwestcorner_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_northwestcorner_sendbuf;
   message_noreschange_northwestcorner_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_northwestcorner_recbuf[i]; 
      message_noreschange_northwestcorner_recbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_northwestcorner_recbuf;
   message_noreschange_northwestcorner_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_northeastcorner_sendbuf[i]; 
      message_noreschange_northeastcorner_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_northeastcorner_sendbuf;
   message_noreschange_northeastcorner_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_northeastcorner_recbuf[i]; 
      message_noreschange_northeastcorner_recbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_northeastcorner_recbuf;
   message_noreschange_northeastcorner_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_southeastcorner_sendbuf[i]; 
      message_noreschange_southeastcorner_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_southeastcorner_sendbuf;
   message_noreschange_southeastcorner_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_southeastcorner_recbuf[i]; 
      message_noreschange_southeastcorner_recbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_southeastcorner_recbuf;
   message_noreschange_southeastcorner_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_southwestcorner_sendbuf[i]; 
      message_noreschange_southwestcorner_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_southwestcorner_sendbuf;
   message_noreschange_southwestcorner_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_noreschange_southwestcorner_recbuf[i]; 
      message_noreschange_southwestcorner_recbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_southwestcorner_recbuf;
   message_noreschange_southwestcorner_recbuf = NULL;
   //
   // Resolution change message buffers.
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      for (j = 0; j < 2 ; ++j) {
         delete []message_reschange_northface_sendbuf[i][j]; 
         message_reschange_northface_sendbuf[i][j] = NULL;
      } /* endfor */
      delete []message_reschange_northface_sendbuf[i]; 
      message_reschange_northface_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_northface_sendbuf;
   message_reschange_northface_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      for (j = 0; j < 2 ; ++j) {
         delete []message_reschange_northface_recbuf[i][j]; 
         message_reschange_northface_recbuf[i][j] = NULL;
      } /* endfor */
      delete []message_reschange_northface_recbuf[i]; 
      message_reschange_northface_recbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_northface_recbuf;
   message_reschange_northface_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      for (j = 0; j < 2 ; ++j) {
         delete []message_reschange_southface_sendbuf[i][j]; 
         message_reschange_southface_sendbuf[i][j] = NULL;
      } /* endfor */
      delete []message_reschange_southface_sendbuf[i]; 
      message_reschange_southface_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_southface_sendbuf;
   message_reschange_southface_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      for (j = 0; j < 2 ; ++j) {
         delete []message_reschange_southface_recbuf[i][j]; 
         message_reschange_southface_recbuf[i][j] = NULL;
      } /* endfor */
      delete []message_reschange_southface_recbuf[i]; 
      message_reschange_southface_recbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_southface_recbuf;
   message_reschange_southface_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      for (j = 0; j < 2 ; ++j) {
         delete []message_reschange_eastface_sendbuf[i][j]; 
         message_reschange_eastface_sendbuf[i][j] = NULL;
      } /* endfor */
      delete []message_reschange_eastface_sendbuf[i]; 
      message_reschange_eastface_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_eastface_sendbuf;
   message_reschange_eastface_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      for (j = 0; j < 2 ; ++j) {
         delete []message_reschange_eastface_recbuf[i][j]; 
         message_reschange_eastface_recbuf[i][j] = NULL;
      } /* endfor */
      delete []message_reschange_eastface_recbuf[i]; 
      message_reschange_eastface_recbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_eastface_recbuf;
   message_reschange_eastface_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      for (j = 0; j < 2 ; ++j) {
         delete []message_reschange_westface_sendbuf[i][j]; 
         message_reschange_westface_sendbuf[i][j] = NULL;
      } /* endfor */
      delete []message_reschange_westface_sendbuf[i]; 
      message_reschange_westface_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_westface_sendbuf;
   message_reschange_westface_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      for (j = 0; j < 2 ; ++j) {
         delete []message_reschange_westface_recbuf[i][j];
         message_reschange_westface_recbuf[i][j] = NULL;
      } /* endfor */
      delete []message_reschange_westface_recbuf[i]; 
      message_reschange_westface_recbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_westface_recbuf;
   message_reschange_westface_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_reschange_northwestcorner_sendbuf[i]; 
      message_reschange_northwestcorner_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_northwestcorner_sendbuf;
   message_reschange_northwestcorner_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_reschange_northwestcorner_recbuf[i]; 
      message_reschange_northwestcorner_recbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_northwestcorner_recbuf;
   message_reschange_northwestcorner_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_reschange_northeastcorner_sendbuf[i]; 
      message_reschange_northeastcorner_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_northeastcorner_sendbuf;
   message_reschange_northeastcorner_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_reschange_northeastcorner_recbuf[i]; 
      message_reschange_northeastcorner_recbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_northeastcorner_recbuf;
   message_reschange_northeastcorner_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_reschange_southeastcorner_sendbuf[i]; 
      message_reschange_southeastcorner_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_southeastcorner_sendbuf;
   message_reschange_southeastcorner_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_reschange_southeastcorner_recbuf[i]; 
      message_reschange_southeastcorner_recbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_southeastcorner_recbuf;
   message_reschange_southeastcorner_recbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_reschange_southwestcorner_sendbuf[i]; 
      message_reschange_southwestcorner_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_southwestcorner_sendbuf;
   message_reschange_southwestcorner_sendbuf = NULL;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      delete []message_reschange_southwestcorner_recbuf[i]; 
      message_reschange_southwestcorner_recbuf[i] = NULL;
   } /* endfor */
   delete []message_reschange_southwestcorner_recbuf;
   message_reschange_southwestcorner_recbuf = NULL;
   //
   // Reset the number of solution blocks.
   Nblk = 0;
}

/*****************************************************************
 * AdaptiveBlock2D_List::nochangeAll -- Set no refinement flags. *
 *****************************************************************/
inline void AdaptiveBlock2D_List::nochangeAll(void) {
   int i; 
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      RefineFlag[i] = ADAPTIVEBLOCK2D_NOCHANGE;
   } /* endfor */
}

/*****************************************************************
 * AdaptiveBlock2D_List::refineAll -- Set refinement flags.      *
 *****************************************************************/
inline void AdaptiveBlock2D_List::refineAll(void) {
   int i; 
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      RefineFlag[i] = ADAPTIVEBLOCK2D_REFINE;
   } /* endfor */
}

/*****************************************************************
 * AdaptiveBlock2D_List::coarsenAll -- Set coarsening flags.     *
 *****************************************************************/
inline void AdaptiveBlock2D_List::coarsenAll(void) {
   int i; 
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      RefineFlag[i] = ADAPTIVEBLOCK2D_COARSEN;
   } /* endfor */
}

/*****************************************************************
 * AdaptiveBlock2D_List::setRefineAll -- Set refinement flags.   *
 *****************************************************************/
inline void AdaptiveBlock2D_List::setRefineAll(const int Flag) {
   int i; 
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      RefineFlag[i] = Flag;
   } /* endfor */
}

/********************************************************************
 * AdaptiveBlock2D_List::Nused -- Return the number of blocks used. *
 ********************************************************************/
inline int AdaptiveBlock2D_List::Nused(void) {
   int i, n = 0;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      if (Block[i].used) n++;
   } /* endfor */
   return n;
}

/*************************************************************
 * AdaptiveBlock2D_List -- Input-output operators.           *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const AdaptiveBlock2D_List &Blk_List) {
  int i;
  out_file << Blk_List.Nblk << "\n";
  for ( i = 0 ; i <= Blk_List.Nblk-1; ++i ) {
     out_file << Blk_List.Block[i] << "\n";
  } /* endfor */
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     AdaptiveBlock2D_List &Blk_List) {
  int i, n;
  in_file.setf(ios::skipws); in_file >> n; in_file.unsetf(ios::skipws);
  if (Blk_List.Block == NULL || Blk_List.Nblk != n) {
      Blk_List.deallocate(); Blk_List.allocate(n);
  } /* endif */
  for ( i = 0 ; i <= Blk_List.Nblk-1; ++i ) {
     in_file >> Blk_List.Block[i];
  } /* endfor */
  return (in_file);
}

/*************************************************************
 * AdaptiveBlock2D_List -- External subroutines.             *
 *************************************************************/

extern void Allocate_Message_Buffers(AdaptiveBlock2D_List &Blk_List,
                                     const int Number_of_Solution_Variables);

extern void Allocate_Message_Buffers_NoResChange(AdaptiveBlock2D_List &Blk_List,
                                                 const int Number_of_Solution_Variables);

extern void Allocate_Message_Buffers_ResChange(AdaptiveBlock2D_List &Blk_List,
                                               const int Number_of_Solution_Variables);

extern void Deallocate_Message_Buffers(AdaptiveBlock2D_List &Blk_List);

extern void Deallocate_Message_Buffers_NoResChange(AdaptiveBlock2D_List &Blk_List);

extern void Deallocate_Message_Buffers_ResChange(AdaptiveBlock2D_List &Blk_List);

extern int Exchange_Messages(AdaptiveBlock2D_List &Blk_List,
                             const int Number_of_Solution_Variables);

extern int Exchange_Messages_NoResChange(AdaptiveBlock2D_List &Blk_List,
                                         const int Number_of_Solution_Variables);

extern int Exchange_Messages_ResChange_FineToCoarse(AdaptiveBlock2D_List &Blk_List,
                                                    const int Number_of_Solution_Variables);

extern int Exchange_Messages_ResChange_CoarseToFine(AdaptiveBlock2D_List &Blk_List,
                                                    const int Number_of_Solution_Variables);

extern int Exchange_Messages_Fine_Grid_Solution_Information(AdaptiveBlock2D_List &Blk_List,
							    const int Number_of_Solution_Variables);

extern void Copy_Refinement_Flags(AdaptiveBlock2D_List &Blk_List_1,
                                  AdaptiveBlock2D_List &Blk_List_2);

/*******************************************************************
 * AdaptiveBlock2D -- Include templated message passing rountines. *
 *******************************************************************/

#ifndef _ADAPTIVEBLOCK2D_MESSAGEPASSING_INCLUDED
#include "AdaptiveBlock2D_MessagePassing.h"
#endif // _ADAPTIVEBLOCK2D_MESSAGEPASSING_INCLUDED

#endif /* _ADAPTIVEBLOCK_INCLUDED  */
