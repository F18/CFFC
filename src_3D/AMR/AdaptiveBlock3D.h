/* AdaptiveBlock.h:  Header file defining adaptive block classes. */

#ifndef _ADAPTIVEBLOCK3D_INCLUDED
#define _ADAPTIVEBLOCK3D_INCLUDED

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

#ifndef _BLOCK_ORIENTATION_INFO_INCLUDED
#include "BlockOrientationInfo.h"
#endif // _BLOCK_ORIENTATION_INFO_INCLUDED

/* Define the block in-use indicators. */

#define	ADAPTIVEBLOCK3D_USED                            1
#define	ADAPTIVEBLOCK3D_NOT_USED                        0

/* Define the maximum number of possible neighbouring
   adaptive blocks. */

#define	ADAPTIVEBLOCK3D_NUMBER_FACE_NEIGHBOURS_MAX      4
#define	ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX      2  
#define	ADAPTIVEBLOCK3D_NUMBER_CORNER_NEIGHBOURS_MAX    1

/* Define the block sector indicators for 
   block refinement and coarsening. */

#define	ADAPTIVEBLOCK3D_SECTOR_NONE                    -1
#define	ADAPTIVEBLOCK3D_SECTOR_SW                       0
#define	ADAPTIVEBLOCK3D_SECTOR_SE                       1
#define	ADAPTIVEBLOCK3D_SECTOR_NW                       2
#define	ADAPTIVEBLOCK3D_SECTOR_NE                       3
#define	ADAPTIVEBLOCK3D_SECTOR_TS                       4
#define	ADAPTIVEBLOCK3D_SECTOR_TW                       5
#define	ADAPTIVEBLOCK3D_SECTOR_TE                       6
#define	ADAPTIVEBLOCK3D_SECTOR_TN                       7
#define	ADAPTIVEBLOCK3D_SECTOR_TSW                      8
#define	ADAPTIVEBLOCK3D_SECTOR_TSE                      9
#define	ADAPTIVEBLOCK3D_SECTOR_TNW                      10
#define	ADAPTIVEBLOCK3D_SECTOR_TNE                      11
#define	ADAPTIVEBLOCK3D_SECTOR_BS                       12
#define	ADAPTIVEBLOCK3D_SECTOR_BW                       13
#define	ADAPTIVEBLOCK3D_SECTOR_BE                       14
#define	ADAPTIVEBLOCK3D_SECTOR_BN                       15
#define	ADAPTIVEBLOCK3D_SECTOR_BSW                      16
#define	ADAPTIVEBLOCK3D_SECTOR_BSE                      17
#define	ADAPTIVEBLOCK3D_SECTOR_BNW                      18
#define	ADAPTIVEBLOCK3D_SECTOR_BNE                      19

/* Define the block refinement indicators. */

#define	ADAPTIVEBLOCK3D_REFINE                          1
#define	ADAPTIVEBLOCK3D_NOCHANGE                        0
#define	ADAPTIVEBLOCK3D_COARSEN                        -1




/*******************************************************************************
 *
 * Special types
 *
 ******************************************************************************/
namespace BE
{
   enum Boundary_Elements{
      
      BSW = 0,
      SW = 1,
      TSW = 2,
      BW = 3,
      W =4,
      TW = 5,
      BNW = 6,
      NW = 7,
      TNW = 8,
      BS = 9,
      S = 10,
      TS = 11,
      B = 12,
      T = 14,
      BN = 15,
      N = 16,
      TN = 17,
      BSE = 18,
      SE = 19,
      TSE = 20,
      BE = 21,
      E = 22,
      TE = 23,
      BNE = 24,
      NE = 25,
      TNE = 26
   
   };
}


/* Define the classes. */

/********************************************************
 * Class: AdaptiveBlock3D_ResourceList                  *
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
class AdaptiveBlock3D_ResourceList{
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
    AdaptiveBlock3D_ResourceList(void) {
       ThisCPU = 0; Ncpu = 0; Nblk = 0; Ntotal = 0; Nused = 0; Nfree = 0;
       CPU = NULL; Block = NULL; 
    }

    AdaptiveBlock3D_ResourceList(const AdaptiveBlock3D_ResourceList &R) {
       ThisCPU = R.ThisCPU; Ncpu = R.Ncpu; Nblk = R.Nblk; Ntotal = R.Ntotal; 
       Nused = R.Nused; Nfree = R.Nfree; CPU = R.CPU; Block = R.Block;
    }

    /* Destructor. */
    // ~AdaptiveBlock3D_ResourceList(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // AdaptiveBlock3D_ResourceList operator = (const AdaptiveBlock3D_ResourceList &R);
    // Use automatically generated assignment operator.

    /* Allocate memory for the resource list. */
    void allocate(const int ncpu, const int nblk);

    /* Deallocate memory for the resource list. */
    void deallocate(void);

    /* Initialize resource list. */
    void initialize(void);

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

    static void Create_Block_Resource_List(AdaptiveBlock3D_ResourceList &List_of_Available_Blocks,
                                           const int Number_of_Processors,
                                           const int Number_of_Blocks_per_Processor);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, 
                                 const AdaptiveBlock3D_ResourceList &R);
    friend istream &operator >> (istream &in_file, 
                                 AdaptiveBlock3D_ResourceList &R);

};

/**********************************************************************
 * AdaptiveBlock3D_ResourceList::allocate -- Allocate memory.         *
 **********************************************************************/
inline void AdaptiveBlock3D_ResourceList::allocate(const int ncpu, 
                                                   const int nblk) {
   assert( ncpu > 0 && nblk > 0); 
   Ncpu = ncpu; Nblk = nblk; Ntotal = ncpu*nblk; Nused = 0; Nfree = Ntotal;
   CPU = new int[Ntotal]; Block = new int[Ntotal];
}

/**********************************************************************
 * AdaptiveBlock3D_ResourceList::deallocate -- Deallocate memory.     *
 **********************************************************************/
inline void AdaptiveBlock3D_ResourceList::deallocate(void) {
   Ncpu = 0; Nblk = 0; Ntotal = 0; Nused = 0; Nfree = 0;
   delete []CPU; CPU = NULL; delete []Block; Block = NULL; 
}

/**********************************************************************
 * AdaptiveBlock3D_ResourceList::initialize -- Initialization.        *
 **********************************************************************/
inline void AdaptiveBlock3D_ResourceList::initialize(void) {
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

/**********************************************************************
 * AdaptiveBlock3D_ResourceList::nextCPU -- Next CPU number.          *
 **********************************************************************/
inline int AdaptiveBlock3D_ResourceList::nextCPU(void) {
   return(CPU[Nused]);
}

/**********************************************************************
 * AdaptiveBlock3D_ResourceList::nextBlock -- Next Block number.      *
 **********************************************************************/
inline int AdaptiveBlock3D_ResourceList::nextBlock(void) {
   return(Block[Nused]);
}

/**********************************************************************
 * AdaptiveBlock3D_ResourceList::update_next -- Update next.          *
 **********************************************************************/
inline void AdaptiveBlock3D_ResourceList::update_next(void) {
   Nused += 1; Nfree -= 1;
}

/**********************************************************************
 * AdaptiveBlock3D_ResourceList::returnCPU -- Return CPU number.      *
 **********************************************************************/
inline void AdaptiveBlock3D_ResourceList::returnCPU(const int cpu) {
   CPU[Nused-1] = cpu;
}

/**********************************************************************
 * AdaptiveBlock3D_ResourceList::returnBlock -- Return Block number.  *
 **********************************************************************/
inline void AdaptiveBlock3D_ResourceList::returnBlock(const int block) {
   Block[Nused-1] = block;
}

/*************************************************************************
 * AdaptiveBlock3D_ResourceList::update_return -- Update returned block. *
 *************************************************************************/
inline void AdaptiveBlock3D_ResourceList::update_return(void) {
   Nused -= 1; Nfree += 1;
}

/*******************************************************************
 * AdaptiveBlock3D_ResourceList -- Input-output operators.         *
 *******************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const AdaptiveBlock3D_ResourceList &R) {
  int i;
  for ( i = 0; i <= R.Ntotal-1; ++i ) {
      out_file << " " << R.CPU[i] << " " << R.Block[i] << "\n";
  } /* endfor */
  return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             AdaptiveBlock3D_ResourceList &R) {
  int i;
  for ( i = 0; i <= R.Ntotal-1; ++i ) {
      in_file >> R.CPU[i] >> R.Block[i];
  } /* endfor */
  return (in_file);
}


/********************************************************
 * Class: AdaptiveBlock3D_Dimensions                    *
 *                                                      *
 * Member functions                                     *
 *         i    -- Return adaptive block dimension in   *
 *                 the i-direction.                     *
 *         j    -- Return adaptive block dimension in   *
 *                 the j-direction.                     *
 *         k    -- Return adaptive block dimension in   *
 *                 the k-direction.                     *
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
class AdaptiveBlock3D_Dimensions{
  private:
  public:
    int         i,j,k;  // Adaptive block dimensions.
    int         ghost;  // Adaptive block number of ghost 
                        // (halo or overlap) cells.
	                // Made public so can access them.
    
    /* Creation, copy, and assignment constructors. */
    AdaptiveBlock3D_Dimensions(void) {
       i = 0; j = 0; k = 0; ghost = 0;
    }

    AdaptiveBlock3D_Dimensions(const AdaptiveBlock3D_Dimensions &Blk_Dimen) {
       i = Blk_Dimen.i; j = Blk_Dimen.j; k = Blk_Dimen.k; ghost = Blk_Dimen.ghost;
    }

    AdaptiveBlock3D_Dimensions(const int i_dimen,
	                       const int j_dimen,
	                       const int k_dimen,
                               const int ghost_dimen) {
       i = i_dimen; j = j_dimen; k = k_dimen; ghost = ghost_dimen;
    }

    /* Destructor. */
    // ~AdaptiveBlock3D_Dimensions(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // AdaptiveBlock3D_Dimensions operator = (const AdaptiveBlock3D_Dimensions &Blk_Dimen);
    // Use automatically generated assignment operator.

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const AdaptiveBlock3D_Dimensions &Blk_Dimen);
    friend istream &operator >> (istream &in_file,
				 AdaptiveBlock3D_Dimensions &Blk_Dimen);

};

/*************************************************************
 * AdaptiveBlock3D_Dimensions -- Input-output operators.     *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const AdaptiveBlock3D_Dimensions &Blk_Dimen) {
  out_file << " " << Blk_Dimen.i << " " << Blk_Dimen.j << " " << Blk_Dimen.k<< " " << Blk_Dimen.ghost;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     AdaptiveBlock3D_Dimensions &Blk_Dimen) {
  in_file.setf(ios::skipws);
  in_file >> Blk_Dimen.i >> Blk_Dimen.j >> Blk_Dimen.k >> Blk_Dimen.ghost;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Class: AdaptiveBlock3D_Info                          *
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
class AdaptiveBlock3D_Info{
  private:
  public:
    int                            cpu; // Adaptive block processor number.
    int                         blknum; // Local adaptive block number.
    AdaptiveBlock3D_Dimensions   dimen; // Adaptive block dimensions. 
    int                         sector; // Adaptive block sector number.
    int                          level; // Adaptive block refinement level.
    Block_Orientation_Info       blkorient; // block orientation information
    
	                                // Made public so can access them.
    
    /* Creation, copy, and assignment constructors. */
    AdaptiveBlock3D_Info(void) {
       cpu = 0; blknum = 0; dimen.i = 0; dimen.j = 0; dimen.k = 0; dimen.ghost = 0; 
       sector = ADAPTIVEBLOCK3D_SECTOR_NONE; level = 0; 
    }

    AdaptiveBlock3D_Info(const AdaptiveBlock3D_Info &Blk_Info) {
       cpu = Blk_Info.cpu; blknum = Blk_Info.blknum;
       dimen = Blk_Info.dimen; sector = Blk_Info.sector; level = Blk_Info.level;
       blkorient = Blk_Info.blkorient;
       
    }

    AdaptiveBlock3D_Info(const int i_processor,
	                 const int i_block,
                         const int i_dimen,
			 const int j_dimen,
 			 const int k_dimen,
                         const int ghost_dimen,
                         const int i_sector,
	                 const int i_level,
                         const Block_Orientation_Info i_blkorient) {
       cpu = i_processor; blknum = i_block;
       dimen.i = i_dimen; dimen.j = j_dimen;dimen.k = k_dimen;  dimen.ghost = ghost_dimen;
       sector = i_sector; level = i_level; blkorient = i_blkorient;
       
    }

    /* Destructor. */
    // ~AdaptiveBlock3D_Info(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // AdaptiveBlock3D_Info operator = (const AdaptiveBlock3D_Info &Blk_Info);
    // Use automatically generated assignment operator.

    /* Reset the block information. */
    void reset(void);


    
    static void Broadcast_Adaptive_Block_Info(AdaptiveBlock3D_Info &Blk_Info);
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const AdaptiveBlock3D_Info &Blk_Info);
    friend istream &operator >> (istream &in_file,
				 AdaptiveBlock3D_Info &Blk_Info);

};

/*************************************************************
 * AdaptiveBlock3D_Info -- Reset the block information.      *
 *************************************************************/
inline void AdaptiveBlock3D_Info::reset(void) {
  cpu = 0; 
  blknum = 0; 
  dimen.i = 0; dimen.j = 0; dimen.k = 0; dimen.ghost = 0; 
  sector = ADAPTIVEBLOCK3D_SECTOR_NONE;
  level = 0;
  
}

/*************************************************************
 * AdaptiveBlock3D_Info -- Input-output operators.           *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const AdaptiveBlock3D_Info &Blk_Info) {
  out_file << " " << Blk_Info.cpu 
           << " " << Blk_Info.blknum << Blk_Info.dimen 
           << " " << Blk_Info.sector
           << " " << Blk_Info.level;
  out_file << " " << Blk_Info.blkorient;
  
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     AdaptiveBlock3D_Info &Blk_Info) {
  in_file.setf(ios::skipws);
  in_file >> Blk_Info.cpu >> Blk_Info.blknum;
  in_file.unsetf(ios::skipws);
  in_file >> Blk_Info.dimen;
  in_file.setf(ios::skipws);
  in_file >> Blk_Info.sector >> Blk_Info.level;
  in_file.unsetf(ios::skipws);
  in_file.setf(ios::skipws);
  in_file >> Blk_Info.blkorient;
  return (in_file);
}


/********************************************************
 * Class: AdaptiveBlock3D                               *
 *                                                      *
 * Member functions                                     *
 *      used     -- Return adaptive block usage.        *
 *   gblknum     -- Return global block number.         *
 *      info     -- Return adaptive block info.         *
 *        nT     -- Return number of Top neighbour      *
 *                  blocks.                             *
 *        nB     -- Return number of Bottom neighbour   *
 *                  blocks.                             *
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
 *       nTN     -- Return number of Top-north          *
 *                  corner neighbour blocks.            *
 *       nTS     -- Return number of Top-south          *
 *                  corner neighbour blocks.            *
 *       nTE     -- Return number of Top-east           *
 *                  corner neighbour blocks.            *
 *       nTW     -- Return number of Top-west           *
 *                  corner neighbour blocks.            *
 *       nBN     -- Return number of bottom-north       *
 *                  corner neighbour blocks.            *
 *       nBS     -- Return number of bottom-south       *
 *                  corner neighbour blocks.            *
 *       nBE     -- Return number of bottom-east        *
 *                  corner neighbour blocks.            *
 *       nBW     -- Return number of bottom-west        *
 *                  corner neighbour blocks.            *
 *       nTNW     -- Return number of Top- north-west   *
 *                  corner neighbour blocks.            *
 *       nTNE     -- Return number of Top- north-east   *
 *                  corner neighbour blocks.            *
 *       nTSE     -- Return number of Top- south-east   *
 *                  corner neighbour blocks.            *
 *       nTSW     -- Return number of Top- south-west   *
 *                  corner neighbour blocks.            *
 *       nBNW     -- Return number of Bottom-north-west *
 *                  corner neighbour blocks.            *
 *       nBNE     -- Return number of Bottom-north-east *
 *                  corner neighbour blocks.            *
 *       nBSE     -- Return number of Bottom-south-east *
 *                  corner neighbour blocks.            *
 *       nBSW     -- Return number of Bottom-south-west *
 *                  corner neighbour blocks.            *
 *     infoT     -- Return info for top neighbour       *
 *                  blocks.                             *
 *     infoB     -- Return info for bottom neighbour    *
 *                  blocks.                             *
 *     infoN     -- Return info for north neighbour     *
 *                  blocks.                             *
 *     infoS     -- Return info for south neighbour     *
 *                  blocks.                             *
 *     infoE     -- Return info for east neighbour      *
 *                  blocks.                             *
 *     infoW     -- Return info for west neighbour      *
 *                  blocks.                             *
 *    infoTN     -- Return info for top  - north        *
 *                  corner neighbour blocks.            *
 *    infoTS     -- Return info for top  - south        *
 *                  corner neighbour blocks.            *
 *    infoTE     -- Return info for top  - east         *
 *                  corner neighbour blocks.            *
 *    infoTW     -- Return info for top  - west         *
 *                  corner neighbour blocks.            *
 *    infoTNW     -- Return info for Top-north-west     *
 *                  corner neighbour blocks.            *
 *    infoTNE     -- Return info for Top-north-east     *
 *                  corner neighbour blocks.            *
 *    infoTSE     -- Return info for Top-south-east     *
 *                  corner neighbour blocks.            *
 *    infoTSW     -- Return info for Top-south-west     *
 *                  corner neighbour blocks.            *
 *    infoNW     -- Return info for north-west corner   *
 *                  neighbour blocks.                   *
 *    infoNE     -- Return info for north-east corner   *
 *                  neighbour blocks.                   *
 *    infoSE     -- Return info for south-east corner   *
 *                  neighbour blocks.                   *
 *    infoSW     -- Return info for south-west corner   *
 *                  neighbour blocks.                   * 
 *    infoBN     -- Return info for bottom-north        *
 *                  corner neighbour blocks.            *
 *    infoBS     -- Return info for bottom-south        *
 *                  corner neighbour blocks.            *
 *    infoBE     -- Return info for bottom-east         *
 *                  corner neighbour blocks.            *
 *    infoBW     -- Return info for bottom-west         *
 *                  corner neighbour blocks.            *
 *    infoBNW     -- Return info for bottom-north-west  *
 *                  corner neighbour blocks.            *
 *    infoBNE     -- Return info for bottom-north-east  *
 *                  corner neighbour blocks.            *
 *    infoBSE     -- Return info for bottom-south-east  *
 *                  corner neighbour blocks.            *
 *    infoBSW     -- Return info for bottom-south-west  *
 *                  corner neighbour blocks.            *
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
class AdaptiveBlock3D{
  private:
  public:
    int                    used;  // Adaptive block usage indicator.
    int                 gblknum;  // Global adaptive block number.
    AdaptiveBlock3D_Info   info;  // Adaptive block Info.
    int       nT,nB,nN,nS,nE,nW;  // Number of neighbouring blocks in each of the 26 directions
    int         nNW,nNE,nSE,nSW;  
    int         nTN,nTS,nTE,nTW;  
    int     nTNW,nTSW,nTNE,nTSE;  
    int         nBN,nBS,nBE,nBW;  
    int     nBNW,nBSW,nBNE,nBSE;  
    AdaptiveBlock3D_Info infoT[ADAPTIVEBLOCK3D_NUMBER_FACE_NEIGHBOURS_MAX],
                         infoB[ADAPTIVEBLOCK3D_NUMBER_FACE_NEIGHBOURS_MAX],
                         infoN[ADAPTIVEBLOCK3D_NUMBER_FACE_NEIGHBOURS_MAX],
                         infoS[ADAPTIVEBLOCK3D_NUMBER_FACE_NEIGHBOURS_MAX],
                         infoE[ADAPTIVEBLOCK3D_NUMBER_FACE_NEIGHBOURS_MAX],
                         infoW[ADAPTIVEBLOCK3D_NUMBER_FACE_NEIGHBOURS_MAX];
    AdaptiveBlock3D_Info infoNW[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX],
                         infoNE[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX],
                         infoSE[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX],
                         infoSW[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX];
    AdaptiveBlock3D_Info infoTN[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX],
                         infoTS[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX],
                         infoTE[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX],
                         infoTW[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX];
    AdaptiveBlock3D_Info infoTNW[ADAPTIVEBLOCK3D_NUMBER_CORNER_NEIGHBOURS_MAX],
                         infoTSW[ADAPTIVEBLOCK3D_NUMBER_CORNER_NEIGHBOURS_MAX],
                         infoTNE[ADAPTIVEBLOCK3D_NUMBER_CORNER_NEIGHBOURS_MAX],
                         infoTSE[ADAPTIVEBLOCK3D_NUMBER_CORNER_NEIGHBOURS_MAX];
    AdaptiveBlock3D_Info infoBN[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX],
                         infoBS[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX],
                         infoBE[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX],
                         infoBW[ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX];
    AdaptiveBlock3D_Info infoBNW[ADAPTIVEBLOCK3D_NUMBER_CORNER_NEIGHBOURS_MAX],
                         infoBSW[ADAPTIVEBLOCK3D_NUMBER_CORNER_NEIGHBOURS_MAX],
                         infoBNE[ADAPTIVEBLOCK3D_NUMBER_CORNER_NEIGHBOURS_MAX],
                         infoBSE[ADAPTIVEBLOCK3D_NUMBER_CORNER_NEIGHBOURS_MAX];
                                 // Block info of neighbouring blocks
 	                          // Made public so can access them.
    
    /* Creation, copy, and assignment constructors. */
    AdaptiveBlock3D(void) {
       used = ADAPTIVEBLOCK3D_NOT_USED; gblknum = 0;
       nN = 0; nS = 0; nE = 0; nW = 0; nT=0;nB=0;
       nNW = 0; nNE = 0; nSE = 0; nSW = 0;
       nTN = 0; nTS = 0; nTE = 0; nTW = 0;
       nTNW = 0; nTNE = 0; nTSE = 0; nTSW = 0;
       nBN = 0; nBS = 0; nBE = 0; nBW = 0;
       nBNW = 0; nBNE = 0; nBSE = 0; nBSW = 0;
    }

    AdaptiveBlock3D(const AdaptiveBlock3D &Blk) {
       int i;
       used = Blk.used; gblknum = Blk.gblknum; info = Blk.info;
       nT = Blk.nT; nB = Blk.nB; nN = Blk.nN; nS = Blk.nS; nE = Blk.nE; nW = Blk.nW;
       for (i = 0; i <= ADAPTIVEBLOCK3D_NUMBER_FACE_NEIGHBOURS_MAX-1; ++i) {
	 infoT[i] = Blk.infoT[i]; infoB[i] = Blk.infoB[i]; 
	 infoN[i] = Blk.infoN[i]; infoS[i] = Blk.infoS[i]; 
	 infoE[i] = Blk.infoE[i]; infoW[i] = Blk.infoW[i];
       } /* endfor */

       nNW = Blk.nNW; nNE = Blk.nNE; nSE = Blk.nSE; nSW = Blk.nSW;
       nTN = Blk.nTN; nTS = Blk.nTS; nTE = Blk.nTE; nTW = Blk.nTW;
       nTNW = Blk.nTNW; nTNE = Blk.nTNE; nTSE = Blk.nTSE; nTSW = Blk.nTSW;
       nBN = Blk.nBN; nBS = Blk.nBS; nBE = Blk.nBE; nBW = Blk.nBW;
       nBNW = Blk.nBNW; nBNE = Blk.nBNE; nBSE = Blk.nBSE; nBSW = Blk.nBSW;

       for (i = 0; i <= ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX-1; ++i) {
	 infoNW[i]=Blk.infoNW[i]; infoNE[i]=Blk.infoNE[i]; infoSE[i]=Blk.infoSE[i]; infoSW[i]=Blk.infoSW[i];
	 infoTN[i]=Blk.infoTN[i]; infoTS[i]=Blk.infoTS[i]; infoTE[i]=Blk.infoTE[i]; infoTW[i]=Blk.infoTW[i];
	 infoBN[i]=Blk.infoBN[i]; infoBS[i]=Blk.infoBS[i]; infoBE[i]=Blk.infoBE[i]; infoBW[i]=Blk.infoBW[i];
       }

       for (i = 0; i <= ADAPTIVEBLOCK3D_NUMBER_CORNER_NEIGHBOURS_MAX-1; ++i) {
	 infoTNW[i]=Blk.infoTNW[i];infoTNE[i]=Blk.infoTNE[i]; infoTSE[i]=Blk.infoTSE[i]; infoTSW[i]=Blk.infoTSW[i];
	 infoBNW[i]=Blk.infoBNW[i];infoBNE[i]=Blk.infoBNE[i]; infoBSE[i]=Blk.infoBSE[i]; infoBSW[i]=Blk.infoBSW[i];
       } /* endfor */
    }

    AdaptiveBlock3D(const AdaptiveBlock3D *Blk) {
       int i;
       used = Blk->used; gblknum = Blk->gblknum; info = Blk->info;
        nT = Blk->nT; nB = Blk->nB; nN = Blk->nN; nS = Blk->nS; nE = Blk->nE; nW = Blk->nW;
       for (i = 0; i <= ADAPTIVEBLOCK3D_NUMBER_FACE_NEIGHBOURS_MAX-1; ++i) {
	 infoT[i] = Blk->infoT[i]; infoB[i] = Blk->infoB[i];
	 infoN[i] = Blk->infoN[i]; infoS[i] = Blk->infoS[i];
	 infoE[i] = Blk->infoE[i]; infoW[i] = Blk->infoW[i];
       } /* endfor */
       nNW = Blk->nNW; nNE = Blk->nNE; nSE = Blk->nSE; nSW = Blk->nSW;
       for (i = 0; i <= ADAPTIVEBLOCK3D_NUMBER_EDGE_NEIGHBOURS_MAX-1; ++i) {
	 infoNW[i]=Blk->infoNW[i]; infoNE[i]=Blk->infoNE[i]; infoSE[i]=Blk->infoSE[i]; infoSW[i]=Blk->infoSW[i];
	 infoTN[i]=Blk->infoTN[i]; infoTS[i]=Blk->infoTS[i]; infoTE[i]=Blk->infoTE[i]; infoTW[i]=Blk->infoTW[i];
	 infoBN[i]=Blk->infoBN[i]; infoBS[i]=Blk->infoBS[i]; infoBE[i]=Blk->infoBE[i]; infoBW[i]=Blk->infoBW[i];
       }
       for (i = 0; i <= ADAPTIVEBLOCK3D_NUMBER_CORNER_NEIGHBOURS_MAX-1; ++i) {
	 infoTNW[i]=Blk->infoTNW[i];infoTNE[i]=Blk->infoTNE[i];infoTSE[i]=Blk->infoTSE[i];infoTSW[i]=Blk->infoTSW[i];
	 infoTNW[i]=Blk->infoTNW[i];infoTNE[i]=Blk->infoTNE[i];infoTSE[i]=Blk->infoTSE[i];infoTSW[i]=Blk->infoTSW[i];
       } /* endfor */
    }

    /* Destructor. */
    // ~AdaptiveBlock3D(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // AdaptiveBlock3D operator =(const AdaptiveBlock3D &Blk);
    // Use automatically generated assignment operator.

    static void Broadcast_Adaptive_Block(AdaptiveBlock3D &Blk);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const AdaptiveBlock3D &Blk);
    friend istream &operator >> (istream &in_file,
				 AdaptiveBlock3D &Blk);

};

/*************************************************************
 * AdaptiveBlock3D -- Input-output operators.                *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const AdaptiveBlock3D &Blk) {
  int i;
  out_file << " " << Blk.used << " " << Blk.gblknum << Blk.info;
  out_file << "\n " << Blk.nT;
  for (i = 0; i <= Blk.nT-1; ++i) {
     out_file << " " << Blk.infoT[i];
  } /* endfor */
  out_file << "\n " << Blk.nB;
  for (i = 0; i <= Blk.nB-1; ++i) {
     out_file << " " << Blk.infoB[i];
  } /* endfor */
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

  out_file << "\n " << Blk.nTN;
  for (i = 0; i <= Blk.nTN-1; ++i) {
     out_file << " " << Blk.infoTN[i];
  } /* endfor */
  out_file << "\n " << Blk.nTS;
  for (i = 0; i <= Blk.nTS-1; ++i) {
     out_file << " " << Blk.infoTS[i];
  } /* endfor */
  out_file << "\n " << Blk.nTE;
  for (i = 0; i <= Blk.nTE-1; ++i) {
     out_file << " " << Blk.infoTE[i];
  } /* endfor */
  out_file << "\n " << Blk.nTW;
  for (i = 0; i <= Blk.nTW-1; ++i) {
     out_file << " " << Blk.infoTW[i];
  } /* endfor */


  out_file << "\n " << Blk.nTNW;
  for (i = 0; i <= Blk.nTNW-1; ++i) {
     out_file << " " << Blk.infoTNW[i];
  } /* endfor */
  out_file << "\n " << Blk.nTNE;
  for (i = 0; i <= Blk.nTNE-1; ++i) {
     out_file << " " << Blk.infoTNE[i];
  } /* endfor */
  out_file << "\n " << Blk.nTSE;
  for (i = 0; i <= Blk.nTSE-1; ++i) {
     out_file << " " << Blk.infoTSE[i];
  } /* endfor */
  out_file << "\n " << Blk.nTSW;
  for (i = 0; i <= Blk.nTSW-1; ++i) {
     out_file << " " << Blk.infoTSW[i];
  } /* endfor */


  out_file << "\n " << Blk.nBN;
  for (i = 0; i <= Blk.nBN-1; ++i) {
     out_file << " " << Blk.infoBN[i];
  } /* endfor */
  out_file << "\n " << Blk.nBS;
  for (i = 0; i <= Blk.nBS-1; ++i) {
     out_file << " " << Blk.infoBS[i];
  } /* endfor */
  out_file << "\n " << Blk.nBE;
  for (i = 0; i <= Blk.nBE-1; ++i) {
     out_file << " " << Blk.infoBE[i];
  } /* endfor */
  out_file << "\n " << Blk.nBW;
  for (i = 0; i <= Blk.nBW-1; ++i) {
     out_file << " " << Blk.infoBW[i];
  } /* endfor */


  out_file << "\n " << Blk.nBNW;
  for (i = 0; i <= Blk.nBNW-1; ++i) {
     out_file << " " << Blk.infoBNW[i];
  } /* endfor */
  out_file << "\n " << Blk.nBNE;
  for (i = 0; i <= Blk.nBNE-1; ++i) {
     out_file << " " << Blk.infoBNE[i];
  } /* endfor */
  out_file << "\n " << Blk.nBSE;
  for (i = 0; i <= Blk.nBSE-1; ++i) {
     out_file << " " << Blk.infoBSE[i];
  } /* endfor */
  out_file << "\n " << Blk.nBSW;
  for (i = 0; i <= Blk.nBSW-1; ++i) {
     out_file << " " << Blk.infoBSW[i];
  } /* endfor */

  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     AdaptiveBlock3D &Blk) {
  int i;
  in_file.setf(ios::skipws); in_file >> Blk.used >> Blk.gblknum; 
  in_file.unsetf(ios::skipws);
  in_file >> Blk.info;
  in_file.setf(ios::skipws); in_file >> Blk.nT; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nT-1; ++i) {
     in_file >> Blk.infoT[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nB; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nB-1; ++i) {
     in_file >> Blk.infoB[i];
  } /* endfor */
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

 in_file.setf(ios::skipws); in_file >> Blk.nTN; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nTN-1; ++i) {
     in_file >> Blk.infoTN[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nTS; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nTS-1; ++i) {
     in_file >> Blk.infoTS[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nTE; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nTE-1; ++i) {
     in_file >> Blk.infoTE[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nTW; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nTW-1; ++i) {
     in_file >> Blk.infoTW[i];
  } /* endfor */

 in_file.setf(ios::skipws); in_file >> Blk.nTNW; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nTNW-1; ++i) {
     in_file >> Blk.infoTNW[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nTNE; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nTNE-1; ++i) {
     in_file >> Blk.infoTNE[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nTSE; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nTSE-1; ++i) {
     in_file >> Blk.infoTSE[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nTSW; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nTSW-1; ++i) {
     in_file >> Blk.infoTSW[i];
  } /* endfor */


 in_file.setf(ios::skipws); in_file >> Blk.nBN; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nBN-1; ++i) {
     in_file >> Blk.infoBN[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nBS; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nBS-1; ++i) {
     in_file >> Blk.infoBS[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nBE; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nBE-1; ++i) {
     in_file >> Blk.infoBE[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nBW; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nBW-1; ++i) {
     in_file >> Blk.infoBW[i];
  } /* endfor */


 in_file.setf(ios::skipws); in_file >> Blk.nBNW; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nBNW-1; ++i) {
     in_file >> Blk.infoBNW[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nBNE; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nBNE-1; ++i) {
     in_file >> Blk.infoBNE[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nBSE; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nBSE-1; ++i) {
     in_file >> Blk.infoBSE[i];
  } /* endfor */
  in_file.setf(ios::skipws); in_file >> Blk.nBSW; in_file.unsetf(ios::skipws);
  for (i = 0; i <= Blk.nBSW-1; ++i) {
     in_file >> Blk.infoBSW[i];
  } /* endfor */


  return (in_file);
}

/**************************************************************
 * Class: AdaptiveBlock3D_List                                *
 *                                                            *
 * Member functions                                           *
 *     ThisCPU -- Return the global MPI number or             *
 *                rank for the CPU that is executing          *
 *                the task or process.                        *
 *        Nblk -- Return number of local adaptive             *
 *                blocks in the list.                         *
 *       Block -- Return list of adaptive blocks.             *
 *       RefineFlag -- Return list of mesh refinement         *
 *                     flags for the adaptive blocks.         *
 *       message_noreschange_sendbuf                          *
 *       message_noreschange_recbuf                           *
 *        -- Return message passing send and receive          *
 *           buffers for block faces and corners with         *
 *           no block resolution (refinement) change.         *
 *    allocate -- Allocate memory for block list.             *
 *  deallocate -- Deallocate memory for block list.           *
 *  nochangeAll-- Sets the mesh refinement flags to           *
 *                force no refinement or coarsening of        *
 *                adaptive blocks (default).                  *
 *  refineAll  -- Sets the mesh refinement flags to           *
 *                force refinement of all adaptive            *
 *                blocks.                                     *
 *  coarsenAll -- Sets the mesh refinement flags to           *
 *                force coasening of all adaptive             *
 *                blocks.                                     *
 *  setRefineAll -- Sets the mesh refinement flags of         *
 *                  all blocks to specified value.            *
 *                                                            *
 * Member operators                                           *
 *      B -- adaptive block list                              *
 *                                                            *
 * B = B;                                                     *
 * cout << B; (output function)                               *
 * cin  >> B; (input function)                                *
 *                                                            *
 **************************************************************/
class AdaptiveBlock3D_List{
  private:
  public:
    int              ThisCPU; // Number or rank of CPU executing the
                              // current process.
    int                 Nblk; // Number of local adaptive
                              // blocks in the list.
    AdaptiveBlock3D   *Block; // List of adaptive blocks.
    int          *RefineFlag; // Block refinement (coarsening) flag.
   
    // Message passing send and receive buffers
    // for block faces and corners with no block
    // resolution (refinement) change.
    
    double   ***message_noreschange_sendbuf,
       ***message_noreschange_recbuf;

 
                              // Message passing send and receive buffers
                              // for block faces and corners with block
                              // resolution (refinement) change.
                              // Made public so can access them.
    
    /* Creation, copy, and assignment constructors. */
    AdaptiveBlock3D_List(void) {
       ThisCPU = 0; Nblk = 0; Block = NULL; RefineFlag = NULL;
       // No resolution change message buffers.
       
       message_noreschange_sendbuf = NULL;
       message_noreschange_recbuf = NULL;
       
    }

    AdaptiveBlock3D_List(const AdaptiveBlock3D_List &Blk_List) {
       ThisCPU = Blk_List.ThisCPU; Nblk = Blk_List.Nblk; 
       Block = Blk_List.Block; RefineFlag = Blk_List.RefineFlag;
       // No resolution change message buffers.
       //******************************************************
       message_noreschange_sendbuf = Blk_List.message_noreschange_sendbuf;
       message_noreschange_recbuf = Blk_List.message_noreschange_recbuf;
       
    }

    /* Destructor. */
    // ~AdaptiveBlock3D_List(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // AdaptiveBlock3D_List operator = (const AdaptiveBlock3D_List &Blk_List);
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

    static void Allocate_Message_Buffers(AdaptiveBlock3D_List &Blk_List,
                                         const int Number_of_Solution_Variables);

    static void Allocate_Message_Buffers_NoResChange(AdaptiveBlock3D_List &Blk_List,
                                                     const int Number_of_Solution_Variables);
    
    static void Allocate_Message_Buffers_ResChange(AdaptiveBlock3D_List &Blk_List,
                                                   const int Number_of_Solution_Variables);
    
    static void Deallocate_Message_Buffers(AdaptiveBlock3D_List &Blk_List);
    
    static void Deallocate_Message_Buffers_NoResChange(AdaptiveBlock3D_List &Blk_List);
    
    static void Deallocate_Message_Buffers_ResChange(AdaptiveBlock3D_List &Blk_List);
    
    static int Exchange_Messages(AdaptiveBlock3D_List &Blk_List,
                                 const int Number_of_Solution_Variables);
    
    static int Exchange_Messages_NoResChange(AdaptiveBlock3D_List &Blk_List,
                                             const int Number_of_Solution_Variables);
    
    static int Exchange_Messages_ResChange_FineToCoarse(AdaptiveBlock3D_List &Blk_List,
                                                        const int Number_of_Solution_Variables);
    
    static int Exchange_Messages_ResChange_CoarseToFine(AdaptiveBlock3D_List &Blk_List,
                                                        const int Number_of_Solution_Variables);
    
    static  void Copy_Refinement_Flags(AdaptiveBlock3D_List &Blk_List_1,
                                       AdaptiveBlock3D_List &Blk_List_2);
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const AdaptiveBlock3D_List &Blk_List);
    friend istream &operator >> (istream &in_file,
				 AdaptiveBlock3D_List &Blk_List);

};

/*************************************************************
 * AdaptiveBlock3D_List::allocate -- Allocate memory.        *
 *************************************************************/
inline void AdaptiveBlock3D_List::allocate(const int N) {
   int i, j; assert( N > 0 ); Nblk = N;
   Block = new AdaptiveBlock3D[Nblk]; 
   RefineFlag = new int[Nblk]; nochangeAll();
 
   // No resolution change message buffers.
   
   message_noreschange_sendbuf = new double**[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
       message_noreschange_sendbuf[i] = new double*[27];
      for (j = 0; j < 27 ; ++j) {
         message_noreschange_sendbuf[i][j] = new double[1];
         message_noreschange_sendbuf[i][j][0] = ZERO;
      } /* endfor */
   } /* endfor */

   message_noreschange_recbuf = new double**[Nblk];
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      message_noreschange_recbuf[i] = new double*[27];
      for (j = 0; j < 27 ; ++j) {
         message_noreschange_recbuf[i][j] = new double[1];
         message_noreschange_recbuf[i][j][0] = ZERO;
      } /* endfor */
   } /* endfor */

}

/*************************************************************
 * AdaptiveBlock3D_List::deallocate -- Deallocate memory.    *
 *************************************************************/
inline void AdaptiveBlock3D_List::deallocate(void) {
   int i, j; delete []Block; Block = NULL;
   delete []RefineFlag; RefineFlag = NULL;
 
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      for (j = 0; j < 27 ; ++j) {
         delete []message_noreschange_sendbuf[i][j];
         message_noreschange_sendbuf[i][j] = NULL;
      } /* endfor */
      delete []message_noreschange_sendbuf[i];
      message_noreschange_sendbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_sendbuf;
   message_noreschange_sendbuf = NULL;

   for ( i = 0; i <= Nblk-1 ; ++i ) {
      for (j = 0; j < 27 ; ++j) {
         delete []message_noreschange_recbuf[i][j];
         message_noreschange_recbuf[i][j] = NULL;
      } /* endfor */
      delete []message_noreschange_recbuf[i];
      message_noreschange_recbuf[i] = NULL;
   } /* endfor */
   delete []message_noreschange_recbuf;
   message_noreschange_recbuf = NULL;


   // Reset the number of solution blocks.
   Nblk = 0;
}

/*****************************************************************
 * AdaptiveBlock3D_List::nochangeAll -- Set no refinement flags. *
 *****************************************************************/
inline void AdaptiveBlock3D_List::nochangeAll(void) {
   int i; 
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      RefineFlag[i] = ADAPTIVEBLOCK3D_NOCHANGE;
   } /* endfor */
}

/*****************************************************************
 * AdaptiveBlock3D_List::refineAll -- Set refinement flags.      *
 *****************************************************************/
inline void AdaptiveBlock3D_List::refineAll(void) {
   int i; 
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      RefineFlag[i] = ADAPTIVEBLOCK3D_REFINE;
   } /* endfor */
}

/*****************************************************************
 * AdaptiveBlock3D_List::coarsenAll -- Set coarsening flags.     *
 *****************************************************************/
inline void AdaptiveBlock3D_List::coarsenAll(void) {
   int i; 
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      RefineFlag[i] = ADAPTIVEBLOCK3D_COARSEN;
   } /* endfor */
}

/*****************************************************************
 * AdaptiveBlock3D_List::setRefineAll -- Set refinement flags.   *
 *****************************************************************/
inline void AdaptiveBlock3D_List::setRefineAll(const int Flag) {
   int i; 
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      RefineFlag[i] = Flag;
   } /* endfor */
}

/********************************************************************
 * AdaptiveBlock3D_List::Nused -- Return the number of blocks used. *
 ********************************************************************/
inline int AdaptiveBlock3D_List::Nused(void) {
   int i, n = 0;
   for ( i = 0; i <= Nblk-1 ; ++i ) {
      if (Block[i].used) n++;
   } /* endfor */
   return n;
}

/*************************************************************
 * AdaptiveBlock3D_List -- Input-output operators.           *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const AdaptiveBlock3D_List &Blk_List) {
  int i;
  out_file << Blk_List.Nblk << "\n";
  for ( i = 0 ; i <= Blk_List.Nblk-1; ++i ) {
     out_file << Blk_List.Block[i] << "\n";
  } /* endfor */
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     AdaptiveBlock3D_List &Blk_List) {
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


/*******************************************************************
 * AdaptiveBlock3D -- Include templated message passing rountines. *
 *******************************************************************/

#ifndef _ADAPTIVEBLOCK3D_MESSAGEPASSING_INCLUDED
#include "AdaptiveBlock3D_MessagePassing.h"
#endif // _ADAPTIVEBLOCK3D_MESSAGEPASSING_INCLUDED

#endif /* _ADAPTIVEBLOCK_INCLUDED  */
