/* Octree.h:  Header file defining octree adaptive blocks 
              hierarchical data structure. */

#ifndef _OCTREE_INCLUDED
#define _OCTREE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include math macro and adaptive block header files. */

#ifndef _ADAPTIVEBLOCK_INCLUDED
#include "AdaptiveBlock3D.h"
#endif // _ADAPTIVEBLOCK_INCLUDED

/* Define the classes. */

#define OCTREE_SECTOR_NONE                        -1000
#define	OCTREE_SECTOR_BSW                             0
#define	OCTREE_SECTOR_BSE                             1
#define	OCTREE_SECTOR_BNW                             2
#define	OCTREE_SECTOR_BNE                             3
#define	OCTREE_SECTOR_TSW                             4
#define	OCTREE_SECTOR_TSE                             5
#define	OCTREE_SECTOR_TNW                             6
#define	OCTREE_SECTOR_TNE                             7

/* Define the octree search direction masks for finding
   nearest neighbours. */

#define	OCTREE_DIRECTION_MASK_TOP                    1
#define	OCTREE_DIRECTION_MASK_BOTTOM                 2
#define	OCTREE_DIRECTION_MASK_NORTH                  3
#define	OCTREE_DIRECTION_MASK_SOUTH                  4
#define	OCTREE_DIRECTION_MASK_EAST                   5
#define	OCTREE_DIRECTION_MASK_WEST                   6
#define	OCTREE_DIRECTION_MASK_NORTHWEST              7
#define	OCTREE_DIRECTION_MASK_NORTHEAST              8
#define	OCTREE_DIRECTION_MASK_SOUTHEAST              9
#define	OCTREE_DIRECTION_MASK_SOUTHWEST              10
#define	OCTREE_DIRECTION_MASK_TOPNORTH               11
#define	OCTREE_DIRECTION_MASK_TOPSOUTH               12
#define	OCTREE_DIRECTION_MASK_TOPEAST                13
#define	OCTREE_DIRECTION_MASK_TOPWEST                14
#define	OCTREE_DIRECTION_MASK_BOTTOMNORTH            15
#define	OCTREE_DIRECTION_MASK_BOTTOMSOUTH            16
#define	OCTREE_DIRECTION_MASK_BOTTOMEAST             17
#define	OCTREE_DIRECTION_MASK_BOTTOMWEST             18
#define	OCTREE_DIRECTION_MASK_TOPNORTHWEST           19
#define	OCTREE_DIRECTION_MASK_TOPNORTHEAST           20
#define	OCTREE_DIRECTION_MASK_TOPSOUTHEAST           21
#define	OCTREE_DIRECTION_MASK_TOPSOUTHWEST           22
#define	OCTREE_DIRECTION_MASK_BOTTOMNORTHWEST        23
#define	OCTREE_DIRECTION_MASK_BOTTOMNORTHEAST        24
#define	OCTREE_DIRECTION_MASK_BOTTOMSOUTHEAST        25
#define	OCTREE_DIRECTION_MASK_BOTTOMSOUTHWEST        26

/* Define the classes. */

/********************************************************
 * Class: OctreeBlock                                   *
 *                                                      *
 * Member functions                                     *
 *  block       -- Return octree adaptive block.        *
 *  parent_ptr  -- Return pointer to parent of octree   *
 *                 adaptive block.                      *
 *  childNW_ptr -- Return pointer to north-west child   *
 *                 of octree adaptive block.            *
 *  childNE_ptr -- Return pointer to north-east child   *
 *                 of octree adaptive block.            *
 *  childSW_ptr -- Return pointer to south-west child   *
 *                 of octree adaptive block.            *
 *  childSE_ptr -- Return pointer to south-east child   *
 *                 of octree adaptive block.            *
 *  child_ptr   -- Return pointer to specified child of *
 *                 octree adaptive block.               *
 *  search_dir  -- Returns the search direction given   *
 *                 the specified search direction mask. *
 *  sibling     -- Returns flag indicating whether or   *
 *                 not there is a sibling in the        *
 *                 specified search direction.          *
 *  neighbour_ptr -- Returns pointer to neighbouring    *
 *                 octree adaptive block found in the   *
 *                 specified search direction.          *
 *  read        -- Reads in the octree block and then   *
 *                 recursively descends the subtree and *
 *                 reads siblings and their siblings,   *
 *                 etc...                               *
 *  write       -- Writes out octree block and then     *
 *                 recursively descends the subtree and *
 *                 writes out siblings and their        *
 *                 siblings, etc...                     *
 *  broadcast   -- Broadcasts the octree block in a     *
 *                 recursive manner.                    *
 *  maxRefinementLevel -- Return the maximum refinement *
 *                        level of all blocks in        *
 *                        subtree.                      *
 *                                                      *
 * Member operators                                     *
 *      B -- element in a octree hierarchical data      *
 *           structure                                  *
 *                                                      *
 * B = B;                                               *
 * cout << B; (output function)                         *
 * cin  >> B; (input function)                          *
 *                                                      *
 ********************************************************/
class OctreeBlock{
  private:
  public:
    AdaptiveBlock3D        block;  // Pointer to octree adaptive block.
    OctreeBlock      *parent_ptr;  // Pointer to parent of adaptive block.
    OctreeBlock     *childTSW_ptr, // Pointers to children of
                    *childTSE_ptr, // adaptive block.
                    *childTNW_ptr, //
                    *childTNE_ptr, //
                    *childBSW_ptr, // Pointers to children of
                    *childBSE_ptr, // adaptive block.
                    *childBNW_ptr, //
                    *childBNE_ptr; // 

    /* Creation, copy, and assignment constructors. */
    OctreeBlock(void) {
      parent_ptr = NULL; 
      childTSW_ptr = NULL;
      childTSE_ptr = NULL;
      childTNW_ptr = NULL;
      childTNE_ptr = NULL;
      childBSW_ptr = NULL; 
      childBSE_ptr = NULL; 
      childBNW_ptr = NULL;
      childBNE_ptr = NULL;	
    }

    OctreeBlock(const OctreeBlock &Block) {
       block = Block.block; parent_ptr = Block.parent_ptr; 
       childTSW_ptr = Block.childTSW_ptr;
       childTSE_ptr = Block.childTSE_ptr;
       childTNW_ptr = Block.childTNW_ptr;
       childTNE_ptr = Block.childTNE_ptr;
       childBSW_ptr = Block.childBSW_ptr; 
       childBSE_ptr = Block.childBSE_ptr; 
       childBNW_ptr = Block.childBNW_ptr;
       childBNE_ptr = Block.childBNE_ptr;	
    }

    /* Destructor. */
    // ~OctreeBlock(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // OctreeBlock operator = (const OctreeBlock &Block);
    // Use automatically generated assignment operator.

    /* Pointer to child. */
    OctreeBlock *child_ptr(const int Sector);

    /* Search direction. */
    AdaptiveBlock3D_Dimensions search_dir(const int Search_Dir_Mask);

    /* Sibling in specified search direction? */
    int sibling(const int Search_Dir_Mask);

    /* Pointer to neighbour. */
    OctreeBlock *neighbour_ptr(const int Search_Dir_Mask);

    /* Read octree block (recursive). */
    void read(istream &in_file);
    void read(istream &in_file,
              AdaptiveBlock3D_ResourceList &List_of_Available_Blocks);

    /* Write octree block (recursive). */
    void write(ostream &out_file);

   /* Broadcast octree block (recursive). */
    void broadcast(void);
    void broadcast(AdaptiveBlock3D_ResourceList &List_of_Available_Blocks);

    /* Determine maximum refinement level. */
    int maxRefinementLevel(int &Max_Level);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const OctreeBlock &Block);
    friend istream &operator >> (istream &in_file,
				 OctreeBlock &Block);
};

/*************************************************************
 * OctreeBlock::child_ptr -- Return children of block.     *
 *************************************************************/
inline OctreeBlock *OctreeBlock::child_ptr(const int Sector) {
  switch(Sector) {
    case OCTREE_SECTOR_NONE :
      return (NULL);
      break;
    case OCTREE_SECTOR_TNW :
      return (childTNW_ptr);
      break;
    case OCTREE_SECTOR_TNE :
      return (childTNE_ptr);
      break;
    case OCTREE_SECTOR_TSE :
      return (childTSE_ptr);
      break;
    case OCTREE_SECTOR_TSW :
      return (childTSW_ptr);
      break;
    case OCTREE_SECTOR_BNW :
      return (childBNW_ptr);
      break;
    case OCTREE_SECTOR_BNE :
      return (childBNE_ptr);
      break;
    case OCTREE_SECTOR_BSE :
      return (childBSE_ptr);
      break;
    case OCTREE_SECTOR_BSW :
    default:
      return (childBSW_ptr);
      break;
  } /* endswitch */
}

/*************************************************************
 * OctreeBlock::search_dir -- Determine search direction.    *
 *************************************************************/
inline AdaptiveBlock3D_Dimensions OctreeBlock::search_dir(const int Search_Dir_Mask) {
  int i, j, k;
  switch(Search_Dir_Mask) {
    case OCTREE_DIRECTION_MASK_EAST :
      i = 1; j = 0; k= 0;
      break;
    case OCTREE_DIRECTION_MASK_WEST :
      i = -1; j = 0; k= 0;
      break;
    case OCTREE_DIRECTION_MASK_NORTH :
      i = 0; j = 2; k= 0;
      break;
    case OCTREE_DIRECTION_MASK_SOUTH :
      i = 0; j = -2; k= 0;
      break;
    case OCTREE_DIRECTION_MASK_TOP :
      i = 0; j = 0; k= 4;
      break;
    case OCTREE_DIRECTION_MASK_BOTTOM :
      i = 0; j = 0; k= -4;
      break;
    case OCTREE_DIRECTION_MASK_SOUTHWEST :
      i = -1; j = -2; k= 0;
      break;
    case OCTREE_DIRECTION_MASK_SOUTHEAST :
      i = 1; j = -2; k= 0;
      break;
    case OCTREE_DIRECTION_MASK_NORTHWEST :
      i = -1; j = 2; k= 0;
      break;
    case OCTREE_DIRECTION_MASK_NORTHEAST :
      i = 1; j = 2; k= 0;
      break;

     case OCTREE_DIRECTION_MASK_TOPNORTH :
     i = 0; j = 2; k= 4;
      break;
    case OCTREE_DIRECTION_MASK_TOPSOUTH :
     i = 0; j = -2; k= 4;
      break;
    case OCTREE_DIRECTION_MASK_TOPEAST :
      i = 1; j = 0; k= 4;
      break;
   case OCTREE_DIRECTION_MASK_TOPWEST :
      i = -1; j = 0; k= 4;
      break;
    case OCTREE_DIRECTION_MASK_BOTTOMNORTH :
      i = 0; j = 2; k= -4;
      break;
    case OCTREE_DIRECTION_MASK_BOTTOMSOUTH :
      i = 0; j = -2; k= -4;
      break;
   case OCTREE_DIRECTION_MASK_BOTTOMEAST :
      i = 1; j = 0; k= -4;
      break;
    case OCTREE_DIRECTION_MASK_BOTTOMWEST :
      i = -1; j = 0; k= -4;
      break;
    case OCTREE_DIRECTION_MASK_TOPNORTHEAST :
      i = 1; j = 2; k= 4;
      break;
    case OCTREE_DIRECTION_MASK_TOPNORTHWEST :
      i = 1; j = -2; k= 4;
      break;
    case OCTREE_DIRECTION_MASK_TOPSOUTHEAST :
      i = -1; j = 2; k= 4;
      break;
    case OCTREE_DIRECTION_MASK_TOPSOUTHWEST :
      i = -1; j = -2; k= 4;
      break;
    case OCTREE_DIRECTION_MASK_BOTTOMNORTHEAST :
      i = 1; j = 2; k= -4;
      break;
    case OCTREE_DIRECTION_MASK_BOTTOMNORTHWEST :
      i = 1; j = -2; k= -4;
      break;
    case OCTREE_DIRECTION_MASK_BOTTOMSOUTHEAST :
      i = -1; j = 2; k= -4;
      break;
    case OCTREE_DIRECTION_MASK_BOTTOMSOUTHWEST :
      i = -1; j = -2; k= -4;
      break;
    default:
      i = 0; j = 0; k= 0;
      break;
  } /* endswitch */
  return (AdaptiveBlock3D_Dimensions(i, j, k,0));
}

/*************************************************************
 * OctreeBlock::sibling -- Sibling in direction?           *
 *************************************************************/
inline int OctreeBlock::sibling(const int Search_Dir_Mask) {
  int i_sibling, j_sibling, k_sibling;
  AdaptiveBlock3D_Dimensions search_direction;
  if (block.info.sector != ADAPTIVEBLOCK3D_SECTOR_NONE) {
     search_direction = search_dir(Search_Dir_Mask);
     if (search_direction.i >= 0) {
        i_sibling = !(block.info.sector & search_direction.i);
     } else {
        i_sibling = !(!(block.info.sector & (-search_direction.i)));
     } /* endif */
     if (search_direction.j >= 0) {
        j_sibling = !(block.info.sector & search_direction.j);
     } else {
        j_sibling = !(!(block.info.sector & (-search_direction.j)));
     } /* endif */
     if (search_direction.k >= 0) {
        k_sibling = !(block.info.sector & search_direction.k);
     } else {
        k_sibling = !(!(block.info.sector & (-search_direction.k)));
     } /* endif */
  } else {
     i_sibling = 0;
     j_sibling = 0;
     k_sibling = 0;
  } /* endif */
  return (i_sibling && j_sibling && k_sibling);  //I'm not sure if this logic is correct
}


/*************************************************************
 * OctreeBlock::neighbour_ptr -- Return neighbouring block.*
 *************************************************************/
inline OctreeBlock *OctreeBlock::neighbour_ptr(const int Search_Dir_Mask) {
  OctreeBlock *block_ptr; AdaptiveBlock3D_Dimensions search_direction;
  return (block_ptr);
}

/*************************************************************
 * OctreeBlock::read -- Read octree block (recursive).       *
 *************************************************************/
inline void OctreeBlock::read(istream &in_file) {
  int number_of_children; 
  in_file >> block;
  in_file.setf(ios::skipws);  in_file >> number_of_children;  
  in_file.unsetf(ios::skipws);
  if (!block.used && number_of_children > 0) {
     if (childBSW_ptr != NULL) { delete childBSW_ptr; childBSW_ptr = NULL; }
     childBSW_ptr = new OctreeBlock;
     childBSW_ptr->parent_ptr = this;
     childBSW_ptr->read(in_file);
     if (childBSE_ptr != NULL) { delete childBSE_ptr; childBSE_ptr = NULL; }
     childBSE_ptr = new OctreeBlock;
     childBSE_ptr->parent_ptr = this;
     childBSE_ptr->read(in_file);
     if (childBNW_ptr != NULL) { delete childBNW_ptr; childBNW_ptr = NULL; }
     childBNW_ptr = new OctreeBlock;
     childBNW_ptr->parent_ptr = this;
     childBNW_ptr->read(in_file);
     if (childBNE_ptr != NULL) { delete childBNE_ptr; childBNE_ptr = NULL; }
     childBNE_ptr = new OctreeBlock;
     childBNE_ptr->parent_ptr = this;
     childBNE_ptr->read(in_file);
     if (childTSW_ptr != NULL) { delete childTSW_ptr; childTSW_ptr = NULL; }
     childTSW_ptr = new OctreeBlock;
     childTSW_ptr->parent_ptr = this;
     childTSW_ptr->read(in_file);
     if (childTSE_ptr != NULL) { delete childTSE_ptr; childTSE_ptr = NULL; }
     childTSE_ptr = new OctreeBlock;
     childTSE_ptr->parent_ptr = this;
     childTSE_ptr->read(in_file);
     if (childTNW_ptr != NULL) { delete childTNW_ptr; childTNW_ptr = NULL; }
     childTNW_ptr = new OctreeBlock;
     childTNW_ptr->parent_ptr = this;
     childTNW_ptr->read(in_file);
     if (childTNE_ptr != NULL) { delete childTNE_ptr; childTNE_ptr = NULL; }
     childTNE_ptr = new OctreeBlock;
     childTNE_ptr->parent_ptr = this;
     childTNE_ptr->read(in_file);
  } else {
     if (childBSW_ptr != NULL) { delete childBSW_ptr; childBSW_ptr = NULL; }
     if (childBSE_ptr != NULL) { delete childBSE_ptr; childBSE_ptr = NULL; }
     if (childBNW_ptr != NULL) { delete childBNW_ptr; childBNW_ptr = NULL; }
     if (childBNE_ptr != NULL) { delete childBNE_ptr; childBNE_ptr = NULL; }
     if (childTSW_ptr != NULL) { delete childTSW_ptr; childTSW_ptr = NULL; }
     if (childTSE_ptr != NULL) { delete childTSE_ptr; childTSE_ptr = NULL; }
     if (childTNW_ptr != NULL) { delete childTNW_ptr; childTNW_ptr = NULL; }
     if (childTNE_ptr != NULL) { delete childTNE_ptr; childTNE_ptr = NULL; }
  } /* endif */
}

inline void OctreeBlock::read(istream &in_file,
                              AdaptiveBlock3D_ResourceList &List_of_Available_Blocks) {
  int number_of_children; 
  in_file >> block;
  if (block.used) {
     if (List_of_Available_Blocks.Nfree > 0) {
        block.info.cpu = List_of_Available_Blocks.nextCPU();
        block.info.blknum = List_of_Available_Blocks.nextBlock();
        List_of_Available_Blocks.update_next();
     } else {
        cout << "\n " << CFFC_Version() 
             << " Octree Read Error: Insufficient number of hexahedral solution blocks.\n";
     } /* endif */
  } /* endif */
  in_file.setf(ios::skipws);  in_file >> number_of_children;  
  in_file.unsetf(ios::skipws);
  if (!block.used && number_of_children > 0) {
     if (childBSW_ptr != NULL) { delete childBSW_ptr; childBSW_ptr = NULL; }
     childBSW_ptr = new OctreeBlock;
     childBSW_ptr->parent_ptr = this;
     childBSW_ptr->read(in_file, List_of_Available_Blocks);
     if (childBSE_ptr != NULL) { delete childBSE_ptr; childBSE_ptr = NULL; }
     childBSE_ptr = new OctreeBlock;
     childBSE_ptr->parent_ptr = this;
     childBSE_ptr->read(in_file, List_of_Available_Blocks);
     if (childBNW_ptr != NULL) { delete childBNW_ptr; childBNW_ptr = NULL; }
     childBNW_ptr = new OctreeBlock;
     childBNW_ptr->parent_ptr = this;
     childBNW_ptr->read(in_file, List_of_Available_Blocks);
     if (childBNE_ptr != NULL) { delete childBNE_ptr; childBNE_ptr = NULL; }
     childBNE_ptr = new OctreeBlock;
     childBNE_ptr->parent_ptr = this;
     childBNE_ptr->read(in_file, List_of_Available_Blocks);
     if (childTSW_ptr != NULL) { delete childTSW_ptr; childTSW_ptr = NULL; }
     childTSW_ptr = new OctreeBlock;
     childTSW_ptr->parent_ptr = this;
     childTSW_ptr->read(in_file, List_of_Available_Blocks);
     if (childTSE_ptr != NULL) { delete childTSE_ptr; childTSE_ptr = NULL; }
     childTSE_ptr = new OctreeBlock;
     childTSE_ptr->parent_ptr = this;
     childTSE_ptr->read(in_file, List_of_Available_Blocks);
     if (childTNW_ptr != NULL) { delete childTNW_ptr; childTNW_ptr = NULL; }
     childTNW_ptr = new OctreeBlock;
     childTNW_ptr->parent_ptr = this;
     childTNW_ptr->read(in_file, List_of_Available_Blocks);
     if (childTNE_ptr != NULL) { delete childTNE_ptr; childTNE_ptr = NULL; }
     childTNE_ptr = new OctreeBlock;
     childTNE_ptr->parent_ptr = this;
     childTNE_ptr->read(in_file, List_of_Available_Blocks);
   } else {
     if (childBSW_ptr != NULL) { delete childBSW_ptr; childBSW_ptr = NULL; }
     if (childBSE_ptr != NULL) { delete childBSE_ptr; childBSE_ptr = NULL; }
     if (childBNW_ptr != NULL) { delete childBNW_ptr; childBNW_ptr = NULL; }
     if (childBNE_ptr != NULL) { delete childBNE_ptr; childBNE_ptr = NULL; }
     if (childTSW_ptr != NULL) { delete childTSW_ptr; childTSW_ptr = NULL; }
     if (childTSE_ptr != NULL) { delete childTSE_ptr; childTSE_ptr = NULL; }
     if (childTNW_ptr != NULL) { delete childTNW_ptr; childTNW_ptr = NULL; }
     if (childTNE_ptr != NULL) { delete childTNE_ptr; childTNE_ptr = NULL; }
  } /* endif */
}

/*************************************************************
 * OctreeBlock::write -- Write octree block (recursive).     *
 *************************************************************/
inline void OctreeBlock::write(ostream &out_file) {
  int number_of_children; 
  out_file << block << "\n";
  if (!block.used && childBSW_ptr != NULL) {
    number_of_children = 8;
  } else {
    number_of_children = 0;
  } /* endif */
  out_file << " " << number_of_children << "\n";
  if (childBSW_ptr != NULL) childBSW_ptr->write(out_file);
  if (childBSE_ptr != NULL) childBSE_ptr->write(out_file);
  if (childBNW_ptr != NULL) childBNW_ptr->write(out_file);
  if (childBNE_ptr != NULL) childBNE_ptr->write(out_file);
  if (childTSW_ptr != NULL) childTSW_ptr->write(out_file);
  if (childTSE_ptr != NULL) childTSE_ptr->write(out_file);
  if (childTNW_ptr != NULL) childTNW_ptr->write(out_file);
  if (childTNE_ptr != NULL) childTNE_ptr->write(out_file);
}

/*************************************************************
 * OctreeBlock::broadcast -- Broadcast octree block.         *
 *************************************************************/
inline void OctreeBlock::broadcast(void) {
#ifdef _MPI_VERSION
  int number_of_children; 
   AdaptiveBlock3D::Broadcast_Adaptive_Block(block);
  if (CFFC_Primary_MPI_Processor()) {
     if (!block.used && childBSW_ptr != NULL) {
       number_of_children = 8;
     } else {
       number_of_children = 0;
     } /* endif */
  } /* endif */
  MPI::COMM_WORLD.Bcast(&number_of_children, 1, MPI::INT, 0);
  if (!block.used && number_of_children > 0) {
     if (!CFFC_Primary_MPI_Processor() && childBSW_ptr != NULL) { delete childBSW_ptr; childBSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childBSW_ptr = new OctreeBlock; childBSW_ptr->parent_ptr = this; }
     childBSW_ptr->broadcast();
     if (!CFFC_Primary_MPI_Processor() && childBSE_ptr != NULL) { delete childBSE_ptr; childBSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childBSE_ptr = new OctreeBlock; childBSE_ptr->parent_ptr = this; }
     childBSE_ptr->broadcast();
     if (!CFFC_Primary_MPI_Processor() && childBNW_ptr != NULL) { delete childBNW_ptr; childBNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childBNW_ptr = new OctreeBlock; childBNW_ptr->parent_ptr = this; }
     childBNW_ptr->broadcast();
     if (!CFFC_Primary_MPI_Processor() && childBNE_ptr != NULL) { delete childBNE_ptr; childBNE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childBNE_ptr = new OctreeBlock; childBNE_ptr->parent_ptr = this; }
     childBNE_ptr->broadcast();
     if (!CFFC_Primary_MPI_Processor() && childTSW_ptr != NULL) { delete childTSW_ptr; childTSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childTSW_ptr = new OctreeBlock; childTSW_ptr->parent_ptr = this; }
     childTSW_ptr->broadcast();
     if (!CFFC_Primary_MPI_Processor() && childTSE_ptr != NULL) { delete childTSE_ptr; childTSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childTSE_ptr = new OctreeBlock; childTSE_ptr->parent_ptr = this; }
     childTSE_ptr->broadcast();
     if (!CFFC_Primary_MPI_Processor() && childTNW_ptr != NULL) { delete childTNW_ptr; childTNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childTNW_ptr = new OctreeBlock; childTNW_ptr->parent_ptr = this; }
     childTNW_ptr->broadcast();
     if (!CFFC_Primary_MPI_Processor() && childTNE_ptr != NULL) { delete childTNE_ptr; childTNE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childTNE_ptr = new OctreeBlock; childTNE_ptr->parent_ptr = this; }
     childTNE_ptr->broadcast();
  } else {
     if (!CFFC_Primary_MPI_Processor() && childBSW_ptr != NULL) { delete childBSW_ptr; childBSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childBSE_ptr != NULL) { delete childBSE_ptr; childBSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childBNW_ptr != NULL) { delete childBNW_ptr; childBNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childBNE_ptr != NULL) { delete childBNE_ptr; childBNE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childTSW_ptr != NULL) { delete childTSW_ptr; childTSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childTSE_ptr != NULL) { delete childTSE_ptr; childTSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childTNW_ptr != NULL) { delete childTNW_ptr; childTNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childTNE_ptr != NULL) { delete childTNE_ptr; childTNE_ptr = NULL; }
  } /* endif */
#endif
}

inline void OctreeBlock::broadcast(AdaptiveBlock3D_ResourceList &List_of_Available_Blocks) {
#ifdef _MPI_VERSION
  int number_of_children; 
  AdaptiveBlock3D::Broadcast_Adaptive_Block(block);
  if (!CFFC_Primary_MPI_Processor()) {
     if (block.used) {
        if (List_of_Available_Blocks.Nfree > 0) {
           block.info.cpu = List_of_Available_Blocks.nextCPU();
           block.info.blknum = List_of_Available_Blocks.nextBlock();
           List_of_Available_Blocks.update_next();
        } else {
           cout << "\n " << CFFC_Version() 
                << " Octree Broadcast Error: Insufficient number of ocrilateral solution blocks.\n";
        } /* endif */
     } /* endif */
  } /* endif */
  if (CFFC_Primary_MPI_Processor()) {
     if (!block.used && childBSW_ptr != NULL) {
       number_of_children = 8;
     } else {
       number_of_children = 0;
     } /* endif */
  } /* endif */
  MPI::COMM_WORLD.Bcast(&number_of_children, 1, MPI::INT, 0);
  if (!block.used && number_of_children > 0) {
     if (!CFFC_Primary_MPI_Processor() && childBSW_ptr != NULL) { delete childBSW_ptr; childBSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childBSW_ptr = new OctreeBlock; childBSW_ptr->parent_ptr = this; }
     childBSW_ptr->broadcast(List_of_Available_Blocks);
     if (!CFFC_Primary_MPI_Processor() && childBSE_ptr != NULL) { delete childBSE_ptr; childBSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childBSE_ptr = new OctreeBlock; childBSE_ptr->parent_ptr = this; }
     childBSE_ptr->broadcast(List_of_Available_Blocks);
     if (!CFFC_Primary_MPI_Processor() && childBNW_ptr != NULL) { delete childBNW_ptr; childBNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childBNW_ptr = new OctreeBlock; childBNW_ptr->parent_ptr = this; }
     childBNW_ptr->broadcast(List_of_Available_Blocks);
     if (!CFFC_Primary_MPI_Processor() && childBNE_ptr != NULL) { delete childBNE_ptr; childBNE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childBNE_ptr = new OctreeBlock; childBNE_ptr->parent_ptr = this; }
     childBNE_ptr->broadcast(List_of_Available_Blocks);
     if (!CFFC_Primary_MPI_Processor() && childTSW_ptr != NULL) { delete childTSW_ptr; childTSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childTSW_ptr = new OctreeBlock; childTSW_ptr->parent_ptr = this; }
     childTSW_ptr->broadcast(List_of_Available_Blocks);
     if (!CFFC_Primary_MPI_Processor() && childTSE_ptr != NULL) { delete childTSE_ptr; childTSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childTSE_ptr = new OctreeBlock; childTSE_ptr->parent_ptr = this; }
     childTSE_ptr->broadcast(List_of_Available_Blocks);
     if (!CFFC_Primary_MPI_Processor() && childTNW_ptr != NULL) { delete childTNW_ptr; childTNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childTNW_ptr = new OctreeBlock; childTNW_ptr->parent_ptr = this; }
     childTNW_ptr->broadcast(List_of_Available_Blocks);
     if (!CFFC_Primary_MPI_Processor() && childTNE_ptr != NULL) { delete childTNE_ptr; childTNE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childTNE_ptr = new OctreeBlock; childTNE_ptr->parent_ptr = this; }
     childTNE_ptr->broadcast(List_of_Available_Blocks);
  } else {
     if (!CFFC_Primary_MPI_Processor() && childTSW_ptr != NULL) { delete childTSW_ptr; childTSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childTSE_ptr != NULL) { delete childTSE_ptr; childTSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childTNW_ptr != NULL) { delete childTNW_ptr; childTNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childTNE_ptr != NULL) { delete childTNE_ptr; childTNE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childBSW_ptr != NULL) { delete childBSW_ptr; childBSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childBSE_ptr != NULL) { delete childBSE_ptr; childBSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childBNW_ptr != NULL) { delete childBNW_ptr; childBNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childBNE_ptr != NULL) { delete childBNE_ptr; childBNE_ptr = NULL; }
  } /* endif */
#endif
}

/*************************************************************
 * OctreeBlock::maxRefinementLevel -- Maximum refinement     *
 *                                    level (recursive).     *
 *************************************************************/
inline int OctreeBlock::maxRefinementLevel(int &Max_Level) {
  if (block.used) {
     Max_Level = max(Max_Level, block.info.level);
  } else if (!block.used && childBSW_ptr != NULL) {
     if (childBSW_ptr != NULL) childBSW_ptr->maxRefinementLevel(Max_Level);
     if (childBSE_ptr != NULL) childBSE_ptr->maxRefinementLevel(Max_Level);
     if (childBNW_ptr != NULL) childBNW_ptr->maxRefinementLevel(Max_Level);
     if (childBNE_ptr != NULL) childBNE_ptr->maxRefinementLevel(Max_Level);
     if (childTSW_ptr != NULL) childTSW_ptr->maxRefinementLevel(Max_Level);
     if (childTSE_ptr != NULL) childTSE_ptr->maxRefinementLevel(Max_Level);
     if (childTNW_ptr != NULL) childTNW_ptr->maxRefinementLevel(Max_Level);
     if (childTNE_ptr != NULL) childTNE_ptr->maxRefinementLevel(Max_Level);
  } /* endif */
  return (Max_Level);
}

/*************************************************************
 * OctreeBlock -- Input-output operators.                    *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const OctreeBlock &Block) {
    out_file << Block.block;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     OctreeBlock &Block) {
  in_file >> Block.block;
  return (in_file);
}

/********************************************************
 * Class: Octree_DataStructure                          *
 *                                                      *
 * Member functions                                     *
 *   Roots -- Return roots of octree data structure.    *
 *   NRi   -- Return number of root blocks in           *
 *            the i-direction (zeta-direction).         *
 *   NRj   -- Return number of root blocks in           *
 *            the j-direction (eta-direction).          *
 *   NRk   -- Return number of root blocks in           *
 *            the k-direction (???-direction).          *
 *   Blocks -- Return global list of pointers to the    *
 *             blocks in use in the octree data         *
 *             structure.                               *
 *   RefineFlags -- Return global list of mesh          *
 *                  refinement flags for blocks in data *
 *                  structure.                          *
 *   Ncpu  -- Return number of CPUs available for the   *
 *            calculation.                              *
 *   Nblk  -- Return the number of local blocks per     *
 *            CPU available for the calculation.        * 
 *   MaximumRefinementLevel -- Return maximum           *
 *            allowable level of refinement for         *
 *            adaptive blocks.                          *
 *   MinimumRefinementLevel -- Return minimum           *
 *            allowable level of refinement for         *
 *            adaptive blocks.                          *
 *   allocate -- Allocate memory for octree data        *
 *               structure.                             *
 *   deallocate -- Deallocate memory for octree data    *
 *                 structure.                           *
 *   allocateRoots -- Allocate memory for roots of      *
 *                    octree data structure.            *
 *   deallocateRoots -- Deallocate memory for roots of  *
 *                      octree data structure.          *
 *   allocateBlocks -- Allocate memory for global block *
 *                     list of octree data structure.   *
 *   deallocateBlocks -- Deallocate memory for global   *
 *                       block list of octree data      *
 *                       structure.                     *
 *   assign_block_pointers, assign_block_ptr            *
 *                    -- Assign block pointers to       *
 *                       all used blocks.               *
 *   renumber -- Renumbers octree adaptive blocks,      *
 *               assigning global block numbers.        *
 *   countUsedBlocks -- Returns the number of used      *
 *                      blocks in the octree data       *
 *                      structure.                      *
 *   countUsedCells -- Returns the number of            *
 *                     computational cells in the       *
 *                     octree data structure.           *
 *   efficiencyRefinement -- Returns the ratio of the   *
 *                           current number of cells    *
 *                           used to the number of cells*
 *                           that would have been used  *
 *                           with a uniform mesh.       *
 *   getRoot -- Return indices of root octree block.    *
 *   getNeighbour -- Return neighbour of octree block   *
 *                   in specified search direction of   *
 *                   interest.                          *
 *   findNeighbours -- Find neighbouring blocks of      *
 *                     octree block.                    *
 *   refineBlock   -- Refines (divides) a octree        *
 *                    block into four offspring.        *
 *   coarsenBlocks -- Coarsens (contracts) four ocree   *
 *                    blocks into a single parent.      *
 *   nochangeAll-- Sets the mesh refinement flags to    *
 *                 force no refinement or coarsening of *
 *                 octree blocks (default).             *
 *   refineAll  -- Sets the mesh refinement flags to    *
 *                 force refinement of all octree       *
 *                 blocks.                              *
 *   coarsenAll -- Sets the mesh refinement flags to    *
 *                 force coasening of all octree        *
 *                 blocks.                              *
 *   setRefineAll -- Sets the mesh refinement flags of  *
 *                   all octree blocks to specified     *
 *                   value.                             *
 *   highestRefinementLevel -- Return highest level     *
 *                of block refinement.                  *
 *   numberToBeRefined -- Return number of blocks to be *
 *                refined.                              *
 *   numberToBeCoarsened -- Return number of blocks to  *
 *                be coarsened.                         * 
 *                                                      *
 * Member operators                                     *
 *      QT -- octree hierarchical data structure        *
 *                                                      *
 * QT = QT;                                             *
 * cout << QT; (output function)                        *
 * cin  >> QT; (input function)                         *
 *                                                      *
 ********************************************************/
class Octree_DataStructure{
  private:
  public:
    int                     NR; // Total number of root blocks.
    int                    NRi; // Number of roots in i-direction.
    int                    NRj; // Number of roots in j-direction.
    int                    NRk; // Number of roots in k-direction.
    OctreeBlock         *Roots; // Roots of octree data structure.
    int                   Ncpu; // Number of CPUs available.
    int                   Nblk; // Number of local blocks per CPU.
    OctreeBlock      ***Blocks; // Global list of pointers to blocks
                                // in use in octree data structure.
    int          **RefineFlags; // Global list of refinement flags.
    double     RefineThreshold, // Thresholds for refinment and coarsening.
              CoarsenThreshold; //
    int MaximumRefinementLevel; // Maximum allowable refinement level.
    int MinimumRefinementLevel; // Maximum allowable refinement level.
                                // for blocks in data structures.

    /* Creation, copy, and assignment constructors. */
    Octree_DataStructure(void) {
       NR = 0; NRi = 0; NRj = 0; NRk = 0; Roots = NULL;
       Ncpu = 0; Nblk = 0; Blocks = NULL; RefineFlags = NULL;
       MaximumRefinementLevel = 99; MinimumRefinementLevel = 0;
       RefineThreshold = 0.50; CoarsenThreshold = 0.10;
    }

    Octree_DataStructure(const Octree_DataStructure &QT) {
       NR = QT.NR; NRi = QT.NRi; NRj = QT.NRj; NRk = QT.NRk;
       Roots = QT.Roots; 
       Ncpu = QT.Ncpu; Nblk = QT.Nblk; Blocks = QT.Blocks;
       RefineFlags = QT.RefineFlags;
       MaximumRefinementLevel = QT.MaximumRefinementLevel; 
       MinimumRefinementLevel = QT.MinimumRefinementLevel; 
       RefineThreshold = QT.RefineThreshold; 
       CoarsenThreshold = QT.CoarsenThreshold;
    }

    /* Destructor. */
    // ~Octree_DataStructure(void);
    // Use automatically generated destructor.
    ~Octree_DataStructure(void) { deallocate(); }

    /* Assignment operator. */
    // Octree_DataStructure operator = (const Octree_DataStructure &QT);
    // Use automatically generated assignment operator.

    /* Allocate memory for octree data structure. */
    void allocate(const int n, const int ncpu, const int nblk);
    void allocate(const int ni, const int nj, const int nk, 
                  const int ncpu, const int nblk);

    /* Deallocate memory for octree data structure. */
    void deallocate(void);

    /* Allocate memory for octree data structure roots. */
    void allocateRoots(const int n);
    void allocateRoots(const int ni, const int nj, const int nk);

    /* Deallocate memory for octree data structure roots. */
    void deallocateRoots(void);

    /* Allocate memory for octree data structure global block list. */
    void allocateBlocks(const int ncpu, const int nblk);

    /* Deallocate memory for octree data structure global block list. */
    void deallocateBlocks(void);

    /* Assign block pointers to all used blocks. */
    void assign_block_ptr(OctreeBlock *Block_Ptr);
    void assign_block_pointers(void);

    /* Renumber all octree blocks. */
    void renumber(void);

    /* Count number of used blocks and cells in octree data structure 
       and other useful mesh statistics. */
    int countUsedBlocks(void);
    int countUsedCells(void);
    double efficiencyRefinement(void);

    /* Octree block root indices. */
    AdaptiveBlock3D_Dimensions getRoot(OctreeBlock *Block_Ptr);

    /* Reassign neighbour information for root blocks. */
    void reassignRootNeighbours(void);

    /* Find neighbour of octree block in specified search direction. */
    OctreeBlock *getNeighbour(OctreeBlock *Block_Ptr,
                              const int Search_Dir_Mask);

    /* Find neighbouring blocks. */
    void findNeighbours(void);

    /* Block refinement and coarsening. */
    void refineBlock(int *new_blocks_CPU, 
                     int *new_blocks_BLK, 
                     int *new_blocks_SECTOR);
    void coarsenBlocks(int *old_blocks_CPU, 
                       int *old_blocks_BLK, 
                       int *old_blocks_SECTOR);

    /* Set mesh refinement and coarsening flags. */
    void nochangeAll(void);
    void refineAll(void);
    void coarsenAll(void);
    void setRefineAll(const int Flag);

    /* Mesh refinement information. */
    int highestRefinementLevel(void);
    int highestRefinementLevel(const int use_tree);
    int numberToBeRefined(void);
    int numberToBeCoarsened(void);

    static void Create_Octree_Data_Structure(Octree_DataStructure &Octree,
                                             const int Number_of_Roots_Idir,
                                             const int Number_of_Roots_Jdir,
                                             const int Number_of_Roots_Kdir,
                                             const int Number_of_Processors,
                                             const int Number_of_Blocks_per_Processor);

    static void Broadcast_Octree_Data_Structure(Octree_DataStructure &Octree,
                                                AdaptiveBlock3D_ResourceList  &List_of_Available_Blocks);
    
    static void Renumber_Solution_Blocks(Octree_DataStructure &Octree);

    static void Renumber_Solution_Blocks(Octree_DataStructure &Octree,
                                         AdaptiveBlock3D_List &LocalSolnBlockList);

    static void Reassign_Root_Neighbours(Octree_DataStructure &Octree);

    static void Reassign_Root_Neighbours(Octree_DataStructure &Octree,
                                         AdaptiveBlock3D_List &LocalSolnBlockList);

    static void Find_Neighbours(Octree_DataStructure &Octree);

    static void Find_Neighbours(Octree_DataStructure &Octree,
                                AdaptiveBlock3D_List &LocalSolnBlockList);

    static void Get_Refinement_List(Octree_DataStructure &Octree,
                                    AdaptiveBlock3D_List &LocalSolnBlockList);    

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Octree_DataStructure &QT);
    friend istream &operator >> (istream &in_file,
				 Octree_DataStructure &QT);
    
};

/*********************************************************************
 * Octree_DataStructure::allocate -- Allocate memory.                *
 *********************************************************************/
inline void Octree_DataStructure::allocate(const int n, 
                                           const int ncpu, 
                                           const int nblk) {
   assert( n > 0 && ncpu > 0 && nblk > 0 ); 
   NR = n; NRi = n; NRj = 1; NRk = 1; Ncpu = ncpu; Nblk = nblk;
   Roots = new OctreeBlock[NR];
   Blocks = new OctreeBlock**[Ncpu];
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      Blocks[i] = new OctreeBlock*[Nblk];
      for (int j = 0; j <= Nblk-1; ++j) Blocks[i][j] = NULL; 
   } /* endfor */
   RefineFlags = new int*[Ncpu];
   for (int i = 0; i <= Ncpu-1 ; ++i ) RefineFlags[i] = new int[Nblk];
   nochangeAll();
}

inline void Octree_DataStructure::allocate(const int ni, 
                                           const int nj, 
                                           const int nk, 
                                           const int ncpu, 
                                           const int nblk) {
   assert( ni > 0 && nj > 0 && nk > 0 && ncpu > 0 && nblk > 0 );
   NR = ni*nj*nk; NRi = ni; NRj = nj; NRk = nk; Ncpu = ncpu; Nblk = nblk;
   Roots = new OctreeBlock[NR];
   Blocks = new OctreeBlock**[Ncpu];
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      Blocks[i] = new OctreeBlock*[Nblk];
      for (int j = 0; j <= Nblk-1; ++j) Blocks[i][j] = NULL; 
   } /* endfor */
   RefineFlags = new int*[Ncpu];
   for (int i = 0; i <= Ncpu-1 ; ++i) RefineFlags[i] = new int[Nblk];
   nochangeAll();
}

inline void Octree_DataStructure::allocateRoots(const int n) {
   assert( n > 0 );
   NR = n; NRi = n; NRj = 1; NRk = 1;
   Roots = new OctreeBlock[NR];
}

inline void Octree_DataStructure::allocateRoots(const int ni, 
						const int nk,
                                                const int nj) {
   assert( ni > 0 && nj > 0 && nk > 0 );
   NR = ni*nj*nk; NRi = ni; NRj = nj; NRk = nk;
   Roots = new OctreeBlock[NR];
}

inline void Octree_DataStructure::allocateBlocks(const int ncpu, 
                                                 const int nblk) {
   assert( ncpu > 0 && nblk > 0 );
   Ncpu = ncpu; Nblk = nblk; Blocks = new OctreeBlock**[Ncpu];
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      Blocks[i] = new OctreeBlock*[Nblk];
      for (int j = 0; j <= Nblk-1; ++j) Blocks[i][j] = NULL; 
   } /* endfor */
   RefineFlags = new int*[Ncpu];
   for (int i = 0; i <= Ncpu-1 ; ++i) RefineFlags[i] = new int[Nblk];
   nochangeAll();
}

/*********************************************************************
 * Octree_DataStructure::deallocate -- Deallocate memory.            *
 *********************************************************************/
inline void Octree_DataStructure::deallocate(void) {
   delete []Roots; Roots = NULL;
   NR = 0; NRi = 0; NRj = 0; NRk = 0;
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      delete []Blocks[i]; Blocks[i] = NULL;
   } /* endfor */
   delete []Blocks; Blocks = NULL;
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      delete []RefineFlags[i]; RefineFlags[i] = NULL;
   } /* endfor */
   delete []RefineFlags; RefineFlags = NULL;
   Ncpu = 0; Nblk = 0;
}

inline void Octree_DataStructure::deallocateRoots(void) {
  delete []Roots; Roots = NULL;
  NR = 0; NRi = 0; NRj = 0; NRk = 0;
}

inline void Octree_DataStructure::deallocateBlocks(void) {
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      delete []Blocks[i]; Blocks[i] = NULL;
   } /* endfor */
   delete []Blocks; Blocks = NULL;
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      delete []RefineFlags[i]; RefineFlags[i] = NULL;
   } /* endfor */
   delete []RefineFlags; RefineFlags = NULL;
   Ncpu = 0; Nblk = 0;
}

/****************************************************************************
 * Octree_DataStructure::assign_block_ptr -- Assign block ptr.              *
 ****************************************************************************/
inline void Octree_DataStructure::assign_block_ptr(OctreeBlock *Block_Ptr) {
   int icpu, iblk;
   if (Block_Ptr != NULL) {
      if (Block_Ptr->block.used) {
	 icpu = Block_Ptr->block.info.cpu;
         iblk = Block_Ptr->block.info.blknum;
         if (icpu < 0 || icpu >= Ncpu ||
             iblk < 0 || iblk >= Nblk) {
            cout << "\n " << CFFC_Version() 
                 << " Octree Block Error: Invalid range for block indices.\n";
         } else {
            Blocks[icpu][iblk] = Block_Ptr;
         } /* endif */
      } else if (Block_Ptr->childBSW_ptr != NULL) {
	 assign_block_ptr(Block_Ptr->childBSW_ptr);
	 assign_block_ptr(Block_Ptr->childBSE_ptr);
	 assign_block_ptr(Block_Ptr->childBNW_ptr);
	 assign_block_ptr(Block_Ptr->childBNE_ptr);
	 assign_block_ptr(Block_Ptr->childTSW_ptr);
	 assign_block_ptr(Block_Ptr->childTSE_ptr);
	 assign_block_ptr(Block_Ptr->childTNW_ptr);
	 assign_block_ptr(Block_Ptr->childTNE_ptr);
       } else {
	 // Not sure what to do!  This shouldn't happen.
      } /* endif */
   } /* endif */
}

inline void Octree_DataStructure::assign_block_pointers(void) {
   for (int jBLK = 0; jBLK <= Nblk-1; ++jBLK) {
      for (int iBLK = 0; iBLK <= Ncpu-1; ++iBLK) {
         Blocks[iBLK][jBLK] = NULL;
      } /* endfor */
   } /* endfor */
   for (int iBLK = 0 ; iBLK <= NR-1; ++iBLK ) {
      assign_block_ptr(&(Roots[iBLK]));
   } /* endfor */
}

/****************************************************************************
 * Octree_DataStructure::renumber -- Renumber octree blocks.                *
 ****************************************************************************/
inline void Octree_DataStructure::renumber(void) {
   int global_block_number(0);
   for (int nBLK = 0 ; nBLK <= Nblk-1 ; ++nBLK) {
      for (int iCPU = 0 ; iCPU <= Ncpu-1 ; ++iCPU) {
         if (Blocks[iCPU][nBLK] != NULL) {
            if (Blocks[iCPU][nBLK]->block.used) {
               Blocks[iCPU][nBLK]->block.info.gblknum = global_block_number;
               ++global_block_number;
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */
}

/****************************************************************************
 * Octree_DataStructure::countUsedBlocks -- Number of used blocks.          *
 ****************************************************************************/
inline int Octree_DataStructure::countUsedBlocks(void) {
  int number_of_used_blocks(0);
  for (int nBLK = 0 ; nBLK <= Nblk-1 ; ++nBLK ) {
      for (int iCPU = 0 ; iCPU <= Ncpu-1 ; ++iCPU ) {
         if (Blocks[iCPU][nBLK] != NULL) {
            if (Blocks[iCPU][nBLK]->block.used) {
               ++number_of_used_blocks;
            } /* endif */
         } /* endif */
      } /* endfor */
  } /* endfor */
  return (number_of_used_blocks);
}

/****************************************************************************
 * Octree_DataStructure::countUsedCells -- Number of used cells.            *
 ****************************************************************************/
inline int Octree_DataStructure::countUsedCells(void) {
  int number_of_used_cells(0);
  for (int nBLK = 0 ; nBLK <= Nblk-1 ; ++nBLK) {
      for (int iCPU = 0 ; iCPU <= Ncpu-1 ; ++iCPU) {
         if (Blocks[iCPU][nBLK] != NULL) {
            if (Blocks[iCPU][nBLK]->block.used) {
               number_of_used_cells += Blocks[iCPU][nBLK]->block.info.dimen.i*
                                       Blocks[iCPU][nBLK]->block.info.dimen.j*
                                       Blocks[iCPU][nBLK]->block.info.dimen.k;
            } /* endif */
         } /* endif */
      } /* endfor */
  } /* endfor */
  return (number_of_used_cells);
}

/******************************************************************************
 * Octree_DataStructure::efficiencyRefinement -- Refinment efficiency.        *
 ******************************************************************************/
inline double Octree_DataStructure::efficiencyRefinement(void) {
  int number_of_used_cells(0), number_of_cells_with_uniform_mesh(0), max_level;
  max_level = highestRefinementLevel();
  for (int nBLK = 0 ; nBLK <= Nblk-1 ; ++nBLK) {
      for (int iCPU = 0 ; iCPU <= Ncpu-1 ; ++iCPU) {
         if (Blocks[iCPU][nBLK] != NULL) {
            if (Blocks[iCPU][nBLK]->block.used) {
               number_of_used_cells += Blocks[iCPU][nBLK]->block.info.dimen.i*
                                       Blocks[iCPU][nBLK]->block.info.dimen.j*
                                       Blocks[iCPU][nBLK]->block.info.dimen.k;
               number_of_cells_with_uniform_mesh +=
                                       int(pow(double(8), double(max_level-Blocks[iCPU][nBLK]->block.info.level)))*
                                       (Blocks[iCPU][nBLK]->block.info.dimen.i*
                                        Blocks[iCPU][nBLK]->block.info.dimen.j*
                                        Blocks[iCPU][nBLK]->block.info.dimen.k);
            } /* endif */
         } /* endif */
      } /* endfor */
  } /* endfor */
  return (ONE-double(number_of_used_cells)/double(number_of_cells_with_uniform_mesh));
}

/****************************************************************************
 * Octree_DataStructure::getRoot -- Get root indices.                       *
 ****************************************************************************/
inline AdaptiveBlock3D_Dimensions Octree_DataStructure::getRoot(OctreeBlock *Block_Ptr) {
  for (int kBLK = 0 ; kBLK <= NRk-1 ; ++kBLK) {
    for (int jBLK = 0 ; jBLK <= NRj-1 ; ++jBLK) {
       for (int iBLK = 0 ; iBLK <= NRi-1 ; ++iBLK) {
          if (Block_Ptr == &(Roots[iBLK*jBLK*kBLK])) {
             return (AdaptiveBlock3D_Dimensions(iBLK, jBLK, kBLK, 0));
          } /* endif */
       } /* endfor */
    } /* endfor */
  } /* endfor */
  return (AdaptiveBlock3D_Dimensions(-1, -1, -1, 0));
}

/****************************************************************************
 * Octree_DataStructure::reassignRootNeighbours -- Reassign root neighbours.*
 ****************************************************************************/
inline void Octree_DataStructure::reassignRootNeighbours(void) {

   int nr_neighbour;

   for (int nr = 0; nr <= NR-1; ++nr) {
      if (Roots[nr].block.nT > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoT[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoT[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoT[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nB > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoB[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoB[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoB[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */
   
      if (Roots[nr].block.nN > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoN[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoN[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoN[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */
   
      if (Roots[nr].block.nS > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoS[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoS[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoS[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */
 
      if (Roots[nr].block.nE > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoE[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoE[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoE[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nW > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoW[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoW[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoW[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */
      
      if (Roots[nr].block.nNW > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoNW[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoNW[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoNW[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nNE > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoNE[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoNE[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoNE[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */
   
      if (Roots[nr].block.nSE > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoSE[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoSE[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoSE[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nSW > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoSW[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoSW[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoSW[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nTN > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoTN[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoTN[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoTN[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nTS > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoTS[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoTS[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoTS[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */
   
      if (Roots[nr].block.nTE > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoTE[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoTE[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoTE[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nTW > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoTW[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoTW[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoTW[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nBN > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoBN[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoBN[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoBN[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nBS > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoBS[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoBS[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoBS[0].cpu = Roots[nr_neighbour].block.info.cpu;;
      } /* endif */

      if (Roots[nr].block.nBE > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoBE[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoBE[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoBE[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */
 
      if (Roots[nr].block.nBW > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoBW[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoBW[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoBW[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nTNW > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoTNW[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoTNW[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoTNW[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nTSW > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoTSW[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoTSW[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoTSW[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nTNE > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoTNE[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoTNE[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoTNE[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nTSE > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoTSE[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoTSE[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoTSE[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nBNW > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoBNW[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoBNW[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoBNW[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nBSW > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoBSW[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoBSW[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoBSW[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nBNE > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoBNE[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoBNE[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoBNE[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */

      if (Roots[nr].block.nBSE > 0) {
         for (int nr_search = 0; nr_search <= NR-1; ++nr_search) {
            if (Roots[nr_search].block.info.gblknum == Roots[nr].block.infoBSE[0].gblknum) {
	       nr_neighbour = nr_search;
	       break;
	    } /* endif */
         } /* endfor */
         Roots[nr].block.infoBSE[0].blknum = Roots[nr_neighbour].block.info.blknum;
         Roots[nr].block.infoBSE[0].cpu = Roots[nr_neighbour].block.info.cpu;
      } /* endif */
   } /* endfor */

}

/****************************************************************************
 * Octree_DataStructure::getNeighbour -- Return neighbour.                  *
 ****************************************************************************/
inline OctreeBlock *Octree_DataStructure::getNeighbour(OctreeBlock *Block_Ptr,
                                                       const int Search_Dir_Mask) {
  assert(Block_Ptr!=NULL);//used for debugging
  OctreeBlock *neighbour_block_ptr;
  return (neighbour_block_ptr);
}


/****************************************************************************
 * Octree_DataStructure::findNeighbours -- Find block neighbours.           *
 ****************************************************************************/
inline void Octree_DataStructure::findNeighbours(void) {

  cout << "\nOctree_DataStructure::findNeighbours() NOT YET IMPLEMENTED\n";

}
 
/****************************************************************************
 * Octree_DataStructure::refineBlock -- Refine (divide) block.              *
 ****************************************************************************/
inline void Octree_DataStructure::refineBlock(int *new_blocks_CPU, 
                                              int *new_blocks_BLK,
                                              int *new_blocks_SECTOR) {

  cout << "\nOctree_DataStructure::refineBlock() NOT YET IMPLEMENTED\n";

}

/****************************************************************************
 * Octree_DataStructure::coarsenBlock -- Coarsen (contract) blocks.         *
 ****************************************************************************/
inline void Octree_DataStructure::coarsenBlocks(int *old_blocks_CPU, 
                                                int *old_blocks_BLK,
                                                int *old_blocks_SECTOR) {

  cout << "\nOctree_DataStructure::coarsenBlock() NOT YET IMPLEMENTED\n";

}

/****************************************************************************
 * Octree_DataStructure::nochangeAll -- Set no refinement flags.            *
 ****************************************************************************/
inline void Octree_DataStructure::nochangeAll(void) {
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      for (int j = 0; j <= Nblk-1 ; ++j) {
         RefineFlags[i][j] = ADAPTIVEBLOCK3D_NOCHANGE;
      } /* endfor */
   } /* endfor */
}

/****************************************************************************
 * Octree_DataStructure::refineAll -- Set refinement flags.                 *
 ****************************************************************************/
inline void Octree_DataStructure::refineAll(void) {
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      for (int j = 0; j <= Nblk-1 ; ++j) {
         if (Blocks[i][j] != NULL) {
            if (Blocks[i][j]->block.used) {
               RefineFlags[i][j] = ADAPTIVEBLOCK3D_REFINE;
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */
}

/****************************************************************************
 * Octree_DataStructure::coarsenAll -- Set coarsening flags.                *
 ****************************************************************************/
inline void Octree_DataStructure::coarsenAll(void) {
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      for (int j = 0; j <= Nblk-1 ; ++j) {
         if (Blocks[i][j] != NULL) {
            if (Blocks[i][j]->block.used) {
               RefineFlags[i][j] = ADAPTIVEBLOCK3D_COARSEN;
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */
}

/****************************************************************************
 * Octree_DataStructure::setRefineAll -- Set refinement flags.              *
 ****************************************************************************/
inline void Octree_DataStructure::setRefineAll(const int Flag) {
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      for (int j = 0; j <= Nblk-1 ; ++j) {
         if (Blocks[i][j] != NULL) {
            if (Blocks[i][j]->block.used) {
               RefineFlags[i][j] = Flag;
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */
}

/****************************************************************************
 * Octree_DataStructure::highestRefinementLevel -- Return highest           *
 *                                 refinement level of all solution blocks. *
 ****************************************************************************/
inline int Octree_DataStructure::highestRefinementLevel(void) {
  int level(0);
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      for (int j = 0; j <= Nblk-1 ; ++j) {
         if (Blocks[i][j] != NULL) {
            if (Blocks[i][j]->block.used) {
               level = max(level, Blocks[i][j]->block.info.level);
            } /* endif */
         } /* endif */             
      } /* endfor */
   } /* endfor */
   return (level);
}

inline int Octree_DataStructure::highestRefinementLevel(const int use_tree) {
   int level(0);
   if (use_tree) {
      for (int n = 0 ; n <= NR-1 ; ++n) {
         level = Roots[n].maxRefinementLevel(level);
      } /* endfor */
   } else {
      for (int i = 0; i <= Ncpu-1 ; ++i) {
         for (int j = 0; j <= Nblk-1 ; ++j) {
            if (Blocks[i][j] != NULL) {
	       if (Blocks[i][j]->block.used) {
                  level = max(level, Blocks[i][j]->block.info.level);
               } /* endif */
            } /* endif */             
         } /* endfor */
      } /* endfor */
   } /* endif */
   return (level);
}

/****************************************************************************
 * Octree_DataStructure::numberToBeRefined -- Return number of              *
 *                                    octree solution blocks to be refined. *
 ****************************************************************************/
inline int Octree_DataStructure::numberToBeRefined(void) {
   int number_to_be_refined(0);
   for (int i = 0; i <= Ncpu-1 ; ++i) {
      for (int j = 0; j <= Nblk-1 ; ++j) {
         if (Blocks[i][j] != NULL) {
   	    if (Blocks[i][j]->block.used) {
               if (RefineFlags[i][j] == ADAPTIVEBLOCK3D_REFINE) {
                  number_to_be_refined += 1;
               } /* endif */
            } /* endif */
         } /* endif */             
      } /* endfor */
   } /* endfor */
   return (number_to_be_refined);
}

/****************************************************************************
 * Octree_DataStructure::numberToBeCoarsened -- Return number of            *
 *                                  octree solution blocks to be coarsened. *
 ****************************************************************************/
inline int Octree_DataStructure::numberToBeCoarsened(void) {
   int number_to_be_coarsened(0);
   for (int i = 0; i <= Ncpu-1 ; ++i ) {
      for (int j = 0; j <= Nblk-1 ; ++j ) {
         if (Blocks[i][j] != NULL) {
   	    if (Blocks[i][j]->block.used) {
               if (RefineFlags[i][j] == ADAPTIVEBLOCK3D_COARSEN) {
                  number_to_be_coarsened += 1;
               } /* endif */
            } /* endif */
         } /* endif */             
      } /* endfor */
   } /* endfor */
   return (number_to_be_coarsened);
}

/****************************************************************************
 * Octree_DataStructure -- Input-output operators.                          *
 ****************************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Octree_DataStructure &QT) {
   out_file << QT.NR << " " 
            << QT.NRi << " " 
            << QT.NRj << " " 
            << QT.NRk << " " 
            << QT.Ncpu << " " 
            << QT.Nblk << "\n";
   for (int n = 0 ; n <= QT.NR-1 ; ++n) {
      QT.Roots[n].write(out_file);
   } /* endfor */
   return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Octree_DataStructure &QT) {
   int nr, nri, nrj, nrk, ncpu, nblk;
   if (QT.Roots != NULL && QT.Blocks != NULL) {
     QT.deallocate();
   } else if (QT.Roots != NULL) {
     QT.deallocateRoots();
   } else if (QT.Blocks != NULL) {
      QT.deallocateBlocks();
   } /* endif */
   in_file.setf(ios::skipws);
   in_file >> nr >> nri >> nrj >> nrk >> ncpu >> nblk;
   in_file.setf(ios::skipws);
   QT.allocate(nri, nrj, nrk, ncpu, nblk);
   for (int n = 0 ; n <= QT.NR-1 ; ++n) {
      QT.Roots[n].read(in_file);
   } /* endfor */
   return (in_file);
}

#endif /* _OCTREE_INCLUDED  */

