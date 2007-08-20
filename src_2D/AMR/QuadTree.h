/* QuadTree.h:  Header file defining quadtree adaptive blocks 
                hierarchical data structure. */

#ifndef _QUADTREE_INCLUDED
#define _QUADTREE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include adaptive block header file. */

#ifndef _ADAPTIVEBLOCK_INCLUDED
#include "AdaptiveBlock.h"
#endif // _ADAPTIVEBLOCK_INCLUDED

/* Define the quadtree sector indicators for finding nearest
   neighbours and performing block refinement and coarsening. */

#define QUADTREE_SECTOR_NONE                          -1000
#define	QUADTREE_SECTOR_SW                             0
#define	QUADTREE_SECTOR_SE                             1
#define	QUADTREE_SECTOR_NW                             2
#define	QUADTREE_SECTOR_NE                             3

/* Define the quadtree search direction masks for finding
   nearest neighbours. */

#define	QUADTREE_DIRECTION_MASK_NORTH                  1
#define	QUADTREE_DIRECTION_MASK_SOUTH                  2
#define	QUADTREE_DIRECTION_MASK_EAST                   3
#define	QUADTREE_DIRECTION_MASK_WEST                   4
#define	QUADTREE_DIRECTION_MASK_NORTHWEST              5
#define	QUADTREE_DIRECTION_MASK_NORTHEAST              6
#define	QUADTREE_DIRECTION_MASK_SOUTHEAST              7
#define	QUADTREE_DIRECTION_MASK_SOUTHWEST              8

/* Define the classes. */

/********************************************************
 * Class: QuadTreeBlock                                 *
 *                                                      *
 * Member functions                                     *
 *  block       -- Return quadtree adaptive block.      *
 *  parent_ptr  -- Return pointer to parent of quadtree *
 *                 adaptive block.                      *
 *  childNW_ptr -- Return pointer to north-west child   *
 *                 of quadtree adaptive block.          *
 *  childNE_ptr -- Return pointer to north-east child   *
 *                 of quadtree adaptive block.          *
 *  childSW_ptr -- Return pointer to south-west child   *
 *                 of quadtree adaptive block.          *
 *  childSE_ptr -- Return pointer to south-east child   *
 *                 of quadtree adaptive block.          *
 *  child_ptr   -- Return pointer to specified child of *
 *                 quadtree adaptive block.             *
 *  search_dir  -- Returns the search direction given   *
 *                 the specified search direction mask. *
 *  sibling     -- Returns flag indicating whether or   *
 *                 not there is a sibling in the        *
 *                 specified search direction.          *
 *  neighbour_ptr -- Returns pointer to neighbouring    *
 *                 quadtree adaptive block found in the *
 *                 specified search direction.          *
 *  read        -- Reads in the quadtree block and then *
 *                 recursively descends the subtree and *
 *                 reads siblings and their siblings,   *
 *                 etc...                               *
 *  write       -- Writes out quadtree block and then   *
 *                 recursively descends the subtree and *
 *                 writes out siblings and their        *
 *                 siblings, etc...                     *
 *  broadcast   -- Broadcasts the quadtree block in a   *
 *                 recursive manner.                    *
 *  maxRefinementLevel -- Return the maximum refinement *
 *                        level of all blocks in        *
 *                        subtree.                      *
 *                                                      *
 * Member operators                                     *
 *      B -- element in a quadtree hierarchical data    *
 *           structure                                  *
 *                                                      *
 * B = B;                                               *
 * cout << B; (output function)                         *
 * cin  >> B; (input function)                          *
 *                                                      *
 ********************************************************/
class QuadTreeBlock{
  private:
  public:
    AdaptiveBlock2D        block;  // Pointer to quadtree adaptive block.
    QuadTreeBlock    *parent_ptr;  // Pointer to parent of adaptive block.
    QuadTreeBlock   *childNW_ptr,  // Pointers to children of
                    *childNE_ptr,  // adaptive block.
                    *childSE_ptr,  //
                    *childSW_ptr;  //
	                           // Made public so can access them.

    /* Creation, copy, and assignment constructors. */
    QuadTreeBlock(void) {
       parent_ptr = NULL; childNW_ptr = NULL; childNE_ptr = NULL; 
       childSE_ptr = NULL; childSW_ptr = NULL;
    }

    QuadTreeBlock(const QuadTreeBlock &Block) {
       block = Block.block; parent_ptr = Block.parent_ptr; 
       childNW_ptr = Block.childNW_ptr; childNE_ptr = Block.childNE_ptr; 
       childSE_ptr = Block.childSE_ptr; childSW_ptr = Block.childSW_ptr;
    }

    /* Destructor. */
    // ~QuadTreeBlock(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // QuadTreeBlock operator = (const QuadTreeBlock &Block);
    // Use automatically generated assignment operator.

    /* Pointer to child. */
    QuadTreeBlock *child_ptr(const int Sector);

    /* Search direction. */
    AdaptiveBlock2D_Dimensions search_dir(const int Search_Dir_Mask);

    /* Sibling in specified search direction? */
    int sibling(const int Search_Dir_Mask);

    /* Pointer to neighbour. */
    QuadTreeBlock *neighbour_ptr(const int Search_Dir_Mask);

    /* Read quadtree block (recursive). */
    void read(istream &in_file);
    void read(istream &in_file,
              AdaptiveBlockResourceList &List_of_Available_Blocks);

    /* Write quadtree block (recursive). */
    void write(ostream &out_file);

    /* Broadcast quadtree block (recursive). */
    void broadcast(void);
    void broadcast(AdaptiveBlockResourceList &List_of_Available_Blocks);

    /* Determine maximum refinement level. */
    int maxRefinementLevel(int &Max_Level);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const QuadTreeBlock &Block);
    friend istream &operator >> (istream &in_file,
				 QuadTreeBlock &Block);
    
};

/*************************************************************
 * QuadTreeBlock::child_ptr -- Return children of block.     *
 *************************************************************/
inline QuadTreeBlock *QuadTreeBlock::child_ptr(const int Sector) {
  switch(Sector) {
    case QUADTREE_SECTOR_NONE :
      return (NULL);
      break;
    case QUADTREE_SECTOR_NW :
      return (childNW_ptr);
      break;
    case QUADTREE_SECTOR_NE :
      return (childNE_ptr);
      break;
    case QUADTREE_SECTOR_SE :
      return (childSE_ptr);
      break;
    case QUADTREE_SECTOR_SW :
    default:
      return (childSW_ptr);
      break;
  } /* endswitch */
}

/*************************************************************
 * QuadTreeBlock::search_dir -- Determine search direction.  *
 *************************************************************/
inline AdaptiveBlock2D_Dimensions QuadTreeBlock::search_dir(const int Search_Dir_Mask) {
  int i, j;
  switch(Search_Dir_Mask) {
    case QUADTREE_DIRECTION_MASK_EAST :
      i = 1; j = 0;
      break;
    case QUADTREE_DIRECTION_MASK_WEST :
      i = -1; j = 0;
      break;
    case QUADTREE_DIRECTION_MASK_NORTH :
      i = 0; j = 2;
      break;
    case QUADTREE_DIRECTION_MASK_SOUTH :
      i = 0; j = -2;
      break;
    case QUADTREE_DIRECTION_MASK_SOUTHWEST :
      i = -1; j = -2;
      break;
    case QUADTREE_DIRECTION_MASK_SOUTHEAST :
      i = 1; j = -2;
      break;
    case QUADTREE_DIRECTION_MASK_NORTHWEST :
      i = -1; j = 2;
      break;
    case QUADTREE_DIRECTION_MASK_NORTHEAST :
      i = 1; j = 2;
      break;
    default:
      i = 0; j = 0;
      break;
  } /* endswitch */
  return (AdaptiveBlock2D_Dimensions(i, j, 0));
}

/*************************************************************
 * QuadTreeBlock::sibling -- Sibling in direction?           *
 *************************************************************/
inline int QuadTreeBlock::sibling(const int Search_Dir_Mask) {
  int i_sibling, j_sibling;
  AdaptiveBlock2D_Dimensions search_direction;
  if (block.info.sector != ADAPTIVEBLOCK2D_SECTOR_NONE) {
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
  } else {
     i_sibling = 0;
     j_sibling = 0;
  } /* endif */
  return (i_sibling && j_sibling);
}

/*************************************************************
 * QuadTreeBlock::neighbour_ptr -- Return neighbouring block.*
 *************************************************************/
inline QuadTreeBlock *QuadTreeBlock::neighbour_ptr(const int Search_Dir_Mask) {
  QuadTreeBlock *block_ptr; AdaptiveBlock2D_Dimensions search_direction;
  if (parent_ptr == NULL || block.info.sector == QUADTREE_SECTOR_NONE) {
     block_ptr = NULL; // Block has no parent.  Block is a root.  Return null pointer.
  } else {
     if (sibling(Search_Dir_Mask)) { // Neighbour is a sibling.  Get it.
        search_direction = search_dir(Search_Dir_Mask);
        block_ptr = parent_ptr->child_ptr(block.info.sector+
                                          search_direction.i+
                                          search_direction.j);
        if ((block_ptr->child_ptr(QUADTREE_SECTOR_SW-
                                  search_direction.i-
                                  search_direction.j+
                                  (max(-search_direction.i*search_direction.j, 0)/2)*
                                  max(search_direction.i, search_direction.j)) == NULL) ||
            (block_ptr->block.used)) { // Sibling is used and has no children.
	   block_ptr = block_ptr; // Return sibling.
	} else { // Return child of sibling.
           block_ptr = block_ptr->child_ptr(QUADTREE_SECTOR_SW-
                                            search_direction.i-
                                            search_direction.j+
                                            (max(-search_direction.i*search_direction.j, 0)/2)*
                                            max(search_direction.i, search_direction.j));
        } /* endif */
     } else if (Search_Dir_Mask <= QUADTREE_DIRECTION_MASK_WEST) { // Neighbour is not a sibling.
        block_ptr = parent_ptr->neighbour_ptr(Search_Dir_Mask); // Move up the tree.
        if (block_ptr == NULL) { // Parent of neighbour does not exist.
	   block_ptr = NULL; // Return null pointer.
        } else {
           search_direction = search_dir(Search_Dir_Mask);
           if ((block_ptr->child_ptr(block.info.sector-
                                     search_direction.i-
                                     search_direction.j) == NULL) ||
               (block_ptr->block.used)) { // Neighbour is used and has no children.
	      block_ptr = block_ptr; // Return neighbour.
	   } else { // Return child of neighbour.
              block_ptr = block_ptr->child_ptr(block.info.sector-
                                               search_direction.i-
                                               search_direction.j);
           } /* endif */
        } /* endif */
     } else { // Neighbour is not a sibling and is a corner neighbour.  Perform two-direction search.
        switch(Search_Dir_Mask) {
          case QUADTREE_DIRECTION_MASK_SOUTHWEST :
            if (sibling(QUADTREE_DIRECTION_MASK_WEST)) {
               block_ptr = neighbour_ptr(QUADTREE_DIRECTION_MASK_WEST);
               if (block_ptr != NULL) { 
                  block_ptr = block_ptr->neighbour_ptr(QUADTREE_DIRECTION_MASK_SOUTH);
	       } /* endif */
            } else {
               block_ptr = neighbour_ptr(QUADTREE_DIRECTION_MASK_SOUTH);
               if (block_ptr != NULL) { 
                  block_ptr = block_ptr->neighbour_ptr(QUADTREE_DIRECTION_MASK_WEST);
	       } /* endif */
            } /* endif */
            break;
          case QUADTREE_DIRECTION_MASK_SOUTHEAST :
            if (sibling(QUADTREE_DIRECTION_MASK_EAST)) {
               block_ptr = neighbour_ptr(QUADTREE_DIRECTION_MASK_EAST);
               if (block_ptr != NULL) { 
                  block_ptr = block_ptr->neighbour_ptr(QUADTREE_DIRECTION_MASK_SOUTH);
               } /* endif */
            } else {
               block_ptr = neighbour_ptr(QUADTREE_DIRECTION_MASK_SOUTH);
               if (block_ptr != NULL) { 
                  block_ptr = block_ptr->neighbour_ptr(QUADTREE_DIRECTION_MASK_EAST);
               } /* endif */
            } /* endif */
            break;
          case QUADTREE_DIRECTION_MASK_NORTHWEST :
            if (sibling(QUADTREE_DIRECTION_MASK_WEST)) {
               block_ptr = neighbour_ptr(QUADTREE_DIRECTION_MASK_WEST);
               if (block_ptr != NULL) { 
                  block_ptr = block_ptr->neighbour_ptr(QUADTREE_DIRECTION_MASK_NORTH);
               } /* endif */
            } else {
               block_ptr = neighbour_ptr(QUADTREE_DIRECTION_MASK_NORTH);
               if (block_ptr != NULL) { 
                  block_ptr = block_ptr->neighbour_ptr(QUADTREE_DIRECTION_MASK_WEST);
               } /* endif */
            } /* endif */
            break;
          case QUADTREE_DIRECTION_MASK_NORTHEAST :
            if (sibling(QUADTREE_DIRECTION_MASK_EAST)) {
               block_ptr = neighbour_ptr(QUADTREE_DIRECTION_MASK_EAST);
               if (block_ptr != NULL) { 
                  block_ptr = block_ptr->neighbour_ptr(QUADTREE_DIRECTION_MASK_NORTH);
               } /* endif */
            } else {
               block_ptr = neighbour_ptr(QUADTREE_DIRECTION_MASK_NORTH);
               if (block_ptr != NULL) { 
                  block_ptr = block_ptr->neighbour_ptr(QUADTREE_DIRECTION_MASK_EAST);
               } /* endif */
            } /* endif */
            break;
          default:
            block_ptr = NULL;
            break;
        } /* endswitch */
     } /* endif */
  } /* endif */
  return (block_ptr);
}

/*************************************************************
 * QuadTreeBlock::read -- Read quadtree block (recursive).   *
 *************************************************************/
inline void QuadTreeBlock::read(istream &in_file) {
  int number_of_children; 
  in_file >> block;
  in_file.setf(ios::skipws);  in_file >> number_of_children;  
  in_file.unsetf(ios::skipws);
  if (!block.used && number_of_children > 0) {
     if (childSW_ptr != NULL) { delete childSW_ptr; childSW_ptr = NULL; }
     childSW_ptr = new QuadTreeBlock;
     childSW_ptr->parent_ptr = this;
     childSW_ptr->read(in_file);
     if (childSE_ptr != NULL) { delete childSE_ptr; childSE_ptr = NULL; }
     childSE_ptr = new QuadTreeBlock;
     childSE_ptr->parent_ptr = this;
     childSE_ptr->read(in_file);
     if (childNW_ptr != NULL) { delete childNW_ptr; childNW_ptr = NULL; }
     childNW_ptr = new QuadTreeBlock;
     childNW_ptr->parent_ptr = this;
     childNW_ptr->read(in_file);
     if (childNE_ptr != NULL) { delete childNE_ptr; childNE_ptr = NULL; }
     childNE_ptr = new QuadTreeBlock;
     childNE_ptr->parent_ptr = this;
     childNE_ptr->read(in_file);
  } else {
     if (childSW_ptr != NULL) { delete childSW_ptr; childSW_ptr = NULL; }
     if (childSE_ptr != NULL) { delete childSE_ptr; childSE_ptr = NULL; }
     if (childNW_ptr != NULL) { delete childNW_ptr; childNW_ptr = NULL; }
     if (childNE_ptr != NULL) { delete childNE_ptr; childNE_ptr = NULL; }
  } /* endif */
}

inline void QuadTreeBlock::read(istream &in_file,
                                AdaptiveBlockResourceList &List_of_Available_Blocks) {
  int number_of_children; 
  in_file >> block;
  if (block.used) {
     if (List_of_Available_Blocks.Nfree > 0) {
        block.info.cpu = List_of_Available_Blocks.nextCPU();
        block.info.blknum = List_of_Available_Blocks.nextBlock();
        List_of_Available_Blocks.update_next();
     } else {
        cout << "\n " << CFFC_Version() 
             << " Quadtree Read Error: Insufficient number of quadrilateral solution blocks.\n";
     } /* endif */
  } /* endif */
  in_file.setf(ios::skipws);  in_file >> number_of_children;  
  in_file.unsetf(ios::skipws);
  if (!block.used && number_of_children > 0) {
     if (childSW_ptr != NULL) { delete childSW_ptr; childSW_ptr = NULL; }
     childSW_ptr = new QuadTreeBlock;
     childSW_ptr->parent_ptr = this;
     childSW_ptr->read(in_file, List_of_Available_Blocks);
     if (childSE_ptr != NULL) { delete childSE_ptr; childSE_ptr = NULL; }
     childSE_ptr = new QuadTreeBlock;
     childSE_ptr->parent_ptr = this;
     childSE_ptr->read(in_file, List_of_Available_Blocks);
     if (childNW_ptr != NULL) { delete childNW_ptr; childNW_ptr = NULL; }
     childNW_ptr = new QuadTreeBlock;
     childNW_ptr->parent_ptr = this;
     childNW_ptr->read(in_file, List_of_Available_Blocks);
     if (childNE_ptr != NULL) { delete childNE_ptr; childNE_ptr = NULL; }
     childNE_ptr = new QuadTreeBlock;
     childNE_ptr->parent_ptr = this;
     childNE_ptr->read(in_file, List_of_Available_Blocks);
  } else {
     if (childSW_ptr != NULL) { delete childSW_ptr; childSW_ptr = NULL; }
     if (childSE_ptr != NULL) { delete childSE_ptr; childSE_ptr = NULL; }
     if (childNW_ptr != NULL) { delete childNW_ptr; childNW_ptr = NULL; }
     if (childNE_ptr != NULL) { delete childNE_ptr; childNE_ptr = NULL; }
  } /* endif */
}

/*************************************************************
 * QuadTreeBlock::write -- Write quadtree block (recursive). *
 *************************************************************/
inline void QuadTreeBlock::write(ostream &out_file) {
  int number_of_children; 
  out_file << block << "\n";
  if (!block.used && childSW_ptr != NULL) {
    number_of_children = 4;
  } else {
    number_of_children = 0;
  } /* endif */
  out_file << " " << number_of_children << "\n";
  if (childSW_ptr != NULL) childSW_ptr->write(out_file);
  if (childSE_ptr != NULL) childSE_ptr->write(out_file);
  if (childNW_ptr != NULL) childNW_ptr->write(out_file);
  if (childNE_ptr != NULL) childNE_ptr->write(out_file);
}

/*************************************************************
 * QuadTreeBlock::broadcast -- Broadcast quadtree block.     *
 *************************************************************/
inline void QuadTreeBlock::broadcast(void) {
#ifdef _MPI_VERSION
  int number_of_children; 
  Broadcast_Adaptive_Block(block);
  if (CFFC_Primary_MPI_Processor()) {
     if (!block.used && childSW_ptr != NULL) {
       number_of_children = 4;
     } else {
       number_of_children = 0;
     } /* endif */
  } /* endif */
  MPI::COMM_WORLD.Bcast(&number_of_children, 1, MPI::INT, 0);
  if (!block.used && number_of_children > 0) {
     if (!CFFC_Primary_MPI_Processor() && childSW_ptr != NULL) { delete childSW_ptr; childSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childSW_ptr = new QuadTreeBlock; childSW_ptr->parent_ptr = this; }
     childSW_ptr->broadcast();
     if (!CFFC_Primary_MPI_Processor() && childSE_ptr != NULL) { delete childSE_ptr; childSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childSE_ptr = new QuadTreeBlock; childSE_ptr->parent_ptr = this; }
     childSE_ptr->broadcast();
     if (!CFFC_Primary_MPI_Processor() && childNW_ptr != NULL) { delete childNW_ptr; childNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childNW_ptr = new QuadTreeBlock; childNW_ptr->parent_ptr = this; }
     childNW_ptr->broadcast();
     if (!CFFC_Primary_MPI_Processor() && childNE_ptr != NULL) { delete childNE_ptr; childNE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childNE_ptr = new QuadTreeBlock; childNE_ptr->parent_ptr = this; }
     childNE_ptr->broadcast();
  } else {
     if (!CFFC_Primary_MPI_Processor() && childSW_ptr != NULL) { delete childSW_ptr; childSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childSE_ptr != NULL) { delete childSE_ptr; childSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childNW_ptr != NULL) { delete childNW_ptr; childNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childNE_ptr != NULL) { delete childNE_ptr; childNE_ptr = NULL; }
  } /* endif */
#endif
}

inline void QuadTreeBlock::broadcast(AdaptiveBlockResourceList &List_of_Available_Blocks) {
#ifdef _MPI_VERSION
  int number_of_children; 
  Broadcast_Adaptive_Block(block);
  if (!CFFC_Primary_MPI_Processor()) {
     if (block.used) {
        if (List_of_Available_Blocks.Nfree > 0) {
           block.info.cpu = List_of_Available_Blocks.nextCPU();
           block.info.blknum = List_of_Available_Blocks.nextBlock();
           List_of_Available_Blocks.update_next();
        } else {
           cout << "\n " << CFFC_Version() 
                << " Quadtree Broadcast Error: Insufficient number of quadrilateral solution blocks.\n";
        } /* endif */
     } /* endif */
  } /* endif */
  if (CFFC_Primary_MPI_Processor()) {
     if (!block.used && childSW_ptr != NULL) {
       number_of_children = 4;
     } else {
       number_of_children = 0;
     } /* endif */
  } /* endif */
  MPI::COMM_WORLD.Bcast(&number_of_children, 1, MPI::INT, 0);
  if (!block.used && number_of_children > 0) {
     if (!CFFC_Primary_MPI_Processor() && childSW_ptr != NULL) { delete childSW_ptr; childSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childSW_ptr = new QuadTreeBlock; childSW_ptr->parent_ptr = this; }
     childSW_ptr->broadcast(List_of_Available_Blocks);
     if (!CFFC_Primary_MPI_Processor() && childSE_ptr != NULL) { delete childSE_ptr; childSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childSE_ptr = new QuadTreeBlock; childSE_ptr->parent_ptr = this; }
     childSE_ptr->broadcast(List_of_Available_Blocks);
     if (!CFFC_Primary_MPI_Processor() && childNW_ptr != NULL) { delete childNW_ptr; childNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childNW_ptr = new QuadTreeBlock; childNW_ptr->parent_ptr = this; }
     childNW_ptr->broadcast(List_of_Available_Blocks);
     if (!CFFC_Primary_MPI_Processor() && childNE_ptr != NULL) { delete childNE_ptr; childNE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor()) { childNE_ptr = new QuadTreeBlock; childNE_ptr->parent_ptr = this; }
     childNE_ptr->broadcast(List_of_Available_Blocks);
  } else {
     if (!CFFC_Primary_MPI_Processor() && childSW_ptr != NULL) { delete childSW_ptr; childSW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childSE_ptr != NULL) { delete childSE_ptr; childSE_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childNW_ptr != NULL) { delete childNW_ptr; childNW_ptr = NULL; }
     if (!CFFC_Primary_MPI_Processor() && childNE_ptr != NULL) { delete childNE_ptr; childNE_ptr = NULL; }
  } /* endif */
#endif
}

/*************************************************************
 * QuadTreeBlock::maxRefinementLevel -- Maximum refinement   *
 *                                      level (recursive).   *
 *************************************************************/
inline int QuadTreeBlock::maxRefinementLevel(int &Max_Level) {
  if (block.used) {
     Max_Level = max(Max_Level, block.info.level);
  } else if (!block.used && childSW_ptr != NULL) {
     if (childSW_ptr != NULL) childSW_ptr->maxRefinementLevel(Max_Level);
     if (childSE_ptr != NULL) childSE_ptr->maxRefinementLevel(Max_Level);
     if (childNW_ptr != NULL) childNW_ptr->maxRefinementLevel(Max_Level);
     if (childNE_ptr != NULL) childNE_ptr->maxRefinementLevel(Max_Level);
  } /* endif */
  return (Max_Level);
}

/*************************************************************
 * QuadTreeBlock -- Input-output operators.                  *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const QuadTreeBlock &Block) {
  out_file << Block.block;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     QuadTreeBlock &Block) {
  in_file >> Block.block;
  return (in_file);
}

/*************************************************************
 * QuadTreeBlock -- External subroutines.                    *
 *************************************************************/


/********************************************************
 * Class: QuadTreeBlock_DataStructure                   *
 *                                                      *
 * Member functions                                     *
 *   Roots -- Return roots of quadtree data structure.  *
 *   NRi   -- Return number of root blocks in           *
 *            the i-direction (zeta-direction).         *
 *   NRj   -- Return number of root blocks in           *
 *            the j-direction (eta-direction).          *
 *   Blocks -- Return global list of pointers to the    *
 *             blocks in use in the quadtree data       *
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
 *   allocate -- Allocate memory for quadtree data      *
 *               structure.                             *
 *   deallocate -- Deallocate memory for quadtree data  *
 *                 structure.                           *
 *   allocateRoots -- Allocate memory for roots of      *
 *                    quadtree data structure.          *
 *   deallocateRoots -- Deallocate memory for roots of  *
 *                      quadtree data structure.        *
 *   allocateBlocks -- Allocate memory for global block *
 *                     list of quadtree data structure. *
 *   deallocateBlocks -- Deallocate memory for global   *
 *                       block list of quadtree data    *
 *                       structure.                     *
 *   assign_block_pointers, assign_block_ptr            *
 *                    -- Assign block pointers to       *
 *                       all used blocks.               *
 *   renumber -- Renumbers quadtree adaptive blocks,    *
 *               assigning global block numbers.        *
 *   countUsedBlocks -- Returns the number of used      *
 *                      blocks in the quadtree data     *
 *                      structure.                      *
 *   countUsedCells -- Returns the number of            *
 *                     computational cells in the       *
 *                     quadtree data structure.         *
 *   efficiencyRefinement -- Returns the ratio of the   *
 *                           current number of cells    *
 *                           used to the number of cells*
 *                           that would have been used  *
 *                           with a uniform mesh.       *
 *   getRoot -- Return indices of root quadtree block.  *
 *   getNeighbour -- Return neighbour of quadtree block *
 *                   in specified search direction of   *
 *                   interest.                          *
 *   findRootNeighbours -- Find neighbouring blocks of  *
 *                         all quadtree root blocks.    *
 *   findNeighbours -- Find neighbouring blocks of      *
 *                     quadtree block.                  *
 *   refineBlock   -- Refines (divides) a quadtree      *
 *                    block into four offspring.        *
 *   coarsenBlocks -- Coarsens (contracts) four quadree *
 *                    blocks into a single parent.      *
 *   nochangeAll-- Sets the mesh refinement flags to    *
 *                 force no refinement or coarsening of *
 *                 quadtree blocks (default).           *
 *   refineAll  -- Sets the mesh refinement flags to    *
 *                 force refinement of all quadtree     *
 *                 blocks.                              *
 *   coarsenAll -- Sets the mesh refinement flags to    *
 *                 force coasening of all quadtree      *
 *                 blocks.                              *
 *   setRefineAll -- Sets the mesh refinement flags of  *
 *                   all quadtree blocks to specified   *
 *                   value.                             *
 *   highestRefinementLevel -- Return highest level     *
 *                of block refinement.                  *
 *   numberToBeRefined -- Return number of blocks to be *
 *                refined.                              *
 *   numberToBeCoarsened -- Return number of blocks to  *
 *                be coarsened.                         * 
 *                                                      *
 * Member operators                                     *
 *      QT -- quadtree hierarchical data structure      *
 *                                                      *
 * QT = QT;                                             *
 * cout << QT; (output function)                        *
 * cin  >> QT; (input function)                         *
 *                                                      *
 ********************************************************/
class QuadTreeBlock_DataStructure{
  private:
  public:
    int                    NRi; // Number of roots in i-direction.
    int                    NRj; // Number of roots in j-direction.
    QuadTreeBlock      **Roots; // Roots of quadtree data structure.
    int                   Ncpu; // Number of CPUs available.
    int                   Nblk; // Number of local blocks per CPU.
    QuadTreeBlock    ***Blocks; // Global list of pointers to blocks
                                // in use in quadtree data structure.
    int          **RefineFlags; // Global list of refinement flags.
    double     RefineThreshold, // Thresholds for refinment and coarsening.
              CoarsenThreshold; //
    int MaximumRefinementLevel; // Maximum allowable refinement level.
    int MinimumRefinementLevel; // Maximum allowable refinement level.
                                // for blocks in data structures.
	                        // Made public so can access them.

    /* Creation, copy, and assignment constructors. */
    QuadTreeBlock_DataStructure(void) {
       NRi = 0; NRj = 0; Roots = NULL;
       Ncpu = 0; Nblk = 0; Blocks = NULL; RefineFlags = NULL;
       MaximumRefinementLevel = 99; MinimumRefinementLevel = 0;
       RefineThreshold = 0.50; CoarsenThreshold = 0.10;
    }

    QuadTreeBlock_DataStructure(const QuadTreeBlock_DataStructure &QT) {
       NRi = QT.NRi; NRj = QT.NRj; Roots = QT.Roots;
       Ncpu = QT.Ncpu; Nblk = QT.Nblk; Blocks = QT.Blocks;
       RefineFlags = QT.RefineFlags;
       MaximumRefinementLevel = QT.MaximumRefinementLevel; 
       MinimumRefinementLevel = QT.MinimumRefinementLevel; 
       RefineThreshold = QT.RefineThreshold; 
       CoarsenThreshold = QT.CoarsenThreshold;
    }

    /* Destructor. */
    // ~QuadTreeBlock_DataStructure(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // QuadTreeBlock_DataStructure operator = (const QuadTreeBlock_DataStructure &QT);
    // Use automatically generated assignment operator.

    /* Allocate memory for quadtree data structure. */
    void allocate(const int ni, const int nj, 
                  const int ncpu, const int nblk);

    /* Deallocate memory for quadtree data structure. */
    void deallocate(void);

    /* Allocate memory for quadtree data structure roots. */
    void allocateRoots(const int ni, const int nj);

    /* Deallocate memory for quadtree data structure roots. */
    void deallocateRoots(void);

    /* Allocate memory for quadtree data structure global block list. */
    void allocateBlocks(const int ncpu, const int nblk);

    /* Deallocate memory for quadtree data structure global block list. */
    void deallocateBlocks(void);

    /* Assign block pointers to all used blocks. */
    void assign_block_ptr(QuadTreeBlock *Block_Ptr);
    void assign_block_pointers(void);

    /* Renumber all quadtree blocks. */
    void renumber(void);

    /* Count number of used blocks and cells in quadtree data structure 
       and other useful mesh statistics. */
    int countUsedBlocks(void);
    int countUsedCells(void);
    double efficiencyRefinement(void);

    /* Quadtree block root indices. */
    AdaptiveBlock2D_Dimensions getRoot(QuadTreeBlock *Block_Ptr);
    AdaptiveBlock2D_Dimensions getRootIndices(QuadTreeBlock *Block_Ptr);

    /* Quadtree block indices. */
    AdaptiveBlock2D_Dimensions getBlockIndices(QuadTreeBlock *Block_Ptr);

    /* Find neighbour of quadtree block in specified search direction. */
    QuadTreeBlock *getNeighbour(QuadTreeBlock *Block_Ptr,
                                const int Search_Dir_Mask);

    /* Find neighbouring blocks of all quadtree root blocks. */
    void findRootNeighbours(void);

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

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const QuadTreeBlock_DataStructure &QT);
    friend istream &operator >> (istream &in_file,
				 QuadTreeBlock_DataStructure &QT);
    
};

/*********************************************************************
 * QuadTreeBlock_DataStructure::allocate -- Allocate memory.         *
 *********************************************************************/
inline void QuadTreeBlock_DataStructure::allocate(const int ni, 
                                                  const int nj, 
                                                  const int ncpu, 
                                                  const int nblk) {
   int i, j; assert( ni > 0 && nj > 0 && ncpu > 0 && nblk > 0 );
   NRi = ni; NRj = nj; Ncpu = ncpu; Nblk = nblk;
   Roots = new QuadTreeBlock*[NRi];
   for ( i = 0; i <= NRi-1 ; ++i ) Roots[i] = new QuadTreeBlock[NRj];
   Blocks = new QuadTreeBlock**[Ncpu];
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      Blocks[i] = new QuadTreeBlock*[Nblk];
      for ( j = 0; j <= Nblk-1; ++j) Blocks[i][j] = NULL; 
   } /* endfor */
   RefineFlags = new int*[Ncpu];
   for ( i = 0; i <= Ncpu-1 ; ++i ) RefineFlags[i] = new int[Nblk];
   nochangeAll();
}

inline void QuadTreeBlock_DataStructure::allocateRoots(const int ni, 
                                                       const int nj) {
   int i; assert( ni > 0 && nj > 0 );
   NRi = ni; NRj = nj; Roots = new QuadTreeBlock*[NRi];
   for ( i = 0; i <= NRi-1 ; ++i ) Roots[i] = new QuadTreeBlock[NRj];
}

inline void QuadTreeBlock_DataStructure::allocateBlocks(const int ncpu, 
                                                        const int nblk) {
   int i, j; assert( ncpu > 0 && nblk > 0 );
   Ncpu = ncpu; Nblk = nblk; Blocks = new QuadTreeBlock**[Ncpu];
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      Blocks[i] = new QuadTreeBlock*[Nblk];
      for ( j = 0; j <= Nblk-1; ++j) Blocks[i][j] = NULL; 
   } /* endfor */
   RefineFlags = new int*[Ncpu];
   for ( i = 0; i <= Ncpu-1 ; ++i ) RefineFlags[i] = new int[Nblk];
   nochangeAll();
}

/*********************************************************************
 * QuadTreeBlock_DataStructure::deallocate -- Deallocate memory.     *
 *********************************************************************/
inline void QuadTreeBlock_DataStructure::deallocate(void) {
   int i;
   for ( i = 0; i <= NRi-1 ; ++i ) {
      delete []Roots[i]; Roots[i] = NULL;
   } /* endfor */
   delete []Roots; Roots = NULL;
   NRi = 0; NRj = 0;
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      delete []Blocks[i]; Blocks[i] = NULL;
   } /* endfor */
   delete []Blocks; Blocks = NULL;
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      delete []RefineFlags[i]; RefineFlags[i] = NULL;
   } /* endfor */
   delete []RefineFlags; RefineFlags = NULL;
   Ncpu = 0; Nblk = 0;
}

inline void QuadTreeBlock_DataStructure::deallocateRoots(void) {
   int i;
   for ( i = 0; i <= NRi-1 ; ++i ) {
      delete []Roots[i]; Roots[i] = NULL;
   } /* endfor */
   delete []Roots; Roots = NULL;
   NRi = 0; NRj = 0;
}

inline void QuadTreeBlock_DataStructure::deallocateBlocks(void) {
   int i;
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      delete []Blocks[i]; Blocks[i] = NULL;
   } /* endfor */
   delete []Blocks; Blocks = NULL;
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      delete []RefineFlags[i]; RefineFlags[i] = NULL;
   } /* endfor */
   delete []RefineFlags; RefineFlags = NULL;
   Ncpu = 0; Nblk = 0;
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::assign_block_ptr -- Assign block ptr.       *
 ****************************************************************************/
inline void QuadTreeBlock_DataStructure::assign_block_ptr(QuadTreeBlock *Block_Ptr) {
   int icpu, iblk;
   if (Block_Ptr != NULL) {
      if (Block_Ptr->block.used) {
	 icpu = Block_Ptr->block.info.cpu;
         iblk = Block_Ptr->block.info.blknum;
         if (icpu < 0 || icpu >= Ncpu ||
             iblk < 0 || iblk >= Nblk) {
            cout << "\n " << CFFC_Version() 
                 << " Quadtree Block Error: Invalid range for block indices.\n";
         } else {
            Blocks[icpu][iblk] = Block_Ptr;
         } /* endif */
      } else if (Block_Ptr->childSW_ptr != NULL) {
	 assign_block_ptr(Block_Ptr->childSW_ptr);
	 assign_block_ptr(Block_Ptr->childSE_ptr);
	 assign_block_ptr(Block_Ptr->childNW_ptr);
	 assign_block_ptr(Block_Ptr->childNE_ptr);
      } else {
	 // Not sure what to do!  This shouldn't happen.
      } /* endif */
   } /* endif */
}

inline void QuadTreeBlock_DataStructure::assign_block_pointers(void) {
   int iBLK, jBLK;
   for ( jBLK = 0; jBLK <= Nblk-1; ++jBLK) {
      for ( iBLK = 0; iBLK <= Ncpu-1 ; ++iBLK ) {
         Blocks[iBLK][jBLK] = NULL;
      } /* endfor */
   } /* endfor */
   for ( jBLK = 0 ; jBLK <= NRj-1 ; ++jBLK ) {
      for ( iBLK = 0 ; iBLK <= NRi-1 ; ++iBLK ) {
 	 assign_block_ptr(&(Roots[iBLK][jBLK]));
      } /* endfor */
   } /* endfor */
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::renumber -- Renumber quadtree blocks.       *
 ****************************************************************************/
inline void QuadTreeBlock_DataStructure::renumber(void) {
   int iCPU, nBLK, global_block_number;
   global_block_number = 0;
   for ( nBLK = 0 ; nBLK <= Nblk-1 ; ++nBLK ) {
      for ( iCPU = 0 ; iCPU <= Ncpu-1 ; ++iCPU ) {
         if (Blocks[iCPU][nBLK] != NULL) {
            if (Blocks[iCPU][nBLK]->block.used) {
               Blocks[iCPU][nBLK]->block.gblknum = global_block_number;
               ++global_block_number;
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::countUsedBlocks -- Number of used blocks.   *
 ****************************************************************************/
inline int QuadTreeBlock_DataStructure::countUsedBlocks(void) {
  int iCPU, nBLK, number_of_used_blocks;
  number_of_used_blocks = 0;
  for ( nBLK = 0 ; nBLK <= Nblk-1 ; ++nBLK ) {
      for ( iCPU = 0 ; iCPU <= Ncpu-1 ; ++iCPU ) {
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
 * QuadTreeBlock_DataStructure::countUsedCells -- Number of used cells.     *
 ****************************************************************************/
inline int QuadTreeBlock_DataStructure::countUsedCells(void) {
  int iCPU, nBLK, number_of_used_cells;
  number_of_used_cells = 0;
  for ( nBLK = 0 ; nBLK <= Nblk-1 ; ++nBLK ) {
      for ( iCPU = 0 ; iCPU <= Ncpu-1 ; ++iCPU ) {
         if (Blocks[iCPU][nBLK] != NULL) {
            if (Blocks[iCPU][nBLK]->block.used) {
               number_of_used_cells += Blocks[iCPU][nBLK]->block.info.dimen.i*
                                       Blocks[iCPU][nBLK]->block.info.dimen.j;
            } /* endif */
         } /* endif */
      } /* endfor */
  } /* endfor */
  return (number_of_used_cells);
}

/******************************************************************************
 * QuadTreeBlock_DataStructure::efficiencyRefinement -- Refinment efficiency. *
 ******************************************************************************/
inline double QuadTreeBlock_DataStructure::efficiencyRefinement(void) {
  int iCPU, nBLK, number_of_used_cells, number_of_cells_with_uniform_mesh, max_level;
  number_of_used_cells = 0; number_of_cells_with_uniform_mesh = 0;
  max_level = highestRefinementLevel();
  for ( nBLK = 0 ; nBLK <= Nblk-1 ; ++nBLK ) {
      for ( iCPU = 0 ; iCPU <= Ncpu-1 ; ++iCPU ) {
         if (Blocks[iCPU][nBLK] != NULL) {
            if (Blocks[iCPU][nBLK]->block.used) {
               number_of_used_cells += Blocks[iCPU][nBLK]->block.info.dimen.i*
                                       Blocks[iCPU][nBLK]->block.info.dimen.j;
               number_of_cells_with_uniform_mesh +=
                                       int(pow(double(4), double(max_level-Blocks[iCPU][nBLK]->block.info.level)))*
                                       (Blocks[iCPU][nBLK]->block.info.dimen.i*
                                        Blocks[iCPU][nBLK]->block.info.dimen.j);
            } /* endif */
         } /* endif */
      } /* endfor */
  } /* endfor */
  return (ONE-double(number_of_used_cells)/double(number_of_cells_with_uniform_mesh));
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::getRoot -- Get root indices.                *
 ****************************************************************************/
inline AdaptiveBlock2D_Dimensions QuadTreeBlock_DataStructure::getRoot(QuadTreeBlock *Block_Ptr) {
  int iBLK, jBLK;
  for ( jBLK = 0 ; jBLK <= NRj-1 ; ++jBLK ) {
      for ( iBLK = 0 ; iBLK <= NRi-1 ; ++iBLK ) {
         if (Block_Ptr == &(Roots[iBLK][jBLK])) {
            return (AdaptiveBlock2D_Dimensions(iBLK, jBLK, 0));
         } /* endif */
      } /* endfor */
  } /* endfor */
  return (AdaptiveBlock2D_Dimensions(-1, -1, 0));
}

inline AdaptiveBlock2D_Dimensions QuadTreeBlock_DataStructure::getRootIndices(QuadTreeBlock *Block_Ptr) {
  QuadTreeBlock *Block_Ptr2 = Block_Ptr;
  if (Block_Ptr2->parent_ptr != NULL &&
      Block_Ptr2->block.info.sector != QUADTREE_SECTOR_NONE) {
    // The block is not a root block.  Descend tree.
    return getRootIndices(Block_Ptr2->parent_ptr);
  }
  return getRoot(Block_Ptr2); 
}

/**********************************************************************
 * QuadTreeBlock_DataStructure::getBlockIndices -- Get block indices. *
 **********************************************************************/
inline AdaptiveBlock2D_Dimensions QuadTreeBlock_DataStructure::getBlockIndices(QuadTreeBlock *Block_Ptr) {
  if (Block_Ptr->parent_ptr == NULL) return AdaptiveBlock2D_Dimensions(0,0,0);
  AdaptiveBlock2D_Dimensions dimensions;
  dimensions.i = max(Block_Ptr->block.info.sector,0)%2;
  dimensions.j = max(Block_Ptr->block.info.sector,0)/2;
  if (Block_Ptr->parent_ptr == NULL) return dimensions;
  dimensions += 2*getBlockIndices(Block_Ptr->parent_ptr);
  return dimensions;
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::getNeighbour -- Return neighbour.           *
 ****************************************************************************/
inline QuadTreeBlock *QuadTreeBlock_DataStructure::getNeighbour(QuadTreeBlock *Block_Ptr,
                                                                const int Search_Dir_Mask) {
  int iBLK, jBLK;
  QuadTreeBlock *neighbour_block_ptr; AdaptiveBlock2D_Dimensions index; 
  AdaptiveBlock2D_Dimensions search_direction;
  if (Block_Ptr->parent_ptr == NULL || // Block has no parent.  Block is a root.  
      Block_Ptr->block.info.sector == QUADTREE_SECTOR_NONE) {
     index = getRoot(Block_Ptr);  // Get indices of root block.
     switch(Search_Dir_Mask) { // Find neighbouring root block.  Check each direction.
       case QUADTREE_DIRECTION_MASK_EAST :
         if (index.i < NRi-1) {
            if (Block_Ptr->block.nE > 0) {
               neighbour_block_ptr = &(Roots[index.i+1][index.j]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else {
	    if (Block_Ptr->block.nE > 0) {
               neighbour_block_ptr = &(Roots[0][index.j]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } /* endif */
         break;
       case QUADTREE_DIRECTION_MASK_WEST :
         if (index.i > 0) {
            if (Block_Ptr->block.nW > 0) {
               neighbour_block_ptr = &(Roots[index.i-1][index.j]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else {
	    if (Block_Ptr->block.nW > 0) {
               neighbour_block_ptr = &(Roots[NRi-1][index.j]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } /* endif */
         break;
       case QUADTREE_DIRECTION_MASK_NORTH :
         if (index.j < NRj-1) {
            if (Block_Ptr->block.nN > 0) {
               neighbour_block_ptr = &(Roots[index.i][index.j+1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else {
            if (Block_Ptr->block.nN > 0) {
               neighbour_block_ptr = &(Roots[index.i][0]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } /* endif */
         break;
       case QUADTREE_DIRECTION_MASK_SOUTH :
         if (index.j > 0) {
            if (Block_Ptr->block.nS > 0) {
               neighbour_block_ptr = &(Roots[index.i][index.j-1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else {
            if (Block_Ptr->block.nS > 0) {
               neighbour_block_ptr = &(Roots[index.i][NRj-1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } /* endif */
         break;
       case QUADTREE_DIRECTION_MASK_SOUTHWEST :
         if (index.i > 0 && index.j > 0) {
            if (Block_Ptr->block.nSW > 0) {
               neighbour_block_ptr = &(Roots[index.i-1][index.j-1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else if (index.j > 0) {
            if (Block_Ptr->block.nSW > 0) {
               neighbour_block_ptr = &(Roots[NRi-1][index.j-1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else if (index.i > 0) {
            if (Block_Ptr->block.nSW > 0) {
               neighbour_block_ptr = &(Roots[index.i-1][NRj-1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else {
            if (Block_Ptr->block.nSW > 0) {
               neighbour_block_ptr = &(Roots[NRi-1][NRj-1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } /* endif */
         break;
       case QUADTREE_DIRECTION_MASK_SOUTHEAST :
         if (index.i < NRi-1 && index.j > 0) {
           if (Block_Ptr->block.nSE > 0) {
               neighbour_block_ptr = &(Roots[index.i+1][index.j-1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else if (index.j > 0) {
           if (Block_Ptr->block.nSE > 0) {
               neighbour_block_ptr = &(Roots[0][index.j-1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else if (index.i < NRi-1) {
           if (Block_Ptr->block.nSE > 0) {
               neighbour_block_ptr = &(Roots[index.i+1][NRj-1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else {
           if (Block_Ptr->block.nSE > 0) {
               neighbour_block_ptr = &(Roots[0][NRj-1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } /* endif */
         break;
       case QUADTREE_DIRECTION_MASK_NORTHWEST :
         if (index.i > 0 && index.j < NRj-1) {
           if (Block_Ptr->block.nNW > 0) {
               neighbour_block_ptr = &(Roots[index.i-1][index.j+1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else if (index.j < NRj-1) {
           if (Block_Ptr->block.nNW > 0) {
               neighbour_block_ptr = &(Roots[NRi-1][index.j+1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else if (index.i > 0) {
           if (Block_Ptr->block.nNW > 0) {
               neighbour_block_ptr = &(Roots[index.i-1][0]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else {
           if (Block_Ptr->block.nNW > 0) {
               neighbour_block_ptr = &(Roots[NRi-1][0]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } /* endif */
         break;
       case QUADTREE_DIRECTION_MASK_NORTHEAST :
         if (index.i < NRi-1 && index.j < NRj-1) {
           if (Block_Ptr->block.nNE > 0) {
               neighbour_block_ptr = &(Roots[index.i+1][index.j+1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else if (index.j < NRj-1) {
           if (Block_Ptr->block.nNE > 0) {
               neighbour_block_ptr = &(Roots[0][index.j+1]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
	 } else if (index.i < NRi-1) {
           if (Block_Ptr->block.nNE > 0) {
               neighbour_block_ptr = &(Roots[index.i+1][0]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } else {
           if (Block_Ptr->block.nNE > 0) {
               neighbour_block_ptr = &(Roots[0][0]);
            } else {
               neighbour_block_ptr = NULL;
            } /* endif */
         } /* endif */
         break;
       default:
         neighbour_block_ptr = NULL;
         break;
     } /* endswitch */
     if (Block_Ptr->block.used &&
         neighbour_block_ptr != NULL) { // Descend neigbouring root block until used children are encountered.
        search_direction = Block_Ptr->search_dir(Search_Dir_Mask);
        while ((neighbour_block_ptr->child_ptr(QUADTREE_SECTOR_SW-
                                               search_direction.i-
                                               search_direction.j+
                                               (max(-search_direction.i*search_direction.j, 0)/2)*
                                               max(search_direction.i, search_direction.j)) != NULL) &&
	       !(neighbour_block_ptr->block.used)) {
           neighbour_block_ptr = neighbour_block_ptr->child_ptr(QUADTREE_SECTOR_SW-
                                                                search_direction.i-
                                                                search_direction.j+
                                                                (max(-search_direction.i*search_direction.j, 0)/2)*
                                                                max(search_direction.i, search_direction.j));
        } /* endwhile */
     } /* endif */ // Finally, return resulting neighbour.
  } else { // Block has a parent.  Block is not a root.
     if (Block_Ptr->sibling(Search_Dir_Mask)) { // Neighbour is a sibling.  Get it.
        search_direction = Block_Ptr->search_dir(Search_Dir_Mask);
        neighbour_block_ptr = Block_Ptr->parent_ptr->child_ptr(Block_Ptr->block.info.sector+
                                                               search_direction.i+
                                                               search_direction.j);
        if (Block_Ptr->block.used) { // If block is not used, don't descend any further!
           if ((neighbour_block_ptr->child_ptr(QUADTREE_SECTOR_SW-
                                               search_direction.i-
                                               search_direction.j+
                                               (max(-search_direction.i*search_direction.j, 0)/2)*
                                               max(search_direction.i, search_direction.j)) == NULL) ||
               (neighbour_block_ptr->block.used)) { // Sibling is used and has no children.
	      neighbour_block_ptr = neighbour_block_ptr; // Return sibling.
	   } else { // Return child of sibling.
              neighbour_block_ptr = neighbour_block_ptr->child_ptr(QUADTREE_SECTOR_SW-
                                                                   search_direction.i-
                                                                   search_direction.j+
                                                                   (max(-search_direction.i*search_direction.j, 0)/2)*
                                                                   max(search_direction.i, search_direction.j));
           } /* endif */
        } /* endif */
     } else if (Search_Dir_Mask <= QUADTREE_DIRECTION_MASK_WEST) { // Neighbour is not a sibling.
        neighbour_block_ptr = getNeighbour(Block_Ptr->parent_ptr, Search_Dir_Mask); // Move up the tree.
        if (neighbour_block_ptr == NULL) { // Parent of neighbour does not exist.
	   neighbour_block_ptr = NULL; // Return null pointer.
        } else {
 	   if (neighbour_block_ptr->block.info.level > Block_Ptr->parent_ptr->block.info.level) { // If neighbour is more refined than parent, get parent of neighbour.
              neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
           } /* endif */
           search_direction = Block_Ptr->search_dir(Search_Dir_Mask);
           if ((neighbour_block_ptr->child_ptr(Block_Ptr->block.info.sector-
                                               search_direction.i-
                                               search_direction.j) == NULL) ||
               (neighbour_block_ptr->block.used)) { // Neighbour is used and has no children.
	      neighbour_block_ptr = neighbour_block_ptr; // Return neighbour.
           } else { // Return child of neighbour.
              neighbour_block_ptr = neighbour_block_ptr->child_ptr(Block_Ptr->block.info.sector-
                                                                   search_direction.i-
                                                                   search_direction.j);
              if ((!neighbour_block_ptr->block.used) &&
                  (neighbour_block_ptr->child_ptr(QUADTREE_SECTOR_SW-
                                                  search_direction.i-
                                                  search_direction.j) != NULL)) { // Child of neighbour is not used, return child of child.
                 neighbour_block_ptr = neighbour_block_ptr->child_ptr(QUADTREE_SECTOR_SW-
                                                                      search_direction.i-
                                                                      search_direction.j);
              } /* endif */
           } /* endif */
        } /* endif */
     } else { // Neighbour is not a sibling and is a corner neighbour.  Perform two-direction search.
        switch(Search_Dir_Mask) {
	  case QUADTREE_DIRECTION_MASK_SOUTHWEST : // SOUTHWEST CORNER NEIGHBOUR
            if (Block_Ptr->sibling(QUADTREE_DIRECTION_MASK_WEST)) {
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_WEST);
               if (neighbour_block_ptr != NULL) { 
		  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_SOUTH);
                  if (neighbour_block_ptr != NULL) {
                     if ((neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_EAST))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_EAST);
                     } /* endif */
                  } /* endif */
               } /* endif */
            } else if (Block_Ptr->sibling(QUADTREE_DIRECTION_MASK_SOUTH)) {
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_SOUTH);
               if (neighbour_block_ptr != NULL) { 
                  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_WEST);
                  if (neighbour_block_ptr != NULL) {
                     if ((neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_NORTH))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_NORTH);
                     } /* endif */
                  } /* endif */
               } /* endif */
            } else {
	       // Usual approach
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_SOUTH);
               if (neighbour_block_ptr != NULL) {
                  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_WEST);
                  if (neighbour_block_ptr != NULL) {
                     if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_NORTH))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_NORTH);
                     } /* endif */
                  } /* endif */
               } else {
                  neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_WEST);
                  if (neighbour_block_ptr != NULL) {
                     if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                        neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                     } /* endif */
                     neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_SOUTH);
                     if (neighbour_block_ptr != NULL) {
                        if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) &&
                            (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_EAST))) {
                           neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_EAST);
                        } /* endif */
                     } /* endif */
                  } /* endif */
               } /* endif */
	       // NASA Rotor grid corner neighbours
/*                neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_SOUTH); */
/*                if (neighbour_block_ptr != NULL) { */
/*                   if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent. */
/*                      neighbour_block_ptr = neighbour_block_ptr->parent_ptr; */
/*                   } /\* endif *\/ */
/*                   neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_WEST); */
/*                   if (neighbour_block_ptr != NULL) { */
/*                      neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_SOUTH); */
/*                      if (neighbour_block_ptr != NULL) { */
/*                         if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) && */
/*                             (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_EAST))) { */
/*                            neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_EAST); */
/*                         } /\* endif *\/ */
/*                      } /\* endif *\/ */
/*                   } /\* endif *\/ */
/*                } /\* endif *\/ */
            } /* endif */
            break;
	  case QUADTREE_DIRECTION_MASK_SOUTHEAST : // SOUTHEAST CORNER NEIGHBOUR
            if (Block_Ptr->sibling(QUADTREE_DIRECTION_MASK_EAST)) {
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_EAST);
               if (neighbour_block_ptr != NULL) { 
                  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_SOUTH);
                  if (neighbour_block_ptr != NULL) {
                     if ((neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_WEST))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_WEST);
                     } /* endif */
                  } /* endif */
               } /* endif */
            } else if (Block_Ptr->sibling(QUADTREE_DIRECTION_MASK_SOUTH)) {
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_SOUTH);
               if (neighbour_block_ptr != NULL) { 
                  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_EAST);
                  if (neighbour_block_ptr != NULL) {
                     if ((neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_NORTH))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_NORTH);
                     } /* endif */
                  } /* endif */
               } /* endif */
            } else {
	       // Usual approach
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_SOUTH);
               if (neighbour_block_ptr != NULL) {
                  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_EAST);
                  if (neighbour_block_ptr != NULL) {
                     if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_NORTH))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_NORTH);
                     } /* endif */
                  } /* endif */
               } else {
                  neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_EAST);
                  if (neighbour_block_ptr != NULL) {
                     if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                        neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                     } /* endif */
                     neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_SOUTH);
                     if (neighbour_block_ptr != NULL) {
                        if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) &&
                            (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_WEST))) {
                           neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_WEST);
                        } /* endif */
                     } /* endif */
                  } /* endif */
               } /* endif */
	       // NASA Rotor grid corner neighbours
/*                neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_SOUTH); */
/*                if (neighbour_block_ptr != NULL) { */
/*                   if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent. */
/*                      neighbour_block_ptr = neighbour_block_ptr->parent_ptr; */
/*                   } /\* endif *\/ */
/*                   neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_EAST); */
/*                   if (neighbour_block_ptr != NULL) { */
/*                      neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_SOUTH); */
/*                      if (neighbour_block_ptr != NULL) { */
/*                         if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) && */
/*                             (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_WEST))) { */
/*                            neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_WEST); */
/*                         } /\* endif *\/ */
/*                      } /\* endif *\/ */
/*                   } /\* endif *\/ */
/*                } /\* endif *\/ */
            } /* endif */
            break;
	  case QUADTREE_DIRECTION_MASK_NORTHWEST : // NORTHWEST CORNER NEIGHBOUR
            if (Block_Ptr->sibling(QUADTREE_DIRECTION_MASK_WEST)) {
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_WEST);
               if (neighbour_block_ptr != NULL) { 
                  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_NORTH);
                  if (neighbour_block_ptr != NULL) {
                     if ((neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_EAST))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_EAST);
                     } /* endif */
                  } /* endif */
               } /* endif */
            } else if (Block_Ptr->sibling(QUADTREE_DIRECTION_MASK_NORTH)) {
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_NORTH);
               if (neighbour_block_ptr != NULL) { 
                  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_WEST);
                  if (neighbour_block_ptr != NULL) {
                     if ((neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_SOUTH))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_SOUTH);
                     } /* endif */
                  } /* endif */
               } /* endif */
            } else {
	       // Usual approach
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_NORTH);
               if (neighbour_block_ptr != NULL) {
                  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_WEST);
                  if (neighbour_block_ptr != NULL) {
                     if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_SOUTH))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_SOUTH);
                     } /* endif */
                  } /* endif */
               } else {
                  neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_WEST);
                  if (neighbour_block_ptr != NULL) {
                     if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                        neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                     } /* endif */
                     neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_NORTH);
                     if (neighbour_block_ptr != NULL) {
                        if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) &&
                            (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_EAST))) {
                           neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_EAST);
                        } /* endif */
                     } /* endif */
                  } /* endif */
               } /* endif */
	       // NASA Rotor grid corner neighbours
/*                neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_NORTH); */
/*                if (neighbour_block_ptr != NULL) { */
/*                   if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent. */
/*                      neighbour_block_ptr = neighbour_block_ptr->parent_ptr; */
/*                   } /\* endif *\/ */
/*                   neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_WEST); */
/*                   if (neighbour_block_ptr != NULL) { */
/*                      neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_NORTH); */
/*                      if (neighbour_block_ptr != NULL) { */
/*                         if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) && */
/*                             (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_EAST))) { */
/*                            neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_EAST); */
/*                         } /\* endif *\/ */
/*                      } /\* endif *\/ */
/*                   } /\* endif *\/ */
/*                } /\* endif *\/ */
            } /* endif */
            break;
 	  case QUADTREE_DIRECTION_MASK_NORTHEAST : // NORTHEAST CORNER NEIGHBOUR
            if (Block_Ptr->sibling(QUADTREE_DIRECTION_MASK_EAST)) {
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_EAST);
               if (neighbour_block_ptr != NULL) {
                  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_NORTH);
                  if (neighbour_block_ptr != NULL) { 
                     if ((neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_WEST))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_WEST);
                     } /* endif */
                  } /* endif */
               } /* endif */
            } else if (Block_Ptr->sibling(QUADTREE_DIRECTION_MASK_NORTH)) {
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_NORTH);
               if (neighbour_block_ptr != NULL) { 
                  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_EAST);
                  if (neighbour_block_ptr != NULL) { 
                     if ((neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_SOUTH))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_SOUTH);
                     } /* endif */
                  } /* endif */
               } /* endif */
            } else {
	       // Usual approach
               neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_NORTH);
               if (neighbour_block_ptr != NULL) {
                  if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                     neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                  } /* endif */
                  neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_EAST);
                  if (neighbour_block_ptr != NULL) {
                     if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) &&
                         (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_SOUTH))) {
                        neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_SOUTH);
                     } /* endif */
                  } /* endif */
               } else {
                  neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_EAST);
                  if (neighbour_block_ptr != NULL) {
                     if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent.
                        neighbour_block_ptr = neighbour_block_ptr->parent_ptr;
                     } /* endif */
                     neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_NORTH);
                     if (neighbour_block_ptr != NULL) {
                        if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) &&
                            (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_WEST))) {
                           neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_WEST);
                        } /* endif */
                     } /* endif */
                  } /* endif */
               } /* endif */
	       // NASA Rotor grid corner neighbours
/*                neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_NORTH); */
/*                if (neighbour_block_ptr != NULL) { */
/*                   if (neighbour_block_ptr->block.info.level > Block_Ptr->block.info.level) { // If more refined, get parent. */
/*                      neighbour_block_ptr = neighbour_block_ptr->parent_ptr; */
/*                   } /\* endif *\/ */
/*                   neighbour_block_ptr = getNeighbour(Block_Ptr, QUADTREE_DIRECTION_MASK_EAST); */
/*                   if (neighbour_block_ptr != NULL) { */
/*                      neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_NORTH); */
/*                      if (neighbour_block_ptr != NULL) {  */
/*                         if ((neighbour_block_ptr->block.info.level >= Block_Ptr->block.info.level) && */
/*                             (neighbour_block_ptr->sibling(QUADTREE_DIRECTION_MASK_WEST))) { */
/*                            neighbour_block_ptr = getNeighbour(neighbour_block_ptr, QUADTREE_DIRECTION_MASK_WEST); */
/*                         } /\* endif *\/ */
/*                      } /\* endif *\/ */
/*                   } /\* endif *\/ */
/*                } /\* endif *\/ */
            } /* endif */
            break;
          default:
            neighbour_block_ptr = NULL;
            break;
        } /* endswitch */
     } /* endif */
  } /* endif */
  return (neighbour_block_ptr);
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::findRootNeighbours -- Find root neighbours. *
 ****************************************************************************/
inline void QuadTreeBlock_DataStructure::findRootNeighbours(void) {
  int iBLK, jBLK;
  for ( jBLK = 0 ; jBLK <= NRj-1 ; ++jBLK ) {
      for ( iBLK = 0 ; iBLK <= NRi-1 ; ++iBLK ) {
	  // First determine the number of neighbours based
	  // on (i,j) index location of the block.
          if (iBLK == 0 && iBLK == NRi-1) {
             Roots[iBLK][jBLK].block.nW = 0;
             Roots[iBLK][jBLK].block.nNW = 0;
             Roots[iBLK][jBLK].block.nSW = 0;
             Roots[iBLK][jBLK].block.nE = 0;
             Roots[iBLK][jBLK].block.nNE = 0;
             Roots[iBLK][jBLK].block.nSE = 0;
          } else if (iBLK == 0) {
             Roots[iBLK][jBLK].block.nW = 0;
             Roots[iBLK][jBLK].block.nNW = 0;
             Roots[iBLK][jBLK].block.nSW = 0;
             Roots[iBLK][jBLK].block.nE = 1;
             Roots[iBLK][jBLK].block.nNE = 1;
             Roots[iBLK][jBLK].block.nSE = 1;
          } else if (iBLK == NRi-1) {
             Roots[iBLK][jBLK].block.nW = 1;
             Roots[iBLK][jBLK].block.nNW = 1;
             Roots[iBLK][jBLK].block.nSW = 1;
             Roots[iBLK][jBLK].block.nE = 0;
             Roots[iBLK][jBLK].block.nNE = 0;
             Roots[iBLK][jBLK].block.nSE = 0;
          } else {
             Roots[iBLK][jBLK].block.nW = 1;
             Roots[iBLK][jBLK].block.nNW = 1;
             Roots[iBLK][jBLK].block.nSW = 1;
             Roots[iBLK][jBLK].block.nE = 1;
             Roots[iBLK][jBLK].block.nNE = 1;
             Roots[iBLK][jBLK].block.nSE = 1;
          } /* endif */
          if (jBLK == 0 && jBLK == NRj-1) {
             Roots[iBLK][jBLK].block.nS = 0;
             Roots[iBLK][jBLK].block.nSW = 0;
             Roots[iBLK][jBLK].block.nSE = 0;
             Roots[iBLK][jBLK].block.nN = 0;
             Roots[iBLK][jBLK].block.nNW = 0;
             Roots[iBLK][jBLK].block.nNE = 0;
          } else if (jBLK == 0) {
             Roots[iBLK][jBLK].block.nS = 0;
             Roots[iBLK][jBLK].block.nSW = 0;
             Roots[iBLK][jBLK].block.nSE = 0;
             Roots[iBLK][jBLK].block.nN = 1;
             Roots[iBLK][jBLK].block.nNW = 1 * Roots[iBLK][jBLK].block.nNW;
             Roots[iBLK][jBLK].block.nNE = 1 * Roots[iBLK][jBLK].block.nNE;
	    } else if (jBLK == NRj-1) {
             Roots[iBLK][jBLK].block.nS = 1;
             Roots[iBLK][jBLK].block.nSW = 1 * Roots[iBLK][jBLK].block.nSW;
             Roots[iBLK][jBLK].block.nSE = 1 * Roots[iBLK][jBLK].block.nSE;
             Roots[iBLK][jBLK].block.nN = 0;
             Roots[iBLK][jBLK].block.nNW = 0;
             Roots[iBLK][jBLK].block.nNE = 0;
          } else {
             Roots[iBLK][jBLK].block.nS = 1;
             Roots[iBLK][jBLK].block.nSW = 1 * Roots[iBLK][jBLK].block.nSW;
             Roots[iBLK][jBLK].block.nSE = 1 * Roots[iBLK][jBLK].block.nSE;
             Roots[iBLK][jBLK].block.nN = 1;
             Roots[iBLK][jBLK].block.nNW = 1 * Roots[iBLK][jBLK].block.nNW;
             Roots[iBLK][jBLK].block.nNE = 1 * Roots[iBLK][jBLK].block.nNE;
          } /* endif */
	  
          // Check to see if the neighbour numbers are valid.
	  if (Roots[iBLK][jBLK].block.nW == 1 &&
	      !Roots[iBLK-1][jBLK].block.used && 
              Roots[iBLK-1][jBLK].childNW_ptr == NULL &&
              Roots[iBLK-1][jBLK].childNE_ptr == NULL &&
              Roots[iBLK-1][jBLK].childSE_ptr == NULL &&
              Roots[iBLK-1][jBLK].childSW_ptr == NULL) {
             Roots[iBLK][jBLK].block.nW = 0;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nE == 1 &&
              !Roots[iBLK+1][jBLK].block.used && 
              Roots[iBLK+1][jBLK].childNW_ptr == NULL &&
              Roots[iBLK+1][jBLK].childNE_ptr == NULL &&
              Roots[iBLK+1][jBLK].childSE_ptr == NULL &&
              Roots[iBLK+1][jBLK].childSW_ptr == NULL) { 
             Roots[iBLK][jBLK].block.nE = 0;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nS == 1 &&
              !Roots[iBLK][jBLK-1].block.used && 
              Roots[iBLK][jBLK-1].childNW_ptr == NULL &&
              Roots[iBLK][jBLK-1].childNE_ptr == NULL &&
              Roots[iBLK][jBLK-1].childSE_ptr == NULL &&
              Roots[iBLK][jBLK-1].childSW_ptr == NULL) { 
             Roots[iBLK][jBLK].block.nS = 0;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nN == 1 &&
              !Roots[iBLK][jBLK+1].block.used && 
              Roots[iBLK][jBLK+1].childNW_ptr == NULL &&
              Roots[iBLK][jBLK+1].childNE_ptr == NULL &&
              Roots[iBLK][jBLK+1].childSE_ptr == NULL &&
              Roots[iBLK][jBLK+1].childSW_ptr == NULL) { 
             Roots[iBLK][jBLK].block.nN = 0;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nNW == 1 &&
              !Roots[iBLK-1][jBLK+1].block.used && 
              Roots[iBLK-1][jBLK+1].childNW_ptr == NULL &&
              Roots[iBLK-1][jBLK+1].childNE_ptr == NULL &&
              Roots[iBLK-1][jBLK+1].childSE_ptr == NULL &&
              Roots[iBLK-1][jBLK+1].childSW_ptr == NULL) { 
             Roots[iBLK][jBLK].block.nNW = 0;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nSW == 1 &&
              !Roots[iBLK-1][jBLK-1].block.used && 
              Roots[iBLK-1][jBLK-1].childNW_ptr == NULL &&
              Roots[iBLK-1][jBLK-1].childNE_ptr == NULL &&
              Roots[iBLK-1][jBLK-1].childSE_ptr == NULL &&
              Roots[iBLK-1][jBLK-1].childSW_ptr == NULL) { 
             Roots[iBLK][jBLK].block.nSW = 0;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nNE == 1 &&
              !Roots[iBLK+1][jBLK+1].block.used && 
              Roots[iBLK+1][jBLK+1].childNW_ptr == NULL &&
              Roots[iBLK+1][jBLK+1].childNE_ptr == NULL &&
              Roots[iBLK+1][jBLK+1].childSE_ptr == NULL &&
              Roots[iBLK+1][jBLK+1].childSW_ptr == NULL) { 
             Roots[iBLK][jBLK].block.nNE = 0;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nSE == 1 &&
              !Roots[iBLK+1][jBLK-1].block.used && 
              Roots[iBLK+1][jBLK-1].childNW_ptr == NULL &&
              Roots[iBLK+1][jBLK-1].childNE_ptr == NULL &&
              Roots[iBLK+1][jBLK-1].childSE_ptr == NULL &&
              Roots[iBLK+1][jBLK-1].childSW_ptr == NULL) {
             Roots[iBLK][jBLK].block.nSE = 0;
          } /* endif */
	  
          // Finally, assign neighbour information.
	  if (Roots[iBLK][jBLK].block.nW == 1) {
             Roots[iBLK][jBLK].block.infoW[0] = Roots[iBLK-1][jBLK].block.info;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nE == 1) {
             Roots[iBLK][jBLK].block.infoE[0] = Roots[iBLK+1][jBLK].block.info;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nS == 1) {
             Roots[iBLK][jBLK].block.infoS[0] = Roots[iBLK][jBLK-1].block.info;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nN == 1) {
             Roots[iBLK][jBLK].block.infoN[0] = Roots[iBLK][jBLK+1].block.info;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nNW == 1) {
             Roots[iBLK][jBLK].block.infoNW[0] = Roots[iBLK-1][jBLK+1].block.info;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nSW == 1) {
             Roots[iBLK][jBLK].block.infoSW[0] = Roots[iBLK-1][jBLK-1].block.info;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nNE == 1) {
             Roots[iBLK][jBLK].block.infoNE[0] = Roots[iBLK+1][jBLK+1].block.info;
          } /* endif */
	  
	  if (Roots[iBLK][jBLK].block.nSE == 1) {
             Roots[iBLK][jBLK].block.infoSE[0] = Roots[iBLK+1][jBLK-1].block.info;
          } /* endif */
      } /* endfor */
  } /* endfor */
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::findNeighbours -- Find block neighbours.    *
 ****************************************************************************/
inline void QuadTreeBlock_DataStructure::findNeighbours(void) {
  int i_blk, j_blk; 
  QuadTreeBlock *neighbour_block_ptr;
  for ( j_blk = 0 ; j_blk <= Nblk-1 ; ++j_blk ) {
     for ( i_blk = 0 ; i_blk <= Ncpu-1 ; ++i_blk ) {
        if (Blocks[i_blk][j_blk] != NULL) {
	   if (Blocks[i_blk][j_blk]->block.used) {
	      // Find neighbours to the north. 
              neighbour_block_ptr = getNeighbour(Blocks[i_blk][j_blk], 
                                                 QUADTREE_DIRECTION_MASK_NORTH);
              if (neighbour_block_ptr != NULL) {
                 if (neighbour_block_ptr->block.used) {
                    Blocks[i_blk][j_blk]->block.nN = 1;
                    Blocks[i_blk][j_blk]->block.infoN[0] = neighbour_block_ptr->block.info;
                    if (neighbour_block_ptr->block.info.level > 
                        Blocks[i_blk][j_blk]->block.info.level) {
                       Blocks[i_blk][j_blk]->block.nN += 1;
                       Blocks[i_blk][j_blk]->block.infoN[1] =
                          getNeighbour(neighbour_block_ptr, 
                                       QUADTREE_DIRECTION_MASK_EAST)->block.info;
                    } /* endif */
		 } else {
                    Blocks[i_blk][j_blk]->block.nN = 0;
                 } /* endif */
              } else {
                 Blocks[i_blk][j_blk]->block.nN = 0;
              } /* endif */

	      // Find neighbours to the south. 
              neighbour_block_ptr = getNeighbour(Blocks[i_blk][j_blk], 
                                                 QUADTREE_DIRECTION_MASK_SOUTH);
              if (neighbour_block_ptr != NULL) {
                 if (neighbour_block_ptr->block.used) {
                    Blocks[i_blk][j_blk]->block.nS = 1;
                    Blocks[i_blk][j_blk]->block.infoS[0] = neighbour_block_ptr->block.info;
                    if (neighbour_block_ptr->block.info.level > 
                        Blocks[i_blk][j_blk]->block.info.level) {
                       Blocks[i_blk][j_blk]->block.nS += 1;
                       Blocks[i_blk][j_blk]->block.infoS[1] =
                          getNeighbour(neighbour_block_ptr, 
                                       QUADTREE_DIRECTION_MASK_EAST)->block.info;
                    } /* endif */
		 } else {
                    Blocks[i_blk][j_blk]->block.nS = 0;
                 } /* endif */
              } else {
                 Blocks[i_blk][j_blk]->block.nS = 0;
              } /* endif */

	      // Find neighbours to the east. 
              neighbour_block_ptr = getNeighbour(Blocks[i_blk][j_blk], 
                                                 QUADTREE_DIRECTION_MASK_EAST);
              if (neighbour_block_ptr != NULL) {
                 if (neighbour_block_ptr->block.used) {
                    Blocks[i_blk][j_blk]->block.nE = 1;
                    Blocks[i_blk][j_blk]->block.infoE[0] = neighbour_block_ptr->block.info;
                    if (neighbour_block_ptr->block.info.level > 
                        Blocks[i_blk][j_blk]->block.info.level) {
                       Blocks[i_blk][j_blk]->block.nE += 1;
                       Blocks[i_blk][j_blk]->block.infoE[1] =
                          getNeighbour(neighbour_block_ptr, 
                                       QUADTREE_DIRECTION_MASK_NORTH)->block.info;
                    } /* endif */
		 } else {
                    Blocks[i_blk][j_blk]->block.nE = 0;
                 } /* endif */
              } else {
                 Blocks[i_blk][j_blk]->block.nE = 0;
              } /* endif */

	      // Find neighbours to the west. 
              neighbour_block_ptr = getNeighbour(Blocks[i_blk][j_blk], 
                                                 QUADTREE_DIRECTION_MASK_WEST);
              if (neighbour_block_ptr != NULL) {
                 if (neighbour_block_ptr->block.used) {
                    Blocks[i_blk][j_blk]->block.nW = 1;
                    Blocks[i_blk][j_blk]->block.infoW[0] = neighbour_block_ptr->block.info;
                    if (neighbour_block_ptr->block.info.level > 
                        Blocks[i_blk][j_blk]->block.info.level) {
                       Blocks[i_blk][j_blk]->block.nW += 1;
                       Blocks[i_blk][j_blk]->block.infoW[1] =
                          getNeighbour(neighbour_block_ptr, 
                                       QUADTREE_DIRECTION_MASK_NORTH)->block.info;
                    } /* endif */
		 } else {
                    Blocks[i_blk][j_blk]->block.nW = 0;
                 } /* endif */
              } else {
                 Blocks[i_blk][j_blk]->block.nW = 0;
              } /* endif */

	      // Find neighbours to the north west. 
              neighbour_block_ptr = getNeighbour(Blocks[i_blk][j_blk], 
                                                 QUADTREE_DIRECTION_MASK_NORTHWEST);
              if (neighbour_block_ptr != NULL) {
                 if (neighbour_block_ptr->block.used) {
                    if (neighbour_block_ptr->block.info.level >= 
                        Blocks[i_blk][j_blk]->block.info.level) {
                       Blocks[i_blk][j_blk]->block.nNW = 1;
                       Blocks[i_blk][j_blk]->block.infoNW[0] = neighbour_block_ptr->block.info;
                    } else if ((neighbour_block_ptr != 
                                getNeighbour(Blocks[i_blk][j_blk]->parent_ptr, 
                                             QUADTREE_DIRECTION_MASK_NORTH)) &&
                               (neighbour_block_ptr != 
                                getNeighbour(Blocks[i_blk][j_blk]->parent_ptr, 
                                             QUADTREE_DIRECTION_MASK_WEST))) {
                       Blocks[i_blk][j_blk]->block.nNW = 1;
                       Blocks[i_blk][j_blk]->block.infoNW[0] = neighbour_block_ptr->block.info;
                    } else {
                       Blocks[i_blk][j_blk]->block.nNW = 0;
                    } /* endif */
		 } else {
                    Blocks[i_blk][j_blk]->block.nNW = 0;
                 } /* endif */
              } else {
                 Blocks[i_blk][j_blk]->block.nNW = 0;
              } /* endif */

	      // Find neighbours to the north east. 
              neighbour_block_ptr = getNeighbour(Blocks[i_blk][j_blk], 
                                                 QUADTREE_DIRECTION_MASK_NORTHEAST);
              if (neighbour_block_ptr != NULL) {
                 if (neighbour_block_ptr->block.used) {
                    if (neighbour_block_ptr->block.info.level >= 
                        Blocks[i_blk][j_blk]->block.info.level) {
                       Blocks[i_blk][j_blk]->block.nNE = 1;
                       Blocks[i_blk][j_blk]->block.infoNE[0] = neighbour_block_ptr->block.info;
                    } else if ((neighbour_block_ptr != 
                                getNeighbour(Blocks[i_blk][j_blk]->parent_ptr, 
                                             QUADTREE_DIRECTION_MASK_NORTH)) &&
                               (neighbour_block_ptr != 
                                getNeighbour(Blocks[i_blk][j_blk]->parent_ptr, 
                                             QUADTREE_DIRECTION_MASK_EAST))) {
                       Blocks[i_blk][j_blk]->block.nNE = 1;
                       Blocks[i_blk][j_blk]->block.infoNE[0] = neighbour_block_ptr->block.info;
                    } else {
                       Blocks[i_blk][j_blk]->block.nNE = 0;
                    } /* endif */
		 } else {
                    Blocks[i_blk][j_blk]->block.nNE = 0;
                 } /* endif */
              } else {
                 Blocks[i_blk][j_blk]->block.nNE = 0;
              } /* endif */

	      // Find neighbours to the south east. 
              neighbour_block_ptr = getNeighbour(Blocks[i_blk][j_blk], 
                                                 QUADTREE_DIRECTION_MASK_SOUTHEAST);
              if (neighbour_block_ptr != NULL) {
                 if (neighbour_block_ptr->block.used) {
                    if (neighbour_block_ptr->block.info.level >= 
                        Blocks[i_blk][j_blk]->block.info.level) {
                       Blocks[i_blk][j_blk]->block.nSE = 1;
                       Blocks[i_blk][j_blk]->block.infoSE[0] = neighbour_block_ptr->block.info;
                    } else if ((neighbour_block_ptr != 
                                getNeighbour(Blocks[i_blk][j_blk]->parent_ptr, 
                                             QUADTREE_DIRECTION_MASK_SOUTH)) &&
                               (neighbour_block_ptr != 
                                getNeighbour(Blocks[i_blk][j_blk]->parent_ptr, 
                                             QUADTREE_DIRECTION_MASK_EAST))) {
                       Blocks[i_blk][j_blk]->block.nSE = 1;
                       Blocks[i_blk][j_blk]->block.infoSE[0] = neighbour_block_ptr->block.info;
                    } else {
                       Blocks[i_blk][j_blk]->block.nSE = 0;
                    } /* endif */
		 } else {
                    Blocks[i_blk][j_blk]->block.nSE = 0;
                 } /* endif */
              } else {
                 Blocks[i_blk][j_blk]->block.nSE = 0;
              } /* endif */

	      // Find neighbours to the south west. 
              neighbour_block_ptr = getNeighbour(Blocks[i_blk][j_blk], 
                                                 QUADTREE_DIRECTION_MASK_SOUTHWEST);
              if (neighbour_block_ptr != NULL) {
                 if (neighbour_block_ptr->block.used) {
                    if (neighbour_block_ptr->block.info.level >= 
                        Blocks[i_blk][j_blk]->block.info.level) {
                       Blocks[i_blk][j_blk]->block.nSW = 1;
                       Blocks[i_blk][j_blk]->block.infoSW[0] = neighbour_block_ptr->block.info;
                    } else if ((neighbour_block_ptr != 
                                getNeighbour(Blocks[i_blk][j_blk]->parent_ptr, 
                                             QUADTREE_DIRECTION_MASK_SOUTH)) &&
                               (neighbour_block_ptr != 
                                getNeighbour(Blocks[i_blk][j_blk]->parent_ptr, 
                                             QUADTREE_DIRECTION_MASK_WEST))) {
                       Blocks[i_blk][j_blk]->block.nSW = 1;
                       Blocks[i_blk][j_blk]->block.infoSW[0] = neighbour_block_ptr->block.info;
                    } else {
                       Blocks[i_blk][j_blk]->block.nSW = 0;
                    } /* endif */
		 } else {
                    Blocks[i_blk][j_blk]->block.nSW = 0;
                 } /* endif */
              } else {
                 Blocks[i_blk][j_blk]->block.nSW = 0;
              } /* endif */
           } /* endif */
        } /* endif */
     } /* endfor */
  } /* endfor */
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::refineBlock -- Refine (divide) block.       *
 ****************************************************************************/
inline void QuadTreeBlock_DataStructure::refineBlock(int *new_blocks_CPU, 
                                                     int *new_blocks_BLK,
                                                     int *new_blocks_SECTOR) {
   int iNEW;
   QuadTreeBlock *block_to_be_refined_ptr;

   // Get pointer to block to be refined.
   block_to_be_refined_ptr = Blocks[new_blocks_CPU[0]][new_blocks_BLK[0]];

   // Set block not used.
   block_to_be_refined_ptr->block.used = ADAPTIVEBLOCK2D_NOT_USED;

   for ( iNEW = 0 ; iNEW <= 3 ; ++iNEW ) {
      switch(new_blocks_SECTOR[iNEW]) {
        case QUADTREE_SECTOR_NW :
          // Create NW offspring.
          if (block_to_be_refined_ptr->childNW_ptr != NULL) {
             delete block_to_be_refined_ptr->childNW_ptr;
             block_to_be_refined_ptr->childNW_ptr = NULL;
          } /* endif */
          block_to_be_refined_ptr->childNW_ptr = new QuadTreeBlock;
          block_to_be_refined_ptr->childNW_ptr->block.used = ADAPTIVEBLOCK2D_USED;
          block_to_be_refined_ptr->childNW_ptr->block.gblknum = 
             block_to_be_refined_ptr->block.gblknum;
          block_to_be_refined_ptr->childNW_ptr->block.info.cpu = new_blocks_CPU[iNEW];
          block_to_be_refined_ptr->childNW_ptr->block.info.blknum = new_blocks_BLK[iNEW];
          block_to_be_refined_ptr->childNW_ptr->block.info.dimen.i = 
             block_to_be_refined_ptr->block.info.dimen.i;
          block_to_be_refined_ptr->childNW_ptr->block.info.dimen.j = 
             block_to_be_refined_ptr->block.info.dimen.j;
          block_to_be_refined_ptr->childNW_ptr->block.info.dimen.ghost = 
             block_to_be_refined_ptr->block.info.dimen.ghost;
          block_to_be_refined_ptr->childNW_ptr->block.info.sector = new_blocks_SECTOR[iNEW];
          block_to_be_refined_ptr->childNW_ptr->block.info.level = 
             block_to_be_refined_ptr->block.info.level + 1;
          block_to_be_refined_ptr->childNW_ptr->parent_ptr = block_to_be_refined_ptr;
          block_to_be_refined_ptr->childNW_ptr->childNW_ptr = NULL;
          block_to_be_refined_ptr->childNW_ptr->childNE_ptr = NULL;
          block_to_be_refined_ptr->childNW_ptr->childSE_ptr = NULL;
          block_to_be_refined_ptr->childNW_ptr->childSW_ptr = NULL;
          Blocks[new_blocks_CPU[iNEW]][new_blocks_BLK[iNEW]] = 
             block_to_be_refined_ptr->childNW_ptr;
          break;
        case QUADTREE_SECTOR_NE :
          // Create NE offspring.
          if (block_to_be_refined_ptr->childNE_ptr != NULL) {
             delete block_to_be_refined_ptr->childNE_ptr;
             block_to_be_refined_ptr->childNE_ptr = NULL;
          } /* endif */
          block_to_be_refined_ptr->childNE_ptr = new QuadTreeBlock;
          block_to_be_refined_ptr->childNE_ptr->block.used = ADAPTIVEBLOCK2D_USED;
          block_to_be_refined_ptr->childNE_ptr->block.gblknum = 
             block_to_be_refined_ptr->block.gblknum;
          block_to_be_refined_ptr->childNE_ptr->block.info.cpu = new_blocks_CPU[iNEW];
          block_to_be_refined_ptr->childNE_ptr->block.info.blknum = new_blocks_BLK[iNEW];
          block_to_be_refined_ptr->childNE_ptr->block.info.dimen.i = 
             block_to_be_refined_ptr->block.info.dimen.i;
          block_to_be_refined_ptr->childNE_ptr->block.info.dimen.j = 
             block_to_be_refined_ptr->block.info.dimen.j;
          block_to_be_refined_ptr->childNE_ptr->block.info.dimen.ghost = 
             block_to_be_refined_ptr->block.info.dimen.ghost;
          block_to_be_refined_ptr->childNE_ptr->block.info.sector = new_blocks_SECTOR[iNEW];
          block_to_be_refined_ptr->childNE_ptr->block.info.level = 
             block_to_be_refined_ptr->block.info.level + 1;
          block_to_be_refined_ptr->childNE_ptr->parent_ptr = block_to_be_refined_ptr;
          block_to_be_refined_ptr->childNE_ptr->childNW_ptr = NULL;
          block_to_be_refined_ptr->childNE_ptr->childNE_ptr = NULL;
          block_to_be_refined_ptr->childNE_ptr->childSE_ptr = NULL;
          block_to_be_refined_ptr->childNE_ptr->childSW_ptr = NULL;
          Blocks[new_blocks_CPU[iNEW]][new_blocks_BLK[iNEW]] = 
             block_to_be_refined_ptr->childNE_ptr;
          break;
        case QUADTREE_SECTOR_SE :
          // Create SE offspring.
          if (block_to_be_refined_ptr->childSE_ptr != NULL) {
             delete block_to_be_refined_ptr->childSE_ptr;
             block_to_be_refined_ptr->childSE_ptr = NULL;
          } /* endif */
          block_to_be_refined_ptr->childSE_ptr = new QuadTreeBlock;
          block_to_be_refined_ptr->childSE_ptr->block.used = ADAPTIVEBLOCK2D_USED;
          block_to_be_refined_ptr->childSE_ptr->block.gblknum = 
             block_to_be_refined_ptr->block.gblknum;
          block_to_be_refined_ptr->childSE_ptr->block.info.cpu = new_blocks_CPU[iNEW];
          block_to_be_refined_ptr->childSE_ptr->block.info.blknum = new_blocks_BLK[iNEW];
          block_to_be_refined_ptr->childSE_ptr->block.info.dimen.i = 
             block_to_be_refined_ptr->block.info.dimen.i;
          block_to_be_refined_ptr->childSE_ptr->block.info.dimen.j = 
             block_to_be_refined_ptr->block.info.dimen.j;
          block_to_be_refined_ptr->childSE_ptr->block.info.dimen.ghost = 
             block_to_be_refined_ptr->block.info.dimen.ghost;
          block_to_be_refined_ptr->childSE_ptr->block.info.sector = new_blocks_SECTOR[iNEW];
          block_to_be_refined_ptr->childSE_ptr->block.info.level = 
             block_to_be_refined_ptr->block.info.level + 1;
          block_to_be_refined_ptr->childSE_ptr->parent_ptr = block_to_be_refined_ptr;
          block_to_be_refined_ptr->childSE_ptr->childNW_ptr = NULL;
          block_to_be_refined_ptr->childSE_ptr->childNE_ptr = NULL;
          block_to_be_refined_ptr->childSE_ptr->childSE_ptr = NULL;
          block_to_be_refined_ptr->childSE_ptr->childSW_ptr = NULL;
          Blocks[new_blocks_CPU[iNEW]][new_blocks_BLK[iNEW]] = 
             block_to_be_refined_ptr->childSE_ptr;
          break;
        case QUADTREE_SECTOR_SW :
        default:
          // Create SW offspring.
          if (block_to_be_refined_ptr->childSW_ptr != NULL) {
             delete block_to_be_refined_ptr->childSW_ptr;
             block_to_be_refined_ptr->childSW_ptr = NULL;
          } /* endif */
          block_to_be_refined_ptr->childSW_ptr = new QuadTreeBlock;
          block_to_be_refined_ptr->childSW_ptr->block.used = ADAPTIVEBLOCK2D_USED;
          block_to_be_refined_ptr->childSW_ptr->block.gblknum = 
             block_to_be_refined_ptr->block.gblknum;
          block_to_be_refined_ptr->childSW_ptr->block.info.cpu = new_blocks_CPU[iNEW];
          block_to_be_refined_ptr->childSW_ptr->block.info.blknum = new_blocks_BLK[iNEW];
          block_to_be_refined_ptr->childSW_ptr->block.info.dimen.i = 
             block_to_be_refined_ptr->block.info.dimen.i;
          block_to_be_refined_ptr->childSW_ptr->block.info.dimen.j = 
             block_to_be_refined_ptr->block.info.dimen.j;
          block_to_be_refined_ptr->childSW_ptr->block.info.dimen.ghost = 
             block_to_be_refined_ptr->block.info.dimen.ghost;
          block_to_be_refined_ptr->childSW_ptr->block.info.sector = new_blocks_SECTOR[iNEW];
          block_to_be_refined_ptr->childSW_ptr->block.info.level = 
             block_to_be_refined_ptr->block.info.level + 1;
          block_to_be_refined_ptr->childSW_ptr->parent_ptr = block_to_be_refined_ptr;
          block_to_be_refined_ptr->childSW_ptr->childNW_ptr = NULL;
          block_to_be_refined_ptr->childSW_ptr->childNE_ptr = NULL;
          block_to_be_refined_ptr->childSW_ptr->childSE_ptr = NULL;
          block_to_be_refined_ptr->childSW_ptr->childSW_ptr = NULL;
          Blocks[new_blocks_CPU[iNEW]][new_blocks_BLK[iNEW]] = 
             block_to_be_refined_ptr->childSW_ptr;
          break;
      } /* endswitch */
   } /* endfor */

}

/****************************************************************************
 * QuadTreeBlock_DataStructure::coarsenBlock -- Coarsen (contract) blocks.  *
 ****************************************************************************/
inline void QuadTreeBlock_DataStructure::coarsenBlocks(int *old_blocks_CPU, 
                                                       int *old_blocks_BLK,
                                                       int *old_blocks_SECTOR) {

   int iOLD;
   QuadTreeBlock *coarsened_block_ptr;

   // Get pointer to block which will result from coarsening process.
   coarsened_block_ptr = Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->parent_ptr;

   // Assign block information to coarsened block.
   coarsened_block_ptr->block.used = ADAPTIVEBLOCK2D_USED;
   coarsened_block_ptr->block.gblknum = 
      Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.gblknum;
   coarsened_block_ptr->block.info.cpu = old_blocks_CPU[0];
   coarsened_block_ptr->block.info.blknum = old_blocks_BLK[0];
   coarsened_block_ptr->block.info.dimen.i = 
      Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.dimen.i;
   coarsened_block_ptr->block.info.dimen.j = 
      Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.dimen.j;
   coarsened_block_ptr->block.info.dimen.ghost = 
      Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.dimen.ghost;
   //coarsened_block_ptr->info.sector = coarsened_block_ptr->info.sector;
   coarsened_block_ptr->block.info.level = 
      Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.level - 1;
   //coarsened_block_ptr->parent_ptr = coarsened_block_ptr->parent_ptr;
   coarsened_block_ptr->childNW_ptr = NULL;
   coarsened_block_ptr->childNE_ptr = NULL;
   coarsened_block_ptr->childSE_ptr = NULL;
   coarsened_block_ptr->childSW_ptr = NULL;
   delete Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]];
   Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]] = coarsened_block_ptr;

   // De-allocate old refined blocks.
   for ( iOLD = 1 ; iOLD <= 3 ; ++iOLD ) {
      if (Blocks[old_blocks_CPU[iOLD]][old_blocks_BLK[iOLD]] != NULL) {
         delete Blocks[old_blocks_CPU[iOLD]][old_blocks_BLK[iOLD]];
         Blocks[old_blocks_CPU[iOLD]][old_blocks_BLK[iOLD]] = NULL;
      } /* endif */
   } /* endfor */
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::nochangeAll -- Set no refinement flags.     *
 ****************************************************************************/
inline void QuadTreeBlock_DataStructure::nochangeAll(void) {
   int i, j; 
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      for ( j = 0; j <= Nblk-1 ; ++j ) {
         RefineFlags[i][j] = ADAPTIVEBLOCK2D_NOCHANGE;
      } /* endfor */
   } /* endfor */
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::refineAll -- Set refinement flags.          *
 ****************************************************************************/
inline void QuadTreeBlock_DataStructure::refineAll(void) {
   int i, j; 
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      for ( j = 0; j <= Nblk-1 ; ++j ) {
         if (Blocks[i][j] != NULL) {
            if (Blocks[i][j]->block.used) {
               RefineFlags[i][j] = ADAPTIVEBLOCK2D_REFINE;
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::coarsenAll -- Set coarsening flags.         *
 ****************************************************************************/
inline void QuadTreeBlock_DataStructure::coarsenAll(void) {
   int i, j; 
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      for ( j = 0; j <= Nblk-1 ; ++j ) {
         if (Blocks[i][j] != NULL) {
            if (Blocks[i][j]->block.used) {
               RefineFlags[i][j] = ADAPTIVEBLOCK2D_COARSEN;
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::setRefineAll -- Set refinement flags.       *
 ****************************************************************************/
inline void QuadTreeBlock_DataStructure::setRefineAll(const int Flag) {
   int i, j; 
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      for ( j = 0; j <= Nblk-1 ; ++j ) {
         if (Blocks[i][j] != NULL) {
            if (Blocks[i][j]->block.used) {
               RefineFlags[i][j] = Flag;
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::highestRefinementLevel -- Return highest    *
 *                      refinement level of all solution blocks.            *
 ****************************************************************************/
inline int QuadTreeBlock_DataStructure::highestRefinementLevel(void) {
   int i, j, level = 0;
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      for ( j = 0; j <= Nblk-1 ; ++j ) {
         if (Blocks[i][j] != NULL) {
            if (Blocks[i][j]->block.used) {
               level = max(level, Blocks[i][j]->block.info.level);
            } /* endif */
         } /* endif */             
      } /* endfor */
   } /* endfor */
   return (level);
}

inline int QuadTreeBlock_DataStructure::highestRefinementLevel(const int use_tree) {
   int i, j, level = 0;
   if (use_tree) {
      for ( j = 0 ; j <= NRj-1 ; ++j ) {
         for ( i = 0 ; i <= NRi-1 ; ++i ) {
	    level = Roots[i][j].maxRefinementLevel(level);
         } /* endfor */
      } /* endfor */
   } else {
      for ( i = 0; i <= Ncpu-1 ; ++i ) {
         for ( j = 0; j <= Nblk-1 ; ++j ) {
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
 * QuadTreeBlock_DataStructure::numberToBeRefined -- Return number of       *
 *                      quadtree solution blocks to be refined.             *
 ****************************************************************************/
inline int QuadTreeBlock_DataStructure::numberToBeRefined(void) {
   int i, j, number_to_be_refined = 0;
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      for ( j = 0; j <= Nblk-1 ; ++j ) {
         if (Blocks[i][j] != NULL) {
   	    if (Blocks[i][j]->block.used) {
               if (RefineFlags[i][j] == ADAPTIVEBLOCK2D_REFINE) {
                  number_to_be_refined += 1;
               } /* endif */
            } /* endif */
         } /* endif */             
      } /* endfor */
   } /* endfor */
   return (number_to_be_refined);
}

/****************************************************************************
 * QuadTreeBlock_DataStructure::numberToBeCoarsened -- Return number of     *
 *                      quadtree solution blocks to be coarsened.           *
 ****************************************************************************/
inline int QuadTreeBlock_DataStructure::numberToBeCoarsened(void) {
   int i, j, number_to_be_coarsened = 0;
   for ( i = 0; i <= Ncpu-1 ; ++i ) {
      for ( j = 0; j <= Nblk-1 ; ++j ) {
         if (Blocks[i][j] != NULL) {
   	    if (Blocks[i][j]->block.used) {
               if (RefineFlags[i][j] == ADAPTIVEBLOCK2D_COARSEN) {
                  number_to_be_coarsened += 1;
               } /* endif */
            } /* endif */
         } /* endif */             
      } /* endfor */
   } /* endfor */
   return (number_to_be_coarsened);
}

/****************************************************************************
 * QuadTreeBlock_DataStructure -- Input-output operators.                   *
 ****************************************************************************/
inline ostream &operator << (ostream &out_file,
			     const QuadTreeBlock_DataStructure &QT) {
   int iBLK, jBLK;
   out_file << QT.NRi << " " << QT.NRj << " " << QT.Ncpu << " " << QT.Nblk << "\n";
   for ( jBLK = 0 ; jBLK <= QT.NRj-1 ; ++jBLK ) {
       for ( iBLK = 0 ; iBLK <= QT.NRi-1 ; ++iBLK ) {
	  QT.Roots[iBLK][jBLK].write(out_file);
       } /* endfor */
   } /* endfor */
   return (out_file);
}

inline istream &operator >> (istream &in_file,
			     QuadTreeBlock_DataStructure &QT) {
   int iBLK, jBLK, nri, nrj, ncpu, nblk;
   if (QT.Roots != NULL && QT.Blocks != NULL) {
     QT.deallocate();
   } else if (QT.Roots != NULL) {
     QT.deallocateRoots();
   } else if (QT.Blocks != NULL) {
      QT.deallocateBlocks();
   } /* endif */
   in_file.setf(ios::skipws);
   in_file >> nri >> nrj >> ncpu >> nblk;
   in_file.setf(ios::skipws);
   QT.allocate(nri, nrj, ncpu, nblk);
   for ( jBLK = 0 ; jBLK <= QT.NRj-1 ; ++jBLK ) {
       for ( iBLK = 0 ; iBLK <= QT.NRi-1 ; ++iBLK ) {
	  QT.Roots[iBLK][jBLK].read(in_file);
       } /* endfor */
   } /* endfor */
   return (in_file);
}

/*********************************************************************
 * QuadTreeBlock_Data_Structure -- External subroutines.             *
 *********************************************************************/

extern void Create_QuadTree_Data_Structure(QuadTreeBlock_DataStructure &QuadTree,
		                           const int Number_of_Roots_Idir,
		                           const int Number_of_Roots_Jdir,
                                           const int Number_of_Processors,
                                           const int Number_of_Blocks_per_Processor);

extern void Broadcast_QuadTree_Data_Structure(QuadTreeBlock_DataStructure &QuadTree,
                                              AdaptiveBlockResourceList   &List_of_Available_Blocks);

extern void Renumber_Solution_Blocks(QuadTreeBlock_DataStructure &QuadTree);

extern void Renumber_Solution_Blocks(QuadTreeBlock_DataStructure &QuadTree,
                                     AdaptiveBlock2D_List &LocalSolnBlockList);

extern void Find_Neighbours_of_Root_Solution_Blocks(QuadTreeBlock_DataStructure &QuadTree);

extern void Find_Neighbours_of_Root_Solution_Blocks(QuadTreeBlock_DataStructure &QuadTree,
                                                    AdaptiveBlock2D_List &LocalSolnBlockList);

extern void Modify_Neighbours_of_Root_Solution_Blocks(QuadTreeBlock_DataStructure &QuadTree,
		                                      const int Grid_Type);

extern void Modify_Neighbours_of_Root_Solution_Blocks(QuadTreeBlock_DataStructure &QuadTree,
                                                      AdaptiveBlock2D_List &LocalSolnBlockList,
		                                      const int Grid_Type);

extern void Find_Neighbours(QuadTreeBlock_DataStructure &QuadTree);

extern void Find_Neighbours(QuadTreeBlock_DataStructure &QuadTree,
                            AdaptiveBlock2D_List &LocalSolnBlockList);

extern void Get_Refinement_List(QuadTreeBlock_DataStructure &QuadTree,
                                AdaptiveBlock2D_List &LocalSolnBlockList);

/* Include Morton re-ordering header file. */

#ifndef _QUADTREE_MORTON_INCLUDED
#include "QuadTree_Morton.h"
#endif // _QUADTREE_MORTON_INCLUDED

#endif /* _QUADTREE_INCLUDED  */
