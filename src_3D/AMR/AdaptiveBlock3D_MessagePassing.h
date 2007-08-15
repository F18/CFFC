/* AdaptiveBlock3D_MessagePassing.h:  Templated subroutines for adaptive blocks classes
                                      which load and unload message passing buffers. */

#ifndef _ADAPTIVEBLOCK3D_MESSAGEPASSING_INCLUDED
#define _ADAPTIVEBLOCK3D_MESSAGEPASSING_INCLUDED

/* Include adaptive block and 3D vector header files. */

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "../Grid/Grid3DHexaBlock.h"
#endif  //_GRID3D_HEXA_BLOCK_INCLUDE

#ifndef _ADAPTIVBLOCK_INCLUDED
#include "AdaptiveBlock3D.h"
#endif // _ADAPTIVEBLOCK_INCLUDED

/*************************************************************
 * AdaptiveBlock3D -- Templated message passing subroutines. *
 *************************************************************/

/********************************************************
 * Routine: Load_Send_Message_Buffers_NoResChange       *
 *                                                      *
 * Loads the send message passing buffers for sending   *
 * solution information to neighbouring 3D              *
 * hexahedrial multi-block solution blocks with no      *
 * mesh resolution changes (works on entire 1D array of *
 * 3D solution blocks).                                 *
 *                                                      *
 ********************************************************/
template <class Hexa_Soln_Block>
int Load_Send_Message_Buffers_NoResChange(Hexa_Soln_Block *Soln_Blks,
                                          AdaptiveBlock3D_List &Soln_Block_List,
                                          const int Number_of_Solution_Variables,
                                          const int Send_Mesh_Geometry_Only) {

  
 int i_blk, buffer_size_neighbour, 
    i_min, i_max, i_inc, i, 
    j_min, j_max, j_inc, j, 
    k_min, k_max, k_inc, k, 
    m, l;
  int i_ref, j_ref, k_ref;
  Vector3D x_ref;

  /* Check to see if the number of solution variables specified
     in the calling argument is correct. */

  if (Number_of_Solution_Variables != NUM_COMP_VECTOR3D &&
      Number_of_Solution_Variables != Soln_Blks[0].NumVar() &&
      Number_of_Solution_Variables != Soln_Blks[0].NumVar() + NUM_COMP_VECTOR3D) {
    return(1);
  } /* endif */
   
    /* Load the send buffers of each solution block. */

  for ( i_blk = 0 ; i_blk < Soln_Block_List.Nblk ; ++i_blk ) {

     ////////////////////////////////////////////////////////////////////
    // Load send buffer for Top face neighbours of solution blocks. //
    ////////////////////////////////////////////////////////////////////
    if (Soln_Block_List.Block[i_blk].used &&
	(Soln_Block_List.Block[i_blk].nT == 1) &&
	(Soln_Block_List.Block[i_blk].info.level ==
	 Soln_Block_List.Block[i_blk].infoT[0].level)) {
       
      if (!Send_Mesh_Geometry_Only) {
	buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost*
	  abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.i)*
	  abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.j)*
	  Soln_Blks[i_blk].NumVar();

       
	if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
	  buffer_size_neighbour = buffer_size_neighbour +
	    Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost*
	    (abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.i)+1)*
	    (abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.j)+1)*
	    NUM_COMP_VECTOR3D+
	    4*Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.i)+
	    4*Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.j);
      
      } else {
	buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost*
	  (abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.i)+1)*
	  (abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.j)+1)*
	  NUM_COMP_VECTOR3D+
	  4*Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.i)+
	  4*Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.j);
      } /* endif */
      l = -1;
      // Load ghost cell solution variable information as required.
      if (!Send_Mesh_Geometry_Only) {
	if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0 &&
	    Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0 &&
	    Soln_Block_List.Block[i_blk].infoT[0].dimen.k > 0) {
	  i_min = Soln_Blks[i_blk].ICl;
	  i_max = Soln_Blks[i_blk].ICu;
	  i_inc = 1;
	  j_min = Soln_Blks[i_blk].JCl;
	  j_max = Soln_Blks[i_blk].JCu;
	  j_inc = 1;
	  k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
	  k_max = Soln_Blks[i_blk].KCu;
	  k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0 &&
		   Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0) {
	  i_min = Soln_Blks[i_blk].ICl;
	  i_max = Soln_Blks[i_blk].ICu;
	  i_inc = 1;
	  j_min = Soln_Blks[i_blk].JCl;
	  j_max = Soln_Blks[i_blk].JCu;
	  j_inc = 1;
	  k_min = Soln_Blks[i_blk].KCu;
	  k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
	  k_inc = -1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0 &&
		   Soln_Block_List.Block[i_blk].infoT[0].dimen.k > 0) {
	  i_min = Soln_Blks[i_blk].ICu;
	  i_max = Soln_Blks[i_blk].ICl;
	  i_inc = -1;
	  j_min = Soln_Blks[i_blk].JCl;
	  j_max = Soln_Blks[i_blk].JCu;
	  j_inc = 1;
	  k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
	  k_max = Soln_Blks[i_blk].KCu;
	  k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0 &&
		   Soln_Block_List.Block[i_blk].infoT[0].dimen.k > 0) {
	  i_min = Soln_Blks[i_blk].ICl;
	  i_max = Soln_Blks[i_blk].ICu;
	  i_inc = 1;
	  j_min = Soln_Blks[i_blk].JCu;
	  j_max = Soln_Blks[i_blk].JCl;
	  j_inc = -1;
	  k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
	  k_max = Soln_Blks[i_blk].KCu;
	  k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.k > 0) {
	  i_min = Soln_Blks[i_blk].ICu;
	  i_max = Soln_Blks[i_blk].ICl;
	  i_inc = -1;
	  j_min = Soln_Blks[i_blk].JCu;
	  j_max = Soln_Blks[i_blk].JCl;
	  j_inc = -1;
	  k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
	  k_max = Soln_Blks[i_blk].KCu;
	  k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0) {
	  i_min = Soln_Blks[i_blk].ICu;
	  i_max = Soln_Blks[i_blk].ICl;
	  i_inc = -1;
	  j_min = Soln_Blks[i_blk].JCl;
	  j_max = Soln_Blks[i_blk].JCu;
	  j_inc = 1;
	  k_min = Soln_Blks[i_blk].KCu;
	  k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
	  k_inc = -1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0) {
	  i_min = Soln_Blks[i_blk].ICl;
	  i_max = Soln_Blks[i_blk].ICu;
	  i_inc = 1;
	  j_min = Soln_Blks[i_blk].JCu;
	  j_max = Soln_Blks[i_blk].JCl;
	  j_inc = -1;
	  k_min = Soln_Blks[i_blk].KCu;
	  k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
	  k_inc = -1;
	} else {
	  i_min = Soln_Blks[i_blk].ICu;
	  i_max = Soln_Blks[i_blk].ICl;
	  i_inc = -1;
	  j_min = Soln_Blks[i_blk].JCu;
	  j_max = Soln_Blks[i_blk].JCl;
	  j_inc = -1;
	  k_min = Soln_Blks[i_blk].KCu;
	  k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
	  k_inc = -1;
	} /* endif */


  
        
	i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_topface_sendbuf[i_blk],
					   l,buffer_size_neighbour,
					   i_min,i_max,i_inc,
					   j_min,j_max,j_inc,
					   k_min,k_max,k_inc);
  
	if (i != 0) return(2200);
      } /* endif */
      // Load ghost cell mesh information as required.
      if (Send_Mesh_Geometry_Only || Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0 &&
	    Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0 &&
	    Soln_Block_List.Block[i_blk].infoT[0].dimen.k > 0) {
	  i_min = Soln_Blks[i_blk].Grid.INl;
	  i_max = Soln_Blks[i_blk].Grid.INu;
	  i_inc = 1;
	  j_min = Soln_Blks[i_blk].Grid.JNl;
	  j_max = Soln_Blks[i_blk].Grid.JNu;
	  j_inc = 1;
	  k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost;
	  k_max = Soln_Blks[i_blk].Grid.KNu-1;
	  k_inc = 1;
	  i_ref = i_min;
	  j_ref = j_min;
	  k_ref = k_max+1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0 &&
	      Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0 ) {
	    i_min = Soln_Blks[i_blk].Grid.INl;
	    i_max = Soln_Blks[i_blk].Grid.INu;
	    i_inc = 1;
	    j_min = Soln_Blks[i_blk].Grid.JNl;
	    j_max = Soln_Blks[i_blk].Grid.JNu;
	    j_inc = 1;
	    k_min = Soln_Blks[i_blk].Grid.KNu-1;
	    k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost;
	    k_inc = -1;
	    i_ref = i_min;
	    j_ref = j_min;
	    k_ref = k_min+1;
	} else if ( Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoT[0].dimen.k > 0) {
	      i_min = Soln_Blks[i_blk].Grid.INu;
	      i_max = Soln_Blks[i_blk].Grid.INl;
	      i_inc = -1;
	      j_min = Soln_Blks[i_blk].Grid.JNl;
	      j_max = Soln_Blks[i_blk].Grid.JNu;
	      j_inc = 1;
	      k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost;
	      k_max = Soln_Blks[i_blk].Grid.KNu-1;
	      k_inc = 1;
	      i_ref = i_min;
	      j_ref = j_min;
	      k_ref = k_max+1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0 &&
		  Soln_Block_List.Block[i_blk].infoT[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_max+1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_max+1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_min+1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_min+1;
	} else {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_min+1;
	} /* endif */

	x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
	for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
	      for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(2201);
                      Soln_Block_List.message_noreschange_topface_sendbuf[i_blk][l] =
			Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
	      } /* endfor */
	if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0 &&
                  Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0 &&
		  Soln_Block_List.Block[i_blk].infoT[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0 &&
		    Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0) {
		  i_min = Soln_Blks[i_blk].ICl;
		  i_max = Soln_Blks[i_blk].ICu;
		  i_inc = 1;
		  j_min = Soln_Blks[i_blk].JCl;
		  j_max = Soln_Blks[i_blk].JCu;
		  j_inc = 1;
		  k_min = Soln_Blks[i_blk].KCu;
		  k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
		  k_inc = -1;
	} else if ( Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0 &&
		       Soln_Block_List.Block[i_blk].infoT[0].dimen.k > 0) {
		    i_min = Soln_Blks[i_blk].ICu;
		    i_max = Soln_Blks[i_blk].ICl;
		    i_inc = -1;
		    j_min = Soln_Blks[i_blk].JCl;
		    j_max = Soln_Blks[i_blk].JCu;
		    j_inc = 1;
		    k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
		    k_max = Soln_Blks[i_blk].KCu;
		    k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0 &&
			Soln_Block_List.Block[i_blk].infoT[0].dimen.k > 0) {
		      i_min = Soln_Blks[i_blk].ICl;
		      i_max = Soln_Blks[i_blk].ICu;
		      i_inc = 1;
		      j_min = Soln_Blks[i_blk].JCu;
		      j_max = Soln_Blks[i_blk].JCl;
		      j_inc = -1;
		      k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
		      k_max = Soln_Blks[i_blk].KCu;
		      k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.k> 0) {
		      i_min = Soln_Blks[i_blk].ICu;
		      i_max = Soln_Blks[i_blk].ICl;
		      i_inc = -1;
		      j_min = Soln_Blks[i_blk].JCu;
		      j_max = Soln_Blks[i_blk].JCl;
		      j_inc = -1;
		      k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
		      k_max = Soln_Blks[i_blk].KCu;
		      k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.j > 0) {
		      i_min = Soln_Blks[i_blk].ICu;
		      i_max = Soln_Blks[i_blk].ICl;
		      i_inc = -1;
		      j_min = Soln_Blks[i_blk].JCl;
		      j_max = Soln_Blks[i_blk].JCu;
		      j_inc = 1;
		      k_min = Soln_Blks[i_blk].KCu;
		      k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
		      k_inc = -1;
	} else if (Soln_Block_List.Block[i_blk].infoT[0].dimen.i > 0) {
		      i_min = Soln_Blks[i_blk].ICl;
		      i_max = Soln_Blks[i_blk].ICu;
		      i_inc = 1;
		      j_min = Soln_Blks[i_blk].JCu;
		      j_max = Soln_Blks[i_blk].JCl;
		      j_inc = -1;
		      k_min = Soln_Blks[i_blk].KCu;
		      k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
		      k_inc = -1;
	} else {
		      i_min = Soln_Blks[i_blk].ICu;
		      i_max = Soln_Blks[i_blk].ICl;
		      i_inc = -1;
		      j_min = Soln_Blks[i_blk].JCu;
		      j_max = Soln_Blks[i_blk].JCl;
		      j_inc = -1;
		      k_min = Soln_Blks[i_blk].KCu;
		      k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost+1;
		      k_inc = -1;
	} /* endif */
	for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ){
	    l = l + 1;
	    if (l >= buffer_size_neighbour) return(2202);
	    Soln_Block_List.message_noreschange_topface_sendbuf[i_blk][l] =
	      double(Soln_Blks[i_blk].Grid.BCtypeW[j][k]);
	  } /* endfor */
	for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ){
		      l = l + 1;
		  if (l >= buffer_size_neighbour) return(2203);
		  Soln_Block_List.message_noreschange_topface_sendbuf[i_blk][l] =
		    double(Soln_Blks[i_blk].Grid.BCtypeE[j][k]);
	  } /* endfor */
	for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	  for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		    l = l + 1;
		      if (l >= buffer_size_neighbour) return(2204);
		    Soln_Block_List.message_noreschange_topface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeS[i][k]);
	  } /* endfor */
	for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	  for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		    l = l + 1;
		      if (l >= buffer_size_neighbour) return(2205);
		    Soln_Block_List.message_noreschange_topface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeN[i][k]);
	  } /* endfor */

      } /* endif */
    } /* endif */
  
 
    
       ////////////////////////////////////////////////////////////////////
       // Load send buffer for Bottom face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////
    if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nB == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoB[0].level)) {
      if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost*
                                     abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.i)*
                                    abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.j)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.i)+1)*
                                        (abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR3D+
	    4*Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.i)+
	    4*Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.j);
      } else {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.i)+1)*
                                        (abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR3D+
	    4*Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.i)+
	    4*Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.j);
      } /* endif */
      l = -1;
          // Load ghost cell solution variable information as required.
      if (!Send_Mesh_Geometry_Only) {
	if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
 	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
	        k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
 	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
                i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
 	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
	        k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
	        k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
	        k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
 	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	} else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	} /* endif */

	i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_bottomface_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
						j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	if (i != 0) return(2100);
      } /* endif */
          // Load ghost cell mesh information as required.
      if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost;;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_min-1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost;;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_max-1;
	} else if ( Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost;;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_min-1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost;;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_min-1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost;;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_min-1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost;;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_max-1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost;;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_max-1;
	} else {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost;;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min;
                k_ref = k_max-1;
	} /* endif */

	x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
	for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
	      for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(2101);
                      Soln_Block_List.message_noreschange_bottomface_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
	      } /* endfor */

	if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0 &&
	    Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0 &&
	    Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	  i_min = Soln_Blks[i_blk].ICl;
	  i_max = Soln_Blks[i_blk].ICu;
	  i_inc = 1;
	  j_min = Soln_Blks[i_blk].JCl;
	  j_max = Soln_Blks[i_blk].JCu;
	  j_inc = 1;
	  k_min = Soln_Blks[i_blk].KCl;
	  k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
	  k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0 &&
                  Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	} else if ( Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
                i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
	        k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
	        k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
	        k_inc = 1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	} else if (Soln_Block_List.Block[i_blk].infoB[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	} else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost-1;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	} /* endif */

	for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
		  l = l + 1;
		  if (l >= buffer_size_neighbour) return(2102);
                   Soln_Block_List.message_noreschange_bottomface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeW[j][k]);
	  } /* endfor */
	for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ){
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2103);
                   Soln_Block_List.message_noreschange_bottomface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeE[j][k]);
	  } /* endfor */
	for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	  for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2104);
                   Soln_Block_List.message_noreschange_bottomface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeS[i][k]);
	  } /* endfor */
	for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	  for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2105);
                   Soln_Block_List.message_noreschange_bottomface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeN[i][k]);
	  } /* endfor */

      } /* endif */
    } /* endif */
  

       ////////////////////////////////////////////////////////////////////
       // Load send buffer for North face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nN == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoN[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                     abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)*
                                     abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.k)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+1)*
                                        (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+
		  4*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.k);
          } else {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+1)*
                                        (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+
		  4*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.k);
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
	    if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
                i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	    } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	    } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_northface_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
						j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(2300);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	    if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_min;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_min;
	    } else if ( Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_min;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_min;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_min;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_min;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_min;
	    } else {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_min;
	    } /* endif */
	    x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
	    for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	      for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
		  for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(2301);
                      Soln_Block_List.message_noreschange_northface_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
		  } /* endfor */
	    if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                  Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                  Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	    } else if ( Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
                i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.k> 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	    } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
	    } /* endif */

	    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	      for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
		  l = l + 1;
                   if (l >= buffer_size_neighbour) return(2302);
                   Soln_Block_List.message_noreschange_northface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeW[j][k]);
	      } /* endfor */
	    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	      for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2303);
                   Soln_Block_List.message_noreschange_northface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeE[j][k]);
	      } /* endfor */
	    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	      for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2304);
                   Soln_Block_List.message_noreschange_northface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeB[i][j]);
	      } /* endfor */
	    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	      for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2305);
                   Soln_Block_List.message_noreschange_northface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeT[i][j]);
	      } /* endfor */

          } /* endif */
       } /* endif */


       ////////////////////////////////////////////////////////////////////
       // Load send buffer for South face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                     abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)*
                                     abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.k)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+1)*
                                        (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+
		  4*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.k);
          } else {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+1)*
                                        (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+
		  4*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.k);
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
         if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
 	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
 	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
                i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
 	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
 	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */


             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_southface_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
						j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(2400);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if ( Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_min;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_min;
             } /* endif */

             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(2401);
                      Soln_Block_List.message_noreschange_southface_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */

              if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                  Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                  Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
              } else if ( Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
                i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCl;
                k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
                k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
		  for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
		  l = l + 1;
                   if (l >= buffer_size_neighbour) return(2402);
                   Soln_Block_List.message_noreschange_southface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeW[j][k]);
                } /* endfor */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
		  for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2403);
                   Soln_Block_List.message_noreschange_southface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeE[j][k]);
                } /* endfor */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
		  for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2404);
                   Soln_Block_List.message_noreschange_southface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeB[i][j]);
                } /* endfor */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
		  for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2405);
                   Soln_Block_List.message_noreschange_southface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeT[i][j]);
		  } /* endfor */

          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////
       // Load send buffer for East face neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                     abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)*
                                     abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.k)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+1)*
                                        (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+
		  4*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.k);
           } else {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+1)*
                                        (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+
		  4*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.k);
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else  if ( Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(1200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_min;
             } else  if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
               i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
		i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_min;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_min;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(1201);
                      Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
             if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } /* endif */
                for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
		  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
		  l = l + 1;
                   if (l >= buffer_size_neighbour) return(1202);
                   Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeB[i][j]);
                } /* endfor */
                for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
		  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1203);
                   Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeT[i][j]);
                } /* endfor */
                for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
		  for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1204);
                   Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeS[i][k]);
                } /* endfor */
                for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
		  for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1205);
                   Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeN[i][k]);
                } /* endfor */



          } /* endif */
       } /* endif */
       
       
       ///////////////////////////////////////////////////////////////////
       // Load send buffer for West face neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                     abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)*
                                     abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.k)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+1)*
                                        (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+
		  4*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.k);
          } else {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+1)*
                                        (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+
		  4*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.k);
          } /* endif */
          l = -1;
         // Load ghost cell solution variable information as required.
	  if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else  if ( Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_westface_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(1300);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_min;
             } else  if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
               i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
		i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min;
                k_ref = k_min;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_min;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(1301);
                      Soln_Block_List.message_noreschange_westface_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
             if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } /* endif */

                 for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
		   for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
		  l = l + 1;
                   if (l >= buffer_size_neighbour) return(1302);
                   Soln_Block_List.message_noreschange_westface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeB[i][j]);
                } /* endfor */
                for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
		  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1303);
                   Soln_Block_List.message_noreschange_westface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeT[i][j]);
                } /* endfor */
                for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
		  for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ){
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1304);
                   Soln_Block_List.message_noreschange_westface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeS[i][k]);
                } /* endfor */
                for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
		  for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1305);
                   Soln_Block_List.message_noreschange_westface_sendbuf[i_blk][l] =
                      double(Soln_Blks[i_blk].Grid.BCtypeN[i][k]);
                } /* endfor */

   
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for North West corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoNW[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoNW[0].dimen.k)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoNW[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoNW[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if ( Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */

           i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_northwestcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(4200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_max+1;
                k_ref = k_min  ;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
			Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_max+1;
                k_ref = k_min  ;
             } else if ( Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_max+1;
                k_ref = k_min  ;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
			Soln_Block_List.Block[i_blk].infoNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min+1;
                k_ref = k_min  ;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min+1;
                k_ref = k_min  ;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_max+1;
                k_ref = k_min  ;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min+1;
		k_ref = k_min  ;
            } else {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min+1;
                k_ref = k_min  ;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4201);
                      Soln_Block_List.message_noreschange_northwestcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for North East corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoNE[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoNE[0].dimen.k)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoNE[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoNE[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
              } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
            } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */

             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_northeastcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(5200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_max+1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_max+1;
                k_ref = k_min;
              } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_max+1;
                k_ref = k_min;
              } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min+1;
                k_ref = k_min;
              } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min+1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_max+1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min+1;
                k_ref = k_min;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min+1;
                k_ref = k_min;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5201);
                      Soln_Block_List.message_noreschange_northeastcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for South East corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoSE[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoSE[0].dimen.k)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoSE[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoSE[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
   	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	     } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	     } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_southeastcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(4100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_max-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_max-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_max-1;
                k_ref = k_min;
              } else {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_max-1;
                k_ref = k_min;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4101);
                      Soln_Block_List.message_noreschange_southeastcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for South West corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoSW[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoSW[0].dimen.k)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoSW[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour =sqr(Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoSW[0].dimen.k)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
           // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
          } else
             if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
          } else
             if ( Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
          } else
             if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_southwestcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(5100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_max-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl;
	        k_max = Soln_Blks[i_blk].Grid.KNu;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_max-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min-1;
                k_ref = k_min;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_max-1;
                k_ref = k_min;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_max-1;
                k_ref = k_min;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5101);
                      Soln_Block_List.message_noreschange_southwestcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
		   } /* endfor */
          } /* endif */
       } /* endif */

   ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Top North corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTN == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTN[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoTN[0].dimen.i)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoTN[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoTN[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */

          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTN[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        k_inc = -1;
              } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
            } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost+1;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_topnorthcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(5300);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTN[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_min+1;
              } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_max+1;
              } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_max+1;
              } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_min+1;
             } else if (Soln_Block_List.Block[i_blk].infoTN[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_min+1;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTN[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_min+1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5301);
                      Soln_Block_List.message_noreschange_topnorthcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

        ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Top East corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTE[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoTE[0].dimen.j)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoTE[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoTE[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        k_inc = -1;
              } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
            } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost+1;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_topeastcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(5400);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_min+1;
              } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_max+1;
              } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_max+1;
              } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_min+1;
             } else if (Soln_Block_List.Block[i_blk].infoTE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_min+1;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_min+1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5401);
                      Soln_Block_List.message_noreschange_topeastcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */
        ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Top South corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTS[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoTS[0].dimen.i)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoTS[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoTS[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
            } else  if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTS[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost+1;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
   	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	     } else if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	     } else if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost-1;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost+1;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost+1;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost+1;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_topsouthcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(4300);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTS[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_min+1;
             } else if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_min+1;
             } else if (Soln_Block_List.Block[i_blk].infoTS[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        k_inc = -1;
               i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_min+1;
              } else {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTS[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_min+1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4301);
                      Soln_Block_List.message_noreschange_topsouthcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */


       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Top West corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTW[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoTW[0].dimen.j)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoTW[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoTW[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTW[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost+1;
	        k_inc = -1;
             } else if ( Soln_Block_List.Block[i_blk].infoTW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost+1;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost+1;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost+1;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_topwestcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc, k_min,k_max,k_inc);
	     if (i != 0) return(4400);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.i > 0 &&
			Soln_Block_List.Block[i_blk].infoTW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min;
                k_ref = k_min+1  ;
             } else if ( Soln_Block_List.Block[i_blk].infoTW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.i > 0 &&
			Soln_Block_List.Block[i_blk].infoTW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_min+1  ;
             } else if (Soln_Block_List.Block[i_blk].infoTW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min;
		k_ref = k_min+1  ;
            } else {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_min+1  ;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4401);
                      Soln_Block_List.message_noreschange_topwestcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

   ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Bottom North corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBN == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBN[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoBN[0].dimen.i)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoBN[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoBN[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBN[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
              } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost-1;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost-1;
	        k_inc = 1;
            } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost-1;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(5500);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBN[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_max-1;
              } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_min-1;
              } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_min-1;
              } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_max+1;
                k_ref = k_max-1;
             } else if (Soln_Block_List.Block[i_blk].infoBN[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_max-1;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBN[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min+1;
                k_ref = k_max-1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5501);
                      Soln_Block_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

        ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Bottom East corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBE[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoBE[0].dimen.j)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoBE[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoBE[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
              } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost-1;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost-1;
	        k_inc = 1;
            } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost-1;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_bottomeastcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(5600);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_max-1;
              } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_min-1;
              } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_min-1;
              } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_max-1;
             } else if (Soln_Block_List.Block[i_blk].infoBE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min;
                k_ref = k_max-1;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min;
                k_ref = k_max-1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5601);
                      Soln_Block_List.message_noreschange_bottomeastcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */
        ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Bottom South corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBS[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoBS[0].dimen.i)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoBS[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoBS[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBS[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
   	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        k_inc = 1;
	     } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        k_inc = 1;
	     } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(4500);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBS[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_max-1;
             } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_min-1;
                k_ref = k_max-1;
             } else if (Soln_Block_List.Block[i_blk].infoBS[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl;
	        i_max = Soln_Blks[i_blk].Grid.INu;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_max-1;
              } else {
	        i_min = Soln_Blks[i_blk].Grid.INu;
	        i_max = Soln_Blks[i_blk].Grid.INl;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBS[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min;
                j_ref = j_max-1;
                k_ref = k_max-1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4501);
                      Soln_Block_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */


       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Bottom West corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBW[0].level)) {
         if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost)*
                                     abs(Soln_Block_List.Block[i_blk].infoBW[0].dimen.j)*
                                     Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size_neighbour = buffer_size_neighbour +
                                        sqr(Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoBW[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR3D;
          } else {
                buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost)*
                                        (abs(Soln_Block_List.Block[i_blk].infoBW[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBW[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if ( Soln_Block_List.Block[i_blk].infoBW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_bottomwestcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc, k_min,k_max,k_inc);
	     if (i != 0) return(4600);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.i > 0 &&
			Soln_Block_List.Block[i_blk].infoBW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min;
                k_ref = k_max-1 ;
             } else if ( Soln_Block_List.Block[i_blk].infoBW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.i > 0 &&
			Soln_Block_List.Block[i_blk].infoBW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl;
	        j_max = Soln_Blks[i_blk].Grid.JNu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_max-1 ;
             } else if (Soln_Block_List.Block[i_blk].infoBW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min;
		k_ref = k_max-1 ;
            } else {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu;
	        j_max = Soln_Blks[i_blk].Grid.JNl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min;
                k_ref = k_max-1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4601);
                      Soln_Block_List.message_noreschange_bottomwestcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

         ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Topnorth West corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTNW[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost)*
	    Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNW[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        k_inc = -1;
             } else if ( Soln_Block_List.Block[i_blk].infoTNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost+1;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(4700);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_max+1;
                k_ref = k_max+1  ;
             } else if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.i > 0 &&
			Soln_Block_List.Block[i_blk].infoTNW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_max+1;
                k_ref = k_min+1;
             } else if ( Soln_Block_List.Block[i_blk].infoTNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_max+1;
                k_ref = k_max+1  ;
             } else if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.i > 0 &&
			Soln_Block_List.Block[i_blk].infoTNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min+1;
                k_ref = k_max+1  ;
             } else if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min+1;
                k_ref = k_max+1  ;
             } else if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_max+1;
                k_ref = k_min+1  ;
             } else if (Soln_Block_List.Block[i_blk].infoTNW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min+1;
		k_ref = k_min+1  ;
            } else {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min+1;
                k_ref = k_min+1  ;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4701);
                      Soln_Block_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */


       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Topsouth East corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTSE[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost)*
	    Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
   	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	     } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
	     } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost-1;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost+1;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(4800);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min-1;
                k_ref = k_min+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_max-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_max-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min-1;
                k_ref = k_min+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_max-1;
                k_ref = k_min+1;
              } else {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_max-1;
                k_ref = k_min+1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4801);
                      Soln_Block_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Topsouth West corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTSW[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost)*
	    Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
          } else
             if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSW[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost+1;
	        k_inc = -1;
          } else
             if ( Soln_Block_List.Block[i_blk].infoTSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost+1;
	        k_inc = -1;
          } else
             if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost+1;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost+1;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost+1;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost+1;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(5800);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.i > 0 &&
                Soln_Block_List.Block[i_blk].infoTSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.i > 0 &&
                Soln_Block_List.Block[i_blk].infoTSW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min-1;
                k_ref = k_min+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_max-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_max-1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min-1;
                k_ref = k_min+1;
             } else if (Soln_Block_List.Block[i_blk].infoTSW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_max-1;
                k_ref = k_min+1;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTSW[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_max-1;
                k_ref = k_min+1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5801);
                      Soln_Block_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
		   } /* endfor */
          } /* endif */
       } /* endif */

      ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Bottomnorth West corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBNW[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost)*
	    Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNW[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if ( Soln_Block_List.Block[i_blk].infoBNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
	        i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc, k_min,k_max,k_inc);
	     if (i != 0) return(4900);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_max+1;
                k_ref = k_min-1  ;
             } else if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.i > 0 &&
			Soln_Block_List.Block[i_blk].infoBNW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_max+1;
                k_ref = k_max-1;
             } else if ( Soln_Block_List.Block[i_blk].infoBNW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_max+1;
                k_ref = k_min-1  ;
             } else if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.i > 0 &&
			Soln_Block_List.Block[i_blk].infoBNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min+1;
                k_ref = k_min-1  ;
             } else if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min+1;
                k_ref = k_min-1  ;
             } else if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_max+1;
                k_ref = k_max-1  ;
             } else if (Soln_Block_List.Block[i_blk].infoBNW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min+1;
		k_ref = k_max-1  ;
            } else {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min+1;
                k_ref = k_max-1  ;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4901);
                      Soln_Block_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Bottomnorth East corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBNE[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost)*
	    Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
              } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost-1;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost-1;
	        k_inc = 1;
            } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost-1;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(5900);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_max+1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_max+1;
                k_ref = k_max-1;
              } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_max+1;
                k_ref = k_min-1;
              } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min+1;
                k_ref = k_min-1;
              } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min+1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_max+1;
                k_ref = k_max-1;
             } else if (Soln_Block_List.Block[i_blk].infoBNE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min+1;
                k_ref = k_max-1;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBNE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min+1;
                k_ref = k_max-1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5901);
                      Soln_Block_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Bottomsouth East corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBSE[0].level)) {
             buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost)*
	       Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        j_inc = 1;
   	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
   	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        k_inc = 1;
	     } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
 	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        k_inc = 1;
	     } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(5000);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min-1;
                k_ref = k_max-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_max-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_max-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min-1;
                k_ref = k_max-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_max-1;
                k_ref = k_max-1;
              } else {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_max-1;
                k_ref = k_max-1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5001);
                      Soln_Block_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Topnorth East corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTNE[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost)*
	    Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        k_inc = -1;
              } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
                i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
            } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        k_max = Soln_Blks[i_blk].KCu;
	        k_inc = 1;
              } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        j_max = Soln_Blks[i_blk].JCu;
	        j_inc = 1;
 	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        i_max = Soln_Blks[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCu;
	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICu;
	        i_max = Soln_Blks[i_blk].ICu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCu;
 	        j_max = Soln_Blks[i_blk].JCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCu;
	        k_max = Soln_Blks[i_blk].KCu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost+1;
	        k_inc = -1;
             } /* endif */
             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_topnortheastcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(5700);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_max+1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNE[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_max+1;
                k_ref = k_min+1;
              } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_max+1;
                k_ref = k_max+1;
              } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoTNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min+1;
                k_ref = k_max+1;
              } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNu-1;
	        k_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min+1;
                k_ref = k_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNu-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_max+1;
                k_ref = k_min+1;
             } else if (Soln_Block_List.Block[i_blk].infoTNE[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min+1;
                k_ref = k_min+1;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INu-1;
	        i_max = Soln_Blks[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNu-1;
	        j_max = Soln_Blks[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNu-1;
	        k_max = Soln_Blks[i_blk].Grid.KNu-Soln_Block_List.Block[i_blk].infoTNE[0].dimen.ghost;
	        k_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min+1;
                k_ref = k_min+1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5701);
                      Soln_Block_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
                   } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for Bottomsouth West corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBSW[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost)*
	    Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        k_inc = 1;
          } else
             if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSW[0].dimen.j > 0 ) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
          } else
             if ( Soln_Block_List.Block[i_blk].infoBSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
          } else
             if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
            } else if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl;
	        k_max = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        k_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl;
	        j_max = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].ICl;
                i_max = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } else {
	        i_min = Soln_Blks[i_blk].ICl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        i_max = Soln_Blks[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_Blks[i_blk].JCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        j_max = Soln_Blks[i_blk].JCl;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].KCl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost-1;
	        k_max = Soln_Blks[i_blk].KCl;
	        k_inc = -1;
             } /* endif */

             i = Soln_Blks[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                k_min,k_max,k_inc);
	     if (i != 0) return(6000);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.i > 0 &&
                Soln_Block_List.Block[i_blk].infoBSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.i > 0 &&
                Soln_Block_List.Block[i_blk].infoBSW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min-1;
                k_ref = k_max-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.j > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoBSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_min-1;
                j_ref = j_max-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.k > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+1;
	        k_max = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        k_inc = 1;
                i_ref = i_max-1;
                j_ref = j_max-1;
                k_ref = k_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.j > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+1;
	        j_max = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        j_inc = 1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min-1;
                k_ref = k_max-1;
             } else if (Soln_Block_List.Block[i_blk].infoBSW[0].dimen.i > 0) {
	        i_min = Soln_Blks[i_blk].Grid.INl+1;
	        i_max = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_min-1;
                j_ref = j_max-1;
                k_ref = k_max-1;
             } else {
	        i_min = Soln_Blks[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        i_max = Soln_Blks[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_Blks[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        j_max = Soln_Blks[i_blk].Grid.JNl+1;
	        j_inc = -1;
	        k_min = Soln_Blks[i_blk].Grid.KNl+Soln_Block_List.Block[i_blk].infoBSW[0].dimen.ghost;
	        k_max = Soln_Blks[i_blk].Grid.KNl+1;
	        k_inc = -1;
                i_ref = i_max-1;
                j_ref = j_max-1;
                k_ref = k_max-1;
             } /* endif */
             x_ref = Soln_Blks[i_blk].Grid.Node[i_ref][j_ref][k_ref].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
                   for ( m = 1 ; m <= NUM_COMP_VECTOR3D; ++ m) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(6001);
                      Soln_Block_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][l] =
                         Soln_Blks[i_blk].Grid.Node[i][j][k].X[m]-x_ref[m];
		   } /* endfor */
          } /* endif */
       } /* endif */

       //tims end

       }   /* endfor */

    /* Loading of send buffers complete.  Return zero value. */

    return(0);

}

/***************************************************************
 * Routine: Load_Send_Message_Buffers_ResChange_FineToCoarse   *
 *                                                             *
 * Loads the send message passing buffers for sending          *
 * solution information from more refined (higher mesh         *
 * resolution) solution blocks to more coarse neighbouring 3D  *
 * hexahedrial multi-block solution blocks with lower mesh     *
 * resolution (works on entire 1D array of 3D solution         *
 * blocks).                                                    *
 *                                                             *
 ***************************************************************/
template <class Hexa_Soln_Block>
int Load_Send_Message_Buffers_ResChange_FineToCoarse(Hexa_Soln_Block *Soln_Blks,
                                                     AdaptiveBlock3D_List &Soln_Block_List,
                                                     const int Number_of_Solution_Variables,
                                                     const int Send_Mesh_Geometry_Only,
                                                     const int Send_Conservative_Solution_Fluxes) {

  cout << "\nERROR: Load_Send_Message_Buffers_ResChange_FineToCoarse not defined for 3 dimensions\n";
  return(2);

}

/***************************************************************
 * Routine: Load_Send_Message_Buffers_ResChange_CoarseToFine   *
 *                                                             *
 * Loads the send message passing buffers for sending          *
 * solution information from more coarse (lower mesh           *
 * resolution) solution blocks to more refined neighbouring 3D *
 * hexahedrial multi-block solution blocks with higher mesh    *
 * resolution (works on entire 1D array of 3D solution         *
 * blocks).                                                    *
 *                                                             *
 ***************************************************************/
template <class Hexa_Soln_Block>
int Load_Send_Message_Buffers_ResChange_CoarseToFine(Hexa_Soln_Block *Soln_Blks,
                                                     AdaptiveBlock3D_List &Soln_Block_List,
                                                     const int Number_of_Solution_Variables,
                                                     const int Send_Mesh_Geometry_Only) {

  cout << "\nERROR: Load_Send_Message_Buffers_ResChange_CoarseToFine not defined for 3 dimensions\n";
  return(2);

}

/********************************************************
 * Routine: Unload_Receive_Message_Buffers_NoResChange  *
 *                                                      *
 * Unloads the receive message passing buffers for      *
 * receiving solution information from neighbouring 3D  *
 * hexahedrial multi-block solution blocks with no    *
 * mesh resolution changes (works on entire 1D array of *
 * 3D solution blocks).                                 *
 *                                                      *
 ********************************************************/
template <class Hexa_Soln_Block>
int Unload_Receive_Message_Buffers_NoResChange(Hexa_Soln_Block *Soln_Blks,
                                               AdaptiveBlock3D_List &Soln_Block_List,
                                               const int Number_of_Solution_Variables,
                                               const int Send_Mesh_Geometry_Only) {

    int i_blk, buffer_size, 
        i_min, i_max, i_inc, i, 
        j_min, j_max, j_inc, j, 
        k_min, k_max, k_inc, k, 
        l;
    Vector3D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != NUM_COMP_VECTOR3D &&
        Number_of_Solution_Variables != Soln_Blks[0].NumVar() &&
        Number_of_Solution_Variables != Soln_Blks[0].NumVar() + NUM_COMP_VECTOR3D) {
      return(1);
    } /* endif */

    /* Unload the receive buffers for each solution block. */

    for ( i_blk = 0 ; i_blk <= Soln_Block_List.Nblk-1 ; ++i_blk ) {

      /////////////////////////////////////////////////////////////////////////
      // Unload receive buffer for Top face neighbours of solution blocks. //
      /////////////////////////////////////////////////////////////////////////
      if (Soln_Block_List.Block[i_blk].used &&
	  (Soln_Block_List.Block[i_blk].nT == 1) &&
	  (Soln_Block_List.Block[i_blk].info.level ==
	   Soln_Block_List.Block[i_blk].infoT[0].level)) {
	if (!Send_Mesh_Geometry_Only) {
	  buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
	    abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
	    abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
	    Soln_Blks[i_blk].NumVar();
	  if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
	    buffer_size = buffer_size +
	      Soln_Block_List.Block[i_blk].info.dimen.ghost*
	      (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
	      (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
	      NUM_COMP_VECTOR3D+
	      4*Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.i)+
	      4*Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.j);
	} else {
	  buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
	    (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
	    (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
	    NUM_COMP_VECTOR3D+
	    4*Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.i)+
	    4*Soln_Block_List.Block[i_blk].infoT[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoT[0].dimen.j);
	} /* endif */
	l = -1;
	// Unload ghost cell solution information as required.
	if (!Send_Mesh_Geometry_Only) {
	  i_min = Soln_Blks[i_blk].ICl;
	  i_max = Soln_Blks[i_blk].ICu;
	  i_inc = 1;
	  j_min = Soln_Blks[i_blk].JCl;
	  j_max = Soln_Blks[i_blk].JCu;
	  j_inc = 1;
	  k_min = Soln_Blks[i_blk].KCu+1;
	  k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	  k_inc = 1;
	  i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_topface_recbuf[i_blk],
						  l,buffer_size,
						  i_min,i_max,i_inc,j_min,j_max,j_inc,
						  k_min,k_max,k_inc);
	  if (i != 0) return(2200);
	} /* endif */
          // Unload ghost cell mesh information as required.
	if (Send_Mesh_Geometry_Only ||
	    Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	  i_min = Soln_Blks[i_blk].Grid.INl;
	  i_max = Soln_Blks[i_blk].Grid.INu;
	  i_inc = 1;
	  j_min = Soln_Blks[i_blk].Grid.JNl;
	  j_max = Soln_Blks[i_blk].Grid.JNu;
	  j_inc = 1;
	  k_min = Soln_Blks[i_blk].Grid.KNu+1;
	  k_max = Soln_Blks[i_blk].Grid.KNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	  k_inc = 1;

/* 	cout<<"\nat Unload Top buffer_size ="<<buffer_size; */
/* 	cout<<" i_min= "<<i_min<< " i_max= "<<i_max<< " i_inc= "<<i_inc << " j_min= "<<j_min<< " j_max= "<<j_max */
/* 	    << " j_inc= "<<j_inc<< " k_min= "<< k_min << " k_max= "<<k_max<< " k_inc= "<<k_inc <<" l="<<l; */

	  x_ref = Soln_Blks[i_blk].Grid.Node[i_min][j_min][k_min-1].X; // Reference node location.
	  for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	      for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		l = l + NUM_COMP_VECTOR3D;
		if (l >= buffer_size) return(2201);
		Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		  Vector3D(Soln_Block_List.message_noreschange_topface_recbuf[i_blk][l-2],
			   Soln_Block_List.message_noreschange_topface_recbuf[i_blk][l-1],
			   Soln_Block_List.message_noreschange_topface_recbuf[i_blk][l])+x_ref;
	    } /* endfor */
	    i_min = Soln_Blks[i_blk].ICl;
	    i_max = Soln_Blks[i_blk].ICu;
	    i_inc = 1;
	    j_min = Soln_Blks[i_blk].JCl;
	    j_max = Soln_Blks[i_blk].JCu;
	    j_inc = 1;
	    k_min = Soln_Blks[i_blk].KCu+1;
	    k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	    k_inc = 1;
	    for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	      for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		  Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i,j,k);
		  Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
	    for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	      for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2202);
                Soln_Blks[i_blk].Grid.BCtypeW[j][k] =
		  int(Soln_Block_List.message_noreschange_topface_recbuf[i_blk][l]);
	      } /* endfor */
	    for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	      for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ){
                l = l + 1;
                if (l >= buffer_size) return(2203);
                Soln_Blks[i_blk].Grid.BCtypeE[j][k] =
		  int(Soln_Block_List.message_noreschange_topface_recbuf[i_blk][l]);
	      } /* endfor */
	    for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	      for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2204);
                Soln_Blks[i_blk].Grid.BCtypeS[i][k] =
		  int(Soln_Block_List.message_noreschange_topface_recbuf[i_blk][l]);
	      } /* endfor */
	    for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	      for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2205);
                Soln_Blks[i_blk].Grid.BCtypeN[i][k] =
		  int(Soln_Block_List.message_noreschange_topface_recbuf[i_blk][l]);
	      } /* endfor */
          } /* endif */
	} /* endif */


	/////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Bottom face neighbours of solution blocks. //
       /////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nB == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoB[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR3D+
	    4*Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.i)+
	    4*Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.j);
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR3D+
	    4*Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.i)+
	    4*Soln_Block_List.Block[i_blk].infoB[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoB[0].dimen.j);
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
   	     i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl;
	     j_max = Soln_Blks[i_blk].JCu;
	     j_inc = 1;
   	     k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
            i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_bottomface_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                     k_min,k_max,k_inc);
	     if (i != 0) return(2100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl;
	     i_max = Soln_Blks[i_blk].Grid.INu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl;
	     j_max = Soln_Blks[i_blk].Grid.JNu;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_max = Soln_Blks[i_blk].Grid.KNl-1;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min][j_min][k_max+1].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(2101);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_bottomface_recbuf[i_blk][l-2],
			      Soln_Block_List.message_noreschange_bottomface_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_bottomface_recbuf[i_blk][l])+x_ref;
                } /* endfor */
  	     i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl;
	     j_max = Soln_Blks[i_blk].JCu;
	     j_inc = 1;
   	     k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
              for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
 	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
		 for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2102);
                Soln_Blks[i_blk].Grid.BCtypeW[j][k] =
                   int(Soln_Block_List.message_noreschange_bottomface_recbuf[i_blk][l]);
             } /* endfor */
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
		 for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ){
                l = l + 1;
                if (l >= buffer_size) return(2103);
                Soln_Blks[i_blk].Grid.BCtypeE[j][k] =
                   int(Soln_Block_List.message_noreschange_bottomface_recbuf[i_blk][l]);
             } /* endfor */
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2104);
                Soln_Blks[i_blk].Grid.BCtypeS[i][k] =
                   int(Soln_Block_List.message_noreschange_bottomface_recbuf[i_blk][l]);
             } /* endfor */
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
		 for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2105);
                Soln_Blks[i_blk].Grid.BCtypeN[i][k] =
                   int(Soln_Block_List.message_noreschange_bottomface_recbuf[i_blk][l]);
		 } /* endfor */
          } /* endif */
       } /* endif */


       /////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for North face neighbours of solution blocks. //
       /////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nN == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoN[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.k)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+
		  4*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.k);
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+
		  4*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.k);
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCu+1;
             j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl;
             k_max = Soln_Blks[i_blk].KCu;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_northface_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                     k_min,k_max,k_inc);
	     if (i != 0) return(2300);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl;
	     i_max = Soln_Blks[i_blk].Grid.INu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNu+1;
	     j_max = Soln_Blks[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl;
	     k_max = Soln_Blks[i_blk].Grid.KNu;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min][j_min-1][k_min].X; // Reference node location.
	     for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(2301);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_northface_recbuf[i_blk][l-2],
			      Soln_Block_List.message_noreschange_northface_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_northface_recbuf[i_blk][l])+x_ref;
             } /* endfor */
             i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCu+1;
	     j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i,j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2302);
                Soln_Blks[i_blk].Grid.BCtypeW[j][k] =
                   int(Soln_Block_List.message_noreschange_northface_recbuf[i_blk][l]);
             } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2303);
                Soln_Blks[i_blk].Grid.BCtypeE[j][k] =
                   int(Soln_Block_List.message_noreschange_northface_recbuf[i_blk][l]);
             } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	       for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2304);
                Soln_Blks[i_blk].Grid.BCtypeB[i][j] =
                   int(Soln_Block_List.message_noreschange_northface_recbuf[i_blk][l]);
             } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	       for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2305);
                Soln_Blks[i_blk].Grid.BCtypeT[i][j] =
                   int(Soln_Block_List.message_noreschange_northface_recbuf[i_blk][l]);
             } /* endfor */
          } /* endif */
       } /* endif */


       /////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for South face neighbours of solution blocks. //
       /////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.k)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+
		  4*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.k);
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+
		  4*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.k);
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
   	     i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_Blks[i_blk].JCl-1;
	     j_inc = 1;
   	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_southface_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                     k_min,k_max,k_inc);
	     if (i != 0) return(2400);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl;
	     i_max = Soln_Blks[i_blk].Grid.INu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_Blks[i_blk].Grid.JNl-1;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl;
	     k_max = Soln_Blks[i_blk].Grid.KNu;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min][j_max+1][k_min].X; // Reference node location.
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(2401);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_southface_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_southface_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_southface_recbuf[i_blk][l])+x_ref;
                } /* endfor */
  	     i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_Blks[i_blk].JCl-1;
	     j_inc = 1;
  	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
             for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
             } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2402);
                Soln_Blks[i_blk].Grid.BCtypeW[j][k] =
                   int(Soln_Block_List.message_noreschange_southface_recbuf[i_blk][l]);
             } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2403);
                Soln_Blks[i_blk].Grid.BCtypeE[j][k] =
                   int(Soln_Block_List.message_noreschange_southface_recbuf[i_blk][l]);
             } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	       for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2404);
                Soln_Blks[i_blk].Grid.BCtypeB[i][j] =
                   int(Soln_Block_List.message_noreschange_southface_recbuf[i_blk][l]);
             } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
	       for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2405);
                Soln_Blks[i_blk].Grid.BCtypeT[i][j] =
                   int(Soln_Block_List.message_noreschange_southface_recbuf[i_blk][l]);
             } /* endfor */


          } /* endif */
       } /* endif */


       ////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for East face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.k)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+
		  4*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.k);
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+
		  4*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.k);
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
   	     i_min = Soln_Blks[i_blk].ICu+1;
 	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl;
	     j_max = Soln_Blks[i_blk].JCu;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_eastface_recbuf[i_blk],
                                                     l,buffer_size, i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                     k_min,k_max,k_inc);
	     if (i != 0) return(1200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INu+1;
	     i_max = Soln_Blks[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl;
	     j_max = Soln_Blks[i_blk].Grid.JNu;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl;
	     k_max = Soln_Blks[i_blk].Grid.KNu;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min-1][j_min][k_min].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(1201);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_eastface_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_eastface_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_eastface_recbuf[i_blk][l])+x_ref;
                } /* endfor */
 	     i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl;
	     j_max = Soln_Blks[i_blk].JCu;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
             for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1202);
                Soln_Blks[i_blk].Grid.BCtypeS[i][k] =
                   int(Soln_Block_List.message_noreschange_eastface_recbuf[i_blk][l]);
             } /* endfor */
             for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1203);
                Soln_Blks[i_blk].Grid.BCtypeN[i][k] =
                   int(Soln_Block_List.message_noreschange_eastface_recbuf[i_blk][l]);
             } /* endfor */
	       for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
		 for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1204);
                Soln_Blks[i_blk].Grid.BCtypeB[i][j] =
                   int(Soln_Block_List.message_noreschange_eastface_recbuf[i_blk][l]);
             } /* endfor */
	     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1205);
                Soln_Blks[i_blk].Grid.BCtypeT[i][j] =
                   int(Soln_Block_List.message_noreschange_eastface_recbuf[i_blk][l]);
             } /* endfor */

          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for West face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.k)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+
		  4*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.k);
         } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D+
		  4*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+
		  4*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.k);
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
  	     i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl;
	     j_max = Soln_Blks[i_blk].JCu;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_westface_recbuf[i_blk],
                                                     l,buffer_size,
		 i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                     k_min,k_max,k_inc);
	     if (i != 0) return(1100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl;
	     j_max = Soln_Blks[i_blk].Grid.JNu;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl;
	     k_max = Soln_Blks[i_blk].Grid.KNu;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_max+1][j_min][k_min].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(1101);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_westface_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_westface_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_westface_recbuf[i_blk][l])+x_ref;
                } /* endfor */
     	     i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl;
	     j_max = Soln_Blks[i_blk].JCu;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */

             for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1102);
                Soln_Blks[i_blk].Grid.BCtypeS[i][k] =
                   int(Soln_Block_List.message_noreschange_westface_recbuf[i_blk][l]);
             } /* endfor */
             for ( i  = i_min ; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1103);
                Soln_Blks[i_blk].Grid.BCtypeN[i][k] =
                   int(Soln_Block_List.message_noreschange_westface_recbuf[i_blk][l]);
             } /* endfor */
	       for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
		 for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1104);
                Soln_Blks[i_blk].Grid.BCtypeB[i][j] =
                   int(Soln_Block_List.message_noreschange_westface_recbuf[i_blk][l]);
             } /* endfor */
	     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc )
	       for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1105);
                Soln_Blks[i_blk].Grid.BCtypeT[i][j] =
                   int(Soln_Block_List.message_noreschange_westface_recbuf[i_blk][l]);
             } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for North West corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoNW[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.k)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
		  sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
		  (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
		  NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCu+1;
             j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_northwestcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(4200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNu+1;
	     j_max = Soln_Blks[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl;
	     k_max = Soln_Blks[i_blk].Grid.KNu;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_max+1][j_min-1][k_min].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(4201);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_northwestcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_northwestcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_northwestcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCu+1;
	     j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for North East corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoNE[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.k)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCu+1;
             j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_northeastcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(5200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INu+1;
	     i_max = Soln_Blks[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNu+1;
	     j_max = Soln_Blks[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl;
	     k_max = Soln_Blks[i_blk].Grid.KNu;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min-1][j_min-1][k_min].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(5201);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_northeastcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_northeastcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_northeastcorner_recbuf[i_blk][l])+x_ref;
              } /* endfor */
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCu+1;
	     j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for South East corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoSE[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.k)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
             j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_southeastcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(4100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INu+1;
	     i_max = Soln_Blks[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_Blks[i_blk].Grid.JNl-1;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl;
	     k_max = Soln_Blks[i_blk].Grid.KNu;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min-1][j_max+1][k_min].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(4101);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_southeastcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_southeastcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_southeastcorner_recbuf[i_blk][l])+x_ref;
             } /* endfor */
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for South West corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoSW[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.k)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
             j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_southwestcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(5100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_Blks[i_blk].Grid.JNl-1;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl;
	     k_max = Soln_Blks[i_blk].Grid.KNu;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_max+1][j_max+1][k_min].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(5101);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_southwestcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_southwestcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_southwestcorner_recbuf[i_blk][l])+x_ref;
             } /* endfor */
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
	     i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].KCl;
	     k_max = Soln_Blks[i_blk].KCu;
	     k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
             } /* endfor */
          } /* endif */
       } /* endif */


     ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Topnorth  corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTN == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTN[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCu+1;
             j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_topnorthcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(5300);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl;
	     i_max = Soln_Blks[i_blk].Grid.INu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNu+1;
	     j_max = Soln_Blks[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNu+1;
	     k_max = Soln_Blks[i_blk].Grid.KNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min][j_min-1][k_min-1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(5301);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_topnorthcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_topnorthcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_topnorthcorner_recbuf[i_blk][l])+x_ref;
              } /* endfor */
             i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
             i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCu+1;
	     j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */



       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Topsouth  corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTS[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
             i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_topsouthcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(5400);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl;
	     i_max = Soln_Blks[i_blk].Grid.INu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_Blks[i_blk].Grid.JNl-1;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNu+1;
	     k_max = Soln_Blks[i_blk].Grid.KNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min][j_max+1][k_min-1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(5401);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_topsouthcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_topsouthcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_topsouthcorner_recbuf[i_blk][l])+x_ref;
             } /* endfor */
             i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
             i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
             } /* endfor */
          } /* endif */
       } /* endif */


       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Top East corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTE[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl;
             j_max = Soln_Blks[i_blk].JCu;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_topeastcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(5500);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INu+1;
	     i_max = Soln_Blks[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl;
	     j_max = Soln_Blks[i_blk].Grid.JNu;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNu+1;
	     k_max = Soln_Blks[i_blk].Grid.KNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min-1][j_min][k_min-1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(5501);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_topeastcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_topeastcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_topeastcorner_recbuf[i_blk][l])+x_ref;
             } /* endfor */
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl;
             j_max = Soln_Blks[i_blk].JCu;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */


       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Top West corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTW[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl;
             j_max = Soln_Blks[i_blk].JCu;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_topwestcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(5600);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl;
	     j_max = Soln_Blks[i_blk].Grid.JNu;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNu+1;
	     k_max = Soln_Blks[i_blk].Grid.KNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_max+1][j_min][k_min-1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(5601);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_topwestcorner_recbuf[i_blk][l-2],
 		              Soln_Block_List.message_noreschange_topwestcorner_recbuf[i_blk][l-1],
                             Soln_Block_List.message_noreschange_topwestcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
	     i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl;
             j_max = Soln_Blks[i_blk].JCu;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */


       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Bottomnorth  corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBN == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBN[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCu+1;
             j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_bottomnorthcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(5700);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl;
	     i_max = Soln_Blks[i_blk].Grid.INu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNu+1;
	     j_max = Soln_Blks[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_max = Soln_Blks[i_blk].Grid.KNl-1;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min][j_min-1][k_max+1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(5701);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_bottomnorthcorner_recbuf[i_blk][l-2],
 		              Soln_Block_List.message_noreschange_bottomnorthcorner_recbuf[i_blk][l-1],
                             Soln_Block_List.message_noreschange_bottomnorthcorner_recbuf[i_blk][l])+x_ref;
              } /* endfor */
             i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
             i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCu+1;
	     j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */


       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Bottomsouth  corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBS[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
             i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_bottomsouthcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(5800);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl;
	     i_max = Soln_Blks[i_blk].Grid.INu;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_Blks[i_blk].Grid.JNl-1;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_max = Soln_Blks[i_blk].Grid.KNl-1;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min][j_max+1][k_max+1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(5801);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_bottomsouthcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_bottomsouthcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_bottomsouthcorner_recbuf[i_blk][l])+x_ref;
             } /* endfor */
             i_min = Soln_Blks[i_blk].ICl;
	     i_max = Soln_Blks[i_blk].ICu;
             i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
             } /* endfor */
          } /* endif */
       } /* endif */


       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Bottom East corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBE[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl;
             j_max = Soln_Blks[i_blk].JCu;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_bottomeastcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(5900);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INu+1;
	     i_max = Soln_Blks[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl;
	     j_max = Soln_Blks[i_blk].Grid.JNu;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_max = Soln_Blks[i_blk].Grid.KNl-1;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min-1][j_min][k_max+1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(5901);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_bottomeastcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_bottomeastcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_bottomeastcorner_recbuf[i_blk][l])+x_ref;
             } /* endfor */
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl;
             j_max = Soln_Blks[i_blk].JCu;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */


       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Bottom West corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBW[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
                           Soln_Blks[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                buffer_size = buffer_size +
                              sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR3D;
          } else {
                buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR3D;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl;
             j_max = Soln_Blks[i_blk].JCu;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_bottomwestcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(6000);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl;
	     j_max = Soln_Blks[i_blk].Grid.JNu;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_max = Soln_Blks[i_blk].Grid.KNl-1;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_max+1][j_min][k_max+1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(6001);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_bottomwestcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_bottomwestcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_bottomwestcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
	     i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl;
             j_max = Soln_Blks[i_blk].JCu;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */


       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Topnorth West corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTNW[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
	    Soln_Block_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCu+1;
             j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_topnorthwestcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(6100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNu+1;
	     j_max = Soln_Blks[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNu+1;
	     k_max = Soln_Blks[i_blk].Grid.KNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_max+1][j_min-1][k_min-1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(6101);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_topnorthwestcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_topnorthwestcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_topnorthwestcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCu+1;
	     j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */



       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Topsouth East corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTSE[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
	    Soln_Block_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_topsoutheastcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(6300);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INu+1;
	     i_max = Soln_Blks[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_Blks[i_blk].Grid.JNl-1;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNu+1;
	     k_max = Soln_Blks[i_blk].Grid.KNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min-1][j_max+1][k_min-1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(6301);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_topsoutheastcorner_recbuf[i_blk][l-2],
 		              Soln_Block_List.message_noreschange_topsoutheastcorner_recbuf[i_blk][l-1],
                             Soln_Block_List.message_noreschange_topsoutheastcorner_recbuf[i_blk][l])+x_ref;
             } /* endfor */
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Topnorth East corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTNE[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
	    Soln_Block_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCu+1;
             j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_topnortheastcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(6200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INu+1;
	     i_max = Soln_Blks[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNu+1;
	     j_max = Soln_Blks[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNu+1;
	     k_max = Soln_Blks[i_blk].Grid.KNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min-1][j_min-1][k_min-1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(6201);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_topnortheastcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_topnortheastcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_topnortheastcorner_recbuf[i_blk][l])+x_ref;
              } /* endfor */
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCu+1;
	     j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */


       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Topsouth West corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nTSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoTSW[0].level)) {
	 buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
	   Soln_Block_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_topsouthwestcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(6400);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_Blks[i_blk].Grid.JNl-1;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNu+1;
	     k_max = Soln_Blks[i_blk].Grid.KNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_max+1][j_max+1][k_min-1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(6401);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_topsouthwestcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_topsouthwestcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_topsouthwestcorner_recbuf[i_blk][l])+x_ref;
             } /* endfor */
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
	     i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCu+1;
             k_max = Soln_Blks[i_blk].KCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
             } /* endfor */
          } /* endif */
       } /* endif */



       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Bottomnorth West corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBNW[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
	    Soln_Block_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCu+1;
             j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(6500);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNu+1;
	     j_max = Soln_Blks[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_max = Soln_Blks[i_blk].Grid.KNl-1;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_max+1][j_min-1][k_max+1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(6501);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCu+1;
	     j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Bottomnorth East corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBNE[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
	    Soln_Block_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCu+1;
             j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(6600);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INu+1;
	     i_max = Soln_Blks[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNu+1;
	     j_max = Soln_Blks[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_max = Soln_Blks[i_blk].Grid.KNl-1;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min-1][j_min-1][k_max+1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(6601);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk][l])+x_ref;
              } /* endfor */
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].JCu+1;
	     j_max = Soln_Blks[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Bottomsouth East corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBSE[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
	    Soln_Block_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(6700);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INu+1;
	     i_max = Soln_Blks[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_Blks[i_blk].Grid.JNl-1;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_max = Soln_Blks[i_blk].Grid.KNl-1;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_min-1][j_max+1][k_max+1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(6701);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk][l])+x_ref;
             } /* endfor */
             i_min = Soln_Blks[i_blk].ICu+1;
	     i_max = Soln_Blks[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
                } /* endfor */
          } /* endif */
       } /* endif */


       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for Bottomsouth West corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nBSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level ==
            Soln_Block_List.Block[i_blk].infoBSW[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
	    Soln_Block_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
             j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
             i = Soln_Blks[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc,k_min,k_max,k_inc);
	     if (i != 0) return(6800);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
	     i_min = Soln_Blks[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_Blks[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_Blks[i_blk].Grid.JNl-1;
	     j_inc = 1;
	     k_min = Soln_Blks[i_blk].Grid.KNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     k_max = Soln_Blks[i_blk].Grid.KNl-1;
	     k_inc = 1;
             x_ref = Soln_Blks[i_blk].Grid.Node[i_max+1][j_max+1][k_max+1].X; // Reference node location.
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR3D;
                   if (l >= buffer_size) return(6801);
                   Soln_Blks[i_blk].Grid.Node[i][j][k].X =
		     Vector3D(Soln_Block_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk][l-2],
		              Soln_Block_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk][l])+x_ref;
             } /* endfor */
             i_min = Soln_Blks[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_Blks[i_blk].ICl-1;
	     i_inc = 1;
             j_min = Soln_Blks[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_Blks[i_blk].JCl-1;
	     j_inc = 1;
             k_min = Soln_Blks[i_blk].KCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             k_max = Soln_Blks[i_blk].KCl-1;
             k_inc = 1;
	       for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc )
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc )
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].Xc = Soln_Blks[i_blk].Grid.centroid(i, j,k);
	           Soln_Blks[i_blk].Grid.Cell[i][j][k].V = Soln_Blks[i_blk].Grid.volume(i, j,k);
             } /* endfor */
          } /* endif */
       } /* endif */




    }  /* endfor */

    /* Unloading of receive buffers complete.  Return zero value. */

    return(0);

}

/******************************************************************
 * Routine: Unload_Receive_Message_Buffers_ResChange_FineToCoarse *
 *                                                                *
 * Unloads the receive message passing buffers as sent from more  *
 * refined (higher mesh resolution) solution blocks to more       *
 * coarse neighbouring 3D hexahedrial multi-block solution      *
 * blocks with lower mesh resolution (works on entire 1D array of *
 * 3D solution blocks).                                           *
 *                                                                *
 ******************************************************************/
template <class Hexa_Soln_Block>
int Unload_Receive_Message_Buffers_ResChange_FineToCoarse(Hexa_Soln_Block *Soln_Blks,
                                                          AdaptiveBlock3D_List &Soln_Block_List,
                                                          const int Number_of_Solution_Variables,
                                                          const int Send_Mesh_Geometry_Only,
                                                          const int Send_Conservative_Solution_Fluxes) {

  cout << "\nERROR: Unload_Receive_Message_Buffers_ResChange_FineToCoarse() not defined for 3 dimensions\n";
  return(2);

}

/******************************************************************
 * Routine: Unload_Receive_Message_Buffers_ResChange_CoarseToFine *
 *                                                                *
 * Unloads the receive message passing buffers as sent from more  *
 * coarse (lower mesh resolution) solution blocks to more         *
 * refined neighbouring 3D hexahedrial multi-block solution     *
 * blocks with higher mesh resolution (works on entire 1D array   *
 * of 3D solution blocks).                                        *
 *                                                                *
 ******************************************************************/
template <class Hexa_Soln_Block>
int Unload_Receive_Message_Buffers_ResChange_CoarseToFine(Hexa_Soln_Block *Soln_Blks,
                                                          AdaptiveBlock3D_List &Soln_Block_List,
                                                          const int Number_of_Solution_Variables,
                                                          const int Send_Mesh_Geometry_Only) {

  cout << "\nERROR: Unload_Receive_Message_Buffers_ResChange_CoarseToFine() not defined for 3 dimensions\n";
  return(2);
 
}

/********************************************************
 * Routine: Send_All_Messages                           *
 *                                                      *
 * Loads, sends and exchanges, and then unloads all     *
 * messages for sharing solution information between    *
 * neighbouring 3D hexahedrial multi-block solution   *
 * blocks (works on entire 1D array of 3D solution      *
 * blocks).                                             *
 *                                                      *
 ********************************************************/
template <class Hexa_Soln_Block>
int Send_All_Messages(Hexa_Soln_Block *Soln_Blks,
                      AdaptiveBlock3D_List &Soln_Block_List,
                      const int Number_of_Solution_Variables,
                      const int Send_Mesh_Geometry_Only) {

    int error_flag;

    /* Load message buffers at block interfaces with no cell resolution change. */


    error_flag = Load_Send_Message_Buffers_NoResChange(Soln_Blks,
                                                       Soln_Block_List,
                                                       Number_of_Solution_Variables,
                                                       Send_Mesh_Geometry_Only);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Message Passing Error: Load_Send_Message_Buffers_NoResChange, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */


    /* Exchange message buffers at block interfaces with no cell resolution change. */

    error_flag = AdaptiveBlock3D_List::Exchange_Messages_NoResChange(Soln_Block_List,
                                                                     Number_of_Solution_Variables);
    if (error_flag) {
       cout << "\n " << CFFC_Version()
            << " Message Passing Error: Exchange_Messages_NoResChange. "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */


    /* Unload message buffers at block interfaces with no cell resolution change. */

    error_flag = Unload_Receive_Message_Buffers_NoResChange(Soln_Blks,
                                                            Soln_Block_List,
                                                            Number_of_Solution_Variables,
                                                            Send_Mesh_Geometry_Only);
    if (error_flag) {
       cout << "\n " << CFFC_Version()
            << " Message Passing Error: Unload_Receive_Message_Buffers_NoResChange, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */


    /* Update corner ghost cell information for cases where there are no corner neighbours. */

/*    if (Send_Mesh_Geometry_Only) { */
/*       for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) { */
/*         if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK3D_USED) Update_Corner_Ghost_Nodes(Soln_Blks[nb].Grid); */
/*       } /\* endfor *\/ */
/*     } /\* endif *\/ */

    /* Return error flag. */

    return(error_flag);

}

/********************************************************
 * Routine: Send_Conservative_Flux_Corrections          *
 *                                                      *
 * Loads, sends and exchanges, and then unloads the     *
 * conservative flux corrections for neighbouring 3D    *
 * hexahedrial multi-block solution blocks (works on  *
 * entire 1D array of 3D solution blocks).              *
 *                                                      *
 ********************************************************/
template <class Hexa_Soln_Block>
int Send_Conservative_Flux_Corrections(Hexa_Soln_Block *Soln_Blks,
                                       AdaptiveBlock3D_List &Soln_Block_List,
                                       const int Number_of_Solution_Variables) {

    int error_flag;

    /* Load message buffers at block interfaces for sending from fine to coarse blocks. */

    error_flag = Load_Send_Message_Buffers_ResChange_FineToCoarse(Soln_Blks,
                                                                  Soln_Block_List,
                                                                  Number_of_Solution_Variables,
                                                                  OFF,
                                                                  ON);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Flux Correction Message Passing Error: Load_Send_Message_Buffers_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Exchange message buffers at block interfaces, 
       sending messages from fine to coarse blocks. */

    error_flag = AdaptiveBlock3D_List::Exchange_Messages_ResChange_FineToCoarse(Soln_Block_List,
                                                                                Number_of_Solution_Variables);
   if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Flux Correction Message Passing Error: Exchange_Messages_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Unload message buffers at block interfaces as sent from fine to coarse blocks. */

    error_flag = Unload_Receive_Message_Buffers_ResChange_FineToCoarse(Soln_Blks,
                                                                       Soln_Block_List,
                                                                       Number_of_Solution_Variables,
                                                                       OFF,
                                                                       ON);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Flux Correction Message Passing Error: Unload_Receive_Message_Buffers_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
    } /* endif */
    return(error_flag);

}

#endif // _ADAPTIVEBLOCK3D_MESSAGEPASSING_INCLUDED
