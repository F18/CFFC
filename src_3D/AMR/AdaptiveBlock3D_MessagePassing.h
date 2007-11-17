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

#ifndef  _SYSTEM_LINUX_INCLUDED
#include "../System/System_Linux.h"
#endif //_SYSTEM_LINUX_INCLUDED


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
 
   int i_bound_elem; // index for boundary element, face edge or vertex
   int n_bound_elem[27];
   int n_imin, n_imax, n_jmin, n_jmax, n_kmin, n_kmax;
   int i_temp;
   
   AdaptiveBlock3D_Info info_bound_elem[27];
   int compact_trans_matrix[3]; // orientation of the block that receives message.
   int ti[3]; // ti[0] -i, ti[1] -j, ti[2] - k; computational orientation of the block that receives message.
   int ts[3]; // sign of the computational orientation of the block that receives message.
   int id_min[3], id_max[3], inc[3];
 
  /* Check to see if the number of solution variables specified
     in the calling argument is correct. */

  if (Number_of_Solution_Variables != NUM_COMP_VECTOR3D &&
      Number_of_Solution_Variables != Soln_Blks[0].NumVar() &&
      Number_of_Solution_Variables != Soln_Blks[0].NumVar() + NUM_COMP_VECTOR3D) {
    return(1);
  } /* endif */
   
    /* Load the send buffers of each solution block. */
  for ( i_blk = 0 ; i_blk <= Soln_Block_List.Nblk-1 ; ++i_blk ) {
     
  
     if(Soln_Block_List.Block[i_blk].used){
        // assign the boudnary element info
        n_bound_elem[BE::BSW] = Soln_Block_List.Block[i_blk].nBSW;
        info_bound_elem[BE::BSW] =  Soln_Block_List.Block[i_blk].infoBSW[0];
        
        n_bound_elem[BE::SW] = Soln_Block_List.Block[i_blk].nSW;
        info_bound_elem[BE::SW] =  Soln_Block_List.Block[i_blk].infoSW[0];
        
        n_bound_elem[BE::TSW] = Soln_Block_List.Block[i_blk].nTSW;
        info_bound_elem[BE::TSW] =  Soln_Block_List.Block[i_blk].infoTSW[0];
        
        n_bound_elem[BE::BW] = Soln_Block_List.Block[i_blk].nBW;
        info_bound_elem[BE::BW] =  Soln_Block_List.Block[i_blk].infoBW[0];
        
        n_bound_elem[BE::W] = Soln_Block_List.Block[i_blk].nW;
        info_bound_elem[BE::W] =  Soln_Block_List.Block[i_blk].infoW[0];
        
        n_bound_elem[BE::TW] = Soln_Block_List.Block[i_blk].nTW;
        info_bound_elem[BE::TW] =  Soln_Block_List.Block[i_blk].infoTW[0];
        
        n_bound_elem[BE::BNW] = Soln_Block_List.Block[i_blk].nBNW;
        info_bound_elem[BE::BNW] =  Soln_Block_List.Block[i_blk].infoBNW[0];
        
        n_bound_elem[BE::NW] = Soln_Block_List.Block[i_blk].nNW;
        info_bound_elem[BE::NW] =  Soln_Block_List.Block[i_blk].infoNW[0];
        
        n_bound_elem[BE::TNW] = Soln_Block_List.Block[i_blk].nTNW;
        info_bound_elem[BE::TNW] =  Soln_Block_List.Block[i_blk].infoTNW[0];
        
        n_bound_elem[BE::BS] = Soln_Block_List.Block[i_blk].nBS;
        info_bound_elem[BE::BS] =  Soln_Block_List.Block[i_blk].infoBS[0];
        
        n_bound_elem[BE::S] = Soln_Block_List.Block[i_blk].nS;
        info_bound_elem[BE::S] =  Soln_Block_List.Block[i_blk].infoS[0];
        
        n_bound_elem[BE::TS] = Soln_Block_List.Block[i_blk].nTS;
        info_bound_elem[BE::TS] =  Soln_Block_List.Block[i_blk].infoTS[0];
        
        n_bound_elem[BE::B] = Soln_Block_List.Block[i_blk].nB;
        info_bound_elem[BE::B] =  Soln_Block_List.Block[i_blk].infoB[0];
        
        n_bound_elem[BE::T] = Soln_Block_List.Block[i_blk].nT;
        info_bound_elem[BE::T] =  Soln_Block_List.Block[i_blk].infoT[0];
        
        n_bound_elem[BE::BN] = Soln_Block_List.Block[i_blk].nBN;
        info_bound_elem[BE::BN] =  Soln_Block_List.Block[i_blk].infoBN[0];
        
        n_bound_elem[BE::N] = Soln_Block_List.Block[i_blk].nN;
        info_bound_elem[BE::N] =  Soln_Block_List.Block[i_blk].infoN[0];
        
        n_bound_elem[BE::TN] = Soln_Block_List.Block[i_blk].nTN;
        info_bound_elem[BE::TN] =  Soln_Block_List.Block[i_blk].infoTN[0];
        
        n_bound_elem[BE::BSE] = Soln_Block_List.Block[i_blk].nBSE;
        info_bound_elem[BE::BSE] =  Soln_Block_List.Block[i_blk].infoBSE[0];
        
        n_bound_elem[BE::SE] = Soln_Block_List.Block[i_blk].nSE;
        info_bound_elem[BE::SE] =  Soln_Block_List.Block[i_blk].infoSE[0];
        
        n_bound_elem[BE::TSE] = Soln_Block_List.Block[i_blk].nTSE;
        info_bound_elem[BE::TSE] =  Soln_Block_List.Block[i_blk].infoTSE[0];
        
        n_bound_elem[BE::BE] = Soln_Block_List.Block[i_blk].nBE;
        info_bound_elem[BE::BE] =  Soln_Block_List.Block[i_blk].infoBE[0];
        
        n_bound_elem[BE::E] = Soln_Block_List.Block[i_blk].nE;
        info_bound_elem[BE::E] =  Soln_Block_List.Block[i_blk].infoE[0];
        
        n_bound_elem[BE::TE] = Soln_Block_List.Block[i_blk].nTE;
        info_bound_elem[BE::TE] =  Soln_Block_List.Block[i_blk].infoTE[0]; 
        
        n_bound_elem[BE::BNE] = Soln_Block_List.Block[i_blk].nBNE;
        info_bound_elem[BE::BNE] =  Soln_Block_List.Block[i_blk].infoBNE[0]; 
        
        n_bound_elem[BE::NE] = Soln_Block_List.Block[i_blk].nNE;
        info_bound_elem[BE::NE] =  Soln_Block_List.Block[i_blk].infoNE[0];  
        
        n_bound_elem[BE::TNE] = Soln_Block_List.Block[i_blk].nTNE;
        info_bound_elem[BE::TNE] =  Soln_Block_List.Block[i_blk].infoTNE[0]; 

     }
     
     for (int ii = -1; ii<2; ii++){
        for (int jj = -1; jj<2; jj++){
           for (int kk = -1; kk<2; kk++){
              i_bound_elem = 9*(ii+1) + 3*(jj+1) + (kk+1);
              
              if (Soln_Block_List.Block[i_blk].used  && (n_bound_elem[i_bound_elem] == 1) && (i_bound_elem != 13) &&
                  (Soln_Block_List.Block[i_blk].info.level == info_bound_elem[i_bound_elem].level)) {
                 if (!Send_Mesh_Geometry_Only) {
                    buffer_size_neighbour = ((abs(ii)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!ii)*abs(Soln_Block_List.Block[i_blk].info.dimen.i)))*
                       ((abs(jj)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!jj)*abs(Soln_Block_List.Block[i_blk].info.dimen.j)))*
                       ((abs(kk)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!kk)*abs(Soln_Block_List.Block[i_blk].info.dimen.k)))*(Number_of_Solution_Variables);
                    
                    if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                       buffer_size_neighbour = buffer_size_neighbour +
                          ((abs(ii)*Soln_Block_List.Block[i_blk].info.dimen.ghost)+ ((!ii)*abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1))*
                          ((abs(jj)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!jj)*abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1))*
                          ((abs(kk)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!kk)*abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1))*(NUM_COMP_VECTOR3D);
                    
                 } else {
                    buffer_size_neighbour = ((abs(ii)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!ii)*abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1))*
                       ((abs(jj)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!jj)*abs(Soln_Block_List.Block[i_blk].info.dimen.j) +1))*
                       ((abs(kk)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!kk)*abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1))*(NUM_COMP_VECTOR3D)+
                         (((!ii)*abs(Soln_Block_List.Block[i_blk].info.dimen.i)))*(((!jj)*abs(Soln_Block_List.Block[i_blk].info.dimen.j)))*
                         (((!kk)*abs(Soln_Block_List.Block[i_blk].info.dimen.k)));
                    
                    
                 } /* endif */
                 l = -1;
                 
                 // Load ghost cell solution variable information as required.
                 if (!Send_Mesh_Geometry_Only) {
                    if( ii == -1 ){
                       n_imin = Soln_Blks[i_blk].Nghost;
                       n_imax = Soln_Blks[i_blk].ICl +1;
                       
                    }else if( ii == 1){
                       n_imin =  Soln_Blks[i_blk].ICu-1;
                       n_imax = Soln_Blks[i_blk].ICu;
                       
                    }else{
                       n_imin = 0;
                       n_imax = 0;
                    }
                    if( jj == -1 ){
                       n_jmin = Soln_Blks[i_blk].Nghost;
                       n_jmax = Soln_Blks[i_blk].JCl +1;
                       
                    }else if( jj == 1){
                       n_jmin =  Soln_Blks[i_blk].JCu-1;
                       n_jmax = Soln_Blks[i_blk].JCu;
                       
                    }else{
                       n_jmin = 0;
                       n_jmax = 0;
                    }
                    if( kk == -1 ){
                       n_kmin = Soln_Blks[i_blk].Nghost;
                       n_kmax = Soln_Blks[i_blk].KCl +1;
                       
                    }else if( kk == 1){
                       n_kmin =  Soln_Blks[i_blk].KCu-1;
                       n_kmax = Soln_Blks[i_blk].KCu;
                       
                    }else{
                       n_kmin = 0;
                       n_kmax = 0;
                    }
   
                    i_min = (!ii)*Soln_Blks[i_blk].ICl + n_imin;
                    i_max = (!ii)*Soln_Blks[i_blk].ICu + n_imax;
                    i_inc = 1;
                    j_min = (!jj)*Soln_Blks[i_blk].JCl + n_jmin;
                    j_max = (!jj)*Soln_Blks[i_blk].JCu + n_jmax;
                    j_inc = 1;
                    k_min = (!kk)*Soln_Blks[i_blk].KCl + n_kmin;
                    k_max = (!kk)*Soln_Blks[i_blk].KCu + n_kmax;
                    k_inc = 1;
                    
                    // Using transformation matrices to load the send buffer in the order of unloading.
                    compact_trans_matrix[0] = info_bound_elem[i_bound_elem].blkorient.ctm_offsets[0];
                    compact_trans_matrix[1] = info_bound_elem[i_bound_elem].blkorient.ctm_offsets[1];
                    compact_trans_matrix[2] = info_bound_elem[i_bound_elem].blkorient.ctm_offsets[2];
                   
                    ts[0] = sgn(compact_trans_matrix[0]);
                    ts[1] = sgn(compact_trans_matrix[1]);
                    ts[2] = sgn(compact_trans_matrix[2]);
                    
                    ti[0] = abs(compact_trans_matrix[0]) -1;
                    ti[1] = abs(compact_trans_matrix[1]) -1;
                    ti[2] = abs(compact_trans_matrix[2]) -1;

                    if(ts[0] <0 ){
                       i_inc = -1;
                       i_temp = i_min;
                       i_min = i_max;
                       i_max = i_temp;
                    }
                    if(ts[1] <0 ){
                       j_inc = -1;
                       i_temp = j_min;
                       j_min = j_max;
                       j_max = i_temp;
                    }
                    if(ts[2] <0 ){
                       k_inc = -1;
                       i_temp = k_min;
                       k_min = k_max;
                       k_max = i_temp;
                    }
                    
                    id_min[0] = i_min;
                    id_min[1] = j_min;
                    id_min[2] = k_min;
                    id_max[0] = i_max;
                    id_max[1] = j_max;
                    id_max[2] = k_max;
                    inc[0] = i_inc;
                    inc[1] = j_inc;
                    inc[2] = k_inc;
                    
                    i = Soln_Blks[i_blk].LoadSendBuffer_Solution(Soln_Block_List.message_noreschange_sendbuf[i_blk][i_bound_elem],
                                                        l,buffer_size_neighbour, id_min, id_max, inc, ti);
                     
                    if (i != 0) return(2200);
                 } /* endif */
                 // Load ghost cell mesh information as required.
                 if (Send_Mesh_Geometry_Only || Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
                    
                    if( ii == -1 ){
                       n_imin = Soln_Blks[i_blk].Grid.INl +1;
                       n_imax = Soln_Blks[i_blk].Grid.INl + Soln_Blks[i_blk].Nghost;
                       
                    }else if( ii == 1){
                       n_imin =  Soln_Blks[i_blk].Grid.INu-Soln_Blks[i_blk].Nghost;
                       n_imax = Soln_Blks[i_blk].Grid.INu-1;
                       
                    }else{
                       n_imin = 0;
                       n_imax = 0;
                    }
                    if( jj == -1 ){
                       n_jmin = Soln_Blks[i_blk].Grid.JNl +1;
                       n_jmax = Soln_Blks[i_blk].Grid.JNl +Soln_Blks[i_blk].Nghost;
                       
                    }else if( jj == 1){
                       n_jmin =  Soln_Blks[i_blk].Grid.JNu-Soln_Blks[i_blk].Nghost;
                       n_jmax = Soln_Blks[i_blk].Grid.JNu-1;
                       
                    }else{
                       n_jmin = 0;
                       n_jmax = 0;
                    }
                    if( kk == -1 ){
                       n_kmin = Soln_Blks[i_blk].Grid.KNl +1;
                       n_kmax = Soln_Blks[i_blk].Grid.KNl +Soln_Blks[i_blk].Nghost;
                       
                    }else if( kk == 1){
                       n_kmin =  Soln_Blks[i_blk].Grid.KNu-Soln_Blks[i_blk].Nghost;
                       n_kmax = Soln_Blks[i_blk].Grid.KNu -1;
                       
                    }else{
                       n_kmin = 0;
                       n_kmax = 0;
                    }
                    
                  
                    i_min = (!ii)*Soln_Blks[i_blk].Grid.INl + n_imin;
                    i_max = (!ii)*Soln_Blks[i_blk].Grid.INu + n_imax;
                    i_inc = 1;
                    j_min = (!jj)*Soln_Blks[i_blk].Grid.JNl + n_jmin;
                    j_max = (!jj)*Soln_Blks[i_blk].Grid.JNu + n_jmax;
                    j_inc = 1;
                    k_min = (!kk)*Soln_Blks[i_blk].Grid.KNl + n_kmin;
                    k_max = (!kk)*Soln_Blks[i_blk].Grid.KNu + n_kmax;
                    k_inc = 1;
           
                    // Using transformation matrices to load the send buffer in the order of unloading.
                    compact_trans_matrix[0] = info_bound_elem[i_bound_elem].blkorient.ctm_offsets[0];
                    compact_trans_matrix[1] = info_bound_elem[i_bound_elem].blkorient.ctm_offsets[1];
                    compact_trans_matrix[2] = info_bound_elem[i_bound_elem].blkorient.ctm_offsets[2];
                   
                    ts[0] = sgn(compact_trans_matrix[0]);
                    ts[1] = sgn(compact_trans_matrix[1]);
                    ts[2] = sgn(compact_trans_matrix[2]);
                    
                    ti[0] = abs(compact_trans_matrix[0]) -1;
                    ti[1] = abs(compact_trans_matrix[1]) -1;
                    ti[2] = abs(compact_trans_matrix[2]) -1;

                    if(ts[0] <0 ){
                       i_inc = -1;
                       i_temp = i_min;
                       i_min = i_max;
                       i_max = i_temp;
                    }
                    if(ts[1] <0 ){
                       j_inc = -1;
                       i_temp = j_min;
                       j_min = j_max;
                       j_max = i_temp;
                    }
                    if(ts[2] <0 ){
                       k_inc = -1;
                       i_temp = k_min;
                       k_min = k_max;
                       k_max = i_temp;
                    }
                    
                    id_min[0] = i_min;
                    id_min[1] = j_min;
                    id_min[2] = k_min;
                    id_max[0] = i_max;
                    id_max[1] = j_max;
                    id_max[2] = k_max;
                    inc[0] = i_inc;
                    inc[1] = j_inc;
                    inc[2] = k_inc;
              
                    i = Soln_Blks[i_blk].LoadSendBuffer_Geometry(
                       Soln_Block_List.message_noreschange_sendbuf[i_blk][i_bound_elem],
                       l,buffer_size_neighbour, id_min, id_max, inc, ti);
    
   
                    if (i != 0) return(2201);
                    
                    if( ii == -1 ){
                       n_imin = Soln_Blks[i_blk].Nghost;
                       n_imax = Soln_Blks[i_blk].ICl +1;
                       
                    }else if( ii == 1){
                       n_imin =  Soln_Blks[i_blk].ICu-1;
                       n_imax = Soln_Blks[i_blk].ICu;
                       
                    }else{
                       n_imin = 0;
                       n_imax = 0;
                    }
                    if( jj == -1 ){
                       n_jmin = Soln_Blks[i_blk].Nghost;
                       n_jmax = Soln_Blks[i_blk].JCl +1;
                       
                    }else if( jj == 1){
                       n_jmin =  Soln_Blks[i_blk].JCu-1;
                       n_jmax = Soln_Blks[i_blk].JCu;
                       
                    }else{
                       n_jmin = 0;
                       n_jmax = 0;
                    }
                    if( kk == -1 ){
                       n_kmin = Soln_Blks[i_blk].Nghost;
                       n_kmax = Soln_Blks[i_blk].KCl +1;
                       
                    }else if( kk == 1){
                       n_kmin =  Soln_Blks[i_blk].KCu-1;
                       n_kmax = Soln_Blks[i_blk].KCu;
                       
                    }else{
                       n_kmin = 0;
                       n_kmax = 0;
                    }
                    
                    i_min = (!ii)*Soln_Blks[i_blk].ICl + n_imin;
                    i_max = (!ii)*Soln_Blks[i_blk].ICu + n_imax;
                    i_inc = 1;
                    j_min = (!jj)*Soln_Blks[i_blk].JCl + n_jmin;
                    j_max = (!jj)*Soln_Blks[i_blk].JCu + n_jmax;
                    j_inc = 1;
                    k_min = (!kk)*Soln_Blks[i_blk].KCl + n_kmin;
                    k_max = (!kk)*Soln_Blks[i_blk].KCu + n_kmax;
                    k_inc = 1;
           
                    // Using transformation matrices to load the send buffer in the order of unloading.
                    compact_trans_matrix[0] = info_bound_elem[i_bound_elem].blkorient.ctm_offsets[0];
                    compact_trans_matrix[1] = info_bound_elem[i_bound_elem].blkorient.ctm_offsets[1];
                    compact_trans_matrix[2] = info_bound_elem[i_bound_elem].blkorient.ctm_offsets[2];
                    
                    ts[0] = sgn(compact_trans_matrix[0]);
                    ts[1] = sgn(compact_trans_matrix[1]);
                    ts[2] = sgn(compact_trans_matrix[2]);
                    
                    ti[0] = abs(compact_trans_matrix[0]) -1;
                    ti[1] = abs(compact_trans_matrix[1]) -1;
                    ti[2] = abs(compact_trans_matrix[2]) -1;

                    if(ts[0] <0 ){
                       i_inc = -1;
                       i_temp = i_min;
                       i_min = i_max;
                       i_max = i_temp;
                    }
                    if(ts[1] <0 ){
                       j_inc = -1;
                       i_temp = j_min;
                       j_min = j_max;
                       j_max = i_temp;
                    }
                    if(ts[2] <0 ){
                       k_inc = -1;
                       i_temp = k_min;
                       k_min = k_max;
                       k_max = i_temp;
                    }
                    
                    id_min[0] = i_min;
                    id_min[1] = j_min;
                    id_min[2] = k_min;
                    id_max[0] = i_max;
                    id_max[1] = j_max;
                    id_max[2] = k_max;
                    inc[0] = i_inc;
                    inc[1] = j_inc;
                    inc[2] = k_inc;

/*                     i = Soln_Blks[i_blk].LoadSendBuffer_BCs( Soln_Block_List.message_noreschange_sendbuf[i_blk][i_bound_elem], */
/*                                                              l,buffer_size_neighbour, id_min, id_max, inc, ti, ii, jj, kk); */
  
 
                 } /* endif */
              } /* endif */
           }/* end of for k*/
        }/* end of for j*/
     }/* end of for i*/
  

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

    int i_bound_elem; // index for boundary element, face edge or vertex
    int n_bound_elem[27];
    AdaptiveBlock3D_Info info_bound_elem[27];
    int n_imin, n_imax, n_jmin, n_jmax, n_kmin, n_kmax;
    int recv_bound_elem, recv_blknum;
    

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != NUM_COMP_VECTOR3D &&
        Number_of_Solution_Variables != Soln_Blks[0].NumVar() &&
        Number_of_Solution_Variables != Soln_Blks[0].NumVar() + NUM_COMP_VECTOR3D) {
      return(1);
    } /* endif */

    /* Unload the receive buffers for each solution block. */

    for ( i_blk = 0 ; i_blk <= Soln_Block_List.Nblk-1 ; ++i_blk ) {

       if (Soln_Block_List.Block[i_blk].used ){
          
          // assign the boudnary element info
          n_bound_elem[BE::BSW] = Soln_Block_List.Block[i_blk].nBSW;
          info_bound_elem[BE::BSW] =  Soln_Block_List.Block[i_blk].infoBSW[0];
          
          n_bound_elem[BE::SW] = Soln_Block_List.Block[i_blk].nSW;
          info_bound_elem[BE::SW] =  Soln_Block_List.Block[i_blk].infoSW[0];
          
          n_bound_elem[BE::TSW] = Soln_Block_List.Block[i_blk].nTSW;
          info_bound_elem[BE::TSW] =  Soln_Block_List.Block[i_blk].infoTSW[0];
          
          n_bound_elem[BE::BW] = Soln_Block_List.Block[i_blk].nBW;
          info_bound_elem[BE::BW] =  Soln_Block_List.Block[i_blk].infoBW[0];
          
          n_bound_elem[BE::W] = Soln_Block_List.Block[i_blk].nW;
          info_bound_elem[BE::W] =  Soln_Block_List.Block[i_blk].infoW[0];
          
          n_bound_elem[BE::TW] = Soln_Block_List.Block[i_blk].nTW;
          info_bound_elem[BE::TW] =  Soln_Block_List.Block[i_blk].infoTW[0];
          
          n_bound_elem[BE::BNW] = Soln_Block_List.Block[i_blk].nBNW;
          info_bound_elem[BE::BNW] =  Soln_Block_List.Block[i_blk].infoBNW[0];
          
          n_bound_elem[BE::NW] = Soln_Block_List.Block[i_blk].nNW;
          info_bound_elem[BE::NW] =  Soln_Block_List.Block[i_blk].infoNW[0];
          
          n_bound_elem[BE::TNW] = Soln_Block_List.Block[i_blk].nTNW;
          info_bound_elem[BE::TNW] =  Soln_Block_List.Block[i_blk].infoTNW[0];
          
          n_bound_elem[BE::BS] = Soln_Block_List.Block[i_blk].nBS;
          info_bound_elem[BE::BS] =  Soln_Block_List.Block[i_blk].infoBS[0];
          
          n_bound_elem[BE::S] = Soln_Block_List.Block[i_blk].nS;
          info_bound_elem[BE::S] =  Soln_Block_List.Block[i_blk].infoS[0];
          
          n_bound_elem[BE::TS] = Soln_Block_List.Block[i_blk].nTS;
          info_bound_elem[BE::TS] =  Soln_Block_List.Block[i_blk].infoTS[0];
          
          n_bound_elem[BE::B] = Soln_Block_List.Block[i_blk].nB;
          info_bound_elem[BE::B] =  Soln_Block_List.Block[i_blk].infoB[0];
          
          n_bound_elem[BE::T] = Soln_Block_List.Block[i_blk].nT;
          info_bound_elem[BE::T] =  Soln_Block_List.Block[i_blk].infoT[0];
          
          n_bound_elem[BE::BN] = Soln_Block_List.Block[i_blk].nBN;
          info_bound_elem[BE::BN] =  Soln_Block_List.Block[i_blk].infoBN[0];
          
          n_bound_elem[BE::N] = Soln_Block_List.Block[i_blk].nN;
          info_bound_elem[BE::N] =  Soln_Block_List.Block[i_blk].infoN[0];
          
          n_bound_elem[BE::TN] = Soln_Block_List.Block[i_blk].nTN;
          info_bound_elem[BE::TN] =  Soln_Block_List.Block[i_blk].infoTN[0];
          
          n_bound_elem[BE::BSE] = Soln_Block_List.Block[i_blk].nBSE;
          info_bound_elem[BE::BSE] =  Soln_Block_List.Block[i_blk].infoBSE[0];
          
          n_bound_elem[BE::SE] = Soln_Block_List.Block[i_blk].nSE;
          info_bound_elem[BE::SE] =  Soln_Block_List.Block[i_blk].infoSE[0];
          
          n_bound_elem[BE::TSE] = Soln_Block_List.Block[i_blk].nTSE;
          info_bound_elem[BE::TSE] =  Soln_Block_List.Block[i_blk].infoTSE[0];
          
          n_bound_elem[BE::BE] = Soln_Block_List.Block[i_blk].nBE;
          info_bound_elem[BE::BE] =  Soln_Block_List.Block[i_blk].infoBE[0];
          
          n_bound_elem[BE::E] = Soln_Block_List.Block[i_blk].nE;
          info_bound_elem[BE::E] =  Soln_Block_List.Block[i_blk].infoE[0];

          n_bound_elem[BE::TE] = Soln_Block_List.Block[i_blk].nTE;
          info_bound_elem[BE::TE] =  Soln_Block_List.Block[i_blk].infoTE[0]; 
          
          n_bound_elem[BE::BNE] = Soln_Block_List.Block[i_blk].nBNE;
          info_bound_elem[BE::BNE] =  Soln_Block_List.Block[i_blk].infoBNE[0]; 
          
          n_bound_elem[BE::NE] = Soln_Block_List.Block[i_blk].nNE;
          info_bound_elem[BE::NE] =  Soln_Block_List.Block[i_blk].infoNE[0];  
          
          n_bound_elem[BE::TNE] = Soln_Block_List.Block[i_blk].nTNE;
          info_bound_elem[BE::TNE] =  Soln_Block_List.Block[i_blk].infoTNE[0]; 
          
       }
       
       for (int ii = -1; ii<2; ii++){
          for (int jj = -1; jj<2; jj++){
             for (int kk = -1; kk<2; kk++){
                
                i_bound_elem = 9*(ii+1) + 3*(jj+1) + (kk+1);
                
                if (Soln_Block_List.Block[i_blk].used  && (n_bound_elem[i_bound_elem] == 1) && (i_bound_elem != 13) &&
                    (Soln_Block_List.Block[i_blk].info.level == info_bound_elem[i_bound_elem].level)) {
                   
                   if (!Send_Mesh_Geometry_Only) {
                      
                      buffer_size =  ((abs(ii)*Soln_Block_List.Block[i_blk].info.dimen.ghost)+ ((!ii)*abs(Soln_Block_List.Block[i_blk].info.dimen.i)))*
                         ((abs(jj)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!jj)*abs(Soln_Block_List.Block[i_blk].info.dimen.j)))*
                         ((abs(kk)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!kk)*abs(Soln_Block_List.Block[i_blk].info.dimen.k)))*( Soln_Blks[i_blk].NumVar());
                      
                      if (Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar())
                         
                         buffer_size = buffer_size + 
                            ((abs(ii)*Soln_Block_List.Block[i_blk].info.dimen.ghost)+ ((!ii)*abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1))*
                            ((abs(jj)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!jj)*abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1))*
                            ((abs(kk)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!kk)*abs(Soln_Block_List.Block[i_blk].info.dimen.k)+1))*(NUM_COMP_VECTOR3D);
                      
                   } else {
                      buffer_size = ((abs(ii)*Soln_Block_List.Block[i_blk].info.dimen.ghost)+ ((!ii)*abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1))*
                         ((abs(jj)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!jj)*abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1))*
                         ((abs(kk)*Soln_Block_List.Block[i_blk].info.dimen.ghost) + ((!kk)*abs(Soln_Block_List.Block[i_blk].info.dimen.k) +1))*(NUM_COMP_VECTOR3D) +
                         (((!ii)*abs(Soln_Block_List.Block[i_blk].info.dimen.i)))*(((!jj)*abs(Soln_Block_List.Block[i_blk].info.dimen.j)))*
                         (((!kk)*abs(Soln_Block_List.Block[i_blk].info.dimen.k)));
                      
                      
                   } /* endif */
                   l = -1;
                   
                   // Unload ghost cell solution information as required.
                   if (!Send_Mesh_Geometry_Only) {
                      
                      if( ii == -1 ){
                         n_imin = Soln_Blks[i_blk].Nghost-2;
                         n_imax = Soln_Blks[i_blk].Nghost-1;
                         
                      }else if( ii == 1){
                         n_imin =  Soln_Blks[i_blk].ICu+ 1;
                         n_imax = Soln_Blks[i_blk].ICu+ Soln_Blks[i_blk].Nghost;
                         
                      }else{
                         n_imin = 0;
                         n_imax = 0;
                      }
                      if( jj == -1 ){
                         n_jmin = Soln_Blks[i_blk].Nghost-2;
                         n_jmax = Soln_Blks[i_blk].Nghost-1;
                       
                      }else if( jj == 1){
                         n_jmin =  Soln_Blks[i_blk].JCu+1;
                         n_jmax = Soln_Blks[i_blk].JCu + Soln_Blks[i_blk].Nghost;
                         
                      }else{
                         n_jmin = 0;
                         n_jmax = 0;
                      }
                      if( kk == -1 ){
                         n_kmin = Soln_Blks[i_blk].Nghost-2;
                         n_kmax = Soln_Blks[i_blk].Nghost-1;
                         
                      }else if( kk == 1){
                         n_kmin =  Soln_Blks[i_blk].KCu + 1;
                         n_kmax = Soln_Blks[i_blk].KCu +  Soln_Blks[i_blk].Nghost;
                         
                      }else{
                         n_kmin = 0;
                         n_kmax = 0;
                      }
   
                      i_min = (!ii)*Soln_Blks[i_blk].ICl + n_imin;
                      i_max = (!ii)*Soln_Blks[i_blk].ICu + n_imax;
                      i_inc = 1;
                      j_min = (!jj)*Soln_Blks[i_blk].JCl + n_jmin;
                      j_max = (!jj)*Soln_Blks[i_blk].JCu + n_jmax;
                      j_inc = 1;
                      k_min = (!kk)*Soln_Blks[i_blk].KCl + n_kmin;
                      k_max = (!kk)*Soln_Blks[i_blk].KCu + n_kmax;
                      k_inc = 1;
                    
                      i = Soln_Blks[i_blk].UnloadReceiveBuffer_Solution(Soln_Block_List.message_noreschange_recbuf[i_blk][i_bound_elem],
                                                               l,buffer_size,
                                                               i_min,i_max,i_inc,j_min,j_max,j_inc,
                                                               k_min,k_max,k_inc);
                    
        

                 
                    if (i != 0) return(2202);
                   } /* endif */
                 // Unload ghost cell mesh information as required.
                 if (Send_Mesh_Geometry_Only ||
                     Number_of_Solution_Variables > Soln_Blks[i_blk].NumVar()) {
                    if( ii == -1 ){
                       n_imin = Soln_Blks[i_blk].Nghost-2;
                       n_imax = Soln_Blks[i_blk].Nghost-1;
                       
                    }else if( ii == 1){
                       n_imin =  Soln_Blks[i_blk].Grid.INu+1;
                       n_imax = Soln_Blks[i_blk].Grid.INu+Soln_Blks[i_blk].Nghost;
                       
                    }else{
                       n_imin = 0;
                       n_imax = 0;
                    }
                    if( jj == -1 ){
                       n_jmin = Soln_Blks[i_blk].Nghost - 2;
                       n_jmax = Soln_Blks[i_blk].Nghost - 1;
                       
                    }else if( jj == 1){
                       n_jmin =  Soln_Blks[i_blk].Grid.JNu+1;
                       n_jmax = Soln_Blks[i_blk].Grid.JNu + Soln_Blks[i_blk].Nghost;
                       
                    }else{
                       n_jmin = 0;
                       n_jmax = 0;
                    }
                    if( kk == -1 ){
                       n_kmin = Soln_Blks[i_blk].Nghost - 2;
                       n_kmax = Soln_Blks[i_blk].Nghost - 1;
                       
                    }else if( kk == 1){
                       n_kmin =  Soln_Blks[i_blk].Grid.KNu+1;
                       n_kmax = Soln_Blks[i_blk].Grid.KNu+ Soln_Blks[i_blk].Nghost;
                       
                    }else{
                       n_kmin = 0;
                       n_kmax = 0;
                    }
                    
                    i_min = (!ii)*Soln_Blks[i_blk].Grid.INl + n_imin;
                    i_max = (!ii)*Soln_Blks[i_blk].Grid.INu + n_imax;
                    i_inc = 1;
                    j_min = (!jj)*Soln_Blks[i_blk].Grid.JNl + n_jmin;
                    j_max = (!jj)*Soln_Blks[i_blk].Grid.JNu + n_jmax;
                    j_inc = 1;
                    k_min = (!kk)*Soln_Blks[i_blk].Grid.KNl + n_kmin;
                    k_max = (!kk)*Soln_Blks[i_blk].Grid.KNu + n_kmax;
                    k_inc = 1;
                            


                 
               
                      i = Soln_Blks[i_blk].UnloadReceiveBuffer_Geometry(
                         Soln_Block_List.message_noreschange_recbuf[i_blk][i_bound_elem],
                         l,buffer_size,
                         i_min,i_max,i_inc,j_min,j_max,j_inc,
                         k_min,k_max, k_inc);
                      
                      
                    
                    if (i != 0) return(2203);
                    
                    if( ii == -1 ){
                       n_imin = Soln_Blks[i_blk].Nghost -2 ;
                       n_imax = Soln_Blks[i_blk].Nghost -1 ;
                       
                    }else if( ii == 1){
                       n_imin = Soln_Blks[i_blk].ICu+1;
                       n_imax = Soln_Blks[i_blk].ICu + Soln_Blks[i_blk].Nghost;
                       
                    }else{
                       n_imin = 0;
                       n_imax = 0;
                    }
                    if( jj == -1 ){
                       n_jmin = Soln_Blks[i_blk].Nghost - 2;
                       n_jmax = Soln_Blks[i_blk].Nghost - 1;
                       
                    }else if( jj == 1){
                       n_jmin =  Soln_Blks[i_blk].JCu+1;
                       n_jmax = Soln_Blks[i_blk].JCu+ Soln_Blks[i_blk].Nghost;
                       
                    }else{
                       n_jmin = 0;
                       n_jmax = 0;
                    }
                    if( kk == -1 ){
                       n_kmin = Soln_Blks[i_blk].Nghost - 2;
                       n_kmax =  Soln_Blks[i_blk].Nghost - 1;
                       
                    }else if( kk == 1){
                       n_kmin =  Soln_Blks[i_blk].KCu+1;
                       n_kmax = Soln_Blks[i_blk].KCu+ Soln_Blks[i_blk].Nghost;
                       
                    }else{
                       n_kmin = 0;
                       n_kmax = 0;
                    }
                    
                    i_min = (!ii)*Soln_Blks[i_blk].ICl + n_imin;
                    i_max = (!ii)*Soln_Blks[i_blk].ICu + n_imax;
                    i_inc = 1;
                    j_min = (!jj)*Soln_Blks[i_blk].JCl + n_jmin;
                    j_max = (!jj)*Soln_Blks[i_blk].JCu + n_jmax;
                    j_inc = 1;
                    k_min = (!kk)*Soln_Blks[i_blk].KCl + n_kmin;
                    k_max = (!kk)*Soln_Blks[i_blk].KCu + n_kmax;
                    k_inc = 1;
                    
                            
/*                     i = Soln_Blks[i_blk].UnloadReceiveBuffer_BCs( */
/*                        Soln_Block_List.message_noreschange_recbuf[i_blk][i_bound_elem], */
/*                        l,buffer_size, i_min,i_max,i_inc,j_min,j_max,j_inc, */
/*                        k_min,k_max,k_inc, ii, jj, kk); */
                    
        



                 } /* endif */
              } /* endif */
           }  /* endfor */ 
        }
     }
    }
 
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

  
    
/*     for ( int iProc = 0; iProc !=  CFFC_MPI::Number_of_Processors; ++iProc ) { */
/*        if (  CFFC_MPI::This_Processor_Number == iProc ) { */
/*           cout<<"\n  after loading buffer function received location "<<endl; */
/*           cout<<"\n CFFC_MPI::This_Processor_Number = "<< CFFC_MPI::This_Processor_Number<<endl; */
     
      
/*           System::sleep(0.1); */
/*        } */
/*        MPI::COMM_WORLD.Barrier(); */
/*     } */

    /* Exchange message buffers at block interfaces with no cell resolution change. */

    error_flag = AdaptiveBlock3D_List::Exchange_Messages_NoResChange(Soln_Block_List,
                                                                     Number_of_Solution_Variables,
                                                                     Send_Mesh_Geometry_Only);
    if (error_flag) {
       cout << "\n " << CFFC_Version()
            << " Message Passing Error: Exchange_Messages_NoResChange. "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

   /*  for ( int iProc = 0; iProc !=  CFFC_MPI::Number_of_Processors; ++iProc ) { */
/*        if (  CFFC_MPI::This_Processor_Number == iProc ) { */
/*           cout<<"\n  after exchange  buffer function received location "<<endl; */
/*           cout<<"\n CFFC_MPI::This_Processor_Number = "<< CFFC_MPI::This_Processor_Number<<endl; */
      
      
/*           System::sleep(0.1); */
/*        } */
/*        MPI::COMM_WORLD.Barrier(); */
/*     } */
    
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
 
 /*    for ( int iProc = 0; iProc !=  CFFC_MPI::Number_of_Processors; ++iProc ) { */
/*        if (  CFFC_MPI::This_Processor_Number == iProc ) { */
/*           cout<<"\n  after unloading  buffer function received location "<<endl; */
/*           cout<<"\n CFFC_MPI::This_Processor_Number = "<< CFFC_MPI::This_Processor_Number<<endl; */
      
      
/*           System::sleep(0.1); */
/*        } */
/*        MPI::COMM_WORLD.Barrier(); */
/*     } */

   if (Send_Mesh_Geometry_Only) {
      for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
        if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK3D_USED) Soln_Blks[nb].Grid.Update_Ghost_Cells();
      } /* endfor */
    } /* endif */

 
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
