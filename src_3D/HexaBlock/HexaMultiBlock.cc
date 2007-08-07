/* Euler3DHexaMultiBlock.cc:  Multi-Block Versions of Subroutines for 3D Euler 
                              Multi-Block Hexarilateral Mesh 
                              Solution Classes. */

/* Include 3D Euler hexahedrial mesh solution header file. */

#ifndef _HEXA_MULTIBLOCK_INCLUDED
#include "HexaMultiBlock.h"
#endif // _HEXA_MULTIBLOCK_INCLUDED

// int  Local_Block_List::Send_All_Boundary_Info
// (MultiBlk_Connectivity &Multiblock_conn,
//  Input_Parameters &Input_Parameters){

//    int error_flag;
//    int size = Local_SolnBlk_List.size();
//    int imssg = 0;
   
//    error_flag = 0;
//    MPI::Request *SendRecv_requests;
   
//    //allocate using new    
//    SendRecv_requests = new MPI::Request [size*52];
   
//    /* Update the solution for each solution block. */
//    for(SolnBlk_Itr SBit = Local_SolnBlk_List.begin(); 
//        SBit!= Local_SolnBlk_List.end(); ++SBit){
      
//       error_flag = Multiblock_conn.Send_MultiBlk_Boundary_Info(
//          SendRecv_requests,
//          imssg,
//          Input_Parameters, (*SBit));

//       if (error_flag) return (error_flag);
      
//    }  /* endfor */
   
//    // call wait for all the messages on this process being sent and received 
//    // to be completed
 
//    if(imssg > 0){
        
//       MPI::Request::Waitall(imssg, SendRecv_requests);
//      //if( CFFC_MPI::This_Processor_Number == 1 ){
// 	//for (int im = 0; im<imssg; ++im){
//    //		cout<<"\n imssg  = "<<im<<" Test=   "<<SendRecv_requests[im].Test()<<endl;
// 	//}       
//       //}
//       // after passing the conservative solution vector
//       // should convert to primitive solution vector
//       WtoU();
      
//    }
   
//    //deallocate the memory allocated by using new
//    delete []SendRecv_requests; SendRecv_requests = NULL;
   
//    return(error_flag);
      
// }

// void Local_Block_List::Message_Passing_Datatype(
//    MultiBlk_Connectivity &MultiBlock_Connectivity,
//    Input_Parameters &Input_Parameters){

//    // only need use one solution block in this block list
//    // to setup the message passing interface datatype   
//    SolnBlk_Itr SBit = Local_SolnBlk_List.begin();
//    MultiBlock_Connectivity.Assign_SendRecv_Datatype_SolutionBlk
//       (Input_Parameters, *SBit);
   
// }


   
