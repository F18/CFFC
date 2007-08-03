#ifndef _SYSTEM_LINUX_INCLUDED
#define _SYSTEM_LINUX_INCLUDED

/*******************************************************************************
 *
 * Routines specific to the Linux operating system.
 *
 ******************************************************************************/

//=====Included Files=====//

//-----Standard Library-----//

#include <string>
#include <iostream>
#include <fstream>

//=====End of Includes=====//

namespace System
{

   // Create or replace a directory
   int Replace_Directory(const std::string &dirName);

   // Compress old restart (do this before
   // over-writing).
   void Compress_Restart();

   // Uncompress old restart
   void Uncompress_Restart();

   // Set flag to indicate restart files are
   // being written
   void Set_Restart_Flag();

   // Remove restart flag
   void Remove_Restart_Flag();

   // check for Restart flag
   int Restart_In_Progress();
   
}  // End of namespace System

#endif  // _SYSTEM_LINUX_INCLUDED
