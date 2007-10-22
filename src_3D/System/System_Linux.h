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
#include <sys/stat.h> 
#include <stdexcept>

//=====End of Includes=====//

#define _MAX_PATH_  260

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
  
  // Change working directory
  int Change_Current_Working_Directory(const std::string &dirName);
  
  // Check if file exists
  bool Check_If_File_Exists(char * filename);
  bool Check_If_File_Exists(const std::string& filename);

  // Get current path
  void Get_Current_Path(char * buffer) throw(std::runtime_error);


}  // End of namespace System

#endif  // _SYSTEM_LINUX_INCLUDED
