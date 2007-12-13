/*******************************************************************************
 *
 * This file contains routines specific to the Linux operating system.
 *
 ******************************************************************************/

//=====Included Files=====//

//-----Standard Library-----//

//#include <cstdlib>    //included in .h file
//#include <string>

//-----Internal-----//

#include "System_Linux.h"

//=====End of Includes=====//


/*------------------------------------------------------------------------------
 *
 * Routine: Replace_Directory
 *
 * Purpose
 * =======
 *
 *   Creates a directory overwriting any that may already exists (including
 *   any contents).
 *
 *----------------------------------------------------------------------------*/

int System::Replace_Directory(const std::string &dirName)
   
{

   std::system(("rm -f -r " + dirName).c_str());
   std::system(("mkdir -p " + dirName).c_str());
   return 0;
   
}

/*------------------------------------------------------------------------------
 *
 * Routine: Compress_Restart
 *
 * Purpose
 * =======
 *
 *   This will compress restart data
 *   (*.tree and *.soln) into backup.tar.gz
 *  
 *   This is done using the "find" command just in case there are
 *   more *.soln files that can be handles simultaneously.
 *
 *----------------------------------------------------------------------------*/

void System::Compress_Restart()

{

  std::system("tar -cf backup.tar *.tree");
  std::system("find ./ -name '*.soln' -print0 | xargs -0 tar -rf backup.tar");
  std::system("gzip -f backup.tar");
  std::system("rm -f *.tree");
  std::system("find ./ -name '*.soln' -print0 | xargs -0 rm -f");
  return;

}

/*------------------------------------------------------------------------------
 *
 * Routine: Uncompress_Restart
 *
 * Purpose
 * =======
 *
 *   This will uncompress restart data
 *   (*.tree and *.soln) from backup.tar.gz
 *  
 *----------------------------------------------------------------------------*/

void System::Uncompress_Restart()

{

  std::system("rm -f *.tree");
  std::system("find ./ -name '*.soln' -print0 | xargs -0 rm -f");
  std::system("tar -xzf backup.tar.gz");
  return;

}

/*------------------------------------------------------------------------------
 *
 * Routine: Set_Restart_Flag
 *
 * Purpose
 * =======
 *
 *  This will create a file ".restart_being_written" to indicate
 *  that a restart is being saved
 *  
 *----------------------------------------------------------------------------*/

void System::Set_Restart_Flag()

{
  std::system("touch .restart_being_written");
  return;
}

/*------------------------------------------------------------------------------
 *
 * Routine: Remove_Restart_Flag
 *
 * Purpose
 * =======
 *
 *  This will delete the restart flag
 *  
 *----------------------------------------------------------------------------*/

void System::Remove_Restart_Flag()

{
  std::system("rm -f .restart_being_written");
  return;
}

/*------------------------------------------------------------------------------
 *
 * Routine: Restart_in_Progress
 *
 * Purpose
 * =======
 *
 *  This will check if ".restart_being_written" exists.  If so
 *  a restart is being written or the process writting the
 *  restarts was killed during the process, in which case
 *  the restart files should be considered corrupt.
 *  
 *----------------------------------------------------------------------------*/

int System::Restart_In_Progress()

{
  std::ifstream test;

  test.open(".restart_being_written");

  if(test.fail()) {
    test.close();
    return 0;
  }

  test.close();
  return 1;

}

/*------------------------------------------------------------------------------
 *
 * Routine: Change_Current_Working_Directory
 *
 * Purpose
 * =======
 *
 *  This will change the current working directory to 
 *  the path 'dirName'.
 *  If the operation is successful the returned value is 0,
 *  otherwise is -1.  
 *
 *----------------------------------------------------------------------------*/

int System::Change_Current_Working_Directory(const std::string &dirName)

{
  return chdir(dirName.c_str());
}

/*------------------------------------------------------------------------------
 *
 * Routine: Check_If_File_Exists
 *
 * Purpose
 * =======
 *
 *  This will check whether the file with the 'filename' exists on the disk.
 *  This subroutine doesn't check to see what permissions the current user
 *  has on the file.
 *
 *----------------------------------------------------------------------------*/

bool System::Check_If_File_Exists(char * filename)

{
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(filename,&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  
  return(blnReturn);
}

/*------------------------------------------------------------------------------
 *
 * Routine: Check_If_File_Exists
 *
 * Purpose
 * =======
 *
 *  Subroutine with 'filename' of type string instead of char*
 *----------------------------------------------------------------------------*/

bool System::Check_If_File_Exists(const std::string& filename)

{
  return Check_If_File_Exists(filename.c_str());
}

/*------------------------------------------------------------------------------
 *
 * Routine: Get_Current_Path
 *
 * Purpose
 * =======
 *
 *  This will write the current path in 'buffer' variable.
 *  If the operation is successful 'Result' is different than NULL.
 *  If the operation is unsuccessful 'Result' is NULL,
 *  and the subroutine throws an error.
 *  Un unsuccessful operation is most likely due to the buffer size (too small!).
 *
 *----------------------------------------------------------------------------*/

void System::Get_Current_Path(char * buffer) throw(std::runtime_error)

{
  char * Result(NULL);
  Result = getcwd(buffer, _MAX_PATH_);

  if (Result == NULL){
    throw std::runtime_error("Get_Current_Path() ERROR: failed to identify the current directory.");
  }

}
