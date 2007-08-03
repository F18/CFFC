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

