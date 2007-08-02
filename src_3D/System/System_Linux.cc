/*******************************************************************************
 *
 * This file contains routines specific to the Linux operating system.
 *
 ******************************************************************************/

//=====Included Files=====//

//-----Standard Library-----//

#include <cstdlib>
#include <string>

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
