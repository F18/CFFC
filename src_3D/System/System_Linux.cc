/*******************************************************************************
 *
 * This file contains routines specific to the Linux operating system.
 *
 ******************************************************************************/

//=====Included Files=====//

#define _POSIX_C_SOURCE 199309

//-----Standard Library-----//

#include <cstdlib>
#include <string>

//-----System-----//

#include <time.h>

//-----Internal-----//

#include "System_Linux.h"

//=====End of Includes=====//


/*==============================================================================
 *
 * Routine: Replace_Directory
 *
 * Purpose
 * =======
 *
 *   Creates a directory overwriting any that may already exists (including
 *   any contents).
 *
 *============================================================================*/

int System::Replace_Directory(const std::string &dirName)
   
{

   std::system(("rm -f -r " + dirName).c_str());
   std::system(("mkdir -p " + dirName).c_str());
   return 0;
   
}


/*==============================================================================
 *
 * Routine: sleep
 *
 * Purpose
 * =======
 *
 *   Let the process sleep for a while
 *
 * I/O
 * ===
 *
 *   s                  - (I) Time to sleep in seconds
 *   returns            -  0: successful sleep
 *                        -1: interupted
 *
 *============================================================================*/

int System::sleep(const double s)
{

   const unsigned int sec = static_cast<unsigned int>(std::abs(s));
   timespec req, rem;
   req.tv_sec = sec;
   req.tv_nsec = static_cast<long>((std::abs(s) - sec)*1.E9);
   return nanosleep(&req, &rem);

}
