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

//=====End of Includes=====//

namespace System
{

   // Create or replace a directory
   int Replace_Directory(const std::string &dirName);
   
   // Let the process sleep for a while
   int sleep(const double s);

}  // End of namespace System

#endif  // _SYSTEM_LINUX_INCLUDED
