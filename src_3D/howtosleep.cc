#define _POSIX_C_SOURCE 199309
#include <time.h>

namespace System
{

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

int sleep(const double s)
{

   const unsigned int sec = static_cast<unsigned int>(std::fabs(s));
   timespec req, rem;
   req.tv_sec = sec;
   req.tv_nsec = static_cast<long>((std::fabs(s) - sec)*1.E9);
   return nanosleep(&req, &rem);

}

}


/*==============================================================================
 *
 * Then in your code ...
 *
 *============================================================================*/

for ( int iProc = 0; iProc != numProc; ++iProc ) {
   if ( thisProc == iProc ) {
      std::cout << "=== Processsor " << iProc << " ===\n";
      System::sleep(0.1);
   }
   MPI::COMM_WORLD.Barrier();
}
