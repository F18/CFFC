/*!\file TestRun.h
  \brief Header file defining the elements needed for running the testing framework.*/

#ifndef _TEST_RUN_INCLUDED
#define _TEST_RUN_INCLUDED

#include <iostream>
#include "tut.h"
#include "tut_reporter.h"
#include <cstring>

using tut::reporter;
using tut::groupnames;

void Test_Run(std::string, int);

#endif // _TEST_RUN_INCLUDED
