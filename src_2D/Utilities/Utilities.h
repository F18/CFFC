/*!\file Utilities.h
 \brief Header file defining useful macros and functions. */

// Original source from Thinking in C++, 2nd Edition
// Available at http://www.BruceEckel.com
// (c) Bruce Eckel 2000
// Copyright notice in Copyright.txt
// Test for error conditions in programs
// Local "using namespace std" for old compilers

#ifndef _UTILITIES_INCLUDED
#define _UTILITIES_INCLUDED

/* Include required C++ libraries. */

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <iostream>

/* Using std namespace functions */
using std::cout;
using std::setprecision;
using std::endl;
using std::runtime_error;

#define Print_(x) cout << setprecision(16) << #x " = " << x << endl;
#define Print_2(x,y) cout<< setprecision(16)<< #x " = " << x << ",  " #y " = " << y << endl;
#define Print_3(x,y,z) cout<< setprecision(16)<< #x " = " << x << ",  " #y " = " << y << ",  " #z " = " << z << endl;
#define Print_File(x,output_stream) output_stream << setprecision(16) << #x " = " << x << endl;
#define Print_File_2(x,y,output_stream) output_stream << setprecision(16)<< #x " = " << x << ",  " #y " = " << y << endl;
#define Print_File_3(x,y,z,output_stream) output_stream << setprecision(16)<< #x " = " << x << ",  " #y " = " << y << ",  " \
  #z " = " << z << endl;
#define Print_Digits(x,digits) cout << setprecision(digits) << #x " = " << x << endl;
#define TurnOff(x) /* #x */ cout << "\n Code: \n \"" #x << "\"\n was TURNED OFF!!!" << endl;

inline void require(bool requirement, 
                    const std::string& msg = "Requirement failed") throw(runtime_error) {
#ifdef _DO_CHECKS
  if (!requirement) {
    throw runtime_error(msg.c_str());
  }
#endif
}

inline void requireArgs(int argc, int args, 
			const std::string& msg = "Must use %d arguments"){
#ifdef _DO_CHECKS
  if (argc != args + 1) {
    fprintf(stderr, msg.c_str(), args);
    fputs("\n", stderr);
    exit(1);
  }
#endif
}

inline void requireMinArgs(int argc, int minArgs,
			   const std::string& msg = "Must use at least %d arguments"){
#ifdef _DO_CHECKS
  if(argc < minArgs + 1) {
    fprintf(stderr, msg.c_str(), minArgs);
    fputs("\n", stderr);
    exit(1);
  }
#endif
}
  
inline void assure(std::ifstream& in, 
		   const std::string& filename = ""){
#ifdef _DO_CHECKS
  if(!in) {
    fprintf(stderr, "Could not open file %s\n",
	    filename.c_str());
    exit(1);
  }
#endif
}

inline void assure(std::ofstream& out, 
		   const std::string& filename = ""){
#ifdef _DO_CHECKS
  if(!out) {
    fprintf(stderr, "Could not open file %s\n", 
	    filename.c_str());
    exit(1);
  }
#endif
}

#endif // _UTILITIES_INCLUDED
