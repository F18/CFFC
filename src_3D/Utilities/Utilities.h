//: :Utilities.h

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
#include <string> 
#include <sstream>


/* Using std namespace functions */
using std::cout;
using std::setprecision;
using std::setw;
using std::endl;

#define Print_(x) cout << setprecision(16) << #x " = " << x << endl;
#define Print_2(x,y) cout<< setprecision(16)<< #x " = " << x << ",  " #y " = " << y << endl;
#define Print_3(x,y,z) cout<< setprecision(16)<< #x " = " << x << ",  " #y " = " << y << ",  " #z " = " << z << endl;
#define Print_File(x,output_stream) output_stream << setprecision(16) << #x " = " << x << endl;
#define Print_File_2(x,y,output_stream) output_stream << setprecision(16)<< #x " = " << x << ",  " #y " = " << y << endl;
#define Print_File_3(x,y,z,output_stream) output_stream << setprecision(16)<< #x " = " << x << ",  " #y " = " << y << ",  " \
  #z " = " << z << endl;
#define Print_Digits(x,digits) cout << setprecision(digits) << #x " = " << x << endl;
#define TurnOff(x) /* #x */ cout << "\nElement " #x << " was turned off" << endl;

inline void require(bool requirement, 
                    const std::string& msg = "Requirement failed"){
#ifndef __No_Checking__
  using namespace std;
  if (!requirement) {
    fputs(msg.c_str(), stderr);
    fputs("\n", stderr);
    exit(1);
  }
#endif
}

inline void requireArgs(int argc, int args, 
			const std::string& msg = "Must use %d arguments") {
  using namespace std;
  if (argc != args + 1) {
    fprintf(stderr, msg.c_str(), args);
    fputs("\n", stderr);
    exit(1);
  }
}

inline void requireMinArgs(int argc, int minArgs,
			   const std::string& msg = "Must use at least %d arguments") {
  using namespace std;
  if(argc < minArgs + 1) {
    fprintf(stderr, msg.c_str(), minArgs);
    fputs("\n", stderr);
    exit(1);
  }
}
  
inline void assure(std::ifstream& in, 
		   const std::string& filename = "") {
  using namespace std;
  if(!in) {
    fprintf(stderr, "Could not open file %s\n",
	    filename.c_str());
    exit(1);
  }
}

inline void assure(std::ofstream& out, 
		   const std::string& filename = "") {
  using namespace std;
  if(!out) {
    fprintf(stderr, "Could not open file %s\n", 
	    filename.c_str());
    exit(1);
  }
}

/// Converts to std::string
/// Don't use this to convert to a char, use c_str for that.
/// Typical use is to convert to numbers.
/// @param str string to convert from
/// @return converter type
template <class T>
static std::string to_str (T v)
{
    std::ostringstream oss;
    oss << v;
    return oss.str();
}

/// Converts from std::string
/// Don't use this to convert to a char, use c_str for that.
/// Typical use is to convert to numbers.
/// @param str string to convert from
/// @return converter type
template <class T>
static T from_str (const std::string& str)
{
    T v;
    if (str.length() > 0)
    {
        std::istringstream iss(str.c_str());
        iss >> v;
    }
    else
    {
        // pretty much everything has an empty constuctor
        v = T();
    }
    return v;
  }

enum ProgressMode {PROGRESS_MODE_SILENT,PROGRESS_MODE_MESSAGE,PROGRESS_MODE_FILE,PROGRESS_MODE_TERMINAL};

inline void ShowProgress(std::string message, int numIn, int maximum, int mode) {

    switch (mode) {
        case PROGRESS_MODE_SILENT :
            break;
        case PROGRESS_MODE_TERMINAL : 
        {
            char barspin[16] = {'\\','\\','\\','\\',
                '|', '|','|','|',
                '/','/', '/', '/',
            '-','-','-','-'};
            
            int whichOne = numIn % 16;
            
            int last_index = maximum;
            int percent = int(100*numIn/double(last_index));
            
            std::cout << '\r'
            << message << setw(3)  <<  percent << " %"
            << "  " << barspin[whichOne] << " ";
            std::cout.flush();
            if (percent == 100) {
                std::cout << '\r'
                << message << setw(3)  <<  percent << " %      " << std::endl;
            }  
            break;
        }
        case PROGRESS_MODE_FILE : 
        {
            int first_index = 1;
            int last_index = maximum;
            int percent = int(100*numIn/double(last_index));
            if (numIn == first_index) {
                std::cout << message << endl;
            }
            int previous_percent = int(100*(numIn-1)/double(last_index));
            if (percent != previous_percent || numIn == first_index) {
                std::cout << " " << percent << "%";
                std::cout.flush();
            }
            if (percent == 100) {
                std::cout << std::endl;
            }
            break;
        }
        case PROGRESS_MODE_MESSAGE :
        {
            // no progress, just message
            int first_index = 1;
            if (numIn == first_index) {
                std::cout << message << endl;
            }
        }
        default :
            break;
    }
    return;
}



#endif // _UTILITIES_INCLUDED
