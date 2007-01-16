//: :require.h
// From Thinking in C++, 2nd Edition
// Available at http://www.BruceEckel.com
// (c) Bruce Eckel 2000
// Copyright notice in Copyright.txt
// Test for error conditions in programs
// Local "using namespace std" for old compilers
#ifndef REQUIRE_H
#define REQUIRE_H
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

#define Print(x) std::cout << #x " = " << x << std::endl;
#define Print_(x) std::cout << setprecision(16) << #x " = " << x << std::endl;
#define Print_2(x,y) std::cout<<setprecision(16)<< #x " = " << x << ",  " #y " = " << y << std::endl;
#define Print_3(x,y,z) std::cout<<setprecision(16)<< #x " = " << x << ",  " #y " = " << y << ",  " #z " = " << z << std::endl;
#define Print_Digits(x,digits) std::cout << setprecision(digits) << #x " = " << x << std::endl;
#define TurnOff(x) /* #x */ std::cout << "\nElement " #x << " was turned off" << std::endl;

inline void require(bool requirement, 
  const std::string& msg = "Requirement failed"){
  using namespace std;
  if (!requirement) {
    fputs(msg.c_str(), stderr);
    fputs("\n", stderr);
    exit(1);
  }
}

inline void requireArgs(int argc, int args, 
  const std::string& msg = 
    "Must use %d arguments") {
  using namespace std;
   if (argc != args + 1) {
     fprintf(stderr, msg.c_str(), args);
     fputs("\n", stderr);
     exit(1);
   }
}

inline void requireMinArgs(int argc, int minArgs,
  const std::string& msg =
    "Must use at least %d arguments") {
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

#endif // REQUIRE_H ///:~
