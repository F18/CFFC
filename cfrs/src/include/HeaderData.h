// HeaderData.h defines a class for storing the names of the variables in structure
#include <cstring>
#include <fstream>
#include <vector>

#ifndef _HEADERDATA_INCLUDED
#define _HEADERDATA_INCLUDED

using namespace std;

class HeaderData : public std::vector<std::string> {
 public:

  HeaderData(){};
  HeaderData(const std::string & Var){
    push_back(Var);
  }

  void add(const std::string & Name){
    push_back(Name);
  }
  void PrintHeaderTecplot(std::ofstream &output_file){
    int i;
    for(i=0; i< size(); ++i)
      output_file << " \" " << at(i) << "\" \\ \n";
  }

  void PrintDataTypeTecplot(std::ofstream &output_file, int SpaceDimension){
    int i;
    output_file << "DT = (";
    for(i=1; i<=SpaceDimension + size(); ++i){
      output_file << " DOUBLE ";
    }
    output_file << ") \n";
  }

  friend std::ostream& operator << (std::ostream &os, const HeaderData & Obj){
    for(int i=0; i< Obj.size(); ++i)
      os << Obj.at(i) << "\t";
    return os;
  }

};

class State{

 public:
  static HeaderData PrintedVariables;


};


#endif
