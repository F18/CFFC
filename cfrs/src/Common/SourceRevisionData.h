/* SourceRevisionData.h:  Header file defining the source revision parameters
                          which are characteristic for the executable that is built. */

#ifndef _SOURCE_REVISION_DATA_INCLUDED
#define _SOURCE_REVISION_DATA_INCLUDED

/* Include header files. */
#include <fstream>

class SourceCode{

 public:
  static char* CodeName;
  static char* CompilationTime;
  static char* LastChanged_Author;
  static char* LastCommitted_Hash;
  static char* LastChanged_Date;

  static const char* ProgramName(void){return CodeName;}
  static const char* TimeAtCompilation(void){return CompilationTime;}
  static const char* LastCommitted_Author(void){return LastChanged_Author;}
  static const char* LastCommitted_HashID(void){return LastCommitted_Hash;}
  static const char* LastCommitted_Date(void){return LastChanged_Date;}

  static void PrintTecplotAuxData(std::ofstream &output_file);

 protected:
  SourceCode(){};
};

/***************************************************
 * PrintTecplotAuxData()                           *
 * Print the source code data as auxiliary dataset *
 *                                                 *
 **************************************************/
inline void SourceCode::PrintTecplotAuxData(std::ofstream &output_file){

  output_file << "DATASETAUXDATA Code_Used = \"  " << SourceCode::ProgramName() << "     \" \n "
	      << "DATASETAUXDATA Compiled_On = \"  " << SourceCode::TimeAtCompilation() << " \" \n "
	      << "DATASETAUXDATA Last_Committed_Date = \"  " << SourceCode::LastCommitted_Date() << "    \" \n "
	      << "DATASETAUXDATA Last_Committed_Revision = \"  " << SourceCode::LastCommitted_HashID() << "     \" \n ";
}


#endif
