/* SourceRevisionData.h:  Header file defining the source revision parameters
                          which are characteristic for the executable that is built. */

#ifndef _SOURCE_REVISION_DATA_INCLUDED
#define _SOURCE_REVISION_DATA_INCLUDED

/* Include header files. */

class SourceCode{

 public:
  static char* CodeName;
  static char* Revision;
  static char* CompilationTime;
  static char* LastChanged_Author;
  static char* LastChanged_Revision;
  static char* LastChanged_Date;

  static const char* ProgramName(void){return CodeName;}
  static const char* RevisionAtCompileTime(void){return Revision;}
  static const char* TimeAtCompilation(void){return CompilationTime;}
  static const char* LastCommitted_Author(void){return LastChanged_Author;}
  static const char* LastCommitted_Revision(void){return LastChanged_Revision;}
  static const char* LastCommitted_Date(void){return LastChanged_Date;}

 private:
  SourceCode(){};
};


#endif
