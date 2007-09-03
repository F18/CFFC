// Define the basic test data class.

/* Include required C++ libraries. */
#include <cmath>
#include <sstream>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <dirent.h>
#include <tut.h>
#include <limits>

/* Using std namespace functions */
using std::invalid_argument;
using std::runtime_error;
using std::ostringstream;

/* Include CFFC header files */

#include "../Utilities/Utilities.h"
#include "../CFD/CFD.h"

namespace tut {

  // TestData class definition
  class TestData{
  public:
    int digits;			// Number of exact imposed digits --> gives the precision
    double tol;	                // Numerical tolerance for comparing the solutions
    double EpsMachine;          // Machine accuracy
    char output_file_name[256];	// Output file name
    char input_file_name[256];	// Input file name
    int  RunRegression;		// Flag for running regression

    char *MasterFile, *CurrentFile; /* File name variables used for regression tests */
    char* Msg;			/* message variable */
    std::string Message;		/* message variable */

    /* Relative path to the "/src_2D" directory for outputing and inputing data */
    char* Global_TestSuit_Path;

    /* Local output and input path --> specific for each test */
    char *Local_Output_Path, *Local_Input_Path;

    // Constructor
    TestData(void);

    // Destructor
    ~TestData(void);

    // Fiels Access
    ofstream & out(void){ return output_file; }
    ifstream & in(void){return in_file; }
    ostringstream & ostm(void){ return output_string_stream;}

    // Open the output file
    void Open_Output_File(const char * file_name);
    void Open_Output_File(char * file_name){const char * ptr = file_name; return Open_Output_File(ptr);}
    void Close_Output_File(void);
    // Initialize output file name
    void InitializeOutputFileName(void);
    // Access to OutputFileReady
    int OFR(void) {return output_file.is_open();}
    // Open the file used for output for reading
    void Open_Output_File_For_Reading(void);
    // Open output stream
    void open_output_stream(ofstream & output_stream, const char * output_stream_file_name);

    // Open the input file
    void Open_Input_File(const char * file_name);
    void Open_Input_File(char * file_name){const char * ptr = file_name; return Open_Input_File(ptr);}
    void Close_Input_File(void);
    // Initialize input file name
    void InitializeInputFileName(void);
    // Access to InputFileReady
    int IFR(void) {return in_file.is_open();}
    // Open input stream
    void open_input_stream(ifstream & input_stream, const char * input_stream_file_name);

    void set_test_suit_path(char * dir_name);
    void set_local_output_path(char * dir_name);
    void set_local_input_path(char * dir_name);

    /* ostmClear() --> clears output_string_stream */
    void ostmClear(void){ ostm().str(""); ostm().clear(); }

    /* Run Regression Test --> compares the CurrentFile against the MasterFile, ignoring numeric differences within the 
                               given absolute or relative tolerances.
                               The files are compared starting at line 'CompareStartsAtLine' 
                               (this argument can be used to jump over text that is changes at every runtime) */
    void RunRegressionTest(const std::string & FailureMsg, char * CurrentFile,
			   char * MasterFile, double AbsoluteTol, double RelativeTol,
			   int CompareStartsAtLine = 1) throw(std::runtime_error,
							     std::invalid_argument,
							     tut::failure);
    void RunRegressionTest(const std::string & FailureMsg, char * CurrentFile,
			   char * MasterFile, double AbsoluteTol, int CompareStartsAtLine = 1) throw(std::runtime_error,
												    std::invalid_argument,
												    tut::failure){
      RunRegressionTest(FailureMsg,CurrentFile,MasterFile,AbsoluteTol, 0.0, CompareStartsAtLine);
    }

    void CompareFilesWithNumericalTolerance(const std::string & FailureMsg, char * CurrentFile,
					    char * MasterFile, double AbsoluteTol, double RelativeTol,
					    int CompareStartsAtLine = 1) throw(std::runtime_error,
									       std::invalid_argument,
									       tut::failure);

    /* Check_Input_Output_Operator --> outputs the initial Obj to a file and then reads the file in order
                                       to initialize an object of the same type as Obj. 
				       Outputs the second object to a different file and then compares
				       them with a tight tolerance using numdiff. 
                     Prerequisites:  Class T must have defined the following operators:
                                       --> operator << (output)
                                       --> operator >> (input)
      NB: You should check the output operator independently, to make sure that what you output
          is what you want!!!
    */
    template<class T>
      void Check_Input_Output_Operator(const char * msg, const T &Obj);
    template<class T>
      void Check_Input_Output_Operator(const T &Obj){ return Check_Input_Output_Operator("operator << & >>", Obj); }
    template<class T>
      void Check_Input_Output_Operator(const std::string & msg, const T &Obj);


    /* Local class */
    class __LocalTest__{
    public:
      double x,y,z;

      __LocalTest__(void):x(0.0),y(0.0),z(0.0){};
      friend bool operator== (const __LocalTest__ &left, const __LocalTest__ &right);
      friend bool operator!= (const __LocalTest__ &left, const __LocalTest__ &right){ return !(left==right); }
      friend ostream & operator<< (ostream & os, const __LocalTest__ &Obj);
      friend istream & operator>> (istream & os,  __LocalTest__ &Obj);
    };

  protected:
    ofstream output_file;	// Output file stream
    ifstream in_file;		// Input file stream
    ostringstream output_string_stream; // Stream for formating to strings 
                                        // (To retrieve the resulting string (of type std::string),
                                        //  you call the member function str() )

  private:
    /* Check the existance of the Global paths */
    void check_global_paths(void);

  };

  /* Constructor */
  inline TestData::TestData(void){

    // Initialize tolerance
    #ifdef _GNU_GCC_V3
    digits = numeric_limits<double>::digits10;
    EpsMachine = numeric_limits<double>::epsilon();
    #else
    digits = 14;
    EpsMachine = 1.0e-14;
    #endif
    tol = 0.5*pow(10.0,-digits);

    /* Set pointers to NULL */
    MasterFile = NULL; CurrentFile = NULL;
    Local_Output_Path = NULL; Local_Input_Path = NULL;
    Global_TestSuit_Path = NULL;
    Msg = NULL;

    /* Set the Global_TestSuit_Path to the current directory */
    Global_TestSuit_Path = new char [3];
    strcpy(Global_TestSuit_Path, "./");

    // Initialize output_file_name to the proper output directory
    InitializeOutputFileName();
    // Initialize input_file_name to the proper input directory
    InitializeInputFileName();

    RunRegression = ON;		/* run regression tests */
  }

  /* Destructor */
  inline TestData::~TestData(void){ 
    // deallocate memory
    delete [] Local_Output_Path; Local_Output_Path = NULL;
    delete [] Local_Input_Path; Local_Input_Path = NULL;
    delete [] Msg; Msg = NULL;

    if (OFR()){
      // close the output stream
      output_file.close();
    }
    if (IFR()){
      // close the input stream
      in_file.close();
    }
  }
  
  /* open_output_stream() */
  inline void TestData::open_output_stream(ofstream & output_stream, const char * output_stream_file_name){
    /* Associate the output_stream with the output_stream_file_name */
    output_stream.open(output_stream_file_name, ios::out);

    /* Check for successful operation */
    if (!output_stream.is_open()){
      string msg;
      msg = "Unable to open output file " + string(output_stream_file_name);
      throw runtime_error(msg);
    }
  }

  /* open_input_stream() */
  inline void TestData::open_input_stream(ifstream & input_stream, const char * input_stream_file_name){
    /* Associate the input_stream with the input_stream_file_name */
    input_stream.open(input_stream_file_name, ios::in);

    /* Check for successful operation */
    if (!input_stream.is_open()){
      string msg;
      msg = "Unable to open input file " + string(input_stream_file_name);
      throw runtime_error(msg);
    }
  }

  /* Initialize output file name */
  inline void TestData::InitializeOutputFileName(void){
    if (Local_Output_Path != NULL){
      strcpy(output_file_name, Local_Output_Path);
    } else {
      strcpy(output_file_name, Global_TestSuit_Path);
    }
  }

  /* Open_Output_File */
  inline void TestData::Open_Output_File(const char * file_name){
    /* close the output stream if open and reset the path */
    Close_Output_File();

    /* set the output_file_name */
    strcat(output_file_name, file_name);

    /* open the output stream */
    open_output_stream(output_file, output_file_name);
  }

  /* Open_Output_File_For_Reading */
  inline void TestData::Open_Output_File_For_Reading(void){
    /* close the in_file */
    if(IFR()){
      in_file.close();
    }

    /* open the in_file stream with 'output_file_name' */
    open_input_stream(in_file, output_file_name);
  }

  /* Close_Output_File */
  inline void TestData::Close_Output_File(void){
    if (OFR()){
      output_file.close();
      InitializeOutputFileName();
    }
  }

  /* Initialize input file name */
  inline void TestData::InitializeInputFileName(void){
    if (Local_Input_Path != NULL){
      strcpy(input_file_name, Local_Input_Path);
    } else {
      strcpy(input_file_name, Global_TestSuit_Path);
    }
  }
  
  /* Open_Input_File */
  inline void TestData::Open_Input_File(const char * file_name){
    /* close the input stream if open and reset the path */
    Close_Input_File();

    /* set the input_file_name */
    strcat(input_file_name, file_name);

    /* open the input stream */
    open_input_stream(in_file, input_file_name);
  }

  /* Close_Input_File */
  inline void TestData::Close_Input_File(void){
    if(IFR()){
      in_file.close();
      InitializeInputFileName();
    }
  }

  /* __LocalTest__ operator<< */
  inline ostream & operator<< (ostream & os, const TestData::__LocalTest__ &Obj){
    os.setf(ios::skipws);
    os << setprecision(16) << Obj.x << "\t" << Obj.y << "\t" << Obj.z << endl;
    return os;
  }

  /* __LocalTest__ operator>> */
  inline istream & operator>> (istream & os, TestData::__LocalTest__ &Obj){
    os.setf(ios::skipws);
    os >> Obj.x >> Obj.y >> Obj.z ;
    os.unsetf(ios::skipws);
    return os;
  }

  /* __LocalTest__ operator== */
  inline bool operator== (const TestData::__LocalTest__ &left, const TestData::__LocalTest__ &right){
    return (left.x==right.x && left.y==right.y && left.z==right.z );
  }

  // Check_Input_Output_Operator()
  template<class T>
    void TestData::Check_Input_Output_Operator(const char * msg, const T &Obj){
    
    ofstream out_CIOO;	// Output file stream
    ifstream in_CIOO;      // Input file stream
    char file_name_CIOO1[100], file_name_CIOO2[100]; 

    /* create a temporary object of class T */
    T Temp;

    strcpy(file_name_CIOO1, "/tmp/RegressionTest_Check_Input_Output_Operator_I");
    strcpy(file_name_CIOO2, "/tmp/RegressionTest_Check_Input_Output_Operator_II");
    open_output_stream(out_CIOO, file_name_CIOO1);
    open_input_stream(in_CIOO, file_name_CIOO1);

    /* output the input object */
    out_CIOO << Obj << endl;

    /* read the temporary object from the file */
    in_CIOO >> Temp;

    /* close the streams */
    out_CIOO.close();
    in_CIOO.close();

    /* output the temporary object */
    open_output_stream(out_CIOO, file_name_CIOO2);
    out_CIOO << Temp << endl;

    /* check equality by comparing the two files */
    CompareFilesWithNumericalTolerance(msg,file_name_CIOO2,file_name_CIOO1,1.0e-14,1.0e-14);

    /* remove the temporary files */
    int stat;
    string FailureMsg;

    stat = remove(file_name_CIOO1);
    if (stat !=0){
      FailureMsg = "Check_Input_Output_Operator() ERROR: The file " + string(file_name_CIOO1) + " couldn't be deleted\n";
      throw tut::warning(FailureMsg);
    }

    stat = remove(file_name_CIOO2);
    if (stat !=0){
      FailureMsg = "Check_Input_Output_Operator() ERROR: The file " + string(file_name_CIOO2) + " couldn't be deleted\n";
      throw tut::warning(FailureMsg);
    }

  }

  template<class T> inline
    void TestData::Check_Input_Output_Operator(const std::string & msg, const T &Obj){
    char * Message = new char [msg.length() + 1];
    strcpy(Message,msg.c_str());

    try{

      Check_Input_Output_Operator(Message, Obj);
      delete [] Message; Message = NULL;

    } catch(...){
      delete [] Message; Message = NULL;
      throw;
    }
  }



}
