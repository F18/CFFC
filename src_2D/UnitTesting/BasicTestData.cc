#include "BasicTestData.h"

namespace tut{

  void TestData::check_global_paths(void){
    DIR *dp;
    std::string FailureMsg;

    if ((dp = opendir(Global_TestSuit_Path)) == NULL){
      FailureMsg = "Critical ERROR: Cannot open directory: " + std::string(Global_TestSuit_Path);
      FailureMsg += "\n\t       Most likely you are NOT running the executable under the '/src_2D' directory!";
      throw tut::bad_ctor(FailureMsg);
    } else {
      // close the directory stream
      closedir(dp);
    }
  }
 
  void TestData::set_test_suit_path(char * dir_name){
    // sets the path to the test suit directory relative to "/src_2D"
    // this path is used to determine where the input and output files are located
    DIR *dp;
    std::string Path;
    std::string FailureMsg;

    Path = "./" + std::string(dir_name) + "/";

    // deallocate memory
    if (Global_TestSuit_Path != NULL){
      delete [] Global_TestSuit_Path; Global_TestSuit_Path = NULL;
    }

    // Set the Global_TestSuit_Path
    Global_TestSuit_Path = new char [Path.capacity() + 1];
    strcpy(Global_TestSuit_Path,Path.c_str());

    check_global_paths();

    if ((dp = opendir(Global_TestSuit_Path)) == NULL){
      try {
	FailureMsg = "Critical ERROR: Cannot open the test-suit directory: " + std::string(Global_TestSuit_Path);
	FailureMsg += "\n\t       Most likely you are NOT running the executable under the '/src_2D' directory!";
	throw tut::bad_ctor(FailureMsg);
      } 
      catch (...){
	delete [] Global_TestSuit_Path; Global_TestSuit_Path = NULL;
	throw;
      }
    } else {
      // close the directory stream
      closedir(dp);

      // Initialize input_file_name to the proper input directory
      InitializeInputFileName();

      // Initialize output_file_name to the proper input directory
      InitializeOutputFileName();
    }

  }


  void TestData::set_local_output_path(char * dir_name){ 
    DIR *dp;
    std::string Path;
    int stat;
    std::string FailureMsg;

    Path = std::string(Global_TestSuit_Path) + std::string(dir_name) + "/";

    // deallocate memory
    if (Local_Output_Path != NULL){
      delete [] Local_Output_Path; Local_Output_Path = NULL;
    }
    // Set the Local_Output_Path
    Local_Output_Path = new char [Path.capacity() + 1];
    strcpy(Local_Output_Path,Path.c_str());

    if ((dp = opendir(Local_Output_Path)) == NULL){
      try{ 

	stat = mkdir(Local_Output_Path, S_IRWXU|S_IRGRP|S_IROTH );
	
	if (stat != 0){
	  FailureMsg = "Critical ERROR: Cannot create directory: " + std::string(Local_Output_Path);
	  throw tut::bad_ctor(FailureMsg);
	}

	// Initialize output_file_name to the proper output directory
	InitializeOutputFileName();

      } catch(...){
	// deallocate memory
	delete [] Local_Output_Path; Local_Output_Path = NULL;
	throw;
      }
    } else {
      // close the directory stream
      closedir(dp);

      // Initialize output_file_name to the proper output directory
      InitializeOutputFileName();
    }

  }
  
  void TestData::set_local_input_path(char * dir_name){
    DIR *dp;
    std::string Path;
    std::string FailureMsg;

    Path = std::string(Global_TestSuit_Path) + std::string(dir_name) + "/";

    // deallocate memory
    if (Local_Input_Path != NULL){
      delete [] Local_Input_Path; Local_Input_Path = NULL;
    }

    Local_Input_Path = NULL;

    // Set the Local_Input_Path
    Local_Input_Path = new char [Path.capacity() + 1];
    strcpy(Local_Input_Path,Path.c_str());

    if ((dp = opendir(Local_Input_Path)) == NULL){
      try {
	FailureMsg = "Critical ERROR: Cannot open directory: " + std::string(Local_Input_Path);
	throw tut::bad_ctor(FailureMsg);
      } 
      catch (...){
	delete [] Local_Input_Path; Local_Input_Path = NULL;
	throw;
      }
    } else {
      // close the directory stream
      closedir(dp);

      // Initialize input_file_name to the proper input directory
      InitializeInputFileName();
    }

  }

  void TestData::RunRegressionTest(const std::string & FailureMsg, char* _CurrentFile_,
				   char* _MasterFile_, double AbsoluteTol, double RelativeTol,
				   int CompareStartsAtLine) throw(std::runtime_error,
								  std::invalid_argument,
								  tut::failure) 
  {

    /* Local Variables */
    char *File1, *File2;


    try {
      // Build File1 and File2 by adding the Local_Output_Path to _CurrentFile_ and _MasterFile_
      if (Local_Output_Path != NULL){
	File1 = new char [strlen(Local_Output_Path) + strlen(_CurrentFile_) + 1];
	strcpy(File1,Local_Output_Path);
	strcat(File1,_CurrentFile_);

	File2 = new char [strlen(Local_Output_Path) + strlen(_MasterFile_) + 1];
	strcpy(File2,Local_Output_Path);
	strcat(File2,_MasterFile_);

      } else {

	File1 = new char [strlen(Global_TestSuit_Path) + strlen(_CurrentFile_) + 1];
	strcpy(File1,Global_TestSuit_Path);
	strcat(File1,_CurrentFile_);

	File2 = new char [strlen(Global_TestSuit_Path) + strlen(_MasterFile_) + 1];
	strcpy(File2,Global_TestSuit_Path);
	strcat(File2,_MasterFile_);
      } // endif

      // Check if the output stream is still open and associated with either File1 or File2
      if ((strcmp(File1,output_file_name) == 0) || (strcmp(File2,output_file_name) == 0)){
	Close_Output_File();
      }

      // Compare the two files
      CompareFilesWithNumericalTolerance(FailureMsg,File1,File2,
					 AbsoluteTol,RelativeTol,CompareStartsAtLine);

      // remove File1(_CurrentFile_)
      int stat;
      string FailureMsg;

      stat = remove(File1);
      if (stat !=0){
	FailureMsg = "RunRegressionTest() WARNING: The file " + string(File1) + " couldn't be deleted\n";
	throw tut::warning(FailureMsg);
      }

      delete [] File1; File1 = NULL;
      delete [] File2; File2 = NULL;
    }
    catch(...) {
      // Try to remove File1(_CurrentFile_)
      remove(File1);

      // Deallocate memory
      delete [] File1; File1 = NULL;
      delete [] File2; File2 = NULL;

      throw;			// re-throw the error
    }
  }


  void TestData::CompareFilesWithNumericalTolerance(const std::string & FailureMsg, char * CurrentFile,
						    char * MasterFile, double AbsoluteTol, double RelativeTol,
						    int CompareStartsAtLine) throw(std::runtime_error,
										   std::invalid_argument,
										   tut::failure){

    /* Local Variables */
    pid_t pid, child_pid;	/* process ID */
    int stat_val(0);		/* status of the child process */
    int val(0);
    char AbsoluteError[50], RelativeError[50];	/* absolute and relative error as a string */
    char StartLine[10];
    string msg;
    char* FailMsg;

    /* Convert the absolute and relative tolerances */
    val = sprintf (AbsoluteError, "%12.12e", AbsoluteTol);
    if (val < 0){
      throw invalid_argument("RunRegressionTest() ERROR: invalid 'AbsoluteTol'.");
    }
    val = sprintf (RelativeError, "%12.12e", RelativeTol);
    if (val < 0){
      throw invalid_argument("RunRegressionTest() ERROR: invalid 'RelativeTol'.");
    }

    /* Build the StartLine string */
    val = sprintf(StartLine, "%d", CompareStartsAtLine);
    if (val < 0){
      throw invalid_argument("RunRegressionTest() ERROR: invalid 'CompareStartsAtLine'.");
    }
    strcat(StartLine,"-");

    /* Build the command line arguments for 'numdiff' */
    char * const numdiff_argv[] = {"numdiff", "-q", "-a", AbsoluteError, "-r", RelativeError ,
				   "-L", StartLine, CurrentFile , MasterFile , 0};
    char * const numdiff_argv_NoRelTol[] = {"numdiff", "-q", "-a", AbsoluteError, 
					    "-L", StartLine, CurrentFile , MasterFile , 0};

    /* Fork the program */
    pid = fork();

    switch(pid){
      /* Error from fork() */
    case -1:
      throw runtime_error("CompareFilesWithNumericalTolerance() ERROR: failed to fork the process.");

    case 0:
      /* Child process */

      /* Run 'numdiff' with 'numdiff_argv' */
      if (fabs(RelativeTol) > 0.0){
	/* Use RelativeTol */
	stat_val = execvp("numdiff", numdiff_argv);
      } else {
	/* Don't use RelativeTol */
	stat_val = execvp("numdiff", numdiff_argv_NoRelTol);
      }
      if (stat_val == -1){
	/* Error in launching numdiff */
	exit(-2);		/* the process returns either (-2) or (254) */
      }
      break;

    default:
      /* Parent process */

      /* Wait for the child process to finish */
      child_pid = wait(&stat_val);
      
      if (WIFEXITED(stat_val)){

	switch(WEXITSTATUS(stat_val)){

	case -1:
	  /* error handled by 255 */
	case 255:
	  /* numdiff terminated abnormally. Most likely wrong input parameters */
	  msg = "CompareFilesWithNumericalTolerance() ERROR: numdiff ERROR.\n";
	  msg += "\t       Most likely wrong input parameters.";
	  throw runtime_error(msg);
	  break;

	case -2:
	  /* error handled by 254 */
	case 254:
	  /* Error in launching numdiff */
	  msg = "CompareFilesWithNumericalTolerance() ERROR: failed to execute 'numdiff'.\n";
	  msg += "\t       Either numdiff is not installed or it is not in your path.";
	  throw runtime_error(msg);
	  break;

	default:
	  try{
	    /* Build the message for test failure */
	    msg = "CompareFilesWithNumericalTolerance() for " + FailureMsg + " failed.\n\t\t\t";
	    FailMsg = new char[msg.capacity()];
	    strcpy(FailMsg,msg.c_str());
	      
	    /* Check if the test succeeded --> i.e. numdiff exited with 0 */
	    ensure_equals(FailMsg, WEXITSTATUS(stat_val), 0);

	    delete [] FailMsg; FailMsg = NULL;
	  }
	  catch(...){
	    /* deallocate FailMsg */
	    delete [] FailMsg; FailMsg = NULL;
	    /* re-throw the error */
	    throw;
	  }

	} /* endswitch (WIFEXITED(stat_val)) */
      } else {
	throw runtime_error("CompareFilesWithNumericalTolerance() ERROR: numdiff terminated abnormally.");
      } /* endif */
    }	/* endswitch */

  }

}
