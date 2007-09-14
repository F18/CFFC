#include <tut.h>
#include "include/HeaderData.h"

using namespace std;

namespace tut
{
  struct Data_HeaderData{
    HeaderData Names;
  };

  typedef test_group<Data_HeaderData> HeaderData_TestGroup;
  typedef HeaderData_TestGroup::object HeaderData_TestObject;
}

namespace tut
{

  /*Test1: Add elements */
  template<>
  template<>
  void HeaderData_TestObject::test<1>()
  {
    // Set Path
    ofstream output_file;
    char output_file_name[256];
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"HeaderData.dat");
    output_file.open(output_file_name, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    output_file << "VARIABLES = ";

    Names.add("Var1");
    Names.add("Var2");
    Names.add("p");
    Names.PrintHeaderTecplot(output_file);
    cout << Names[1] << endl;

    output_file.close();



  }

}

namespace tut
{
  HeaderData_TestGroup HeaderData_Test("HeaderData_Test");
}
