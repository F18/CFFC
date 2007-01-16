#include <tut.h>
#include "Grid/SubGrid/SubGridMesh.h"
#include "include/require.h"
#include "Common/Containers/PointWiseSolution.h"

using namespace std;
using namespace std;

namespace tut
{
  struct Data_SubGridMesh{
    SubGridMesh<Node2D,OneD,double> SubGrid1D;
    SubGridMesh<> SubGrid2D;
    SubGridMesh<Node3D,ThreeD,double> SubGrid3D;
    double val;			// the same type as the SolutionType declared for the SubGridMeshes
    vector<int> XYZ;
    int Nx,Ny,Nz;
  };

  typedef test_group<Data_SubGridMesh> SubGridMesh_TestGroup;
  typedef SubGridMesh_TestGroup::object SubGridMesh_TestObject;
}

namespace tut
{

//   /*Test1: Constructor 1D*/
  template<>
  template<>
  void SubGridMesh_TestObject::test<1>()
  {
    Nx = 5;			// number of nodes
    XYZ.push_back(Nx);
    SubGrid1D.SetMesh(XYZ);
    val=0;
    Node2D node(1,1);
    for (int i=0; i<Nx; i++){
      ensure("Constructor1D solution", val == SubGrid1D(i).GetValue());
      ensure("Constructor1D node", node == SubGrid1D(i).GetNode());
    }
  }

  /*Test2: Constructor 2D*/
  template<>
  template<>
  void SubGridMesh_TestObject::test<2>()
  {
    Nx = 5; Ny = 10;
    // 2D vector
    XYZ.push_back(Nx); XYZ.push_back(Ny);
    SubGrid2D.SetMesh(XYZ);
    val=0;
    Node2D node(1,1);
    for (int i=0; i<Nx; i++)
      for (int j=0; j<Ny; j++){
	ensure("Constructor2D solution", val == SubGrid2D(i,j).GetValue());
	ensure("Constructor2D node", node == SubGrid2D(i,j).GetNode());
    }
  }

  /*Test3: Constructor 3D*/
  template<>
  template<>
  void SubGridMesh_TestObject::test<3>()
  {
    Nx = 5; Ny = 10; Nz = 2;
    XYZ.push_back(Nx); XYZ.push_back(Ny); XYZ.push_back(Ny);
    SubGrid3D.SetMesh(XYZ);
    val=0.0;
    Node3D node(1,1,1);
    for (int i=0; i<Nx; i++)
      for (int j=0; j<Ny; j++)
	for (int k=0; k<Nz; k++){
	  ensure("Constructor3D solution", val == SubGrid3D(i,j,k).GetValue());
	  ensure("Constructor3D node", node == SubGrid3D(i,j,k).GetNode());
	}

  }

  /*Test4: Assign Values to the nodes of the Grid*/
  template<>
  template<>
  void SubGridMesh_TestObject::test<4>()
  {
    Nx = 5; Ny = 10; Nz = 2;
    XYZ.push_back(Nx); XYZ.push_back(Ny); XYZ.push_back(Ny);
    SubGrid3D.SetMesh(XYZ);
    val=0.0;
    Node3D node(1,1,1);
    for (int i=0; i<Nx; i++)
      for (int j=0; j<Ny; j++)
	for (int k=0; k<Nz; k++){
	  ensure("Constructor3D solution", val == SubGrid3D(i,j,k).GetValue());
	  ensure("Constructor3D node", node == SubGrid3D(i,j,k).GetNode());
	}

  }

  /*Test5: SetMesh(1D) */
  template<>
  template<>
  void SubGridMesh_TestObject::test<5>()
  {
    Nx = 5;			// number of nodes
    XYZ.push_back(Nx);
    // use first method for allocation
    SubGrid1D.SetMesh(XYZ);

    SubGridMesh<Node2D,OneD,double> SubGrid1D_Second; // second mesh
    // use second method for allocation
    SubGrid1D_Second.SetMesh(Nx);
    Node2D node(1,1);
    ensure("SetMesh() method", SubGrid1D == SubGrid1D_Second);

    // Change value
    Node2D node2(1.3, 3.1);
    SubGrid1D_Second(4).SetNode(node2);
    SubGrid1D_Second(4).SetValue(3.45545);
    ensure("SetMesh() method", !(SubGrid1D == SubGrid1D_Second));
  }


  /*Test6: SetMesh(2D) */
  template<>
  template<>
  void SubGridMesh_TestObject::test<6>()
  {
    Nx = 5; Ny=26;			// number of nodes
    XYZ.push_back(Nx); XYZ.push_back(Ny);
    // use first method for allocation
    SubGrid2D.SetMesh(XYZ);

    SubGridMesh<Node2D,TwoD,double> SubGrid2D_Second; // second mesh
    // use second method for allocation
    SubGrid2D_Second.SetMesh(Nx,Ny);
    ensure("SetMesh() method", SubGrid2D == SubGrid2D_Second);

    // Change value
    Node2D node2(1.3, 3.1);
    SubGrid2D_Second(3,17).SetNode(node2);
    SubGrid2D_Second(4,14).SetValue(3.45545);
    ensure("SetMesh() method", !(SubGrid2D == SubGrid2D_Second));
  }

  /*Test7: SetMesh(2D) */
  template<>
  template<>
  void SubGridMesh_TestObject::test<7>()
  {
    Nx = 5; Ny=26;			// number of nodes
    XYZ.push_back(Nx); XYZ.push_back(Ny);
    // use first method for allocation
    SubGrid2D.SetMesh(XYZ);

    SubGridMesh<Node2D,TwoD,double> SubGrid2D_Second; // second mesh
    // use second method for allocation
    SubGrid2D_Second.SetMesh(Nx,Ny);
    ensure("SetMesh() method", SubGrid2D == SubGrid2D_Second);

    // Change value using SetParam()
    Node2D node2(1.3, 3.1);
    SubGrid2D_Second(4,25).SetParam(node2,2.4545);
    SubGrid2D_Second(4,0).SetParam(node2,10.1);
    ensure("SetMesh() method", !(SubGrid2D == SubGrid2D_Second));
  }

  /*Test8: SetMesh(3D) */
  template<>
  template<>
  void SubGridMesh_TestObject::test<8>()
  {
    Nx = 2; Ny=3; Nz=2;	// number of nodes
    SubGrid3D.SetMesh(Nx,Ny,Nz);

    SubGridMesh<Node3D,ThreeD,double> SubGrid3D_Second; // second mesh
    XYZ.push_back(Nx); XYZ.push_back(Ny);  XYZ.push_back(Nz);
    SubGrid3D_Second.SetMesh(XYZ);

    //check method of allocation
    ensure("SetMesh() method", SubGrid3D == SubGrid3D_Second);

    // Change value using SetParam()
    Node3D node(1.3, 3.1, 3.4);
    SubGrid3D_Second(1,2,1).SetParam(node,0);
    SubGrid3D_Second(0,1,1).SetParam(node,10.1);
    SubGrid3D(0,1,1).SetParam(node,22.34);
    ensure("SetMesh() method", (SubGrid3D != SubGrid3D_Second));
  }

  /*Test9: cout */
  template<>
  template<>
  void SubGridMesh_TestObject::test<9>()
  {
    Nx = 2; Ny=3; Nz=2;	// number of nodes
    SubGrid3D.SetMesh(Nx,Ny,Nz);

    SubGridMesh<Node3D,ThreeD,double> SubGrid3D_Second; // second mesh
    XYZ.push_back(Nx); XYZ.push_back(Ny);  XYZ.push_back(Nz);
    SubGrid3D_Second.SetMesh(XYZ);

    // Change value using SetParam()
    Node3D node(1.3, 3.1, 3.4);
    SubGrid3D_Second(1,0,1).SetParam(node,0);
    SubGrid3D_Second(1,2,0).SetParam(node,10.1);
    //    cout << endl << SubGrid3D_Second << endl;
  }

  /*Test10: Copy Constructor */
  template<>
  template<>
  void SubGridMesh_TestObject::test<10>()
  {
    Nx = 2; Ny=3; Nz=2;	// number of nodes
    SubGrid3D.SetMesh(Nx,Ny,Nz);
    Node3D node(1.3, 3.1, 3.4);
    SubGrid3D(1,0,1).SetParam(node,0);
    SubGrid3D(1,2,0).SetParam(node,10.1);

    SubGridMesh<Node3D,ThreeD,double> SubGrid3D_Second = SubGrid3D; // second mesh
    ensure("Copy Constructor", SubGrid3D_Second == SubGrid3D);
    SubGrid3D(1,2,0).SetParam(node,1.1);
    ensure("Copy Constructor 2", SubGrid3D_Second!=SubGrid3D);
  }

  /*Test11: Assignment operator */
  template<>
  template<>
  void SubGridMesh_TestObject::test<11>()
  {
    Nx = 2; Ny=3; Nz=2;	// number of nodes
    SubGrid3D.SetMesh(Nx,Ny,Nz);
    Node3D node(1.3, 3.1, 3.4);
    SubGrid3D(1,0,1).SetParam(node,0);
    SubGrid3D(1,2,0).SetParam(node,10.1);

    SubGridMesh<Node3D,ThreeD,double> SubGrid3D_Second; // second mesh

    SubGrid3D_Second = SubGrid3D;
    ensure("Assignment Operator", SubGrid3D_Second == SubGrid3D);
  }

  /*Test12: OutputGeometryTecplot() */
  template<>
  template<>
  void SubGridMesh_TestObject::test<12>()
  {
    Nx = 5; Ny=2;			// number of nodes
    XYZ.push_back(Nx); XYZ.push_back(Ny);
    // use first method for allocation
    SubGrid2D.SetMesh(XYZ);

    // Change value using SetParam()
    Node2D node2(1.3, 3.1);
    SubGrid2D(4,1).SetParam(node2,2.4545);
    SubGrid2D(4,0).SetParam(node2,10.1);

    // output file

    // Set Path
    ofstream output_file;
    char *output_file_name_ptr;
    char output_file_name[256];
    output_file_name_ptr = output_file_name;

    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"SubgridTecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    SubGrid2D.OutputGeometryTecplot(output_file,true);
    output_file.close();
  }


  /*Test13: OutputGeometryTecplot() */
  template<>
  template<>
  void SubGridMesh_TestObject::test<13>()
  {

    Nx = 2; Ny=3; Nz=2;	// number of nodes
    SubGrid3D.SetMesh(Nx,Ny,Nz);
    Node3D node(1.3, 3.1, 3.4);
    SubGrid3D(1,0,1).SetParam(node,0);
    SubGrid3D(1,2,0).SetParam(node,10.1);

    // output file

    // Set Path
    ofstream output_file;
    char *output_file_name_ptr;
    char output_file_name[256];
    output_file_name_ptr = output_file_name;
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"SubgridTecplot3D.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    SubGrid3D.OutputGeometryTecplot(output_file);
    output_file.close();

    HeaderData Names;
    Names.add("Var1");
    // SubGrid3D.DefineHeader(Names);
    //    SubGrid3D.VarNames.add("Var1");
    SubGrid3D.VarNames = State::PrintedVariables;
    cout << SubGrid3D.VarNames << endl;

  }

}

namespace tut
{
  SubGridMesh_TestGroup SubGridMesh_Test("SubGridMesh_Test");
}
