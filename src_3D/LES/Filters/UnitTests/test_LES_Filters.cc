/* Include required CFFC header files. */


#include "TestData.h"
#include "../LES_Filters.h"

// Change this file to your
#include "../../LES_Polytropic/LES3DPolytropic.h"


/* --------------- copied from HexaBlock.h ---------------- */
#ifndef _HEXA_MULTIBLOCK_INCLUDED
#include "../../../HexaBlock/HexaMultiBlock.h"
#endif //_HEXA_MULTIBLOCK_INCLUDED

#ifndef _AMR_INCLUDED
#include "../../../AMR/AMR.h"
#endif // _AMR_INCLUDED

#ifndef _ADAPTIVEBLOCK3D_MESSAGEPASSING_INCLUDED
#include "../../../AMR/AdaptiveBlock3D_MessagePassing.h"
#endif // _ADAPTIVEBLOCK3D_MESSAGEPASSING_INCLUDED

#ifndef _HEXA_SOLVER_CLASSES_INCLUDED
#include  "../../../HexaBlock/HexaSolverClasses.h"
#endif  //_HEXA_SOLVER_CLASSES_INCLUDED

#ifndef _HEXA_PRE_PROCESSING_INCLUDED
#include "../../../HexaBlock/HexaPreProcessing.h"
#endif //_HEXA_PRE_PROCESSING_INCLUDED

#ifndef _HEXA_POST_PROCESSING_INCLUDED
#include "../../../HexaBlock/HexaPostProcessing.h"
#endif //_HEXA_POST_PROCESSING_INCLUDED

#ifndef _HEXA_EXPLICIT_SOLVER
#include "../../../HexaBlock/HexaExplicitSolver.h"
#endif //_HEXA_EXPLICIT_SOLVER

#ifndef _NKS_INCLUDED
#include "../../../NewtonKrylovSchwarz/NKS.h"
#endif //_NKS_INCLUDED
/* -------------------------------------------------------- */

namespace tut
{
    
    /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
    class Data_LES_Filters : public TestData {
        
        // Local variables
    public:
        typedef LES3D_Polytropic_pState Soln_pState;
        typedef LES3D_Polytropic_cState Soln_cState;
        
        int error_flag;
        int batch_flag;
        char* Input_File_Name_ptr;
        // Non Solution Specific
        HexaSolver_Data                                     Data;  
        // Solution Specific 
        HexaSolver_Solution_Data<Soln_pState, Soln_cState>  Solution_Data; 
        
        // Constructor
        Data_LES_Filters(){
            /* Set the global path for this test suite. 
             It's a good practice to put it in the constructor of the data class in order to be set
             automatically for each individual test. Declare it relative to the /src_3D directory,
             otherwise the framework might not find the input and output files. */
            set_test_suite_path("LES/Filters/UnitTests");
            
            
            
            batch_flag = 1;
            error_flag = false;
            Input_File_Name_ptr = strcat(Global_TestSuite_Path, "LES_Filters_input.in");
            
            

        }
        
        
        
    private:
        
    };
    
    
    
    /**
     * This group of declarations is just to register
     * test group in test-application-wide singleton.
     * Name of test group object (e.g. 'NAME'_TestSuite) shall
     * be unique in tut:: namespace. Alternatively, you
     * you may put it into anonymous namespace.
     */
    typedef test_group<Data_LES_Filters> LES_Filters_TestSuite;
    typedef LES_Filters_TestSuite::object LES_Filters_object;
    
    
    /******************************************************************************************************
     ******************************************************************************************************
     *                      ********    *******        ***       ********        ***                      *
     *                         **       **           **   **        **         **   **                    *
     *                         **       **          **              **        **                          *
     *                         **       *******      *****          **          *****                     *
     *                         **       **                **        **               **                   *
     *                         **       **            *    **       **          *   **                    *
     *                         **       *******        ****         **           ****                     *
     ******************************************************************************************************
     ******************************************************************************************************/
    
    
    /*********************************************
     TESTING OPERATIONS:
     -->  ensure("ConditionName", test_condition);
     -->  ensure_not("ConditionName", test_condition);
     -->  ensure_equals("ConditionName", expected_value, got_value);
     // operators '<<' and '!=' must be defined!!!
     -->  ensure_distance("ConditionName", expected_value, got_value, tolerance);
     // operators '<<', '>=' and '<=' must be defined!!!
     -->  fail("Message");
     
     Obs: "ConditionName" is optional
     */
    
    

    /* Test 1:*/
    template<>
    template<>
    void LES_Filters_object::test<1>()
    {
        
        set_test_name("Test Haselbacher functions");
        
        Haselbacher_Filter<Soln_pState,Soln_cState> filter;
        
        cout << endl << endl << endl;
        ensure("fac(0)",filter.fac(0)==1);
        ensure("fac(1)",filter.fac(1)==1);
        ensure("fac(2)",filter.fac(2)==2);
        ensure("fac(3)",filter.fac(3)==6);
        
        ensure("trinomial(0,0,0)",filter.trinomial_coefficient(0,0,0)==1);
        ensure("trinomial(0,0,0)",filter.trinomial_coefficient(1,2,1)==12);
        ensure("trinomial(0,0,0)",filter.trinomial_coefficient(3,2,1)==60);


        
    }


    /* Test 2:*/
    template<>
    template<>
    void LES_Filters_object::test<2>()
    {
        
        set_test_name("Test transfer_function");
        
        
        Data.batch_flag=batch_flag;
        error_flag = Solution_Data.Get_Input_Parameters(Input_File_Name_ptr, batch_flag);  
        ensure("Get Input Parameters",error_flag==false);
        CFFC_Barrier_MPI();
        error_flag = Initialize_Solution_Blocks(Data,Solution_Data);
        ensure("Initialize_Solution_Blocks",error_flag==false);
        error_flag = Initial_Conditions(Data,Solution_Data);      
        ensure("Initial_Conditions",error_flag==false);
        
        
        
        LES_Filter<Soln_pState,Soln_cState> myfilter(Data,Solution_Data,LES_FILTER_HASELBACHER);
        
        myfilter.transfer_function();

        
    }
    
    
    
    /* Test 4:*/
    template<>
    template<>
    void LES_Filters_object::test<4>()
    {
        
        set_test_name("Test DiagonalMatrix");
        
        
        DiagonalMatrix D(3);
        DenseMatrix M(3,3);
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                M(i,j) = (i+1)*(j+1);
            }
            D(i)=i+1;
        }
        cout << "D = " << endl << D << endl;
        cout << "M = " << endl << M << endl;
        cout << "M*D = " << endl << M*D << endl;
        cout << "D*M = " << endl << D*M << endl;
        
    }
    
    /* Test 3:*/
    template<>
    template<>
    void LES_Filters_object::test<3>()
    {
        
        set_test_name("Test Filter");
        
        
        Data.batch_flag=batch_flag;
        error_flag = Solution_Data.Get_Input_Parameters(Input_File_Name_ptr, batch_flag);  
        ensure("Get Input Parameters",error_flag==false);
        CFFC_Barrier_MPI();
        error_flag = Initialize_Solution_Blocks(Data,Solution_Data);
        ensure("Initialize_Solution_Blocks",error_flag==false);
        error_flag = Initial_Conditions(Data,Solution_Data);      
        ensure("Initial_Conditions",error_flag==false);
        
        
        
        
        cout << endl<< endl << endl << "FILTERING..." << endl;
        LES_Filter<Soln_pState,Soln_cState> myfilter(Data,Solution_Data,LES_FILTER_HASELBACHER);
        myfilter.filter();
        
        error_flag = Hexa_Post_Processing(Data,Solution_Data);
        ensure("Post_Processing",error_flag==false);
    
    }
    


}





// Test suite constructor
tut::LES_Filters_TestSuite LES_Filters_TestSuite("Class:LES_Filters");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"
 
 e.g. "Class:NameOfTheClass"
 "Template Class:NameOfTheTemplateClass [&& type used for testing]"
 "Integration:NameOfTheMethod"
 "Linear Systems Solvers:NameOfTheMethod"
 
 */
