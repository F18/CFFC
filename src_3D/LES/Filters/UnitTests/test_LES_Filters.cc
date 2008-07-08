/* Include required CFFC header files. */


#include "TestData.h"
#include "../LES_Filters.h"
#include "../../../TurbulenceModelling/SpectralAnalysis.h"

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
        
        void Initialize(void) {
            Data.batch_flag=batch_flag;
            error_flag = Solution_Data.Get_Input_Parameters(Input_File_Name_ptr, batch_flag);  
            ensure("Get Input Parameters",error_flag==false);
            CFFC_Barrier_MPI();
            error_flag = Initialize_Solution_Blocks(Data,Solution_Data);
            ensure("Initialize_Solution_Blocks",error_flag==false);
            error_flag = Initial_Conditions(Data,Solution_Data);      
            ensure("Initial_Conditions",error_flag==false);
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
    void LES_Filters_object::test<3>()
    {
        
        set_test_name("Test Filter");
        
        

        Initialize();
        
        SpectralAnalysis<Soln_pState,Soln_cState> Spectrum(Data,Solution_Data);
        
        double J = Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.jacobian(5,8,8,2);
        cout << "J = " << J << "   Delta = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[5][8][8].dXc << endl;
         J = Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.jacobian(6,8,8,2);
        cout << "J = " << J << "   Delta = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[6][8][8].dXc<< endl;
         J = Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.jacobian(7,8,8,2);
        cout << "J = " << J << "   Delta = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[7][8][8].dXc<< endl;
         J = Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.jacobian(8,8,8,2);
        cout << "J = " << J << "   Delta = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[8][8][8].dXc<< endl;
         J = Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.jacobian(9,8,8,2);
        cout << "J = " << J << "   Delta = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[9][8][8].dXc<< endl;
         J = Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.jacobian(10,8,8,2);
        cout << "J = " << J << "   Delta = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[10][8][8].dXc<< endl;
         J = Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.jacobian(11,8,8,2);
        cout << "J = " << J << "   Delta = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[11][8][8].dXc<< endl;
        //cout << "J = " << J << endl;
            
        
        typedef double (Soln_pState::*member_ptr);
        typedef Soln_cState *** (Hexa_Block<Soln_pState,Soln_cState>::*Soln_cState_3D_ptr_type);

        Soln_cState_3D_ptr_type U_ptr = &Hexa_Block<Soln_pState,Soln_cState>::U; 

        member_ptr rho_member = &Soln_pState::rho;
//        member_ptr p_member = &Soln_pState::p;
//
//        Spectrum.Set_Spectrum(p_member);
        Spectrum.Set_Spectrum(rho_member);
//                
//        //cout << endl<< endl << endl << "FILTERING..." << endl;
        LES_Filter<Soln_pState,Soln_cState> myfilter(Data,Solution_Data);
        //LES_Filter<Soln_pState,Soln_cState> myfilter(Data,Solution_Data,LES_FILTER_HASELBACHER);

//        Soln_cState **** (Hexa_Block<Soln_pState,Soln_cState>::*dUdt_ptr) = NULL;
//        dUdt_ptr = &Hexa_Block<Soln_pState,Soln_cState>::dUdt;
//        myfilter.Set_Filter_Variables(dUdt_ptr,0);
//        myfilter.test();
       // myfilter.reset();   // make sure that we will use the new filter
        
       // Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Update_Grid_Cells();
        cout << "J(0,0,0) = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[0][0][0].Jacobian << endl;
        cout << "J(0,0,10) = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[0][0][10].Jacobian << endl;
        cout << "J(0,10,10) = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[0][10][10].Jacobian << endl;
        cout << "J(10,10,10) = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[10][10][10].Jacobian << endl;

        cout << "J(5,8,8) = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[5][8][8].Jacobian << endl;
        cout << "J(6,8,8) = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[6][5][8].Jacobian << endl;
        cout << "J(7,8,8) = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[7][8][8].Jacobian << endl;
        cout << "J(8,8,8) = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[8][8][8].Jacobian << endl;
        cout << "J(9,8,8) = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[9][8][8].Jacobian << endl;
        cout << "J(10,8,8) = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[10][8][8].Jacobian << endl;
        cout << "J(11,8,8) = " << Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Grid.Cell[11][8][8].Jacobian << endl;

//        myfilter.test();
//        myfilter.transfer_function(FILTER_CORNER_CELL);
//        myfilter.transfer_function(FILTER_FACE_CELL);
//        myfilter.transfer_function(FILTER_EDGE_CELL);
//        myfilter.transfer_function(FILTER_INNER_CELL);
//        myfilter.transfer_function(FILTER_MIDDLE_CELL);

        
//        cout << "start filtering" << endl;
//        myfilter.filter(rho_member);
        

//        cout << "done filtering" << endl;
//        Spectrum.Get_Spectrum(rho_member,"density");
        
        //if (Solution_Data.Input.Turbulence_IP.i_filter_type != FILTER_TYPE_RESTART)
        //    myfilter.Write_to_file();

        

        
        cout << " Calculating commutation error " << endl;

        myfilter.Commutation_Error(rho_member);
        
//        Spectrum.Get_Spectrum(p_member,"pressure");
//
        error_flag = Hexa_Post_Processing(Data,Solution_Data);
//        ensure("Post_Processing",error_flag==false);
    
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
