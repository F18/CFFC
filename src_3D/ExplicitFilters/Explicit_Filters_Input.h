/*
 *  Explicit_Filters_Input.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 18/04/09.
 *
 */

//        MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);


/* TurbulenceModellingInput.h:  Header file declaring turbulence model
 input class. */

#ifndef _EXPLICIT_FILTERS_INPUT_INCLUDED
#define _EXPLICIT_FILTERS_INPUT_INCLUDED

// Include required C++ libraries.

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <vector>

using namespace std;

// Include CFD and MPI header files.

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _EXPLICIT_FILTER_CONSTANTS_INCLUDED
#include "../ExplicitFilters/Explicit_Filter_Constants.h"
#endif // _EXPLICIT_FILTER_CONSTANTS_INCLUDED

#define EXPLICIT_FILTERSL_INPUT_PARAMETER_LENGTH 256
#define NUMBER_OF_FILTERS 2

using namespace Explicit_Filter_Constants;

/*!
 * Class: Turbulence_Modelling_Input_Parameters
 *
 * @brief Input Parameters for the modelling of turbulence flows.
 *
 * This class defines and handles the input variables related to the
 * modelling and treatment of turbulent flows.
 *
 */
class Explicit_Filters_Input_Parameters{
private:
public:
    // operations
    int Filter_Initial_Condition;           //!< Flag defines to filter the initial condition
    int Solution_Filtering_Frequency;       //!< The solution is filtered (with secondary filter) each time after this number of timesteps has passed.
    int Filter_Solution_Before_Execution;   //!< Define if the solution is filtered before execution with secondary filter
    
    // general input
    int Commutation_Order;              //!< Commutation order of explicit filter
    int Filter_Method;                  //!< Chooses either to filter residuals or solutions
    int Finite_Differencing_Order;      //!< Finite differencing order in commutation error calculations

    // type independent
    std::vector<int>    Filter_Type;                 //!< Filter type : default = implicit filtering
    std::vector<double> FGR;                         //!< Filter width to mesh size ratio, default : $f \bar{\Delta} = 2 \Delta x $f
    std::vector<double> Filter_Width;                //!< Constant filter width
    std::vector<int>    Number_Of_Rings;             //!< Number of rings used in the explicit filter
    std::vector<int>    Filter_Memory_Efficient;     //!< Determines whether to store filter weights or not (not storing slows down tremendously)
    std::vector<double> Filter_Strength;             //!< Defines how strong the filter will act (from 0 to 1)
    std::vector<int>    Filter_Relative;             //!< Defines if the filtering happens relative to the mesh spacing
    std::vector<int>    Generate_At_Startup;         //!< Sets whether the filter weights will be calculated on the go or during Initialization

    // Vasilyev
    std::vector<double> Target_Filter_Sharpness;     //!< Sharpness of target filter in case of Least squares constraint (vasilyev)
    std::vector<int>    Filter_Width_Strict;         //!< This will strictly set the FGR and won't allow Least squares to approximate : default = false
    std::vector<int>    LS_Constraints;              //!< This will turn on or off Least squares constraints for Vasilyev's filter : default = true
    std::vector<int>    Derivative_Constraints;      //!< Determines the number of derivative constraints for Vasilyev's filter : default = true
    
    // Haselbacher
    std::vector<double> Relaxation_Factor;                      //!< This coefficient is used in the least-squares reconstruction filter and should be left "DEFAULT".
    std::vector<int>    Least_Squares_Filter_Weighting;         //!< This coefficient controls if a weighted least-squares is used in Haselbacher filter.
    std::vector<double> Least_Squares_Filter_Weighting_Factor;  //!< This coefficient controls if a weighted least-squares is used in Haselbacher filter.
    std::vector<int>    Reconstruction_Type;                    //!< This flag defines whether to use the CENO geometric coefficients or standard LS
    //@{ @name Constructors and desctructors:
    //! Constructor (assign default values)
    Explicit_Filters_Input_Parameters();    
    //! Destructor
    ~Explicit_Filters_Input_Parameters(void){ }
    //@}
    
    //@{ @name Other Member functions:
    //! Broadcast input parameters to all processors
    void Broadcast(void);
    //! Parse next input line
    int Parse_Next_Input_Control_Parameter(char *code, stringstream &value);
    //! Check validity of specified input parameters
    int Check_Inputs(void);
    //@}
    
    //@{ @name Input-output operators:
    friend ostream &operator << (ostream &out_file,
                                 const Explicit_Filters_Input_Parameters &IP);
    friend istream &operator >> (istream &in_file,
                                 Explicit_Filters_Input_Parameters &IP);
    void Output(ostream &out_file) const;
    //@}
    
};






#endif // _EXPLICIT_FILTERS_INPUT_INCLUDED