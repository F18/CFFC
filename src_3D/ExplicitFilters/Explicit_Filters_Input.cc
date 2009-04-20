/*
 *  Explicit_Filters_Input.cc
 *  CFFC
 *
 *  Created by Willem Deconinck on 18/04/09.
 *
 */


#include "Explicit_Filters_Input.h"

using namespace Explicit_Filter_Constants;

/************************************************************************************
 *                       Constructor                                                *
 ************************************************************************************/
Explicit_Filters_Input_Parameters::Explicit_Filters_Input_Parameters() :
    Filter_Type(NUMBER_OF_FILTERS),
    FGR(NUMBER_OF_FILTERS),
    Filter_Width(NUMBER_OF_FILTERS),
    Number_Of_Rings(NUMBER_OF_FILTERS),
    Filter_Memory_Efficient(NUMBER_OF_FILTERS),
    Filter_Relative(NUMBER_OF_FILTERS),
    Filter_Strength(NUMBER_OF_FILTERS),
    Relaxation_Factor(NUMBER_OF_FILTERS),
    Least_Squares_Filter_Weighting(NUMBER_OF_FILTERS),
    Least_Squares_Filter_Weighting_Factor(NUMBER_OF_FILTERS),
    Target_Filter_Sharpness(NUMBER_OF_FILTERS),
    Filter_Width_Strict(NUMBER_OF_FILTERS),
    LS_Constraints(NUMBER_OF_FILTERS),
    Derivative_Constraints(NUMBER_OF_FILTERS)
{
    
    // ============================================================== //
    // Explicit Filter Operations
    // ============================================================== //
    
    Filter_Initial_Condition = ON;
    Solution_Filtering_Frequency = 0;
    Filter_Solution_Before_Execution = OFF;
    
    // ============================================================== //
    // General Explicit Filter Input (doesn't depend on which one)
    // ============================================================== //
    
    Commutation_Order = 2;
    Filter_Method = FILTER_RESIDUALS;
    Finite_Differencing_Order = Commutation_Order + 1;
    
    // ============================================================== //
    // Type Independent Input
    // ============================================================== //
    
    Filter_Type[PRIMARY_FILTER]   = IMPLICIT_FILTER;
    Filter_Type[SECONDARY_FILTER] = IMPLICIT_FILTER;
    
    FGR[PRIMARY_FILTER]   = TWO;
    FGR[SECONDARY_FILTER] = TWO;
    
    Filter_Width[PRIMARY_FILTER]   = ZERO;
    Filter_Width[SECONDARY_FILTER] = ZERO;
    
    Number_Of_Rings[PRIMARY_FILTER]   = 2;
    Number_Of_Rings[SECONDARY_FILTER] = 2;
    
    Filter_Memory_Efficient[PRIMARY_FILTER]   = OFF;
    Filter_Memory_Efficient[SECONDARY_FILTER] = OFF;
    
    Filter_Relative[PRIMARY_FILTER] = ON;
    Filter_Relative[SECONDARY_FILTER] = ON;
    
    Filter_Strength[PRIMARY_FILTER]   = 1.0;
    Filter_Strength[SECONDARY_FILTER] = 1.0;        
    
    // ============================================================== //
    // Haselbacher Filter Input
    // ============================================================== //
    
    Relaxation_Factor[PRIMARY_FILTER]   = DEFAULT;
    Relaxation_Factor[SECONDARY_FILTER] = DEFAULT;
    
    Least_Squares_Filter_Weighting_Factor[PRIMARY_FILTER]   = DEFAULT;
    Least_Squares_Filter_Weighting_Factor[SECONDARY_FILTER] = DEFAULT;
    
    Least_Squares_Filter_Weighting[PRIMARY_FILTER]   = ON;
    Least_Squares_Filter_Weighting[SECONDARY_FILTER] = ON;
    
    // ============================================================== //
    // Vasilyev Filter Input
    // ============================================================== //
    
    Target_Filter_Sharpness[PRIMARY_FILTER]   = DEFAULT;
    Target_Filter_Sharpness[SECONDARY_FILTER] = DEFAULT;
    
    Filter_Width_Strict[PRIMARY_FILTER]   = OFF;
    Filter_Width_Strict[SECONDARY_FILTER] = OFF;
    
    LS_Constraints[PRIMARY_FILTER]   = ON;
    LS_Constraints[SECONDARY_FILTER] = ON;
    
    Derivative_Constraints[PRIMARY_FILTER]   = DEFAULT; 
    Derivative_Constraints[SECONDARY_FILTER] = DEFAULT;
}


/********************************************************************************************
 * Explicit_Filters_Input_Parameters::Parse_Next_Input_Control_Parameter - Parse input. *
 ********************************************************************************************/
int Explicit_Filters_Input_Parameters::Parse_Next_Input_Control_Parameter(char *code, 
                                                                              stringstream &value) {
    
    // Returns:
    //  - INVALID_INPUT_VALUE if code is valid but value is invalid
    //  - INVALID_INPUT_CODE  if unknown code
    
    int i_command = INVALID_INPUT_CODE;
    string value_string;
            
    // ============================================================== //
    // Explicit Filter Operations
    // ============================================================== //
    
    if (strcmp(code, "ExplicitFilter[1].Filter_Initial_Condition") == 0) {
        i_command = 132;
        value >> value_string;
        if (strcmp(value_string.c_str(), "ON") == 0) {
            Filter_Initial_Condition = ON;
        } else if(strcmp(value_string.c_str(), "OFF") == 0) {
            Filter_Initial_Condition = OFF; 
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    } 
    
    else if (strcmp(code, "ExplicitFilter[2].Solution_Filtering_Frequency") == 0) {
        i_command = 132;
        value >> Solution_Filtering_Frequency;
    } 
    
    else if (strcmp(code, "ExplicitFilter[2].Filter_Solution_Before_Execution") == 0) {
        i_command = 143;
        value >> value_string;
        if      (value_string == "ON")   Filter_Solution_Before_Execution = ON;
        else if (value_string == "OFF")  Filter_Solution_Before_Execution = OFF;
        else 
            i_command = INVALID_INPUT_VALUE;
    } 
    
    // ============================================================== //
    // General Explicit Filter Input (doesn't depend on which one)
    // ============================================================== //
    
    else if (strcmp(code, "ExplicitFilter.Commutation_Order") == 0) {
        i_command = 131;
        value >> Commutation_Order;
        if (Commutation_Order < 1)
            i_command = INVALID_INPUT_VALUE;
    }
    
    else if (strcmp(code,"ExplicitFilter.Method") == 0) {
        i_command = 143;
        value >> value_string;
        if      (value_string == "Variables"){
            Filter_Method = FILTER_VARIABLES;
            Filter_Strength[PRIMARY_FILTER] = 0.2;
            Filter_Strength[SECONDARY_FILTER] = 0.0;
        }
        else if (value_string == "Residuals") {
            Filter_Method = FILTER_RESIDUALS;
            Filter_Strength[PRIMARY_FILTER] = 1.0;
            Filter_Strength[SECONDARY_FILTER] = 1.0;
        }
        else
            i_command = INVALID_INPUT_VALUE;
    } 
    
    else if (strcmp(code, "ExplicitFilter.Finite_Differencing_Order") == 0) {
        i_command = 131;
        value >> Finite_Differencing_Order;
        if (Finite_Differencing_Order < 1)
            i_command = INVALID_INPUT_VALUE;
    }
    
    // ============================================================== //
    // Type Independent Input
    // ============================================================== //
    
    else if (strcmp(code, "ExplicitFilter[1].Type") == 0) {
        i_command = 130;
        value >> value_string;
        if (value_string == "Implicit") {
            Filter_Type[PRIMARY_FILTER] = IMPLICIT_FILTER;
        } else if (value_string == "Haselbacher") {
            Filter_Type[PRIMARY_FILTER] = HASELBACHER_FILTER;
        } else if (value_string == "Vasilyev") {
            Filter_Type[PRIMARY_FILTER] = VASILYEV_FILTER;
        } else if (value_string == "Restart") {
            Filter_Type[PRIMARY_FILTER] = RESTART_FILTER;
        } else {
            i_command = INVALID_INPUT_VALUE;
        } /* endif */
    } 
    else if (strcmp(code, "ExplicitFilter[2].Type") == 0) {
        i_command = 130;
        value >> value_string;
        if (value_string == "Implicit") {
            Filter_Type[SECONDARY_FILTER] = IMPLICIT_FILTER;
        } else if (value_string == "Haselbacher") {
            Filter_Type[SECONDARY_FILTER] = HASELBACHER_FILTER;
        } else if (value_string == "Vasilyev") {
            Filter_Type[SECONDARY_FILTER] = VASILYEV_FILTER;
        } else if (value_string == "Restart") {
            Filter_Type[SECONDARY_FILTER] = RESTART_FILTER;
        } else {
            i_command = INVALID_INPUT_VALUE;
        } /* endif */
    } 
    
    else if (strcmp(code, "ExplicitFilter[1].Filter_Memory_Efficient") == 0) {
        i_command = 132;
        value >> value_string;
        if (strcmp(value_string.c_str(), "ON") == 0) {
            Filter_Memory_Efficient[PRIMARY_FILTER] = ON;
        } else if(strcmp(value_string.c_str(), "OFF") == 0) {
            Filter_Memory_Efficient[PRIMARY_FILTER] = OFF; 
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    } 
    else if (strcmp(code, "ExplicitFilter[2].Filter_Memory_Efficient") == 0) {
        i_command = 132;
        value >> value_string;
        if (strcmp(value_string.c_str(), "ON") == 0) {
            Filter_Memory_Efficient[SECONDARY_FILTER] = ON;
        } else if(strcmp(value_string.c_str(), "OFF") == 0) {
            Filter_Memory_Efficient[SECONDARY_FILTER] = OFF; 
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    }
    
    else if (strcmp(code, "ExplicitFilter[1].FGR") == 0) {
        i_command = 131;
        value >> FGR[PRIMARY_FILTER];
        Filter_Relative[PRIMARY_FILTER] = true;
        if (FGR[PRIMARY_FILTER] < 1)
            i_command = INVALID_INPUT_VALUE;
    }
    else if (strcmp(code, "ExplicitFilter[2].FGR") == 0) {
        i_command = 131;
        value >> FGR[SECONDARY_FILTER];
        Filter_Relative[SECONDARY_FILTER] = true;
        if (FGR[SECONDARY_FILTER] < 1)
            i_command = INVALID_INPUT_VALUE;
    }
    
    else if (strcmp(code, "ExplicitFilter[1].Filter_Width") == 0) {
        i_command = 132;
        value >> Filter_Width[PRIMARY_FILTER];
        Filter_Relative[PRIMARY_FILTER] = false;
        if ( Filter_Width[PRIMARY_FILTER] < 0.0 )
            i_command = INVALID_INPUT_VALUE;
    } 
    else if (strcmp(code, "ExplicitFilter[2].Filter_Width") == 0) {
        i_command = 132;
        value >> Filter_Width[SECONDARY_FILTER];
        Filter_Relative[SECONDARY_FILTER] = false;
        if ( Filter_Width[SECONDARY_FILTER] < 0.0 )
            i_command = INVALID_INPUT_VALUE;
    } 
    
    else if (strcmp(code, "ExplicitFilter[1].Filter_Relative") == 0) {
        i_command = 131;
        value >> value_string;
        if (strcmp(value_string.c_str(), "ON") == 0) {
            Filter_Relative[PRIMARY_FILTER] = ON;
        } else if(strcmp(value_string.c_str(), "OFF") == 0) {
            Filter_Relative[PRIMARY_FILTER] = OFF; 
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    }
    else if (strcmp(code, "ExplicitFilter[2].Filter_Relative") == 0) {
        i_command = 131;
        value >> value_string;
        if (strcmp(value_string.c_str(), "ON") == 0) {
            Filter_Relative[SECONDARY_FILTER] = ON;
        } else if(strcmp(value_string.c_str(), "OFF") == 0) {
            Filter_Relative[SECONDARY_FILTER] = OFF; 
        } else {
            i_command = INVALID_INPUT_VALUE;
        }    }
    
    else if (strcmp(code, "ExplicitFilter[1].Number_Of_Rings") == 0) {
        i_command = 131;
        value >> Number_Of_Rings[PRIMARY_FILTER];
        if (Number_Of_Rings[PRIMARY_FILTER] < 1)
            i_command = INVALID_INPUT_VALUE;
    }
    else if (strcmp(code, "ExplicitFilter[2].Number_Of_Rings") == 0) {
        i_command = 131;
        value >> Number_Of_Rings[SECONDARY_FILTER];
        if (Number_Of_Rings[SECONDARY_FILTER] < 1)
            i_command = INVALID_INPUT_VALUE;                    
    } 
    
    else if (strcmp(code,"ExplicitFilter[1].Strength") == 0) {
        i_command = 143;
        value >> Filter_Strength[PRIMARY_FILTER];
    } 
    else if (strcmp(code,"ExplicitFilter[2].Strength") == 0) {
        i_command = 143;
        value >> Filter_Strength[PRIMARY_FILTER];        
    } 
    
    // ============================================================== //
    // Vasilyev Filter Input
    // ============================================================== //
    
    else if (strcmp(code, "ExplicitFilter[1].Vasilyev.Target_Filter_Sharpness") == 0) {
        i_command = 132;
        value >> Target_Filter_Sharpness[PRIMARY_FILTER];
    } 
    else if (strcmp(code, "ExplicitFilter[2].Vasilyev.Target_Filter_Sharpness") == 0) {
        i_command = 132;
        value >> Target_Filter_Sharpness[SECONDARY_FILTER];
    }
    
    else if (strcmp(code, "ExplicitFilter[1].Vasilyev.Filter_Width_Strict") == 0) {
        i_command = 132;
        value >> value_string;
        if (strcmp(value_string.c_str(), "ON") == 0) {
            Filter_Width_Strict[PRIMARY_FILTER] = ON;
        } else if(strcmp(value_string.c_str(), "OFF") == 0) {
            Filter_Width_Strict[PRIMARY_FILTER] = OFF; 
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    } 
    else if (strcmp(code, "ExplicitFilter[2].Vasilyev.Filter_Width_Strict") == 0) {
        i_command = 132;
        value >> value_string;
        if (strcmp(value_string.c_str(), "ON") == 0) {
            Filter_Width_Strict[SECONDARY_FILTER] = ON;
        } else if(strcmp(value_string.c_str(), "OFF") == 0) {
            Filter_Width_Strict[SECONDARY_FILTER] = OFF; 
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    } 
    
    else if (strcmp(code, "ExplicitFilter[1].Vasilyev.LS_Constraints") == 0) {
        i_command = 132;
        value >> value_string;
        if (strcmp(value_string.c_str(), "ON") == 0) {
            LS_Constraints[PRIMARY_FILTER] = ON;
        } else if(strcmp(value_string.c_str(), "OFF") == 0) {
            LS_Constraints[PRIMARY_FILTER] = OFF; 
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    } 
    else if (strcmp(code, "ExplicitFilter[2].Vasilyev.LS_Constraints") == 0) {
        i_command = 132;
        value >> value_string;
        if (strcmp(value_string.c_str(), "ON") == 0) {
            LS_Constraints[SECONDARY_FILTER] = ON;
        } else if(strcmp(value_string.c_str(), "OFF") == 0) {
            LS_Constraints[SECONDARY_FILTER] = OFF; 
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    } 
    
    else if (strcmp(code, "ExplicitFilter[1].Vasilyev.Derivative_Constraints") == 0) {
        i_command = 132;
        value >> Derivative_Constraints[PRIMARY_FILTER];
    } 
    else if (strcmp(code, "ExplicitFilter[2].Vasilyev.Derivative_Constraints") == 0) {
        i_command = 132;
        value >> Derivative_Constraints[SECONDARY_FILTER];
    } 
    
    // ============================================================== //
    // Haselbacher Filter Input
    // ============================================================== //
    
    else if (strcmp(code, "ExplicitFilter[1].Haselbacher.Least_Squares_Filter_Weighting") == 0) {
        i_command = 132;
        value >> value_string;
        if (strcmp(value_string.c_str(), "ON") == 0) {
            Least_Squares_Filter_Weighting[PRIMARY_FILTER] = ON;
        } else if(strcmp(value_string.c_str(), "OFF") == 0) {
            Least_Squares_Filter_Weighting[PRIMARY_FILTER] = OFF; 
        } else {
            i_command = INVALID_INPUT_VALUE; 
        }
    } 
    else if (strcmp(code, "ExplicitFilter[2].Haselbacher.Least_Squares_Filter_Weighting") == 0) {
        i_command = 132;
        value >> value_string;
        if (strcmp(value_string.c_str(), "ON") == 0) {
            Least_Squares_Filter_Weighting[PRIMARY_FILTER] = ON;
        } else if(strcmp(value_string.c_str(), "OFF") == 0) {
            Least_Squares_Filter_Weighting[PRIMARY_FILTER] = OFF; 
        } else {
            i_command = INVALID_INPUT_VALUE; 
        }
    } 
    
    else if (strcmp(code, "ExplicitFilter[1].Haselbacher.Least_Squares_Filter_Weighting_Factor") == 0) {
        i_command = 131;
        value >> value_string;
        if (strcmp(value_string.c_str(), "DEFAULT") == 0) {
            Least_Squares_Filter_Weighting_Factor[PRIMARY_FILTER] = DEFAULT;
        } else {
            std::stringstream valuestream;
            valuestream << value_string;
            valuestream >> Least_Squares_Filter_Weighting_Factor[PRIMARY_FILTER];
        }
    } 
    else if (strcmp(code, "ExplicitFilter[2].Haselbacher.Least_Squares_Filter_Weighting_Factor") == 0) {
        i_command = 131;
        value >> value_string;
        if (strcmp(value_string.c_str(), "DEFAULT") == 0) {
            Least_Squares_Filter_Weighting_Factor[SECONDARY_FILTER] = DEFAULT;
        } else {
            std::stringstream valuestream;
            valuestream << value_string;
            valuestream >> Least_Squares_Filter_Weighting_Factor[SECONDARY_FILTER];
        }
    } 
    
    else if (strcmp(code, "ExplicitFilter[1].Haselbacher.Relaxation_Factor") == 0) {
        i_command = 131;
        value >> value_string;
        if (strcmp(value_string.c_str(), "DEFAULT") == 0) {
            Relaxation_Factor[PRIMARY_FILTER] = DEFAULT;
        } else {
            std::stringstream valuestream;
            valuestream << value_string;
            valuestream >> Relaxation_Factor[PRIMARY_FILTER];
        }
    } 
    else if (strcmp(code, "ExplicitFilter[2].Haselbacher.Relaxation_Factor") == 0) {
        i_command = 131;
        value >> value_string;
        if (strcmp(value_string.c_str(), "DEFAULT") == 0) {
            Relaxation_Factor[SECONDARY_FILTER] = DEFAULT;
        } else {
            std::stringstream valuestream;
            valuestream << value_string;
            valuestream >> Relaxation_Factor[SECONDARY_FILTER];
        }
    } 
    
    else {
        i_command = INVALID_INPUT_CODE;
    } /* endif */
    
    return i_command;
    
}

/******************************************************************************
 * Explicit_Filters_Input_Parameters::Check_Inputs -- Check input values.     *
 ******************************************************************************/
int Explicit_Filters_Input_Parameters::Check_Inputs(void) {
    
    // Input parameters are consistent.  Exit successfully.
    return 0;
    
}

/***************************************************************************
 * Explicit_Filters_Input_Parameters -- Input-output operators.            *
 ***************************************************************************/
ostream &operator << (ostream &out_file,
                      const Explicit_Filters_Input_Parameters &IP) {
    
    IP.Output(out_file);
    return (out_file);
    
}

istream &operator >> (istream &in_file,
                      Explicit_Filters_Input_Parameters &IP) {
    
    return in_file;
    
}

void Explicit_Filters_Input_Parameters::Output(ostream &out_file) const {
    

    out_file << "\n    -> Filter type: ";
    if (Filter_Type[PRIMARY_FILTER] == IMPLICIT_FILTER){
        out_file << "\n       -> Filter Grid Ratio: " << FGR[PRIMARY_FILTER];
    }
    
}


/************************************************************************************
 * Explicit_Filters_Input_Parameters::Broadcast -- Broadcast to all processors.     *
 ************************************************************************************/
void Explicit_Filters_Input_Parameters::Broadcast(void) {
    
#ifdef _MPI_VERSION
    
    MPI::COMM_WORLD.Bcast(&(Filter_Initial_Condition), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Solution_Filtering_Frequency), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Filter_Solution_Before_Execution), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Commutation_Order), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Filter_Method), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Finite_Differencing_Order), 1, MPI::INT, 0);

    for (int filter_number = PRIMARY_FILTER, filter_number <= SECONDARY_FILTER, filter_number++) {
        MPI::COMM_WORLD.Bcast(&(Filter_Type[filter_number]), 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&(FGR[filter_number]), 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&(Filter_Width[filter_number]), 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&(Number_Of_Rings[filter_number]), 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&(Filter_Memory_Efficient[filter_number]), 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&(Filter_Relative[filter_number]), 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&(Filter_Strength[filter_number]), 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&(Relaxation_Factor[filter_number]), 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&(Least_Squares_Filter_Weighting_Factor[filter_number]), 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&(Least_Squares_Filter_Weighting[filter_number]), 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&(Target_Filter_Sharpness[filter_number]), 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&(Filter_Width_Strict[filter_number]), 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&(LS_Constraints[filter_number]), 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&(Derivative_Constraints[filter_number]), 1, MPI::INT, 0);        
    }
    

#endif
    
}