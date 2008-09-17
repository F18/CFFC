/*
 *  Explicit_Filter.cc
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/03/08. 
 */
#include "Explicit_Filter_Helpers.h"
bool   Explicit_Filter_Properties::debug_flag = true;
bool   Explicit_Filter_Properties::Memory_Efficient = false;
int    Explicit_Filter_Properties::commutation_order = 2;
int    Explicit_Filter_Properties::finite_differencing_order = 4;
int    Explicit_Filter_Properties::number_of_rings = 2;
bool   Explicit_Filter_Properties::LS_constraints = ON;
int    Explicit_Filter_Properties::Derivative_constraints = DEFAULT;
bool   Explicit_Filter_Properties::Filter_Width_strict = OFF;
double Explicit_Filter_Properties::FGR = 2.0;
double Explicit_Filter_Properties::target_filter_sharpness = -1.0;
bool   Explicit_Filter_Properties::restarted = false;
int    Explicit_Filter_Properties::progress_mode = PROGRESS_MODE_MESSAGE;
int    Explicit_Filter_Properties::filter_type = FILTER_TYPE_HASELBACHER;
char   *(Explicit_Filter_Properties::output_file_name) = NULL;
int    Explicit_Filter_Properties::batch_flag = 0;
int    Explicit_Filter_Properties::number_of_rings_increased = 3;
int    Explicit_Filter_Properties::derivative_accuracy = 3;
double Explicit_Filter_Properties::G_cutoff = exp(-sqr(PI)/(4.*6.));
