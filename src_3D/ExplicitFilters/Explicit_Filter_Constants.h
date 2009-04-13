/*
 *  Explicit_Filter_Constants.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 11/04/09.
 *  Copyright 2009. All rights reserved.
 *
 */

#ifndef _EXPLICIT_FILTER_CONSTANTS_INCLUDED
#define _EXPLICIT_FILTER_CONSTANTS_INCLUDED

namespace Explicit_Filter_Constants {
    
    enum Filter_Method         {FILTER_VARIABLES, FILTER_RESIDUALS};
    enum Filter_Cell           {CORNER_CELL, FACE_CELL, EDGE_CELL, INNER_CELL, MIDDLE_CELL};
    enum Filter_Type           {IMPLICIT_FILTER, HASELBACHER_FILTER, VASILYEV_FILTER, RESTART_FILTER};
    enum Filter_Number         {PRIMARY_FILTER,SECONDARY_FILTER};
}

#endif
