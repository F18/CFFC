/*
 *  General_Filter.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/04/08.
 *
 */


#ifndef _GENERAL_FILTER_INCLUDED
#define _GENERAL_FILTER_INCLUDED

#include "Neighbours.h"
#include "LES_Filters.h"

/**
 * CLASS: General Filter
 * Abstract class
 */
template <typename Soln_pState, typename Soln_cState>
class General_Filter {
public:
    virtual RowVector filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) = 0;
    virtual void transfer_function(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) = 0;
    virtual void filter_tests(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) = 0;
    virtual void Reset_Filter_Weights(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk) = 0;
};

#endif