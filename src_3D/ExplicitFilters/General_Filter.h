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
//#include "Explicit_Filter.h"

template <typename Soln_pState, typename Soln_cState> class Explicit_Filters;
/**
 * CLASS: General Filter
 * Abstract class
 */
template <typename Soln_pState, typename Soln_cState>
class General_Filter {
public:
    virtual RowVector filter(Explicit_Filters<Soln_pState,Soln_cState> &explicit_filter, Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) =0;
    virtual RowVector filter_1D(Explicit_Filters<Soln_pState,Soln_cState> &explicit_filter, Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, int direction) = 0;
    virtual void transfer_function(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) = 0;
    virtual void filter_tests(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) = 0;
    virtual int filter_type(void) = 0;
    virtual string filter_name(void) = 0;
//    virtual void Write_to_file(Grid3D_Hexa_Block &Grid_Blk, ofstream &out_file) = 0;
//    virtual void Read_from_file(Grid3D_Hexa_Block &Grid_Blk, ifstream &in_file) = 0;

};

#endif
