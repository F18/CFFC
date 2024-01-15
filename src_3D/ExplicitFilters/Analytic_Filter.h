/*
 *  Analytic_Filter.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 04/05/08.
 *
 */


#ifndef _ANALYTIC_FILTER_INCLUDED
#define _ANALYTIC_FILTER_INCLUDED

#include "General_Filter.h"
#include "Explicit_Filter_Helpers.h"
#include "Neighbours.h"

#include <complex>

#include "../Math/Math.h"
#include "../Math/LinearSystems.h"
#include "../Math/Matrix.h"

#include "../Utilities/Utilities.h"

#ifndef _GNUPLOT_INCLUDED
#include "../system/gnuplot.h"
#endif

#include <cstdlib>

typedef complex<double> Complex;

/**
 * CLASS: Analytic_Filter
 */

template <typename Soln_pState, typename Soln_cState>
class Analytic_Filter : public General_Filter<Soln_pState,Soln_cState> {
public:
    
    Analytic_Filter(void) {
        FGR = Explicit_Filter_Properties::FGR;
        number_of_rings = Explicit_Filter_Properties::number_of_rings;
        theNeighbours.allocate(number_of_rings);
        I = Complex(0,1);
    }
    
    Complex I;
    Neighbours theNeighbours;
    Cell3D theCell;
    double FGR;
    int number_of_rings;

    RowVector filter(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell);
    RowVector filter_1D(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, int direction);

    virtual double Filter_Function(Vector3D &X, Vector3D &Delta) = 0;
    virtual double Filter_Function(double dx, double dy, double dz, double FGR) = 0;
    virtual double Filter_Function(Vector3D &X, double &Delta) = 0;

    virtual void Get_Neighbours(Cell3D &theCell) = 0;

    /* --------------- Transfer function calculations ------------- */
    void transfer_function(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell);
    virtual string filter_name(void) = 0;
    virtual int filter_type(void) = 0;
    
    void Write_to_file(Grid3D_Hexa_Block &Grid_Blk, ofstream &out_file);
    void Read_from_file(Grid3D_Hexa_Block &Grid_Blk, ifstream &in_file);
    
    void filter_tests(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell);
    void Reset_Filter_Weights(Grid3D_Hexa_Block &Grid_Blk);
};

template <typename Soln_pState, typename Soln_cState>
RowVector Analytic_Filter<Soln_pState,Soln_cState>::filter(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) {
    
    theNeighbours.set_grid(Grid_Blk);
    Get_Neighbours(theCell);

    double Delta = Explicit_Filter_Properties::FGR*theCell.dXc.abs();
    
    RowVector phi_bar = Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::FilterVariable(theCell);
    phi_bar.zero();
    Vector3D dX;

    
    /* quick and dirty piecewise constant convolution -- only good for tophat filter */
    for (int i=0; i<theNeighbours.number_of_neighbours; i++) {
        dX = theNeighbours.neighbour[i].Xc-theCell.Xc;
        phi_bar += theNeighbours.neighbour[i].V * Filter_Function(dX,Delta) * Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::FilterVariable(theNeighbours.neighbour[i]);
    }
    
    return phi_bar;
}


template <typename Soln_pState, typename Soln_cState>
RowVector Analytic_Filter<Soln_pState,Soln_cState>::filter_1D(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, int direction) {
    
    // Not implemented yet
    
    return RowVector(0);
}

template <typename Soln_pState, typename Soln_cState>
void Analytic_Filter<Soln_pState,Soln_cState>::transfer_function(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) {
    // Not implemented yet

}


template <typename Soln_pState, typename Soln_cState>
void Analytic_Filter<Soln_pState,Soln_cState>::Write_to_file(Grid3D_Hexa_Block &Grid_Blk, ofstream &out_file) {
    // Not implemented yet

}

template <typename Soln_pState, typename Soln_cState>
void Analytic_Filter<Soln_pState,Soln_cState>::Read_from_file(Grid3D_Hexa_Block &Grid_Blk, ifstream &in_file) {
    // Not implemented yet

}

template <typename Soln_pState, typename Soln_cState>
void Analytic_Filter<Soln_pState,Soln_cState>::filter_tests(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) {
    // Not implemented yet

}

template <typename Soln_pState, typename Soln_cState>
void Analytic_Filter<Soln_pState,Soln_cState>::Reset_Filter_Weights(Grid3D_Hexa_Block &Grid_Blk) {
    // Not implemented yet

}


#endif
