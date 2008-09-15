/*
 *  Gaussian_Filter.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 04/05/08.
 *
 */


#ifndef _GAUSSIAN_FILTER_INCLUDED
#define _GAUSSIAN_FILTER_INCLUDED

#include "Analytic_Filter.h"

template <typename Soln_pState, typename Soln_cState>
class Gaussian_Filter : public Analytic_Filter<Soln_pState,Soln_cState> {
public:
    using Analytic_Filter<Soln_pState,Soln_cState>::theNeighbours;
    using Analytic_Filter<Soln_pState,Soln_cState>::FGR;
    
    Gaussian_Filter(void) : Analytic_Filter<Soln_pState,Soln_cState>() { 
    }
    
    
    void Get_Neighbours(Cell3D &theCell);
    
    string filter_name(void);
    int filter_type(void);
        
    double Filter_Function(Vector3D &X, Vector3D &Delta) ;
    double Filter_Function(double dx, double dy, double dz, double FGR) ;
    double Filter_Function(Vector3D &X, double &Delta);

};



template <typename Soln_pState, typename Soln_cState>
void Gaussian_Filter<Soln_pState,Soln_cState>::Get_Neighbours(Cell3D &theCell) {
    int number_of_rings = int(ceil(FGR));
    theNeighbours.GetNeighbours(theCell,int(ceil(FGR))+1,FILTER_TYPE_VASILYEV); //temporary choose same method as for Vasilyev
}

template <typename Soln_pState, typename Soln_cState>
string Gaussian_Filter<Soln_pState,Soln_cState>::filter_name(void) {
    return "Gaussian";
}

template <typename Soln_pState, typename Soln_cState>
int Gaussian_Filter<Soln_pState,Soln_cState>::filter_type(void) {
    return FILTER_TYPE_GAUSSIAN;
}
    
template <typename Soln_pState, typename Soln_cState>
double Gaussian_Filter<Soln_pState,Soln_cState>::Filter_Function(Vector3D &X, Vector3D &Delta) {
    
    Vector3D absX(fabs(X.x),fabs(X.y),fabs(X.z));
    if(absX<=Delta/2)
		return (1/Delta.abs());
	else
		return (0.0);
}

template <typename Soln_pState, typename Soln_cState>
double Gaussian_Filter<Soln_pState,Soln_cState>::Filter_Function(Vector3D &X, double &Delta) {
    
    //Vector3D Delta(theNeighbours.Delta);
    double dr;
    
    //Delta *= weight_factor;
    
        dr=X.sqr();
        return  sqrt(SIX/(PI*sqr(Delta))*exp(- SIX * dr/sqr(Delta)));
    
}

/* This one is for computational uniform grid : dx = dy = dz = 1 */
template <typename Soln_pState, typename Soln_cState>
double Gaussian_Filter<Soln_pState,Soln_cState>::Filter_Function(double dx, double dy, double dz, double FGR)  {
    
    double dr = sqrt( sqr(dx) + sqr(dy) + sqr(dz) );
    if(dr<=FGR/2.0)
		return (1.0/FGR);
	else
		return (0.0);
}
#endif
