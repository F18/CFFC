/*
 *  Tophat_Filter.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 04/05/08.
 *
 */


#ifndef _TOPHAT_FILTER_INCLUDED
#define _TOPHAT_FILTER_INCLUDED

#include "Analytic_Filter.h"

template <typename Soln_pState, typename Soln_cState>
class Tophat_Filter : public Analytic_Filter<Soln_pState,Soln_cState> {
public:
    using Analytic_Filter<Soln_pState,Soln_cState>::theNeighbours;
    using Analytic_Filter<Soln_pState,Soln_cState>::FGR;
    
    Tophat_Filter(void) : Analytic_Filter<Soln_pState,Soln_cState>() { 
    }
    
    
    void Get_Neighbours(Cell3D &theCell);
    
    string filter_name(void);
    int filter_type(void);
        
    double Filter_Function(Vector3D &X, Vector3D &Delta) ;
    double Filter_Function(double dx, double dy, double dz, double FGR) ;
    double Filter_Function(Vector3D &X, double &Delta);

};



template <typename Soln_pState, typename Soln_cState>
void Tophat_Filter<Soln_pState,Soln_cState>::Get_Neighbours(Cell3D &theCell) {
    int number_of_rings = int(ceil(FGR));
    theNeighbours.GetNeighbours_Vasilyev(theCell,int(ceil(FGR/2.0)));
}

template <typename Soln_pState, typename Soln_cState>
string Tophat_Filter<Soln_pState,Soln_cState>::filter_name(void) {
    return "Tophat";
}

template <typename Soln_pState, typename Soln_cState>
int Tophat_Filter<Soln_pState,Soln_cState>::filter_type(void) {
    return FILTER_TYPE_TOPHAT;
}
    
template <typename Soln_pState, typename Soln_cState>
double Tophat_Filter<Soln_pState,Soln_cState>::Filter_Function(Vector3D &X, Vector3D &Delta) {
    
    Vector3D absX(fabs(X.x),fabs(X.y),fabs(X.z));
    if(absX<=Delta/2)
		return (1/Delta.abs());
	else
		return (0.0);
}

template <typename Soln_pState, typename Soln_cState>
double Tophat_Filter<Soln_pState,Soln_cState>::Filter_Function(Vector3D &X, double &Delta) {
    
    Vector3D absX(fabs(X.x),fabs(X.y),fabs(X.z));
    if(absX<=Vector3D(Delta/2.0,Delta/2.0,Delta/2.0))
		return (1.0/Delta);
	else
		return (0.0);
}

/* This one is for computational uniform grid : dx = dy = dz = 1 */
template <typename Soln_pState, typename Soln_cState>
double Tophat_Filter<Soln_pState,Soln_cState>::Filter_Function(double dx, double dy, double dz, double FGR)  {
    
    double dr = sqrt( sqr(dx) + sqr(dy) + sqr(dz) );
    if(dr<=FGR/2.0)
		return (1.0/FGR);
	else
		return (0.0);
}
#endif
