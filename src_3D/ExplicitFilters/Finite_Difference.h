/*
 *  NonLinearSystems.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 27/04/08.
 *
 */

#ifndef _FINITE_DIFFERENCE_FOR_EXPLICIT_FILTERS_INCLUDED
#define _FINITE_DIFFERENCE_FOR_EXPLICIT_FILTERS_INCLUDED

#include "../Math/Math.h"
#include "../Math/LinearSystems.h"
#include "../Math/Matrix.h"
#include "../Utilities/Utilities.h"
#include "Explicit_Filter_Helpers.h"

#define DFDI 1
#define DFDJ 2
#define DFDK 3

#define DFDX 4
#define DFDY 5
#define DFDZ 6
#define DFDR 7

template <typename Soln_pState, typename Soln_cState>
class Finite_Difference_Class {
private:
    typedef Explicit_Filter_Adaptor<Soln_pState,Soln_cState> adaptor;
    int order;
    RowVector Central_Finite_Difference(const int i, const int j, const int k, const int derivative, const double &dt, int order);
    RowVector Forward_Finite_Difference(const int i, const int j, const int k, const int derivative, const double &dt, int order);
    RowVector Backward_Finite_Difference(const int i, const int j, const int k, const int derivative, const double &dt, int order);
    RowVector Equally_Spaced_Finite_Difference(Grid3D_Hexa_Block &Grid, const int i, const int j, const int k, const int derivative);
public:
    Finite_Difference_Class(int order_of_accuracy) {
        order = order_of_accuracy;
        cout << "\n\n finite difference order = " << order << "    central_rings = " << Get_central_rings() << endl;
    }
    RowVector Finite_Difference(Grid3D_Hexa_Block &Grid, Cell3D &theCell, const int derivative);
    int Get_central_rings(void) { return int(ceil(order/2)); }

};

/* Finite difference for equally spaced samples */
template <typename Soln_pState, typename Soln_cState>
RowVector Finite_Difference_Class<Soln_pState,Soln_cState>::Central_Finite_Difference(const int i, const int j, const int k, const int derivative, const double &dt, int order) {
    int n=int(ceil(order/2));
    int N=2*n;
    
    RowVector coefficients(N);
    RowVector samplerow = adaptor::FilterVariable(i,j,k);
    DenseMatrix samples(N,samplerow.size());
    
    for(int p=1; p<=n; p++) {
        switch (derivative) {
            case DFDI:
                adaptor::FillMatrixRow(samples, n-p,    i-p,j,k);
                adaptor::FillMatrixRow(samples, n+p-1,  i+p,j,k);
                break;
            case DFDJ:
                adaptor::FillMatrixRow(samples, n-p,    i,j-p,k);
                adaptor::FillMatrixRow(samples, n+p-1,  i,j+p,k);
                break;
            case DFDK:
                adaptor::FillMatrixRow(samples, n-p,    i,j,k-p);
                adaptor::FillMatrixRow(samples, n+p-1,  i,j,k+p);
                break;
        }
    }
    
    
    switch (N) {
        case 2:
            /* 2nd order */
            coefficients(0) = -1.0;
            coefficients(1) = 1.0;
            coefficients /= 2.0;
            break;
        case 4:
            /* 4th order */
            coefficients(0) =  1.0;
            coefficients(1) = -8.0;
            coefficients(2) =  8.0;
            coefficients(3) = -1.0;
            coefficients /= 12.0;
            break;
        case 6:
            /* 6th order */
            coefficients(0) = -1.0;
            coefficients(1) =  9.0;
            coefficients(2) = -45.0;
            coefficients(3) =  45.0;
            coefficients(4) = -9.0;
            coefficients(5) =  1.0;
            coefficients /= 60.0;
            break;
        case 8:
            /* 8th order */
            coefficients(0) =  3.0;
            coefficients(1) = -32.0;
            coefficients(2) =  168.0;
            coefficients(3) = -672.0;
            coefficients(4) =  672.0;
            coefficients(5) = -168.0;
            coefficients(6) =  32.0;
            coefficients(7) = -3.0;
            coefficients /= 840.0;
            break;
        case 10:
            /* 10th order */
            coefficients(0) = -2.0;
            coefficients(1) =  25.0;
            coefficients(2) = -150.0;
            coefficients(3) =  600.0;
            coefficients(4) = -2100.0;
            coefficients(5) =  2100.0;
            coefficients(6) = -600.0;
            coefficients(7) =  150.0;
            coefficients(8) = -25.0;
            coefficients(9) =  2.0;
            coefficients /= 2520.0;
            break;
    }
    
    return (coefficients*samples)/dt;
}

/* Finite difference for equally spaced samples */
template <typename Soln_pState, typename Soln_cState>
RowVector Finite_Difference_Class<Soln_pState,Soln_cState>::Forward_Finite_Difference(const int i, const int j, const int k, const int derivative, const double &dt, int order) {
    
    
    int n = int(ceil(order/2));
    order = 2*n;
    int N=2*n+1;
    
    RowVector coefficients(N);
    RowVector samplerow = adaptor::FilterVariable(i,j,k);
    DenseMatrix samples(N,samplerow.size());
    
    for(int p=1; p<=N; p++) {
        switch (derivative) {
            case DFDI:
                adaptor::FillMatrixRow(samples, N-p,    i-p,j,k);
                break;
            case DFDJ:
                adaptor::FillMatrixRow(samples, N-p,    i,j-p,k);
                break;
            case DFDK:
                adaptor::FillMatrixRow(samples, N-p,    i,j,k-p);
                break;
        }
    }
    
    switch (order) {
        case 2:
            /* 2nd order */
            coefficients(0) =  1.0;
            coefficients(1) = -4.0;
            coefficients(2) =  3.0;
            coefficients /= 2.0;
            break;
        case 4:
            /* 4th order */
            coefficients(0) =  3.0;
            coefficients(1) = -16.0;
            coefficients(2) =  36.0;
            coefficients(3) = -48.0;
            coefficients(4) =  25.0;
            coefficients /= 12.0;
            break;
        case 6:
            /* 6th order */
            coefficients(0) =  10.0;
            coefficients(1) = -72.0;
            coefficients(2) =  225.0;
            coefficients(3) = -400.0;
            coefficients(4) =  450.0;
            coefficients(5) = -360.0;
            coefficients(6) =  147.0;
            coefficients /= 60.0;
            break;
        case 8:
            /* 8th order */
            coefficients(0) =  105.0;
            coefficients(1) = -960.0;
            coefficients(2) =  3920.0;
            coefficients(3) = -9408.0;
            coefficients(4) =  14700.0;
            coefficients(5) = -15680.0;
            coefficients(6) =  11760.0;
            coefficients(7) = -6720.0;
            coefficients(8) =  2283.0;
            coefficients /= 840.0;
            break;
        case 10:
            /* 10th order */
            coefficients(0)  =  252.0;
            coefficients(1)  = -2800.0;
            coefficients(2)  =  14175.0;
            coefficients(3)  = -43200.0;
            coefficients(4)  =  88200.0;
            coefficients(5)  = -127008.0;
            coefficients(6)  =  132300.0;
            coefficients(7)  = -100800.0;
            coefficients(8)  =  56700.0;
            coefficients(9)  = -25200.0;
            coefficients(10) =  7381.0;
            coefficients /= 2520.0;
            break;
    }
    
    return (coefficients*samples)/dt;
}


/* Finite difference for equally spaced samples */
template <typename Soln_pState, typename Soln_cState>
RowVector Finite_Difference_Class<Soln_pState,Soln_cState>::Backward_Finite_Difference(const int i, const int j, const int k, const int derivative, const double &dt, int order) {
    
    int n=int(ceil(order/2));
    order=2*n;
    int N=2*n+1;
    
    RowVector coefficients(N);
    RowVector samplerow = adaptor::FilterVariable(i,j,k);
    DenseMatrix samples(N,samplerow.size());
    
    for(int p=0; p<N; p++) {
        switch (derivative) {
            case DFDI:
                adaptor::FillMatrixRow(samples, p,    i+p,j,k);
                break;
            case DFDJ:
                adaptor::FillMatrixRow(samples, p,    i,j+p,k);
                break;
            case DFDK:
                adaptor::FillMatrixRow(samples, p,    i,j,k+p);
                break;
        }
    }
    
    switch (order) {
        case 2:
            /* 2nd order */
            coefficients(0) = -3.0;
            coefficients(1) =  4.0;
            coefficients(2) = -1.0;
            coefficients /= 2.0;
            break;
        case 4:
            /* 4th order */
            coefficients(0) = -25.0;
            coefficients(1) =  48.0;
            coefficients(2) = -36.0;
            coefficients(3) =  16.0;
            coefficients(4) = -3.0;
            coefficients /= 12.0;
            break;
        case 6:
            /* 6th order */
            coefficients(0) = -147.0;
            coefficients(1) =  360.0;
            coefficients(2) = -450.0;
            coefficients(3) =  400.0;
            coefficients(4) = -225.0;
            coefficients(5) =  72.0;
            coefficients(6) = -10.0;
            coefficients /= 60.0;
            break;
        case 8:
            /* 8th order */
            coefficients(0) =  2283.0;
            coefficients(1) = -6720.0;
            coefficients(2) =  11760.0;
            coefficients(3) = -15680.0;
            coefficients(4) =  14700.0;
            coefficients(5) = -9408.0;
            coefficients(6) =  3920.0;
            coefficients(7) = -960.0;
            coefficients(8) =  105.0;
            coefficients /= 840.0;
            break;
        case 10:
            /* 10th order */
            coefficients(0)  = -7381.0;
            coefficients(1)  =  25200.0;
            coefficients(2)  = -56700.0;
            coefficients(3)  =  100800.0;
            coefficients(4)  = -132300.0;
            coefficients(5)  =  127008.0;
            coefficients(6)  = -88200.0;
            coefficients(7)  =  43200.0;
            coefficients(8)  = -14175.0;
            coefficients(9)  =  2800.0;
            coefficients(10) = -252.0;
            coefficients /= 2520.0;
            break;
    }
    
    return (coefficients*samples)/dt;
}

template <typename Soln_pState, typename Soln_cState>
RowVector Finite_Difference_Class<Soln_pState,Soln_cState>::Equally_Spaced_Finite_Difference(Grid3D_Hexa_Block &Grid, const int i, const int j, const int k, const int derivative) {
    
    double dt = 1.0;
    
    int n = Get_central_rings();
    
    int index, first_index, last_index;
    switch(derivative) {
        case DFDI:
            index = i;
            first_index = Grid.Nghost;
            last_index = Grid.NCi-Grid.Nghost-1;
            break;
        case DFDJ:
            index = j;
            first_index = Grid.Nghost;
            last_index = Grid.NCj-Grid.Nghost-1;
            break;
        case DFDK:
            index = k;
            first_index = Grid.Nghost;
            last_index = Grid.NCk-Grid.Nghost-1;
            break;
    }
        
    if (n > index - first_index) {
        if (last_index - index < n) {
            cout << "backward finite differencing: " << endl;
            cout << "  index = " << index << endl;
            cout << "  first_index = " << first_index << endl;
            cout << "  last_index = " << last_index << endl;
            cout << "  n = " << n << endl;
            cout << "  order = " << order << endl;
            exit(1);
        }
        return Backward_Finite_Difference(i,j,k, derivative, dt, order);
    } else if (n > last_index - index) {
        if (index - first_index < n) {
            cout << "forward finite differencing: " << endl;
            cout << "  index = " << index << endl;
            cout << "  first_index = " << first_index << endl;
            cout << "  last_index = " << last_index << endl;
            cout << "  n = " << n << endl;
            cout << "  order = " << order << endl;
            exit(1);
        }
        return  Forward_Finite_Difference(i,j,k, derivative, dt, order);
    } else {
        if (last_index - index < n && index - first_index < n) {
            cout << "central finite differencing: " << endl;
            cout << "  index = " << index << endl;
            cout << "  first_index = " << first_index << endl;
            cout << "  last_index = " << last_index << endl;
            cout << "  n = " << n << endl;
            cout << "  order = " << order << endl;
            exit(1);
        }
        //assert(last_index - index < n && index - first_index > n);
        return  Central_Finite_Difference(i,j,k, derivative, dt, order);
    }
}

template <typename Soln_pState, typename Soln_cState>
RowVector Finite_Difference_Class<Soln_pState,Soln_cState>::Finite_Difference(Grid3D_Hexa_Block &Grid, Cell3D &theCell, const int derivative) {
    
    int i = theCell.I;
    int j = theCell.J;
    int k = theCell.K;

    // Calculate the Jacobian matrix
    DenseMatrix J(3,3);
    Grid.Jacobian_Matrix(J,i,j,k,order);
        
    
    RowVector dfdi = Equally_Spaced_Finite_Difference(Grid,i,j,k,DFDI);
    RowVector dfdj = Equally_Spaced_Finite_Difference(Grid,i,j,k,DFDJ);
    RowVector dfdk = Equally_Spaced_Finite_Difference(Grid,i,j,k,DFDK);
    
    switch(derivative) {
        case DFDX:
            return J(0,0) * dfdi   +   J(0,1) * dfdj   +   J(0,2) * dfdk;
        case DFDY:
            return J(1,0) * dfdi   +   J(1,1) * dfdj   +   J(1,2) * dfdk;
        case DFDZ:
            return J(2,0) * dfdi   +   J(2,1) * dfdj   +   J(2,2) * dfdk;
        case DFDR: {
            double x = theCell.Xc.x;
            double y = theCell.Xc.y;
            double z = theCell.Xc.z;
            double r = sqrt(sqr(x)+sqr(y)+sqr(z));
            double theta = atan2(y,x);
            double phi = acos(z/r);
            double dxdr = cos(theta)*sin(phi);
            double dydr = sin(theta)*sin(phi);
            double dzdr = cos(phi);
            
            return dxdr * (J(0,0) * dfdi   +   J(0,1) * dfdj   +   J(0,2) * dfdk) +
                   dydr * (J(1,0) * dfdi   +   J(1,1) * dfdj   +   J(1,2) * dfdk) +
                   dzdr * (J(2,0) * dfdi   +   J(2,1) * dfdj   +   J(2,2) * dfdk);
                   
            
        }
    }
}

#endif
