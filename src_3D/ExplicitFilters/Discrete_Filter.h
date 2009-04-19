/*
 *  Discrete_Filter.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 04/05/08.
 *
 */


#ifndef _DISCRETE_FILTER_INCLUDED
#define _DISCRETE_FILTER_INCLUDED

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


template <typename Soln_pState, typename Soln_cState> class Explicit_Filters;

/**
 * CLASS: Discrete_Filter
 */

template <typename Soln_pState, typename Soln_cState>
class Discrete_Filter : public General_Filter<Soln_pState,Soln_cState> {
public:
    
    Discrete_Filter(Explicit_Filter_Properties &explicit_filter_properties) {
        properties = &explicit_filter_properties;
        I = Complex(0,1);
    }
    
    virtual void Read_Basic_Properties(void) {
        properties->Get_Property(debug_flag,"debug_flag");
        properties->Get_Property(batch_flag,"batch_flag");
        properties->Get_Property(FGR,"FGR");
        properties->Get_Property(use_fixed_filter_width,"use_fixed_filter_width");
        properties->Get_Property(fixed_filter_width,"fixed_filter_width");
        properties->Get_Property(commutation_order,"commutation_order");
        properties->Get_Property(number_of_rings,"number_of_rings");
        properties->Get_Property(G_cutoff,"G_cutoff");
        Store_Filter_Weights = !properties->Get_Property_int("memory_efficient");
    }
    
    virtual void Read_Properties(void) = 0;
    
    Explicit_Filter_Properties *properties;
    Complex I;
    Neighbours theNeighbours;
    Cell3D theCell;
    int debug_flag;
    int batch_flag;
    int commutation_order;
    double FGR;
    int number_of_rings;
    int Store_Filter_Weights;
    double G_cutoff;
    int use_fixed_filter_width;
    double fixed_filter_width;

    DenseMatrix Neighbouring_Values;
    void Set_Neighbouring_Values(DenseMatrix &Neighbouring_Values, Neighbours &theNeighbours);

    RowVector filter(Explicit_Filters<Soln_pState,Soln_cState> &explicit_filter, Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell);
    RowVector filter_1D(Explicit_Filters<Soln_pState,Soln_cState> &explicit_filter,Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, int direction);

    
    
    /* --------------- Transfer function calculations ------------- */
    void transfer_function(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell);
    Complex G_function(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &k, RowVector &w);
    Complex G_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    Complex G_function_1D_100(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    Complex G_function_1D_110(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    Complex dG_function(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &k, RowVector &w);
    Complex dG_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    
    int Output_Transfer_Function_tecplot(string prefix_string, Cell3D  ***kCells, Complex ***G, int N);
    int Output_Transfer_Function_gnuplot(string prefix_string, double *k_111, double *G_111, string legend_111, double *k_110, double *G_110, string legend_110, double *K_100, double *G_100, string legend_100, int N, string &title);

    double Filter_Grid_Ratio(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax);
    double Filter_Grid_Ratio_111(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax);
    double Filter_Grid_Ratio_110(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax);
    double Filter_Grid_Ratio_100(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax);
    double Filter_Grid_Ratio_010(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax);
    double Filter_Grid_Ratio_001(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax);

    double Filter_Width(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kdir, RowVector &w, double G_cutoff);

    Vector3D Calculate_wavenumber_of_Gvalue(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w, double G_value);
    Vector3D Calculate_wavenumber_of_dGvalue(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmin, Vector3D &kmax, RowVector &w, double dG_value);

    virtual double filter_moment(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, int q, int r, int s) = 0;
    virtual double filter_moment_1D(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, int q, int direction) = 0;
    void check_filter_moments(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell);
    void check_filter_moments_1D(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, int direction);

    
    virtual void Get_Neighbours(Cell3D &theCell) = 0;
    virtual void Get_Neighbours_1D(Cell3D &theCell, int direction) = 0;
    virtual RowVector Get_Weights(Cell3D &theCell, Neighbours &theNeighbours) = 0;
    virtual RowVector Get_Weights_1D(Cell3D &theCell, Neighbours &theNeighbours, int direction) = 0;

    virtual string filter_name(void) = 0;
    virtual int filter_type(void) = 0;    
};

template <typename Soln_pState, typename Soln_cState>
inline Complex Discrete_Filter<Soln_pState,Soln_cState>::G_function(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &k, RowVector &w) {
    Complex G(0,0);

    Vector3D X0(theCell.Xc);
    Vector3D X;
    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
        X = theNeighbours.neighbour[n]->Xc;
        G += w(n) * exp(-I*(k*(X-X0)));
    }
    return G;
}

template <typename Soln_pState, typename Soln_cState>
inline Complex Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    Vector3D k(wave_number/sqrt(THREE),wave_number/sqrt(THREE),wave_number/sqrt(THREE));    
    return (G_function(theCell, theNeighbours, k, w));
}


template <typename Soln_pState, typename Soln_cState>
inline Complex Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D_100(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    Vector3D k(wave_number,ZERO,ZERO);    
    return (G_function(theCell, theNeighbours, k, w));
}

template <typename Soln_pState, typename Soln_cState>
inline Complex Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D_110(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    Vector3D k(wave_number/sqrt(TWO),wave_number/sqrt(TWO),ZERO);    
    return (G_function(theCell, theNeighbours, k, w));
}



template <typename Soln_pState, typename Soln_cState>
void Discrete_Filter<Soln_pState,Soln_cState>::check_filter_moments(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) {
    
    if (!batch_flag) {
        cout << "\nCalculating filter moments at cell ("<<theCell.I<<","<<theCell.J<<","<<theCell.K<<")";
        cout.flush();
    }
    
    theNeighbours.set_grid(Grid_Blk);
    
    Get_Neighbours(theCell);
    
    RowVector w = Get_Weights(theCell, theNeighbours);
    
    /*
    double ***M = new double **[commutation_order+1];
    for (int q=0; q<=commutation_order; q++) {
        M[q] = new double *[commutation_order+1];
        for (int r=0; r<=commutation_order; r++) {
            M[q][r] = new double [commutation_order+1];
        }
    }
    */
    
    int real_commutation_order(1);
    int error_in_M0=0;
    int previous_order=0;
    for(int order=0; order<=commutation_order+2; order++){
        for( int q=0; q<=commutation_order+2; q++) {
            for (int r=0; r<=commutation_order+2; r++) {
                for (int s=0; s<=commutation_order+2; s++) {
                    if (q + r + s == order) {
                        double M = filter_moment(theCell, theNeighbours, w, q,r,s);
                        if (order==0) {
                            cout << "\n M(" << q << "," << r << "," << s << ") = " << M;
                            if (!(M-PICO<ONE && ONE<M+PICO)) {
                                error_in_M0=1;
                            }
                        }
                        if (order!=0 && order>previous_order && fabs(M)<PICO) {
                            real_commutation_order = max(real_commutation_order,order+1);
                        }
                        if (fabs(M)<PICO) {
                            cout << "\n M(" << q << "," << r << "," << s << ") = " << 0;
                        } else {
                            cout << "\n M(" << q << "," << r << "," << s << ") = " << M;
                        }
                        
                        previous_order=order;
                    }
                }
            }
        }
    }
    cout << endl;
    if (error_in_M0) {
        cout << " transfer function  G(0)!=1   :(" << endl;
    }
    cout << " Commutation order = " << real_commutation_order << endl;
    
    
}

template <typename Soln_pState, typename Soln_cState>
void Discrete_Filter<Soln_pState,Soln_cState>::check_filter_moments_1D(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, int direction) {
    
    if (!batch_flag) {
        cout << "\n\n Calculating filter moments in 1D at cell ("<<theCell.I<<","<<theCell.J<<","<<theCell.K<<")";
    }
    
    theNeighbours.set_grid(Grid_Blk);
    
    Get_Neighbours_1D(theCell,direction);
    
    RowVector w = Get_Weights_1D(theCell, theNeighbours,direction);

    
    int real_commutation_order(1);
    int error_in_M0=0;
    int previous_order=0;
    for( int q=0; q<=commutation_order+1; q++) {
        double M = filter_moment_1D(theCell, theNeighbours, w, q, direction);
        //cout << "\n M(" << q << ") = " << M;
        if (q==0) {
            if (!(M-NANO<ONE && ONE<M+NANO)) {
                error_in_M0=1;
            }
        }
        if (q!=0 && q>previous_order && fabs(M)<NANO) {
            real_commutation_order = max(real_commutation_order,q+1);
        }
        
        previous_order=q;
    }
    
    cout << endl;
    if (error_in_M0) {
        cout << " M0!=1   :(" << endl;
    }
    cout << " Commutation order = " << real_commutation_order << endl;
    
    
}




template <typename Soln_pState, typename Soln_cState>
RowVector Discrete_Filter<Soln_pState,Soln_cState>::filter(Explicit_Filters<Soln_pState,Soln_cState> &explicit_filter, Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) {
    
    if (!explicit_filter.Filter_Weights_Allocated && Store_Filter_Weights) {
        explicit_filter.Allocate_Filter_Weights();
        explicit_filter.Filter_Weights_Allocated = true;
    }
     
     
     theNeighbours.set_grid(Grid_Blk);
     Get_Neighbours(theCell);
     
     

    if (Store_Filter_Weights) {
        int I(theCell.I);
        int J(theCell.J);
        int K(theCell.K);
        if (!explicit_filter.Filter_Weights_Assigned[I][J][K]) {
            explicit_filter.Filter_Weights[I][J][K] = Get_Weights(theCell,theNeighbours);
            explicit_filter.Filter_Weights_Assigned[I][J][K] = true;
        }
        Set_Neighbouring_Values(Neighbouring_Values,theNeighbours);
        return (explicit_filter.Filter_Weights[I][J][K]*Neighbouring_Values);
    } else {
        RowVector W = Get_Weights(theCell,theNeighbours);
        Set_Neighbouring_Values(Neighbouring_Values,theNeighbours);
        return (W*Neighbouring_Values);
    }
    


    
    
}


template <typename Soln_pState, typename Soln_cState>
RowVector Discrete_Filter<Soln_pState,Soln_cState>::filter_1D(Explicit_Filters<Soln_pState,Soln_cState> &explicit_filter, Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, int direction) {
    
    if (!explicit_filter.Filter_Weights_Allocated && Store_Filter_Weights) {
        explicit_filter.Allocate_Filter_Weights();
        explicit_filter.Filter_Weights_Allocated = true;
    }
    
    
    theNeighbours.set_grid(Grid_Blk);
    Get_Neighbours_1D(theCell, direction);
    
    
    
    if (Store_Filter_Weights) {
        int I(theCell.I);
        int J(theCell.J);
        int K(theCell.K);
        if (!explicit_filter.Filter_Weights_Assigned[I][J][K]) {
            explicit_filter.Filter_Weights[I][J][K] = Get_Weights_1D(theCell,theNeighbours,direction);
            explicit_filter.Filter_Weights_Assigned[I][J][K] = true;
        }
        Set_Neighbouring_Values(Neighbouring_Values,theNeighbours);
        return (explicit_filter.Filter_Weights[I][J][K]*Neighbouring_Values);
    } else {
        RowVector W = Get_Weights_1D(theCell,theNeighbours,direction);
        Set_Neighbouring_Values(Neighbouring_Values,theNeighbours);
        return (W*Neighbouring_Values);
    }

}

template <typename Soln_pState, typename Soln_cState>
void Discrete_Filter<Soln_pState,Soln_cState>::Set_Neighbouring_Values(DenseMatrix &Neighbouring_Values, Neighbours &theNeighbours) {
    Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::FillMatrix(Neighbouring_Values, theNeighbours);
}


template <typename Soln_pState, typename Soln_cState>
void Discrete_Filter<Soln_pState,Soln_cState>::transfer_function(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) {
    
    
    theNeighbours.set_grid(Grid_Blk);
    
    Get_Neighbours(theCell);
    
    RowVector w = Get_Weights(theCell, theNeighbours);
    
    Vector3D kmax;
    Vector3D Delta = theNeighbours.Delta;
    kmax.x = fabs(PI/Delta.x);
    kmax.y = fabs(PI/Delta.y);
    kmax.z = fabs(PI/Delta.z);
    
    if (!batch_flag) {
        cout << "\nCalculating Explicit Filter Transfer function at cell ("<<theCell.I<<","<<theCell.J<<","<<theCell.K<<")";
        cout.flush();
    }
    int N=25;
    std::stringstream prefix;
    prefix << properties->Get_Property_string("output_file_name") << "_";
    prefix << "transfer_function_I" << theCell.I << "J" << theCell.J << "K" << theCell.K;
    
    /* -------------- Output gnuplot ----------------- */
    double *k_111 = new double [N];
    double *k_110 = new double [N];
    double *k_100 = new double [N];
    double *k_010 = new double [N];
    double *k_001 = new double [N];
    double *G_111 = new double [N];
    double *G_110 = new double [N];
    double *G_100 = new double [N];
    double *G_010 = new double [N];
    double *G_001 = new double [N];

    
    Vector3D K_111, K_110, K_100, K_010, K_001;
    K_110.zero(); K_100.zero(); K_010.zero(); K_001.zero();
    for (int i=0; i<N; i++) {
        K_111 = i*kmax/(N-1.0);
        k_111[i] = K_111.abs();
        K_110.x = i*kmax.x/(N-1.0);     K_110.y = i*kmax.y/(N-1.0);
        k_110[i] = K_110.abs();
        K_100.x = i*kmax.x/(N-1.0);
        k_100[i] = K_100.abs();
        K_010.y = i*kmax.y/(N-1.0);
        k_010[i] = K_010.abs();
        K_001.z = i*kmax.z/(N-1.0);
        k_001[i] = K_001.abs();
        
        G_111[i]=real(G_function(theCell,theNeighbours,K_111,w));
        G_110[i]=real(G_function(theCell,theNeighbours,K_110,w));
        G_100[i]=real(G_function(theCell,theNeighbours,K_100,w));
        G_010[i]=real(G_function(theCell,theNeighbours,K_010,w));
        G_001[i]=real(G_function(theCell,theNeighbours,K_001,w));
        
        if(!use_fixed_filter_width) {
            k_111[i]/= (kmax.abs()/sqrt(THREE));
            k_110[i]/= sqrt(sqr(kmax.x)+sqr(kmax.y))/sqrt(TWO);
            k_100[i]/= (kmax.x);
            k_010[i]/= (kmax.y);
            k_001[i]/= (kmax.z);            
        }

    }
    string title ;
    std::stringstream Cellstring, legend_111, legend_110, legend_100, legend_010, legend_001;
    Cellstring << "Cell = ("<<theCell.I<<","<<theCell.J<<","<<theCell.K<<")";
    if (use_fixed_filter_width) {
        legend_111 << "111  k_filt = " << fixed << setprecision(2) <<  Filter_Width(theCell,theNeighbours,K_111,w,G_cutoff);
        legend_110 << "110  k_filt = " << fixed << setprecision(2) <<  Filter_Width(theCell,theNeighbours,K_110,w,G_cutoff);
        legend_001 << "001  k_filt = " << fixed << setprecision(2) <<  Filter_Width(theCell,theNeighbours,K_001,w,G_cutoff);
        legend_010 << "010  k_filt = " << fixed << setprecision(2) <<  Filter_Width(theCell,theNeighbours,K_010,w,G_cutoff);
        legend_100 << "100  k_filt = " << fixed << setprecision(2) <<  Filter_Width(theCell,theNeighbours,K_100,w,G_cutoff);
    } else {
        legend_111 << "111  FGR = " << fixed << setprecision(2) <<  Filter_Grid_Ratio_111(theCell,theNeighbours,w,kmax);
        legend_110 << "110  FGR = " << fixed << setprecision(2) <<  Filter_Grid_Ratio_110(theCell,theNeighbours,w,kmax);
        legend_001 << "001  FGR = " << fixed << setprecision(2) <<  Filter_Grid_Ratio_001(theCell,theNeighbours,w,kmax);
        legend_010 << "010  FGR = " << fixed << setprecision(2) <<  Filter_Grid_Ratio_010(theCell,theNeighbours,w,kmax);
        legend_100 << "100  FGR = " << fixed << setprecision(2) <<  Filter_Grid_Ratio_100(theCell,theNeighbours,w,kmax);
    }
    
    title = "Transfer function " + filter_name() + " :   " + Cellstring.str() ;    
    
        
#ifdef _GNUPLOT
    Gnuplot_Control h1;
    h1.gnuplot_init(); 
    h1.gnuplot_setstyle("lines") ;
    h1.gnuplot_cmd("set grid");
    if(use_fixed_filter_width) {
        h1.gnuplot_set_xlabel("k");
    } else {
        h1.gnuplot_set_xlabel("k/kmax");
    }
    h1.gnuplot_set_ylabel("G(k)");
    h1.gnuplot_set_title(title);
    
    h1.gnuplot_plot1d_var2(k_111,G_111,N,legend_111.str().c_str());
    h1.gnuplot_plot1d_var2(k_110,G_110,N,legend_110.str().c_str());
    h1.gnuplot_plot1d_var2(k_001,G_001,N,legend_001.str().c_str());
    h1.gnuplot_plot1d_var2(k_010,G_010,N,legend_010.str().c_str());
    h1.gnuplot_plot1d_var2(k_100,G_100,N,legend_100.str().c_str());
#endif
    
    Output_Transfer_Function_gnuplot(prefix.str(),k_111,G_111,legend_111.str(),k_110,G_110,legend_110.str(),k_100,G_100,legend_100.str(),N,title);
    
    
    delete[] k_111;
    delete[] k_110;
    delete[] k_100;
    delete[] G_111;
    delete[] G_110;
    delete[] G_100;
    
    
    
    /* --------------- allocations ----------------- */
    Cell3D  ***kCells   = new Cell3D  **[N];
    Complex ***G        = new Complex **[N];
    for (int i=0; i<N; i++) {
        kCells[i]   = new Cell3D  *[N];
        G[i]        = new Complex *[N];
        for (int j=0; j<N; j++) {
            kCells[i][j]    = new Cell3D  [N];
            G[i][j]         = new Complex [N];
        }
    }
    
    /* --------------- assign kCells ------------------- */
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            for (int k=0; k<N; k++) {
                kCells[i][j][k].setindex(i,j,k);
                kCells[i][j][k].setloc(i*kmax.x/(N-ONE),j*kmax.y/(N-ONE),k*kmax.z/(N-ONE));
            }
        }
    }
    
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            for (int k=0; k<N; k++) {
                G[i][j][k] = G_function(theCell,theNeighbours,kCells[i][j][k].Xc,w);
            }
        }
    }
    
    Output_Transfer_Function_tecplot(prefix.str(),kCells,G,N);
    
    /* --------- deallocations --------- */
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            delete[] kCells[i][j];
            delete[] G[i][j];
        }
        delete[] kCells[i];
        delete[] G[i];
    }
    delete[] kCells;
    delete[] G;
    
}


template <typename Soln_pState, typename Soln_cState>
int Discrete_Filter<Soln_pState,Soln_cState>::
Output_Transfer_Function_gnuplot(string prefix_string, double *k_111, double *G_111, string legend_111, double *k_110, double *G_110, string legend_110, double *k_100, double *G_100, string legend_100, int N, string &title){
    int i, i_output_title;
    const char *prefix;
    char extension[256], output_file_name[256], gnuplot_file_name[256];
    char *output_file_name_ptr, *gnuplot_file_name_ptr;
    ofstream output_file, gnuplot_file;
    
    prefix = prefix_string.c_str();
    
    /* ------------------- output file ----------------------*/
    
    /* Determine output data file name for this processor. */
    sprintf(extension, "_gnuplot.dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;
    
    /* Open the output data file. */
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.fail()) return (1);
    
    
    output_file << setprecision(6);
    output_file.setf(ios::scientific);
    
    for (int i=0; i<N; i++) {
        output_file
        << " " << k_111[i]  << " " << G_111[i] 
        << " " << k_110[i]  << " " << G_110[i] 
        << " " << k_100[i]  << " " << G_100[i] 
        << "\n";
    } 
    
    output_file.unsetf(ios::scientific);
    output_file.flush();
    
    /* Close the output data file. */
    output_file.close();
    
    /* -------------------- gnuplot command file -------------------- */
    strcpy(extension, ".gplt");
    strcpy(gnuplot_file_name, prefix);
    strcat(gnuplot_file_name, extension);
    
    gnuplot_file_name_ptr = gnuplot_file_name;
    
    gnuplot_file.open(gnuplot_file_name_ptr, ios::out);
    if (gnuplot_file.fail()) return(1);
    
	std::string xlabel;
	
	if(use_fixed_filter_width) {
        xlabel = "k";
    } else {
        xlabel = "k/kmax";
    }
	
    gnuplot_file 
    << "set title \""<< title << "\"\n"
    << "set xlabel \"" << xlabel << "\"\n"
    << "set ylabel \"G(k)\"\n" 
    << "set grid \n"
    //<< "set logscale xy\n"
    << "plot \"" << output_file_name_ptr << "\" using 1:2 \\\n"
    << "     title \"" << legend_111    << "\" with lines , \\\n"
    << "\"" << output_file_name_ptr << "\" using 3:4 \\\n"
    << "     title \"" << legend_110    << "\" with lines , \\\n"
    << "\"" << output_file_name_ptr << "\" using 5:6 \\\n"
    << "     title \"" << legend_100    << "\" with lines \n"
    << "pause -1  \"Hit return to continue\"\n";
    
    gnuplot_file.close();
    return (0);
}


template <typename Soln_pState, typename Soln_cState>
int Discrete_Filter<Soln_pState,Soln_cState>::
Output_Transfer_Function_tecplot(string prefix_string, Cell3D  ***kCells, Complex ***G, int N) {
    
    int i, i_output_title;
    const char *prefix;
    char extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    
    
    prefix = prefix_string.c_str();
    
    /* Determine output data file name for this processor. */
    
    sprintf(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;
    
    /* Open the output data file. */
    
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.fail()) return (1);
    
    /* Write the solution data for each solution block. */
    
    output_file << setprecision(14);
    
    output_file << "TITLE = \"" << "Transfer function, "
    << "\"" << "\n"
    << "VARIABLES = "
    << "\"kx\" \\ \n"
    << "\"ky\" \\ \n"
    << "\"kz\" \\ \n"
    << "\"G real\" \\ \n"
    << "\"G imag\" \\ \n"
    << "\"G abs\" \\ \n";
    
    output_file<< "ZONE T =  \"Block Number = 0" 
    << "\" \\ \n"
    << "I = " << N << " \\ \n"
    << "J = " << N << " \\ \n"
    << "K = " << N << " \\ \n"
    << "DATAPACKING = POINT \n";
    
    
    for (int k=0; k<N; k++) {
        for (int j=0; j<N; j++) {
            for (int i=0; i<N ; i++) {
                output_file << " " << kCells[i][j][k].Xc 
                << " " << real(G[i][j][k]) 
                << " " << imag(G[i][j][k])
                << " " << abs(G[i][j][k])
                << "\n";
            } 
        } 
    }
    
    output_file << setprecision(6);
    
    
    
    /* Close the output data file. */
    
    output_file.close();
    
    /* Writing of output data files complete.  Return zero value. */
    
    return(0);
    
}
template <typename Soln_pState, typename Soln_cState>
Complex Discrete_Filter<Soln_pState,Soln_cState>::dG_function(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &k, RowVector &w) {    
    Complex dG(0,0);
    double abs_k = k.abs();
    if (abs_k == ZERO) return ZERO;
    Vector3D X0(theCell.Xc);
    Vector3D X;
    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
        X = theNeighbours.neighbour[n]->Xc;
        dG += w(n) * k/abs_k * (X-X0) * exp(-I*(k*(X-X0)));
    }
    return dG;
}

template <typename Soln_pState, typename Soln_cState>
Complex Discrete_Filter<Soln_pState,Soln_cState>::dG_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    
    Vector3D k(wave_number/sqrt(THREE),wave_number/sqrt(THREE),wave_number/sqrt(THREE));    
    return (dG_function(theCell, theNeighbours, k, w));
}



/* ------- calculate k_cutoff where G(k_cutoff) = G_cutoff -------- */
template <typename Soln_pState, typename Soln_cState>
Vector3D Discrete_Filter<Soln_pState,Soln_cState>::
Calculate_wavenumber_of_Gvalue(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w, double G_value) {
    
    Vector3D k_dir = kmax/kmax.abs();
    
    double a,b,m,s,fa,fb,fm,fp,p;
    Vector3D ka, kb, km, kp;
    int counter;
    a=ZERO, b=kmax.abs();
    ka = a*k_dir;
    kb = b*k_dir;
    fa = real(G_function(theCell,theNeighbours,ka,w)) - G_value ;
    fb = real(G_function(theCell,theNeighbours,kb,w)) - G_value ;
    fp = 100;
    counter = 0;
    while( fabs(fp) >= 0.001 ) {
        m = a + (b-a)/TWO;
        km = m*k_dir;
        fm = real(G_function(theCell,theNeighbours,km,w)) - G_value ;
        
        if (fa<ZERO)    s = -ONE;
        else            s = ONE;
        p = m + (m-a)*s*fm/(sqrt(fabs(fm*fm - fa*fb))+PICO);
        kp = p*k_dir;
        fp = real(G_function(theCell,theNeighbours,kp,w)) - G_value ;
        
        if (fa*fm < ZERO) { b = m; fb = fm; } else { a = m; fa = fm; }; 
        if (fa*fp < ZERO) { b = p; fb = fp; } else { a = p; fa = fp; }; 
        
        
        counter++;
        if(counter >= 50) {
            if (!batch_flag)
                cout << "max reached for Gvalue" << endl;
            break;
        }
    }
    double k_value = p;
    return k_value*k_dir;
    
}


template <typename Soln_pState, typename Soln_cState>
Vector3D Discrete_Filter<Soln_pState,Soln_cState>::
Calculate_wavenumber_of_dGvalue(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmin, Vector3D &kmax, RowVector &w, double dG_value) {
    Vector3D k_dir = kmax/kmax.abs();
    
    double a,b,m,s,fa,fb,fm,fp,p;
    Vector3D ka, kb, km, kp;
    
    int counter;
    if (kmin.abs()==kmax.abs()) {
        a=MILLI, b=kmax.abs();
    } else {
        a=kmin.abs()+MILLI, b=kmax.abs();
    }
    
    ka = a*k_dir;
    kb = b*k_dir;
    fa = real(dG_function(theCell,theNeighbours,ka,w)) - dG_value ;
    fb = real(dG_function(theCell,theNeighbours,kb,w)) - dG_value ;
    fp = 100;
    counter = 0;
    while( fabs(fp) >= 0.001 ) {
        m = a + (b-a)/TWO;
        km = m*k_dir;
        fm = real(dG_function(theCell,theNeighbours,km,w)) - dG_value ;
        
        if (fa<ZERO)    s = -ONE;
        else            s = ONE;
        
        p = m + (m-a)*s*fm/(sqrt(fabs(fm*fm - fa*fb))+PICO);
        kp = p*k_dir;
        
        fp = real(dG_function(theCell,theNeighbours,kp,w)) - dG_value ;
        
        if (fa*fm < ZERO) { b = m; fb = fm; } else { a = m; fa = fm; }; 
        if (fa*fp < ZERO) { b = p; fb = fp; } else { a = p; fa = fp; }; 
        
        
        counter++;
        if(counter >= 50) {
            if (!batch_flag)
                cout << "max reached for dGvalue" << endl;
            return kmax;
            break;
        }
    }
    double k_value = p;
    return k_value*k_dir;
}


template <typename Soln_pState, typename Soln_cState>
double Discrete_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax) {
    return Filter_Grid_Ratio_111(theCell,theNeighbours,w,kmax);
}


template <typename Soln_pState, typename Soln_cState>
double Discrete_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio_111(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax) {
    Vector3D k_cutoff = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,G_cutoff);
    
    // this k_cutoff is on the diagonal that connects (0,0,0) with kmax
    // the norm of the vector touching the ellipsoid connecting the boundaries of the spectral domain:
    //      kmax.abs()/sqrt(3)
    return (kmax.abs()/sqrt(THREE))/k_cutoff.abs();
}



template <typename Soln_pState, typename Soln_cState>
double Discrete_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio_110(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax) {
    Vector3D kmax_110;
    kmax_110.zero();
    kmax_110.x = kmax.x;
    kmax_110.y = kmax.y;
    
    Vector3D k_cutoff = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax_110,w,G_cutoff);
    
    return (sqrt(sqr(kmax.x)+sqr(kmax.y))/sqrt(TWO))/k_cutoff.abs();
}


template <typename Soln_pState, typename Soln_cState>
double Discrete_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio_100(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax) {
    Vector3D kmax_100;
    kmax_100.zero();
    kmax_100.x = kmax.x;
    
    Vector3D k_cutoff = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax_100,w,G_cutoff);
    
    return (kmax.x/k_cutoff.x);
}

template <typename Soln_pState, typename Soln_cState>
double Discrete_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio_010(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax) {
    Vector3D kmax_010;
    kmax_010.zero();
    kmax_010.y = kmax.y;
    
    Vector3D k_cutoff = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax_010,w,G_cutoff);
    
    return (kmax.y/k_cutoff.y);
}

template <typename Soln_pState, typename Soln_cState>
double Discrete_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio_001(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax) {
    Vector3D kmax_001;
    kmax_001.zero();
    kmax_001.z = kmax.z;
    
    Vector3D k_cutoff = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax_001,w,G_cutoff);
    
    return (kmax.z/k_cutoff.z);
}

template <typename Soln_pState, typename Soln_cState>
double Discrete_Filter<Soln_pState,Soln_cState>::Filter_Width(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kdir, RowVector &w, double G_cutoff) {
    Vector3D k_cutoff = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kdir,w,G_cutoff);
    return k_cutoff.abs();
}




#endif
