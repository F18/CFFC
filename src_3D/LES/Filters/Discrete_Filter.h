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

#include <complex>

#include "../../Math/Math.h"
#include "../../Math/LinearSystems.h"
#include "../../Math/Matrix.h"

#ifndef _GNUPLOT
#include "../../System/gnuplot.h"
#endif

#include "../../Utilities/Utilities.h"


#include <cstdlib>

template <typename Soln_pState, typename Soln_cState>
class LES_Filter;

typedef complex<double> Complex;

/**
 * CLASS: Discrete_Filter
 */
template <typename Soln_pState, typename Soln_cState>
class Discrete_Filter : public General_Filter<Soln_pState,Soln_cState> {
public:
    
    Discrete_Filter(void) {
        FGR = LES_Filter<Soln_pState,Soln_cState>::FGR;
        commutation_order = LES_Filter<Soln_pState,Soln_cState>::commutation_order;
        number_of_rings = LES_Filter<Soln_pState,Soln_cState>::number_of_rings;
        target_filter_sharpness = LES_Filter<Soln_pState,Soln_cState>::target_filter_sharpness;
        theNeighbours.allocate(number_of_rings+1);   // shouldn't be +1!!!! look up why
    }
    
    Neighbours theNeighbours;
    int commutation_order;
    double FGR;
    int number_of_rings;
    double target_filter_sharpness;

    void Allocate_Filter_Weights(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk);
    void Reset_Filter_Weights(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk);

    DenseMatrix Neighbouring_Values;
    void Set_Neighbouring_Values(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Neighbours &theNeighbours, DenseMatrix &Neighbouring_Values);

    RowVector filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);

    
    
    /* --------------- Transfer function calculations ------------- */
    void transfer_function(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);
    Complex G_function(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &k, RowVector &w);
    Complex G_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    Complex G_function_1D_100(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    Complex G_function_1D_110(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    Complex dG_function(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &k, RowVector &w);
    Complex dG_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    
    int Output_Transfer_Function_tecplot(Cell3D  ***kCells, Complex ***G, int N);
    int Output_Transfer_Function_gnuplot(double *k_111, double *G_111, double *k_110, double *G_110, double *K_100, double *G_100, int N, string &title);

    double Filter_Grid_Ratio(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax);
    Vector3D Calculate_wavenumber_of_Gvalue(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w, double G_value);
    Vector3D Calculate_wavenumber_of_dGvalue(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmin, Vector3D &kmax, RowVector &w, double dG_value);

    
    virtual void Get_Neighbours(Cell3D &theCell) = 0;
    virtual RowVector Get_Weights(Cell3D &theCell, Neighbours &theNeighbours) = 0;
    virtual string filter_name(void) = 0;
    
    void Write_to_file(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, ofstream &out_file);
    void Read_from_file(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, ifstream &in_file);
    
    
};

template<typename Soln_pState, typename Soln_cState>
inline RowVector Discrete_Filter<Soln_pState,Soln_cState>::filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {

    if (!SolnBlk.Filter_Weights_Allocated) {
        Allocate_Filter_Weights(SolnBlk);
        SolnBlk.Filter_Weights_Allocated = true;
    }
    
    
    theNeighbours.set_grid(SolnBlk.Grid);
    Get_Neighbours(theCell);
    
    int I(theCell.I);
    int J(theCell.J);
    int K(theCell.K);
    if (!SolnBlk.Filter_Weights_Assigned[I][J][K]) {
        SolnBlk.Filter_Weights[I][J][K] = Get_Weights(theCell,theNeighbours);
        SolnBlk.Filter_Weights_Assigned[I][J][K] = true;
    }
    Set_Neighbouring_Values(SolnBlk,theNeighbours,Neighbouring_Values);
    
    return SolnBlk.Filter_Weights[I][J][K]*Neighbouring_Values;
    
}


template <typename Soln_pState, typename Soln_cState>
void Discrete_Filter<Soln_pState,Soln_cState>::Set_Neighbouring_Values(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Neighbours &theNeighbours, DenseMatrix &Neighbouring_Values) {
    RowVector RowVec;
    LES_Filter<Soln_pState,Soln_cState>::what_to_filter(SolnBlk,0,0,0,RowVec);
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    int the_number_of_variables = RowVec.size();
    
    if (Neighbouring_Values.get_n()!=the_number_of_neighbours || Neighbouring_Values.get_m()!=the_number_of_variables) {
        Neighbouring_Values.newsize(the_number_of_neighbours,the_number_of_variables);
    }
    
    int I,J,K;
    for (int i=0; i<the_number_of_neighbours; i++) {
        I = theNeighbours.neighbour[i].I;
        J = theNeighbours.neighbour[i].J;
        K = theNeighbours.neighbour[i].K;
        LES_Filter<Soln_pState,Soln_cState>::what_to_filter(SolnBlk,I,J,K,Neighbouring_Values,i);
    }    
}

template<typename Soln_pState, typename Soln_cState>
inline void Discrete_Filter<Soln_pState,Soln_cState>::Allocate_Filter_Weights(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk) {
    SolnBlk.Filter_Weights = new RowVector **[SolnBlk.NCi];
    SolnBlk.Filter_Weights_Assigned = new bool **[SolnBlk.NCi];
    for (int i=0; i<SolnBlk.NCi; i++) {
        SolnBlk.Filter_Weights[i] = new RowVector *[SolnBlk.NCj];
        SolnBlk.Filter_Weights_Assigned[i] = new bool *[SolnBlk.NCj];
        for (int j=0; j<SolnBlk.NCj; j++) {
            SolnBlk.Filter_Weights[i][j] = new RowVector [SolnBlk.NCk];
            SolnBlk.Filter_Weights_Assigned[i][j] = new bool [SolnBlk.NCk];
            
            for (int k=0; k<SolnBlk.NCk; k++) {
                SolnBlk.Filter_Weights_Assigned[i][j][k] = false;
            }
        }
    }
}

template<typename Soln_pState, typename Soln_cState>
inline void Discrete_Filter<Soln_pState,Soln_cState>::Reset_Filter_Weights(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk) {
    SolnBlk.Filter_Weights_Assigned = new bool **[SolnBlk.NCi];
    for (int i=0; i<SolnBlk.NCi; i++) {
        SolnBlk.Filter_Weights_Assigned[i] = new bool *[SolnBlk.NCj];
        for (int j=0; j<SolnBlk.NCj; j++) {
            SolnBlk.Filter_Weights_Assigned[i][j] = new bool [SolnBlk.NCk];
            
            for (int k=0; k<SolnBlk.NCk; k++) {
                SolnBlk.Filter_Weights_Assigned[i][j][k] = false;
            }
        }
    }
}


template<typename Soln_pState, typename Soln_cState>
inline void Discrete_Filter<Soln_pState,Soln_cState>::transfer_function(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {
    
    
    theNeighbours.set_grid(SolnBlk.Grid);
    
    Get_Neighbours(theCell);
    
    RowVector w = Get_Weights(theCell, theNeighbours);
    
    Vector3D kmax;
    Vector3D Delta = theNeighbours.Delta;
    kmax.x = PI/Delta.x;
    kmax.y = PI/Delta.y;
    kmax.z = PI/Delta.z;
    
    cout << "\nkmax = " << kmax << endl;
    
    
    /* --------------- allocations ----------------- */
    int N=50;
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
    

    double *k_111 = new double [N];
    double *k_110 = new double [N];
    double *k_100 = new double [N];
    double *G_111 = new double [N];
    double *G_110 = new double [N];
    double *G_100 = new double [N];
    
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            for (int k=0; k<N; k++) {
                G[i][j][k] = G_function(theCell,theNeighbours,kCells[i][j][k].Xc,w);
            }
        }
    }
    
    
    Output_Transfer_Function_tecplot(kCells,G,N);
    
    
    Vector3D K_111, K_110, K_100;
    K_110.zero(); K_100.zero();
    for (int i=0; i<N; i++) {
        K_111 = i*kmax/(N-1.0);
        k_111[i] = K_111.abs();
        K_110.x = i*kmax.x/(N-1.0);     K_110.y = i*kmax.y/(N-1.0);
        k_110[i] = K_110.abs();
        K_100.x = i*kmax.x/(N-1.0);
        k_100[i] = K_100.abs();
        G_111[i]=real(G_function(theCell,theNeighbours,K_111,w));
        G_110[i]=real(G_function(theCell,theNeighbours,K_110,w));
        G_100[i]=real(G_function(theCell,theNeighbours,K_100,w));
        k_111[i]/= (kmax.abs()/sqrt(THREE));
        k_110[i]/= (kmax.abs()/sqrt(THREE));
        k_100[i]/= (kmax.abs()/sqrt(THREE));
    }
    double FGR = Filter_Grid_Ratio(theCell,theNeighbours,w,kmax);
    string title;
    std::stringstream ss ;
    ss << fixed << setprecision(2) << FGR ;
    title = "Transfer function " + filter_name() + " :  FGR = " + ss.str();    
    
#ifdef _GNUPLOT
    Gnuplot_Control h1;
    h1.gnuplot_init(); 
    h1.gnuplot_setstyle("lines") ;
    h1.gnuplot_cmd("set grid");
    h1.gnuplot_set_xlabel("k");
    h1.gnuplot_set_ylabel("G(k)");
    h1.gnuplot_set_title(title);
    
    h1.gnuplot_plot1d_var2(k_111,G_111,N,"111");
    h1.gnuplot_plot1d_var2(k_110,G_110,N,"110");
    h1.gnuplot_plot1d_var2(k_100,G_100,N,"100");
#endif
    
    Output_Transfer_Function_gnuplot(k_111,G_111,k_110,G_110,k_100,G_100,N,title);

}


template<typename Soln_pState, typename Soln_cState>
int Discrete_Filter<Soln_pState,Soln_cState>::
Output_Transfer_Function_gnuplot(double *k_111, double *G_111, double *k_110, double *G_110, double *k_100, double *G_100, int N, string &title){
    int i, i_output_title;
    char prefix[256], extension[256], output_file_name[256], gnuplot_file_name[256];
    char *output_file_name_ptr, *gnuplot_file_name_ptr;
    ofstream output_file, gnuplot_file;
    
    /* Determine prefix of output data file names. */
    
    //    i = 0;
    //    while (1) {
    //        if (Input.Output_File_Name[i] == ' ' ||
    //            Input.Output_File_Name[i] == '.') break;
    //        prefix[i]=Input.Output_File_Name[i];
    //        i = i + 1;
    //        if (i > strlen(Input.Output_File_Name) ) break;
    //    } /* endwhile */
    //    prefix[i] = '\0';
    sprintf(prefix,"transfer_function");
    
    
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
    
    gnuplot_file 
    << "set title \""<< title << "\"\n"
    << "set xlabel \"k Delta / pi \"\n"
    << "set ylabel \"G(k)\"\n" 
    << "set grid \n"
    //<< "set logscale xy\n"
    << "plot \"" << output_file_name_ptr << "\" using 1:2 \\\n"
    << "     title \"" << "111"    << "\" with lines , \\\n"
    << "\"" << output_file_name_ptr << "\" using 3:4 \\\n"
    << "     title \"" << "110"    << "\" with lines , \\\n"
    << "\"" << output_file_name_ptr << "\" using 5:6 \\\n"
    << "     title \"" << "100"    << "\" with lines \n"
    << "pause -1  \"Hit return to continue\"\n";
    
    gnuplot_file.close();
    return (0);
}


template<typename Soln_pState, typename Soln_cState>
int Discrete_Filter<Soln_pState,Soln_cState>::
Output_Transfer_Function_tecplot(Cell3D  ***kCells, Complex ***G, int N) {
    
    int i, i_output_title;
    char prefix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    
    
    /* Determine prefix of output data file names. */
    
    //    i = 0;
    //    while (1) {
    //        if (Input.Output_File_Name[i] == ' ' ||
    //            Input.Output_File_Name[i] == '.') break;
    //        prefix[i]=Input.Output_File_Name[i];
    //        i = i + 1;
    //        if (i > strlen(Input.Output_File_Name) ) break;
    //    } /* endwhile */
    //    prefix[i] = '\0';
    sprintf(prefix,"transfer_function");
    
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



template<typename Soln_pState, typename Soln_cState>
inline Complex Discrete_Filter<Soln_pState,Soln_cState>::G_function(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &k, RowVector &w) {
    Complex G(0,0);

    Vector3D X0(theCell.Xc);
    Vector3D X;
    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
        X = theNeighbours.neighbour[n].Xc;
        G += w(n) * exp(-I*(k*(X-X0)));
    }
    return G;
}

template<typename Soln_pState, typename Soln_cState>
inline Complex Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    Vector3D k(wave_number/sqrt(THREE),wave_number/sqrt(THREE),wave_number/sqrt(THREE));    
    return (G_function(theCell, theNeighbours, k, w));
}


template<typename Soln_pState, typename Soln_cState>
inline Complex Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D_100(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    Vector3D k(wave_number,ZERO,ZERO);    
    return (G_function(theCell, theNeighbours, k, w));
}

template<typename Soln_pState, typename Soln_cState>
inline Complex Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D_110(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    Vector3D k(wave_number/sqrt(TWO),wave_number/sqrt(TWO),ZERO);    
    return (G_function(theCell, theNeighbours, k, w));
}


template<typename Soln_pState, typename Soln_cState>
inline Complex Discrete_Filter<Soln_pState,Soln_cState>::dG_function(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &k, RowVector &w) {
//    Complex G(0,0);
//    double x,y,z,x0,y0,z0;
//    double abs_k = kCell.Xc.abs();
//    if (abs_k == ZERO) return ZERO;
//    x0 = theCell.Xc.x;
//    y0 = theCell.Xc.y;
//    z0 = theCell.Xc.z;
//    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
//        x=theNeighbours.neighbour[n].Xc.x;
//        y=theNeighbours.neighbour[n].Xc.y;
//        z=theNeighbours.neighbour[n].Xc.z;
//        G += w(n) * (kCell.Xc.x/abs_k*(x0-x) + kCell.Xc.y/abs_k*(y0-y) + kCell.Xc.z/abs_k*(z0-z)) * exp(I*(kCell.Xc.x*(x0-x) + kCell.Xc.y*(y0-y) + kCell.Xc.z*(z0-z)));
//    }
//    return G;
    
    Complex dG(0,0);
    double abs_k = k.abs();
    if (abs_k == ZERO) return ZERO;
    Vector3D X0(theCell.Xc);
    Vector3D X;
    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
        X = theNeighbours.neighbour[n].Xc;
        dG += w(n) * k/abs_k * (X-X0) * exp(-I*(k*(X-X0)));
    }
    return dG;
}

template<typename Soln_pState, typename Soln_cState>
inline Complex Discrete_Filter<Soln_pState,Soln_cState>::dG_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    
    Vector3D k(wave_number/sqrt(THREE),wave_number/sqrt(THREE),wave_number/sqrt(THREE));    
    return (dG_function(theCell, theNeighbours, k, w));
}



/* ------- calculate k_HALF where G(k_HALF) = 0.5 -------- */
template<typename Soln_pState, typename Soln_cState>
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
            cout << "max reached for Gvalue" << endl;
            break;
        }
    }
    double k_value = p;
    return k_value*k_dir;
    
}


template<typename Soln_pState, typename Soln_cState>
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
            cout << "max reached for dGvalue" << endl;
            return kmax;
            break;
        }
    }
    double k_value = p;
    
    cout << "dG_k_value = " << p << endl;
    return k_value*k_dir;
}


template<typename Soln_pState, typename Soln_cState>
double Discrete_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, Vector3D &kmax) {
    Vector3D k_HALF = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,HALF);
    
    // this k_HALF is on the diagonal that connects (0,0,0) with kmax
    // the norm of the vector touching the ellipsoid connecting the boundaries of the spectral domain:
    //      kmax.abs()/sqrt(3)
    return (kmax.abs()/sqrt(THREE))/k_HALF.abs();
}


template<typename Soln_pState, typename Soln_cState>
void Discrete_Filter<Soln_pState,Soln_cState>::Write_to_file(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, ofstream &out_file) {
    for (int i=0; i<SolnBlk.NCi; i++) {
        for (int j=0; j<SolnBlk.NCj; j++) {                
            for (int k=0; k<SolnBlk.NCk; k++) {
                out_file << SolnBlk.Filter_Weights_Assigned[i][j][k] << " ";
                if (SolnBlk.Filter_Weights_Assigned[i][j][k]) {
                    SolnBlk.Filter_Weights[i][j][k].write(out_file);
                } else {
                    out_file << "\n";
                }
            }
        }
    }
    out_file << "\n"; // extra line to separate SolnBlks
}

template<typename Soln_pState, typename Soln_cState>
void Discrete_Filter<Soln_pState,Soln_cState>::Read_from_file(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, ifstream &in_file) {
    Allocate_Filter_Weights(SolnBlk);
    for (int i=0; i<SolnBlk.NCi; i++) {
        for (int j=0; j<SolnBlk.NCj; j++) {                
            for (int k=0; k<SolnBlk.NCk; k++) {
                in_file.setf(ios::skipws);
                in_file >> SolnBlk.Filter_Weights_Assigned[i][j][k];
                in_file.unsetf(ios::skipws);
                if (SolnBlk.Filter_Weights_Assigned[i][j][k]) {
                    SolnBlk.Filter_Weights[i][j][k].read(in_file);
                }
            }
        }
    }
    SolnBlk.Filter_Weights_Allocated = true;
}


#endif