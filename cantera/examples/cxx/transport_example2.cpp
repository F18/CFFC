/////////////////////////////////////////////////////////////
//
//  mixture-averaged transport properties
//
//  $Author: dggoodwin $
//  $Revision: 1.5 $
//  $Date: 2005/07/22 12:50:27 $
//
//  copyright California Institute of Technology 2002
//
/////////////////////////////////////////////////////////////

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <cantera/Cantera.h>
#include <cantera/transport.h>
#include <time.h>
#include "example_utils.h"
#include <cantera/equilibrium.h>
#include <cantera/IdealGasMix.h>

template<class G, class V>
void makeTransportDataLabels(const G& gas, V& names) {
    int nsp = gas.nSpecies();
    names.resize(nsp + 3);
    names[0]  = "Temperature (K)";
    names[1]  = "Viscosity ()";
    names[2]  = "Thermal Conductivity (W/m-K)";
    int k;
    for (k = 0; k < nsp; k++) names[3+k] = gas.speciesName(k);
}

template<class G, class A>
void plotTransportSoln(string fname, string fmt, string title, const G& gas, 
    const A& soln) {
    vector<string> names;
    makeTransportDataLabels(gas, names);
    writePlotFile(fname, fmt, title, names, soln);
}


int transport_example2(int job) {

    try {

        cout << "Multicomponent transport properties." << endl;
        if (job > 0) {
            cout << "Viscosity, thermal conductivity, and thermal diffusion\n"
                " coefficients at 2 atm for a "
                 << "range of temperatures" << endl;
        }
        if (job <= 1) return 0;

        // header
        writeCanteraHeader(cout);


        // create a gas mixture, and set its state

        IdealGasMix gas("gri30.cti", "gri30");
        doublereal temp = 2000.0;
        doublereal pres = 2.0*OneAtm;
        gas.setState_TPX(temp, pres, "H2:1.0, O2:0.5, CH4:0.1, N2:0.2");
        equilibrate(gas,"TP");


        // create a transport manager that implements
        // multicomponent transport properties

        Transport* tr = newTransportMgr("Multi",  &gas);
        int nsp = gas.nSpecies();


        // create a 2D array to hold the outputs
        int ntemps = 20;
        Array2D output(nsp+3, ntemps);

        // main loop
        clock_t t0 = clock();
        for (int i = 0; i < ntemps; i++) {
            temp = 500.0 + 100.0*i;
            gas.setState_TP(temp, pres);
            output(0,i) = temp;
            output(1,i) = tr->viscosity();
            output(2,i) = tr->thermalConductivity();
            tr->getThermalDiffCoeffs(&output(3,i));
        }
        clock_t t1 = clock();

        // make a Tecplot data file and an Excel spreadsheet
        string plotTitle = "transport example 2: "
                           "multicomponent transport properties";
        plotTransportSoln("tr2.dat", "TEC", plotTitle, gas, output);
        plotTransportSoln("tr2.csv", "XL", plotTitle, gas, output);

        // print final temperature and timing data
        doublereal tmm = 1.0*(t1 - t0)/CLOCKS_PER_SEC;
        cout << " time = " << tmm << endl << endl;

        cout << "Output files:" << endl
             << "  tr2.csv    (Excel CSV file)" << endl
             << "  tr2.dat    (Tecplot data file)" << endl;

        return 0;
    }

    // handle exceptions thrown by Cantera
    catch (CanteraError) {
        showErrors(cout);
        cout << " terminating... " << endl;
        appdelete();
        return -1;
    }
}
