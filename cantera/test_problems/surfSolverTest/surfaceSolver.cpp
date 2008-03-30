/**
 *  @file surfaceSolver.cpp
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
//  Example 
//
//  Read a surface growth mechanism and calculate the solution
//  using Placid.
//

#include <iostream>
#include <string>
#include <vector>

#define MSSIZE 200
using namespace std;

#ifdef DEBUG_HKM
int iDebug_HKM = 0;
#endif

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
static void printUsage()
{

}

#include "Cantera.h"
#include "Interface.h"

#include "kinetics.h"
#include "kernel/ImplicitSurfChem.h"
#include "kernel/solveSP.h"

using namespace Cantera;


void printGas(ThermoPhase *gasTP, InterfaceKinetics * iKin_ptr, double *src) {
  double x[MSSIZE];
  double C[MSSIZE];
  string gasPhaseName          = "gas";
  gasTP->getMoleFractions(x);
  gasTP->getConcentrations(C);
  double Temp = gasTP->temperature();
  double p = gasTP->pressure();
  cout << "Gas Temperature = " << Temp << endl;
  cout << "Gas Pressure    = " << p << endl;
  int kstart = iKin_ptr->kineticsSpeciesIndex(0, 0);
  cout << "Gas Phase:  " << gasPhaseName << "   " 
       << "(" << kstart << ")" << endl;
  cout << "                       Name      " 
       << "     Conc              MoleF       SrcRate " << endl;
  cout << "                                 "
       << "   (kmol/m^3)                   (kmol/m^2/s) " << endl;
  double sum = 0.0;
  int nspGas = gasTP->nSpecies();
  for (int k = 0; k < nspGas; k++) {
    kstart = iKin_ptr->kineticsSpeciesIndex(k, 0);
    printf("%4d %24s   %14g %14g  %14e\n", 
	   k, gasTP->speciesName(k).c_str(),
	   C[k], x[k], src[kstart]);
    sum += x[k];
  }
  cout << "Sum of gas mole fractions= " << sum << endl;
  cout << endl;
}

void printBulk(ThermoPhase *bulkPhaseTP, InterfaceKinetics * iKin_ptr, double *src) {
  double x[MSSIZE];
  double C[MSSIZE];
  string bulkParticlePhaseName = bulkPhaseTP->id();
  bulkPhaseTP->getMoleFractions(x);
  bulkPhaseTP->getConcentrations(C);
  int kstart = iKin_ptr->kineticsSpeciesIndex(0, 1);
  double dens = bulkPhaseTP->density();
  cout << "Bulk Phase:  " << bulkParticlePhaseName << "   " 
       << "(" << kstart << ")" << endl;
  double Temp = bulkPhaseTP->temperature();
  double p = bulkPhaseTP->pressure();
  cout << "Bulk Temperature = " << Temp << endl;
  cout << "Bulk Pressure    = " << p << endl;
  cout << "                       Name      " 
       << "     Conc              MoleF       SrcRate " << endl;
  cout << "                                 "
       << "   (kmol/m^3)                   (kmol/m^2/s) " << endl;
  double sum = 0.0;
  double Wsum = 0.0;
  const array_fp& molecW = bulkPhaseTP->molecularWeights();
  int nspBulk = bulkPhaseTP->nSpecies();
  for (int k = 0; k < nspBulk; k++) {
    kstart = iKin_ptr->kineticsSpeciesIndex(k, 1);
    printf("%4d %24s   %14g %14g  %14e\n", 
	   k, bulkPhaseTP->speciesName(k).c_str(),
	   C[k], x[k], src[kstart]);
    sum += x[k];
    Wsum += src[kstart] * molecW[k];
  }
  cout << "Bulk Weight Growth Rate = " << Wsum << " kg/m^2/s" << endl;
  double gr = Wsum / dens;
  cout << "Bulk Growth Rate = " << gr << " m/s" << endl;
  cout << "Bulk Growth Rate = " << gr * 1.0E6 * 3600.
       << " microns / hour" << endl;
  cout << "Density of bulk phase = " << dens << " kg / m^3 "<< endl;
  cout << "                      = " << dens / 1.0E3 
       <<" gm / cm^3 " << endl;
  cout << "Sum of bulk mole fractions= " << sum << endl;
  cout << endl;
}

void printSurf(ThermoPhase *surfPhaseTP, InterfaceKinetics * iKin_ptr, double *src) {
  double x[MSSIZE];
  string surfParticlePhaseName = surfPhaseTP->id();
  surfPhaseTP->getMoleFractions(x);
  int kstart = iKin_ptr->kineticsSpeciesIndex(0, 2);
  cout << "Surface Phase:  " << surfParticlePhaseName 
       << " (" << kstart << ")" << endl;
  double Temp = surfPhaseTP->temperature();
  double p = surfPhaseTP->pressure();
  cout << "Surface Temperature = " << Temp << endl;
  cout << "Surface Pressure    = " << p << endl;
  cout << "                       Name      " 
       << "   Coverage         SrcRate " << endl;
  double sum = 0.0;
 int nspSurf = surfPhaseTP->nSpecies();
  for (int k = 0; k < nspSurf; k++) {
    kstart = iKin_ptr->kineticsSpeciesIndex(0, 2);
    printf("%4d %24s   %14g   %14e\n", 
	   k, surfPhaseTP->speciesName(k).c_str(),
	   x[k], src[kstart]);
    sum += x[k];
  }
  cout << "Sum of coverages = " << sum << endl;
}

int main(int argc, char** argv) {
    string infile;
    int ioflag = 1;
    int i, k;
    // look for command-line options
    if (argc > 1) {
      string tok;
      for (int j = 1; j < argc; j++) {
	tok = string(argv[j]);
	if (tok[0] == '-') {
	  int nopt = static_cast<int>(tok.size());
	  for (int n = 1; n < nopt; n++) {
	    if (tok[n] == 'h') {
	      printUsage();
	      exit(0);
            } else if (tok[n] == 'd') {
              int lvl = 0;
              if (j < (argc - 1)) {
                string tokla = string(argv[j+1]);
                if (strlen(tokla.c_str()) > 0) {
                  lvl = atoi(tokla.c_str());
                  n = nopt - 1;
                  j += 1;
                  ioflag = lvl;
               }
             }
	    } else {
	      printUsage();
	      exit(1);
	    }
	  }
	} else if (infile == "") {
	  infile = tok;
	}
	else {
	  printUsage();
	  exit(1);
	}
      }
    }
    if (infile == "") {
      infile = "diamond.cti";
    }
 
    try {
      /*************************************************************/

      /*
       *  FILL IN THESE NAMES FOR EACH PROBLEM
       */
      /*
       * ProblemNumber = 0 : diamond.cti
       *               = 1 : haca.cti
       */
      int ProblemNumber = 1;
      string gasPhaseName          = "gas";
      string bulkParticlePhaseName = "diamond";
      string surfParticlePhaseName = "diamond_100";
      if (ProblemNumber == 1) { 
	gasPhaseName          = "gas";
	bulkParticlePhaseName = "soot";
	surfParticlePhaseName = "soot_interface";
      }

      /************************************************************/
      XML_Node *xc = new XML_Node();
      string path = findInputFile(infile);
      ctml::get_CTML_Tree(xc, path);

      XML_Node * const xg = (XML_Node *) findXMLPhase(xc, gasPhaseName);
      if (!xg) {
	printf("ERROR: Could not find gas phase named, %s, in file\n",
	       gasPhaseName.c_str());
        exit(-1);
      }
      ThermoPhase *gasTP = newPhase(*xg);
      int nspGas = gasTP->nSpecies();
      cout << "Number of species = " << nspGas << endl;

      XML_Node * const xd = 
	  (XML_Node *) findXMLPhase(xc, bulkParticlePhaseName);
      if (!xd) {
	printf("ERROR: Could not find bulk phase named, %s, in file\n",
	       bulkParticlePhaseName.c_str());
        exit(-1);
      }
      ThermoPhase *bulkPhaseTP = newPhase(*xd);
      int nspBulk = bulkPhaseTP->nSpecies();
      cout << "Number of species in bulk phase named " <<
	  bulkParticlePhaseName << " = " << nspBulk << endl;


      XML_Node * const xs =
	  (XML_Node *) findXMLPhase(xc, surfParticlePhaseName);
      if (!xs) {
	printf("ERROR: Could not find surf Particle phase named, %s, in file\n",
	       surfParticlePhaseName.c_str());
        exit(-1);
      }
      ThermoPhase *surfPhaseTP = newPhase(*xs);
      int nsp_d100 = surfPhaseTP->nSpecies();
      cout << "Number of species in surface phase, " << surfParticlePhaseName 
	   << " = " << nsp_d100 << endl;

      vector<ThermoPhase *> phaseList;
      phaseList.push_back(gasTP);     
      phaseList.push_back(bulkPhaseTP);
      phaseList.push_back(surfPhaseTP);

      InterfaceKinetics *iKin_ptr = new InterfaceKinetics();
      importKinetics(*xs, phaseList, iKin_ptr);
      int nr = iKin_ptr->nReactions();
      cout << "Number of reactions = " << nr << endl;

      double x[MSSIZE], p = OneAtm;

      /*
       * Set the Gas State: 
       * -> note that the states are set in the xml files too
       */
      for (i = 0; i < MSSIZE; i++) x[i] = 0.0;
      if (ProblemNumber == 0) {
	x[0] = 0.0010;
	x[1] = 0.9888;
	x[2] = 0.0002;
	x[3] = 0.0100;
	p = 20.0*OneAtm/760.0;
	gasTP->setState_TPX(1200., p, x);
      }

      /*
       * Set the surface initial state
       */
      for (i = 0; i < MSSIZE; i++) x[i] = 0.0;
      if (ProblemNumber == 0) {
	int i0 = surfPhaseTP->speciesIndex("c6H*");
	if (i0 >= 0) {
	  x[i0] = 0.1;
	}
	int i1 = surfPhaseTP->speciesIndex("c6HH");
	if (i1 >= 0) {
	  x[i1] = 0.9;
	}
	surfPhaseTP->setState_TX(1200., x);
      }

      /*
       * Set the bulk Phase State
       */
      for (i = 0; i < MSSIZE; i++) x[i] = 0.0;
      if (ProblemNumber == 0) { 
	x[0] = 1.0;
	bulkPhaseTP->setState_TPX(1200., p, x);
      }

      iKin_ptr->setIOFlag(ioflag);
      /*
       *  Solve the Equation system
       */
      //iKin_ptr->advanceCoverages(100.);
      iKin_ptr->solvePseudoSteadyStateProblem();
   
      /*
       * Download the source terms for the rate equations
       */
      double src[MSSIZE];
      iKin_ptr->getNetProductionRates(src);

      double sum = 0.0;
      if (ProblemNumber == 0) {
	double naH;
	for (k = 0; k < 13; k++) {
	  if (k < 4) {
	    naH = gasTP->nAtoms(k, 0);
	  } else if (k == 4) {
	    naH = 0;
	  } else if (k > 4) {
	    int itp = k - 5;
	    naH = surfPhaseTP->nAtoms(itp, 0);
	  }
	  cout << k << "  " << naH << "  " ;
	  if (fabs(src[k]) < 2.0E-17) {
	    cout << " nil" << endl;
	  } else {
	    cout << src[k] << endl;
	  }
	  sum += naH * src[k];
	}
	cout << "sum = " << sum << endl;
      }
  
   
      printGas(gasTP, iKin_ptr, src);
      printBulk(bulkPhaseTP, iKin_ptr, src);
      printSurf(surfPhaseTP, iKin_ptr, src) ;

      /*****************************************************************************/
      /*  Now Tweak the inputs and do a quick calculation */
      /****************************************************************************/

      /*
       * Set the Gas State: 
       * -> note that the states are set in the xml files too
       */
      double pres = gasTP->pressure();
      gasTP->getMoleFractions(x);
      double tmp = 0.3 * x[0];
      double tmp2 = 0.3 * x[1];
      if (tmp2 < tmp) tmp = tmp2;
      x[0] += tmp;
      x[1] -= tmp;
      gasTP->setState_PX(pres, x);
      
      iKin_ptr->solvePseudoSteadyStateProblem();
      iKin_ptr->getNetProductionRates(src);

      printGas(gasTP, iKin_ptr, src);
      printBulk(bulkPhaseTP, iKin_ptr, src);
      printSurf(surfPhaseTP, iKin_ptr, src) ;


      /*****************************************************************************/
      /*  Now Tweak the inputs and do a quick calculation */
      /****************************************************************************/

      /*
       * Set the Gas State: 
       * -> note that the states are set in the xml files too
       */

      /*
       * Set the Gas State: 
       * -> note that the states are set in the xml files too
       */
      pres = gasTP->pressure();
      double temp = gasTP->temperature();
      temp += 95;
      gasTP->setState_TP(temp, pres);
      
      iKin_ptr->solvePseudoSteadyStateProblem();
      iKin_ptr->getNetProductionRates(src);

      printGas(gasTP, iKin_ptr, src);
      printBulk(bulkPhaseTP, iKin_ptr, src);
      printSurf(surfPhaseTP, iKin_ptr, src) ;


      /*****************************************************************************/
      /*  Now Don't Tweak the inputs at all */
      /****************************************************************************/
      gasTP->setState_TP(temp, pres);
  
      iKin_ptr->solvePseudoSteadyStateProblem();
      iKin_ptr->getNetProductionRates(src);

      printGas(gasTP, iKin_ptr, src);
      printBulk(bulkPhaseTP, iKin_ptr, src);
      printSurf(surfPhaseTP, iKin_ptr, src) ;

      delete iKin_ptr; 
      delete gasTP; gasTP = 0;
      delete bulkPhaseTP; bulkPhaseTP = 0;
      delete surfPhaseTP; surfPhaseTP = 0;
      delete xc;
      appdelete();
    }
    catch (CanteraError) {
      showErrors(cout);
    }

    return 0;
}
/***********************************************************/