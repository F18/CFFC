/////////////////////////////////////////////////////////////////////
///
/// \file RadiatingGas.h
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the struct definition for 
///        describing the state of a participating 
///        (radiation absorbing/emitting) gas.
///
/////////////////////////////////////////////////////////////////////
#ifndef _RADIATING_GAS_INCLUDED 
#define _RADIATING_GAS_INCLUDED


class RadiatingGas {

public:
  double p;        //!< pressure [atm]
  double T;        //!< temperature [K]
  double xco;      //!< mole fraction of CO
  double xh2o;     //!< mole fraction of H2O
  double xco2;     //!< mole fraction of CO2
  double xo2;      //!< mole fraction of O2
  double fsoot;    //!< volume fraction of soot 

  //! default constructors
  RadiatingGas () : p(0.0), T(0.0), xco(0.0), xh2o(0.0), 
		    xco2(0.0), xo2(0.0), fsoot(0.0) {};

  RadiatingGas ( const double &pp,
		 const double &tt,
		 const double &xCO,
		 const double &xH2O,
		 const double &xCO2,
		 const double &xO2,
		 const double &xSOOT ) : p(pp), T(tt), xco(xCO), 
					 xh2o(xH2O), xco2(xCO2), 
					 xo2(xO2), fsoot(xSOOT) {};

  // MPI broadcast operator
#ifdef _MPI_VERSION
  void Broadcast( MPI::Intracomm &Communicator = MPI::COMM_WORLD,
		  const int Source_Rank = 0) {
    Communicator.Bcast(&(p),     1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(T),     1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(xco),   1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(xh2o),  1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(xco2),  1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(xo2),   1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(fsoot), 1, MPI::DOUBLE, Source_Rank);
  }
#endif


  // output operator
  friend ostream& operator << (ostream &out, const RadiatingGas &g) {
    out << "p = " << g.p << endl;
    out << "T = " << g.T << endl;
    out << "xCO = " << g.xco << endl;
    out << "xH2O = " << g.xh2o << endl;
    out << "xCO2 = " << g.xco2 << endl;
    out << "xO2 = " << g.xo2 << endl;
    out << "xSOOT = " << g.fsoot << endl;
  }

};

#endif // _RADIATING_GAS_INCLUDED
