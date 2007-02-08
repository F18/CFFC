/* TecplotAuxData.h:  Header file defining the data that are included in 
                      the tecplot files as auxiliary data. */

#ifndef _TECPLOT_AUX_DATA_INCLUDED
#define _TECPLOT_AUX_DATA_INCLUDED

/* Include header files. */
#include "CFD/CFD.h"
#include "SourceRevisionData.h"

class TecplotDatasetAuxData: SourceCode{

 public:

  static void PrintAuxData(std::ofstream &output_file);

 protected:
  TecplotDatasetAuxData(){};
};

/***************************************************
 * PrintAuxData()                                  *
 * Print the tecplot auxiliary dataset             *
 ************************************************/
inline void TecplotDatasetAuxData::PrintAuxData(std::ofstream &output_file){

  /* Output SourceCode data */
  SourceCode::PrintTecplotAuxData(output_file);

  /* Output Simulation data */
  output_file << "DATASETAUXDATA Simulation_Date = \"  " << Date_And_Time() << "   \" \n ";
}

#endif
