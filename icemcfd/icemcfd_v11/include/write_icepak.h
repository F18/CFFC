
/* Interface for writing out result files in ICEM CFD format.
 *
 * Copyright (c) 1999 ICEM CFD Engineering.
 * Author:  Wayne A. Christopher
 *          ICEM CFD Engineering, Berkeley, California
 *          wayne@icemcfd.com
 */

#ifndef WRITE_ICEPAK_H
#define WRITE_ICEPAK_H

extern int df_write_icepak(char* filename, int n_nodes,
			   int timestep, double timeval, double timeincr,
			   int n_faces3, int* face3info,
			   int n_faces4, int* face4info,
			   int n_fresults, char** fresult_names,
			   double** fresult_vals,
			   int n_nresults, char** nresult_names,
			   double** nresult_vals,
			   int n_resids, int niters, char** resid_names,
			   double** resid_vals);

#endif
