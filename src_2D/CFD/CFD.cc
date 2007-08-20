/* CFD.cc:  Basic CFD Subroutines. */

/* Include the CFD header file. */

#ifndef _CFD_INCLUDED
#include "CFD.h"
#endif // _CFD_INCLUDED

/********************************************************
 * CFD -- I/O Routines.                                 *
 ********************************************************/

/********************************************************
 * Routine: Output_Progress                             *
 *                                                      *
 * This routine writes out progress information for a   *
 * CFD calculation, including iteration level, time,    *
 * CPU time, and residual norms to the standard output  *
 * device.                                              *
 *                                                      *
 ********************************************************/
void Output_Progress(const int Number_of_Time_Steps,
                     const double &Time,
                     const CPUTime &CPU_Time,
                     const double &Residual_L1_Norm,
                     const int First_Step,
                     const int Frequency) {

    cout << setprecision(6);
    if (First_Step) {
        cout << "  n = " << Number_of_Time_Steps
             << "   t = " << Time
             << "   CPU t = " << CPU_Time.min()
             << "   L1-norm = ";
        cout.setf(ios::scientific);
        cout << Residual_L1_Norm << "\n  .";
        cout.unsetf(ios::scientific);
        cout.flush();
    } else if (Number_of_Time_Steps-Frequency*
               (Number_of_Time_Steps/Frequency) == 0 ) {
        cout << "\n" 
             << "  n = " << Number_of_Time_Steps
             << "   t = " << Time 
             << "   CPU t = " << CPU_Time.min()
             << "   L1-norm = ";
        cout.setf(ios::scientific); 
        cout << Residual_L1_Norm << "\n  .";
        cout.unsetf(ios::scientific);
        cout.flush();
//    } else if (Number_of_Time_Steps-(Frequency/2)*
//               (Number_of_Time_Steps/(Frequency/2)) == 0 ) {
//        cout << "\n  .";
    } else {
        cout << "."; 
        cout.flush();
    } /* endif */

}

/********************************************************
 * Routine: Output_Progress_L2norm                      *
 *                                                      *
 * This routine writes out progress information for a   *
 * CFD calculation, including iteration level, time,    *
 * CPU time, and residual norms to the standard output  *
 * device.                                              *
 *                                                      *
 ********************************************************/
void Output_Progress_L2norm(const int Number_of_Time_Steps,
			    const double &Time,
			    const CPUTime &CPU_Time,
			    const double &Residual_L2_Norm,
			    const int First_Step,
			    const int Frequency) {

   static const int progress_character = 0;
   static const double Ratio_Residual_L2_Norm = -1.0;

   Output_Progress_L2norm(Number_of_Time_Steps,
       			  Time,
			  CPU_Time,
			  Residual_L2_Norm,
			  Ratio_Residual_L2_Norm,
			  First_Step,
			  Frequency,
			  progress_character);
}

void Output_Progress_L2norm(const int Number_of_Time_Steps,
			    const double &Time,
			    const CPUTime &CPU_Time,
			    const double &Residual_L2_Norm,
			    const int First_Step,
			    const int Frequency,
			    const int progress_character) {

   static const double Ratio_Residual_L2_Norm = -1.0;

   Output_Progress_L2norm(Number_of_Time_Steps,
			  Time,
	 		  CPU_Time,
			  Residual_L2_Norm,
			  Ratio_Residual_L2_Norm,
			  First_Step,
			  Frequency,
			  progress_character);
}

void Output_Progress_L2norm(const int Number_of_Time_Steps,
			    const double &Time,
			    const CPUTime &CPU_Time,
			    const double &Residual_L2_Norm,
			    const double &Ratio_Residual_L2_Norm,
			    const int First_Step,
			    const int Frequency) {
   static const int progress_character = 0;
   Output_Progress_L2norm(Number_of_Time_Steps,
	  		  Time,
			  CPU_Time,
			  Residual_L2_Norm,
			  Ratio_Residual_L2_Norm,
			  First_Step,
			  Frequency,
			  progress_character);
}

void Output_Progress_L2norm(const int Number_of_Time_Steps,
			    const double &Time,
			    const CPUTime &CPU_Time,
			    const double &Residual_L2_Norm,
			    const double &Ratio_Residual_L2_Norm,
			    const int First_Step,
			    const int Frequency,
			    const int progress_character) {

    cout << setprecision(6);
    if (First_Step || Number_of_Time_Steps % Frequency == 0) {
      cout << endl
	   << "  n = " << Number_of_Time_Steps
	   << "   t = " << Time 
	   << "   CPU t = " << CPU_Time.min();
      cout.setf(ios::scientific);
      cout << "   L2-norm = " << Residual_L2_Norm;
      if (Ratio_Residual_L2_Norm > 0) {
         cout << "   L2-norm ratio = " << Ratio_Residual_L2_Norm;
      } /* endif */
      cout.unsetf(ios::scientific);
      cout << endl << "  ";
    } /* endif */

    switch(progress_character) {
       case 0: 
         cout << "."; 
         break;
       case 1: 
         cout << "+"; 
         break;
       case 2: 
         cout << "-"; 
         break;
       case 3: 
         cout << "o"; 
         break;
       case 4: 
         cout << ":"; 
         break;
       case 5: 
         cout << "~"; 
         break;
       default: 
         cout << "."; 
         break;
    };
    cout.flush();

}

/********************************************************
 * Routine: Open_Progress_File                          *
 *                                                      *
 * This routine opens the progress file for a CFD       *
 * calculation and prepares it for I/O.                 *
 *                                                      *
 ********************************************************/
int Open_Progress_File(ofstream &Progress_File,
                       char *File_Name,
                       const int Append_to_File) {

    char prefix[256], 
         progress_file_name[256], 
         gnuplot_file_name[256],
	 eps_output_file_name[256];

    ofstream gnuplot_file;

    /* Determine the name of the progress file. */

    strcpy(prefix, File_Name);
    for (int i = 0; i < strlen(prefix); i++) {
       if (prefix[i] == ' ' || prefix[i] == '.') {
	   prefix[i] = '\0'; break;
       } /* endif */
    } /* endfor */

    strcpy(progress_file_name, prefix);
    strcat(progress_file_name, "_residual.dat");

    /* Open the progress file. */

    if (Append_to_File) {
       Progress_File.open(progress_file_name, ios::out|ios::app);
    } else {
       Progress_File.open(progress_file_name, ios::out);
       int out_precision = 5;
       int out_width = out_precision + 8;
       int swidth = 6;
       Progress_File << left;
       Progress_File << "# | Number of Time Steps\n";
       Progress_File << setw(swidth) << "# |" << "  | Solution Time\n";
       Progress_File << setw(swidth) << "# |" << setw(out_width) << "  |" << " ";
       Progress_File << "  | Total CPU Time (minutes) (from all CPUs)\n";
       Progress_File << setw(swidth) << "# |" << setw(out_width) << "  |" << " " << setw(out_width) << "  |" << " ";
       Progress_File << "  | Residual L1-Norm (search on residual_variable in the code for more info)\n";
       Progress_File << setw(swidth) << "# |" << setw(out_width) << "  |" << " " << setw(out_width) << "  |" << " "
	             << setw(out_width) << "  |" << " ";
       Progress_File << "  | Residual L2-Norm\n";
       Progress_File << setw(swidth) << "# |" << setw(out_width) << "  |" << " " << setw(out_width) << "  |" << " "
	             << setw(out_width) << "  |" << " " << setw(out_width) << "  |" << " ";
       Progress_File << "  | Residual Max-Norm\n";
       Progress_File << right;
       Progress_File.flush();
    } /* endif */

    if (!Progress_File.good()) return (1);

    /* Write the appropriate GNUPLOT command file for 
       plotting progress file information. */

    strcpy(gnuplot_file_name, prefix);
    strcat(gnuplot_file_name, "_residual.gplt");

    gnuplot_file.open(gnuplot_file_name, ios::out);
    if (!gnuplot_file.good()) return(1);

    strcpy(eps_output_file_name, prefix);
    strcat(eps_output_file_name, "_plot_residual_iter.eps");

    gnuplot_file << "set title \"Solution Convergence\"\n"
                 << "set xlabel \"N (iterations)\"\n"
                 << "set ylabel \"Residual\"\n" 
                 << "set grid\n"
                 << "set logscale y\n"
                 << "plot \"" << progress_file_name << "\""
                 << " using 1:2 \"%lf%*lf%*lf%lf%*lf%*lf\" \\\n"
                 << "     title \"L1-norm\" with lines, \\\n"
                 << "\"" << progress_file_name << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%lf%*lf\" \\\n"
	         << "     title \"L2-norm\" with lines, \\\n"
                 << "\"" << progress_file_name << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%*lf%lf\" \\\n"
	         << "     title \"Max-norm\" with lines\n"
                 << "set terminal postscript eps enhanced color\n"
                 << "set output \"" << eps_output_file_name << "\"\n"
                 << "replot\n"
                 << "pause -1  \"Press return to continue\"\n";

    gnuplot_file.close();

    strcpy(gnuplot_file_name, prefix);
    strcat(gnuplot_file_name, "_view_residual.gplt");
    gnuplot_file.open(gnuplot_file_name, ios::out);
    if (!gnuplot_file.good()) return(1);

    gnuplot_file << "set title \"Solution Convergence\"\n"
                 << "set xlabel \"N (iterations)\"\n"
                 << "set ylabel \"Residual\"\n" 
                 << "set grid\n"
                 << "set logscale y\n"
                 << "plot \"" << progress_file_name << "\""
                 << " using 1:2 \"%lf%*lf%*lf%lf%*lf%*lf\" \\\n"
                 << "     title \"L1-norm\" with lines, \\\n"
                 << "\"" << progress_file_name << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%lf%*lf\" \\\n"
	         << "     title \"L2-norm\" with lines, \\\n"
                 << "\"" << progress_file_name << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%*lf%lf\" \\\n"
	         << "     title \"Max-norm\" with lines\n"
                 << "pause 1\n"
                 << "reread\n";

    gnuplot_file.close();

    strcpy(gnuplot_file_name, prefix);
    strcat(gnuplot_file_name, "_residual_cpu.gplt");

    gnuplot_file.open(gnuplot_file_name, ios::out);
    if (!gnuplot_file.good()) return(1);

    strcpy(eps_output_file_name, prefix);
    strcat(eps_output_file_name, "_plot_residual_cpu.eps");

    gnuplot_file << "set title \"Solution Convergence\"\n"
                 << "set xlabel \"Total CPU time (minutes)\"\n"
                 << "set ylabel \"Residual\"\n" 
                 << "set grid\n"
                 << "set logscale y\n"
                 << "plot \"" << progress_file_name << "\""
                 << " using 1:2 \"%*lf%*lf%lf%lf%*lf%*lf\" \\\n"
                 << "     title \"L1-norm\" with linespoints, \\\n"
                 << "\"" << progress_file_name << "\""
                 << " using 1:2 \"%*lf%*lf%lf%*lf%lf%*lf\" \\\n"
	         << "     title \"L2-norm\" with linespoints, \\\n"
                 << "\"" << progress_file_name << "\""
                 << " using 1:2 \"%*lf%*lf%lf%*lf%*lf%lf\" \\\n"
	         << "     title \"Max-norm\" with linespoints\n"
                 << "set terminal postscript eps enhanced color\n"
                 << "set output \"" << eps_output_file_name << "\"\n"
                 << "replot\n"
                 << "pause -1  \"Press return to continue\"\n";

    gnuplot_file.close();

    /* Preparation of progress file complete.
       Return zero value. */

    return(0);

}

/**********************************************************
 * Routine: Output_Progress_to_File                       *
 *                                                        *
 * This routine writes out progress information for       *
 * a CFD calculation to a progress file, including        *
 * iteration level, time, CPU time, and residual norm(s). *
 *                                                        *
 **********************************************************/
void Output_Progress_to_File(ostream &Progress_File,
                             const int Number_of_Time_Steps,
                             const double &Time,
                             const CPUTime &CPU_Time,
                             const double *Residual_L1_Norm,
                             const double *Residual_L2_Norm,
                             const double *Residual_Max_Norm,
			     const int &Residual_Norm,
			     const int &Number_of_Residual_Norms) {

    Progress_File << setprecision(6);
    Progress_File << Number_of_Time_Steps
                  << " " << Time
                  << " " << CPU_Time.min();
    Progress_File.setf(ios::scientific);

    for(int q=0; q < Number_of_Residual_Norms; q++) {
	Progress_File << " " << Residual_L1_Norm[q]
		      << " " << Residual_L2_Norm[q] 
		      << " " << Residual_Max_Norm[q];
    } /* endfor */

    Progress_File << "\n";
    Progress_File.unsetf(ios::scientific);
    Progress_File.flush();

}

/********************************************************
 * Routine: Open_Progress_File                          *
 *                                                      *
 * This routine opens the progress file for a CFD       *
 * calculation and prepares it for I/O.                 *
 *                                                      *
 ********************************************************/
 //CURRENTLY SET TO SHOW FIRST 4 SOLN VARS NORMS
int Open_Progress_File(ofstream &Progress_File,
                       char *File_Name,
                       const int Append_to_File,
		       const int &Residual_Norm,
		       const int &Number_of_Residual_Norms) {

    char prefix[256], extension[256], 
         progress_file_name[256], 
         gnuplot_file_name[256];

    char *progress_file_name_ptr, 
         *gnuplot_file_name_ptr;

    ofstream gnuplot_file;

    /* Determine the name of the progress file. */
    int i = 0;
    while (1) {
       if (File_Name[i] == ' ' ||
           File_Name[i] == '.') break;
       prefix[i] = File_Name[i];
       i = i + 1;
       if (i > strlen(File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_residual");

    strcpy(extension, ".dat");
    strcpy(progress_file_name, prefix);
    strcat(progress_file_name, extension);

    progress_file_name_ptr = progress_file_name;

    /* Open the progress file. */
    if (Append_to_File) {
       Progress_File.open(progress_file_name_ptr, ios::out|ios::app);
    } else {
       Progress_File.open(progress_file_name_ptr, ios::out);
    } 
    if (Progress_File.bad()) return (1);

    /* Write the appropriate GNUPLOT command file for 
       plotting progress file information. */
    strcpy(extension, ".gplt");
    strcpy(gnuplot_file_name, prefix);
    strcat(gnuplot_file_name, extension);

    gnuplot_file_name_ptr = gnuplot_file_name;

    gnuplot_file.open(gnuplot_file_name_ptr, ios::out);
    if (gnuplot_file.bad()) return(1);

    // Write gnuplot script
    gnuplot_file << "set term x11 enhanced font \",10\"\n"
//                  << "set xlabel \"N (iterations)\"\n"
//                  << "set ylabel \"residual\"\n" 
                 << "set logscale y \n"
		 << "set style data lines \n"
		 << "set key spacing 0.5 \n"
                 << "set size 1.0,1.0 \n"
                 << "set origin 0.0,0.0 \n"
		 << "set title 0,-1 \n"
                 << "set multiplot \n";
    
    // Format  time, cput, {L1, L2, Max}xNumber_of_Residual_Norms for 4
    for(int q=4; q < Number_of_Residual_Norms*3+3; q=q+3){
      gnuplot_file << "set size 0.5,0.5 \n";
      switch(q){
      case 4:
	gnuplot_file << "set origin 0.0,0.5 \n"
		     << "set title \"Density convergence\"\n";
	break;
      case 7:
	gnuplot_file << "set origin 0.5,0.5 \n"
		     << "set title \"Momentum-x convergence\"\n";
	break;
      case 10:
	gnuplot_file << "set origin 0.0,0.0 \n"
		     << "set title \"Momentum-y convergence\"\n";
	break;
      case 13:
	gnuplot_file << "set origin 0.5,0.0 \n"
		     << "set title \"Energy convergence\"\n";
	break;
      default:
	gnuplot_file << "ISSUES IN Open_Progress_File \n";
	break;
      }            
      gnuplot_file<<"plot \""<< progress_file_name_ptr << "\" using 1:"<<q<<" title \"L1-norm\", \\\n"
		  <<"\""<< progress_file_name_ptr <<"\" using 1:"<<q+1<<" title \"L2-norm\", \\\n"
		  <<"\""<< progress_file_name_ptr <<"\" using 1:"<<q+2<<" title \"Max-norm\" \n";
    }
    gnuplot_file << "unset multiplot \n";
    gnuplot_file<< "pause -1  \"Hit return to continue\"\n";

    gnuplot_file.close();

    /* Preparation of progress file complete.
       Return zero value. */

    return(0);

}
/********************************************************
 * Routine: Close_Progress_File                         *
 *                                                      *
 * This routine closes the progress file for a CFD      *
 * calculation.                                         *
 *                                                      *
 ********************************************************/
int Close_Progress_File(ofstream &Progress_File) {

    Progress_File.close();

    return(0);

}

/********************************************************
 * Routine: Output_Progress_to_File                     *
 *                                                      *
 * This routine writes out progress information for     *
 * a CFD calculation to a progress file, including      *
 * iteration level, time, CPU time, and residual norms. *
 *                                                      *
 ********************************************************/
void Output_Progress_to_File(ostream &Progress_File,
                             const int Number_of_Time_Steps,
                             const double &Time,
                             const CPUTime &CPU_Time,
                             const double &Residual_L1_Norm,
                             const double &Residual_L2_Norm,
                             const double &Residual_Max_Norm) {

    int out_precision = 5;
    int out_width = out_precision + 8;
    int swidth = 6;
    Progress_File << setprecision(out_precision);
    Progress_File << setw(swidth) << Number_of_Time_Steps
                  << " " << setw(out_width) << Time
                  << " " << setw(out_width) << CPU_Time.min();
    Progress_File.setf(ios::scientific);
    Progress_File << " " << setw(out_width) << Residual_L1_Norm
                  << " " << setw(out_width) << Residual_L2_Norm 
                  << " " << setw(out_width) << Residual_Max_Norm
                  << "\n";
    Progress_File.unsetf(ios::scientific);
    Progress_File.flush();

}

/********************************************************
 * CFD -- Algebraic Grid Point Distribution Stretching  *
 *        Functions.                                    *
 ********************************************************/

/********************************************************
 * Routine: StretchingFcn (Grid Point Distribution      *
 *                         Stretching Function)         *
 *                                                      *
 * This function returns the value of the various       *
 * algebraic stretching functions used to control       *
 * grid point distributions in CFD mesh.  See Anderson, *
 * Tannehill, and Pletcher (1984).                      *
 *                                                      *
 ********************************************************/
double StretchingFcn(const double &xx,
                     const double &beta,
                     const double &tau,
                     const int i_stretch) {
    double a, zz;
  
    switch (i_stretch) {
      case STRETCHING_FCN_LINEAR:
        zz = xx;
        break;
      case STRETCHING_FCN_MIN_CLUSTERING:
        zz = pow((beta+ONE) / (beta-ONE), ONE - xx);
        zz = ((beta+ONE) - (beta-ONE) * zz) /
             (zz + ONE);
        break;
      case STRETCHING_FCN_MAX_CLUSTERING:
        a = ZERO;
        zz = pow((beta+ONE) / (beta-ONE), (xx - a) /
                 (ONE - a));
        zz = (zz * (beta + TWO * a) - beta + TWO * a) /
             ((ONE + TWO * a) * (zz + ONE));
        break;
      case STRETCHING_FCN_MINMAX_CLUSTERING:
        a = HALF;
        zz = pow((beta+ONE) / (beta-ONE), (xx - a) /
                 (ONE - a));
        zz = (zz * (beta + TWO * a) - beta + TWO * a) /
             ((ONE + TWO * a) * (zz + ONE));
        break;
      case STRETCHING_FCN_MIDPT_CLUSTERING:
        a = HALF * log ( (ONE + (exp(tau) - ONE) * beta) /
            (ONE + (exp(-tau) - ONE) * beta) ) / tau;
        zz = beta * (ONE + sinh(tau * (xx - a)) / sinh(tau * a));
        break;
      case STRETCHING_FCN_SINE:
        zz = sin(HALF*PI*xx);
        break;
      case STRETCHING_FCN_COSINE:
        zz = ONE - cos(HALF*PI*xx);
        break;
      default:
        zz = xx;
        break;
    } /* endswitch */

    return (zz);

}

/********************************************************
 * CFD -- Entropy Fixes.                                *
 ********************************************************/

/********************************************************
 * Routine: HartenFixPos (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the positive part of a         *
 * corrected elemental wave speed or eigenvalue         *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
double HartenFixPos(const double &lambda_a,
                    const double &lambda_l,
                    const double &lambda_r) {
    double delta;

    delta = max(ZERO, TWO*(lambda_r-lambda_l));
    if (fabs(lambda_a) > delta || delta < TOLER) {
        delta = ZERO;
    } else {
        delta = HALF*(sqr(lambda_a)/delta+delta)
                - fabs(lambda_a);
    } /* endif */

    return (HALF*(lambda_a+fabs(lambda_a)+delta));

}

/********************************************************
 * Routine: HartenFixNeg (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the negative part of a         *
 * corrected elemental wave speed or eigenvalue         *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
double HartenFixNeg(const double &lambda_a,
                    const double &lambda_l,
                    const double &lambda_r) {
    double delta;

    delta = max(ZERO, TWO*(lambda_r-lambda_l));
    if (fabs(lambda_a) > delta || delta < TOLER) {
        delta = ZERO;
    } else {
        delta = HALF*(sqr(lambda_a)/delta+delta)
                -fabs(lambda_a);
    } /* endif */
    
    return (HALF*(lambda_a-fabs(lambda_a)-delta));

}

/********************************************************
 * Routine: HartenFixAbs (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the absolute value of a        *
 * corrected elemental wave speed or eigenvalue         *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
double HartenFixAbs(const double &lambda_a,
                    const double &lambda_l,
                    const double &lambda_r) {
    double delta;

    delta = max(ZERO, TWO*(lambda_r-lambda_l));
    if (fabs(lambda_a) > delta || delta < TOLER) {
        delta = ZERO;
    } else {
        delta = HALF*(sqr(lambda_a)/delta+delta)
                -fabs(lambda_a);
    } /* endif */

    return (fabs(lambda_a)+delta);

}

/********************************************************
 * CFD -- Slope Limiters.                               *
 ********************************************************/

/********************************************************
 * Routine: Limiter_BarthJespersen (Slope Limiter of    *
 *          Barth and Jespersen)                        *
 *                                                      *
 * This function returns the value of the slope or      *
 * gradient limiter according to the formulation        *
 * proposed by Barth and Jespersen (1989).  Given the   *
 * minimum and maximum values of all cell centered      *
 * values used in the reconstruction and the unlimited  *
 * values of the solution at the flux quadrature        *
 * points, the routine returns the computed value of    *
 * the limiter.                                         *
 *                                                      *
 ********************************************************/
double Limiter_BarthJespersen(double *uQuad,
                              const double &u0,
                              const double &u0Min,
	      	              const double &u0Max,
			      const int nQuad) {
    int i;
    double phi, y;

    phi = ONE;

    for (i = 0; i <= nQuad-1; i++) {
       if (uQuad[i] - u0 > TOLER) {
           phi = min(phi, min(ONE,(u0Max-u0)/(uQuad[i]-u0)));
       } else if (uQuad[i] - u0 < -TOLER) {
           phi = min(phi, min(ONE,(u0Min-u0)/(uQuad[i]-u0)));
       } else {
           phi = min(phi, ONE);
       } /* endif */
    } /* endfor */

    return(phi);

}

/********************************************************
 * Routine: Limiter_Venkatakrishnan (Slope Limiter of   *
 *          Venkatakrishnan)                            *
 *                                                      *
 * This function returns the value of the slope or      *
 * gradient limiter according to the formulation        *
 * proposed by Venkatakrishnan (1993).  Given the       *
 * minimum and maximum values of all cell centered      *
 * values used in the reconstruction and the unlimited  *
 * values of the solution at the flux quadrature        *
 * points, the routine returns the computed value of    *
 * the limiter.                                         *
 *                                                      *
 ********************************************************/
double Limiter_Venkatakrishnan(double *uQuad,
                               const double &u0,
                               const double &u0Min,
	      	               const double &u0Max,
			       const int nQuad) {
    int i;
    double phi, y;

    phi = ONE;

    for (i = 0; i <= nQuad-1; i++) {
       if (uQuad[i] - u0 > TOLER) {
           y = (u0Max-u0)/(uQuad[i]-u0);
           phi = min(phi, (sqr(y)+TWO*y)/(sqr(y)+y+TWO));
       } else if (uQuad[i] - u0 < -TOLER) {
           y = (u0Min-u0)/(uQuad[i]-u0);
           phi = min(phi, (sqr(y)+TWO*y)/(sqr(y)+y+TWO));
       } else {
           phi = min(phi, ONE);
       } /* endif */
    } /* endfor */

    return(phi);

}

/********************************************************
 * Routine: Limiter_VanLeer (Slope Limiter of Van Leer) *
 *                                                      *
 * This function returns the value of the slope or      *
 * gradient limiter according to the formulation        *
 * proposed by Van Leer (1978).  Given the              *
 * minimum and maximum values of all cell centered      *
 * values used in the reconstruction and the unlimited  *
 * values of the solution at the flux quadrature        *
 * points, the routine returns the computed value of    *
 * the limiter.                                         *
 *                                                      *
 ********************************************************/
double Limiter_VanLeer(double *uQuad,
                       const double &u0,
                       const double &u0Min,
	      	       const double &u0Max,
		       const int nQuad) {
    int i;
    double phi, y;

    phi = ONE;

    for (i = 0; i <= nQuad-1; i++) {
       if (uQuad[i] - u0 > TOLER) {
           y = (u0Max-u0)/(uQuad[i]-u0);
           phi = min(phi, HALF*vanleer(ONE, y));
       } else if (uQuad[i] - u0 < -TOLER) {
           y = (u0Min-u0)/(uQuad[i]-u0);
           phi = min(phi, HALF*vanleer(ONE, y));
       } else {
           phi = min(phi, ONE);
       } /* endif */
    } /* endfor */ 

    return(phi);

}

/********************************************************
 * Routine: Limiter_VanAlbada (Slope Limiter of         *
 *          Van Albada)                                 *
 *                                                      *
 * This function returns the value of the slope or      *
 * gradient limiter according to the formulation        *
 * proposed by Van Albada (1982).  Given the            *
 * minimum and maximum values of all cell centered      *
 * values used in the reconstruction and the unlimited  *
 * values of the solution at the flux quadrature        *
 * points, the routine returns the computed value of    *
 * the limiter.                                         *
 *                                                      *
 ********************************************************/
double Limiter_VanAlbada(double *uQuad,
                         const double &u0,
                         const double &u0Min,
	      	         const double &u0Max,
			 const int nQuad) {
    int i;
    double phi, y;

    phi = ONE;

    for (i = 0; i <= nQuad-1; i++) {
       if (uQuad[i] - u0 > TOLER) {
           y = (u0Max-u0)/(uQuad[i]-u0);
           phi = min(phi, HALF*vanalbada(ONE, y, 0.10));
       } else if (uQuad[i] - u0 < -TOLER) {
           y = (u0Min-u0)/(uQuad[i]-u0);
           phi = min(phi, HALF*vanalbada(ONE, y, 0.10));
       } else {
           phi = min(phi, ONE);
       } /* endif */
    } /* endfor */

    return(phi);

}

/********************************************************
 * CFD -- Coefficients of Runge-Kutta and Multistage    *
 *        Explicit Time Stepping Schemes.               *
 ********************************************************/

/********************************************************
 * Routine: Runge_Kutta(Coefficients for Runge-Kutta    *
 *          Multistage Explicit Time Stepping Schemes)  *
 *                                                      *
 * This function returns the coefficient for the        *
 * I_stage of the general N_stage Runge-Kutta scheme.   *
 *                                                      *
 ********************************************************/
double Runge_Kutta(const int I_stage,
                   const int N_stage) {

    double beta;
   
    switch(N_stage) {
      case 1:
        beta = ONE;
        break;
      case 2 :
        beta = ONE;
	if (I_stage == 2) beta = HALF*beta;
        break;
      case 4 :
        beta = ONE;
	if (I_stage == 1) beta = HALF*beta;
	if (I_stage == 2) beta = HALF*beta;
	if (I_stage == 4) beta = beta/SIX;
        break;
      case 5 :
// 	beta = ONE;
// 	if (I_stage == 1) beta *= QUARTER;
// 	if (I_stage == 2) beta /= SIX;
// 	if (I_stage == 3) beta *= THREE/EIGHT;
// 	if (I_stage == 4) beta *= HALF;
	beta = ONE;
	if (I_stage == 1) beta *= 0.059;
	if (I_stage == 2) beta *= 0.14;
	if (I_stage == 3) beta *= 0.273;
	if (I_stage == 4) beta *= HALF;
        break;
      default:
        beta = ONE;
        break;
    } /* endswitch */
      
    return (beta);

}

/********************************************************
 * Routine: MultiStage_Optimally_Smoothing(Coefficients *
 *          for Optimally-Smoothing Multistage Explicit *
 *          Time Stepping Schemes of Van Leer, Tai, and *
 *          Powell)                                     *
 *                                                      *
 * This function returns the coefficient for the        *
 * I_stage of the general N_stage optimally-smoothing   *
 * scheme developed by Van Leer, Tai, and Powell        *
 * (1989).                                              *
 *                                                      *
 ********************************************************/
double MultiStage_Optimally_Smoothing(const int I_stage,
                                      const int N_stage,
                                      const int Limiter_Type) {

    double beta;
   
    switch(N_stage) {
      case 1:
        beta = ONE;
        break;
      case 2 :
	if (Limiter_Type == LIMITER_ZERO) {
          beta = ONE;
	  if (I_stage == 1) beta = 0.3333*beta;
	} else {
          beta = 0.4693;
	  if (I_stage == 1) beta = 0.4242*beta;
	} /* endif */
        break;
      case 3 :
	if (Limiter_Type == LIMITER_ZERO) {
          beta = 1.50;
	  if (I_stage == 1) beta = 0.1481*beta;
	  if (I_stage == 2) beta = 0.4000*beta;
	} else {
          beta = 0.6936;
	  if (I_stage == 1) beta = 0.1918*beta;
	  if (I_stage == 2) beta = 0.4929*beta;
	} /* endif */
        break;
      case 4 :
	if (Limiter_Type == LIMITER_ZERO) {
          beta = 2.00;
	  if (I_stage == 1) beta = 0.0833*beta;
	  if (I_stage == 2) beta = 0.2069*beta;
	  if (I_stage == 3) beta = 0.4265*beta;
	} else {
          beta = 0.9214;
	  if (I_stage == 1) beta = 0.1084*beta;
	  if (I_stage == 2) beta = 0.2602*beta;
	  if (I_stage == 3) beta = 0.5052*beta;
	} /* endif */
        break;
      case 5 :
	if (Limiter_Type == LIMITER_ZERO) {
          beta = 2.50;
	  if (I_stage == 1) beta = 0.0533*beta;
	  if (I_stage == 2) beta = 0.1263*beta;
	  if (I_stage == 3) beta = 0.2375*beta;
	  if (I_stage == 4) beta = 0.4414*beta;
	} else {
          beta = 1.1508;
	  if (I_stage == 1) beta = 0.0695*beta;
	  if (I_stage == 2) beta = 0.1602*beta;
	  if (I_stage == 3) beta = 0.2898*beta;
	  if (I_stage == 4) beta = 0.5060*beta;
	} /* endif */
        break;
      case 6 :
	if (Limiter_Type == LIMITER_ZERO) {
          beta = 3.00;
	  if (I_stage == 1) beta = 0.0370*beta;
	  if (I_stage == 2) beta = 0.0851*beta;
	  if (I_stage == 3) beta = 0.1521*beta;
	  if (I_stage == 4) beta = 0.2562*beta;
	  if (I_stage == 5) beta = 0.4512*beta;
	} else {
          beta = 1.3805;
	  if (I_stage == 1) beta = 0.0482*beta;
	  if (I_stage == 2) beta = 0.1085*beta;
	  if (I_stage == 3) beta = 0.1885*beta;
	  if (I_stage == 4) beta = 0.3050*beta;
	  if (I_stage == 5) beta = 0.5063*beta;
	} /* endif */
        break;
      default:
        beta = ONE;
        break;
    } /* endswitch */
      
    return (beta);

}

/**********************************************************************
 * Routine: A_Stable_Implicit_Method_Coefficients                     *
 *                                                                    *
 * This routine assigns the coefficients theta, xi, and phi for a     *
 * specified unconditionally stable (A-stable) implicit method.  See  *
 * the text book by Lomax, Pulliam, and Zingg (2001) pages 128-130    *
 * for details.                                                       *
 *                                                                    *
 **********************************************************************/
extern void A_Stable_Implicit_Method_Coefficients(double &theta,
						  double &xi,
						  double &phi,
						  const int Time_Integration_Scheme) {

  switch(Time_Integration_Scheme) {
  case TIME_STEPPING_IMPLICIT_EULER :
    theta = 1.000;
    xi = 0.000;
    phi = 0.000;
    break;
  case TIME_STEPPING_IMPLICIT_TRAPEZOIDAL :
    theta = 0.500;
    xi = 0.000;
    phi = 0.000;
    break;
  case TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD :
    theta = 1.000;
    xi = 0.500;
    phi = 0.000;
    break;
  case TIME_STEPPING_IMPLICIT_ADAMS_TYPE :
    theta = 0.750;
    xi = 0.000;
    phi = -0.250;
    break;
  case TIME_STEPPING_IMPLICIT_LEES :
    theta = 0.333333333333;
    xi = -0.500;
    phi = -0.333333333333;
    break;
  case TIME_STEPPING_IMPLICIT_TWO_STEP_TRAPEZOIDAL :
    theta = 0.500;
    xi = -0.500;
    phi = -0.500;
    break;
  case TIME_STEPPING_IMPLICIT_A_CONTRACTIVE :
    theta = 0.625;
    xi = -0.166666666667;
    phi = -0.222222222222;
    break;
  default :
    // Second-order backward.
    theta = 1.000;
    xi = 0.500;
    phi = 0.000;
  };

}

/*************************************************************
 * CFD1D_Input_Parameters -- External subroutines.           *
 *************************************************************/

/********************************************************
 * Routine: Open_Input_File                             *
 *                                                      *
 * Opens the appropriate input data file.               *
 *                                                      *
 ********************************************************/
void Open_Input_File(CFD1D_Input_Parameters &IP) {

    IP.Input_File.open(IP.Input_File_Name, ios::in);
    if (! IP.Input_File.bad()) {
       IP.Line_Number = 0;
       IP.Input_File.setf(ios::skipws);
    } /* endif */

}

/********************************************************
 * Routine: Close_Input_File                            *
 *                                                      *
 * Closes the appropriate input data file.              *
 *                                                      *
 ********************************************************/
void Close_Input_File(CFD1D_Input_Parameters &IP) {

    IP.Input_File.unsetf(ios::skipws);
    IP.Input_File.close();

}

/********************************************************
 * Routine: Set_Default_Input_Parameters                *
 *                                                      *
 * Assigns default values to the input parameters.      *
 *                                                      *
 ********************************************************/
void Set_Default_Input_Parameters(CFD1D_Input_Parameters &IP) {

    char *string_ptr;

    string_ptr = "CFD1D.in";
    strcpy(IP.Input_File_Name, string_ptr);

    string_ptr = "Explicit_Euler";
    strcpy(IP.Time_Integration_Type, string_ptr);
    IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
    IP.Time_Accurate = 1;
    IP.Local_Time_Stepping = 0;
    IP.Maximum_Number_of_Time_Steps = 100;
    IP.N_Stage = 1;
    IP.Time_Max = ZERO;

    string_ptr = "MUSCL";
    strcpy(IP.Reconstruction_Type, string_ptr);
    IP.i_Reconstruction = RECONSTRUCTION_MUSCL;

    string_ptr = "VanLeer";
    strcpy(IP.Limiter_Type, string_ptr);
    IP.i_Limiter = LIMITER_VANLEER;

    string_ptr = "Roe";
    strcpy(IP.Flux_Function_Type, string_ptr);
    IP.i_Flux_Function = FLUX_FUNCTION_ROE;

    string_ptr = "Uniform";
    strcpy(IP.ICs_Type, string_ptr);
    IP.i_ICs = IC_UNIFORM;

    IP.Kappa = 0.02;
    IP.a = ONE;
    IP.Tau = ONE;

    string_ptr = "Uniform";
    strcpy(IP.Grid_Type, string_ptr);
    IP.i_Grid = GRID_CARTESIAN_UNIFORM;
    IP.X_Min = -ONE;
    IP.X_Max = ONE;
    IP.Number_of_Cells = 50;
    IP.Number_of_Nodes = 51;

    string_ptr = "outputfile.dat";
    strcpy(IP.Output_File_Name, string_ptr);

    string_ptr = "gnuplotfile.gplt";
    strcpy(IP.Gnuplot_File_Name, string_ptr);

    string_ptr = "Gnuplot";
    strcpy(IP.Output_Format_Type, string_ptr);
    IP.i_Output_Format = IO_GNUPLOT;

    string_ptr = " ";
    strcpy(IP.Next_Control_Parameter, string_ptr);

    IP.Line_Number = 0;

}

/********************************************************
 * Routine: Get_Next_Input_Control_Parameter            *
 *                                                      *
 * Get the next input control parameter from the input  *
 * file.                                                *
 *                                                      *
 ********************************************************/
void Get_Next_Input_Control_Parameter(CFD1D_Input_Parameters &IP) {

    int i;
    char buffer[256];

    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File.getline(buffer, sizeof(buffer));
    i = 0;
    if (buffer[0] != '#') {
       while (1) {
          if (buffer[i] == ' ' || buffer[i] == '=') break;
          i = i + 1;
          if (i > strlen(buffer) ) break;
       } /* endwhile */
       buffer[i] = '\0';
    } /* endif */
    strcpy(IP.Next_Control_Parameter, buffer);

}

/********************************************************
 * Routine: Parse_Next_Input_Control_Parameter          *
 *                                                      *
 * Parses and executes the next input control parameter *
 * from the input file.                                 *
 *                                                      *
 ********************************************************/
int Parse_Next_Input_Control_Parameter(CFD1D_Input_Parameters &IP) {

    int i_command;
    char buffer[256];

    i_command = 0;

    if (strcmp(IP.Next_Control_Parameter, "Time_Integration_Type") == 0) {
       i_command = 1;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Time_Integration_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Time_Integration_Type, "Explicit_Euler") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Implicit_Euler") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Implicit_Trapezoidal") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_IMPLICIT_TRAPEZOIDAL;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Explicit_Predictor_Corrector") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
           IP.N_Stage = 2;
       } else if (strcmp(IP.Time_Integration_Type, "Semi_Implicit_Euler") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_SEMI_IMPLICIT_EULER;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Semi_Implicit_Trapezoidal") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_SEMI_IMPLICIT_TRAPEZOIDAL;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Semi_Implicit_Predictor_Corrector") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_SEMI_IMPLICIT_PREDICTOR_CORRECTOR;
           IP.N_Stage = 2;
       } else if (strcmp(IP.Time_Integration_Type, "Multistage_Optimal_Smoothing") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
           IP.N_Stage = 4;
       } else if (strcmp(IP.Time_Integration_Type, "Lax_Friedrichs") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_LAX_FRIEDRICHS;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Lax_Wendroff") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_LAX_WENDROFF;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "MacCormack") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_MACCORMACK;
           IP.N_Stage = 2;
       } else if (strcmp(IP.Time_Integration_Type, "Hancock") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_HANCOCK;
           IP.N_Stage = 2;
       } else if (strcmp(IP.Time_Integration_Type, "Simple_Explicit") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_SIMPLE_EXPLICIT;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Simple_Implicit") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_SIMPLE_IMPLICIT;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Crank_Nicolson") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_CRANK_NICOLSON;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "ADE") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_ADE;
           IP.N_Stage = 2;
       } else {
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
           IP.N_Stage = 1;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Reconstruction_Type") == 0) {
       i_command = 2;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Reconstruction_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Reconstruction_Type, "MUSCL") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_MUSCL;
       } else if (strcmp(IP.Reconstruction_Type, "Green_Gauss") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
       } else if (strcmp(IP.Reconstruction_Type, "Least_Squares") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
       } else if (strcmp(IP.Reconstruction_Type, "Characteristic") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_CHARACTERISTIC;
       } else {
          IP.i_Reconstruction = RECONSTRUCTION_MUSCL;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Limiter_Type") == 0) {
       i_command = 3;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Limiter_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Limiter_Type, "One") == 0) {
          IP.i_Limiter = LIMITER_ONE;
       } else if (strcmp(IP.Limiter_Type, "Zero") == 0) {
          IP.i_Limiter = LIMITER_ZERO;
       } else if (strcmp(IP.Limiter_Type, "Minmod") == 0) {
          IP.i_Limiter = LIMITER_MINMOD;
       } else if (strcmp(IP.Limiter_Type, "UMIST") == 0) {
          IP.i_Limiter = LIMITER_UMIST;
       } else if (strcmp(IP.Limiter_Type, "Double_Minmod") == 0) {
          IP.i_Limiter = LIMITER_DOUBLE_MINMOD;
       } else if (strcmp(IP.Limiter_Type, "Superbee") == 0) {
          IP.i_Limiter = LIMITER_SUPERBEE;
       } else if (strcmp(IP.Limiter_Type, "Phi") == 0) {
          IP.i_Limiter = LIMITER_PHI;
       } else if (strcmp(IP.Limiter_Type, "VanLeer") == 0) {
          IP.i_Limiter = LIMITER_VANLEER;
       } else if (strcmp(IP.Limiter_Type, "VanAlbada") == 0) {
          IP.i_Limiter = LIMITER_VANALBADA;
       } else if (strcmp(IP.Limiter_Type, "Barth_Jespersen") == 0) {
          IP.i_Limiter = LIMITER_BARTH_JESPERSEN;
       } else if (strcmp(IP.Limiter_Type, "Venkatakrishnan") == 0) {
          IP.i_Limiter = LIMITER_VENKATAKRISHNAN;
       } else {
          IP.i_Limiter = LIMITER_VANLEER ;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Flux_Function_Type") == 0) {
       i_command = 4;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Flux_Function_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Flux_Function_Type, "Godunov") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV;
       } else if (strcmp(IP.Flux_Function_Type, "Roe") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_ROE;
       } else if (strcmp(IP.Flux_Function_Type, "Rusanov") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_RUSANOV;
       } else if (strcmp(IP.Flux_Function_Type, "HLLE") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_HLLE;
       } else if (strcmp(IP.Flux_Function_Type, "Linde") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_LINDE;
       } else if (strcmp(IP.Flux_Function_Type, "HLLC") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_HLLC;
       } else if (strcmp(IP.Flux_Function_Type, "Osher") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_OSHER;
       } else {
          IP.i_Flux_Function = FLUX_FUNCTION_ROE;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "ICs_Type") == 0) {
       i_command = 5;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICs_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.ICs_Type, "Constant") == 0) {
          IP.i_ICs = IC_CONSTANT;
       } else if (strcmp(IP.ICs_Type, "Uniform") == 0) {
          IP.i_ICs = IC_UNIFORM;
       } else if (strcmp(IP.ICs_Type, "Sod") == 0) {
          IP.i_ICs = IC_SOD;
       } else if (strcmp(IP.ICs_Type, "Groth") == 0) {
          IP.i_ICs = IC_GROTH;
       } else if (strcmp(IP.ICs_Type, "Einfeldt") == 0) {
          IP.i_ICs = IC_EINFELDT;
       } else if (strcmp(IP.ICs_Type, "Shock_Wave") == 0) {
          IP.i_ICs = IC_SHOCK_WAVE;
       } else if (strcmp(IP.ICs_Type, "Contact_Surface") == 0) {
          IP.i_ICs = IC_CONTACT_SURFACE;
       } else if (strcmp(IP.ICs_Type, "Rarefaction_Wave") == 0) {
          IP.i_ICs = IC_RAREFACTION_WAVE;
       } else if (strcmp(IP.ICs_Type, "Brio_Wu") == 0) {
          IP.i_ICs = IC_BRIO_WU;
       } else if (strcmp(IP.ICs_Type, "Riemann") == 0) {
          IP.i_ICs = IC_RIEMANN;
       } else if (strcmp(IP.ICs_Type, "Square_Wave") == 0) {
          IP.i_ICs = IC_SQUARE_WAVE;
       } else if (strcmp(IP.ICs_Type, "Sinx2_Wave") == 0) {
          IP.i_ICs = IC_SINX2_WAVE;
       } else if (strcmp(IP.ICs_Type, "Impulsive_Rod") == 0) {
          IP.i_ICs = IC_IMPULSIVE_ROD;
          IP.X_Min = ZERO;
          IP.X_Max = TEN;
       } else if (strcmp(IP.ICs_Type, "Sinusoidal_Rod1") == 0) {
          IP.i_ICs = IC_SINUSOIDAL_ROD1;
          IP.X_Min = ZERO;
          IP.X_Max = TEN;
       } else if (strcmp(IP.ICs_Type, "Sinusoidal_Rod4") == 0) {
          IP.i_ICs = IC_SINUSOIDAL_ROD4;
          IP.X_Min = ZERO;
          IP.X_Max = TEN;
       } else if (strcmp(IP.ICs_Type, "Riemann_IVP_qx=0") == 0) {
          IP.i_ICs = IC_RIEMANN_IVP_QX0;
          IP.X_Min = -ONE;
          IP.X_Max = ONE;
       } else if (strcmp(IP.ICs_Type, "Riemann_IVP_T=0") == 0) {
          IP.i_ICs = IC_RIEMANN_IVP_T0;
          IP.X_Min = -ONE;
          IP.X_Max = ONE;
       } else if (strcmp(IP.ICs_Type, "Riemann_IVP") == 0) {
          IP.i_ICs = IC_RIEMANN_IVP;
          IP.X_Min = -ONE;
          IP.X_Max = ONE;
       } else {
          IP.i_ICs = IC_UNIFORM;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Grid_Type") == 0) {
       i_command = 6;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Grid_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Grid_Type, "Uniform") == 0) {
          IP.i_Grid = GRID_CARTESIAN_UNIFORM;
       } else {
          IP.i_Grid = GRID_CARTESIAN_UNIFORM;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Output_File_Name") == 0) {
       i_command = 7;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Output_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Output_File_Name, ".dat");
       strcpy(IP.Gnuplot_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Gnuplot_File_Name, ".gplt");

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells") == 0) {
       i_command = 8;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Cells;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Cells < 1) i_command = INVALID_INPUT_VALUE;
       IP.Number_of_Nodes = IP.Number_of_Cells + 1;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Nodes") == 0) {
       i_command = 9;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Nodes;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Nodes < 2) i_command = INVALID_INPUT_VALUE;
       IP.Number_of_Cells = IP.Number_of_Nodes - 1;

    } else if (strcmp(IP.Next_Control_Parameter, "Time_Accurate") == 0) {
       i_command = 10;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Accurate;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Time_Accurate != 0 &&
           IP.Time_Accurate != 1) IP.Time_Accurate = 0;
       if (IP.Time_Accurate) {
          IP.Local_Time_Stepping = 0;
       } else {
          IP.Local_Time_Stepping = 1;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Local_Time_Stepping") == 0) {
       i_command = 11;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Local_Time_Stepping;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Local_Time_Stepping != 0 &&
           IP.Local_Time_Stepping != 1) IP.Local_Time_Stepping = 1;

    } else if (strcmp(IP.Next_Control_Parameter, "Maximum_Number_of_Time_Steps") == 0) {
       i_command = 12;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Maximum_Number_of_Time_Steps;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "N_Stage") == 0) {
       i_command = 13;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.N_Stage;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.N_Stage < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "CFL_Number") == 0) {
       i_command = 14;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.CFL_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "X_Min") == 0) {
       i_command = 15;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Min;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.X_Min <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "X_Max") == 0) {
       i_command = 16;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Max;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.X_Max <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Time_Max") == 0) {
       i_command = 17;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Max;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Kappa") == 0) {
       i_command = 18;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Kappa;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Kappa < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "a") == 0) {
       i_command = 19;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.a;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Tau") == 0) {
       i_command = 20;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Tau;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Tau < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Output_Format_Type") == 0) {
       i_command = 21;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Output_Format_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Output_Format_Type, "Gnuplot") == 0  ||
           strcmp(IP.Output_Format_Type, "gnuplot") == 0  ||
           strcmp(IP.Output_Format_Type, "GNUPLOT") == 0) {
          IP.i_Output_Format = IO_GNUPLOT;
       } else if (strcmp(IP.Output_Format_Type, "Tecplot") == 0  ||
                  strcmp(IP.Output_Format_Type, "tecplot") == 0  ||
                  strcmp(IP.Output_Format_Type, "TECPLOT") == 0) {
          IP.i_Output_Format = IO_TECPLOT;
       } else if (strcmp(IP.Output_Format_Type, "Matlab") == 0  ||
                  strcmp(IP.Output_Format_Type, "matlab") == 0  ||
                  strcmp(IP.Output_Format_Type, "MATLAB") == 0) {
          IP.i_Output_Format = IO_MATLAB;
       } else if (strcmp(IP.Output_Format_Type, "Octave") == 0  ||
                  strcmp(IP.Output_Format_Type, "octave") == 0  ||
                  strcmp(IP.Output_Format_Type, "OCTAVE") == 0) {
          IP.i_Output_Format = IO_OCTAVE;
       } else {
          IP.i_Output_Format = IO_TECPLOT;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Execute") == 0) {
       i_command = EXECUTE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Terminate") == 0) {
       i_command = TERMINATE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Continue") == 0) {
       i_command = CONTINUE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output") == 0) {
       i_command = WRITE_OUTPUT_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Restart") == 0) {
       i_command = WRITE_RESTART_CODE;

    } else if (IP.Next_Control_Parameter[0] == '#') {
       i_command = COMMENT_CODE;

    } else {
       i_command = INVALID_INPUT_CODE;

    } /* endif */

    /* Return the parser command type indicator. */

    return (i_command);

}

// for use with qsort
int Compare_Ints(const void *p, const void *q) {
    return *(int *)p - *(int *)q;
}

double Least_Squares_Slope(double *values, int n, int position) {

    double ssxx = 0.0, ssxy = 0.0;
    double ymean = 0.0, xmean = (n-1.0)/2.0;

    for (int x = 0; x < n; x++) { 
	ymean += values[x];
	ssxx += x * x;
	ssxy += x * values[ (position+x) % n ];
    }

    ymean *= 1.0 / n;

    ssxx -= n * xmean * xmean;
    ssxy -= n * xmean * ymean;

    return ssxy / ssxx;
}

