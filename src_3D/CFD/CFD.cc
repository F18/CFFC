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

    cout << setprecision(6);
    if (First_Step) {
        cout << "  n = " << Number_of_Time_Steps
             << "   t = " << Time
             << "   CPU t = " << CPU_Time.min()
             << "   L2-norm = ";
        cout.setf(ios::scientific);
        cout << Residual_L2_Norm << "\n  .";
        cout.unsetf(ios::scientific);
        cout.flush();
    } else if (Number_of_Time_Steps-Frequency*
               (Number_of_Time_Steps/Frequency) == 0 ) {
        cout << "\n" 
             << "  n = " << Number_of_Time_Steps
             << "   t = " << Time 
             << "   CPU t = " << CPU_Time.min()
             << "   L2-norm = ";
        cout.setf(ios::scientific); 
        cout << Residual_L2_Norm << "\n  .";
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
			    const int Frequency,
			    const int progress_character) {

  cout << setprecision(6);
  if (First_Step) {
    cout << "  n = " << Number_of_Time_Steps
	 << "   t = " << Time
	 << "   CPU t = " << CPU_Time.min()
	 << "   L2-norm = ";
    cout.setf(ios::scientific);
    cout << Residual_L2_Norm << "\n  .";
    cout.unsetf(ios::scientific);
    cout.flush();
  } else if (Number_of_Time_Steps-Frequency*
	     (Number_of_Time_Steps/Frequency) == 0) {
    cout << "\n" 
	 << "  n = " << Number_of_Time_Steps
	 << "   t = " << Time 
	 << "   CPU t = " << CPU_Time.min()
	 << "   L2-norm = ";
    cout.setf(ios::scientific); 
    cout << Residual_L2_Norm << "\n  .";
    cout.unsetf(ios::scientific);
    cout.flush();
  }
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

    int i;
    char prefix[256], prefix2[256], extension[256], 
         progress_file_name[256], gnuplot_file_name[256];
    char *progress_file_name_ptr, *gnuplot_file_name_ptr;
    ofstream gnuplot_file;

    /* Determine the name of the progress file. */

    i = 0;
    while (1) {
       if (File_Name[i] == ' ' ||
           File_Name[i] == '.') break;
       prefix[i] = File_Name[i];
       prefix2[i] = File_Name[i];
       i = i + 1;
       if (i > strlen(File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    prefix2[i] = '\0';
    strcat(prefix, "_residual");
    strcat(prefix2, "_view_residual");

    strcpy(extension, ".dat");
    strcpy(progress_file_name, prefix);
    strcat(progress_file_name, extension);

    progress_file_name_ptr = progress_file_name;

    /* Open the progress file. */

    if (Append_to_File) {
       Progress_File.open(progress_file_name_ptr, ios::out|ios::app);
    } else {
       Progress_File.open(progress_file_name_ptr, ios::out);
    } /* endif */
    if (Progress_File.bad()) return (1);

    /* Write the appropriate GNUPLOT command file for 
       plotting progress file information. */

    strcpy(extension, ".gplt");
    strcpy(gnuplot_file_name, prefix);
    strcat(gnuplot_file_name, extension);

    gnuplot_file_name_ptr = gnuplot_file_name;

    gnuplot_file.open(gnuplot_file_name_ptr, ios::out);
    if (gnuplot_file.bad()) return(1);

    gnuplot_file << "set title \"Solution Convergence\"\n"
                 << "set xlabel \"N (iterations)\"\n"
                 << "set ylabel \"residual\"\n" 
                 << "set grid\n"
                 << "set logscale y\n"
                 << "plot \"" << progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%lf%*lf%*lf\" \\\n"
                 << "     title \"L1-norm\" with lines, \\\n"
                 << "\"" << progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%lf%*lf\" \\\n"
	         << "     title \"L2-norm\" with lines, \\\n"
                 << "\"" << progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%*lf%lf\" \\\n"
	         << "     title \"Max-norm\" with lines\n"
                 << "pause -1  \"Hit return to continue\"\n";

    gnuplot_file.close();

    strcpy(gnuplot_file_name, prefix2);
    strcat(gnuplot_file_name, extension);

    gnuplot_file_name_ptr = gnuplot_file_name;

    gnuplot_file.open(gnuplot_file_name_ptr, ios::out);
    if (gnuplot_file.bad()) return(1);

    gnuplot_file << "set title \"Solution Convergence\"\n"
                 << "set xlabel \"N (iterations)\"\n"
                 << "set ylabel \"residual\"\n" 
                 << "set grid\n"
                 << "set logscale y\n"
                 << "plot \"" << progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%lf%*lf%*lf\" \\\n"
                 << "     title \"L1-norm\" with lines, \\\n"
                 << "\"" << progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%lf%*lf\" \\\n"
	         << "     title \"L2-norm\" with lines, \\\n"
                 << "\"" << progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%*lf%lf\" \\\n"
	         << "     title \"Max-norm\" with lines\n"
                 << "pause 1\n"
                 << "reread\n";

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

    Progress_File << setprecision(6);
    Progress_File << Number_of_Time_Steps
                  << " " << Time
                  << " " << CPU_Time.min();
    Progress_File.setf(ios::scientific);
    Progress_File << " " << Residual_L1_Norm
                  << " " << Residual_L2_Norm 
                  << " " << Residual_Max_Norm
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


