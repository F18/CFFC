#ifndef _TURBULENT_VELOCITY_FIELD_INCLUDED 
#define _TURBULENT_VELOCITY_FIELD_INCLUDED

/* Include required C++ header files. */

#include <cmath> 
#include <cassert>
#include <cstdlib>     // defines the drand48() function
#include <ctime>       // defines the time() function
#include <limits>
#include <complex>

/* Include required CFFC header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif // _VECTOR3D_INCLUDED

#ifndef _TENSOR3D_INCLUDED
#include "../Math/Tensor3D.h"
#endif // _TENSOR3D_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // _INPUT_INCLUDED

#ifndef _GRID3D_HEXA_MULTIBLOCK_INCLUDED
#include "../Grid/Grid3DHexaMultiBlock.h"
#endif // _GRID3D_HEXA_MULTIBLOCK_INCLUDED

#ifndef _ADAPTIVEBLOCK3D_INCLUDED
#include "../AMR/AdaptiveBlock3D.h"
#endif // _ADAPTIVEBLOCK3D_INCLUDED

#ifndef _FFTW_INCLUDED
#include "fftw3.h"
#endif // _FFTW_INCLUDED

// Constants
const complex<double>  I(0.0, 1.0);      // sqrt(-1.0)

#define TURBULENT_VELOCITY_FIELD_DATA_USED        1

#define TURBULENT_VELOCITY_FIELD_DATA_NOT_USED    0

/*!
 * Class: Turbulent_Velocity_Field_Block
 *
 * \brief Class used to store turbulent velocity field for a single block.
 *
 */
class Turbulent_Velocity_Field_Block {
  public:
    int           NCi,ICl,ICu; // i-direction turbulent velocity field cell counters
    int           NCj,JCl,JCu; // j-direction turbulent velocity field cell counters
    int           NCk,KCl,KCu; // k-direction turbulent velocity field cell counters
    int                Nghost; // number of ghost cells
    Vector3D      ***Velocity; // array of turbulent velocity field vectors

    int Allocated; // Indicates whether or not the turbulent velocity field data has been allocated.
    
    /* Constructors. */
    Turbulent_Velocity_Field_Block(void) {
       NCi = 0; ICl = 0; ICu = 0; 
       NCj = 0; JCl = 0; JCu = 0;
       NCk = 0; KCl = 0; KCu = 0;
       Nghost = 0;
       Allocated = TURBULENT_VELOCITY_FIELD_DATA_NOT_USED;
       Velocity = NULL; 
    }

    Turbulent_Velocity_Field_Block(const int Ni, 
                                   const int Nj, 
                                   const int Nk,
                                   const int Ng) {
      allocate(Ni, Nj, Nk, Ng);
    }

    /* Destructor. */
    ~Turbulent_Velocity_Field_Block(void) {
       deallocate();
    }

    /* Allocate memory for velocity field data. */
    void allocate(const int Ni, 
                  const int Nj, 
                  const int Nk,
                  const int Ng);
    
    /* Deallocate memory for velocity field data. */
    void deallocate(void);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, 
                                 const Turbulent_Velocity_Field_Block &V);

    friend istream &operator >> (istream &in_file, 
                                 Turbulent_Velocity_Field_Block &V);
  
    /* Other useful member functions. */

  private:
    //copy and assignment are not permitted
    Turbulent_Velocity_Field_Block(const Turbulent_Velocity_Field_Block &V);
    Turbulent_Velocity_Field_Block &operator =(const Turbulent_Velocity_Field_Block &V);
};

/*************************************************************************
 * Turbulent_Velocity_Field_Block::allocate -- Allocate memory.          *
 *************************************************************************/
inline void Turbulent_Velocity_Field_Block::allocate(const int Ni, 
                                                     const int Nj, 
                                                     const int Nk,
                                                     const int Ng) {
   assert( Ni >= 1 && Nj >= 1 && Nk >= 1 && Ng >=1 && !Allocated);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
   NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1;
   NCk = Nk+2*Ng; KCl = Ng; KCu = Nk+Ng-1;
   Nghost = Ng;
   Allocated = TURBULENT_VELOCITY_FIELD_DATA_USED;

   Velocity = new Vector3D**[NCi];
   for (int i = 0; i <= NCi-1; ++i ){
      Velocity[i] = new Vector3D*[NCj];
      for (int j = 0; j <= NCj-1; ++j ){
         Velocity[i][j] = new Vector3D[NCk];
      } /* endfor */
   } /* endfor */
}

/*************************************************************************
 * Turbulent_Velocity_Field_Block::deallocate -- Deallocate memory.      *
 *************************************************************************/
inline void Turbulent_Velocity_Field_Block::deallocate(void) {
   if (Allocated) {
      assert(NCi >= 1 && NCj >= 1 && NCk >= 1);
      for (int i = 0; i <= NCi-1 ; ++i ) {
         for ( int j = 0 ; j <= NCj-1 ; ++j) {
            delete []Velocity[i][j]; Velocity[i][j] = NULL;
         } /* endfor */
         delete []Velocity[i]; Velocity[i] = NULL;
      }/*endfor*/
      delete []Velocity; Velocity = NULL;
  
      NCi = 0; ICl = 0; ICu = 0; 
      NCj = 0; JCl = 0; JCu = 0; 
      NCk = 0; KCl = 0; KCu = 0; 
      Nghost = 0;
      Allocated = TURBULENT_VELOCITY_FIELD_DATA_NOT_USED;
   } /* endif */
}

/*!
 * Class: Turbulent_Velocity_Field_Multi_Block
 *
 * \brief Class used to store turbulent velocity field for a
 *        1D array of turbulent velocity field solution blocks.
 *
 */
class Turbulent_Velocity_Field_Multi_Block_List {
  public:
    Turbulent_Velocity_Field_Block  *Vel_Blks; // one dimensional array of grid block.
    int                                  NBlk;
    int                             NBlk_Idir, 
                                    NBlk_Jdir, 
                                    NBlk_Kdir; // Number of blocks in i, j and k directions.
    int                             Allocated; // Indicates if the grid blocks have been allocated or not.
 
    /* Creation constructors. */
    Turbulent_Velocity_Field_Multi_Block_List(void) : 
       NBlk(0), NBlk_Idir(0), NBlk_Jdir(0), NBlk_Kdir(0), Vel_Blks(NULL), Allocated(0) { }

    Turbulent_Velocity_Field_Multi_Block_List(const int N) {
       Allocate(N);
    }

    Turbulent_Velocity_Field_Multi_Block_List(const int Ni, 
                                              const int Nj, 
                                              const int Nk) {
       Allocate(Ni, Nj, Nk);
    }

    /* Destructor. */
    ~Turbulent_Velocity_Field_Multi_Block_List(void) {
        Deallocate();
    }

    /* Other member functions  */

    void Allocate(const int Ni, const int Nj, const int Nk);

    void Allocate(const int N);

    void Deallocate(void);

    void Create(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                const Grid3D_Input_Parameters &Input);

  private:
    //copy and assignment are not permitted
    Turbulent_Velocity_Field_Multi_Block_List(const Turbulent_Velocity_Field_Multi_Block_List &V);
    Turbulent_Velocity_Field_Multi_Block_List &operator = (const Turbulent_Velocity_Field_Multi_Block_List &V);
};

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Allocate -- Allocate memory.   *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Allocate(const int Ni, 
                                                                const int Nj, 
                                                                const int Nk) {
   if (Ni >= 1 && Nj >= 1 && Nk >= 1 && !Allocated) {
      NBlk_Idir = Ni; 
      NBlk_Jdir = Nj; 
      NBlk_Kdir = Nk; 
      NBlk = Ni*Nj*Nk;
      Vel_Blks = new Turbulent_Velocity_Field_Block[NBlk];
      Allocated = 1;
   } /* endif */
}

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Allocate -- Allocate memory.   *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Allocate(const int N) {
   if (N >= 1 && !Allocated) {
      NBlk_Idir = N; 
      NBlk_Jdir = 1; 
      NBlk_Kdir = 1; 
      NBlk = N;
      Vel_Blks = new Turbulent_Velocity_Field_Block[NBlk];
      Allocated = 1;
   } /* endif */
}

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Allocate -- Deallocate memory.   *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Deallocate(void) {
   if (NBlk >= 1 && Allocated) {
       delete []Vel_Blks;
       Vel_Blks = NULL;
       NBlk_Idir = 0; 
       NBlk_Jdir = 0; 
       NBlk_Kdir = 0;
       NBlk = 0;
       Allocated = 0;
   } /* endif */
}

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Create --                      *
 *       Create memory containers for velocity field data.                   *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Create(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                                              const Grid3D_Input_Parameters &Input) {
   if (Initial_Mesh.NBlk >= 1 && Initial_Mesh.Allocated) {
      Allocate(Initial_Mesh.NBlk_Idir, Initial_Mesh.NBlk_Jdir, Initial_Mesh.NBlk_Kdir);
      for (int nBlk = 0; nBlk <= NBlk-1; ++nBlk ) {
	if (Initial_Mesh.Grid_Blks[nBlk].Allocated) Vel_Blks[nBlk].allocate(Initial_Mesh.Grid_Blks[nBlk].NCi-
                                                                            2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                            Initial_Mesh.Grid_Blks[nBlk].NCj-
                                                                            2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                            Initial_Mesh.Grid_Blks[nBlk].NCk-
                                                                            2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                            Initial_Mesh.Grid_Blks[nBlk].Nghost);
      } /* endfor */
   } /* endif */
}

/*!
 * Class: RandomFieldRogallo
 *
 * \brief Class defined to generate random fluctuations using
 * Rogallo's procedure.
 *
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
class RandomFieldRogallo {
  private:
    int spectrum_flag;  //!< Turbulence kinetic energy spectrum flag.

  public:

    //! Default constructor.
    RandomFieldRogallo() : spectrum_flag(SPECTRUM_VON_KARMAN_PAO) { }

    //! Another constructor.
    RandomFieldRogallo(const int &SPECTRUM) : spectrum_flag(SPECTRUM) { }  


    double k_1(const int &n1, const double &L1) const { return double(n1); }//2.0*PI*n1/L1; }//double(n1); }

    double k_2(const int &n2, const double &L2) const { return double(n2); }//2.0*PI*n2/L2; }//double(n2); }

    double k_3(const int &n3, const double &L3) const { return double(n3); }//2.0*PI*n3/L3; }//double(n3); }

    double random_double() const { return  drand48(); }

    // alpha and beta as defined by Rogallo, 1981
    complex<double> alpha(const double &abs_wave_num, 
			  const double &theta1, 
			  const double &phi) const;

    complex<double> beta(const double &abs_wave_num, 
			 const double &theta2, 
			 const double &phi) const;

    double Energy_Spectrum_Value(const double &abs_wave_num) const;

    int Create_Homogeneous_Turbulence_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                                     Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
			                             const Grid3D_Input_Parameters &IPs);

    void Write_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                        const Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
					const Grid3D_Input_Parameters &IPs) const;
       
    void Read_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                       Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
				       Grid3D_Input_Parameters &IPs);

};

//-----------------------------------------------------//
//         Members of RandomFieldRogallo class         //
//-----------------------------------------------------//

// alpha
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
alpha(const double &abs_wave_num, const double &theta1, const double &phi) const {

  double E, k;
  k = abs_wave_num;  E = Energy_Spectrum_Value(k);

  return (k == 0.0)  ?  0.0 : sqrt(E/(4.0*PI*k*k)) * exp(I*theta1) * cos(phi);
}

// beta
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
beta(const double &abs_wave_num, const double &theta2, const double &phi) const {

  double E, k;
  k = abs_wave_num;  E = Energy_Spectrum_Value(k);

  return (k == 0.0)  ?  0.0 : sqrt(E/(4.0*PI*k*k)) * exp(I*theta2) * sin(phi);
}

// Prescribed energy spectrum
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
double RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Energy_Spectrum_Value(const double &abs_wave_num) const {

  double k, kp, kd, eps, u, Lp, EE;
  double C, A, alpha, a_s;
  int s;
 
  k = abs_wave_num;
  
  switch (spectrum_flag) {
    /*****  Lee and Reynolds  *****/
  case SPECTRUM_LEE_REYNOLDS :
    C = 20.0;  kp = 80.0;   
    EE = (k <= kp)  ?  C*k*k : C*kp*kp*pow(k/kp, -5.0/3.0);
    break;

    /*****  Laval and Nazarenko paper  *****/
  case SPECTRUM_LAVAL_NAZARENKO :
    //double EE, kp, C;
    C = 1.0;  kp = 4.0;
    EE = C*k*exp(-pow(k/kp, 2.0));
    break;

    /*****   von Karman-Pao   *****/
  case SPECTRUM_VON_KARMAN_PAO :
    A = 1.5;    alpha = 1.5;
    //Lp = TWO*PI/ONE;
    kp =1.0;   kd = 300.0 /*362.0*/;
    u = 0.8;    eps = 2.73E-3; //2.73E-03;
    EE = (A/eps)*pow(u, 5.0);
    EE = EE*pow(k/kp, 4.0) * exp(-3.0*alpha*pow(k/kd, 4.0/3.0)/2.0);
    EE = EE/pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
    break;

    /*****  Haworth and Poinsot paper  *****/
  case SPECTRUM_HAWORTH_POINSOT :
    //double EE, kp, 
    u = 4.5;  Lp = TWO*PI/4.0;  // kp=4.0, u=2.5  kp=8
    kp = TWO*PI/Lp;
    EE = (32.0/3.0) * sqrt(2.0/PI)* (u*u/kp) * pow(k/kp, 4.0) * exp(-2.0*(k/kp)*(k/kp));
    break;

    /*****  Chasnov paper 1996  *****/
  case SPECTRUM_CHASNOV :  
    //double a_s, EE, kp = 4.0 /*8.0 20.0 4.0*/, u = 0.095; /*0.1 28.3 21.21  0.001*/  
    //int s = 3;
    s = 3;
    kp = 4.0;
    u = 0.009; 
    // u = 0.001  ->  Re_lambda = 10 
    // u = 0.009  ->  Re_lambda = 98
    // u = 0.095  ->  Re_lambda = 950
   
    a_s = pow(2.0*double(s)+1.0, double(s)+1.0)/(factorial(s)*pow(2.0, double(s)));
    EE = (HALF*a_s*u*u/kp)*pow(k/kp, 2.0*double(s)+1.0);
    EE = EE*exp(-(double(s)+HALF)*pow(k/kp, 2.0));
    break;

    /*****   Bell & Day report   *****/
  case SPECTRUM_BELL_DAY :
    //  kd = 1/(2*dx)
    kp = 3.0;   kd = ONE/0.576E-3;
    EE = pow(k/kp, 4.0) * exp(-9.0*pow(k/kd, 4.0/3.0)/4.0);
    EE /= pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
    break;

    /*****   von Karman-Pao   *****/
  default :
    A = 1.5;    alpha = 1.5;
    //Lp = TWO*PI/ONE;
    kp = 3.0;   kd = 1000.0 /*362.0*/;
    u = 0.107;    eps = 2.73E-3; //2.73E-03;
    EE = (A/eps)*pow(u, 5.0);
    EE = EE*pow(k/kp, 4.0) * exp(-3.0*alpha*pow(k/kd, 4.0/3.0)/2.0);
    EE = EE/pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
    break;    

  }  // end switch

  return (k == 0.0)  ?  0.0 : EE;
 
}
 
// Create_Homogeneous_Turbulence_Velocity_Field
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Create_Homogeneous_Turbulence_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                             Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
                                             const Grid3D_Input_Parameters &IPs) {

  double L1, L2, L3;
  L1 = IPs.Box_Length;
  L2 = IPs.Box_Width;
  L3 = IPs.Box_Height;

  int Nx, Ny, Nz;
  Nx = IPs.NCells_Idir * Initial_Mesh.NBlk_Idir;
  Ny = IPs.NCells_Jdir * Initial_Mesh.NBlk_Jdir;
  Nz = IPs.NCells_Kdir * Initial_Mesh.NBlk_Kdir;

  double scaling_factor = 1.0/double(Nx*Ny*Nz);  // Scaling factor for the complex to real transform

  double        *u, *v, *w;        // Arrays to store the velocity fluctuations in physical space
  fftw_complex  *uu, *vv, *ww;     // Arrays to store the velocity fluctuations in Fourier space
  fftw_plan      physical;

  int index;
  int nz = Nz/2+1;

  // Allocation of arrays used in the transforms
  u = (double *) malloc(Nx*Ny*Nz * sizeof(double));
  v = (double *) malloc(Nx*Ny*Nz * sizeof(double));
  w = (double *) malloc(Nx*Ny*Nz * sizeof(double));

  uu = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
  vv = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
  ww = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
  
  int seed = 1; 
  //int seed = time(NULL);   // assigns the current time to the seed
  srand48(seed);             // changes the seed for drand48()
  int iconj, jconj;          // Position of the conjugate complex for the i index
  double k1, k2, k3;         // Wave numbers
  
  double theta1, theta2, phi;
  complex<double> aa, bb;
  double deno;

  if (CFFC_Primary_MPI_Processor()) {
    cout << "\n\n ==========================================================================\n"; 
    cout << " Generating Homegeneous Turbulent Velocity Field"<<endl;
    cout << " ==========================================================================" << endl;
  }

  for (int i=0; i<Nx; ++i) {
    iconj = (i==0  ?  0 : Nx-i);
    for (int j=0; j<Ny; ++j) {
      jconj = (j==0  ?  0 : Ny-j);
      for (int l=0; l<Nz/2+1; ++l) {
	//Components of the wave number vector
	if( i<=Nx/2) {
	  k1 = k_1(i, L1);
	} else {
	  k1 = k_1(i-Nx, L1);
	}

	if( j<=Ny/2) {
	  k2 = k_2(j, L2);
	} else {
	  k2 = k_2(j-Ny, L2);
	}

	k3 = k_3(l, L3);

	// Wave number magnitude
	double abs_k = sqrt(k1*k1 + k2*k2 + k3*k3);

	//Rogallo's function  
	theta1 = 2.0*PI*random_double();  // Random number (0, 2*PI)
	theta2 = 2.0*PI*random_double();  // Random number (0, 2*PI)
	phi = 2.0*PI*random_double();      // Random number (0, 2*PI)

	if ( theta1 == theta2  && theta2 == phi ) {
	  cerr << "\n theta1, theta2 and phi are all equal.";
	} /* endif */

	aa = alpha(abs_k, theta1, phi);
	bb = beta(abs_k, theta2, phi);

	deno = abs_k * sqrt(k1*k1 + k2*k2);

	index = l + nz*(j+Ny*i);

	if (deno != 0.0) {
	  uu[index][0] = real( (aa*abs_k*k2 + bb*k1*k3)/deno );
	  uu[index][1] = imag( (aa*abs_k*k2 + bb*k1*k3)/deno );
	  
	  vv[index][0] = real( (bb*k2*k3 - aa*abs_k*k1)/deno );
	  vv[index][1] = imag( (bb*k2*k3 - aa*abs_k*k1)/deno );

	  ww[index][0] = real( -( bb*(k1*k1 + k2*k2) )/deno );
	  ww[index][1] = imag( -( bb*(k1*k1 + k2*k2) )/deno );

	} else {
	  uu[index][0] = 0.0;
	  uu[index][1] = 0.0;

	  vv[index][0] = 0.0;
	  vv[index][1] = 0.0;

	  ww[index][0] = 0.0;
	  ww[index][1] = 0.0;
	} /* endif */

 	if ( l==0  ||  l==Nz/2) {
	  // complex conjugates
	  if ( j>Ny/2  ||  ( i>Nx/2  &&  (j==0 || j==Ny/2) ) ) {
	    int index_conj = l + nz*(jconj+Ny*iconj);

	    uu[index][0] =  uu[index_conj][0];  
	    uu[index][1] = -uu[index_conj][1];

	    vv[index][0] =  vv[index_conj][0];
	    vv[index][1] = -vv[index_conj][1];

	    ww[index][0] =  ww[index_conj][0];
	    ww[index][1] = -ww[index_conj][1];

	  // real values at 8 corners
	  } else if ( (i==0 || i==Nx/2)  &&  (j==0 || j==Ny/2) ) {
	    uu[index][1] = 0.0;
	    vv[index][1] = 0.0;
	    ww[index][1] = 0.0; 
	  }

	} // end if

      } /* end for */
    } /* end for */
  } /* end for */
  
  physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, uu, u, /*FFTW_BACKWARD,*/ FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);
 
  physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, vv, v, /*FFTW_BACKWARD,*/ FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);
  
  physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, ww, w, /*FFTW_BACKWARD,*/ FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);
  
  for (int i=0; i<Nx; ++i) {
    for (int j=0; j<Ny; ++j) {
      for (int l=0; l<Nz; ++l) {
	index = l + Nz*(j+Ny*i);
	u[index] =  u[index] /*scaling_factor*/;
	v[index] =  v[index] /*scaling_factor*/;
	w[index] =  w[index] /*scaling_factor*/;
      } /* endfor */
    } /* endfor */
  } /* endfor */

  // Assign turbulent velocity field
  int nBlk, ix, iy, iz;
  for (int kBlk = 0; kBlk <= Initial_Mesh.NBlk_Kdir-1; ++kBlk) {
     for (int jBlk = 0; jBlk <= Initial_Mesh.NBlk_Jdir-1; ++jBlk) {
        for (int iBlk = 0; iBlk <= Initial_Mesh.NBlk_Idir-1; ++iBlk) {
            nBlk = iBlk + 
                   jBlk*Initial_Mesh.NBlk_Idir + 
                   kBlk*Initial_Mesh.NBlk_Idir*Initial_Mesh.NBlk_Jdir;
            for (int i = Initial_Mesh.Grid_Blks[nBlk].ICl; i <= Initial_Mesh.Grid_Blks[nBlk].ICu; ++i) {
               for (int j = Initial_Mesh.Grid_Blks[nBlk].JCl; j <= Initial_Mesh.Grid_Blks[nBlk].JCu; ++j) {
                  for (int k = Initial_Mesh.Grid_Blks[nBlk].KCl; k <= Initial_Mesh.Grid_Blks[nBlk].KCu; ++k) {
		     ix = iBlk*IPs.NCells_Idir+(i-Initial_Mesh.Grid_Blks[nBlk].ICl);
		     iy = jBlk*IPs.NCells_Jdir+(j-Initial_Mesh.Grid_Blks[nBlk].JCl);
		     iz = kBlk*IPs.NCells_Kdir+(k-Initial_Mesh.Grid_Blks[nBlk].KCl);
	             index = iz + 
                             iy*Nz + 
                             ix*Ny*Nz;
                     Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k].x = u[index];
                     Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k].y = v[index];
                     Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k].z = w[index];
		  } /* endfor */
	       } /* endfor */
	    } /* endfor */
	} /* endfor */
     } /* endfor */
  } /* endfor */

  // Deallocations   
  fftw_free(u);
  fftw_free(v);
  fftw_free(w);
  fftw_free(uu);
  fftw_free(vv);
  fftw_free(ww);

  return (0);

}

// Write_Turbulent_Velocity_Field
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Write_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
			       const Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
                               const Grid3D_Input_Parameters &IPs) const {

    ofstream out_file;
    out_file.open("Initial_Turbulence_Fluctuations.dat", ios::out);
    if (out_file.fail()) {
      cerr<<"\n Error opening file: Initial_Turbulence_Fluctuations.dat to write" << endl;
      exit(1);
    } /* endif */
  
    out_file.setf(ios::scientific);

    for (int nBlk = 0; nBlk <= Initial_Mesh.NBlk-1; ++nBlk ) {
       for (int i = Initial_Mesh.Grid_Blks[nBlk].ICl; i <= Initial_Mesh.Grid_Blks[nBlk].ICu; ++i) {
          for (int j = Initial_Mesh.Grid_Blks[nBlk].JCl; j <= Initial_Mesh.Grid_Blks[nBlk].JCu; ++j) {
             for (int k = Initial_Mesh.Grid_Blks[nBlk].KCl; k <= Initial_Mesh.Grid_Blks[nBlk].KCu; ++k) {
	      out_file << setprecision(10)
		       << Initial_Mesh.Grid_Blks[nBlk].Cell[i][j][k].Xc << " " 
                       << Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k];
	     } /* endfor */
	  } /* endfor */
       } /* endfor */
       out_file.setf(ios::scientific);
       out_file << endl;
    } /* endfor */

    out_file.close();

}

// Read_Turbulent_Velocity_Field
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Read_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
			      Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
                              Grid3D_Input_Parameters &IPs) {

   ifstream InFile;

   // Open data file for reading
   InFile.open("Initial_Turbulence_Fluctuations.dat", ios::in); 
   // Check to see if successful
   if(InFile.fail()){ 
     cerr<<"\n Error opening file: Initial_Turbulence.dat to read" <<endl;
     exit(1); 
   } 

   bool interpolate_flag;
   double xx, yy, zz, dx, dy, dz, dd, uprime_x, uprime_y, uprime_z;

   InFile.setf(ios::skipws);
 
   InFile.unsetf(ios::skipws);
   InFile.close();

}

// Assign_Homogeneous_Turbulence_Velocity_Field
template<typename HEXA_BLOCK>
void Assign_Homogeneous_Turbulence_Velocity_Field(HEXA_BLOCK *Solution_Block,
                                                  const AdaptiveBlock3D_List &LocalSolnBlockList,
                                                  const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field) {

   /* Assign initial turbulent velocity field to each solution block. */

   for (int nBlk = 0 ; nBlk <= LocalSolnBlockList.Nblk-1 ; nBlk++) {
      if (LocalSolnBlockList.Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
	 if (Velocity_Field.Vel_Blks[LocalSolnBlockList.Block[nBlk].gblknum].Allocated) {
  	    Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Block[nBlk],
                                                         Velocity_Field.Vel_Blks[LocalSolnBlockList.Block[nBlk].gblknum]);
	 } /* endif */
      } /* endif */
   }  /* endfor */

}

template<typename HEXA_BLOCK>
void Assign_Homogeneous_Turbulence_Velocity_Field(HEXA_BLOCK &Solution_Block,
                                                  const Turbulent_Velocity_Field_Block &Velocity_Field) {

  /* Assign initial turbulent velocity field to solution block. */

  for (int i = Solution_Block.ICl ; i <= Solution_Block.ICu ; i++) {
     for (int j = Solution_Block.JCl ; j <= Solution_Block.JCu ; j++) {
        for (int k = Solution_Block.KCl ; k <= Solution_Block.KCu ; k++) {
           Solution_Block.W[i][j][k].v += Velocity_Field.Velocity[i][j][k];
 	   Solution_Block.U[i][j][k] = Solution_Block.W[i][j][k].U();
	} /* endfor */
     } /* endfor */
  } /* endfor */

}

template<typename HEXA_BLOCK>
void Time_Averaging_of_Velocity_Field(HEXA_BLOCK *Solution_Block,
                                      AdaptiveBlock3D_List &LocalSolnBlockList,
                                      Grid3D_Input_Parameters &IPs,
                                      double &u_average,
                                      double &v_average,
                                      double &w_average) {

  double local_vol, Volume = ZERO;
  double Yfuel_conditional = ZERO;
  Vector3D vel;  vel.zero();
  //Conditional average on fresh gas
  Yfuel_conditional = 0.95*0.05518;//IPs.Fresh_Fuel_Mass_Fraction;
  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
          if (Solution_Block[p].W[i][j][k].spec[0].c >= Yfuel_conditional) {
	    local_vol = Solution_Block[p].Grid.volume(i,j,k);
	    vel += Solution_Block[p].W[i][j][k].v * local_vol;
            Volume += Solution_Block[p].Grid.volume(i,j,k);//Total_Block_Volume(Solution_Block[p]);
	  }
        }
      }
     }
   }
 }
  Volume = CFFC_Summation_MPI(Volume);
  vel.x = CFFC_Summation_MPI(vel.x);
  vel.y = CFFC_Summation_MPI(vel.y);
  vel.z = CFFC_Summation_MPI(vel.z);
  u_average = vel.x/Volume;
  v_average = vel.y/Volume;
  w_average = vel.z/Volume;
}

template <typename HEXA_BLOCK>
double Time_Averaging_of_Turbulent_Burning_Rate(HEXA_BLOCK *Solution_Block,
                                                AdaptiveBlock3D_List &LocalSolnBlockList,
                                                Grid3D_Input_Parameters &IPs){

  double local_vol, Yf_u, rho_u, Ly, Lz, burning_rate = ZERO;
  Yf_u = 0.05518;//IPs.Fresh_Fuel_Mass_Fraction;
  rho_u = 1.13;//IPs.Fresh_Density;
  Ly = IPs.Box_Width;
  Lz = IPs.Box_Height;
  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	    local_vol = Solution_Block[p].Grid.volume(i,j,k);
	    burning_rate +=  Solution_Block[p].W[i][j][k].Fsd*local_vol*Solution_Block[p].W[i][j][k].rho;
        }
      }
    }
  }
}
  burning_rate = CFFC_Summation_MPI(burning_rate);
  burning_rate = burning_rate*0.3837/(Ly*Lz);//IPs.laminar_flame_speed/Ly;//(rho_u*Ly);  //(rho_u*Yf_u*Ly);
  if (CFFC_Primary_MPI_Processor()) {
    cout << "\n\n ==========================================================================\n"; 
    cout << " Turbulent Burning Rate = " << burning_rate<<endl;
    cout << " ==========================================================================" << endl;
  } 
  return burning_rate;
}

template<typename HEXA_BLOCK>
void Time_Averaging_of_Solution(HEXA_BLOCK *Solution_Block,
                                AdaptiveBlock3D_List &LocalSolnBlockList,
                                Grid3D_Input_Parameters &IPs,
                                const double &u_average, 
                                const double &v_average,
                                const double &w_average,
                                double &sqr_u) {

  double vis, u_ave, v_ave, w_ave, local_vol, total_vol = ZERO, vis_ave = ZERO;
  double u_p=ZERO, v_p=ZERO, w_p=ZERO, ens=ZERO, eps_w=ZERO, eps_ss=ZERO;
  double Yfuel_conditional = ZERO;
  
  //Conditional average on fresh gas
  Yfuel_conditional = 0.95*0.05518;//IPs.Fresh_Fuel_Mass_Fraction;
  u_ave = u_average;
  v_ave = v_average;  
  w_ave = w_average;  
    
  for (int p = 0; p < LocalSolnBlockList.Nblk; p++) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl; i <= Solution_Block[p].ICu; ++i) {
        for (int j  = Solution_Block[p].JCl; j <= Solution_Block[p].JCu; ++j) {
           for (int k  = Solution_Block[p].KCl; k <= Solution_Block[p].KCu; ++k) {
              if (Solution_Block[p].W[i][j][k].spec[0].c >= Yfuel_conditional) {
                local_vol = Solution_Block[p].Grid.volume(i,j,k);
                total_vol += local_vol;
                u_p += sqr(Solution_Block[p].W[i][j][k].v.x - u_ave) * local_vol;
                v_p += sqr(Solution_Block[p].W[i][j][k].v.y - v_ave) * local_vol;
                w_p += sqr(Solution_Block[p].W[i][j][k].v.z - w_ave) * local_vol;
                vis = Solution_Block[p].W[i][j][k].mu()/Solution_Block[p].W[i][j][k].rho;
/*              vis = Solution_Block[p].W[i][j][k].mu_t(Solution_Block[p].dWdx[i][j][k],
                                                        Solution_Block[p].dWdy[i][j][k],
                                                        Solution_Block[p].dWdz[i][j][k],
                                                        Solution_Block[p].Flow_Type,Solution_Block[p].Grid.volume(i,j,k))/
                      (Solution_Block[p].W[i][j][k].rho); */
                vis_ave += vis*local_vol;
                ens += Solution_Block[p].W[i][j][k].Enstrophy(Solution_Block[p].dWdx[i][j][k],
                                                              Solution_Block[p].dWdy[i][j][k],
                                                              Solution_Block[p].dWdz[i][j][k]) * local_vol;
                eps_w += 2.0*vis* Solution_Block[p].W[i][j][k].Enstrophy(Solution_Block[p].dWdx[i][j][k],
                                                                         Solution_Block[p].dWdy[i][j][k],
                                                                         Solution_Block[p].dWdz[i][j][k])*local_vol;
                eps_ss += 2.0*vis*(sqr(Solution_Block[p].W[i][j][k].abs_strain_rate(Solution_Block[p].dWdx[i][j][k],
                                                                                    Solution_Block[p].dWdy[i][j][k],
                                                                                    Solution_Block[p].dWdz[i][j][k]))/ 
                          2.0)*local_vol;
	      } /* endif */
	   } /* endfor */
	} /* endfor */
      } /* endfor*/
    } /* endif */
  } /* endfor */

  total_vol = CFFC_Summation_MPI(total_vol);
  vis_ave = CFFC_Summation_MPI(vis_ave);
  u_p = CFFC_Summation_MPI(u_p);
  v_p = CFFC_Summation_MPI(v_p);
  w_p = CFFC_Summation_MPI(w_p);
  ens = CFFC_Summation_MPI(ens);
  eps_w = CFFC_Summation_MPI(eps_w);
  eps_ss = CFFC_Summation_MPI(eps_ss);

  sqr_u = u_p/total_vol;
  vis_ave = vis_ave/total_vol;
  u_p = u_p/total_vol;
  v_p = v_p/total_vol;
  w_p = w_p/total_vol;
  ens = ens/total_vol;
  eps_w = eps_w/total_vol;
  eps_ss = eps_ss/total_vol;
  
  double u_rms = sqrt((u_p + v_p + w_p)/3.0);
  double Taylor_scale, Kolmogorov_scale, Re_Taylor, L11; 
  double l_1, l_2;

  Kolmogorov_scale = pow(pow(vis_ave, THREE)/eps_w, 0.25);

  if (ens == ZERO) {
    Taylor_scale = ZERO;
  } else {
    Taylor_scale = u_rms*sqrt(15.0)*Kolmogorov_scale;
    //sqrt(10.0*vis_ave*1.5*u_rms*u_rms/ens);//sqrt(TWO*u_rms*u_rms/(TWO*ens));//?????????
  }

  Re_Taylor = u_rms*Taylor_scale/vis_ave;
  L11= 0.09*pow(1.5*u_rms*u_rms, 1.5)/eps_ss;

  if (eps_w > 0.0) {
    l_1 = 0.42*pow(u_rms, 3.0)/eps_w;
  } else {
    l_1 = 0.0;
  }

  if (eps_ss > 0.0) {
    l_2 = 0.42*pow(u_rms, 3.0)/eps_ss;
  } else {
    l_2 = 0.0;
  }

  if (CFFC_Primary_MPI_Processor()) {
    cout << "\n\n ==========================================================================\n"; 
    cout << " In physical space:\n";
    cout << "\n <u^2> = "<< u_p <<"  "<< "<v^2> = "<< v_p <<"  "<< "<v^2> = "<< w_p <<"  "
	 << "u_rms  = " << u_rms <<"  "
	 << "\n <u> = " << u_ave <<"  "<< "<v> = " << v_ave <<"  "<< "<w> = " << w_ave <<"  "
	 << "ens = "<< ens <<"  " 
	 << "\n eps_w = "<< eps_w <<"  "<< "eps_ss = "<< eps_ss <<"  "
	 << "l_1 = "<< l_1 <<"  "<< "l_2 = " << l_2 <<"  "
	 << "\n vis = "<< vis_ave << "  Re_Taylor = " << Re_Taylor <<"  "
	 << "\n Taylor_scale = " << Taylor_scale <<"  " 
	 << "Kolmogorov_scale = " << Kolmogorov_scale <<" "
         << "\n L11 = " << L11 <<  endl;
    cout << " ==========================================================================" << endl;
  } /* endif */

}

template <typename HEXA_BLOCK>
int Longitudinal_Correlation(Octree_DataStructure &OcTree,
                             AdaptiveBlock3D_ResourceList &Global_Soln_Block_List,
			     AdaptiveBlock3D_List &LocalSolnBlockList,
			     HEXA_BLOCK *Solution_Block,
			     Grid3D_Input_Parameters &IPs,
			     const double &u_ave, 
                             const double &sqr_u) {

  ofstream Out_Corr_Function;
  int error_flag = 0;
  int Nblks = Global_Soln_Block_List.Nused;
  double xref, yref, zref, uref, ucorr;
  double Yfuel_conditional = ZERO;
 
  int gblknum;
  double r, dr, max_r, fr, Lx, volume, local_vol, total_vol, R11, L11;
  Lx = IPs.Box_Width;
  dr = MILLION;

  //Conditional longitudinal correlation on fresh gas
  //  if (IPs.react_name != "NO_REACTIONS") {
  Yfuel_conditional = 0.95*0.05518;//IPs.Fresh_Fuel_Mass_Fraction;
    //  }

  for (int k = 0; k < LocalSolnBlockList.Nblk; ++k) {
    if (LocalSolnBlockList.Block[k].used == ADAPTIVEBLOCK3D_USED) {
      for (int ii = Solution_Block[k].ICl; ii <= Solution_Block[k].ICu; ++ii) {
        for (int jj = Solution_Block[k].JCl; jj <= Solution_Block[k].JCu; ++jj) {
           for (int kk = Solution_Block[k].KCl; kk <= Solution_Block[k].KCu; ++kk) {
	  if (Solution_Block[k].W[ii][jj][kk].spec[0].c >= Yfuel_conditional) {
	    dr = min(Solution_Block[k].Grid.Cell[ii][jj][kk].Xc.x - Solution_Block[k].Grid.Cell[ii-1][jj][kk].Xc.x, dr);
	  }
	}
      }
    }
   }
  }

  dr = CFFC_Minimum_MPI(dr);

  if (CFFC_Primary_MPI_Processor()) {
    Out_Corr_Function.open("Correlation_Function.dat", ios::out);
    if(Out_Corr_Function.fail()){
      cerr<<"\n Error opening file: Correlation_Function.dat to write" << endl;
      exit(1);
    }
  }
  
  int *new_blocks_CPU, *CPUs_in_new_blocks;
  int my_rank, undefined_rank, new_blocks_BLK;
  CPUs_in_new_blocks = new int[Global_Soln_Block_List.Ncpu];
  new_blocks_CPU = new int[Global_Soln_Block_List.Ncpu];
  
#ifdef _MPI_VERSION
    MPI::Intracomm new_comm;
    MPI::Group     big_group = MPI::COMM_WORLD.Get_group();
    MPI::Group     new_group;
    undefined_rank = MPI::UNDEFINED;
#else
    undefined_rank = -1;
#endif
  
  HEXA_BLOCK  SolnBlk_duplicated;
  OctreeBlock *octree_block_duplicated_ptr;
  int SolnBlk_duplicated_info_level;

  Vector3D Vcorr, Xcorr;
  Vcorr.zero();  Xcorr.zero();
  bool correlated_flag, flag = false;
  int count = 0, count1 = 0;
  r = ZERO;
  L11 = ZERO;

/*   if (IPs.i_Grid == GRID_TURBULENT_PREMIXED_FLAME && */
/*       IPs.react_name == "NO_REACTIONS") { */
/*     max_r = HALF*(Lx-dr); */
/*   } else  { */
    max_r = Lx-dr;
/*   } */
  
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
 
  while( r <= max_r) {
    total_vol = ZERO;
    R11 = ZERO;
    correlated_flag = false;
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

    for (int iCPU = 0; iCPU <= OcTree.Ncpu-1; ++iCPU ) { // Loop over available processors.
      for (int iBLK = 0; iBLK <= OcTree.Nblk-1; ++iBLK ) { // Loop over available blocks.
        if (OcTree.Blocks[iCPU][iBLK] != NULL) {
          if (OcTree.Blocks[iCPU][iBLK]->block.used) { // Check if the solution block is used.
        
            octree_block_duplicated_ptr = OcTree.Blocks[iCPU][iBLK];
            new_blocks_CPU[0] = octree_block_duplicated_ptr->block.info.cpu;
            new_blocks_BLK = octree_block_duplicated_ptr->block.info.blknum;
            CPUs_in_new_blocks[0] = new_blocks_CPU[0];

	    if (Global_Soln_Block_List.Ncpu > 1) {
	      for (int iNEW = 1 ; iNEW < Global_Soln_Block_List.Ncpu; ++iNEW ) {
		if (iNEW != new_blocks_CPU[0]) {
		  new_blocks_CPU[iNEW] = iNEW;
		} else {
		  new_blocks_CPU[iNEW] = 0;
		}
		CPUs_in_new_blocks[iNEW] = new_blocks_CPU[iNEW];
	      }
	    }
  
#ifdef _MPI_VERSION
	    new_group = big_group.Incl(Global_Soln_Block_List.Ncpu, CPUs_in_new_blocks);
	    new_comm = MPI::COMM_WORLD.Create(new_group);
#endif
/* 	    if (LocalSolnBlockList.ThisCPU == new_blocks_CPU[0]) { */
/* 	      Copy_Solution_Block(SolnBlk_duplicated, Solution_Block[new_blocks_BLK]); */
/* 	      SolnBlk_duplicated_info_level = LocalSolnBlockList.Block[new_blocks_BLK].info.level; */
/* 	    } */
/* #ifdef _MPI_VERSION */
/* 	    if (my_rank != undefined_rank) { */
/* 	      Broadcast_Solution_Block(SolnBlk_duplicated, new_comm, new_blocks_CPU[0]); */
/* 	      new_comm.Bcast(&SolnBlk_duplicated_info_level, 1, MPI::INT, 0); */
/* 	    } */
/* #endif */

#ifdef _MPI_VERSION
	    if (new_comm != MPI::COMM_NULL) new_comm.Free();
	    new_group.Free();
#endif
	    for (int i_ref = SolnBlk_duplicated.ICl; i_ref <= SolnBlk_duplicated.ICu; ++i_ref) {
	      for (int j_ref = SolnBlk_duplicated.JCl; j_ref <= SolnBlk_duplicated.JCu; ++j_ref) {
  	         for (int k_ref = SolnBlk_duplicated.KCl; k_ref <= SolnBlk_duplicated.KCu; ++k_ref) {
/*                 if (SolnBlk_duplicated.W[i_ref][j_ref][k_ref].spec[0].c >= Yfuel_conditional) { */
/* 		  xref = SolnBlk_duplicated.Grid.Cell[i_ref][j_ref][k_ref].Xc.x; */
/* 		  yref = SolnBlk_duplicated.Grid.Cell[i_ref][j_ref][k_ref].Xc.y; */
/* 		  zref = SolnBlk_duplicated.Grid.Cell[i_ref][j_ref][k_ref].Xc.z; */
/* 		  uref = SolnBlk_duplicated.W[i_ref][j_ref][k_ref].v.x; */
/* 		  flag = false; */
/* 		} else { */
/*                   continue; */
/* 		} */
	       
		for (int q = 0; q < LocalSolnBlockList.Nblk; ++q) {
		  if (LocalSolnBlockList.Block[q].used == ADAPTIVEBLOCK3D_USED) {
		    for (int i = Solution_Block[q].ICl; i <= Solution_Block[q].ICu; ++i) {
		      for (int j = Solution_Block[q].JCl; j <= Solution_Block[q].JCu; ++j) {
		         for (int k = Solution_Block[q].KCl; k <= Solution_Block[q].KCu; ++k) {

			if ((Solution_Block[q].W[i][j][k].spec[0].c >= Yfuel_conditional) &&
                            (Solution_Block[q].Grid.xfaceW(i,j,k).x < (xref + r) &&
			     Solution_Block[q].Grid.xfaceE(i,j,k).x > (xref + r))) {

			  // Finer reference block than local block
			  if (SolnBlk_duplicated_info_level >= LocalSolnBlockList.Block[q].info.level  &&
			      (Solution_Block[q].Grid.xfaceS(i,j,k).y < yref  &&  Solution_Block[q].Grid.xfaceN(i,j,k).y > yref)) {

			    if (Solution_Block[q].Grid.Cell[i][j][k].Xc.y == yref &&
				Solution_Block[q].Grid.Cell[i][j][k].Xc.x == (xref + r)) {
			      Vcorr.x = Solution_Block[q].W[i][j][k].v.x;
			    } else {
			      Vector2D dX;
			      dX.x = (xref + r) - Solution_Block[q].Grid.Cell[i][j][k].Xc.x;
			      dX.y = yref - Solution_Block[q].Grid.Cell[i][j][k].Xc.y;
			      Vcorr.x = Solution_Block[q].W[i][j][k].v.x + Solution_Block[q].phi[i][j][k].v.x
				* (Solution_Block[q].dWdx[i][j][k].v.x * dX.x + Solution_Block[q].dWdy[i][j][k].v.x * dX.y);
			    }
			    // correlated u
			    ucorr = Vcorr.x;
			    local_vol = Solution_Block[q].Grid.volume(i,j,k);
			    count1 ++;
			    flag = true;

			  // Coarser reference block than local block
			  } else if (SolnBlk_duplicated_info_level < LocalSolnBlockList.Block[q].info.level  &&
				     (Solution_Block[q].Grid.Cell[i][j][k].Xc.y < yref  &&  Solution_Block[q].Grid.Cell[i][j+1][k].Xc.y > yref)) {
			    Xcorr.x = xref + r;
			    Xcorr.y = yref;
                            Trilinear_Interpolation(Solution_Block[q].Grid.Cell[i-1][j][k].Xc, Solution_Block[q].W[i-1][j][k],
                                                    Solution_Block[q].Grid.Cell[i][j][k].Xc, Solution_Block[q].W[i][j][k],
                                                    Solution_Block[q].Grid.Cell[i][j-1][k].Xc, Solution_Block[q].W[i][j-1][k],
                                                    Solution_Block[q].Grid.Cell[i-1][j-1][k].Xc, Solution_Block[q].W[i-1][j-1][k],
                                                    Solution_Block[q].Grid.Cell[i-1][j][k-1].Xc, Solution_Block[q].W[i-1][j][k-1],
                                                    Solution_Block[q].Grid.Cell[i][j][k-1].Xc, Solution_Block[q].W[i][j][k-1],
                                                    Solution_Block[q].Grid.Cell[i][j-1][k-1].Xc, Solution_Block[q].W[i][j-1][k-1],
                                                    Solution_Block[q].Grid.Cell[i-1][j-1][k-1].Xc, Solution_Block[q].W[i-1][j-1][k-1],
                                                    Solution_Block[q].Grid.Node[i][j][k].X);
/* 			    Bilinear_Interpolation(Solution_Block[q].W[i][j][k].v, Solution_Block[q].Grid.Cell[i][j][k].Xc, */
/* 						   Solution_Block[q].W[i][j+1][k].v, Solution_Block[q].Grid.Cell[i][j+1][k].Xc, */
/* 						   Solution_Block[q].W[i+1][j+1][k].v, Solution_Block[q].Grid.Cell[i+1][j+1][k].Xc, */
/* 						   Solution_Block[q].W[i+1][j][k].v, Solution_Block[q].Grid.Cell[i+1][j][k].Xc, */
/* 						   Xcorr, Vcorr); */
			    // correlated u
			    ucorr = Vcorr.x;
			    local_vol = Solution_Block[q].Grid.volume(i,j,k);
			    count1 ++;
			    flag = true;
			  } // end if
                                                                                          
// 			  if (Soln_ptr[q].Grid.Cell[i][j].Xc.y == yref &&
// 			      Soln_ptr[q].Grid.Cell[i][j].Xc.x == (xref + r)) {
// 			    Vcorr.x = Soln_ptr[q].W[i][j].v.x;
// 			  } else {
// 			    Xcorr.x = xref + r;
// 			    Xcorr.y = yref;
// 			    // Bilinear_Interpolation(Soln_ptr[q].WnNW(i,j).v, Soln_ptr[q].Grid.nodeNW(i,j).X,
// 			    // 						 Soln_ptr[q].WnNE(i,j).v, Soln_ptr[q].Grid.nodeNE(i,j).X,
// 			    // 						 Soln_ptr[q].WnSE(i,j).v, Soln_ptr[q].Grid.nodeSE(i,j).X,
// 			    // 						 Soln_ptr[q].WnSW(i,j).v, Soln_ptr[q].Grid.nodeSW(i,j).X,
// 			    // 						 Xcorr, Vcorr);
// 			    Bilinear_Interpolation(Soln_ptr[q].W[i-1][j].v, Soln_ptr[q].Grid.Cell[i-1][j].Xc,
// 						   Soln_ptr[q].W[i][j+1].v, Soln_ptr[q].Grid.Cell[i][j+1].Xc,
// 						   Soln_ptr[q].W[i+1][j].v, Soln_ptr[q].Grid.Cell[i+1][j].Xc,
// 						   Soln_ptr[q].W[i][j-1].v, Soln_ptr[q].Grid.Cell[i][j-1].Xc,
// 						   Xcorr, Vcorr);
// 			  }
		  
			} // end if
		      }
		    }
		  } //end if
		  if (flag) break;
		} // end for
        
		if (flag) {
		  volume = SolnBlk_duplicated.Grid.volume(i_ref,j_ref,k_ref) + local_vol;
		  total_vol += volume;
		  R11 += (ucorr - u_ave)*(uref - u_ave)*volume;
		  correlated_flag = true;
		  count++;
		}
           
		} // end for
		 } // end for
	      }
	    }
	  }
	}
      }
    }
    if (CFFC_OR_MPI(correlated_flag)) {
      total_vol = CFFC_Summation_MPI(total_vol);
      R11 = CFFC_Summation_MPI(R11);
       
      R11 = R11/total_vol;     // Area-averaged two-point correlation
      fr = R11/sqr_u;        // Longitudinal autocorrelation function
      L11 += fr*dr;          // Integrate the above funtion to obtain L11
      if (CFFC_Primary_MPI_Processor()) {
        Out_Corr_Function << r << "  " << fr << endl;
      }
    }
       
    //cout << "\n->" << r;
    r += dr;
  } // end while
  
  count = CFFC_Summation_MPI(count);
  count1 = CFFC_Summation_MPI(count1);
    
  if (CFFC_Primary_MPI_Processor()) {
    Out_Corr_Function.close();
    cout << "\n\n *** L11 = " << L11 << " ***" << endl;
  }

  delete[] CPUs_in_new_blocks;   CPUs_in_new_blocks = NULL;
  delete[] new_blocks_CPU;       new_blocks_CPU = NULL;
  
  if (CPUs_in_new_blocks != NULL || new_blocks_CPU != NULL) error_flag = 1;
       
  return (error_flag);
}

#endif // _TURBULENT_VELOCITY_FIELD_INCLUDED 
