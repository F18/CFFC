/* HexaBlock.h: Header file defining 
                3D Hexahedral Mesh Solution Classes. */

#ifndef _HEXA_BLOCK_INCLUDED
#define _HEXA_BLOCK_INCLUDED

/* Include required CFFC header files. */

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // _INPUT_INCLUDED

#ifndef _CELL3D_INCLUDED
#include "../Grid/Cell3D.h"
#endif // _CELL3D_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "../Grid/Grid3DHexaBlock.h"
#endif // _GRID3D_HEXA_BLOCK_INCLUDED

#ifndef _GASCONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GASCONSTANTS_INCLUDED

#ifndef _TURBULENCE_MODELLING_INCLUDED
#include "../Turbulence/TurbulenceModelling.h"
#endif // TURBULENCE_MODELLING_INCLUDED   

/* Define the solution block in-use indicators. */

#define	HEXA_BLOCK_USED                            1
#define	HEXA_BLOCK_NOT_USED                        0

/* Define the structures and classes. */

#define	NUMBER_OF_RESIDUAL_VECTORS    3

template<class SOLN_pSTATE, class SOLN_cSTATE>
class Hexa_Block{
  protected:
 
  public:

   typedef SOLN_pSTATE Soln_pState;
   typedef SOLN_cSTATE Soln_cState;

   int              NCi,ICl,ICu;  // i-direction cell counters.
   int              NCj,JCl,JCu;  // j-direction cell counters.
   int              NCk,KCl,KCu;  // k-direction cell counters.
   int                   Nghost;  // Number of ghost cells.
   int           Freeze_Limiter;  // Limiter freezing indicator.
   int                Flow_Type;  // Flow type indicator.
   double                 ***dt;  // Local time step.
   
   SOLN_pSTATE             ***W;  // Primitive solution state.
   SOLN_cSTATE             ***U;  // Conserved solution state.
   Grid3D_Hexa_Block       Grid;  // pointer to the 3D Hexa grid.
   
   SOLN_cSTATE         ****dUdt;  // Solution residual.
   SOLN_pSTATE          ***dWdx;  // Unlimited solution gradient
                                  // (x-direction).
   SOLN_pSTATE          ***dWdy;  // Unlimited solution gradient
                                  // (y-direction).
   SOLN_pSTATE          ***dWdz;  // Unlimited solution gradient
                                  // (z-direction).
   SOLN_pSTATE           ***phi;  // Solution slope limiter.
   SOLN_cSTATE            ***Uo;  // Initial solution state.
   //   SOLN_cSTATE **FluxN,**FluxS,  // Boundary solution fluxes.
   //               **FluxE,**FluxW,
   //               **FluxT,**FluxB;

   SOLN_pSTATE     **WoN, **WoS,  // Boundary condition reference states for
                   **WoE, **WoW,  // north, south, east, west, top
                   **WoT, **WoB;  // bottom boundaries.
   
   // Only allocate for turbulent flow (depending on flow type indicator)
   Turbulent3DWallData ***WallData;
		      
   int Allocated; // Indicates whether or not the solution block has been allocated.

   /* Constructors. */

   Hexa_Block(void) {
      NCi=0; ICl=0; ICu=0; NCj=0; JCl=0; JCu=0;
      NCk=0; KCl=0; KCu=0; Nghost=0; 
      Flow_Type = FLOWTYPE_INVISCID; Freeze_Limiter = OFF;
      Allocated = HEXA_BLOCK_NOT_USED;
      dt = NULL; W = NULL; U = NULL; dUdt = NULL;
      dWdx = NULL; dWdy = NULL; dWdz = NULL;
      phi = NULL; Uo = NULL;
      //FluxN = NULL; FluxS = NULL;
      //FluxE = NULL; FluxW = NULL;
      //FluxT = NULL; FluxB = NULL;
      WoN = NULL; WoS = NULL; 
      WoE = NULL; WoW = NULL;
      WoT = NULL; WoB = NULL; WallData = NULL;
   }

   Hexa_Block(const int Ni, const int Nj, const int Nk, 
              const int Ng) {
      allocate(Ni, Nj, Nk, Ng);
   }

   Hexa_Block(const int Ni, const int Nj, const int Nk, 
              const int Ng, const int flow_type,  
              Grid3D_Hexa_Block &Grid2) {       
      allocate(Ni, Nj, Nk, Ng);
      Grid.Copy(Grid2);
   }
   
   /* Destructors. */
  ~Hexa_Block() {
      deallocate(); 
   }

   /* Allocate memory for structured hexahedrial solution block. */
   void allocate(void);
   void allocate(const int Ni, const int Nj, const int Nk, const int Ng);
 
   /* Deallocate memory for structured hexahedral solution block. */
   void deallocate(void);
   
   /* Return primitive solution state at specified node. */
   SOLN_pSTATE Wn(const int ii, const int jj, const int kk);
   
   /* Retern conserverd solution state at specified node. */
   SOLN_cSTATE Un(const int ii, const int jj, const int kk);
   
   /* Return primitive solution state at cell nodes. */
   SOLN_pSTATE WnNW(const int ii, const int jj, const int kk);
   SOLN_pSTATE WnNE(const int ii, const int jj, const int kk);
   SOLN_pSTATE WnSE(const int ii, const int jj, const int kk);
   SOLN_pSTATE WnSW(const int ii, const int jj, const int kk);
   
   /* Return conserved solution state at cell nodes. */
   SOLN_cSTATE UnNW(const int ii, const int jj, const int kk);
   SOLN_cSTATE UnNE(const int ii, const int jj, const int kk);
   SOLN_cSTATE UnSE(const int ii, const int jj, const int kk);
   SOLN_cSTATE UnSW(const int ii, const int jj, const int kk);
   
   /* Set flags for limiter evaluation. */
   void evaluate_limiters(void);
   void freeze_limiters(void);
   
   /* Number of solution state variables. */
   int NumVar(void);

   /* Other important member functions. */

   void Create_Block(Grid3D_Hexa_Block &Grid2);

   void Copy(Hexa_Block &Block2);

   void Broadcast(void);

#ifdef _MPI_VERSION
   void Broadcast(MPI::Intracomm &Communicator);
#endif

   void Output_Tecplot(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs,
                       const int Number_of_Time_Steps,
                       const double &Time,
                       const int Block_Number,
                       const int Output_Title,
                       ostream &Out_File);

   void Output_Cells_Tecplot(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs,
                             const int Number_of_Time_Steps,
                             const double &Time,
                             const int Block_Number,
                             const int Output_Title,
                             ostream &Out_File);

   void Output_Nodes_Tecplot(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs,
                             const int Number_of_Time_Steps,
                             const double &Time,
                             const int Block_Number,
                             const int Output_Title,
                             ostream &Out_File);

   int ICs(const int i_ICtype, 
           Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs);

   void BCs(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs);

   double CFL(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs);

   void Set_Global_TimeStep(const double &Dt_min);

   double L1_Norm_Residual(void);

   double L2_Norm_Residual(void);

   double Max_Norm_Residual(void);

   int WtoU(void);

   void Linear_Reconstruction_LeastSquares(const int i, 
                                           const int j, 
                                           const int k, 
                                           const int Limiter);

   void Linear_Reconstruction_LeastSquares(const int Limiter);

   int Wall_Shear(void);
   
   double Wall_Friction_Velocity(const int i, 
                                 const int j, 
                                 const int k);

   int dUdt_Multistage_Explicit(const int i_stage, 
                                Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs);

   int Update_Solution_Multistage_Explicit(const int i_stage, 
                                           Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs);

   /* MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING. */

   /* Load send message passing buffer. */
   int LoadSendBuffer(double *buffer,
                      int &buffer_count,
                      const int buffer_size,
                      const int i_min,
                      const int i_max,
                      const int i_inc,
                      const int j_min,
                      const int j_max,
                      const int j_inc,
                      const int k_min,
                      const int k_max,
                      const int k_inc);

   int LoadSendBuffer_F2C(double *buffer,
                          int &buffer_count,
                          const int buffer_size,
                          const int i_min,
                          const int i_max,
                          const int i_inc,
                          const int j_min,
                          const int j_max,
			  const int j_inc,
			  const int k_min,
			  const int k_max,
			  const int k_inc);

   int LoadSendBuffer_C2F(double *buffer,
                          int &buffer_count,
                          const int buffer_size,
                          const int i_min,
                          const int i_max,
                          const int i_inc,
                          const int j_min,
                          const int j_max,
			  const int j_inc,
			  const int k_min,
			  const int k_max,
			  const int k_inc);

   /* Unload receive message passing buffer. */

   int UnloadReceiveBuffer(double *buffer,
                           int &buffer_count,
                           const int buffer_size,
                           const int i_min,
                           const int i_max,
                           const int i_inc,
                           const int j_min,
                           const int j_max,
			   const int j_inc,
			   const int k_min,
			   const int k_max,
			   const int k_inc);

   int UnloadReceiveBuffer_F2C(double *buffer,
                               int &buffer_count,
                               const int buffer_size,
                               const int i_min,
                               const int i_max,
                               const int i_inc,
                               const int j_min,
                               const int j_max,
                               const int j_inc,
                               const int k_min,
			       const int k_max,
			       const int k_inc);

   int UnloadReceiveBuffer_C2F(double *buffer,
                               int &buffer_count,
                               const int buffer_size,
                               const int i_min,
                               const int i_max,
                               const int i_inc,
                               const int j_min,
                               const int j_max,
	  		       const int j_inc,
			       const int k_min,
			       const int k_max,
			       const int k_inc);

   /* Subcell solution reconstruction within given computational cell. */

   void SubcellReconstruction(const int i,
                              const int j,
                              const int k,
                              const int Limiter);

   /* Load and unload conservative flux message passing buffer. */

   int LoadSendBuffer_Flux_F2C(double *buffer,
                               int &buffer_count,
                               const int buffer_size,
                               const int i_min,
                               const int i_max,
                               const int i_inc,
                               const int j_min,
                               const int j_max,
			       const int j_inc,
			       const int k_min,
			       const int k_max,
			       const int k_inc);

   int UnloadReceiveBuffer_Flux_F2C(double *buffer,
                                    int &buffer_count,
                                    const int buffer_size,
                                    const int i_min,
                                    const int i_max,
                                    const int i_inc,
                                    const int j_min,
                                    const int j_max,
				    const int j_inc,
				    const int k_min,
				    const int k_max,
				    const int k_inc);
   
  private:
   //copy and assignment are not permitted ...
   Hexa_Block(const Hexa_Block &Soln);
   Hexa_Block &operator =(const Hexa_Block &Soln);
   
};

 /* Input-output operators. */

template<class SOLN_pSTATE, class SOLN_cSTATE>
ostream &operator << (ostream &out_file,
                      const Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> &Soln);

template<class SOLN_pSTATE, class SOLN_cSTATE>
istream &operator >> (istream &in_file,
                      Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> &Soln);

/******************************************************************
 * Hexa_Block::allocate -- Allocate memory.                       *
 ******************************************************************/

template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::allocate(void){
     
   W = new SOLN_pSTATE**[NCi]; 
   U = new SOLN_cSTATE**[NCi];
   Uo = new SOLN_cSTATE**[NCi];
   dUdt = new SOLN_cSTATE***[NCi];
   dWdx = new SOLN_pSTATE**[NCi];
   dWdy = new SOLN_pSTATE**[NCi];
   dWdz = new SOLN_pSTATE**[NCi];
   phi = new SOLN_pSTATE**[NCi]; 
   dt = new double**[NCi];
   
   for (int i = 0; i <= NCi-1 ; ++i) {
      W[i] = new SOLN_pSTATE*[NCj]; 
      U[i] = new SOLN_cSTATE*[NCj];
      Uo[i] = new SOLN_cSTATE*[NCj];
      dUdt[i] = new SOLN_cSTATE**[NCj];
      dWdx[i] = new SOLN_pSTATE*[NCj];
      dWdy[i] = new SOLN_pSTATE*[NCj];
      dWdz[i] = new SOLN_pSTATE*[NCj];
      phi[i] = new SOLN_pSTATE*[NCj];
      dt[i] = new double*[NCj];
      for (int j = 0; j <= NCj-1 ; ++j) {
         W[i][j] = new SOLN_pSTATE[NCk]; 
         U[i][j] = new SOLN_cSTATE[NCk];
         Uo[i][j] = new SOLN_cSTATE[NCk];
         dUdt[i][j] = new SOLN_cSTATE*[NCk];
         dWdx[i][j] = new SOLN_pSTATE[NCk]; 
         dWdy[i][j] = new SOLN_pSTATE[NCk];
         dWdz[i][j] = new SOLN_pSTATE[NCk];
         phi[i][j] = new SOLN_pSTATE[NCk]; 
         dt[i][j] = new double [NCk];
         for (int k = 0; k <= NCk-1 ; ++k) {
            dUdt[i][j][k]=new SOLN_cSTATE[NUMBER_OF_RESIDUAL_VECTORS];
         } /* endfor */
      } /* endfor */
   } /* endfor */
   
   SOLN_cSTATE U_VACUUM;
   U_VACUUM.Vacuum();
   SOLN_pSTATE W_VACUUM;
   W_VACUUM.Vacuum();

   // Set the solution residuals, gradients, limiters, and other values to zero.
   for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
      for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
         for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
            Uo[i][j][k] = U_VACUUM; 
            for (int n = 0 ; n <= NUMBER_OF_RESIDUAL_VECTORS-1 ; ++n) {
               dUdt[i][j][k][n] = U_VACUUM;
            } /* endfor */
            dWdx[i][j][k] = W_VACUUM; 
            dWdy[i][j][k] = W_VACUUM;
            dWdz[i][j][k] = W_VACUUM; 
            phi[i][j][k] =  W_VACUUM; 
            dt[i][j][k] = ZERO;
         } /* endfor */
      } /* endfor */
   } /* endfor */
 
   //Boundary References
   WoW = new SOLN_pSTATE *[NCj]; WoE = new SOLN_pSTATE *[NCj];
   for (int j = 0; j < NCj; ++j) {
      WoW[j] = new SOLN_pSTATE[NCk]; WoE[j] = new SOLN_pSTATE[NCk];
   } /* endfor */

   WoN = new SOLN_pSTATE *[NCi]; WoS = new SOLN_pSTATE *[NCi];
   WoT = new SOLN_pSTATE *[NCi]; WoB = new SOLN_pSTATE *[NCi];
   for (int i = 0; i < NCi; ++i) {
      WoN[i] = new SOLN_pSTATE[NCk]; WoS[i] = new SOLN_pSTATE[NCk];
      WoT[i] = new SOLN_pSTATE[NCj]; WoB[i] = new SOLN_pSTATE[NCj];
   } /* endfor */

   // Only allocate for turbulent flows
   if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
       Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      WallData = new Turbulent3DWallData**[NCi]; 
      for (int i = 0; i <= NCi-1 ; ++i) {
         WallData[i] = new Turbulent3DWallData*[NCj]; 
         for (int j = 0; j <= NCj-1 ; ++j) {
            WallData[i][j] = new Turbulent3DWallData[NCk]; 
         } /* endfor */
      } /* endfor */
   } /* endif */

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::allocate(const int Ni, 
                                                    const int Nj, 
                                                    const int Nk, 
                                                    const int Ng) {

    assert(Ni > 1 && Nj > 1  && Nk > 1 && Ng > 1 && !Allocated); 
    NCi=Ni+2*Ng; ICl=Ng; ICu=Ni+Ng-1; 
    NCj=Nj+2*Ng; JCl=Ng; JCu=Nj+Ng-1;
    NCk=Nk+2*Ng; KCl=Ng; KCu=Nk+Ng-1; Nghost=Ng;
    Flow_Type = FLOWTYPE_INVISCID; Freeze_Limiter = OFF;
    Allocated = HEXA_BLOCK_USED;
    Grid.allocate(Ni, Nj, Nk, Ng);
    allocate();

}

/**************************************************************************
 * Hexa_Block::deallocate -- Deallocate memory.                           *
 **************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::deallocate(void) {

   if (Allocated) {   
      Grid.deallocate(); 
    
      for (int i = 0; i <= NCi-1; ++i) {
         for (int j = 0; j <= NCj-1 ; ++j) {
            delete []W[i][j];  W[i][j] = NULL; 
            delete []U[i][j];  U[i][j] = NULL;
            delete []Uo[i][j]; Uo[i][j] = NULL;
            for (int k = 0; k <= NCk-1 ; ++k) {
               delete []dUdt[i][j][k]; dUdt[i][j][k] = NULL; 
            } /* endfor */
            delete []dUdt[i][j]; dUdt[i][j] = NULL;
            delete []dWdx[i][j]; dWdx[i][j] = NULL; 
            delete []dWdy[i][j]; dWdy[i][j] = NULL;
            delete []dWdz[i][j]; dWdz[i][j] = NULL;
            delete []phi[i][j]; phi[i][j] = NULL; 
            delete []dt[i][j]; dt[i][j] = NULL; 
         } /* endfor */
         delete []W[i]; W[i] = NULL; 
         delete []U[i]; U[i] = NULL;
         delete []Uo[i]; Uo[i] = NULL;
         delete []dUdt[i]; dUdt[i] = NULL;
         delete []dWdx[i]; dWdx[i] = NULL; 
         delete []dWdy[i]; dWdy[i] = NULL;
         delete []dWdz[i]; dWdz[i] = NULL;
         delete []phi[i]; phi[i] = NULL; 
         delete []dt[i]; dt[i] = NULL; 
      } /* endfor */
   
      delete []W; W = NULL; delete []U; U = NULL;
      delete []Uo; Uo = NULL; delete []dUdt; dUdt = NULL;
      delete []dWdx; dWdx = NULL; delete []dWdy; dWdy = NULL; 
      delete []dWdz; dWdz = NULL; 
      delete []phi; phi = NULL; delete []dt; dt = NULL; 

      //Boundary references
      for (int j = 0; j <= NCj-1; ++j) {
         delete []WoW[j]; WoW[j] = NULL;
         delete []WoE[j]; WoE[j] = NULL;
      } /* endfor */
      delete []WoW; WoW = NULL;
      delete []WoE; WoE = NULL;
   
      for (int i = 0; i <= NCi-1; ++i) {
         delete []WoN[i]; WoN[i] = NULL;
         delete []WoS[i]; WoS[i] = NULL;
         delete []WoT[i]; WoT[i] = NULL;
         delete []WoB[i]; WoB[i] = NULL;
      } /* endfor */
   
      delete []WoN; WoN = NULL; delete []WoS; WoS = NULL; 
      delete []WoE; WoE = NULL; delete []WoW; WoW = NULL;
      delete []WoT; WoT = NULL; delete []WoB; WoB = NULL;
   
      if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
          Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
         for (int i = 0; i <= NCi-1 ; ++i) {
            for (int j = 0; j <= NCj-1 ; ++j) {
               delete []WallData[i][j];  WallData[i][j] = NULL; 
            } /* endfor */
            delete []WallData[i]; WallData[i] = NULL; 
         } /* endfor */
         delete []WallData; WallData = NULL; 
      } /* endif */

      NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0;
      NCk = 0; KCl = 0; KCu = 0;  Nghost = 0;
      Allocated = HEXA_BLOCK_NOT_USED;
   } /* endif */

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::NumVar(void) {
   
   int num=0;
   int NUM_VAR_3D = W[num][num][num].NUM_VAR_3D;    

   return (int(NUM_VAR_3D));

}

/* Return primitive solution state at specified node. */
/* Trilinear interpolation --- in process */
template<class SOLN_pSTATE, class SOLN_cSTATE>
SOLN_pSTATE Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Wn(const int ii, const int jj, 
                                                     const int kk){
   return (Trilinear_Interpolation(
           Grid.Cell[ii-1][jj][kk].Xc, W[ii-1][jj][kk],
           Grid.Cell[ii][jj][kk].Xc, W[ii][jj][kk],
           Grid.Cell[ii][jj-1][kk].Xc, W[ii][jj-1][kk],
           Grid.Cell[ii-1][jj-1][kk].Xc, W[ii-1][jj-1][kk],
           Grid.Cell[ii-1][jj][kk-1].Xc, W[ii-1][jj][kk-1],
           Grid.Cell[ii][jj][kk-1].Xc, W[ii][jj][kk-1],
           Grid.Cell[ii][jj-1][kk-1].Xc, W[ii][jj-1][kk-1],
           Grid.Cell[ii-1][jj-1][kk-1].Xc, W[ii-1][jj-1][kk-1],
           Grid.Node[ii][jj][kk].X));

}

/********************************************************
 * Routine: Create_Block                                *
 *                                                      *
 * Create a new hexahedral solution block.              *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Create_Block(Grid3D_Hexa_Block &Grid2) {

   if (Grid2.Allocated) {
      if (NCi != Grid2.NCi ||
          NCj != Grid2.NCj ||
          NCk != Grid2.NCk ||
          Nghost != Grid2.Nghost) {
	 if (Allocated) {
           deallocate();
         } /*endif */

         allocate(Grid2.NCi-2*Grid2.Nghost, 
                  Grid2.NCj-2*Grid2.Nghost, 
                  Grid2.NCk-2*Grid2.Nghost, 
                  Grid2.Nghost);
      } /* endif */

      Grid.Copy(Grid2);
   } /* endif */
 
}

/********************************************************
 * Routine: Copy                                        *
 *                                                      *
 * Make a copy of solution block contents contained in  *
 * the hexahedral solution block Block2.                *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Copy(Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> &Block2) {

   if (Block2.Allocated) {
      //  Allocate memory as required.
      if (NCi != Block2.NCi ||
          NCj != Block2.NCj ||
          NCk != Block2.NCk ||
          Nghost != Block2.Nghost) {
	 if (Allocated) {
           deallocate();
         } /*endif */

         allocate(Block2.NCi-2*Block2.Nghost, 
                  Block2.NCj-2*Block2.Nghost, 
                  Block2.NCk-2*Block2.Nghost, 
                  Block2.Nghost);
      } /* endif */

      // Copy the grid.
      Grid.Copy(Block2.Grid);

      // Assign flow type.
      Flow_Type = Block2.Flow_Type;

      // Copy the freeze limite indicator.
      Freeze_Limiter = Block2.Freeze_Limter;

      // Copy the solution, solutioon residuals, gradients, limiters, and 
      // other stored values.
      for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
         for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
            for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
               U[i][j][k] = Block2.U[i][j][k];
               W[i][j][k] = Block2.W[i][j][k];
               Uo[i][j][k] = Block2.Uo[i][j][k];
               for (int n = 0 ; n <= NUMBER_OF_RESIDUAL_VECTORS-1 ; ++n) {
                  dUdt[i][j][k][n] = Block2.dUdt[i][j][k][n];
               } /* endfor */
               dWdx[i][j][k] = Block2.dWdx[i][j][k]; 
               dWdy[i][j][k] = Block2.dWdy[i][j][k];
               dWdz[i][j][k] = Block2.dWdz[i][j][k]; 
               phi[i][j][k] = Block2.phi[i][j][k]; 
               dt[i][j][k] = Block2.dt[i][j][k];
            } /* endfor */
         } /* endfor */
      } /* endfor */
 
      // Copy boundary reference states.
      for (int k = 0; k < NCk; ++k) {
         for (int j = 0; j < NCj; ++j) {
            WoW[j][k] = Block2.WoW[j][k]; 
            WoE[j][k] = Block2.WoE[j][k];
         } /* endfor */
      } /* endfor */

      for (int k = 0; k < NCk; ++k) {
         for (int i = 0; i < NCi; ++i) {
            WoN[i][k] = Block2.WoN[i][k]; 
            WoS[i][k] = Block2.WoS[i][k];
         } /* endfor */
      } /* endfor */

      for (int j = 0; j < NCj; ++j) {
         for (int i = 0; i < NCi; ++i) {
            WoT[i][j] = Block2.WoT[i][j]; 
            WoB[i][j] = Block2.WoB[i][j];
         } /* endfor */
      } /* endfor */

      // Copy wall data for turbulent flows only.
      if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
          Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
         for (int k = 0; k < NCk; ++k) {
            for (int j = 0; j < NCj; ++j) {
               for (int i = 0; i < NCi; ++i) {
                  WallData[i][j][k] = Block2.WallData[i][j][k]; 
               } /* endfor */
            } /* endfor */
         } /* endfor */
      } /* endif */

   } /* endif */

}

/********************************************************
 * Routine: Broadcast                                   *
 *                                                      *
 * Broadcast the solution information for hexahedral    * 
 * solution block.                                      *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Broadcast(void) {

#ifdef _MPI_VERSION

#endif

}

#ifdef _MPI_VERSION
/********************************************************
 * Routine: Broadcast                                   *
 *                                                      *
 * Broadcast the solution information for hexahedral    * 
 * solution block.                                      *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Broadcast(MPI::Intracomm &Communicator) {

}
#endif

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution values at the nodes of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Output_Tecplot(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs,
                                                          const int Number_of_Time_Steps,
                                                          const double &Time,
                                                          const int Block_Number,
                                                          const int Output_Title,
                                                          ostream &Out_File) {

   SOLN_pSTATE W_node;
   
   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
   
   BCs(IPs);
   
   /* Output nodal solution data. */
   
   Out_File << setprecision(14);
   if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 3D Solution, "
               << "Time Step/Iteration Level = " << Number_of_Time_Steps
               << ", Time = " << Time
               << "\"" << "\n"
               << "VARIABLES = \"x\" \\ \n"
               << "\"y\" \\ \n"
               << "\"z\" \\ \n"
               << "\"rho\" \\ \n"
               << "\"u\" \\ \n"
               << "\"v\" \\ \n"
               << "\"w\" \\ \n"
               << "\"p\" \\ \n";
      //n species mass fractions names
      for (int i =0 ;i<W[0][0][0].ns ;i++) {
         Out_File <<"\"c_"<<W[0][0][0].specdata[i].Speciesname()<<"\" \\ \n";
      }
      
      Out_File <<"\"T\" \\ \n";
      
      Out_File<< "ZONE T =  \"Block Number = " << Block_Number
              << "\" \\ \n"
              << "I = " << Grid.INu -  Grid.INl + 1 << " \\ \n"
              << "J = " << Grid.JNu -  Grid.JNl + 1 << " \\ \n"
              << "K = " << Grid.KNu -  Grid.KNl + 1 << " \\ \n"
              << "DATAPACKING = POINT \n";
      
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
               << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
               << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";              
   } /* endif */
   
   
   for (int k  =  Grid.KNl ; k <=  Grid.KNu ; ++k) {
      for (int j  =  Grid.JNl ; j <=  Grid.JNu ; ++j) {
         for (int i =  Grid.INl ; i <=  Grid.INu ; ++i) {
            W_node = Wn(i, j, k);
            Out_File << " "  << Grid.Node[i][j][k].X << W_node;
            Out_File.setf(ios::scientific);
            Out_File << " " << W_node.T() << "\n";
            Out_File.unsetf(ios::scientific);
         } /* endfor */
      } /* endfor */
   } /* endfor */

   Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Output_Cells_Tecplot(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs,
                                                                const int Number_of_Time_Steps,
                                                                const double &Time,
                                                                const int Block_Number,
                                                                const int Output_Title,
                                                                ostream &Out_File) {

   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
  
   BCs(IPs);

   /* Output cell-centred solution data. */

   Out_File << setprecision(14);

   if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 3D Solution, "
               << "Time Step/Iteration Level = " << Number_of_Time_Steps
               << ", Time = " << Time
               << "\"" << "\n"
               << "VARIABLES = \"x\" \\ \n"
               << "\"y\" \\ \n"
               << "\"z\" \\ \n"
               << "\"rho\" \\ \n"
               << "\"u\" \\ \n"
               << "\"v\" \\ \n"
               << "\"w\" \\ \n"
               << "\"p\" \\ \n";
                
      //n species mass fractions names
      for (int i =0;i<W[0][0][0].ns ;i++) {
         Out_File <<"\"c"<<W[0][0][0].specdata[i].Speciesname()<<"\" \\ \n";
      }
     
      Out_File <<"\"T\" \\ \n";
      Out_File <<"\"R\" \\ \n";

      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " <<  ICu -  ICl + 2*Nghost +1 << " \\ \n"
               << "J = " <<  JCu -  JCl + 2*Nghost +1 << " \\ \n"
               << "K = " <<  KCu -  KCl + 2*Nghost +1 << " \\ \n"
               << "DATAPACKING = POINT \n";
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " <<  Grid.ICu -  Grid.ICl + 2* Nghost + 1 << " \\ \n"
               << "J = " <<  Grid.JCu -  Grid.JCl + 2* Nghost + 1 << " \\ \n"
               << "K = " <<  Grid.KCu -  Grid.KCl + 2* Nghost + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";
      
   } /* endif */

   for (int k =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k) {
      for (int j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j) {
         for (int i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i) {
            Out_File << " "  <<  Grid.Cell[i][j][k].Xc << W[i][j][k];
            Out_File.setf(ios::scientific);
            Out_File << " " <<  W[i][j][k].T() << " " <<W[i][j][k].Rtot()<< "\n";
            Out_File.unsetf(ios::scientific);
         } /* endfor */
      } /* endfor */
   } /* endfor */
   
   Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the solution values at the nodes of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Output_Nodes_Tecplot(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs,
                                                                const int Number_of_Time_Steps,
                                                                const double &Time,
                                                                const int Block_Number,
                                                                const int Output_Title,
                                                                ostream &Out_File) {

   SOLN_pSTATE W_node;
   
   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
   
   BCs(IPs);
   
   /* Output nodal solution data. */
   
   Out_File << setprecision(14);
   if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 3D Solution, "
               << "Time Step/Iteration Level = " << Number_of_Time_Steps
               << ", Time = " << Time
               << "\"" << "\n"
               << "VARIABLES = \"x\" \\ \n"
               << "\"y\" \\ \n"
               << "\"z\" \\ \n"
               << "\"rho\" \\ \n"
               << "\"u\" \\ \n"
               << "\"v\" \\ \n"
               << "\"w\" \\ \n"
               << "\"p\" \\ \n";
      //n species mass fractions names
      for (int i =0 ;i<W[0][0][0].ns ;i++) {
         Out_File <<"\"c_"<<W[0][0][0].specdata[i].Speciesname()<<"\" \\ \n";
      }
      
      Out_File <<"\"T\" \\ \n";
      
      Out_File<< "ZONE T =  \"Block Number = " << Block_Number
              << "\" \\ \n"
              << "I = " << Grid.INu - Grid.INl + 1 + 2 << " \\ \n"
              << "J = " << Grid.JNu - Grid.JNl + 1 + 2 << " \\ \n"
              << "K = " << Grid.KNu - Grid.KNl + 1 + 2 << " \\ \n"
              << "DATAPACKING = POINT \n";
      
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " << Grid.INu - Grid.INl + 1 + 2 << " \\ \n"
               << "J = " << Grid.JNu - Grid.JNl + 1 + 2 << " \\ \n"
               << "K = " << Grid.KNu - Grid.KNl + 1 + 2 << " \\ \n"
               << "DATAPACKING = POINT \n";
   } /* endif */
   
   
   for (int k  = Grid.KNl - 1; k <= Grid.KNu + 1; ++k) {
      for (int j  = Grid.JNl - 1; j <= Grid.JNu + 1; ++j) {
         for (int i = Grid.INl - 1; i <= Grid.INu + 1; ++i) {
	    W_node = Wn(i, j, k);
            Out_File << " "  << Grid.Node[i][j][k].X << W_node;
            Out_File.setf(ios::scientific);
            Out_File << " " << W_node.T() << "\n";
            Out_File.unsetf(ios::scientific);
         } /* endfor */
      } /* endfor */
   } /* endfor */

   Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Apply initial conditions for the hexahedral          *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::ICs(const int i_ICtype,
                                              Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs){
   
   int i, j, k;
   
   double dpdx, dpdy, dpdz, delta_pres;
   double di, U_axi, Um;
      
   Flow_Type = IPs.i_Flow_Type;
   
   SOLN_pSTATE Wl, Wr;
   
   switch(i_ICtype) {
      case IC_VISCOUS_COUETTE :
         dpdx = IPs.Pressure_Gradient.x;  
         delta_pres = dpdx*IPs.IP_Grid.Box_Length;
         //starts with linear pressure gradient 
         for (  k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for (  j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  /* velocity is function of z only (2D flow) */
                  W[i][j][k].v.x = HALF*(-dpdx)/W[i][j][k].mu()*(Grid.Cell[i][j][k].Xc.z + 0.5*IPs.IP_Grid.Box_Height)*
                     (Grid.Cell[i][j][k].Xc.z - 0.5*IPs.IP_Grid.Box_Height) +
                     IPs.Moving_wall_velocity.x*(Grid.Cell[i][j][k].Xc.z/IPs.IP_Grid.Box_Height + HALF);
                  // assuming the wall velocity is parallel to x direction
                  W[i][j][k].p = IPs.Wo.p - (i-ICl-Nghost)*delta_pres/(ICu-ICl);	 
                  //start
                  if( i == ICl-Nghost || i == ICl-Nghost+1 ){
                    W[i][j][k].p = IPs.Wo.p;
                  }
                  //end
                  if( i == ICu+Nghost || i == ICu+Nghost-1 ){
                     W[i][j][k].p = IPs.Wo.p - delta_pres; 
                  }
                  U[i][j][k] = W[i][j][k].U( );
	       } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break; 

      case IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_X :
         dpdx = IPs.Pressure_Gradient.x;  
         delta_pres = dpdx*IPs.IP_Grid.Box_Length;
         //starts with linear pressure gradient 
         for (  k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for (  j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (  i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  /* velocity is function of z only (2D flow) */
                  W[i][j][k].v.x = HALF*(-dpdx)/W[i][j][k].mu()*(Grid.Cell[i][j][k].Xc.z + 0.5*IPs.IP_Grid.Box_Height)*
                     (Grid.Cell[i][j][k].Xc.z - 0.5*IPs.IP_Grid.Box_Height) +
                     IPs.Moving_wall_velocity.x*(Grid.Cell[i][j][k].Xc.z/IPs.IP_Grid.Box_Height + HALF);
                  // assuming the wall velocity is parallel to x direction
                  W[i][j][k].p = IPs.Wo.p - (Grid.Cell[i][j][k].Xc.x)*  delta_pres/IPs.IP_Grid.Box_Length;	 
                  U[i][j][k] = W[i][j][k].U( );
	       } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break; 

      case IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Y :
         dpdy = IPs.Pressure_Gradient.y;  
         delta_pres = dpdy*IPs.IP_Grid.Box_Width;
         //starts with linear pressure gradient 
         for (  k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for (  j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (  i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  /* velocity is function of z only (2D flow) */
                  W[i][j][k].v.y = HALF*(-dpdy)/W[i][j][k].mu()*(Grid.Cell[i][j][k].Xc.z + 0.5*IPs.IP_Grid.Box_Height)*
                     (Grid.Cell[i][j][k].Xc.z - 0.5*IPs.IP_Grid.Box_Height) +
                     IPs.Moving_wall_velocity.y*(Grid.Cell[i][j][k].Xc.z/IPs.IP_Grid.Box_Height + HALF);
                  //assuming the wall velocity is parallel to x direction
                  W[i][j][k].p = IPs.Wo.p - (Grid.Cell[i][j][k].Xc.y)*delta_pres/IPs.IP_Grid.Box_Width;	 
                  //start
                  if( j == JCl-Nghost || j == JCl-Nghost+1 ){
                     W[i][j][k].p = IPs.Wo.p;
                  }
                  //end
                  if( j == JCu+Nghost || j == JCu+Nghost-1 ){
                     W[i][j][k].p = IPs.Wo.p - delta_pres; 
                  }
                  U[i][j][k] = W[i][j][k].U( );
 	       } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break; 

      case IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Z :
         dpdz = IPs.Pressure_Gradient.z;  
         delta_pres = dpdz*IPs.IP_Grid.Box_Height;
         //starts with linear pressure gradient 
         for (  k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for (  j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (  i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  /* velocity is function of y only (2D flow) */
                  W[i][j][k].v.z = HALF*(-dpdz)/W[i][j][k].mu()*(Grid.Cell[i][j][k].Xc.y + 0.5*IPs.IP_Grid.Box_Width)*
                     (Grid.Cell[i][j][k].Xc.y - 0.5*IPs.IP_Grid.Box_Width) +
                     IPs.Moving_wall_velocity.z*(Grid.Cell[i][j][k].Xc.y/IPs.IP_Grid.Box_Width + HALF);
                  W[i][j][k].p = IPs.Wo.p - (Grid.Cell[i][j][k].Xc.z)*delta_pres/IPs.IP_Grid.Box_Height;	 
                  //start
                 if( k == KCl-Nghost || k == KCl-Nghost+1 ){
                     W[i][j][k].p = IPs.Wo.p;
                  }
                  //end
                  if( k == KCu+Nghost || k == KCu+Nghost-1 ){
                     W[i][j][k].p = IPs.Wo.p - delta_pres; 
                  }
                  U[i][j][k] = W[i][j][k].U( );
	       } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break; 
      
      case IC_PRESSURE_GRADIENT_X :
         dpdx = IPs.Pressure_Gradient.x;  
         delta_pres = dpdx*IPs.IP_Grid.Box_Length;
         di =  1.0/sqr(IPs.IP_Grid.Box_Height/2.0)+ 1.0/sqr(IPs.IP_Grid.Box_Width/2.0);
         U_axi = delta_pres/(2.0*W[0][0][0].mu())/di;
         //starts with linear pressure gradient 
         for (  k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for (  j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (  i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  /* velocity is function of z only (2D flow) */
                  W[i][j][k].v.x = 0.5*(-dpdx)/W[i][j][k].mu()*(Grid.Cell[i][j][k].Xc.z + 0.5*IPs.IP_Grid.Box_Height)*
                     (Grid.Cell[i][j][k].Xc.z - 0.5*IPs.IP_Grid.Box_Height);
                  W[i][j][k].p = IPs.Wo.p - (i-ICl-Nghost)*delta_pres/(ICu-ICl);	 
                  //start
                  if( i == ICl-Nghost || i == ICl-Nghost+1 ){
                     W[i][j][k].p = IPs.Wo.p;
                  }
                  //end
                  if( i == ICu+Nghost || i == ICu+Nghost-1){
                     W[i][j][k].p = IPs.Wo.p - delta_pres; 
                  }
                  U[i][j][k] = W[i][j][k].U( );
	       } /* endfor */ 
	    } /* endfor */       
	 } /* endfor */
         break; 

      case IC_CHANNEL_FLOW :
         // this case is from John Laufer
         // Investigation of Turbulent Flow in a Two-Dimensional Channel
         // debugging purpose here (laminar flow)
         dpdx = IPs.Pressure_Gradient.x;  
         delta_pres = dpdx*IPs.IP_Grid.Box_Length;
         Um = 22.3;
         //Um = 7.07; // computed from Reynolds number 61600
         for (  k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
	    for (  j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (  i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  W[i][j][k].p = IPs.Wo.p - (i-ICl-Nghost)*delta_pres/(ICu-ICl);	 
                  //start
                  if( i == ICl-Nghost || i == ICl-Nghost+1 ){
                     W[i][j][k].p = IPs.Wo.p;
                  }
                  //end
                  if( i == ICu+Nghost || i == ICu+Nghost-1){
                     W[i][j][k].p = IPs.Wo.p - delta_pres; 
                  }
                  /* velocity is function of y only (2D flow) */
                  W[i][j][k].v.x = 0.5*(-dpdx)/W[i][j][k].mu()*(Grid.Cell[i][j][k].Xc.y + 0.5*IPs.IP_Grid.Box_Width)*
                     (Grid.Cell[i][j][k].Xc.y - 0.5*IPs.IP_Grid.Box_Width);
                  //conservative solution state
                  U[i][j][k] = W[i][j][k].U();
 	       } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break; 

      case IC_SHOCK_BOX :
         Wl = SOLN_pSTATE(DENSITY_STDATM, Vector3D_ZERO,
                          PRESSURE_STDATM,
                          IPs.mass_fractions);
         Wr = SOLN_pSTATE(DENSITY_STDATM*8.0, Vector3D_ZERO,
                          PRESSURE_STDATM*10.0, 
                          IPs.mass_fractions);
         for (k = KCl-Nghost; k <= KCu+Nghost; ++k) {
	    for (j = JCl-Nghost; j <= JCu+Nghost; ++j) {
               for (i = ICl-Nghost; i <= ICu+Nghost; ++i) {
		  if (Grid.Cell[i][j][k].Xc.x <= ZERO &&
                      Grid.Cell[i][j][k].Xc.y <= ZERO &&
                      Grid.Cell[i][j][k].Xc.z <= ZERO) {
                      W[i][j][k] = Wl; 
                  } else {
                      W[i][j][k] = Wr;
                  } /* end if */
                  U[i][j][k] = W[i][j][k].U();
	       } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break;

      case IC_UNIFORM :
      default:
         // Set the solution state to the initial state Wo[0].
	 for (k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
	    for ( j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  U[i][j][k] = W[i][j][k].U( );
               } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break;
      
   } //end of switch

   /* Set default values for the boundary conditions
      reference states. */

   for ( k = KCl-Nghost ; k<= KCu+Nghost; ++k ) {
      for ( j = JCl-Nghost ; j<= JCu+Nghost; ++j ){
         if ((k >= KCl && k <= KCu) && (j >= JCl && j <= JCu)) {
            WoW[j][k] = W[ICl][j][k];
            WoE[j][k] = W[ICu][j][k];
         } else if (j < JCl && k < KCl ) {
            WoW[j][k] = W[ICl][JCl][KCl];
            WoE[j][k] = W[ICu][JCl][KCl];
         } else if (j > JCu && k> KCu) {
            WoW[j][k] = W[ICl][JCu][KCu];
            WoE[j][k] = W[ICu][JCu][KCu];
         } else if(j < JCl &&(k >= KCl && k <= KCu)){
            WoW[j][k] = W[ICl][JCl][k];
            WoE[j][k] = W[ICu][JCl][k];
         } else if(j > JCu && (k >= KCl && k <= KCu)){
            WoW[j][k] = W[ICl][JCu][k];
            WoE[j][k] = W[ICu][JCu][k];
         } else if(k < KCl &&(j >= JCl && j <= JCu)){
            WoW[j][k] = W[ICl][j][KCl];
            WoE[j][k] = W[ICu][j][KCl];
         } else if(k > KCu && (j >= JCl && j <= JCu)){
            WoW[j][k] = W[ICl][j][KCu];
            WoE[j][k] = W[ICu][j][KCu];
         } else if(k > KCu && j < JCl ){
            WoW[j][k] = W[ICl][JCl][KCu];
            WoE[j][k] = W[ICu][JCl][KCu];
         } else if(k < KCl && j > JCu){
            WoW[j][k] = W[ICl][JCu][KCl];
            WoE[j][k] = W[ICu][JCu][KCl];
         }
      } /* endfor */ 
   } /* endfor */
    
   for (  k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
      for (  i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
         if ((k >= KCl && k <= KCu) && (i >= ICl && i <= ICu)) {
            WoS[i][k] = W[i][JCl][k];
            WoN[i][k] = W[i][JCu][k];
         } else if (i < ICl && k< KCl) {
            WoS[i][k] = W[ICl][JCl][KCl];
            WoN[i][k] = W[ICl][JCu][KCl];
         } else if (i > ICu && k > KCu) {
            WoS[i][k] = W[ICu][JCl][KCu];
            WoN[i][k] = W[ICu][JCu][KCu];
         } else if (i<ICl && (k >= KCl && k <= KCu)){
            WoS[i][k] = W[ICl][JCl][k];
            WoN[i][k] = W[ICl][JCu][k];
         } else if (i>ICu && (k >= KCl && k <= KCu)){
            WoS[i][k] = W[ICu][JCl][k];
            WoN[i][k] = W[ICu][JCu][k];
         } else if ((i >= ICl && i <= ICu) && k< KCl) {
            WoS[i][k] = W[i][JCl][KCl];
            WoN[i][k] = W[i][JCu][KCl];
         } else if ((i >= ICl && i <= ICu) && k > KCu) {
            WoS[i][k] = W[i][JCl][KCu];
            WoN[i][k] = W[i][JCu][KCu];
         } else if (i < ICl  && k > KCu) {
            WoS[i][k] = W[ICl][JCl][KCu];
            WoN[i][k] = W[ICl][JCu][KCu];
         } else if (i >ICu  && k < KCl) {
            WoS[i][k] = W[ICu][JCl][KCl];
            WoN[i][k] = W[ICu][JCu][KCl];
         }
      } /* endfor */
   } /* endfor */

   for (  j = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for (  i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
          if ((j >= JCl && j <= JCu) && (i >= ICl && i <= ICu)) {
             WoT[i][j] = W[i][j][KCu];
             WoB[i][j] = W[i][j][KCl];
          } else if (i < ICl &&  j< JCl) {
             WoT[i][j] = W[ICl][JCl][KCu];
             WoB[i][j] = W[ICl][JCl][KCl];
          } else if(i > ICu &&  j > JCu) {
             WoT[i][j] = W[ICu][JCu][KCu];
             WoB[i][j] = W[ICu][JCu][KCl];
          }else if (i < ICl && (j >= JCl && j <= JCu)) {
             WoT[i][j] = W[ICl][j][KCu];
             WoB[i][j] = W[ICl][j][KCl];
          }else if (i > ICu && (j >= JCl && j <= JCu)) {
             WoT[i][j] = W[ICu][j][KCu];
             WoB[i][j] = W[ICu][j][KCl];
          } else if ((i >= ICl && i <= ICu) &&  j< JCl) {
             WoT[i][j] = W[i][JCl][KCu];
             WoB[i][j] = W[i][JCl][KCl];
          } else if ((i >= ICl && i <= ICu) &&  j> JCu) {
             WoT[i][j] = W[i][JCu][KCu];
             WoB[i][j] = W[i][JCu][KCl];
          } else if (i > ICu && j < JCl) {
             WoT[i][j] = W[ICu][JCl][KCu];
             WoB[i][j] = W[ICu][JCl][KCl];
          } else if (i < ICl && j > JCu) {
             WoT[i][j] = W[ICl][JCu][KCu];
             WoB[i][j] = W[ICl][JCu][KCl];
          }
       } /* endfor */
   } /* endfor */
      
   return (0);
    
}

/********************************************************
 * Routine: BCs                                         *
 *                                                      *
 * Apply boundary conditions at boundaries of the       *
 * specified hexa solution block.                       *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::BCs(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs) {
   
   int i, j, k;
   
   Vector3D MOVING_WALL_VELOCITY = IPs.Moving_wall_velocity, dX;
   double dpdx;
   
 
     
   for ( k =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k) 
      for ( j =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) {
         
         // Prescribe West boundary conditions.
         switch( Grid.BCtypeW[j][k]) {
            
         case BC_REFLECTION :
            W[ ICl-1][j][k] = SOLN_pSTATE::Reflect( W[ ICl][j][k],
                                                    Grid.nfaceW( ICl,j,k));
            U[ ICl-1][j][k] =W[ ICl-1][j][k].U();
            W[ ICl-2][j][k] = SOLN_pSTATE::Reflect( W[ ICl+1][j][k],
                                       Grid.nfaceW( ICl,j,k));
            U[ ICl-2][j][k] =  W[ ICl-2][j][k].U();
     
            
            break;

         case BC_FIXED_PRESSURE :
            W[ ICl-1][j][k] =  W[ ICl][j][k];
            W[ ICl-1][j][k].p = WoW[j][k].p;
            U[ ICl-1][j][k] = W[ ICl-1][j][k].U();
            W[ ICl-2][j][k] = W[ ICl][j][k] ;
            W[ ICl-2][j][k].p = WoW[j][k].p;
            U[ ICl-2][j][k] =  W[ ICl-2][j][k].U();
                   
            break;

         case BC_CONSTANT_EXTRAPOLATION :
            W[ ICl-1][j][k] =  W[ ICl][j][k];
            U[ ICl-1][j][k] = W[ ICl-1][j][k].U();
            W[ ICl-2][j][k] = W[ ICl][j][k] ;
            U[ ICl-2][j][k] =  W[ ICl-2][j][k].U();
          

            break;
            
         case BC_PERIODIC :  
            W[ICl-1][j][k] = W[ICu-1][j][k];
            U[ICl-1][j][k] = U[ICu-1][j][k];
            W[ICl-2][j][k] = W[ICu-2][j][k];
            U[ICl-2][j][k] = U[ICu-2][j][k];

            break;

         case BC_NO_SLIP :
            W[ICl-1][j][k] = SOLN_pSTATE::No_Slip(W[ICl][j][k],WoW[j][k], Grid.nfaceW(ICl,j,k),
                                                  IPs.Pressure_Gradient,
                                     FIXED_TEMPERATURE_WALL);
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = SOLN_pSTATE::No_Slip(W[ICl+1][j][k], WoW[j][k], Grid.nfaceW(ICl,j,k),
                                                  IPs.Pressure_Gradient,
                                                  FIXED_TEMPERATURE_WALL);
            U[ICl-2][j][k] = W[ICl-2][j][k].U();

        
            break;

         case BC_MOVING_WALL :
            W[ICl-1][j][k] = SOLN_pSTATE::Moving_Wall(W[ICl][j][k],WoW[j][k],
                                         Grid.nfaceW(ICl,j,k),
                                         MOVING_WALL_VELOCITY,
                                                      IPs.Pressure_Gradient,
                                         FIXED_TEMPERATURE_WALL);
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = SOLN_pSTATE::Moving_Wall(W[ICl+1][j][k], WoW[j][k], 
                                         Grid.nfaceW(ICl,j,k),
                                         MOVING_WALL_VELOCITY, 
                                                      IPs.Pressure_Gradient,
                                         FIXED_TEMPERATURE_WALL);
            U[ICl-2][j][k] = W[ICl-2][j][k].U();

            break;
        
          
         case BC_CHANNEL_INFLOW:
            dpdx = IPs.Pressure_Gradient.x; 
            //for turbulent channel flow
            // p linearly varys based on constant pressure gradient 
            dX = Grid.Cell[ICl-1][j][k].Xc - Grid.Cell[ICl][j][k].Xc; 
            W[ICl-1][j][k] = WoW[j][k];
            W[ICl-1][j][k].v.x = W[ICl][j][k].v.x;
            W[ICl-1][j][k].p = WoW[j][k].p  ;//-dpdx*dX.x;
            U[ICl-1][j][k] = W[ICl-1][j][k].U( );
            
            dX = Grid.Cell[ICl-2][j][k].Xc - Grid.Cell[ICl][j][k].Xc;
            W[ICl-2][j][k] = WoW[j][k];
            W[ICl-2][j][k].v.x = W[ICl][j][k].v.x;
            W[ICl-2][j][k].p =  W[ICl-1][j][k].p;// (WoW[j][k].p ) - dpdx*dX.x;
            U[ICl-2][j][k] = W[ICl-2][j][k].U( );
            break;  
            
         } /* endswitch */
         
         // Prescribe East boundary conditions.
         switch( Grid.BCtypeE[j][k]) {
         case BC_REFLECTION :
            W[ ICu+1][j][k] = SOLN_pSTATE::Reflect( W[ ICu][j][k],
                                       Grid.nfaceE( ICu,j,k));
            U[ ICu+1][j][k] = W[ ICu+1][j][k].U();
            W[ ICu+2][j][k] = SOLN_pSTATE::Reflect( W[ ICu-1][j][k],
                                       Grid.nfaceE( ICu,j,k));
            U[ ICu+2][j][k] = W[ ICu+2][j][k].U();
            

            break;
         case BC_FIXED_PRESSURE :
            W[ ICu+1][j][k] =  W[ ICu][j][k];
            W[ ICu+1][j][k].p = WoE[j][k].p;
            U[ ICu+1][j][k] = W[ ICu+1][j][k].U();
            
            W[ ICu+2][j][k] =  W[ ICu-1][j][k];
            W[ ICu+2][j][k].p = WoE[j][k].p; 
            U[ ICu+2][j][k] = W[ ICu+2][j][k].U();
          
           
             break;

         case BC_CONSTANT_EXTRAPOLATION :
            W[ ICu+1][j][k] =  W[ ICu][j][k];
            U[ ICu+1][j][k] = W[ ICu+1][j][k].U();
            W[ ICu+2][j][k] =  W[ ICu][j][k];
            U[ ICu+2][j][k] = W[ ICu+2][j][k].U();


            break;

         case BC_PERIODIC :
            W[ICu+1][j][k] = W[ICl+1][j][k];
            U[ICu+1][j][k] = U[ICl+1][j][k];
            W[ICu+2][j][k] = W[ICl+2][j][k];
            U[ICu+2][j][k] = U[ICl+2][j][k];


            break;

         case BC_NO_SLIP :
            W[ICu+1][j][k] = SOLN_pSTATE::No_Slip(W[ICu][j][k], WoE[j][k], Grid.nfaceE(ICu,j,k),
                                                  IPs.Pressure_Gradient,
                                     FIXED_TEMPERATURE_WALL);
            U[ICu+1][j][k] = W[ICu+1][j][k].U( );
            W[ICu+2][j][k] = SOLN_pSTATE::No_Slip(W[ICu-1][j][k], WoE[j][k], 
                                     Grid.nfaceE(ICu,j,k),
                                                  IPs.Pressure_Gradient,
                                     FIXED_TEMPERATURE_WALL);
            U[ICu+2][j][k] = W[ICu+2][j][k].U( );

            break;
            
         case BC_MOVING_WALL :
            W[ICu+1][j][k] = SOLN_pSTATE::Moving_Wall(W[ICu][j][k], WoE[j][k], Grid.nfaceE(ICu,j,k),
                                                      MOVING_WALL_VELOCITY,     IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
            U[ICu+1][j][k] = W[ICu+1][j][k].U( );
            W[ICu+2][j][k] = SOLN_pSTATE::Moving_Wall(W[ICu-1][j][k], WoE[j][k], 
                                     Grid.nfaceE(ICu,j,k),
                                                      MOVING_WALL_VELOCITY,
                                                      IPs.Pressure_Gradient,
                                     FIXED_TEMPERATURE_WALL);
            U[ICu+2][j][k] = W[ICu+2][j][k].U( );

            break;

         case BC_CHANNEL_OUTFLOW:
            dpdx = IPs.Pressure_Gradient.x; 
            // all constant extrapolation except pressure specified which linearly varys if there is pressure gradient
            dX = Grid.Cell[ICu+1][j][k].Xc - Grid.Cell[ICu][j][k].Xc; 
            W[ICu+1][j][k] =  W[ICu][j][k];
            W[ICu+1][j][k].p = (WoE[j][k].p);//-dpdx*dX.x;
            U[ICu+1][j][k] = W[ ICu+1][j][k].U();
            
            dX = Grid.Cell[ICu+2][j][k].Xc - Grid.Cell[ICu][j][k].Xc; 
            W[ICu+2][j][k] = W[ICu][j][k];
            W[ICu+2][j][k].p =  W[ICu+1][j][k].p;// (WoE[j][k].p)-dpdx*dX.x; 	
            U[ICu+2][j][k] = W[ICu+2][j][k].U( );
            
 
            break;
            
         }//endofeastface
         
      } /* endfor */
   
   for ( k =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k )
      for ( i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
              
         // Prescribe South boundary conditions.
         switch( Grid.BCtypeS[i][k]) {
         case BC_REFLECTION :
            W[i][ JCl-1][k] = SOLN_pSTATE::Reflect( W[i][ JCl][k],
                                       Grid.nfaceS(i,  JCl,k));
            U[i][ JCl-1][k] =  W[i][ JCl-1][k].U();
            W[i][ JCl-2][k] = SOLN_pSTATE::Reflect( W[i][ JCl+1][k],
                                       Grid.nfaceS(i,  JCl,k));
            U[i][ JCl-2][k] =  W[i][ JCl-2][k].U();

            break;
         
         case BC_FIXED_PRESSURE :
            W[i][ JCl-1][k] = W[i][ JCl][k];
            W[i][ JCl-1][k].p = WoS[i][k].p;
            U[i][ JCl-1][k] =  W[i][ JCl-1][k].U();
            W[i][ JCl-2][k] = W[i][ JCl+1][k];
            W[i][ JCl-2][k].p = WoS[i][k].p;
            U[i][ JCl-2][k] =  W[i][ JCl-2][k].U();
  
	    break;

         case BC_CONSTANT_EXTRAPOLATION :
            W[i][ JCl-1][k] = W[i][ JCl][k];
            U[i][ JCl-1][k] =  W[i][ JCl-1][k].U();
            W[i][ JCl-2][k] = W[i][ JCl][k];
            U[i][ JCl-2][k] =  W[i][ JCl-2][k].U();

	    break;

         case BC_PERIODIC :
            W[i][JCl-1][k] = W[i][JCu-1][k];
            U[i][JCl-1][k] = U[i][JCu-1][k];
            W[i][JCl-2][k] = W[i][JCu-2][k];
            U[i][JCl-2][k] = U[i][JCu-2][k];
           
            break;

         case BC_NO_SLIP :

            W[i][JCl-1][k] = SOLN_pSTATE::No_Slip(W[i][JCl][k], WoS[i][k], 
                                                  Grid.nfaceS(i, JCl,k),
                                                  IPs.Pressure_Gradient,
                                                  FIXED_TEMPERATURE_WALL);
            U[i][JCl-1][k] =  W[i][JCl-1][k].U( );
            W[i][JCl-2][k] = SOLN_pSTATE::No_Slip(W[i][JCl+1][k], WoS[i][k],
                                                  Grid.nfaceS(i, JCl,k),
                                                  IPs.Pressure_Gradient,
                                                  FIXED_TEMPERATURE_WALL);
            U[i][JCl-2][k] =  W[i][JCl-2][k].U( );

            break;
         case BC_MOVING_WALL :
            W[i][JCl-1][k] = SOLN_pSTATE::Moving_Wall(W[i][JCl][k], WoS[i][k], 
                                                      Grid.nfaceS(i, JCl,k),
                                                      MOVING_WALL_VELOCITY,
                                                      IPs.Pressure_Gradient, 
                                                      FIXED_TEMPERATURE_WALL);
            U[i][JCl-1][k] =  W[i][JCl-1][k].U( );
            W[i][JCl-2][k] = SOLN_pSTATE::Moving_Wall(W[i][JCl+1][k], WoS[i][k],
                                                      Grid.nfaceS(i, JCl,k),
                                                      MOVING_WALL_VELOCITY, 
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
            U[i][JCl-2][k] =  W[i][JCl-2][k].U( );
            
            break;
     
            
         } /* endswitch */
          
         // Prescribe North boundary conditions.
         switch( Grid.BCtypeN[i][k]) {
            
         case BC_REFLECTION :
            W[i][ JCu+1][k] = SOLN_pSTATE::Reflect( W[i][ JCu][k],
                                       Grid.nfaceN(i,  JCu,k));
            U[i][ JCu+1][k] =  W[i][ JCu+1][k].U();
            W[i][ JCu+2][k] = SOLN_pSTATE::Reflect( W[i][ JCu-1][k],
                                       Grid.nfaceN(i,  JCu, k));
            U[i][ JCu+2][k] =  W[i][ JCu+2][k].U();

            
            break;
            
         case BC_FIXED_PRESSURE :
            
            W[i][ JCu+1][k] = W[i][ JCu][k];
            W[i][ JCu+1][k].p = WoN[i][k].p;
            U[i][ JCu+1][k] =  W[i][ JCu+1][k].U();
          
            W[i][ JCu+2][k] =  W[i][ JCu-1][k];
            W[i][ JCu+2][k].p = WoN[i][k].p;
            U[i][ JCu+2][k] =  W[i][ JCu+2][k].U();

            break;
            
         case BC_CONSTANT_EXTRAPOLATION :
            W[i][ JCu+1][k] = W[i][ JCu][k];
            U[i][ JCu+1][k] =  W[i][ JCu+1][k].U();
            W[i][ JCu+2][k] =  W[i][ JCu][k];
            U[i][ JCu+2][k] =  W[i][ JCu+2][k].U();

            break;

         case BC_PERIODIC :
            W[i][JCu+1][k] = W[i][JCl+1][k];
            U[i][JCu+1][k] = U[i][JCl+1][k];
            W[i][JCu+2][k] = W[i][JCl+2][k];
            U[i][JCu+2][k] = U[i][JCl+2][k];
          
           
            break;

         case BC_NO_SLIP :
            W[i][JCu+1][k] = SOLN_pSTATE::No_Slip(W[i][JCu][k], WoN[i][k],
                                     Grid.nfaceN(i, JCu, k),
                                                  IPs.Pressure_Gradient,
                                                  FIXED_TEMPERATURE_WALL);
            U[i][JCu+1][k] = W[i][JCu+1][k].U( );
            W[i][JCu+2][k] = SOLN_pSTATE::No_Slip(W[i][JCu-1][k], WoN[i][k],
                                     Grid.nfaceN(i, JCu, k),
                                                  IPs.Pressure_Gradient,
                                     FIXED_TEMPERATURE_WALL);
            U[i][JCu+2][k] = W[i][JCu+2][k].U( );


         break;

         case BC_MOVING_WALL :
            W[i][JCu+1][k] = SOLN_pSTATE::Moving_Wall(W[i][JCu][k], WoN[i][k],
                                         Grid.nfaceN(i, JCu, k),
                                         MOVING_WALL_VELOCITY,
                                                      IPs.Pressure_Gradient,
                                         FIXED_TEMPERATURE_WALL);
            U[i][JCu+1][k] = W[i][JCu+1][k].U( );
            W[i][JCu+2][k] = SOLN_pSTATE::Moving_Wall(W[i][JCu-1][k], WoN[i][k],
                                         Grid.nfaceN(i, JCu, k),
                                         MOVING_WALL_VELOCITY,
                                                      IPs.Pressure_Gradient,
                                         FIXED_TEMPERATURE_WALL);
            U[i][JCu+2][k] = W[i][JCu+2][k].U( );

         break;

         } /* endswitch */
      } /* endfor */
   
   
   for ( j =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j )
      for ( i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
         
         // Prescribe Bottom boundary conditions.
         switch( Grid.BCtypeB[i][j]) {
         case BC_REFLECTION :
            W[i][j][ KCl-1 ] = SOLN_pSTATE::Reflect( W[i][j][ KCl],
                                        Grid.nfaceBot(i, j,  KCl));
            U[i][j][ KCl-1] =  W[i][j][ KCl-1].U();
            W[i][j][ KCl-2] = SOLN_pSTATE::Reflect( W[i][j][ KCl+1],
                                       Grid.nfaceBot(i,j,  KCl));
            U[i][j][ KCl-2] = W[i][j][ KCl-2].U();

            
            break;
            
         case BC_FIXED_PRESSURE :
            W[i][j][ KCl-1 ] = W[i][j][ KCl];
            W[i][j][ KCl-1 ].p = WoB[i][j].p;
            U[i][j][ KCl-1] =  W[i][j][ KCl-1].U();
            W[i][j][ KCl-2] =  W[i][j][ KCl+1];
            W[i][j][ KCl-2].p =  WoB[i][j].p;
            U[i][j][ KCl-2] = W[i][j][ KCl-2].U();

            break;

         case BC_CONSTANT_EXTRAPOLATION :
            W[i][j][ KCl-1 ] = W[i][j][ KCl];
            U[i][j][ KCl-1] =  W[i][j][ KCl-1].U();
            W[i][j][ KCl-2] =  W[i][j][ KCl];
            U[i][j][ KCl-2] = W[i][j][ KCl-2].U();

            break;

         case BC_PERIODIC :
            W[i][j][KCl- 1] = W[i][j][KCu -1];
            U[i][j][KCl- 1] = U[i][j][KCu -1];
            W[i][j][KCl- 2] = W[i][j][KCu -2];
            U[i][j][KCl -2] = U[i][j][KCu -2];

            break;

         case BC_NO_SLIP :

            W[i][j][KCl-1] = SOLN_pSTATE::No_Slip(W[i][j][KCl], WoB[i][j],
                                     Grid.nfaceBot(i, j, KCl),
                                                  IPs.Pressure_Gradient,
                                     FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-1] =  W[i][j][KCl-1].U( ); 
            W[i][j][KCl-2] =  SOLN_pSTATE::No_Slip(W[i][j][KCl+1], WoB[i][j],
                                     Grid.nfaceBot(i, j, KCl),
                                                   IPs.Pressure_Gradient,
                                     FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-2] =  W[i][j][KCl-2].U( );
            
            break;
            
            
         case BC_MOVING_WALL :
            W[i][j][KCl-1] =  SOLN_pSTATE::Moving_Wall(W[i][j][KCl], WoB[i][j],
                                         Grid.nfaceBot(i, j,KCl),
                                         MOVING_WALL_VELOCITY,
                                                       IPs.Pressure_Gradient,
                                         FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-1] =  W[i][j][KCl-1].U( );
            W[i][j][KCl-2] =  SOLN_pSTATE::Moving_Wall(W[i][j][KCl+1], WoB[i][j],
                                         Grid.nfaceBot(i, j ,KCl),
                                         MOVING_WALL_VELOCITY,
                                                       IPs.Pressure_Gradient,
                                         FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-2] =  W[i][j][KCl-2].U( );

            
            break;

            
         } /* endswitch */
         
         // Prescribe Top boundary conditions.
         switch( Grid.BCtypeT[i][j]) {
         case BC_REFLECTION :
            W[i][j][ KCu+1] =  SOLN_pSTATE::Reflect( W[i][j][ KCu],
                                       Grid.nfaceTop(i,j,  KCu));
            U[i][j][ KCu+1] =  W[i][j][ KCu+1].U();
            W[i][j][ KCu+2] =  SOLN_pSTATE::Reflect( W[i][j][ KCu-1],
                                       Grid.nfaceTop(i, j,  KCu));
            U[i][j][ KCu+2] =  W[i][j][ KCu+2].U();

                        
            break;
         case BC_FIXED_PRESSURE :
            W[i][j][ KCu+1] = W[i][j][ KCu];
            W[i][j][ KCu+1].p = WoT[i][j].p;
            U[i][j][ KCu+1] =  W[i][j][ KCu+1].U();
            W[i][j][ KCu+2] =  W[i][j][ KCu-1];
            W[i][j][ KCu+2].p = WoT[i][j].p;
            U[i][j][ KCu+2] =  W[i][j][ KCu+2].U();

            break;

         case BC_CONSTANT_EXTRAPOLATION :
            W[i][j][ KCu+1] = W[i][j][ KCu];
            U[i][j][ KCu+1] =  W[i][j][ KCu+1].U();
            W[i][j][ KCu+2] =  W[i][j][ KCu];
            U[i][j][ KCu+2] =  W[i][j][ KCu+2].U();

            break;

         case BC_PERIODIC :
            W[i][j][KCu + 1] = W[i][j][KCl +1];
            U[i][j][KCu + 1] = U[i][j][KCl +1];
            W[i][j][KCu + 2] = W[i][j][KCl +2];
            U[i][j][KCu + 2] = U[i][j][KCl +2];

            break;

         case BC_NO_SLIP :
            W[i][j][KCu +1] =  SOLN_pSTATE::No_Slip(W[i][j][KCu], WoT[i][j],
                                      Grid.nfaceTop(i, j , KCu),
                                                    IPs.Pressure_Gradient,
                                      FIXED_TEMPERATURE_WALL);
            U[i][j][KCu +1] =  W[i][j][KCu+1].U( );
            W[i][j][KCu +2] =  SOLN_pSTATE::No_Slip(W[i][j][KCu -1], WoT[i][j],
                                      Grid.nfaceTop(i, j , KCu),
                                                    IPs.Pressure_Gradient, 
                                     FIXED_TEMPERATURE_WALL);
            U[i][j][KCu +2] =  W[i][j][KCu +2].U( );

            
            break;
            
            
         case BC_MOVING_WALL :
            W[i][j][KCu +1] =  SOLN_pSTATE::Moving_Wall(W[i][j][KCu], WoT[i][j],
                                          Grid.nfaceTop(i, j , KCu),
                                          MOVING_WALL_VELOCITY,
                                                        IPs.Pressure_Gradient,
                                          FIXED_TEMPERATURE_WALL);
            
            U[i][j][KCu +1] =  W[i][j][KCu+1].U( );

            W[i][j][KCu +2] =  SOLN_pSTATE::Moving_Wall(W[i][j][KCu -1], WoT[i][j],
                                          Grid.nfaceTop(i, j , KCu),
                                          MOVING_WALL_VELOCITY,
                                                        IPs.Pressure_Gradient,
                                          FIXED_TEMPERATURE_WALL);
            U[i][j][KCu +2] =  W[i][j][KCu +2].U( );


            break;
            
         } /* endswitch */
      } /* endfor */

}

/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Euler time stepping scheme) for the    *
 * specified hexa solution block according to           *
 * the Courant-Friedrichs-Lewy condition.               *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
double Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::CFL(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs) {
   
   int i, j, k;
   double dtMin, d_i, d_j, d_k, v_i, v_j, v_k, a, dt_vis, nv, dt_chem;
   double mr, aa_i, aa_j, aa_k;
   
   dtMin = MILLION;
   
   for ( k  =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k ) 
      for ( j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) 
         for ( i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
            if (i <  ICl || i >  ICu ||
                j <  JCl || j >  JCu || 
                k <  KCl || k >  KCu) {
               dt[i][j][k] = ZERO;
            } else {
               d_i = TWO*(Grid.volume(i,j,k)/
                          (Grid.AfaceE(i, j, k)+ Grid.AfaceW(i, j, k)));
               d_j = TWO*( Grid.volume(i, j, k)/
                           (Grid.AfaceN(i, j, k)+ Grid.AfaceS(i, j, k)));
               d_k = TWO*( Grid.volume(i, j, k)/
                           (Grid.AfaceTop(i, j, k)+ Grid.AfaceBot(i, j, k)));
               v_i = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                           (Grid.nfaceE(i, j, k)- Grid.nfaceW(i, j, k)));
               v_j = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                           ( Grid.nfaceN(i, j, k)- Grid.nfaceS(i, j, k)));
               v_k = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                            (Grid.nfaceTop(i, j, k)- Grid.nfaceBot(i, j, k)));
               a =  W[i][j][k].a();
               dt[i][j][k] = min(min(d_i/(a+fabs(v_i)), d_j/(a+fabs(v_j))),
                                 d_k/(a+fabs(v_k)));
               
               if (IPs.i_Flow_Type != FLOWTYPE_INVISCID) {  
                  nv = W[i][j][k].mu()/W[i][j][k].rho;
                  dt_vis = min(min((d_i*d_i)/(3.0*nv), (d_j*d_j)/(3.0*nv)), (d_k*d_k)/(3.0*nv)); 
                  dt[i][j][k]  = min(dt_vis, dt[i][j][k]);
               }

               /******** Chemical Source Term deltat calculation ************/   
               if (W[i][j][k].React.reactset_flag != NO_REACTIONS){
                  dt_chem = HALF/W[i][j][k].dSwdU_max_diagonal();
                  dt[i][j][k] = min(dt_chem, dt[i][j][k]);
               } /* endif */
               
               dtMin = min(dtMin,  dt[i][j][k]);
            } /* endif */
         } /* endfor */
   
   for ( k  =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k ) 
      for ( j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) 
         for ( i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
            if (i <  ICl || i >  ICu ||
                j <  JCl || j >  JCu || 
                k <  KCl || k >  KCu) {
               dt[i][j][k] = dtMin;
            } /* endif */
         } /* endfor */
   
   /* Return the global time step. */
   
   return (dtMin);
}


/********************************************************
 * Routine: Set_Global_TimeStep                         *
 *                                                      *
 * Assigns global time step to specified solution block *
 * for time-accurate calculations.                      *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Set_Global_TimeStep(const double &Dt_min) {
   
   int i, j, k;

   for ( k  = KCl ; k <= KCu ; ++k ) {
      for ( j  = JCl ; j <= JCu ; ++j ) {
         for ( i = ICl ; i <= ICu ; ++i ) {
            dt[i][j][k] = Dt_min;
         } /* endfor */
      } /* endfor */
   } /*endfor */

}

/********************************************************
 * Routine: L1_Norm_Residual                            *
 *                                                      *
 * Determines the L1-norm of the solution residual for  *
 * the specified hexahedral solution block.             *
 * Useful for monitoring convergence of the solution    *
 * for steady state problems.                           *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
double Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::L1_Norm_Residual(void) {
   
    int i, j, k;
    double l1_norm;

    l1_norm = ZERO;
    for ( k  =  KCl ; k <=  KCu ; ++k ) 
       for ( j  =  JCl ; j <=  JCu ; ++j ) 
          for ( i =  ICl ; i <=  ICu ; ++i ) {
             l1_norm += fabs( abs(dUdt[i][j][k][0].rhov));
          } /* endfor */

    return (l1_norm);
    
}

/********************************************************
 * Routine: L2_Norm_Residual                            *
 *                                                      *
 * Determines the L2-norm of the solution residual for  *
 * the specified hexahedral solution block.             *
 * Useful for monitoring convergence of the solution    *
 * for steady state problems.                           *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
double Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::L2_Norm_Residual(void) {

   int i, j, k;
   double l2_norm;
   
   l2_norm = ZERO;
   
   for ( k  =  KCl ; k <=  KCu ; ++k ) 
      for ( j  =  JCl ; j <=  JCu ; ++j ) 
         for ( i =  ICl ; i <=  ICu ; ++i ) {
            l2_norm += sqr( abs(dUdt[i][j][k][0].rhov));
         } /* endfor */
   
   l2_norm = sqrt(l2_norm);
   
   return (l2_norm);
   
}

/********************************************************
 * Routine: Max_Norm_Residual                           *
 *                                                      *
 * Determines the maximum norm of the solution residual *
 * for the specified hexahedral solution block.         *
 * Useful for monitoring convergence of the solution    *
 * for steady state problems.                           *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
double Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Max_Norm_Residual(void) {

   int i, j, k;
   double max_norm;

   max_norm = ZERO;
   for ( k  =  KCl ; k <=  KCu ; ++k ) 
      for ( j  =  JCl ; j <=  JCu ; ++j ) 
         for ( i =  ICl ; i <=  ICu ; ++i ) {
            max_norm = max(max_norm, fabs( abs(dUdt[i][j][k][0].rhov)));
            
         } /* endfor */
   
   return (max_norm);
   
}

/********************************************************
 * Routine: WtoU                                        *
 *                                                      *
 * Convert primitive solution vector to conservative    *
 * solution vector.                                     *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::WtoU(void) {
   for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) 
      for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) 
         for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
            U[i][j][k]= W[i][j][k].U();
         }// convert W to U 
   
   return (0);
   
}

/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares          *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within each cell of the        *
 * computational mesh of the specified Hexahedral       *
 * solution block.  A least squares approach is         *
 * used in the evaluation of the unlimited solution     *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Linear_Reconstruction_LeastSquares(const int Limiter) {
   
   int i, j, k;
   
   /* Carry out the limited solution reconstruction in
      each cell of the computational mesh. */

   for ( k  =  KCl- Nghost+1 ; k <=  KCu+ Nghost-1 ; ++k ) 
      for ( j  =  JCl- Nghost+1 ; j <=  JCu+ Nghost-1 ; ++j ) 
         for ( i =  ICl- Nghost+1 ; i <=  ICu+ Nghost-1 ; ++i ) {
            Linear_Reconstruction_LeastSquares(i, j, k, Limiter);
         } /* endfor */
   
}

/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares          *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within a given cell (i,j,k) of *
 * the computational mesh for the specified             *
 * Hexahedral solution block.  A least squares          *
 * approach is used in the evaluation of the unlimited  *
 * solution gradients.  Several slope limiters may be   *
 * used.                                                *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Linear_Reconstruction_LeastSquares(const int i, 
									      const int j,
									      const int k,
									      const int Limiter) {
   
   int n, n2, n_pts, i_index[26], j_index[26], k_index[26];
   double u0Min, u0Max, uHexa[6], PHI;
   double DxDx_ave, DxDy_ave, DyDy_ave, DxDz_ave, DyDz_ave, DzDz_ave;
   double D;;
   Vector3D dX;
   
   /* solnvec in  DU (DUDx_ave, DUDy_ave, DUDz_ave, D1, D2, D3) 
      is allocated using new  */ 
   SOLN_pSTATE DU, DUDx_ave, DUDy_ave, DUDz_ave;
   SOLN_pSTATE D1, D2, D3;
   SOLN_pSTATE Temp1, Temp2, Temp3;
   
   int NUM_VAR_3D = NumVar();
   
   Vector3D dX_neigbor;

   int num=0;
   
   if (i ==  ICl- Nghost || i ==  ICu+ Nghost ||
       j ==  JCl- Nghost || j ==  JCu+ Nghost || 
       k ==  KCl- Nghost || k ==  KCu+ Nghost) {
      n_pts = 0;
   } else {
      n_pts = 26;
      // k plane
      i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
      i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
      i_index[2] = i+1; j_index[2] = j-1; k_index[2] = k;
      i_index[3] = i-1; j_index[3] = j  ; k_index[3] = k;
      i_index[4] = i+1; j_index[4] = j  ; k_index[4] = k;
      i_index[5] = i-1; j_index[5] = j+1; k_index[5] = k;
      i_index[6] = i  ; j_index[6] = j+1; k_index[6] = k;
      i_index[7] = i+1; j_index[7] = j+1; k_index[7] = k;
      //k-1 plane
      i_index[8] = i-1; j_index[8] = j-1; k_index[8] = k-1;
      i_index[9] = i  ; j_index[9] = j-1; k_index[9] = k-1;
      i_index[10] = i+1; j_index[10] = j-1; k_index[10] = k-1;
      i_index[11] = i-1; j_index[11] = j  ; k_index[11] = k-1;
      i_index[12] = i  ; j_index[12] = j  ; k_index[12] = k-1;
      i_index[13] = i+1; j_index[13] = j  ; k_index[13] = k-1;
      i_index[14] = i-1; j_index[14] = j+1; k_index[14] = k-1;
      i_index[15] = i  ; j_index[15] = j+1; k_index[15] = k-1;
      i_index[16] = i+1; j_index[16] = j+1; k_index[16] = k-1;
      //k+1 plane
      i_index[17] = i-1; j_index[17] = j-1; k_index[17] = k+1;
      i_index[18] = i  ; j_index[18] = j-1; k_index[18] = k+1;
      i_index[19] = i+1; j_index[19] = j-1; k_index[19] = k+1;
      i_index[20] = i-1; j_index[20] = j  ; k_index[20] = k+1;
      i_index[21] = i  ; j_index[21] = j  ; k_index[21] = k+1;
      i_index[22] = i+1; j_index[22] = j  ; k_index[22] = k+1;
      i_index[23] = i-1; j_index[23] = j+1; k_index[23] = k+1;
      i_index[24] = i  ; j_index[24] = j+1; k_index[24] = k+1;
      i_index[25] = i+1; j_index[25] = j+1; k_index[25] = k+1;
   } /* endif */
     
   if (n_pts > 0) {
      DUDx_ave.Vacuum();
      DUDy_ave.Vacuum();
      DUDz_ave.Vacuum();
      D1.Vacuum();
      D2.Vacuum();
      D3.Vacuum();
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DxDz_ave = ZERO;
      DyDy_ave = ZERO;
      DyDz_ave = ZERO;
      DzDz_ave = ZERO;
      D = ZERO;
      
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
         dX =  Grid.Cell[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ].Xc - 
            Grid.Cell[i][j][k].Xc;
         DU =  W[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ] -  W[i][j][k];
         
         DUDx_ave += DU*dX.x;
         DUDy_ave += DU*dX.y;
         DUDz_ave += DU*dX.z;
         DxDx_ave += dX.x*dX.x;
         DxDy_ave += dX.x*dX.y;
         DxDz_ave += dX.x*dX.z;
         DyDy_ave += dX.y*dX.y;
         DyDz_ave += dX.y*dX.z;
         DzDz_ave += dX.z*dX.z;
         
      } /* endfor */
      
      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DUDz_ave = DUDz_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DxDz_ave = DxDz_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);
      DyDz_ave = DyDz_ave/double(n_pts);
      DzDz_ave = DzDz_ave/double(n_pts);
      
      // (1) Either write a linear solver for 3x3 linear system
      // (2) Or simplely use cramer's rule for this simple system

      D = DxDx_ave*(DyDy_ave* DzDz_ave - DyDz_ave*DyDz_ave) + 
          DxDy_ave*(DxDz_ave*DyDz_ave - DxDy_ave*DzDz_ave)+
          DxDz_ave*(DxDy_ave*DyDz_ave - DxDz_ave*DyDy_ave);
      
      D1 = DUDx_ave*(DyDy_ave* DzDz_ave - DyDz_ave*DyDz_ave) + 
           DUDy_ave*(DxDz_ave*DyDz_ave - DxDy_ave*DzDz_ave)+
           DUDz_ave*(DxDy_ave*DyDz_ave - DxDz_ave*DyDy_ave);
      
      D2 =DxDx_ave*(DUDy_ave* DzDz_ave - DUDz_ave*DyDz_ave) + 
          DxDy_ave*(DxDz_ave*DUDz_ave - DUDx_ave*DzDz_ave)+
          DxDz_ave*(DUDx_ave*DyDz_ave - DxDz_ave*DUDy_ave);

      D3 =DxDx_ave*(DyDy_ave* DUDz_ave - DyDz_ave*DUDy_ave) + 
          DxDy_ave*(DUDx_ave*DyDz_ave - DxDy_ave*DUDz_ave)+
          DxDz_ave*(DxDy_ave*DUDy_ave - DUDx_ave*DyDy_ave);

      dWdx[i][j][k] = D1/D;
      dWdy[i][j][k] = D2/D;
      dWdz[i][j][k] = D3/D;
      
      if (! Freeze_Limiter) {
         for ( n = 1 ; n <= NUM_VAR_3D ; ++n ) {
            
            u0Min =  W[i][j][k][n];
            u0Max = u0Min;
            for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
               u0Min = min(u0Min,  
                           W[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ][n]);
               u0Max = max(u0Max,  
                           W[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ][n]);
            } /* endfor */
            
            dX =  Grid.xfaceE(i, j, k)- Grid.Cell[i][j][k].Xc;
            uHexa[0] =  W[i][j][k][n] + 
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
            dX =  Grid.xfaceW(i, j, k)- Grid.Cell[i][j][k].Xc;
            uHexa[1] =  W[i][j][k][n] + 
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
            dX =  Grid.xfaceN(i, j, k)- Grid.Cell[i][j][k].Xc;
            uHexa[2] =  W[i][j][k][n] + 
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
	    dX =  Grid.xfaceS(i, j, k)- Grid.Cell[i][j][k].Xc;
	    uHexa[3] =  W[i][j][k][n] + 
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
            dX =  Grid.xfaceTop(i, j, k)- Grid.Cell[i][j][k].Xc;
	    uHexa[4] =  W[i][j][k][n] + 
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
            dX =  Grid.xfaceBot(i, j, k)- Grid.Cell[i][j][k].Xc;
	    uHexa[5] =  W[i][j][k][n] + 
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
	    
	    switch(Limiter) {
	    case LIMITER_ONE :
               PHI = ONE;
               break;
	    case LIMITER_ZERO :
               PHI = ZERO;
               break;
	    case LIMITER_BARTH_JESPERSEN :
               PHI = Limiter_BarthJespersen(uHexa,  W[i][j][k][n], 
                                            u0Min, u0Max, 6);
               break;
	    case LIMITER_VENKATAKRISHNAN :
               PHI = Limiter_Venkatakrishnan(uHexa,  W[i][j][k][n], 
                                             u0Min, u0Max, 6);
               break;
	    case LIMITER_VANLEER :
               PHI = Limiter_VanLeer(uHexa,  W[i][j][k][n], 
                                     u0Min, u0Max, 6);
               break;
	    case LIMITER_VANALBADA :
               PHI = Limiter_VanAlbada(uHexa,  W[i][j][k][n], 
                                       u0Min, u0Max, 6);
               break;
	    default:
               PHI = Limiter_BarthJespersen(uHexa,  W[i][j][k][n], 
                                            u0Min, u0Max, 6);
               break;
 	    }// endswitch
	    
	     phi[i][j][k][n] = PHI;
                        
         } /* endfor */
      } /* endif */
   } else {
       dWdx[i][j][k].Vacuum();
       dWdy[i][j][k].Vacuum();
       dWdz[i][j][k].Vacuum();
       phi[i][j][k].Vacuum(); 

   } /* endif */
   
    
}

/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine determines the solution residuals for a *
 * given stage of a variety of multi-stage explicit     *
 * time integration schemes for a given solution block. *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::dUdt_Multistage_Explicit(const int i_stage,
                                                                   Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs) {

   int i, j, k,  k_residual;
   double omega; 
   Vector3D dX;
   
   SOLN_pSTATE Wl, Wr;
   SOLN_cSTATE Flux;

   SOLN_cSTATE U_VACUUM;
   U_VACUUM.Vacuum();
   SOLN_pSTATE W_VACUUM;
   W_VACUUM.Vacuum();
 

   /* Evaluate the solution residual for stage 
    *       i_stage of an N stage scheme. */
   /* Evaluate the time step fraction and residual storage location 
    *       for the stage. */
   switch(IPs.i_Time_Integration) {
   case TIME_STEPPING_EXPLICIT_EULER :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      
      break;
   case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      break;
   case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0; 
      if (IPs.N_Stage == 4) { 
         if (i_stage == 4) { 
            k_residual = 0; 
         } else { 
            k_residual = i_stage - 1; 
         } // endif  
      } // endif 
      break;
   case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
      omega = MultiStage_Optimally_Smoothing(i_stage, 
                                             IPs.N_Stage,
                                             IPs.i_Limiter);
      k_residual = 0;
      break;
   default: 
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      break;
   } /* endswitch */

   /* Perform the linear reconstruction within each cell
      of the computational grid for this stage. */
   
   switch(IPs.i_Reconstruction) {
      
   case RECONSTRUCTION_LEAST_SQUARES :
      
      Linear_Reconstruction_LeastSquares(IPs.i_Limiter);
      
      break;
   default:
      Linear_Reconstruction_LeastSquares(IPs.i_Limiter);
      break;
   } /* endswitch */

   /* Evaluate the time rate of change of the solution
      (i.e., the solution residuals) using a second-order
      limited upwind scheme with a variety of flux functions. */
 
   // Add i-direction (zeta-direction) fluxes.
   for ( k  =  KCl-1 ; k <=  KCu+1 ; ++k ){
      for ( j  =  JCl-1 ; j <=  JCu+1 ; ++j ) {
         if ( i_stage == 1 ) {
            Uo[ ICl-1][j][k] =  U[ ICl-1][j][k];
            dUdt[ ICl-1][j][k][k_residual] = U_VACUUM;
            
            
         } else {
            dUdt[ ICl-1][j][k][k_residual] = U_VACUUM;
         } /* endif */
         
         for ( i =  ICl-1 ; i <=  ICu ; ++i ) {
            if ( i_stage == 1 ) {
               Uo[i+1][j][k] =  U[i+1][j][k];
               
               dUdt[i+1][j][k][k_residual] = U_VACUUM;
            } else if (( j >  JCl-1 && j <  JCu+1 )
                       && (k >  KCl-1 && k <  KCu+1) ){
               switch(IPs.i_Time_Integration) {

               case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
                  //dUdt[i+1][j][k][k_residual] =
                  //  dUdt[i+1][j][k][k_residual];
                  break;
               case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
		  if (IPs.N_Stage == 2) {
                     //dUdt[i+1][j][k_residual] =
                     // dUdt[i+1][j][k_residual];
                  } else if (IPs.N_Stage == 4 && i_stage == 4) {
                     dUdt[i+1][j][k][k_residual] =
                        dUdt[i+1][j][k][0] +
                        TWO*dUdt[i+1][j][k][1] +TWO*dUdt[i+1][j][k][2];
                  } else {
         	     dUdt[i+1][j][k][k_residual].Vacuum();
                  } /* endif */
                  break;
               case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
                  dUdt[i+1][j][k][k_residual] = U_VACUUM;
                  break;
               default:
                  dUdt[i+1][j][k][k_residual] = U_VACUUM;
                  
                  break;
               } /* endswitch */
            } /* endif */
            if (( j >  JCl-1 && j <  JCu+1 ) &&
                ( k >  KCl-1 && k <  KCu+1 )) {
               
               /* Evaluate the cell interface i-direction fluxes. */
               if (i ==  ICl-1 &&  (Grid.BCtypeW[j][k] == BC_REFLECTION ||
                                    Grid.BCtypeW[j][k] == BC_NO_SLIP ||
                                    Grid.BCtypeW[j][k] == BC_MOVING_WALL) ) {
                  
                  dX =  Grid.xfaceW(i+1, j, k)- Grid.Cell[i+1][j][k].Xc;
                  
                  Wr =  W[i+1][j][k] +
                     ( phi[i+1][j][k]^dWdx[i+1][j][k])*dX.x +
                     ( phi[i+1][j][k]^dWdy[i+1][j][k])*dX.y +
                     ( phi[i+1][j][k]^dWdz[i+1][j][k])*dX.z;
                  
                  if ( Grid.BCtypeW[j][k] == BC_REFLECTION) {
                     Wl =  SOLN_pSTATE::Reflect(Wr,  Grid.nfaceW(i+1, j, k));
                  }
                  if ( Grid.BCtypeW[j][k] == BC_NO_SLIP) {
                     Wl =  SOLN_pSTATE::No_Slip(Wr, WoW[j][k], Grid.nfaceW(i+1, j, k),
                                                IPs.Pressure_Gradient,
                                                FIXED_TEMPERATURE_WALL);
                  }
                  if ( Grid.BCtypeW[j][k] == BC_MOVING_WALL) {
                     Wl =  SOLN_pSTATE::Moving_Wall(Wr, WoW[j][k], Grid.nfaceW(i+1, j, k),
                                                    IPs.Moving_wall_velocity,
                                                    IPs.Pressure_Gradient,
                                                    FIXED_TEMPERATURE_WALL);
                  }
               } else if (i ==  ICu &&
                          ( Grid.BCtypeE[j][k] == BC_REFLECTION ||
                            Grid.BCtypeE[j][k] == BC_NO_SLIP ||
                            Grid.BCtypeE[j][k] == BC_MOVING_WALL )) {
                  
                  dX =  Grid.xfaceE(i, j, k)- Grid.Cell[i][j][k].Xc;
                  Wl =  W[i][j][k] +
                     ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                     ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                     ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                  if ( Grid.BCtypeE[j][k] == BC_REFLECTION) {
                     Wr =  SOLN_pSTATE::Reflect(Wl,  Grid.nfaceE(i, j, k));
                  }
                  if ( Grid.BCtypeE[j][k] == BC_NO_SLIP) {
                     Wr =  SOLN_pSTATE::No_Slip(Wl, WoE[j][k], Grid.nfaceE(i, j, k),
                                                IPs.Pressure_Gradient,
                                                FIXED_TEMPERATURE_WALL);
                  }
                  if ( Grid.BCtypeE[j][k] == BC_MOVING_WALL) {
                     Wr =  SOLN_pSTATE::Moving_Wall(Wl, WoE[j][k], Grid.nfaceE(i, j, k),
                                                    IPs.Moving_wall_velocity,
                                                    IPs.Pressure_Gradient,
                                                    FIXED_TEMPERATURE_WALL);
                  }
                  
               } else {
                  dX =  Grid.xfaceE(i, j, k)- Grid.Cell[i][j][k].Xc;
                  Wl =  W[i][j][k] +
                     ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                     ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                     ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                  dX =  Grid.xfaceW(i+1, j, k)- Grid.Cell[i+1][j][k].Xc;
                  Wr =  W[i+1][j][k] +
                     ( phi[i+1][j][k]^ dWdx[i+1][j][k])*dX.x +
                     ( phi[i+1][j][k]^ dWdy[i+1][j][k])*dX.y +
                     ( phi[i+1][j][k]^ dWdz[i+1][j][k])*dX.z;                  
               } /* endif */
               
               
               switch(IPs.i_Flux_Function) {
               case FLUX_FUNCTION_HLLE :
                  
                  Flux =  SOLN_pSTATE::FluxHLLE_n(Wl, Wr,  Grid.nfaceE(i, j, k));
                  
                  break;
                  
               case FLUX_FUNCTION_ROE :
                  
                  Flux =  SOLN_pSTATE::FluxRoe_n(Wl, Wr,  Grid.nfaceE(i, j, k));
                  break;
                  
               } /* endswitch */
               
               /* Evaluate cell-averaged solution changes. */
               
               dUdt[i][j][k][k_residual] -=
                  (IPs.CFL_Number* dt[i][j][k])*
                  Flux* Grid.AfaceE(i, j, k)/
                  Grid.volume(i, j, k);
               dUdt[i+1][j][k][k_residual] +=
                  (IPs.CFL_Number* dt[i+1][j][k])*
                  Flux* Grid.AfaceW(i+1, j, k)/
                  Grid.volume(i+1, j ,k);
               
               
               /* Include source terms associated with the finite-rate chemistry and
                  turbulence/chemistry interactions */
               
               if (W[i][j][k].React.reactset_flag != NO_REACTIONS) {
                  
                  dUdt[i][j][k][k_residual] += IPs.CFL_Number*dt[i][j][k]*W[i][j][k].Sw(
                     W[i][j][k].React.reactset_flag);
                  
               } /* endif */
              
               
               /* Save west and east face boundary flux. */
               
               //    if (i ==  ICl-1) {
//                     FluxW[j] = -Flux* Grid.lfaceW(i+1, j);
//                 } else if (i ==  ICu) {
//                     FluxE[j] = Flux* Grid.lfaceE(i, j);
//                 } /* endif */
               
            } /* endif */
         } /* endfor */
         
         if (( j >  JCl-1 && j <  JCu+1 ) &&
             ( k >  KCl-1 && k <  KCu+1 ) ) {
            dUdt[ ICl-1][j][k][k_residual] = U_VACUUM;
            dUdt[ ICu+1][j][k][k_residual] = U_VACUUM;
         } /* endif */
      } /* endfor */
   } /* endfor */

   // Add j-direction (eta-direction) fluxes.
   for ( k =  KCl ; k <=  KCu ; ++k ){
      for ( i =  ICl ; i <=  ICu ; ++i ) {
         for ( j  =  JCl-1 ; j <=  JCu ; ++j ) {
            
            /* Evaluate the cell interface j-direction fluxes. */
            if (j ==  JCl-1 &&
                ( Grid.BCtypeS[i][k] == BC_REFLECTION ||
                  Grid.BCtypeS[i][k] == BC_NO_SLIP||
                  Grid.BCtypeS[i][k] == BC_MOVING_WALL)) {
               
               dX =  Grid.xfaceS(i, j+1, k)- Grid.Cell[i][j+1][k].Xc;
               Wr =  W[i][j+1][k] +
                  ( phi[i][j+1][k]^ dWdx[i][j+1][k])*dX.x +
                  ( phi[i][j+1][k]^ dWdy[i][j+1][k])*dX.y+
                  ( phi[i][j+1][k]^ dWdz[i][j+1][k])*dX.z;
               if ( Grid.BCtypeS[i][k] == BC_REFLECTION) {
                  Wl =  SOLN_pSTATE::Reflect(Wr,  Grid.nfaceS(i, j+1, k));
               }
               if ( Grid.BCtypeS[i][k] == BC_NO_SLIP) {
                  Wl =  SOLN_pSTATE::No_Slip(Wr, WoS[i][k], Grid.nfaceS(i, j+1, k),
                                             IPs.Pressure_Gradient,
                                             FIXED_TEMPERATURE_WALL);
               }
               if ( Grid.BCtypeS[i][k] == BC_MOVING_WALL) {
                  Wl =  SOLN_pSTATE::Moving_Wall(Wr, WoS[i][k], Grid.nfaceS(i, j+1, k),
                                                 IPs.Moving_wall_velocity,
                                                 IPs.Pressure_Gradient,
                                                 FIXED_TEMPERATURE_WALL);
               }
               
            } else if (j ==  JCu &&
                       ( Grid.BCtypeN[i][k] == BC_REFLECTION ||
                         Grid.BCtypeN[i][k] == BC_NO_SLIP||
                         Grid.BCtypeN[i][k] == BC_MOVING_WALL)) {
               dX =  Grid.xfaceN(i, j, k)- Grid.Cell[i][j][k].Xc;
               Wl =  W[i][j][k] +
                  ( phi[i][j][k]^ dWdx[i][j][k])*dX.x+
                  ( phi[i][j][k]^ dWdy[i][j][k])*dX.y+
                  ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
               if ( Grid.BCtypeN[i][k] == BC_REFLECTION) {
                  Wr =  SOLN_pSTATE::Reflect(Wl,  Grid.nfaceN(i, j, k));
               }
               if ( Grid.BCtypeN[i][k] == BC_NO_SLIP) {
                  Wr =  SOLN_pSTATE::No_Slip(Wl, WoN[i][k], Grid.nfaceN(i, j, k),
                                             IPs.Pressure_Gradient,
                                             FIXED_TEMPERATURE_WALL );
               }
               if ( Grid.BCtypeN[i][k] == BC_MOVING_WALL) {
                  Wr =  SOLN_pSTATE::Moving_Wall(Wl, WoN[i][k], Grid.nfaceN(i, j, k),
                                                 IPs.Moving_wall_velocity,
                                                 IPs.Pressure_Gradient,
                                                 FIXED_TEMPERATURE_WALL );
               }
               
            } else {
               dX =  Grid.xfaceN(i, j, k)- Grid.Cell[i][j][k].Xc;
               Wl =  W[i][j][k] +
                  ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                  ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                  ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
               dX =  Grid.xfaceS(i, j+1, k)- Grid.Cell[i][j+1][k].Xc;
               Wr =  W[i][j+1][k] +
                  ( phi[i][j+1][k]^ dWdx[i][j+1][k])*dX.x +
                  ( phi[i][j+1][k]^ dWdy[i][j+1][k])*dX.y +
                  ( phi[i][j+1][k]^ dWdz[i][j+1][k])*dX.z;
            } /* endif */
            
            switch(IPs.i_Flux_Function) {
            case FLUX_FUNCTION_HLLE :
               Flux =  SOLN_pSTATE::FluxHLLE_n(Wl, Wr,  Grid.nfaceN(i, j, k));
               break;
            case FLUX_FUNCTION_ROE :
               Flux =  SOLN_pSTATE::FluxRoe_n(Wl, Wr,  Grid.nfaceN(i, j, k));
               break;
               
            } /* endswitch */
            
            // add viscous flux in j direction
            
            /* Evaluate cell-averaged solution changes. */

            dUdt[i][j][k][k_residual] -=
               (IPs.CFL_Number* dt[i][j][k])*
               Flux* Grid.AfaceN(i, j, k)/
               Grid.volume(i, j, k);
            dUdt[i][j+1][k][k_residual] +=
               (IPs.CFL_Number* dt[i][j+1][k])*
               Flux* Grid.AfaceS(i, j+1, k)/
               Grid.volume(i, j+1, k);

            
            /* Save south and north face boundary flux. */
            
            // if (j ==  JCl-1) {
//               FluxS[i] = -Flux* Grid.lfaceS(i, j+1);
//           } else if (j ==  JCu) {
//               FluxN[i] = Flux* Grid.lfaceN(i, j);
//           } /* endif */
            
         } /* endfor */
         if (( i >  ICl-1 && i <  ICu+1 ) &&
             ( k >  KCl-1 && k <  KCu+1 ) ) {
            dUdt[i][ JCl-1][k][k_residual] = U_VACUUM;
            dUdt[i][ JCu+1][k][k_residual] = U_VACUUM;
         }
 
      } /* endfor */
   }/* endfor */
           
   // Add k-direction (gamma-direction) fluxes.
   for ( i =  ICl; i <=  ICu ; ++i ){
      for ( j  =  JCl ; j <=  JCu ; ++j ){
         for ( k =  KCl-1 ; k <=  KCu ; ++k )  {
            
            /* Evaluate the cell interface j-direction fluxes. */
            if (k ==  KCl-1 &&
                ( Grid.BCtypeB[i][j] == BC_REFLECTION ||
                  Grid.BCtypeB[i][j] == BC_NO_SLIP ||
                  Grid.BCtypeB[i][j] == BC_MOVING_WALL)) {
               
               dX =  Grid.xfaceBot(i, j, k+1)- Grid.Cell[i][j][k+1].Xc;
               Wr =  W[i][j][k+1] +
                  ( phi[i][j][k+1]^ dWdx[i][j][k+1])*dX.x +
                  ( phi[i][j][k+1]^ dWdy[i][j][k+1])*dX.y+
                  ( phi[i][j][k+1]^ dWdz[i][j][k+1])*dX.z;
               if ( Grid.BCtypeB[i][j] == BC_REFLECTION) {
                  Wl =  SOLN_pSTATE::Reflect(Wr,  Grid.nfaceBot(i, j, k+1));
               }
               if ( Grid.BCtypeB[i][j] == BC_NO_SLIP) {
                  Wl =  SOLN_pSTATE::No_Slip(Wr, WoB[i][j], Grid.nfaceBot(i, j, k+1),
                                             IPs.Pressure_Gradient,
                                             FIXED_TEMPERATURE_WALL);
               }
               if ( Grid.BCtypeB[i][j] == BC_MOVING_WALL) {
                  Wl =  SOLN_pSTATE::Moving_Wall(Wr, WoB[i][j], Grid.nfaceBot(i, j, k+1),
                                                 IPs.Moving_wall_velocity,
                                                 IPs.Pressure_Gradient,
                                                 FIXED_TEMPERATURE_WALL);
               }
               
               
            } else if (k ==  KCu &&
                       ( Grid.BCtypeT[i][j] == BC_REFLECTION ||
                         Grid.BCtypeT[i][j] == BC_NO_SLIP||
                         Grid.BCtypeT[i][j] == BC_MOVING_WALL)) {
               
               dX =  Grid.xfaceTop(i, j, k)- Grid.Cell[i][j][k].Xc;
               Wl =  W[i][j][k] +
                  ( phi[i][j][k]^ dWdx[i][j][k])*dX.x+
                  ( phi[i][j][k]^ dWdy[i][j][k])*dX.y+
                  ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
               if ( Grid.BCtypeT[i][j] == BC_REFLECTION) {
                  Wr =  SOLN_pSTATE::Reflect(Wl,  Grid.nfaceTop(i, j, k));
               }
               if ( Grid.BCtypeT[i][j] == BC_NO_SLIP) {
                  Wr =  SOLN_pSTATE::No_Slip(Wl, WoT[i][j], Grid.nfaceTop(i, j, k),
                                             IPs.Pressure_Gradient,
                                             FIXED_TEMPERATURE_WALL );
               }
               if ( Grid.BCtypeT[i][j] == BC_MOVING_WALL) {
                  Wr =  SOLN_pSTATE::Moving_Wall(Wl, WoT[i][j], Grid.nfaceTop(i, j, k),
                                                 IPs.Moving_wall_velocity,
                                                 IPs.Pressure_Gradient,
                                                 FIXED_TEMPERATURE_WALL );
               }
            } else {
               dX =  Grid.xfaceTop(i, j, k)- Grid.Cell[i][j][k].Xc;
               Wl =  W[i][j][k] +
                  ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                  ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                  ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
               dX =  Grid.xfaceBot(i, j, k+1)- Grid.Cell[i][j][k+1].Xc;
               Wr =  W[i][j][k+1] +
                  ( phi[i][j][k+1]^ dWdx[i][j][k+1])*dX.x +
                  ( phi[i][j][k+1]^ dWdy[i][j][k+1])*dX.y +
                  ( phi[i][j][k+1]^ dWdz[i][j][k+1])*dX.z;
               
               
            } /* endif */
             
            switch(IPs.i_Flux_Function) {
            case FLUX_FUNCTION_HLLE :
               Flux =  SOLN_pSTATE::FluxHLLE_n(Wl, Wr,  Grid.nfaceTop(i, j, k));
               break;
            case FLUX_FUNCTION_ROE :
               Flux =  SOLN_pSTATE::FluxRoe_n(Wl, Wr,  Grid.nfaceTop(i, j, k));
               break;
               
            } /* endswitch */

            
            /* Evaluate cell-averaged solution changes. */
            
            dUdt[i][j][k][k_residual] -=
               (IPs.CFL_Number* dt[i][j][k])*
               Flux* Grid.AfaceTop(i, j, k)/
               Grid.volume(i, j, k);
            dUdt[i][j][k+1][k_residual] +=
               (IPs.CFL_Number* dt[i][j][k+1])*
               Flux* Grid.AfaceBot(i, j, k+1)/
               Grid.volume(i, j, k+1);
            
            /* Save top and bottom face boundary flux. */
            
            // if (j ==  JCl-1) {
//               FluxS[i] = -Flux* Grid.lfaceS(i, j+1);
//           } else if (j ==  JCu) {
//               FluxN[i] = Flux* Grid.lfaceN(i, j);
//           } /* endif */
            
         } /* endfor */
         dUdt[i][j][ KCl-1][k_residual] = U_VACUUM;
         dUdt[i][j][ KCu+1][k_residual] = U_VACUUM;
      } /* endfor */
   }/* endfor */
   
   return (0);
    
}

/********************************************************
 * Routine: Update_Solution_Multistage_Explicit         *
 *                                                      *
 * This routine updates solution states of the given    *
 * solution block for a variety of multi-stage explicit *
 * time integration schemes.                            *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Update_Solution_Multistage_Explicit(const int i_stage,
                                                                              Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs) {
   
   int i, j, k,  k_residual;
   double omega;
   
   // Memory for linear system solver.
   
   SOLN_cSTATE dU_precon;
   int NUM_VAR_3D = NumVar();
   /* Allocate memory for linear system solver. */
   /* Perform update of solution variables for stage 
      i_stage of an N stage scheme. */
   /* Evaluate the time step fraction and residual storage 
      location for the stage. */
  
   switch(IPs.i_Time_Integration) {
   case TIME_STEPPING_EXPLICIT_EULER :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      break;
   case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      break;
   case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
        if (IPs.N_Stage == 4) { 
          if (i_stage == 4) { 
             k_residual = 0;
          } else { 
             k_residual = i_stage - 1;
          } // endif  
       } // endif 
      break;
   case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
        
      omega = MultiStage_Optimally_Smoothing(i_stage, 
                                             IPs.N_Stage,
                                             IPs.i_Limiter);
      k_residual = 0;
      break;
   default:
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      break;
   } /* endswitch */
    

   /* Update solution variables for this stage. */
   for ( k  =  KCl ; k <=  KCu ; ++k ) {
      for ( j  =  JCl ; j <=  JCu ; ++j ) {
         for ( i =  ICl ; i <=  ICu ; ++i ) {
            
            if (IPs.Local_Time_Stepping == 
                GLOBAL_TIME_STEPPING || 
                IPs.Local_Time_Stepping == 
                SCALAR_LOCAL_TIME_STEPPING) {
               
               U[i][j][k] = Uo[i][j][k] + omega* dUdt[i][j][k][k_residual];
                             
               //N-1 species
               U[i][j][k][NUM_VAR_3D] = U[i][j][k].rho*(ONE - U[i][j][k].sum_species());
            }
            
            Uo[i][j][k].negative_speccheck( Uo[i][j][k], W[i][j][k].React.reactset_flag );
            if (IPs.Local_Time_Stepping == GLOBAL_TIME_STEPPING && 
                ( U[i][j][k].rho  <= ZERO ||  
                  U[i][j][k].es() <= ZERO)) {
               
               cout << "\n " << CFFC_Name() << 
                  " ERROR: Negative Density, Mass Fractions, and/or Sensible Energy: \n"
                    << " cell = (" << i << ", " << j <<", "<< k << ") " 
                    << " X = " <<  Grid.Cell[i][j][k].Xc << "\n U = " 
                    <<  U[i][j][k] << "\n dUdt = " 
                    <<  dUdt[i][j][k][k_residual] << "\n";
               return (i);
               
            } 
                          
            W[i][j][k] = U[i][j][k].W();
            
         } /* endfor */    	 
      } /* endfor */    
   } /* endfor */
    
    /* Deallocate memory for linear system solver. */

   /* Solution successfully updated. */
    
   return (0);   

}

/**************************************************************************
 * Hexa_Block -- Input-output operators.                                  *
 **************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
ostream &operator << (ostream &out_file,
                      const Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> &SolnBlk){
   
   out_file << SolnBlk.Grid;
   out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " << SolnBlk.ICu << "\n";
   out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << "\n";
   out_file << SolnBlk.NCk << " " << SolnBlk.KCl << " " << SolnBlk.KCu << "\n";
   out_file << SolnBlk.Nghost << "\n";
  
   if (SolnBlk.NCi == 0 || SolnBlk.NCj == 0 || SolnBlk.NCk == 0 || SolnBlk.Nghost == 0) 
      return(out_file);
   
   for (int k=SolnBlk.KCl-SolnBlk.Nghost; k<= SolnBlk.KCu+SolnBlk.Nghost; ++k) {
      for(int j= SolnBlk.JCl-SolnBlk.Nghost; j<= SolnBlk.JCu+SolnBlk.Nghost; ++j) {
         for(int i=SolnBlk.ICl-SolnBlk.Nghost; i<= SolnBlk.ICu+SolnBlk.Nghost; ++i) {
            out_file << SolnBlk.U[i][j][k] << "\n";
         } /* endfor */
      } /* endfor */
   } /* endfor */

   // boundary values
   for (int k = SolnBlk.KCl-SolnBlk.Nghost ; k<= SolnBlk.KCu+SolnBlk.Nghost; ++k ) {
      for (int j = SolnBlk.JCl-SolnBlk.Nghost ; j<= SolnBlk.JCu+SolnBlk.Nghost; ++j ) {
         out_file << SolnBlk.WoW[j][k] << "\n";
         out_file << SolnBlk.WoE[j][k] << "\n";
      } /* endfor */
   } /* endfor */
    
   for (int k = SolnBlk.KCl-SolnBlk.Nghost ; k <= SolnBlk.KCu+SolnBlk.Nghost ; ++k ) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
         out_file << SolnBlk.WoS[i][k] << "\n";
         out_file << SolnBlk.WoN[i][k] << "\n";
      } /* endfor */
   } /* endfor */

   for ( int j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
         out_file << SolnBlk.WoB[i][j] << "\n";
         out_file << SolnBlk.WoT[i][j] << "\n";
      } /* endfor */
   } /* endfor */
    
   return (out_file);

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
istream &operator >> (istream &in_file,
                      Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> &SolnBlk) {
   
   int ni, il, iu, nj, jl, ju, nk, kl, ku, ng;
    
   // Grid3D_Hexa_Block New_Grid; in_file >> New_Grid;
   in_file>> SolnBlk.Grid;
   in_file.setf(ios::skipws);
   in_file >> ni >> il >> iu;
   in_file >> nj >> jl >> ju;
   in_file >> nk >> kl >> ku;
   in_file >> ng;
   in_file.unsetf(ios::skipws);
   
   if (ni == 0 || nj == 0 || nk == 0 || ng == 0) {
      if (SolnBlk.Allocated) SolnBlk.deallocate(); 
      return(in_file);
   } /* endif */
   
   if (!SolnBlk.Allocated || SolnBlk.NCi != ni ||SolnBlk.NCj != nj
       || SolnBlk.NCk != nk || SolnBlk.Nghost != ng) {
      if (SolnBlk.Allocated) SolnBlk.deallocate();
      SolnBlk.allocate(ni, nj, nk, ng);
   } /* endif */

   SOLN_cSTATE U_VACUUM;
   U_VACUUM.Vacuum();
   SOLN_pSTATE W_VACUUM;
   W_VACUUM.Vacuum();

   // Copy_Hexa_Block(SolnBlk.Grid, New_Grid); New_Grid.deallocate();
   for (int k=SolnBlk.KCl-SolnBlk.Nghost; k<= SolnBlk.KCu+SolnBlk.Nghost ; ++k ) {
      for (int j=SolnBlk.JCl-SolnBlk.Nghost; j<= SolnBlk.JCu+SolnBlk.Nghost; ++j) {
         for (int i=SolnBlk.ICl-SolnBlk.Nghost; i<= SolnBlk.ICu+SolnBlk.Nghost; ++i) {
            in_file >> SolnBlk.U[i][j][k];
            SolnBlk.W[i][j][k] = SolnBlk.U[i][j][k].W();
            for (int m = 0 ; m <= NUMBER_OF_RESIDUAL_VECTORS-1 ; ++m ) {
               SolnBlk.dUdt[i][j][k][m] = U_VACUUM;
            } /* endfor */
            SolnBlk.dWdx[i][j][k] = W_VACUUM;
            SolnBlk.dWdy[i][j][k] = W_VACUUM;
            SolnBlk.dWdz[i][j][k] = W_VACUUM;
            
            SolnBlk.phi[i][j][k] = W_VACUUM;
            SolnBlk.Uo[i][j][k] =  U_VACUUM;
            SolnBlk.dt[i][j][k] = ZERO;
         } /* endfor */
      } /* endfor */
   } /* endfor */

   // boundary values
   for (int k = SolnBlk.KCl-SolnBlk.Nghost ; k<= SolnBlk.KCu+SolnBlk.Nghost; ++k ) {
      for (int j = SolnBlk.JCl-SolnBlk.Nghost ; j<= SolnBlk.JCu+SolnBlk.Nghost; ++j ) {
         in_file >> SolnBlk.WoW[j][k];
         in_file >> SolnBlk.WoE[j][k];
      } /* endfor */ 
   } /* endfor */
    
   for (int k = SolnBlk.KCl-SolnBlk.Nghost ; k <= SolnBlk.KCu+SolnBlk.Nghost ; ++k) {
      for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
         in_file >> SolnBlk.WoS[i][k];
         in_file >> SolnBlk.WoN[i][k];
      } /* endfor */
   } /* endfor */

   for ( int j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
         in_file >> SolnBlk.WoB[i][j] ;
         in_file >> SolnBlk.WoT[i][j] ;   
      } /* endfor */
   } /* endfor */
   
   return (in_file);

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::Wall_Shear(void) {

   // for Euler and NavierStokes ... 
   // do nothing
   return (0);
   
}

/*******************************************************************************
 * Hexa_Block::LoadSendBuffer -- Loads send message buffer.                    *
 *******************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::LoadSendBuffer(double *buffer,
                                                         int &buffer_count,
                                                         const int buffer_size,
                                                         const int i_min, 
                                                         const int i_max,
                                                         const int i_inc,
                                                         const int j_min, 
                                                         const int j_max,
                                                         const int j_inc,
						         const int k_min, 
                                                         const int k_max,
                                                         const int k_inc) {
  int i, j, k,t;
  for ( k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc ) 
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for ( t = 1 ; t <= NumVar(); ++ t) {
  	   buffer_count = buffer_count + 1;
	 
	   if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = U[i][j][k][t];
        } /* endfor */
     } /* endfor */
  } /* endfor */

  return(0);

}

/*******************************************************************************
 * Hexa_Block::LoadSendBuffer_F2C -- Loads send message buffer for             *
 *                                   fine to coarse block message passing.     *
 *******************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::LoadSendBuffer_F2C(double *buffer,
                                                             int &buffer_count,
                                                             const int buffer_size,
                                                             const int i_min, 
                                                             const int i_max,
                                                             const int i_inc,
                                                             const int j_min, 
                                                             const int j_max,
                                                             const int j_inc,
					               	     const int k_min, 
                                                             const int k_max,
                                                             const int k_inc) {

   cout << "\nError: LoadSendBuffer_F2C() is not written for Hexa"; cout.flush();

   return(2); 

}

/*******************************************************************************
 * Hexa_Block::LoadSendBuffer_C2F -- Loads send message buffer for             *
 *                                   coarse to fine block message passing.     *
 *******************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::LoadSendBuffer_C2F(double *buffer,
                                                             int &buffer_count,
                                                             const int buffer_size,
                                                             const int i_min, 
                                                             const int i_max,
                                                             const int i_inc,
                                                             const int j_min, 
                                                             const int j_max,
                                                             const int j_inc,
						             const int k_min, 
                                                             const int k_max,
                                                             const int k_inc) {

   cout << "\nError: LoadSendBuffer_C2F() is not written for Hexa"; cout.flush();

   return(2); 

}

/*******************************************************************************
 * Hexa_Block::UnloadReceiveBuffer -- Unloads receive message buffer.          *
 *******************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::UnloadReceiveBuffer(double *buffer,
                                                   int &buffer_count,
                                                   const int buffer_size,
                                                   const int i_min, 
                                                   const int i_max,
                                                   const int i_inc,
                                                   const int j_min, 
                                                   const int j_max,
                                                       const int j_inc,
						       const int k_min, 
                                                       const int k_max,
                                                       const int k_inc) {

  int i, j,k;
  for (k  = k_min ; ((k_inc+1)/2) ? (k <= k_max):(k >= k_max) ; k += k_inc) {
     for (j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc) {
        for (i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc) {
           for (int nV = 1 ; nV <=NumVar() ; ++ nV) {
              buffer_count = buffer_count + 1;
              if (buffer_count >= buffer_size) return(1);    
              U[i][j][k][nV] = buffer[buffer_count];
           } /* endfor */
           W[i][j][k] = U[i][j][k].W();
        } /* endfor */
     } /* endfor */
  } /* endfor */ 

  return(0);

}

/*******************************************************************************
 * Hexa_Block::UnloadReceiveBuffer_F2C -- Unloads receive message buffer for   *
 *                                        fine to coarse block message passing.*
 *******************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::UnloadReceiveBuffer_F2C(double *buffer,
                                                                  int &buffer_count,
                                                                  const int buffer_size,
                                                                  const int i_min, 
                                                                  const int i_max,
                                                                  const int i_inc,
                                                                  const int j_min, 
                                                                  const int j_max,
                                                                  const int j_inc,
						                  const int k_min, 
                                                                  const int k_max,
                                                                  const int k_inc) {

   cout << "\nError: UnloadReceiveBuffer_F2C() is not written for Hexa"; cout.flush();

   return(2); 

}

/*******************************************************************************
 * Hexa_Block::UnloadReceiveBuffer_C2F -- Unloads receive message buffer for   *
 *                                        coarse to fine block message passing.*
 *******************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::UnloadReceiveBuffer_C2F(double *buffer,
                                                                  int &buffer_count,
                                                                  const int buffer_size,
                                                                  const int i_min, 
                                                                  const int i_max,
                                                                  const int i_inc,
                                                                  const int j_min, 
                                                                  const int j_max,
                                                                  const int j_inc,
						                  const int k_min, 
                                                                  const int k_max,
                                                                  const int k_inc) {

   cout << "\nError: UnloadReceiveBuffer_C2F() is not written for Hexa"; cout.flush();

   return(2);

}

/**************************************************************************
 * Hexa_Block::SubcellReconstruction --                                   *
 *             Performs the subcell reconstruction of solution stat  e    *
 *             within a given cell (i,j,k) of the computational mesh for  *
 *             the specified hexahedral solution block.                   *
 **************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::SubcellReconstruction(const int i, 
                                                                 const int j,
                                                                 const int k,
                                                                 const int Limiter) {

  cout << "\nError: SubcellReconstruction() is not written for Hexa"; cout.flush();

  return;
    
}

/*******************************************************************************
 * Hexa_Block::LoadSendBuffer_Flux_F2C -- Loads send message buffer for        *
 *                                        fine to coarse block message passing *
 *                                        of conservative solution fluxes.     *
 *******************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::LoadSendBuffer_Flux_F2C(double *buffer,
                                                                  int &buffer_count,
                                                                  const int buffer_size,
                                                                  const int i_min, 
                                                                  const int i_max,
                                                                  const int i_inc,
                                                                  const int j_min, 
                                                                  const int j_max,
                                                                  const int j_inc,
	       					                  const int k_min, 
                                                                  const int k_max,
                                                                  const int k_inc) {

   cout << "\nError: LoadSendBuffer_Flux_F2C() is not written for Hexa"; cout.flush();

   return(2);


}

/*******************************************************************************
 * Hexa_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message *
 *                                                buffer for fine to coarse    *
 *                                                block message passing of     *
 *                                                conservative solution fluxes.*
 *******************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::UnloadReceiveBuffer_Flux_F2C(double *buffer,
                                                                       int &buffer_count,
                                                                       const int buffer_size,
                                                                       const int i_min, 
                                                                       const int i_max,
                                                                       const int i_inc,
                                                                       const int j_min, 
                                                                       const int j_max,
                                                                       const int j_inc,
						                       const int k_min, 
                                                                       const int k_max,
                                                                       const int k_inc) {

   cout << "\nError: UnloadReceiveBuffer_Flux_F2C() is not written for Hexa"; cout.flush();

   return(2);

}

#endif // _HEXA_BLOCK_INCLUDED
