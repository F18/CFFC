/* AdvectDiffuse2DState.cc:  Subroutines for 
                             2D Advection Diffusion Equation Solution State Class. */

/* Include 2D advection diffusion equation solution state header file. */

#ifndef _ADVECTDIFFUSE2D_STATE_INCLUDED
#include "AdvectDiffuse2DState.h"
#endif // _ADVECTDIFFUSE2D_STATE_INCLUDED


/**********************************************************
 * Routine: Riemann (Riemann solver, x-direction)         *
 *                                                        *
 * This function returns a solution to Riemann problem    *
 * for the 2D advection-diffusion equations in the        *
 * x-direction, returning the intermediate state          *
 * variables along the ray x/t=0.                         *
 *                                                        *
 **********************************************************/
AdvectDiffuse2D_State Riemann(const AdvectDiffuse2D_State &Ul,
	      	              const AdvectDiffuse2D_State &Ur) {
       
    AdvectDiffuse2D_State Um;

    Um.V = HALF*(Ul.V+Ur.V);
    Um.k = HALF*(Ul.k+Ur.k);
    Um.T = HALF*(Ul.T+Ur.T);
    Um.Fd = HALF*(Ul.Fd+Ur.Fd);

    if (Um.V.x >= ZERO) {
      Um.u = Ul.u;
    } else {
      Um.u = Ur.u;
    } /* endif */

    return (Um);

}

/**********************************************************
 * Routine: Riemann_x (Riemann solver, x-direction)       *
 *                                                        *
 * This function returns a solution to Riemann problem    *
 * for the 2D advection-diffusion equations in the        *
 * x-direction, returning the intermediate state          *
 * variables along the ray x/t=0.                         *
 *                                                        *
 **********************************************************/
AdvectDiffuse2D_State Riemann_x(const AdvectDiffuse2D_State &Ul,
	      	                const AdvectDiffuse2D_State &Ur) {
       
    AdvectDiffuse2D_State Um;

    Um.V = HALF*(Ul.V+Ur.V);
    Um.k = HALF*(Ul.k+Ur.k);
    Um.T = HALF*(Ul.T+Ur.T);
    Um.Fd = HALF*(Ul.Fd+Ur.Fd);

    if (Um.V.x >= ZERO) {
      Um.u = Ul.u;
    } else {
      Um.u = Ur.u;
    } /* endif */

    return (Um);

}

/**********************************************************
 * Routine: Riemann_y (Riemann solver, y-direction)       *
 *                                                        *
 * This function returns a solution to Riemann problem    *
 * for the 2D advection-diffusion equations in the        *
 * y-direction, returning the intermediate state          *
 * variables along the ray y/t=0.                         *
 *                                                        *
 **********************************************************/
AdvectDiffuse2D_State Riemann_y(const AdvectDiffuse2D_State &Ul,
	      	                const AdvectDiffuse2D_State &Ur) {

    AdvectDiffuse2D_State Um;

    Um.V = HALF*(Ul.V+Ur.V);
    Um.k = HALF*(Ul.k+Ur.k);
    Um.T = HALF*(Ul.T+Ur.T);
    Um.Fd = HALF*(Ul.Fd+Ur.Fd);

    if (Um.V.y >= ZERO) {
      Um.u = Ul.u;
    } else {
      Um.u = Ur.u;
    } /* endif */

    return (Um);

}

/*********************************************************
 * Routine: Flux (solution flux function, x-direction)   *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states.                                      *
 *                                                       *
 *********************************************************/
double Flux(const AdvectDiffuse2D_State &Ul,
	    const AdvectDiffuse2D_State &Ur) {

    AdvectDiffuse2D_State Um;
    Um = Riemann(Ul, Ur);
    return (Um.F().x+Um.Fd.x);

}

/*********************************************************
 * Routine: Flux_x (solution flux function, x-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states.                                      *
 *                                                       *
 *********************************************************/
double Flux_x(const AdvectDiffuse2D_State &Ul,
	      	const AdvectDiffuse2D_State &Ur) {

    AdvectDiffuse2D_State Um;
    Um = Riemann_x(Ul, Ur);
    return (Um.F().x+Um.Fd.x);

}

/*********************************************************
 * Routine: Flux_y (solution flux function, y-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the y-direction given left and right         *
 * solution states.                                      *
 *                                                       *
 *********************************************************/
double Flux_y(const AdvectDiffuse2D_State &Ul,
	      const AdvectDiffuse2D_State &Ur) {

    AdvectDiffuse2D_State Um;
    Um = Riemann_y(Ul, Ur);
    return (Um.F().y+Um.Fd.y);

}

/*********************************************************
 * Routine: Flux_n (solution flux function, n-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.                       *
 *                                                       *
 *********************************************************/
double Flux_n(const AdvectDiffuse2D_State &Ul,
	      const AdvectDiffuse2D_State &Ur,
              const Vector2D &norm_dir) {

    double cos_angle, sin_angle;
    AdvectDiffuse2D_State Ul_rotated, Ur_rotated, Um_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */

    Ul_rotated.u = Ul.u;
    Ul_rotated.V.x = Ul.V.x*cos_angle +
                     Ul.V.y*sin_angle;
    Ul_rotated.V.y = - Ul.V.x*sin_angle +
                       Ul.V.y*cos_angle;
    Ul_rotated.k = Ul.k;
    Ul_rotated.T = Ul.T;
    Ul_rotated.Fd.x = Ul.Fd.x*cos_angle +
                      Ul.Fd.y*sin_angle;
    Ul_rotated.Fd.y = - Ul.Fd.x*sin_angle +
                        Ul.Fd.y*cos_angle;

    Ur_rotated.u = Ur.u;
    Ur_rotated.V.x = Ur.V.x*cos_angle +
                     Ur.V.y*sin_angle;
    Ur_rotated.V.y = - Ur.V.x*sin_angle +
                       Ur.V.y*cos_angle;
    Ur_rotated.k = Ur.k;
    Ur_rotated.T = Ur.T;
    Ur_rotated.Fd.x = Ur.Fd.x*cos_angle +
                      Ur.Fd.y*sin_angle;
    Ur_rotated.Fd.y = - Ur.Fd.x*sin_angle +
                        Ur.Fd.y*cos_angle;

    /* Solve the Riemann problem in the rotated frame */

    Um_rotated = Riemann(Ul_rotated, Ur_rotated);

    /* Return the solution flux. */

    return (Um_rotated.F().x+Um_rotated.Fd.x);

}

/********************************************************
 * Routine: BC_Dirichlet (Dirichlet Boundary Condition) *
 *                                                      *
 * This function returns the Dirichlet boundary data    *
 * for a given direction given the solution on the      *
 * interior of the boundary, Ui, the desired flow state *
 * to be imposed at the boundary, Uo, and the unit      *
 * normal vector in the direction of interest.          *
 *                                                      *
 ********************************************************/
AdvectDiffuse2D_State BC_Dirichlet(const AdvectDiffuse2D_State &Ui,
                                   const AdvectDiffuse2D_State &Uo,
	      	                   const Vector2D &norm_dir) {

    AdvectDiffuse2D_State Ub;
    return (Ub);
       
}

/********************************************************
 * Routine: BC_Dirichlet (Neumann Boundary Condition)   *
 *                                                      *
 * This function returns the Neumann boundary data      *
 * for a given direction given the solution on the      *
 * interior of the boundary, Ui, the desired flow state *
 * to be imposed at the boundary, Uo, and the unit      *
 * normal vector in the direction of interest.          *
 *                                                      *
 ********************************************************/
AdvectDiffuse2D_State BC_Neumann(const AdvectDiffuse2D_State &Ui,
                                 const AdvectDiffuse2D_State &Uo,
	      	                 const Vector2D &norm_dir) {

    AdvectDiffuse2D_State Ub;
    return (Ub);

}

/********************************************************
 * Routine: BC_Robin (Robin Boundary Condition)         *
 *                                                      *
 * This function returns the Robin boundary data        *
 * for a given direction given the solution on the      *
 * interior of the boundary, Ui, the desired flow state *
 * to be imposed at the boundary, Uo, and the unit      *
 * normal vector in the direction of interest.          *
 *                                                      *
 ********************************************************/
AdvectDiffuse2D_State BC_Robin(const AdvectDiffuse2D_State &Ui,
                               const AdvectDiffuse2D_State &Uo,
	      	               const Vector2D &norm_dir) {

    AdvectDiffuse2D_State Ub;
    return (Ub);
       
}
