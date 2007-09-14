/* Spline2D.cc:  Subroutines for various 2D spline classes. */

/* Include the 2D Spline header file. */

#ifndef _SPLINE2D_INCLUDED
#include "Spline2D.h"
#endif // _SPLINE2D_INCLUDED

/********************************************************
 * Spline2D -- Spline fits in 2 space dimensions.       *
 ********************************************************/

/********************************************************
 * Routine: Spline                                      *
 *                                                      *
 *    This routine performs a piecewise linear,         *
 * quadratic, cubic, or quintic blended spline          *
 * interpolation to return the 2D position vector on    *
 * the splined surface (body or boundary) of interest   *
 * given the path length s along the spline and the     *
 * set of np discrete points (Xp.x,Xp.y) that define    *
 * the spline geometry.                                 *
 *                                                      *
 ********************************************************/
Vector2D Spline(const double &s,
                const Spline2D &S) {

    int i, il, ir, subinterval_found, npts_used;
    double sp_used[4];
    double ds0, ds1, ds2, ds3,
           ds01, ds02, ds12, ds13, ds23,
           a0, a1, a2, a3, a_blend;
    Vector2D Xp_used[4], Xl, Xr; 

    /* Determine the subinterval of the spline
       containing the point (position) of interest. */

    if (s <= S.sp[0]) return(S.Xp[0]);
    if (s >= S.sp[S.np-1]) return(S.Xp[S.np-1]);      

    subinterval_found = 0;
    i = 0;
    while (subinterval_found == 0 && i < S.np-1) {
      if ((s-S.sp[i])*(s-S.sp[i+1]) <= ZERO) {
        subinterval_found = 1;
      } else {
        i += 1;
      } /* endif */
    } /* endwhile */
    il = i;
    ir = i+1;

    /* Determine number and which points to use in the 
       blended spline interpolation. */

    if (S.type == SPLINE2D_CONSTANT) return(S.Xp[il]);
    if (S.type == SPLINE2D_LINEAR) {
       npts_used = 2;
    } else if (S.tp[il] == SPLINE2D_POINT_SHARP_CORNER &&
               S.tp[ir] == SPLINE2D_POINT_SHARP_CORNER) {
       npts_used = 2;
    } else if (S.tp[il] == SPLINE2D_POINT_NORMAL &&
               S.tp[ir] == SPLINE2D_POINT_SHARP_CORNER) {
       npts_used = 3;
       if (il == 0) {
          Xp_used[0]=S.Xp[S.np-2];
          sp_used[0]=S.sp[0]-(S.sp[S.np-1]-S.sp[S.np-2]);
          Xp_used[1]=S.Xp[il];
          sp_used[1]=S.sp[il];
          Xp_used[2]=S.Xp[ir];
          sp_used[2]=S.sp[ir];
       } else {
          Xp_used[0]=S.Xp[il-1];
          sp_used[0]=S.sp[il-1];
          Xp_used[1]=S.Xp[il];
          sp_used[1]=S.sp[il];
          Xp_used[2]=S.Xp[ir];
          sp_used[2]=S.sp[ir];
       } /* endif */
    } else if (S.tp[il] == SPLINE2D_POINT_SHARP_CORNER &&
               S.tp[ir] == SPLINE2D_POINT_NORMAL) {
       npts_used = 3;
       if (il == S.np-2) {
          Xp_used[0]=S.Xp[il];
          sp_used[0]=S.sp[il];
          Xp_used[1]=S.Xp[ir];
          sp_used[1]=S.sp[ir];
          Xp_used[2]=S.Xp[1];
          sp_used[2]=S.sp[S.np-1]+(S.sp[1]-S.sp[0]);
       } else {
          Xp_used[0]=S.Xp[il];
          sp_used[0]=S.sp[il];
          Xp_used[1]=S.Xp[ir];
          sp_used[1]=S.sp[ir];
          Xp_used[2]=S.Xp[ir+1];
          sp_used[2]=S.sp[ir+1];
       } /* endif */
    } else if (S.tp[il] == SPLINE2D_POINT_NORMAL &&
               S.tp[ir] == SPLINE2D_POINT_NORMAL) {
       npts_used = 4;
       if (il == 0) {
          Xp_used[0]=S.Xp[S.np-2];
          sp_used[0]=S.sp[0]-(S.sp[S.np-1]-S.sp[S.np-2]);
          Xp_used[1]=S.Xp[il];
          sp_used[1]=S.sp[il];
          Xp_used[2]=S.Xp[ir];
          sp_used[2]=S.sp[ir];
          Xp_used[3]=S.Xp[ir+1];
          sp_used[3]=S.sp[ir+1];
       } else if (il == S.np-2) {
          Xp_used[0]=S.Xp[il-1];
          sp_used[0]=S.sp[il-1];
          Xp_used[1]=S.Xp[il];
          sp_used[1]=S.sp[il];
          Xp_used[2]=S.Xp[ir];
          sp_used[2]=S.sp[ir];
          Xp_used[3]=S.Xp[1];
          sp_used[3]=S.sp[S.np-1]+(S.sp[1]-S.sp[0]);
       } else {
          Xp_used[0]=S.Xp[il-1];
          sp_used[0]=S.sp[il-1];
          Xp_used[1]=S.Xp[il];
          sp_used[1]=S.sp[il];
          Xp_used[2]=S.Xp[ir];
          sp_used[2]=S.sp[ir];
          Xp_used[3]=S.Xp[ir+1];
          sp_used[3]=S.sp[ir+1];
       } /* endif */
    } /* endif */

    /* Perform the interpolation. */

    if (npts_used == 2) {
       Xl = S.Xp[il]+(s-S.sp[il])*(S.Xp[ir]-S.Xp[il])/
	    (S.sp[ir]-S.sp[il]);
       return(Xl);
    } else if (npts_used == 3) {
       ds0 = s - sp_used[0];
       ds1 = s - sp_used[1];
       ds2 = s - sp_used[2];
       
       ds01 = sp_used[0] - sp_used[1];
       ds02 = sp_used[0] - sp_used[2];
       ds12 = sp_used[1] - sp_used[2];
       
       a0 =  ds1*ds2/(ds01*ds02);
       a1 = -ds0*ds2/(ds01*ds12);
       a2 =  ds0*ds1/(ds02*ds12);
         
       Xl = a0*Xp_used[0] + a1*Xp_used[1] +
	    a2*Xp_used[2];
       return(Xl);
    } else if (npts_used == 4) {
       ds0 = s - sp_used[0];
       ds1 = s - sp_used[1];
       ds2 = s - sp_used[2];
       ds3 = s - sp_used[3];
       
       ds01 = sp_used[0] - sp_used[1];
       ds02 = sp_used[0] - sp_used[2];
       ds12 = sp_used[1] - sp_used[2];
       ds13 = sp_used[1] - sp_used[3];
       ds23 = sp_used[2] - sp_used[3];
       
       a0 =  ds1*ds2/(ds01*ds02);
       a1 = -ds0*ds2/(ds01*ds12);
       a2 =  ds0*ds1/(ds02*ds12);

       Xl = a0*Xp_used[0] + a1*Xp_used[1] +
	    a2*Xp_used[2];

       a1 =  ds2*ds3/(ds12*ds13);
       a2 = -ds1*ds3/(ds12*ds23);
       a3 =  ds1*ds2/(ds13*ds23);

       Xr = a1*Xp_used[1] + a2*Xp_used[2] +
	    a3*Xp_used[3];

       switch(S.type) {
         case SPLINE2D_QUADRATIC :
           Xl = HALF*(Xl+Xr);
           break; 
         case SPLINE2D_CUBIC :
           a_blend = ds1/ds12;
	   Xl = (ONE+a_blend)*Xl - a_blend*Xr;
           break;
         case SPLINE2D_QUINTIC :
           a_blend = ds1/ds12;
           a_blend =-a_blend*a_blend*(THREE+TWO*a_blend);
	   Xl = (ONE+a_blend)*Xl - a_blend*Xr;
           break;
         default:
           Xl = HALF*(Xl+Xr);
           break;
       } /* endswitch */

       return(Xl);

    } /* endif */

}

/********************************************************
 * Routine: BCtype                                      *
 *                                                      *
 *    This routine returns the boundary condition type  *
 * given the path length s along the spline and the     *
 * the set of np boundary types (bc) that have been     *
 * defined for each spline point.                       *
 *                                                      *
 ********************************************************/
int BCtype(const double &s,
           const Spline2D &S) {

    int i, il, ir, subinterval_found;

    /* Determine the subinterval of the spline
       containing the point (position) of interest. */

    if (s <= S.sp[0]) return(S.bc[0]);
    if (s >= S.sp[S.np-1]) return(S.bc[S.np-1]);      

    subinterval_found = 0;
    i = 0;
    while (subinterval_found == 0 && i < S.np-1) {
      if ((s-S.sp[i])*(s-S.sp[i+1]) <= ZERO) {
        subinterval_found = 1;
      } else {
        i += 1;
      } /* endif */
    } /* endwhile */
    il = i;
    ir = i+1;

    /* Determine the boundary condition type. */

    if (S.bc[il] == S.bc[ir]) {
        return(S.bc[il]);
    } else {
        return(S.bc[il]);
    } /* endif */
 
}

/********************************************************
 * Routine: Copy_Spline                                 *
 *                                                      *
 * Copies the spline S2 to spline S1.                   *
 *                                                      *
 ********************************************************/
void Copy_Spline(Spline2D &S1,
	      	 Spline2D &S2) {

    int i;
 
    /* Allocate (re-allocate) memory for the spline S1 
       as necessary. */

    if (S1.np != S2.np) {
       if (S1.np != 0) S1.deallocate();
       if (S2.np >= 2) S1.allocate(S2.np);
    } /* endif */

    /* Set spline type for spline S1. */

    S1.settype(S2.type);

    /* Copy the spline coordinates, pathlength, point type, and boundary 
       condition information from spline S2 to spline S1. */

    if (S2.np >= 2) { 
       for ( i = 0; i <= S1.np-1; ++i ) {
           S1.Xp[i] = S2.Xp[i];
           S1.sp[i] = S2.sp[i];
           S1.tp[i] = S2.tp[i];
           S1.bc[i] = S2.bc[i];
       } /* endfor */
    } /* endif */

}

/********************************************************
 * Routine: Broadcast_Spline                            *
 *                                                      *
 * Broadcasts a spline to all processors involved in    *
 * the calculation from the primary processor using     *
 * the MPI broadcast routine.                           *
 *                                                      *
 ********************************************************/
void Broadcast_Spline(Spline2D &S) {

#ifdef _MPI_VERSION
    int i, npts, buffer_size, i_buffer_size;
    int *i_buffer;
    double *buffer;
 
    /* Broadcast the number of spline points. */

    if (CFDkit_Primary_MPI_Processor()) {
      npts = S.np;
    } /* endif */

    MPI::COMM_WORLD.Bcast(&npts, 1, MPI::INT, 0);

    /* On non-primary MPI processors, allocate (re-allocate) 
       memory for the spline as necessary. */

    if (!CFDkit_Primary_MPI_Processor()) {
       if (S.np != npts) {
          if (S.np != 0) S.deallocate();
          if (npts >= 2) S.allocate(npts);
       } /* endif */
    } /* endif */

    /* Broadcast the spline type. */

    MPI::COMM_WORLD.Bcast(&(S.type), 1, MPI::INT, 0);

    /* Broadcast the the spline coordinates, pathlength, 
       point type, and boundary condition information. */

    if (npts >= 2) {
       buffer = new double[3*npts];
       i_buffer = new int[2*npts];

       if (CFDkit_Primary_MPI_Processor()) {
          buffer_size = 0;
          i_buffer_size = 0;
          for ( i = 0; i <= S.np-1; ++i ) {
              buffer[buffer_size] = S.Xp[i].x;
              buffer[buffer_size+1] = S.Xp[i].y;
              buffer[buffer_size+2] = S.sp[i];
              i_buffer[i_buffer_size] = S.tp[i];
              i_buffer[i_buffer_size+1] = S.bc[i];
              buffer_size = buffer_size + 3;
              i_buffer_size = i_buffer_size + 2;
          } /* endfor */
       } /* endif */

       buffer_size = 3*npts;
       i_buffer_size = 2*npts;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);
       MPI::COMM_WORLD.Bcast(i_buffer, i_buffer_size, MPI::INT, 0);

       if (!CFDkit_Primary_MPI_Processor()) {
          buffer_size = 0;
          i_buffer_size = 0;
          for ( i = 0; i <= S.np-1; ++i ) {
              S.Xp[i].x = buffer[buffer_size];
              S.Xp[i].y = buffer[buffer_size+1];
              S.sp[i] = buffer[buffer_size+2];
              S.tp[i] = i_buffer[i_buffer_size];
              S.bc[i] = i_buffer[i_buffer_size+1];
              buffer_size = buffer_size + 3;
              i_buffer_size = i_buffer_size + 2;
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;
       delete []i_buffer; 
       i_buffer = NULL;
    } /* endif */
#endif

}

#ifdef _MPI_VERSION
/********************************************************
 * Routine: Broadcast_Spline                            *
 *                                                      *
 * Broadcasts a spline to all processors associated     *
 * with the specified communicator from the specified   *
 * processor using the MPI broadcast routine.           *
 *                                                      *
 ********************************************************/
void Broadcast_Spline(Spline2D &S,
                      MPI::Intracomm &Communicator, 
                      const int Source_CPU) {

    int Source_Rank = 0;
    int i, npts, buffer_size, i_buffer_size;
    int *i_buffer;
    double *buffer;
 
    /* Broadcast the number of spline points. */

    if (CFDkit_MPI::This_Processor_Number == Source_CPU) {
      npts = S.np;
    } /* endif */

    Communicator.Bcast(&npts, 1, MPI::INT, Source_Rank);

    /* On non-source MPI processors, allocate (re-allocate) 
       memory for the spline as necessary. */

    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
       if (S.np != npts) {
          if (S.np != 0) S.deallocate();
          if (npts >= 2) S.allocate(npts);
       } /* endif */
    } /* endif */

    /* Broadcast the spline type. */

    Communicator.Bcast(&(S.type), 1, MPI::INT, Source_Rank);

    /* Broadcast the the spline coordinates, pathlength, 
       point type, and boundary condition information. */

    if (npts >= 2) {
       buffer = new double[3*npts];
       i_buffer = new int[2*npts];

       if (CFDkit_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          i_buffer_size = 0;
          for ( i = 0; i <= S.np-1; ++i ) {
              buffer[buffer_size] = S.Xp[i].x;
              buffer[buffer_size+1] = S.Xp[i].y;
              buffer[buffer_size+2] = S.sp[i];
              i_buffer[i_buffer_size] = S.tp[i];
              i_buffer[i_buffer_size+1] = S.bc[i];
              buffer_size = buffer_size + 3;
              i_buffer_size = i_buffer_size + 2;
          } /* endfor */
       } /* endif */

       buffer_size = 3*npts;
       i_buffer_size = 2*npts;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);
       Communicator.Bcast(i_buffer, i_buffer_size, MPI::INT, Source_Rank);

       if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          i_buffer_size = 0;
          for ( i = 0; i <= S.np-1; ++i ) {
              S.Xp[i].x = buffer[buffer_size];
              S.Xp[i].y = buffer[buffer_size+1];
              S.sp[i] = buffer[buffer_size+2];
              S.tp[i] = i_buffer[i_buffer_size];
              S.bc[i] = i_buffer[i_buffer_size+1];
              buffer_size = buffer_size + 3;
              i_buffer_size = i_buffer_size + 2;
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;
       delete []i_buffer; 
       i_buffer = NULL;
    } /* endif */

}
#endif

/********************************************************
 * Routine: Translate_Spline                            *
 *                                                      *
 * Translates or shifts the positions of all of the     *
 * points defining the spline.                          *
 *                                                      *
 ********************************************************/
void Translate_Spline(Spline2D &S,
	              const Vector2D &V) {

    int i;
 
    /* Translate each point defining the spline. */

    for ( i = 0; i <= S.np-1; ++i ) {
        S.Xp[i] += V;
    } /* endfor */

}

/********************************************************
 * Routine: Scale_Spline                                *
 *                                                      *
 * Scale the positions of all of the points defining    *
 * the spline.                                          *
 *                                                      *
 ********************************************************/
void Scale_Spline(Spline2D &S,
	          const double &Scaling_Factor) {

    int i;
 
    /* Scale each point defining the spline. */

    for ( i = 0; i <= S.np-1; ++i ) {
        S.Xp[i] = S.Xp[i]*Scaling_Factor;
        S.sp[i] = S.sp[i]*Scaling_Factor;
    } /* endfor */

}

/********************************************************
 * Routine: Rotate_Spline                               *
 *                                                      *
 * Apply a solid body rotation about the origin and     *
 * recompute the positions of all of the points in      *
 * the spline accordingly.                              *
 *                                                      *
 ********************************************************/
void Rotate_Spline(Spline2D &S,
	           const double &Angle) {

    int i;
    double cos_angle, sin_angle;
    Vector2D X;
 
    /* Apply the rotation to each point defining the spline. */

    cos_angle = cos(-Angle); 
    sin_angle = sin(-Angle);

    for ( i = 0; i <= S.np-1; ++i ) {
        X.x = S.Xp[i].x*cos_angle +
              S.Xp[i].y*sin_angle;
        X.y = - S.Xp[i].x*sin_angle +
                S.Xp[i].y*cos_angle;

        S.Xp[i] = X;
    } /* endfor */

}

/********************************************************
 * Routine: Reflect_Spline                              *
 *                                                      *
 * Apply a mirror reflection about the y=0 axis and     *
 * recompute the positions of all of the points in      *
 * the spline accordingly.                              *
 *                                                      *
 ********************************************************/
void Reflect_Spline(Spline2D &S) {

    int i;
    Vector2D X;
 
    /* Apply a mirror reflection about the y=0 axis. */

    for ( i = 0; i <= S.np-1; ++i ) {
        X.x = S.Xp[i].x;
        X.y = - S.Xp[i].y;
        S.Xp[i] = X;
    } /* endfor */

}

/********************************************************
 * Routine: Reverse_Spline                              *
 *                                                      *
 * Reverses the order of the spline points.             *
 *                                                      *
 ********************************************************/
void Reverse_Spline(Spline2D &S) {

   int i, tp, bc;
   double sp;
   Vector2D Xp;

   /* Reverse the order of the spline points. */

   for ( i = 0; i <= (S.np-1)/2; ++i ) {
      Xp = S.Xp[i];  S.Xp[i] = S.Xp[S.np-1-i];  S.Xp[S.np-1-i] = Xp;
      tp = S.tp[i];  S.tp[i] = S.tp[S.np-1-i];  S.tp[S.np-1-i] = tp;
      bc = S.bc[i];  S.bc[i] = S.bc[S.np-1-i];  S.bc[S.np-1-i] = bc;
      sp = S.sp[i];  S.sp[i] = S.sp[S.np-1-i];  S.sp[S.np-1-i] = sp;
   } /* endfor */

}

/********************************************************
 * Routine: Concatenate_Splines                         *
 *                                                      *
 * Concatenate (combines) spline S1 and S2 and returns  *
 * a new spline which is a combination of both.         *
 *                                                      *
 ********************************************************/
Spline2D Concatenate_Splines(const Spline2D &S1,
	      	             const Spline2D &S2) {

    int i, npts;
    Spline2D Sc;

    /* Set number of points for concatenated spline. */

    npts = S1.np + S2.np - 1;
 
    /* Allocate memory for concatenated spline. */ 

    Sc.allocate(npts);

    /* Set spline type for the concatenated spline. */

    Sc.settype(S1.type);

    /* Copy the spline coordinates, point type, and boundary 
       condition information from spline S1 to the new 
       concatenated spline Sc. */

    for ( i = 0; i <= S1.np-1; ++i ) {
        Sc.Xp[i] = S1.Xp[i];
        Sc.tp[i] = S1.tp[i];
        Sc.bc[i] = S1.bc[i];
    } /* endfor */

    Sc.tp[S1.np-1] = SPLINE2D_POINT_SHARP_CORNER;

    /* Copy the spline coordinates, point type, and boundary 
       condition information from spline S2 to the new 
       concatenated spline Sc. A translation is applied to
       ensure the first point of spline S2 is the same as
       the last point of spline S1. */

    for ( i = 1; i <= S2.np-1; ++i ) {
        Sc.Xp[i+S1.np-1] = S2.Xp[i] + (S1.Xp[S1.np-1]-S2.Xp[0]);
        Sc.tp[i+S1.np-1] = S2.tp[i];
        Sc.bc[i+S1.np-1] = S2.bc[i];
    } /* endfor */

    /* Calculate the concatenated spline pathlengths. */

    Sc.pathlength();

    /* Return the concatenated spline. */

    return(Sc);

}

/********************************************************
 * Routine: Create_Spline_Line                          *
 *                                                      *
 * This routine calculates and returns a 2D spline      *
 * representing a straight line between two points.     *
 *                                                      *
 ********************************************************/
void Create_Spline_Line(Spline2D &Line_Spline,
                        const Vector2D &V1,
			const Vector2D &V2,
  	                const int Number_of_Spline_Points) {

    int i;
    Vector2D dX;

    /* Allocate memory for the straight line spline. */

    if (Line_Spline.np != 0) Line_Spline.deallocate();
    Line_Spline.allocate(Number_of_Spline_Points);

    /* Set the spline type. */

    Line_Spline.settype(SPLINE2D_LINEAR);

    /* Compute the locations of the spline points on the 
       straight line between points V1 and V2. */

    dX = (V2-V1)/double(Number_of_Spline_Points-1);

    for (i = 0; i <= Number_of_Spline_Points-1; i++) {
        Line_Spline.Xp[i] = V1 + double(i)*dX;

        if (i == 0 || i == Number_of_Spline_Points - 1) {
           Line_Spline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
        } else {
           Line_Spline.tp[i] = SPLINE2D_POINT_NORMAL;
        } /* endif */

        Line_Spline.bc[i] = BC_NONE;
    } /* endfor */

    /* Calculate the spline pathlengths. */

    Line_Spline.pathlength();

}

/********************************************************
 * Routine: Create_Spline_Circular_Arc                  *
 *                                                      *
 * This routine calculates and returns a 2D spline      *
 * representing a segment of a circular arc.            *
 *                                                      *
 ********************************************************/
void Create_Spline_Circular_Arc(Spline2D &Circle_Spline,
			        const Vector2D &Origin,
				const double &Radius,
                                const double &Angle1,
			        const double &Angle2,
  	                        const int Number_of_Spline_Points) {

    int i;
    double theta;

    /* Allocate memory for the circular arc spline. */

    if (Circle_Spline.np != 0) Circle_Spline.deallocate();
    Circle_Spline.allocate(Number_of_Spline_Points);

    /* Set the spline type. */

    Circle_Spline.settype(SPLINE2D_QUINTIC);

    /* Compute the locations of the spline points on the 
       circular arc. */

    for (i = 0; i <= Number_of_Spline_Points-1; i++) {
        theta = Angle1+
                (Angle2-Angle1)*double(i)/double(Number_of_Spline_Points-1);
        theta = TWO*PI*theta/360.0;

        Circle_Spline.Xp[i].x = Radius*cos(theta);
        Circle_Spline.Xp[i].y = Radius*sin(theta);

        Circle_Spline.Xp[i] += Origin;        

        if (i == 0 || i == Number_of_Spline_Points - 1) {
           Circle_Spline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
        } else {
           Circle_Spline.tp[i] = SPLINE2D_POINT_NORMAL;
        } /* endif */

        Circle_Spline.bc[i] = BC_NONE;
    } /* endfor */

    /* Ensure that arc closes on itself if it is a
       complete circle. */

    if ( fabs(fabs(Angle2-Angle1)-360.00) < TOLER ) {
       Circle_Spline.Xp[Circle_Spline.np-1] = Circle_Spline.Xp[0];
       Circle_Spline.tp[0] = SPLINE2D_POINT_NORMAL;
       Circle_Spline.tp[Circle_Spline.np-1] = SPLINE2D_POINT_NORMAL;
    } /* endif */

    /* Calculate the spline pathlengths. */

    Circle_Spline.pathlength();

}

/********************************************************
 * Routine: Create_Spline_Ellipsoidal_Arc               *
 *                                                      *
 * This routine calculates and returns a 2D spline      *
 * representing a segment of an ellipsoidal arc.        *
 *                                                      *
 ********************************************************/
void Create_Spline_Ellipsoidal_Arc(Spline2D &Ellipse_Spline,
			           const Vector2D &Origin,
				   const double &A,
				   const double &B,
                                   const double &Angle1,
			           const double &Angle2,
  	                           const int Number_of_Spline_Points) {

    int i;
    double theta;

    /* Allocate memory for the ellipsoidal arc spline. */

    if (Ellipse_Spline.np != 0) Ellipse_Spline.deallocate();
    Ellipse_Spline.allocate(Number_of_Spline_Points);

    /* Set the spline type. */

    Ellipse_Spline.settype(SPLINE2D_QUINTIC);

    /* Compute the locations of the spline points on the 
       ellipsoidal arc. */

    for (i = 0; i <= Number_of_Spline_Points-1; i++) {
        theta = Angle1+
                (Angle2-Angle1)*double(i)/double(Number_of_Spline_Points-1);
        theta = TWO*PI*theta/360.0;

        Ellipse_Spline.Xp[i].x = A*cos(theta);
        Ellipse_Spline.Xp[i].y = B*sin(theta);

        Ellipse_Spline.Xp[i] += Origin;        

        if (i == 0 || i == Number_of_Spline_Points - 1) {
           Ellipse_Spline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
        } else {
           Ellipse_Spline.tp[i] = SPLINE2D_POINT_NORMAL;
        } /* endif */

        Ellipse_Spline.bc[i] = BC_NONE;
    } /* endfor */

    /* Ensure that arc closes on itself if it is a
       complete ellipse. */

    if ( fabs(fabs(Angle2-Angle1)-360.00) < TOLER ) {
       Ellipse_Spline.Xp[Ellipse_Spline.np-1] = Ellipse_Spline.Xp[0];
       Ellipse_Spline.tp[0] = SPLINE2D_POINT_NORMAL;
       Ellipse_Spline.tp[Ellipse_Spline.np-1] = SPLINE2D_POINT_NORMAL;
    } /* endif */

    /* Calculate the spline pathlengths. */

    Ellipse_Spline.pathlength();

}

/********************************************************
 * Routine: Create_Spline_NACA_Aerofoil                 *
 *                                                      *
 * This routine calculates and returns a 2D spline for  *
 * both NACA 4- and 5-digit aerofoil sections.  The     *
 * algorithm for calculating the splines is straight    *
 * out of Abbott and von Doenhoff, with a little        *
 * Reigels mixed in.                                    *
 *                                                      *
 * Basically, a NACA aerofoil is composed of a camber   *
 * line and a thickness distribution.  The thickness    *
 * distribution is a single equation, while the camber  *
 * is usually two joined quadratics.                    *
 *                                                      *
 * The equations for the upper and lower coordinates    *
 * are:                                                 *
 *                                                      *
 *  x(upper) =                                          *
 *    x - yt*sin(theta)  y(upper) = yc + yt*cos(theta)  *
 *  x(lower) =                                          *
 *    x + yt*sin(theta)  y(lower) = yc - yt*cos(theta)  *
 *                                                      *
 * where tan(theta) = d(yc)/dx.  In these equations, yc *
 * is the camber line, yt is the thickness distribution.*
 * A common approximation (small-angle) is to assume    *
 * theta is small, so that sin(theta) is approx. 0 and  *
 * cos(theta) is approx. 1.  The equations become:      *
 *                                                      *
 * 	  x(upper) = x     y(upper) = yc + yt           *
 * 	  x(lower) = x     y(lower) = yc - yt           *
 *                                                      *
 * For 4-digit aerofoils, the camber lines and          *
 * thickness:                                           *
 *                                                      *
 * (yc/c) = (f/c)*(1/(x1^2))*(2*x1*(x/c) - (x/c)^2)     *
 *              for 0<=(x/c)<=x1                        *
 *                                                      *
 * and                                                  *
 *                                                      *
 * (yc/c) = (f/c)*(1/(1-x1)^2)*((1-2x1)+2x1*(x/c)-      *
 *          (x/c)^2)                                    *
 *              for x1<=(x/c)<=1 with x1=(xf/c)         *
 *                                                      *
 * (yt/c) = 5t*(0.29690*x^0.5 - 0.12600X - 0.35160*x^2  *
 *         + 0.28430*x^3 - 0.10150*x^4)                 *
 *                                                      *
 * where t = thickness/chord,                           *
 * x = position along x-axis,                           *
 * xf = position of maximum camber,                     *
 * f = maximum camber.                                  *
 *                                                      *
 * For 5-digit aerofoils the thickness distribution is  *
 * the same, only the camber line is changed.  There are*
 * two types, based on the third digit.  The majority   *
 * of 5-digit aerofoils are 'type 0' (i.e.,             *
 * NACA 23015) hand have a camber line given by:        *
 *                                                      *
 * (yc/c) = (k1/6)*(x^3 - 3*x1*x^2 + x1^2*(3-x1)*x)     *
 *               for 0<=x<=x1                           *
 *                                                      *
 * and                                                  *
 *                                                      *
 * (yc/c) = (k1*x1^3 / 6) * (1 - x)                     *
 *        for x1<=x<=1 with x1 = (x1/c) and x = (x/c),  *
 * and x1 is related to xf, which is the position of    *
 * maximum camber.
 *                                                      * 
 * The constants x1 and k1 are determined from the      *
 * following table:                                     *
 *                                                      *
 * xf           0.05    0.10    0.15     0.20    0.25   *
 * x1           0.0580  0.1260  0.2025   0.2900  0.3910 *
 * (Cl* /(f/c)) 26.9    19.6    16.4     14.5    11.3   *
 * (k1/Cl*)     1205    172.1   53.2     22.13   10.77  *
 *                                                      *
 * The 'type 1' aerofoil camber line is not provided    *
 * here.                                                *
 *                                                      *
 * Breakdown of the NACA designations:                  *
 *                                                      *
 * In a 4-digit aerofoil, the first digit is the value  *
 * of the maximum camber (in percent of the chord), the *
 * second digit is the position of the maximum camber   *
 * from the leading edge in tenths of the chord, and    *
 * the last two digits denote the maximum thickness of  *
 * the aerofoil in percent.  For the NACA 2415 aerofoil,*
 * the maximum camber is 2%, the position of the        *
 * maximum camber is 0.4c, and the thickness is 15%.    *
 * This example is displayed in below.                  *
 *                                                      *
 * NACA    2415                                         *
 *                                                      *
 * 24      mean line (24)                               *
 * 2       maximum camber, (2%)                         *
 * 4       10 * position of maximum camber, (0.4c)      *
 * 15      thickness, (15%)                             *
 *                                                      *
 * The NACA 5-digit aerofoils are set up in a similar   *
 * manner to the 4-digit aerofoils.  The primary        *
 * difference is the use of a different camber line.    *
 * In a 5-digit aerofoil, 1.5 times the first digit is  *
 * the design lift coefficient in tenths, the second    *
 * and third digits are one-half the distance from the  *
 * leading edge to the location of maximum camber in    *
 * percent of the chord, and the fourth and fifth       *
 * digits are the thickness in percent of the chord.    *
 * For example, a NACA 23015 aerofoil has a design lift *
 * coefficient of 0.3, has the maximum camber at 0.15c, *
 * and is 15% thick.  Additionally, the first three     *
 * digits indicate the mean line used.  In this case,   *
 * the mean line designation is 230.  The 5-digit       *
 * aerofoils use the same thickness distribution as the *
 * 4-digit aerofoils.  This example is displayed below. *
 *                                                      *
 * NACA    23015                                        *
 *                                                      *
 * 230     mean line (230)                              *
 * 2       (design lift coefficient * 10) / 1.5, (0.3)  *
 * 30      2 * (position of maximum camber),            *
 *         (0.30 / 2 = 0.15c)                           *
 * 0       type of camber line used                     *
 * 15      thickness, (15%)                             *
 *                                                      *
 * References:                                          *
 *                                                      *
 * I. H. Abbot and von A. E. Doenhoff, "Theory of       *
 * Wing Sections", Dover Publications, 1959.            *
 *                                                      *
 * F. W. Reigels, "Aerofoil Sections", Butterworth &    *
 * Co., 1961.                                           *
 *                                                      *
 ********************************************************/
void Create_Spline_NACA_Aerofoil(Spline2D &NACA_Spline,
                                 char *NACA_Aerofoil_Type_ptr,
                                 const double &Chord_Length,
                                 const int i_Up_All_Low,
  	                         const int Number_of_Spline_Points) {

    int i, number_of_digits;
    char camber_max_indicator[2], x_camber_max_indicator[3], 
         thickness_max_indicator[3], design_lift_indicator[2];
    double camber_max, x_camber_max, thickness_max, design_lift,
           x1, k1, x_thickness_scaling_factor=1.008930411356;
    double theta, x, y, xt, yt, yc;

    /* Allocate memory for the NACA aerofoil spline. */

    if (NACA_Spline.np != 0) NACA_Spline.deallocate();
    NACA_Spline.allocate(Number_of_Spline_Points);

    /* Set the spline type. */

    NACA_Spline.settype(SPLINE2D_QUINTIC);

    /* Determine the number of digits in the NACA Aerofoil
       type designator. */

    number_of_digits = strlen(NACA_Aerofoil_Type_ptr);

    /* Determine the aerofoil chamber line and thickness 
       distribution parameters. */

    if (number_of_digits == 4) {
      camber_max_indicator[0] = NACA_Aerofoil_Type_ptr[0];
      camber_max_indicator[1] = '\n';
      camber_max = double(atof(camber_max_indicator)/HUNDRED);
      x_camber_max_indicator[0] = NACA_Aerofoil_Type_ptr[1];
      x_camber_max_indicator[1] = '\n';
      x_camber_max = double(atof(x_camber_max_indicator)/TEN);
      thickness_max_indicator[0] = NACA_Aerofoil_Type_ptr[2];
      thickness_max_indicator[1] = NACA_Aerofoil_Type_ptr[3];
      thickness_max_indicator[2] = '\n';
      thickness_max = double(atof(thickness_max_indicator)/HUNDRED);
    } else if (number_of_digits == 5) {
      design_lift_indicator[0] = NACA_Aerofoil_Type_ptr[0];
      design_lift_indicator[1] = '\n';
      design_lift = double(1.50*atof(design_lift_indicator)/TEN);
      x_camber_max_indicator[0] = NACA_Aerofoil_Type_ptr[1];
      x_camber_max_indicator[1] = NACA_Aerofoil_Type_ptr[2];
      x_camber_max_indicator[2] = '\n';
      x_camber_max = double(HALF*atof(x_camber_max_indicator)/HUNDRED);
      thickness_max_indicator[0] = NACA_Aerofoil_Type_ptr[3];
      thickness_max_indicator[1] = NACA_Aerofoil_Type_ptr[4];
      thickness_max_indicator[2] = '\n';
      thickness_max = double(atof(thickness_max_indicator)/HUNDRED);
      if (fabs(x_camber_max-0.05) < TOLER) {
        x1 = 0.0580;
	k1 = 1205.0*design_lift;
      } else if (fabs(x_camber_max-0.10) < TOLER) {
        x1 = 0.1260;
	k1 = 172.1*design_lift;
      } else if (fabs(x_camber_max-0.15) < TOLER) {
        x1 = 0.2025;
	k1 = 53.2*design_lift;
      } else if (fabs(x_camber_max-0.20) < TOLER) {
        x1 = 0.2900;
	k1 = 22.13*design_lift;
      } else if (fabs(x_camber_max-0.25) < TOLER) {
        x1 = 0.3910;
	k1 = 10.77*design_lift;
      } else {
        x1 = 0.3910;
	k1 = 10.77*design_lift;
      } /* endif */
    } else {
    } /* end if */

    /* Compute the locations of the spline points on the surface
       of the aerofoil. */

    for (i = 0; i <= Number_of_Spline_Points-1; i++) {
        // Calculate theta.
        if ( i_Up_All_Low == 0 ) {
           theta = TWO*PI*double(i)/double(Number_of_Spline_Points-1);
	} else if ( i_Up_All_Low > 0 ) {
           theta = PI + PI*double(i)/double(Number_of_Spline_Points-1);
        } else {
           theta = PI*double(i)/double(Number_of_Spline_Points-1);
        } /* endif */

        // Calculate x location.
        x = HALF*(ONE+cos(theta));

        // Determine the y position of the cord line, yc.
        if (number_of_digits == 4) {
           if (x < x_camber_max && x_camber_max != ZERO) {
	      yc = camber_max*x*(TWO*x_camber_max-x)/sqr(x_camber_max);
           } else {
              yc = camber_max*(ONE-TWO*x_camber_max+x*(TWO*x_camber_max-x))/
	           sqr(ONE-x_camber_max);
           } /* endif */
        } else if (number_of_digits == 5) {
           if (x < x1) {
	      yc = (k1/SIX)*(x*x*x - THREE*x1*x*x + 
		   x1*x1*(THREE-x1)*x); 
           } else {
	      yc = (k1/SIX)*x1*x1*x1*(ONE-x); 
           } /* endif */
        } else {
        } /* end if */

        // Determine the local aerofoil thickness, yt.
        xt = x_thickness_scaling_factor*x;
        yt = FIVE*thickness_max*(0.29690*sqrt(xt) - 
                                 0.12600*xt -
                                 0.35160*xt*xt + 
	                         0.28430*xt*xt*xt - 
                                 0.10150*xt*xt*xt*xt);

        // Evaluate the coordinates of the points on the aerofoil surface.
        if ( theta >= PI ) {
	   y = yc + yt;
	} else {
           y = yc - yt;
        } /* endif */

        // Ensure closure on trailing edge.
        if ( ( i_Up_All_Low == 0 && 
             (i == 0 || i == Number_of_Spline_Points - 1) ) ||
             ( i_Up_All_Low > 0 && i == Number_of_Spline_Points - 1 ) ||
             ( i_Up_All_Low < 0 && i == 0 ) ) {
           x = ONE;
           y = ZERO; 
        } /* endif */

        // Scale the resulting aerofoil coordinates.
        x = Chord_Length*x;
        y = Chord_Length*y;

        // Assign values to the spline variables.
        NACA_Spline.Xp[i].x = x;
        NACA_Spline.Xp[i].y = y;

        if (i == 0 || i == Number_of_Spline_Points - 1) {
           NACA_Spline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
        } else {
           NACA_Spline.tp[i] = SPLINE2D_POINT_NORMAL;
        } /* endif */

        NACA_Spline.bc[i] = BC_NONE;
    } /* endfor */

    /* Calculate the spline pathlengths. */

    NACA_Spline.pathlength();

}

/********************************************************
 * Routine: Create_Spline_Bow_Shock                     *
 *                                                      *
 * This routine calculates and returns a 2D spline      *
 * representing the approximate position of the bow     *
 * shock for supersonic flow over a circular cylinder.  *
 *                                                      *
 ********************************************************/
void Create_Spline_Bow_Shock(Spline2D &Shock_Spline,
                             const double &Radius,
			     const double &Mach_Number,
			     const int i_Up_All_Low,
  	                     const int Number_of_Spline_Points) {

    int i;
    double theta, x, y, d, rs, b, y_shock_max;

    /* Allocate memory for the bow shock spline. */

    if (Shock_Spline.np != 0) Shock_Spline.deallocate();
    Shock_Spline.allocate(Number_of_Spline_Points);

    /* Set the spline type. */

    Shock_Spline.settype(SPLINE2D_QUINTIC);

    /* Compute the locations of the spline points which
       define the position of the bow shock. */

    d  = 0.386*Radius*exp(4.67/(Mach_Number*Mach_Number));
    rs = 1.386*Radius*exp(1.89/pow(Mach_Number-ONE, 0.75));
    b  = asin(ONE/Mach_Number);
    y_shock_max = (rs/tan(b))*sqrt(sqr(ONE+(Radius+d)*tan(b)*tan(b)/rs)-ONE);

    for (i = 0; i <= Number_of_Spline_Points-1; i++) {
        // Calculate theta.
        if ( i_Up_All_Low == 0 ) {
           theta = -HALF*PI + PI*double(i)/double(Number_of_Spline_Points-1);
	} else if ( i_Up_All_Low > 0 ) {
           theta = HALF*PI*double(i)/double(Number_of_Spline_Points-1);
        } else {
           theta = -HALF*PI + HALF*PI*double(i)/double(Number_of_Spline_Points-1);
        } /* endif */

        y = y_shock_max*sin(theta);
        x = -(Radius+d-(rs/(tan(b)*tan(b)))*(sqrt(ONE+sqr(y*tan(b)/rs))-ONE));

        Shock_Spline.Xp[i].x = x;
        Shock_Spline.Xp[i].y = y;

        if (i == 0 || i == Number_of_Spline_Points - 1) {
           Shock_Spline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
        } else {
           Shock_Spline.tp[i] = SPLINE2D_POINT_NORMAL;
        } /* endif */

        Shock_Spline.bc[i] = BC_NONE;
    } /* endfor */

    /* Calculate the spline pathlengths. */

    Shock_Spline.pathlength();

}

/********************************************************
 * Routine: Create_Spline_Area_Variation                *
 *                                                      *
 * This routine calculates and returns a 2D spline      *
 * representing the radius as a function of axial       *
 * position for a duct with a smoothly varying change   *
 * in cross-sectional area.                             *
 *                                                      *
 ********************************************************/
void Create_Spline_Area_Variation(Spline2D &Radius_Spline,
                                  const double &Xup,
				  const double &Xthroat,
                                  const double &Xdown,
				  const double &Rup,
				  const double &Rthroat,
				  const double &Rdown,
  	                          const int Number_of_Spline_Points) {

    int i;
    double x, radius, Aup, Athroat, Adown;

    /* Allocate memory for the area variation spline. */

    if (Radius_Spline.np != 0) Radius_Spline.deallocate();
    Radius_Spline.allocate(Number_of_Spline_Points);

    /* Set the spline type. */

    Radius_Spline.settype(SPLINE2D_QUINTIC);

    /* Compute the locations of the spline points which
       define the radius of the duct as a function of the
       axial postion. */

    Aup = PI*Rup*Rup;
    Athroat = PI*Rthroat*Rthroat;
    Adown = PI*Rdown*Rdown;

    for (i = 0; i <= Number_of_Spline_Points-1; i++) {
        x = (Xdown-Xup)*double(i)/double(Number_of_Spline_Points-1);

        if (x <= Xthroat) {
           radius = Aup*exp(HALF*log(Athroat/Aup)*
                    (ONE-cos(PI*(x-Xup)/(Xthroat-Xup))));
           radius = sqrt(radius/PI);
        } else {
           radius = Athroat*exp(0.50*log(Adown/Athroat)*
                    (ONE-cos(PI*(x-Xthroat)/(Xdown-Xthroat))));
	   radius = sqrt(radius/PI);
        } /* endif */

        Radius_Spline.Xp[i].x = x;
        Radius_Spline.Xp[i].y = radius;

        if (i == 0 || i == Number_of_Spline_Points - 1) {
           Radius_Spline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
        } else {
           Radius_Spline.tp[i] = SPLINE2D_POINT_NORMAL;
        } /* endif */

        Radius_Spline.bc[i] = BC_NONE;
    } /* endfor */

    /* Calculate the spline pathlengths. */

    Radius_Spline.pathlength();

}

/********************************************************
 * Routine: Create_Spline_Rectangle                     *
 *                                                      *
 * This routine calculates and returns a 2D spline      *
 * representing a rectangle.                            *
 *                                                      *
 ********************************************************/
void Create_Spline_Rectangle(Spline2D &Rectangle_Spline,
			     const Vector2D &Origin,
			     const double &Length,
			     const double &Width) {

    int i;

    /* Allocate memory for the rectangular spline. */

    if (Rectangle_Spline.np != 0) Rectangle_Spline.deallocate();
    Rectangle_Spline.allocate(5);

    /* Set the spline type. */

    Rectangle_Spline.settype(SPLINE2D_LINEAR);

    /* Compute the locations of the spline points on the 
       rectangle. */

    Rectangle_Spline.Xp[0].x = Origin.x - HALF*Width;
    Rectangle_Spline.Xp[0].y = Origin.y - HALF*Length;
    Rectangle_Spline.Xp[1].x = Origin.x + HALF*Width;
    Rectangle_Spline.Xp[1].y = Origin.y - HALF*Length;
    Rectangle_Spline.Xp[2].x = Origin.x + HALF*Width;
    Rectangle_Spline.Xp[2].y = Origin.y + HALF*Length;
    Rectangle_Spline.Xp[3].x = Origin.x - HALF*Width;
    Rectangle_Spline.Xp[3].y = Origin.y + HALF*Length;
    Rectangle_Spline.Xp[4]   = Rectangle_Spline.Xp[0];

    /* Set the point and boundary condition type. */

    for (i = 0; i < 5; i++) {
      Rectangle_Spline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
      Rectangle_Spline.bc[i] = BC_NONE;
    }

    /* Calculate the spline pathlengths. */

    Rectangle_Spline.pathlength();

}

/********************************************************
 * Routine: Create_Spline_Zalesaks_Disk                 *
 *                                                      *
 * This routine calculates and returns a 2D spline      *
 * representing Zalesak's disk.                         *
 *                                                      *
 ********************************************************/
void Create_Spline_Zalesaks_Disk(Spline2D &Zalesak_Spline,
				 const Vector2D &Origin,
				 const double &Radius) {

    int npts;
    int i;
    double theta, theta1, theta2;

    /* Set the number of spline points. */

    npts = 341 + 3;

    /* Set the angles. */

    theta1 = 10.0;
    theta2 = 350.0;

    /* Allocate memory for the circular arc spline. */

    if (Zalesak_Spline.np != 0) Zalesak_Spline.deallocate();
    Zalesak_Spline.allocate(npts);

    /* Set the spline type. */

    Zalesak_Spline.settype(SPLINE2D_LINEAR);

    /* Compute the locations of the spline points on the 
       spline and the point and boundary condition types. */

    for (i = 0; i < npts-3; i++) {
      theta = theta1 + (theta2 - theta1)*double(i)/double(npts-4);
      theta = TWO*PI*theta/360.0;
      Zalesak_Spline.Xp[i].x = Origin.x + Radius*cos(theta);
      Zalesak_Spline.Xp[i].y = Origin.y + Radius*sin(theta);
      Zalesak_Spline.bc[i]   = BC_NONE;
      Zalesak_Spline.tp[i]   = SPLINE2D_POINT_NORMAL;
      if (i == 0 || i == npts - 4)
	Zalesak_Spline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    }
    Zalesak_Spline.Xp[npts-3] = Zalesak_Spline.Xp[npts-4] - Vector2D(1.5*Radius,ZERO);
    Zalesak_Spline.bc[npts-3] = BC_NONE;
    Zalesak_Spline.tp[npts-3] = SPLINE2D_POINT_SHARP_CORNER;
    Zalesak_Spline.Xp[npts-2] = Zalesak_Spline.Xp[0] - Vector2D(1.5*Radius,ZERO);
    Zalesak_Spline.bc[npts-2] = BC_NONE;
    Zalesak_Spline.tp[npts-2] = SPLINE2D_POINT_SHARP_CORNER;

    /* Last point must be the same as the first point. */

    Zalesak_Spline.Xp[npts-1] = Zalesak_Spline.Xp[0];
    Zalesak_Spline.bc[npts-1] = Zalesak_Spline.bc[0];
    Zalesak_Spline.tp[npts-1] = Zalesak_Spline.tp[0];

    /* Calculate the spline pathlengths. */

    Zalesak_Spline.pathlength();

    /* Rotate the spline. */
    Rotate_Spline(Zalesak_Spline,-HALF*PI);

}

/**********************************************************************
 * Routine: Create_Spline_Ringleb_Flow                                *
 *                                                                    *
 * This routine calculates and returns a 2D spline representing       *
 * Ringleb's flow domain.                                             *
 *                                                                    *
 **********************************************************************/
void Create_Spline_Ringleb_Flow(Spline2D &Ringleb_Spline) {

  int nk, nq;
  double **rho;
  double delta_k, delta_q;
  double k_init, q_init, q_final;
  double **q, **k, qo, ko, c, J;
  double g = 1.40;
  Spline2D Spline_North, Spline_East, Spline_West, Spline_Reflection;

  // Allocate memory.
  nk = 32;
  nq = 50;
  k = new double*[nk];
  for (int i = 0; i < nk; i++) k[i] = new double[nq];
  q = new double*[nk];
  for (int i = 0; i < nk; i++) q[i] = new double[nq];
  rho = new double*[nk];
  for (int i = 0; i < nk; i++) rho[i] = new double[nq];

  // Determine q, k.
  delta_k = (1.50-0.75)/double(nk-1);
  k_init = 0.75;
  q_init = 0.50;
  for (int i = 0; i < nk; i++) {
    ko = k_init + double(i)*delta_k;
    q_final = ko; // condition y = 0
    delta_q = (q_final - q_init)/double(nq-1);
    for (int j = 0; j < nq; j++) {
      qo = q_init + double(j)*delta_q;
      k[i][j] = ko;
      q[i][j] = qo;
    }
  }

  // Allocate memory for the interface splines.
  Spline_North.allocate(nk);  Spline_North.settype(SPLINE2D_QUINTIC);
  Spline_East.allocate(nq);   Spline_East.settype(SPLINE2D_QUINTIC);
  Spline_West.allocate(nq);   Spline_West.settype(SPLINE2D_QUINTIC);

  // Determine the spline points.
  for (int i = 0; i < nk; i++) {
    for (int j = 0; j < nq; j++) {

      c = sqrt(ONE - ((g-ONE)/TWO)*q[i][j]*q[i][j]);
      rho[i][j] = pow(c,TWO/(g-ONE));
      J = ONE/c + ONE/(THREE*pow(c,THREE)) + ONE/(FIVE*pow(c,FIVE)) - HALF*log((ONE+c)/(ONE-c));
      // NORTH spline.
      if (j == 0) {
	Spline_North.Xp[nk-1-i].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Spline_North.Xp[nk-1-i].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Spline_North.bc[nk-1-i] = BC_RINGLEB_FLOW;
	if (i == 0 || i == nk-1) Spline_North.tp[nk-1-i] = SPLINE2D_POINT_SHARP_CORNER;
	else Spline_North.tp[nk-1-i] = SPLINE2D_POINT_NORMAL;
      }
      // EAST spline.
      if (i == 0) {
	Spline_East.Xp[j].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Spline_East.Xp[j].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Spline_East.bc[j] = BC_FIXED;//REFLECTION;
	if (j == 0 || j == nq-1) Spline_East.tp[j] = SPLINE2D_POINT_SHARP_CORNER;
	else Spline_East.tp[j] = SPLINE2D_POINT_NORMAL;
      }
      // WEST spline.
      if (i == nk-1) {
	Spline_West.Xp[nq-1-j].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Spline_West.Xp[nq-1-j].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Spline_West.bc[nq-1-j] = BC_REFLECTION;
	if (j == 0 || j == nq-1) Spline_West.tp[nq-1-j] = SPLINE2D_POINT_SHARP_CORNER;
	else Spline_West.tp[nq-1-j] = SPLINE2D_POINT_NORMAL;
      }

    }
  }

  // Calculate spline path lengths.
  Spline_North.pathlength();
  Spline_East.pathlength();
  Spline_West.pathlength();

  // Set the boundary condition types for each of the interface splines.
  Spline_North.setBCtype(BC_RINGLEB_FLOW);
  Spline_East.setBCtype(BC_FIXED);//REFLECTION);
  Spline_West.setBCtype(BC_REFLECTION);
  //Spline_East.setBCtype(BC_RINGLEB_FLOW);
  //Spline_West.setBCtype(BC_RINGLEB_FLOW);

  // Reflect Spline_West and concatenate to Spline_West.
  Copy_Spline(Spline_Reflection,Spline_West);
  Reflect_Spline(Spline_Reflection);
  Reverse_Spline(Spline_Reflection);
  Spline_Reflection = Concatenate_Splines(Spline_Reflection,Spline_West);
  Copy_Spline(Spline_West,Spline_Reflection);

  // Reflect Spline_East and concatenate to Spline_East.
  Copy_Spline(Spline_Reflection,Spline_East);
  Reflect_Spline(Spline_Reflection);
  Reverse_Spline(Spline_Reflection);
  Spline_East = Concatenate_Splines(Spline_East,Spline_Reflection);

  // Concatenate sub-splines to the Ringleb spline.
  Ringleb_Spline = Concatenate_Splines(Spline_West,Spline_North);
  Ringleb_Spline = Concatenate_Splines(Ringleb_Spline,Spline_East);

  // Deallocate the memory for the interface splines.
  Spline_North.deallocate();
  Spline_East.deallocate();
  Spline_West.deallocate();

  // Deallocate memory for k, q, and rho.
  for (int i = 0; i < nk; i++) {
    delete []k[i];   k[i]   = NULL;
    delete []q[i];   q[i]   = NULL;
    delete []rho[i]; rho[i] = NULL;
  }
  delete []k;   k   = NULL;
  delete []q;   q   = NULL;
  delete []rho; rho = NULL;

}

/****************************************************************
 * Routine: getY                                                *
 *                                                              *
 *   This routine returns the splined value of Y given X.       *
 *                                                              *
 *                                                              *
 *   Returns a VECTOR2D LINKED LIST containing all the vectors  *
 *   on the spline with the desired X value.                    *
 *                                                              *
 ****************************************************************/
Vector2D_LL getY(const double &x, const Spline2D &S) {
    int i, iNext, found, icount;
    double s1, s2, s3, s3_old, x1, x2, x3;
    Vector2D_LL LL;

    iNext=0;

    while (0<1) {
		  
      /* Do incremental search to find interval. */

      found=0;
      for (i=iNext; i<=S.np-2; ++i) {
        x1=S.Xp[i].x;
        x2=S.Xp[i+1].x;
        s1=S.sp[i];
        s2=S.sp[i+1];
        if ((x>=x1 && x<=x2) || (x<=x1 && x>=x2)) {
	  iNext=i+1;
	  found=1;
  	  break;
        } /* endif */
      } /* endfor */
    
      if (found==0) break;

      /* Use method of false position to find 
         value of 's' that will give 'x'. */

      s3_old = s1;
      s3 = s2;
      icount = 0;
      while ( fabs(s3-s3_old) > TOLER*TOLER*s3) {
        s3_old = s3;
        s3 = (s1*(x2-x) - s2*(x1-x))/(x2-x1);
        x3=Spline(s3,S).x;
     
        if ((x1-x)*(x3-x) < ZERO) {
	  x2=x3;
  	  s2=s3;
        } else {
	  x1=x3;
	  s1=s3;
        } /* endif */
        icount += 1;
        if (icount > 1000) break;
      } /* endwhile */
    
      /* Check for duplicate. */

      if (LL.np==0) {
        LL.addNode( Vector2D(x,Spline(s3,S).y) );
      } else if (LL[LL.np-1].y!=Spline(s3,S).y) {
        LL.addNode( Vector2D(x,Spline(s3,S).y) );
      } /* endif */ 

    } /* endwhile */
  
    return(LL);
  
}

/****************************************************************
 * Routine: getX                                                *
 *                                                              *
 *   This routine returns the splined value of X given Y.       *
 *                                                              *
 *                                                              *
 *   Returns a VECTOR2D LINKED LIST containing all the vectors  *
 *   on the spline with the desired Y value.                    *
 *                                                              *
 ****************************************************************/
Vector2D_LL getX(const double &y, const Spline2D &S) {
    int i, iNext, found, icount;
    double s1, s2, s3, s3_old, y1, y2, y3;
    Vector2D_LL LL;

    iNext=0;

    while (0<1) {
		  
      /* Do incremental search to find interval. */

      found=0;
      for (i=iNext; i<=S.np-2; ++i) {
        y1=S.Xp[i].y;
        y2=S.Xp[i+1].y;
        s1=S.sp[i];
        s2=S.sp[i+1];
        if ((y>=y1 && y<=y2) || (y<=y1 && y>=y2)) {
	  iNext=i+1;
	  found=1;
	  break;
        } /* endif */
      } /* endfor */
    
      if (found==0) break;

      /* Use method of false position to find 
         value of 's' that will give 'y'. */

      s3_old = s1;
      s3 = s2;
      icount = 0;
      while ( fabs(s3-s3_old) > TOLER*TOLER*s3 ) {
        s3_old = s3;
        s3 = (s1*(y2-y) - s2*(y1-y))/(y2-y1);
        y3=Spline(s3,S).y;
     
        if ((y1-y)*(y3-y) < ZERO) {
	  y2=y3;
	  s2=s3;
        } else {
	  y1=y3;
	  s1=s3;
        } /* endif */
        icount += 1;
        if (icount > 1000) break;
      } /* endwhile */

      /* Check for duplicates. */

      if (LL.np==0) {
        LL.addNode( Vector2D(Spline(s3,S).x, y) ); 
      } else if(LL[LL.np-1].x!=Spline(s3,S).x) {
        LL.addNode( Vector2D(Spline(s3,S).x, y) );
      } /* endif */

    } /* endwhile */

    return(LL);
  
}

/****************************************************************
 * Routine: getS                                                *
 *                                                              *
 *   This routine returns the pathlength of vector location X   *
 *   on spline S.                                               *
 *                                                              *
 ****************************************************************/
double getS(const Vector2D &X, const Spline2D &S) {
    int i, subinterval_found, icount;
    double s1, s2, s3, s3_old, x1, x2, x3, y1, y2, y3;

    /* Determine the subinterval of the spline
       containing the point, X, of interest. */

    subinterval_found = 0;
    i = 0;
    while (subinterval_found == 0 && i < S.np-1) {
      x1=S.Xp[i].x;
      x2=S.Xp[i+1].x;
      y1=S.Xp[i].y;
      y2=S.Xp[i+1].y;
      s1=S.sp[i];
      s2=S.sp[i+1];
      if ( ((X.x-x1)*(X.x-x2) <= ZERO) &&
           ((X.y-y1)*(X.y-y2) <= ZERO) ) {
        subinterval_found = 1;
      } else {
        i += 1;
      } /* endif */
    } /* endwhile */

    if (!subinterval_found) {
       subinterval_found = 0;
       i = 0;
       while (subinterval_found == 0 && i < S.np-1) {
          x1=S.Xp[i].x;
          x2=S.Xp[i+1].x;
          y1=S.Xp[i].y;
          y2=S.Xp[i+1].y;
          s1=S.sp[i];
          s2=S.sp[i+1];
          if ( ((X.x-x1)*(X.x-x2) <= sqr(TOLER*max(fabs(X.x),TOLER))) &&
               ((X.y-y1)*(X.y-y2) <= sqr(TOLER*max(fabs(X.y),TOLER))) ) {
            subinterval_found = 1;
          } else {
            i += 1;
          } /* endif */
       } /* endwhile */
       if (!subinterval_found) {
          subinterval_found = 0;
          i = 0;
          while (subinterval_found == 0 && i < S.np-1) {
             x1=S.Xp[i].x;
             x2=S.Xp[i+1].x;
             y1=S.Xp[i].y;
             y2=S.Xp[i+1].y;
             s1=S.sp[i];
             s2=S.sp[i+1];
             if ( ((X.x-x1)*(X.x-x2) <= sqr(TOLER*max(fabs(X.x),TOLER))) ||
                  ((X.y-y1)*(X.y-y2) <= sqr(TOLER*max(fabs(X.y),TOLER))) )  {
               subinterval_found = 1;
             } else {
               i += 1;
             } /* endif */
          } /* endwhile */
          if (!subinterval_found) {
	     //cout << X << "\n" << S << "\n"; cout.flush();
             return(S.sp[0]);
          } /* endif */
       } /* endif */
    } /* endif */


    /* Use method of false position to find value of S
       that will yield vector position X. */

    if ( (fabs(y2-y1) > fabs(x2-x1)) &&
         ((X.y-y1)*(X.y-y2) <= sqr(TOLER*X.y)) ) { 
       s3_old = s1;
       s3 = s2;
       icount = 0;
       while ( fabs(s3-s3_old) > TOLER*TOLER*s3 ) {
         s3_old = s3;
         s3 = (s1*(y2-X.y) - s2*(y1-X.y))/(y2-y1);
         y3=Spline(s3,S).y;
     
         if ((y1-X.y)*(y3-X.y) < ZERO) {
	   y2=y3;
   	   s2=s3;
         } else {
	   y1=y3;
	   s1=s3;
         } /* endif */

         icount += 1;
         if (icount > 1000) break;
       } /* endwhile */
    } else {
       s3_old = s1;
       s3 = s2;
       icount = 0;
       while ( fabs(s3-s3_old) > TOLER*TOLER*s3 ) {
         s3_old = s3;
         s3 = (s1*(x2-X.x) - s2*(x1-X.x))/(x2-x1);
         x3=Spline(s3,S).x;
     
         if ((x1-X.x)*(x3-X.x) < ZERO) {
	   x2=x3;
   	   s2=s3;
         } else {
	   x1=x3;
	   s1=s3;
         } /* endif */

         icount += 1;
         if (icount > 1000) break;
       } /* endwhile */
    } /* endif */

    /* Return pathlength. */

    return(s3);
  
}

/****************************************************************
 * Routine: getBCtype                                           *
 *                                                              *
 *   This routine returns the boundary condition of vector      *
 *   location X on spline S.                                    *
 *                                                              *
 ****************************************************************/
int getBCtype(const Vector2D &X, const Spline2D &S) {

  int i, subinterval_found, icount, bctype, bc1, bc2;
  Vector2D x1, x2, Xs, Xx, Xy;

  // Determine the subinterval of the spline containing the point X.

  Xy = getminY(X,S);
  Xx = getminX(X,S);
  if (abs(X-Xy) <= abs(X-Xx)) Xs = Xy;
  else Xs = Xx;

  subinterval_found = 0;
  i = 0;
  while (subinterval_found == 0 && i < S.np-1) {
    x1 = S.Xp[i];
    x2 = S.Xp[i+1];
    bc1 = S.bc[i];
    bc2 = S.bc[i+1];
    if (((Xs.x-x1.x)*(Xs.x-x2.x) <= ZERO) && ((Xs.y-x1.y)*(Xs.y-x2.y) <= ZERO)) {
      subinterval_found = 1;
    } else {
      i += 1;
    }
  }

  if (!subinterval_found) {
    subinterval_found = 0;
    i = 0;
    while (subinterval_found == 0 && i < S.np-1) {
      x1 = S.Xp[i];
      x2 = S.Xp[i+1];
      bc1 = S.bc[i];
      bc2 = S.bc[i+1];
      if (((Xs.x-x1.x)*(Xs.x-x2.x) <= sqr(TOLER*max(fabs(Xs.x),TOLER))) &&
	  ((Xs.y-x1.y)*(Xs.y-x2.y) <= sqr(TOLER*max(fabs(Xs.y),TOLER)))) {
	subinterval_found = 1;
      } else {
	i += 1;
      }
    }
    if (!subinterval_found) {
      subinterval_found = 0;
      i = 0;
      while (subinterval_found == 0 && i < S.np-1) {
	x1 = S.Xp[i];
	x2 = S.Xp[i+1];
	bc1 = S.bc[i];
        bc2 = S.bc[i+1];
	if (((X.x-x1.x)*(X.x-x2.x) <= sqr(TOLER*max(fabs(X.x),TOLER))) ||
	    ((X.y-x1.y)*(X.y-x2.y) <= sqr(TOLER*max(fabs(X.y),TOLER))))  {
	  subinterval_found = 1;
	} else {
	  i += 1;
	}
      }
      if (!subinterval_found) {
	cout << "ERROR: Can't find spline segment for X =" << Xs << endl; cout.flush();
	return BC_REFLECTION;
      }
    }
  }

  if (!subinterval_found) {
    cout << "ERROR: Can't find spline segment for X =" << Xs << endl; cout.flush();
    return BC_REFLECTION;
  }

  // Determine the boundary condition type on that spline segment.
  if (bc1 == bc2) bctype = bc1;
  else bctype = max(bc1,bc2);

  // Return the boundary condition type.
  return bctype;
  
}

/**********************************************************************
 * Routine: getminY                                                   *
 *                                                                    *
 * This routine returns the nearest splined value of Y given X.       *
 *                                                                    *
 **********************************************************************/
Vector2D getminY(const Vector2D &X, const Spline2D &S) {

  int iNext, found, icount;
  double s1, s2, s3, s3_old, x1, x2, x3;
  Vector2D Xp;

  iNext = 0;
  Xp = Vector2D(MILLION,MILLION);

  while (0 < 1) {
		  
    // Do incremental search to find interval.
    found = 0;
    for (int i = iNext; i <= S.np-2; i++) {
      x1 = S.Xp[i].x;
      x2 = S.Xp[i+1].x;
      s1 = S.sp[i];
      s2 = S.sp[i+1];
      if ((X.x >= x1 && X.x <= x2) || (X.x <= x1 && X.x >= x2)) {
	iNext=i+1;
	found=1;
	break;
      }
    }

    if (found==0) break;

    // Use method of false position to find value of 's' that will give 'x'.
    s3_old = s1;
    s3 = s2;
    icount = 0;
    while (fabs(s3-s3_old) > TOLER*TOLER*s3) {
      s3_old = s3;
      s3 = (s1*(x2-X.x) - s2*(x1-X.x))/(x2-x1);
      x3=Spline(s3,S).x;
      if ((x1-X.x)*(x3-X.x) < ZERO) {
	x2 = x3;
	s2 = s3;
      } else {
	x1 = x3;
	s1 = s3;
      }
      icount += 1;
      if (icount > 1000) break;
    }

    // Choose nearest point..
    if (abs(X-Vector2D(X.x,Spline(s3,S).y)) < abs(X-Xp)) 
      Xp = Vector2D(X.x,Spline(s3,S).x);

  }

  // Return the nearest spline point.
  return Xp;
  
}

/**********************************************************************
 * Routine: getminX                                                   *
 *                                                                    *
 * This routine returns the nearest splined value of X given Y.       *
 *                                                                    *
 **********************************************************************/
Vector2D getminX(const Vector2D &X, const Spline2D &S) {

  int iNext, found, icount;
  double s1, s2, s3, s3_old, y1, y2, y3;
  Vector2D Xp;

  iNext = 0;
  Xp = Vector2D(MILLION,MILLION);

  while (0 < 1) {
		  
    // Do incremental search to find interval.
    found = 0;
    for (int i = iNext; i < S.np-1; i++) {
      y1 = S.Xp[i].y;
      y2 = S.Xp[i+1].y;
      s1 = S.sp[i];
      s2 = S.sp[i+1];
      if ((X.y >= y1 && X.y <= y2) || (X.y <= y1 && X.y >= y2)) {
	iNext=i+1;
	found=1;
	break;
      }
    }
    
    if (found == 0) break;

    // Use method of false position to find value of 's' that will give 'y'.
    s3_old = s1;
    s3 = s2;
    icount = 0;
    while (fabs(s3-s3_old) > TOLER*TOLER*s3) {
      s3_old = s3;
      s3 = (s1*(y2-X.y) - s2*(y1-X.y))/(y2-y1);
      y3 = Spline(s3,S).y;
      if ((y1-X.y)*(y3-X.y) < ZERO) {
	y2 = y3;
	s2 = s3;
      } else {
	y1 = y3;
	s1 = s3;
      }
      icount += 1;
      if (icount > 1000) break;
    }

    // Choose nearest point..
    if (abs(X-Vector2D(Spline(s3,S).x,X.y)) < abs(X-Xp)) 
      Xp = Vector2D(Spline(s3,S).x,X.y);

  }

  // Return the nearest spline point.
  return Xp;
  
}
