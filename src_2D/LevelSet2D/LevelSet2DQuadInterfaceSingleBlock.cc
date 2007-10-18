/******************************************************************//**
 * \file LevelSet2DQuadInterfaceSingleBlock.cc:                       
 *              Single-block versions of subroutines for 2D Level Set 
 *              multi-block quadrilateral mesh solution classes.      
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.
#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Single Block External Subroutines.        *
 **********************************************************************/

/******************************************************************//**
 * Routine: Set_Interface_List                                        
 *                                                                    
 * Sets the interface list to the given list.                         
 *                                                                    
 **********************************************************************/
int Set_Interface_List(LevelSet2D_Quad_Block &SolnBlk,
		       const Interface2D_List &Interface_List) {

  // Copy the given interface list.
  SolnBlk.Interface_List.Copy(Interface_List);

  // Interface list set successfully.
  return 0;

}

/******************************************************************//**
 * Routine: Initialize_Interfaces                                     
 *                                                                    
 * Constructs the initial interfaces: assigns appropriate flags and   
 * reference values and constructs the actual interface spline.       
 *                                                                    
 **********************************************************************/
int Initialize_Interfaces(LevelSet2D_Quad_Block &SolnBlk,
			  LevelSet2D_Input_Parameters &IP) {

  int error_flag, NP;
  double theta, Angle1, Angle2, Radius;

  // Allocate memory for the interface data.
  SolnBlk.Interface_List.allocate(IP.Interface_IP.Component_List.Ni);

  // Construct each interface.
  for (int ni = 1; ni <= SolnBlk.Interface_List.Ni; ni++) {

    // Copy interface type.
    SolnBlk.Interface_List[ni].Type = IP.Interface_IP.Component_List[ni].Type;
    // Copy interface lengths.
    SolnBlk.Interface_List[ni].Length1 = IP.Interface_IP.Component_List[ni].Length1;
    SolnBlk.Interface_List[ni].Length2 = IP.Interface_IP.Component_List[ni].Length2;
    // Copy interface motion type.
    SolnBlk.Interface_List[ni].Motion = IP.Interface_IP.Component_List[ni].Motion;
    // Copy interface motion speed.
    SolnBlk.Interface_List[ni].Speed = IP.Interface_IP.Component_List[ni].Speed;
    // Copy interface reference point.
    SolnBlk.Interface_List[ni].Xref = IP.Interface_IP.Component_List[ni].Xref;

    // Deallocate memory for interface spline if already allocated.
    if (SolnBlk.Interface_List[ni].Spline.np != 0) 
      SolnBlk.Interface_List[ni].Spline.deallocate();

    // Construct the interface.
    switch(SolnBlk.Interface_List[ni].Type) {
    case INTERFACE_CIRCLE :
      // Create circular spline.
      Create_Spline_Circular_Arc(SolnBlk.Interface_List[ni].Spline,
				 IP.Interface_IP.Component_List[ni].Spline.Xp[0],
				 SolnBlk.Interface_List[ni].Length1,
				 ZERO,
				 360.0,
				 361);
      // Determine the reference point (centroid) of the interface.
      SolnBlk.Interface_List[ni].Centroid();
      break;
    case INTERFACE_ELLIPSE :
      // Create circular spline.
      Create_Spline_Ellipsoidal_Arc(SolnBlk.Interface_List[ni].Spline,
				    IP.Interface_IP.Component_List[ni].Spline.Xp[0],
				    SolnBlk.Interface_List[ni].Length1,
				    SolnBlk.Interface_List[ni].Length2,
				    ZERO,
				    360.0,
				    361);
      // Determine the reference point (centroid) of the interface.
      SolnBlk.Interface_List[ni].Centroid();
      break;
    case INTERFACE_SQUARE :
      // Create square spline.
      Create_Spline_Rectangle(SolnBlk.Interface_List[ni].Spline,
			      IP.Interface_IP.Component_List[ni].Spline.Xp[0],
			      SolnBlk.Interface_List[ni].Length1,
			      SolnBlk.Interface_List[ni].Length1);
      // Determine the reference point (centroid) of the interface.
      SolnBlk.Interface_List[ni].Centroid();
      break;
    case INTERFACE_RECTANGLE :
      // Create square spline.
      Create_Spline_Rectangle(SolnBlk.Interface_List[ni].Spline,
			      IP.Interface_IP.Component_List[ni].Spline.Xp[0],
			      SolnBlk.Interface_List[ni].Length1,
			      SolnBlk.Interface_List[ni].Length2);
      // Determine the reference point (centroid) of the interface.
      SolnBlk.Interface_List[ni].Centroid();
      break;
    case INTERFACE_ROCKET_PROPELLANT_GRAIN :
      // Create solid propellant spline.
      // Allocate memory for the circular arc spline. 
      SolnBlk.Interface_List[ni].Spline.allocate(3);
      // Set the spline type.
      SolnBlk.Interface_List[ni].Spline.settype(SPLINE2D_LINEAR);
      // Set the spline points.
      SolnBlk.Interface_List[ni].Spline.Xp[0] = IP.Interface_IP.Component_List[ni].Spline.Xp[0];
      SolnBlk.Interface_List[ni].Spline.Xp[1] = IP.Interface_IP.Component_List[ni].Spline.Xp[1];
      SolnBlk.Interface_List[ni].Spline.Xp[2] = IP.Interface_IP.Component_List[ni].Spline.Xp[2];
      // Set spline tp and bc.
      SolnBlk.Interface_List[ni].Spline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
      SolnBlk.Interface_List[ni].Spline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
      SolnBlk.Interface_List[ni].Spline.tp[2] = SPLINE2D_POINT_SHARP_CORNER;
      for (int j = 0; j < 3; j++) SolnBlk.Interface_List[ni].Spline.bc[j] = BC_NONE;
      // Calculate the spline pathlengths.
      SolnBlk.Interface_List[ni].Spline.pathlength();
      // Determine the reference point (centroid) of the interface.
      SolnBlk.Interface_List[ni].Xref = IP.Interface_IP.Component_List[ni].Xref;
      break;
    case INTERFACE_NACA0012_AEROFOIL :
      // Create NACA aerofoil spline.
      Create_Spline_NACA_Aerofoil(SolnBlk.Interface_List[ni].Spline,
				  "0012",//IP.NACA_Aerofoil_Type,
				  SolnBlk.Interface_List[ni].Length1,
				  0,
				  41);
      SolnBlk.Interface_List[ni].Translate(IP.Interface_IP.Component_List[ni].Spline.Xp[0]);
      // Determine the reference point (centroid) of the interface.
      SolnBlk.Interface_List[ni].Centroid();
      break;
    case INTERFACE_NACA0015_AEROFOIL :
      // Create NACA aerofoil spline.
      Create_Spline_NACA_Aerofoil(SolnBlk.Interface_List[ni].Spline,
				  "0015",//IP.NACA_Aerofoil_Type,
				  SolnBlk.Interface_List[ni].Length1,
				  0,
				  41);
      SolnBlk.Interface_List[ni].Translate(IP.Interface_IP.Component_List[ni].Spline.Xp[0]);
      // Determine the reference point (centroid) of the interface.
      SolnBlk.Interface_List[ni].Centroid();
      break;
    case INTERFACE_ZALESAK :
      // Create Zalesak's Disk.
      SolnBlk.Interface_List[ni].Zalesak(IP.Interface_IP.Component_List[ni].Spline.Xp[0],
					 SolnBlk.Interface_List[ni].Length1);
      // Determine the reference point (centroid) of the interface.
      SolnBlk.Interface_List[ni].Xref = IP.Interface_IP.Component_List[ni].Xref;
      break;
    case INTERFACE_STAR:
      // Create multi-point star.
      SolnBlk.Interface_List[ni].Star(IP.Interface_IP.Component_List[ni].Spline.Xp[0],
				      SolnBlk.Interface_List[ni].Length1,
				      7);
      // Determine the reference point (centroid) of the interface.
      SolnBlk.Interface_List[ni].Xref = IP.Interface_IP.Component_List[ni].Xref;
      break;
    case INTERFACE_USER_SPECIFIED :
      // Create user specified spline.
      // Allocate memory for the circular arc spline. 
      SolnBlk.Interface_List[ni].Spline.allocate(IP.Interface_IP.Component_List[ni].Spline.np);
      // Set the spline type.
      SolnBlk.Interface_List[ni].Spline.settype(SPLINE2D_LINEAR);
      // Set the spline points, tp, and bc.
      for (int np = 0; np < SolnBlk.Interface_List[ni].Spline.np; np++) {
	SolnBlk.Interface_List[ni].Spline.Xp[np] = IP.Interface_IP.Component_List[ni].Spline.Xp[np];
	SolnBlk.Interface_List[ni].Spline.tp[np] = IP.Interface_IP.Component_List[ni].Spline.tp[np];
	SolnBlk.Interface_List[ni].Spline.bc[np] = IP.Interface_IP.Component_List[ni].Spline.bc[np];
      }
      // Calculate the spline pathlengths.
      SolnBlk.Interface_List[ni].Spline.pathlength();
      // Determine the reference point (centroid) of the interface.
      SolnBlk.Interface_List[ni].Xref = IP.Interface_IP.Component_List[ni].Xref;
      break;
    };
    // Initialize the interface velocity function.
    SolnBlk.Interface_List[ni].Initialize_Velocity_Function(SolnBlk.Interface_List[ni].Motion);
    // Set the interface velocity function.
    SolnBlk.Interface_List[ni].Set_Velocity_Function(ZERO);
    // Sort the interface into counter-clockwise order.
    SolnBlk.Interface_List[ni].Sort();
    // Determine the bounding box for each interface.
    SolnBlk.Interface_List[ni].BoundingBox(IP.Extension_Distance);
  }

  // Interfaces successfully created.
  return 0;

}

int Initialize_Interfaces(LevelSet2D_Quad_Block &SolnBlk,
			  LevelSet2D_Input_Parameters &IP,
			  const Interface2D_List &Interface_List) {

  int error_flag;

  // Copy the input interface list to the local interface list.
  SolnBlk.Interface_List.Copy(Interface_List);

  // Recompute the bounding box for each interface.
  for (int ni = 1; ni <= SolnBlk.Interface_List.Ni; ni++) {
    SolnBlk.Interface_List[ni].BoundingBox(IP.Extension_Distance);
  }

  // Interfaces successfully created.
  return 0;

}

/******************************************************************//**
 * Routine: Exact_Initial_Extension                                   
 *                                                                    
 * Exact initialization of the level set function and the extended    
 * front speeds.  Only available for line, circle, ellipse, and       
 * square interface.                                                  
 *                                                                    
 **********************************************************************/
int Exact_Initial_Extension(LevelSet2D_Quad_Block &SolnBlk,
			    LevelSet2D_Input_Parameters &IP) {

  int error_flag;

  // Compute the signed distances and front speed exactly.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {

      switch(SolnBlk.Interface_List[1].Type) {
      case INTERFACE_LINE :
	break;
      case INTERFACE_CIRCLE :
	SolnBlk.U[i][j].psi = abs(SolnBlk.Grid.Cell[i][j].Xc-SolnBlk.Interface_List[1].Xref) - SolnBlk.Interface_List[1].Length1;
	break;
      case INTERFACE_ELLIPSE :
	break;
      case INTERFACE_SQUARE :
	break;
      case INTERFACE_RECTANGLE :
	break;
      case INTERFACE_USER_SPECIFIED :
	// NOTE THIS IS HARD-CODED FOR A SPECIFIC CASE... constant volume bomb.
	SolnBlk.U[i][j].psi = SolnBlk.Interface_List[1].Spline.Xp[0].y-SolnBlk.Grid.Cell[i][j].Xc.y;
	break;
      };

      // Intentionally destroy the signed distance function to test the Eikonal equation.
      if (IP.Perturb_Distance_Function) SolnBlk.U[i][j].psi *= 0.1 + sqr(SolnBlk.Grid.Cell[i][j].Xc.x+3.5) + sqr(SolnBlk.Grid.Cell[i][j].Xc.y+2.0);

    }
  }

  // Exact extension problem computed successfully.
  return 0;

}

/******************************************************************//**
 * Routine: Geometric_Extension_Problem                               
 *                                                                    
 * Geometric solution of the extension problems involving the level   
 * set function and the normal front speed.  This algorithm based on  
 * triangles constructed out of the point of contention and the two   
 * nearest interface nodes.  Note, that this algorithm does not, in   
 * general, initialize the level set function as a signed distance    
 * funtion.  It does, however, produce a decent initial estimate.     
 * Iterative solution of the Eikonal equations (with unit from speed) 
 * is required to produce ensure a signed distance function.          
 *                                                                    
 **********************************************************************/
int Geometric_Extension_Problem(LevelSet2D_Quad_Block &SolnBlk,
				LevelSet2D_Input_Parameters &IP) {

  int error_flag;
  double d = MILLION, h = MILLION;
  double d12, d23, d31, costheta, dt;
//   double f, f2, f3, ft, g, A, A2, A3;
  Vector2D p1, p2, p3;
  Vector2D Xs1, Xs2, Xm1, Xm2, Xp;
  double epsilon = TOLER*min(SolnBlk.Grid.lfaceN(SolnBlk.Grid.ICl,SolnBlk.Grid.JCl),
			     SolnBlk.Grid.lfaceE(SolnBlk.Grid.ICl,SolnBlk.Grid.JCl));

  // Initialize the signed distance function to a million and the speed
  // function to zero.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      SolnBlk.U[i][j].psi = MILLION;
      SolnBlk.U[i][j].F   = ZERO;
    }
  }

  // Ensure that the bounding box has been computed for all of the
  // interfaces.
  for (int ni = 1; ni <= SolnBlk.Interface_List.Ni; ni++) {
    SolnBlk.Interface_List[ni].BoundingBox(IP.Extension_Distance);
  }

  // Determine the signed distance at each cell center.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {

      d = MILLION; // acute/right calculation.
      h = MILLION; // obtuse calcuation.

      // Determine the mesh cell-center.
      p1 = SolnBlk.Grid.Cell[i][j].Xc;

      // The geometric extension problem is solved for each interface
      // and the solution for the minimum signed distance is used.
      for (int ni = 1; ni <= SolnBlk.Interface_List.Ni; ni++) {

	// Only solve the geometric extension if the current point, p, is
	// contained within the bounding box of the interface.
	if (p1.x >= SolnBlk.Interface_List[ni].Xmin.x &&
	    p1.y >= SolnBlk.Interface_List[ni].Xmin.y &&
	    p1.x <= SolnBlk.Interface_List[ni].Xmax.x &&
	    p1.y <= SolnBlk.Interface_List[ni].Xmax.y) {

	  // For each interface segment.
	  for (int np = 0; np < SolnBlk.Interface_List[ni].Spline.np-1; np++) {
	  
	    // Determine the interface points.
	    p2 = SolnBlk.Interface_List[ni].Spline.Xp[np  ];
	    p3 = SolnBlk.Interface_List[ni].Spline.Xp[np+1];
	    if (abs(p2 - p3) < epsilon) break;

// 	    // Determine the front speeds at the interface points.
// 	    f2 = SolnBlk.Interface_List[ni].Fn(np);
// 	    f3 = SolnBlk.Interface_List[ni].Fn(np+1);

	    // Calculate the lengths of the sides of the triangle.
	    d12 = abs(p1-p2);
	    d23 = abs(p2-p3);
	    d31 = abs(p3-p1);

	    // If the two angles of the triangle located opposite the node
	    // in question are acute or right angles then determine the 
	    // distance from the segment to the node that is normal to the 
	    // interface segment.  Choose the minimum distance.  Also 
	    // track the nearest interface point in case no such minimum 
	    // distance can be determined.  The distance between that 
	    // interface point and the node will be used as the distance 
	    // value.
	    if (fabs(d12) < epsilon) {
	      // p1 = p2 which requires d = zero and set f to the interface 
	      // speed at p2.
	      d = ZERO;
// 	      f = f2;

	    } else if (fabs(d31) < epsilon) {
	      // p1 = p3 which requires d = zero and set f to the interface 
	      // speed at p3.
	      d = ZERO;
// 	      f = f3;

	    } else if ((d12*d12 + d23*d23 > d31*d31 || fabs(d12*d12 + d23*d23 - d31*d31) < epsilon) &&
		       (d23*d23 + d31*d31 > d12*d12 || fabs(d23*d23 + d31*d31 - d12*d12) < epsilon)) {
	      // The angles at p2 and p3 are either acute or right angles.
	      // Calculate the value of d by using the Pythagorean 
	      // theorem.  Determine f by taking a weighted average of the
	      // interface speeds at p2 and p3.
	      if (d12 <= d31) {
		assert(d12*d23 > ZERO);
		costheta = (d12*d12 + d23*d23 - d31*d31)/(TWO*d12*d23);
		if (fabs(ONE - costheta*costheta) < epsilon) dt = ZERO;
		else dt = d12*sqrt(ONE - costheta*costheta);
	      } else {
		assert(d23*d31 > ZERO);
		costheta = (d23*d23 + d31*d31 - d12*d12)/(TWO*d23*d31);
		if (fabs(ONE - costheta*costheta) < epsilon) dt = ZERO;
		else dt = d31*sqrt(ONE - costheta*costheta);
	      }
	      // Determine the speed value.
// 	      if (dt > ZERO) {
// 		// Determine the areas.
// 		A  = fabs(HALF*((p1^p2) + (p2^p3) + (p3^p1)));
// 		A2 = HALF*dt*sqrt(d12*d12 - dt*dt);
// 		A3 = fabs(A - A2);
// 		if (A3 < ZERO) { cout << endl << " -------------> A = " << A << " A2 = " << A2 << " A3 = " << A3; cout.flush(); }
// 		//assert(A3 > ZERO);
// 		ft = (A3/A)*f2 + (A2/A)*f3;
// 	      } else {
// 		assert(d23 > ZERO);
// 		ft = (d31/d23)*f2 + (d12/d23)*f3;
// 	      }
	      // Set d and f.
	      if (dt < d) {
		d = dt;
// 		f = ft;
	      }

	    } else {
	      // One of the angles is obtuse.  Set d to the distance to 
	      // the nearest interface node and set f to the 
	      // corresponding speed,
	      if (d12 <= d && d12 <= d31) {
		// d2 is the closest interface node.
		d = d12;
// 		f = f2;
	      } else if (d31 <= d && d31 <= d12) {
		//  d3 is the nearest interface node.
		d = d31;
// 		f = f3;
	      }

	    }

	  }

	  // Enforce the signed distance to zero if need be.
	  if (fabs(d) < epsilon) d = ZERO;

	  // Choose the minimum of the previously computed signed distance
	  // or the newly determined value.  This is important for when 
	  // multiple interfaces exist.  
	  if (d < fabs(SolnBlk.U[i][j].psi)) {
	    // Set the signed distance value, determine signs afterward.
	    SolnBlk.U[i][j].psi = d;
// 	    // Set the appropriate front speed.
// 	    SolnBlk.U[i][j].F = f;
// 	    if (fabs(SolnBlk.U[i][j].F) < epsilon) SolnBlk.U[i][j].F = ZERO;
	  }

	}

      }

      if (fabs(SolnBlk.U[i][j].psi) > IP.Extension_Distance) {
	SolnBlk.U[i][j].psi = IP.Extension_Distance;
// 	SolnBlk.U[i][j].F = ZERO;
      }

      // Determine the appropriate sign of the signed distance function
      // by counting the number of times it crosses the interface.  An
      // even number corresponds to a negative sign and a odd number 
      // corresponds to a positive sign.
      for (int ni = 1; ni <= SolnBlk.Interface_List.Ni; ni++) {
	if (SolnBlk.U[i][j].psi > ZERO) {
 	  if (SolnBlk.Interface_List[ni].Point_In_Interface(SolnBlk.Grid.Cell[i][j].Xc)) {
 	    SolnBlk.U[i][j].psi = -SolnBlk.U[i][j].psi;
 	  }
	}
      }

      // Special case for Zalesak's disk.
      if (SolnBlk.Interface_List[1].Type == INTERFACE_ZALESAK) {
	//if (SolnBlk.Interface_List[1].Point_In_Interface(SolnBlk.Grid.Cell[i][j].Xc))
	SolnBlk.U[i][j].psi = -SolnBlk.U[i][j].psi;
      }

      // Intentionally destroy the signed distance function to test the Eikonal equation.
      if (IP.Perturb_Distance_Function) SolnBlk.U[i][j].psi *= 0.1 + sqr(SolnBlk.Grid.Cell[i][j].Xc.x+3.5) + sqr(SolnBlk.Grid.Cell[i][j].Xc.y+2.0);

    }
  }

  // Perform the scalar geometric extension problem.
  error_flag = Scalar_Geometric_Extension_Problem(SolnBlk,IP);
  if (error_flag) return error_flag;

  // Extension problem solved successfully.
  return 0;

}

/******************************************************************//**
 * Routine: Scalar_Geometric_Extension_Problem                        
 *                                                                    
 * Geometric solution of the scalar extension problem (eg, the normal 
 * front speed).  This algorithm based on triangles constructed out   
 * of the point of contention and the two nearest interface nodes.    
 *                                                                    
 **********************************************************************/
int Scalar_Geometric_Extension_Problem(LevelSet2D_Quad_Block &SolnBlk,
				       LevelSet2D_Input_Parameters &IP) {

  double d = MILLION, h = MILLION;
  double d12, d23, d31, costheta, dt;
  double f, f2, f3, ft, A, A2, A3;
  Vector2D p1, p2, p3, norm_dir;
  double epsilon = TOLER*min(SolnBlk.Grid.lfaceN(SolnBlk.Grid.ICl,SolnBlk.Grid.JCl),
			     SolnBlk.Grid.lfaceE(SolnBlk.Grid.ICl,SolnBlk.Grid.JCl));

  // Initialize the signed distance function to a million and the speed
  // function to zero.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      SolnBlk.U[i][j].F = ZERO;
    }
  }

  // Ensure that the bounding box has been computed for all of the
  // interfaces.
  for (int ni = 1; ni <= SolnBlk.Interface_List.Ni; ni++) {
    SolnBlk.Interface_List[ni].BoundingBox(IP.Extension_Distance);
  }

  // Determine the signed distance at each cell center.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {

      d = MILLION; // acute/right calculation.
      h = MILLION; // obtuse calcuation.

      // Determine the local normal direction defined by the level set function.
      Linear_Reconstruction_LeastSquares(SolnBlk,i,j,LIMITER_ONE);
      if (sqrt(sqr(SolnBlk.dUdx[i][j].psi) + sqr(SolnBlk.dUdy[i][j].psi)) > TOLER*TOLER) {
	norm_dir = Vector2D(SolnBlk.dUdx[i][j].psi,SolnBlk.dUdy[i][j].psi)/
   	           sqrt(sqr(SolnBlk.dUdx[i][j].psi) + sqr(SolnBlk.dUdy[i][j].psi));
      } else {
	norm_dir = Vector2D_ZERO;
      }

      // Determine the mesh cell-center.
      p1 = SolnBlk.Grid.Cell[i][j].Xc;

      // The geometric extension problem is solved for each interface
      // and the solution for the minimum signed distance is used.
      for (int ni = 1; ni <= SolnBlk.Interface_List.Ni; ni++) {

	// Only solve the geometric extension if the current point, p, is
	// contained within the bounding box of the interface.
	if (p1.x >= SolnBlk.Interface_List[ni].Xmin.x &&
	    p1.y >= SolnBlk.Interface_List[ni].Xmin.y &&
	    p1.x <= SolnBlk.Interface_List[ni].Xmax.x &&
	    p1.y <= SolnBlk.Interface_List[ni].Xmax.y) {

	  // For each interface segment.
	  for (int np = 0; np < SolnBlk.Interface_List[ni].Spline.np-1; np++) {

	    // Determine the interface points.
	    p2 = SolnBlk.Interface_List[ni].Spline.Xp[np  ];
	    p3 = SolnBlk.Interface_List[ni].Spline.Xp[np+1];
	    if (abs(p2 - p3) < epsilon) break;

	    // Determine the front speeds at the interface points.
	    f2 = SolnBlk.Interface_List[ni].Speed.x;
	    f3 = SolnBlk.Interface_List[ni].Speed.x;

	    //  f2 = SolnBlk.Interface_List[ni].Fn(np);//,norm_dir);
	    //  f3 = SolnBlk.Interface_List[ni].Fn(np+1);//,norm_dir);

	    // Calculate the lengths of the sides of the triangle.
	    d12 = abs(p1-p2);
	    d23 = abs(p2-p3);
	    d31 = abs(p3-p1);

	    // If the two angles of the triangle located opposite the node
	    // in question are acute or right angles then	determine the 
	    // distance from the segment to the node that is normal to the 
	    // interface segment.  Choose the minimum distance.  Also 
	    // track the nearest interface point in case no such minimum 
	    // distance can be determined.  The distance between that 
	    // interface point and the node will be used as the distance 
	    // value.
	    if (fabs(d12) < epsilon) {
	      // p1 = p2 which requires d = zero and set f to the interface 
	      // speed at p2.
	      d = ZERO;
	      f = f2;

	    } else if (fabs(d31) < epsilon) {
	      // p1 = p3 which requires d = zero and set f to the inetrface 
	      // speed at p3.
	      d = ZERO;
	      f = f3;

	    } else if ((d12*d12 + d23*d23 > d31*d31 || fabs(d12*d12 + d23*d23 - d31*d31) < epsilon*d23) &&
		       (d23*d23 + d31*d31 > d12*d12 || fabs(d23*d23 + d31*d31 - d12*d12) < epsilon*d23)) {
	      // The angles at p2 and p3 are either acute or right angles.
	      // Calculate the value of d by using the Pythagorean 
	      // theorem.  Determine f by taking a weighted average of the
	      // interface speeds at p2 and p3.
	      if (d12 <= d31) {
		assert(d12*d23 > ZERO);
		costheta = (d12*d12 + d23*d23 - d31*d31)/(TWO*d12*d23);
		if (fabs(ONE - costheta*costheta) < epsilon) dt = ZERO;
		else dt = d12*sqrt(ONE - costheta*costheta);
	      } else {
		assert(d23*d31 > ZERO);
		costheta = (d23*d23 + d31*d31 - d12*d12)/(TWO*d23*d31);
		if (fabs(ONE - costheta*costheta) < epsilon) dt = ZERO;
		else dt = d31*sqrt(ONE - costheta*costheta);
	      }
	      // Determine the speed value.
	      if (dt > ZERO) {
		// Determine the areas.
		A  = fabs(HALF*((p1^p2) + (p2^p3) + (p3^p1)));
		A2 = HALF*dt*sqrt(d12*d12 - dt*dt);
		A3 = A - A2;
		assert(A > ZERO);
		ft = (A3/A)*f2 + (A2/A)*f3;
	      } else {
		assert(d23 > ZERO);
		ft = (d31/d23)*f2 + (d12/d23)*f3;
	      }
	      // Set d and f.
	      if (dt < d) {
		d = dt;
		f = ft;
	      }

	    } else {
	      // One of the angles is obtuse.  Set d to the distance to 
	      // the nearest interface node and set f to the 
	      // corresponding speed,
	      if (d12 <= d && d12 <= d31) {
		// d2 is the closest interface node.
		d = d12;
		f = f2;
	      } else if (d31 <= d && d31 <= d12) {
		//  d3 is the nearest interface node.
		d = d31;
		f = f3;
	      }

	    }

	  }

	}

      }

      // Enforce the signed distance to zero if need be.
      if (fabs(d) < epsilon) d = ZERO;
      if (fabs(f) < epsilon) f = ZERO;

      // Assign the solution to the scalar geometric extension problem.
      if (fabs(SolnBlk.U[i][j].psi) < IP.Extension_Distance) SolnBlk.U[i][j].F = f;
      else SolnBlk.U[i][j].F = ZERO;

    }
  }

  // Scalar geometric extension problem solved successfully.
  return 0;

}

/******************************************************************//**
 * Routine: Retrieve_Interface_Spline                                 
 *                                                                    
 * This routine locates the zero level set contained within the       
 * specified solution block. The interface location is saved as a     
 * spline(s).                                                         
 *                                                                    
 **********************************************************************/
#ifndef _RETRIEVE_DEBUG_
int Retrieve_Interface_Spline(LevelSet2D_Quad_Block &SolnBlk,
			      const int &gblknum) {
#endif
#ifdef _RETRIEVE_DEBUG_
int Retrieve_Interface_Spline(LevelSet2D_Quad_Block &SolnBlk,
			      const int &gblknum,
			      ofstream &dout) {
#endif

  int found;

  //  int numInfectedFaces;   // number of infected faces/centers

  int ncells;     // number of infected cells.

  double epsilon = TOLER*min(SolnBlk.Grid.lfaceN(SolnBlk.Grid.ICl,SolnBlk.Grid.JCl),
			     SolnBlk.Grid.lfaceE(SolnBlk.Grid.ICl,SolnBlk.Grid.JCl));
#ifdef _RETRIEVE_DEBUG_
  dout << endl << " epsilon = " << epsilon; dout.flush();
#endif

  // Deallocate current level set spline(s).
  if (SolnBlk.Interface_List.Ni) SolnBlk.Interface_List.deallocate();

  // Initialize cut type matrix.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      SolnBlk.cut[i][j].clean();
    }
  }

  // Flag infected cells and identify start points for spline tracing.
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      Flag_Infected_Cell(SolnBlk,i,j,ncells,epsilon);
    }
  }

  // Exit immediately if no interface points have been located.
  if (!ncells) return 0;

#ifdef _RETRIEVE_DEBUG_
  dout << endl << "Cells infected in CENTER:" << endl;
  for (int j = SolnBlk.JCu+1; j >= SolnBlk.JCl-1; j--) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      if (SolnBlk.cut[i][j].center == INFECTED) dout << "c ";
      else dout << "o ";
      dout.flush();
    }
    dout << endl;
  }
  dout << endl << "Cells infected in NORTH:" << endl;
  for (int j = SolnBlk.JCu+1; j >= SolnBlk.JCl-1; j--) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      if (SolnBlk.cut[i][j].north == INFECTED) dout << "n ";
      else dout << "o ";
      dout.flush();
    }
    dout << endl;
  }
  dout << endl << "Cells infected in SOUTH:" << endl;
  for (int j = SolnBlk.JCu+1; j >= SolnBlk.JCl-1; j--) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      if (SolnBlk.cut[i][j].south == INFECTED) dout << "s ";
      else dout << "o ";
      dout.flush();
    }
    dout << endl;
  }
  dout << endl << "Cells infected in EAST:" << endl;
  for (int j = SolnBlk.JCu+1; j >= SolnBlk.JCl-1; j--) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      if (SolnBlk.cut[i][j].east == INFECTED) dout << "e ";
      else dout << "o ";
      dout.flush();
    }
    dout << endl;
  }
  dout << endl << "Cells infected in WEST:" << endl;
  for (int j = SolnBlk.JCu+1; j >= SolnBlk.JCl-1; j--) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      if (SolnBlk.cut[i][j].west == INFECTED) dout << "w ";
      else dout << "o ";
      dout.flush();
    }
    dout << endl;
  }
  dout << endl << "Cells infected in NORTHEAST:" << endl;
  for (int j = SolnBlk.JCu+1; j >= SolnBlk.JCl-1; j--) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      if (SolnBlk.cut[i][j].north_east == INFECTED) dout << "a ";
      else dout << "o ";
      dout.flush();
    }
    dout << endl;
  }
  dout << endl << "Cells infected in SOUTHWEST:" << endl;
  for (int j = SolnBlk.JCu+1; j >= SolnBlk.JCl-1; j--) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      if (SolnBlk.cut[i][j].south_west == INFECTED) dout << "d ";
      else dout << "o ";
      dout.flush();
    }
    dout << endl;
  }
  dout << endl << "Cells that are infected somewhere:" << endl;
  for (int j = SolnBlk.JCu+1; j >= SolnBlk.JCl-1; j--) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      if (SolnBlk.cut[i][j].center == INFECTED ||
 	  SolnBlk.cut[i][j].north == INFECTED ||
 	  SolnBlk.cut[i][j].east  == INFECTED ||
 	  SolnBlk.cut[i][j].south == INFECTED ||
 	  SolnBlk.cut[i][j].west  == INFECTED ||
 	  SolnBlk.cut[i][j].north_east == INFECTED ||
 	  SolnBlk.cut[i][j].south_west == INFECTED) {
 	dout << "+ ";
      } else {
 	dout << "o ";
      }
      dout.flush();
    }
    dout << endl;
  }
  dout << endl;
  dout << endl;
  dout << endl << " Number of infected faces/centers = " << ncells;
  dout << endl << " Number of starting points        = " << SolnBlk.Trace.np;
  dout << endl;
  dout << endl << " Starting data: ";
//   for (int nspts = 0; nspts < start.np; nspts++) {
//     dout << endl << nspts
// 	 << " " << start[nspts]
// 	 << " " << start[nspts].i
// 	 << " " << start[nspts].j
// 	 << " " << start[nspts].face;
//   }
//   dout << endl << " Trace data: ";
  for (int nifs = 0; nifs < SolnBlk.Trace.np; nifs++) {
    dout << endl
	 << " " << nifs
	 << " " << SolnBlk.Trace[nifs].i
	 << " " << SolnBlk.Trace[nifs].j
	 << " " << SolnBlk.Trace[nifs].face
 	 << SolnBlk.Grid.Cell[SolnBlk.Trace[nifs].i][SolnBlk.Trace[nifs].j].Xc;
  }

  dout << endl;
#endif



  // Declare required integer variables.
  int error_flag;

  LinkedList<Vector2D> p,  // spline points for all splines
                       F;  // spline velocities for all splines

  LinkedList<int> npts;    // number of points in each spline

  int numpts;              // number of points in this spline
  int icell, jcell;
  int start = 0;           // Index of the 'Trace' LinkedList at which we should
                           // begin tracing the interface.

  int doneTracing = 0;     // completion flag

  // Trace all interfaces until there are no more starting points.
  while (!doneTracing) {

    // Trace the interface at specified start point.
#ifndef _RETRIEVE_DEBUG_
    error_flag = Trace_Interface_Spline(SolnBlk,p,F,epsilon,numpts,start);
#endif
#ifdef _RETRIEVE_DEBUG_
    error_flag = Trace_Interface_Spline(SolnBlk,p,F,epsilon,numpts,start,dout);
#endif
    if (error_flag) return error_flag;

    // Add to linked list.
    npts.add(numpts);

    // Next starting point index.
    start++;

    // Find the next infected starting point to begin tracing.
    while (start <= SolnBlk.Trace.np) {
      icell = SolnBlk.Trace[start].i;
      jcell = SolnBlk.Trace[start].j;
      if (SolnBlk.cut[icell][jcell].center == INFECTED ||
	  SolnBlk.cut[icell][jcell].north  == INFECTED ||
	  SolnBlk.cut[icell][jcell].east   == INFECTED ||
	  SolnBlk.cut[icell][jcell].south  == INFECTED ||
	  SolnBlk.cut[icell][jcell].west   == INFECTED) {
	break;
      } else {
	if (start == SolnBlk.Trace.np) {
	  doneTracing = 1;   // All starting points have been traced.
	  break;
	} else {
	  start++;
	}
      }
    }

  }

  // Look for any more infected cells that have not been traced.
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      if (SolnBlk.cut[i][j].center == INFECTED ||
	  SolnBlk.cut[i][j].north  == INFECTED ||
	  SolnBlk.cut[i][j].east   == INFECTED ||
	  SolnBlk.cut[i][j].south  == INFECTED ||
	  SolnBlk.cut[i][j].west   == INFECTED) {
	// Trace using this as the starting point.
#ifndef _RETRIEVE_DEBUG_
	error_flag = Trace_Interface_Spline(SolnBlk,p,F,epsilon,numpts,0);
#endif
#ifdef _RETRIEVE_DEBUG_
	error_flag = Trace_Interface_Spline(SolnBlk,p,F,epsilon,numpts,0,dout);
#endif
	if (error_flag) return error_flag;

	// Add to linked list.
	npts.add(numpts);
      }
    }
  }

#ifdef _RETRIEVE_DEBUG_
  dout << endl << "==========================================";
  dout << endl << "Number of interfaces traced = " << npts.np-1;
  for (int ni = 0; ni < npts.np; ni++) {
    dout << endl << "Interface " << ni;
    dout << endl << "Spline points:";
    for (int pts = 0; pts <= npts[ni]; pts++) {
      dout << endl << p[pts] << " " << F[pts];
    }
    dout << endl << "---------";
  }
#endif

  // Clip interfaces as necessary.



// //   for (int m = 0; m < ncells; m++) {
// //     if (((SolnBlk.Trace[m].i == SolnBlk.ICl-1 || SolnBlk.Trace[m].i == SolnBlk.ICu+1) && SolnBlk.Trace[m].face == NORTH && SolnBlk.cut[SolnBlk.Trace[m].i][SolnBlk.Trace[m].j].north == INFECTED) ||
// //   	((SolnBlk.Trace[m].j == SolnBlk.JCl-1 || SolnBlk.Trace[m].j == SolnBlk.JCu+1) && SolnBlk.Trace[m].face == EAST  && SolnBlk.cut[SolnBlk.Trace[m].i][SolnBlk.Trace[m].j].east  == INFECTED) ||
// //   	((SolnBlk.Trace[m].i == SolnBlk.ICl-1 || SolnBlk.Trace[m].i == SolnBlk.ICu+1) && SolnBlk.Trace[m].face == SOUTH && SolnBlk.cut[SolnBlk.Trace[m].i][SolnBlk.Trace[m].j].south == INFECTED) ||
// //   	((SolnBlk.Trace[m].j == SolnBlk.JCl-1 || SolnBlk.Trace[m].j == SolnBlk.JCu+1) && SolnBlk.Trace[m].face == WEST  && SolnBlk.cut[SolnBlk.Trace[m].i][SolnBlk.Trace[m].j].west  == INFECTED)) {
//   nsp++;
//   numpts = 0;

// #ifndef _RETRIEVE_DEBUG_
//   //error_flag = Trace_Interface_Spline(SolnBlk,p,F,numpts,m);
//   error_flag = Trace_Interface_Spline(SolnBlk,p,F,epsilon,numpts,start);
// #endif
// #ifdef _RETRIEVE_DEBUG_
//   dout << endl << " Begin Trace"; dout.flush();
//   //error_flag = Trace_Interface_Spline(SolnBlk,p,F,numpts,m,dout);
//   error_flag = Trace_Interface_Spline(SolnBlk,p,F,epsilon,numpts,start,dout);
//   dout << endl << "Number of points traced: " << numpts;
//   dout << endl << " End Trace"; dout.flush();
// #endif
//   if (error_flag) return error_flag;
//   npts.add(numpts);
// //     }
// //   }

//   // Exit immediately if no interfaces have been found.
//   if (!nsp) { 
//     p.deallocate();
//     F.deallocate();
//     npts.deallocate();
//     return 0;
//   }

//   // Parse the list of interfaces and clip interfaces as required.
//   LinkedList<Vector2D> pn, Fn;
//   LinkedList<int> nptsn;
//   int added_flag, in_counter;
//   Vector2D Xmin = SolnBlk.Grid.Node[SolnBlk.Grid.INl][SolnBlk.Grid.JNl].X;
//   Vector2D Xmax = SolnBlk.Grid.Node[SolnBlk.Grid.INu][SolnBlk.Grid.JNu].X;
//   Vector2D Xp, Xp1, Xp2, Xp3, Xp4, fm;

//   start = 0;

//   for (int ni = 0; ni < nsp; ni++) {
//     nptsn.add(0);
//     if (npts[ni] > 1) {
//       // If none of the traced nodes are within the block then disregard
//       // the trace.
//       in_counter = 0;
//       for (int np = start; np < npts[ni]; np++) {
//  	if (p[np].x >= Xmin.x && p[np].x <= Xmax.x && p[np].y >= Xmin.y && p[np].y <= Xmax.y) {
//  	  in_counter++;
//  	}
//       }
//       // None of the traced nodes are within the block but if the trace
//       // interesects the block (can occur at corners) then interpolate
//       // an internal point.
//       if (!in_counter) {
// 	for (int np = start; np < npts[ni]-1; np++) {
// 	  found = 0;
//  	  if ((p[np].x < Xmin.x || p[np].x > Xmax.x || p[np].y < Xmin.y || p[np].y > Xmax.y) &&
//  	      (p[np+1].x < Xmin.x || p[np+1].x > Xmax.x || p[np+1].y < Xmin.y || p[np+1].y > Xmax.y)) {
//   	    if (Line_Intersection(p[np],p[np+1],Vector2D(Xmin.x,Xmin.y),Vector2D(Xmax.x,Xmin.y),Xp1)) in_counter++;
//  	    if (Line_Intersection(p[np],p[np+1],Vector2D(Xmax.x,Xmin.y),Vector2D(Xmax.x,Xmax.y),Xp2)) in_counter++;
//  	    if (Line_Intersection(p[np],p[np+1],Vector2D(Xmax.x,Xmax.y),Vector2D(Xmin.x,Xmax.y),Xp3)) in_counter++;
//  	    if (Line_Intersection(p[np],p[np+1],Vector2D(Xmin.x,Xmax.y),Vector2D(Xmin.x,Xmin.y),Xp4)) in_counter++;
// 	    if (in_counter == 2) {
// 	      Xp = HALF*(Xp1 + Xp2 + Xp3 + Xp4);
//  	      fm = Linear_Interpolation(p[np],F[np],p[np+1],F[np+1],Xp);
// 	      p.insert(Xp,np);
// 	      F.insert(fm,np);
// 	      found = 1;
// 	    }
//  	  }
// 	  if (found) break;
// 	}
//       }
//       // Include the traced nodes.
//       if (in_counter) {
//         added_flag = 0;
// 	for (int np = start; np < npts[ni]; np++) {
// 	  if (!added_flag && np < npts[ni]-1) {
// 	    if ((p[np].x < Xmin.x || p[np].x > Xmax.x || p[np].y < Xmin.y || p[np].y > Xmax.y) &&
// 		(p[np+1].x >= Xmin.x && p[np+1].x <= Xmax.x && p[np+1].y >= Xmin.y && p[np+1].y <= Xmax.y)) {
// 	      added_flag = 1;
// 	    }
// 	  }
// 	  if (added_flag) {
// 	    pn.add(p[np]);
// 	    Fn.add(F[np]);
// 	    nptsn.put(nptsn[ni]+1,ni);
// 	    if (np+1 < npts[ni]) {
// 	      if ((p[np].x < Xmin.x || p[np].x > Xmax.x || p[np].y < Xmin.y || p[np].y > Xmax.y) &&
// 		  (p[np+1].x < Xmin.x || p[np+1].x > Xmax.x || p[np+1].y < Xmin.y || p[np+1].y > Xmax.y)) {
// 		break;
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//     start = npts[ni];
//   }

//   // Count the new number of interfaces (ignore all 1-point interfaces).
//   int nnsp = 0;
//   for (int ni = 0; ni < nsp; ni++) if (nptsn[ni] > 1) nnsp++;

//   cout << endl << "###########";
//   cout << endl << "New number of interfaces: " << nnsp;
//   cout << endl << "###########";

//   // Exit immediately if no interfaces exist.
//   if (!nnsp) {
//     p.deallocate(); F.deallocate(); npts.deallocate();
//     pn.deallocate(); Fn.deallocate(); nptsn.deallocate();
//     return 0;
//   }

//   // Assemble the interface list.
//   int nnii = 0; start = 0;
//   double fx;

//   // Allocate the number of interfaces.
//   SolnBlk.Interface_List.allocate(nnsp);
//   for (int ni = 0; ni < nsp; ni++) {
//     if (nptsn[ni] > 1) {
//       // Increment the interface number.
//       nnii++;
//       // Allocate interface.
//       SolnBlk.Interface_List[nnii].allocate(nptsn[ni]);
//       // Set the interface type to the global block number.
//       SolnBlk.Interface_List[nnii].Type = gblknum;
//       // Set interface spline type.
//       SolnBlk.Interface_List[nnii].Spline.settype(SPLINE2D_LINEAR);
//       // Set the interface data.
//       for (int np = 0; np < nptsn[ni]; np++) {
//  	SolnBlk.Interface_List[nnii].Spline.Xp[np] = pn[start + np];
//  	if ((np == 0) || (np == npts[ni-1]-1)) {
//  	  SolnBlk.Interface_List[nnii].Spline.tp[np] = SPLINE2D_POINT_SHARP_CORNER;
//  	} else {
//  	  SolnBlk.Interface_List[nnii].Spline.tp[np] = SPLINE2D_POINT_NORMAL;
//  	}
// 	SolnBlk.Interface_List[nnii].Spline.tp[np] = SPLINE2D_POINT_NORMAL;
//  	SolnBlk.Interface_List[nnii].Spline.bc[np] = BC_NONE;
//  	SolnBlk.Interface_List[nnii].F[np] = Fn[start + np];
//       }
//       SolnBlk.Interface_List[nnii].Spline.pathlength();
//     }
//     start = nptsn[ni];
//   }

//   // Clip interfaces to block boundaries.

//   int clipped_flag;

//   for (int ni = 1; ni <= SolnBlk.Interface_List.Ni; ni++) {

//     // Determine if the first spline point needs to be clipped.
//     if (SolnBlk.Interface_List[ni].Spline.Xp[0].x < Xmin.x ||
// 	SolnBlk.Interface_List[ni].Spline.Xp[0].x > Xmax.x ||
// 	SolnBlk.Interface_List[ni].Spline.Xp[0].y < Xmin.y ||
// 	SolnBlk.Interface_List[ni].Spline.Xp[0].y > Xmax.y) {
//       clipped_flag = 0;
//       // Try clipping against the north boundary:
//       if (SolnBlk.Interface_List[ni].Spline.Xp[0].y > Xmax.y)
// 	clipped_flag = Line_Intersection(Vector2D(Xmin.x,Xmax.y),Vector2D(Xmax.x,Xmax.y),
// 					 SolnBlk.Interface_List[ni].Spline.Xp[0],SolnBlk.Interface_List[ni].Spline.Xp[1],Xp);
//       // Try clipping against the south boundary if not already clipped:
//       if (SolnBlk.Interface_List[ni].Spline.Xp[0].y < Xmin.y && !clipped_flag)
// 	clipped_flag = Line_Intersection(Vector2D(Xmin.x,Xmin.y),Vector2D(Xmax.x,Xmin.y),
// 					 SolnBlk.Interface_List[ni].Spline.Xp[0],SolnBlk.Interface_List[ni].Spline.Xp[1],Xp);
//       // Try clipping against the east boundary if not already clipped:
//       if (SolnBlk.Interface_List[ni].Spline.Xp[0].x > Xmax.x && !clipped_flag)
// 	clipped_flag = Line_Intersection(Vector2D(Xmax.x,Xmin.y),Vector2D(Xmax.x,Xmax.y),
// 					 SolnBlk.Interface_List[ni].Spline.Xp[0],SolnBlk.Interface_List[ni].Spline.Xp[1],Xp);
//       // Try clipping against the west boundary if not already clipped:
//       if (SolnBlk.Interface_List[ni].Spline.Xp[0].x < Xmin.x && !clipped_flag)
// 	clipped_flag = Line_Intersection(Vector2D(Xmin.x,Xmin.y),Vector2D(Xmin.x,Xmax.y),
// 					 SolnBlk.Interface_List[ni].Spline.Xp[0],SolnBlk.Interface_List[ni].Spline.Xp[1],Xp);
//       if (!clipped_flag) {
// #ifdef _RETRIEVE_DEBUG_
// 	dout << endl << " " << 7771 << SolnBlk.Interface_List[ni].Spline.Xp[0]
// 	     << SolnBlk.Interface_List[ni].Spline.Xp[1] << Xmin << Xmax;
// #endif
// 	return 7771;
//       }
//       // Determine the front speed at the clipped point if found
//       // and add the information to the interface list:
//       fm = Linear_Interpolation(SolnBlk.Interface_List[ni].Spline.Xp[0],SolnBlk.Interface_List[ni].F[0],
// 				SolnBlk.Interface_List[ni].Spline.Xp[1],SolnBlk.Interface_List[ni].F[1],Xp);
//       SolnBlk.Interface_List[ni].Spline.Xp[0] = Xp;
//       SolnBlk.Interface_List[ni].F[0] = fm;
//     }
//     // Determine if the last spline point needs to be clipped.
//     for (int np = SolnBlk.Interface_List[ni].Spline.np-1; np >= SolnBlk.Interface_List[ni].Spline.np-2; np--) {
//       if (SolnBlk.Interface_List[ni].Spline.Xp[np].x >= Xmin.x &&
// 	  SolnBlk.Interface_List[ni].Spline.Xp[np].x <= Xmax.x &&
// 	  SolnBlk.Interface_List[ni].Spline.Xp[np].y >= Xmin.y &&
// 	  SolnBlk.Interface_List[ni].Spline.Xp[np].y <= Xmax.y) {
// 	// Do nothing.
// 	break;
//       } else if ((SolnBlk.Interface_List[ni].Spline.Xp[np].x < Xmin.x ||
// 		  SolnBlk.Interface_List[ni].Spline.Xp[np].x > Xmax.x ||
// 		  SolnBlk.Interface_List[ni].Spline.Xp[np].y < Xmin.y ||
// 		  SolnBlk.Interface_List[ni].Spline.Xp[np].y > Xmax.y) &&
// 		 (SolnBlk.Interface_List[ni].Spline.Xp[np-1].x >= Xmin.x &&
// 		  SolnBlk.Interface_List[ni].Spline.Xp[np-1].x <= Xmax.x &&
// 		  SolnBlk.Interface_List[ni].Spline.Xp[np-1].y >= Xmin.y &&
// 		  SolnBlk.Interface_List[ni].Spline.Xp[np-1].y <= Xmax.y)) {
// 	clipped_flag = 0;
// 	// Try clipping against the north boundary:
// 	if (SolnBlk.Interface_List[ni].Spline.Xp[np].y > Xmax.y)
// 	  clipped_flag = Line_Intersection(Vector2D(Xmin.x,Xmax.y),Vector2D(Xmax.x,Xmax.y),
// 					   SolnBlk.Interface_List[ni].Spline.Xp[np],SolnBlk.Interface_List[ni].Spline.Xp[np-1],Xp);
// 	// Try clipping against the south boundary if not already clipped:
// 	if (SolnBlk.Interface_List[ni].Spline.Xp[np].y < Xmin.y && !clipped_flag)
// 	  clipped_flag = Line_Intersection(Vector2D(Xmin.x,Xmin.y),Vector2D(Xmax.x,Xmin.y),
// 					   SolnBlk.Interface_List[ni].Spline.Xp[np],SolnBlk.Interface_List[ni].Spline.Xp[np-1],Xp);
// 	// Try clipping against the east boundary if not already clipped:
// 	if (SolnBlk.Interface_List[ni].Spline.Xp[np].x > Xmax.x && !clipped_flag)
// 	  clipped_flag = Line_Intersection(Vector2D(Xmax.x,Xmin.y),Vector2D(Xmax.x,Xmax.y),
// 					   SolnBlk.Interface_List[ni].Spline.Xp[np],SolnBlk.Interface_List[ni].Spline.Xp[np-1],Xp);
// 	// Try clipping against the west boundary if not already clipped:
// 	if (SolnBlk.Interface_List[ni].Spline.Xp[np].x < Xmin.x && !clipped_flag)
// 	  clipped_flag = Line_Intersection(Vector2D(Xmin.x,Xmin.y),Vector2D(Xmin.x,Xmax.y),
// 					   SolnBlk.Interface_List[ni].Spline.Xp[np],SolnBlk.Interface_List[ni].Spline.Xp[np-1],Xp);
// 	if (!clipped_flag) {
// #ifdef _RETRIEVE_DEBUG_
// 	  dout << endl << " " << 7772;
// #endif
// 	  return 7772;
// 	}
// 	// Determine the front speed at the clipped point if found
// 	// and add the information to the interface list:
// 	fm = Linear_Interpolation(SolnBlk.Interface_List[ni].Spline.Xp[np],SolnBlk.Interface_List[ni].F[np],
// 				  SolnBlk.Interface_List[ni].Spline.Xp[np-1],SolnBlk.Interface_List[ni].F[np-1],Xp);
// 	SolnBlk.Interface_List[ni].Spline.Xp[np] = Xp;
// 	SolnBlk.Interface_List[ni].F[np] = fm;
// 	break;
//       } else {
// 	SolnBlk.Interface_List[ni].Spline.np--;
//       }
//     }
//   }

//   p.deallocate(); F.deallocate(); npts.deallocate();
//   pn.deallocate(); Fn.deallocate(); nptsn.deallocate();

  // Interface spline retrieval completed.
  return 0;

}

/******************************************************************//**
 * Routine: Trace_Interface_Spline                                    
 *                                                                    
 * This routine will trace the interface throughout a domain given a  
 * starting location.  The interface is saved as a linked list of     
 * spline points.                                                     
 *                                                                    
 **********************************************************************/
#ifndef _RETRIEVE_DEBUG_
int Trace_Interface_Spline(LevelSet2D_Quad_Block &SolnBlk,
			   LinkedList<Vector2D> &p,
			   LinkedList<Vector2D> &F,
			   const double &epsilon,
			   int &numpts,
			   const int &start) {
#endif
#ifdef _RETRIEVE_DEBUG_
int Trace_Interface_Spline(LevelSet2D_Quad_Block &SolnBlk,
			   LinkedList<Vector2D> &p,
			   LinkedList<Vector2D> &F,
			   const double &epsilon,
			   int &numpts,
			   const int &start,
			   ofstream &dout) {
#endif

  int done = 0;
  int face = SolnBlk.Trace[start].face, faceo = face;
  int ic = SolnBlk.Trace[start].i, ico = ic;
  int jc = SolnBlk.Trace[start].j, jco = jc;
  double fn = ZERO;
  Vector2D pn;
  Vector2D Xmin = SolnBlk.Grid.Node[SolnBlk.Grid.INl][SolnBlk.Grid.JNl].X;
  Vector2D Xmax = SolnBlk.Grid.Node[SolnBlk.Grid.INu][SolnBlk.Grid.JNu].X;

  // Clean the infected cell.
  Clean_Infected_Cell(SolnBlk,ic,jc,face);

  // Trace spline.
  while (!done) {

    ico = ic;  jco = jc;  faceo = face;

    // Determine interface spline point.
    if (face == CENTER) {
      pn = SolnBlk.Grid.Cell[ico][jco].Xc;
      fn = SolnBlk.U[ico][jco].F;
    } else if (face == NORTH) {
      pn = Interpolate(SolnBlk.Grid.Cell[ico][jco  ].Xc,SolnBlk.U[ico][jco  ].psi,
		       SolnBlk.Grid.Cell[ico][jco+1].Xc,SolnBlk.U[ico][jco+1].psi,
		       epsilon);
      fn = Linear_Interpolation(SolnBlk.Grid.Cell[ico][jco  ].Xc,SolnBlk.U[ico][jco  ].F,
				SolnBlk.Grid.Cell[ico][jco+1].Xc,SolnBlk.U[ico][jco+1].F,pn);
    } else if (face == SOUTH) {
      pn = Interpolate(SolnBlk.Grid.Cell[ico][jco  ].Xc,SolnBlk.U[ico][jco  ].psi,
		       SolnBlk.Grid.Cell[ico][jco-1].Xc,SolnBlk.U[ico][jco-1].psi,
		       epsilon);
      fn = Linear_Interpolation(SolnBlk.Grid.Cell[ico][jco  ].Xc,SolnBlk.U[ico][jco  ].F,
				SolnBlk.Grid.Cell[ico][jco-1].Xc,SolnBlk.U[ico][jco-1].F,pn);
    } else if (face == EAST) {
      pn = Interpolate(SolnBlk.Grid.Cell[ico  ][jco].Xc,SolnBlk.U[ico  ][jco].psi,
		       SolnBlk.Grid.Cell[ico+1][jco].Xc,SolnBlk.U[ico+1][jco].psi,
		       epsilon);
      fn = Linear_Interpolation(SolnBlk.Grid.Cell[ico  ][jco].Xc,SolnBlk.U[ico  ][jco].F,
				SolnBlk.Grid.Cell[ico+1][jco].Xc,SolnBlk.U[ico+1][jco].F,pn);
    } else if (face == WEST) {
      pn = Interpolate(SolnBlk.Grid.Cell[ico  ][jco].Xc,SolnBlk.U[ico  ][jco].psi,
		       SolnBlk.Grid.Cell[ico-1][jco].Xc,SolnBlk.U[ico-1][jco].psi,
		       epsilon);
      fn = Linear_Interpolation(SolnBlk.Grid.Cell[ico  ][jco].Xc,SolnBlk.U[ico  ][jco].F,
				SolnBlk.Grid.Cell[ico-1][jco].Xc,SolnBlk.U[ico-1][jco].F,pn);
    } else if (face == NORTH_EAST) {
      pn = Interpolate(SolnBlk.Grid.Cell[ico+1][jco].Xc,SolnBlk.U[ico+1][jco].psi,
		       SolnBlk.Grid.Cell[ico][jco+1].Xc,SolnBlk.U[ico][jco+1].psi,
		       epsilon);
      fn = Linear_Interpolation(SolnBlk.Grid.Cell[ico+1][jco].Xc,SolnBlk.U[ico+1][jco].F,
				SolnBlk.Grid.Cell[ico][jco+1].Xc,SolnBlk.U[ico][jco+1].F,pn);
    } else if (face == SOUTH_WEST) {
      pn = Interpolate(SolnBlk.Grid.Cell[ico-1][jco].Xc,SolnBlk.U[ico-1][jco].psi,
		       SolnBlk.Grid.Cell[ico][jco-1].Xc,SolnBlk.U[ico][jco-1].psi,
		       epsilon);
      fn = Linear_Interpolation(SolnBlk.Grid.Cell[ico-1][jco].Xc,SolnBlk.U[ico-1][jco].F,
				SolnBlk.Grid.Cell[ico][jco-1].Xc,SolnBlk.U[ico][jco-1].F,pn);
    }

#ifdef _RETRIEVE_DEBUG_
//     dout << endl << ic << " " << jc << " " << face << SolnBlk.Grid.Cell[ic][jc].Xc;
//     if (face == CENTER) {
//       dout << endl << " C:" << SolnBlk.Grid.Cell[ico][jco].Xc << " " << SolnBlk.U[ico][jco].psi;
//     } else if (face == NORTH) {
//       dout << endl << " N:" << SolnBlk.Grid.Cell[ico][jco].Xc << SolnBlk.Grid.Cell[ico][jco+1].Xc << " " << SolnBlk.U[ico][jco  ].psi << " " << SolnBlk.U[ico][jco+1].psi;
//     } else if (face == SOUTH) {
//       dout << endl << " S:" << SolnBlk.Grid.Cell[ico][jco].Xc << SolnBlk.Grid.Cell[ico][jco-1].Xc << " " << SolnBlk.U[ico][jco  ].psi << " " << SolnBlk.U[ico][jco-1].psi;
//     } else if (face == EAST) {
//       dout << endl << " E:" << SolnBlk.Grid.Cell[ico][jco].Xc << SolnBlk.Grid.Cell[ico+1][jco].Xc << " " << SolnBlk.U[ico][jco  ].psi << " " << SolnBlk.U[ico+1][jco].psi;
//     } else if (face == WEST) {
//       dout << endl << " W:" << SolnBlk.Grid.Cell[ico][jco].Xc << SolnBlk.Grid.Cell[ico-1][jco].Xc << " " << SolnBlk.U[ico][jco  ].psi << " " << SolnBlk.U[ico-1][jco].psi;
//     } else if (face == NORTH_EAST) {
//       dout << endl << " NE:" << SolnBlk.Grid.Cell[ico+1][jco].Xc << SolnBlk.Grid.Cell[ico][jco+1].Xc << " " << SolnBlk.U[ico+1][jco].psi << " " << SolnBlk.U[ico][jco+1].psi;
//     } else if (face == SOUTH_WEST) {
//       dout << endl << " SW:" << SolnBlk.Grid.Cell[ico-1][jco].Xc << SolnBlk.Grid.Cell[ico][jco-1].Xc << " " << SolnBlk.U[ico-1][jco].psi << " " << SolnBlk.U[ico][jco-1].psi;
//     }
    dout << endl << pn;// << " " << fn;
#endif

    // Add the computed point to the list.
    if (face != SOUTH_WEST && face != NORTH_EAST) {
      p.add(pn); F.add(Vector2D(fn,fn)); numpts++;
    }

    // Search for next face.
    if (face == CENTER) {
      // CENTRE.
      if (SolnBlk.cut[ico][jco].north == INFECTED) {
	ic = ico; jc = jco; face = NORTH;
      } else if (SolnBlk.cut[ico][jco+1].center == INFECTED) {
	ic = ico; jc = jco+1; face = CENTER;
      } else if (SolnBlk.cut[ico][jco].north_east == INFECTED) {
	ic = ico; jc = jco; face = NORTH_EAST;
      } else if (SolnBlk.cut[ico][jco].east == INFECTED) {
	ic = ico; jc = jco; face = EAST;
      } else if (SolnBlk.cut[ico+1][jco].center == INFECTED) {
	ic = ico+1; jc = jco; face = CENTER;
      } else if (SolnBlk.cut[ico+1][jco].south == INFECTED) {
	ic = ico+1; jc = jco; face = SOUTH;
      } else if (SolnBlk.cut[ico][jco-1].north_east == INFECTED) {
	ic = ico; jc = jco-1; face = NORTH_EAST;
      } else if (SolnBlk.cut[ico][jco-1].east == INFECTED) {
	ic = ico; jc = jco-1; face = EAST;
      } else if (SolnBlk.cut[ico][jco-1].north == INFECTED) {
	ic = ico; jc = jco-1; face = NORTH;
      } else if (SolnBlk.cut[ico][jco-1].center == INFECTED) {
	ic = ico; jc = jco-1; face = CENTER;
      } else if (SolnBlk.cut[ico][jco].south_west == INFECTED) {
	ic = ico; jc = jco; face = SOUTH_WEST;
      } else if (SolnBlk.cut[ico][jco].west == INFECTED) {
	ic = ico; jc = jco; face = WEST;
      } else if (SolnBlk.cut[ico-1][jco].center == INFECTED) {
	ic = ico-1; jc = jco; face = CENTER;
      } else if (SolnBlk.cut[ico-1][jco].north == INFECTED) {
	ic = ico-1; jc = jco; face = NORTH;
      } else if (SolnBlk.cut[ico-1][jco].north_east == INFECTED) {
	ic = ico-1; jc = jco; face = NORTH_EAST;
      } else if (SolnBlk.cut[ico][jco+1].west == INFECTED) {
	ic = ico; jc = jco+1; face = WEST;
      }
    } else if (face == NORTH) {
      // NORTH face.
      if (SolnBlk.cut[ico][jco].east == INFECTED) {
	ic = ico; jc = jco; face = EAST;
      } else if (SolnBlk.cut[ico][jco].north_east == INFECTED) {
	ic = ico; jc = jco; face = NORTH_EAST;
      } else if (SolnBlk.cut[ico][jco+1].west == INFECTED) {
	ic = ico; jc = jco+1; face = WEST;
      } else if (SolnBlk.cut[ico][jco+1].south_west == INFECTED) {
	ic = ico; jc = jco+1; face = SOUTH_WEST;
      } else if (SolnBlk.cut[ico+1][jco].center == INFECTED) {
	ic = ico+1; jc = jco; face = CENTER;
      } else if (SolnBlk.cut[ico-1][jco+1].center == INFECTED) {
	ic = ico-1; jc = jco+1; face = CENTER;
      }
    } else if (face == SOUTH) {
      // SOUTH face.
      if (SolnBlk.cut[ico][jco].west == INFECTED) {
	ic = ico; jc = jco; face = WEST;
      } else if (SolnBlk.cut[ico][jco].south_west == INFECTED) {
	ic = ico; jc = jco; face = SOUTH_WEST;
      } else if (SolnBlk.cut[ico][jco-1].east == INFECTED) {
	ic = ico; jc = jco-1; face = EAST;
      } else if (SolnBlk.cut[ico][jco-1].north_east == INFECTED) {
	ic = ico; jc = jco-1; face = NORTH_EAST;
      } else if (SolnBlk.cut[ico+1][jco-1].center == INFECTED) {
	ic = ico+1; jc = jco-1; face = CENTER;
      } else if (SolnBlk.cut[ico-1][jco].center == INFECTED) {
	ic = ico-1; jc = jco; face = CENTER;
      }
    } else if (face == EAST) {
      // EAST face.
      if (SolnBlk.cut[ico][jco].north == INFECTED) {
	ic = ico; jc = jco; face = NORTH;
      } else if (SolnBlk.cut[ico][jco].north_east == INFECTED) {
	ic = ico; jc = jco; face = NORTH_EAST;
      } else if (SolnBlk.cut[ico+1][jco].south == INFECTED) {
	ic = ico+1; jc = jco; face = SOUTH;
      } else if (SolnBlk.cut[ico+1][jco].south_west == INFECTED) {
	ic = ico+1; jc = jco; face = SOUTH_WEST;
      } else if (SolnBlk.cut[ico][jco+1].center == INFECTED) {
	ic = ico; jc = jco+1; face = CENTER;
      } else if (SolnBlk.cut[ico+1][jco-1].center == INFECTED) {
	ic = ico+1; jc = jco-1; face = CENTER;
      }
    } else if (face == WEST) {
      // WEST face.
      if (SolnBlk.cut[ico][jco].south == INFECTED) {
	ic = ico; jc = jco; face = SOUTH;
      } else if (SolnBlk.cut[ico][jco].south_west == INFECTED) {
	ic = ico; jc = jco; face = SOUTH_WEST;
      } else if (SolnBlk.cut[ico-1][jco].north == INFECTED) {
	ic = ico-1; jc = jco; face = NORTH;
      } else if (SolnBlk.cut[ico-1][jco].north_east == INFECTED) {
	ic = ico-1; jc = jco; face = NORTH_EAST;
      } else if (SolnBlk.cut[ico][jco-1].center == INFECTED) {
	ic = ico; jc = jco-1; face = CENTER;
      } else if (SolnBlk.cut[ico-1][jco+1].center == INFECTED) {
	ic = ico-1; jc = jco+1; face = CENTER;
      }
    } else if (face == NORTH_EAST) {
      // NORTH-EAST face.
      if (SolnBlk.cut[ico][jco].north == INFECTED) {
	ic = ico; jc = jco; face = NORTH;
      } else if (SolnBlk.cut[ico][jco].east == INFECTED) {
	ic = ico; jc = jco; face = EAST;
      } else if (SolnBlk.cut[ico+1][jco+1].south == INFECTED) {
	ic = ico+1; jc = jco+1; face = SOUTH;
      } else if (SolnBlk.cut[ico+1][jco+1].west == INFECTED) {
	ic = ico+1; jc = jco+1; face = WEST;
      } else if (SolnBlk.cut[ico+1][jco+1].center == INFECTED) {
	ic = ico+1; jc = jco+1; face = CENTER;
      }
    } else if (face == SOUTH_WEST) {
      // SOUTH-WEST face.
      if (SolnBlk.cut[ico][jco].south == INFECTED) {
	ic = ico; jc = jco; face = SOUTH;
      } else if (SolnBlk.cut[ico][jco].west == INFECTED) {
	ic = ico; jc = jco; face = WEST;
      } else if (SolnBlk.cut[ico-1][jco-1].north == INFECTED) {
	ic = ico-1; jc = jco-1; face = NORTH;
      } else if (SolnBlk.cut[ico-1][jco-1].east == INFECTED) {
	ic = ico-1; jc = jco-1; face = EAST;
      } else if (SolnBlk.cut[ico-1][jco-1].center == INFECTED) {
	ic = ico-1; jc = jco-1; face = CENTER;
      }
    }

    // Now clean new face.
    Clean_Infected_Cell(SolnBlk,ic,jc,face);

    // If no new point has been found (the indices have not changed)
    // then determine exit criteria.
    if (ic == ico && jc == jco && face == faceo) {
//       // Close the curve if the starting point can be associated with
//       // the end point.
// //       if ((ic > SolnBlk.ICl || (ic == SolnBlk.ICl && face != WEST)) &&
// // 	  (ic < SolnBlk.ICu || (ic == SolnBlk.ICu && face != EAST)) &&
// // 	  (jc > SolnBlk.JCl || (jc == SolnBlk.JCl && face != SOUTH)) &&
// // 	  (jc < SolnBlk.JCu || (jc == SolnBlk.JCu && face != NORTH))) {
// // 	p.add(p[0]); F.add(F[0]); npts++;
// //  	dout << endl << " -> " << p[0];
// //       }
      done = 1;
    }

  }

  // Spline trace succesful.
  return 0;

}


/******************************************************************//**
 * Routine: Flag_Infected_Cell                                        
 *                                                                    
 * Flag infected cells and compile a list of all faces that the       
 * interfaces intersects with. Also identifies possible starting points
 * for the Trace_Interface_Spline routine.                            
 **********************************************************************/
void Flag_Infected_Cell(LevelSet2D_Quad_Block &SolnBlk,
			const int &ic,
			const int &jc,
			int &ncells,
			const double &epsilon) {

  if (fabs(SolnBlk.U[ic][jc].psi) < epsilon) {

    // Flag cell CENTER.
    SolnBlk.cut[ic][jc].center = INFECTED;
    ncells++;
    //  SolnBlk.Trace.add(Trace_Data(ic,jc,CENTER));
    if (ic == SolnBlk.ICu+1 || ic == SolnBlk.ICl-1 || jc == SolnBlk.JCu+1 || jc == SolnBlk.JCl-1) {
      SolnBlk.Trace.add(Trace_Data(ic,jc,CENTER));
    }
  } else {

    // Flag NORTH face.
//     if (jc <= SolnBlk.JCu && (SolnBlk.U[ic][jc].psi*SolnBlk.U[ic][jc+1].psi <= ZERO ||
// 			      fabs(SolnBlk.U[ic][jc+1].psi) < epsilon)) {	
    if (jc <= SolnBlk.JCu && SolnBlk.U[ic][jc].psi*SolnBlk.U[ic][jc+1].psi <= ZERO &&
	fabs(SolnBlk.U[ic][jc+1].psi) > epsilon) {	
      SolnBlk.cut[ic][jc].north = INFECTED;
      ncells++;
      //    SolnBlk.Trace.add(Trace_Data(ic,jc,NORTH));
      if (ic == SolnBlk.ICu+1) {
	SolnBlk.Trace.add(Trace_Data(ic,jc,NORTH));
      }
    }

    // Flag SOUTH face.
//     if (jc >= SolnBlk.JCl && (SolnBlk.U[ic][jc].psi*SolnBlk.U[ic][jc-1].psi <= ZERO ||
// 			      fabs(SolnBlk.U[ic][jc-1].psi) < epsilon)) {
    if (jc >= SolnBlk.JCl && SolnBlk.U[ic][jc].psi*SolnBlk.U[ic][jc-1].psi <= ZERO &&
	fabs(SolnBlk.U[ic][jc-1].psi) > epsilon) {
      SolnBlk.cut[ic][jc].south = INFECTED;
      ncells++;
      //    SolnBlk.Trace.add(Trace_Data(ic,jc,SOUTH));
      if (ic == SolnBlk.ICl-1) {
	SolnBlk.Trace.add(Trace_Data(ic,jc,SOUTH));
      }
    }

    // Flag EAST face.
//     if (ic <= SolnBlk.ICu && (SolnBlk.U[ic][jc].psi*SolnBlk.U[ic+1][jc].psi <= ZERO ||
// 			      fabs(SolnBlk.U[ic+1][jc].psi) < epsilon)) {
    if (ic <= SolnBlk.ICu && SolnBlk.U[ic][jc].psi*SolnBlk.U[ic+1][jc].psi <= ZERO &&
	fabs(SolnBlk.U[ic+1][jc].psi) > epsilon) {
      SolnBlk.cut[ic][jc].east = INFECTED;
      ncells++;
      //      SolnBlk.Trace.add(Trace_Data(ic,jc,EAST));
      if (jc == SolnBlk.JCl-1) {
	SolnBlk.Trace.add(Trace_Data(ic,jc,EAST));
      }
    }

    // Flag WEST face.
//     if (ic >= SolnBlk.ICl && (SolnBlk.U[ic][jc].psi*SolnBlk.U[ic-1][jc].psi <= ZERO ||
// 			      fabs(SolnBlk.U[ic-1][jc].psi) < epsilon)) {
    if (ic >= SolnBlk.ICl && SolnBlk.U[ic][jc].psi*SolnBlk.U[ic-1][jc].psi <= ZERO &&
	fabs(SolnBlk.U[ic-1][jc].psi) > epsilon) {
      SolnBlk.cut[ic][jc].west = INFECTED;
      ncells++;
      //      SolnBlk.Trace.add(Trace_Data(ic,jc,WEST));
      if (jc == SolnBlk.JCu+1) {
	SolnBlk.Trace.add(Trace_Data(ic,jc,WEST));
      }
    }

    // Flag NORTH-EAST face.
    if (ic <= SolnBlk.ICu && jc <= SolnBlk.JCu &&
	fabs(SolnBlk.U[ic+1][jc].psi) > epsilon && fabs(SolnBlk.U[ic][jc+1].psi) > epsilon &&
	SolnBlk.U[ic+1][jc].psi*SolnBlk.U[ic][jc+1].psi <= ZERO) {
      SolnBlk.cut[ic][jc].north_east = INFECTED;
      ncells++;
      //      SolnBlk.Trace.add(Trace_Data(ic,jc,NORTH_EAST));
    }

    // Flag SOUTH-WEST face.
    if (ic >= SolnBlk.ICl && jc >= SolnBlk.JCl &&
	fabs(SolnBlk.U[ic-1][jc].psi) > epsilon && fabs(SolnBlk.U[ic][jc-1].psi) > epsilon &&
	SolnBlk.U[ic-1][jc].psi*SolnBlk.U[ic][jc-1].psi <= ZERO) {
      SolnBlk.cut[ic][jc].south_west = INFECTED;
      ncells++;
      //      SolnBlk.Trace.add(Trace_Data(ic,jc,SOUTH_WEST));
    }

  }

}

/******************************************************************//**
 * Routine: Clean_Infected_Cell                                       
 *                                                                    
 * This routine turns off the appropriate flags of the adjacent cells 
 * for nodes that have been used to locate an interface point.        
 *                                                                    
 **********************************************************************/
void Clean_Infected_Cell(LevelSet2D_Quad_Block &SolnBlk,
			 const int &ic,
			 const int &jc,
			 const int &face) {

  if (face == CENTER) {
    SolnBlk.cut[ic][jc  ].center = CLEAN;
  } else if (face == NORTH) {
    SolnBlk.cut[ic][jc  ].north = CLEAN;
    SolnBlk.cut[ic][jc+1].south = CLEAN;
  } else if (face == SOUTH) {
    SolnBlk.cut[ic][jc  ].south = CLEAN;
    SolnBlk.cut[ic][jc-1].north = CLEAN;
  } else if (face == EAST) {
    SolnBlk.cut[ic  ][jc].east = CLEAN;
    SolnBlk.cut[ic+1][jc].west = CLEAN;
  } else if (face == WEST) {
    SolnBlk.cut[ic  ][jc].west = CLEAN;
    SolnBlk.cut[ic-1][jc].east = CLEAN;
  } else if (face == NORTH_EAST) {
    SolnBlk.cut[ic  ][jc  ].north_east = CLEAN;
    SolnBlk.cut[ic+1][jc+1].south_west = CLEAN;
  } else if (face == SOUTH_WEST) {
    SolnBlk.cut[ic  ][jc  ].south_west = CLEAN;
    SolnBlk.cut[ic-1][jc-1].north_east = CLEAN;
  }

}
