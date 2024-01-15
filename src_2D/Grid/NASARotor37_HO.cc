/*!\file NASARotor37_HO.cc
   \brief Header file implementing member functions of NASARotor37 class related to high-order grid type. */

#include "NASARotor37.h" /* Include NASARotor37 header file */


/**********************************************************************
 * Routine: genMeshH_3x2_AUTO                                         *
 *                                                                    *
 * Calls genMeshH_3x2 with preset parameters for high-order geometry. *
 *                                                                    *
 **********************************************************************/
void NASARotor37::genMeshH_3x2_AUTO(Grid2D_Quad_MultiBlock_HO &mesh,
				    double pspan, double zMin, double zMax,  
				    const int Number_of_Cells_Idir,
				    const int Number_of_Cells_Jdir,
				    const int Number_of_Ghost_Cells,
				    const int & Highest_Order_of_Reconstruction,
				    int n) {

  double stretch_i, stretch_j;
  
  cout << "\n  NASA Rotor 37 automatic 3x2 multi-block mesh generation...";

  //'hard-wire' stretching parameters
  if (Number_of_Cells_Idir/2 < 10) {
     stretch_i=1.05;
  } else if (Number_of_Cells_Idir/2 < 25) {
     stretch_i=1.01;
  } else if (Number_of_Cells_Idir/2 < 50) {
     stretch_i=1.005;
  } else if (Number_of_Cells_Idir/2 < 100) {
     stretch_i=1.0025;
  } else {
     stretch_i=1.001;
  } /* endif */
  if (Number_of_Cells_Jdir < 10) {
     stretch_j=1.05;
  } else if (Number_of_Cells_Jdir < 25) {
     stretch_j=1.01;
  } else if (Number_of_Cells_Jdir < 50) {
     stretch_j=1.005;
  } else if (Number_of_Cells_Jdir < 100) {
     stretch_j=1.0025;
  } else {
     stretch_j=1.001;
  } /* endif */

  //generate mesh
  genMeshH_3x2(mesh,
	       Highest_Order_of_Reconstruction,

	       pspan, zMin, zMax, n,
	   
	       Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	       STRETCHING_FCN_MAX_CLUSTERING, stretch_i, 0,
	       STRETCHING_FCN_MAX_CLUSTERING, stretch_j, 0,
	       BC_NONE, BC_NONE, BC_NONE, BC_FIXED,
	       ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, NOT_ORTHOGONAL, 
	   
	       Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	       STRETCHING_FCN_MINMAX_CLUSTERING, ONE+(stretch_i-ONE)/ONE, 0,
	       STRETCHING_FCN_MAX_CLUSTERING, stretch_j, 0,
	       BC_REFLECTION, BC_NONE, BC_NONE, BC_NONE,
	       ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, 
	   
	       Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	       STRETCHING_FCN_MIN_CLUSTERING, stretch_i, 0,
	       STRETCHING_FCN_MAX_CLUSTERING, stretch_j, 0,
	       BC_NONE, BC_CHARACTERISTIC, BC_NONE, BC_NONE,
	       NOT_ORTHOGONAL, NOT_ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, 
	   
	       Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	       STRETCHING_FCN_MAX_CLUSTERING, stretch_i, 0,
	       STRETCHING_FCN_MIN_CLUSTERING, stretch_j, 0,
	       BC_NONE, BC_NONE, BC_NONE, BC_FIXED,
	       ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, NOT_ORTHOGONAL, 
	   
	       Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	       STRETCHING_FCN_MINMAX_CLUSTERING, ONE+(stretch_i-ONE)/ONE, 0,
	       STRETCHING_FCN_MIN_CLUSTERING, stretch_j, 0,
	       BC_NONE, BC_NONE, BC_REFLECTION, BC_NONE,
	       ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, ORTHOGONAL, 
	   
	       Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Cells_Jdir,
	       STRETCHING_FCN_MIN_CLUSTERING, stretch_i, 0,
	       STRETCHING_FCN_MIN_CLUSTERING, stretch_j, 0,
	       BC_NONE, BC_CHARACTERISTIC, BC_NONE, BC_NONE,
	       ORTHOGONAL, NOT_ORTHOGONAL, ORTHOGONAL, ORTHOGONAL);   
}

/**********************************************************************
 * Routine: genMeshH_3x2                                              *
 *                                                                    *
 * -> generate a mesh for the cross-section of the blade at any %span *
 * -> H GRID (multi-block)                                            *
 * -> PARAMETERS MUST BE SUPPLIED FOR EVERY BLOCK                     *
 *                                                                    *
 * __________________________________________________________________ *
 * |                /                              \                | *
 * |               /                                \               | *
 * | Block (0,1)  |           Block (1,1)            |  Block (2,1) | *
 * |               \                                /               | *
 * |                \                              /                | *
 * |----------------- <<<<<<<<<<<<BLADE>>>>>>>>>>> -----------------| *
 * |                /                              \                | *
 * |               /                                \               | *
 * | Block (0,0)  |           Block (1,0)            |  Block (2,0) | *
 * |               \                                /               | *
 * |                \                              /                | *
 * |________________________________________________________________| *
 *                                                                    *
 **********************************************************************/
void NASARotor37::genMeshH_3x2(Grid2D_Quad_MultiBlock_HO &mesh,
			       const int & Highest_Order_of_Reconstruction,
			       
			       double pspan, double zMin, double zMax, int n, 
			       
			       //Parameters for block (0,0) - (BOTTOM LEFT CORNER)
			       const int B00_Number_of_Cells_Idir, const int B00_Number_of_Cells_Jdir,
			       const int B00_Number_of_Ghost_Cells,
			       int B00_StretchFcnI, double B00_BetaI, double B00_TauI,
			       int B00_StretchFcnJ, double B00_BetaJ, double B00_TauJ,
			       int B00_NbndType, int B00_EbndType, int B00_SbndType, int B00_WbndType,
			       int B00_OrthN, int B00_OrthE, int B00_OrthS, int B00_OrthW,
			       
			       //Parameters for block (1,0)
			       const int B10_Number_of_Cells_Idir, const int B10_Number_of_Cells_Jdir,
			       const int B10_Number_of_Ghost_Cells,
			       int B10_StretchFcnI, double B10_BetaI, double B10_TauI,
			       int B10_StretchFcnJ, double B10_BetaJ, double B10_TauJ,
			       int B10_NbndType, int B10_EbndType, int B10_SbndType, int B10_WbndType,
			       int B10_OrthN, int B10_OrthE, int B10_OrthS, int B10_OrthW,

			       //Parameters for block (2,0)
			       const int B20_Number_of_Cells_Idir, const int B20_Number_of_Cells_Jdir,
			       const int B20_Number_of_Ghost_Cells,
			       int B20_StretchFcnI, double B20_BetaI, double B20_TauI,
			       int B20_StretchFcnJ, double B20_BetaJ, double B20_TauJ,
			       int B20_NbndType, int B20_EbndType, int B20_SbndType, int B20_WbndType,
			       int B20_OrthN, int B20_OrthE, int B20_OrthS, int B20_OrthW,

			       //Parameters for block (0,1)
			       const int B01_Number_of_Cells_Idir, const int B01_Number_of_Cells_Jdir,
			       const int B01_Number_of_Ghost_Cells,
			       int B01_StretchFcnI, double B01_BetaI, double B01_TauI,
			       int B01_StretchFcnJ, double B01_BetaJ, double B01_TauJ,
			       int B01_NbndType, int B01_EbndType, int B01_SbndType, int B01_WbndType,
			       int B01_OrthN, int B01_OrthE, int B01_OrthS, int B01_OrthW,
			
			       //Parameters for block (1,1)
			       const int B11_Number_of_Cells_Idir, const int B11_Number_of_Cells_Jdir,
			       const int B11_Number_of_Ghost_Cells,
			       int B11_StretchFcnI, double B11_BetaI, double B11_TauI,
			       int B11_StretchFcnJ, double B11_BetaJ, double B11_TauJ,
			       int B11_NbndType, int B11_EbndType, int B11_SbndType, int B11_WbndType,
			       int B11_OrthN, int B11_OrthE, int B11_OrthS, int B11_OrthW,

			       //Parameters for block (2,1)
			       const int B21_Number_of_Cells_Idir, const int B21_Number_of_Cells_Jdir,
			       const int B21_Number_of_Ghost_Cells,
			       int B21_StretchFcnI, double B21_BetaI, double B21_TauI,
			       int B21_StretchFcnJ, double B21_BetaJ, double B21_TauJ,
			       int B21_NbndType, int B21_EbndType, int B21_SbndType, int B21_WbndType,
			       int B21_OrthN, int B21_OrthE, int B21_OrthS, int B21_OrthW) {


  Spline2D_HO upperB, lowerB, upperMiddleB, lowerMiddleB, camberTrail, camberLead, camberBlade, cs;
  int i, iTrail, iLead, pos, mIndex;
  Vector2D swap, leadV, trailV, x_temp;
  double zlU, zlL, zrU, zrL, mLead, mTrail, mLeadc, mTrailc, 
         mTop, mBot, m1, m2, z, dz, A, m, zLead, zTrail, dm, swapm;

  //*******************************************************************************
  //*******************************************************************************
  //*******************************  INITIALIZATION *******************************
  //*******************************************************************************
  //*******************************************************************************

  cout << "\n  Initializing Mesh Variables...";
  cout.flush();

  //Allocate memory for 3x2 multiblock grid.
  mesh.allocate(3,2);

  //find the index for the leading and trailing edge
  findLT(pspan, iTrail, iLead);

  //get cross section spline
  cs=getBladeCS(pspan);

  //re-organize pts so they ascend in the clockwise direction
  for ( i=0; i<R37GEOM_NUMPTS; ++i ) {
    swap=cs.Xp[2*R37GEOM_NUMPTS-2-i];
    cs.Xp[2*R37GEOM_NUMPTS-2-i]=cs.Xp[i];
    cs.Xp[i]=swap;
  }

  //change iTrail and iLead accordingly
  iTrail=cs.np-1-iTrail;
  iLead=cs.np-1-iLead;

  //save leading and trailing edge vectors
  leadV=cs.Xp[iLead];
  trailV=cs.Xp[iTrail];

  //get camber line splines - 40 pts
  camberTrail=getCamberTrail(pspan, 40, zMax); 
  camberLead=getCamberLead(pspan, 40, zMin);
  camberBlade=getCamberBlade(pspan, 40);
  
  //get upper, lower, upperMiddle, and lowerMiddle boundaries
  upperB=getCamberAndShift(pspan, 40, zMin, zMax, M_PI/num_blades);
  lowerB=getCamberAndShift(pspan, 40, zMin, zMax, -M_PI/num_blades);
  upperMiddleB=getCamberAndShift(pspan, 40, zMin, zMax, M_PI/(2*num_blades));
  lowerMiddleB=getCamberAndShift(pspan, 40, zMin, zMax, -M_PI/(2*num_blades));

  // Locate boundary corners to ensure that boundaries meeting at the leading
  // and trailing edges are separated by approximately 45deg
   
  //calculate slope of leading and trailing camber lines
  mLeadc=(camberLead.Xp[camberLead.np-1].y-camberLead.Xp[0].y)/
    (camberLead.Xp[camberLead.np-1].x-camberLead.Xp[0].x);

  mTrailc=(camberTrail.Xp[camberTrail.np-1].y-camberTrail.Xp[0].y)/
    (camberTrail.Xp[camberTrail.np-1].x-camberTrail.Xp[0].x);

  
  //calculate slope of tangent to cross-section at the leading and 
  //trailing edge using a central-difference
  z=cs.Xp[iLead].x;
  dz=(cs.Xp[iLead].x-cs.Xp[iLead+1].x)/10.0;
  mLead=(-getY(z+2*dz, cs)[0].y+8*getY(z+dz, cs)[0].y-8*getY(z-dz, cs)[0].y+
        getY(z-2*dz, cs)[0].y)/(12*dz);

  z=cs.Xp[iTrail].x;
  dz=(cs.Xp[iTrail].x-cs.Xp[iTrail-1].x)/10.0;
  mTrail=(-getY(z+2*dz, cs)[1].y+8*getY(z+dz, cs)[1].y-8*getY(z-dz, cs)[1].y+
         getY(z-2*dz, cs)[1].y)/(12*dz);

  //get slopes of lines approx 45deg to camber line at the LEADING EDGE
  A=(mLead+mLeadc)/(1-mLead*mLeadc);
  mTop=(-1+sqrt(1+A*A))/A;
  mBot=(-1-sqrt(1+A*A))/A;

  //find point on upperMiddle spline - leading edge
  //do search
  dm=1e15;
  mIndex=-1;
  for ( i=0; i<upperMiddleB.np; ++i ) {
    if (fabs(upperMiddleB.Xp[i].x-leadV.x) > TOLER) {
       m=(upperMiddleB.Xp[i].y-leadV.y)/(upperMiddleB.Xp[i].x-leadV.x);
      
       if( fabs(m-mTop)<dm) {
         dm=fabs(m-mTop);
         mIndex=i;
       } /* endif */
    } /* endif */
  }
  assert(mIndex!=-1);
  zlU=upperMiddleB.Xp[mIndex].x;  

  //find point on lowerMiddle spline - leading edge
  //do search
  dm=1e15;
  mIndex=-1;
  for ( i=0; i<lowerMiddleB.np; ++i ) {
    if (fabs(lowerMiddleB.Xp[i].x-leadV.x) > TOLER) {
       m=(lowerMiddleB.Xp[i].y-leadV.y)/(lowerMiddleB.Xp[i].x-leadV.x);
      
       if( fabs(m-mBot)<dm) {
         dm=fabs(m-mBot);
         mIndex=i;
       } /* endif */
    } /* endif */
  }
  assert(mIndex!=-1);
  zlL=lowerMiddleB.Xp[mIndex].x;  

  //get slopes of lines approx 45deg to camber line at the TRAILING EDGE
  A=(mTrail+mTrailc)/(1-mTrail*mTrailc);
  mBot=(-1+sqrt(1+A*A))/A;
  mTop=(-1-sqrt(1+A*A))/A;
  
  //find point on upperMiddle spline - trailing edge
  //do search
  dm=1e15;
  mIndex=-1;
  for ( i=0; i<upperMiddleB.np; ++i ) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
       m=(upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      
       if( fabs(m-mTop)<dm) {
         dm=fabs(m-mTop);
         mIndex=i;
       } /* endif */
    } /* endif */
  }
  assert(mIndex!=-1);
  zrU=upperMiddleB.Xp[mIndex].x;  

  dm=1e15;
  mIndex=-1;
  for ( i=0; i<upperMiddleB.np; ++i ) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
       m=(upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      
       if( fabs(m-mBot)<dm) {
         dm=fabs(m-mBot);
         mIndex=i;
       } /* endif */
    } /* endif */
  }
  assert(mIndex!=-1);
  
  //determine which to use
  if(upperMiddleB.Xp[mIndex].x>zrU) {
    zrU=upperMiddleB.Xp[mIndex].x;  
    swapm=mBot;
    mBot=mTop;
    mTop=swapm;
  };  

  //find point on lowerMiddle spline - trailing edge
  //do search
  dm=1e15;
  mIndex=-1;
  for ( i=0; i<lowerMiddleB.np; ++i ) {
    if (fabs(lowerMiddleB.Xp[i].x-trailV.x) > TOLER) {
       m=(lowerMiddleB.Xp[i].y-trailV.y)/(lowerMiddleB.Xp[i].x-trailV.x);
      
       if( fabs(m-mBot)<dm) {
         dm=fabs(m-mBot);
         mIndex=i;
       } /* endif */
    } /* endif */
  }
  assert(mIndex!=-1);
  zrL=lowerMiddleB.Xp[mIndex].x;  

  //*******************************************************************************
  //*******************************************************************************
  //****************************  GENERATE BOUNDARIES  ****************************
  //*******************************************************************************
  //*******************************************************************************  

  cout << "\n  Generating Boundaries...";
  cout.flush();

  // **********************
  // **** BLOCK (1,1) *****
  // **********************

  // ** NORTH SPLINE **
  //allocate space
  mesh(1,1).BndNorthSpline.allocate(70);

  //set type
  mesh(1,1).BndNorthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=(zlL+zlU)/2;
  zTrail=(zrL+zrU)/2;
  dz=(zTrail-zLead)/69;
  for ( i=0; i<70; ++i ) {
    mesh(1,1).BndNorthSpline.Xp[i]=getY(zLead+i*dz, upperB)[0];
  }

  //set point types
  mesh(1,1).BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,1).BndNorthSpline.tp[mesh(1,1).BndNorthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i< mesh(1,1).BndNorthSpline.np-1; ++i )
    mesh(1,1).BndNorthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh(1,1).BndNorthSpline.pathlength();

  // ** SOUTH SPLINE **
  //allocate space
  mesh(1,1).BndSouthSpline.allocate(iTrail-iLead+1);
  
  //set type
  mesh(1,1).BndSouthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  for ( i=iLead; i<=iTrail; ++i )
    mesh(1,1).BndSouthSpline.Xp[i-iLead]=cs.Xp[i];

  //set point types
  mesh(1,1).BndSouthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,1).BndSouthSpline.tp[mesh(1,1).BndSouthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh(1,1).BndSouthSpline.np-1; ++i )
    mesh(1,1).BndSouthSpline.tp[i]=SPLINE2D_POINT_NORMAL;
  
  //calc pathlength
  mesh(1,1).BndSouthSpline.pathlength();

  // ** EAST SPLINE **
  //allocate space
  mesh(1,1).BndEastSpline.allocate(5);
  
  //set type
  mesh(1,1).BndEastSpline.settype(SPLINE2D_CUBIC);

  //set points
  mesh(1,1).BndEastSpline.Xp[0]=cs.Xp[iTrail];
  x_temp = getY(zrU, upperMiddleB)[0];
  mesh(1,1).BndEastSpline.Xp[1]=cs.Xp[iTrail]+HALF*(x_temp-cs.Xp[iTrail]);
  mesh(1,1).BndEastSpline.Xp[2]=x_temp;
  mesh(1,1).BndEastSpline.Xp[3]=x_temp+
     HALF*(mesh(1,1).BndNorthSpline.Xp[mesh(1,1).BndNorthSpline.np-1]-x_temp);
  mesh(1,1).BndEastSpline.Xp[4]=mesh(1,1).BndNorthSpline.Xp[mesh(1,1).BndNorthSpline.np-1];
  
  //set point types
  mesh(1,1).BndEastSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,1).BndEastSpline.tp[1]=SPLINE2D_POINT_NORMAL;
  mesh(1,1).BndEastSpline.tp[2]=SPLINE2D_POINT_NORMAL;
  mesh(1,1).BndEastSpline.tp[3]=SPLINE2D_POINT_NORMAL;
  mesh(1,1).BndEastSpline.tp[4]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh(1,1).BndEastSpline.pathlength();

  // ** WEST SPLINE **
  //allocate space
  mesh(1,1).BndWestSpline.allocate(6);
  
  //set type
  mesh(1,1).BndWestSpline.settype(SPLINE2D_CUBIC);

  //set points
  mesh(1,1).BndWestSpline.Xp[0]=cs.Xp[iLead];
  x_temp = getY(zlU, upperMiddleB)[0];
  mesh(1,1).BndWestSpline.Xp[1]=cs.Xp[iLead]+0.3333*(x_temp-cs.Xp[iLead]);
  mesh(1,1).BndWestSpline.Xp[2]=HALF*(cs.Xp[iLead]+HALF*(x_temp-cs.Xp[iLead])+
                                       x_temp+0.20*(mesh(1,1).BndNorthSpline.Xp[0]-x_temp));
  mesh(1,1).BndWestSpline.Xp[3]=HALF*(cs.Xp[iLead]+0.6666*(x_temp-cs.Xp[iLead])+
                                       x_temp+0.3333*(mesh(1,1).BndNorthSpline.Xp[0]-x_temp));
  mesh(1,1).BndWestSpline.Xp[4]=x_temp+HALF*(mesh(1,1).BndNorthSpline.Xp[0]-x_temp);
  mesh(1,1).BndWestSpline.Xp[5]=mesh(1,1).BndNorthSpline.Xp[0];
  
  //set point types
  mesh(1,1).BndWestSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,1).BndWestSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,1).BndWestSpline.tp[2]=SPLINE2D_POINT_NORMAL;
  mesh(1,1).BndWestSpline.tp[3]=SPLINE2D_POINT_NORMAL;
  mesh(1,1).BndWestSpline.tp[4]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,1).BndWestSpline.tp[5]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh(1,1).BndWestSpline.pathlength();

  // **********************
  // **** BLOCK (1,0) *****
  // **********************
 
  // ** NORTH SPLINE **
  //allocate space
  mesh(1,0).BndNorthSpline.allocate(cs.np-iTrail+iLead);
  
  //set type
  mesh(1,0).BndNorthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  for ( i=iLead; i>=0; --i )
    mesh(1,0).BndNorthSpline.Xp[iLead-i]=cs.Xp[i];
  for ( i=cs.np-2; i>=iTrail; --i )
    mesh(1,0).BndNorthSpline.Xp[iLead+1+(cs.np-2)-i]=cs.Xp[i];

  //set point types
  mesh(1,0).BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,0).BndNorthSpline.tp[mesh(1,0).BndNorthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh(1,0).BndNorthSpline.np-1; ++i )
    mesh(1,0).BndNorthSpline.tp[i]=SPLINE2D_POINT_NORMAL;
  
  //calc pathlength
  mesh(1,0).BndNorthSpline.pathlength();  

  // ** SOUTH SPLINE **
  //allocate space
  mesh(1,0).BndSouthSpline.allocate(70);
  
  //set type
  mesh(1,0).BndSouthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=(zlL+zlU)/2;;
  zTrail=(zrL+zrU)/2;;
  dz=(zTrail-zLead)/69;
  for ( i=0; i<70; ++i ) {
    mesh(1,0).BndSouthSpline.Xp[i]=getY(zLead+i*dz, lowerB)[0];
  }

  //set point types
  mesh(1,0).BndSouthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,0).BndSouthSpline.tp[mesh(1,0).BndSouthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh(1,0).BndSouthSpline.np-1; ++i )
    mesh(1,0).BndSouthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh(1,0).BndSouthSpline.pathlength();

  // ** WEST SPLINE **
  //allocate space
  mesh(1,0).BndWestSpline.allocate(5);
  
  //set type
  mesh(1,0).BndWestSpline.settype(SPLINE2D_CUBIC);

  //set points
  mesh(1,0).BndWestSpline.Xp[0]=mesh(1,0).BndSouthSpline.Xp[0];
  x_temp = getY(zlL, lowerMiddleB)[0];
  mesh(1,0).BndWestSpline.Xp[1]=mesh(1,0).BndSouthSpline.Xp[0]+
                                 HALF*(x_temp-mesh(1,0).BndSouthSpline.Xp[0]);
  mesh(1,0).BndWestSpline.Xp[2]=x_temp;
  mesh(1,0).BndWestSpline.Xp[3]=x_temp+HALF*(cs.Xp[iLead]-x_temp);
  mesh(1,0).BndWestSpline.Xp[4]=cs.Xp[iLead];
  
  //set point types
  mesh(1,0).BndWestSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,0).BndWestSpline.tp[1]=SPLINE2D_POINT_NORMAL;
  mesh(1,0).BndWestSpline.tp[2]=SPLINE2D_POINT_NORMAL;
  mesh(1,0).BndWestSpline.tp[3]=SPLINE2D_POINT_NORMAL;
  mesh(1,0).BndWestSpline.tp[4]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh(1,0).BndWestSpline.pathlength();

  // ** EAST SPLINE **
  //allocate space
  mesh(1,0).BndEastSpline.allocate(6);
  
  //set type
  mesh(1,0).BndEastSpline.settype(SPLINE2D_CUBIC);

  //set points
  mesh(1,0).BndEastSpline.Xp[0]=mesh(1,0).BndSouthSpline.Xp[mesh(1,0).BndSouthSpline.np-1];
  x_temp = getY(zrL, lowerMiddleB)[0];
  mesh(1,0).BndEastSpline.Xp[1]=x_temp + 
     HALF*(mesh(1,0).BndSouthSpline.Xp[mesh(1,0).BndSouthSpline.np-1]-x_temp);
  mesh(1,0).BndEastSpline.Xp[2]=HALF*(cs.Xp[iTrail]+0.6666*(x_temp-cs.Xp[iTrail])+
     x_temp+0.3333*(mesh(1,0).BndSouthSpline.Xp[mesh(1,0).BndSouthSpline.np-1]-x_temp));
  mesh(1,0).BndEastSpline.Xp[3]=0.95*(cs.Xp[iTrail]+HALF*(x_temp-cs.Xp[iTrail]))+
     0.05*(x_temp+0.80*(mesh(1,0).BndSouthSpline.Xp[mesh(1,0).BndSouthSpline.np-1]-x_temp));
  mesh(1,0).BndEastSpline.Xp[4]=cs.Xp[iTrail]+0.3333*(x_temp-cs.Xp[iTrail]);
  mesh(1,0).BndEastSpline.Xp[5]=cs.Xp[iTrail];
   
  //set point types
  mesh(1,0).BndEastSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,0).BndEastSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,0).BndEastSpline.tp[2]=SPLINE2D_POINT_NORMAL;
  mesh(1,0).BndEastSpline.tp[3]=SPLINE2D_POINT_NORMAL;
  mesh(1,0).BndEastSpline.tp[4]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(1,0).BndEastSpline.tp[5]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh(1,0).BndEastSpline.pathlength();

  // **********************
  // **** BLOCK (0,0) *****
  // **********************
  
  // ** NORTH SPLINE **
  //allocate space
  mesh(0,0).BndNorthSpline.allocate(2);
  
  //set type
  mesh(0,0).BndNorthSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh(0,0).BndNorthSpline.Xp[0]=camberLead.Xp[0];
  mesh(0,0).BndNorthSpline.Xp[1]=cs.Xp[iLead];
  
  //set point types
  mesh(0,0).BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(0,0).BndNorthSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh(0,0).BndNorthSpline.pathlength();

  // ** SOUTH SPLINE ** 
  //allocate space
  mesh(0,0).BndSouthSpline.allocate(30);
  
  //set type
  mesh(0,0).BndSouthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=lowerB.Xp[0].x;
  zTrail=mesh(1,0).BndSouthSpline.Xp[0].x;
  dz=(zTrail-zLead)/29;
  for ( i=0; i<30; ++i ) {
    mesh(0,0).BndSouthSpline.Xp[i]=getY(zLead+i*dz, lowerB)[0];
  }
  mesh(0,0).BndSouthSpline.Xp[0]=lowerB.Xp[0];

  //set point types
  mesh(0,0).BndSouthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(0,0).BndSouthSpline.tp[mesh(0,0).BndSouthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh(0,0).BndSouthSpline.np-1; ++i )
    mesh(0,0).BndSouthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh(0,0).BndSouthSpline.pathlength();

  // ** EAST SPLINE **
  mesh(0,0).BndEastSpline = mesh(1,0).BndWestSpline;  

  // ** WEST SPLINE **
  //allocate space
  mesh(0,0).BndWestSpline.allocate(2);
  
  //set type
  mesh(0,0).BndWestSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh(0,0).BndWestSpline.Xp[0]=mesh(0,0).BndSouthSpline.Xp[0];
  mesh(0,0).BndWestSpline.Xp[1]=mesh(0,0).BndNorthSpline.Xp[0];
  
  //set point types
  mesh(0,0).BndWestSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(0,0).BndWestSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh(0,0).BndWestSpline.pathlength();

  // **********************
  // **** BLOCK (2,0) *****
  // **********************
  
  // ** NORTH SPLINE **
  //allocate space
  mesh(2,0).BndNorthSpline.allocate(2);
  
  //set type
  mesh(2,0).BndNorthSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh(2,0).BndNorthSpline.Xp[0]=cs.Xp[iTrail];
  mesh(2,0).BndNorthSpline.Xp[1]=camberTrail.Xp[camberTrail.np-1];
  
  //set point types
  mesh(2,0).BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(2,0).BndNorthSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh(2,0).BndNorthSpline.pathlength();

  // ** SOUTH SPLINE ** 
  //allocate space
  mesh(2,0).BndSouthSpline.allocate(30);
  
  //set type
  mesh(2,0).BndSouthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=(zrL+zrU)/2;;
  zTrail=lowerB.Xp[lowerB.np-1].x;
  dz=(zTrail-zLead)/29;

  for ( i=0; i<30; ++i ) {
    mesh(2,0).BndSouthSpline.Xp[i].x=zLead+i*dz;
    //prevent round-off error
    if(mesh(2,0).BndSouthSpline.Xp[i].x>zTrail) mesh(2,0).BndSouthSpline.Xp[i].x=zTrail;
    mesh(2,0).BndSouthSpline.Xp[i].y=getY(mesh(2,0).BndSouthSpline.Xp[i].x, lowerB)[0].y;
  }
  mesh(2,0).BndSouthSpline.Xp[mesh(2,0).BndSouthSpline.np-1]=lowerB.Xp[lowerB.np-1];

  //set point types
  mesh(2,0).BndSouthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(2,0).BndSouthSpline.tp[mesh(2,0).BndSouthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh(2,0).BndSouthSpline.np-1; ++i )
    mesh(2,0).BndSouthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh(2,0).BndSouthSpline.pathlength();

  // ** WEST SPLINE **
  mesh(2,0).BndWestSpline = mesh(1,0).BndEastSpline; 

  // ** EAST SPLINE **
  //allocate space
  mesh(2,0).BndEastSpline.allocate(2);
  
  //set type
  mesh(2,0).BndEastSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh(2,0).BndEastSpline.Xp[0]=mesh(2,0).BndSouthSpline.Xp[mesh(2,0).BndSouthSpline.np-1];
  mesh(2,0).BndEastSpline.Xp[1]=mesh(2,0).BndNorthSpline.Xp[mesh(2,0).BndNorthSpline.np-1];
  
  //set point types
  mesh(2,0).BndEastSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(2,0).BndEastSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh(2,0).BndEastSpline.pathlength();

  // **********************
  // **** BLOCK (0,1) *****
  // **********************
  
  // ** NORTH SPLINE **
  //allocate space
  mesh(0,1).BndNorthSpline.allocate(30);
  
  //set type
  mesh(0,1).BndNorthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=upperB.Xp[0].x;
  zTrail=(zlL+zlU)/2;
  dz=(zTrail-zLead)/29;
  for ( i=0; i<30; ++i ) {
    mesh(0,1).BndNorthSpline.Xp[i]=getY(zLead+i*dz, upperB)[0];
  }
  mesh(0,1).BndNorthSpline.Xp[0]=upperB.Xp[0];

  //set point types
  mesh(0,1).BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(0,1).BndNorthSpline.tp[mesh(0,1).BndNorthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh(0,1).BndNorthSpline.np-1; ++i )
    mesh(0,1).BndNorthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh(0,1).BndNorthSpline.pathlength();

  // ** SOUTH SPLINE **
  mesh(0,1).BndSouthSpline = mesh(0,0).BndNorthSpline;

  // ** EAST SPLINE **
  mesh(0,1).BndEastSpline = mesh(1,1).BndWestSpline;

  // ** WEST SPLINE **
  //allocate space
  mesh(0,1).BndWestSpline.allocate(2);
  
  //set type
  mesh(0,1).BndWestSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh(0,1).BndWestSpline.Xp[0]=mesh(0,1).BndSouthSpline.Xp[0];
  mesh(0,1).BndWestSpline.Xp[1]=mesh(0,1).BndNorthSpline.Xp[0];
  
  //set point types
  mesh(0,1).BndWestSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(0,1).BndWestSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh(0,1).BndWestSpline.pathlength();

  // **********************
  // **** BLOCK (2,1) *****
  // **********************
  
  // ** NORTH SPLINE **
  //allocate space
  mesh(2,1).BndNorthSpline.allocate(30);
  
  //set type
  mesh(2,1).BndNorthSpline.settype(SPLINE2D_QUINTIC);

  //set points
  zLead=(zrL+zrU)/2;
  zTrail=upperB.Xp[upperB.np-1].x;
  dz=(zTrail-zLead)/29;

  for ( i=0; i<30; ++i ) {
    mesh(2,1).BndNorthSpline.Xp[i].x=zLead+i*dz;
    //prevent round-off error
    if(mesh(2,1).BndNorthSpline.Xp[i].x>zTrail)
      mesh(2,1).BndNorthSpline.Xp[i].x=zTrail;
    mesh(2,1).BndNorthSpline.Xp[i].y=getY(mesh(2,1).BndNorthSpline.Xp[i].x, upperB)[0].y;
  }
  mesh(2,1).BndNorthSpline.Xp[mesh(2,1).BndNorthSpline.np-1]=upperB.Xp[upperB.np-1];

  //set point types
  mesh(2,1).BndNorthSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(2,1).BndNorthSpline.tp[mesh(2,1).BndNorthSpline.np-1]=SPLINE2D_POINT_SHARP_CORNER;
  for ( i=1; i<mesh(2,1).BndNorthSpline.np-1; ++i )
    mesh(2,1).BndNorthSpline.tp[i]=SPLINE2D_POINT_NORMAL;

  //calc pathlength
  mesh(2,1).BndNorthSpline.pathlength();

  // ** SOUTH SPLINE **
  mesh(2,1).BndSouthSpline = mesh(2,0).BndNorthSpline;

  // ** WEST SPLINE **
  mesh(2,1).BndWestSpline = mesh(1,1).BndEastSpline;

  // ** EAST SPLINE **
  //allocate space
  mesh(2,1).BndEastSpline.allocate(2);
  
  //set type
  mesh(2,1).BndEastSpline.settype(SPLINE2D_LINEAR);

  //set points
  mesh(2,1).BndEastSpline.Xp[0]=mesh(2,1).BndSouthSpline.Xp[mesh(2,1).BndSouthSpline.np-1];
  mesh(2,1).BndEastSpline.Xp[1]=mesh(2,1).BndNorthSpline.Xp[mesh(2,1).BndNorthSpline.np-1];
  
  //set point types
  mesh(2,1).BndEastSpline.tp[0]=SPLINE2D_POINT_SHARP_CORNER;
  mesh(2,1).BndEastSpline.tp[1]=SPLINE2D_POINT_SHARP_CORNER;
  
  //calc pathlength
  mesh(2,1).BndEastSpline.pathlength();

  //*******************************************************************************
  //*******************************************************************************
  //*******************************  GENERATE MESH ********************************
  //*******************************************************************************
  //*******************************************************************************


  cout << "\n  Generating Mesh - Block (0,0)...";
  cout.flush();

  genMeshBlock(mesh(0,0), Highest_Order_of_Reconstruction,
	       n, B00_Number_of_Cells_Idir, B00_Number_of_Cells_Jdir, B00_Number_of_Ghost_Cells,
	       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
	       B00_StretchFcnI, B00_BetaI, B00_TauI,
	       B00_StretchFcnJ, B00_BetaJ, B00_TauJ,
	       B00_NbndType, B00_EbndType, B00_SbndType, B00_WbndType,
	       B00_OrthN, B00_OrthE, B00_OrthS, B00_OrthW);

  cout << "\n  Generating Mesh - Block (1,0)...";
  cout.flush();

  genMeshBlock(mesh(1,0), Highest_Order_of_Reconstruction,
	       n, B10_Number_of_Cells_Idir, B10_Number_of_Cells_Jdir, B10_Number_of_Ghost_Cells,
	       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY,
	       B10_StretchFcnI, B10_BetaI, B10_TauI,
	       B10_StretchFcnJ, B10_BetaJ, B10_TauJ,
	       B10_NbndType, B10_EbndType, B10_SbndType, B10_WbndType,
	       B10_OrthN, B10_OrthE, B10_OrthS, B10_OrthW);

  cout << "\n  Generating Mesh - Block (2,0)...";
  cout.flush();

  genMeshBlock(mesh(2,0), Highest_Order_of_Reconstruction,
	       n, B20_Number_of_Cells_Idir, B20_Number_of_Cells_Jdir, B20_Number_of_Ghost_Cells,
	       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
	       B20_StretchFcnI, B20_BetaI, B20_TauI,
	       B20_StretchFcnJ, B20_BetaJ, B20_TauJ,
	       B20_NbndType, B20_EbndType, B20_SbndType, B20_WbndType,
	       B20_OrthN, B20_OrthE, B20_OrthS, B20_OrthW);

  cout << "\n  Generating Mesh - Block (0,1)...";
  cout.flush();

  genMeshBlock(mesh(0,1), Highest_Order_of_Reconstruction,
	       n, B01_Number_of_Cells_Idir, B01_Number_of_Cells_Jdir, B01_Number_of_Ghost_Cells,
               GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
	       B01_StretchFcnI, B01_BetaI, B01_TauI,
	       B01_StretchFcnJ, B01_BetaJ, B01_TauJ,
	       B01_NbndType, B01_EbndType, B01_SbndType, B01_WbndType,
	       B01_OrthN, B01_OrthE, B01_OrthS, B01_OrthW);       

  cout << "\n  Generating Mesh - Block (1,1)...";
  cout.flush();

  genMeshBlock(mesh(1,1), Highest_Order_of_Reconstruction,
	       n, B11_Number_of_Cells_Idir, B11_Number_of_Cells_Jdir, B11_Number_of_Ghost_Cells,
	       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY,
	       B11_StretchFcnI, B11_BetaI, B11_TauI,
	       B11_StretchFcnJ, B11_BetaJ, B11_TauJ,
	       B11_NbndType, B11_EbndType, B11_SbndType, B11_WbndType,
	       B11_OrthN, B11_OrthE, B11_OrthS, B11_OrthW);

  cout << "\n  Generating Mesh - Block (2,1)...";
  cout.flush();

  genMeshBlock(mesh(2,1), Highest_Order_of_Reconstruction,
	       n, B21_Number_of_Cells_Idir, B21_Number_of_Cells_Jdir, B21_Number_of_Ghost_Cells,
	       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
	       B21_StretchFcnI, B21_BetaI, B21_TauI,
	       B21_StretchFcnJ, B21_BetaJ, B21_TauJ,
	       B21_NbndType, B21_EbndType, B21_SbndType, B21_WbndType,
	       B21_OrthN, B21_OrthE, B21_OrthS, B21_OrthW);
}


/**********************************************************************
 * Routine: genMeshBlock                                              *
 **********************************************************************/
void NASARotor37::genMeshBlock(Grid2D_Quad_Block_HO &block,
			       const int & Highest_Order_of_Reconstruction,
			       int n,
			       const int Number_of_Cells_Idir,
			       const int Number_of_Cells_Jdir,
			       const int Number_of_Ghost_Cells,
			       int node_init_procedure,
			       int StretchFcnI, double  BetaI, double  TauI,
			       int StretchFcnJ, double  BetaJ, double  TauJ,
			       int NbndType, int  EbndType, int  SbndType, int WbndType,
			       int OrthN, int  OrthE, int  OrthS, int  OrthW) {

  int i;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South, Bnd_Spline_East, Bnd_Spline_West;

  //set north, south, east and west boundary types for each boundary point
  for ( i=0; i<block.BndNorthSpline.np; ++i )
    block.BndNorthSpline.bc[i]= NbndType;

  for ( i=0; i<block.BndSouthSpline.np; ++i )
    block.BndSouthSpline.bc[i]= SbndType;
  
  for ( i=0; i<block.BndEastSpline.np; ++i )
    block.BndEastSpline.bc[i]= EbndType;
  
  for ( i=0; i<block.BndWestSpline.np; ++i )
    block.BndWestSpline.bc[i]= WbndType;

  //copy boundary splines to temp variables
  Bnd_Spline_North =  block.BndNorthSpline;
  Bnd_Spline_South = block.BndSouthSpline;
  Bnd_Spline_East = block.BndEastSpline;
  Bnd_Spline_West = block.BndWestSpline;
  
  /* Create the 2D quadrilateral grid block. */

  block.Create_Quad_Block_Without_Update(Bnd_Spline_North,
					 Bnd_Spline_South,
					 Bnd_Spline_East,
					 Bnd_Spline_West,
					 Number_of_Cells_Idir,
					 Number_of_Cells_Jdir,
					 Number_of_Ghost_Cells,
					 Highest_Order_of_Reconstruction,
					 node_init_procedure,
					 StretchFcnI, 
					 BetaI, 
					 TauI,
					 StretchFcnJ, 
					 BetaJ, 
					 TauJ,
					 OrthN, 
					 OrthE, 
					 OrthS, 
					 OrthW);

  /* Smooth the 2D quadrilateral grid block. */

  Smooth_Quad_Block(block, n);

}
