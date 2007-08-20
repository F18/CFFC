/* AdvectDiffuse2DQuad_NKS.cc:  Functions for NKS Solver to solve 2D 
                                Advection Diffusion Equation. */

/* Include 2D advection diffusion equation quadrilateral NKS header file. */

#ifndef _ADVECTDIFFUSE2D_QUAD_NKS_INCLUDED 
#include "AdvectDiffuse2DQuad_NKS.h" 
#endif // _ADVECTDIFFUSE2D_QUAD_NKS_INCLUDED 
 
/********************************************************
 * Routine: Jacobian                                    *
 *                                                      *
 * This routine returns the Jacobian matrix for the     *
 * specified solution block.                            *
 *                                                      *
 ********************************************************/
CompRow_Mat_double Jacobian(AdvectDiffuse2D_Quad_Block &SolnBlk, 
			    int JaS, 
			    int JdS, 
			    int JsS, 
			    int overlap) {

  int i,j,k;
  int ncol      = 0;
  int nzcount   = 0;
  int count     = 0;
  int rownum    = 1;
  int xpts = (SolnBlk.ICu-SolnBlk.ICl+1)+SolnBlk.Nghost*2;
  int ypts = (SolnBlk.JCu-SolnBlk.JCl+1)+SolnBlk.Nghost*2;
  
  short int ** model = new short int * [xpts*ypts];
  for (int i = 0; i < xpts*ypts; i++) {
      model[i]= new short int [13];
  }

  /**************************
   * Construct Model Matrix *
   **************************/
  for (i = 0; i < xpts * ypts; i++) {
    for (j = 0; j < 13 ; j++) {
      model[i][j] = 0;
    }
  }

  for ( i = 0; i < xpts*ypts; i++){
    if ( i >= xpts * 2) {
      model[i][0]=1; 
    }
    if ( i >= xpts){
      if ( i % xpts != 0){
	model[i][1]=2;}
      model[i][2]=3;
      if ( i % xpts != xpts - 1){
	model[i][3]=4;}
    }
    if (( i % xpts != 0) & ( i % xpts != 1)){
      model[i][4]=5;}
    if ( i % xpts != 0){
      model[i][5]=6;}
    model[i][6]=7;
    if ( i % xpts != xpts - 1){
      model[i][7]=8;}
    if (( i % xpts != xpts - 1) && ( i % xpts != xpts - 2)){
      model[i][8]=9;}
    if( i < xpts * ypts - xpts){
      if ( i % xpts != 0){
	model[i][9]=10;}
      model[i][10]=11;
      if ( i % xpts != xpts - 1){
	model[i][11]=12;}
    }
    if ( i < xpts * ypts - xpts * 2){
      model[i][12]=13;}
  }

  for ( i = 0; i < xpts * ypts; i++){
    for ( j = 0; j < 13; j++){ 
      if (model[i][j] != 0){
	nzcount++;}
    }
  }

  /*******************************************************************
   * Initialize and assign vectors for "CompRow_Mat_Double" function *
   *******************************************************************/
  VECTOR_double     val(nzcount,0.0);  
  VECTOR_int     colind(nzcount,0);
  VECTOR_int   wholerow(nzcount,0);
  VECTOR_int     rowind(xpts*ypts+1,0); 
  
  for ( i = 0; i < xpts * ypts; i++){
    for ( j = 0; j < 13; j++){
      if (model[i][j] != 0){
	if ( j == 0){
	  ncol= i-xpts*2;}
	else if ( j == 1){
	  ncol=i-xpts-1;}
	else if ( j == 2){
	  ncol=i-xpts;}
	else if ( j == 3){
	  ncol=i+1-xpts;}
	else if ( j == 4){
	  ncol=i-2;}
	else if ( j == 5){
	  ncol=i-1;}
	else if ( j == 6){
	  ncol=i;}
	else if ( j == 7){
	  ncol=i+1;}
	else if ( j == 8){
	  ncol=i+2;}
	else if ( j == 9){
	  ncol=i+xpts-1;}
	else if ( j == 10){
	  ncol=i+xpts;}
	else if ( j == 11){
	  ncol=i+1+xpts;}
	else if ( j == 12){
	  ncol=i+xpts*2;}
	else {ncol=0;
	}
	wholerow[count] = i;
	colind[count] = ncol;
	val[count] =   Ja(SolnBlk, ncol, i, JaS, overlap) 
	             + Jd(SolnBlk, ncol, i, JdS, JaS, overlap) 
	             + Js(SolnBlk, ncol, i, JsS, JdS, JaS, overlap);
	count++;
      }
    }
  }
  
  rowind[0]=0;
  for (i=0; i<=count-2; i++){
    if (wholerow[i] != wholerow[i+1]) {
      rowind[rownum] = i+1;
      rownum++;
    }
  }
  rowind[rownum]=count;

  CompRow_Mat_double JM(xpts*ypts,xpts*ypts,nzcount,val,rowind,colind);
  
  for (int i = 0; i < xpts*ypts; i++) {
      delete [] model[i];
  }
  delete [] model;

  return(JM);

} /* End of Jacobian Function */ 

/********************************************************
 * Routine: Ja                                          *
 *                                                      *
 * This routine returns the Jacobian matrix cefficient  *
 * based on the advection equation for the specified    *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
double Ja(AdvectDiffuse2D_Quad_Block &SolnBlk, 
	  int ncol, 
	  int nrow, 
	  int JaS, 
	  int overlap) {
  
  double Jadv; 
  
  int xpts = (SolnBlk.ICu-SolnBlk.ICl+1)+SolnBlk.Nghost*2;
  int ypts = (SolnBlk.JCu-SolnBlk.JCl+1)+SolnBlk.Nghost*2; 
  int a_xpts = SolnBlk.ICu-SolnBlk.ICl+1;
  int a_ypts = SolnBlk.JCu-SolnBlk.JCl+1;
  int row = int (double (nrow / xpts));
  int col = nrow % xpts;
  int i = col;
  int j = row;
  
  int Fswitch = ncol - nrow; 

  if (JaS == 1) { 
    
    if ((col >= a_xpts + SolnBlk.Nghost + overlap) || 
	(col < SolnBlk.Nghost - overlap) || 
	(row >= a_ypts + SolnBlk.Nghost + overlap) || 
	(row < SolnBlk.Nghost - overlap)) 
      {
	
	// GHOST CELLS
	if (Fswitch == 0) {
	  Jadv = 1;
	} else {
	  Jadv = 0;
	}
	
      } else {
	
	// INTERIOR CELLS
	double Cn, Cs, Ce, Cw;
	double Ns, Ss, Es, Ws;
	double Na, Sa, Ea, Wa;
	
	Vector2D VN, VS, VE, VW;
	Vector2D NN, NS, NORTHEAST, NORTHWEST;
	
	MV_Vector_double Vn(2,0.0);
	MV_Vector_double Vs(2,0.0);      
	MV_Vector_double Ve(2,0.0);
	MV_Vector_double Vw(2,0.0);
	
	MV_Vector_double nn(2,0.0);
	MV_Vector_double ns(2,0.0);      
	MV_Vector_double ne(2,0.0);
	MV_Vector_double nw(2,0.0);
	
	/**************************************/
	/** CALCULATE ADVECTION TERMS *********/
	/**************************************/
	
	/* Advection velocities *******/
	VN = HALF*(  SolnBlk.U[i][j].V + SolnBlk.U[i][j+1].V);
	VS = HALF*(SolnBlk.U[i][j-1].V +   SolnBlk.U[i][j].V);
	VE = HALF*(SolnBlk.U[i+1][j].V +   SolnBlk.U[i][j].V);
	VW = HALF*(  SolnBlk.U[i][j].V + SolnBlk.U[i-1][j].V);
	
	Vn[0]=VN.x;Vn[1]=VN.y;
	Vs[0]=VS.x;Vs[1]=VS.y;
	Ve[0]=VE.x;Ve[1]=VE.y;
	Vw[0]=VW.x;Vw[1]=VW.y;
	/******************************/
	
	/* Normal vectors *************/
	NN = SolnBlk.Grid.nfaceN(i,j);
	NS = SolnBlk.Grid.nfaceS(i,j); 
	NORTHEAST = SolnBlk.Grid.nfaceE(i,j); 
	NORTHWEST = SolnBlk.Grid.nfaceW(i,j); 
	
	nn[0]=NN.x;nn[1]=NN.y;
	ns[0]=NS.x;ns[1]=NS.y;
	ne[0]=NORTHEAST.x;ne[1]=NORTHEAST.y;
	nw[0]=NORTHWEST.x;nw[1]=NORTHWEST.y;
	/******************************/
	
	Cn = dot(Vn, nn);
	Cs = dot(Vs, ns);
	Ce = dot(Ve, ne);
	Cw = dot(Vw, nw);
	
	Ns = Cn - abs(Cn);
	Ss = Cs - abs(Cs);
	Es = Ce - abs(Ce);
	Ws = Cw - abs(Cw);
	
	Na = abs(Cn) + Cn;
	Sa = abs(Cs) + Cs;
	Ea = abs(Ce) + Ce;
	Wa = abs(Cw) + Cw;
	
	
	if (Fswitch == - xpts) {
	  // Construct Coefficient X
	  Jadv = -HALF*SolnBlk.Grid.lfaceS(i,j)*Ss/SolnBlk.Grid.area(i,j);
	} else if (Fswitch == -1) {
	  // Construct Coefficient B
	  Jadv = -HALF*SolnBlk.Grid.lfaceW(i,j)*Ws/SolnBlk.Grid.area(i,j);
	} else if (Fswitch == 0) {
	  // Construct Coefficient C
	  Jadv = -HALF*(SolnBlk.Grid.lfaceN(i,j)*Na+SolnBlk.Grid.lfaceS(i,j)*Sa+
                 SolnBlk.Grid.lfaceE(i,j)*Ea+SolnBlk.Grid.lfaceW(i,j)*Wa)/SolnBlk.Grid.area(i,j);
	} else if (Fswitch == 1) {
	  // Construct Coefficient D
	  Jadv = -HALF*SolnBlk.Grid.lfaceE(i,j)*Es/SolnBlk.Grid.area(i,j);
	} else if (Fswitch == xpts) {
	  // Construct Coefficient Y 
	  Jadv = -HALF*SolnBlk.Grid.lfaceN(i,j)*Ns/SolnBlk.Grid.area(i,j);
	} else { // Others
	  Jadv = 0;
	}
      }
  } else {
    Jadv = 0;
  } 
  
  return (Jadv);
  
}/* End of Ja Function */ 
 
/********************************************************
 * Routine: Jd                                          *
 *                                                      *
 * This routine returns Jacobian matrix coefficient     *
 * based on the diffusion equation for the specified    *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
double Jd(AdvectDiffuse2D_Quad_Block &SolnBlk, 
	  int ncol, 
	  int nrow, 
	  int JdS, 
	  int JaS, 
	  int overlap) {

  double Jdiff;

  int xpts = (SolnBlk.ICu-SolnBlk.ICl+1)+SolnBlk.Nghost*2;
  int ypts = (SolnBlk.JCu-SolnBlk.JCl+1)+SolnBlk.Nghost*2; 
  int a_xpts = SolnBlk.ICu-SolnBlk.ICl+1;
  int a_ypts = SolnBlk.JCu-SolnBlk.JCl+1;  
  int row = int (double (nrow / xpts));
  int col = nrow % xpts;
  int i = col;
  int j = row;
  int Fswitch = ncol - nrow; 

  if (JdS == 1) { 
    
    if ((col >= a_xpts + SolnBlk.Nghost + overlap) || 
	(col < SolnBlk.Nghost - overlap) || 
	(row >= a_ypts + SolnBlk.Nghost + overlap) || 
	(row < SolnBlk.Nghost - overlap)) 
      {
	
	// GHOST CELLS
	if (JaS == 1) {
	  Jdiff = 0;
	} else if (Fswitch == 0) {
	  Jdiff = 1;
	} else {
	  Jdiff = 0;	
	}
      } else {
	
	// INTERIOR CELLS
	
	double Cn, Cs, Ce, Cw;
	double Ns, Ss, Es, Ws;
	double Na, Sa, Ea, Wa;
	double kappaN, kappaS, kappaE, kappaW;
	
	Vector2D VN, VS, VE, VW;
	Vector2D NN, NS, NORTHEAST, NORTHWEST;

	MV_Vector_double Vn(2,0.0);
	MV_Vector_double Vs(2,0.0);      
	MV_Vector_double Ve(2,0.0);
	MV_Vector_double Vw(2,0.0);

	// Normal vectors for block (i,j) 
	MV_Vector_double nn_ij(2,0.0);
	MV_Vector_double ns_ij(2,0.0);      
	MV_Vector_double ne_ij(2,0.0);
	MV_Vector_double nw_ij(2,0.0);
	// Normal vectors for block (i-1,j) 
	MV_Vector_double nn_mij(2,0.0);
	MV_Vector_double ns_mij(2,0.0);      
	MV_Vector_double ne_mij(2,0.0);
	MV_Vector_double nw_mij(2,0.0);
	// Normal vectors for block (i+1,j) 
	MV_Vector_double nn_pij(2,0.0);
	MV_Vector_double ns_pij(2,0.0);      
	MV_Vector_double ne_pij(2,0.0);
	MV_Vector_double nw_pij(2,0.0);
	// Normal vectors for block (i,j-1) 
	MV_Vector_double nn_ijm(2,0.0);
	MV_Vector_double ns_ijm(2,0.0);      
	MV_Vector_double ne_ijm(2,0.0);
	MV_Vector_double nw_ijm(2,0.0);
	// Normal vectors for block (i,j+1) 
	MV_Vector_double nn_ijp(2,0.0);
	MV_Vector_double ns_ijp(2,0.0);      
	MV_Vector_double ne_ijp(2,0.0);
	MV_Vector_double nw_ijp(2,0.0);

	/**************************************/
	/** CALCULATE DIFFUSION TERMS *********/
	/**************************************/
	
	/* Diffusion coefficients ******/
	kappaN = HALF*(SolnBlk.U[i][j+1].k +   SolnBlk.U[i][j].k);
	kappaS = HALF*(  SolnBlk.U[i][j].k + SolnBlk.U[i][j-1].k);
	kappaE = HALF*(SolnBlk.U[i+1][j].k +   SolnBlk.U[i][j].k);
	kappaW = HALF*(  SolnBlk.U[i][j].k + SolnBlk.U[i-1][j].k);
	/******************************/
	
	/* Normal vectors *************/
	// Block (i,j)
	NN = SolnBlk.Grid.nfaceN(i,j);
	NS = SolnBlk.Grid.nfaceS(i,j); 
	NORTHEAST = SolnBlk.Grid.nfaceE(i,j); 
	NORTHWEST = SolnBlk.Grid.nfaceW(i,j); 
	
	nn_ij[0]=NN.x;nn_ij[1]=NN.y;
	ns_ij[0]=NS.x;ns_ij[1]=NS.y;
	ne_ij[0]=NORTHEAST.x;ne_ij[1]=NORTHEAST.y;
	nw_ij[0]=NORTHWEST.x;nw_ij[1]=NORTHWEST.y;

	// Block (i-1,j)
	NN = SolnBlk.Grid.nfaceN(i-1,j);
	NS = SolnBlk.Grid.nfaceS(i-1,j); 
	NORTHEAST = SolnBlk.Grid.nfaceE(i-1,j); 
	NORTHWEST = SolnBlk.Grid.nfaceW(i-1,j); 

	nn_mij[0]=NN.x;nn_mij[1]=NN.y;
	ns_mij[0]=NS.x;ns_mij[1]=NS.y;
	ne_mij[0]=NORTHEAST.x;ne_mij[1]=NORTHEAST.y;
	nw_mij[0]=NORTHWEST.x;nw_mij[1]=NORTHWEST.y;

	// Block (i+1,j)
	NN = SolnBlk.Grid.nfaceN(i+1,j);
	NS = SolnBlk.Grid.nfaceS(i+1,j); 
	NORTHEAST = SolnBlk.Grid.nfaceE(i+1,j); 
	NORTHWEST = SolnBlk.Grid.nfaceW(i+1,j); 

	nn_pij[0]=NN.x;nn_pij[1]=NN.y;
	ns_pij[0]=NS.x;ns_pij[1]=NS.y;
	ne_pij[0]=NORTHEAST.x;ne_pij[1]=NORTHEAST.y;
	nw_pij[0]=NORTHWEST.x;nw_pij[1]=NORTHWEST.y;

	// Block (i,j-1)
	NN = SolnBlk.Grid.nfaceN(i,j-1);
	NS = SolnBlk.Grid.nfaceS(i,j-1); 
	NORTHEAST = SolnBlk.Grid.nfaceE(i,j-1); 
	NORTHWEST = SolnBlk.Grid.nfaceW(i,j-1); 
	
	nn_ijm[0]=NN.x;nn_ijm[1]=NN.y;
	ns_ijm[0]=NS.x;ns_ijm[1]=NS.y;
	ne_ijm[0]=NORTHEAST.x;ne_ijm[1]=NORTHEAST.y;
	nw_ijm[0]=NORTHWEST.x;nw_ijm[1]=NORTHWEST.y;
	
	// Block (i,j+1)
	NN = SolnBlk.Grid.nfaceN(i,j+1);
	NS = SolnBlk.Grid.nfaceS(i,j+1); 
	NORTHEAST = SolnBlk.Grid.nfaceE(i,j+1); 
	NORTHWEST = SolnBlk.Grid.nfaceW(i,j+1); 
	
	nn_ijp[0]=NN.x;nn_ijp[1]=NN.y;
	ns_ijp[0]=NS.x;ns_ijp[1]=NS.y;
	ne_ijp[0]=NORTHEAST.x;ne_ijp[1]=NORTHEAST.y;
	nw_ijp[0]=NORTHWEST.x;nw_ijp[1]=NORTHWEST.y;
	/******************************/
	
	if (Fswitch == -(xpts*2)) {
	  // Construct Coefficient W
	  Jdiff = QUARTER*(dot(ns_ij,ns_ijm)*SolnBlk.Grid.lfaceS(i,j)*SolnBlk.Grid.lfaceS(i,j-1)*kappaS)/
                  (SolnBlk.Grid.area(i,j)*SolnBlk.Grid.area(i,j-1));
	} else if (Fswitch == -(xpts+1)) {
	  // Construct Coefficient P
	  Jdiff = (QUARTER/SolnBlk.Grid.area(i,j))*
	    ((dot(ns_ij,nw_ijm)*SolnBlk.Grid.lfaceS(i,j)*SolnBlk.Grid.lfaceW(i,j-1)*
	      kappaS)/SolnBlk.Grid.area(i,j-1)+
	     (dot(nw_ij,ns_mij)*SolnBlk.Grid.lfaceW(i,j)*SolnBlk.Grid.lfaceS(i-1,j)*
	      kappaW)/SolnBlk.Grid.area(i-1,j));
	} else if (Fswitch == -xpts) {
	  // Construct Coefficient X
	  Jdiff = (QUARTER*SolnBlk.Grid.lfaceS(i,j)/SolnBlk.Grid.area(i,j))*
	    ((kappaS/SolnBlk.Grid.area(i,j-1))*
	     (dot(ns_ij,nn_ijm)*SolnBlk.Grid.lfaceN(i,j-1)+
	      dot(ns_ij,ns_ijm)*SolnBlk.Grid.lfaceS(i,j-1)+
	      dot(ns_ij,ne_ijm)*SolnBlk.Grid.lfaceE(i,j-1)+
	      dot(ns_ij,nw_ijm)*SolnBlk.Grid.lfaceW(i,j-1))+
	     (1.0/SolnBlk.Grid.area(i,j))*
	     (dot(ns_ij,nn_ij)*SolnBlk.Grid.lfaceN(i,j)*kappaN+
	      SolnBlk.Grid.lfaceS(i,j)*kappaS+
	      dot(ns_ij,ne_ij)*SolnBlk.Grid.lfaceE(i,j)*kappaE+
	      dot(ns_ij,nw_ij)*SolnBlk.Grid.lfaceW(i,j)*kappaW));	
	} else if (Fswitch == -(xpts-1)) {
	  // Construct Coefficient Q
	  Jdiff = (QUARTER/SolnBlk.Grid.area(i,j))*
	    ((dot(ns_ij,ne_ijm)*SolnBlk.Grid.lfaceS(i,j)*SolnBlk.Grid.lfaceE(i,j-1)*
	      kappaS)/SolnBlk.Grid.area(i,j-1)+
	     (dot(ne_ij,ns_pij)*SolnBlk.Grid.lfaceE(i,j)*SolnBlk.Grid.lfaceS(i+1,j)*
	      kappaE)/SolnBlk.Grid.area(i+1,j));
	} else if (Fswitch == -2) {
	  // Construct Coefficient A
	  Jdiff = QUARTER*(dot(nw_ij,nw_mij)*SolnBlk.Grid.lfaceW(i,j)*
			   SolnBlk.Grid.lfaceW(i-1,j)*kappaW)/
	    (SolnBlk.Grid.area(i,j)*SolnBlk.Grid.area(i-1,j));
	} else if (Fswitch == -1) {
	  // Construct Coefficient B
	  Jdiff = (QUARTER*SolnBlk.Grid.lfaceW(i,j)/SolnBlk.Grid.area(i,j))*
	    ((kappaW/SolnBlk.Grid.area(i-1,j))*
	     (dot(nw_ij,nn_mij)*SolnBlk.Grid.lfaceN(i-1,j)+
	      dot(nw_ij,ns_mij)*SolnBlk.Grid.lfaceS(i-1,j)+
	      dot(nw_ij,ne_mij)*SolnBlk.Grid.lfaceE(i-1,j)+
	      dot(nw_ij,nw_mij)*SolnBlk.Grid.lfaceW(i-1,j))+
	     (1.0/SolnBlk.Grid.area(i,j))*
	     (dot(nw_ij,nn_ij)*SolnBlk.Grid.lfaceN(i,j)*kappaN+
	      dot(nw_ij,ns_ij)*SolnBlk.Grid.lfaceS(i,j)*kappaS+
	      dot(nw_ij,ne_ij)*SolnBlk.Grid.lfaceE(i,j)*kappaE+
	      SolnBlk.Grid.lfaceW(i,j)*kappaW));
	} else if (Fswitch == 0) {	   				    
	  // Construct Coefficient C
	  Jdiff = (QUARTER/SolnBlk.Grid.area(i,j))*
	    ((dot(nn_ij,ns_ijp)*SolnBlk.Grid.lfaceN(i,j)*SolnBlk.Grid.lfaceS(i,j+1)*
	      kappaN)/SolnBlk.Grid.area(i,j+1)+
	     (dot(ns_ij,nn_ijm)*SolnBlk.Grid.lfaceS(i,j)*SolnBlk.Grid.lfaceN(i,j-1)*
	      kappaS)/SolnBlk.Grid.area(i,j-1)+
	     (dot(ne_ij,nw_pij)*SolnBlk.Grid.lfaceE(i,j)*SolnBlk.Grid.lfaceW(i+1,j)*
	      kappaE)/SolnBlk.Grid.area(i+1,j)+
	     (dot(nw_ij,ne_mij)*SolnBlk.Grid.lfaceW(i,j)*SolnBlk.Grid.lfaceE(i-1,j)*
	      kappaW)/SolnBlk.Grid.area(i-1,j)+
	   
	     (SolnBlk.Grid.lfaceN(i,j)*kappaN/SolnBlk.Grid.area(i,j))*
	     (SolnBlk.Grid.lfaceN(i,j)+
	      dot(nn_ij,ns_ij)*SolnBlk.Grid.lfaceS(i,j)+
	      dot(nn_ij,ne_ij)*SolnBlk.Grid.lfaceE(i,j)+
	      dot(nn_ij,nw_ij)*SolnBlk.Grid.lfaceW(i,j))+
	   
	     (SolnBlk.Grid.lfaceS(i,j)*kappaS/SolnBlk.Grid.area(i,j))*
	     (dot(ns_ij,nn_ij)*SolnBlk.Grid.lfaceN(i,j)+
	      SolnBlk.Grid.lfaceS(i,j)+
	      dot(ns_ij,ne_ij)*SolnBlk.Grid.lfaceE(i,j)+
	      dot(ns_ij,nw_ij)*SolnBlk.Grid.lfaceW(i,j))+
	   
	     (SolnBlk.Grid.lfaceE(i,j)*kappaE/SolnBlk.Grid.area(i,j))*
	     (dot(ne_ij,nn_ij)*SolnBlk.Grid.lfaceN(i,j)+
	      dot(ne_ij,ns_ij)*SolnBlk.Grid.lfaceS(i,j)+
	      SolnBlk.Grid.lfaceE(i,j)+
	      dot(ne_ij,nw_ij)*SolnBlk.Grid.lfaceW(i,j))+
	   
	     (SolnBlk.Grid.lfaceW(i,j)*kappaW/SolnBlk.Grid.area(i,j))*
	     (dot(nw_ij,nn_ij)*SolnBlk.Grid.lfaceN(i,j)+
	      dot(nw_ij,ns_ij)*SolnBlk.Grid.lfaceS(i,j)+
	      dot(nw_ij,ne_ij)*SolnBlk.Grid.lfaceE(i,j)+
	      SolnBlk.Grid.lfaceW(i,j))
	     );
	} else if (Fswitch == 1) {	
	  // Construct Coefficient D
	  Jdiff = (QUARTER*SolnBlk.Grid.lfaceE(i,j)/SolnBlk.Grid.area(i,j))*
	    ((kappaE/SolnBlk.Grid.area(i+1,j))*
	     (dot(ne_ij,nn_pij)*SolnBlk.Grid.lfaceN(i+1,j)+
	      dot(ne_ij,ns_pij)*SolnBlk.Grid.lfaceS(i+1,j)+
	      dot(ne_ij,ne_pij)*SolnBlk.Grid.lfaceE(i+1,j)+
	      dot(ne_ij,nw_pij)*SolnBlk.Grid.lfaceW(i+1,j))+
	     (1.0/SolnBlk.Grid.area(i,j))*
	     (dot(ne_ij,nn_ij)*SolnBlk.Grid.lfaceN(i,j)*kappaN+
	      dot(ne_ij,ns_ij)*SolnBlk.Grid.lfaceS(i,j)*kappaS+
	      SolnBlk.Grid.lfaceE(i,j)*kappaE+
	      dot(ne_ij,nw_ij)*SolnBlk.Grid.lfaceW(i,j)*kappaW));
	} else if (Fswitch == 2) {	
	  // Construct Coefficient E
	  Jdiff = QUARTER*(dot(ne_pij,ne_ij)*SolnBlk.Grid.lfaceE(i+1,j)*
			   SolnBlk.Grid.lfaceE(i,j)*kappaE)/
	    (SolnBlk.Grid.area(i+1,j)*SolnBlk.Grid.area(i,j));
	} else if (Fswitch == (xpts-1)) {	
	  // Construct Coefficient R
	  Jdiff = (QUARTER/SolnBlk.Grid.area(i,j))*
	    ((dot(nn_ij,nw_ijp)*SolnBlk.Grid.lfaceN(i,j)*SolnBlk.Grid.lfaceW(i,j+1)*
	      kappaN)/SolnBlk.Grid.area(i,j+1)+
	     (dot(nw_ij,nn_mij)*SolnBlk.Grid.lfaceW(i,j)*SolnBlk.Grid.lfaceN(i-1,j)*
	      kappaE)/SolnBlk.Grid.area(i-1,j));
	} else if (Fswitch == xpts) {	
	  // Construct Coefficient Y
	  Jdiff = (QUARTER*SolnBlk.Grid.lfaceN(i,j)/SolnBlk.Grid.area(i,j))*
	    ((kappaN/SolnBlk.Grid.area(i,j+1))*
	     (dot(nn_ij,nn_ijp)*SolnBlk.Grid.lfaceN(i,j+1)+
	      dot(nn_ij,ns_ijp)*SolnBlk.Grid.lfaceS(i,j+1)+
	      dot(nn_ij,ne_ijp)*SolnBlk.Grid.lfaceE(i,j+1)+
	      dot(nn_ij,nw_ijp)*SolnBlk.Grid.lfaceW(i,j+1))+
	     (1.0/SolnBlk.Grid.area(i,j))*
	     (                 SolnBlk.Grid.lfaceN(i,j)*kappaN+
			       dot(nn_ij,ns_ij)*SolnBlk.Grid.lfaceS(i,j)*kappaS+
			       dot(nn_ij,ne_ij)*SolnBlk.Grid.lfaceE(i,j)*kappaE+
			       dot(nn_ij,nw_ij)*SolnBlk.Grid.lfaceW(i,j)*kappaW));
	} else if (Fswitch == (xpts+1)) {
	  // Construct Coefficient S
	  Jdiff = (QUARTER/SolnBlk.Grid.area(i,j))*
	    ((dot(nn_ij,ne_ijp)*SolnBlk.Grid.lfaceN(i,j)*SolnBlk.Grid.lfaceE(i,j+1)*
	      kappaN)/SolnBlk.Grid.area(i,j+1)+
	     (dot(ne_ij,nn_pij)*SolnBlk.Grid.lfaceE(i,j)*SolnBlk.Grid.lfaceN(i+1,j)*
	      kappaE)/SolnBlk.Grid.area(i+1,j));
	} else if (Fswitch == (xpts*2)) {
	  // Construct Coefficient Z
	  Jdiff = QUARTER*(dot(nn_ij,nn_ijp)*SolnBlk.Grid.lfaceN(i,j)*
			   SolnBlk.Grid.lfaceN(i,j+1)*kappaN)/
	    (SolnBlk.Grid.area(i,j)*SolnBlk.Grid.area(i,j+1));
	} else {
	  Jdiff = 5e20;
	  cout << "************* DIFFUSION : PROBLEM !!!!!! **************" << endl;
	}
      }

  } else {
    Jdiff = 0;
  }
  
  return (Jdiff); 
  
} /* End of Jd Function */ 

/********************************************************
 * Routine: Js                                          *
 *                                                      *
 * This routine returns Jacobian matrix coefficient     *
 * generated from the source term for the specified     *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
double Js(AdvectDiffuse2D_Quad_Block &SolnBlk, 
	  int ncol, 
	  int nrow, 
	  int JsS, 
	  int JdS, 
	  int JaS, 
	  int overlap) {
  
  double Jsrc; 
  
  int xpts = (SolnBlk.ICu-SolnBlk.ICl+1)+SolnBlk.Nghost*2;
  int ypts = (SolnBlk.JCu-SolnBlk.JCl+1)+SolnBlk.Nghost*2; 
  int a_xpts = SolnBlk.ICu-SolnBlk.ICl+1;
  int a_ypts = SolnBlk.JCu-SolnBlk.JCl+1;  
  int row = int (double (nrow / xpts));
  int col = nrow % xpts;
  int i = col;
  int j = row;
  int Fswitch = ncol - nrow; 
  
  if (JsS == 1) { 
    
    if ((col >= a_xpts + SolnBlk.Nghost + overlap ) || 
	(col < SolnBlk.Nghost - overlap) || 
	(row >= a_ypts + SolnBlk.Nghost + overlap) || 
	(row < SolnBlk.Nghost - overlap)) 
      {
	
	// GHOST CELLS
	Jsrc =0;
	
      } else {
	
	// INTERIOR CELLS    
	
	// Construct Coefficient C
	if (Fswitch == 0) {
	  Jsrc = - 1 / SolnBlk.U[i][j].T;
	} else {
	  Jsrc = 0;
	}
      }
    
  } else {
    Jsrc = 0;
  }
  
  return (Jsrc);   
  
}/* End of Js Function */ 

/********************************************************
 * Routine: Newton_Krylov_Schwarz_Solver                *
 *                                                      *
 * This routine updates the specified solution block    *
 * using Newton-Krylov-Schwarz method.                  *
 *                                                      *
 ********************************************************/ 
int Newton_Krylov_Schwarz_Solver(ostream &Progress_File,
                                 int &Number_of_Startup_Interations,
				 AdvectDiffuse2D_Quad_Block *Soln_ptr,
  			         AdaptiveBlock2D_List &Soln_Block_List,
  			         AdvectDiffuse2D_Input_Parameters &Input_Parameters) { 
  
  CPUTime Time;

  // Overall Convergence tolerance
  double tol      = 1e-10;              
  // GMRES Convergence tolerance
  double gmrestol = Input_Parameters.GMRES_Toler;

  double L2norm_current     = 0.0, 
         L1norm_current     = 0.0, 
         Max_norm_current   = 0.0; 
  double L2norm_first       = 0.0, 
         L1norm_first       = 0.0,
         Max_norm_first     = 0.0;
  double L2norm_current_n   = 0.0, 
         L1norm_current_n   = 0.0,
         Max_norm_current_n = 0.0;
  double previous_L2norm    = 0.0;
  double norm_ratio         = 0.0;
  double dTime; 

  int count           = 0;
  int stop            = 1;
  int max_newton_iter = Input_Parameters.Maximum_Number_of_NKS_Iterations;
  int max_gmres_iter  = Input_Parameters.Maximum_Number_of_GMRES_Iterations;
  int restart         = Input_Parameters.GMRES_Restart;
  int overlap         = Input_Parameters.GMRES_Overlap;
  int NBLK            = Soln_Block_List.Nblk;

  int i, j, Bcount, error_flag, xpts, ypts, xcount, ycount;
  int Adv_Switch, Diff_Switch, Src_Switch;
  int P_Switch;

  int    limiter_check    = ON;
  int    Freeze_Limiter = Input_Parameters.Freeze_Limiter;
  double Freeze_Limiter_Residual_Level = Input_Parameters.Freeze_Limiter_Residual_Level;
  int    finite_time_step = Input_Parameters.Finite_Time_Step;

  /* Switch for advection terms */
  Adv_Switch = (Input_Parameters.a == 0 && Input_Parameters.b == 0)?0:1;

  /* Switch for diffusion terms */
  Diff_Switch = (Input_Parameters.Kappa == 0)?0:1;

  /* Switch for source terms */
  Src_Switch = (Input_Parameters.Tau > 1e7)?0:1;

  /* Switch for preconditioner */
  P_Switch = Input_Parameters.GMRES_P_Switch;

  if (CFFC_Primary_MPI_Processor()) {
    for (int star=0;star<75;star++){cout <<"*";}
    cout << "\nAdvection  = " << Adv_Switch << "; Diffusion = " << Diff_Switch 
         << "; Source = " << Src_Switch << endl;
    cout << "(1 = ON; 0 = OFF)" << endl;  
    
    for (int star=0;star<75;star++){cout <<"*";}
    cout << "\n********     SparseLib++ SparseLib++ SparseLib++ SparseLib++     **********" << endl;
    cout << "********                     Newton-Krylov-Schwarz               **********" << endl;   
    for (int star=0;star<75;star++){cout <<"*";}
    cout << "\nLimiter Type          ====> ";
    if (Input_Parameters.i_Limiter == LIMITER_ZERO)   { 
      cout << "ZERO" << endl;
    } else if (Input_Parameters.i_Limiter == LIMITER_BARTH_JESPERSEN) { 
      cout <<"BARTH_JESPERSEN" << endl;
    } else if (Input_Parameters.i_Limiter == LIMITER_VENKATAKRISHNAN) { 
      cout << "VENKATAKRISHNAN" << endl;
    } else { 
      cout << "??" << endl;
    } /* endif */
    if (P_Switch == 1) {         // ILU(0) 
      cout << "Local Preconditioner  ====> ILU(0)" << endl;
    } else if (P_Switch == 2){   // Diagonal
      cout << "Local Preconditioner  ====> Diagonal" << endl; 
    } else {                     // Identity
      cout << "Local Preconditioner  ====> Identity" << endl; 
    } /* endif */  
    cout <<     "Overall Tolerance     ====> " << tol << endl;
    cout <<     "GMRES Tolerance       ====> " << gmrestol << endl;
    cout <<     "GMRES Restart         ====> " << restart << endl;
    cout <<     "Level of Overlap      ====> " << overlap << endl;
    cout <<     "Number of Processors  ====> " << Input_Parameters.Number_of_Processors << endl;
    cout <<     "Block Size            ====> " << Input_Parameters.Number_of_Blocks_Idir 
         << " x " << Input_Parameters.Number_of_Blocks_Jdir << endl;
    cout <<     "Overall Cell Size     ====> " << Input_Parameters.Number_of_Cells_Idir 
         << " x " << Input_Parameters.Number_of_Cells_Jdir << endl;
    if (finite_time_step == ON) {
      cout << "Finite Time Step      ====> ON" << endl;
      cout << "Initial_CFL           ====> " 
	   << Input_Parameters.Finite_Time_Step_Initial_CFL << endl;
    } else {
      cout << "Finite Time Step      ====> OFF" << endl; 
    } /* endif */ 
    if (Freeze_Limiter == ON) {
      cout << "Freeze_Limiter        ====> ON" << endl;
      cout << "Residual_Level        ====> "
	   << Input_Parameters.Freeze_Limiter_Residual_Level << endl;
    } else {
      cout << "Freeze_Limiter      ====> OFF" << endl; 
    } /* endif */
    for (int star=0;star<75;star++){cout <<"*";}
    cout << "  " << endl;

  } /* endif (CFFC_Primary_MPI_Processor())  */ 
  
  /*******************************************/
  /* BEGIN NEWTON-KRYLOV-SCHWARZ CALCULATION */
  /*******************************************/
  
  /* Create GMRES vector for each block. */
  GMRES *G = new GMRES[NBLK];
  for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      xpts = (Soln_ptr[Bcount].ICu-Soln_ptr[Bcount].ICl+1) +
	      Soln_ptr[Bcount].Nghost*2;
      ypts = (Soln_ptr[Bcount].JCu-Soln_ptr[Bcount].JCl+1) +
	      Soln_ptr[Bcount].Nghost*2;
      G[Bcount].allocate(restart, xpts, ypts);
      G[Bcount].xpts           = xpts;
      G[Bcount].ypts           = ypts;
      G[Bcount].Nghost         = Soln_ptr[Bcount].Nghost;
      G[Bcount].m              = restart;
      G[Bcount].overlap        = overlap;
      G[Bcount].P_Switch       = P_Switch;
      G[Bcount].max_gmres_iter = max_gmres_iter;

      G[Bcount].NCi = Soln_ptr[Bcount].NCi;
      G[Bcount].ICl = Soln_ptr[Bcount].ICl;
      G[Bcount].ICu = Soln_ptr[Bcount].ICu;
      G[Bcount].NCj = Soln_ptr[Bcount].NCj;
      G[Bcount].JCl = Soln_ptr[Bcount].JCl;
      G[Bcount].JCu = Soln_ptr[Bcount].JCu;
      G[Bcount].Nghost = Soln_ptr[Bcount].Nghost;

      G[Bcount].SolBlk_ptr = &(Soln_ptr[Bcount]);

    } /* endif */
  } /* endfor */  

  double total_norm_b;
  while (stop && count < max_newton_iter) {
    
    count=count+1;

    /* Calculate residual : use dudt[i][j][0] for residual calculations */
    for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	dUdt_Residual_Evaluation(Soln_ptr[Bcount],
				 Input_Parameters);
      } /* endif */
    } /* endfor */

    // Send boundary flux corrections at block interfaces with resolution changes.
    error_flag = Send_Conservative_Flux_Corrections(Soln_ptr, 
		                		    Soln_Block_List,
						    NUM_VAR_ADVECTDIFFUSE2D);
    if (error_flag) {
       cout << "\n AdvectDiffuse2D NKS ERROR: AdvectDiffuse2D flux correction message passing error on processor "
            << Soln_Block_List.ThisCPU
            << ".\n";
       cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
	  
    // Apply boundary flux corrections to residual to ensure that method is conservative.
    Apply_Boundary_Flux_Corrections(Soln_ptr, 
		                    Soln_Block_List);

    // Copy residual to b vector.
    for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	G[Bcount].copy_residual(G[Bcount].b, Soln_ptr[Bcount]);
	
	G[Bcount].b = -1.0 * G[Bcount].b;
      } /* endif */
    } /* endfor */

    /* Calculate norms for all blocks */    
    L2norm_current   = sqr(L2_Norm_Residual(Soln_ptr, Soln_Block_List));
    L2norm_current   = CFFC_Summation_MPI(L2norm_current);
    L2norm_current   = sqrt(L2norm_current);
    L1norm_current   = L1_Norm_Residual(Soln_ptr, Soln_Block_List);
    L1norm_current   = CFFC_Summation_MPI(L1norm_current);
    Max_norm_current = Max_Norm_Residual(Soln_ptr, Soln_Block_List);
    Max_norm_current = CFFC_Summation_MPI(Max_norm_current);

    if (count == 1) {
      L2norm_first   = L2norm_current;
      L1norm_first   = L1norm_current;
      Max_norm_first = Max_norm_current;
    }
    
    L2norm_current_n   = L2norm_current / L2norm_first; 
    L1norm_current_n   = L1norm_current / L1norm_first; 
    Max_norm_current_n = Max_norm_current / Max_norm_first;

    previous_L2norm = L2norm_current;
    
    /* Calculate L2norm Ratio */
    norm_ratio = L2norm_first / L2norm_current;
    
    Time.zero();
    Output_Progress_to_File(Progress_File,
	 		    Number_of_Startup_Interations+count-1,
			    ZERO,
			    Time,
			    L1norm_current_n,
			    L2norm_current_n,
			    Max_norm_current_n);

    /* Freeze Limiter */
    if (Freeze_Limiter == ON) {
     if (L2norm_current_n <= Freeze_Limiter_Residual_Level &&
	 limiter_check == ON)  {
       if (CFFC_Primary_MPI_Processor()) cout << "********** Apply Limiter Freezing **********" << endl;
       
       Freeze_Limiters(Soln_ptr, Soln_Block_List);

       limiter_check = OFF;
     } /* endif */
   } /* endif */
   
   /* Calculate delta t */
   dTime = CFL(Soln_ptr, Soln_Block_List, Input_Parameters);
   dTime = CFFC_Minimum_MPI(dTime); 
   if (CFFC_Primary_MPI_Processor()) cout << "dTime = " << dTime; 
   if (finite_time_step == ON) {
     dTime = norm_ratio > 1e6 ? 1e200 : min(Input_Parameters.Finite_Time_Step_Initial_CFL* 
	                                    dTime*sqr(max(ONE,norm_ratio)), 1e200);
   } else {
     dTime = 1e200;
   }

   if (CFFC_Primary_MPI_Processor()) cout << "->  modified dTime = " << dTime << endl;        
   Set_Global_TimeStep(Soln_ptr, Soln_Block_List, dTime);

    if (count <= 1) {

      if (CFFC_Primary_MPI_Processor()) cout << "  **** Create Jacobian Matrix **** " << endl;  

      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {

	  G[Bcount].dTime = dTime;

	  /* Update Jacobian Matrix */  
	  G[Bcount].A = Jacobian(Soln_ptr[Bcount], Adv_Switch, 
				 Diff_Switch, Src_Switch, overlap);

  	  /* Apply finite time step on Jacobian Matrix (Diagonally) */
	  for (j=0;j<G[Bcount].xpts*G[Bcount].ypts;j++) {
	    for (i=0;i<G[Bcount].xpts*G[Bcount].ypts;i++) {
	      if ((fabs(G[Bcount].A(j,i)) > 1e-10) && 
		  (G[Bcount].A(j,i) != 1.0) && (j==i)) {
		  G[Bcount].A.set(j,i) = G[Bcount].A(j,i) - (1.0/dTime);
		} /* endif */
	    } /* endfor */ 
	  } /* endfor */
	
	} /* endif */
      } /* endfor */
   
    } /* endif */
    
    if (L2norm_current_n > tol){   
      
      if (CFFC_Primary_MPI_Processor()) {      
	cout << "\nBegin Newton Step (Outer Iterations) = "<< setw(3) << count 
            << "   L2norm_normalized = "<< setw(12)<< L2norm_current_n << endl;
      } /* endif */
      
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {

	  /* Choose Preconditioner */
	  if (P_Switch == 1) {        
	    /* ILU */
	    CompRow_ILUPreconditioner_double M(G[Bcount].A); 
	    G[Bcount].BM = M;
	  } else { 
	    /* Diagonal */
	    DiagPreconditioner_double D(G[Bcount].A);
	    G[Bcount].BD = D;
	    
	    /* Identity */
	    if (P_Switch == 3) {
	      G[Bcount].BD.edit(1);
	    } /* endif */
	  } /* endif */
	  
	  /* Copy solution u to uo for all blocks */
	  for (j  = Soln_ptr[Bcount].JCl-Soln_ptr[Bcount].Nghost; 
	       j <= Soln_ptr[Bcount].JCu+Soln_ptr[Bcount].Nghost; ++j ) {
           for (i = Soln_ptr[Bcount].ICl-Soln_ptr[Bcount].Nghost; 
		i <= Soln_ptr[Bcount].ICu+Soln_ptr[Bcount].Nghost; ++i ) {
	      Soln_ptr[Bcount].uo[i][j] = Soln_ptr[Bcount].U[i][j].u;

	    } /* endfor */
	  } /* endfor */
	  
	} /* endif */
      } /* endfor */
      
      /* Solve system with Right Preconditioned Matrix Free GMRES */ 
      error_flag = GMRES_Algorithm(Soln_ptr,
				   Soln_Block_List,
                                   Input_Parameters,
				   G, gmrestol);
      if (error_flag) {
	delete []G;
	return 1;
      } else {
	// Update Solution : use dudt[i][j][1] to store delta u
	for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	  if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	    for (i = Soln_ptr[Bcount].ICl-Soln_ptr[Bcount].Nghost; 
		 i <= Soln_ptr[Bcount].ICu+Soln_ptr[Bcount].Nghost; ++i) {
	      for (j  = Soln_ptr[Bcount].JCl-Soln_ptr[Bcount].Nghost; 
		   j <= Soln_ptr[Bcount].JCu+Soln_ptr[Bcount].Nghost; ++j) {
		Soln_ptr[Bcount].U[i][j].u = Soln_ptr[Bcount].uo[i][j] + 
		  Soln_ptr[Bcount].dudt[i][j][1];
		
	      } /* endfor */
	    } /* endfor */
	    
	  } /* endif */
	} /* endfor */

        // Exchange solution information between neighbouring blocks.
	error_flag = Send_All_Messages(Soln_ptr, 
		        	       Soln_Block_List,
				       NUM_VAR_ADVECTDIFFUSE2D,
				       OFF);
	if (error_flag) {
	  cout << "\n AdvectDiffuse2D ERROR: AdvectDiffuse2D message passing error on processor "
	       << Soln_Block_List.ThisCPU
	       << ".\n";
	  cout.flush();
	} /* endif */
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return (error_flag);
	  
	// Apply boundary conditions for Newton step.
	BCs(Soln_ptr, 
	    Soln_Block_List,
	    Input_Parameters);

      } /* endif */
    } else {  
      stop = 0;
    } /* endif */
  } /* endwhile */

  /* Calculate Final Residual for all blocks  : */
  /*   use dudt[i][j][0] for residual calculations */
  for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      dUdt_Residual_Evaluation(Soln_ptr[Bcount],
			       Input_Parameters);
    } /* endif */
  } /* endfor */
  
  // Send boundary flux corrections at block interfaces with resolution changes.
  error_flag = Send_Conservative_Flux_Corrections(Soln_ptr, 
       	                		          Soln_Block_List,
						  NUM_VAR_ADVECTDIFFUSE2D);
  if (error_flag) {
     cout << "\n AdvectDiffuse2D NKS ERROR: AdvectDiffuse2D flux correction message passing error on processor "
          << Soln_Block_List.ThisCPU
          << ".\n";
     cout.flush();
  } /* endif */
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return (error_flag);
  	  
  // Apply boundary flux corrections to residual to ensure that method is conservative.
  Apply_Boundary_Flux_Corrections(Soln_ptr, 
		                  Soln_Block_List);

  // Copy residual to b.
  for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      G[Bcount].copy_residual(G[Bcount].b, Soln_ptr[Bcount]);
      
      G[Bcount].b = -1.0 * G[Bcount].b;
    } /* endif */
  } /* endfor */

  /* Calculate L2norm for all blocks */
  double L2norm;   
  L2norm = sqr(L2_Norm_Residual(Soln_ptr, Soln_Block_List));
  L2norm = sqrt(CFFC_Summation_MPI(L2norm));

  if (CFFC_Primary_MPI_Processor()) {  
    cout << " " << endl;
    for (int star=0;star<75;star++){cout <<"*";}
    cout.precision(20);
    if (count == 1) count = 2;
    cout << "\nEnd of Newton Steps = "<< setw(3) << count-1 << "   L2norm_normalized = "
         << setw(12)<< L2norm / L2norm_first << endl;
    
    for (int star=0;star<75;star++){cout <<"*";}
  } /* endif */

  // Update iteration counter.
  Number_of_Startup_Interations = Number_of_Startup_Interations+count-1;

  // Deallocate memories.
  for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      G[Bcount].deallocate();
    }
  }
  delete [] G;

  return (0);
  
} /* End of Newton_Krylov_Schwarz_Solver Function */ 
