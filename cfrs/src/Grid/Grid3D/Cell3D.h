/* Cell3D.h:  Header file defining 3D cell types. */

#ifndef _CELL3D_INCLUDED
#define _CELL3D_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

/* Include math macro and 3D vector header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../../../src_3D/Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../../../src_3D/Math/Vector3D.h"
#endif //_VECTOR3D_INCLUDED

/* Define the basic 3D cell class. */

/*********************************************************
 * Class: Cell3D                                         *
 *                                                       *
 * Member functions                                      *
 *     Xc       -- Return 3D vector containing location  *
 *                 of cell center (centroid).            *
 *     I        -- Return I index for cell.              *
 *     J        -- Return J index for cell.              *
 *     K        -- Return K index for cell.              *
 *     V        -- Return the volume of the cell.        *
 *     setloc   -- Set cell center location (position).  *
 *     setindex -- Set cell (i,j,k) indices.             * 
 *                                                       *
 * Member operators                                      *
 *      C -- a cell                                      *
 *                                                       *
 * C = C;                                                *
 * C == C;                                               *
 * C != C;                                               *
 * cout << C; (output function)                          *
 * cin  >> C; (input function)                           *
 *                                                       *
 *********************************************************/
class Cell3D{
  private:
  public:
    Vector3D         Xc;   // Location of cell center.
    int             I,J,K; // (i,j,k) indices for cell.
    double            V;   // Cell volume
	                   // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Cell3D(void) {
       Xc.x = ONE; Xc.y = ONE; Xc.z=ONE; I = 0; J =0; K=0; V = ONE;
    }

    Cell3D(const Cell3D &Cell) {
       Xc = Cell.Xc; I = Cell.I; J = Cell.J; K=Cell.K; V = Cell.V;
    }

    Cell3D(const Vector3D &V) {
       Xc = V;
    }

    Cell3D(const double &xx, const double &yy, const double &zz) {
      Xc.x = xx; Xc.y = yy; Xc.z=zz;
    }
    
    /* Destructor. */
    // ~Cell2D(void);
    // Use automatically generated destructor.

    /* Set cell center location. */
    void setloc(void);
    void setloc(const Cell3D &Cell);
    void setloc(const Vector3D &V);
    void setloc(const double &xx, const double &yy,const double &zz);
    
    /* Set cell (i,j,k) indices. */
    void setindex(void);
    void setindex(const Cell3D &Cell);
    void setindex(const int II, const int JJ, const int KK);

    const double & Volume(void) const { return V;}
    
    /* Assignment operator. */
    // Cell2D operator = (const Cell2D &Cell);
    // Use automatically generated assignment operator.

    /* Relational operators. */
    friend int operator ==(const Cell3D &Cell1,
			   const Cell3D &Cell2);
    friend int operator !=(const Cell3D &Cell1,
			   const Cell3D &Cell2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Cell3D &Cell);
    friend istream &operator >> (istream &in_file,
				 Cell3D &Cell);
    
};

/******************************************************************
 * Cell3D::setloc -- Set cell center location.                    *
 ******************************************************************/
inline void Cell3D::setloc(void) {
  Xc.x = ONE; Xc.y = ONE; Xc.z=ONE;
}

inline void Cell3D::setloc(const Cell3D &Cell) {
  Xc = Cell.Xc;
}

inline void Cell3D::setloc(const Vector3D &V) {
  Xc = V;
}

inline void Cell3D::setloc(const double &xx, const double &yy, const double &zz) {
  Xc.x = xx; Xc.y = yy; Xc.z=zz;
}

/******************************************************************
 * Cell3D::setindex -- Set cell (i,j) indices.                    *
 ******************************************************************/
inline void Cell3D::setindex(void) {
  I = 0; J = 0; K = 0;
}

inline void Cell3D::setindex(const Cell3D &Cell) {
  I = Cell.I; J = Cell.J; K=Cell.K;
}

inline void Cell3D::setindex(const int II, const int JJ, const int KK) {
  I = II; J = JJ; K=KK;
}

/******************************************************************
 * Cell3D -- Relational operators.                                *
 ******************************************************************/
inline int operator ==(const Cell3D &Cell1,
		       const Cell3D &Cell2) {
    return (Cell1.Xc == Cell2.Xc && Cell1.V == Cell2.V);
}

inline int operator !=(const Cell3D &Cell1,
		       const Cell3D &Cell2) {
    return (Cell1.Xc != Cell2.Xc || Cell1.V != Cell2.V);
}

/******************************************************************
 * Cell3D -- Input-output operators.                              *
 ******************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Cell3D &Cell) {
    out_file << " " << Cell.I << " " << Cell.J << " "<<Cell.K<<" "<<Cell.Xc;
    out_file.setf(ios::scientific);
    out_file << " " << Cell.V;
    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Cell3D &Cell) {
    in_file.setf(ios::skipws);
    in_file >> Cell.I >> Cell.J >>Cell.K; 
    in_file.unsetf(ios::skipws);
    in_file >> Cell.Xc;
    in_file.setf(ios::skipws);
    in_file >> Cell.V;
    in_file.unsetf(ios::skipws);
    return (in_file);
}

/* Define the basic 3D node class. */

/*********************************************************
 * Class: Node3D                                         *
 *                                                       *
 * Member functions                                      *
 *      X       -- Return 3D vector containing location  *
 *                 of node.                              *
 *      setloc  -- Set node location (position).         * 
 *                                                       *
 * Member operators                                      *
 *      N -- a node                                      *
 *                                                       *
 * N = N;                                                *
 * N == N;                                               *
 * N != N;                                               *
 * cout << N; (output function)                          *
 * cin  >> N; (input function)                           *
 *                                                       *
 *********************************************************/
class Node3D{
  private:
  public:
    Vector3D         X;   // Node location.
	                  // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Node3D(void) {
      X.x = ONE; X.y = ONE; X.z=ONE;
    }

    Node3D(const Node3D &Node) {
       X = Node.X;
    }

    Node3D(const Vector3D &V) {
       X = V;
    }

    Node3D(const double &xx, const double &yy, const double &zz) {
      X.x = xx; X.y = yy; X.z=zz;
    }
    
    /* Destructor. */
    // ~Node3D(void);
    // Use automatically generated destructor.

    /* Set cell center location. */
    void setloc(void);
    void setloc(const Node3D &Node);
    void setloc(const Vector3D &V);
    void setloc(const double &xx, const double &yy, const double &zz);
    double & x(void) { return X.x;}
    double & y(void) { return X.y;}
    double & z(void) { return X.z;}

    double tetvolume(const Vector3D V1, const Node3D N2, const Node3D N3, const Vector3D V4);

    /* Assignment operator. */
    // Node2D operator = (const Node2D &Node);
    // Use automatically generated assignment operator.

    /* Relational operators. */
    friend int operator ==(const Node3D &Node1,
			   const Node3D &Node2);
    friend int operator !=(const Node3D &Node1,
			   const Node3D &Node2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Node3D &Node);
    friend istream &operator >> (istream &in_file,
				 Node3D &Node);
    
};

/******************************************************************
 * Node3D::setloc -- Set cell center location.  *
 ******************************************************************/
inline void Node3D::setloc(void) {
  X.x = ONE; X.y = ONE; X.z = ONE;
}

inline void Node3D::setloc(const Node3D &Node) {
  X = Node.X;
}

inline void Node3D::setloc(const Vector3D &V) {
  X = V;
}

inline void Node3D::setloc(const double &xx, const double &yy, const double &zz) {
  X.x = xx; X.y = yy; X.z=zz;
}

inline double tetvolume(const Vector3D V1,const  Node3D N2,const  Node3D N3,const  Vector3D V4){
  Vector3D A,B,C; 
  A = (N2.X-V1);
  B = (N3.X-V1);
  C = (V4-V1);
  return ((fabs(A*(B^C)))/6); 
}

/******************************************************************
 * Node3D -- Relational operators.              *
 ******************************************************************/
inline int operator ==(const Node3D &Node1,
		       const Node3D &Node2) {
    return (Node1.X == Node2.X);
}

inline int operator !=(const Node3D &Node1,
		       const Node3D &Node2) {
    return (Node1.X != Node2.X);
}

/******************************************************************
 * Node2D -- Input-output operators.                *
 ******************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Node3D &Node) {
    out_file << Node.X;
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Node3D &Node) {
    in_file >> Node.X; 
    return (in_file);
}

/* Define the cell class for a 3D hexahedral mesh. */

/*********************************************************
 * Class: Cell3D_Hexa                                    *
 *                                                       *
 * Member functions                                      *
 *      xc       -- Return 3D vector containing location *
 *                  of cell center.                      *
 *      xnNWTop  -- Return 3D vector containing location *
 *                  of north-west top node.              *
 *      xnNETop  -- Return 3D vector containing location *
 *                  of north-east top node.              *
 *      xnSETop  -- Return 3D vector containing location *
 *                  of south-east top node.              *
 *      xnSWTop  -- Return 3D vector containing location *
 *                  of south-west top node.              *
 *      xnNWBot  -- Return 3D vector containing location *
 *                  of north-west bottom node.           *
 *      xnNEBot  -- Return 3D vector containing location *
 *                  of north-east bottom node.           *
 *      xnSEBot  -- Return 3D vector containing location *
 *                  of south-east bottom node.           *
 *      xnSWBot  -- Return 3D vector containing location *
 *                  of south-west bottom node.           *
 *                                                       *
 *      setnodes -- Set cell node locations.             * 
 *                                                       *
 * Member operators                                      *
 *      C -- a cell                                      *
 *                                                       *
 * C = C;                                                *
 * C == C;                                               *
 * C != C;                                               *
 * cout << C; (output function)                          *
 * cin  >> C; (input function)                           *
 *                                                       *
 *********************************************************/
class Cell3D_Hexa{
private:
public:
  Vector3D            xc;   // Location of cell center.
  Vector3D       xnNWTop;   // Location of top north-west node.
  Vector3D       xnNETop;   // Location of top north-east node.
  Vector3D       xnSETop;   // Location of top south-east node.
  Vector3D       xnSWTop;   // Location of top south-west node.
  Vector3D       xnNWBot;   // Location of bot north-west node.
  Vector3D       xnNEBot;   // Location of bot north-east node.
  Vector3D       xnSEBot;   // Location of bot south-east node.
  Vector3D       xnSWBot;   // Location of bot south-west node.
  // Made public so can access them.
		      
  /* Creation, copy, and assignment constructors. */
  Cell3D_Hexa(void);
  Cell3D_Hexa(const Cell3D_Hexa &Cell);
  //  Cell3D_Hexa(const Vector3D &V1, const Vector3D &V2, const Vector3D &V3,
  //	      const Vector3D &V4, const Vector3D &V5, const Vector3D &V6);
  Cell3D_Hexa(const double &x1, const double &y1, const double &z1,
	      const double &x2, const double &y2, const double &z2,
	      const double &x3, const double &y3, const double &z3,
	      const double &x4, const double &y4, const double &z4,
	      const double &x5, const double &y5, const double &z5,
	      const double &x6, const double &y6, const double &z6,
	      const double &x7, const double &y7, const double &z7,
	      const double &x8, const double &y8, const double &z8);
    
  /* Destructor. */
  // ~Cell3D_Hexa(void);
  // Use automatically generated destructor.

  /* Cell Volume. */
  const double V(void) const;
  
  /* Cell centroid. */
  void centroid(void);

  /* Cell face-center */
  const Vector3D xfaceN(void) const;
  const Vector3D xfaceE(void) const;
  const Vector3D xfaceS(void) const;
  const Vector3D xfaceW(void) const;
  const Vector3D xfaceTop(void) const;
  const Vector3D xfaceBot(void) const;

  /* Node location */
  Vector3D NodeNWTop(void){return xnNWTop;}
  Vector3D NodeNWTop(void) const {return xnNWTop;}

  Vector3D NodeNETop(void){return xnNETop;}
  Vector3D NodeNETop(void) const {return xnNETop;}

  Vector3D NodeSWTop(void){return xnSWTop;}
  Vector3D NodeSWTop(void) const {return xnSWTop;}

  Vector3D NodeSETop(void){return xnSETop;}
  Vector3D NodeSETop(void) const {return xnSETop;}

  Vector3D NodeNWBot(void){return xnNWBot;}
  Vector3D NodeNWBot(void) const {return xnNWBot;}

  Vector3D NodeNEBot(void){return xnNEBot;}
  Vector3D NodeNEBot(void) const {return xnNEBot;}

  Vector3D NodeSWBot(void){return xnSWBot;}
  Vector3D NodeSWBot(void) const {return xnSWBot;}

  Vector3D NodeSEBot(void){return xnSEBot;}
  Vector3D NodeSEBot(void) const {return xnSEBot;}


  /* Set cell node locations. */
  void setnodes(void);
  void setnodes(const Cell3D_Hexa &Cell);
  //void setnodes(const Vector2D &V1, const Vector2D &V2,
  //		const Vector2D &V3, const Vector2D &V4);
  void setnodes(const double &x1, const double &y1, const double &z1,
		const double &x2, const double &y2, const double &z2,
		const double &x3, const double &y3, const double &z3,
		const double &x4, const double &y4, const double &z4,
		const double &x5, const double &y5, const double &z5,
		const double &x6, const double &y6, const double &z6,
		const double &x7, const double &y7, const double &z7,
		const double &x8, const double &y8, const double &z8);

  /* Assignment operator. */
  // Cell3D_Hexa operator = (const Cell3D_Hexa &Cell);
  // Use automatically generated assignment operator.

  /* Relational operators. */
  friend int operator ==(const Cell3D_Hexa &Cell1,
			 const Cell3D_Hexa &Cell2);
  friend int operator !=(const Cell3D_Hexa &Cell1,
			 const Cell3D_Hexa &Cell2);
    
  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file,
			       const Cell3D_Hexa &Cell);
  friend istream &operator >> (istream &in_file,
			       Cell3D_Hexa &Cell);
    
};

/* Creation, copy, and assignment constructors. */
inline Cell3D_Hexa::Cell3D_Hexa(void) {
  xnNWTop.x = ZERO; xnNWTop.y = ZERO; xnNWTop.z =  ONE;
  xnNETop.x =  ONE; xnNETop.y = ZERO; xnNETop.z =  ONE;  
  xnSETop.x =  ONE; xnSETop.y =  ONE; xnSETop.z =  ONE; 
  xnSWTop.x = ZERO; xnSWTop.y =  ONE; xnSWTop.z =  ONE;

  xnNWBot.x = ZERO; xnNWBot.y = ZERO; xnNWBot.z = ZERO;
  xnNEBot.x =  ONE; xnNEBot.y = ZERO; xnNEBot.z = ZERO;  
  xnSEBot.x =  ONE; xnSEBot.y =  ONE; xnSEBot.z = ZERO; 
  xnSWBot.x = ZERO; xnSWBot.y =  ONE; xnSWBot.z = ZERO;

  centroid();
}

inline Cell3D_Hexa::Cell3D_Hexa(const Cell3D_Hexa &Cell) {
  xnNWTop = Cell.xnNWTop; xnNETop = Cell.xnNETop; xnSETop = Cell.xnSETop; xnSWTop = Cell.xnSWTop; 
  xnNWBot = Cell.xnNWBot; xnNEBot = Cell.xnNEBot; xnSEBot = Cell.xnSEBot; xnSWBot = Cell.xnSWBot; 
  xc = Cell.xc;
}

//inline Cell3D_Hexa::Cell3D_Hexa(const Vector3D &V1, const Vector3D &V2, const Vector3D &V3,
//				const Vector3D &V4, const Vector3D &V5, const Vector3D &V6) {
//  xnNW = V4; xnNE = V3; xnSE = V2; xnSW = V1; 
//  centroid();
//}

inline Cell3D_Hexa::Cell3D_Hexa(const double &x1, const double &y1, const double &z1,
				const double &x2, const double &y2, const double &z2,
				const double &x3, const double &y3, const double &z3,
				const double &x4, const double &y4, const double &z4,
				const double &x5, const double &y5, const double &z5,
				const double &x6, const double &y6, const double &z6,
				const double &x7, const double &y7, const double &z7,
				const double &x8, const double &y8, const double &z8) {

  xnNWTop.x = x4; xnNWTop.y = y4; xnNWTop.z = z4;
  xnNETop.x = x3; xnNETop.y = y3; xnNETop.z = z3;  
  xnSETop.x = x2; xnSETop.y = y2; xnSETop.z = z2; 
  xnSWTop.x = x1; xnSWTop.y = y1; xnSWTop.z = z1;

  xnNWBot.x = x8; xnNWBot.y = y8; xnNWBot.z = z8;
  xnNEBot.x = x7; xnNEBot.y = y7; xnNEBot.z = z7;  
  xnSEBot.x = x6; xnSEBot.y = y6; xnSEBot.z = z6; 
  xnSWBot.x = x5; xnSWBot.y = y5; xnSWBot.z = z5;
  
  centroid();  

  //xnNW.x = x4; xnNW.y = y4; xnNE.x = x3; xnNE.y = y3;
  //xnSE.x = x2; xnSE.y = y2; xnSW.x = x1; xnSW.y = y1;
  
}


/******************************************************************
 * Cell3D_Hexa::setnodes -- Set cell node locations.              *
 ******************************************************************/
inline void Cell3D_Hexa::setnodes(void) {
  xnNWTop.x = ZERO; xnNWTop.y = ZERO; xnNWTop.z =  ONE;
  xnNETop.x =  ONE; xnNETop.y = ZERO; xnNETop.z =  ONE;  
  xnSETop.x =  ONE; xnSETop.y =  ONE; xnSETop.z =  ONE; 
  xnSWTop.x = ZERO; xnSWTop.y =  ONE; xnSWTop.z =  ONE;

  xnNWBot.x = ZERO; xnNWBot.y = ZERO; xnNWBot.z = ZERO;
  xnNEBot.x =  ONE; xnNEBot.y = ZERO; xnNEBot.z = ZERO;  
  xnSEBot.x =  ONE; xnSEBot.y =  ONE; xnSEBot.z = ZERO; 
  xnSWBot.x = ZERO; xnSWBot.y =  ONE; xnSWBot.z = ZERO;

  centroid();

}

inline void Cell3D_Hexa::setnodes(const Cell3D_Hexa &Cell) {
  xnNWTop = Cell.xnNWTop; xnNETop = Cell.xnNETop; xnSETop = Cell.xnSETop; xnSWTop = Cell.xnSWTop; 
  xnNWBot = Cell.xnNWBot; xnNEBot = Cell.xnNEBot; xnSEBot = Cell.xnSEBot; xnSWBot = Cell.xnSWBot; 
  xc = Cell.xc;
}

//inline void Cell2D_Quad::setnodes(const Vector2D &V1, const Vector2D &V2,
//                                 const Vector2D &V3, const Vector2D &V4) {
//  xnNW = V4; xnNE = V3; xnSE = V2; xnSW = V1; 
//  centroid();
//}

inline void Cell3D_Hexa::setnodes(const double &x1, const double &y1, const double &z1,
				  const double &x2, const double &y2, const double &z2,
				  const double &x3, const double &y3, const double &z3,
				  const double &x4, const double &y4, const double &z4,
				  const double &x5, const double &y5, const double &z5,
				  const double &x6, const double &y6, const double &z6,
				  const double &x7, const double &y7, const double &z7,
				  const double &x8, const double &y8, const double &z8) {

  xnNWTop.x = x4; xnNWTop.y = y4; xnNWTop.z = z4;
  xnNETop.x = x3; xnNETop.y = y3; xnNETop.z = z3;  
  xnSETop.x = x2; xnSETop.y = y2; xnSETop.z = z2; 
  xnSWTop.x = x1; xnSWTop.y = y1; xnSWTop.z = z1;

  xnNWBot.x = x8; xnNWBot.y = y8; xnNWBot.z = z8;
  xnNEBot.x = x7; xnNEBot.y = y7; xnNEBot.z = z7;  
  xnSEBot.x = x6; xnSEBot.y = y6; xnSEBot.z = z6; 
  xnSWBot.x = x5; xnSWBot.y = y5; xnSWBot.z = z5;
  
  centroid();  

  //xnNW.x = x4; xnNW.y = y4; xnNE.x = x3; xnNE.y = y3;
  //xnSE.x = x2; xnSE.y = y2; xnSW.x = x1; xnSW.y = y1;

}

/*************************************************************************
 * Cell3D_Hexa::centroid -- Cell centre.                           *
 *************************************************************************/
inline void Cell3D_Hexa::centroid(void){

  xc = (xnNWTop+xnNETop+xnSETop+xnSWTop+xnNWBot+xnNEBot+xnSEBot+xnSWBot)/8.0;
}

/********************************************************
 * Cell3D_Hexa -- Get midpoints of faces                *
 ********************************************************/
inline const Vector3D Cell3D_Hexa::xfaceN(void) const{
  return ((xnNWTop+xnNETop+xnNEBot+xnNWBot)/FOUR);
}


inline const Vector3D Cell3D_Hexa::xfaceS(void) const{
  return ((xnSWTop+xnSETop+xnSEBot+xnSWBot)/FOUR);
}

inline const Vector3D Cell3D_Hexa::xfaceE(void) const{
  return ((xnNETop+xnNEBot+xnSEBot+xnSETop)/FOUR);
}

inline const Vector3D Cell3D_Hexa::xfaceW(void) const{
  return ((xnNWTop+xnNWBot+xnSWBot+xnSWTop)/FOUR);
}

inline const Vector3D Cell3D_Hexa::xfaceTop(void) const{
  return ((xnNETop+xnSETop+xnSWTop+xnNWTop)/FOUR);
}

inline const Vector3D Cell3D_Hexa::xfaceBot(void) const{
  return ((xnNEBot+xnSEBot+xnSWBot+xnSWBot)/FOUR);
}


/*************************************************************************
 * Cell3D_Hexa::centroid -- Cell Volume.                           *
 *************************************************************************/
inline const double Cell3D_Hexa::V(void) const {

  return (tetvolume(xc, NodeNWBot(), NodeNEBot(), xfaceBot())+
	  tetvolume(xc, NodeNWBot(), NodeSWBot(), xfaceBot())+
	  tetvolume(xc, NodeSWBot(), NodeSEBot(), xfaceBot())+
	  tetvolume(xc, NodeSEBot(), NodeNEBot(), xfaceBot())+
	  tetvolume(xc, NodeSEBot(), NodeNEBot(), xfaceE())+
	  tetvolume(xc, NodeNEBot(), NodeNETop(), xfaceE())+
	  tetvolume(xc, NodeNETop(), NodeSETop(), xfaceE())+
	  tetvolume(xc, NodeSETop(), NodeSEBot(), xfaceE())+
	  tetvolume(xc, NodeSETop(), NodeSEBot(), xfaceS())+
	  tetvolume(xc, NodeSEBot(), NodeSWBot(), xfaceS())+
	  tetvolume(xc, NodeSWBot(), NodeSWTop(), xfaceS())+
	  tetvolume(xc, NodeSWTop(), NodeSETop(), xfaceS())+
	  tetvolume(xc, NodeSWTop(), NodeSETop(), xfaceTop())+
	  tetvolume(xc, NodeSETop(), NodeNETop(), xfaceTop())+
	  tetvolume(xc, NodeNETop(), NodeNWTop(), xfaceTop())+
	  tetvolume(xc, NodeNWTop(), NodeSWTop(), xfaceTop())+
	  tetvolume(xc, NodeNWTop(), NodeSWTop(), xfaceW())+
	  tetvolume(xc, NodeSWTop(), NodeSWBot(), xfaceW())+
	  tetvolume(xc, NodeSWBot(), NodeNWBot(), xfaceW())+
	  tetvolume(xc, NodeNWBot(), NodeNWTop(), xfaceW())+
	  tetvolume(xc, NodeNWBot(), NodeNWTop(), xfaceN())+
	  tetvolume(xc, NodeNWTop(), NodeNETop(), xfaceN())+
	  tetvolume(xc, NodeNETop(), NodeNEBot(), xfaceN())+
	  tetvolume(xc, NodeNEBot(), NodeNWBot(), xfaceN())); 
}

/******************************************************************
 * Cell3D_Hexa -- Relational operators.                           *
 ******************************************************************/
inline int operator ==(const Cell3D_Hexa &Cell1,
		       const Cell3D_Hexa &Cell2) {
  return (Cell1.xnNETop == Cell2.xnNETop && Cell1.xnNWTop == Cell2.xnNWTop &&
	  Cell1.xnSWTop == Cell2.xnSWTop && Cell1.xnSETop == Cell2.xnSETop &&
	  Cell1.xnNEBot == Cell2.xnNEBot && Cell1.xnNWBot == Cell2.xnNWBot &&
	  Cell1.xnSWBot == Cell2.xnSWBot && Cell1.xnSEBot == Cell2.xnSEBot);
}

inline int operator !=(const Cell3D_Hexa &Cell1,
		       const Cell3D_Hexa &Cell2) {
  return (Cell1.xnNETop != Cell2.xnNETop || Cell1.xnNWTop != Cell2.xnNWTop ||
	  Cell1.xnSWTop != Cell2.xnSWTop || Cell1.xnSETop != Cell2.xnSETop ||
	  Cell1.xnNEBot != Cell2.xnNEBot || Cell1.xnNWBot != Cell2.xnNWBot ||
	  Cell1.xnSWBot != Cell2.xnSWBot || Cell1.xnSEBot != Cell2.xnSEBot);
}

/******************************************************************
 * Cell3D_Hexa -- Input-output operators.                         *
 ******************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Cell3D_Hexa &Cell) {
  out_file << Cell.xc;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Cell3D_Hexa &Cell) {
  in_file >> Cell.xc; 
  return (in_file);
}

#endif /* _CELL3D_INCLUDED  */

