

#ifndef _BLUFFBODYBURNER_INCLUDED
#include "BluffBodyBurner.h"
#endif //


//! Routine: getU -- Conduct the interpolation of the data spline.
double getU(const Vector2D &Xt, const Spline2D &U) {
  if (Xt.y < U.Xp[0].y) {
    // If the radial component of the target point is less than the
    // first point of the spline then use the first point of the spline.
    return U.Xp[0].x;
  } else if (Xt.y > U.Xp[U.np-1].y) {
    // If the radial component of the target point is greater than the
    // last point of the spline then use the last point of the spline.
    return U.Xp[U.np-1].x;
  } else {
    // Conduct spline interpolation based on the prespecified order.
    return getminX(Xt,U).x;
  }
}


//! NonreactiveVelocityField:: interpolation
//!
//! Given a point Xt this function returns the flow velocity at that
//! point interpolated from the splines representing the data at the
//! measurement stations.
Vector2D  NonreactiveVelocityField::interpolation(const Vector2D &Xt) {

  Vector2D Uflow, U1, U2;
  int ns1, ns2;

  // Determine the neighbour stations.
  if (Xt.x <= x[0]) {
    // If the target point is before the first station then just use
    // the first station.
    ns1 = 0; ns2 = 0;
  } else if (Xt.x >= x[Ns-1]) { 
    // If the target point is after the last station then just use
    // the last station.
    ns1 = Ns-1; ns2 = Ns-1;
  } else {
    // Determine the neighbour stations.
    for (int ns = 0; ns < Ns-1; ns++) {
      if (Xt.x > x[ns] && Xt.x <= x[ns+1]) {
	ns1 = ns; ns2 = ns+1;
	break;
      }
    }
  }

  // // cout << "Neighbour stations: " << ns1 << " " << ns2 << endl;

  // Conduct the interpolation.
  if (ns1 == ns2) {
    // Either before the first station or after the last station.
    Uflow.x = getU(Xt, U[ns1]);
    Uflow.y = getU(Xt, V[ns1]);
    
  } else {
    U1.x = getU(Xt, U[ns1]);
    U2.x = getU(Xt, U[ns2]);
    U1.y = getU(Xt,V[ns1]);
    U2.y = getU(Xt, V[ns2]);
    
    // Use linear interpolation to get the flow velocity.
    Uflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),U1,
				 Vector2D(x[ns2],Xt.y),U2,Xt);
  
  }

  // cout << Uflow.x <<"  "<< Uflow.y<< endl;

  // Return the interpolated velocity.
  return Uflow;

}
//! NonreactiveScalarField:: interpolation
//!
//! Given a point Xt this function returns the flow velocity at that
//! point interpolated from the splines representing the data at the
//! measurement stations.
Vector2D  NonreactiveScalarField::interpolation(const Vector2D &Xt) {

  Vector2D Mflow, M1, M2;
  int ns1, ns2;

  // Determine the neighbour stations.
  if (Xt.x <= x[0]) {
    // If the target point is before the first station then just use
    // the first station.
    ns1 = 0; ns2 = 0;
  } else if (Xt.x >= x[Ns-1]) { 
    // If the target point is after the last station then just use
    // the last station.
    ns1 = Ns-1; ns2 = Ns-1;
  } else {
    // Determine the neighbour stations.
    for (int ns = 0; ns < Ns-1; ns++) {
      if (Xt.x > x[ns] && Xt.x <= x[ns+1]) {
	ns1 = ns; ns2 = ns+1;
	break;
      }
    }
  }

  // cout << "scalar field station  " << ns1 << " " << ns2 << endl;

  // Conduct the interpolation.
  if (ns1 == ns2) {
    // Either before the first station or after the last station.
    Mflow.x = getU(Xt, M[ns1]); 
    
  } else {
    M1.x = getU(Xt,M[ns1]);
    M2.x = getU(Xt,M[ns2]);
    // Use linear interpolation to get mixture fraction.
    Mflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),M1,
				 Vector2D(x[ns2],Xt.y),M2,Xt);
  }

  // cout << Mflow.x << endl;


  return Mflow;

}
//! ReactiveVelocityField:: interpolation
//!
//! Given a point Xt this function returns the flow velocity at that
//! point interpolated from the splines representing the data at the
//! measurement stations.
Vector2D  ReactiveVelocityField_Fuel_CH4H2::interpolation(const Vector2D &Xt) {

  Vector2D Uflow, U1, U2;
  int ns1, ns2;

  // Determine the neighbour stations.
  if (Xt.x <= x[0]) {
    // If the target point is before the first station then just use
    // the first station.
    ns1 = 0; ns2 = 0;
  } else if (Xt.x >= x[Ns-1]) { 
    // If the target point is after the last station then just use
    // the last station.
    ns1 = Ns-1; ns2 = Ns-1;
  } else {
    // Determine the neighbour stations.
    for (int ns = 0; ns < Ns-1; ns++) {
      if (Xt.x > x[ns] && Xt.x <= x[ns+1]) {
	ns1 = ns; ns2 = ns+1;
	break;
      }
    }
  }

  // cout << "Neighbour stations: " << ns1 << " " << ns2 << endl;

  // Conduct the interpolation.
  if (ns1 == ns2) {
    // Either before the first station or after the last station.
    Uflow.x = getU(Xt,U[ns1]);
    Uflow.y = getU(Xt,V[ns1]);
  } else {
    U1.x = getU(Xt,U[ns1]);
    U2.x = getU(Xt,U[ns2]);
    U1.y = getU(Xt,V[ns1]);
    U2.y = getU(Xt,V[ns2]);
    
    // Use linear interpolation to get the flow velocity.
    Uflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),U1,
				 Vector2D(x[ns2],Xt.y),U2,Xt);
  }

  // cout << Uflow.x <<"  "<< Uflow.y<< endl;

  // Return the interpolated velocity.
  return Uflow;

}
//! ReactiveScalarField:: interpolation
//!
//! Given a point Xt this function returns the flow velocity at that
//! point interpolated from the splines representing the data at the
//! measurement stations.
Vector2D  ReactiveScalarField_Fuel_CH4H2::interpolation(const Vector2D &Xt, double c[]) {

  Vector2D Mflow, M1, M2;
  Vector2D Tflow, T1, T2;
  Vector2D O2flow, O21, O22;
  Vector2D N2flow, N21, N22;
  Vector2D H2flow, H21, H22;
  Vector2D H2Oflow, H2O1, H2O2;
  Vector2D Coflow, Co1, Co2;
  Vector2D CO2flow, CO21, CO22;
  Vector2D HCflow, HC1, HC2;
  Vector2D OHflow, OH1, OH2;
  int ns1, ns2;

  // Determine the neighbour stations.
  if (Xt.x <= x[0]) {
    // If the target point is before the first station then just use
    // the first station.
    ns1 = 0; ns2 = 0;
  } else if (Xt.x >= x[Ns-1]) { 
    // If the target point is after the last station then just use
    // the last station.
    ns1 = Ns-1; ns2 = Ns-1;
  } else {
    // Determine the neighbour stations.
    for (int ns = 0; ns < Ns-1; ns++) {
      if (Xt.x > x[ns] && Xt.x <= x[ns+1]) {
	ns1 = ns; ns2 = ns+1;
	break;
      }
    }
  }

  // cout << "scalar field station  " << ns1 << " " << ns2 << endl;
  
  // Conduct the interpolation.
  if (ns1 == ns2) {
    // Either before the first station or after the last station.
    
    Mflow.x  = getU(Xt,M[ns1]);
    Tflow.x  = getU(Xt,T[ns1]);
    O2flow.x = getU(Xt,O2[ns1]);
    N2flow.x = getU(Xt,N2[ns1]);
    H2flow.x = getU(Xt,H2[ns1]);
    H2Oflow.x= getU(Xt,H2O[ns1]);
    Coflow.x = getU(Xt,CO[ns1]);
    CO2flow.x= getU(Xt,CO2[ns1]);
    HCflow.x = getU(Xt,HC[ns1]);
    OHflow.x = getU(Xt,OH[ns1]);
    
  } else {
    M1.x = getU(Xt,M[ns1]);
    M2.x = getU(Xt,M[ns2]);
    T1.x = getU(Xt,T[ns1]);
    T2.x = getU(Xt,T[ns2]);
    O22.x = getU(Xt,O2[ns2]);
    O21.x = getU(Xt,O2[ns1]);
    N22.x = getU(Xt,N2[ns2]);
    N21.x = getU(Xt,N2[ns1]);
    H22.x = getU(Xt,H2[ns2]);
    H21.x = getU(Xt,H2[ns1]);
    H2O2.x = getU(Xt,H2O[ns2]);
    H2O1.x = getU(Xt,H2O[ns1]);
    Co2.x = getU(Xt,CO[ns2]);
    Co1.x = getU(Xt,CO[ns1]);
    CO22.x = getU(Xt,CO2[ns2]);
    CO21.x = getU(Xt,CO2[ns1]);
    HC2.x = getU(Xt,HC[ns2]);
    HC1.x = getU(Xt,HC[ns1]);
    OH2.x = getU(Xt,OH[ns2]);
    OH1.x = getU(Xt,OH[ns1]);

    // Use linear interpolation to get mixture fraction.
    Mflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),M1,
				 Vector2D(x[ns2],Xt.y),M2,Xt);
    Tflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),T1,
				 Vector2D(x[ns2],Xt.y),T2,Xt);
    O2flow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),O21,
				 Vector2D(x[ns2],Xt.y),O22,Xt);

    N2flow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),N21,
				  Vector2D(x[ns2],Xt.y),N22,Xt);


    H2flow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),H21,
				 Vector2D(x[ns2],Xt.y),H22,Xt);

    H2Oflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),H2O1,
				 Vector2D(x[ns2],Xt.y),H2O2,Xt);

    Coflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),Co1,
			       Vector2D(x[ns2],Xt.y),Co2,Xt);

    CO2flow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),CO21,
			       Vector2D(x[ns2],Xt.y),CO22,Xt);

    HCflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),HC1,
				  Vector2D(x[ns2],Xt.y),HC2,Xt);


    OHflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),OH1,
				  Vector2D(x[ns2],Xt.y),OH2,Xt);
  }

  c[0] =  Tflow.x;
  c[1] =  O2flow.x;
  c[2] =  N2flow.x;
  c[3] =  H2flow.x;
  c[4] =  H2Oflow.x;
  c[5] =  Coflow.x;
  c[6] =  CO2flow.x;
  c[7] =  HCflow.x;
  c[8] =  OHflow.x;
  c[9] =  Mflow.x;


  // cout << O2flow.x << endl;
//   cout << N2flow.x << endl;
//   cout << H2flow.x << endl;
//   cout << H2Oflow.x << endl;
//   cout << Coflow.x << endl;
//   cout << CO2flow.x << endl;
//   cout << HCflow.x << endl;
//   cout << OHflow.x << endl;
//   cout << Mflow.x << endl;

  
  return Mflow;

}
// //! ReactiveScalarField:: interpolation
// //!
// //! Given a point Xt this function returns the flow velocity at that
// //! point interpolated from the splines representing the data at the
// //! measurement stations.
// Vector2D  ReactiveScalarField_Fuel_CH4::interpolation(const Vector2D &Xt, double c[]) {
   
//    Vector2D Mflow, M1, M2;
//    Vector2D Tflow, T1, T2;
//    Vector2D O2flow, O21, O22;
//    Vector2D N2flow, N21, N22;
//    Vector2D H2flow, H21, H22;
//    Vector2D H2Oflow, H2O1, H2O2;
//    Vector2D Coflow, Co1, Co2;
//    Vector2D CO2flow, CO21, CO22;
//    Vector2D HCflow, HC1, HC2;
//    Vector2D OHflow, OH1, OH2;
//    int ns1, ns2;

//   // Determine the neighbour stations.
//   if (Xt.x <= x[0]) {
//     // If the target point is before the first station then just use
//     // the first station.
//     ns1 = 0; ns2 = 0;
//   } else if (Xt.x >= x[Ns-1]) { 
//     // If the target point is after the last station then just use
//     // the last station.
//     ns1 = Ns-1; ns2 = Ns-1;
//   } else {
//     // Determine the neighbour stations.
//     for (int ns = 0; ns < Ns-1; ns++) {
//       if (Xt.x > x[ns] && Xt.x <= x[ns+1]) {
// 	ns1 = ns; ns2 = ns+1;
// 	break;
//       }
//     }
//   }

//   // cout << "scalar field station  " << ns1 << " " << ns2 << endl;
  
//   // Conduct the interpolation.
//   if (ns1 == ns2) {
//     // Either before the first station or after the last station.
    
//     Mflow.x  = getU(Xt,M[ns1]);
//     Tflow.x  = getU(Xt,T[ns1]);
//     O2flow.x = getU(Xt,O2[ns1]);
//     N2flow.x = getU(Xt,N2[ns1]);
//     H2flow.x = getU(Xt,H2[ns1]);
//     H2Oflow.x= getU(Xt,H2O[ns1]);
//     Coflow.x = getU(Xt,CO[ns1]);
//     CO2flow.x= getU(Xt,CO2[ns1]);
//     HCflow.x = getU(Xt,HC[ns1]);
//     OHflow.x = getU(Xt,OH[ns1]);
    
//   } else {
//     M1.x = getU(Xt,M[ns1]);
//     M2.x = getU(Xt,M[ns2]);
//     T1.x = getU(Xt,T[ns1]);
//     T2.x = getU(Xt,T[ns2]);
//     O22.x = getU(Xt,O2[ns2]);
//     O21.x = getU(Xt,O2[ns1]);
//     N22.x = getU(Xt,N2[ns2]);
//     N21.x = getU(Xt,N2[ns1]);
//     H22.x = getU(Xt,H2[ns2]);
//     H21.x = getU(Xt,H2[ns1]);
//     H2O2.x = getU(Xt,H2O[ns2]);
//     H2O1.x = getU(Xt,H2O[ns1]);
//     Co2.x = getU(Xt,CO[ns2]);
//     Co1.x = getU(Xt,CO[ns1]);
//     CO22.x = getU(Xt,CO2[ns2]);
//     CO21.x = getU(Xt,CO2[ns1]);
//     HC2.x = getU(Xt,HC[ns2]);
//     HC1.x = getU(Xt,HC[ns1]);
//     OH2.x = getU(Xt,OH[ns2]);
//     OH1.x = getU(Xt,OH[ns1]);

//     // Use linear interpolation to get mixture fraction.
//     Mflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),M1,
// 				 Vector2D(x[ns2],Xt.y),M2,Xt);
//     Tflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),T1,
// 				 Vector2D(x[ns2],Xt.y),T2,Xt);
//     O2flow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),O21,
// 				 Vector2D(x[ns2],Xt.y),O22,Xt);

//     N2flow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),N21,
// 				  Vector2D(x[ns2],Xt.y),N22,Xt);


//     H2flow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),H21,
// 				 Vector2D(x[ns2],Xt.y),H22,Xt);

//     H2Oflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),H2O1,
// 				 Vector2D(x[ns2],Xt.y),H2O2,Xt);

//     Coflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),Co1,
// 			       Vector2D(x[ns2],Xt.y),Co2,Xt);

//     CO2flow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),CO21,
// 			       Vector2D(x[ns2],Xt.y),CO22,Xt);

//     HCflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),HC1,
// 				  Vector2D(x[ns2],Xt.y),HC2,Xt);


//     OHflow = Linear_Interpolation(Vector2D(x[ns1],Xt.y),OH1,
// 				  Vector2D(x[ns2],Xt.y),OH2,Xt);
//   }

//   c[0] =  Tflow.x;
//   c[1] =  O2flow.x;
//   c[2] =  N2flow.x;
//   c[3] =  H2flow.x;
//   c[4] =  H2Oflow.x;
//   c[5] =  Coflow.x;
//   c[6] =  CO2flow.x;
//   c[7] =  HCflow.x;
//   c[8] =  OHflow.x;
//   c[9] =  Mflow.x;


//   // cout << O2flow.x << endl;
// //   cout << N2flow.x << endl;
// //   cout << H2flow.x << endl;
// //   cout << H2Oflow.x << endl;
// //   cout << Coflow.x << endl;
// //   cout << CO2flow.x << endl;
// //   cout << HCflow.x << endl;
// //   cout << OHflow.x << endl;
// //   cout << Mflow.x << endl;

  
//   return Mflow;

// }
