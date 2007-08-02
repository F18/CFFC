/* Ion5Moment2DState.cc:  Subroutines for 
                          2D 5-Moment Ion Transport Model
                          Solution State Classes. */

/* Include 2D 5-moment ion transport model solution state header file. */

#ifndef _ION5MOMENT2D_STATE_INCLUDED
#include "Ion5Moment2DState.h"
#endif // _ION5MOMENT2D_STATE_INCLUDED

/*******************************************************************
 * Ion5Moment2D_pState -- Create storage and assign ion constants. *
 *******************************************************************/
double Ion5Moment2D_pState::g = GAMMA_H;
double Ion5Moment2D_pState::gm1 = GAMMA_H-ONE;
double Ion5Moment2D_pState::gm1i = ONE/(GAMMA_H-ONE);
double Ion5Moment2D_pState::R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
double Ion5Moment2D_pState::q = E_CHARGE;
double Ion5Moment2D_pState::m = MOLE_WT_H/(AVOGADRO*THOUSAND);
int    Ion5Moment2D_pState::ion = ION_H;
int    Ion5Moment2D_pState::gas = GAS_AIR;

/*******************************************************************
 * Ion5Moment2D_cState -- Create storage and assign ion constants. *
 *******************************************************************/
double Ion5Moment2D_cState::g = GAMMA_H;
double Ion5Moment2D_cState::gm1 = GAMMA_H-ONE;
double Ion5Moment2D_cState::gm1i = ONE/(GAMMA_H-ONE);
double Ion5Moment2D_cState::R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
double Ion5Moment2D_cState::q = E_CHARGE;
double Ion5Moment2D_cState::m = MOLE_WT_H/(AVOGADRO*THOUSAND);
int    Ion5Moment2D_cState::ion = ION_H;
int    Ion5Moment2D_cState::gas = GAS_AIR;

/**********************************************************
 * Routine: Riemann (Exact Riemann solver, x-direction)   *
 *                                                        *
 * This function uses a Newton-Raphson interative         *
 * procedure to obtain the exact solution to the          *
 * Riemann problem for the 2D 5-moment equations in the   *
 * x-direction, returning the intermediate state          *
 * variables along the ray x/t=0.  See Gottlieb and       *
 * Groth (1987).                                          *
 *                                                        *
 **********************************************************/
Ion5Moment2D_pState Riemann(const Ion5Moment2D_pState &Wl,
                            const Ion5Moment2D_pState &Wr) {

    int number_of_iterations;
    
    double al, ar, CL, CR, Z;
    double dml, dmr, um, vml, vmr, pm, aml, amr;
    double msl, pml, dpmldum, msr, pmr, dpmrdum;
    double vsl, vhl, vtl, vsr, vhr, vtr;

    /* Determine the left and right state sound speeds. */

    al = Wl.a();
    ar = Wr.a();

    /* Compute the left and right state Riemann invariants. */

    CL=Wl.v.x+TWO*al/Wl.gm1;
    CR=Wr.v.x-TWO*ar/Wr.gm1;

    /* Check for vacuum state. */

    if ( CL-CR <= ZERO ) {
        return (Ion5Moment2D_W_VACUUM);
    } /* endif */

    /* Make an initial estimate of the intermediate state flow
       velocity to begin the Newton-Raphson iterative solution
       procedure.  The initial guess tate velocity is made
       based on isentropic flow theory. */

    Z = (ar/al)*pow((Wl.p/Wr.p), HALF*Wl.gm1/Wl.g);
    um = (CL*Z+CR)/(ONE+Z);

    /* In the case that two rarefaction waves are present,
       then an exact solution has been found and the iterative
       procedure is not required.  Check for this. */

    if ( um >= Wl.v.x && um <= Wr.v.x ) {
        if (um >= ZERO) {
            aml = al-HALF*Wl.gm1*(um-Wl.v.x);
	    pm = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
            vhl = Wl.v.x-al;
	    vtl = um-aml;
            if (vhl >= ZERO) {
                return (Ion5Moment2D_pState(Wl));
            } else if (vtl <= ZERO) {
                vml = Wl.v.y;
 	        dml = Wl.g*pm/sqr(aml);
                return (Ion5Moment2D_pState(dml, um, vml, pm));
            } else {
	        um = (Wl.gm1*Wl.v.x+TWO*al)/(Wl.g+ONE);
                pm = Wl.p*pow((um/al), TWO*Wl.g/Wl.gm1);
                vml = Wl.v.y;
		dml = Wl.g*pm/sqr(um);
                return (Ion5Moment2D_pState(dml, um, vml, pm));
            } /* endif */
        } else {
            amr = ar+HALF*Wr.gm1*(um-Wr.v.x);
	    pm = Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1);
            vhr = Wr.v.x+ar;
	    vtr = um+amr;
            if (vhr <= ZERO) {
                return (Ion5Moment2D_pState(Wr));
            } else if (vtr >= ZERO) {
                vmr = Wr.v.y;
 	        dmr = Wr.g*pm/sqr(amr);
                return (Ion5Moment2D_pState(dmr, um, vmr, pm));
            } else {
	        um = (Wr.gm1*Wr.v.x-TWO*ar)/(Wr.g+ONE);
                pm = Wr.p*pow((-um/ar), TWO*Wr.g/Wr.gm1);
                vmr = Wr.v.y;
 	        dmr = Wr.g*pm/sqr(um);
                return (Ion5Moment2D_pState(dmr, um, vmr, pm));
            } /* endif */
        } /* end if */
    } /* endif */
    
    /* Perform the Newton-Raphson iterative procedure and solve for
       the velocity in the intermediate state.  During this iterative
       process the pressure in the intermediate state is also found. */

    number_of_iterations = 0;
    
    while (1) {
        /* Update the iteration counter. */
        number_of_iterations = number_of_iterations + 1;
 
        /* Determine solution changes for left wave. */
        if ( um < Wl.v.x ) {
            msl = (Wl.g+ONE)*(um-Wl.v.x)/(FOUR*al);
            msl = msl-sqrt(ONE+sqr(msl));
            pml = Wl.p*(ONE+Wl.g*(um-Wl.v.x)*msl/al);
            dpmldum = TWO*Wl.g*Wl.p*cube(msl)/(al*(ONE+sqr(msl)));
	} else {
            aml = al-HALF*Wl.gm1*(um-Wl.v.x);
            pml = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
	    dpmldum = -Wl.g*pml/aml;
        } /* end if */
	  
        /* Determine solution changes for right wave. */
        if ( um > Wr.v.x ) {
            msr = (Wr.g+ONE)*(um-Wr.v.x)/(FOUR*ar);
            msr = msr+sqrt(ONE+sqr(msr));
            pmr = Wr.p*(ONE+Wr.g*(um-Wr.v.x)*msr/ar);
            dpmrdum = TWO*Wr.g*Wr.p*cube(msr)/(ar*(ONE+sqr(msr)));
	} else {
            amr = ar+HALF*Wr.gm1*(um-Wr.v.x);
            pmr = Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1);
	    dpmrdum = Wr.g*pmr/amr;
        } /* end if */

	/* Check for convergence (i.e., pml=pmr). */

	if ( fabs(ONE-pml/pmr) <= TOLER) break;

	/* Compute next estimate for the intermediate
	   state velocity, um. */

	um = um-(pml-pmr)/(dpmldum-dpmrdum);
	
    } /* endwhile */

    pm = HALF*(pml+pmr);

    /* Return the intermediate state solution. */

    if ( um >= ZERO ) {
        if ( um < Wl.v.x ) {
            aml = al * sqrt(((Wl.g+ONE)+Wl.gm1*pm/Wl.p) /
			    ((Wl.g+ONE)+Wl.gm1*Wl.p/pm) );
            vsl = Wl.v.x+msl*al;
            if (vsl >= ZERO) {
                return (Ion5Moment2D_pState(Wl));                
            } else {
                vml = Wl.v.y;
 	        dml = Wl.g*pm/sqr(aml);
                return (Ion5Moment2D_pState(dml, um, vml, pm));
            } /* endif */
        } else {
            vhl = Wl.v.x-al;
	    vtl = um-aml;
            if (vhl >= ZERO) {
                return (Ion5Moment2D_pState(Wl));
            } else if (vtl <= ZERO) {
                vml = Wl.v.y;
 	        dml = Wl.g*pm/sqr(aml);
                return (Ion5Moment2D_pState(dml, um, vml, pm));
            } else {
	        um = (Wl.gm1*Wl.v.x+TWO*al)/(Wl.g+ONE);
                pm = Wl.p*pow((um/al), TWO*Wl.g/Wl.gm1);
                vml = Wl.v.y;
		dml = Wl.g*pm/sqr(um);
                return (Ion5Moment2D_pState(dml, um, vml, pm));
            } /* endif */
        } /* endif */
    } else {
        if ( um > Wr.v.x ) {
            amr = ar * sqrt(((Wr.g+ONE)+Wr.gm1*pm/Wr.p) /
			    ((Wr.g+ONE)+Wr.gm1*Wr.p/pm) );
            vsr = Wr.v.x+msr*ar;
            if (vsr <= ZERO) {
                return (Ion5Moment2D_pState(Wr));                
            } else {
                vmr = Wr.v.y;
 	        dmr = Wr.g*pm/sqr(amr);
                return (Ion5Moment2D_pState(dmr, um, vmr, pm));
            } /* endif */
        } else {
            vhr = Wr.v.x+ar;
	    vtr = um+amr;
            if (vhr <= ZERO) {
                return (Ion5Moment2D_pState(Wr));
            } else if (vtr >= ZERO) {
                vmr = Wr.v.y;
 	        dmr = Wr.g*pm/sqr(amr);
                return (Ion5Moment2D_pState(dmr, um, vmr, pm));
            } else {
	        um = (Wr.gm1*Wr.v.x-TWO*ar)/(Wr.g+ONE);
                pm = Wr.p*pow((-um/ar), TWO*Wr.g/Wr.gm1);
                vmr = Wr.v.y;
 	        dmr = Wr.g*pm/sqr(um);
                return (Ion5Moment2D_pState(dmr, um, vmr, pm));
            } /* endif */
        } /* endif */
    } /* end if */
       
}

/**********************************************************
 * Routine: Riemann_x (Exact Riemann solver, x-direction) *
 *                                                        *
 * This function uses a Newton-Raphson interative         *
 * procedure to obtain the exact solution to the          *
 * Riemann problem for the 2D 5-moment equations in the   *
 * x-direction, returning the intermediate state          *
 * variables along the ray x/t=0.  See Gottlieb and       *
 * Groth (1987).                                          *
 *                                                        *
 **********************************************************/
Ion5Moment2D_pState Riemann_x(const Ion5Moment2D_pState &Wl,
                              const Ion5Moment2D_pState &Wr) {

    int number_of_iterations;
    
    double al, ar, CL, CR, Z;
    double dml, dmr, um, vml, vmr, pm, aml, amr;
    double msl, pml, dpmldum, msr, pmr, dpmrdum;
    double vsl, vhl, vtl, vsr, vhr, vtr;

    /* Determine the left and right state sound speeds. */

    al = Wl.a();
    ar = Wr.a();

    /* Compute the left and right state Riemann invariants. */

    CL=Wl.v.x+TWO*al/Wl.gm1;
    CR=Wr.v.x-TWO*ar/Wr.gm1;

    /* Check for vacuum state. */

    if ( CL-CR <= ZERO ) {
        return (Ion5Moment2D_W_VACUUM);
    } /* endif */

    /* Make an initial estimate of the intermediate state flow
       velocity to begin the Newton-Raphson iterative solution
       procedure.  The initial guess tate velocity is made
       based on isentropic flow theory. */

    Z = (ar/al)*pow((Wl.p/Wr.p), HALF*Wl.gm1/Wl.g);
    um = (CL*Z+CR)/(ONE+Z);

    /* In the case that two rarefaction waves are present,
       then an exact solution has been found and the iterative
       procedure is not required.  Check for this. */

    if ( um >= Wl.v.x && um <= Wr.v.x ) {
        if (um >= ZERO) {
            aml = al-HALF*Wl.gm1*(um-Wl.v.x);
	    pm = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
            vhl = Wl.v.x-al;
	    vtl = um-aml;
            if (vhl >= ZERO) {
                return (Ion5Moment2D_pState(Wl));
            } else if (vtl <= ZERO) {
                vml = Wl.v.y;
 	        dml = Wl.g*pm/sqr(aml);
                return (Ion5Moment2D_pState(dml, um, vml, pm));
            } else {
	        um = (Wl.gm1*Wl.v.x+TWO*al)/(Wl.g+ONE);
                pm = Wl.p*pow((um/al), TWO*Wl.g/Wl.gm1);
                vml = Wl.v.y;
		dml = Wl.g*pm/sqr(um);
                return (Ion5Moment2D_pState(dml, um, vml, pm));
            } /* endif */
        } else {
            amr = ar+HALF*Wr.gm1*(um-Wr.v.x);
	    pm = Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1);
            vhr = Wr.v.x+ar;
	    vtr = um+amr;
            if (vhr <= ZERO) {
                return (Ion5Moment2D_pState(Wr));
            } else if (vtr >= ZERO) {
                vmr = Wr.v.y;
 	        dmr = Wr.g*pm/sqr(amr);
                return (Ion5Moment2D_pState(dmr, um, vmr, pm));
            } else {
	        um = (Wr.gm1*Wr.v.x-TWO*ar)/(Wr.g+ONE);
                pm = Wr.p*pow((-um/ar), TWO*Wr.g/Wr.gm1);
                vmr = Wr.v.y;
 	        dmr = Wr.g*pm/sqr(um);
                return (Ion5Moment2D_pState(dmr, um, vmr, pm));
            } /* endif */
        } /* end if */
    } /* endif */
    
    /* Perform the Newton-Raphson iterative procedure and solve for
       the velocity in the intermediate state.  During this iterative
       process the pressure in the intermediate state is also found. */

    number_of_iterations = 0;
    
    while (1) {
        /* Update the iteration counter. */
        number_of_iterations = number_of_iterations + 1;
 
        /* Determine solution changes for left wave. */
        if ( um < Wl.v.x ) {
            msl = (Wl.g+ONE)*(um-Wl.v.x)/(FOUR*al);
            msl = msl-sqrt(ONE+sqr(msl));
            pml = Wl.p*(ONE+Wl.g*(um-Wl.v.x)*msl/al);
            dpmldum = TWO*Wl.g*Wl.p*cube(msl)/(al*(ONE+sqr(msl)));
	} else {
            aml = al-HALF*Wl.gm1*(um-Wl.v.x);
            pml = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
	    dpmldum = -Wl.g*pml/aml;
        } /* end if */
	  
        /* Determine solution changes for right wave. */
        if ( um > Wr.v.x ) {
            msr = (Wr.g+ONE)*(um-Wr.v.x)/(FOUR*ar);
            msr = msr+sqrt(ONE+sqr(msr));
            pmr = Wr.p*(ONE+Wr.g*(um-Wr.v.x)*msr/ar);
            dpmrdum = TWO*Wr.g*Wr.p*cube(msr)/(ar*(ONE+sqr(msr)));
	} else {
            amr = ar+HALF*Wr.gm1*(um-Wr.v.x);
            pmr = Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1);
	    dpmrdum = Wr.g*pmr/amr;
        } /* end if */

	/* Check for convergence (i.e., pml=pmr). */

	if ( fabs(ONE-pml/pmr) <= TOLER) break;

	/* Compute next estimate for the intermediate
	   state velocity, um. */

	um = um-(pml-pmr)/(dpmldum-dpmrdum);
	
    } /* endwhile */

    pm = HALF*(pml+pmr);

    /* Return the intermediate state solution. */

    if ( um >= ZERO ) {
        if ( um < Wl.v.x ) {
            aml = al * sqrt(((Wl.g+ONE)+Wl.gm1*pm/Wl.p) /
			    ((Wl.g+ONE)+Wl.gm1*Wl.p/pm) );
            vsl = Wl.v.x+msl*al;
            if (vsl >= ZERO) {
                return (Ion5Moment2D_pState(Wl));                
            } else {
                vml = Wl.v.y;
 	        dml = Wl.g*pm/sqr(aml);
                return (Ion5Moment2D_pState(dml, um, vml, pm));
            } /* endif */
        } else {
            vhl = Wl.v.x-al;
	    vtl = um-aml;
            if (vhl >= ZERO) {
                return (Ion5Moment2D_pState(Wl));
            } else if (vtl <= ZERO) {
                vml = Wl.v.y;
 	        dml = Wl.g*pm/sqr(aml);
                return (Ion5Moment2D_pState(dml, um, vml, pm));
            } else {
	        um = (Wl.gm1*Wl.v.x+TWO*al)/(Wl.g+ONE);
                pm = Wl.p*pow((um/al), TWO*Wl.g/Wl.gm1);
                vml = Wl.v.y;
		dml = Wl.g*pm/sqr(um);
                return (Ion5Moment2D_pState(dml, um, vml, pm));
            } /* endif */
        } /* endif */
    } else {
        if ( um > Wr.v.x ) {
            amr = ar * sqrt(((Wr.g+ONE)+Wr.gm1*pm/Wr.p) /
			    ((Wr.g+ONE)+Wr.gm1*Wr.p/pm) );
            vsr = Wr.v.x+msr*ar;
            if (vsr <= ZERO) {
                return (Ion5Moment2D_pState(Wr));                
            } else {
                vmr = Wr.v.y;
 	        dmr = Wr.g*pm/sqr(amr);
                return (Ion5Moment2D_pState(dmr, um, vmr, pm));
            } /* endif */
        } else {
            vhr = Wr.v.x+ar;
	    vtr = um+amr;
            if (vhr <= ZERO) {
                return (Ion5Moment2D_pState(Wr));
            } else if (vtr >= ZERO) {
                vmr = Wr.v.y;
 	        dmr = Wr.g*pm/sqr(amr);
                return (Ion5Moment2D_pState(dmr, um, vmr, pm));
            } else {
	        um = (Wr.gm1*Wr.v.x-TWO*ar)/(Wr.g+ONE);
                pm = Wr.p*pow((-um/ar), TWO*Wr.g/Wr.gm1);
                vmr = Wr.v.y;
 	        dmr = Wr.g*pm/sqr(um);
                return (Ion5Moment2D_pState(dmr, um, vmr, pm));
            } /* endif */
        } /* endif */
    } /* end if */
       
}

/**********************************************************
 * Routine: Riemann_y (Exact Riemann solver, y-direction) *
 *                                                        *
 * This function uses a Newton-Raphson interative         *
 * procedure to obtain the exact solution to the          *
 * Riemann problem for the 2D 5-moment equations in the   *
 * y-direction, returning the intermediate state          *
 * variables along the ray y/t=0.  See Gottlieb and       *
 * Groth (1987).                                          *
 *                                                        *
 **********************************************************/
Ion5Moment2D_pState Riemann_y(const Ion5Moment2D_pState &Wl,
                              const Ion5Moment2D_pState &Wr) {

    int number_of_iterations;
    
    double al, ar, CL, CR, Z;
    double dml, dmr, uml, umr, vm, pm, aml, amr;
    double msl, pml, dpmldvm, msr, pmr, dpmrdvm;
    double vsl, vhl, vtl, vsr, vhr, vtr;

    /* Determine the left and right state sound speeds. */

    al = Wl.a();
    ar = Wr.a();

    /* Compute the left and right state Riemann invariants. */

    CL=Wl.v.y+TWO*al/Wl.gm1;
    CR=Wr.v.y-TWO*ar/Wr.gm1;

    /* Check for vacuum state. */

    if ( CL-CR <= ZERO ) {
        return (Ion5Moment2D_W_VACUUM);
    } /* endif */

    /* Make an initial estimate of the intermediate state flow
       velocity to begin the Newton-Raphson iterative solution
       procedure.  The initial guess tate velocity is made
       based on isentropic flow theory. */

    Z = (ar/al)*pow((Wl.p/Wr.p), HALF*Wl.gm1/Wl.g);
    vm = (CL*Z+CR)/(ONE+Z);

    /* In the case that two rarefaction waves are present,
       then an exact solution has been found and the iterative
       procedure is not required.  Check for this. */

    if ( vm >= Wl.v.y && vm <= Wr.v.y ) {
        if (vm >= ZERO) {
            aml = al-HALF*Wl.gm1*(vm-Wl.v.y);
	    pm = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
            vhl = Wl.v.y-al;
	    vtl = vm-aml;
            if (vhl >= ZERO) {
                return (Ion5Moment2D_pState(Wl));
            } else if (vtl <= ZERO) {
                uml = Wl.v.x;
 	        dml = Wl.g*pm/sqr(aml);
                return (Ion5Moment2D_pState(dml, uml, vm, pm));
            } else {
	        vm = (Wl.gm1*Wl.v.y+TWO*al)/(Wl.g+ONE);
                pm = Wl.p*pow((vm/al), TWO*Wl.g/Wl.gm1);
                uml = Wl.v.x;
		dml = Wl.g*pm/sqr(vm);
                return (Ion5Moment2D_pState(dml, uml, vm, pm));
            } /* endif */
        } else {
            amr = ar+HALF*Wr.gm1*(vm-Wr.v.y);
	    pm = Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1);
            vhr = Wr.v.y+ar;
	    vtr = vm+amr;
            if (vhr <= ZERO) {
                return (Ion5Moment2D_pState(Wr));
            } else if (vtr >= ZERO) {
                umr = Wr.v.x;
 	        dmr = Wr.g*pm/sqr(amr);
                return (Ion5Moment2D_pState(dmr, umr, vm, pm));
            } else {
	        vm = (Wr.gm1*Wr.v.y-TWO*ar)/(Wr.g+ONE);
                pm = Wr.p*pow((-vm/ar), TWO*Wr.g/Wr.gm1);
                umr = Wr.v.x;
 	        dmr = Wr.g*pm/sqr(vm);
                return (Ion5Moment2D_pState(dmr, umr, vm, pm));
            } /* endif */
        } /* end if */
    } /* endif */
    
    /* Perform the Newton-Raphson iterative procedure and solve for
       the velocity in the intermediate state.  During this iterative
       process the pressure in the intermediate state is also found. */

    number_of_iterations = 0;
    
    while (1) {
        /* Update the iteration counter. */
        number_of_iterations = number_of_iterations + 1;
 
        /* Determine solution changes for left wave. */
        if ( vm < Wl.v.y ) {
            msl = (Wl.g+ONE)*(vm-Wl.v.y)/(FOUR*al);
            msl = msl-sqrt(ONE+sqr(msl));
            pml = Wl.p*(ONE+Wl.g*(vm-Wl.v.y)*msl/al);
            dpmldvm = TWO*Wl.g*Wl.p*cube(msl)/(al*(ONE+sqr(msl)));
	} else {
            aml = al-HALF*Wl.gm1*(vm-Wl.v.y);
            pml = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
	    dpmldvm = -Wl.g*pml/aml;
        } /* end if */
	  
        /* Determine solution changes for right wave. */
        if ( vm > Wr.v.y ) {
            msr = (Wr.g+ONE)*(vm-Wr.v.y)/(FOUR*ar);
            msr = msr+sqrt(ONE+sqr(msr));
            pmr = Wr.p*(ONE+Wr.g*(vm-Wr.v.y)*msr/ar);
            dpmrdvm = TWO*Wr.g*Wr.p*cube(msr)/(ar*(ONE+sqr(msr)));
	} else {
            amr = ar+HALF*Wr.gm1*(vm-Wr.v.y);
            pmr = Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1);
	    dpmrdvm = Wr.g*pmr/amr;
        } /* end if */

	/* Check for convergence (i.e., pml=pmr). */

	if ( fabs(ONE-pml/pmr) <= TOLER) break;

	/* Compute next estimate for the intermediate
	   state velocity, vm. */

	vm = vm-(pml-pmr)/(dpmldvm-dpmrdvm);
	
    } /* endwhile */

    pm = HALF*(pml+pmr);

    /* Return the intermediate state solution. */

    if ( vm >= ZERO ) {
        if ( vm < Wl.v.y ) {
            aml = al * sqrt(((Wl.g+ONE)+Wl.gm1*pm/Wl.p) /
			    ((Wl.g+ONE)+Wl.gm1*Wl.p/pm) );
            vsl = Wl.v.y+msl*al;
            if (vsl >= ZERO) {
                return (Ion5Moment2D_pState(Wl));                
            } else {
                uml = Wl.v.x;
 	        dml = Wl.g*pm/sqr(aml);
                return (Ion5Moment2D_pState(dml, uml, vm, pm));
            } /* endif */
        } else {
            vhl = Wl.v.y-al;
	    vtl = vm-aml;
            if (vhl >= ZERO) {
                return (Ion5Moment2D_pState(Wl));
            } else if (vtl <= ZERO) {
                uml = Wl.v.x;
 	        dml = Wl.g*pm/sqr(aml);
                return (Ion5Moment2D_pState(dml, uml, vm, pm));
            } else {
	        vm = (Wl.gm1*Wl.v.y+TWO*al)/(Wl.g+ONE);
                pm = Wl.p*pow((vm/al), TWO*Wl.g/Wl.gm1);
                uml = Wl.v.x;
		dml = Wl.g*pm/sqr(vm);
                return (Ion5Moment2D_pState(dml, uml, vm, pm));
            } /* endif */
        } /* endif */
    } else {
        if ( vm > Wr.v.y ) {
            amr = ar * sqrt(((Wr.g+ONE)+Wr.gm1*pm/Wr.p) /
			    ((Wr.g+ONE)+Wr.gm1*Wr.p/pm) );
            vsr = Wr.v.y+msr*ar;
            if (vsr <= ZERO) {
                return (Ion5Moment2D_pState(Wr));                
            } else {
                umr = Wr.v.x;
 	        dmr = Wr.g*pm/sqr(amr);
                return (Ion5Moment2D_pState(dmr, umr, vm, pm));
            } /* endif */
        } else {
            vhr = Wr.v.y+ar;
	    vtr = vm+amr;
            if (vhr <= ZERO) {
                return (Ion5Moment2D_pState(Wr));
            } else if (vtr >= ZERO) {
                umr = Wr.v.x;
 	        dmr = Wr.g*pm/sqr(amr);
                return (Ion5Moment2D_pState(dmr, umr, vm, pm));
            } else {
	        vm = (Wr.gm1*Wr.v.y-TWO*ar)/(Wr.g+ONE);
                pm = Wr.p*pow((-vm/ar), TWO*Wr.g/Wr.gm1);
                umr = Wr.v.x;
 	        dmr = Wr.g*pm/sqr(vm);
                return (Ion5Moment2D_pState(dmr, umr, vm, pm));
            } /* endif */
        } /* endif */
    } /* end if */

}

/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged primitive     *
 * solution state given left and right primitive        *
 * solution variables.  See Roe (1981).                 *
 *                                                      *
 ********************************************************/
Ion5Moment2D_pState RoeAverage(const Ion5Moment2D_pState &Wl,
	      	               const Ion5Moment2D_pState &Wr) {

    double hl, hr, sdl, sdr;
    double da, ua, va, pa, aa2, ha, ga, gam1;

    /* Determine the left and right state specific enthalpies
       and square roots of the density. */

    hl = Wl.h();
    hr = Wr.h();
    sdl = sqrt(Wl.d);
    sdr = sqrt(Wr.d);

    /* Determine the appropriate Roe averages. */

    ga = Wl.g;
    gam1 = Wl.gm1;
    da = sdl*sdr;
    ua = (sdl*Wl.v.x+sdr*Wr.v.x)/(sdl+sdr);
    va = (sdl*Wl.v.y+sdr*Wr.v.y)/(sdl+sdr);
    ha = (sdl*hl+sdr*hr)/(sdl+sdr);
    aa2 = gam1*(ha-HALF*(sqr(ua)+sqr(va)));
    pa = da*aa2/ga;

    /* Return the Roe-averged state. */

    return (Ion5Moment2D_pState(da, ua, va, pa));
       
}

/********************************************************
 * Routine: Reflect                                     *
 *                                                      *
 * This function returns the reflected solution state   *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 ********************************************************/
Ion5Moment2D_pState Reflect(const Ion5Moment2D_pState &W,
	      	            const Vector2D &norm_dir) {

    double dr, ur, vr, pr, u, v;
    double cos_angle, sin_angle;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and calculate the primitive
       solution state variables in the local rotated frame
       defined by the unit normal vector. */

    dr = W.d;
    ur = W.v.x*cos_angle +
         W.v.y*sin_angle;
    vr = - W.v.x*sin_angle +
           W.v.y*cos_angle;
    pr = W.p;

    /* Reflect the normal velocity in the rotated frame. */

    ur = -ur;

    /* Rotate back to the original Cartesian reference frame. */

    u = ur*cos_angle - vr*sin_angle;
    v = ur*sin_angle + vr*cos_angle;

    /* Return the reflected state. */

    return (Ion5Moment2D_pState(dr, u, v, pr));
       
}

/********************************************************
 * Routine: BC_Characteristic (Characteristic-Based     *
 *                             Boundary Condition)      *
 *                                                      *
 * This function returns the boundary solution state    *
 * for a given direction given the primitive solution   *
 * state on the interior of the boundary, Wi, the       *
 * desired flow state to be imposed at the boundary,    *
 * Wo, and the unit normal vector in the direction of   *
 * interest. The boundary flow state is based on the    *
 * solution of a simplified Riemann problem as          *
 * described by Gottlieb and Groth (1999).  The         *
 * imposition of the boundary-data respects the         *
 * directions of propogation for the solution           *
 * characteristics at the boundary and thereby ensures  *
 * that the boundary conditions are not ill-posed (i.e.,*
 * the boundary data is not under- or over-determined). *
 *                                                      *
 ********************************************************/
Ion5Moment2D_pState BC_Characteristic(const Ion5Moment2D_pState &Wi,
                                      const Ion5Moment2D_pState &Wo,
	      	                      const Vector2D &norm_dir) {

    Ion5Moment2D_pState Wi_rotated, Wo_rotated;
    char pattern;
    double mi, poi, pai, pbi, mab, mac1, mac2, mbc1,
           mc1c2, mbd, mad;
    double de, ue, ve, pe, ae, ue_rotated, ve_rotated;
    double cos_angle, sin_angle;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate interior and 
       imposed boundary solution states in the local rotated 
       frame defined by the unit normal vector. */

    Wi_rotated.d = Wi.d;
    Wi_rotated.v.x = Wi.v.x*cos_angle +
                     Wi.v.y*sin_angle;
    Wi_rotated.v.y = - Wi.v.x*sin_angle +
                       Wi.v.y*cos_angle;
    Wi_rotated.p = Wi.p;

    Wo_rotated.d = Wo.d;
    Wo_rotated.v.x = Wo.v.x*cos_angle +
                     Wo.v.y*sin_angle;
    Wo_rotated.v.y = - Wo.v.x*sin_angle +
                       Wo.v.y*cos_angle;
    Wo_rotated.p = Wo.p;

    /* Determine the Mach number at the interior node and the
       imposed boundary to interior node pressure ratio in
       the rotated frame. */

    mi = Wi_rotated.v.x/Wi_rotated.a();
    poi = Wo_rotated.p/Wi_rotated.p;

    /* Identify the wave pattern between the interior and 
       boundary node in the rotated frame. */

    if (poi >= ONE) {
       pai = pow(HALF*(Wi_rotated.g+ONE), 
                 TWO*Wi_rotated.g*Wi_rotated.gm1i);
       pbi = pow(HALF*(Wi_rotated.g+ONE)+
                 HALF*Wi_rotated.gm1*(Wo_rotated.a()/Wi_rotated.a()), 
                 TWO*Wi_rotated.g*Wi_rotated.gm1i);
       if (poi <= pai) {
          mab = ONE;
          // Pattern A, pi <= po <= pa.
          if (mi >= mab) {
	     pattern='a';
	     goto impose_boundary_conditions;
          } /* endif */
          mbc1 = TWO*Wi_rotated.gm1i*(pow(poi, 
                 HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
          // Pattern B, pi <= po <= pa.
          if (mi >= mbc1) {
	     pattern='b';
	     goto impose_boundary_conditions;
          } /* endif */
          mc1c2 = mbc1 - Wo_rotated.a()/Wi_rotated.a();
          // Pattern C1, pi <= po <= pa.
          if (mi >= mc1c2) {
	     pattern='c';
	     goto impose_boundary_conditions;
          // Pattern C2, pi <= po <= pa.
          } else {
	     pattern='c';
	     goto impose_boundary_conditions;
          } /* endif */
       } else if (poi <= pbi) {
          mac1 = ONE;
          // Pattern A, pa < po <= pb.
          if (mi >= mac1) {
	     pattern='a';
	     goto impose_boundary_conditions;
          } /* endif */
          mbc1 = TWO*Wi_rotated.gm1i*(pow(poi, 
                 HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
          mc1c2 = mbc1 - Wo_rotated.a()/Wi_rotated.a();
          // Pattern C1, pa < po <= pb.
          if (mi >= mc1c2) {
	     pattern='c';
	     goto impose_boundary_conditions;
          // Pattern C2, pa < po <= pb.
          } else {
	     pattern='c';
	     goto impose_boundary_conditions;
          } /* endif */
       } else {
          mac2 = ONE;
          // Pattern A, po > pb.
          if (mi >= mac2) {
	     pattern='a';
	     goto impose_boundary_conditions;
          // Pattern C2, po > pb.
          } else {
	     pattern='c';
	     goto impose_boundary_conditions;
          } /* endif */
       } /* endif */
    } else {
      mad = ONE;
      // Pattern A, po < pi.
      if (mi >= mad) {
         pattern='a';
         goto impose_boundary_conditions;
      } /* endif */
      mbd = (Wi_rotated.g+ONE)*Wi_rotated.gm1i*pow(poi, 
            HALF*Wi_rotated.gm1/Wi_rotated.g)-
            TWO*Wi_rotated.gm1i;
      // Pattern D, po < pi.
      if (mi >= mbd) {
         pattern='d';
         goto impose_boundary_conditions;
      } /* endif */
      mbc1 = TWO*Wi_rotated.gm1i*(pow(poi, 
             HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
      // Pattern B, po < pi.
      if (mi >= mbc1) {
         pattern='b';
         goto impose_boundary_conditions;
      } /* endif */
      mc1c2 = mbc1 - Wo_rotated.a()/Wi_rotated.a();
      // Pattern C1, po < pi.
      if (mi >= mc1c2) {
         pattern='c';
         goto impose_boundary_conditions;
      // Pattern C2, po < pi.
      } else {
         pattern='c';
         goto impose_boundary_conditions;
      } /* endif */
    } /* endif */

    /* Prescribe the appropriate boundary conditions
       in the rotated frame. */

    impose_boundary_conditions: ;
    switch(pattern) {
      case 'a' :
        de = Wi_rotated.d;
        ue_rotated = Wi_rotated.v.x;
        ve_rotated = Wi_rotated.v.y;
        pe = Wi_rotated.p;
        break;        
      case 'b' :
        pe = Wo_rotated.p;
        ue_rotated = Wi_rotated.v.x -
             TWO*Wi_rotated.a()*Wi_rotated.gm1i*
             (pow(poi, HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
        ve_rotated = Wi_rotated.v.y;
        de = Wi_rotated.d*pow(poi, ONE/Wi_rotated.g);
        break;
      case 'c' :
        pe = Wo_rotated.p;
        ue_rotated = Wi_rotated.v.x -
             TWO*Wi_rotated.a()*Wi_rotated.gm1i*
             (pow(poi, HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
        ve_rotated = Wo_rotated.v.y;
        de = Wo_rotated.d;
        break;
      case 'd' :
        pe = Wi_rotated.p*pow((TWO/(Wi_rotated.g+ONE)+
                      (Wi_rotated.gm1/(Wi_rotated.g+ONE))*mi), 
                      TWO*Wi_rotated.g*Wi_rotated.gm1i);
        ue_rotated = Wi_rotated.v.x -
             TWO*Wi_rotated.a()*Wi_rotated.gm1i*
             (pow(pe/Wi_rotated.p, HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
        ve_rotated = Wi_rotated.v.y;
        de = Wi_rotated.d*pow(pe/Wi_rotated.p, ONE/Wi_rotated.g);
        break;
      default:
        de = Wo_rotated.d;
        ue_rotated = Wo_rotated.v.x;
        ve_rotated = Wo_rotated.v.y;
        pe = Wo_rotated.p;
        break;
    } /* endswitch */

    /* Rotate the resulting boundary state back to the original 
       Cartesian reference frame. */

    ue = ue_rotated*cos_angle - ve_rotated*sin_angle;
    ve = ue_rotated*sin_angle + ve_rotated*cos_angle;

    /* Return boundary solution state. */

    return (Ion5Moment2D_pState(de, ue, ve, pe));
       
}

/********************************************************
 * Routine: BC_Characteristic_Pressure                  *
 *   (Characteristic-Based Boundary Condition with      *
 *    Static Pressure Specified Whenever Possible)      *
 *                                                      *
 * This function returns the boundary solution state    *
 * for a given direction given the primitive solution   *
 * state on the interior of the boundary, Wi, the       *
 * desired flow state to be imposed at the boundary,    *
 * Wo, and the unit normal vector in the direction of   *
 * interest. A simplified characteristic analysis is    *
 * used to specify the boundary flow state in which the *
 * static pressure is specified whenever possible.  The *
 * imposition of the boundary-data respects the         *
 * directions of propogation for the solution           *
 * characteristics at the boundary and thereby ensures  *
 * that the boundary conditions are not ill-posed (i.e.,*
 * the boundary data is not under- or over-determined). *
 * The following procedure is adopted:                  *
 *                                                      *
 * 1) for supersonic outflow: constant extrapolation    *
 *    is employed to specify boundary conditions        *
 *    using known interior solution state,              *
 * 2) for subsonic outflow: boundary conditions         *
 *    specified by employing a 1D unsteady isentropic   *
 *    wave approximation to match boundary pressure,    *
 * 3) for subsonic inflow: boundary conditions          *
 *    specified by employing a 1D unsteady isentropic   *
 *    wave approximation to match both boundary state   *
 *    pressure and sound speed,                         *
 * 4) for supersonic inflow: the known boundary state   *
 *    is used to specify the boundary state.            *
 *                                                      *
 ********************************************************/
Ion5Moment2D_pState BC_Characteristic_Pressure(const Ion5Moment2D_pState &Wi,
                                               const Ion5Moment2D_pState &Wo,
	      	                               const Vector2D &norm_dir) {

    Ion5Moment2D_pState Wi_rotated, Wo_rotated;
    double mi, db, ub, vb, pb, ab, ub_rotated, vb_rotated;
    double cos_angle, sin_angle;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate interior and 
       imposed boundary solution states in the local rotated 
       frame defined by the unit normal vector. */

    Wi_rotated.d = Wi.d;
    Wi_rotated.v.x = Wi.v.x*cos_angle +
                     Wi.v.y*sin_angle;
    Wi_rotated.v.y = - Wi.v.x*sin_angle +
                       Wi.v.y*cos_angle;
    Wi_rotated.p = Wi.p;

    Wo_rotated.d = Wo.d;
    Wo_rotated.v.x = Wo.v.x*cos_angle +
                     Wo.v.y*sin_angle;
    Wo_rotated.v.y = - Wo.v.x*sin_angle +
                       Wo.v.y*cos_angle;
    Wo_rotated.p = Wo.p;

    /* Determine the Mach number at the interior node. */

    mi = Wi_rotated.v.x/Wi_rotated.a();

    /* Boundary condition for supersonic outflow. */

    if (mi >= ONE) {
       db = Wi_rotated.d;
       ub_rotated = Wi_rotated.v.x;
       vb_rotated = Wi_rotated.v.y;
       pb = Wi_rotated.p;
   
    /* Boundary condition for subsonic outflow. 
       Pressure specified. */

    } else if (mi >= ZERO) {
       pb = Wo_rotated.p;
       db = Wi_rotated.d*pow(pb/Wi_rotated.p, ONE/Wi_rotated.g);
       ab = sqrt(Wi_rotated.g*pb/db);
       ub_rotated = Wi_rotated.v.x + TWO*(Wi_rotated.a()-ab)*Wi_rotated.gm1i;
       vb_rotated = Wi_rotated.v.y;

    /* Boundary condition for subsonic inflow. 
       Pressure specified. */

    } else if (mi >= -ONE) {
       pb = Wo_rotated.p;
       db = Wo_rotated.d;
       vb_rotated = Wo_rotated.v.y;
       ab = sqrt(Wo_rotated.g*pb/db);
       ub_rotated = Wi_rotated.v.x + TWO*(Wi_rotated.a()-ab)*Wo_rotated.gm1i;

    /* Boundary condition for supersonic inflow.  */

    } else {
       db = Wo_rotated.d;
       ub_rotated = Wo_rotated.v.x;
       vb_rotated = Wo_rotated.v.y;
       pb = Wo_rotated.p;
    } /* endif */

    /* Rotate the resulting boundary state back to the original 
       Cartesian reference frame. */

    ub = ub_rotated*cos_angle - vb_rotated*sin_angle;
    vb = ub_rotated*sin_angle + vb_rotated*cos_angle;

    /* Return boundary solution state. */

    return (Ion5Moment2D_pState(db, ub, vb, pb));
       
}

/********************************************************
 * Routine: BC_Characteristic_Mach_Number               *
 *   (Characteristic-Based Boundary Condition with      *
 *    Flow Mach Number Specified Whenever Possible)     *
 *                                                      *
 * This function returns the boundary solution state    *
 * for a given direction given the primitive solution   *
 * state on the interior of the boundary, Wi, the       *
 * desired flow state to be imposed at the boundary,    *
 * Wo, and the unit normal vector in the direction of   *
 * interest. A simplified characteristic analysis is    *
 * used to specify the boundary flow state in which the *
 * flow Mach number is specified whenever possible.     *
 * The imposition of the boundary-data respects the     *
 * directions of propogation for the solution           *
 * characteristics at the boundary and thereby ensures  *
 * that the boundary conditions are not ill-posed (i.e.,*
 * the boundary data is not under- or over-determined). *
 * The following procedure is adopted:                  *
 *                                                      *
 * 1) for supersonic outflow: constant extrapolation    *
 *    is employed to specify boundary conditions        *
 *    using known interior solution state,              *
 * 2) for subsonic outflow: boundary conditions         *
 *    specified by employing a 1D unsteady isentropic   *
 *    wave approximation to match boundary Mach number, *
 * 3) for subsonic inflow: boundary conditions          *
 *    specified by employing a 1D unsteady isentropic   *
 *    wave approximation to match both boundary state   *
 *    Mach number and density,                          *
 * 4) for supersonic inflow: the known boundary state   *
 *    is used to specify the boundary state.            *
 *                                                      *
 ********************************************************/
Ion5Moment2D_pState BC_Characteristic_Mach_Number(const Ion5Moment2D_pState &Wi,
                                                  const Ion5Moment2D_pState &Wo,
	      	                                  const Vector2D &norm_dir) {

    Ion5Moment2D_pState Wi_rotated, Wo_rotated;
    double mi, mb, db, ub, vb, pb, ab, ub_rotated, vb_rotated;
    double cos_angle, sin_angle;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate interior and 
       imposed boundary solution states in the local rotated 
       frame defined by the unit normal vector. */

    Wi_rotated.d = Wi.d;
    Wi_rotated.v.x = Wi.v.x*cos_angle +
                     Wi.v.y*sin_angle;
    Wi_rotated.v.y = - Wi.v.x*sin_angle +
                       Wi.v.y*cos_angle;
    Wi_rotated.p = Wi.p;

    Wo_rotated.d = Wo.d;
    Wo_rotated.v.x = Wo.v.x*cos_angle +
                     Wo.v.y*sin_angle;
    Wo_rotated.v.y = - Wo.v.x*sin_angle +
                       Wo.v.y*cos_angle;
    Wo_rotated.p = Wo.p;

    /* Determine the Mach number at the interior node. */

    mi = Wi_rotated.v.x/Wi_rotated.a();

    /* Boundary condition for supersonic outflow. */

    if (mi >= ONE) {
       db = Wi_rotated.d;
       ub_rotated = Wi_rotated.v.x;
       vb_rotated = Wi_rotated.v.y;
       pb = Wi_rotated.p;
   
    /* Boundary condition for subsonic outflow. 
       Mach number specified. */

    } else if (mi >= ZERO) {
       db = Wi_rotated.d;
       mb = Wo_rotated.v.x/Wo_rotated.a();
       ab= (Wi_rotated.v.x+TWO*Wi_rotated.a()*Wi_rotated.gm1i)/
           (mb+TWO*Wi_rotated.gm1i);
       pb = db*ab*ab/Wi_rotated.g;
       ub_rotated = mb*ab;
       vb_rotated = Wi_rotated.v.y;

    /* Boundary condition for subsonic inflow. 
       Pressure specified. */

    } else if (mi >= -ONE) {
       db = Wo_rotated.d;
       mb = Wo_rotated.v.x/Wo_rotated.a();
       ab= (Wi_rotated.v.x+TWO*Wi_rotated.a()*Wi_rotated.gm1i)/
           (mb+TWO*Wi_rotated.gm1i);
       pb = db*ab*ab/Wo_rotated.g;
       ub_rotated = mb*ab;
       vb_rotated = Wo_rotated.v.y;

    /* Boundary condition for supersonic inflow.  */

    } else {
       db = Wo_rotated.d;
       ub_rotated = Wo_rotated.v.x;
       vb_rotated = Wo_rotated.v.y;
       pb = Wo_rotated.p;
    } /* endif */

    /* Rotate the resulting boundary state back to the original 
       Cartesian reference frame. */

    ub = ub_rotated*cos_angle - vb_rotated*sin_angle;
    vb = ub_rotated*sin_angle + vb_rotated*cos_angle;

    /* Return boundary solution state. */

    return (Ion5Moment2D_pState(db, ub, vb, pb));
       
}

/********************************************************
 * Routine: WaveSpeedPos                                *
 *                                                      *
 * This function returns the positive parts of the      *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
Ion5Moment2D_pState WaveSpeedPos(const Ion5Moment2D_pState &lambdas_a,
                                 const Ion5Moment2D_pState &lambdas_l,
                                 const Ion5Moment2D_pState &lambdas_r) {
  return (Ion5Moment2D_pState(HALF*(lambdas_a[1]+fabs(lambdas_a[1])),
                              HALF*(lambdas_a[2]+fabs(lambdas_a[2])),
                              HALF*(lambdas_a[3]+fabs(lambdas_a[3])),
                              HALF*(lambdas_a[4]+fabs(lambdas_a[4]))));
}

/********************************************************
 * Routine: WaveSpeedNeg                                *
 *                                                      *
 * This function returns the negative parts of the      *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
Ion5Moment2D_pState WaveSpeedNeg(const Ion5Moment2D_pState &lambdas_a,
                                 const Ion5Moment2D_pState &lambdas_l,
                                 const Ion5Moment2D_pState &lambdas_r) {
  return (Ion5Moment2D_pState(HALF*(lambdas_a[1]-fabs(lambdas_a[1])),
                              HALF*(lambdas_a[2]-fabs(lambdas_a[2])),
                              HALF*(lambdas_a[3]-fabs(lambdas_a[3])),
                              HALF*(lambdas_a[4]-fabs(lambdas_a[4]))));
}

/********************************************************
 * Routine: WaveSpeedAbs                                *
 *                                                      *
 * This function returns the absolute values of the     *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
Ion5Moment2D_pState WaveSpeedAbs(const Ion5Moment2D_pState &lambdas_a,
                                 const Ion5Moment2D_pState &lambdas_l,
                                 const Ion5Moment2D_pState &lambdas_r) {
  return (Ion5Moment2D_pState(fabs(lambdas_a[1]),
                              fabs(lambdas_a[2]),
                              fabs(lambdas_a[3]),
                              fabs(lambdas_a[4])));
}

/********************************************************
 * Routine: HartenFixPos (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the positive parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
Ion5Moment2D_pState HartenFixPos(const Ion5Moment2D_pState &lambdas_a,
                                 const Ion5Moment2D_pState &lambdas_l,
                                 const Ion5Moment2D_pState &lambdas_r) {
  return (Ion5Moment2D_pState(HartenFixPos(lambdas_a[1],
                                           lambdas_l[1],
                                           lambdas_r[1]),
                              HALF*(lambdas_a[2]+fabs(lambdas_a[2])),
                              HALF*(lambdas_a[3]+fabs(lambdas_a[3])),
                              HartenFixPos(lambdas_a[4],
                                           lambdas_l[4],
                                           lambdas_r[4])));
}

/********************************************************
 * Routine: HartenFixNeg (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the negative parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
Ion5Moment2D_pState HartenFixNeg(const Ion5Moment2D_pState &lambdas_a,
                            const Ion5Moment2D_pState &lambdas_l,
                            const Ion5Moment2D_pState &lambdas_r) {
  return (Ion5Moment2D_pState(HartenFixNeg(lambdas_a[1],
                                           lambdas_l[1],
                                           lambdas_r[1]),
                              HALF*(lambdas_a[2]-fabs(lambdas_a[2])),
                              HALF*(lambdas_a[3]-fabs(lambdas_a[3])),
                              HartenFixNeg(lambdas_a[4],
                                           lambdas_l[4],
                                           lambdas_r[4])));
}

/********************************************************
 * Routine: HartenFixAbs (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the absolute values of the     *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
Ion5Moment2D_pState HartenFixAbs(const Ion5Moment2D_pState &lambdas_a,
                            const Ion5Moment2D_pState &lambdas_l,
                            const Ion5Moment2D_pState &lambdas_r) {
  return (Ion5Moment2D_pState(HartenFixAbs(lambdas_a[1],
                                           lambdas_l[1],
                                           lambdas_r[1]),
                              fabs(lambdas_a[2]),
                              fabs(lambdas_a[3]),
                              HartenFixAbs(lambdas_a[4],
                                           lambdas_l[4],
                                           lambdas_r[4])));
}

/*********************************************************
 * Routine: FluxGodunov (Godunov's flux function,        *
 *                       x-direction)                    *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by solving exactly the Riemann        *
 * problem associated with the two states.  See Godunov  *
 * (1959).                                               *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxGodunov(const Ion5Moment2D_pState &Wl,
	      	                const Ion5Moment2D_pState &Wr) {
    return (F(Riemann(Wl, Wr)));
}

Ion5Moment2D_cState FluxGodunov(const Ion5Moment2D_cState &Ul,
	      	                const Ion5Moment2D_cState &Ur) {
    return (F(Riemann(Ul.W(), Ur.W())));
}

/*********************************************************
 * Routine: FluxGodunov_x (Godunov's flux function,      *
 *                         x-direction)                  *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by solving exactly the Riemann        *
 * problem associated with the two states.  See Godunov  *
 * (1959).                                               *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxGodunov_x(const Ion5Moment2D_pState &Wl,
	      	                  const Ion5Moment2D_pState &Wr) {
    return (Fx(Riemann_x(Wl, Wr)));
}

Ion5Moment2D_cState FluxGodunov_x(const Ion5Moment2D_cState &Ul,
	      	                  const Ion5Moment2D_cState &Ur) {
    return (Fx(Riemann_x(Ul.W(), Ur.W())));
}

/*********************************************************
 * Routine: FluxGodunov_y (Godunov's flux function,      *
 *                         y-direction)                  *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the y-direction given left and right         *
 * solution states by solving exactly the Riemann        *
 * problem associated with the two states.  See Godunov  *
 * (1959).                                               *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxGodunov_y(const Ion5Moment2D_pState &Wl,
	      	                  const Ion5Moment2D_pState &Wr) {
    return (Fy(Riemann_y(Wl, Wr)));
}

Ion5Moment2D_cState FluxGodunov_y(const Ion5Moment2D_cState &Ul,
	      	                  const Ion5Moment2D_cState &Ur) {
    return (Fy(Riemann_y(Ul.W(), Ur.W())));
}

/*********************************************************
 * Routine: FluxGodunov_n (Godunov's flux function,      *
 *                         n-direction)                  *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then solving exactly the Riemann    *
 * problem associated with the two states in the rotated *
 * frame.  See Godunov (1959).                           *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxGodunov_n(const Ion5Moment2D_pState &Wl,
	      	                  const Ion5Moment2D_pState &Wr,
                                  const Vector2D &norm_dir) {

    double cos_angle, sin_angle;
    Ion5Moment2D_pState Wl_rotated, Wr_rotated;
    Ion5Moment2D_cState Flux, Flux_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */

    Wl_rotated.d = Wl.d;
    Wl_rotated.v.x = Wl.v.x*cos_angle +
                     Wl.v.y*sin_angle;
    Wl_rotated.v.y = - Wl.v.x*sin_angle +
                       Wl.v.y*cos_angle;
    Wl_rotated.p = Wl.p;

    Wr_rotated.d = Wr.d;
    Wr_rotated.v.x = Wr.v.x*cos_angle +
                     Wr.v.y*sin_angle;
    Wr_rotated.v.y = - Wr.v.x*sin_angle +
                       Wr.v.y*cos_angle;
    Wr_rotated.p = Wr.p;

    /* Solve the Riemann problem in the rotated frame
       and evaluate the intermediate state solution 
       flux. */

    Flux_rotated = Fx(Riemann_x(Wl_rotated, Wr_rotated));

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */

    Flux.d = Flux_rotated.d;
    Flux.dv.x = Flux_rotated.dv.x*cos_angle -
                Flux_rotated.dv.y*sin_angle;
    Flux.dv.y = Flux_rotated.dv.x*sin_angle +
                Flux_rotated.dv.y*cos_angle;
    Flux.E = Flux_rotated.E;

    return (Flux);

}

Ion5Moment2D_cState FluxGodunov_n(const Ion5Moment2D_cState &Ul,
	      	                  const Ion5Moment2D_cState &Ur,
                                  const Vector2D &norm_dir) {
    return (FluxGodunov_n(Ul.W(), Ur.W(), norm_dir));
}

/*********************************************************
 * Routine: FluxRoe (Roe's flux function, x-direction)   *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the "linearized" approximate *
 * Riemann solver of Roe for the two states.  See Roe    *
 * (1981).                                               *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxRoe(const Ion5Moment2D_pState &Wl,
	      	            const Ion5Moment2D_pState &Wr) {

    int i;
    Ion5Moment2D_pState Wa, dWrl, 
                   wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux;

    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the primitive solution states. */

    dWrl = Wr-Wl;

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda();
    lambdas_r = Wr.lambda();
    lambdas_a = Wa.lambda();

    /* Determine the intermediate state flux. */

    if (Wa.v.x >= ZERO) {
        Flux = Wl.F();
        wavespeeds = HartenFixNeg(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_ION5MOMENT2D ; ++i ) {
           if (wavespeeds[i] < ZERO)
              Flux += wavespeeds[i]*(Wa.lp(i)*dWrl)*Wa.rc(i);
        } /* endfor */
    } else {
        Flux = Wr.F();
        wavespeeds = HartenFixPos(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_ION5MOMENT2D ; ++i ) {
           if (wavespeeds[i] > ZERO)
              Flux -= wavespeeds[i]*(Wa.lp(i)*dWrl)*Wa.rc(i);
        } /* endfor */
    } /* endif */

    /* Return solution flux. */
    
    return (Flux);
    
}

Ion5Moment2D_cState FluxRoe(const Ion5Moment2D_cState &Ul,
	      	            const Ion5Moment2D_cState &Ur) {
   return (FluxRoe(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxRoe_x (Roe's flux function, x-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the "linearized" approximate *
 * Riemann solver of Roe for the two states.  See Roe    *
 * (1981).                                               *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxRoe_x(const Ion5Moment2D_pState &Wl,
	      	              const Ion5Moment2D_pState &Wr) {

    int i;
    Ion5Moment2D_pState Wa, dWrl, 
                   wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux;

    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the primitive solution states. */

    dWrl = Wr-Wl;

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_x();
    lambdas_r = Wr.lambda_x();
    lambdas_a = Wa.lambda_x();

    /* Determine the intermediate state flux. */

    if (Wa.v.x >= ZERO) {
        Flux = Wl.Fx();
        wavespeeds = HartenFixNeg(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_ION5MOMENT2D ; ++i ) {
           if (wavespeeds[i] < ZERO)
              Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
        } /* endfor */
    } else {
        Flux = Wr.Fx();
        wavespeeds = HartenFixPos(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_ION5MOMENT2D ; ++i ) {
           if (wavespeeds[i] > ZERO)
              Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
        } /* endfor */
    } /* endif */

    /* Return solution flux. */
    
    return (Flux);
    
}

Ion5Moment2D_cState FluxRoe_x(const Ion5Moment2D_cState &Ul,
	      	              const Ion5Moment2D_cState &Ur) {
   return (FluxRoe_x(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxRoe_y (Roe's flux function, y-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the y-direction given left and right         *
 * solution states by using the "linearized" approximate *
 * Riemann solver of Roe for the two states.  See Roe    *
 * (1981).                                               *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxRoe_y(const Ion5Moment2D_pState &Wl,
	      	              const Ion5Moment2D_pState &Wr) {

    int i;
    Ion5Moment2D_pState Wa, dWrl, 
                   wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux;

    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the primitive solution states. */

    dWrl = Wr-Wl;

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_y();
    lambdas_r = Wr.lambda_y();
    lambdas_a = Wa.lambda_y();

    /* Determine the intermediate state flux. */

    if (Wa.v.y >= ZERO) {
        Flux = Wl.Fy();
        wavespeeds = HartenFixNeg(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_ION5MOMENT2D ; ++i ) {
           if (wavespeeds[i] < ZERO)
              Flux += wavespeeds[i]*(Wa.lp_y(i)*dWrl)*Wa.rc_y(i);
        } /* endfor */
    } else {
        Flux = Wr.Fy();
        wavespeeds = HartenFixPos(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_ION5MOMENT2D ; ++i ) {
           if (wavespeeds[i] > ZERO)
              Flux -= wavespeeds[i]*(Wa.lp_y(i)*dWrl)*Wa.rc_y(i);
        } /* endfor */
    } /* endif */

    /* Return solution flux. */
    
    return (Flux);
    
}

Ion5Moment2D_cState FluxRoe_y(const Ion5Moment2D_cState &Ul,
	      	              const Ion5Moment2D_cState &Ur) {
   return (FluxRoe_y(Ul.W(), Ur.W()));
}  

/*********************************************************
 * Routine: FluxRoe_n (Roe's flux function, n-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then by using the "linearized"      *
 * approximate Riemann solver of Roe to specify the flux *
 * in terms of the rotated solution states.  See Roe     *
 * (1981).                                               *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxRoe_n(const Ion5Moment2D_pState &Wl,
	      	              const Ion5Moment2D_pState &Wr,
                              const Vector2D &norm_dir) {

    double cos_angle, sin_angle;
    Ion5Moment2D_pState Wl_rotated, Wr_rotated;
    Ion5Moment2D_cState Flux, Flux_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */

    Wl_rotated.d = Wl.d;
    Wl_rotated.v.x = Wl.v.x*cos_angle +
                     Wl.v.y*sin_angle;
    Wl_rotated.v.y = - Wl.v.x*sin_angle +
                       Wl.v.y*cos_angle;
    Wl_rotated.p = Wl.p;

    Wr_rotated.d = Wr.d;
    Wr_rotated.v.x = Wr.v.x*cos_angle +
                     Wr.v.y*sin_angle;
    Wr_rotated.v.y = - Wr.v.x*sin_angle +
                       Wr.v.y*cos_angle;
    Wr_rotated.p = Wr.p;

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

    Flux_rotated = FluxRoe_x(Wl_rotated, Wr_rotated);

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */

    Flux.d = Flux_rotated.d;
    Flux.dv.x = Flux_rotated.dv.x*cos_angle -
                Flux_rotated.dv.y*sin_angle;
    Flux.dv.y = Flux_rotated.dv.x*sin_angle +
                Flux_rotated.dv.y*cos_angle;
    Flux.E = Flux_rotated.E;

    return (Flux);

}

Ion5Moment2D_cState FluxRoe_n(const Ion5Moment2D_cState &Ul,
	      	              const Ion5Moment2D_cState &Ur,
                              const Vector2D &norm_dir) {
    return (FluxRoe_n(Ul.W(), Ur.W(), norm_dir));
}

/*********************************************************
 * Routine: FluxRusanov (Rusanov flux function,          *
 *                       x-direction)                    *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the Rusanov approximation    *
 * for the fluxes.  See Rusanov (1964).                  *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxRusanov(const Ion5Moment2D_pState &Wl,
	      	                const Ion5Moment2D_pState &Wr) {

    int i;
    double wavespeed_max;
    Ion5Moment2D_pState Wa, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda();
    lambdas_r = Wr.lambda();
    lambdas_a = Wa.lambda();

    /* Determine the intermediate state flux. */

    Flux = HALF*(Wl.F()+Wr.F());
    wavespeeds = HartenFixAbs(lambdas_a,
                              lambdas_l,
                              lambdas_r);

    wavespeed_max = wavespeeds[1];
    for ( i = 2 ; i <= NUM_VAR_ION5MOMENT2D ; ++i ) {
        wavespeed_max = max(wavespeed_max, wavespeeds[i]);
    } /* endfor */

    Flux -= HALF*wavespeed_max*dUrl;

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxRusanov(const Ion5Moment2D_cState &Ul,
	      	                const Ion5Moment2D_cState &Ur) {
   return (FluxRusanov(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxRusanov_x (Rusanov flux function,        *
 *                         x-direction)                  *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the Rusanov approximation    *
 * for the fluxes.  See Rusanov (1964).                  *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxRusanov_x(const Ion5Moment2D_pState &Wl,
	      	                  const Ion5Moment2D_pState &Wr) {

    int i;
    double wavespeed_max;
    Ion5Moment2D_pState Wa, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_x();
    lambdas_r = Wr.lambda_x();
    lambdas_a = Wa.lambda_x();

    /* Determine the intermediate state flux. */

    Flux = HALF*(Wl.Fx()+Wr.Fx());
    wavespeeds = HartenFixAbs(lambdas_a,
                              lambdas_l,
                              lambdas_r);

    wavespeed_max = wavespeeds[1];
    for ( i = 2 ; i <= NUM_VAR_ION5MOMENT2D ; ++i ) {
        wavespeed_max = max(wavespeed_max, wavespeeds[i]);
    } /* endfor */

    Flux -= HALF*wavespeed_max*dUrl;

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxRusanov_x(const Ion5Moment2D_cState &Ul,
	      	                  const Ion5Moment2D_cState &Ur) {
   return (FluxRusanov_x(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxRusanov_y (Rusanov flux function,        *
 *                         y-direction)                  *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the y-direction given left and right         *
 * solution states by using the Rusanov approximation    *
 * for the fluxes.  See Rusanov (1964).                  *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxRusanov_y(const Ion5Moment2D_pState &Wl,
	      	                  const Ion5Moment2D_pState &Wr) {

    int i;
    double wavespeed_max;
    Ion5Moment2D_pState Wa, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_y();
    lambdas_r = Wr.lambda_y();
    lambdas_a = Wa.lambda_y();

    /* Determine the intermediate state flux. */

    Flux = HALF*(Wl.Fy()+Wr.Fy());
    wavespeeds = HartenFixAbs(lambdas_a,
                              lambdas_l,
                              lambdas_r);

    wavespeed_max = wavespeeds[1];
    for ( i = 2 ; i <= NUM_VAR_ION5MOMENT2D ; ++i ) {
        wavespeed_max = max(wavespeed_max, wavespeeds[i]);
    } /* endfor */

    Flux -= HALF*wavespeed_max*dUrl;

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxRusanov_y(const Ion5Moment2D_cState &Ul,
	      	                  const Ion5Moment2D_cState &Ur) {
   return (FluxRusanov_y(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxRusanov_n (Rusanov flux function,        *
 *                         n-direction)                  *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then by using the Rusanov           *
 * approximation for the intermediate state flux in      *
 * in of the rotated solution states.  See Rusanov       *
 * (1964).                                               *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxRusanov_n(const Ion5Moment2D_pState &Wl,
	      	                  const Ion5Moment2D_pState &Wr,
                                  const Vector2D &norm_dir) {

    double cos_angle, sin_angle;
    Ion5Moment2D_pState Wl_rotated, Wr_rotated;
    Ion5Moment2D_cState Flux, Flux_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */

    Wl_rotated.d = Wl.d;
    Wl_rotated.v.x = Wl.v.x*cos_angle +
                     Wl.v.y*sin_angle;
    Wl_rotated.v.y = - Wl.v.x*sin_angle +
                       Wl.v.y*cos_angle;
    Wl_rotated.p = Wl.p;

    Wr_rotated.d = Wr.d;
    Wr_rotated.v.x = Wr.v.x*cos_angle +
                     Wr.v.y*sin_angle;
    Wr_rotated.v.y = - Wr.v.x*sin_angle +
                       Wr.v.y*cos_angle;
    Wr_rotated.p = Wr.p;

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

    Flux_rotated = FluxRusanov_x(Wl_rotated, Wr_rotated);

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */

    Flux.d = Flux_rotated.d;
    Flux.dv.x = Flux_rotated.dv.x*cos_angle -
                Flux_rotated.dv.y*sin_angle;
    Flux.dv.y = Flux_rotated.dv.x*sin_angle +
                Flux_rotated.dv.y*cos_angle;
    Flux.E = Flux_rotated.E;

    return (Flux);

}

Ion5Moment2D_cState FluxRusanov_n(const Ion5Moment2D_cState &Ul,
	      	                  const Ion5Moment2D_cState &Ur,
                                  const Vector2D &norm_dir) {
    return (FluxRusanov_n(Ul.W(), Ur.W(), norm_dir));
}

/*********************************************************
 * Routine: FluxHLLE (Harten-Lax-van Leer flux function, *
 *                    x-direction)                       *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the so-called Harten-Lax-    *
 * van Leer approximation for the fluxes.  See Harten,   *
 * Lax, van Leer (1983).                                 *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxHLLE(const Ion5Moment2D_pState &Wl,
	      	             const Ion5Moment2D_pState &Wr) {

    double wavespeed_l, wavespeed_r;
    Ion5Moment2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

   /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda();
    lambdas_r = Wr.lambda();
    lambdas_a = Wa.lambda();

    /* Determine the intermediate state flux. */

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_ION5MOMENT2D],
                      lambdas_a[NUM_VAR_ION5MOMENT2D]);
 
    wavespeed_l = min(wavespeed_l, ZERO);
    wavespeed_r = max(wavespeed_r, ZERO);

    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else {
        Flux =   ( (wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
                  +(wavespeed_l*wavespeed_r)*dUrl)/
                 (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxHLLE(const Ion5Moment2D_cState &Ul,
	      	             const Ion5Moment2D_cState &Ur) {
   return (FluxHLLE(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLE_x (Harten-Lax-van Leer flux         *
 *                      function, x-direction)           *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the so-called Harten-Lax-    *
 * van Leer approximation for the fluxes.  See Harten,   *
 * Lax, van Leer (1983).                                 *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxHLLE_x(const Ion5Moment2D_pState &Wl,
	      	               const Ion5Moment2D_pState &Wr) {

    double wavespeed_l, wavespeed_r;
    Ion5Moment2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

   /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_x();
    lambdas_r = Wr.lambda_x();
    lambdas_a = Wa.lambda_x();

    /* Determine the intermediate state flux. */

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_ION5MOMENT2D],
                      lambdas_a[NUM_VAR_ION5MOMENT2D]);
 
    wavespeed_l = min(wavespeed_l, ZERO);
    wavespeed_r = max(wavespeed_r, ZERO);

    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else {
        Flux =   ( (wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx())
                  +(wavespeed_l*wavespeed_r)*dUrl)/
                 (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxHLLE_x(const Ion5Moment2D_cState &Ul,
	      	               const Ion5Moment2D_cState &Ur) {
   return (FluxHLLE_x(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLE_y (Harten-Lax-van Leer flux         *
 *                      function, y-direction)           *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the y-direction given left and right         *
 * solution states by using the so-called Harten-Lax-    *
 * van Leer approximation for the fluxes.  See Harten,   *
 * Lax, van Leer (1983).                                 *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxHLLE_y(const Ion5Moment2D_pState &Wl,
	      	               const Ion5Moment2D_pState &Wr) {

    double wavespeed_l, wavespeed_r;
    Ion5Moment2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

   /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_y();
    lambdas_r = Wr.lambda_y();
    lambdas_a = Wa.lambda_y();

    /* Determine the intermediate state flux. */

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_ION5MOMENT2D],
                      lambdas_a[NUM_VAR_ION5MOMENT2D]);
 
    wavespeed_l = min(wavespeed_l, ZERO);
    wavespeed_r = max(wavespeed_r, ZERO);

    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else {
        Flux =   ( (wavespeed_r*Wl.Fy()-wavespeed_l*Wr.Fy())
                  +(wavespeed_l*wavespeed_r)*dUrl)/
                 (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxHLLE_y(const Ion5Moment2D_cState &Ul,
	      	               const Ion5Moment2D_cState &Ur) {
   return (FluxHLLE_y(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLE_n (Harten-Lax-van Leer flux         *
 *                      function, n-direction)           *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then by using the so-called         *
 * Harten-Lax-van Leer approximation to specify the      *
 * intermediate state fluxes in terms of the rotated     *
 * solution states.  See Harten, Lax, van Leer (1983).   *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxHLLE_n(const Ion5Moment2D_pState &Wl,
	      	               const Ion5Moment2D_pState &Wr,
                               const Vector2D &norm_dir) {

    double cos_angle, sin_angle;
    Ion5Moment2D_pState Wl_rotated, Wr_rotated;
    Ion5Moment2D_cState Flux, Flux_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */

    Wl_rotated.d = Wl.d;
    Wl_rotated.v.x = Wl.v.x*cos_angle +
                     Wl.v.y*sin_angle;
    Wl_rotated.v.y = - Wl.v.x*sin_angle +
                       Wl.v.y*cos_angle;
    Wl_rotated.p = Wl.p;

    Wr_rotated.d = Wr.d;
    Wr_rotated.v.x = Wr.v.x*cos_angle +
                     Wr.v.y*sin_angle;
    Wr_rotated.v.y = - Wr.v.x*sin_angle +
                       Wr.v.y*cos_angle;
    Wr_rotated.p = Wr.p;

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

    Flux_rotated = FluxHLLE_x(Wl_rotated, Wr_rotated);

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */

    Flux.d = Flux_rotated.d;
    Flux.dv.x = Flux_rotated.dv.x*cos_angle -
                Flux_rotated.dv.y*sin_angle;
    Flux.dv.y = Flux_rotated.dv.x*sin_angle +
                Flux_rotated.dv.y*cos_angle;
    Flux.E = Flux_rotated.E;

    return (Flux);

}

Ion5Moment2D_cState FluxHLLE_n(const Ion5Moment2D_cState &Ul,
	      	               const Ion5Moment2D_cState &Ur,
                               const Vector2D &norm_dir) {
    return (FluxHLLE_n(Ul.W(), Ur.W(), norm_dir));
}

/*********************************************************
 * Routine: FluxLinde (Timur Linde's flux function,      *
 *                     x-direction)                      *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the Linde approximation for  *
 * the fluxes.  See Linde (1998).                        *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxLinde(const Ion5Moment2D_pState &Wl,
                              const Ion5Moment2D_pState &Wr) {

    double wavespeed_l, wavespeed_r, wavespeed_m, da, ca, dU, alpha;
    Ion5Moment2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux, dFrl, dUrl, dFwave;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda();
    lambdas_r = Wr.lambda();
    lambdas_a = Wa.lambda();

    /* Determine the intermediate state flux. */

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_ION5MOMENT2D],
                      lambdas_a[NUM_VAR_ION5MOMENT2D]);

    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else {
        dUrl = Wr.U()-Wl.U();
        dFrl = Wr.F()-Wl.F();
        wavespeed_m = Wa.v.x;
        da = Wa.d;
        ca = Wa.a();
        dU = fabs(dUrl.d)/da+
             fabs(dUrl.dv.x)/(da*ca)+
             fabs(dUrl.dv.y)/(da*ca)+
             fabs(dUrl.E)/(da*ca*ca);
        if (dU <= TOLER) {
            alpha = ZERO;
        } else {
            dU = ONE/dU;
            dFwave = dFrl - wavespeed_m*dUrl;
            alpha = ONE - (fabs(dFwave.d)/(da*ca)+
                           fabs(dFwave.dv.x)/(da*ca*ca)+
                           fabs(dFwave.dv.y)/(da*ca*ca)+
                           fabs(dFwave.E)/(da*ca*ca*ca))*dU;
 	    alpha = max(ZERO, alpha);
        } /* endif */
 
        Flux =   ( (wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
                  +(wavespeed_l*wavespeed_r)*
		   (ONE-(ONE-max(wavespeed_m/wavespeed_r,
                                 wavespeed_m/wavespeed_l))*alpha)*dUrl)/
                   (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxLinde(const Ion5Moment2D_cState &Ul,
	      	              const Ion5Moment2D_cState &Ur) {
   return (FluxLinde(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxLinde_x (Timur Linde's flux function,    *
 *                       x-direction)                    *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the Linde approximation for  *
 * the fluxes.  See Linde (1998).                        *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxLinde_x(const Ion5Moment2D_pState &Wl,
                                const Ion5Moment2D_pState &Wr) {

    double wavespeed_l, wavespeed_r, wavespeed_m, da, ca, dU, alpha;
    Ion5Moment2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux, dFrl, dUrl, dFwave;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_x();
    lambdas_r = Wr.lambda_x();
    lambdas_a = Wa.lambda_x();

    /* Determine the intermediate state flux. */

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_ION5MOMENT2D],
                      lambdas_a[NUM_VAR_ION5MOMENT2D]);

    if (wavespeed_l >= ZERO) {
        Flux = Wl.Fx();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.Fx();
    } else {
        dUrl = Wr.U()-Wl.U();
        dFrl = Wr.Fx()-Wl.Fx();
        wavespeed_m = Wa.v.x;
        da = Wa.d;
        ca = Wa.a();
        dU = fabs(dUrl.d)/da+
             fabs(dUrl.dv.x)/(da*ca)+
             fabs(dUrl.dv.y)/(da*ca)+
             fabs(dUrl.E)/(da*ca*ca);
        if (dU <= TOLER) {
            alpha = ZERO;
        } else {
            dU = ONE/dU;
            dFwave = dFrl - wavespeed_m*dUrl;
            alpha = ONE - (fabs(dFwave.d)/(da*ca)+
                           fabs(dFwave.dv.x)/(da*ca*ca)+
                           fabs(dFwave.dv.y)/(da*ca*ca)+
                           fabs(dFwave.E)/(da*ca*ca*ca))*dU;
 	    alpha = max(ZERO, alpha);
        } /* endif */
 
        Flux =   ( (wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx())
                  +(wavespeed_l*wavespeed_r)*
		   (ONE-(ONE-max(wavespeed_m/wavespeed_r,
                                 wavespeed_m/wavespeed_l))*alpha)*dUrl)/
                   (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxLinde_x(const Ion5Moment2D_cState &Ul,
	      	                const Ion5Moment2D_cState &Ur) {
   return (FluxLinde_x(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxLinde_y (Timur Linde's flux function,    *
 *                       y-direction)                    *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the y-direction given left and right         *
 * solution states by using the Linde approximation for  *
 * the fluxes.  See Linde (1998).                        *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxLinde_y(const Ion5Moment2D_pState &Wl,
                                const Ion5Moment2D_pState &Wr) {

    double wavespeed_l, wavespeed_r, wavespeed_m, da, ca, dU, alpha;
    Ion5Moment2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Ion5Moment2D_cState Flux, dFrl, dUrl, dFwave;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_y();
    lambdas_r = Wr.lambda_y();
    lambdas_a = Wa.lambda_y();

    /* Determine the intermediate state flux. */

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_ION5MOMENT2D],
                      lambdas_a[NUM_VAR_ION5MOMENT2D]);

    if (wavespeed_l >= ZERO) {
        Flux = Wl.Fy();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.Fy();
    } else {
        dUrl = Wr.U()-Wl.U();
        dFrl = Wr.Fy()-Wl.Fy();
        wavespeed_m = Wa.v.y;
        da = Wa.d;
        ca = Wa.a();
        dU = fabs(dUrl.d)/da+
             fabs(dUrl.dv.x)/(da*ca)+
             fabs(dUrl.dv.y)/(da*ca)+
             fabs(dUrl.E)/(da*ca*ca);
        if (dU <= TOLER) {
            alpha = ZERO;
        } else {
            dU = ONE/dU;
            dFwave = dFrl - wavespeed_m*dUrl;
            alpha = ONE - (fabs(dFwave.d)/(da*ca)+
                           fabs(dFwave.dv.x)/(da*ca*ca)+
                           fabs(dFwave.dv.y)/(da*ca*ca)+
                           fabs(dFwave.E)/(da*ca*ca*ca))*dU;
 	    alpha = max(ZERO, alpha);
        } /* endif */
 
        Flux =   ( (wavespeed_r*Wl.Fy()-wavespeed_l*Wr.Fy())
                  +(wavespeed_l*wavespeed_r)*
		   (ONE-(ONE-max(wavespeed_m/wavespeed_r,
                                 wavespeed_m/wavespeed_l))*alpha)*dUrl)/
                   (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxLinde_y(const Ion5Moment2D_cState &Ul,
	      	                const Ion5Moment2D_cState &Ur) {
   return (FluxLinde_y(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxLinde_n (Timur Linde's flux function,    *
 *                       n-direction)                    *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then by using the Linde             *
 * approximation to specify the intermediate state flux  *
 * in terms of the rotated solution states.  See Linde   *
 * (1998).                                               *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxLinde_n(const Ion5Moment2D_pState &Wl,
	      	                const Ion5Moment2D_pState &Wr,
                                const Vector2D &norm_dir) {

    double cos_angle, sin_angle;
    Ion5Moment2D_pState Wl_rotated, Wr_rotated;
    Ion5Moment2D_cState Flux, Flux_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */

    Wl_rotated.d = Wl.d;
    Wl_rotated.v.x = Wl.v.x*cos_angle +
                     Wl.v.y*sin_angle;
    Wl_rotated.v.y = - Wl.v.x*sin_angle +
                       Wl.v.y*cos_angle;
    Wl_rotated.p = Wl.p;

    Wr_rotated.d = Wr.d;
    Wr_rotated.v.x = Wr.v.x*cos_angle +
                     Wr.v.y*sin_angle;
    Wr_rotated.v.y = - Wr.v.x*sin_angle +
                       Wr.v.y*cos_angle;
    Wr_rotated.p = Wr.p;

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

    Flux_rotated = FluxLinde_x(Wl_rotated, Wr_rotated);

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */

    Flux.d = Flux_rotated.d;
    Flux.dv.x = Flux_rotated.dv.x*cos_angle -
                Flux_rotated.dv.y*sin_angle;
    Flux.dv.y = Flux_rotated.dv.x*sin_angle +
                Flux_rotated.dv.y*cos_angle;
    Flux.E = Flux_rotated.E;

    return (Flux);

}

Ion5Moment2D_cState FluxLinde_n(const Ion5Moment2D_cState &Ul,
	      	                const Ion5Moment2D_cState &Ur,
                                const Vector2D &norm_dir) {
    return (FluxLinde_n(Ul.W(), Ur.W(), norm_dir));
}

/*********************************************************
 * Routine: FluxHLLC (Tito Toro's flux function,         *
 *                    x-direction)                       *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the HLLC approximation for   *
 * the fluxes.  See Toro, Spruce, and Speares (1994).    *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxHLLC(const Ion5Moment2D_pState &Wl,
                             const Ion5Moment2D_pState &Wr) {

    double wavespeed_l, wavespeed_r, wavespeed_m;
    double al, ar, CL, CR, Z, ql, qr;
    double um, pm, aml, amr;

    Ion5Moment2D_cState Flux, Uml, Umr;
    
    /* Determine the left, intermediate, and right
       wave speeds. */

    al = Wl.a();
    ar = Wr.a();
    CL=Wl.v.x+TWO*al/Wl.gm1;
    CR=Wr.v.x-TWO*ar/Wr.gm1;
    Z = (ar/al)*pow((Wl.p/Wr.p), HALF*Wl.gm1/Wl.g);
    um = (CL*Z+CR)/(ONE+Z);
    aml = al-HALF*Wl.gm1*(um-Wl.v.x);
    pm = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
    amr = ar+HALF*Wr.gm1*(um-Wr.v.x);
    pm = HALF*(pm + Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1));

    if (pm/Wl.p <= ONE) {
      ql = ONE;
    } else {
      ql = sqrt(ONE+HALF*((Wl.g+ONE)/Wl.g)*(pm/Wl.p));
    } /* endif */
    wavespeed_l = Wl.v.x - ql*al;

    if (pm/Wr.p <= ONE) {
      qr = ONE;
    } else {
      qr = sqrt(ONE+HALF*((Wr.g+ONE)/Wr.g)*(pm/Wr.p));
    } /* endif */
    wavespeed_r = Wr.v.x + qr*ar;

    wavespeed_m = um;

    /* Determine the intermediate state flux. */

    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else  if (wavespeed_m >= ZERO) {
        Uml = Ion5Moment2D_cState(Wl.d, Wl.d*wavespeed_m, Wl.d*Wl.v.y,
          		     Wl.E()+Wl.d*(wavespeed_m-Wl.v.x)*(wavespeed_m+
        		     Wl.p/(Wl.d*(wavespeed_l-Wl.v.x))))*
              ((wavespeed_l-Wl.v.x)/(wavespeed_l-wavespeed_m));
        Flux = Wl.F()+wavespeed_l*(Uml-Wl.U());
    } else {
        Umr = Ion5Moment2D_cState(Wr.d, Wr.d*wavespeed_m, Wr.d*Wr.v.y,
          		     Wr.E()+Wr.d*(wavespeed_m-Wr.v.x)*(wavespeed_m+
        		     Wr.p/(Wr.d*(wavespeed_r-Wr.v.x))))*
              ((wavespeed_r-Wr.v.x)/(wavespeed_r-wavespeed_m));
        Flux = Wr.F()+wavespeed_r*(Umr-Wr.U());
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxHLLC(const Ion5Moment2D_cState &Ul,
	      	             const Ion5Moment2D_cState &Ur) {
   return (FluxHLLC(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLC_x (Tito Toro's flux function,       *
 *                      x-direction)                     *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the HLLC approximation for   *
 * the fluxes.  See Toro, Spruce, and Speares (1994).    *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxHLLC_x(const Ion5Moment2D_pState &Wl,
                               const Ion5Moment2D_pState &Wr) {

    double wavespeed_l, wavespeed_r, wavespeed_m;
    double al, ar, CL, CR, Z, ql, qr;
    double um, pm, aml, amr;

    Ion5Moment2D_cState Flux, Uml, Umr;
    
    /* Determine the left, intermediate, and right
       wave speeds. */

    al = Wl.a();
    ar = Wr.a();
    CL=Wl.v.x+TWO*al/Wl.gm1;
    CR=Wr.v.x-TWO*ar/Wr.gm1;
    Z = (ar/al)*pow((Wl.p/Wr.p), HALF*Wl.gm1/Wl.g);
    um = (CL*Z+CR)/(ONE+Z);
    aml = al-HALF*Wl.gm1*(um-Wl.v.x);
    pm = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
    amr = ar+HALF*Wr.gm1*(um-Wr.v.x);
    pm = HALF*(pm + Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1));

    if (pm/Wl.p <= ONE) {
      ql = ONE;
    } else {
      ql = sqrt(ONE+HALF*((Wl.g+ONE)/Wl.g)*(pm/Wl.p));
    } /* endif */
    wavespeed_l = Wl.v.x - ql*al;

    if (pm/Wr.p <= ONE) {
      qr = ONE;
    } else {
      qr = sqrt(ONE+HALF*((Wr.g+ONE)/Wr.g)*(pm/Wr.p));
    } /* endif */
    wavespeed_r = Wr.v.x + qr*ar;

    wavespeed_m = um;

    /* Determine the intermediate state flux. */

    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else  if (wavespeed_m >= ZERO) {
        Uml = Ion5Moment2D_cState(Wl.d, Wl.d*wavespeed_m, Wl.d*Wl.v.y,
          		     Wl.E()+Wl.d*(wavespeed_m-Wl.v.x)*(wavespeed_m+
        		     Wl.p/(Wl.d*(wavespeed_l-Wl.v.x))))*
              ((wavespeed_l-Wl.v.x)/(wavespeed_l-wavespeed_m));
        Flux = Wl.F()+wavespeed_l*(Uml-Wl.U());
    } else {
        Umr = Ion5Moment2D_cState(Wr.d, Wr.d*wavespeed_m, Wr.d*Wr.v.y,
          		     Wr.E()+Wr.d*(wavespeed_m-Wr.v.x)*(wavespeed_m+
        		     Wr.p/(Wr.d*(wavespeed_r-Wr.v.x))))*
              ((wavespeed_r-Wr.v.x)/(wavespeed_r-wavespeed_m));
        Flux = Wr.F()+wavespeed_r*(Umr-Wr.U());
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxHLLC_x(const Ion5Moment2D_cState &Ul,
	      	               const Ion5Moment2D_cState &Ur) {
   return (FluxHLLC_x(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLC_y (Tito Toro's flux function,       *
 *                      y-direction)                     *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the y-direction given left and right         *
 * solution states by using the HLLC approximation for   *
 * the fluxes.  See Toro, Spruce, and Speares (1994).    *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxHLLC_y(const Ion5Moment2D_pState &Wl,
                               const Ion5Moment2D_pState &Wr) {

    double wavespeed_l, wavespeed_r, wavespeed_m;
    double al, ar, CL, CR, Z, ql, qr;
    double vm, pm, aml, amr;

    Ion5Moment2D_cState Flux, Uml, Umr;
    
    /* Determine the left, intermediate, and right
       wave speeds. */

    al = Wl.a();
    ar = Wr.a();
    CL=Wl.v.y+TWO*al/Wl.gm1;
    CR=Wr.v.y-TWO*ar/Wr.gm1;
    Z = (ar/al)*pow((Wl.p/Wr.p), HALF*Wl.gm1/Wl.g);
    vm = (CL*Z+CR)/(ONE+Z);
    aml = al-HALF*Wl.gm1*(vm-Wl.v.y);
    pm = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
    amr = ar+HALF*Wr.gm1*(vm-Wr.v.y);
    pm = HALF*(pm + Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1));

    if (pm/Wl.p <= ONE) {
      ql = ONE;
    } else {
      ql = sqrt(ONE+HALF*((Wl.g+ONE)/Wl.g)*(pm/Wl.p));
    } /* endif */
    wavespeed_l = Wl.v.y - ql*al;

    if (pm/Wr.p <= ONE) {
      qr = ONE;
    } else {
      qr = sqrt(ONE+HALF*((Wr.g+ONE)/Wr.g)*(pm/Wr.p));
    } /* endif */
    wavespeed_r = Wr.v.y + qr*ar;

    wavespeed_m = vm;

    /* Determine the intermediate state flux. */

    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else  if (wavespeed_m >= ZERO) {
        Uml = Ion5Moment2D_cState(Wl.d,  Wl.d*Wl.v.x, Wl.d*wavespeed_m,
          		     Wl.E()+Wl.d*(wavespeed_m-Wl.v.y)*(wavespeed_m+
        		     Wl.p/(Wl.d*(wavespeed_l-Wl.v.y))))*
              ((wavespeed_l-Wl.v.y)/(wavespeed_l-wavespeed_m));
        Flux = Wl.F()+wavespeed_l*(Uml-Wl.U());
    } else {
        Umr = Ion5Moment2D_cState(Wr.d,  Wr.d*Wr.v.x, Wr.d*wavespeed_m,
          		     Wr.E()+Wr.d*(wavespeed_m-Wr.v.y)*(wavespeed_m+
        		     Wr.p/(Wr.d*(wavespeed_r-Wr.v.y))))*
              ((wavespeed_r-Wr.v.y)/(wavespeed_r-wavespeed_m));
        Flux = Wr.F()+wavespeed_r*(Umr-Wr.U());
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Ion5Moment2D_cState FluxHLLC_y(const Ion5Moment2D_cState &Ul,
	      	               const Ion5Moment2D_cState &Ur) {
   return (FluxHLLC_y(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLC_n (Tito Toro's flux function,       *
 *                      n-direction)                     *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then by using the HLLC              *
 * approximation to specify the intermediate state in    *
 * terms of the rotated solution states.  See Toro,      *
 * Spruce, and Speares (1994).                           *
 *                                                       *
 *********************************************************/
Ion5Moment2D_cState FluxHLLC_n(const Ion5Moment2D_pState &Wl,
	      	               const Ion5Moment2D_pState &Wr,
                               const Vector2D &norm_dir) {

    double cos_angle, sin_angle;
    Ion5Moment2D_pState Wl_rotated, Wr_rotated;
    Ion5Moment2D_cState Flux, Flux_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */

    Wl_rotated.d = Wl.d;
    Wl_rotated.v.x = Wl.v.x*cos_angle +
                     Wl.v.y*sin_angle;
    Wl_rotated.v.y = - Wl.v.x*sin_angle +
                       Wl.v.y*cos_angle;
    Wl_rotated.p = Wl.p;

    Wr_rotated.d = Wr.d;
    Wr_rotated.v.x = Wr.v.x*cos_angle +
                     Wr.v.y*sin_angle;
    Wr_rotated.v.y = - Wr.v.x*sin_angle +
                       Wr.v.y*cos_angle;
    Wr_rotated.p = Wr.p;

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

    Flux_rotated = FluxHLLC_x(Wl_rotated, Wr_rotated);

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */

    Flux.d = Flux_rotated.d;
    Flux.dv.x = Flux_rotated.dv.x*cos_angle -
                Flux_rotated.dv.y*sin_angle;
    Flux.dv.y = Flux_rotated.dv.x*sin_angle +
                Flux_rotated.dv.y*cos_angle;
    Flux.E = Flux_rotated.E;

    return (Flux);

}

Ion5Moment2D_cState FluxHLLC_n(const Ion5Moment2D_cState &Ul,
	      	               const Ion5Moment2D_cState &Ur,
                               const Vector2D &norm_dir) {
    return (FluxHLLC_n(Ul.W(), Ur.W(), norm_dir));
}
