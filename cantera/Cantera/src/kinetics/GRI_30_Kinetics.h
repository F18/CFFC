/**
 *
 * @file GRI_30_Kinetics.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_GRI30_KINETICS_H
#define CT_GRI30_KINETICS_H

#include "GasKinetics.h"

namespace Cantera {

    const int cGRI_30_Kinetics = cGasKinetics + 1;

    /**
     *  Kinetics manager implementing reaction mechanism GRI-Mech 3.0
     */    
    class GRI_30_Kinetics : public GasKinetics {

    public:

        /// Default constructor.
        GRI_30_Kinetics(thermo_t* th=0);

        /// Destructor.
        virtual ~GRI_30_Kinetics(){}

        virtual int ID() { return cGRI_30_Kinetics; }

        virtual void getNetProductionRates(doublereal* net) {
            gri30_updateROP();
            get_wdot(&m_kdata->m_ropnet[0], net);
        }

    private:
        void gri30_update_rates_T();
        void gri30_updateROP();
        void gri30_updateKc();
        void get_wdot(const doublereal* rop, doublereal* wdot);
        void update_kc(const double* grt, double c0, double* rkc);
        void update_rates(double t, double tlog, double* rf);
        void eval_ropnet(const double* c, const double* rf, const double* rkc, double* r);
    };
}

#endif
