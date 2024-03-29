/**
 *
 *  @file LatticePhase.cpp
 *
 * $Id: LatticePhase.cpp,v 1.3 2007/06/05 15:17:14 hkmoffa Exp $
 */

#include "config.h"
#ifdef WITH_LATTICE_SOLID

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "mix_defs.h"
#include "LatticePhase.h"
#include "SpeciesThermo.h"
#include <math.h>

namespace Cantera {

  //! Base Empty constructor
  LatticePhase::LatticePhase() :
    m_tlast(0.0) 
  {
  }
  
  //! Copy Constructor
  /*!
   * @param right Object to be copied
   */
  LatticePhase::LatticePhase(const LatticePhase &right) :
    m_tlast(0.0)
  {
    *this = operator=(right);
  }
  
  //! Assignment operator
  /*!
   * @param right Object to be copied
   */
  LatticePhase& LatticePhase::operator=(const LatticePhase& right) {
    if (&right != this) {
      ThermoPhase::operator=(right);
      m_mm         = right.m_mm;
      m_tmin       = right.m_tmin;
      m_tmax       = right.m_tmax;
      m_p0         = right.m_p0;
      m_tlast      = right.m_tlast;
      m_h0_RT      = right.m_h0_RT;
      m_cp0_R      = right.m_cp0_R;
      m_g0_RT      = right.m_g0_RT;
      m_s0_R       = right.m_s0_R;
      m_press      = right.m_press;
      m_vacancy    = right.m_vacancy;
      m_molar_density = right.m_molar_density;
    }
    return *this;
  }
  
  //! Destructor
  LatticePhase::~LatticePhase() {
  }
  
  //! Duplication function
  /*!
   * This virtual function is used to create a duplicate of the
   * current phase. It's used to duplicate the phase when given
   * a ThermoPhase pointer to the phase.
   *
   * @return It returns a ThermoPhase pointer.
   */
  ThermoPhase *LatticePhase::duplMyselfAsThermoPhase() const {
    LatticePhase *igp = new LatticePhase(*this);
    return (ThermoPhase *) igp;
  }
  
  
    doublereal LatticePhase::
    enthalpy_mole() const {
        doublereal p0 = m_spthermo->refPressure();
        return GasConstant * temperature() * 
            mean_X(&enthalpy_RT()[0]) 
            + (pressure() - p0)/molarDensity();
    }

    doublereal LatticePhase::intEnergy_mole() const {
        doublereal p0 = m_spthermo->refPressure();
        return GasConstant * temperature() * 
            mean_X(&enthalpy_RT()[0]) 
            - p0/molarDensity();
    }

    doublereal LatticePhase::entropy_mole() const {
        return GasConstant * (mean_X(&entropy_R()[0]) -
            sum_xlogx());
    }

    doublereal LatticePhase::gibbs_mole() const {
        return enthalpy_mole() - temperature() * entropy_mole();
    }

    doublereal LatticePhase::cp_mole() const {
        return GasConstant * mean_X(&cp_R()[0]);
    }

    void LatticePhase::getActivityConcentrations(doublereal* c) const {
        getMoleFractions(c);
    }

    void LatticePhase::getActivityCoefficients(doublereal* ac) const {
        for (int k = 0; k < m_kk; k++) {
	  ac[k] = 1.0;
	}
    }

    doublereal LatticePhase::standardConcentration(int k) const {
        return 1.0;
    }

    doublereal LatticePhase::logStandardConc(int k) const {
        return 0.0;
    }

    void LatticePhase::getChemPotentials(doublereal* mu) const {
        doublereal vdp = (pressure() - m_spthermo->refPressure())/
                         molarDensity();
        doublereal xx;
        doublereal rt = temperature() * GasConstant;
        const array_fp& g_RT = gibbs_RT();
        for (int k = 0; k < m_kk; k++) {
            xx = fmaxx(SmallNumber, moleFraction(k));
            mu[k] = rt*(g_RT[k] + log(xx)) + vdp;
        }
    }

    void LatticePhase::getStandardChemPotentials(doublereal* mu0) const {
        const array_fp& gibbsrt = gibbs_RT();
        scale(gibbsrt.begin(), gibbsrt.end(), mu0, _RT());
    }

    void LatticePhase::initThermo() {
        m_kk = nSpecies();
        m_mm = nElements();
        doublereal tmin = m_spthermo->minTemp();
        doublereal tmax = m_spthermo->maxTemp();
        if (tmin > 0.0) m_tmin = tmin;
        if (tmax > 0.0) m_tmax = tmax;
        m_p0 = refPressure();

        int leng = m_kk;
        m_h0_RT.resize(leng);
        m_g0_RT.resize(leng);
        m_cp0_R.resize(leng);
        m_s0_R.resize(leng);
        setMolarDensity(m_molar_density);
    }


    void LatticePhase::_updateThermo() const {
        doublereal tnow = temperature();
        if (fabs(molarDensity() - m_molar_density)/m_molar_density > 0.0001) {
            throw CanteraError("_updateThermo","molar density changed from "
                +fp2str(m_molar_density)+" to "+fp2str(molarDensity()));
        }
        if (m_tlast != tnow) {
            m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0], 
                &m_s0_R[0]);
            m_tlast = tnow;
            int k;
            for (k = 0; k < m_kk; k++) {
                m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
            }
            m_tlast = tnow;
        }
    }

    void LatticePhase::setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","Lattice");
        m_molar_density = getFloat(eosdata, "site_density", "-");
        m_vacancy = getString(eosdata, "vacancy_species");
    }
}

#endif
