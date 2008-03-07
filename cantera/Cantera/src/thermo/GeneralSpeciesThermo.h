/**
 * @file GeneralSpeciesThermo.h
 *  Headers for a completely general species thermodynamic property
 *  manager for a phase (see \ref spthermo and
 * \link Cantera::GeneralSpeciesThermo GeneralSpeciesThermo\endlink).
 *
 *  Because it is general, it is slow.
 */

/*
 * $Author: hkmoffa $
 * $Revision: 1.3 $
 * $Date: 2007/09/13 15:05:39 $
 */

#ifndef CT_GENERALSPECIESTHERMO_H
#define CT_GENERALSPECIESTHERMO_H
#include <string>
#include "ct_defs.h"
#include "SpeciesThermoMgr.h"
#include "NasaPoly1.h"
#include "Nasa9Poly1.h"
#include "speciesThermoTypes.h"
//#include "polyfit.h"

namespace Cantera {

  
  //! A species thermodynamic property manager for a phase.
  /*!
   * This is a general manager that can handle a wide variety
   * of species thermodynamic polynomials for individual species.
   * It is slow, however, because it recomputes the functions of
   * temperature needed for each species. What it does is to create
   * a vector of SpeciesThermoInterpType objects.
   *
   * @ingroup spthermo
   */
  class GeneralSpeciesThermo : public SpeciesThermo {
    
  public:

    //! Constructor
    GeneralSpeciesThermo();

    //! Copy constructor
    GeneralSpeciesThermo(const GeneralSpeciesThermo &);

    //! Assignment operator
    GeneralSpeciesThermo & operator=(const GeneralSpeciesThermo &);

    //! destructor
    virtual ~GeneralSpeciesThermo();

    //! Duplicator
    virtual SpeciesThermo *duplMyselfAsSpeciesThermo() const ;

    //! Install a new species thermodynamic property
    //! parameterization for one species.  
    /*!
     * Install a SpeciesThermoInterpType object for the species, index.
     * This routine contains an internal list of  SpeciesThermoInterpType
     * objects that it knows about. A factory-type lookup is done
     * to create the object.
     *
     * @param name      Name of the species
     * @param index     The 'update' method will update the property 
     *                  values for this species 
     *                  at position i index in the property arrays.  
     * @param type      int flag specifying the type of parameterization to be
     *                 installed. 
     * @param c        vector of coefficients for the parameterization. 
     *                 This vector is simply passed through to the
     *                 parameterization constructor. It's length depends upon
     *                 the parameterization.
     * @param minTemp  minimum temperature for which this parameterization
     *                 is valid.
     * @param maxTemp  maximum temperature for which this parameterization
     *                 is valid.
     * @param refPressure standard-state pressure for this 
     *                    parameterization. 
     * @see speciesThermoTypes.h 
     *
     * @todo Create a factory method for SpeciesThermoInterpType.
     *       That's basically what we are doing here.
     */
    virtual void install(std::string name, int index, int type, 
			 const doublereal* c, 
			 doublereal minTemp, doublereal maxTemp,
			 doublereal refPressure);

    //! Install a new species thermodynamic property
    //! parameterization for one species.
    /*!
     * @param stit_ptr Pointer to the SpeciesThermoInterpType object
     *          This will set up the thermo for one species
     */
    virtual void install_STIT(SpeciesThermoInterpType *stit_ptr);

    //! Like update(), but only updates the single species k.
    /*!
     * @param k       species index
     * @param T       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    virtual void update_one(int k, doublereal T, doublereal* cp_R, 
			    doublereal* h_RT,
			    doublereal* s_R) const;

    //! Compute the reference-state properties for all species.
    /*!
     * Given temperature T in K, this method updates the values of
     * the non-dimensional heat capacity at constant pressure,
     * enthalpy, and entropy, at the reference pressure, Pref
     * of each of the standard states.
     *
     * @param T       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */	
    virtual void update(doublereal T, doublereal* cp_R, 
			doublereal* h_RT, doublereal* s_R) const;
                
    //! Minimum temperature.
    /*!
     * If no argument is supplied, this
     * method returns the minimum temperature for which \e all
     * parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the minimum
     * temperature for species k in the phase.
     *
     * @param k    Species index
     */ 
    virtual doublereal minTemp(int k=-1) const;

    //! Maximum temperature.
    /*!
     * If no argument is supplied, this
     * method returns the maximum temperature for which \e all
     * parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the maximum
     * temperature for parameterization k.
     *
     * @param k  Species Index
     */
    virtual doublereal maxTemp(int k=-1) const;

    //! The reference-state pressure for species k.
    /*!
     *
     * returns the reference state pressure in Pascals for
     * species k. If k is left out of the argument list,
     * it returns the reference state pressure for the first
     * species.
     * Note that some SpeciesThermo implementations, such
     * as those for ideal gases, require that all species
     * in the same phase have the same reference state pressures.
     *
     * @param k Species Index
     */
    virtual doublereal refPressure(int k = -1) const;

    //! This utility function reports the type of parameterization
    //! used for the species with index number index.
    /*!
     *
     * @param index  Species index
     */
    virtual int reportType(int index) const;

    //! This utility function reports back the type of 
    //! parameterization and all of the parameters for the species, index.
    /*!
     * @param index     Species index
     * @param type      Integer type of the standard type
     * @param c         Vector of coefficients used to set the
     *                  parameters for the standard state.
     * @param minTemp   output - Minimum temperature
     * @param maxTemp   output - Maximum temperature
     * @param refPressure output - reference pressure (Pa).
     */
    virtual void reportParams(int index, int &type, 
			      doublereal * const c, 
			      doublereal &minTemp, 
			      doublereal &maxTemp,
			      doublereal &refPressure) const;

    //! Modify parameters for the standard state
    /*!
     * @param index Species index
     * @param c     Vector of coefficients used to set the
     *              parameters for the standard state.
     */
    virtual void modifyParams(int index, doublereal *c);

  protected:

    /**
     * This is the main unknown in the object. It is 
     * a list of pointers to type SpeciesThermoInterpType.
     * Note, this object owns the objects, so they are deleted
     * in the destructor of this object.
     */
    std::vector<SpeciesThermoInterpType *>  m_sp;

    //! Maximum value of the lowest temperature
    doublereal                         m_tlow_max;

    //! Minimum value of the highest temperature
    doublereal                         m_thigh_min;

    //! reference pressure (Pa)
    doublereal                         m_p0;

    /**
     * Internal variable indicating the length of the 
     * number of species in the phase.
     */
    int m_kk;

  private:

  

  };

}

#endif

