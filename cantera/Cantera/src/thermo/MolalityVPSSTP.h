/**
 *  @file MolalityVPSSTP.h
 *   Header for intermediate ThermoPhase object for phases which
 *   employ molality based activity coefficient formulations
 *  (see \ref thermoprops 
 * and class \link Cantera::MolalityVPSSTP MolalityVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon activities
 * based on the molality scale.  These include most of the methods for
 * calculating liquid electrolyte thermodynamics.
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Id: MolalityVPSSTP.h,v 1.14 2007/06/13 00:08:35 hkmoffa Exp $
 */

#ifndef CT_MOLALITYVPSSTP_H
#define CT_MOLALITYVPSSTP_H

#include "VPStandardStateTP.h"

namespace Cantera {

  /**
   * @ingroup thermoprops
   */

  /*!
   * MolalityVPSSTP is a derived class of ThermoPhase that handles
   * variable pressure standard state methods for calculating
   * thermodynamic properties that are further based on
   * molality-scaled activities.
   * This category incorporates most of the methods
   * for calculating liquid electrolyte thermodynamics that have been
   * developed since the 1970's.
   *
   * This class adds additional functions onto the %ThermoPhase interface
   * that handle molality based standard states. The %ThermoPhase
   * class includes a member function, ThermoPhase::activityConvention() 
   * that indicates which convention the activities are based on. The
   * default is to assume activities are based on the molar convention. 
   * However, classes which derive from the MolalityVPSSTP class return
   * <b>cAC_CONVENTION_MOLALITY</b> from this member function.
   *
   * The molality of a solute, \f$ m_i \f$, is defined as
   *
   * \f[
   *     m_i = \frac{n_i}{\tilde{M}_o n_o} 
   * \f]
   * where
   * \f[
   *     \tilde{M}_o = \frac{M_o}{1000} 
   * \f]
   *    
   * where \f$ M_o \f$ is the molecular weight of the solvent. The molality
   * has units of gmol kg<SUP>-1</SUP>. For the solute, the molality may be
   * considered as the amount of gmol's of solute per kg of solvent, a natural
   * experimental quantity.
   *
   * The formulas for calculating mole fractions if given the molalities of
   * the solutes is stated below. First calculate \f$ L^{sum} \f$, an intermediate
   * quantity.
   *
   *   \f[
   *       L^{sum} = \frac{1}{\tilde{M}_o X_o} = \frac{1}{\tilde{M}_o} + \sum_{i\ne o} m_i
   *  \f]
   *  Then,
   *  \f[
   *       X_o =   \frac{1}{\tilde{M}_o L^{sum}}  
   *  \f]
   *  \f[
   *       X_i =   \frac{m_i}{L^{sum}}  
   *  \f]
   * where \f$ X_o \f$ is the mole fraction of solvent, and \f$ X_o \f$ is the
   * mole fraction of solute <I>i</I>. Thus, the molality scale and the mole fraction
   * scale offer a one-to-one mapping between each other, except in the limit
   * of a zero solvent mole fraction.
   *
   * The standard states for thermodynamic objects that derive from <b>MolalityVPSSTP</b>
   * are on the unit molality basis. Chemical potentials
   * of the solutes,  \f$ \mu_k \f$, and the solvent, \f$ \mu_o \f$, which are based 
   * on the molality form, have the following general format:
   *
   * \f[
   *    \mu_k = \mu^{\triangle}_k(T,P) + R T ln(\gamma_k^{\triangle} \frac{m_k}{m^\triangle}) 
   * \f]
   * \f[
   *    \mu_o = \mu^o_o(T,P) + RT ln(a_o) 
   * \f]
   *
   * where \f$ \gamma_k^{\triangle} \f$ is the molality based activity coefficient for species
   * \f$k\f$. 
   *
   * The chemical potential of the solvent is thus expressed in a different format
   * than the chemical potential of the solutes. Additionally, the activity of the
   * solvent, \f$ a_o \f$, is further reexpressed in terms of an osmotic coefficient,
   * \f$ \phi \f$.
   *   \f[
   *       \phi = \frac{- ln(a_o)}{\tilde{M}_o \sum_{i \ne o} m_i}
   *   \f]
   * 
   *  MolalityVPSSTP::osmoticCoefficient() returns the value of \f$ \phi \f$. 
   *  Note there  are a few of definitions of the osmotic coefficient floating
   *  around. We use the one defined in 
   *  (Activity Coefficients in Electrolyte Solutions, K. S. Pitzer
   *   CRC Press, Boca Raton, 1991, p. 85, Eqn. 28). This definition is most clearly
   *  related to theoretical calculation.
   *
   * The molar-based activity coefficients \f$ \gamma_k \f$ may be calculated 
   * from the molality-based
   * activity coefficients, \f$ \gamma_k^\triangle \f$ by the following 
   * formula.
   * \f[
   *     \gamma_k = \frac{\gamma_k^\triangle}{X_o}
   * \f]
   * For purposes of establishing a convention, the molar activity coefficient of the
   * solvent is set equal to the molality-based activity coefficient of the
   * solvent: 
   * \f[
   *     \gamma_o = \gamma_o^\triangle
   * \f]
   *
   *  The molality-based and molarity-based standard states may be related to one
   *  another by the following formula.
   *
   * \f[
   *   \mu_k^\triangle(T,P) = \mu_k^o(T,P) + R T \ln(\tilde{M}_o m^\triangle)
   * \f]
   *
   *  An important convention is followed in all routines that derive from <b>%MolalityVPSSTP</b>.
   *  Standard state thermodynamic functions and reference state thermodynamic functions
   *  return the molality-based quantities. Also all functions which return 
   *  activities return the molality-based activities. The reason for this convention
   *  has been discussed in supporting memos. However, it's important because the
   *  term in the equation above is non-trivial. For example it's equal to 2.38 kcal gmol<SUP>-1</SUP>
   *  for water at 298 K. 
   *
   *
   * In order to prevent a singularity, this class includes the concept of a minimum
   * value for the solvent mole fraction. All calculations involving the formulation
   * of activity coefficients and other non-ideal solution behavior adhere to 
   * this concept of a minimul value for the solvent mole fraction. This makes sense
   * because these solution behavior were all designed and measured far away from 
   * the zero solvent singularity condition and are not applicable in that limit. 
   *   
   *
   * This objects add a layer that supports molality. It inherits from VPStandardStateTP.
   *
   *  All objects that derive from this are assumed to have molality based standard states.
   *
   * @todo Make two solvent minimum fractions. One would be for calculation of the non-ideal
   *       factors. The other one would be for purposes of stoichiometry evaluation. the
   *       stoichiometry evaluation one would be a 1E-13 limit. Anything less would create
   *       problems with roundoff error.
   *
   */
  class MolalityVPSSTP : public VPStandardStateTP  {

  public:
        
    /// Constructors 
    /*!
     * This doesn't do much more than initialize constants with
     * default values for water at 25C. Water molecular weight 
     * comes from the default elements.xml file. It actually
     * differs slightly from the IAPWS95 value of 18.015268. However,
     * density conservation and therefore element conservation
     * is the more important principle to follow.
     */
    MolalityVPSSTP();

    //! Copy constructor
    /*!
     *  Note this stuff will not work until the underlying phase
     *  has a working copy constructor
     *
     * @param b class to be copied
     */
    MolalityVPSSTP(const MolalityVPSSTP &b);

    /// Assignment operator
    /*!
     *  Note this stuff will not work until the underlying phase
     *  has a working assignment operator
     *
     * @param b class to be copied.
     */
    MolalityVPSSTP& operator=(const	MolalityVPSSTP&b);

    /// Destructor. 
    virtual ~MolalityVPSSTP();

    //! Duplication routine for objects which inherit from  ThermoPhase.
    /*!
     *  This virtual routine can be used to duplicate thermophase objects
     *  inherited from ThermoPhase even if the application only has
     *  a pointer to ThermoPhase to work with.
     */
    virtual ThermoPhase *duplMyselfAsThermoPhase() const;
    
    /**
     *   
     * @name  Utilities  
     * @{
     */

   
    //! Equation of state type flag.
    /*!
     * The ThermoPhase base class returns
     * zero. Subclasses should define this to return a unique
     * non-zero value. Known constants defined for this purpose are
     * listed in mix_defs.h. The MolalityVPSSTP class also returns
     * zero, as it is a non-complete class.
     */
    virtual int eosType() const { return 0; }


    /**
     * @} 
     * @name  Molar Thermodynamic Properties 
     * @{
     */


    /**
     * @}
     * @name Utilities for Solvent ID and Molality
     * @{
     */

    /**
     * This routine sets the index number of the solvent for
     * the phase.
     *  
     *  Note, having a solvent 
     *   is a precursor to many things having to do with molality.
     *
     * @param k the solvent index number
     */
    void setSolvent(int k);

    /**
     * Sets the minimum mole fraction in the molality formulation.
     * Note the molality formulation is singular in the limit that
     * the solvent mole fraction goes to zero. Numerically, how
     * this limit is treated and resolved is an ongoing issue within
     * Cantera.
     *
     * @param xmolSolventMIN  Input double containing the minimum mole fraction
     */
    void setMoleFSolventMin(doublereal xmolSolventMIN);

    //! Returns the solvent index.
    int solventIndex() const;

    /**
     * Returns the minimum mole fraction in the molality 
     * formulation.
     */
    doublereal moleFSolventMin() const;

    //! Calculates the molality of all species and stores the result internally.
    /*!
     *   We calculate the vector of molalities of the species
     *   in the phase and store the result internally:
     *   \f[
     *     m_i = \frac{X_i}{1000 * M_o * X_{o,p}}
     *   \f]
     *    where 
     *    - \f$ M_o \f$ is the molecular weight of the solvent
     *    - \f$ X_o \f$ is the mole fraction of the solvent
     *    - \f$ X_i \f$ is the mole fraction of the solute.
     *    - \f$ X_{o,p} = max (X_{o}^{min}, X_o) \f$
     *    - \f$ X_{o}^{min} \f$ = minimum mole fraction of solvent allowed
     *              in the denominator.
     */
    void calcMolalities() const;

    //!  This function will return the molalities of the species.
    /*!
     *   We calculate the vector of molalities of the species
     *   in the phase
     * \f[
     *     m_i = \frac{X_i}{1000 * M_o * X_{o,p}}
     * \f]
     *    where 
     *    - \f$ M_o \f$ is the molecular weight of the solvent
     *    - \f$ X_o \f$ is the mole fraction of the solvent
     *    - \f$ X_i \f$ is the mole fraction of the solute.
     *    - \f$ X_{o,p} = \max (X_{o}^{min}, X_o) \f$
     *    - \f$ X_{o}^{min} \f$ = minimum mole fraction of solvent allowed
     *              in the denominator.
     *
     * @param molal       Output vector of molalities. Length: m_kk.
     */
    void getMolalities(doublereal * const molal) const;

    //! Set the molalities of the solutes in a phase
    /*!
     * Note, the entry for the solvent is not used.
     *   We are supplied with the molalities of all of the
     *   solute species. We then calculate the mole fractions of all
     *   species and update the %ThermoPhase object.
     *  \f[
     *     m_i = \frac{X_i}{M_o/1000 * X_{o,p}}
     *  \f]
     *    where 
     *    -  \f$M_o\f$ is the molecular weight of the solvent
     *    -  \f$X_o\f$ is the mole fraction of the solvent
     *    -  \f$X_i\f$ is the mole fraction of the solute.
     *    -  \f$X_{o,p} = \max(X_o^{min}, X_o)\f$
     *    -  \f$X_o^{min}\f$ = minimum mole fraction of solvent allowed
     *                     in the denominator.
     *
     * The formulas for calculating mole fractions are
     *  \f[
     *       L^{sum} = \frac{1}{\tilde{M}_o X_o} = \frac{1}{\tilde{M}_o} + \sum_{i\ne o} m_i
     *  \f]
     *  Then,
     *  \f[
     *       X_o =   \frac{1}{\tilde{M}_o L^{sum}}  
     *  \f]
     *  \f[
     *       X_i =   \frac{m_i}{L^{sum}}  
     *  \f]
     *  It is currently an error if the solvent mole fraction is attempted to be set 
     *  to a value lower than  \f$X_o^{min}\f$.
     *
     * @param molal   Input vector of molalities. Length: m_kk.
     */
    void setMolalities(const doublereal * const molal);

    //! Set the molalities of a phase
    /*!
     * Set the molalities of the solutes in a phase. Note, the entry for the
     * solvent is not used.
     *
     * @param xMap  Composition Map containing the molalities.
     */
    void setMolalitiesByName(compositionMap& xMap);

    //! Set the molalities of a phase
    /*!
     * Set the molalities of the solutes in a phase. Note, the entry for the
     * solvent is not used.
     *
     * @param name  String containing the information for a composition map.
     */
    void setMolalitiesByName(const std::string& name);

    /**
     * @}
     * @name Mechanical Properties
     * @{
     */

    /**
     * @} 
     * @name Potential Energy
     * 
     * Species may have an additional potential energy due to the
     * presence of external gravitation or electric fields. These
     * methods allow specifying a potential energy for individual
     * species.
     * @{
     */

    /**
     * @}
     * @name Activities, Standard States, and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by \f[ \mu_k = \mu_k^0(T)
     * + \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T,P)\f$ is
     * the chemical potential at unit activity, which depends only
     * on temperature and pressure.
     * @{
     */

    /**
     * This method returns the activity convention.
     * Currently, there are two activity conventions
     *  Molar-based activities
     *       Unit activity of species at either a hypothetical pure
     *       solution of the species or at a hypothetical
     *       pure ideal solution at infinite dilution
     *   cAC_CONVENTION_MOLAR 0
     *      - default
     *  
     *  Molality based acvtivities
     *       (unit activity of solutes at a hypothetical 1 molal
     *        solution referenced to infinite dilution at all
     *        pressures and temperatures).
     *   cAC_CONVENTION_MOLALITY 1
     *
     *  We set the convention to molality here.
     */
    int activityConvention() const;

    /**
     * This method returns an array of generalized concentrations
     * \f$ C_k\f$ that are defined such that 
     * \f$ a_k = C_k / C^0_k, \f$ where \f$ C^0_k \f$ 
     * is a standard concentration
     * defined below.  These generalized concentrations are used
     * by kinetics manager classes to compute the forward and
     * reverse rates of elementary reactions. 
     *
     * @param c Array of generalized concentrations. The 
     *           units depend upon the implementation of the
     *           reaction rate expressions within the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const {
      err("getActivityConcentrations");
    }

    /**
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the generalized concentration. In many cases, this quantity
     * will be the same for all species in a phase - for example,
     * for an ideal gas \f$ C^0_k = P/\hat R T \f$. For this
     * reason, this method returns a single value, instead of an
     * array.  However, for phases in which the standard
     * concentration is species-specific (e.g. surface species of
     * different sizes), this method may be called with an
     * optional parameter indicating the species.
     *
     * @param k species index. Defaults to zero.
     */
    virtual doublereal standardConcentration(int k=0) const {
      err("standardConcentration");
      return -1.0;
    }

    /**
     * Returns the natural logarithm of the standard 
     * concentration of the kth species
     *
     * @param k  species index
     */
    virtual doublereal logStandardConc(int k=0) const {
      err("logStandardConc");
      return -1.0;
    }

    /**
     * Returns the units of the standard and generalized
     * concentrations Note they have the same units, as their
     * ratio is defined to be equal to the activity of the kth
     * species in the solution, which is unitless.
     *
     * This routine is used in print out applications where the
     * units are needed. Usually, MKS units are assumed throughout
     * the program and in the XML input files.
     *
     * @param uA Output vector containing the units
     *  uA[0] = kmol units - default  = 1
     *  uA[1] = m    units - default  = -nDim(), the number of spatial
     *                                dimensions in the Phase class.
     *  uA[2] = kg   units - default  = 0;
     *  uA[3] = Pa(pressure) units - default = 0;
     *  uA[4] = Temperature units - default = 0;
     *  uA[5] = time units - default = 0
     * @param k species index. Defaults to 0.
     * @param sizeUA output int containing the size of the vector.
     *        Currently, this is equal to 6.
     */
    virtual void getUnitsStandardConc(double *uA, int k = 0,
				      int sizeUA = 6) const;

    
    //! Get the array of non-dimensional activities (molality
    //! based for this class and classes that derive from it) at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * All standard state properties for molality-based phases are
     * evaluated consistent with the molality scale. Therefore, this function
     * must return molality-based activities.
     *
     * \f[
     *  a_i^\triangle = \gamma_k^{\triangle} \frac{m_k}{m^\triangle}
     * \f]
     *
     * This function must be implemented in derived classes.
     *
     * @param ac     Output vector of molality-based activities. Length: m_kk.
     */
    virtual void getActivities(doublereal* ac) const {
      err("getActivities");
    }

    //! Get the array of non-dimensional activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * These are mole-fraction based activity coefficients. In this
     * object, their calculation is based on translating the values
     * of the molality-based activity coefficients.
     *  See Denbigh p. 278 for a thorough discussion.
     *
     * The molar-based activity coefficients \f$ \gamma_k \f$ may be calculated from the 
     * molality-based
     * activity coefficients, \f$ \gamma_k^\triangle \f$ by the following 
     * formula.
     * \f[
     *     \gamma_k = \frac{\gamma_k^\triangle}{X_o}
     * \f]
     *
     * For purposes of establishing a convention, the molar activity coefficient of the
     * solvent is set equal to the molality-based activity coefficient of the
     * solvent:
     * 
     * \f[
     *     \gamma_o = \gamma_o^\triangle
     * \f]
     *
     * Derived classes don't need to overload this function. This function is
     * handled at this level.
     *
     * @param ac  Output vector containing the mole-fraction based activity coefficients.
     *            length: m_kk.
     */
    void getActivityCoefficients(doublereal* ac) const;

    //! Get the array of non-dimensional molality based 
    //!  activity coefficients at the current solution temperature, 
    //!  pressure, and  solution concentration.
    /*!
     *  See Denbigh p. 278 for a thorough discussion. This class must be overwritten in
     *  classes which derive from %MolalityVPSSTP. This function takes over from the
     *  molar-based activity coefficient calculation, getActivityCoefficients(), in
     *  derived classes.
     *
     * @param acMolality Output vector containing the molality based activity coefficients.
     *                   length: m_kk.
     */
    virtual void getMolalityActivityCoefficients(doublereal *acMolality) const {
      err("getMolalityActivityCoefficients");
    }
   
    //! Calculate the osmotic coefficient
    /*!
     *   \f[
     *       \phi = \frac{- ln(a_o)}{\tilde{M}_o \sum_{i \ne o} m_i}
     *   \f]
     * 
     *  Note there  are a few of definitions of the osmotic coefficient floating
     *  around. We use the one defined in 
     *  (Activity Coefficients in Electrolyte Solutions, K. S. Pitzer
     *   CRC Press, Boca Raton, 1991, p. 85, Eqn. 28). This definition is most clearly
     *  related to theoretical calculation.
     *
     *  units = dimensionless
     */
    virtual double osmoticCoefficient() const;

    //@}
    /// @name  Partial Molar Properties of the Solution 
    //@{


    /**
     * Get the species electrochemical potentials. 
     * These are partial molar quantities.
     * This method adds a term \f$ Fz_k \phi_k \f$ to the 
     * to each chemical potential.
     *
     * Units: J/kmol
     *
     * @param mu     output vector containing the species electrochemical potentials.
     *               Length: m_kk.
     */
    void getElectrochemPotentials(doublereal* mu) const {
      getChemPotentials(mu);
      double ve = Faraday * electricPotential();
      for (int k = 0; k < m_kk; k++) {
	mu[k] += ve*charge(k);
      }
    }

 
    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    //@{

     

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{


    ///////////////////////////////////////////////////////
    //
    //  The methods below are not virtual, and should not
    //  be overloaded.
    //
    //////////////////////////////////////////////////////

    /**
     * @name Specific Properties
     * @{
     */


    /**
     * @name Setting the State
     *
     * These methods set all or part of the thermodynamic
     * state.
     * @{
     */

    //@}

    /**
     * @name Chemical Equilibrium
     * Routines that implement the Chemical equilibrium capability
     * for a single phase, based on the element-potential method.
     * @{
     */

    /**
     * This method is used by the ChemEquil element-potential
     * based equilibrium solver.
     * It sets the state such that the chemical potentials of the 
     * species within the current phase satisfy
     * \f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
     * \left(\frac{\lambda_m} {\hat R T}\right) \f] where 
     * \f$ \lambda_m \f$ is the element potential of element m. The
     * temperature is unchanged.  Any phase (ideal or not) that
     * implements this method can be equilibrated by ChemEquil.
     *
     * @param lambda_RT Input vector containing the dimensionless 
     *                  element potentials.
     */
    virtual void setToEquilState(const doublereal* lambda_RT) {
      err("setToEquilState");
    }


    //@}


    //! Set equation of state parameter values from XML entries.
    /*!
     * This method is called by function importPhase() in
     * file importCTML.cpp when processing a phase definition in
     * an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. 
     *
     * The MolalityVPSSTP object defines a new method for setting
     * the concentrations of a phase. The new method is defined by a
     * block called "soluteMolalities". If this block
     * is found, the concentrations within that phase are
     * set to the "name":"molalities pairs found within that
     * XML block. The solvent concentration is then set
     * to everything else.
     *
     * The function first calls the overloaded function ,
     * VPStandardStateTP::setStateFromXML(), to pick up the parent class
     * behavior.
     *
     * usage: Overloaded functions should call this function
     *        before carrying out their own behavior.
     *   
     * @param state An XML_Node object corresponding to
     *              the "state" entry for this phase in the input file.
     */
    virtual void setStateFromXML(const XML_Node& state);

    /// The following methods are used in the process of constructing
    /// the phase and setting its parameters from a specification in an 
    /// input file. They are not normally used in application programs.
    /// To see how they are used, see files importCTML.cpp and 
    /// ThermoFactory.cpp.


    /*!
     * @internal Initialize. This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase.
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();


    /**
     *   Import and initialize a ThermoPhase object
     *
     * @param phaseNode This object must be the phase node of a
     *             complete XML tree
     *             description of the phase, including all of the
     *             species data. In other words while "phase" must
     *             point to an XML phase object, it must have
     *             sibling nodes "speciesData" that describe
     *             the species in the phase.
     * @param id   ID of the phase. If nonnull, a check is done
     *             to see if phaseNode is pointing to the phase
     *             with the correct id. 
     */
    void initThermoXML(XML_Node& phaseNode, std::string id);

    
    //! Set the temperature (K), pressure (Pa), and molalities
    //!(gmol kg-1) of the solutes
    /*!
     * @param t           Temperature (K)
     * @param p           Pressure (Pa)
     * @param molalities  Input vector of molalities of the solutes.
     *                    Length: m_kk.
     */
    void setState_TPM(doublereal t, doublereal p, 
		      const doublereal * const molalities);

    //! Set the temperature (K), pressure (Pa), and molalities.
    /*!
     * @param t           Temperature (K)
     * @param p           Pressure (Pa)
     * @param m           compositionMap containing the molalities
     */
    void setState_TPM(doublereal t, doublereal p, compositionMap& m);

    //! Set the temperature (K), pressure (Pa), and molalities. 
    /*!
     * @param t           Temperature (K)
     * @param p           Pressure (Pa)
     * @param m           String which gets translated into a composition
     *                    map for the molalities of the solutes.
     */
    void setState_TPM(doublereal t, doublereal p, const std::string& m);

  private:
    void initLengths();
            
  protected:

    //! Index of the solvent
    /*!
     * Currently the index of the solvent is hard-coded to the value 0
     */
    int        m_indexSolvent;

    //! Molecular weight of the Solvent
    doublereal m_weightSolvent;

    /*!
     * In any molality implementation, it makes sense to have
     * a minimum solvent mole fraction requirement, since the 
     * implementation becomes singular in the xmolSolvent=0
     * limit. The default is to set it to 0.01.
     * We then modify the molality definition to ensure that
     * molal_solvent = 0 when xmol_solvent = 0. 
     */
    doublereal m_xmolSolventMIN;

    //! This is the multiplication factor that goes inside 
    //! log expressions involving the molalities of species.
    /*!
     * Its equal to Wt_0 / 1000,
     *     where Wt_0 = weight of solvent (kg/kmol)
     */
    doublereal m_Mnaught;

    //! Current value of the molalities of the species in the phase.
    /*!
     * Note this vector is a mutable quantity.
     * units are (kg/kmol)
     */
    mutable vector_fp  m_molalities;
  private:
    doublereal err(std::string msg) const;

  };

}
        
#endif





