#include "MultiPhaseEquil.h"
#include "MultiPhase.h"
#include "sort.h"
#include "global.h"

#include <math.h>
#include <iostream>
using namespace std;


namespace Cantera {

    const doublereal TINY = 1.0e-20;

#if defined(WITH_HTML_LOGS)
    /// Used to print reaction equations. Given a stoichiometric
    /// coefficient 'nu' and a chemical symbol 'sym', return a string
    /// for this species in the reaction. 
    /// @param first  if this is false, then a " + " string will be
    /// added to the beginning of the string.
    /// @param nu  Stoichiometric coefficient. May be positive or negative. The
    /// absolute value will be used in the string. 
    /// @param sym Species chemical symbol.
    /// 
    static string coeffString(bool first, doublereal nu, string sym) {
        if (nu == 0.0) return "";
        string strt = " + ";
        if (first) strt = "";
        if (nu == 1.0 || nu == -1.0) 
            return strt + sym;
        string s = fp2str(fabs(nu));
        return strt + s + " " + sym;
    }
#endif

    /// Constructor. Construct a multiphase equilibrium manager for a
    /// multiphase mixture.
    /// @param mix Pointer to a multiphase mixture object.
    /// @param start If true, the initial composition will be
    /// determined by a linear Gibbs minimization, otherwise the
    /// initial mixture composition will be used.
    MultiPhaseEquil::MultiPhaseEquil(mix_t* mix, bool start, int loglevel) : m_mix(mix)
    {
        // the multi-phase mixture
        //        m_mix = mix;

        // store some mixture parameters locally
        m_nel_mix = mix->nElements();
        m_nsp_mix = mix->nSpecies();
        m_np = mix->nPhases();
        m_press = mix->pressure();
        m_temp = mix->temperature();

        index_t m, k;
        m_nel = 0;
        m_nsp = 0;
        m_eloc = 1000;
        m_incl_species.resize(m_nsp_mix,1);
        m_incl_element.resize(m_nel_mix,1);
        for (m = 0; m < m_nel_mix; m++) {
            string enm = mix->elementName(m);
            // element 'E' or 'e' represents an electron; this
            // requires special handling, so save its index
            // for later use
            if (enm == "E" || enm == "e") {
                m_eloc = m;
            }
            // if an element other than electrons is not present in
            // the mixture, then exclude it and all species containing
            // it from the calculation. Electrons are a special case,
            // since a species can have a negative number of 'atoms'
            // of electrons (positive ions).
            if (m_mix->elementMoles(m) <= 0.0) {
                if (m != m_eloc) {
                    m_incl_element[m] = 0;
                    for (k = 0; k < m_nsp_mix; k++) {
                        if (m_mix->nAtoms(k,m) != 0.0) {
                            m_incl_species[k] = 0;
                        }
                    }
                }
            }
        }
        
        // Now build the list of elements to be included, starting with
        // electrons, if they are present.
        if (m_eloc < m_nel_mix) {
            m_element.push_back(m_eloc);
            m_nel++;
        }
        // add the included elements other than electrons
        for (m = 0; m < m_nel_mix; m++) {
            if (m_incl_element[m] == 1 && m != m_eloc) {
                m_nel++;
                m_element.push_back(m);
            }
        }

        // include pure single-constituent phases only if their thermo
        // data are valid for this temperature. This is necessary,
        // since some thermo polynomial fits are done only for a
        // limited temperature range. For example, using the NASA
        // polynomial fits for solid ice and liquid water, if this
        // were not done the calculation would predict solid ice to be
        // present far above its melting point, since the thermo
        // polynomial fits only extend to 273.15 K, and give
        // unphysical results above this temperature, leading
        // (incorrectly) to Gibbs free energies at high temperature
        // lower than for liquid water.
        index_t ip;
        for (k = 0; k < m_nsp_mix; k++) {
            ip = m_mix->speciesPhaseIndex(k);
            if (!m_mix->solutionSpecies(k) && 
                !m_mix->tempOK(ip)) {
                m_incl_species[k] = 0;
                if (m_mix->speciesMoles(k) > 0.0) {
                    throw CanteraError("MultiPhaseEquil",
                        "condensed-phase species"+ m_mix->speciesName(k) 
                        + " is excluded since its thermo properties are \n"
                        "not valid at this temperature, but it has "
                        "non-zero moles in the initial state.");
                }
            }
        }

        // Now build the list of all species to be included in the
        // calculation.
        for (k = 0; k < m_nsp_mix; k++) {
            if (m_incl_species[k] ==1) {
                m_nsp++;
                m_species.push_back(k);
            }
        }

        // some work arrays for internal use
        m_work.resize(m_nsp);
        m_work2.resize(m_nsp);
        m_work3.resize(m_nsp_mix);
        m_mu.resize(m_nsp_mix);

        // number of moles of each species
        m_moles.resize(m_nsp);
        m_lastmoles.resize(m_nsp);
        m_dxi.resize(m_nsp - m_nel);

        // initialize the mole numbers to the mixture composition
        index_t ik;
        for (ik = 0; ik < m_nsp; ik++) {
            m_moles[ik] = m_mix->speciesMoles(m_species[ik]);
        }

        // Delta G / RT for each reaction
        m_deltaG_RT.resize(m_nsp - m_nel, 0.0);

        m_majorsp.resize(m_nsp);
        m_sortindex.resize(m_nsp,0);
        m_lastsort.resize(m_nel);
        m_solnrxn.resize(m_nsp - m_nel);
        m_A.resize(m_nel, m_nsp, 0.0);
        m_N.resize(m_nsp, m_nsp - m_nel);
        m_order.resize(m_nsp, 0);

        // if the 'start' flag is set, estimate the initial mole
        // numbers by doing a linear Gibbs minimization. In this case,
        // only the elemental composition of the initial mixture state
        // matters.
        if (start) {
            setInitialMoles(loglevel-1);
        }
        computeN();

        // Take a very small step in composition space, so that no
        // species has precisely zero moles.
        vector_fp dxi(m_nsp - m_nel, 1.0e-20);
        multiply(m_N, DATA_PTR(dxi), DATA_PTR(m_work));
        unsort(m_work);

        for (k = 0; k < m_nsp; k++) {
            m_moles[k] += m_work[k];
            m_lastmoles[k] = m_moles[k];
            if (m_mix->solutionSpecies(m_species[k])) 
                m_dsoln.push_back(1);
            else
                m_dsoln.push_back(0);
        }
        m_force = false; 
        updateMixMoles();

        // At this point, the instance has been created, the species
        // to be included have been determined, and an initial
        // composition has been selected that has all non-zero mole
        // numbers for the included species.
    }


    doublereal MultiPhaseEquil::equilibrate(int XY, doublereal err, 
        int maxsteps, int loglevel) {
        int i;
        m_iter = 0;
        string iterstr;
        if (loglevel > 0)
            beginLogGroup("MultiPhaseEquil::equilibrate", loglevel);

        for (i = 0; i < maxsteps; i++) {
            if (loglevel > 0) {
                iterstr = "iteration "+int2str(i);
                beginLogGroup(iterstr);
            }
            stepComposition(loglevel-1);
            if (loglevel > 0) {
                addLogEntry("error",fp2str(error()));
                endLogGroup(iterstr);
            }
            if (error() < err) break;
        }
        if (i >= maxsteps) {
            if (loglevel > 0) {
                addLogEntry("Error","no convergence in "+int2str(maxsteps)
                    +" iterations");
                endLogGroup("MultiPhaseEquil::equilibrate");
            }
            throw CanteraError("MultiPhaseEquil::equilibrate",
                "no convergence in " + int2str(maxsteps) + 
                " iterations. Error = " + fp2str(error()));
        }
        if (loglevel > 0) {
            addLogEntry("iterations",int2str(iterations()));
            addLogEntry("error tolerance",fp2str(err));
            addLogEntry("error",fp2str(error()));
            endLogGroup("MultiPhaseEquil::equilibrate");
        }
        finish();
        return error();
    }

    void MultiPhaseEquil::updateMixMoles() {
        fill(m_work3.begin(), m_work3.end(), 0.0);
        index_t k;
        for (k = 0; k < m_nsp; k++) {
            m_work3[m_species[k]] = m_moles[k];
        }
        m_mix->setMoles(DATA_PTR(m_work3));
    }

    /// Clean up the composition. The solution algorithm can leave
    /// some species in stoichiometric condensed phases with very
    /// small negative mole numbers. This method simply sets these to
    /// zero.
    void MultiPhaseEquil::finish() {
        fill(m_work3.begin(), m_work3.end(), 0.0);
        index_t k;
        for (k = 0; k < m_nsp; k++) {
            m_work3[m_species[k]] = (m_moles[k] > 0.0 ? m_moles[k] : 0.0);
        }
        m_mix->setMoles(DATA_PTR(m_work3));
    }


    /// Extimate the initial mole numbers. This is done by running
    /// each reaction as far forward or backward as possible, subject
    /// to the constraint that all mole numbers remain
    /// non-negative. Reactions for which \f$ \Delta \mu^0 \f$ are
    /// positive are run in reverse, and ones for which it is negative
    /// are run in the forward direction. The end result is equivalent
    /// to solving the linear programming problem of minimizing the
    /// linear Gibbs function subject to the element and
    /// non-negativity constraints.
    int MultiPhaseEquil::setInitialMoles(int loglevel) {
        index_t ik, j;

        double not_mu = 1.0e12;
        if (loglevel > 0)
            beginLogGroup("MultiPhaseEquil::setInitialMoles");

        m_mix->getValidChemPotentials(not_mu, DATA_PTR(m_mu), true);
        doublereal dg_rt;

        int idir;
        double nu;
        double delta_xi, dxi_min = 1.0e10;
        bool redo = true;
        int iter = 0;
        while (redo) {

            // choose a set of components based on the current
            // composition
            computeN();
            if (loglevel > 0)
                addLogEntry("iteration",iter);
            redo = false;
            iter++;
            if (iter > 4) break;

            // loop over all reactions
            for (j = 0; j < m_nsp - m_nel; j++) {
                dg_rt = 0.0;
                dxi_min = 1.0e10;
                for (ik = 0; ik < m_nsp; ik++) {
                    dg_rt += mu(ik) * m_N(ik,j);
                }
                // fwd or rev direction
                idir = (dg_rt < 0.0 ? 1 : -1);

                for (ik = 0; ik < m_nsp; ik++) {
                    nu = m_N(ik, j);
                    
                    // set max change in progress variable by
                    // non-negativity requirement
                    if (nu*idir < 0) {
                        delta_xi = fabs(moles(ik)/nu);
                        // if a component has nearly zero moles, redo
                        // with a new set of components
                        if (!redo && delta_xi < 1.0e-10 && ik < m_nel) {
                            if (loglevel > 0)
                                addLogEntry("component too small",speciesName(ik));
                            redo = true;
                        }
                        if (delta_xi < dxi_min) dxi_min = delta_xi;
                    }
                }
                // step the composition by dxi_min
                for (ik = 0; ik < m_nsp; ik++) {
                    moles(ik) += m_N(ik, j) * idir*dxi_min;
                }
            }
            // set the moles of the phase objects to match
            updateMixMoles();
        }
        for (ik = 0; ik < m_nsp; ik++) 
            if (moles(ik) != 0.0) addLogEntry(speciesName(ik), moles(ik));
        if (loglevel > 0)
            endLogGroup("MultiPhaseEquil::setInitialMoles");
        return 0;
    }


    ///  This method finds a set of component species and a complete
    ///  set of formation reactions for the non-components in terms of
    ///  the components. Note that in most cases, many different
    ///  component sets are possible, and therefore neither the
    ///  components returned by this method nor the formation
    ///  reactions are unique. The algorithm used here is described in
    ///  Smith and Missen, Chemical Reaction Equilibrium Analysis.
    /// 
    ///  The component species are taken to be the first M species
    ///  in array 'species' that have linearly-independent compositions.
    ///
    ///  @param order On entry, vector \a order should contain species
    ///  index numbers in the order of decreasing desirability as a
    ///  component. For example, if it is desired to choose the
    ///  components from among the major species, this array might
    ///  list species index numbers in decreasing order of mole
    ///  fraction. If array 'species' does not have length =
    ///  nSpecies(), then the species will be considered as candidates
    ///  to be components in declaration order, beginning with the
    ///  first phase added.
    ///
    void MultiPhaseEquil::getComponents(const vector_int& order) {
        index_t m, k, j;
        int n;

        // if the input species array has the wrong size, ignore it
        // and consider the species for components in declaration order.
        if (order.size() != m_nsp) {
            for (k = 0; k < m_nsp; k++) m_order[k] = k;
        }
        else {
            for (k = 0; k < m_nsp; k++) m_order[k] = order[k];
        }
        doublereal tmp;
        index_t itmp;

        index_t nRows = m_nel;
        index_t nColumns = m_nsp;
        doublereal fctr;

        // set up the atomic composition matrix
        for (m = 0; m < nRows; m++) {
            for (k = 0; k < nColumns; k++) {
                m_A(m, k) = m_mix->nAtoms(m_species[m_order[k]], m_element[m]);
            }
        }

        // Do Gauss elimination  
        for (m = 0; m < nRows; m++) {

            // If a pivot is zero, exchange columns.  This occurs when
            // a species has an elemental composition that is not
            // linearly independent of the component species that have
            // already been assigned
            if (m_A(m,m) == 0.0) {

                // First, we need to find a good candidate for a
                // component species to swap in for the one that has
                // zero pivot. It must contain element m, be linearly
                // independent of the components processed so far
                // (m_A(m,k) != 0), and should be a major species if
                // possible. We'll choose the species with greatest
                // mole fraction that satisfies these criteria.
                doublereal maxmoles = -999.0;
                index_t kmax = 0;
                for (k = m+1; k < nColumns; k++) {
                    if (m_A(m,k) != 0.0) {
                        if (fabs(m_moles[m_order[k]]) > maxmoles) {
                            kmax = k;
                            maxmoles = fabs(m_moles[m_order[k]]);
                        }
                    }
                }

                // Now exchange the column with zero pivot with the 
                // column for this major species
                for (n = 0; n < int(nRows); n++) {
                    tmp = m_A(n,m);
                    m_A(n, m) = m_A(n, kmax);
                    m_A(n, kmax) = tmp;
                }

                // exchange the species labels on the columns
                itmp = m_order[m];
                m_order[m] = m_order[kmax];
                m_order[kmax] = itmp;
            
            }

            // scale row m so that the diagonal element is unity
            fctr = 1.0/m_A(m,m);
            for (k = 0; k < nColumns; k++) {
                m_A(m,k) *= fctr;
            }

            // For all rows below the diagonal, subtract A(n,m)/A(m,m)
            // * (row m) from row n, so that A(n,m) = 0.
            for (n = int(m+1); n < int(m_nel); n++) {
                fctr = m_A(n,m)/m_A(m,m);
                for (k = 0; k < m_nsp; k++) {
                    m_A(n,k) -= m_A(m,k)*fctr;
                }
            }
        }


        // The left m_nel columns of A are now upper-diagonal.  Now
        // reduce the m_nel columns to diagonal form by back-solving
        for (m = nRows-1; m > 0; m--) {
            for (n = m-1; n>= 0; n--) {
                if (m_A(n,m) != 0.0) {
                    fctr = m_A(n,m);
                    for (k = m; k < m_nsp; k++) {
                        m_A(n,k) -= fctr*m_A(m,k);
                    }
                }
            }
        }

        // create stoichometric coefficient matrix. 
        for (n = 0; n < int(m_nsp); n++) {
            if (n < int(m_nel)) 
                for (k = 0; k < m_nsp - m_nel; k++) 
                    m_N(n, k) = -m_A(n, k + m_nel);
            else {
                for (k = 0; k < m_nsp - m_nel; k++) m_N(n, k) = 0.0;
                m_N(n, n - m_nel) = 1.0;
            }
        }

        // find reactions involving solution phase species
        for (j = 0; j < m_nsp - m_nel; j++) {
            m_solnrxn[j] = false;
            for (k = 0; k < m_nsp; k++) {
                if (m_N(k, j) != 0)
                    if (m_mix->solutionSpecies(m_species[m_order[k]])) 
                        m_solnrxn[j] = true;
            }
        }
    }




    /// Re-arrange a vector of species properties in sorted form
    /// (components first) into unsorted, sequential form.
    void MultiPhaseEquil::unsort(vector_fp& x) {
        copy(x.begin(), x.end(), m_work2.begin());
        index_t k;
        for (k = 0; k < m_nsp; k++) {
            x[m_order[k]] = m_work2[k];
        }
    }

#if defined(WITH_HTML_LOGS)
    void MultiPhaseEquil::printInfo(int loglevel) {
        index_t m, ik, k;
        if (loglevel > 0) {
            beginLogGroup("info");
            beginLogGroup("components");
        }
        for (m = 0; m < m_nel; m++) {
            ik = m_order[m];
            k = m_species[ik];
            if (loglevel > 0)
                addLogEntry(m_mix->speciesName(k), fp2str(m_moles[ik]));
        }
        if (loglevel > 0) {
            endLogGroup("components");
            beginLogGroup("non-components");
        }
        for (m = m_nel; m < m_nsp; m++) {
            ik = m_order[m];
            k = m_species[ik];
            if (loglevel > 0) 
                addLogEntry(m_mix->speciesName(k), fp2str(m_moles[ik]));
        }
        if (loglevel > 0) {
            endLogGroup("non-components");
            addLogEntry("Error",fp2str(error()));
            beginLogGroup("Delta G / RT");
        }
        for (k = 0; k < m_nsp - m_nel; k++) {
            if (loglevel > 0)
                addLogEntry(reactionString(k), fp2str(m_deltaG_RT[k]));
        }
        if (loglevel > 0) {
            endLogGroup("Delta G / RT");
            endLogGroup("info");
        }
    }

    /// Return a string specifying the jth reaction. 
    string MultiPhaseEquil::reactionString(index_t j) {
        string sr = "", sp = "";
        index_t i, k;
        bool rstrt = true;
        bool pstrt = true;
        doublereal nu;
        for (i = 0; i < m_nsp; i++) {
            nu = m_N(i, j);
            k = m_species[m_order[i]];
            if (nu < 0.0) {
                sr += coeffString(rstrt, nu, m_mix->speciesName(k));
                rstrt = false;
            }
            if (nu > 0.0) {
                sp += coeffString(pstrt, nu, m_mix->speciesName(k));
                pstrt = false;
            }
        }
        return sr + " <=> " + sp;
    }
#endif

    void MultiPhaseEquil::step(doublereal omega, vector_fp& deltaN, 
        int loglevel) {
        index_t k, ik;
        if (loglevel > 0)
            beginLogGroup("MultiPhaseEquil::step");
        if (omega < 0.0) 
            throw CanteraError("step","negative omega");

        for (ik = 0; ik < m_nel; ik++) {
            k = m_order[ik];
            m_lastmoles[k] = m_moles[k];
            if (loglevel > 0) {
                addLogEntry("component "+m_mix->speciesName(m_species[k])+" moles",
                m_moles[k]);
                addLogEntry("component "+m_mix->speciesName(m_species[k])+" step",
                omega*deltaN[k]);
            }
            m_moles[k] += omega * deltaN[k];
        }

        for (ik = m_nel; ik < m_nsp; ik++) {
            k = m_order[ik];
            m_lastmoles[k] = m_moles[k];
            if (m_majorsp[k]) {
                m_moles[k] += omega * deltaN[k];
            }
            else {
                m_moles[k] = fabs(m_moles[k])*fminn(10.0, 
                    exp(-m_deltaG_RT[ik - m_nel]));
            }
        }
        updateMixMoles();
        if (loglevel > 0)
            endLogGroup("MultiPhaseEquil::step");
    }


    /// Take one step in composition, given the gradient of G at the
    /// starting point, and a vector of reaction steps dxi. 
    doublereal MultiPhaseEquil::
    stepComposition(int loglevel) {
        if (loglevel > 0)
            beginLogGroup("MultiPhaseEquil::stepComposition");

        m_iter++;
        index_t ik, k = 0;
        doublereal grad0 = computeReactionSteps(m_dxi);

        // compute the mole fraction changes. 
        multiply(m_N, DATA_PTR(m_dxi), DATA_PTR(m_work));

        // change to sequential form
        unsort(m_work);

        // scale omega to keep the major species non-negative
        doublereal FCTR = 0.99;
        const doublereal MAJOR_THRESHOLD = 1.0e-12;

        doublereal omega = 1.0, omax, omegamax = 1.0;
        for (ik = 0; ik < m_nsp; ik++) {
            k = m_order[ik];
            if (ik < m_nel) {
                FCTR = 0.99;
                if (m_moles[k] < MAJOR_THRESHOLD) m_force = true;
            }
            else FCTR = 0.9;
            // if species k is in a multi-species solution phase, then its 
            // mole number must remain positive, unless the entire phase 
            // goes away. First we'll determine an upper bound on omega,
            // such that all 
            if (m_dsoln[k] == 1) {

                if ((m_moles[k] > MAJOR_THRESHOLD)  || (ik < m_nel)) {
                    if (m_moles[k] < MAJOR_THRESHOLD) m_force = true;
                    omax = m_moles[k]*FCTR/(fabs(m_work[k]) + TINY);
                    if (m_work[k] < 0.0 && omax < omegamax) {
                        omegamax = omax;
                        if (omegamax < 1.0e-5) {
                            m_force = true;
                        }
                    }
                    m_majorsp[k] = true;
                }
                else {
                    m_majorsp[k] = false;
                }
            }
            else {
                if (m_work[k] < 0.0 && m_moles[k] > 0.0) {
                    omax = -m_moles[k]/m_work[k];
                    if (omax < omegamax) {
                        omegamax = omax; //*1.000001;
                        if (omegamax < 1.0e-5) {
                            m_force = true;
                        }
                    }
                }
                if (m_moles[k] < -Tiny) {
                    if (loglevel > 0)
                        addLogEntry("Negative moles for "
                            +m_mix->speciesName(m_species[k]), fp2str(m_moles[k]));
                }
                m_majorsp[k] = true;
            }
        }

        // now take a step with this scaled omega
        if (loglevel > 0)
            addLogEntry("Stepping by ", fp2str(omegamax));
        step(omegamax, m_work);
        // compute the gradient of G at this new position in the
        // current direction. If it is positive, then we have overshot
        // the minimum. In this case, interpolate back.
        doublereal not_mu = 1.0e12;
        m_mix->getValidChemPotentials(not_mu, DATA_PTR(m_mu)); 
        doublereal grad1 = 0.0;
        for (k = 0; k < m_nsp; k++) {
            grad1 += m_work[k] * m_mu[m_species[k]];
        }

        omega = omegamax;
        if (grad1 > 0.0) {
            omega *= fabs(grad0) / (grad1 + fabs(grad0));
            for (k = 0; k < m_nsp; k++) m_moles[k] = m_lastmoles[k];
            if (loglevel > 0)
                addLogEntry("Stepped over minimum. Take smaller step ", fp2str(omega));
            step(omega, m_work);
        }
        printInfo(loglevel); 
        if (loglevel > 0)
            endLogGroup("MultiPhaseEquil::stepComposition");
        return omega;
    }

    
    /// Compute the change in extent of reaction for each reaction.

    doublereal MultiPhaseEquil::computeReactionSteps(vector_fp& dxi) {

        index_t j, k, ik, kc, ip;
        doublereal stoich, nmoles, csum, term1, fctr, rfctr;
        vector_fp nu;
        const doublereal TINY = 1.0e-20;
        doublereal grad = 0.0;

        dxi.resize(m_nsp - m_nel);
        computeN();
        doublereal not_mu = 1.0e12;
        m_mix->getValidChemPotentials(not_mu, DATA_PTR(m_mu));

        for (j = 0; j < m_nsp - m_nel; j++) {

            // get stoichiometric vector
            getStoichVector(j, nu);

            // compute Delta G
            doublereal dg_rt = 0.0;
            for (k = 0; k < m_nsp; k++) {
                dg_rt += m_mu[m_species[k]] * nu[k];
            }
            dg_rt /= (m_temp * GasConstant);

            m_deltaG_RT[j] = dg_rt;
            fctr = 1.0;

            // if this is a formation reaction for a single-component phase, 
            // check whether reaction should be included
            ik = j + m_nel;
            k = m_order[ik];
            if (!m_dsoln[k]) {
                if (m_moles[k] <= 0.0 && dg_rt > 0.0) {
                    fctr = 0.0;
                }
                else {
                    fctr = 0.5;
                }
            }
            else if (!m_solnrxn[j]) {
                fctr = 1.0;
            }
            else {

                // component sum
                csum = 0.0;
                for (k = 0; k < m_nel; k++) {
                    kc = m_order[k];
                    stoich = nu[kc];
                    nmoles = fabs(m_mix->speciesMoles(m_species[kc])) + TINY;
                    csum += stoich*stoich*m_dsoln[kc]/nmoles;
                }

                // noncomponent term
                kc = m_order[j + m_nel];  
                nmoles = fabs(m_mix->speciesMoles(m_species[kc])) + TINY;
                term1 = m_dsoln[kc]/nmoles;

                // sum over solution phases
                doublereal sum = 0.0, psum;
                for (ip = 0; ip < m_np; ip++) {
                    phase_t& p = m_mix->phase(ip);
                    if (p.nSpecies() > 1) {
                        psum = 0.0;
                        for (k = 0; k < m_nsp; k++) {
                            kc = m_species[k];
                            if (m_mix->speciesPhaseIndex(kc) == ip) {
                                // bug fixed 7/12/06 DGG
                                stoich = nu[k]; // nu[kc];
                                psum += stoich * stoich;
                            }
                        }
                        sum -= psum / (fabs(m_mix->phaseMoles(ip)) + TINY);
                    }
                }
                rfctr = term1 + csum + sum;
                if (fabs(rfctr) < TINY) 
                    fctr = 1.0;
                else
                    fctr = 1.0/(term1 + csum + sum);
            }
            dxi[j] = -fctr*dg_rt;
                
            index_t m;
             for (m = 0; m < m_nel; m++) {
                if (m_moles[m_order[m]] <= 0.0 && (m_N(m, j)*dxi[j] < 0.0))
                    dxi[j] = 0.0;
             }
            grad += dxi[j]*dg_rt;

        }
        return grad*GasConstant*m_temp;
    }

    void MultiPhaseEquil::computeN() {
        index_t m, k;

        // get the species moles

        // sort mole fractions
        doublereal molesum = 0.0;
        for (k = 0; k < m_nsp; k++) {
            m_work[k] = m_mix->speciesMoles(m_species[k]);
            m_sortindex[k] = k;
            molesum += m_work[k];
        }
        heapsort(m_work, m_sortindex);

        // reverse order in sort index
        index_t itmp;
        for (k = 0; k < m_nsp/2; k++) {
            itmp = m_sortindex[m_nsp-k-1];
            m_sortindex[m_nsp-k-1] = m_sortindex[k];
            m_sortindex[k] = itmp;
        }
        index_t ik, ij;
        bool ok;
        for (m = 0; m < m_nel; m++) {
            for (ik = 0; ik < m_nsp; ik++) {
                k = m_sortindex[ik];
                if (m_mix->nAtoms(m_species[k],m_element[m]) != 0) break;
            }
            ok = false;
            for (ij = 0; ij < m_nel; ij++) {
                if (int(k) == m_order[ij]) ok = true;
            }
            if (!ok || m_force) {
                getComponents(m_sortindex);
                m_force = true;
                break;
            }
        }
    }

    doublereal MultiPhaseEquil::error() {
        index_t j, ik, k;
        doublereal err, maxerr = 0.0;

        // examine every reaction
        for (j = 0; j < m_nsp - m_nel; j++) {
            ik = j + m_nel;
            k = m_order[ik];

            // don't require formation reactions for solution species
            // present in trace amounts to be equilibrated
            if (!isStoichPhase(ik) && fabs(moles(ik)) <= SmallNumber) {
                err = 0.0;
            }

            // for stoichiometric phase species, no error if not present and 
            // delta G for the formation reaction is positive 
             if (isStoichPhase(ik) && moles(ik) <= 0.0 && 
                m_deltaG_RT[j] >= 0.0) {
                err = 0.0;
            }
            else {
                err = fabs(m_deltaG_RT[j]);
            }
            if (err > maxerr) {
                maxerr = err;
            }
        }
        return maxerr;
    }
}
