/**
 *  @file Group.h
 *
 * $Author: dggoodwin $
 * $Revision: 1.1 $
 * $Date: 2007/05/04 14:27:23 $
 */


// Copyright 2001  California Institute of Technology

#ifndef CT_RXNPATH_GROUP
#define CT_RXNPATH_GROUP

#include "ct_defs.h"

//using namespace std;
namespace Cantera {

    /**
     * Class Group is an internal class used by class ReactionPath. It
     * represents some subset of the atoms of a molecule.
     */
    class Group {
    public:
        Group() : m_sign(-999) { }
        Group(int n) : m_sign(0) { m_comp.resize(n,0);}
        Group(const vector_int& elnumbers) :
            m_comp(elnumbers), m_sign(0) {
                validate();
            }
        Group(const Group& g) :
            m_comp(g.m_comp), m_sign(g.m_sign) { }
        Group& operator=(const Group& g) {
            if (&g != this) {
                m_comp = g.m_comp;
                m_sign = g.m_sign;
            }
	    return *this;
        }
        virtual ~Group(){}

        /**
         * Decrement the atom numbers by those in group 'other'.
         */
        void operator-=(const Group& other) {
            verifyInputs(*this, other);
            int n = m_comp.size();
            for (int m = 0; m < n; m++)
                m_comp[m] -= other.m_comp[m];
            validate();
        }
        void operator+=(const Group& other) {
            verifyInputs(*this, other);
            int n = m_comp.size();
            for (int m = 0; m < n; m++)
                m_comp[m] += other.m_comp[m];
            validate();
        }
        void operator*=(int a) {
            int n = m_comp.size();
            for (int m = 0; m < n; m++)
                m_comp[m] *= a;
            validate();
        }
        bool operator==(const Group& other) const {
            verifyInputs(*this, other);
            int n = m_comp.size();
            for (int m = 0; m < n; m++) {
                if (m_comp[m] != other.m_comp[m]) return false;
            }
            return true;
        }
        friend Group operator-(const Group& g1, const Group& g2) {
            verifyInputs(g1, g2);
            Group diff(g1);
            diff -= g2;
            return diff;
        }
        friend Group operator+(const Group& g1, const Group& g2) {
            verifyInputs(g1, g2);
            Group sum(g1);
            sum += g2;
            return sum;
        }
        friend void verifyInputs(const Group& g1, const Group& g2) {
//             if (Debug::on) {
//                 if (g1.size() != g2.size()) {
//                     cerr << "Group: size mismatch!" << std::endl;
//                     cerr << " group 1 = " << g1 << std::endl;
//                     cerr << " group 2 = " << g2 << std::endl;
//                 }
//             }
        }

        void validate();

        /**
         * True if all non-zero atom numbers have the same sign.
         */
        bool valid() const { return (m_sign != -999); }
        bool operator!() const { return (m_sign == -999); }
        int sign() const { return m_sign; }
        int size() const { return m_comp.size(); }

        /// Number of atoms in the group (>= 0)
        int nAtoms() const {
            int n = m_comp.size();
            int sum = 0;
            for (int m = 0; m < n; m++) sum += std::abs(m_comp[m]);
            return sum;
        }
        /// Number of atoms of element m (positive or negative)
        int nAtoms(int m) const {
            if (m_comp.empty()) return 0;
            return m_comp[m];
        }

        std::ostream& fmt(std::ostream& s, const std::vector<std::string>& esymbols) const;

        friend std::ostream& operator<<(std::ostream& s,
                                        const Group& g);

    private:
        vector_int m_comp;
        int m_sign;
    };

}

#endif
