/**
 * @file EdgeKinetics.h
 *
 * $Author: dggoodwin $
 * $Revision: 1.2 $
 * $Date: 2007/05/18 10:54:54 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_EDGEKINETICS_H
#define CT_EDGEKINETICS_H

#include "InterfaceKinetics.h"

namespace Cantera {

    class EdgeKinetics : public InterfaceKinetics {

    public:

        /**
	 * Constructor
	 *
	 */
        EdgeKinetics() : InterfaceKinetics() {}

        /// Destructor.
        virtual ~EdgeKinetics() {}

	/**
	 * Identifies the subclass of the Kinetics manager type.
	 * These are listed in mix_defs.h.
	 */
        virtual int ID() { return cEdgeKinetics; }

	/**
	 * Identifies the subclass of the Kinetics manager type.
	 * These are listed in mix_defs.h.
	 */
        virtual int type() { return cEdgeKinetics; }

        // defined in InterfaceKinetics.cpp
        virtual void finalize();
    };
}

#endif
