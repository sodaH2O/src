#ifndef OPAL_Wake_HH
#define OPAL_Wake_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalWake.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalWake
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"

class ElementBase;
class WakeFunction;

// Class OpalWake
// ------------------------------------------------------------------------
/// The WAKE definition.
//  A WAKE definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.

class OpalWake: public Definition {

public:

    /// Exemplar constructor.
    OpalWake();

    virtual ~OpalWake();

    /// Test if replacement is allowed.
    //  Can replace only by another WAKE.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual OpalWake *clone(const std::string &name);

    /// Check the WAKE data.
    virtual void execute();

    /// Find named WAKE.
    static OpalWake *find(const std::string &name);

    /// Update the WAKE data.
    virtual void update();

    void print(std::ostream &os) const;

    int getNumberOfBins();

    void initWakefunction(ElementBase &element);

    WakeFunction *wf_m;

private:

    // Not implemented.
    OpalWake(const OpalWake &);
    void operator=(const OpalWake &);

    // Clone constructor.
    OpalWake(const std::string &name, OpalWake *parent);

    // The particle reference data.
    PartData reference;

    // The conversion from GeV to eV.
    static const double energy_scale;

    // the element the wake is attached to
    ElementBase *itsElement_m;

};

inline std::ostream &operator<<(std::ostream &os, const OpalWake &b) {
    b.print(os);
    return os;
}

#endif // OPAL_Wake_HH
