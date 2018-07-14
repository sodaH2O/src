#ifndef OPAL_OpalSBend_HH
#define OPAL_OpalSBend_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalSBend.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSBend
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalBend.h"

class OpalWake;
class ParticleMatterInteraction;

// Class OpalSBend
// ------------------------------------------------------------------------
/// The SBEND element.

class OpalSBend: public OpalBend {

public:

    /// Exemplar constructor.
    OpalSBend();

    virtual ~OpalSBend();

    /// Make clone.
    virtual OpalSBend *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC bend.
    virtual void update();

private:

    // Not implemented.
    OpalSBend(const OpalSBend &);
    void operator=(const OpalSBend &);

    // Clone constructor.
    OpalSBend(const std::string &name, OpalSBend *parent);

    OpalWake *owk_m;
    ParticleMatterInteraction *parmatint_m;
};

#endif // OPAL_OpalSBend_HH
