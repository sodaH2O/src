#ifndef OPAL_OpalSlit_HH
#define OPAL_OpalSlit_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalSlit.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSlit
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;

// Class OpalSlit
// ------------------------------------------------------------------------
/// The ECOLLIMATOR element.

class OpalSlit: public OpalElement {

public:

    /// The attributes of class OpalSlit.
    enum {
        XSIZE = COMMON,  // The horizontal half-size.
        YSIZE,           // The vertical half-size.
        OUTFN,
        SIZE
    };

    /// Exemplar constructor.
    OpalSlit();

    virtual ~OpalSlit();

    /// Make clone.
    virtual OpalSlit *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalSlit(const OpalSlit &);
    void operator=(const OpalSlit &);

    // Clone constructor.
    OpalSlit(const std::string &name, OpalSlit *parent);

    ParticleMatterInteraction *parmatint_m;
};

#endif // OPAL_OpalSlit_HH
