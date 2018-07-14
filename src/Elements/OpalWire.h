#ifndef OPAL_OpalWire_HH
#define OPAL_OpalWire_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalSlit.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalWire
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;
// Class OpalWire
// ------------------------------------------------------------------------
/// The ECOLLIMATOR element.

class OpalWire: public OpalElement {

public:

    /// The attributes of class OpalWire.
    enum {
        XSIZE = COMMON,  // The horizontal half-size.
        YSIZE,
        XPOS,  // The horizontal position.
        YPOS,           // The vertical position.
        OUTFN,
        SIZE
    };

    /// Exemplar constructor.
    OpalWire();

    virtual ~OpalWire();

    /// Make clone.
    virtual OpalWire *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalWire(const OpalWire &);
    void operator=(const OpalWire &);

    // Clone constructor.
    OpalWire(const std::string &name, OpalWire *parent);
    ParticleMatterInteraction *parmatint_m;
};

#endif // OPAL_OpalWire_HH
