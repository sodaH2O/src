#ifndef OPAL_OpalCyclotronValley_HH
#define OPAL_OpalCyclotronValley_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalCyclotronValley.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalCyclotronValley
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

// Class OpalCyclotronValley
// ------------------------------------------------------------------------
/// The RFCAVITY element.

class OpalWake;

class OpalCyclotronValley: public OpalElement {

public:

    /// The attributes of class OpalCyclotronValley.
    enum {
        FMAPFN = COMMON,         // The filename of the fieldmap
        BFLG,
        SIZE
    };

    /// Exemplar constructor.
    OpalCyclotronValley();

    virtual ~OpalCyclotronValley();

    /// Make clone.
    virtual OpalCyclotronValley *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC cavity.
    virtual void update();

private:

    // Not implemented.
    OpalCyclotronValley(const OpalCyclotronValley &);
    void operator=(const OpalCyclotronValley &);

    // Clone constructor.
    OpalCyclotronValley(const std::string &name, OpalCyclotronValley *parent);


};

#endif // OPAL_OpalCyclotronValley_HH
