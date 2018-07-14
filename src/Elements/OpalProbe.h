#ifndef OPAL_OpalProbe_HH
#define OPAL_OpalProbe_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalProbe.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalProbe
//
// ------------------------------------------------------------------------
//
// $Date: 2009/10/07 10:08:06 $
// $Author: bi $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class OpalWake;

// Class OpalProbe
// ------------------------------------------------------------------------
/// The Probe element.


class OpalProbe: public OpalElement {

public:

    /// The attributes of class OpalProbe.
    enum {
        XSTART = COMMON, // Start of x coordinate
        XEND,            // End of x coordinate
        YSTART,          // Start of y coordinate
        YEND1,           // Not used now
        YEND,            // End of y coordinate
        WIDTH,           // Width of the probe
        SIZE
    };
    ///YEND1 is not used since it is treated as a string in the input file which should be a real argument.
    /// Exemplar constructor.
    OpalProbe();

    virtual ~OpalProbe();

    /// Make clone.
    virtual OpalProbe *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC septum.
    virtual void update();

private:

    // Not implemented.
    OpalProbe(const OpalProbe &);
    void operator=(const OpalProbe &);

    // Clone constructor.
    OpalProbe(const std::string &name, OpalProbe *parent);

    OpalWake *owk_m;
};

#endif // OPAL_OpalProbe_HH
