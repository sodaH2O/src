#ifndef OPAL_Insertion_HH
#define OPAL_Insertion_HH 1

// ------------------------------------------------------------------------
// $RCSfile: Insertion.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Insertion
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Tables/Twiss.h"


// Class Insertion
// ------------------------------------------------------------------------
/// The TWISSTRACK command.

class Insertion: public Twiss {

public:

    /// Exemplar construction.
    Insertion();

    virtual ~Insertion();

    /// Make clone.
    virtual Insertion *clone(const std::string &name);

    /// Fill the buffer using the defined algorithm.
    virtual void fill();

    /// Print the table on an ASCII stream.
    virtual void printTable(std::ostream &, const CellArray &) const;

private:

    // Not implemented.
    Insertion(const Insertion &);
    void operator=(const Insertion &);

    // Clone constructor.
    Insertion(const std::string &name, Insertion *parent);

    // The attribute of class Insertion.
    enum {
        // Input values.
        INIT = Twiss::SIZE,  // A table row defining the initial conditions.
        BETX,                // Initial lattice functions.
        ALFX,
        BETY,
        ALFY,
        DX,                  // Initial dispersion.
        DPX,
        DY,
        DPY,
        XC,                  // Initial orbit.
        PXC,
        YC,
        PYC,
        TC,
        PTC,
        // Read-only values.
        LENGTH,              // Total length.
        MU1,                 // Phases.
        MU2,
        MU3,
        DELTAP,              // Differential momentum variation
        SIZE
    };
};

#endif // OPAL_Insertion_HH
