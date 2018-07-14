#ifndef CLASSIC_Lambertson_HH
#define CLASSIC_Lambertson_HH

// ------------------------------------------------------------------------
// $RCSfile: Lambertson.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Lambertson
//   *** MISSING *** Lambertson interface is still incomplete.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"


// Class Lambertson
// ------------------------------------------------------------------------
/// Interface for a Lambertson septum.
//  Class Lambertson defines the abstract interface for a Lambertson
//  septum magnet.

class Lambertson: public Component {

public:

    /// Constructor with given name.
    explicit Lambertson(const std::string &name);

    Lambertson();
    Lambertson(const Lambertson &);
    virtual ~Lambertson();

    /// Apply visitor to Lambertson.
    virtual void accept(BeamlineVisitor &) const;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);

    virtual void finalise();

    virtual bool bends() const;

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

private:

    // Not implemented.
    void operator=(const Lambertson &);
};

#endif // CLASSIC_Lambertson_HH