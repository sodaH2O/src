#ifndef OPAL_TRIM_COIL_H
#define OPAL_TRIM_COIL_H

#include <string>
#include <memory>
#include "AbstractObjects/Definition.h"

class TrimCoil;

// Class OpalTrimCoil
// ------------------------------------------------------------------------
/// The TRIMCOIL definition.
//  A TRIMCOIL definition is used to define a trim coil which can be applied
//  to a Cyclotron

class OpalTrimCoil: public Definition {

public:

    /// Exemplar constructor.
    OpalTrimCoil();

    virtual ~OpalTrimCoil();

    /// Test if replacement is allowed.
    //  Can replace only by another OpalTrimCoil
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual OpalTrimCoil *clone(const std::string &name);

    /// Check the OpalTrimCoil data.
    virtual void execute();

    /// Find named trim coil.
    static OpalTrimCoil *find(const std::string &name);

    /// Update the OpalTrimCoil data.
    virtual void update();

    /// Print method, called at initialisation
    Inform& print(Inform& os) const;

    /// Initialise implementation
    void initOpalTrimCoil();

    /// Actual implementation
    std::unique_ptr<TrimCoil> trimcoil_m;

private:

    ///@{ Not implemented.
    OpalTrimCoil  (const OpalTrimCoil &) = delete;
    void operator=(const OpalTrimCoil &) = delete;
    ///@}
    /// Private copy constructor, called by clone
    OpalTrimCoil(const std::string &name, OpalTrimCoil *parent);

};

inline Inform &operator<<(Inform &os, const OpalTrimCoil &b) {
    return b.print(os);
}

#endif // OPAL_TRIM_COIL_H