#ifndef CLASSIC_RBendRep_HH
#define CLASSIC_RBendRep_HH

// ------------------------------------------------------------------------
// $RCSfile: RBendRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBendRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RBend.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Fields/BMultipoleField.h"


// Class RBendRep
// ------------------------------------------------------------------------
/// Representation for a rectangular bend magnet.
//  A rectangular bend magnet has a rectilinear geometry about which its
//  multipole components are specified.

class RBendRep: public RBend {

public:

    /// Constructor with given name.
    explicit RBendRep(const std::string &name);

    RBendRep();
    RBendRep(const RBendRep &);
    virtual ~RBendRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

    /// Get field.
    //  Version for non-constant object.
    virtual BMultipoleField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const BMultipoleField &getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual RBendGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const RBendGeometry &getGeometry() const;

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

    /// Get field.
    //  Return the vertical component of the field in Teslas.
    virtual double getB() const;

    /// Set vertical component.
    //  Assign the vertical component of the field in Teslas.
    virtual void setB(double By);

    /// Set field.
    //  Assign the multipole expansion.
    virtual void setField(const BMultipoleField &field);

    /// Allow field errors.
    //  Build a FieldWrapper pointing to the bend and return a pointer to
    //  that wrapper.
    virtual ElementBase *makeFieldWrapper();

    /// Get pole entry face rotation.
    //  Return the rotation of the entry pole face with respect to the x-axis.
    //  A positive angle rotates the pole face normal away from the centre
    //  of the machine.
    virtual double getEntryFaceRotation() const;

    /// Get exit pole face rotation.
    //  Return the rotation of the exit pole face with respect to the x-axis.
    //  A positive angle rotates the pole face normal away from the centre
    //  of the machine.
    virtual double getExitFaceRotation() const;

    /// Get entry pole face curvature.
    //  Return the curvature of the entry pole face.
    //  A positive curvature creates a convex pole face.
    virtual double getEntryFaceCurvature() const;

    /// Get exit pole face curvature.
    //  Return the curvature of the exit pole face.
    //  A positive curvature creates a convex pole face.
    virtual double getExitFaceCurvature() const;

    /// Set pole entry face rotation.
    //  Return the rotation of the entry pole face with respect to the x-axis.
    //  A positive angle rotates the pole face normal away from the centre
    //  of the machine.
    virtual void setEntryFaceRotation(double e1);

    /// Set exit pole face rotation.
    //  Return the rotation of the exit pole face with respect to the x-axis.
    //  A positive angle rotates the pole face normal away from the centre
    //  of the machine.
    virtual void setExitFaceRotation(double e2);

    /// Set entry pole face curvature.
    //  Return the curvature of the entry pole face.
    //  A positive curvature creates a convex pole face.
    virtual void setEntryFaceCurvature(double h1);

    /// Set exit pole face curvature.
    //  Return the curvature of the exit pole face.
    //  A positive curvature creates a convex pole face.
    virtual void setExitFaceCurvature(double h2);

    /// Get number of slices.
    virtual double getSlices() const;

    /// Get stepsize.
    virtual double getStepsize() const;

    /// Set number of slices.
    virtual void setSlices(double sl);

    /// Set stepsize.
    virtual void setStepsize(double ds);

private:

    // Not implemented.
    void operator=(const RBendRep &);

    /// The bend geometry.
    RBendGeometry geometry;

    /// The multipole expansion.
    BMultipoleField field;

    // The pole face angles and curvatures.
    double rEntry;
    double rExit;
    double hEntry;
    double hExit;

    // Parameters that determine integration step-size.
    double slices;
    double stepsize;
};

#endif // CLASSIC_RBendRep_HH
