#ifndef CLASSIC_RBendWrapper_HH
#define CLASSIC_RBendWrapper_HH

// ------------------------------------------------------------------------
// $RCSfile: RBendWrapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBendWrapper
//
// ------------------------------------------------------------------------
// Class category: ComponentWrappers
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:54 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RBend.h"
#include "Fields/BMultipoleField.h"
#include "MemoryManagement/Pointer.h"

class BMultipoleField;


// Class RBendWrapper
// ------------------------------------------------------------------------
/// Representation of a perturbed rectangular bend.
//  A RBendWrapper represents a unique instance of a bend magnet
//  in the accelerator model. It defines imperfections of the field,
//  related to an ``ideal'' magnet contained in the wrapper.

class RBendWrapper: public RBend {

public:

    /// Constructor.
    //  Constructs a wrapper for its argument.
    //  The wrapper is not sharable by default.
    explicit RBendWrapper(RBend *);

    RBendWrapper(const RBendWrapper &right);
    virtual ~RBendWrapper();


    /// Apply visitor to modified bend.
    virtual void accept(BeamlineVisitor &) const;

    /// Make clone.
    virtual ElementBase *clone() const;

    /// Make structural copy.
    virtual ElementBase *copyStructure();

    /// Get multipole field error.
    //  This method can be used to get or set the error field. The error field
    //  offset is mutable, so as to allow changing it in a constant structure.
    virtual BMultipoleField &errorField() const;

    /// Set sharable flag.
    //  The whole structure depending on [b]this[/b] is marked as sharable.
    //  After this call a [b]copyStructure()[/b] call reuses the element.
    virtual void makeSharable();

    /// Get field.
    //  Return the perturbed field.
    //  Version for non-constant object.
    virtual BMultipoleField &getField();

    /// Get field.
    //  Return the perturbed field.
    //  Version for constant object.
    virtual const BMultipoleField &getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual RBendGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const RBendGeometry &getGeometry() const;

    /// Get dipole component.
    virtual double getB() const;

    /// Get pole face rotation.
    //  Return the rotation of the entry pole face with respect to the x-axis.
    //  A positive angle rotates the pole face normal away from the centre
    //  of the machine.
    virtual double getEntryFaceRotation() const;

    /// Get pole face rotation.
    //  Return the rotation of the exit pole face with respect to the x-axis.
    //  A positive angle rotates the pole face normal away from the centre
    //  of the machine.
    virtual double getExitFaceRotation() const;

    /// Get pole face curvature.
    //  Return the curvature of the entry pole face.
    //  A positive curvature creates a convex pole face.
    virtual double getEntryFaceCurvature() const;

    /// Get pole face curvature.
    //  Return the curvature of the exit pole face.
    //  A positive curvature creates a convex pole face.
    virtual double getExitFaceCurvature() const;

    /// Get number of slices.
    virtual double getSlices() const;

    /// Get stepsize.
    virtual double getStepsize() const;

    /// Get element type string.
    virtual ElementBase::ElementType getType() const;

    /// Get design RBend.
    //  Version for constant object.
    virtual const RBend &getDesign() const;

    /// Get design RBend.
    //  Version for non-constant object.
    virtual RBend &getDesign();

    /// Make wrapper for this bend.
    //  Return [b]this[/b], since this is already a field wrapper.
    virtual ElementBase *makeFieldWrapper();

    /// Remove field wrapper.
    virtual ElementBase *removeFieldWrapper();

    /// Remove field wrapper.
    virtual const ElementBase *removeFieldWrapper() const;

    /// Remove all wrappers.
    virtual ElementBase *removeWrappers();

    /// Remove all wrappers.
    virtual const ElementBase *removeWrappers() const;

private:

    // Not implemented.
    RBendWrapper();
    void operator=(const RBendWrapper &);

    // The ``design'' magnet.
    Pointer<RBend> itsDesign;

    /// The error field.
    mutable BMultipoleField itsError;

    // The magnetic field after application of all modifiers.
    mutable BMultipoleField tempField;
};

#endif // CLASSIC_RBendWrapper_HH