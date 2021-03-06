#ifndef CLASSIC_OpalParticle_HH
#define CLASSIC_OpalParticle_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalParticle.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalParticle
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------


// Class OpalParticle
// ------------------------------------------------------------------------
/// OpalParticle position.
//  This class represents the canonical coordinates of a particle.
//  {P}
//  NOTE. The order of phase space coordinates is,
//  as opposed to the original CLASSIC definition:
//  {CENTER}
//  X = 0, PX = 1, Y = 2, PY = 3, T = 4, PT = 5.
//  {/CENTER}
//  The copy constructor, destructor, and assignment operator generated
//  by the compiler perform the correct operation.  For speed reasons
//  they are not implemented.

class OpalParticle {

public:

    // Particle coordinate numbers.
    enum { X, PX, Y, PY, T, PT };

    /// Constructor.
    //  Construct particle with the given coordinates.
    OpalParticle(double x, double px, double y, double py, double t, double pt);

    OpalParticle();

    /// Get coordinate.
    //  Access coordinate by index.
    //  Note above order of phase space coordinates!
    double &operator[](int);

    /// Get reference to horizontal position in m.
    double &x() ;

    /// Get reference to horizontal momentum (no dimension).
    double &px();

    /// Get reference to vertical displacement in m.
    double &y() ;

    /// Get reference to vertical momentum (no dimension).
    double &py();

    /// Get reference to longitudinal displacement c*t in m.
    double &t() ;

    /// Get reference to relative momentum error (no dimension).
    double &pt();

    /// Get coordinate.
    //  Access coordinate by index for constant particle.
    double operator[](int) const;

    /// Get horizontal position in m.
    double x() const;

    /// Get horizontal momentum (no dimension).
    double px() const;

    /// Get vertical displacement in m.
    double y() const;

    /// Get vertical momentum (no dimension).
    double py() const;

    /// Get longitudinal displacement c*t in m.
    double t() const;

    /// Get relative momentum error (no dimension).
    double pt() const;

private:

    // The particle's phase space coordinates:
    double phase[6];
};


// Inline member functions.
// ------------------------------------------------------------------------

inline double &OpalParticle::operator[](int i) {
    return phase[i];
}

inline double &OpalParticle::x()
{ return phase[X]; }

inline double &OpalParticle::y()
{ return phase[Y]; }

inline double &OpalParticle::t()
{ return phase[T]; }

inline double &OpalParticle::px()
{ return phase[PX]; }

inline double &OpalParticle::py()
{ return phase[PY]; }

inline double &OpalParticle::pt() {
    return phase[PT];
}


inline double OpalParticle::operator[](int i) const {
    return phase[i];
}

inline double OpalParticle::x() const
{ return phase[X]; }

inline double OpalParticle::y() const
{ return phase[Y]; }

inline double OpalParticle::t() const
{ return phase[T]; }

inline double OpalParticle::px() const
{ return phase[PX]; }

inline double OpalParticle::py() const
{ return phase[PY]; }

inline double OpalParticle::pt() const {
    return phase[PT];
}

#endif // CLASSIC_OpalParticle_HH
