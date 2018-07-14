// ------------------------------------------------------------------------
// $RCSfile: SBendRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SBendRep
//   Defines a concrete representation for a sector (curved) bend.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/SBendRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndexedChannel.h"
#include "Channels/IndirectChannel.h"
#include "ComponentWrappers/SBendWrapper.h"
#include <cctype>

// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(SBendRep::*get)() const;
        void (SBendRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "L",
            &SBendRep::getElementLength,
            &SBendRep::setElementLength
        },
        {
            "BY",
            &SBendRep::getB,
            &SBendRep::setB
        },
        {
            "E1",
            &SBendRep::getEntryFaceRotation,
            &SBendRep::setEntryFaceRotation
        },
        {
            "E2",
            &SBendRep::getExitFaceRotation,
            &SBendRep::setExitFaceRotation
        },
        {
            "H1",
            &SBendRep::getEntryFaceCurvature,
            &SBendRep::setEntryFaceCurvature
        },
        {
            "H2",
            &SBendRep::getExitFaceCurvature,
            &SBendRep::setExitFaceCurvature
        },
        { 0, 0, 0 }
    };
}


// Class SBendRep
// ------------------------------------------------------------------------

SBendRep::SBendRep():
    SBend(),
    geometry(0.0, 0.0),
    field() {
    rEntry = rExit = hEntry = hExit = 0.0;
}


SBendRep::SBendRep(const SBendRep &rhs):
    SBend(rhs),
    geometry(rhs.geometry),
    field(rhs.field) {
    rEntry = rhs.rEntry;
    rExit  = rhs.hExit;
    hEntry = rhs.rEntry;
    hExit  = rhs.hExit;
}


SBendRep::SBendRep(const std::string &name):
    SBend(name),
    geometry(0.0, 0.0),
    field() {
    rEntry = rExit = hEntry = hExit = 0.0;
}


SBendRep::~SBendRep()
{}


ElementBase *SBendRep::clone() const {
    return new SBendRep(*this);
}


Channel *SBendRep::getChannel(const std::string &aKey, bool create) {
    if(aKey[0] == 'a'  ||  aKey[0] == 'b') {
        int n = 0;

        for(std::string::size_type k = 1; k < aKey.length(); k++) {
            if(isdigit(aKey[k])) {
                n = 10 * n + aKey[k] - '0';
            } else {
                return 0;
            }
        }

        if(aKey[0] == 'b') {
            return new IndexedChannel<SBendRep>
                   (*this, &SBendRep::getNormalComponent,
                    &SBendRep::setNormalComponent, n);
        } else {
            return new IndexedChannel<SBendRep>
                   (*this, &SBendRep::getSkewComponent,
                    &SBendRep::setSkewComponent, n);
        }
    } else {
        for(const Entry *table = entries; table->name != 0; ++table) {
            if(aKey == table->name) {
                return new IndirectChannel<SBendRep>(*this, table->get, table->set);
            }
        }

        return ElementBase::getChannel(aKey, create);
    }
}


BMultipoleField &SBendRep::getField() {
    return field;
}

const BMultipoleField &SBendRep::getField() const {
    return field;
}


PlanarArcGeometry &SBendRep::getGeometry() {
    return geometry;
}

const PlanarArcGeometry &SBendRep::getGeometry() const {
    return geometry;
}


ElementImage *SBendRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *table = entries; table->name != 0; ++table) {
        image->setAttribute(table->name, (this->*(table->get))());
    }

    for(int n = 1; n <= field.order(); n++) {
        char buffer[20];
        char *p = buffer;
        int k = n;

        while(k != 0) {
            *p++ = k % 10 + '0';
            k /= 10;
        }

        std::string name(" ");
        while(p > buffer) name += *--p;

        double b = field.getNormalComponent(n);
        if(b != 0.0) {
            name[0] = 'b';
            image->setAttribute(name, b);
        }

        double a = field.getSkewComponent(n);
        if(a != 0.0) {
            name[0] = 'a';
            image->setAttribute(name, a);
        }
    }

    return image;
}


double SBendRep::getB() const {
    return field.getNormalComponent(1);
}

void SBendRep::setB(double B) {
    field.setNormalComponent(1, B);
}


double SBendRep::getEntryFaceRotation() const {
    return rEntry;
}


double SBendRep::getExitFaceRotation() const {
    return rExit;
}

void SBendRep::setEntryFaceRotation(double e1) {
    rEntry = e1;
}

void SBendRep::setExitFaceRotation(double e2) {
    rExit = e2;
}

double SBendRep::getEntryFaceCurvature() const {
    return hEntry;
}

double SBendRep::getExitFaceCurvature() const {
    return hExit;
}

void SBendRep::setEntryFaceCurvature(double h1) {
    hEntry = h1;
}

void SBendRep::setExitFaceCurvature(double h2) {
    hExit = h2;
}


double SBendRep::getSlices() const {
    return slices;
}

double SBendRep::getStepsize() const {
    return stepsize;
}

void SBendRep::setSlices(double sl) {
    slices = sl;
}

void SBendRep::setStepsize(double ds) {
    stepsize = ds;
}


void SBendRep::setField(const BMultipoleField &f) {
    field = f;
}


ElementBase *SBendRep::makeFieldWrapper() {
    ElementBase *wrap = new SBendWrapper(this);
    wrap->setName(getName());
    return wrap;
}
