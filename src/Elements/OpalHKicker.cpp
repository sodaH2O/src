// ------------------------------------------------------------------------
// $RCSfile: OpalHKicker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalHKicker
//   The class of OPAL horizontal orbit correctors.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalHKicker.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/XCorrectorRep.h"
#include "ComponentWrappers/CorrectorWrapper.h"
#include "Physics/Physics.h"


// Class OpalHKicker
// ------------------------------------------------------------------------

OpalHKicker::OpalHKicker():
    OpalElement(SIZE, "HKICKER",
                "The \"HKICKER\" element defines a closed orbit corrector "
                "acting on the horizontal plane.") {
    itsAttr[KICK] = Attributes::makeReal
                    ("KICK", "Horizontal deflection in rad");
    itsAttr[DESIGNENERGY] = Attributes::makeReal
                           ("DESIGNENERGY", "the mean energy of the particles");
    itsAttr[K0] = Attributes::makeReal
                  ("K0", "Normal dipole field in T");

    registerRealAttribute("HKICK");
    registerRealAttribute("DESIGNENERGY");
    registerRealAttribute("K0");

    registerOwnership();

    setElement((new XCorrectorRep("HKICKER"))->makeWrappers());
}


OpalHKicker::OpalHKicker(const std::string &name, OpalHKicker *parent):
    OpalElement(name, parent) {
    setElement((new XCorrectorRep(name))->makeWrappers());
}


OpalHKicker::~OpalHKicker()
{}


OpalHKicker *OpalHKicker::clone(const std::string &name) {
    return new OpalHKicker(name, this);
}


void OpalHKicker::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
    const CorrectorWrapper *corr =
        dynamic_cast<const CorrectorWrapper *>(base.removeAlignWrapper());
    BDipoleField field;

    if(flag == ERROR_FLAG) {
        field = corr->errorField();
    } else if(flag == ACTUAL_FLAG) {
        field = corr->getField();
    } else if(flag == IDEAL_FLAG) {
        field = corr->getDesign().getField();
    }

    double scale = Physics::c / OpalData::getInstance()->getP0();
    attributeRegistry["HKICK"]->setReal(- field.getBy() * scale);
    attributeRegistry["VKICK"]->setReal(+ field.getBx() * scale);
}


void OpalHKicker::update() {
    OpalElement::update();

    XCorrectorRep *corr =
        dynamic_cast<XCorrectorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    double kick = Attributes::getReal(itsAttr[KICK]);

    corr->setElementLength(length);
    corr->setBy(- kick * factor);

    corr->setKickX(kick);
    if(itsAttr[DESIGNENERGY]) {
        double kineticEnergy = Attributes::getReal(itsAttr[DESIGNENERGY]) * 1e6;
        corr->setDesignEnergy(kineticEnergy, false);
    }

    if (itsAttr[K0]) {
        corr->setKickField(Vector_t(0, Attributes::getReal(itsAttr[K0]), 0));
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(corr);
}