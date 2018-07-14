// ------------------------------------------------------------------------
// $RCSfile: OpalWire.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalWire
//   The class of OPAL wire collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalWire.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"

// Class OpalWire
// ------------------------------------------------------------------------

OpalWire::OpalWire():
    OpalElement(SIZE, "WIRE",
                "The \"WIRE\" element defines a wire."),
    parmatint_m(NULL) {
    itsAttr[XSIZE] = Attributes::makeReal
                     ("XSIZE", "Horizontal half-aperture in m");
    itsAttr[YSIZE] = Attributes::makeReal
                     ("YSIZE", "Vertical half-aperture in m");
    itsAttr[XPOS] = Attributes::makeReal
                    ("XPOS", "Horizontal position in m");
    itsAttr[YPOS] = Attributes::makeReal
                    ("YPOS", "Vertical position in m");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "Wire output filename");

    registerStringAttribute("OUTFN");
    registerRealAttribute("XSIZE");
    registerRealAttribute("YSIZE");
    registerRealAttribute("XPOS");
    registerRealAttribute("YPOS");

    registerOwnership();

    setElement((new CCollimatorRep("WIRE"))->makeAlignWrapper());
}


OpalWire::OpalWire(const std::string &name, OpalWire *parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement((new CCollimatorRep(name))->makeAlignWrapper());
}


OpalWire::~OpalWire() {
    if(parmatint_m)
        delete parmatint_m;
}


OpalWire *OpalWire::clone(const std::string &name) {
    return new OpalWire(name, this);
}


void OpalWire::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    const CCollimatorRep *coll =
        dynamic_cast<const CCollimatorRep *>(base.removeWrappers());
    attributeRegistry["XSIZE"]->setReal(coll->getXsize());
    attributeRegistry["YSIZE"]->setReal(coll->getYsize());
}


void OpalWire::update() {
    OpalElement::update();

    CCollimatorRep *coll =
        dynamic_cast<CCollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    coll->setElementLength(length);
    coll->setXsize(Attributes::getReal(itsAttr[XSIZE]));
    coll->setYsize(Attributes::getReal(itsAttr[YSIZE]));
    coll->setXpos(Attributes::getReal(itsAttr[XPOS]));
    coll->setYpos(Attributes::getReal(itsAttr[YPOS]));
    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    coll->setWire();

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*coll);
        coll->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}