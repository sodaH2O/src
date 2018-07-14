/*
 *  Copyright (c) 2017, Titus Dascalu
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */


#include "Elements/OpalMultipoleT.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include "Expressions/SValue.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include <iostream>
#if defined(__GNUC__) && __GNUC__ < 3
#include <strstream>
#else
#include <sstream>
#endif
#include <vector>


// Class OpalMultipoleT
// ------------------------------------------------------------------------

OpalMultipoleT::OpalMultipoleT():
    OpalElement(SIZE, "MULTIPOLET",
    "The \"MULTIPOLET\" element defines a combined function multipole.") {
    itsAttr[TP] = Attributes::makeRealArray
                  ("TP", "Transverse Profile derivatives in m^(-k)");
    itsAttr[LFRINGE] = Attributes::makeReal
                  ("LFRINGE", "The length of the left end field in m");
    itsAttr[RFRINGE] = Attributes::makeReal
                  ("RFRINGE", "The length of the right end field in m");
    itsAttr[HAPERT] = Attributes::makeReal
                  ("HAPERT", "The aperture width in m");
    itsAttr[VAPERT] = Attributes::makeReal
                  ("VAPERT", "The aperture height in m");
    itsAttr[ANGLE] = Attributes::makeReal
                  ("ANGLE", "The azimuthal angle of the magnet in ring (rad)");
    itsAttr[EANGLE] = Attributes::makeReal
                  ("EANGLE", "The entrance angle (rad)");
    itsAttr[MAXFORDER] = Attributes::makeReal
                  ("MAXFORDER", 
                   "Number of terms used in each field component");
    itsAttr[ROTATION] = Attributes::makeReal
                  ("ROTATION", 
                   "Rotation angle about its axis for skew elements (rad)");
    itsAttr[VARRADIUS] = Attributes::makeBool
                  ("VARRADIUS",
                   "Set true if radius of magnet is variable");
    itsAttr[VARSTEP] = Attributes::makeReal
                  ("VARSTEP",
                   "Step size used in rotating coords along ref trajectory");
    //registerRealAttribute("FRINGELEN");

    registerOwnership();

    setElement((new MultipoleT("MULTIPOLET"))->makeWrappers());
}


OpalMultipoleT::OpalMultipoleT(const std::string &name, 
			       OpalMultipoleT *parent):
    OpalElement(name, parent) {
    setElement((new MultipoleT(name))->makeWrappers());
}


OpalMultipoleT::~OpalMultipoleT()
{}


OpalMultipoleT *OpalMultipoleT::clone(const std::string &name) {
    return new OpalMultipoleT(name, this);
}


void OpalMultipoleT::print(std::ostream &os) const {
    OpalElement::print(os);
}


void OpalMultipoleT::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);   
    const MultipoleT *multT = 
        dynamic_cast<const MultipoleT*>(base.removeAlignWrapper());
    
    for(unsigned int order = 1; order <= multT->getTransMaxOrder(); order++) {
        std::ostringstream ss;
	ss << order;
	std::string orderString = ss.str();
        std::string attrName = "TP" + orderString;
	registerRealAttribute(attrName)->setReal(multT->getTransProfile(order));
    }

    registerRealAttribute("LFRINGE")->setReal(multT->getFringeLength().at(0));
    registerRealAttribute("RFRINGE")->setReal(multT->getFringeLength().at(1));
    registerRealAttribute("VAPERT")->setReal(multT->getAperture()[0]);
    registerRealAttribute("HAPERT")->setReal(multT->getAperture()[1]);
    registerRealAttribute("MAXFORDER")->setReal(multT->getMaxOrder());
    registerRealAttribute("ROTATION")->setReal(multT->getRotation());
    registerRealAttribute("EANGLE")->setReal(multT->getEntranceAngle());
    registerRealAttribute("VARSTEP")->setReal(multT->getVarStep());
    //registerRealAttribute("VARRADIUS")->setReal(multT->getVarRadius());
    
}


void OpalMultipoleT::update() {
    OpalElement::update();

    // Magnet length.
    MultipoleT *multT =
        dynamic_cast<MultipoleT*>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double angle = Attributes::getReal(itsAttr[ANGLE]);
    multT->setElementLength(length/2);
    multT->setLength(length);
    multT->setBendAngle(angle);
    multT->setAperture(Attributes::getReal(itsAttr[VAPERT]), 
		       Attributes::getReal(itsAttr[HAPERT]));
  
    multT->setFringeField(Attributes::getReal(itsAttr[LENGTH])/2,
                          Attributes::getReal(itsAttr[LFRINGE]),
                          Attributes::getReal(itsAttr[RFRINGE])); 
    if (Attributes::getBool(itsAttr[VARRADIUS])) {
        multT->setVarRadius();
	if(Attributes::getReal(itsAttr[VARSTEP])) {
	    multT->setVarStep(Attributes::getReal(itsAttr[VARSTEP]));
	}
    }
    const std::vector<double> transProfile = 
                              Attributes::getRealArray(itsAttr[TP]);
    int transSize = transProfile.size();

    multT->setTransMaxOrder(transSize - 1);
    multT->setMaxOrder(Attributes::getReal(itsAttr[MAXFORDER]));
    multT->setRotation(Attributes::getReal(itsAttr[ROTATION]));
    multT->setEntranceAngle(Attributes::getReal(itsAttr[EANGLE]));

    PlanarArcGeometry &geometry = multT->getGeometry();
    
    if(length) {
        geometry = PlanarArcGeometry(length/2, angle / length);
    } else {
        geometry = PlanarArcGeometry(angle);
    }
    
    for(int comp = 0; comp < transSize; comp++) {
        multT->setTransProfile(comp, transProfile[comp]);
    }
    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(multT);

    setElement(multT->makeWrappers());
}
