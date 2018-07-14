/*
 *  Copyright (c) 2012, Chris Rogers
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

#include "Elements/OpalOffset/OpalLocalCartesianOffset.h"

#include <string>

#include "AbsBeamline/Offset.h"
#include "Utilities/OpalException.h"
#include "Elements/OpalElement.h"

#include "Attributes/Attributes.h"

namespace OpalOffset {

const std::string OpalLocalCartesianOffset::doc_string =
    std::string("The \"LOCAL_CARTESIAN_OFFSET\" element defines an offset")+
    std::string("in cartesian coordinates, relative to the last placed ")+
    std::string("element.");

OpalLocalCartesianOffset::OpalLocalCartesianOffset()
       : OpalElement(int(SIZE),
                     "LOCAL_CARTESIAN_OFFSET",
                     doc_string.c_str()) {
    itsAttr[END_POSITION_X] = Attributes::makeReal("END_POSITION_X",
             "x component of position of end of the offset in coordinate system of the end of the upstream element.");
    itsAttr[END_POSITION_Y] = Attributes::makeReal("END_POSITION_Y",
             "y component of position of end of the offset in coordinate system of the end of the upstream element.");
    itsAttr[END_NORMAL_X] = Attributes::makeReal("END_NORMAL_X",
             "x component of normal of end of the offset in coordinate system of the end of the upstream element.");
    itsAttr[END_NORMAL_Y] = Attributes::makeReal("END_NORMAL_Y",
             "y component of normal of end of the offset in coordinate system of the end of the upstream element.");
    registerRealAttribute("END_POSITION_X");
    registerRealAttribute("END_POSITION_Y");
    registerRealAttribute("END_NORMAL_X");
    registerRealAttribute("END_NORMAL_Y");

    registerOwnership();

    setElement((new Offset("LOCAL_CARTESIAN_OFFSET"))->makeAlignWrapper());
}

OpalLocalCartesianOffset* OpalLocalCartesianOffset::clone() {
    return new OpalLocalCartesianOffset(this->getOpalName(), this);
}


OpalLocalCartesianOffset* OpalLocalCartesianOffset::clone(const std::string &name) {
    return new OpalLocalCartesianOffset(name, this);
}

void OpalLocalCartesianOffset::print(std::ostream& out) const {
    OpalElement::print(out);
}

OpalLocalCartesianOffset::OpalLocalCartesianOffset(const std::string &name, OpalLocalCartesianOffset *parent):
    OpalElement(name, parent) {
    setElement(new Offset(name));
}

OpalLocalCartesianOffset::~OpalLocalCartesianOffset() {}

void OpalLocalCartesianOffset::fillRegisteredAttributes
                                     (const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
    const Offset* offset = dynamic_cast<const Offset*>(&base);
    if (offset == NULL) {
        throw OpalException("OpalVariableRFCavity::fillRegisteredAttributes",
                            "Failed to cast ElementBase to a VariableRFCavity");
    }

    Euclid3D trans = offset->getGeometry().getTotalTransform();
    double rot = trans.getRotation().getAxis()(1);
    attributeRegistry["END_POSITION_X"]->setReal(trans.getVector()(2));
    attributeRegistry["END_POSITION_Y"]->setReal(trans.getVector()(0));
    attributeRegistry["END_NORMAL_X"]->setReal(cos(rot));
    attributeRegistry["END_NORMAL_Y"]->setReal(sin(rot));
}

void OpalLocalCartesianOffset::update() {
    // getOpalName() comes from AbstractObjects/Object.h
    Offset *offset = dynamic_cast<Offset*>(getElement()->removeWrappers());
    std::string name = getOpalName();
    Vector_t pos(Attributes::getReal(itsAttr[END_POSITION_X]),
                 Attributes::getReal(itsAttr[END_POSITION_Y]), 0.);
    Vector_t norm(Attributes::getReal(itsAttr[END_NORMAL_X]),
                  Attributes::getReal(itsAttr[END_NORMAL_Y]), 0.);
    Offset off = Offset(Offset::localCartesianOffset(name, pos, norm));
    *offset = off;
    setElement(offset->makeAlignWrapper());
}
}