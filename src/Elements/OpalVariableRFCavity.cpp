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

#include "Elements/OpalVariableRFCavity.h"

#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Attributes/Attributes.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "AbsBeamline/VariableRFCavity.h"

extern Inform *gmsg;

const std::string OpalVariableRFCavity::doc_string =
      std::string("The \"VARIABLE_RF_CAVITY\" element defines an RF cavity ")+
      std::string("with time dependent frequency, phase and amplitude.");

OpalVariableRFCavity::OpalVariableRFCavity():
    OpalElement(SIZE, "VARIABLE_RF_CAVITY", doc_string.c_str()) {
    itsAttr[PHASE_MODEL] = Attributes::makeString("PHASE_MODEL",
                "The name of the phase time dependence model.");
    itsAttr[AMPLITUDE_MODEL] = Attributes::makeString("AMPLITUDE_MODEL",
                "The name of the amplitude time dependence model.");
    itsAttr[FREQUENCY_MODEL] = Attributes::makeString("FREQUENCY_MODEL",
                "The name of the frequency time dependence model.");
    itsAttr[WIDTH] = Attributes::makeReal("WIDTH",
                "Full width of the cavity.");
    itsAttr[HEIGHT] = Attributes::makeReal("HEIGHT",
                "Full height of the cavity.");
    registerStringAttribute("PHASE_MODEL");
    registerStringAttribute("AMPLITUDE_MODEL");
    registerStringAttribute("FREQUENCY_MODEL");
    registerRealAttribute("WIDTH");
    registerRealAttribute("HEIGHT");

    registerOwnership();

    setElement((new VariableRFCavity("VARIABLE_RF_CAVITY"))->makeAlignWrapper());
}

OpalVariableRFCavity::OpalVariableRFCavity(const std::string &name,
                                           OpalVariableRFCavity *parent) :
          OpalElement(name, parent) {
    VariableRFCavity *cavity = dynamic_cast<VariableRFCavity*>(
                                        parent->getElement()->removeWrappers());
    setElement((new VariableRFCavity(*cavity))->makeAlignWrapper());
}

OpalVariableRFCavity::~OpalVariableRFCavity() {
}

OpalVariableRFCavity *OpalVariableRFCavity::clone(const std::string &name) {
    return new OpalVariableRFCavity(name, this);
}

OpalVariableRFCavity *OpalVariableRFCavity::clone() {
    return new OpalVariableRFCavity(this->getOpalName(), this);
}

void OpalVariableRFCavity::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
    const VariableRFCavity* cavity = dynamic_cast<const VariableRFCavity*>(&base);
    if (cavity == NULL) {
        throw OpalException("OpalVariableRFCavity::fillRegisteredAttributes",
                            "Failed to cast ElementBase to a VariableRFCavity");
    }

    attributeRegistry["L"]->setReal(cavity->getLength());
    std::shared_ptr<AbstractTimeDependence> phase_model = cavity->getPhaseModel();
    std::shared_ptr<AbstractTimeDependence> freq_model = cavity->getFrequencyModel();
    std::shared_ptr<AbstractTimeDependence> amp_model = cavity->getAmplitudeModel();
    std::string phase_name = AbstractTimeDependence::getName(phase_model);
    std::string amp_name = AbstractTimeDependence::getName(amp_model);
    std::string freq_name = AbstractTimeDependence::getName(freq_model);
    attributeRegistry["PHASE_MODEL"]->setString(phase_name);
    attributeRegistry["AMPLITUDE_MODEL"]->setString(amp_name);
    attributeRegistry["FREQUENCY_MODEL"]->setString(freq_name);
    attributeRegistry["WIDTH"]->setReal(cavity->getWidth());
    attributeRegistry["HEIGHT"]->setReal(cavity->getHeight());
}

void OpalVariableRFCavity::update() {
    OpalElement::update();

    VariableRFCavity *cavity = dynamic_cast<VariableRFCavity*>(
                                                getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    cavity->setLength(length);
    std::string phaseName = Attributes::getString(itsAttr[PHASE_MODEL]);
    cavity->setPhaseName(phaseName);
    std::string ampName = Attributes::getString(itsAttr[AMPLITUDE_MODEL]);
    cavity->setAmplitudeName(ampName);
    std::string freqName = Attributes::getString(itsAttr[FREQUENCY_MODEL]);
    cavity->setFrequencyName(freqName);
    double width = Attributes::getReal(itsAttr[WIDTH]);
    cavity->setWidth(width);
    double height = Attributes::getReal(itsAttr[HEIGHT]);
    cavity->setHeight(height);
    setElement(cavity->makeAlignWrapper());
}