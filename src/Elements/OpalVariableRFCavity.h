/* 
 *  Copyright (c) 2014, Chris Rogers
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

#ifndef OPAL_OPALVARIABLERFCAVITY_H
#define OPAL_OPALVARIABLERFCAVITY_H

#include "Elements/OpalElement.h"

/** OpalVariableRFCavity
 */
class OpalVariableRFCavity: public OpalElement {
  public:
    /** enum maps string to integer value for UI definitions */
    enum {
        PHASE_MODEL = COMMON,
        AMPLITUDE_MODEL,
        FREQUENCY_MODEL,
        WIDTH,
        HEIGHT,
        SIZE // size of the enum
    };

    /** Copy constructor **/
    OpalVariableRFCavity(const std::string &name, OpalVariableRFCavity *parent);

    /** Default constructor **/
    OpalVariableRFCavity();

    /** Inherited copy constructor
     *
     *  Call on a base class to instantiate an object of derived class's type
    **/
    OpalVariableRFCavity* clone();

    /** Inherited copy constructor
     *
     *  Call on a base class to instantiate an object of derived class's type
     */
    virtual OpalVariableRFCavity *clone(const std::string &name);

    /** Destructor does nothing */
    virtual ~OpalVariableRFCavity();

    /** Fill in all registered attributes
     *
     *  This updates the registered attributed with values from the ElementBase
     */
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /** Update the OpalVariableRFCavity with new parameters from UI parser */
    virtual void update();

  private:
    // Not implemented.
    OpalVariableRFCavity(const OpalVariableRFCavity &);
    void operator=(const OpalVariableRFCavity &);

    static const std::string doc_string;
};

#endif // OPAL_OPALVARIABLERFCAVITY_H
