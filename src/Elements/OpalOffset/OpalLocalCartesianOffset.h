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


#ifndef OPAL_ELEMENTS_OpalOffset_OpalLocalCartesianOffset_HH
#define OPAL_ELEMENTS_OpalOffset_OpalLocalCartesianOffset_HH

#include "AbsBeamline/Offset.h"
#include "Elements/OpalElement.h"

/** OpalOffset provides classes for making offsets between elements in a ring
 *  geometry; only condition is that they must be in the midplane
 */
namespace OpalOffset {

/** OpalLocalCartesianOffset provides classes for making offsets between
 *  elements in a ring geometry; only condition is that they must be in the
 *  midplane
 */
class OpalLocalCartesianOffset : public OpalElement {
  public:
    /** Enumeration maps to UI parameters */
    enum {
        END_POSITION_X = COMMON,
        END_POSITION_Y,
        END_NORMAL_X,
        END_NORMAL_Y,
        SIZE // size of the enum
    };

    /** Define mapping from enum variables to string UI parameter names */
    OpalLocalCartesianOffset();

    /** No memory allocated so does nothing */
    virtual ~OpalLocalCartesianOffset();

    /** Inherited copy constructor */
    virtual OpalLocalCartesianOffset *clone();
    /** Inherited copy constructor */
    virtual OpalLocalCartesianOffset *clone(const std::string &name);

    /** Calls fillRegisteredAttributes on the OpalElement */
    void fillRegisteredAttributes(const ElementBase &base, ValueFlag flag);

    /** Receive parameters from the parser and hand them off to the
     *  OpalCylindricalOffset
     */
    void update();

    /** Calls print on the OpalElement */
    virtual void print(std::ostream &) const;
  private:
    // Not implemented.
    OpalLocalCartesianOffset(const OpalLocalCartesianOffset &);
    void operator=(const OpalLocalCartesianOffset &);

    // Clone constructor.
    OpalLocalCartesianOffset(const std::string &name,
                                 OpalLocalCartesianOffset *parent);

    static const std::string doc_string;
};

}

#endif // #define OPAL_ELEMENTS_OpalOffset_OpalLocalCartesianOffset_HH
