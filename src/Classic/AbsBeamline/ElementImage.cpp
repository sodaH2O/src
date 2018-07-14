// ------------------------------------------------------------------------
// $RCSfile: ElementImage.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ElementImage
//   Contains two std::strings, representing the element name and type, and a
//   map of its attributes used to represent this element's state.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------
//

#include "AbsBeamline/ElementImage.h"


// Class ElementImage
// ------------------------------------------------------------------------


ElementImage::ElementImage():
    AttributeSet()
{}

ElementImage::ElementImage(const ElementImage &right):
    AttributeSet(right),
    elementName(right.elementName),
    elementType(right.elementType)
{}

ElementImage::ElementImage(const std::string &name, const std::string &type,
                           const AttributeSet &attrib):
    AttributeSet(attrib),
    elementName(name),
    elementType(type)
{}

ElementImage::~ElementImage()
{}


const ElementImage &ElementImage::operator=(const ElementImage &right) {
    AttributeSet::operator=(right);
    elementName = right.elementName;
    elementType = right.elementType;
    return *this;
}


const std::string &ElementImage::getName() const {
    return elementName;
}


void ElementImage::setName(const std::string &name) {
    elementName = name;
}


const std::string &ElementImage::getType() const {
    return elementType;
}


void ElementImage::setType(const std::string &type) {
    elementType = type;
}
