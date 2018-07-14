// ------------------------------------------------------------------------
// $RCSfile: Select.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Select
//   Class for the OPAL SELECT command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "BasicActions/Select.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Table.h"
#include "Algorithms/Flagger.h"
#include "Attributes/Attributes.h"
#include "Tables/Selector.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <iostream>


// Class Select
// ------------------------------------------------------------------------

// The attributes of class Select.

namespace {
    enum {
        LINE,    // The line to be affected.
        FULL,    // If true, all elements are selected.
        CLEAR,   // If true, all selections are cleared.
        RANGE,   // The range to be considered.
        CLASS,   // The class of elements to be selected.
        TYPE,    // The type name of elements to be selected.
        PATTERN, // The regular expression for matching names.
        SIZE
    };
}

extern Inform *gmsg;

Select::Select():
    Action(SIZE, "SELECT",
           "The \"SELECT\" sub-command selects the positions to be affected "
           "by subsequent error sub-commands.") {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE",
                     "Name of the lattice to be affected by selections",
                     "UNNAMED_USE");
    itsAttr[FULL] = Attributes::makeBool
                    ("FULL",
                     "If true, all elements are selected");
    itsAttr[CLEAR] = Attributes::makeBool
                     ("CLEAR",
                      "If true, all selections are cleared");
    itsAttr[RANGE] = Attributes::makeRange
                     ("RANGE",
                      "Range to be considered for selection (default: full range)");
    itsAttr[CLASS] = Attributes::makeString
                     ("CLASS",
                      "Name of class to be selected (default: all classes)");
    itsAttr[TYPE] = Attributes::makeString
                    ("TYPE",
                     "The type name of elements to be selected (default: all types)");
    itsAttr[PATTERN] = Attributes::makeString
                       ("PATTERN",
                        "Regular expression for matching names (default: all names)");

    registerOwnership(AttributeHandler::SUB_COMMAND);
}


Select::Select(const std::string &name, Select *parent):
    Action(name, parent)
{}


Select::~Select()
{}


Select *Select::clone(const std::string &name) {
    return new Select(name, this);
}


void Select::execute() {
    // Find beam sequence  or table definition.
    const std::string name = Attributes::getString(itsAttr[LINE]);

    if(Object *obj = OpalData::getInstance()->find(name)) {
        if(BeamSequence *line = dynamic_cast<BeamSequence *>(obj)) {
            select(*line->fetchLine());
        } else if(Table *table = dynamic_cast<Table *>(obj)) {
            select(*table->getLine());
        } else {
            throw OpalException("Select::execute()",
                                "You cannot do a \"SELECT\" on \"" + name +
                                "\", it is neither a line nor a table.");
        }
    } else {
        throw OpalException("Select::execute()",
                            "Object \"" + name + "\" not found.");
    }
}


void Select::select(const Beamline &bl) {
    if(Attributes::getBool(itsAttr[FULL])) {
        // Select all positions.
        Flagger flagger(bl, true);
        flagger.execute();
        if(Options::info) {
            *gmsg << level3 << "\nAll elements selected.\n" << endl;
        }
    } else if(Attributes::getBool(itsAttr[CLEAR])) {
        // Deselect all selections.
        Flagger flagger(bl, false);
        flagger.execute();
        if(Options::info) {
            *gmsg << level3 << "\nAll elements de-selected.\n" << endl;
        }
    } else {
        Selector sel(bl,
                     Attributes::getRange(itsAttr[RANGE]),
                     Attributes::getString(itsAttr[CLASS]),
                     Attributes::getString(itsAttr[TYPE]),
                     Attributes::getString(itsAttr[PATTERN]));
        sel.execute();

        int count = sel.getCount();

        if(count == 0) {
            *gmsg << level3 << "No elements";
        } else if(count == 1) {
            *gmsg << level3 << "\n1 element";
        } else {
            *gmsg << level3 << '\n' << count << " elements";
        }

        *gmsg << level3 << " selected.\n" << endl;
    }
}