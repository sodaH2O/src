// ------------------------------------------------------------------------
// $RCSfile: Eigen.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Eigen
//   The class for OPAL EIGEN commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/Eigen.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Tables/Twiss.h"
#include "Utilities/OpalException.h"

#include <fstream>
#include <iomanip>
#if defined(__GNUC__) && __GNUC__ < 3
#include <strstream>
#else
#include <sstream>
#endif
#include <iostream>


// Class Eigen
// ------------------------------------------------------------------------

// The attributes of class Eigen.
namespace {
    enum {
        TABLE,       // The name of the table to be listed.
        FNAME,       // The name of the file to be written.
        SIZE
    };
}


Eigen::Eigen():
    Action(SIZE, "EIGEN",
           "The \"EIGEN\" statement lists the eigenvectors for a named "
           "\"TWISS\" table.") {
    itsAttr[TABLE] = Attributes::makeString
                     ("TABLE", "Name of table to be listed");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to receive output", "EIGEN");

    registerOwnership(AttributeHandler::STATEMENT);
}


Eigen::Eigen(const std::string &name, Eigen *parent):
    Action(name, parent)
{}


Eigen::~Eigen()
{}


Eigen *Eigen::clone(const std::string &name) {
    return new Eigen(name, this);
}


void Eigen::execute() {
    std::string tableName = Attributes::getString(itsAttr[TABLE]);
    Twiss *table = dynamic_cast<Twiss *>(OpalData::getInstance()->find(tableName));

    if(table) {
        std::string fileName = Attributes::getString(itsAttr[FNAME]);
        if(fileName == "TERM") {
            format(std::cout, table);
        } else {
            std::ofstream os(fileName.c_str());

            if(os.good()) {
                format(os, table);
            } else {
                throw OpalException("Eigen::execute()",
                                    "Unable to open output stream \"" +
                                    fileName + "\".");
            }
        }
    } else {
        throw OpalException("Eigen::execute()",
                            "Twiss table \"" + tableName + "\" not found.");
    }
}


void Eigen::format(std::ostream &os, const Twiss *table) {
        formatPrint(os, table);
}


void Eigen::formatPrint(std::ostream &os, const Twiss *table) const {
    // Save the formatting flags.
    std::streamsize old_prec = os.precision(6);
    os.setf(std::ios::fixed, std::ios::floatfield);

    // Print table specific header.
    table->printTableTitle(os, "Eigenvectors");
    os << std::string(118, '-') << '\n';
    os << "Element" << std::string(24, ' ') << "S       Orbit |"
       << std::string(25, ' ') << "E i g e n v e c t o r s\n";
    os << std::string(118, '-') << '\n';

    // Print table body.
    for(Twiss::TLine::const_iterator row = table->begin();
        row != table->end(); ++row) {
        if(row->getSelectionFlag()) {
            os << '\n';
            std::string name = row->getElement()->getName();
            if(int occur = row->getCounter()) {
#if defined(__GNUC__) && __GNUC__ < 3
                char buffer[128];
                std::ostrstream tos(buffer, 128);
#else
                std::ostringstream tos;
#endif
                tos << name << '[' << occur << ']' << std::ends;
#if defined(__GNUC__) && __GNUC__ < 3
                name = buffer;
#else
                name = tos.str();
#endif
            }

            if(name.length() > 16) {
                // Truncate the element name.
                os << std::string(name, 0, 13) << ".. ";
            } else {
                // Left adjust the element name.
                os << name << std::string(16 - name.length(), ' ');
            }
            os << std::setw(16) << table->getS(*row);

            FMatrix<double, 6, 6> eigen = table->getCurlyA(*row);
            FVector<double, 6> orbit = table->getOrbit(*row);
            for(int i = 0; i < 6; ++i) {
                if(i != 0) os << std::string(32, ' ');
                os << std::setw(12) << orbit[i] << " |";
                for(int j = 0; j < 6; ++j) {
                    os << std::setw(12)  << eigen[i][j];
                }
                os << '\n';
            }
        }
    }

    os << std::string(118, '-') << std::endl;

    // Restore the formatting flags.
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
}