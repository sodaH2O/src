#include <map>
#include <limits>
#include <iostream>
#include <fstream>

#include "Algorithms/IndexMap.h"
#include "AbstractObjects/OpalData.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Bend.h"
#include "Utilities/Util.h"
#include "Physics/Physics.h"
#include "OPALconfig.h"

#include <boost/filesystem.hpp>

extern Inform *gmsg;

const double IndexMap::oneMinusEpsilon_m = 1.0 - std::numeric_limits<double>::epsilon();

IndexMap::IndexMap():
    mapRange2Element_m(),
    mapElement2Range_m(),
    totalPathLength_m(0.0)
{ }

void IndexMap::print(std::ostream &out) const {
    out << "Size of map " << mapRange2Element_m.size() << " sections " << std::endl;
    out << std::fixed << std::setprecision(6);
    auto mapIti = mapRange2Element_m.begin();
    auto mapItf = mapRange2Element_m.end();

    double totalLength = (*mapRange2Element_m.rbegin()).first.second;
    unsigned int numDigits = std::floor(std::max(0.0, log(totalLength) / log(10.0))) + 1;

    for (; mapIti != mapItf; mapIti++) {
        const auto key = (*mapIti).first;
        const auto val = (*mapIti).second;
        out << "Key: ("
            << std::setw(numDigits + 7) << std::right << key.first
            << " - "
            << std::setw(numDigits + 7) << std::right << key.second
            << ") number of overlapping elements " << val.size() << "\n";

        for (auto element: val) {
            out << std::setw(25 + 2 * numDigits) << " " << element->getName() << "\n";
        }
    }
}

IndexMap::value_t IndexMap::query(key_t::first_type s, key_t::second_type ds) {
    const double lowerLimit = (ds < s? s - ds: 0);
    const double upperLimit = std::min(totalPathLength_m, s + ds);
    value_t elementSet;

    map_t::reverse_iterator rit = mapRange2Element_m.rbegin();
    std::cout << "(rit != mapRange2Element_m.rend()) is :"
	      << (rit != mapRange2Element_m.rend()) << "\n";
    std::cout << "(lowerLimit > (*rit).first.second) is :"
	      << (lowerLimit > (*rit).first.second) << "\n";
    if (rit != mapRange2Element_m.rend() && lowerLimit > (*rit).first.second) {
        throw OutOfBounds("IndexMap::query", "out of bounds");
    }

    map_t::iterator it = mapRange2Element_m.begin();
    const map_t::iterator end = mapRange2Element_m.end();

    for (; it != end; ++ it) {
        const double low = (*it).first.first;
        const double high = (*it).first.second;

        if (lowerLimit < high && upperLimit >= low) break;
    }

    if (it == end) return elementSet;

    map_t::iterator last = std::next(it);
    for (; last != end; ++ last) {
        const double low = (*last).first.first;

        if (upperLimit < low) break;
    }

    for (; it != last; ++ it) {
        const value_t &a = (*it).second;
        elementSet.insert(a.cbegin(), a.cend());
    }

    return elementSet;
}

void IndexMap::add(key_t::first_type initialS, key_t::second_type finalS, const value_t &val) {
    key_t key(initialS, finalS * oneMinusEpsilon_m);
    mapRange2Element_m.insert(std::pair<key_t, value_t>(key, val));
    totalPathLength_m = (*mapRange2Element_m.rbegin()).first.second;

    value_t::iterator setIt = val.begin();
    const value_t::iterator setEnd = val.end();

    for (; setIt != setEnd; ++ setIt) {
        if (mapElement2Range_m.find(*setIt) == mapElement2Range_m.end()) {
            mapElement2Range_m.insert(std::make_pair(*setIt, key));
        } else {
            auto itpair = mapElement2Range_m.equal_range(*setIt);

            bool extendedExisting = false;
            for (auto it = itpair.first; it != itpair.second; ++ it) {
                key_t &currentRange = it->second;

                if (almostEqual(key.first, currentRange.second / oneMinusEpsilon_m)) {
                    currentRange.second = key.second;
                    extendedExisting = true;
                    break;
                }
            }
            if (!extendedExisting) {
                mapElement2Range_m.insert(std::make_pair(*setIt, key));
            }
        }
    }
}

void IndexMap::tidyUp() {
    map_t::reverse_iterator rit = mapRange2Element_m.rbegin();

    if (rit != mapRange2Element_m.rend() && (*rit).second.size() == 0) {
        mapRange2Element_m.erase(std::next(rit).base());
    }
}

enum elements {
    DIPOLE = 0,
    QUADRUPOLE,
    SEXTUPOLE,
    OCTUPOLE,
    DECAPOLE,
    MULTIPOLE,
    SOLENOID,
    RFCAVITY,
    MONITOR,
    SIZE
};

void IndexMap::saveSDDS(double startS) const {

    std::string fileName("data/" + OpalData::getInstance()->getInputBasename() + "_ElementPositions.sdds");
    std::ofstream sdds;
    if (OpalData::getInstance()->hasPriorTrack() && boost::filesystem::exists(fileName)) {
        Util::rewindLinesSDDS(fileName, startS, false);
        sdds.open(fileName, std::ios::app);
    } else {
        sdds.open(fileName);

        std::string indent("        ");

        sdds << "SDDS1\n"
             << "&description \n"
             << indent << "text=\"Element positions '" << OpalData::getInstance()->getInputFn() << "'\", \n"
             << indent << "contents=\"element positions\" \n"
             << "&end\n";
        sdds << "&parameter\n"
             << indent << "name=revision, \n"
             << indent << "type=string,\n"
             << indent << "description=\"git revision of opal\"\n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=s, \n"
             << indent << "type=double, \n"
             << indent << "units=m, \n"
             << indent << "description=\"1 longitudinal position\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=dipole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"2 dipole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=quadrupole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"3 quadrupole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=sextupole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"4 sextupole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=octupole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"5 octupole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=decapole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"6 decapole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=multipole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"7 higher multipole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=solenoid, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"8 solenoid field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=rfcavity, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"9 RF field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=monitor, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"10 monitor present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=element_names, \n"
             << indent << "type=string, \n"
             << indent << "description=\"names of elements\" \n"
             << "&end\n";
        sdds << "&data \n"
             << indent << "mode=ascii, \n"
             << indent << "no_row_counts=1 \n"
             << "&end\n";

        sdds << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " # git rev. " << Util::getGitRevision() << std::endl;
    }

    std::vector<std::vector<int> > allItems(SIZE);
    std::vector<double> allPositions;
    std::vector<std::string> allNames;
    std::vector<double> typeMultipliers = {3.3333e-1, 1.0, 0.5, 0.25, 1.0, 1.0, 1.0, 1.0, 1.0};
    unsigned int counter = 0;

    auto mapIti = mapRange2Element_m.begin();
    auto mapItf = mapRange2Element_m.end();
    for (; mapIti != mapItf; mapIti++) {
        const auto key = (*mapIti).first;
        const auto val = (*mapIti).second;

        std::vector<int> items(SIZE, 0);
        std::string names("");
        for (auto element: val) {
            switch (element->getType()) {
            case ElementBase::RBEND:
            case ElementBase::SBEND:
                {
                    const Bend* bend = static_cast<const Bend*>(element.get());
                    if (bend->getRotationAboutZ() > 0.5 * Physics::pi &&
                        bend->getRotationAboutZ() < 1.5 * Physics::pi) {
                        items[DIPOLE] = -1;
                    } else {
                        items[DIPOLE] = 1;
                    }
                }
                break;
            case ElementBase::MULTIPOLE:
                {
                    const Multipole* mult = static_cast<const Multipole*>(element.get());
                    switch(mult->getMaxNormalComponentIndex()) {
                    case 1:
                        items[DIPOLE] = (mult->isFocusing(0)? 1: -1);
                        break;
                    case 2:
                        items[QUADRUPOLE] = (mult->isFocusing(1)? 1: -1);
                        break;
                    case 3:
                        items[SEXTUPOLE] = (mult->isFocusing(2)? 1: -1);
                        break;
                    case 4:
                        items[OCTUPOLE] = (mult->isFocusing(3)? 1: -1);
                        break;
                    case 5:
                        items[DECAPOLE] = (mult->isFocusing(4)? 1: -1);
                        break;
                    default:
                        items[MULTIPOLE] = 1;
                    }
                }
                break;
            case ElementBase::SOLENOID:
                items[SOLENOID] = 1;
                break;
            case ElementBase::RFCAVITY:
            case ElementBase::TRAVELINGWAVE:
                items[RFCAVITY] = 1;
                break;
            case ElementBase::MONITOR:
                items[MONITOR] = 1;
                break;
            default:
                continue;
            }

            names += element->getName() + ", ";
        }

        if (names.length() == 0) continue;

        allPositions.push_back(key.first);
        allPositions.push_back(key.first);
        allPositions.push_back(key.second);
        allPositions.push_back(key.second);

        allNames.push_back("");
        allNames.push_back(names.substr(0, names.length() - 2));
        allNames.push_back(allNames.back());
        allNames.push_back("");

        for (unsigned int i = 0; i < SIZE; ++ i) {
            allItems[i].push_back(0);
            allItems[i].push_back(items[i]);
            allItems[i].push_back(items[i]);
            allItems[i].push_back(0);
        }

        ++ counter;
    }

    if (counter == 0) return;

    const unsigned int totalSize = counter;
    for (unsigned int i = 0; i < SIZE; ++ i) {
        for (unsigned int j = 0; j < totalSize - 1; ++ j) {
            if (allItems[i][4 * j + 1] == 0) continue;
            if (allItems[i][4 * (j + 1) + 1] == 0) continue;
            if (allItems[i][4 * j + 1] * allItems[i][4 * (j + 1) + 1] < 0) continue;
            if (std::abs(allPositions[4 * j + 3] - allPositions[4 * (j + 1)]) > 1e-6) continue;

            allItems[i][4 * j + 3] = allItems[i][4 * (j + 1)] = allItems[i][4 * j + 1];
        }
    }

    for (unsigned int j = 0; j < totalSize - 1; ++ j) {
        allItems[RFCAVITY][4 * j + 2] = -allItems[RFCAVITY][4 * j + 1];
    }

    for (unsigned int i = 0; i < 4 * totalSize; ++ i) {
        sdds << std::setw(16) << std::setprecision(8) << allPositions[i];

        sdds << std::setw(10) << std::setprecision(5) << allItems[0][i] * typeMultipliers[0];
        for (unsigned int j = 1; j < SIZE; ++ j)
            sdds << std::setw(6) << std::setprecision(3) << allItems[j][i] * typeMultipliers[j];
        sdds << " \"" << allNames[i] << "\"\n";
    }

    sdds.close();
}

std::pair<double, double> IndexMap::getRange(const IndexMap::value_t::value_type &element,
                                             double position) const {
    double minDistance = std::numeric_limits<double>::max();
    std::pair<double, double> range(0.0, 0.0);
    const std::pair<invertedMap_t::const_iterator, invertedMap_t::const_iterator> its = mapElement2Range_m.equal_range(element);
    if (std::distance(its.first, its.second) == 0)
        throw OpalException("IndexMap::getRange()",
                            "Element \"" + element->getName() + "\" not registered");

    for (invertedMap_t::const_iterator it = its.first; it != its.second; ++ it) {
        double distance = std::min(std::abs((*it).second.first - position),
                                   std::abs((*it).second.second - position));
        if (distance < minDistance) {
            minDistance = distance;
            range = (*it).second;
        }
    }

    return range;
}

IndexMap::value_t IndexMap::getTouchingElements(const std::pair<double, double> &range) {
    map_t::iterator it = mapRange2Element_m.begin();
    const map_t::iterator end = mapRange2Element_m.end();
    IndexMap::value_t touchingElements;

    for (; it != end; ++ it) {
        if (almostEqual(it->first.first, range.second) ||
            almostEqual(it->first.second, range.first))
            touchingElements.insert((it->second).begin(), (it->second).end());
    }

    return touchingElements;
}

bool IndexMap::almostEqual(double x, double y) {
    return (std::abs(x - y) < std::numeric_limits<double>::epsilon() * std::abs(x + y) * 2 ||
            std::abs(x - y) < std::numeric_limits<double>::min());
}
