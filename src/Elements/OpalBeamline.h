#ifndef OPAL_BEAMLINE_H
#define OPAL_BEAMLINE_H

#include <list>
#include <limits>
#include <vector>

#include "Algorithms/Tracker.h"
#include "Beamlines/Beamline.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Source.h"

#include "BasicActions/Option.h"
#include "Utilities/Options.h"
#include "Utilities/OpalSection.h"
#include "Utilities/ClassicField.h"

#include "Algorithms/CoordinateSystemTrafo.h"

class Tracker;
template <class T, unsigned Dim>
class PartBunchBase;
class ParticleMaterInteractionHandler;
class BoundaryGeometry;
class WakeFunction;
class Bend;

#define BEAMLINE_EOL  0x80000000   // end of line
#define BEAMLINE_OOB  0x40000000   // out of bounds
#define BEAMLINE_GEOM 0x30000000   // has geometry
#define BEAMLINE_WAKE 0x20000000   // has wake
#define BEAMLINE_BEND 0x10000000   // bends
#define BEAMLINE_PARTICLEMATTERINTERACTION 0x08000000 // has particle mater interaction

class OpalBeamline {

public:
    OpalBeamline();
    OpalBeamline(const Vector_t& origin,
                 const Quaternion& coordTrafoTo);
    ~OpalBeamline();

    OpalSection &getSectionAt(const Vector_t &, long &);
    OpalSection &getSection(const unsigned int &);

    void activateElements();
    std::set<std::shared_ptr<Component>> getElements(const Vector_t &x);
    Vector_t transformTo(const Vector_t &r) const;
    Vector_t transformFrom(const Vector_t &r) const;
    Vector_t rotateTo(const Vector_t &r) const;
    Vector_t rotateFrom(const Vector_t &r) const;

    Vector_t transformToLocalCS(const std::shared_ptr<Component> &comp,
                                const Vector_t &r) const;
    Vector_t transformFromLocalCS(const std::shared_ptr<Component> &comp,
                                  const Vector_t &r) const;
    Vector_t rotateToLocalCS(const std::shared_ptr<Component> &comp,
                             const Vector_t &r) const;
    Vector_t rotateFromLocalCS(const std::shared_ptr<Component> &comp,
                               const Vector_t &r) const;
    CoordinateSystemTrafo getCSTrafoLab2Local(const std::shared_ptr<Component> &comp) const;
    CoordinateSystemTrafo getCSTrafoLab2Local() const;
    CoordinateSystemTrafo getMisalignment(const std::shared_ptr<Component> &comp) const;
    // void getSectionIndexAt(const Vector_t &, long &) const;
    // double getSectionStart(const long &) const;
    // double getSectionEnd(const unsigned int &) const;
    // double getSectionEnd(const Vector_t &, long);

    double getStart(const Vector_t &) const;
    double getEnd(const Vector_t &) const;

    // void setOrientation(const Vector_t &, const Vector_t &);
    // void setOrientation(const Vector_t &, const unsigned int &);
    // void updateOrientation(const Vector_t &, const Vector_t &, const double &, const double &);

    // const Vector_t &getOrientation(const Vector_t &) const;
    // const Vector_t &getOrientation(const long &) const;

    // void resetStatus();
    // void setStatus(const unsigned int &, const bool &);
    // const bool &getStatus(const unsigned int &) const;

    void switchElements(const double &, const double &, const double &kineticEnergy, const bool &nomonitors = false);
    void switchAllElements();

    void switchElementsOff(const double &, ElementBase::ElementType eltype = ElementBase::ANY);
    void switchElementsOff();

    WakeFunction *getWakeFunction(const unsigned int &);
    std::shared_ptr<const ElementBase> getWakeFunctionOwner(const unsigned int &);

    ParticleMatterInteractionHandler *getParticleMatterInteractionHandler(const unsigned int &);

    BoundaryGeometry *getBoundaryGeometry(const unsigned int &);

    void getKFactors(const unsigned int &index, const Vector_t &pos, const long &sindex, const double &t, Vector_t &KR, Vector_t &KT);
    unsigned long getFieldAt(const unsigned int &, const Vector_t &, const long &, const double &, Vector_t &, Vector_t &);
    unsigned long getFieldAt(const Vector_t &, const Vector_t &, const double &, Vector_t &, Vector_t &);

    template<class T>
    void visit(const T &, BeamlineVisitor &, PartBunchBase<double, 3> *);

    void prepareSections();
    void compute3DLattice();
    void plot3DLattice();
    void save3DLattice();
    void save3DInput();
    void print(Inform &) const;

    FieldList getElementByType(ElementBase::ElementType);

    // need this for autophasing in case we have multiple tracks
    double calcBeamlineLength();

    void removeElement(const std::string &ElName);

    void swap(OpalBeamline & rhs);
    void merge(OpalBeamline &rhs);

    bool containsSource();
private:
    FieldList::iterator partiallyInsideDipole(const FieldList::iterator &it,
                                              const FieldList::iterator &begin,
                                              const FieldList::iterator &end,
                                              const unsigned int &minOrder);

    FieldList elements_m;
    bool prepared_m;
    bool containsSource_m;

    CoordinateSystemTrafo coordTransformationTo_m;

    static CompVec dummy_list_m;
    static OpalSection dummy_section_m;
};


// inline WakeFunction *OpalBeamline::getWakeFunction(const unsigned int &index) {
//     if(index < sections_m.size()) {
//         return sections_m[index].getWakeFunction();
//     }
//     return NULL;
// }

// inline std::shared_ptr<const ElementBase> OpalBeamline::getWakeFunctionOwner(const unsigned int &index) {
//     if(index < sections_m.size()) {
//         return sections_m[index].getWakeFunctionOwner();
//     }
//     return NULL;
// }

// inline BoundaryGeometry *OpalBeamline::getBoundaryGeometry(const unsigned int &index) {
//     if(index < sections_m.size()) {
//         return sections_m[index].getBoundaryGeometry();
//     }
//     return NULL;
// }

// inline ParticleMatterInteractionHandler *OpalBeamline::getParticleMatterInteractionHandler(const unsigned int &index) {
//     if(index < sections_m.size()) {
//         return sections_m[index].getParticleMatterInteractionHandler();
//     }
//     return 0;
// }

template<class T> inline
void OpalBeamline::visit(const T &element, BeamlineVisitor &, PartBunchBase<double, 3> *bunch) {
    Inform msg("OPAL ");
    double startField = 0.0;
    double endField = 0.0;
    std::shared_ptr<T> elptr(dynamic_cast<T *>(element.removeWrappers()->clone()));

    if (elptr->isElementPositionSet())
        startField = elptr->getElementPosition();

    elptr->initialise(bunch, startField, endField);
    elements_m.push_back(ClassicField(elptr, startField, endField));
}

template<> inline
void OpalBeamline::visit<Source>(const Source &element, BeamlineVisitor &, PartBunchBase<double, 3> *bunch) {
    containsSource_m = true;
}

template<> inline
void OpalBeamline::visit<AlignWrapper>(const AlignWrapper &wrap, BeamlineVisitor &visitor, PartBunchBase<double, 3> *) {
    wrap.getElement()->accept(visitor);
}

template<> inline
void OpalBeamline::visit<BeamBeam>(const BeamBeam &element, BeamlineVisitor &, PartBunchBase<double, 3> *) {
    WARNMSG(element.getTypeString() << " not implemented yet!" << endl);
}

template<> inline
void OpalBeamline::visit<Diagnostic>(const Diagnostic &element, BeamlineVisitor &, PartBunchBase<double, 3> *) {
    WARNMSG(element.getTypeString() << " not implemented yet!" << endl);
}

template<> inline
void OpalBeamline::visit<Lambertson>(const Lambertson &element, BeamlineVisitor &, PartBunchBase<double, 3> *) {
    WARNMSG(element.getTypeString() << " not implemented yet!" << endl);
}

template<> inline
void OpalBeamline::visit<Marker>(const Marker &element, BeamlineVisitor &, PartBunchBase<double, 3> *) {
}

template<> inline
void OpalBeamline::visit<RFQuadrupole>(const RFQuadrupole &element, BeamlineVisitor &, PartBunchBase<double, 3> *) {
    WARNMSG(element.getTypeString() << " not implemented yet!" << endl);
}

template<> inline
void OpalBeamline::visit<Separator>(const Separator &element, BeamlineVisitor &, PartBunchBase<double, 3> *) {
    WARNMSG(element.getTypeString() << " not implemented yet!" << endl);
}

template<> inline
void OpalBeamline::visit<Septum>(const Septum &element, BeamlineVisitor &, PartBunchBase<double, 3> *) {
    WARNMSG(element.getTypeString() << " not implemented yet!" << endl);
}

inline
void OpalBeamline::removeElement(const std::string &ElName) {
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
        if(flit->getElement()->getName() == ElName) {
            flit->setStart(-flit->getEnd());
            flit->setEnd(-flit->getEnd());
        }
    }
}

inline
Vector_t OpalBeamline::transformTo(const Vector_t &r) const {
    return coordTransformationTo_m.transformTo(r);
}

inline
Vector_t OpalBeamline::transformFrom(const Vector_t &r) const {
    return coordTransformationTo_m.transformFrom(r);
}

inline
Vector_t OpalBeamline::rotateTo(const Vector_t &r) const {
    return coordTransformationTo_m.rotateTo(r);
}

inline
Vector_t OpalBeamline::rotateFrom(const Vector_t &r) const {
    return coordTransformationTo_m.rotateFrom(r);
}

inline
Vector_t OpalBeamline::transformToLocalCS(const std::shared_ptr<Component> &comp,
                                          const Vector_t &r) const {
    return comp->getCSTrafoGlobal2Local().transformTo(r);
}

inline
Vector_t OpalBeamline::transformFromLocalCS(const std::shared_ptr<Component> &comp,
                                            const Vector_t &r) const {
    return comp->getCSTrafoGlobal2Local().transformFrom(r);
}

inline
Vector_t OpalBeamline::rotateToLocalCS(const std::shared_ptr<Component> &comp,
                                       const Vector_t &r) const {
    return comp->getCSTrafoGlobal2Local().rotateTo(r);
}

inline
Vector_t OpalBeamline::rotateFromLocalCS(const std::shared_ptr<Component> &comp,
                                         const Vector_t &r) const {
    return comp->getCSTrafoGlobal2Local().rotateFrom(r);
}

inline
CoordinateSystemTrafo OpalBeamline::getCSTrafoLab2Local(const std::shared_ptr<Component> &comp) const {
    return comp->getCSTrafoGlobal2Local();
}

inline
CoordinateSystemTrafo OpalBeamline::getCSTrafoLab2Local() const {
    return coordTransformationTo_m;
}

inline
CoordinateSystemTrafo OpalBeamline::getMisalignment(const std::shared_ptr<Component> &comp) const {
    return comp->getMisalignment();
}

inline
bool OpalBeamline::containsSource() {
    return containsSource_m;
}
#endif // OPAL_BEAMLINE_H
