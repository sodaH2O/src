#ifndef CLASSIC_CCollimator_HH
#define CLASSIC_CCollimator_HH

// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CCollimator
//   Defines the abstract interface for a beam CCollimator.
//   *** MISSING *** CCollimator interface is still incomplete.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"

#include "gsl/gsl_spline.h"
#include "gsl/gsl_interp.h"

#include <vector>

class BeamlineVisitor;
class LossDataSink;

// Class CCollimator
// ------------------------------------------------------------------------
/// Abstract collimator.
//  Class CCollimator defines the abstract interface for a collimator.

class CCollimator: public Component {

public:

    /// Constructor with given name.
    explicit CCollimator(const std::string &name);

    CCollimator();
    CCollimator(const CCollimator &rhs);
    virtual ~CCollimator();

    /// Apply visitor to CCollimator.
    virtual void accept(BeamlineVisitor &) const;

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B);

    virtual bool checkCollimator(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep);

    virtual bool checkCollimator(Vector_t r, Vector_t rmin, Vector_t rmax);

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);

    virtual void initialise(PartBunchBase<double, 3> *bunch);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline(const double &kineticEnergy);

    virtual void goOffline();

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    void print();

    std::string  getCollimatorShape();
    void setOutputFN(std::string fn);
    std::string getOutputFN();

    unsigned int getLosses() const;

    // void setXsize(double a) ;

    // void setYsize(double b) ;

    // void setXpos(double x0) ;

    // void setYpos(double y0) ;

    // double getXsize(double a) ;

    // double getYsize(double b) ;

    // double getXpos() ;

    // double getYpos() ;

    // --------Cyclotron collimator

    void setXStart(double xstart) ;
    void setYStart(double ystart) ;
    void setZStart(double zstart) ;
    void setXEnd(double xend) ;
    void setYEnd(double yend) ;
    void setZEnd(double zend) ;
    void setWidth(double width) ;

    double getXStart() ;
    double getYStart() ;
    double getZStart() ;
    double getXEnd() ;
    double getYEnd() ;
    double getZEnd() ;
    double getWidth() ;

    int checkPoint(const double & x, const double & y );

private:

    // Not implemented.
    void operator=(const CCollimator &);

    std::string filename_m;               /**< The name of the outputfile*/

    bool informed_m;

    //parameters for CCollimator
    double xstart_m;
    double xend_m;
    double ystart_m;
    double yend_m;
    double zstart_m;
    double zend_m;
    double width_m;

    Point  geom_m[5];
    void setGeom();

    unsigned int losses_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    ParticleMatterInteractionHandler *parmatint_m;
};

inline
unsigned int CCollimator::getLosses() const {
    return losses_m;
}

inline
void CCollimator::setXStart(double xstart) {
    xstart_m = xstart;
}

inline
void CCollimator::setXEnd(double xend) {
    xend_m = xend;
}

inline
void CCollimator::setYStart(double ystart) {
    ystart_m = ystart;
}

inline
void CCollimator::setYEnd(double yend) {
    yend_m = yend;
}

inline
void CCollimator::setZStart(double zstart) {
    zstart_m = zstart;
}

inline
void CCollimator::setZEnd(double zend) {
    zend_m = zend;
}

inline
void CCollimator::setWidth(double width) {
    width_m = width;
}

inline
double CCollimator::getXStart() {
    return xstart_m;
}

inline
double CCollimator::getXEnd() {
    return xend_m;
}

inline
double CCollimator::getYStart() {
    return ystart_m;
}

inline
double CCollimator::getYEnd() {
    return yend_m;
}

inline
double CCollimator::getZStart() {
    return zstart_m;
}

inline
double CCollimator::getZEnd() {
    return zend_m;
}

inline
double CCollimator::getWidth() {
    return width_m;
}
#endif // CLASSIC_CCollimator_HH