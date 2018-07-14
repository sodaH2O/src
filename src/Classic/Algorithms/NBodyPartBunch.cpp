#include "Algorithms/NBodyPartBunch.h"
//let me just include everything else
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include <iostream>
#include <cfloat>
#include <fstream>
#include <iomanip>
#include <memory>
#include <utility>


#include "Distribution/Distribution.h"  // OPAL file
#include "Structure/FieldSolver.h"      // OPAL file
#include "Utilities/GeneralClassicException.h"
#include "Structure/LossDataSink.h"

#include "Algorithms/ListElem.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
//#include <gsl/gsl_sf_erf.h>
//#include <gsl/gsl_qrng.h>
#include <boost/format.hpp>


// Class NBodyPartBunch
// ------------------------------------------------------------------------

NBodyPartBunch::NBodyPartBunch(const PartData *ref):
    PartBunch(ref)
{
    std::cout << "\n!!!DAVID> NBodyPartBunch(const PartData) called\n"
	      << std::endl;
    std::cout << ref->getE() << std::endl;
    std::cout << P.size() << std::endl;
}

NBodyPartBunch::NBodyPartBunch(const std::vector<OpalParticle> &rhs,
                               const PartData *ref):
    PartBunch(rhs, ref)
{
    std::cout << "\n!!!DAVID> "
	      << "NBodyPartBunch(const std::vector<OpalParticles> &rhs, "
	      << "const PartData *ref) called\n"
              << std::endl;
}

NBodyPartBunch::NBodyPartBunch(const PartBunch &rhs):
    PartBunch(rhs)
{
    std::cout << "\n!!!DAVID> NBodyPartBunch(const ParfBunch &rhs) called\n"
	      << std::endl;
}

NBodyPartBunch::~NBodyPartBunch() {

}


void NBodyPartBunch::computeSelfFields(int binNumber) {
    std::cout << "\n!!!DAVID> not yet implemented: "
	      << "NBodyPartBunch::computeSelfFields(int binNumber) called\n"
              << std::endl;
}

void NBodyPartBunch::computeSelfFields() {
    IpplTimings::startTimer(selfFieldTimer_m);
    std::cout << "\n!!!DAVID>Yay, Reached: "
	      << "NBodyPartBunch::computeSelfFields(){\n" << std::endl;
    IpplTimings::stopTimer(selfFieldTimer_m);
}


void NBodyPartBunch::computeSelfFields_cycl(double gamma) {
    std::cout << "\n!!!DAVID> not yet implemented: "
	      << "NBodyPartBunch::computeSelfFields_cycl(double gamma) "
	      << "called\n"
              << std::endl;
}

void NBodyPartBunch::computeSelfFields_cycl(int bin) {
    std::cout << "\n!!!DAVID> not yet implemented: "
	      << "NBodyPartBunch::computeSelfFields_cycl(int bin) called\n"
              << std::endl;
}
