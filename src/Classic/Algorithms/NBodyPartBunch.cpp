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
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_qrng.h>
#include <boost/format.hpp>


// Class NBodyPartBunch
// ------------------------------------------------------------------------

NBodyPartBunch::NBodyPartBunch(const PartData *ref):
    PartBunch(ref)
{
    std::cout << "\n!!!DAVID> NBodyPartBunch(const PartData) called\n"
	      << std::endl;
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
    computeSelfFields();
}

void NBodyPartBunch::computeSelfFields() {
    //Q: C
    //R: m
    //P: [gambet]
    //Ef: [V/m]
    //Bf: ?

#define TURN_ON_SELF_FIELDS
#ifndef TURN_ON_SELF_FIELDS
    std::cout << "self field off\n";
#else
    IpplTimings::startTimer(selfFieldTimer_m);
    /*std::cout << "\nDAVID> Reached: "
                << "NBodyPartBunch::computeSelfFields(){\n" << std::endl;*/
    
    // dump phase space
    /*for (size_t i = 0; i < Q.size(); i++) {
	std::cout << "R[" << i << "]: " << R[i] << std::endl;
	std::cout << "P[" << i << "]: " << P[i] << std::endl;
    }*/

    //get gamma for frame where <pz'> = 0
    Vector_t betaFrame = get_pmean()/get_gamma(); // betaFrame = <pz>/<gamma>
    betaFrame(0) = 0;
    betaFrame(1) = 0; // throw away x, y, component
    double gammaFrame = 1 / sqrt(1-dot(betaFrame, betaFrame));
    //std::cout << "betaFrame is " << betaFrame << std::endl;
    //std::cout << "gammaFrame is " << gammaFrame << std::endl;

    //static calc in lorentz frame
    for(size_t i = 0; i < R.size(); i++) {
	Ef[i] = 0;
	for(size_t j = 0; j < R.size(); j++) {
	    if (j == i){
		continue;
	    }
	    Vector_t r = R[i] - R[j];
	    r(2) *= gammaFrame;
	    r /= pow(euclidean_norm(r), 3);
	    r *= getCouplingConstant();
	    Ef[i] += r*Q[j];
        }
    }

    /*dump e field in lorentz frame
    for(size_t i = 0; i < R.size(); i++) {
	std::cout << "Ef0[" << i << "] is " << Ef[i] << std::endl;
    }*/

    //transform E, B to lab frame
    for (size_t i = 0; i < R.size(); i++) {
	// use original E0 to calc B
	Bf[i] = cross(betaFrame, Ef[i]) * gammaFrame / Physics::c;
	// transform E
	Ef[i](0) *= gammaFrame;
	Ef[i](1) *= gammaFrame;
    }
#endif

#define TURN_ON_DAMPING
#ifdef TURN_ON_DAMPING
    std::cout << "calculating damping e-field\n";
    // this only applys to electrons, as of right now
    // just an linear damping force on momentum
    double eEnergy = 510998.95; // eV
    double lightSpd = 299792458; // m/s
    double qElectron = -1; // e
    double charTime = 10*pow(10,-12); // sec, of damping to 1/e
    for(size_t i = 0; i < R.size(); i++) {
	Ef[i] -= eEnergy / (qElectron * lightSpd * charTime) * P[i];
    }
#endif

    /* dump E fields after transform
    for(size_t i = 0; i < R.size(); i++) {
	std::cout << "Ef[" << i << "] is " << Ef[i] << std::endl;
	std::cout << "Bf[" << i << "] is " << Bf[i] << std::endl;
    }*/
	
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
