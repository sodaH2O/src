// ------------------------------------------------------------------------
// $RCSfile: ParallelTTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelTTracker
//   The visitor class for tracking particles with time as independent
//   variable.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/ParallelTTracker.h"

#include <cfloat>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>

#include "Algorithms/OrbitThreader.h"
#include "Algorithms/CavityAutophaser.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedBeamline.h"
#include "Lines/Sequence.h"

#include "Solvers/CSRWakeFunction.hh"

#include "AbstractObjects/OpalData.h"

#include "BasicActions/Option.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include "Distribution/Distribution.h"
#include "ValueDefinitions/RealVariable.h"
#include "Utilities/Timer.h"
#include "Utilities/OpalException.h"
#include "Solvers/ParticleMatterInteractionHandler.hh"
#include "Structure/BoundaryGeometry.h"
#include "Structure/LossDataSink.h"

class PartData;

using namespace std;

ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack),
    itsDataSink_m(NULL),
    itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getCoordTransformationTo()),
    RefPartR_m(0.0),
    RefPartP_m(0.0),
    globalEOL_m(false),
    wakeStatus_m(false),
    wakeFunction_m(NULL),
    pathLength_m(0.0),
    zstart_m(0.0),
    zStop_m(),
    dtCurrentTrack_m(0.0),
    dtAllTracks_m(),
    localTrackSteps_m(),
    minStepforReBin_m(-1),
    minBinEmitted_m(std::numeric_limits<size_t>::max()),
    repartFreq_m(-1),
    emissionSteps_m(std::numeric_limits<unsigned int>::max()),
    numParticlesInSimulation_m(0),
    timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
    timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
    fieldEvaluationTimer_m(IpplTimings::getTimer("External field eval")),
    BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart")),
    WakeFieldTimer_m(IpplTimings::getTimer("WakeField")),
    particleMaterStatus_m(false),
    totalParticlesInSimulation_m(0)
    // , logger_m("designPath_" + std::to_string(Ippl::myNode()) + ".dat")
{

    CoordinateSystemTrafo labToRef(beamline.getOrigin3D(),
                                   beamline.getCoordTransformationTo().conjugate());
    referenceToLabCSTrafo_m = labToRef.inverted();

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled)
        setupDKS();
#endif
}

ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   PartBunchBase<double, 3> *bunch,
                                   DataSink &ds,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack,
                                   const std::vector<unsigned long long> &maxSteps,
                                   double zstart,
                                   const std::vector<double> &zstop,
                                   const std::vector<double> &dt):
    Tracker(beamline, bunch, reference, revBeam, revTrack),
    itsDataSink_m(&ds),
    itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getCoordTransformationTo()),
    RefPartR_m(0.0),
    RefPartP_m(0.0),
    globalEOL_m(false),
    wakeStatus_m(false),
    wakeFunction_m(NULL),
    pathLength_m(0.0),
    zstart_m(zstart),
    dtCurrentTrack_m(0.0),
    minStepforReBin_m(-1),
    minBinEmitted_m(std::numeric_limits<size_t>::max()),
    repartFreq_m(-1),
    emissionSteps_m(numeric_limits<unsigned int>::max()),
    numParticlesInSimulation_m(0),
    timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
    timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
    fieldEvaluationTimer_m(IpplTimings::getTimer("External field eval")),
    BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart")),
    WakeFieldTimer_m(IpplTimings::getTimer("WakeField")),
    particleMaterStatus_m(false),
    totalParticlesInSimulation_m(0)
    // , logger_m("designPath_" + std::to_string(Ippl::myNode()) + ".dat")
{

    CoordinateSystemTrafo labToRef(beamline.getOrigin3D(),
                                   beamline.getCoordTransformationTo());
    referenceToLabCSTrafo_m = labToRef.inverted();

    for (std::vector<unsigned long long>::const_iterator it = maxSteps.begin(); it != maxSteps.end(); ++ it) {
        localTrackSteps_m.push(*it);
    }
    for (std::vector<double>::const_iterator it = dt.begin(); it != dt.end(); ++ it) {
        dtAllTracks_m.push(*it);
    }
    for (std::vector<double>::const_iterator it = zstop.begin(); it != zstop.end(); ++ it) {
        zStop_m.push(*it);
    }

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled)
        setupDKS();
#endif
}

ParallelTTracker::~ParallelTTracker() {
#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled)
        delete dksbase;
#endif
}

void ParallelTTracker::visitBeamline(const Beamline &bl) {
    const FlaggedBeamline* fbl = static_cast<const FlaggedBeamline*>(&bl);
    if (fbl->getRelativeFlag()) {
        OpalBeamline stash(fbl->getOrigin3D(), fbl->getCoordTransformationTo());
        stash.swap(itsOpalBeamline_m);
        fbl->iterate(*this, false);
        itsOpalBeamline_m.prepareSections();
        itsOpalBeamline_m.compute3DLattice();
        stash.merge(itsOpalBeamline_m);
        stash.swap(itsOpalBeamline_m);
    } else {
        fbl->iterate(*this, false);
    }
}

void ParallelTTracker::updateRFElement(std::string elName, double maxPhase) {
    FieldList cavities = itsOpalBeamline_m.getElementByType(ElementBase::RFCAVITY);
    FieldList travelingwaves = itsOpalBeamline_m.getElementByType(ElementBase::TRAVELINGWAVE);
    cavities.insert(cavities.end(), travelingwaves.begin(), travelingwaves.end());

    for (FieldList::iterator fit = cavities.begin(); fit != cavities.end(); ++ fit) {
        if ((*fit).getElement()->getName() == elName) {

            RFCavity *element = static_cast<RFCavity *>((*fit).getElement().get());

            element->setPhasem(maxPhase);
            element->setAutophaseVeto();

            INFOMSG("Restored cavity phase from the h5 file. Name: " << element->getName() << ", phase: " << maxPhase << " rad" << endl);
            return;
        }
    }
}

void ParallelTTracker::saveCavityPhases() {
    itsDataSink_m->storeCavityInformation();
}

void ParallelTTracker::restoreCavityPhases() {
    typedef std::vector<MaxPhasesT>::iterator iterator_t;

    if (OpalData::getInstance()->hasPriorTrack() ||
        OpalData::getInstance()->inRestartRun()) {
        iterator_t it = OpalData::getInstance()->getFirstMaxPhases();
        iterator_t end = OpalData::getInstance()->getLastMaxPhases();
        for (; it < end; ++ it) {
            updateRFElement((*it).first, (*it).second);
        }
    }
}

void ParallelTTracker::execute() {
    Inform msg("ParallelTTracker ", *gmsg);

    OpalData::getInstance()->setInPrepState(true);

    BorisPusher pusher(itsReference);
    const double globalTimeShift = itsBunch_m->weHaveEnergyBins()? OpalData::getInstance()->getGlobalPhaseShift(): 0.0;
    OpalData::getInstance()->setGlobalPhaseShift(0.0);

    dtCurrentTrack_m = itsBunch_m->getdT();

    evenlyDistributeParticles();

    if (OpalData::getInstance()->hasPriorTrack() || OpalData::getInstance()->inRestartRun()) {
        Options::openMode = Options::APPEND;
    }

    prepareSections();

    itsOpalBeamline_m.compute3DLattice();
    itsOpalBeamline_m.save3DLattice();
    itsOpalBeamline_m.save3DInput();

    std::queue<double> timeStepSizes(dtAllTracks_m);
    std::queue<unsigned long long> numSteps(localTrackSteps_m);
    double minTimeStep = timeStepSizes.front();
    unsigned long long totalNumSteps = 0;
    while (timeStepSizes.size() > 0) {
        if (minTimeStep > timeStepSizes.front()) {
            totalNumSteps = std::ceil(totalNumSteps * minTimeStep / timeStepSizes.front());
            minTimeStep = timeStepSizes.front();
        }
        totalNumSteps += std::ceil(numSteps.front() * timeStepSizes.front() / minTimeStep);

        numSteps.pop();
        timeStepSizes.pop();
    }

    itsOpalBeamline_m.activateElements();

    if ( OpalData::getInstance()->hasPriorTrack() ||
         OpalData::getInstance()->inRestartRun()) {

        referenceToLabCSTrafo_m = itsBunch_m->toLabTrafo_m;
        RefPartR_m = referenceToLabCSTrafo_m.transformFrom(itsBunch_m->RefPartR_m);
        RefPartP_m = referenceToLabCSTrafo_m.rotateFrom(itsBunch_m->RefPartP_m);

        pathLength_m = itsBunch_m->get_sPos();
        zstart_m = pathLength_m;

        restoreCavityPhases();
    } else {
        RefPartR_m = Vector_t(0.0);
        RefPartP_m = euclidean_norm(itsBunch_m->get_pmean_Distribution()) * Vector_t(0, 0, 1);

        if (itsBunch_m->getTotalNum() > 0) {
            if (!itsOpalBeamline_m.containsSource()) {
                RefPartP_m = itsReference.getP() / itsBunch_m->getM() * Vector_t(0, 0, 1);
            }

            if (zstart_m > pathLength_m) {
                findStartPosition(pusher);
            }

            itsBunch_m->set_sPos(pathLength_m);
        }
    }

    Vector_t rmin, rmax;
    itsBunch_m->get_bounds(rmin, rmax);
    OrbitThreader oth(itsReference,
                      referenceToLabCSTrafo_m.transformTo(RefPartR_m),
                      referenceToLabCSTrafo_m.rotateTo(RefPartP_m),
                      pathLength_m,
                      -rmin(2),
                      itsBunch_m->getT(),
                      minTimeStep,
                      totalNumSteps,
                      zStop_m.back() + 2 * rmax(2),
                      itsOpalBeamline_m);

    oth.execute();

    saveCavityPhases();

    numParticlesInSimulation_m = itsBunch_m->getTotalNum();
    totalParticlesInSimulation_m = itsBunch_m->getTotalNum();

    setTime();

    double t = itsBunch_m->getT() - globalTimeShift;
    itsBunch_m->setT(t);

    itsBunch_m->RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
    itsBunch_m->RefPartP_m = referenceToLabCSTrafo_m.rotateTo(RefPartP_m);

    *gmsg << level1 << *itsBunch_m << endl;

    unsigned long long step = itsBunch_m->getGlobalTrackStep();
    OPALTimer::Timer myt1;
    *gmsg << "Track start at: " << myt1.time() << ", t= " << Util::getTimeString(t) << "; "
          << "zstart at: " << Util::getLengthString(pathLength_m)
          << endl;

    prepareEmission();

    *gmsg << level1
          << "Executing ParallelTTracker, initial dt= " << Util::getTimeString(itsBunch_m->getdT()) << ";\n"
          << "max integration steps " << getMaxSteps(localTrackSteps_m) << ", next step= " << step << endl;

    setOptionalVariables();

    itsBunch_m->toLabTrafo_m = referenceToLabCSTrafo_m;

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled)
        allocateDeviceMemory();
#endif

    // loggingFrequency_m = floor(1e-11/itsBunch_m->getdT() + 0.5);
    globalEOL_m = false;
    wakeStatus_m = false;
    deletedParticles_m = false;
    OpalData::getInstance()->setInPrepState(false);
    while (localTrackSteps_m.size() > 0) {
        localTrackSteps_m.front() += step;
        dtCurrentTrack_m = dtAllTracks_m.front();
        changeDT();

        for (; step < localTrackSteps_m.front(); ++step) {

            timeIntegration1(pusher);

            itsBunch_m->Ef = Vector_t(0.0);
            itsBunch_m->Bf = Vector_t(0.0);

            computeSpaceChargeFields(step);

            selectDT();
            emitParticles(step);
            selectDT();

            computeExternalFields(oth);

            timeIntegration2(pusher);

            t += itsBunch_m->getdT();
            itsBunch_m->setT(t);

            if (t > 0.0) {
                updateRefToLabCSTrafo(pusher);
            }

            if (deletedParticles_m) {
                evenlyDistributeParticles();
                deletedParticles_m = false;
            }

            itsBunch_m->toLabTrafo_m = referenceToLabCSTrafo_m;
            itsBunch_m->RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
            itsBunch_m->RefPartP_m = referenceToLabCSTrafo_m.rotateTo(RefPartP_m);
            itsBunch_m->set_sPos(pathLength_m);

            if (hasEndOfLineReached()) break;

            bool const psDump = ((itsBunch_m->getGlobalTrackStep() % Options::psDumpFreq) + 1 == Options::psDumpFreq);
            bool const statDump = ((itsBunch_m->getGlobalTrackStep() % Options::statDumpFreq) + 1 == Options::statDumpFreq);
            dumpStats(step, psDump, statDump);

            itsBunch_m->incTrackSteps();

            double driftPerTimeStep = euclidean_norm(itsBunch_m->getdT() * Physics::c * RefPartP_m / Util::getGamma(RefPartP_m));
            if (std::abs(zStop_m.front() - pathLength_m) < 0.5 * driftPerTimeStep)
                localTrackSteps_m.front() = step;
        }

        if (globalEOL_m)
            break;

        dtAllTracks_m.pop();
        localTrackSteps_m.pop();
        zStop_m.pop();
    }

    itsBunch_m->toLabTrafo_m = referenceToLabCSTrafo_m;
    itsBunch_m->RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
    itsBunch_m->RefPartP_m = referenceToLabCSTrafo_m.rotateTo(RefPartP_m);
    itsBunch_m->set_sPos(pathLength_m);

    if (numParticlesInSimulation_m > minBinEmitted_m) {
        numParticlesInSimulation_m = itsBunch_m->getTotalNum();
    }

    bool const psDump = (((itsBunch_m->getGlobalTrackStep() - 1) % Options::psDumpFreq) + 1 != Options::psDumpFreq);
    bool const statDump = (((itsBunch_m->getGlobalTrackStep() - 1) % Options::statDumpFreq) + 1 != Options::statDumpFreq);
    writePhaseSpace((step + 1), psDump, statDump);

    msg << level2 << "Dump phase space of last step" << endl;

    itsOpalBeamline_m.switchElementsOff();

    OPALTimer::Timer myt3;
    *gmsg << "done executing ParallelTTracker at " << myt3.time() << endl;

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled)
        freeDeviceMemory();
#endif

    LossDataSink::writeStatistics();

    OpalData::getInstance()->setPriorTrack();
}

void ParallelTTracker::prepareSections() {

    itsBeamline_m.accept(*this);
    itsOpalBeamline_m.prepareSections();
}

void ParallelTTracker::timeIntegration1(BorisPusher & pusher) {

    IpplTimings::startTimer(timeIntegrationTimer1_m);
#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled)
        pushParticlesDKS();
    else
        pushParticles(pusher);
#else
    pushParticles(pusher);
#endif

    IpplTimings::stopTimer(timeIntegrationTimer1_m);
}


void ParallelTTracker::timeIntegration2(BorisPusher & pusher) {
    /*
      transport and emit particles
      that passed the cathode in the first
      half-step or that would pass it in the
      second half-step.

      to make IPPL and the field solver happy
      make sure that at least 10 particles are emitted

      also remember that node 0 has
      all the particles to be emitted

      this has to be done *after* the calculation of the
      space charges! thereby we neglect space charge effects
      in the very first step of a new-born particle.

    */

    IpplTimings::startTimer(timeIntegrationTimer2_m);
#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled) {
        kickParticlesDKS();
        pushParticlesDKS(false);
    } else {
        kickParticles(pusher);
        pushParticles(pusher);
    }
#else
    kickParticles(pusher);

    //switchElements();
    pushParticles(pusher);
#endif

    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++ i) {
        itsBunch_m->dt[i] = itsBunch_m->getdT();
    }

    IpplTimings::stopTimer(timeIntegrationTimer2_m);
}

void ParallelTTracker::selectDT() {

    if (itsBunch_m->getIfBeamEmitting()) {
        double dt = itsBunch_m->getEmissionDeltaT();
        itsBunch_m->setdT(dt);
    } else {
        double dt = dtCurrentTrack_m;
        itsBunch_m->setdT(dt);
    }
}

void ParallelTTracker::changeDT() {
    selectDT();
    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++ i) {
        itsBunch_m->dt[i] = itsBunch_m->getdT();
    }
}

void ParallelTTracker::emitParticles(long long step) {
    if (!itsBunch_m->weHaveEnergyBins()) return;

    if (itsBunch_m->getIfBeamEmitting()) {
        CoordinateSystemTrafo gunToLab = itsOpalBeamline_m.getCSTrafoLab2Local().inverted();
        CoordinateSystemTrafo refToGun = itsOpalBeamline_m.getCSTrafoLab2Local() * referenceToLabCSTrafo_m;

        transformBunch(refToGun);

        itsBunch_m->switchToUnitlessPositions(true);

        Vector_t externalE = Vector_t(0.0);
        Vector_t externalB = Vector_t(0.0);
        itsOpalBeamline_m.getFieldAt(gunToLab.transformTo(Vector_t(0.0)),
                                     gunToLab.rotateTo(Vector_t(0.0)),
                                     itsBunch_m->getT() + 0.5 * itsBunch_m->getdT(),
                                     externalE,
                                     externalB);
        itsBunch_m->emitParticles(externalE(2));

        itsBunch_m->updateNumTotal();
        numParticlesInSimulation_m = itsBunch_m->getTotalNum();
        itsBunch_m->switchOffUnitlessPositions(true);

        transformBunch(refToGun.inverted());
    }

    if (step > minStepforReBin_m) {
        itsBunch_m->Rebin();
        itsBunch_m->resetInterpolationCache(true);
    }
}


void ParallelTTracker::computeSpaceChargeFields(unsigned long long step) {
    if (numParticlesInSimulation_m <= minBinEmitted_m) return;

    if (!itsBunch_m->hasFieldSolver()) return;

    itsBunch_m->calcBeamParameters();
    Quaternion alignment = getQuaternion(itsBunch_m->get_pmean(), Vector_t(0, 0, 1));
    CoordinateSystemTrafo beamToReferenceCSTrafo(Vector_t(0, 0, pathLength_m), alignment.conjugate());
    CoordinateSystemTrafo referenceToBeamCSTrafo = beamToReferenceCSTrafo.inverted();

    const unsigned int localNum1 = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum1; ++ i) {
        itsBunch_m->R[i] = referenceToBeamCSTrafo.transformTo(itsBunch_m->R[i]);
        //I added the following
        itsBunch_m->P[i] = referenceToBeamCSTrafo.rotateTo(itsBunch_m->P[i]);
    }

    itsBunch_m->boundp();

    if (step % repartFreq_m + 1 == repartFreq_m) {
        doBinaryRepartition();
    }

    if (itsBunch_m->weHaveEnergyBins()) {
        itsBunch_m->calcGammas();
        itsBunch_m->resetInterpolationCache();
        ParticleAttrib<double> Q_back = itsBunch_m->Q;
        for (int binNumber = 0; binNumber < itsBunch_m->getLastEmittedEnergyBin() &&
                 binNumber < itsBunch_m->getNumberOfEnergyBins(); ++binNumber) {

            itsBunch_m->setBinCharge(binNumber);
            itsBunch_m->setGlobalMeanR(itsBunch_m->get_centroid());
            itsBunch_m->computeSelfFields(binNumber);
            itsBunch_m->Q = Q_back;

        }

    } else {
        itsBunch_m->setGlobalMeanR(itsBunch_m->get_centroid());
        itsBunch_m->computeSelfFields();
    }

    const unsigned int localNum2 = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum2; ++ i) {
        itsBunch_m->R[i] = beamToReferenceCSTrafo.transformTo(itsBunch_m->R[i]);
        itsBunch_m->Ef[i] = beamToReferenceCSTrafo.rotateTo(itsBunch_m->Ef[i]);
        itsBunch_m->Bf[i] = beamToReferenceCSTrafo.rotateTo(itsBunch_m->Bf[i]);
    }
}


void ParallelTTracker::computeExternalFields(OrbitThreader &oth) {
    IpplTimings::startTimer(fieldEvaluationTimer_m);
    Inform msg("ParallelTTracker ", *gmsg);
    const unsigned int localNum = itsBunch_m->getLocalNum();
    bool locPartOutOfBounds = false, globPartOutOfBounds = false;
    Vector_t rmin, rmax;
    itsBunch_m->get_bounds(rmin, rmax);
    IndexMap::value_t elements;

    try {
        elements = oth.query(pathLength_m + 0.5 * (rmax(2) + rmin(2)), rmax(2) - rmin(2));
    } catch(IndexMap::OutOfBounds &e) {
        globalEOL_m = true;
        return;
    }

    IndexMap::value_t::const_iterator it = elements.begin();
    const IndexMap::value_t::const_iterator end = elements.end();

    for (; it != end; ++ it) {
        CoordinateSystemTrafo refToLocalCSTrafo = (itsOpalBeamline_m.getMisalignment((*it)) *
                                                   (itsOpalBeamline_m.getCSTrafoLab2Local((*it)) * referenceToLabCSTrafo_m));
        CoordinateSystemTrafo localToRefCSTrafo = refToLocalCSTrafo.inverted();

        (*it)->setCurrentSCoordinate(pathLength_m + rmin(2));

        for (unsigned int i = 0; i < localNum; ++ i) {
            if (itsBunch_m->Bin[i] < 0) continue;

            itsBunch_m->R[i] = refToLocalCSTrafo.transformTo(itsBunch_m->R[i]);
            itsBunch_m->P[i] = refToLocalCSTrafo.rotateTo(itsBunch_m->P[i]);
            double &dt = itsBunch_m->dt[i];

            Vector_t localE(0.0), localB(0.0);

            if ((*it)->apply(i, itsBunch_m->getT() + 0.5 * dt, localE, localB)) {
                itsBunch_m->R[i] = localToRefCSTrafo.transformTo(itsBunch_m->R[i]);
                itsBunch_m->P[i] = localToRefCSTrafo.rotateTo(itsBunch_m->P[i]);
                itsBunch_m->Bin[i] = -1;
                locPartOutOfBounds = true;

                continue;
            }

            itsBunch_m->R[i] = localToRefCSTrafo.transformTo(itsBunch_m->R[i]);
            itsBunch_m->P[i] = localToRefCSTrafo.rotateTo(itsBunch_m->P[i]);
            itsBunch_m->Ef[i] += localToRefCSTrafo.rotateTo(localE);
            itsBunch_m->Bf[i] += localToRefCSTrafo.rotateTo(localB);
        }
    }

    IpplTimings::stopTimer(fieldEvaluationTimer_m);

    computeWakefield(elements);
    computeParticleMatterInteraction(elements, oth);

    reduce(locPartOutOfBounds, globPartOutOfBounds, OpOrAssign());

    size_t ne = 0;
    if (globPartOutOfBounds) {
        if (itsBunch_m->hasFieldSolver()) {
            ne = itsBunch_m->boundp_destroyT();
        } else {
            ne = itsBunch_m->destroyT();
        }
        numParticlesInSimulation_m  = itsBunch_m->getTotalNum();
        totalParticlesInSimulation_m -= ne;
        deletedParticles_m = true;
    }

    size_t totalNum = itsBunch_m->getTotalNum();
    if (numParticlesInSimulation_m > minBinEmitted_m || totalNum > minBinEmitted_m) {
        numParticlesInSimulation_m = totalNum;
    }

    if (ne > 0) {
        msg << level1 << "* Deleted " << ne << " particles, "
            << "remaining " << itsBunch_m->getTotalNum() << " particles" << endl;
    }
}

void ParallelTTracker::computeWakefield(IndexMap::value_t &elements) {
    bool hasWake = false;
    WakeFunction *wfInstance;

    Inform msg("ParallelTTracker ", *gmsg);

    const unsigned int localNum = itsBunch_m->getLocalNum();

    IndexMap::value_t::const_iterator it = elements.begin();
    const IndexMap::value_t::const_iterator end = elements.end();

    for (; it != end; ++ it) {
        if ((*it)->hasWake() && !hasWake) {
            IpplTimings::startTimer(WakeFieldTimer_m);

            hasWake = true;

            if ((*it)->getWake()->getType() == "CSRWakeFunction" ||
                (*it)->getWake()->getType() == "CSRIGFWakeFunction") {
                if ((*it)->getType() == ElementBase::RBEND ||
                    (*it)->getType() == ElementBase::SBEND) {
                    wfInstance = (*it)->getWake();
                    wakeFunction_m = wfInstance;
                } else {
                    wfInstance = wakeFunction_m;
                }
            } else {
                wfInstance = (*it)->getWake();
            }

            if (!wfInstance) {
                throw OpalException("ParallelTTracker::computeWakefield",
                                    "empty wake function");
            }

            if (!wakeStatus_m) {
                msg << level2 << "============== START WAKE CALCULATION =============" << endl;
                wfInstance->initialize((*it).get());
                wakeStatus_m = true;
            }

            if (!itsBunch_m->hasFieldSolver()) itsBunch_m->calcBeamParameters();

            Quaternion alignment = getQuaternion(itsBunch_m->get_pmean(), Vector_t(0, 0, 1));
            CoordinateSystemTrafo referenceToBeamCSTrafo(Vector_t(0.0), alignment);
            CoordinateSystemTrafo beamToReferenceCSTrafo = referenceToBeamCSTrafo.inverted();

            for (unsigned int i = 0; i < localNum; ++ i) {
                itsBunch_m->R[i] = referenceToBeamCSTrafo.transformTo(itsBunch_m->R[i]);
                itsBunch_m->P[i] = referenceToBeamCSTrafo.rotateTo(itsBunch_m->P[i]);
                itsBunch_m->Ef[i] = referenceToBeamCSTrafo.rotateTo(itsBunch_m->Ef[i]);
            }

            wfInstance->apply(itsBunch_m);

            for (unsigned int i = 0; i < localNum; ++ i) {
                itsBunch_m->R[i] = beamToReferenceCSTrafo.transformTo(itsBunch_m->R[i]);
                itsBunch_m->P[i] = beamToReferenceCSTrafo.rotateTo(itsBunch_m->P[i]);
                itsBunch_m->Ef[i] = beamToReferenceCSTrafo.rotateTo(itsBunch_m->Ef[i]);
            }

            IpplTimings::stopTimer(WakeFieldTimer_m);
        }
    }

    if (wakeStatus_m && !hasWake) {
        msg << level2 << "=============== END WAKE CALCULATION ==============" << endl;
        wakeStatus_m = false;
    }
}

void ParallelTTracker::computeParticleMatterInteraction(IndexMap::value_t elements, OrbitThreader &oth) {
    Inform msg("ParallelTTracker ", *gmsg);
    std::set<IndexMap::value_t::value_type> elementsWithParticleMatterInteraction;
    std::set<ParticleMatterInteractionHandler*> particleMatterinteractionHandlers;
    std::pair<double, double> currentRange(0.0, 0.0);

    while (elements.size() > 0) {
        auto it = elements.begin();
        if ((*it)->hasParticleMatterInteraction()) {
            elementsWithParticleMatterInteraction.insert(*it);
            particleMatterinteractionHandlers.insert((*it)->getParticleMatterInteraction());

            std::pair<double, double> range = oth.getRange(*it, pathLength_m);
            currentRange.first = std::min(currentRange.first, range.first);
            currentRange.second = std::max(currentRange.second, range.second);

            IndexMap::value_t touching = oth.getTouchingElements(range);
            elements.insert(touching.begin(), touching.end());
        }

        elements.erase(it);
    }

    if (elementsWithParticleMatterInteraction.size() > 0) {
        std::set<ParticleMatterInteractionHandler*> oldSPHandlers;
        std::vector<ParticleMatterInteractionHandler*> leftBehindSPHandlers, newSPHandlers;
        for (auto it: activeParticleMatterInteractionHandlers_m) {
            oldSPHandlers.insert(it);
        }

        leftBehindSPHandlers.resize(std::max(oldSPHandlers.size(),
                                             particleMatterinteractionHandlers.size()));
        auto last = std::set_difference(oldSPHandlers.begin(), oldSPHandlers.end(),
                                        particleMatterinteractionHandlers.begin(), particleMatterinteractionHandlers.end(),
                                        leftBehindSPHandlers.begin());
        leftBehindSPHandlers.resize(last - leftBehindSPHandlers.begin());

        for (auto it: leftBehindSPHandlers) {
            if (!it->stillActive()) {
                activeParticleMatterInteractionHandlers_m.erase(it);
            }
        }

        newSPHandlers.resize(std::max(oldSPHandlers.size(),
                                      elementsWithParticleMatterInteraction.size()));
        last = std::set_difference(particleMatterinteractionHandlers.begin(), particleMatterinteractionHandlers.end(),
                                   oldSPHandlers.begin(), oldSPHandlers.end(),
                                   newSPHandlers.begin());
        newSPHandlers.resize(last - newSPHandlers.begin());

        for (auto it: newSPHandlers) {
            activeParticleMatterInteractionHandlers_m.insert(it);
        }

        if(!particleMaterStatus_m) {
            msg << level2 << "============== START PARTICLE MATER INTERACTION CALCULATION =============" << endl;
            particleMaterStatus_m = true;
        }
    }

    if (particleMaterStatus_m) {
        do {
            ///all particles in material if max per node is 2 and other degraders have 0 particles
            //check if more than one degrader has particles
            ParticleMatterInteractionHandler* onlyDegraderWithParticles = NULL;
            int degradersWithParticlesCount = 0;
            for (auto it: activeParticleMatterInteractionHandlers_m) {
                it->setFlagAllParticlesIn(false);
                if (it->getParticlesInMat() > 0) {
                    onlyDegraderWithParticles = it;
                    ++ degradersWithParticlesCount;
                }
            }

            //if max particles per node is 2, and only one degrader has particles set
            //AllParticlesIn for this degrader to true
            int maxPerNode = itsBunch_m->getLocalNum();
            reduce(maxPerNode, maxPerNode, OpMaxAssign());
            bool allParticlesInMat = (maxPerNode == 0 &&
                                      degradersWithParticlesCount == 1);

            if (allParticlesInMat) {
                onlyDegraderWithParticles->setFlagAllParticlesIn(true);
                msg << "All particles in degrader" << endl;
            }

            auto boundingSphere = itsBunch_m->getLocalBoundingSphere();
            unsigned redifusedParticles = 0;
            for (auto it: activeParticleMatterInteractionHandlers_m) {
                ElementBase* element = it->getElement();
                CoordinateSystemTrafo refToLocalCSTrafo = (element->getMisalignment() *
                                                           (element->getCSTrafoGlobal2Local() * referenceToLabCSTrafo_m));
                CoordinateSystemTrafo localToRefCSTrafo = refToLocalCSTrafo.inverted();

                const unsigned int localNum = itsBunch_m->getLocalNum();
                for (unsigned int i = 0; i < localNum; ++i) {
                    itsBunch_m->R[i] = refToLocalCSTrafo.transformTo(itsBunch_m->R[i]);
                    itsBunch_m->P[i] = refToLocalCSTrafo.rotateTo(itsBunch_m->P[i]);
                }
                boundingSphere.first = refToLocalCSTrafo.transformTo(boundingSphere.first);

                it->apply(itsBunch_m, boundingSphere, totalParticlesInSimulation_m);
                it->print(msg);

                boundingSphere.first = localToRefCSTrafo.transformTo(boundingSphere.first);

                const unsigned int newLocalNum = itsBunch_m->getLocalNum();
                for (unsigned int i = 0; i < newLocalNum; ++i) {
                    itsBunch_m->R[i] = localToRefCSTrafo.transformTo(itsBunch_m->R[i]);
                    itsBunch_m->P[i] = localToRefCSTrafo.rotateTo(itsBunch_m->P[i]);
                }

                redifusedParticles += it->getRediffused();
                //if all particles where in material update time to time in degrader
                if (it->getFlagAllParticlesIn()) {
                    double timeDifference = it->getTime() - itsBunch_m->getdT() - itsBunch_m->getT();
                    if (timeDifference > 0.0) {
                        const unsigned int numSteps = std::ceil(timeDifference / itsBunch_m->getdT());
                        const double origdT = itsBunch_m->getdT();
                        itsBunch_m->setdT(timeDifference / numSteps);
                        BorisPusher pusher(itsReference);
                        for (unsigned int i = 0; i < numSteps; ++ i) {
                            updateReferenceParticle(pusher);
                        }
                        itsBunch_m->setdT(origdT);
                    }
                    itsBunch_m->setT(it->getTime() - itsBunch_m->getdT());
                }
            }

            //perform boundp only if there are particles in the bunch, or there are particles
            //comming out of the degrader
            if (maxPerNode > 0 || redifusedParticles > 0) {
                if (itsBunch_m->hasFieldSolver()) {
                    itsBunch_m->boundp();
                } else {
                    itsBunch_m->updateNumTotal();
                }
            }

            //if bunch has no particles update time to time in degrader
            if (itsBunch_m->getTotalNum() == 0)
                itsBunch_m->setT(itsBunch_m->getT() + itsBunch_m->getdT());

        } while (itsBunch_m->getTotalNum() == 0);


        if (activeParticleMatterInteractionHandlers_m.size() == 0) {
            msg << level2 << "============== END PARTICLE MATER INTERACTION CALCULATION =============" << endl;
            particleMaterStatus_m = false;
        }
    }
}

void ParallelTTracker::doBinaryRepartition() {
    size_t particles_or_bins = std::max(minBinEmitted_m, size_t(1000));
    if (itsBunch_m->hasFieldSolver() && numParticlesInSimulation_m > particles_or_bins) {

        INFOMSG("*****************************************************************" << endl);
        INFOMSG("do repartition because of repartFreq_m" << endl);
        INFOMSG("*****************************************************************" << endl);
        IpplTimings::startTimer(BinRepartTimer_m);
        itsBunch_m->do_binaryRepart();
        Ippl::Comm->barrier();
        IpplTimings::stopTimer(BinRepartTimer_m);
        INFOMSG("*****************************************************************" << endl);
        INFOMSG("do repartition done" << endl);
        INFOMSG("*****************************************************************" << endl);
    }
}

void ParallelTTracker::dumpStats(long long step, bool psDump, bool statDump) {
    OPALTimer::Timer myt2;
    Inform msg("ParallelTTracker ", *gmsg);

    if (itsBunch_m->getGlobalTrackStep() % 1000 + 1 == 1000) {
        msg << level1;
    } else if (itsBunch_m->getGlobalTrackStep() % 100 + 1 == 100) {
        msg << level2;
    } else {
        msg << level3;
    }

    if (numParticlesInSimulation_m == 0) {
        msg << myt2.time() << " "
            << "Step " << setw(6) <<  itsBunch_m->getGlobalTrackStep() << "; "
            << "   -- no emission yet --     "
            << "t= "   << Util::getTimeString(itsBunch_m->getT())
            << endl;
        return;
    }

    itsBunch_m->calcEMean();
    size_t totalParticles_f = numParticlesInSimulation_m;
    if (totalParticles_f <= minBinEmitted_m) {
        msg << myt2.time() << " "
            << "Step " << setw(6) << itsBunch_m->getGlobalTrackStep() << "; "
            << "only " << setw(4) << totalParticles_f << " particles emitted; "
            << "t= "   << Util::getTimeString(itsBunch_m->getT()) << ", "
            << "E="    << Util::getEnergyString(itsBunch_m->get_meanKineticEnergy()) << ", "
            << endl;
    } else if (std::isnan(pathLength_m) || std::isinf(pathLength_m)) {
        throw OpalException("ParallelTTracker::dumpStats()",
                            "there seems to be something wrong with the position of the bunch!");
    } else {

        msg << myt2.time() << " "
            << "Step " << setw(6)  <<  itsBunch_m->getGlobalTrackStep() << " "
            << "at "   << Util::getLengthString(pathLength_m) << ", "
            << "t= "   << Util::getTimeString(itsBunch_m->getT()) << ", "
            << "E="    << Util::getEnergyString(itsBunch_m->get_meanKineticEnergy())
            << endl;

        writePhaseSpace(step, psDump, statDump);
    }
}


void ParallelTTracker::setOptionalVariables() {
    Inform msg("ParallelTTracker ", *gmsg);

    minBinEmitted_m  = 10;
    RealVariable *ar = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINBINEMITTED"));
    if (ar)
        minBinEmitted_m = static_cast<size_t>(ar->getReal());
    msg << level2 << "MINBINEMITTED " << minBinEmitted_m << endl;

    minStepforReBin_m  = 200;
    RealVariable *br = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINSTEPFORREBIN"));
    if (br)
        minStepforReBin_m = static_cast<int>(br->getReal());
    msg << level2 << "MINSTEPFORREBIN " << minStepforReBin_m << endl;

    // there is no point to do repartitioning with one node
    if (Ippl::getNodes() == 1) {
        repartFreq_m = numeric_limits<unsigned int>::max();
    } else {
        repartFreq_m = 1000;
        RealVariable *rep = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("REPARTFREQ"));
        if (rep)
            repartFreq_m = static_cast<int>(rep->getReal());
        msg << level2 << "REPARTFREQ " << repartFreq_m << endl;
    }
}


bool ParallelTTracker::hasEndOfLineReached() {
    reduce(&globalEOL_m, &globalEOL_m + 1, &globalEOL_m, OpBitwiseAndAssign());
    return globalEOL_m;
}

void ParallelTTracker::setTime() {
    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++i) {
        itsBunch_m->dt[i] = itsBunch_m->getdT();
    }
}

void ParallelTTracker::prepareEmission() {
    Inform msg("ParallelTTracker ", *gmsg);

    if (!itsBunch_m->doEmission()) return;

    emissionSteps_m = itsBunch_m->getNumberOfEmissionSteps();
    msg << level2 << "Do emission for " << itsBunch_m->getTEmission() << " [s] using "
        << itsBunch_m->getNumberOfEnergyBins() << " energy bins " << endl
        << "Change dT from " <<  itsBunch_m->getdT() << " [s] to "
        <<  itsBunch_m->getEmissionDeltaT() << " [s] during emission " << endl;;

}

void ParallelTTracker::writePhaseSpace(const long long step, bool psDump, bool statDump) {
    extern Inform *gmsg;
    Inform msg("OPAL ", *gmsg);
    Vector_t externalE, externalB;
    Vector_t FDext[2];  // FDext = {BHead, EHead, BRef, ERef, BTail, ETail}.

    // Sample fields at (xmin, ymin, zmin), (xmax, ymax, zmax) and the centroid location. We
    // are sampling the electric and magnetic fields at the back, front and
    // center of the beam.
    Vector_t rmin, rmax;
    itsBunch_m->get_bounds(rmin, rmax);

    if (psDump || statDump) {
        externalB = Vector_t(0.0);
        externalE = Vector_t(0.0);
        itsOpalBeamline_m.getFieldAt(referenceToLabCSTrafo_m.transformTo(RefPartR_m),
                                     referenceToLabCSTrafo_m.rotateTo(RefPartP_m),
                                     itsBunch_m->getT() - 0.5 * itsBunch_m->getdT(),
                                     externalE,
                                     externalB);
        FDext[0] = referenceToLabCSTrafo_m.rotateFrom(externalB);
        FDext[1] = referenceToLabCSTrafo_m.rotateFrom(externalE * 1e-6);
    }

    if (statDump) {
        std::vector<std::pair<std::string, unsigned int> > collimatorLosses;
        FieldList collimators = itsOpalBeamline_m.getElementByType(ElementBase::CCOLLIMATOR);
	if (collimators.size() != 0) {
            for (FieldList::iterator it = collimators.begin(); it != collimators.end(); ++ it) {
                FlexibleCollimator* coll = static_cast<FlexibleCollimator*>(it->getElement().get());
                std::string name = coll->getName();
                unsigned int losses = coll->getLosses();
                collimatorLosses.push_back(std::make_pair(name, losses));
            }
            std::sort(collimatorLosses.begin(), collimatorLosses.end(),
                      [](const std::pair<std::string, unsigned int>& a, const std::pair<std::string, unsigned int>& b) ->bool {
                          return a.first < b.first;
                      });
            std::vector<unsigned int> bareLosses(collimatorLosses.size(),0);
            for (size_t i = 0; i < collimatorLosses.size(); ++ i){
                bareLosses[i] = collimatorLosses[i].second;
            }

            reduce(&bareLosses[0], &bareLosses[0] + bareLosses.size(), &bareLosses[0], OpAddAssign());

            for (size_t i = 0; i < collimatorLosses.size(); ++ i){
                collimatorLosses[i].second = bareLosses[i];
            }
	}
        // Write statistical data.
        itsDataSink_m->writeStatData(itsBunch_m, FDext, collimatorLosses);

        msg << level3 << "* Wrote beam statistics." << endl;
    }

    if (psDump && (itsBunch_m->getTotalNum() > 0)) {
        // Write fields to .h5 file.
        const size_t localNum = itsBunch_m->getLocalNum();
        double distToLastStop = zStop_m.back() - pathLength_m;
        Vector_t driftPerTimeStep = itsBunch_m->getdT() * Physics::c * RefPartP_m / Util::getGamma(RefPartP_m);
        bool driftToCorrectPosition = std::abs(distToLastStop) < 0.5 * euclidean_norm(driftPerTimeStep);
        Ppos_t stashedR;

        if (driftToCorrectPosition) {
            const double tau = distToLastStop / euclidean_norm(driftPerTimeStep) * itsBunch_m->getdT();
            if (localNum > 0) {
                stashedR.create(localNum);
                stashedR = itsBunch_m->R;

                for (size_t i = 0; i < localNum; ++ i) {
                    itsBunch_m->R[i] += tau * (Physics::c * itsBunch_m->P[i] / Util::getGamma(itsBunch_m->P[i]) -
                                               driftPerTimeStep / itsBunch_m->getdT());
                }
            }

            itsBunch_m->RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m + tau * driftPerTimeStep / itsBunch_m->getdT());
            CoordinateSystemTrafo update(tau * driftPerTimeStep / itsBunch_m->getdT(),
                                         Quaternion());
            itsBunch_m->toLabTrafo_m = referenceToLabCSTrafo_m * update.inverted();

            itsBunch_m->set_sPos(zStop_m.back());

            itsBunch_m->calcBeamParameters();
        }
        if (!statDump && !driftToCorrectPosition) itsBunch_m->calcBeamParameters();

        msg << *itsBunch_m << endl;
        itsDataSink_m->writePhaseSpace(itsBunch_m, FDext);

        if (driftToCorrectPosition) {
            if (localNum > 0) {
                itsBunch_m->R = stashedR;
                stashedR.destroy(localNum, 0);
            }

            itsBunch_m->RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
            itsBunch_m->set_sPos(pathLength_m);

            itsBunch_m->calcBeamParameters();
        }

        msg << level2 << "* Wrote beam phase space." << endl;
    }
}

void ParallelTTracker::updateReferenceParticle(const BorisPusher &pusher) {
    //static size_t step = 0;
    const double dt = std::min(itsBunch_m->getT(), itsBunch_m->getdT());
    const double scaleFactor = Physics::c * dt;
    Vector_t Ef(0.0), Bf(0.0);
    // Vector_t oldR = RefPartR_m;

    RefPartR_m /= scaleFactor;
    pusher.push(RefPartR_m, RefPartP_m, dt);
    RefPartR_m *= scaleFactor;

    IndexMap::value_t elements = itsOpalBeamline_m.getElements(referenceToLabCSTrafo_m.transformTo(RefPartR_m));
    IndexMap::value_t::const_iterator it = elements.begin();
    const IndexMap::value_t::const_iterator end = elements.end();

    for (; it != end; ++ it) {
        CoordinateSystemTrafo refToLocalCSTrafo = itsOpalBeamline_m.getCSTrafoLab2Local((*it)) * referenceToLabCSTrafo_m;

        Vector_t localR = refToLocalCSTrafo.transformTo(RefPartR_m);
        Vector_t localP = refToLocalCSTrafo.rotateTo(RefPartP_m);
        Vector_t localE(0.0), localB(0.0);

        if ((*it)->applyToReferenceParticle(localR,
                                            localP,
                                            itsBunch_m->getT() - 0.5 * dt,
                                            localE,
                                            localB)) {
            *gmsg << level1 << "The reference particle hit an element" << endl;
            globalEOL_m = true;
        }

        Ef += refToLocalCSTrafo.rotateFrom(localE);
        Bf += refToLocalCSTrafo.rotateFrom(localB);
    }

    // if (step % loggingFrequency_m == 0) {
    //     Vector_t R = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
    //     Vector_t P = referenceToLabCSTrafo_m.rotateTo(RefPartP_m);
    //     Vector_t lEf = referenceToLabCSTrafo_m.rotateTo(Ef);
    //     Vector_t lBf = referenceToLabCSTrafo_m.rotateTo(Bf);
    //     logger_m << std::setw(18) << std::setprecision(8) << pathLength_m + euclidean_norm(RefPartR_m - oldR)
    //              << std::setw(18) << std::setprecision(8) << R(0)
    //              << std::setw(18) << std::setprecision(8) << R(1)
    //              << std::setw(18) << std::setprecision(8) << R(2)
    //              << std::setw(18) << std::setprecision(8) << P(0)
    //              << std::setw(18) << std::setprecision(8) << P(1)
    //              << std::setw(18) << std::setprecision(8) << P(2)
    //              << std::setw(18) << std::setprecision(8) << lEf(0)
    //              << std::setw(18) << std::setprecision(8) << lEf(1)
    //              << std::setw(18) << std::setprecision(8) << lEf(2)
    //              << std::setw(18) << std::setprecision(8) << lBf(0)
    //              << std::setw(18) << std::setprecision(8) << lBf(1)
    //              << std::setw(18) << std::setprecision(8) << lBf(2)
    //              << std::setw(18) << std::setprecision(8) << Util::getEnergy(RefPartP_m, itsBunch_m->getM() * 1e-6)
    //              << std::setw(18) << std::setprecision(8) << (itsBunch_m->getT() - 0.5 * itsBunch_m->getdT()) * 1e9
    //              << std::endl;
    // }

    pusher.kick(RefPartR_m, RefPartP_m, Ef, Bf, dt);

    RefPartR_m /= scaleFactor;
    pusher.push(RefPartR_m, RefPartP_m, dt);
    RefPartR_m *= scaleFactor;
    //++ step;
}

void ParallelTTracker::transformBunch(const CoordinateSystemTrafo &trafo) {
    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++i) {
        itsBunch_m->R[i] = trafo.transformTo(itsBunch_m->R[i]);
        itsBunch_m->P[i] = trafo.rotateTo(itsBunch_m->P[i]);
        itsBunch_m->Ef[i] = trafo.rotateTo(itsBunch_m->Ef[i]);
        itsBunch_m->Bf[i] = trafo.rotateTo(itsBunch_m->Bf[i]);
    }
}

void ParallelTTracker::updateRefToLabCSTrafo(const BorisPusher &pusher) {
    updateReferenceParticle(pusher);

    pathLength_m += euclidean_norm(RefPartR_m);

    CoordinateSystemTrafo update(RefPartR_m,
                                 getQuaternion(RefPartP_m, Vector_t(0, 0, 1)));

    RefPartR_m = update.transformTo(RefPartR_m);
    RefPartP_m = update.rotateTo(RefPartP_m);

    transformBunch(update);

    referenceToLabCSTrafo_m = referenceToLabCSTrafo_m * update.inverted();
}

void ParallelTTracker::findStartPosition(const BorisPusher &pusher) {

    double t = 0.0;
    itsBunch_m->setT(t);

    dtCurrentTrack_m = dtAllTracks_m.front();
    changeDT();

    if (Util::getEnergy(RefPartP_m, itsBunch_m->getM()) < 1e-3) {
        double gamma = 0.1 / itsBunch_m->getM() + 1.0;
        RefPartP_m = sqrt(std::pow(gamma, 2) - 1) * Vector_t(0, 0, 1);
    }

    while (true) {
        autophaseCavities(pusher);

        t += itsBunch_m->getdT();
        itsBunch_m->setT(t);

        Vector_t oldR = RefPartR_m;
        updateReferenceParticle(pusher);
        pathLength_m += euclidean_norm(RefPartR_m - oldR);

        if (pathLength_m > zStop_m.front()) {
            if (localTrackSteps_m.size() == 0) return;

            dtAllTracks_m.pop();
            localTrackSteps_m.pop();
            zStop_m.pop();

            changeDT();
        }

        double speed = euclidean_norm(RefPartP_m) * Physics::c / sqrt(dot(RefPartP_m, RefPartP_m) + 1);
        if (std::abs(pathLength_m - zstart_m) <=  0.5 * itsBunch_m->getdT() * speed) {
            double tau = (pathLength_m - zstart_m) / speed;

            t += tau;
            itsBunch_m->setT(t);

            RefPartR_m /= (Physics::c * tau);
            pusher.push(RefPartR_m, RefPartP_m, tau);
            RefPartR_m *= (Physics::c * tau);

            pathLength_m = zstart_m;

            CoordinateSystemTrafo update(RefPartR_m,
                                         getQuaternion(RefPartP_m, Vector_t(0, 0, 1)));
            referenceToLabCSTrafo_m = referenceToLabCSTrafo_m * update.inverted();

            RefPartR_m = update.transformTo(RefPartR_m);
            RefPartP_m = update.rotateTo(RefPartP_m);

            return;
        }
    }
}

void ParallelTTracker::autophaseCavities(const BorisPusher &pusher) {

    double t = itsBunch_m->getT();
    Vector_t nextR = RefPartR_m / (Physics::c * itsBunch_m->getdT());
    pusher.push(nextR, RefPartP_m, itsBunch_m->getdT());
    nextR *= Physics::c * itsBunch_m->getdT();

    auto elementSet = itsOpalBeamline_m.getElements(referenceToLabCSTrafo_m.transformTo(nextR));
    for (auto element: elementSet) {
        if (element->getType() == ElementBase::TRAVELINGWAVE) {
            const TravelingWave *TWelement = static_cast<const TravelingWave *>(element.get());
            if (!TWelement->getAutophaseVeto()) {
                RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateTo(RefPartP_m);
                CavityAutophaser ap(itsReference, element);
                ap.getPhaseAtMaxEnergy(itsOpalBeamline_m.transformToLocalCS(element, RefPartR_m),
                                       itsOpalBeamline_m.rotateToLocalCS(element, RefPartP_m),
                                       t, itsBunch_m->getdT());
                RefPartR_m = referenceToLabCSTrafo_m.transformFrom(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateFrom(RefPartP_m);
            }

        } else if (element->getType() == ElementBase::RFCAVITY) {
            const RFCavity *RFelement = static_cast<const RFCavity *>(element.get());
            if (!RFelement->getAutophaseVeto()) {
                RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateTo(RefPartP_m);
                CavityAutophaser ap(itsReference, element);
                ap.getPhaseAtMaxEnergy(itsOpalBeamline_m.transformToLocalCS(element, RefPartR_m),
                                       itsOpalBeamline_m.rotateToLocalCS(element, RefPartP_m),
                                       t, itsBunch_m->getdT());
                RefPartR_m = referenceToLabCSTrafo_m.transformFrom(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateFrom(RefPartP_m);
            }
        }
    }
}

struct DistributionInfo {
    unsigned int who;
    unsigned int whom;
    unsigned int howMany;
};

void ParallelTTracker::evenlyDistributeParticles() {
    const int numNodes = Ippl::getNodes();
    if (itsBunch_m->hasFieldSolver() || numNodes == 1) return;

    long onAverage = itsBunch_m->getTotalNum() / Ippl::getNodes();
    if (itsBunch_m->getTotalNum() % Ippl::getNodes() > 0.5 * Ippl::getNodes())
        ++ onAverage;

    std::vector<long> localParticles(numNodes, 0);
    localParticles[Ippl::myNode()] = itsBunch_m->getLocalNum() - onAverage;
    allreduce(&(localParticles[0]),
              numNodes,
              std::plus<long>());

    std::vector<DistributionInfo> send;
    std::vector<DistributionInfo> receive;

    for (int i = 0; i < Ippl::getNodes(); ++ i) {
        if (localParticles[i] <= 0) continue;

        for (int j = 0; j < Ippl::getNodes(); ++ j) {
            if (j == i || localParticles[j] >= 0) continue;

            long numParts = std::min(localParticles[i], -localParticles[j]);
            localParticles[i] -= numParts;
            localParticles[j] += numParts;

            if (i == Ippl::myNode() || j == Ippl::myNode()) {
                DistributionInfo request;
                request.who = i;
                request.whom = j;
                request.howMany = numParts;

                if (i == Ippl::myNode()) {
                    send.push_back(request);
                } else {
                    receive.push_back(request);
                }
            }

            if (localParticles[i] == 0) break;
        }
    }

    std::vector<MPI_Request> requests;
    const long sizeSingleParticle = 9 * sizeof(double) + sizeof(short) + sizeof(int) + sizeof(PID_t::Return_t);
    long idx = itsBunch_m->getLocalNum() - 1;
    int tag = Ippl::Comm->next_tag(P_SPATIAL_TRANSFER_TAG, P_LAYOUT_CYCLE);

    std::vector<char> send_msgbuf;

    if (send.size() > 0) {
        const char *buffer;

        unsigned int totalSend = 0, startIndex = 0;
        for (DistributionInfo &request: send) {
            totalSend += request.howMany;
        }
        send_msgbuf.reserve(totalSend * sizeSingleParticle);

        for (DistributionInfo &request: send) {
            size_t sizePrior = send_msgbuf.size();
            for (long i = 0; i < request.howMany; ++ i, -- idx) {
                buffer = reinterpret_cast<const char*>(&(itsBunch_m->R[idx](0)));
                send_msgbuf.insert(send_msgbuf.end(), buffer, buffer + 3 * sizeof(double));
                buffer = reinterpret_cast<const char*>(&(itsBunch_m->P[idx](0)));
                send_msgbuf.insert(send_msgbuf.end(), buffer, buffer + 3 * sizeof(double));
                buffer = reinterpret_cast<const char*>(&(itsBunch_m->Q[idx]));
                send_msgbuf.insert(send_msgbuf.end(), buffer, buffer + sizeof(double));
                buffer = reinterpret_cast<const char*>(&(itsBunch_m->M[idx]));
                send_msgbuf.insert(send_msgbuf.end(), buffer, buffer + sizeof(double));
                buffer = reinterpret_cast<const char*>(&(itsBunch_m->dt[idx]));
                send_msgbuf.insert(send_msgbuf.end(), buffer, buffer + sizeof(double));
                buffer = reinterpret_cast<const char*>(&(itsBunch_m->PType[idx]));
                send_msgbuf.insert(send_msgbuf.end(), buffer, buffer + sizeof(short));
                buffer = reinterpret_cast<const char*>(&(itsBunch_m->TriID[idx]));
                send_msgbuf.insert(send_msgbuf.end(), buffer, buffer + sizeof(int));
                buffer = reinterpret_cast<const char*>(&(itsBunch_m->ID[idx]));
                send_msgbuf.insert(send_msgbuf.end(), buffer, buffer + sizeof(PID_t::Return_t));
            }

            size_t sendSizeThis = send_msgbuf.size() - sizePrior;
            MPI_Request req = Ippl::Comm->raw_isend(&(send_msgbuf[startIndex]),
                                                    sendSizeThis,
                                                    request.whom,
                                                    tag);

            requests.push_back(req);

            startIndex += sendSizeThis;
        }

        itsBunch_m->destroy(totalSend, idx + 1, true);
    }

    for (unsigned int i = 0; i < receive.size(); ++ i) {
        int node = Communicate::COMM_ANY_NODE;
        char *recvbuf;
        const int bufsize = Ippl::Comm->raw_probe_receive(recvbuf, node, tag);

        int j = 0;

        while (j < bufsize) {
            ++ idx;
            itsBunch_m->create(1);
            {
                const double *buffer = reinterpret_cast<const double*>(recvbuf + j);
                itsBunch_m->R[idx] = Vector_t(buffer[0], buffer[1], buffer[2]);
                itsBunch_m->P[idx] = Vector_t(buffer[3], buffer[4], buffer[5]);
                itsBunch_m->Q[idx] = buffer[6];
                itsBunch_m->M[idx] = buffer[7];
                itsBunch_m->dt[idx] = buffer[8];
            }
            j += 9 * sizeof(double);

            {
                const short *buffer = reinterpret_cast<const short*>(recvbuf + j);
                itsBunch_m->PType[idx] = buffer[0];
            }
            j += sizeof(short);

            {
                const int *buffer = reinterpret_cast<const int*>(recvbuf + j);
                itsBunch_m->TriID[idx] = buffer[0];
            }
            j += sizeof(int);

            {
                const PID_t::Return_t *buffer = reinterpret_cast<const PID_t::Return_t*>(recvbuf + j);
                itsBunch_m->ID[idx] = buffer[0];
            }
            j += sizeof(PID_t::Return_t);
        }

        delete[] recvbuf;
    }

    if (requests.size() > 0) {
        MPI_Waitall(requests.size(), &(requests[0]), MPI_STATUSES_IGNORE);
    }
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c++
// c-basic-offset: 4
// indent-tabs-mode:nil
// End:
