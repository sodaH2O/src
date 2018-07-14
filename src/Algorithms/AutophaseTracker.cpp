#include "Algorithms/AutophaseTracker.h"
#include "Utilities/ClassicField.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "Utilities/Timer.h"

#include <iostream>
#include <fstream>

#define DTFRACTION 2

extern Inform *gmsg;

AutophaseTracker::AutophaseTracker(const Beamline &beamline,
                                   const PartData& data,
                                   const double &T0,
                                   double initialR,
                                   double initialP):
    DefaultVisitor(beamline, false, false),
    itsBunch_m(&data),
    itsPusher_m(data),
    particleLayout_m(NULL)
{
    particleLayout_m = new Layout_t();
    timeIntegrationTimer_m = IpplTimings::getTimer("AP: time integration");
    fieldEvaluationTimer_m = IpplTimings::getTimer("AP: field evaluation");

    itsBunch_m.initialize(particleLayout_m);
    itsBunch_m.setT(T0);

    if(Ippl::myNode() == 0) {
        itsBunch_m.create(1);
        itsBunch_m.R[0] = Vector_t(0.0, 0.0, initialR);
        itsBunch_m.P[0] = Vector_t(0.0, 0.0, initialP);
        itsBunch_m.Bin[0] = 0;
        itsBunch_m.Q[0] = itsBunch_m.getChargePerParticle();
        itsBunch_m.PType[0] = ParticleType::REGULAR;
    }

    visitBeamline(beamline);
    itsOpalBeamline_m.prepareSections();
}

AutophaseTracker::~AutophaseTracker()
{ }

void AutophaseTracker::execute(const std::queue<double> &dtAllTracks,
                               const std::queue<double> &maxZ,
                               const std::queue<unsigned long long> &maxTrackSteps) {

    if (Ippl::myNode() != 0) {
        receiveCavityPhases();
        return;
    }
    OPALTimer::Timer clock;
    *gmsg << level1 << "\n\nstart autophasing at " << clock.time() << "\n" << endl;

    Vector_t rmin, rmax;
    Vector_t &refR = itsBunch_m.R[0];
    Vector_t &refP = itsBunch_m.P[0];
    std::queue<double> dtAutoPhasing(dtAllTracks);
    std::queue<double> maxSPosAutoPhasing(maxZ);
    std::queue<unsigned long long> maxStepsAutoPhasing(maxTrackSteps);

    maxSPosAutoPhasing.back() = itsOpalBeamline_m.calcBeamlineLength();

    size_t step = 0;
    const int dtfraction = DTFRACTION;
    double newDt = dtAutoPhasing.front() / dtfraction;
    double scaleFactor = newDt * Physics::c;

    itsOpalBeamline_m.print(*gmsg);

    itsOpalBeamline_m.switchAllElements();

    std::shared_ptr<Component> cavity = NULL;
    while (true) {
        std::shared_ptr<Component> next = getNextCavity(cavity);
        if (next == NULL) break;

        double endNext = getEndCavity(next);
        if (refR(2) < endNext) break;

        cavity = next;
    }

    while (maxStepsAutoPhasing.size() > 0) {
        maxStepsAutoPhasing.front() = maxStepsAutoPhasing.front() * dtfraction + step;
        newDt = dtAutoPhasing.front() / dtfraction;
        scaleFactor = newDt * Physics::c;
        Vector_t vscaleFactor = Vector_t(scaleFactor);
        itsBunch_m.setdT(newDt);
        itsBunch_m.dt = newDt;

        INFOMSG(level2 << "\n"
                << "*************************************************\n"
                << "new set of dt, zstop and max number of time steps\n"
                << "*************************************************\n\n"
                << "start at t = " << itsBunch_m.getT() << " [s], "
                << "zstop at z = " << maxSPosAutoPhasing.front() << " [m],\n"
                << "dt = " << itsBunch_m.dt[0] << " [s], "
                << "step = " << step << ", "
                << "R =  " << refR << " [m]" << endl);

        for(; step < maxStepsAutoPhasing.front(); ++ step) {
            std::shared_ptr<Component> next = getNextCavity(cavity);
            if (next == NULL) break;

            double beginNext, endNext;
            next->getDimensions(beginNext, endNext);
            double sStop = beginNext;
            bool reachesNextCavity = true;
            if (sStop > maxSPosAutoPhasing.front()) {
                sStop = maxSPosAutoPhasing.front();
                reachesNextCavity = false;
            }

            advanceTime(step, maxStepsAutoPhasing.front(), sStop);
            if (step >= maxStepsAutoPhasing.front() || !reachesNextCavity) break;

            INFOMSG(level2 << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
            INFOMSG(level2 << "\nFound " << next->getName()
                    << " at " << sStop << " [m], "
                    << "step  " << step << ", "
                    << "t = " << itsBunch_m.getT() << " [s], "
                    << "E = " << getEnergyMeV(refP) << " [MeV]\n\n"
                    << "start phase scan ... " << endl;);

            double guess = guessCavityPhase(next);
            optimizeCavityPhase(next,
                                guess,
                                step,
                                newDt);

            INFOMSG(level2 << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl);

            track(endNext,
                  step,
                  dtAutoPhasing,
                  maxSPosAutoPhasing,
                  maxStepsAutoPhasing);

            cavity = next;
        }

        if (dtAutoPhasing.size() == 0) break;

        dtAutoPhasing.pop();
        maxSPosAutoPhasing.pop();
        maxStepsAutoPhasing.pop();
    }

    INFOMSG(level2 << "step = " << step << ", spos = " << refR(2) << " [m], "
            << "t= " << itsBunch_m.getT() << " [s], " << "E = " << getEnergyMeV(refP) << " [MeV] "
            << endl);

    printCavityPhases();
    sendCavityPhases();

    *gmsg << level1 << "\n\nfinished autophasing " << clock.time() << "\n\n"  << endl;
}

void AutophaseTracker::track(double uptoSPos,
                             size_t &step,
                             std::queue<double> &dtAutoPhasing,
                             std::queue<double> &maxSPosAutoPhasing,
                             std::queue<unsigned long long> &maxStepsAutoPhasing) {
    Vector_t &refR = itsBunch_m.R[0];
    Vector_t &refP = itsBunch_m.P[0];

    double t = itsBunch_m.getT();
    const int dtfraction = DTFRACTION;

    while (maxStepsAutoPhasing.size() > 0) {
        maxStepsAutoPhasing.front() = maxStepsAutoPhasing.front() * dtfraction + step;
        double newDt = dtAutoPhasing.front() / dtfraction;
        double scaleFactor = newDt * Physics::c;
        Vector_t vscaleFactor = Vector_t(scaleFactor);
        itsBunch_m.setdT(newDt);
        itsBunch_m.dt[0] = newDt;

        for(; step < maxStepsAutoPhasing.front(); ++ step) {
            IpplTimings::startTimer(timeIntegrationTimer_m);

            refR /= vscaleFactor;
            itsPusher_m.push(refR, refP, itsBunch_m.dt[0]);
            refR *= vscaleFactor;

            IpplTimings::stopTimer(timeIntegrationTimer_m);

            evaluateField();

            IpplTimings::startTimer(timeIntegrationTimer_m);

            itsPusher_m.kick(refR, refP, itsBunch_m.Ef[0], itsBunch_m.Bf[0], itsBunch_m.dt[0]);

            refR /= vscaleFactor;
            itsPusher_m.push(refR, refP, itsBunch_m.dt[0]);
            refR *= vscaleFactor;

            t += itsBunch_m.getdT();
            itsBunch_m.setT(t);
            IpplTimings::stopTimer(timeIntegrationTimer_m);

            if (refR(2) > uptoSPos) return;

            if(refR(2) > maxSPosAutoPhasing.front())
                maxStepsAutoPhasing.front() = step;

            if(step % 5000 == 0) {
                INFOMSG(level2 << "step = " << step << ", spos = " << refR(2) << " [m], "
                        << "t= " << itsBunch_m.getT() << " [s], " << "E = " << getEnergyMeV(refP) << " [MeV] "
                        << endl);
            }
        }

        dtAutoPhasing.pop();
        maxSPosAutoPhasing.pop();
        maxStepsAutoPhasing.pop();
    }
}

std::shared_ptr<Component> AutophaseTracker::getNextCavity(const std::shared_ptr<Component> &current) {
    auto cavities = itsOpalBeamline_m.getSuccessors(current);

    if (cavities.size() > 0) {
        return cavities.front();
    }

    std::shared_ptr<Component> noCavity;
    return noCavity;
}

void AutophaseTracker::advanceTime(size_t &step, const size_t maxSteps, const double beginNextCavity) {
    Vector_t &refR = itsBunch_m.R[0];
    Vector_t &refP = itsBunch_m.P[0];

    const double distance = beginNextCavity - refR(2);
    if (distance < 0.0) return;

    const double beta = refP(2) / sqrt(dot(refP, refP) + 1.0);
    size_t  numSteps = std::floor(distance / (beta * Physics::c * itsBunch_m.dt[0]));
    numSteps = std::min(numSteps, maxSteps - step);
    double timeDifference = numSteps * itsBunch_m.dt[0];

    itsBunch_m.setT(timeDifference + itsBunch_m.getT());
    step += numSteps;
    refR = Vector_t(0.0, 0.0, refR(2) + beta * Physics::c * timeDifference);
    refP = Vector_t(0.0, 0.0, sqrt(dot(refP, refP)));
}

double AutophaseTracker::guessCavityPhase(const std::shared_ptr<Component> &cavity) {
    Vector_t &refR = itsBunch_m.R[0];
    Vector_t &refP = itsBunch_m.P[0];

    double cavityStart = getBeginCavity(cavity);
    double orig_phi = 0.0, Phimax = 0.0;

    const double beta = std::max(refP(2) / sqrt(dot(refP, refP) + 1.0), 1e-6);
    const double tErr  = (cavityStart - refR(2)) / (Physics::c * beta);

    bool apVeto;
    if(cavity->getType() == ElementBase::TRAVELINGWAVE) {
        TravelingWave *element = static_cast<TravelingWave *>(cavity.get());
        orig_phi = element->getPhasem();
        apVeto = element->getAutophaseVeto();
        if(apVeto) {
            INFOMSG(level2 << " ----> APVETO -----> "
                    << element->getName() <<  endl);
            return orig_phi;
        }
        INFOMSG(level2 << cavity->getName() << ", "
                << "phase = " << orig_phi << " rad" << endl);

        Phimax = element->getAutoPhaseEstimate(getEnergyMeV(refP),
                                               itsBunch_m.getT() + tErr,
                                               itsBunch_m.getQ(),
                                               itsBunch_m.getM() * 1e-6);
    } else {
        RFCavity *element = static_cast<RFCavity *>(cavity.get());
        orig_phi = element->getPhasem();
        apVeto = element->getAutophaseVeto();
        if(apVeto) {
            INFOMSG(level2 << " ----> APVETO -----> "
                    << element->getName() << endl);
            return orig_phi;
        }
        INFOMSG(level2 << cavity->getName() << ", "
                << "phase = " << orig_phi << " rad" << endl);

        Phimax = element->getAutoPhaseEstimate(getEnergyMeV(refP),
                                               itsBunch_m.getT() + tErr,
                                               itsBunch_m.getQ(),
                                               itsBunch_m.getM() * 1e-6);
    }

    return Phimax;
}

double AutophaseTracker::optimizeCavityPhase(const std::shared_ptr<Component> &cavity,
                                             double initialPhase,
                                             size_t currentStep,
                                             double dt) {

    if(cavity->getType() == ElementBase::TRAVELINGWAVE) {
        if (static_cast<TravelingWave *>(cavity.get())->getAutophaseVeto()) return initialPhase;
    } else {
        if (static_cast<RFCavity *>(cavity.get())->getAutophaseVeto()) return initialPhase;
    }

    itsBunch_m.setdT(dt);
    itsBunch_m.dt[0] = dt;
    double Phimax = initialPhase;
    double phi = initialPhase;
    double originalPhase = 0.0, PhiAstra = 0.0;
    double dphi = Physics::pi / 360.0;
    double cavityStart = getBeginCavity(cavity);
    const int numRefinements = Options::autoPhase;

    int j = -1;

    double E = APtrack(cavity, cavityStart, phi);
    double Emax = E;
    do {
        j ++;
        Emax = E;
        initialPhase = phi;
        phi -= dphi;
        E = APtrack(cavity, cavityStart, phi);
    } while(E > Emax);

    if(j == 0) {
        phi = initialPhase;
        E = Emax;
        j = -1;
        do {
            j ++;
            Emax = E;
            initialPhase = phi;
            phi += dphi;
            E = APtrack(cavity, cavityStart, phi);
        } while(E > Emax);
    }

    for(int rl = 0; rl < numRefinements; ++ rl) {
        dphi /= 2.;
        phi = initialPhase - dphi;
        E = APtrack(cavity, cavityStart, phi);
        if(E > Emax) {
            initialPhase = phi;
            Emax = E;
        } else {
            phi = initialPhase + dphi;
            E = APtrack(cavity, cavityStart, phi);
            if(E > Emax) {
                initialPhase = phi;
                Emax = E;
            }
        }
    }
    Phimax = std::fmod(initialPhase, 2 * Physics::pi);

    if(cavity->getType() == ElementBase::TRAVELINGWAVE) {
        TravelingWave *element = static_cast<TravelingWave *>(cavity.get());
        originalPhase = element->getPhasem();
        double newPhase = std::fmod(originalPhase + Phimax, 2 * Physics::pi);
        element->updatePhasem(newPhase);
    } else {
        RFCavity *element = static_cast<RFCavity *>(cavity.get());
        originalPhase = element->getPhasem();
        double newPhase = std::fmod(originalPhase + Phimax, 2 * Physics::pi);
        element->updatePhasem(newPhase);
    }

    PhiAstra = std::fmod((Phimax * Physics::rad2deg) + 90.0, 360.0);

    INFOMSG(cavity->getName() << "_phi = "  << Phimax << " rad / "
            << Phimax *Physics::rad2deg <<  " deg, AstraPhi = " << PhiAstra << " deg,\n"
            << "E = " << Emax << " (MeV), " << "phi_nom = " << originalPhase *Physics::rad2deg << "\n");

    OpalData::getInstance()->setMaxPhase(cavity->getName(), Phimax);

    double cavBegin, cavEnd;
    cavity->getDimensions(cavBegin, cavEnd);
    OpalData::getInstance()->addEnergyData(cavEnd, Emax * 1e6);

    return Phimax + originalPhase;
}

void AutophaseTracker::evaluateField() {
    IpplTimings::startTimer(fieldEvaluationTimer_m);

    Vector_t &refR = itsBunch_m.R[0];

    itsBunch_m.Ef = Vector_t(0.0);
    itsBunch_m.Bf = Vector_t(0.0);

    Vector_t externalE, externalB;
    itsOpalBeamline_m.getFieldAt(0, refR, ls, itsBunch_m.getT() + itsBunch_m.dt[0] / 2., externalE, externalB);

    itsBunch_m.Ef[0] += externalE;
    itsBunch_m.Bf[0] += externalB;

    IpplTimings::stopTimer(fieldEvaluationTimer_m);
}

double AutophaseTracker::APtrack(const std::shared_ptr<Component> &cavity, double cavityStart, const double &phi) const {
    const Vector_t &refR = itsBunch_m.R[0];
    const Vector_t &refP = itsBunch_m.P[0];

    double beta = refP(2) / std::sqrt(dot(refP, refP) + 1.0);
    beta = std::max(beta, 0.0001);
    double tErr  = (cavityStart - refR(2)) / (Physics::c * beta);

    double initialPhase = 0.0;
    double finalMomentum = 0.0;
    if(cavity->getType() == ElementBase::TRAVELINGWAVE) {
        TravelingWave *tws = static_cast<TravelingWave *>(cavity.get());
        initialPhase = tws->getPhasem();
        tws->updatePhasem(phi);
        std::pair<double, double> pe = tws->trackOnAxisParticle(refP(2),
                                                                itsBunch_m.getT() + tErr,
                                                                itsBunch_m.dt[0],
                                                                itsBunch_m.getQ(),
                                                                itsBunch_m.getM() * 1e-6);
        finalMomentum = pe.first;
        tws->updatePhasem(initialPhase);
    } else {
        RFCavity *rfc = static_cast<RFCavity *>(cavity.get());
        initialPhase = rfc->getPhasem();
        rfc->updatePhasem(phi);

        std::pair<double, double> pe = rfc->trackOnAxisParticle(refP(2),
                                                                itsBunch_m.getT() + tErr,
                                                                itsBunch_m.dt[0],
                                                                itsBunch_m.getQ(),
                                                                itsBunch_m.getM() * 1e-6);
        finalMomentum = pe.first;
        rfc->updatePhasem(initialPhase);
    }
    double finalGamma = sqrt(1.0 + finalMomentum * finalMomentum);
    double finalKineticEnergy = (finalGamma - 1.0) * itsBunch_m.getM() * 1e-6;
    return finalKineticEnergy;
}


void AutophaseTracker::sendCavityPhases() {

    int tag = 101;
    Message *mess = new Message();
    putMessage(*mess, OpalData::getInstance()->getNumberOfMaxPhases());

    {
        typedef std::vector<MaxPhasesT>::iterator iterator_t;
        iterator_t it = OpalData::getInstance()->getFirstMaxPhases();
        iterator_t end = OpalData::getInstance()->getLastMaxPhases();
        for(; it != end; ++ it) {
            putMessage(*mess, (*it).first);
            putMessage(*mess, (*it).second);
        }
    }
    {
        typedef energyEvolution_t::iterator iterator_t;
        iterator_t it = OpalData::getInstance()->getFirstEnergyData();
        iterator_t end = OpalData::getInstance()->getLastEnergyData();
        for(; it != end; ++ it) {
            putMessage(*mess, (*it).first);
            putMessage(*mess, (*it).second);
        }
    }
    Ippl::Comm->broadcast_all(mess, tag);

    delete mess;
}

void AutophaseTracker::receiveCavityPhases() {

    int parent = 0;
    int tag = 101;
    int nData = 0;
    Message *mess = Ippl::Comm->receive_block(parent, tag);
    getMessage(*mess, nData);

    {
        typedef std::vector<MaxPhasesT>::iterator iterator_t;
        for(int i = 0; i < nData; i++) {
            std::string elName;
            double maxPhi;
            getMessage(*mess, elName);
            getMessage(*mess, maxPhi);
            OpalData::getInstance()->setMaxPhase(elName, maxPhi);
        }
    }
    {
        typedef energyEvolution_t::iterator iterator_t;
        for(int i = 0; i < nData; i++) {
            double spos;
            double maxEnergy;
            getMessage(*mess, spos);
            getMessage(*mess, maxEnergy);
            OpalData::getInstance()->addEnergyData(spos, maxEnergy);
        }
    }
    delete mess;
}

void AutophaseTracker::printCavityPhases() {
    std::shared_ptr<Component> cavity = NULL;
    const double RADDEG = 180.0 / Physics::pi;
    const double globalTimeShift = OpalData::getInstance()->getGlobalPhaseShift();

    *gmsg << level1 << "\n-------------------------------------------------------------------------------------\n";
    while (true) {
        std::shared_ptr<Component> next = getNextCavity(cavity);
        if (next == NULL) break;

        std::string name = next->getName();
        double frequency;
        double phase;

        if (next->getType() == ElementBase::TRAVELINGWAVE) {
            phase = static_cast<TravelingWave *>(next.get())->getPhasem();
	    frequency = static_cast<TravelingWave *>(next.get())->getFrequencym();
        } else {
            phase = static_cast<RFCavity *>(next.get())->getPhasem();
	    frequency = static_cast<RFCavity *>(next.get())->getFrequencym();
        }

        phase -= globalTimeShift * frequency;

        *gmsg << (cavity == NULL? "": "\n")
            << name
            << ": phi = phi_nom + phi_maxE + global phase shift = " << phase * RADDEG << " degree, "
            << "(global phase shift = " << -globalTimeShift *frequency *RADDEG << " degree) \n";

        cavity = next;
    }
    *gmsg << "-------------------------------------------------------------------------------------\n"
          << endl;
}

void AutophaseTracker::save(const std::string &fname) {
    std::ofstream out(fname);

    std::shared_ptr<Component> cavity = NULL;
    const double globalTimeShift = OpalData::getInstance()->getGlobalPhaseShift();

    while (true) {
        std::shared_ptr<Component> next = getNextCavity(cavity);
        if (next == NULL) break;

        std::string name = next->getName();
        double frequency;
        double phase;

        if (next->getType() == ElementBase::TRAVELINGWAVE) {
            phase = static_cast<TravelingWave *>(next.get())->getPhasem();
	    frequency = static_cast<TravelingWave *>(next.get())->getFrequencym();
        } else {
            phase = static_cast<RFCavity *>(next.get())->getPhasem();
	    frequency = static_cast<RFCavity *>(next.get())->getFrequencym();
        }

        phase -= globalTimeShift * frequency;

        out << name << "_lag = " << phase << ";\n";

        cavity = next;
    }

    out.close();
}
