// ------------------------------------------------------------------------
// $RCSfile: DefaultVisitor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DefaultVisitor
//   Defines the default interface for a BeamlineVisitor.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/05/03 08:16:04 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"

#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Patch.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"

#include "Algorithms/MapIntegrator.h"
#include "Algorithms/TrackIntegrator.h"

#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"

#include "ComponentWrappers/CorrectorWrapper.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include "ComponentWrappers/RBendWrapper.h"
#include "ComponentWrappers/SBendWrapper.h"

// Class DefaultVisitor
// ------------------------------------------------------------------------

DefaultVisitor::DefaultVisitor(const Beamline &beamline,
			       bool backBeam, bool backTrack):
  itsLine(beamline), back_beam(backBeam), back_track(backTrack)
{
  local_flip = back_path =
    back_beam && ! back_track  ||  back_track && ! back_beam;
  flip_B     = back_beam ? -1.0 : 1.0;
  flip_s     = back_path ? -1.0 : 1.0;
}


DefaultVisitor::~DefaultVisitor()
{}


void DefaultVisitor::execute()
{
  local_flip = back_beam && ! back_track  ||  back_track && ! back_beam;
  itsLine.accept(*this);
}


void DefaultVisitor::visitBeamBeam(const BeamBeam &bb)
{
  applyDefault(bb);
}


void DefaultVisitor::visitCollimator(const Collimator &coll)
{
  applyDefault(coll);
}


void DefaultVisitor::visitComponent(const Component &comp)
{
  applyDefault(comp);
}


void DefaultVisitor::visitCorrector(const Corrector &corr)
{
  applyDefault(corr);
}


void DefaultVisitor::visitDiagnostic(const Diagnostic &diag)
{
  applyDefault(diag);
}


void DefaultVisitor::visitDrift(const Drift &drf)
{
  applyDefault(drf);
}


void DefaultVisitor::visitLambertson(const Lambertson &lamb)
{
  applyDefault(lamb);
}


void DefaultVisitor::visitMarker(const Marker &mark)
{
  applyDefault(mark);
}


void DefaultVisitor::visitMonitor(const Monitor &mon)
{
  applyDefault(mon);
}


void DefaultVisitor::visitMultipole(const Multipole &mult)
{
  applyDefault(mult);
}


void DefaultVisitor::visitPatch(const Patch &patch)
{
  applyDefault(patch);
}


void DefaultVisitor::visitRBend(const RBend &bend)
{
  applyDefault(bend);
}


void DefaultVisitor::visitRFCavity(const RFCavity &cav)
{
  applyDefault(cav);
}


void DefaultVisitor::visitRFQuadrupole(const RFQuadrupole &quad)
{
  applyDefault(quad);
}


void DefaultVisitor::visitSBend(const SBend &bend)
{
  applyDefault(bend);
}


void DefaultVisitor::visitSeparator(const Separator &sep)
{
  applyDefault(sep);
}


void DefaultVisitor::visitSeptum(const Septum &sept)
{
  applyDefault(sept);
}


void DefaultVisitor::visitSolenoid(const Solenoid &sol)
{
  applyDefault(sol);
}


void DefaultVisitor::visitBeamline(const Beamline &bl)
{
  // Default behaviour: Apply algorithm to all beamline members.
  // If flip_local is true, track from right to left.
  bl.iterate(*this, local_flip);
}
  

void DefaultVisitor::visitFlaggedElmPtr(const FlaggedElmPtr &fep)
{
  if (fep.getReflectionFlag()) {
    local_flip = ! local_flip;
    fep.getElement()->accept(*this);
    local_flip = ! local_flip;
  } else {
    fep.getElement()->accept(*this);
  }
}


void DefaultVisitor::visitAlignWrapper(const AlignWrapper &wrap)
{
  // Default behaviour: Apply algorithm to the non-offset element.
  wrap.getElement()->accept(*this);
}


void DefaultVisitor::visitCorrectorWrapper(const CorrectorWrapper &wrap)
{
  visitCorrector(wrap);
}


void DefaultVisitor::visitMultipoleWrapper(const MultipoleWrapper &wrap)
{
  visitMultipole(wrap);
}


void DefaultVisitor::visitRBendWrapper(const RBendWrapper &wrap)
{
  visitRBend(wrap);
}


void DefaultVisitor::visitSBendWrapper(const SBendWrapper &wrap)
{
  visitSBend(wrap);
}


void DefaultVisitor::visitIntegrator(const Integrator &i)
{
  // Default: cannot use integrator.
  i.getElement()->accept(*this);
}


void DefaultVisitor::visitMapIntegrator(const MapIntegrator &i)
{
  // Default: cannot use integrator.
  i.getElement()->accept(*this);
}


void DefaultVisitor::visitTrackIntegrator(const TrackIntegrator &i)
{
  // Default: cannot use integrator.
  i.getElement()->accept(*this);
}


void DefaultVisitor::applyDefault(const ElementBase &)
{}
