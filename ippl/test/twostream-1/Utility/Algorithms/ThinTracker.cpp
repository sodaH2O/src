// ------------------------------------------------------------------------
// $RCSfile: ThinTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThinTracker
//   The visitor class for tracking a bunch of particles through a beamline
//   using a thin-lens approximation for all elements.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/ThinTracker.h"
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
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "Algorithms/MapIntegrator.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/PartData.h"
#include "Algorithms/Particle.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/Geometry.h"
#include "Beamlines/Beamline.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Utilities/ComplexErrorFun.h"
#include <cmath>
#include <complex>

using namespace Physics;


// Class ThinTracker
// ------------------------------------------------------------------------

ThinTracker::ThinTracker(const Beamline &beamline, const PartData &reference,
			 bool revBeam, bool revTrack):
  Tracker(beamline, reference, revBeam, revTrack)
{}


ThinTracker::ThinTracker(const Beamline &beamline,
			 const PartBunch &bunch,
			 const PartData &reference,
			 bool revBeam, bool revTrack):
  Tracker(beamline, bunch, reference, revBeam, revTrack)
{}


ThinTracker::~ThinTracker()
{}


void ThinTracker::visitBeamBeam(const BeamBeam &bb)
{
  // If x > explim, exp(-x) is outside machine limits.
  static const double explim = 150.0;
 
  // Parameters of the opposite bunch.
  const double NN = itsReference.getQ()*bb.getBunchCharge();
  const Vector3D &displacement = bb.getBunchDisplacement();
  const Matrix3D &sigma        = bb.getBunchMoment();

  if (NN != 0.0) {
    //
    //               N1 N2     mu_0 c^2 q_e^2     N1 N2
    // fk = 2 r_e * ------- = ---------------- * ------- =
    //               gamma        2 pi m_0        gamma
    //
    //       mu_0 c^2 q_e N1 N2         q_e N1 N2
    //    = -------------------- = -------------------- .
    //            2 pi p_r          2 pi epsilon_0 p_r
    //
    // where m_0 is the particle rest mass in Joule.
    //
    // and p_r = (gamma m_0) / q_e is the momentum in eV.
    //
    const double fk = q_e*NN/
      (epsilon_0*two_pi*itsReference.getP());
    const double dx = displacement(0);
    const double dy = displacement(1);
    const double sx2 = abs(sigma(0,0));
    const double sy2 = abs(sigma(1,1));
    const double sx = sqrt(sx2);
    const double sy = sqrt(sy2);
    
    if (sx2 == sy2) {
      // Limit for sigma(x)^2 = sigma(y)^2.
      const iterator last = itsBunch.end();

      for (iterator part = itsBunch.begin(); part != last; ++part) {
	double xs = part->x() - dx;
	double ys = part->y() - dy;
	double rho2 = xs*xs + ys*ys;
	double tk = rho2/(2.0*sx2);
	double phi = 0.0;

	if (tk > explim) {
	  phi = fk/rho2;
	} else if (tk != 0.0) {
	  phi = fk * (1.0 - exp(- tk))/rho2;
	}

	part->px() += xs*phi;
	part->py() += ys*phi;
      }
    } else {
      // Case sigma(x)^2 != sigma(y)^2.
      const double r = sqrt(2.0 * abs(sx2 - sy2));
      double rk = flip_s * flip_B * fk * sqrt(pi) / r;
      const iterator last = itsBunch.end();

      for (iterator part = itsBunch.begin(); part != last; ++part) {
	double xs = part->x() - dx;
	double ys = part->y() - dy;
	double xr = abs(xs)/r;
	double yr = abs(ys)/r;
	std::complex<double> W = Werrf(std::complex<double>(xr, yr));

	double tk = (xs*xs/sx2 + ys*ys/sy2)/2.0;
	if (tk <= explim) {
	  W -= exp(- tk) * Werrf(std::complex<double>(xr*sy/sx, yr*sx/sy));
	}

	part->px() += rk * ((xs > 0.0) ? std::imag(W) : - std::imag(W));
	part->py() += rk * ((ys > 0.0) ? std::real(W) : - std::real(W));
      }
    }
  }
}


void ThinTracker::visitCollimator(const Collimator &coll)
{
  applyDrift(flip_s * coll.getElementLength());
}


void ThinTracker::visitCorrector(const Corrector &corr)
{
  // Drift through first half of length.
  double length = flip_s * corr.getElementLength();
  if (length) applyDrift(length / 2.0);

  // Apply kick.
  double scale = (flip_s * flip_B * corr.getElementLength() *
		  itsReference.getQ() * c) / itsReference.getP();
  const BDipoleField &field = corr.getField();
  const iterator last = itsBunch.end();

  for (iterator part = itsBunch.begin(); part != last; ++part) {
    part->px() -= field.getBy() * scale;
    part->py() += field.getBx() * scale;
  }

  // Drift through second half of length.
  if (length) applyDrift(length / 2.0);
}


void ThinTracker::visitDiagnostic(const Diagnostic &diag)
{
  // The diagnostic has no effect on particle tracking.
  applyDrift(flip_s * diag.getElementLength());
}


void ThinTracker::visitDrift(const Drift &drift)
{
  applyDrift(flip_s * drift.getElementLength());
}


void ThinTracker::visitLambertson(const Lambertson &lamb)
{
  // Assume the particle go through the magnet's window.
  applyDrift(flip_s * lamb.getElementLength());
}


void ThinTracker::visitMarker(const Marker &)
{
  // Do nothing.
}


void ThinTracker::visitMonitor(const Monitor &corr)
{
  applyDrift(flip_s * corr.getElementLength());
}


void ThinTracker::visitMultipole(const Multipole &mult)
{
  double length = flip_s * mult.getElementLength();
  const BMultipoleField &field = mult.getField();
  double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();

  if (length) {
    // Drift through first half of the length.
    applyDrift(length / 2.0);
    
    // Apply multipole kick, field is per unit length.
    scale *= length;
    applyThinMultipole(field, scale);
    
    // Drift through second half of the length, field is per unit length.
    applyDrift(length / 2.0);
  } else {
    // Thin multipole, field is integral(K*dl).
    scale *= flip_s;
    applyThinMultipole(field, scale);
  }
}


void ThinTracker::visitRBend(const RBend &bend)
{
  // Geometry.
  const RBendGeometry &geometry = bend.getGeometry();
  double length = flip_s * geometry.getElementLength();
  double angle = flip_s * geometry.getBendAngle();

  // Magnetic field.
  const BMultipoleField &field = bend.getField();

  // Drift to mid-plane.
  applyDrift(length / 2.0);

  // Apply multipole kick and linear approximation to geometric bend.
  double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
  if (length) scale *= length;
  iterator last = itsBunch.end();
  for (iterator part = itsBunch.begin(); part != last; ++part) {
    int order = field.order();

    if (order > 0) {
      double x = part->x();
      double y = part->y();
      double kx = + field.normal(order);
      double ky = - field.skew(order);
      
      while (--order > 0) {
	double kxt = x * kx - y * ky;
	double kyt = x * ky + y * kx;
	kx = kxt + field.normal(order);
	ky = kyt - field.skew(order);
      }
      
      part->px() -= kx * scale - angle * (1.0 + part->pt());
      part->py() += ky * scale;
      part->t()  -= angle * x;
    }
  }

  // Drift to out-plane.
  applyDrift(length / 2.0);
}


void ThinTracker::visitRFCavity(const RFCavity &as)
{
  // Drift through half length.
  double length = flip_s * as.getElementLength();
  if (length) applyDrift(length / 2.0);

  // Apply accelerating voltage.
  double freq = as.getFrequency();
  double peak = flip_s * as.getAmplitude() / itsReference.getP();
  double kin = itsReference.getM() / itsReference.getP();
  const iterator last = itsBunch.end();

  for (iterator part = itsBunch.begin(); part != last; ++part) {
    double pt    = (part->pt() + 1.0);
    double speed = (c * pt) / sqrt(pt*pt + kin*kin);
    double phase = as.getPhase() + (freq * part->t()) / speed;
    part->pt() += peak * sin(phase) / pt;
  }

  if (length) applyDrift(length / 2.0);
}


void ThinTracker::visitRFQuadrupole(const RFQuadrupole &rfq)
{
  // *** MISSING *** Tracking for RF Quadrupole.
  applyDrift(flip_s * rfq.getElementLength());
}


void ThinTracker::visitSBend(const SBend &bend)
{
  const PlanarArcGeometry &geometry = bend.getGeometry();
  double length = flip_s * geometry.getElementLength();
  double angle = flip_s * geometry.getBendAngle();

  // Magnetic field.
  const BMultipoleField &field = bend.getField();

  // Drift to mid-plane.
  applyDrift(length / 2.0);

  // Apply multipole kick and linear approximation to geometric bend.
  double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
  if (length) scale *= length;
  iterator last = itsBunch.end();
  for (iterator part = itsBunch.begin(); part != last; ++part) {
    int order = field.order();

    if (order > 0) {
      double x = part->x();
      double y = part->y();
      double kx = + field.normal(order);
      double ky = - field.skew(order);
      
      while (--order > 0) {
	double kxt = x * kx - y * ky;
	double kyt = x * ky + y * kx;
	kx = kxt + field.normal(order);
	ky = kyt - field.skew(order);
      }
      
      part->px() -= kx * scale - angle * (1.0 + part->pt());
      part->py() += ky * scale;
      part->t()  -= angle * x;
    }
  }

  // Drift to out-plane.
  applyDrift(length / 2.0);
}


void ThinTracker::visitSeparator(const Separator &sep)
{
  // Drift through first half of length.
  double length = flip_s * sep.getElementLength();
  if (length) {
    applyDrift(length / 2.0);

    double scale = (length * itsReference.getQ()) / itsReference.getP();
    double Ex = scale * sep.getEx();
    double Ey = scale * sep.getEy();
    const iterator last = itsBunch.end();
    
    for (iterator part = itsBunch.begin(); part != last; ++part) {
      double pt = 1.0 + part->pt();
      part->px() += Ex / pt;
      part->py() += Ey / pt;
    }
    
    applyDrift(length / 2.0);
  }
}


void ThinTracker::visitSeptum(const Septum &sept)
{
  // Assume the particles go through the magnet's window.
  applyDrift(flip_s * sept.getElementLength());
}


void ThinTracker::visitSolenoid(const Solenoid &solenoid)
{
  double length = flip_s * solenoid.getElementLength();

  if (length) {
    double ks = (flip_B*itsReference.getQ()*solenoid.getBz()*c) /
      (2.0*itsReference.getP());

    if (ks) {
      double kin = itsReference.getM() / itsReference.getP();
      double ref = kin * kin;
      const iterator last = itsBunch.end();
      
      for (iterator part = itsBunch.begin(); part != last; ++part) {
	double pt = part->pt() + 1.0;
	double px = part->px() + ks * part->y();
	double py = part->py() - ks * part->x();
	double pz = sqrt(pt*pt - px*px - py*py);
	
	double k = ks / pz;
	double C = cos(k*length);
	double S = sin(k*length);
	
	double xt  = C * part->x()  + S * part->y();
	double yt  = C * part->y()  - S * part->x();
	double pxt = C * part->px() + S * part->py();
	double pyt = C * part->py() - S * part->px();
	
	part->x()  = C * xt  + (S / k) * pxt;
	part->y()  = C * yt  + (S / k) * pyt;
	part->px() = C * pxt - (S * k) * xt;
	part->py() = C * pyt - (S * k) * yt;
	part->t() += length * (pt*ref - (px*px + py*py + 3.0*pt*pt*ref) / 2.0);
      }
    } else {
      applyDrift(length);
    }
  }
}


void ThinTracker::applyDrift(double length)
{
  double   kin  = itsReference.getM() / itsReference.getP();
  double   ref  = kin * kin;
  iterator last = itsBunch.end();

  for (iterator part = itsBunch.begin(); part != last; ++part) {
    double px = part->px();
    double py = part->py();
    double pt = part->pt();
    double lByPz = length / (1.0 + pt);
    part->x() += px * lByPz;
    part->y() += py * lByPz;
    part->t() += length * (pt*ref - (px*px + py*py + 3.0*pt*pt*ref) / 2.0);
  }
}
