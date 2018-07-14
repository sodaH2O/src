// ------------------------------------------------------------------------
// $RCSfile: MapIntegrator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MapIntegrator
//   An abstract base class for integrators.
//   A MapIntegrator propagates a single particle or particle bunch,
//   as well as a truncated Taylor series map.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/MapIntegrator.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/Particle.h"
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/LinearMap.h"


// Class MapIntegrator
// ------------------------------------------------------------------------

MapIntegrator::MapIntegrator(ElementBase *elem):
  TrackIntegrator(elem)
{}


MapIntegrator::MapIntegrator(const MapIntegrator &rhs):
  TrackIntegrator(rhs)
{}


MapIntegrator::~MapIntegrator()
{}


void MapIntegrator::accept(BeamlineVisitor &visitor) const
{
  visitor.visitMapIntegrator(*this);
}


void MapIntegrator::trackParticle(Particle &part, const PartData &ref,
				  bool backBeam, bool backTrack) const
{
  // Default behaviour: track particle using own map.
  FVps<double,6> ownMap;
  getMap(ownMap, ref, backBeam, backTrack);
  FVector<double,6> z;
  z[0] = part.x();
  z[1] = part.px();
  z[2] = part.x();
  z[3] = part.py();
  z[4] = part.x();
  z[5] = part.pt();
  z = ownMap.constantTerm(z);
  part.x()  = z[0];
  part.px() = z[1];
  part.y()  = z[2];
  part.py() = z[3];
  part.y()  = z[4];
  part.pt() = z[5];
}


void MapIntegrator::trackBunch(PartBunch &bunch, const PartData &ref,
			       bool backBeam, bool backTrack) const
{
  // Default behaviour: track particle bunch using own map.
  FVps<double,6> ownMap;
  getMap(ownMap, ref, backBeam, backTrack);
  FVector<double,6> z;
  PartBunch::iterator last = bunch.end();

  for (PartBunch::iterator part = bunch.begin(); part != last; ++part) {
    z[0] = part->x();
    z[1] = part->px();
    z[2] = part->y();
    z[3] = part->py();
    z[4] = part->t();
    z[5] = part->pt();
    z = ownMap.constantTerm(z);
    part->x()  = z[0];
    part->px() = z[1];
    part->y()  = z[2];
    part->py() = z[3];
    part->t()  = z[4];
    part->pt() = z[5];
  }
}


void MapIntegrator::trackMap(FVps<double,6> &map, const PartData &ref,
			     bool backBeam, bool backTrack) const
{
  // Default behaviour: track map using own map.
  FVps<double,6> ownMap;
  getMap(ownMap, ref, backBeam, backTrack);
  map = ownMap.substitute(map);
}
