// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "Particle/IpplParticleBase.h"
#include "Particle/ParticleSpatialLayout.h"
#include "Meshes/UniformCartesian.h"

// 2D UniformCartesian instantiations
template class IpplParticleBase<ParticleSpatialLayout<
   IPPL_PRECISION_TYPE, 2U, UniformCartesian<2U,IPPL_PRECISION_TYPE> > >;


/***************************************************************************
 * $RCSfile: ParticleBaseSpatial_Unif_2D.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:29 $
 * IPPL_VERSION_ID: $Id: ParticleBaseSpatial_Unif_2D.cpp,v 1.1.1.1 2003/01/23 07:40:29 adelmann Exp $
 ***************************************************************************/