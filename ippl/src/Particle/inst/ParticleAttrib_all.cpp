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
#include "Particle/ParticleAttrib.h"
#include "Particle/PAssign.h"
#include "Particle/PAssignDefs.h"
#include "AppTypes/dcomplex.h"

// instantiations
template class ParticleAttrib<unsigned>; // need this for Particle ID
template class ParticleAttrib<int>;
template class ParticleAttrib<IPPL_PRECISION_TYPE>;
template class ParticleAttrib<dcomplex>;


/***************************************************************************
 * $RCSfile: ParticleAttrib_all.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:29 $
 * IPPL_VERSION_ID: $Id: ParticleAttrib_all.cpp,v 1.1.1.1 2003/01/23 07:40:29 adelmann Exp $ 
 ***************************************************************************/
