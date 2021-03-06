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
#include "Field/BCond.h"
#include "Meshes/UniformCartesian.h"
#include "Meshes/Centering.h"
#include "AppTypes/dcomplex.h"

// 2D UniformCartesian Vert instantiations
template class BConds<int,2U,UniformCartesian<2U,IPPL_PRECISION_TYPE>,Vert>;
template class BConds<dcomplex,2U,UniformCartesian<2U,IPPL_PRECISION_TYPE>,Vert>;
template class BConds<IPPL_PRECISION_TYPE,2U,
                      UniformCartesian<2U,IPPL_PRECISION_TYPE>,Vert>;


/***************************************************************************
 * $RCSfile: BCond_Unif_Vert_2D.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:26 $
 * IPPL_VERSION_ID: $Id: BCond_Unif_Vert_2D.cpp,v 1.1.1.1 2003/01/23 07:40:26 adelmann Exp $ 
 ***************************************************************************/
