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
#include "DataSource/FileFieldDataSource.h"
#include "Meshes/Cartesian.h"
#include "Meshes/Centering.h"
#include "AppTypes/Vektor.h"
#include "AppTypes/Tenzor.h"
#include "AppTypes/SymTenzor.h"
#include "AppTypes/AntiSymTenzor.h"

// 3D Cartesian Vert instantiations
template class FileFieldDataSource<Vektor<IPPL_PRECISION_TYPE,3U>, 3U,
                                   Cartesian<3U,IPPL_PRECISION_TYPE>, Vert>;
template class FileFieldDataSource<Tenzor<IPPL_PRECISION_TYPE,3U>, 3U,
                                   Cartesian<3U,IPPL_PRECISION_TYPE>, Vert>;
template class FileFieldDataSource<SymTenzor<IPPL_PRECISION_TYPE,3U>, 3U,
                                   Cartesian<3U,IPPL_PRECISION_TYPE>, Vert>;
template class FileFieldDataSource<AntiSymTenzor<IPPL_PRECISION_TYPE,3U>, 3U,
                                   Cartesian<3U,IPPL_PRECISION_TYPE>, Vert>;


/***************************************************************************
 * $RCSfile: FileFieldAppTypesDS_Cart_Vert_3D.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:25 $
 * IPPL_VERSION_ID: $Id: FileFieldAppTypesDS_Cart_Vert_3D.cpp,v 1.1.1.1 2003/01/23 07:40:25 adelmann Exp $ 
 ***************************************************************************/

