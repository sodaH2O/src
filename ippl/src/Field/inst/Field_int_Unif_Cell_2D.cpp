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
#include "Field/Field.h"
#include "Field/Assign.h"
#include "Field/AssignDefs.h"
#include "Meshes/UniformCartesian.h"
#include "Meshes/Centering.h"
#include "AppTypes/dcomplex.h"

// 2D UniformCartesian Cell instantiations
template class Field<int,2U,UniformCartesian<2U,IPPL_PRECISION_TYPE>,Cell>;


/***************************************************************************
 * $RCSfile: Field_int_Unif_Cell_2D.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:27 $
 * IPPL_VERSION_ID: $Id: Field_int_Unif_Cell_2D.cpp,v 1.1.1.1 2003/01/23 07:40:27 adelmann Exp $ 
 ***************************************************************************/

