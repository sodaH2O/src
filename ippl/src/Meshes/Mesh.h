// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef MESH_H
#define MESH_H

/***********************************************************************
 * 
 * The Mesh base class. Right now, this mainly acts as a standard base
 * class for all meshes so that other objects can register as users of
 * the mesh and can be notified if the mesh changes (e.g., it is rescaled
 * or restructured entirely).
 *
 ***********************************************************************/

// include files
#include "FieldLayout/FieldLayout.h"
#include "FieldLayout/FieldLayoutUser.h"
#include "Utility/UserList.h"

// Enumeration used for specifying mesh boundary conditions. Mesh BC are used
// for things like figuring out how to return the mesh spacing for a cell 
// beyond the edge of the physical mesh, as might arise in stencil operations
// on Field's on the mesh.
enum MeshBC_E { Reflective, Periodic, NoBC }; 

template<unsigned Dim>
class Mesh : private UserList {

public:
  //# public typedefs
  typedef UserList::ID_t             ID_t;
  typedef iterator_user              iterator_if;
  typedef size_type_user             size_type_if;

  //# public enumerations
  enum { Dimension = Dim };

  // static data member
  static std::string MeshBC_E_Names[3];

  // constructor
  Mesh();

  // destructor
  virtual ~Mesh();

  // Mesh geometry queries:
  //-----------------------------------------------------------
  // ....

  // Mesh subsetting functions:
  //----------------------------------------
  // ....
  // Return appropriate SubSetType object for appropriate IndexType argument:
  //  virtual template<class SubSetType, class IndexType> 
  //  SubSetType& getSubSet(IndexType& i);
  // N.B.: maybe above is bad because maybe this:
  // anything except *all* of the mesh can't be generalized to include
  // here; the argument type would have to be appropriate for the Mesh...????

  // Mesh coordinate mapping data and functions:
  //--------------------------------------------
  // Follow Overature's Mapping class design as much as possible.
  // All this is for future implementation, when we need something that is
  // not genuine cartesian geometry.

  // UserList operations
  // -------------------
  // checkin should be called by any objects which want to be informed
  // when the mesh is destroyed or when it is changed (the function
  // 'Repartition' will be called when the mesh changes).  checkout should
  // be called when an object is destroyed or no longer needs to use the
  // Mesh.

  // Return our ID, as generated by UserList.
  ID_t get_Id() const { return getUserListID(); }

  // Tell the Mesh that a FieldLayoutUser has been declared on it.
  // This is just a wrapper around UserList::checkinUser; it only allows
  // FieldLayoutUser's to register.
  void checkin(FieldLayoutUser& f) { checkinUser(f); }

  // Tell the Mesh that a FieldLayoutUser is no longer using it.
  // This is different than the checkoutUser from UserList,
  // for symmetry with checkin and to limit checkout's to FieldLayoutUser's.
  void checkout(FieldLayoutUser& f) { checkoutUser(f); }

  // Accessors for the users accessing this FieldLayout
  size_type_if      size_if() const { return getNumUsers(); }
  iterator_if       begin_if() { return begin_user(); }
  iterator_if       end_if() { return end_user(); }

  // notify all the registered FieldLayoutUser's that this Mesh has
  // changed.  This is done by calling the 'Repartition' virtual function
  // in FieldLayoutUser.
  void notifyOfChange();

};

#include "Meshes/Mesh.hpp"

#endif // MESH_H

/***************************************************************************
 * $RCSfile: Mesh.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:28 $
 * IPPL_VERSION_ID: $Id: Mesh.h,v 1.1.1.1 2003/01/23 07:40:28 adelmann Exp $ 
 ***************************************************************************/
