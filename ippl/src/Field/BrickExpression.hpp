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
 * Visit www.amas.web.psi for more details
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



//////////////////////////////////////////////////////////////////////
// BrickExpLoop::apply member functions

template<class LHS, class RHS, class OP, unsigned Dim>
class BrickExpLoop
{
public: 
  static inline void
  apply(LHS& restrict Lhs, RHS& restrict Rhs, OP Op)
    {
      int n0 = Lhs.size(0);
      int n1 = Lhs.size(1);
      int n2 = Lhs.size(2);
      if ( (n0>0)&&(n1>0)&&(n2>0) )
	{
	  unsigned d;
	  do
	    {
	      for (int i2=0; i2<n2; ++i2)
		for (int i1=0; i1<n1; ++i1)
		  for (int i0=0; i0<n0; ++i0)
		    PETE_apply(Op,Lhs.offset(i0,i1,i2),
			       for_each(Rhs,EvalFunctor_3(i0,i1,i2)));

	      for (d=3; d<Dim; ++d)
		{
		  Lhs.step(d);
		  // #pragma inline here
		  for_each(Rhs,StepFunctor(d),PETE_NullCombiner());
		  if ( ! Lhs.done(d) ) 
		    break;
		  Lhs.rewind(d);
		  // #pragma inline here
		  for_each(Rhs,RewindFunctor(d),PETE_NullCombiner());
		} 
	    } while (d<Dim);
	}
    }
};


//////////////////////////////////////////////////////////////////////
//a specialization of BrickExpLoop::apply for a 1D loop evaluation
template<class LHS, class RHS, class OP>
class BrickExpLoop<LHS,RHS,OP,1U>
{
public:
  static inline void apply(LHS& restrict Lhs, RHS& restrict Rhs, OP Op)
    {
      int n0 = Lhs.size(0);
      for (int i0=0; i0<n0; ++i0) {
	PETE_apply(Op,Lhs.offset(i0),for_each(Rhs,EvalFunctor_1(i0)));
      }
    }
};

//////////////////////////////////////////////////////////////////////
//a specialization of BrickExpLoops::apply for a 2D loop evaluation
template<class LHS, class RHS, class OP>
class BrickExpLoop<LHS,RHS,OP,2U>
{
public:
  static inline void
  apply(LHS& restrict Lhs, RHS& restrict Rhs, OP Op)
    {
      int n0 = Lhs.size(0);
      int n1 = Lhs.size(1);
      for (int i1=0; i1<n1; ++i1)
	for (int i0=0; i0<n0; ++i0)
	  PETE_apply(Op,Lhs.offset(i0,i1),
		     for_each(Rhs,EvalFunctor_2(i0,i1)));
    }
};


//////////////////////////////////////////////////////////////////////
//a specialization of BrickExpLoops::apply for a 3D loop evaluation
template<class LHS, class RHS, class OP>
class BrickExpLoop<LHS,RHS,OP,3U>
{
public:
  static inline void
  apply(LHS& restrict Lhs, RHS& restrict Rhs, OP Op)
    {
      int n0 = Lhs.size(0);
      int n1 = Lhs.size(1);
      int n2 = Lhs.size(2);
      for (int i2=0; i2<n2; ++i2) 
	for (int i1=0; i1<n1; ++i1) 
	  for (int i0=0; i0<n0; ++i0) 
	    PETE_apply(Op,Lhs.offset(i0,i1,i2),
		       for_each(Rhs,EvalFunctor_3(i0,i1,i2)));
    }
};

//////////////////////////////////////////////////////////////////////
// BrickExpression::apply - just use BrickExpLoop
// ada: remove restrict from apply to make  g++ 2.95.3 happy
template<unsigned Dim, class LHS, class RHS, class OP>
#ifdef IPPL_RESTRICT_BUG
 void BrickExpression<Dim,LHS,RHS,OP>::apply()
#else
 void BrickExpression<Dim,LHS,RHS,OP>::apply() restrict
#endif
{
  
  
  
  // #pragma inline here
  BrickExpLoop<LHS,RHS,OP,Dim>::apply(Lhs,Rhs,Op);
}

/***************************************************************************
 * $RCSfile: BrickExpression.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:26 $
 * IPPL_VERSION_ID: $Id: BrickExpression.cpp,v 1.1.1.1 2003/01/23 07:40:26 adelmann Exp $ 
 ***************************************************************************/
