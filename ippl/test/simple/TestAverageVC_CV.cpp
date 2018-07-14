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

// TestAverageVC_CV.cpp , Tim Williams 2/6/1997
// This tests the Average functions in [Uniform]Cartesian class, which do 
// weighted averages between Cell and Vertex centered Field's.

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Field/GuardCellSizes.h"
#include "Meshes/UniformCartesian.h"
#include "AppTypes/Vektor.h"

// set dimensionality and problem size
const unsigned Dim3 = 3;
const unsigned nx = 4, ny = 4, nz = 4;


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);
  bool passed = true; // Pass/fail test

  GuardCellSizes<Dim3> gc(1);
  Index I(nx);
  Index J(ny);
  Index K(nz);
  FieldLayout<Dim3> layoutVert(I,J,K);
  Index Ic(nx-1);
  Index Jc(ny-1);
  Index Kc(nz-1);
  FieldLayout<Dim3> layoutCell(Ic,Jc,Kc);

  // Scalar Field, Cell-Centered
  Field<double,Dim3,UniformCartesian<Dim3>,Cell> C(layoutCell,gc);
  C = 1.0;

  // Scalar weight Field, Cell-Centered
  Field<double,Dim3,UniformCartesian<Dim3>,Cell> wC(layoutCell,gc);
  wC = 2.0;

  // Scalar Field, Vert-Centered
  Field<double,Dim3,UniformCartesian<Dim3>,Vert> V(layoutVert,gc);
  V = 1.0;

  // Scalar weight Field, Vert-Centered
  Field<double,Dim3,UniformCartesian<Dim3>,Vert> wV(layoutVert,gc);
  wV = 2.0;

  // Field's to hold weighted averages:
  Field<double,Dim3,UniformCartesian<Dim3>,Cell> avgToC(layoutCell,gc);
  Field<double,Dim3,UniformCartesian<Dim3>,Vert> avgToV(layoutVert,gc);

  // Weighted average from Cell to Vert:
  assign(avgToV, Average(C, wC, avgToV));

  // Weighted average from Vert to Cell:
  assign(avgToC, Average(V, wV, avgToC));

  // Weight from Field<Vektor,Vert> to Field<Vektor,Cell>, using scalar 
  // weight field wV (Field<double,Vert>)
  // Vector Field, Cell-centered
  Field<Vektor<double,Dim3>,Dim3,UniformCartesian<Dim3>,Vert> 
    vV(layoutVert,gc);
  Field<Vektor<double,Dim3>,Dim3,UniformCartesian<Dim3>,Cell> 
    avgToCSW(layoutCell,gc);
  Vektor<double,Dim3> ones;
  for (int d=0; d < Dim3; d++) ones(d) = 1.0;
  vV = ones;
  assign(avgToCSW, Average(vV, wV, avgToCSW));

  // Check results:
  if (sum(avgToV) != nx*ny*nz) {
    passed = false;
    testmsg << "Test of avgToV failed." << endl;
  }
  if (sum(avgToC) != (nx-1)*(ny-1)*(nz-1)) {
    passed = false;
    testmsg << "Test of avgToC failed." << endl;
  }
  // Following triggers two IPPL bugs:
  // 1) unsigned*Vektor<double,Dim> doesn't work, have to cast the
  //    unsigned to a double
  // 2) sum(Field<Vektor,...>,...>) doesn't work.
  //  if (sum(avgToCSW) != (nx-1)*(ny-1)*(nz-1)*ones) passed = false;
  // WORKAROUND:
  Field<Vektor<double,Dim3>,Dim3,UniformCartesian<Dim3>,Cell>::iterator fi;
  Vektor<double,Dim3> thesum = 0.0;
  Vektor<double,Dim3> globalSum = 0.0;
  for (fi = avgToCSW.begin(); fi != avgToCSW.end(); ++fi) thesum += *fi;
  // This gets the global sum to every PE's copy when multiprocessing:
  globalSum = thesum;
  reduce(globalSum, globalSum, OpAddAssign());
  if (globalSum != ((double)(nx-1)*(ny-1)*(nz-1))*ones) {
    passed = false;
    testmsg << "Test of avgToCSW failed." << endl;
  }

  // Test one of the unweighted (2-argument) Average() functions: 
  // Average from Field<Vektor,Vert> to Field<Vektor,Cell>
  // Vector Field, Cell-centered
  Field<Vektor<double,Dim3>,Dim3,UniformCartesian<Dim3>,Vert> 
    vVu(layoutVert,gc);
  Field<Vektor<double,Dim3>,Dim3,UniformCartesian<Dim3>,Cell> 
    avgToCSu(layoutCell,gc);
  vVu = ones;
  assign(avgToCSu, Average(vV, avgToCSu));
  // Check results:
  if (sum(avgToCSu) != (nx-1)*(ny-1)*(nz-1)) {
    passed = false;
    testmsg << "Test of avgToCSu failed." << endl;
  }


  testmsg << ( (passed) ? "PASSED" : "FAILED" ) << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: TestAverageVC_CV.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: TestAverageVC_CV.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/
