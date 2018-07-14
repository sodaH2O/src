// ------------------------------------------------------------------------
// $RCSfile: Gauss.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Gaussian density function.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/Gauss.h"
#include "Physics/Physics.h"
#include <cmath>


using Physics::two_pi;

double Gauss(double sigx, double sigy, double dx, double dy) {
    double xx = dx / sigx;
    double yy = dy / sigy;
    double arg = (xx * xx + yy * yy) / 2.0;

    if(arg < 72.0) {
        return exp(- arg) / (two_pi * sigx * sigy);
    } else {
        return 0.0;
    }
}

