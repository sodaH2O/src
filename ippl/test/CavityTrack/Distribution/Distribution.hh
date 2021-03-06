#ifndef Distribution_HH
#define Distribution_HH
// ------------------------------------------------------------------------
// $RCSfile: Distribution.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/pcyclint/prog/Distribution/Distribution.hh,v 1.3 2004/10/01 20:33:09 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// $State : $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
//
// ------------------------------------------------------------------------
// Class category: 
// ------------------------------------------------------------------------
//
// $Date: 2004/10/01 20:33:09 $
// $Author: adelmann $
// $Log: Distribution.hh,v $
// Revision 1.3  2004/10/01 20:33:09  adelmann
// add GAUSS and UNIFORM distribution
//
// Revision 1.2  2004/09/24 19:19:13  adelmann
// - Add stuff for Gaussian Distribution
//
// Revision 1.1.1.1  2004/09/22 12:10:44  adelmann
// Imported sources pcyclubt
//
// Revision 1.6  2004/04/06 13:23:19  adelmann
// *** empty log message ***
//
// Revision 1.5  2004/04/05 13:07:52  adelmann
// Many changes to use new IPPL instead of POOMA
//
// Revision 1.4  2004/04/03 12:52:25  adelmann
// New gauss distribution based on cernlib
//
// Revision 1.3  2003/05/02 13:57:22  adelmann
// First electron cloud simulation 3d without SC and probable still with the wrong
// drive beam integration scheme. Otherwise :=) the program runs
// sureprisingly stable and the first inspection of the data looks promissing
//
// Revision 1.2  2003/01/29 13:50:47  adelmann
// Fix READFROMFILE works now
//
// Revision 1.1.1.1  2003/01/23 09:13:58  adelmann
// GenTrackE
//
// ------------------------------------------------------------------------

#include "ChargedParticles.hh"
#include "Configure.h"
#include <string>
#include <cmath>

#define RANLIBX
#include "Distribution/ranlib.h"

#define sqr(x) pow(x,2.)

class Distribution {
    
    typedef ChargedParticles<double,3>::ParticlePos_t ParticlePos_t;
    typedef ChargedParticles<double,3>::Vector_t Vector_t;

public:
  
    // Construction/destruction.
    Distribution(ChargedParticles<double,3> *beam, DistrData *distrData, double len);
    Distribution(ChargedParticles<double,3> *beam, DistrData *distrData, string distrFn);
    
    // recreated distribution without particle generation
    Distribution(ChargedParticles<double,3> *beam);

    ~Distribution();
  
    void createBinom(unsigned long int particles,
        Vector_t emit,
        Vector_t alpha,
        Vector_t beta,
        Vector_t bincoef);

    void regenerateBinom(unsigned long int particles,
        Vector_t emit,
        Vector_t alpha,
        Vector_t beta,
        Vector_t bincoef);
      
    Vector_t generateGaussian(
        double qTot , 
        double mTot, 
        unsigned long particles,
        Vector_t xmean, 
        Vector_t pmean, 
        Vector_t xstddev, 
        Vector_t pstddev, 
        Vector_t angle,
        double ekin, 
        double r, 
        double alpha, 
        double vvertical,
        double PR);
    

    Vector_t generateUniform(
        double qTot, 
        double mTot, 
        unsigned long totalP, 
        Vector_t sigmar, 
        Vector_t sigmav, 
        double ekin, 
        double r, 
        double alpha, 
        double vvertical,
        double PR);
    
    Vector_t readInputDistribution(string fn);
    
    
    void calcError();


    string const convert(double m) {
        if (m>1000.0)
            return "Gauss ";
        else if (m==0.0)
            return "KV    ";
        else
            return "not tabulated distribution";
    }
    
    void print() const
        {
            Inform os("Create distribution  ");
            int pwi = 7;
            if (distrData_m->distrType==BINOMINAL)
                os << "Distribution related informations (binominal distribution   ) " << endl;
            else if (distrData_m->distrType==ELLIPSOIDALUNIFORM)
                os << "Distribution related informations (uniform ellipsoidal distr) " << endl;
            else
                os << "Distribution related informations (distr NOT KNOWN, read in from file) " << endl;
            os << "------------------------------------------------------------------------" << endl;

//    beam_m->calcBeamParameters();

            os << "N=        " << beam_m->getTotalNum() << setw(pwi) << endl;
/*    os << "x (rms) = " << beam_m->getRrms()(0) << setw(pwi)<< endl;
      os << "y (rms) = " << beam_m->getRrms()(1) << setw(pwi) << endl;
      os << "t (rms) = " << beam_m->getRrms()(2)<< setw(pwi) << endl;
      os << "max(r)  = " << beam_m->getRmax() << setw(pwi) << endl;
      os << "min(r)  = " << beam_m->getRmin() << setw(pwi) << endl;
      os << "px (rms)= " << beam_m->getPrms()(0)<< setw(pwi) << endl;
      os << "py (rms)= " << beam_m->getPrms()(1)<< setw(pwi)<< endl;
      os << "pt (rms)= " << beam_m->getPrms()(2)<< setw(pwi) << endl;
      os << "xpx     = " << beam_m->getRPrms()(0)<< setw(pwi)<< endl;
      os << "ypy     = " << beam_m->getRPrms()(1)<< setw(pwi)<< endl;
      os << "tpt     = " << beam_m->getRPrms()(2)<< setw(pwi)<< endl;
      os << "ex (rms)= " << beam_m->getEmitrms()(0)<< setw(pwi)<< endl;
      os << "ey (rms)= " << beam_m->getEmitrms()(1)<< setw(pwi)<< endl;
      os << "et (rms)= " << beam_m->getEmitrms()(2)<< setw(pwi)<< endl;
*/
            if (distrData_m->distrType!=READFROMFILE) {
                os << "mx        " << distrData_m->mx<< endl;
                os << "my        " << distrData_m->my<< endl;
                os << "mt        " << distrData_m->mz << endl;
            }
            os << "----------------------------------------------------------------" << endl;
        }


    double EmittanceUncertainty(double Cxx, double Cxy, double Cyy,double emat[3][3])
        /*-----------------------------------------------------------------------------
         * Calculates the emittance uncertainty from the error matrix
         * ==========================================================
         *
         *---------------------------------------------------------------------------*/
        {
            double tmp = 0.0;
            /*---------------------------------------------------------------------------*/
    
            double b1 = Cxx;
            double b2 = Cyy;
            double b3 = Cxy;
    
            // calculate the derivatives and sum up...
    
            double deb[3] = {b2, -2.0 * b3, b1};
    
            tmp = 0.0;
            for(int k=0; k<3; k++) {
                for(int l=0; l<3; l++) {
                    tmp += emat[k][l] * deb[k] * deb[l];
                } 
            }

            double emittance = sqrt(Cxx*Cyy - sqr(Cxy));
            tmp = tmp / (4.0 * sqr(emittance)); 
    
            return sqrt(tmp);
    
        }
    
    void TwissUncertainty(double Cxx, double Cxy, double Cyy, double emat[3][3],double *dalpha, double *dbeta, double *dgamma)
        {
            double tmp1 = 0.0;
            double tmp2 = 0.0;
            double tmp3 = 0.0;
            /*---------------------------------------------------------------------------*/
    
            double b1 = Cxx;
            double b2 = Cyy;
            double b3 = Cxy;
            double emittance = sqrt(Cxx*Cyy - sqr(Cxy));
    
            // calculate the derivatives and sum up...
    
            double da[3] = {b3*b2, -2.0*(sqr(b3)+sqr(emittance)), b1*b3};
            double db[3] = {(-1.0*b1*b2 + 2.0*sqr(emittance)), 
                            2.0*b1*b3, -1.0*sqr(b1)};
            double dg[3] = {-1.0*sqr(b2), 2.0*b3*b2, 
                            (-1.0*b1*b2 + 2.0*sqr(emittance))};
    
            tmp1 = 0.0;
            tmp2 = 0.0;
            tmp3 = 0.0;
            for(int k=0; k<3; k++) {
                for(int l=0; l<3; l++) {
                    tmp1 += emat[k][l] * da[k] * da[l];
                    tmp2 += emat[k][l] * db[k] * db[l];
                    tmp3 += emat[k][l] * dg[k] * dg[l];
                } 
            }
    
            *dalpha = sqrt(tmp1 / (4.0 * pow(emittance,6.))); 
            *dbeta  = sqrt(tmp2 / (4.0 * pow(emittance,6.))); 
            *dgamma = sqrt(tmp3 / (4.0 * pow(emittance,6.))); 
        }
    
    void RmsUncertainty(double Cxx, double Cxy, double Cyy, double emat[3][3],
        double *dxrms, double *dyrms, double *dxyrms)
        {
 
        }

private:
    // the beam to write
    ChargedParticles<double,3> *beam_m;

    DistrData *distrData_m;

};

// Output operator.
inline ostream &operator<<(ostream &os, const Distribution &data)
{
    //data.print(os);
    return os;
}
#include "Distribution/Distribution.cpp"
#endif // Distribution_HH























