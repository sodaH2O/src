#include "AmrObject.h"

AmrObject::AmrObject()
    : tagging_m(CHARGE_DENSITY),
      scaling_m(0.75),
      chargedensity_m(1.0e-15),
      maxNumPart_m(1),
      minNumPart_m(1),
      refined_m(false),
      amrSolveTimer_m(IpplTimings::getTimer("AMR solve"))
{ }


AmrObject::AmrObject(TaggingCriteria tagging,
                     double scaling,
                     double chargedensity)
    : tagging_m(tagging),
      scaling_m(scaling),
      chargedensity_m(chargedensity),
      maxNumPart_m(1),
      minNumPart_m(1),
      refined_m(false),
      amrSolveTimer_m(IpplTimings::getTimer("AMR solve"))
{ }


AmrObject::~AmrObject()
{ }


void AmrObject::setTagging(TaggingCriteria tagging) {
    tagging_m = tagging;
}


void AmrObject::setTagging(std::string tagging) {
    tagging = Util::toUpper(tagging);
    
    if ( tagging == "POTENTIAL" )
        tagging_m = TaggingCriteria::POTENTIAL;
    else if (tagging == "EFIELD" )
        tagging_m = TaggingCriteria::EFIELD;
    else if ( tagging == "MOMENTA" )
        tagging_m = TaggingCriteria::MOMENTA;
    else if ( tagging == "MAX_NUM_PARTICLES" )
        tagging_m = TaggingCriteria::MAX_NUM_PARTICLES;
    else if ( tagging == "MIN_NUM_PARTICLES" )
        tagging_m = TaggingCriteria::MIN_NUM_PARTICLES;
    else if ( tagging == "CHARGE_DENSITY" )
        tagging_m = TaggingCriteria::CHARGE_DENSITY;
    else
        throw OpalException("AmrObject::setTagging(std::string)",
                            "Not supported refinement criteria "
                            "[CHARGE_DENSITY | POTENTIAL | EFIELD | "
                            "MOMENTA | MAX_NUM_PARTICLES | MIN_NUM_PARTICLES].");
}


void AmrObject::setScalingFactor(double scaling) {
    scaling_m = scaling;
}


void AmrObject::setChargeDensity(double chargedensity) {
    chargedensity_m = chargedensity;
}


void AmrObject::setMaxNumParticles(size_t maxNumPart) {
    maxNumPart_m = maxNumPart;
}


void AmrObject::setMinNumParticles(size_t minNumPart) {
    minNumPart_m = minNumPart;
}


const bool& AmrObject::isRefined() const {
    return refined_m;
}
