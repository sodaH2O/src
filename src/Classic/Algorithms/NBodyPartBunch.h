#ifndef NBODY_PART_BUNCH_H
#define NBODY_PART_BUNCH_H

#include "Algorithms/PartBunch.h"

class NBodyPartBunch: public PartBunch {

public:

    NBodyPartBunch(const PartData *ref);
    NBodyPartBunch(const std::vector<OpalParticle> &,
                   const PartData *ref);
    NBodyPartBunch(const PartBunch &);
    ~NBodyPartBunch();

    void computeSelfFields();
    void computeSelfFields(int b);

    void computeSelfFields_cycl(double gamma);
    void computeSelfFields_cycl(int b);

private:
    bool selfFieldIsOn;
};

#endif
