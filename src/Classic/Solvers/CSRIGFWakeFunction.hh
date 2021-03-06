#ifndef CSRIGFWAKEFUNCTION_HH
#define CSRIGFWAKEFUNCTION_HH

#include "Solvers/WakeFunction.hh"
#include <vector>
#include <string>

class ElementBase;
class Filter;

class CSRIGFWakeFunction: public WakeFunction {
public:
    CSRIGFWakeFunction(const std::string &name, ElementBase *element, std::vector<Filter *> filters, const unsigned int &N);

    void apply(PartBunchBase<double, 3> *bunch);

    void initialize(const ElementBase * ref);

    virtual const std::string getType() const;

private:
    void calculateLineDensity(PartBunchBase<double, 3> * bunch, std::pair<double, double> &meshInfo);

    void calculateContributionInside(size_t sliceNumber, double angleOfSlice, double meshSpacing);
    void calculateContributionAfter(size_t sliceNumber, double angleOfSlice, double meshSpacing);
    void calculateGreenFunction(PartBunchBase<double, 3> *bunch, double meshSpacing);
    double calcPsi(const double &psiInitial, const double &x, const double &Ds) const;

    std::vector<Filter *> filters_m;
    LineDensity lineDensity_m;
    LineDensity dlineDensitydz_m;

    // Longitudinal CSR field.
    std::vector<double> Ez_m;

    // Chi used to calculate the Green's function
    std::vector<double> Chi_m;

    // Green's function
    std::vector<double> Grn_m;

    // Retarded angle Psi; (for after dipole, Stupakov)
    std::vector<double> Psi_m;

    // Start position of CSR wake.
    double Begin_m;

    // Start position of equivalent hard edge dipole that approximates actual
    // dipole.
    double FieldBegin_m;

    // Effective length of dipole.
    double Length_m;

    // Radius of curvature of effective dipole.
    double bendRadius_m;

    std::string bendName_m;

    double totalBendAngle_m;

};

#endif //CSRIGFWAKEFUNCTION_HH
