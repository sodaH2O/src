#ifndef CLASSIC_FIELDMAP1DDYNAMICFAST_HH
#define CLASSIC_FIELDMAP1DDYNAMICFAST_HH

#include "Fields/Fieldmap.h"

class FM1DDynamic_fast: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E,
                                    Vector_t &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd,
                                    double &rBegin, double &rEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal,
                                    double &yIni, double &yFinal,
                                    double &zIni, double &zFinal) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);
    virtual void getOnaxisEz(std::vector<std::pair<double, double>> &eZ);

    virtual bool isInside(const Vector_t &r) const;
private:
    FM1DDynamic_fast(std::string aFilename);
    ~FM1DDynamic_fast();

    virtual void readMap();
    virtual void freeMap();

    bool checkFileData(std::ifstream &fieldFile, bool parsingPassed);
    void computeFieldDerivatives(std::vector<double> fourierCoefs,
                                 double onAxisFieldP[], double onAxisFieldPP[],
                                 double onAxisFieldPPP[]);
    void computeFieldOffAxis(const Vector_t &R, Vector_t &E, Vector_t &B,
                             std::vector<double> fieldComponents) const;
    void computeFieldOnAxis(double z, std::vector<double> &fieldComponents) const;
    std::vector<double> computeFourierCoefficients(double fieldData[]);
    void computeInterpolationVectors(double onAxisFieldP[],
                                     double onAxisFieldPP[],
                                     double onAxisFieldPPP[]);
    void convertHeaderData();
    void normalizeField(double maxEz, std::vector<double> &fourierCoefs);
    double readFileData(std::ifstream &fieldFile, double fieldData[]);
    double readFileData(std::ifstream &fieldFile,
                        std::vector<std::pair<double, double>> &eZ);
    bool readFileHeader(std::ifstream &fieldFile);
    void scaleField(double maxEz, std::vector<std::pair<double, double>> &eZ);
    void stripFileHeader(std::ifstream &fieldFile);

    void prepareForMapCheck(std::vector<double> &fourierCoefs);

    double frequency_m;                     /// Field angular frequency (Hz).
    double twoPiOverLambdaSq_m;             /// 2 Pi divided by the field RF wavelength squared.

    double rBegin_m;                        /// Minimum radius of field.
    double rEnd_m;                          /// Maximum radius of field.
    double zBegin_m;                        /// Longitudinal start of field.
    double zEnd_m;                          /// Longitudinal end of field.
    double length_m;                        /// Field length.
    unsigned int numberOfGridPoints_m;      /// Number of grid points in field input file.
    double deltaZ_m;                        /// Field grid point spacing.
    unsigned int accuracy_m;

    double* onAxisField_m;                      /// On axis field data.
    gsl_spline *onAxisFieldInterpolants_m;      /// On axis field interpolation structure.
    gsl_spline *onAxisFieldPInterpolants_m;     /// On axis field first derivative interpolation structure.
    gsl_spline *onAxisFieldPPInterpolants_m;    /// On axis field second derivative interpolation structure.
    gsl_spline *onAxisFieldPPPInterpolants_m;   /// On axis field third derivative interpolation structure.

    /// Corresponding interpolation evaluation accelerators.
    gsl_interp_accel *onAxisFieldAccel_m;
    gsl_interp_accel *onAxisFieldPAccel_m;
    gsl_interp_accel *onAxisFieldPPAccel_m;
    gsl_interp_accel *onAxisFieldPPPAccel_m;

    friend class Fieldmap;
};

inline bool FM1DDynamic_fast::isInside(const Vector_t &r) const
{
    return r(2) >= zBegin_m && r(2) < zEnd_m;
}

#endif