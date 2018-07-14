#ifndef POISSONPROBLEMS_H
#define POISSONPROBLEMS_H

#include <BoxLib.H>
#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>

#include "Solver.h"
#include "AmrPartBunch.h"
#include "Distribution.h"

#include "writePlotFile.H"

/*!
 * @file PoissonProblems.h
 * @author Matthias Frey
 * @date 25. - 26. October 2016
 * @details Solve \f$ \Delta\phi = \rho \f$ on \f$ [0, 1]\times [0, 1]\times [0, 1]\f$ for different \f$\rho\f$.\n
 * Every function returns the \f$L_{2}\f$-error compared to the solution of a single-level problem.
 * 
 * - doSolveNoParticles:          \f$\rho = -1\f$ everywhere on the domain (no particles)\n
 * - doSolveParticlesUniform:     \f$\rho = -1\f$ everywhere (initialized by particles, idea: we generate particles on the
 *                                cell center and then do an AssignDensity-function call. The charge is scaled such that
 *                                \f$\rho = -1\f$ everywhere on the domain\n
 * - doSolveParticlesGaussian:    \f$\rho\f$ is a Gaussian distribution initialized by particles\n
 * - doSolveParticlesReal:        Read in a H5 cyclotron file and use its particle distribution for \f$\rho\f$\n
 * @brief Several ways of Poisson problems, i.e. different \f$\rho\f$
 */

/// Defines several Poisson problems and solves them on \f$ [0, 1]\times [0, 1]\times [0, 1]\f$
class PoissonProblems {
    
public:
    /*!
     * @param nr is the number of grid cells in x, y, and z of the coarsest level
     * @param maxGridSize is the max. size of a grid
     * @param nLevels is the max. number of nLevels
     * @param lower boundary of physical domain (applied in all dimensions)
     * @param upper boundary of phyiscal domain (applied in all dimensions)
     */
    PoissonProblems(int nr[3], int maxGridSize, int nLevels,
                    const std::vector<double>& lower,
                    const std::vector<double>& upper);
    
    
    /*!
     * Solves \f$\Delta\phi = -1\f$ only on the grid. In case of
     * nLevels > 0, the refinement is performed on the whole domain.
     * @returns l2 error (single-level vs. multi-level solve)
     */
    double doSolveNoParticles();
    
    /*!
     * Solves \f$\Delta\phi = -1\f$ by initializing
     * particles on the finest level. The charge is scaled
     * such that the rhs is -1 everywhere on the domain. In case of
     * nLevels > 0, the refinement is performed on the whole domain.
     * @returns l2 error (single-level vs. multi-level solve)
     */
    double doSolveParticlesUniform();
    
    /*!
     * Solves \f$\Delta\phi = \rho\f$ where the particles are
     * randomly initialized.
     * @param nParticles to be generated
     * @param mean of the Gaussian distribution
     * @param stddev is the standard deviation of the Gaussian distribution
     * @returns l2 error (single-level vs. multi-level solve)
     */
    double doSolveParticlesGaussian(int nParticles, double mean, double stddev);
    
    /*!
     * Solve the Poisson equation with a real particle distribution
     * read in from a H5 file (Cyclotron)
     * @param step specifies where to read in the H5 file
     * @param h5file specifies the path and filename
     * @returns l2 error (single-level vs. multi-level solve)
     */
    double doSolveParticlesReal(int step, std::string h5file);
    
    /*!
     * Generates 3 bunches with nParticles and same standard deviation
     * @param nParticles is the number of particles per bunch
     * @param stddev is the standard deviation of each bunch
     */
    double doSolveMultiGaussians(int nParticles, double stddev);
    
private:
    void refineWholeDomain_m();     ///< Create refined levels (DistributionMapping and BoxArray)
    void initMultiFabs_m();         ///< Initialize the MultiFab's for solving Poisson with MultiGrid solver
    
private:
    RealBox domain_m;                       ///< Physical domain [0, 1] x [0, 1] x [0, 1]
    int nr_m[3];                            ///< Number of grid cells in each dimension (x, y, z)
    int maxGridSize_m;                      ///< Max. grid size of each level
    int nLevels_m;                          ///< Number of levels
    Array<Geometry> geom_m;                 ///< Geometry of every level
    Array<DistributionMapping> dmap_m;      ///< Distribution to cores of each level
    Array<BoxArray> ba_m;                   ///< All boxes of each level
    Array<int> refRatio_m;                  ///< Refinement ratios among levels (here: 2)
    
    Array<std::unique_ptr<MultiFab> > rho_m;                 ///< Density (i.e. rhs)
    Array<std::unique_ptr<MultiFab> > phi_m;                 ///< Potential
    Array<std::unique_ptr<MultiFab> > efield_m;              ///< Electric field
    Array<std::unique_ptr<MultiFab> > phi_single_m;          ///< Potential for single-level solve
    
};

#endif