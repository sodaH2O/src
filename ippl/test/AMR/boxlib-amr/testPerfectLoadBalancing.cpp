/*!
 * @file testPerfectLoadBalancing.cpp
 * @author Matthias Frey
 * @date November 2016
 * @details Compute \f$\Delta\phi = -\rho\f$ where the charges are
 * distributed in such a way that every cell has the exact amount of particles.
 * Write plot files that can be visualized by yt (python visualize.py)
 * or by AmrVis of the CCSE group at LBNL.
 * 
 * Domain:  [0.0, 1.0] x [0.0, 1.0] x [0.0, 1.0]\n
 * BC:      Dirichlet boundary conditions\n
 * Charge/particle: elementary charge\n
 * Uniform particle distribution per cell
 * 
 * Call:\n
 *  mpirun -np [#cores] testPerfectLoadBalancing [#gridpoints x] [#gridpoints y]
 *                                               [#gridpoints z] [#particles per cell]
 *                                               [#levels] [max. box size]
 * 
 * @brief Computes \f$\Delta\phi = -\rho\f$
 */

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>
#include <iomanip>

#include <ParmParse.H>

#include "PartBunch.h"
#include "AmrPartBunch.h"

#include "Distribution.h"
#include "Solver.h"
#include "AmrOpal.h"

#include "helper_functions.h"

#include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

void doSolve(AmrOpal& myAmrOpal, PartBunchBase* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             int nLevels,
             Inform& msg)
{
    static IpplTimings::TimerRef allocTimer = IpplTimings::getTimer("alloc-memory-grid");
    
    static IpplTimings::TimerRef assignTimer = IpplTimings::getTimer("assign-charge");
    
    // =======================================================================                                                                                                                                   
    // 4. prepare for multi-level solve                                                                                                                                                                          
    // =======================================================================
    
    rhs.resize(nLevels);
    phi.resize(nLevels);
    grad_phi.resize(nLevels);
    
    IpplTimings::startTimer(allocTimer);
    
    for (int lev = 0; lev < nLevels; ++lev) {
        initGridData(rhs, phi, grad_phi, myAmrOpal.boxArray()[lev], lev);
    }
    
    IpplTimings::stopTimer(allocTimer);

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();

    IpplTimings::startTimer(assignTimer);
    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rhs, base_level, 1, finest_level);
    IpplTimings::stopTimer(assignTimer);
    
    double totCharge = totalCharge(rhs, finest_level, geom);
    msg << "Total Charge: " << totCharge << " C" << endl
        << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << vol << " m^3" << endl;
    
    // eps in C / (V * m)
    double constant = -1.0 / Physics::epsilon_0;  // in [V m / C]
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(constant, 0, 1);
#endif
    }
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;

    // solve                                                                                                                                                                                                     
    Solver sol;
    sol.solve_for_accel(rhs,
                        phi,
                        grad_phi,
                        geom,
                        base_level,
                        finest_level,
                        offset);
    
    // for plotting unnormalize
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(1.0 / constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(1.0 / constant, 0, 1);
#endif
    }
}

void doBoxLib(const Vektor<std::size_t, 3>& nr, std::size_t nParticlesPerCell,
              int nLevels, std::size_t maxBoxSize, Inform& msg)
{
    static IpplTimings::TimerRef distTimer = IpplTimings::getTimer("init-distribution");
    static IpplTimings::TimerRef regridTimer = IpplTimings::getTimer("regrid");
    static IpplTimings::TimerRef redistTimer = IpplTimings::getTimer("particle-redistr");
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    std::array<double, 3> lower = {{0.0, 0.0, 0.0}}; // m
    std::array<double, 3> upper = {{1.0, 1.0, 1.0}}; // m
    
    RealBox domain;
    
    // in helper_functions.h
    init(domain, nr, lower, upper);
    
    /*
     * create an Amr object
     */
    ParmParse pp("amr");
    pp.add("max_grid_size", int(maxBoxSize));
    
    Array<int> error_buf(nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr[i];
    
    AmrOpal myAmrOpal(&domain, nLevels - 1, nCells, 0 /* cartesian */);
    
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const Array<BoxArray>& ba = myAmrOpal.boxArray();
    const Array<DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    const Array<Geometry>& geom = myAmrOpal.Geom();
    
    Array<int> rr(nLevels);
    for (int i = 0; i < nLevels; ++i)
        rr[i] = 2;
    
    
    PartBunchBase* bunch = new AmrPartBunch(geom, dmap, ba, rr);
    
    dynamic_cast<AmrPartBunch*>(bunch)->SetAllowParticlesNearBoundary(true);
    
    // initialize a particle distribution
    Distribution dist;
    IpplTimings::startTimer(distTimer);
    dist.uniformPerCell(geom, ba, nr, nParticlesPerCell, ParallelDescriptor::MyProc());
    
    // copy particles to the PartBunchBase object.
    dist.injectBeam(*bunch);
    IpplTimings::stopTimer(distTimer);
    
    
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
        bunch->setQM(Physics::q_e /** 1.0e10*/, i);  // in [C]
    
    // redistribute on single-level
    IpplTimings::startTimer(redistTimer);
    bunch->myUpdate();
    IpplTimings::stopTimer(redistTimer);
    
    bunch->gatherStatistics();
    
    std::size_t nParticles = bunch->getTotalNum();
    msg << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->getQM(0) << " C" << endl
        << "Total charge: " << nParticles * bunch->getQM(0) << " C" << endl;
    
    // ========================================================================
    // 2. tagging (i.e. create BoxArray's, DistributionMapping's of all
    //    other levels)
    // ========================================================================
    
    /*
     * do tagging
     */
    myAmrOpal.setBunch(dynamic_cast<AmrPartBunch*>(bunch));
    
    //    dynamic_cast<AmrPartBunch*>(bunch)->python_format(0);
    
    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    IpplTimings::startTimer(regridTimer);
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    IpplTimings::stopTimer(regridTimer);
    
    bunch->gatherStatistics();
    
    dynamic_cast<AmrPartBunch*>(bunch)->python_format(0);
    
    container_t rhs(PArrayManage);
    container_t phi(PArrayManage);
    container_t grad_phi(PArrayManage);
    
    std::string plotsolve = BoxLib::Concatenate("plt", 0, 4);
    doSolve(myAmrOpal, bunch, rhs, phi, grad_phi, geom, rr, nLevels, msg);
    
    msg << "Total field energy: " << totalFieldEnergy(grad_phi, rr) << endl;
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
#ifdef UNIQUE_PTR
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << grad_phi[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << grad_phi[i]->min(0) << endl;
#else
        msg << "Max. potential level " << i << ": "<< phi[i].max(0) << endl
            << "Min. potential level " << i << ": " << phi[i].min(0) << endl
            << "Max. ex-field level " << i << ": " << grad_phi[i].max(0) << endl
            << "Min. ex-field level " << i << ": " << grad_phi[i].min(0) << endl;
#endif
    }
    
    
    //writePlotFile(plotsolve, rhs, phi, grad_phi, rr, geom, 0);
    
//     dynamic_cast<AmrPartBunch*>(bunch)->python_format(0);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z]"
         << " [#particles per cell] [#levels] [max. box size]";
    
    if ( argc < 7 ) {
        msg << call.str() << endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<std::size_t, 3> nr(std::atoi(argv[1]),
                              std::atoi(argv[2]),
                              std::atoi(argv[3]));
    
    
    std::size_t nParticlesPerCell = std::atoi(argv[4]);
    
    
    msg << "Particle test running with" << endl
        << "- #particles per cell = " << nParticlesPerCell << endl
        << "- grid                = " << nr << endl;
        
    BoxLib::Initialize(argc,argv, false);
    std::size_t nLevels = std::atoi(argv[5]) + 1; // i.e. nLevels = 0 --> only single level
    std::size_t maxBoxSize = std::atoi(argv[6]);
    doBoxLib(nr, nParticlesPerCell, nLevels, maxBoxSize, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    std::stringstream timefile;
    timefile << std::string(argv[0]) << "-timing-cores-"
             << std::setfill('0') << std::setw(6) << Ippl::getNodes()
             << "-threads-1.dat";
    
    IpplTimings::print(timefile.str());
    
    return 0;
}
