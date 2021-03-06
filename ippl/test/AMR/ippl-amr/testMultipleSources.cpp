/*!
 * @file testMultipleSources.cpp
 * @author Matthias Frey
 * @date 21. Nov. 2017
 * @brief Solve the electrostatic potential of 3 Gaussian
 * sources.
 */
#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>

#include <AMReX_ParmParse.H>

#include "../Solver.h"
#include "../MGTSolver.h"
#include "../AmrOpal.h"
#include "../H5Reader.h"

#include "../helper_functions.h"

#include "../writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"
#include <random>

#ifdef HAVE_AMR_MG_SOLVER
    #include "../trilinos/AmrMultiGrid.h"
#endif



#include <getopt.h>

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

struct param_t {
    Vektor<size_t, 3> nr;
    size_t nLevels;
    size_t maxBoxSize;
    size_t blocking_factor;
    double length;
    size_t nParticlesPerBunch;
    double pCharge;
    bool isFixedCharge;
    bool isWriteYt;
    bool isWriteCSV;
    bool isWriteParticles;
    bool isHelp;
    bool useMgtSolver;
    size_t nsolve;
    size_t nbunch;
#ifdef HAVE_AMR_MG_SOLVER
    bool useTrilinos;
    size_t nsweeps;
    std::string bcx;
    std::string bcy;
    std::string bcz;
    std::string bs;
    std::string smoother;
    std::string prec;
    bool rebalance;
#endif
    AmrOpal::TaggingCriteria criteria;
    double tagfactor;
};


bool parseProgOptions(int argc, char* argv[], param_t& params, Inform& msg) {
    /* Parsing Command Line Arguments
     * 
     * 26. June 2017
     * https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
     */
    
    params.isWriteYt = false;
    params.isFixedCharge = false;
    params.isWriteCSV = false;
    params.isWriteParticles = false;
    params.isHelp = false;
    params.useMgtSolver = false;
    params.criteria = AmrOpal::kChargeDensity;
    params.tagfactor = 1.0e-14;
    params.nsolve = 1;
    params.nbunch = 3;
    
#ifdef HAVE_AMR_MG_SOLVER
    params.useTrilinos = false;
    params.nsweeps = 12;
    params.bcx = "dirichlet";
    params.bcy = "dirichlet";
    params.bcz = "dirichlet";
    params.bs = "CG";
    params.smoother = "GS";
    params.prec = "NONE";
    params.rebalance = false;
#endif
    
    
    int c = 0;
    
    int cnt = 0;
    
    int required = 7;
    
    while ( true ) {
        static struct option long_options[] = {
            { "gridx",           required_argument, 0, 'x' },
            { "gridy",           required_argument, 0, 'y' },
            { "gridz",           required_argument, 0, 'z' },
            { "level",           required_argument, 0, 'l' },
            { "maxgrid",         required_argument, 0, 'm' },
            { "blocking_factor", required_argument, 0, 'd' },
            { "npartperbunch",   required_argument, 0, 'n' },
            { "writeYt",         no_argument,       0, 'w' },
            { "help",            no_argument,       0, 'h' },
            { "pcharge",         required_argument, 0, 'c' },
            { "writeCSV",        no_argument,       0, 'v' },
            { "writeParticles",  no_argument,       0, 'p' },
            { "use-mgt-solver",  no_argument,       0, 's' },
            { "nsolve",          required_argument, 0, 'r' },
            { "nbunch",          required_argument, 0, 'B' },
#ifdef HAVE_AMR_MG_SOLVER
            { "use-trilinos",    no_argument,       0, 'a' },
            { "nsweeps",         required_argument, 0, 'g' },
            { "smoother",        required_argument, 0, 'q' },
            { "prec",            required_argument, 0, 'o' },
            { "bcx",             required_argument, 0, 'i' },
            { "bcy",             required_argument, 0, 'j' },
            { "bcz",             required_argument, 0, 'k' },
            { "basesolver",      required_argument, 0, 'u' },
            { "rebalance",      no_argument,       0, 'R' },
#endif
            { "tagging",         required_argument, 0, 't' },
            { "tagging-factor",  required_argument, 0, 'f' },
            { 0,                 0,                 0,  0  }
        };
        
        int option_index = 0;
        
#ifdef HAVE_AMR_MG_SOLVER
        c = getopt_long(argc, argv, "a:B:cd:f:g:hi:j:k:l:m:n:o:pq:r:R:st:u:vwx:y:z:", long_options, &option_index);
#else
        c = getopt_long(argc, argv, "x:y:z:l:m:B:n:whcvpst:f:d:r:", long_options, &option_index);
#endif
        
        if ( c == -1 )
            break;
        
        switch ( c ) {
#ifdef HAVE_AMR_MG_SOLVER
            case 'a':
                params.useTrilinos = true; break;
            case 'g':
                params.nsweeps = std::atoi(optarg); break;
            case 'q':
            {
                params.smoother = optarg;
                break;
            }
            case 'o':
            {
                params.prec = optarg;
                break;
            }
            case 'i':
            {
                params.bcx =  optarg;
                break;
            }
            case 'j':
            {
                params.bcy =  optarg;
                break;
            }
            case 'k':
            {
                params.bcz =  optarg;
                break;
            }
            case 'u':
            {
                params.bs = optarg;
                break;
            }
            case 'R':
            {
                params.rebalance = true;
                break;
            }
#endif
            case 'x':
                params.nr[0] = std::atoi(optarg); ++cnt; break;
            case 'y':
                params.nr[1] = std::atoi(optarg); ++cnt; break;
            case 'z':
                params.nr[2] = std::atoi(optarg); ++cnt; break;
            case 'l':
                params.nLevels = std::atoi(optarg) + 1; ++cnt; break;
            case 'm':
                params.maxBoxSize = std::atoi(optarg); ++cnt; break;
            case 'd':
                params.blocking_factor = std::atoi(optarg); ++cnt; break;
            case 'n':
                params.nParticlesPerBunch = std::atoi(optarg); ++cnt; break;
            case 'c':
                params.pCharge = std::atof(optarg);
                params.isFixedCharge = true;
                break;
            case 'w':
                params.isWriteYt = true;
                break;
            case 'v':
                params.isWriteCSV = true;
                break;
            case 'p':
                params.isWriteParticles = true;
                break;
            case 's':
                params.useMgtSolver = true;
                break;
            case 't':
                if ( std::strcmp("efield", optarg) == 0 )
                    params.criteria = AmrOpal::kEfieldStrength;
                else if ( std::strcmp("potential", optarg) == 0 )
                    params.criteria = AmrOpal::kPotentialStrength;
                break;
            case 'f':
                params.tagfactor = std::atof(optarg);
                break;
            case 'r':
                params.nsolve = std::atoi(optarg); break;
            case 'B':
                params.nbunch = std::atoi(optarg); break;
            case 'h':
                msg << "Usage: " << argv[0]
                    << endl
                    << "--gridx [#gridpoints in x]" << endl
                    << "--gridy [#gridpoints in y]" << endl
                    << "--gridz [#gridpoints in z]" << endl
                    << "--level [#levels]" << endl
                    << "--maxgrid [max. grid]" << endl
                    << "--blocking_factor [val] (only grids modulo bf == 0 allowed)" << endl
                    << "--npartperbunch [#particles per bunch]" << endl
                    << "--pcharge [charge per particle] (optional)" << endl
                    << "--writeYt (optional)" << endl
                    << "--writeCSV (optional)" << endl
                    << "--writeParticles (optional)" << endl
                    << "--use-mgt-solver (optional)" << endl
                    << "--nsolve (optional)" << endl
                    << "--nbunch (optional)" << endl
#ifdef HAVE_AMR_MG_SOLVER
                    << "--use-trilinos (optional)" << endl
                    << "--nsweeps (optional, trilinos only, default: 12)" << endl
                    << "--smoother (optional, trilinos only, default: GAUSS_SEIDEL)" << endl
                    << "--prec (optional, trilinos only, default: NONE)" << endl
                    << "--bcx (optional, dirichlet or open, default: dirichlet)" << endl
                    << "--bcy (optional, dirichlet or open, default: dirichlet)" << endl
                    << "--bcz (optional, dirichlet or open, default: dirichlet)" << endl
                    << "--basesolver (optional, trilinos only, default: CG)" << endl
                    << "--rebalance (optional, trilinos only)" << endl
#endif
                    << "--tagging charge (default) / efield / potential (optional)" << endl
                    << "--tagfactor [charge value / 0 ... 1] (optiona)" << endl;
                params.isHelp = true;
                break;
            case '?':
                break;
            
            default:
                break;
            
        }
    }
    
#ifdef HAVE_AMR_MG_SOLVER
    if ( params.useMgtSolver && params.useTrilinos ) {
        params.useMgtSolver = false;
        msg << "Favouring Trilinos over MGT." << endl;
    }
#endif
    
    return ( cnt == required );
}


void writeCSV(const container_t& phi,
              const container_t& efield,
              double lower, double dx)
{
    // Immediate debug output:
    // Output potential and e-field along axis
    std::string outfile = "potential.grid";
    std::ofstream out;
    for (MFIter mfi(*phi[0]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& lhs = (*phi[0])[mfi];
        
        for (int proc = 0; proc < amrex::ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == amrex::ParallelDescriptor::MyProc() ) {
                
                if ( proc == 0 ) {
                    out.open(outfile, std::ios::out);
                    out << "$x$ [m], $\\Phi$ [V]" << std::endl;
                } else
                    out.open(outfile, std::ios::app);
                
                int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
                int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    amrex::IntVect ivec(i, j, k);
                    // add one in order to have same convention as PartBunch::computeSelfField()
                    out << lower + i * dx << ", " << lhs(ivec, 0)  << std::endl;
                }
                out.close();
            }
            amrex::ParallelDescriptor::Barrier();
        }
    }
    
    outfile = "efield.grid";
    for (MFIter mfi(*efield[0]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& lhs = (*efield[0])[mfi];
        
        for (int proc = 0; proc < amrex::ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == amrex::ParallelDescriptor::MyProc() ) {
                
                if ( proc == 0 ) {
                    out.open(outfile, std::ios::out);
                    out << "$x$ [m], $E_x$ [V/m], $E_x$ [V/m], $E_x$ [V/m]" << std::endl;
                } else
                    out.open(outfile, std::ios::app);
                
                int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
                int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    amrex::IntVect ivec(i, j, k);
                    // add one in order to have same convention as PartBunch::computeSelfField()
                    out << lower + i * dx << ", " << lhs(ivec, 0) << ", "
                        << lhs(ivec, 1) << ", " << lhs(ivec, 2) << std::endl;
                }
                out.close();
            }
            amrex::ParallelDescriptor::Barrier();
        }
    }
    
}


void writeYt(container_t& rho,
             const container_t& phi,
             const container_t& efield,
             const amrex::Array<amrex::Geometry>& geom,
             const amrex::Array<int>& rr,
             const double& scalefactor,
             const param_t& params)
{
    std::string dir = "yt-testMultipleSources";
    
    double time = 0.0;
    
    double fac = 1.0;

    if ( !params.isFixedCharge )
        fac = Physics::epsilon_0;

    for (unsigned int i = 0; i < rho.size(); ++i)
        rho[i]->mult(- fac / scalefactor, 0, 1);
    
    writePlotFile(dir, rho, phi, efield, rr, geom, time, scalefactor);
}


void randomMove(amrbunch_t* bunch, int seed, Inform& msg) {
    /* We move every particle in the range of 1 mm in original space
     */
    
    
    std::mt19937_64 mt(seed);
    
    std::uniform_real_distribution<> dist(-0.001, 0.001);
    
    msg << endl << "Move particles randomly" << endl;
    
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
        for (int d = 0; d < 3; ++d)
            bunch->R[i](d) += dist(mt);
    }
    
    double scale = 1.0;
    scale = domainMapping(*bunch, scale);
    
    msg << "Transformed positions (scale = " << scale << ")" << endl
        << "Update particle-to-core" << endl;
    
    bunch->update();
    
    domainMapping(*bunch, scale, true);
    
    msg << "Back to normal positions" << endl;
}


double depositCharge(AmrOpal& myAmrOpal, amrbunch_t* bunch,
                     container_t& rhs, Inform& msg, const param_t& params, double& scale)
{
    scale = 1.0;
    
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();
    
    container_t partMF(params.nLevels);
    for (uint lev = 0; lev < params.nLevels; lev++) {
        rhs[lev]->setVal(0.0);
        
        const amrex::BoxArray& ba = myAmrOpal.boxArray()[lev];
        const amrex::DistributionMapping& dmap = myAmrOpal.DistributionMap(lev);
        partMF[lev].reset(new amrex::MultiFab(ba, dmap, 1, 2));
        partMF[lev]->setVal(0.0, 2);
    }
    
    scale = domainMapping(*bunch, scale);
    
    msg << endl << "Transformed positions (scale = " << scale << ")" << endl
        << "Assign density" << endl;
    
    bunch->AssignDensityFort(bunch->qm, partMF, base_level, 1, finest_level);
    
    domainMapping(*bunch, scale, true);
    
    msg << "Back to normal positions" << endl;
    
    for (uint lev = 0; lev < params.nLevels; ++lev) {
        amrex::MultiFab::Copy(*rhs[lev], *partMF[lev], 0, 0, 1, 0);
    }
    
    for (uint i = 0; i < rhs.size(); ++i)
        if ( rhs[i]->contains_nan() || rhs[i]->contains_inf() )
            throw std::runtime_error("\033[1;31mError: NANs or INFs on charge grid.\033[0m");
    
    const amrex::Array<amrex::Geometry>& geom = myAmrOpal.Geom();
    
    amrex::Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    msg << "Cell volume: " << *(geom[0].CellSize()) << "^3 = " << vol << " m^3" << endl;
    
    // Check charge conservation
    double totCharge = totalCharge(rhs, finest_level, geom);
    
    msg << "Total Charge (computed): " << totCharge << " C" << endl
        << "Vacuum permittivity: " << Physics::epsilon_0 << " F/m (= C/(m V)" << endl;
    
    // eps in C / (V * m)
    double constant = -1.0;

    if ( !params.isFixedCharge )
        constant /= Physics::epsilon_0 ; //* scale;  // in [V m / C]

    for (int i = 0; i <= finest_level; ++i) {
        rhs[i]->mult(constant, 0, 1);       // in [V m]
    }
    
    
    // normalize each level
    double l0norm = rhs[finest_level]->norm0(0);
    msg << "l0norm = " << l0norm << endl;
    for (int i = 0; i <= finest_level; ++i)
        rhs[i]->mult(1.0 / l0norm, 0, 1);
    return l0norm;
}


void prepareSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
                  container_t& rhs,
                  container_t& phi,
                  container_t& efield,
                  const amrex::Array<int>& rr,
                  const param_t& params)
{
    // =======================================================================
    // 4. prepare for multi-level solve               
    // =======================================================================
    rhs.resize(params.nLevels);
    phi.resize(params.nLevels);
    efield.resize(params.nLevels);
    

    // Define the density on level 0 from all particles at all levels
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();
    
    for (uint lev = 0; lev < params.nLevels; ++lev) {
        initGridData(rhs, phi, efield,
                     myAmrOpal.boxArray()[lev], myAmrOpal.DistributionMap(lev), lev);
    }
}


void print(AmrOpal& myAmrOpal, const double& scale,
           container_t& rhs, container_t& phi,
           container_t& efield, Inform& msg,
           const amrex::Array<int>& rr,
           const param_t& params)
{
    const amrex::Array<amrex::Geometry>& geom = myAmrOpal.Geom();
    
    for (int i = 0; i <= myAmrOpal.finestLevel(); ++i) {
        msg << "Max. potential level " << i << ": "<< phi[i]->max(0) << endl
            << "Min. potential level " << i << ": " << phi[i]->min(0) << endl
            << "Max. ex-field level " << i << ": " << efield[i]->max(0) << endl
            << "Min. ex-field level " << i << ": " << efield[i]->min(0) << endl;
    }
    
    double fieldenergy = totalFieldEnergy(efield, rr);
    
    msg << "Total field energy: " << fieldenergy << endl;
    
    amrex::RealBox amr_domain = geom[0].ProbDomain();
    if (params.isWriteCSV && Ippl::getNodes() == 1 && myAmrOpal.maxGridSize(0) == (int)params.nr[0] )
        writeCSV(phi, efield, amr_domain.lo(0) / scale, geom[0].CellSize(0) / scale);
    
    if ( params.isWriteYt )
        writeYt(rhs, phi, efield, geom, rr, scale, params);
}


void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& efield,
             const amrex::Array<int>& rr,
             Inform& msg, const param_t& params)
{
    static IpplTimings::TimerRef solvTimer = IpplTimings::getTimer("solve");
    
    int maxiter = 100;
    int maxiter_b = 100;
    int verbose = 0;
    bool usecg = true;
    double bottom_solver_eps = 1.0e-4;
    int max_nlevel = 1024;
    
    /* MG_SMOOTHER_GS_RB  = 1
     * MG_SMOOTHER_JACOBI = 2
     * MG_SMOOTHER_MINION_CROSS = 5
     * MG_SMOOTHER_MINION_FULL = 6
     * MG_SMOOTHER_EFF_RB = 7
     */
    int smoother = 1;
    
    // #smoothings at each level on the way DOWN the V-cycle
    int nu_1 = 2;
    
    // #smoothings at each level on the way UP the V-cycle
    int nu_2 = 2;
    
    // #smoothings before and after the bottom solver
    int nu_b = 0;
    
    // #smoothings
    int nu_f = 8;

    /* MG_FCycle = 1  (full multigrid)
     * MG_WCycle = 2
     * MG_VCycle = 3
     * MG_FVCycle = 4
     */
    int cycle = 1;
    
    bool cg_solver = true;
    
    /* if cg_solver == true:
     * - BiCG --> 1
     * - CG --> 2
     * - CABiCG --> 3
     * 
     * else if cg_solver == false
     * - CABiCG is taken
     */
    int bottom_solver = 1;
    
    amrex::ParmParse pp("mg");

    pp.add("maxiter", maxiter);
    pp.add("maxiter_b", maxiter_b);
    pp.add("nu_1", nu_1);
    pp.add("nu_2", nu_2);
    pp.add("nu_b", nu_b);
    pp.add("nu_f", nu_f);
    pp.add("v"   , verbose);
    pp.add("usecg", usecg);
    pp.add("cg_solver", cg_solver);

    pp.add("rtol_b", bottom_solver_eps);
    pp.add("numLevelsMAX", max_nlevel);
    pp.add("smoother", smoother);
    pp.add("cycle_type", cycle); // 1 -> F, 2 -> W, 3 -> V, 4 -> F+V
    //
    // The C++ code usually sets CG solver type using cg.cg_solver.
    // We'll allow people to also use mg.cg_solver but pick up the former as well.
    //
    if (!pp.query("cg_solver", cg_solver))
    {
        amrex::ParmParse pp("cg");

        pp.add("cg_solver", cg_solver);
    }

    pp.add("bottom_solver", bottom_solver);
    
    
    prepareSolve(myAmrOpal, bunch, rhs, phi, efield, rr, params);
    
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();
    const amrex::Array<amrex::Geometry>& geom = myAmrOpal.Geom();
    
    double scale = 1.0;
    double l0norm = depositCharge(myAmrOpal, bunch, rhs, msg, params, scale);
    
    // **************************************************************************
    // Compute the total charge of all particles in order to compute the offset
    //     to make the Poisson equations solvable
    // **************************************************************************
    amrex::Real offset = 0.;

    // solve
#ifdef HAVE_AMR_MG_SOLVER
    if ( params.useTrilinos ) {
	std::string interp = "PC";
        std::string norm = "LINF";
        AmrMultiGrid sol(&myAmrOpal, params.bs, params.prec,
                         params.rebalance, params.bcx, params.bcy,
                         params.bcz, params.smoother, params.nsweeps,
                         interp, norm);
        
        for (uint i = 0; i < params.nsolve; ++i) {
            
            IpplTimings::startTimer(solvTimer);
            
            sol.solve(rhs,            // [V m]
                      phi,            // [V m^3]
                      efield,       // [V m^2]
                      base_level,
                      finest_level, i);
            
            IpplTimings::stopTimer(solvTimer);
            
            // undo normalization
            for (int j = 0; j <= finest_level; ++j) {
                phi[j]->mult(scale * l0norm, 0, 1);
            }
            
            // undo scale
            for (int j = 0; j <= finest_level; ++j)
                efield[j]->mult(scale * scale * l0norm, 0, 3);
            
            print(myAmrOpal, scale, rhs, phi, efield, msg, rr, params);
            
            msg << "#iterations: " << sol.getNumIters() << endl;
            
            for (int j = 0; j <= finest_level; ++j) {
                msg << "norm of residual (level " << j << "): "
                    << sol.getLevelResidualNorm(j) << endl;
            }
            
            if ( params.nsolve > 1 ) {
                randomMove(bunch, Ippl::myNode() + i, msg);
                
                l0norm = depositCharge(myAmrOpal, bunch, rhs, msg, params, scale);
            }
        }
    } else if ( params.useMgtSolver ) {
#else
    if ( params.useMgtSolver ) {
#endif
        MGTSolver sol;
        
        for (uint i = 0; i < params.nsolve; ++i) {
            IpplTimings::startTimer(solvTimer);
            
            sol.solve(rhs,            // [V m]
                      phi,            // [V m^3]
                      efield,       // [V m^2]
                      geom);
            
            IpplTimings::stopTimer(solvTimer);
            
            // undo normalization
            for (int i = 0; i <= finest_level; ++i) {
                phi[i]->mult(scale * l0norm, 0, 1);
            }
            
            // undo scale
            for (int i = 0; i <= finest_level; ++i)
                efield[i]->mult(scale * scale * l0norm, 0, 3);
            
            print(myAmrOpal, scale, rhs, phi, efield, msg, rr, params);
            
            if ( params.nsolve > 1 ) {
                randomMove(bunch, Ippl::myNode() + i, msg);
                
                l0norm = depositCharge(myAmrOpal, bunch, rhs, msg, params, scale);
            }
        }
    } else {
        Solver sol;
        
        for (uint i = 0; i < params.nsolve; ++i) {
            IpplTimings::startTimer(solvTimer);
            
            sol.solve_for_accel(rhs,            // [V m]
                                phi,            // [V m^3]
                                efield,       // [V m^2]
                                geom,
                                base_level,
                                finest_level,
                                offset);
            
            IpplTimings::stopTimer(solvTimer);
            
            // undo normalization
            for (int i = 0; i <= finest_level; ++i) {
                phi[i]->mult(scale * l0norm, 0, 1);
            }
            
            // undo scale
            for (int i = 0; i <= finest_level; ++i)
                efield[i]->mult(scale * scale * l0norm, 0, 3);
            
            print(myAmrOpal, scale, rhs, phi, efield, msg, rr, params);
            
            if ( params.nsolve > 1 ) {
                randomMove(bunch, Ippl::myNode() + i, msg);
                
                l0norm = depositCharge(myAmrOpal, bunch, rhs, msg, params, scale);
            }
        }
    }
}

void doAMReX(const param_t& params, Inform& msg)
{
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    /*
     * create an Amr object
     */
    amrex::ParmParse pp("amr");
    pp.add("max_grid_size", int(params.maxBoxSize));
    
    amrex::Array<int> error_buf(params.nLevels, 0);
    
    amrex::Array<int> bf(params.nLevels, int(params.blocking_factor));
    pp.addarr("blocking_factor", bf);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    amrex::ParmParse pgeom("geometry");
    amrex::Array<int> is_per = { 0, 0, 0};
    pgeom.addarr("is_periodic", is_per);
    
    amrex::Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = params.nr[i];
    
    
    std::vector<int> rr(params.nLevels);
    amrex::Array<int> rrr(params.nLevels);
    for (unsigned int i = 0; i < params.nLevels; ++i) {
        rr[i] = 2;
        rrr[i] = 2;
    }
    
    amrex::RealBox amr_domain;
    
    std::array<double, AMREX_SPACEDIM> amr_lower = {{-1.04, -1.04, -1.04}}; // m
    std::array<double, AMREX_SPACEDIM> amr_upper = {{ 1.04,  1.04,  1.04}}; // m
    
    init(amr_domain, params.nr, amr_lower, amr_upper);
    
    AmrOpal myAmrOpal(&amr_domain, params.nLevels - 1, nCells, 0 /* cartesian */, rr);
    
    myAmrOpal.setTagging(params.criteria);
    
    if ( params.criteria == AmrOpal::kChargeDensity )
        myAmrOpal.setCharge(params.tagfactor);
    else
        myAmrOpal.setScalingFactor(params.tagfactor);
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const amrex::Array<amrex::BoxArray>& ba = myAmrOpal.boxArray();
    const amrex::Array<amrex::DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    amrex::Array<amrex::Geometry>& geom = myAmrOpal.Geom();
    
    
    amrplayout_t* playout = new amrplayout_t(geom, dmap, ba, rrr);
    
    std::unique_ptr<amrbunch_t> bunch( new amrbunch_t() );
    bunch->initialize(playout);
    bunch->initializeAmr(); // add attributes: level, grid
    
    bunch->setAllowParticlesNearBoundary(true);
    
    // map particles
    double scale = 1.0;
    
    // initialize particle distribution
    unsigned long int nloc = params.nParticlesPerBunch / Ippl::getNodes();
    Distribution dist;
    
    
    double stddev[] = {0.0025, 0.0025, 0.0025};
    double d = 0.04 /* [m] */ / std::sqrt(2.0);
    for (size_t b = 0; b < params.nbunch; ++b) {
        // first source
        double mean[] = {d * b, d * b, 0};
        
        dist.gaussian(&mean[0], &stddev[0], nloc, Ippl::myNode());
        dist.injectBeam(*bunch);
        
        scale = domainMapping(*bunch, scale);
        
        msg << endl << "Transformed positions (scale = " << scale << ")" << endl << endl;
        
        bunch->update();
        
        domainMapping(*bunch, scale, true);
        
        msg << endl << "Back to normal positions" << endl << endl;
    }
    
    
    if ( !params.isFixedCharge )
        for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
            bunch->qm[i] = Physics::q_e;  // in [C]
    else
        for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
            bunch->qm[i] = params.pCharge;
    
    if ( params.isWriteParticles ) {
        H5Reader h5("testMultipleSources.h5");
        h5.open(0, H5_O_WRONLY);
        h5.write(bunch.get());
        h5.close();
    }
    
    msg << "#Particles: " << bunch->getTotalNum() << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << bunch->getTotalNum() * bunch->qm[0] << " C" << endl;
    
    
    scale = domainMapping(*bunch, scale);
    
    msg << "Scale: " << scale << endl;
    
    myAmrOpal.setBunch(bunch.get());
    
    msg << "//\n//  Process Statistic BEFORE regrid\n//" << endl;
    
    bunch->gatherStatistics();

    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    
    container_t rhs;
    container_t phi;
    container_t efield;
    
    msg << endl << "Transformed positions (scale = " << scale << ")" << endl << endl;
    
    msg << "#Particles: " << bunch->getTotalNum() << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << bunch->getTotalNum() * bunch->qm[0] << " C" << endl;
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, scale/*0.0*/ /*time*/);
    
    
    msg << endl << "Back to normal positions" << endl << endl;
    
    domainMapping(*bunch, scale, true);
    

    if ( Ippl::myNode() == 0 ) {
        std::ofstream out("boxes-per-level-ncores-" + std::to_string(Ippl::getNodes()) + ".dat");
        for (int i = 0; i <= myAmrOpal.finestLevel(); ++i)
            out << i << " " << myAmrOpal.boxArray(i).size() << std::endl;
        out.close();
    }
    
    bunch->gatherLevelStatistics();
    
    msg << "//\n//  Process Statistic AFTER regrid\n//" << endl;
    
    bunch->gatherStatistics();
    
    static IpplTimings::TimerRef statisticsTimer = IpplTimings::getTimer("dump-statistics");
    std::string statistics = "particle-statistics-ncores-" + std::to_string(Ippl::getNodes()) + ".dat";
    IpplTimings::startTimer(statisticsTimer);
    bunch->dumpStatistics(statistics);
    IpplTimings::stopTimer(statisticsTimer);
    
    doSolve(myAmrOpal, bunch.get(), rhs, phi, efield, rrr, msg, params);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg(argv[0]);
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);
    
    param_t params;
    
    amrex::Initialize(argc,argv, false);
    
    try {
        if ( !parseProgOptions(argc, argv, params, msg) && !params.isHelp )
            throw std::runtime_error("\033[1;31mError: Check the program options.\033[0m");
        else if ( params.isHelp )
            return 0;
        
        
        std::string tagging = "charge";
        if ( params.criteria == AmrOpal::kEfieldStrength )
            tagging = "efield";
        else if ( params.criteria == AmrOpal::kPotentialStrength )
            tagging = "potential";
    
        msg << "Particle test running with" << endl
            << "- grid                  = " << params.nr << endl
            << "- max. grid             = " << params.maxBoxSize << endl
            << "- blocking factor       = " << params.blocking_factor << endl
            << "- #level                = " << params.nLevels - 1 << endl
            << "- #particles per bunch  = " << params.nParticlesPerBunch << endl
            << "- tagging               = " << tagging << endl
            << "- tagging factor        = " << params.tagfactor << endl
            << "- #solves               = " << params.nsolve << endl
            << "- #bunches              = " << params.nbunch << endl;

        if ( params.isFixedCharge )
            msg << "- charge per particle   = " << params.pCharge << endl;
        
        if ( params.useMgtSolver )
            msg << "- MGT solver is used" << endl;
        
#ifdef HAVE_AMR_MG_SOLVER
        if ( params.useTrilinos ) {
            msg << "- Trilinos solver is used with: "
                << endl
                << "    - nsweeps:     " << params.nsweeps
                << endl
                << "    - smoother:    " << params.smoother
                << endl;
        }
#endif
        doAMReX(params, msg);
        
    } catch(const std::exception& ex) {
        msg << ex.what() << endl;
    }
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    return 0;
}
