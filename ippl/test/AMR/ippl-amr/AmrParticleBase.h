#ifndef AMRPARTBASE_H
#define AMRPARTBASE_H

#include "Ippl.h"

#include <map>
#include <deque>
#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <array>

#include <AMReX_ParmParse.H>

#include <AMReX_ParGDB.H>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_Particles.H>
#include <AMReX_RealBox.H>


#include <AMReX_BLFort.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFabUtil_F.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>

#include "LevelNumCounter.h"


//AMRPArticleBase class definition. Template parameter is the specific AmrParticleLayout-derived
//class which determines how the particles are distribute amoung processors.
//AmrParticleBase is derived from IpplParticle base and 

typedef double RealType;

typedef long                    SortListIndex_t;
typedef std::vector<SortListIndex_t> SortList_t;
typedef std::vector<ParticleAttribBase *>      attrib_container_t;

typedef std::deque<amrex::ParticleCommData> C;
typedef std::deque<amrex::ParticleCommData> PBox;
typedef amrex::ParticleCommData ParticleType;
typedef typename std::map<int,PBox> PMap;

template<class PLayout>
class AmrParticleBase : public IpplParticleBase<PLayout> {
    
public:
     
    typedef typename PLayout::ParticlePos_t   ParticlePos_t;
    typedef typename PLayout::ParticleIndex_t ParticleIndex_t;
    typedef typename PLayout::SingleParticlePos_t SingleParticlePos_t;
    typedef LevelNumCounter<size_t, size_t>   LevelNumCounter_t;

    ParticleIndex_t m_lev;
    ParticleIndex_t m_grid;

 private:

    IpplTimings::TimerRef AssignDensityTimer_m;
    IpplTimings::TimerRef SortParticlesTimer_m;
    IpplTimings::TimerRef UpdateParticlesTimer_m;

    bool allow_particles_near_boundary;
    
    LevelNumCounter_t LocalNumPerLevel_m;

    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    static void CIC_Cells_Fracs_Basic (const SingleParticlePos_t &R, const amrex::Real* plo, 
                                       const amrex::Real* dx, amrex::Real* fracs,  amrex::IntVect* cells);

    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    static int CIC_Cells_Fracs (const SingleParticlePos_t &R,
                                const amrex::Real*         plo,
                                const amrex::Real*         dx_geom,
                                const amrex::Real*         dx_part,
                                amrex::Array<amrex::Real>&        fracs,
                                amrex::Array<amrex::IntVect>&     cells);
    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    bool FineToCrse (const int ip,
                     int                                flev,
                     const amrex::Array<amrex::IntVect>&              fcells,
                     const amrex::BoxArray&                    fvalid,
                     const amrex::BoxArray&                    compfvalid_grown,
                     amrex::Array<amrex::IntVect>&                    ccells,
                     amrex::Array<amrex::Real>&                       cfracs,
                     amrex::Array<int>&                        which,
                     amrex::Array<int>&                        cgrid,
                     amrex::Array<amrex::IntVect>&                    pshifts,
                     std::vector< std::pair<int,amrex::Box> >& isects);

    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    void FineCellsToUpdateFromCrse (const int ip,
                                    int lev,
                                    const amrex::IntVect& ccell,
                                    const amrex::IntVect& cshift,
                                    amrex::Array<int>& fgrid,
                                    amrex::Array<amrex::Real>& ffrac,
                                    amrex::Array<amrex::IntVect>& fcells,
                                    std::vector< std::pair<int,amrex::Box> >& isects);

    //Function from AMReX adjusted to work with Ippl AmrParticleBase class
    //sends/receivs the particles that are needed by other processes to during AssignDensity
    void AssignDensityDoit(int level, amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, PMap& data,
                           int ncomp, int lev_min = 0);

    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    // Assign values from grid back to particles
    void Interp(const SingleParticlePos_t &R, const amrex::Geometry &geom, const amrex::FArrayBox& fab, 
                const int* idx, amrex::Real* val, int cnt);
    
    
    
    
    
    // amrex repository AMReX_MultiFabUtil.H (missing in AMReX repository)
    void sum_fine_to_coarse(/*const */amrex::MultiFab& S_fine, amrex::MultiFab& S_crse,
                            int scomp, int ncomp, const amrex::IntVect& ratio,
                            const amrex::Geometry& cgeom, const amrex::Geometry& fgeom) const
    {
        BL_ASSERT(S_crse.nComp() == S_fine.nComp());
        BL_ASSERT(ratio == ratio[0]);
        BL_ASSERT(S_fine.nGrow() % ratio[0] == 0);

        const int nGrow = S_fine.nGrow() / ratio[0];

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        amrex::BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);

        amrex::MultiFab crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(), ncomp, nGrow);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (amrex::MFIter mfi(crse_S_fine, true); mfi.isValid(); ++mfi)
        {
            //  NOTE: The tilebox is defined at the coarse level.
            const amrex::Box& tbx = mfi.growntilebox(nGrow);
            
            BL_FORT_PROC_CALL(BL_AVGDOWN, bl_avgdown)
                (tbx.loVect(), tbx.hiVect(),
                 BL_TO_FORTRAN_N(S_fine[mfi] , scomp),
                 BL_TO_FORTRAN_N(crse_S_fine[mfi], 0),
                 ratio.getVect(),&ncomp);
        }
        
        S_crse.copy(crse_S_fine, 0, scomp, ncomp, nGrow, 0,
                    cgeom.periodicity(), amrex::FabArrayBase::ADD);
    }
  
public: 

    //constructor: initializes timers and default variables
    AmrParticleBase() : allow_particles_near_boundary(false), LocalNumPerLevel_m() { 
        AssignDensityTimer_m = IpplTimings::getTimer("AMR AssignDensity");
        SortParticlesTimer_m = IpplTimings::getTimer("AMR sort particles");
        UpdateParticlesTimer_m = IpplTimings::getTimer("AMR update particles");
    }

    // destructor - delete the layout if necessary
    ~AmrParticleBase() { }

    //initialize AmrParticleBase class - add level and grid variables to attribute list
    void initializeAmr() {
        this->addAttribute(m_lev);
        this->addAttribute(m_grid);
    }

    void setAllowParticlesNearBoundary(bool value) {
        allow_particles_near_boundary = value;
    }
    
    void setLocalNumPerLevel(const LevelNumCounter_t& LocalNumPerLevel) {
        LocalNumPerLevel_m = LocalNumPerLevel_m;
    }
    
    LevelNumCounter_t& getLocalNumPerLevel() {
        return LocalNumPerLevel_m;
    }
    
    const LevelNumCounter_t& getLocalNumPerLevel() const {
        return LocalNumPerLevel_m;
    }
    
    void destroy(size_t M, size_t I, bool doNow = false) {
        /* if the particles are deleted directly
         * we need to update the particle level count
         */
        if ( doNow ) {
            for (size_t ip = I; ip < M + I; ++ip)
                --LocalNumPerLevel_m[ m_lev[ip] ];
        }
        IpplParticleBase<PLayout>::destroy(M, I, doNow);
    }

    void performDestroy(bool updateLocalNum = false) {
        // nothing to do if destroy list is empty
        if ( this->DestroyList.empty() )
            return;
        
        if ( updateLocalNum ) {
            typedef std::vector< std::pair<size_t,size_t> > dlist_t;
            dlist_t::const_iterator curr = this->DestroyList.begin();
            const dlist_t::const_iterator last = this->DestroyList.end();
            
            while ( curr != last ) {
                
                for (size_t ip = curr->first;
                     ip < curr->first + curr->second;
                     ++ip)
                {
                    --LocalNumPerLevel_m[ m_lev[ip] ];
                }
            
                ++curr;
            }
        }
        
        IpplParticleBase<PLayout>::performDestroy(updateLocalNum);
    }
    
    void create(size_t M) {
        // particles are created at the coarsest level
        LocalNumPerLevel_m[0] += M;
        
        IpplParticleBase<PLayout>::create(M);
    }
    
    void createWithID(unsigned id) {
        ++LocalNumPerLevel_m[0];
        
        IpplParticleBase<PLayout>::createWithID(id);
    }
    
    
    
    // Update the particle object after a timestep.  This routine will change
    // our local, total, create particle counts properly.
    void update(int lbase = 0, int lfine = -1) {
    
        IpplTimings::startTimer(UpdateParticlesTimer_m);

        // make sure we've been initialized
        PLayout *Layout = &this->getLayout();

        PAssert(Layout != 0);
    
        // ask the layout manager to update our atoms, etc.
        Layout->update(*this, 0, lbase, lfine);
        //sort the particles by grid and level
        sort();
        INCIPPLSTAT(incParticleUpdates);
    
        IpplTimings::stopTimer(UpdateParticlesTimer_m);

    }

    // Update the particle object after a timestep.  This routine will change
    // our local, total, create particle counts properly.
    void update(const ParticleAttrib<char>& canSwap) {

        IpplTimings::startTimer(UpdateParticlesTimer_m);

        // make sure we've been initialized
        PLayout *Layout = &this->getLayout();
        PAssert(Layout != 0);
    
        // ask the layout manager to update our atoms, etc.
        Layout->update(*this, &canSwap);
        //sort the particles by grid and level
        sort();
        INCIPPLSTAT(incParticleUpdates);
    
        IpplTimings::stopTimer(UpdateParticlesTimer_m);
    }

    //sort particles based on the grid and level that they belong to
    void sort() {

        IpplTimings::startTimer(SortParticlesTimer_m);
        size_t LocalNum = this->getLocalNum();
        SortList_t slist1(LocalNum); //slist1 holds the index of where each element should go
        SortList_t slist2(LocalNum); //slist2 holds the index of which element should go in this position

        //sort the lists by grid and level
        //slist1 hold the index of where each element should go in the list
        std::iota(slist1.begin(), slist1.end(), 0);
        std::sort(slist1.begin(), slist1.end(), [this](const SortListIndex_t &i, 
                                                       const SortListIndex_t &j) 
                  {
                      return (this->m_lev[i] < this->m_lev[j] || 
                              (this->m_lev[i] == this->m_lev[j] && this->m_grid[i] < this->m_grid[j]));
                  });

        //slist2 holds the index of which element should go in this position
        for (unsigned int i = 0; i < LocalNum; ++i)
            slist2[slist1[i]] = i;

        //sort the array according to slist2
        this->sort(slist2);

        IpplTimings::stopTimer(SortParticlesTimer_m);
    }

    //sort the particles given a sortlist
    void sort(SortList_t &sortlist) {
        attrib_container_t::iterator abeg = this->begin();
        attrib_container_t::iterator aend = this->end();
        for ( ; abeg != aend; ++abeg )
            (*abeg)->sort(sortlist);
    }
    
    template <class AType>
    void AssignDensityFort (ParticleAttrib<AType> &pa,
                            amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf_to_be_filled, 
                            int lev_min, int ncomp, int finest_level) const;
    
    template <class AType>
    void InterpolateFort (ParticleAttrib<AType> &pa,
                          amrex::Array<std::unique_ptr<amrex::MultiFab> >& mesh_data, 
                          int lev_min, int lev_max);
    
    template <class AType>
    void InterpolateSingleLevelFort (ParticleAttrib<AType> &pa, amrex::MultiFab& mesh_data, int lev);
    
    template <class AType>
    void AssignCellDensitySingleLevelFort (ParticleAttrib<AType> &pa, amrex::MultiFab& mf, int level,
					   int ncomp=1, int particle_lvl_offset = 0) const;


    //Function from AMReX adjusted to work with Ippl AmrParticleBase class
    //Scatter the particle attribute pa on the grid 
    template <class AType>
        void AssignDensity(ParticleAttrib<AType> &pa,
                           bool sub_cycle,
                           amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf_to_be_filled,
                           int lev_min,
                           int finest_level)
        {

            IpplTimings::startTimer(AssignDensityTimer_m);

            PLayout *Layout = &this->getLayout();
            const amrex::ParGDBBase* m_gdb = Layout->GetParGDB();
            size_t LocalNum = this->getLocalNum();

            //TODO: lev min should be > 0 if there are no particles from level 0 in the node
            if (lev_min < 0 || lev_min > finest_level)
                lev_min = 0;

            //if lev_min < m_lev[0] there are no particles from lev_min on the node
            //move the lev_min to the first level that is present on the node
            //if ((unsigned)lev_min < m_lev[0])
            //  lev_min = m_lev[0];

            //while lev_min > m_lev[start_idx] we need to skip these particles since there level is
            //higher than the specified lev_min
            int start_idx = 0;
            while ((unsigned)lev_min > m_lev[start_idx])
                start_idx++;

            if (finest_level == -1)
                finest_level = m_gdb->finestLevel();

            while (!m_gdb->LevelDefined(finest_level))
                finest_level--;

            //
            // The size of the returned multifab is limited by lev_min and 
            // finest_level. In the following code, lev is the real level, 
            // lev_index is the corresponding index for mf. 
            //

            // Create the space for mf_to_be_filled, regardless of whether we'll need a temporary mf
            mf_to_be_filled.resize(finest_level+1-lev_min);
            for (int lev = lev_min; lev <= finest_level; lev++)
            { 
                const int lev_index = lev - lev_min;
                mf_to_be_filled[lev_index].reset(new amrex::MultiFab(m_gdb->boxArray(lev),
                                                              m_gdb->DistributionMap(lev),
                                                              1, 1));
                mf_to_be_filled[lev_index]->setVal(0.0);
            }

            // Test whether the grid structure of the boxArray is the same
            //       as the ParticleBoxArray at all levels 
            bool all_grids_the_same = true; 
            for (int lev = lev_min; lev <= finest_level; lev++) {
                if (!m_gdb->OnSameGrids(lev, *mf_to_be_filled[lev-lev_min])) {
                    all_grids_the_same = false;
                    break;
                }
            }
            
            amrex::Array<std::unique_ptr<amrex::MultiFab> > mf_part;
            if (!all_grids_the_same)
            { 
                // Create the space for the temporary, mf_part
                mf_part.resize(finest_level+1-lev_min);
                for (int lev = lev_min; lev <= finest_level; lev++)
                {
                    const int lev_index = lev - lev_min;
                    mf_part[lev_index].reset(new amrex::MultiFab(m_gdb->ParticleBoxArray(lev),
                                                          m_gdb->ParticleDistributionMap(lev),
                                                          1, 1));
                    mf_part[lev_index]->setVal(0.0);
                }
            }
            
            auto & mf = (all_grids_the_same) ? mf_to_be_filled : mf_part;
            
            if (finest_level == 0)
            {
                //
                // Just use the far simpler single-level version.
                //
                AssignDensitySingleLevel(pa, *mf[0], 0, 1);
                //
                // I believe that we don't need any information in ghost cells so we don't copy those.
                //
                if ( ! all_grids_the_same) {
                    mf_to_be_filled[0]->copy(*mf[0],0,0,1);
                }
                return;
            }

            //
            // This'll hold all the info I need for parallel.
            // // What I'll use: m_lev, m_grid, m_cell & m_data[0..ncomp-1].
            //
            // This is the "data" needed by other MPI procs.
            //
            PMap data;

            //
            // Minimum M required.
            //
            const int M = D_TERM(2,+2,+4);

            amrex::Array<int>     cgrid(M);
            amrex::Array<int>    cwhich(M),  fwhich(M);
            amrex::Array<amrex::Real>    fracs(M),  cfracs(M);
            amrex::Array<amrex::IntVect> cells(M),  ccells(M), cfshifts(M);

            ParticleType pb;

            //
            // I'm going to allocate these badboys here & pass'm into routines that use'm.
            // This should greatly cut down on memory allocation/deallocation.
            //
            amrex::Array<amrex::IntVect>                    pshifts(27);
            std::vector< std::pair<int,amrex::Box> > isects;
            amrex::Array<int>                        fgrid(M);
            amrex::Array<amrex::Real>                       ffracs(M);
            amrex::Array<amrex::IntVect>                    fcells;
            //
            // "fvalid" contains all the valid region of the MultiFab at this level, together
            // with any ghost cells lying outside the domain, that can be periodically shifted into the
            // valid region.  "compfvalid" is the complement of the "fvalid", while "compfvalid_grown" is 
            // "compfvalid" grown by one.  Using these we can figure out whether or not a cell is in the
            // valid region of our MultiFab as well as whether or not we're at a Fine->Crse boundary.
            //   

            for (int lev = lev_min; lev <= finest_level; lev++)
            {
                const amrex::Geometry& gm        = m_gdb->Geom(lev);
                const amrex::Geometry& gm_fine   = (lev < finest_level) ? m_gdb->Geom(lev+1) : gm;
                const amrex::Geometry& gm_coarse = (lev > 0) ? m_gdb->Geom(lev-1) : gm;
                const amrex::Box&      dm        = gm.Domain();
                const amrex::Real*     dx        = gm.CellSize();
                const amrex::Real*     plo       = gm.ProbLo();
                const amrex::Real*     dx_fine   = (lev < finest_level) ? m_gdb->Geom(lev+1).CellSize() : dx;
                const amrex::Real*     dx_coarse = (lev > 0) ? m_gdb->Geom(lev-1).CellSize() : dx;
                const int       lev_index = lev - lev_min;
                const amrex::BoxArray& grids     = mf[lev_index]->boxArray();
                const int       dgrow     = (lev == 0) ? 1 : m_gdb->MaxRefRatio(lev-1);

                amrex::BoxArray compfvalid, compfvalid_grown, fvalid = mf[lev_index]->boxArray();

                //
                // Do we have Fine->Crse overlap on a periodic boundary?
                // We want to add all ghost cells that can be shifted into valid region.
                //
                amrex::BoxList valid;

                for (int i = 0; i < grids.size(); i++)
                {
                    if (gm.isAnyPeriodic())
                    {
                        const amrex::Box& dest = amrex::grow(grids[i],dgrow);

                        if ( ! dm.contains(dest))
                        {
                            for (int j = 0; j < grids.size(); j++)
                            {
                                BL_ASSERT(dm.contains(grids[j]));
            
                                gm.periodicShift(dest, grids[j], pshifts);

                                for (const auto& iv : pshifts)
                                {
                                    const amrex::Box& sbx = grids[j] + iv;
                                    const amrex::Box& dbx = dest & sbx;

                                    BL_ASSERT(dbx.ok());

                                    valid.push_back(dbx);
                                }
                            }
                        }
                    }
                }
                if (valid.isNotEmpty())
                {
                    //
                    // We've got some Fine->Crse periodic overlap.
                    // Don't forget to add the valid boxes too.
                    //
                    for (int i = 0; i < grids.size(); i++) {
                        valid.push_back(grids[i]);
                    }
                    fvalid = amrex::BoxArray(valid);
                    fvalid.removeOverlap();
                }

                //
                // If we're at a lev < finestLevel, this is the coarsened fine BoxArray.
                // We use this for figuring out Crse->Fine issues.
                //
                amrex::BoxArray ccba;
                if (lev > 0)
                {
                    ccba = m_gdb->boxArray(lev);
                    ccba.coarsen(m_gdb->refRatio(lev-1));
                }
                amrex::BoxArray cfba;
                if (lev < finest_level)
                {
                    cfba = m_gdb->boxArray(lev+1);
                    cfba.coarsen(m_gdb->refRatio(lev));

                    BL_ASSERT(mf[lev_index]->boxArray().contains(cfba));
                }

                //
                // This is cfba with any shifted ghost cells.
                //
                amrex::BoxArray cfvalid = cfba;

                if (lev < finest_level)
                {
                    amrex::BoxList cvalid;

                    const amrex::BoxArray& cgrids = mf[lev_index]->boxArray();

                    for (int i = 0; i < cfba.size(); i++)
                    {
                        if (gm.isAnyPeriodic())
                        {
                            const amrex::Box& dest = amrex::grow(cfba[i], mf[lev_index]->nGrow());

                            if ( ! dm.contains(dest))
                            {
                                for (int j = 0; j < cgrids.size(); j++)
                                {
                                    BL_ASSERT(dm.contains(cgrids[j]));

                                    gm.periodicShift(dest, cgrids[j], pshifts);

                                    for (const auto& kiv : pshifts)
                                    {
                                        const amrex::Box& sbx = cfba[i] - kiv;

                                        cvalid.push_back(sbx);
                                    }
                                }
                            }
                        }
                    }
                    if (cvalid.isNotEmpty())
                    {
                        //
                        // We've got some Fine->Crse periodic overlap.
                        // Don't forget to add the valid boxes too.
                        //
                        for (int i = 0; i < cfba.size(); i++) {
                            cvalid.push_back(cfba[i]);
                        }
                        cfvalid = amrex::BoxArray(cvalid);
                        cfvalid.removeOverlap();
                    }
                }

                //
                // The "+1" is so we enclose the valid region together with any
                //  ghost cells that can be periodically shifted into valid.
                //
                compfvalid = amrex::complementIn(amrex::grow(dm,dgrow+1), fvalid);

                compfvalid_grown = compfvalid;
                compfvalid_grown.grow(1);
                compfvalid_grown.removeOverlap();
            
                if (gm.isAnyPeriodic() && ! gm.isAllPeriodic())
                {
                    amrex::Error("AssignDensity: problem must be periodic in no or all directions");
                }
                //
                // If we're at a lev > 0, this is the coarsened BoxArray.
                // We use this for figuring out Fine->Crse issues.
                //
                amrex::BoxArray cba;
                if (lev > 0)
                {
                    cba = m_gdb->boxArray(lev);
                    cba.coarsen(m_gdb->refRatio(lev-1));
                }
                //
                // Do the grids at this level cover the full domain? If they do
                // there can be no Fine->Crse interactions at this level.
                //
                const bool GridsCoverDomain = fvalid.contains(m_gdb->Geom(lev).Domain());

                for (size_t ip = start_idx; ip < LocalNum; ++ip) {
                    //there are no more particles in level lev on this node
                    //exit the loop and move to the next level
                    if (m_lev[ip] != (unsigned)lev) {
                        start_idx = ip;
                        break;
                    }

                    amrex::FArrayBox&  fab = (*mf[lev_index])[m_grid[ip]];
    
                    //
                    // Get "fracs" and "cells" for the particle "p" at this level.
                    //
                    const int M = CIC_Cells_Fracs(this->R[ip], plo, dx, dx, fracs, cells);

                    //
                    // If this is not fully periodic then we have to be careful that no
                    // particle's support leaves the domain. We test this by checking the low
                    // and high corners respectively.
                    //
                    if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
                        if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
                            amrex::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
                        }
                    }

                    //
                    // This section differs based on whether we subcycle.
                    // Without subcycling we use the "stretchy" support for particles.
                    // With subcycling a particles support is strictly defined 
                    // by its resident level.
                    //
                    if (sub_cycle)
                    {
                        bool isFiner = false;
                        bool isBoundary = false;
                        //
                        // First sum the mass in the valid region
                        //
                        for (int i = 0; i < M; i++)
                        {
                            if (cfvalid.contains(cells[i]))
                            {
                                //
                                // Some part of the particle's mass lies in a 
                                // finer region; we'll deal with it shortly.
                                //
                                isFiner    = true;
                                isBoundary = true;
                                continue;
                            }
                            if ( ! fvalid.contains(cells[i]))
                            {
                                //
                                // We're out of the valid region.
                                //
                                isBoundary = true;
                                continue;
                            }
                            //
                            // Sum up mass in first component.
                            //
                            fab(cells[i],0) += pa[ip] * fracs[i];

                            // If the domain is not periodic and we want to let particles
                            //    live near the boundary but "throw away" the contribution that 
                            //    does not fall into the domain ...
                            if ( ! gm.isAllPeriodic() && allow_particles_near_boundary &&
                                ! gm.Domain().contains(cells[i]))
                            {
                                continue;
                            }
                        }
                        //
                        // Deal with mass that doesn't belong at this level.
                        // Here we assume proper nesting so that only one special case can
                        // be true for a given particle.
                        //
                        if (isBoundary)
                        {
                            if (isFiner)
                            {
                                BL_ASSERT(lev < finest_level);
                                //
                                // We're at a coarse->fine interface
                                //
                                // get fine cells/fracs
                                //
                                const int MF = CIC_Cells_Fracs(this->R[ip], plo, dx_fine ,dx, ffracs, fcells);

                                for (int j = 0; j < MF; j++)
                                {
                                    //
                                    // Make sure this fine cell is valid. Check for periodicity.
                                    //
                                    const amrex::Box bx(fcells[j],fcells[j]);
                                    gm_fine.periodicShift(bx, gm_fine.Domain(), pshifts);
                                    if ( !pshifts.empty() )
                                    {
                                        BL_ASSERT(pshifts.size() == 1);
                                        fcells[j] = fcells[j] - pshifts[0];
                                    }
                                    mf[lev_index+1]->boxArray().intersections(amrex::Box(fcells[j], fcells[j]),
                                                                             isects,true,0);
                                    if (isects.size() == 0) {
                                        continue;
                                    }
                                    const int grid = isects[0].first; 
                                    const int who  = mf[lev_index+1]->DistributionMap()[m_grid[ip]];
                                    if (who == amrex::ParallelDescriptor::MyProc())
                                    {
                                        //
                                        // Sum up mass in first component.
                                        //
                                        (*mf[lev_index+1])[m_grid[ip]](fcells[j],0) += pa[ip] * ffracs[j];
                                    }
                                    else
                                    {
                                        pb.m_lev  = lev+1;
                                        pb.m_grid = m_grid[ip];
                                        pb.m_cell = fcells[j];
                                        pb.m_data[0] = pa[ip] *  ffracs[j];
                                        data[who].push_back(pb);
                                    }
                                }
                            }
                            else if (lev_index > 0)
                            {
                                const int MC = CIC_Cells_Fracs(this->R[ip], plo, dx_coarse, dx, cfracs, ccells);
                                for (int j = 0; j < MC; j++)
                                {
                                    //
                                    // Make sure this coarse cell isn't in this level's valid region.
                                    // This may not matter.
                                    //
                                    if (cba.contains(ccells[j]))
                                        continue;
                                    //
                                    // Check for periodicity.
                                    //
                                    const amrex::Box bx(ccells[j],ccells[j]);
                                    gm_coarse.periodicShift(bx, gm_coarse.Domain(), pshifts);
            
                                    if ( ! pshifts.empty())
                                    {
                                        BL_ASSERT(pshifts.size() == 1);
                                        ccells[j] = ccells[j] - pshifts[0]; 
                                    }
                                    //
                                    // Find its resident grid.
                                    //
                                    mf[lev_index - 1]->boxArray().intersections(amrex::Box(ccells[j],ccells[j]),
                                                                                isects,true,0);
                                    if (isects.size() == 0) {
                                        continue;
                                    }
                                    const int grid = isects[0].first;
                                    const int who  = mf[lev_index-1]->DistributionMap()[grid];
                                    if (who == amrex::ParallelDescriptor::MyProc())
                                    {
                                        //
                                        // Sum up mass in first component.
                                        //
                                        (*mf[lev_index-1])[m_grid[ip]](ccells[j],0) += pa[ip] * cfracs[j];
                                    }
                                    else
                                    {
                                        pb.m_lev  = lev-1;
                                        pb.m_grid = m_grid[ip];
                                        pb.m_cell = ccells[j];
                                        //
                                        // Sum up mass in first component.
                                        //
                                        pb.m_data[0] = pa[ip] * cfracs[j];
                                        data[who].push_back(pb);
                                    }
                                }
                            }
                            else
                            {
                                // The mass is below levels we care about. Ignore it.
                            }
                        }
                    }
                    else
                    {
                        bool AnyCrseToFine = false;
                        if (lev < finest_level) {
                            // dummy template values
                            amrex::ParticleContainer<0> p;
                            AnyCrseToFine = p.CrseToFine(cfba,cells,cfshifts,gm,cwhich,pshifts);
                        }
                        //
                        // lev_index > 0 means that we don't do F->C for lower levels
                        // This may mean that the mass fraction is off.
                        //
                        bool AnyFineToCrse = false;
                        if ( lev_index > 0 && !GridsCoverDomain )
                            AnyFineToCrse = FineToCrse(ip,lev,cells,fvalid,compfvalid_grown,
                                                        ccells,cfracs,fwhich,cgrid,pshifts,isects);
	  
                        BL_ASSERT(!(AnyCrseToFine && AnyFineToCrse));
                        if ( ! AnyCrseToFine && ! AnyFineToCrse)
                        {
                            //
                            // By far the most common case.  Just do it!
                            //
                            for (int i = 0; i < M; i++)
                            {
	      
                                // If the domain is not periodic and we want to let particles
                                //    live near the boundary but "throw away" the contribution that 
                                //    does not fall into the domain ...
                                if (! gm.isAllPeriodic() && allow_particles_near_boundary 
                                    && ! gm.Domain().contains(cells[i]))
                                {
                                    continue;
                                }
                                //
                                // Sum up mass in first component.
                                //
                                fab(cells[i],0) += pa[ip] * fracs[i];
                            }
                        }
                        else if (AnyFineToCrse)
                        {
                            amrex::Real sum_crse = 0, sum_fine = 0;

                            for (int i = 0; i < M; i++)
                            {
                                if (fwhich[i])
                                {
                                    //
                                    // We're at a Fine->Crse boundary.
                                    //
                                    BL_ASSERT(cgrid[i] >= 0);
                                    BL_ASSERT(cgrid[i] < mf[lev_index-1]->size());
                                    //
                                    // Here we need to update the crse region.  The coarse
                                    // region is always going to be updated if we have a
                                    // particle in a cell bordering a Fine->Crse boundary.
                                    //
                                    const int who = mf[lev_index-1]->DistributionMap()[cgrid[i]];
                                    if (who == amrex::ParallelDescriptor::MyProc())
                                    {
                                        if ( ! (*mf[lev_index-1])[cgrid[i]].box().contains(ccells[i])) {
                                            continue;
                                        }

                                        // If the domain is not periodic and we want to let particles
                                        //    live near the boundary but "throw away" the contribution that 
                                        //    does not fall into the domain ...
                                        if (! gm_coarse.isAllPeriodic() && allow_particles_near_boundary &&
                                            ! gm_coarse.Domain().contains(ccells[i]))
                                        {
                                            continue;
                                        }

                                        //
                                        // Sum up mass in first component.
                                        //
                                        {
                                            (*mf[lev_index-1])[cgrid[i]](ccells[i],0) += pa[ip] * cfracs[i];
                                        }
                                    }
                                    else
                                    {
                                        pb.m_lev  = lev-1;
                                        pb.m_grid = cgrid[i];
                                        pb.m_cell = ccells[i];
                                        //
                                        // Sum up mass in first component.
                                        //
                                        pb.m_data[0] = pa[ip] * cfracs[i];
                                        data[who].push_back(pb);
                                    }
                                    sum_crse += cfracs[i];
                                }
                            }

                            //
                            // We've updated the Crse cells.  Now we have to update the fine
                            // cells in such a way that the total amount of mass we move
                            // around is precisely p.m_data[0]. In other words, the fractions
                            // we use at crse and fine have to sum to zero.  In the fine
                            // case, we have to account for the case where one or more of the
                            // cell indices is not in the valid region of the box containing 
                            // the particle.
                            //
                            sum_fine = 0;
                            for (int i = 0; i < M; i++) 
                            {
                                //
                                // Reusing "fwhich" to indicate fine cells that need massaging.
                                //
                                fwhich[i] = true;

                                if ( ! compfvalid_grown.contains(cells[i]))
                                {
                                    //
                                    // Go ahead and add the full correct amount to these cells.
                                    // They can't touch a Fine->Crse boundary.
                                    //
                                    sum_fine += fracs[i];
                                    //
                                    // Sum up mass in first component.
                                    //
            
                                    fab(cells[i],0) += pa[ip] * fracs[i];
                                    fwhich[i] = false;
                                }
                                else if (compfvalid.contains(cells[i]))
                                {
                                    fwhich[i] = false;
                                }
                            }

                            const amrex::Real sum_so_far = sum_crse + sum_fine; 
        
                            BL_ASSERT(sum_so_far > 0);
                            BL_ASSERT(sum_so_far < 1);

                            sum_fine = 0;
                            for (int i = 0; i < M; i++) 
                            {       
                                if (fwhich[i])
                                    //
                                    // Got to weight cells in this direction differently.
                                    //
                                    sum_fine += fracs[i];
                            }
        
                            const amrex::Real mult = (1 - sum_so_far) / sum_fine;
                            sum_fine = 0;
                            for (int i = 0; i < M; i++)
                            {
                                if (fwhich[i])
                                {
                                    //
                                    // Sum up mass in first component.
                                    //
                                    fab(cells[i],0) += pa[ip] * fracs[i] * mult;
                
                                    sum_fine += fracs[i] * mult;
                                }
                            }

                            BL_ASSERT(std::abs(1-(sum_fine+sum_so_far)) < 1.e-9);

                        }
                        else if (AnyCrseToFine)
                        {
                            amrex::Real sum = 0;
                            for (int i = 0; i < M; i++)
                            {
                                if (!cwhich[i])
                                {
                                    // If the domain is not periodic and we want to let particles
                                    //    live near the boundary but "throw away" the contribution that 
                                    //    does not fall into the domain ...
                                    if ( ! gm.isAllPeriodic() && allow_particles_near_boundary &&
                                        ! gm.Domain().contains(ccells[i]))
                                    {
                                        continue;
                                    }
                                    //
                                    // Sum up mass in first component.
                                    //
                                    fab(cells[i],0) += pa[ip] * fracs[i];

                                    sum += fracs[i];
                                }
                                else 
                                {
                                    //
                                    // We're at a Crse->Fine boundary.
                                    //
                                    FineCellsToUpdateFromCrse(ip,lev,cells[i],cfshifts[i],
                                                            fgrid,ffracs,fcells,isects);
            
                                    for (int j = 0, nj=fcells.size(); j < nj; ++j)
                                    {
                                        const int who = mf[lev_index+1]->DistributionMap()[fgrid[j]];
                                        if (who == amrex::ParallelDescriptor::MyProc())
                                        {
                                            //
                                            // Sum up mass in first component.
                                            //
                                            (*mf[lev_index+1])[fgrid[j]](fcells[j],0) += pa[ip] * fracs[i] * ffracs[j];
                                        }
                                        else
                                        {
                                            pb.m_lev  = lev+1;
                                            pb.m_grid = fgrid[j];
                                            pb.m_cell = fcells[j];
                                            //
                                            // Sum up mass in first component.
                                            //
                                            pb.m_data[0] = pa[ip] * fracs[i] * ffracs[j];
                
                                            data[who].push_back(pb);
                                        }
                
                                        sum += fracs[i] * ffracs[j];
                                    }
                                }
                            }
                            BL_ASSERT(std::abs(1-sum) < 1.e-9);
                        }
                    }
                }
            }
            
            AssignDensityDoit(0, mf, data, 1, lev_min);
            for (int lev = lev_min; lev <= finest_level; lev++)
            {
                const int       lev_index = lev - lev_min;
                const amrex::Geometry& gm        = m_gdb->Geom(lev);
                const amrex::Real*     dx        = gm.CellSize();
                const amrex::Real      vol       = D_TERM(dx[0], *dx[1], *dx[2]);
    
                mf[lev_index]->SumBoundary(gm.periodicity());

                //
                // Only multiply the first component by (1/vol) because this converts mass
                // to density. If there are additional components (like velocity), we don't
                // want to divide those by volume.
                //
                mf[lev_index]->mult(1/vol,0,1);
            }

            //
            // The size of the returned multifab is limited by lev_min and 
            // finest_level. In the following code, lev is the real level,  
            // lev_index is the corresponding index for mf. 
            //
            // I believe that we don't need any information in ghost cells so we don't copy those.
            //
            if ( ! all_grids_the_same)
            {
                for (int lev = lev_min; lev <= finest_level; lev++)
                {
                    const int lev_index = lev - lev_min;
                    mf_to_be_filled[lev_index]->copy(*mf_part[lev_index],0,0,1);
                }
            }

            IpplTimings::stopTimer(AssignDensityTimer_m);
        }

    //Function from AMReX adjusted to work with Ippl AmrParticleBase class
    //Assign density for a single level
    template <class AType> 
        void AssignDensitySingleLevel (ParticleAttrib<AType> &pa, 
                                       amrex::MultiFab& mf_to_be_filled,
                                       int lev,
                                       int particle_lvl_offset = 0)
        {
            if (mf_to_be_filled.is_nodal()) {
                NodalDepositionSingleLevel(pa, mf_to_be_filled, lev, particle_lvl_offset);
            } else if (mf_to_be_filled.boxArray().ixType().cellCentered()) {
                AssignCellDensitySingleLevel(pa, mf_to_be_filled, lev, particle_lvl_offset);
            } else {
                amrex::Abort("AssignCellDensitySingleLevel: mixed type not supported");
            }
        }

    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    template <class AType>
        void AssignCellDensitySingleLevel(ParticleAttrib<AType> &pa,
                                          amrex::MultiFab& mf_to_be_filled,
                                          int lev,
                                          int particle_lvl_offset = 0)
        {
    
            amrex::MultiFab* mf_pointer;
            PLayout *Layout = &this->getLayout();
            const amrex::ParGDBBase* m_gdb = Layout->GetParGDB();

            if ( m_gdb->OnSameGrids(lev, mf_to_be_filled) ) {
                // If we are already working with the internal mf defined on the 
                // particle_box_array, then we just work with this.
                mf_pointer = &mf_to_be_filled;
            } else {
                // If mf_to_be_filled is not defined on the particle_box_array, then we need 
                // to make a temporary here and copy into mf_to_be_filled at the end.
                mf_pointer = new amrex::MultiFab(m_gdb->ParticleBoxArray(lev), m_gdb->ParticleDistributionMap(lev), 1,
                                          mf_to_be_filled.nGrow());
            }
  
            // We must have ghost cells for each FAB so that a particle in one grid can spread its effect to an
            //    adjacent grid by first putting the value into ghost cells of its own grid.  The mf->sumBoundary call then
            //    adds the value from one grid's ghost cell to another grid's valid region.
            if (mf_pointer->nGrow() < 1) 
                amrex::Error("Must have at least one ghost cell when in AssignDensitySingleLevel");

            const amrex::Geometry& gm          = m_gdb->Geom(lev);
            const amrex::Real*     plo         = gm.ProbLo();
            const amrex::Real*     dx_particle = m_gdb->Geom(lev + particle_lvl_offset).CellSize();
            const amrex::Real*     dx          = gm.CellSize();

            if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) {
                amrex::Error("AssignDensity: problem must be periodic in no or all directions");
            }

            for (amrex::MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
                (*mf_pointer)[mfi].setVal(0);
            }
    
            //loop trough particles and distribute values on the grid
            size_t LocalNum = this->getLocalNum();
            for (size_t ip = 0; ip < LocalNum; ++ip) 
                {
      
                    //if particle doesn't belong on this level exit loop
                    if (m_lev[ip] != (unsigned)lev)
                        break;

                    amrex::FArrayBox& fab = (*mf_pointer)[m_grid[ip]];

                    amrex::Array<amrex::Real> fracs;
                    amrex::Array<amrex::IntVect> cells;

                    const int M = CIC_Cells_Fracs(this->R[ip], plo, dx, dx_particle, fracs, cells);
                    //
                    // If this is not fully periodic then we have to be careful that the
                    // particle's support leaves the domain unless we specifically want to ignore
                    // any contribution outside the boundary (i.e. if allow_particles_near_boundary = true). 
                    // We test this by checking the low and high corners respectively.
                    //
                    if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
                        if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
                            amrex::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
                        }
                    }

                    for (int i = 0; i < M; ++i)
                        {
                            if ( !fab.box().contains(cells[i])) {
                                continue;
                            }

                            // If the domain is not periodic and we want to let particles
                            //    live near the boundary but "throw away" the contribution that 
                            //    does not fall into the domain ...
                            if ( ! gm.isAllPeriodic() && allow_particles_near_boundary && ! gm.Domain().contains(cells[i])) {
                                continue;
                            }

                            //sum up mass in the particle attribute
                            fab(cells[i], 0) += pa[ip] * fracs[i];

                        }
      
                }

            mf_pointer->SumBoundary(gm.periodicity());

            //
            // Only multiply the first component by (1/vol) because this converts mass
            // to density. If there are additional components (like velocity), we don't
            // want to divide those by volume.
            //
            const amrex::Real vol = D_TERM(dx[0], *dx[1], *dx[2]);

            mf_pointer->mult(1/vol,0,1);

            // If mf_to_be_filled is not defined on the particle_box_array, then we need
            // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
            // need any information in ghost cells so we don't copy those.
            if (mf_pointer != &mf_to_be_filled)
                {
                    mf_to_be_filled.copy(*mf_pointer,0,0,0);
                    delete mf_pointer;
                }
    
        }

    //Function from AMReX adjusted to work with Ippl AmrParticleBase class
    template<class AType>
        void NodalDepositionSingleLevel(ParticleAttrib<AType> &pa,
                                        amrex::MultiFab& mf_to_be_filled,
                                        int lev,
                                        int particle_lvl_offset = 0)
        {
            amrex::MultiFab* mf_pointer;
            PLayout *Layout = &this->getLayout();
            const amrex::ParGDBBase* m_gdb = Layout->GetParGDB();

            if ( m_gdb->OnSameGrids(lev, mf_to_be_filled) )
                {
                    // If we are already working with the internal mf defined on the 
                    // particle_box_array, then we just work with this.
                    mf_pointer = &mf_to_be_filled;
                }
            else
                {
                    // If mf_to_be_filled is not defined on the particle_box_array, then we need 
                    // to make a temporary here and copy into mf_to_be_filled at the end.
                    mf_pointer = new amrex::MultiFab(amrex::convert(m_gdb->ParticleBoxArray(lev),
                                                             mf_to_be_filled.boxArray().ixType()),
                                              m_gdb->ParticleDistributionMap(lev),
                                              1, mf_to_be_filled.nGrow());
                }

            const amrex::Geometry& gm          = m_gdb->Geom(lev);
            const amrex::Real*     dx          = gm.CellSize();
    
            if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) 
                amrex::Error("AssignDensity: problem must be periodic in no or all directions");

            mf_pointer->setVal(0.0);

            amrex::Array<amrex::IntVect> cells;
            cells.resize(8);
    
            amrex::Array<amrex::Real> fracs;
            fracs.resize(8);

            amrex::Array<amrex::Real> sx;
            sx.resize(2);
            amrex::Array<amrex::Real> sy;
            sy.resize(2);
            amrex::Array<amrex::Real> sz;
            sz.resize(2);

            //loop trough particles and distribute values on the grid
            size_t LocalNum = this->getLocalNum();
            for (size_t ip = 0; ip < LocalNum; ++ip) 
                {
                    amrex::FArrayBox& fab = (*mf_pointer)[m_grid[ip]];
                    
//                     // FIXME Is "Where" really needed?
//                     Where(this->R, ip);
                    
                    amrex::IntVect m_cell = Layout->Index(this->R[ip], this->m_lev[ip]);
                    cells[0] = m_cell;
                    cells[1] = m_cell+amrex::IntVect(1,0,0);
                    cells[2] = m_cell+amrex::IntVect(0,1,0);
                    cells[3] = m_cell+amrex::IntVect(1,1,0);
                    cells[4] = m_cell+amrex::IntVect(0,0,1);
                    cells[5] = m_cell+amrex::IntVect(1,0,1);
                    cells[6] = m_cell+amrex::IntVect(0,1,1);
                    cells[7] = m_cell+amrex::IntVect(1,1,1);

                    amrex::Real x = this->R[ip][0] / dx[0];
                    amrex::Real y = this->R[ip][1] / dx[1];
                    amrex::Real z = this->R[ip][2] / dx[2];

                    int i = m_cell[0];
                    int j = m_cell[1];
                    int k = m_cell[2];

                    amrex::Real xint = x - i;
                    amrex::Real yint = y - j;
                    amrex::Real zint = z - k;

                    sx[0] = 1.0-xint;
                    sx[1] = xint;
                    sy[0] = 1.0-yint;
                    sy[1] = yint;
                    sz[0] = 1.0-zint;
                    sz[1] = zint;

                    fracs[0] = sx[0] * sy[0] * sz[0];
                    fracs[1] = sx[1] * sy[0] * sz[0];
                    fracs[2] = sx[0] * sy[1] * sz[0];
                    fracs[3] = sx[1] * sy[1] * sz[0];
                    fracs[4] = sx[0] * sy[0] * sz[1];
                    fracs[5] = sx[1] * sy[0] * sz[1];
                    fracs[6] = sx[0] * sy[1] * sz[1];
                    fracs[7] = sx[1] * sy[1] * sz[1];

                    for (int i = 0; i < 8; i++)
                        fab(cells[i],0) += pa[ip] * fracs[i];

                }

            mf_pointer->SumBoundary(gm.periodicity());
            //
            // Only multiply the first component by (1/vol) because this converts mass
            // to density. If there are additional components (like velocity), we don't
            // want to divide those by volume.
            //
            const amrex::Real vol = D_TERM(dx[0], *dx[1], *dx[2]);

            mf_pointer->mult(1/vol,0,1);

            // If mf_to_be_filled is not defined on the particle_box_array, then we need
            // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
            // need any information in ghost cells so we don't copy those.
            if (mf_pointer != &mf_to_be_filled)
                {
                    mf_to_be_filled.copy(*mf_pointer,0,0,0);
                    delete mf_pointer;
                }

        }

    //gather values from grid back to the particles
    //loop trough particles and use BoxLin Interp functions to get the particles value
    //
    template <class AType>
        void GetGravity(ParticleAttrib<AType> &pa,
                        amrex::Array<std::unique_ptr<amrex::MultiFab> > &mf) 
        {

            PLayout *Layout = &this->getLayout();
            const amrex::ParGDBBase* m_gdb = Layout->GetParGDB();
            size_t LocalNum = this->getLocalNum();

            //loop trough all the particles
            for (size_t ip = 0; ip < LocalNum; ++ip) {
                int lev = m_lev[ip];
                int grid = m_grid[ip];

                //get the FArrayBox where this particle is located
                amrex::FArrayBox& fab = (*(mf[lev].get()))[grid];
      
                amrex::Real grav[AMREX_SPACEDIM];
                int idx[AMREX_SPACEDIM] = {  D_DECL(0,1,2) };
      
                //get the value at grid point in grav array
                Interp(this->R[ip], m_gdb->Geom(lev), fab, idx, grav, AMREX_SPACEDIM);

                //assign to particle attribute
                for (int i = 0; i < AMREX_SPACEDIM; ++i)
                    pa[ip][i] = grav[i];
            }

        }

 
};

#include "AmrParticleBase.hpp"

#endif
