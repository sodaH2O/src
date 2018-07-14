#ifndef AMR_TRILINEAR_INTERPOLATER_H
#define AMR_TRILINEAR_INTERPOLATER_H

#include "AmrInterpolater.h"

#include <algorithm>
#include <iterator>

template <class Level>
class AmrTrilinearInterpolater : public AmrInterpolater<Level>
{
public:
    typedef typename Level::go_t        go_t;
    typedef typename Level::lo_t        lo_t;
    typedef typename Level::scalar_t    scalar_t;
    typedef typename Level::umap_t      umap_t;
    typedef typename Level::basefab_t   basefab_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    
public:
    
    AmrTrilinearInterpolater();
    
    void stencil(const AmrIntVect_t& iv,
                 const basefab_t& fab,
                 umap_t& map,
                 const scalar_t& scale,
                 Level* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                umap_t& map,
                const scalar_t& scale,
                lo_t dir, lo_t shift, const basefab_t& rfab,
                const AmrIntVect_t& riv,
                Level* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              umap_t& map,
              const scalar_t& scale,
              lo_t dir, lo_t shift, const basefab_t& fab,
              Level* mglevel);
};


#include "AmrTrilinearInterpolater.hpp"

#endif
