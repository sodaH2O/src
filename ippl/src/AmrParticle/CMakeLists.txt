SET (_SRCS
    )

SET (_HDRS
    ParticleAmrLayout.h
    AmrParticleBase.h
    AmrParticleBase.hpp
    AmrParticleLevelCounter.h
    )

INCLUDE_DIRECTORIES (
    ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_IPPL_SOURCES (${_SRCS})
ADD_IPPL_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/AmrParticle)
