set (CLASSIC_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

add_subdirectory(AbsBeamline)
add_subdirectory(Algebra)
add_subdirectory(Algorithms)
add_subdirectory(BeamlineCore)
add_subdirectory(BeamlineGeometry)
add_subdirectory(Beamlines)
add_subdirectory(Channels)
add_subdirectory(ComponentWrappers)
add_subdirectory(Construction)
add_subdirectory(Fields)
add_subdirectory(Filters)
add_subdirectory(FixedAlgebra)
add_subdirectory(MemoryManagement)
add_subdirectory(Parser)
add_subdirectory(Physics)
add_subdirectory(Solvers)
add_subdirectory(Structure)
add_subdirectory(Utilities)
add_subdirectory(TrimCoils)

install (FILES Utilities/ClassicException.h Utilities/ClassicRandom.h Utilities/Options.h Utilities/OptionTypes.h 
    DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Utilities"
)

set(OPAL_SRCS "${OPAL_SRCS}" PARENT_SCOPE)
