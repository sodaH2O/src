CMAKE_MINIMUM_REQUIRED (VERSION 2.8.10)
PROJECT (SDDS2OPAL)
SET (SDDS2OPAL_VERSION_MAJOR 0)
SET (SDDS2OPAL_VERSION_MINOR 1)

configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

INCLUDE_DIRECTORIES (
    ${SDDS_INCLUDE_DIR}
)

LINK_DIRECTORIES (
    ${SDDS_LIBRARY_DIR}
)

SET (SDDS_LIBS
    mdbcommon
    SDDS1
    mdblib
    lzma
    z
)

ADD_EXECUTABLE( sdds2opal main.cpp )
TARGET_LINK_LIBRARIES( sdds2opal ${SDDS_LIBS} )

INSTALL(TARGETS sdds2opal RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")