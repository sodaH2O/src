set (_SRCS
  SavitzkyGolay.cpp
  Stencil.cpp
  FixedFFTLowPass.cpp
  RelativeFFTLowPass.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    Filter.h
    Filters.h
    FixedFFTLowPass.h
    RelativeFFTLowPass.h
    SavitzkyGolay.h
    Stencil.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Filters")