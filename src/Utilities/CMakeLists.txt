set (_SRCS
  OpalException.cpp
  OpalFilter.cpp
  OpalSection.cpp
  RegularExpression.cpp
  Round.cpp
  Timer.cpp
  TpsWerrf.cpp
  Truncate.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    NumToStr.h
    OpalException.h
    OpalFilter.h
    OpalSection.h
    RegularExpression.h
    Round.h
    Timer.h
    TpsWerrf.h
    Truncate.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Utilities")