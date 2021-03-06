set (_SRCS
  ComplexEigen.cpp
  DoubleEigen.cpp
  DynamicFixedPoint.cpp
  NormalForm.cpp
  StaticFixedPoint.cpp
  TpsComplex.cpp
  Tps.cpp
  TpsData.cpp
  TpsDouble.cpp
  TpsMonomial.cpp
  TpsTps.cpp
  Vps.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    Array1D.h
    Array2D.h
    ComplexEigen.h
    DoubleEigen.h
    DynamicFixedPoint.h
    LieMap.h
    LUMatrix.h
    Matrix.h
    NormalForm.h
    SliceIterator.h
    StaticFixedPoint.h
    TpsData.h
    Tps.h
    TpsMath.h
    TpsMonomial.h
    TpsSubstitution.h
    Vector.h
    Vps.h
    VpsInvMap.h
    VpsMap.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algebra")
install (FILES ${_SRCS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algebra")
