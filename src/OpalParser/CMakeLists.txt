set (_SRCS
  CompoundStatement.cpp
  IfStatement.cpp
  Macro.cpp
  MacroCmd.cpp
  MacroStream.cpp
  OpalParser.cpp
  WhileStatement.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    CompoundStatement.h
    IfStatement.h
    MacroCmd.h
    Macro.h
    MacroStream.h
    OpalParser.h
    WhileStatement.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpalParser")