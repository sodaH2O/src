# this prevents compile errors (due to warnings) since Boost Spirits grammar
# definitions are interpreted as parentheses errors (overloaded operators).
set_source_files_properties(${CMAKE_SOURCE_DIR}/Expression/Parser/expression.cpp PROPERTIES COMPILE_FLAGS "-Wno-error=parentheses" )

set (PILOT_FEVAL_SRCS
    ${CMAKE_SOURCE_DIR}/optimizer/Simulation/FunctionEvaluator.cpp
)

IF (ENABLE_DKS)
   include_directories (${DKS_INCLUDE_DIR})
   link_directories (${DKS_LIBRARY_DIR})
ENDIF (ENABLE_DKS)

add_subdirectory(pisa-standalone)