find_package (GTest REQUIRED)

set (TEST_EXE opal_unit_tests)
set (SRC_DIRS classic_src opal_src)

add_compile_options (${IPPL_CMAKE_CXX_FLAGS} )
add_definitions (-DPARALLEL_IO -DWITH_UNIT_TESTS)

set (TEST_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
macro (add_sources)
    file (RELATIVE_PATH _relPath "${TEST_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND TEST_SRCS_LOCAL "${_relPath}/${_src}")
        else()
            list (APPEND TEST_SRCS_LOCAL "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (TEST_SRCS_LOCAL ${TEST_SRCS_LOCAL} PARENT_SCOPE)
    endif()
endmacro()

add_subdirectory (opal_src)
add_subdirectory (classic_src)
add_subdirectory (ippl_src)
add_subdirectory (opal_test_utilities)

set (TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Main.cpp ${TEST_SRCS_LOCAL})

link_directories (${GTEST_LIB_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${IPPL_LIBRARY_DIR}
    ${OPAL_SOURCE_DIR}/src
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
    ${AMREX_LIBRARY_DIR}
)

include_directories (${GTEST_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${GSL_INCLUDE_DIR}
    ${OPAL_SOURCE_DIR}/src
    ${OPAL_SOURCE_DIR}/src/Classic
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${IPPL_SOURCE_DIRS}
    ${AMREX_INCLUDE_DIRS}
    ${GTEST_INCLUDE_DIRS}
)

# Check to see if cmake finds the test files
# MESSAGE(STATUS "unit test src files: ${TEST_SRCS}")

# Build the test exe. We don't do an install on the unit test exe as it is
# assumed that this is internal to opal
# the opal and classic sources are not needed again if we link against libOPAL and libCLASSIC!
add_executable (${TEST_EXE} ${TEST_SRCS})
target_link_libraries (
	${TEST_EXE}
	libOPAL
	${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
	${AMREX_LIBRARIES}
	${OTHER_CMAKE_EXE_LINKER_FLAGS}
	${GTEST_BOTH_LIBRARIES}
	-lpthread
)
