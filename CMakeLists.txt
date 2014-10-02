set (EXEC "simulate_dos")    # Name of the test source file to compile.

set (BINDIR  "bin" )    # Directory where to store the executables.
set (DOCDIR  "doc" )    # Directory where to store the API documentation.

##-- Do not edit below this line --##

cmake_minimum_required (VERSION 2.8)
project (DOS_Extraction CXX)

############################################################
## Set up some variables.
############################################################
set (SRCDIR  "src" )    # Directory where to find the library sources.
set (TESTDIR "test")    # Directory where to find the test sources.

file (GLOB SRCS ./${SRCDIR}/*.cc ./${TESTDIR}/${EXEC}.cc)
file (GLOB HDRS ./${SRCDIR}/*.h  ./${TESTDIR}/${EXEC}.h )

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set (EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${BINDIR}")
set (CMAKE_CXX_FLAGS "-std=c++11 -Wall --pedantic")

option (DEBUG_MODE "Additional debug informations" OFF)
if (DEBUG_MODE)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
else ()
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif ()

############################################################
## Dependencies and third-party libraries.
############################################################
include_directories (SYSTEM "./" "./include")

find_package (Eigen3 REQUIRED)
include_directories (SYSTEM ${EIGEN3_INCLUDE_DIR})

find_package (Boost REQUIRED COMPONENTS iostreams system filesystem)
include_directories (SYSTEM ${Boost_INCLUDE_DIR})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package (Gnuplot)

############################################################
## Documentation generation settings.
############################################################
set (DOXYFILE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${DOCDIR}")
set (DOXYFILE_OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${DOCDIR}")
set (DOXYFILE_LATEX "ON")

include (RunDoxygen)

############################################################
## Artistic style.
############################################################
add_custom_target (astyle ALL
    COMMAND astyle -q -A10 -t4 -C -S -N -f -p -H -E ${SRCS} ${HDRS}
    COMMENT "Formatting source codes..."
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

############################################################
## Executables definition.
############################################################
add_executable (${EXEC} ${SRCS})
target_link_libraries (${EXEC} ${Boost_LIBRARIES})
