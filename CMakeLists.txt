set(TARGET_NAME "simulate_dos")    # Name of the test source file to compile.

set(BINDIR  "${CMAKE_CURRENT_SOURCE_DIR}/bin" )    # Directory where to store the executables.
set(DOCDIR  "${CMAKE_CURRENT_SOURCE_DIR}/doc" )    # Directory where to store the API documentation.

##--------------- Do not edit below this line ----------------##

cmake_minimum_required(VERSION 2.8)
project(DOS_Extraction CXX)

################################################################
## Set up some variables.
################################################################
set(SRCDIR  "${CMAKE_CURRENT_SOURCE_DIR}/src" )    # Directory where to find the library sources.
set(TESTDIR "${CMAKE_CURRENT_SOURCE_DIR}/test")    # Directory where to find the test sources.

file(GLOB SRCS ${SRCDIR}/*.cc ${TESTDIR}/${TARGET_NAME}.cc)
file(GLOB HDRS ${SRCDIR}/*.h  ${TESTDIR}/${TARGET_NAME}.h )

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(EXECUTABLE_OUTPUT_PATH "${BINDIR}")

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall --pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING    # Default build type: Release.
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

################################################################
## Dependencies and third-party libraries.
################################################################
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/include")

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

find_package(OpenMP)

if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Gnuplot)

################################################################
## Documentation generation settings.
################################################################
set(DOXYFILE_SOURCE_DIR "${DOCDIR}")
set(DOXYFILE_OUTPUT_DIR "${DOCDIR}")
set(DOXYFILE_LATEX "ON")

include(RunDoxygen)

################################################################
## Artistic style.
################################################################
find_package(AStyle)

if(ASTYLE_FOUND)
    add_custom_target(astyle ALL
        COMMAND ${ASTYLE_EXECUTABLE} -q -A7 -s2 -C -S -N -Y -f -p -H -E ${SRCS} ${HDRS}
        COMMENT "Formatting source codes..."
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

################################################################
## Executables definition.
################################################################
add_executable(${TARGET_NAME} ${SRCS})

set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME "${TARGET_NAME}")    # Executable filename.

target_link_libraries(${TARGET_NAME} ${Boost_LIBRARIES})
