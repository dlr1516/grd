cmake_minimum_required(VERSION 2.8)
project(grd)
add_definitions(-std=c++0x)
set(CMAKE_BUILD_TYPE Release)

# Set the output directory for the build executables and libraries
SET(grd_RUNTIME_OUTPUT_DIRECTORY ${grd_SOURCE_DIR}/bin CACHE PATH "Target for the binaries")
SET(grd_LIBRARY_OUTPUT_DIRECTORY ${grd_SOURCE_DIR}/lib CACHE PATH "Target for the libraries")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${grd_LIBRARY_OUTPUT_DIRECTORY})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${grd_LIBRARY_OUTPUT_DIRECTORY})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${grd_RUNTIME_OUTPUT_DIRECTORY})
SET(grd_ENABLE_PCL 0)


include_directories(include)

# Eigen library v. 3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Boost
find_package(Boost COMPONENTS iostreams filesystem regex system REQUIRED)
#find_package(Boost COMPONENTS lexical_cast)
include_directories(${Boost_INCLUDE_DIRS}) 

# Lists of grd library files
set(grd_LIBRARY_SOURCE_LIST
  src/BBOptimizer1d.cpp
  src/CorrespondenceGraph.cpp
  src/DistributionSignature.cpp
  src/functions.cpp
  src/GlarotSignature.cpp
  src/GnuplotVisualizer.cpp
  src/KernelSignature.cpp
  src/ParamMap.cpp
  src/Profiler.cpp
  src/SkipExtractor.cpp
  src/TfUtils.cpp
  src/fileutils.cpp
)

# The "object library" target: compiles the sources only once
add_library(grdobj OBJECT ${grd_LIBRARY_SOURCE_LIST})

# shared libraries need PIC
set_property(TARGET grdobj PROPERTY POSITION_INDEPENDENT_CODE 1)
 
add_library(grd SHARED $<TARGET_OBJECTS:grdobj>)
add_library(grd-static STATIC $<TARGET_OBJECTS:grdobj>)

# Executables:

add_executable(extractSkip app/extractSkip.cpp)
target_link_libraries(extractSkip grd boost_system boost_filesystem boost_iostreams)

add_executable(compareSignatures app/compareSignatures.cpp)
target_link_libraries(compareSignatures grd boost_system boost_filesystem boost_iostreams)

add_executable(computeLoopClosureStats app/computeLoopClosureStats.cpp)
target_link_libraries(computeLoopClosureStats grd boost_system boost_filesystem boost_iostreams)

add_executable(testGlarot app/testGlarot.cpp)
target_link_libraries(testGlarot grd boost_system boost_filesystem boost_iostreams)


# Option "make install": copy binaries
INSTALL(TARGETS grd 
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

# Option "make install": copy headers
FILE(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/include/grd/*.h")
FILE(GLOB headers_3rdparty "${CMAKE_CURRENT_SOURCE_DIR}/include/grd/thirdparty/*.h")
INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/grd/)
INSTALL(FILES ${headers_3rdparty} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/grd/thirdparty/)

# Option "make install": copy cmake script
FILE(GLOB cmake_script "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/grdConfig.cmake")
message(STATUS "cmake_script " ${cmake_script})
INSTALL(FILES ${cmake_script} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/grd/)


