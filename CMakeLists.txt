cmake_minimum_required(VERSION 2.8)
project(grd)
add_definitions(-std=c++0x)
set(CMAKE_BUILD_TYPE Release)


include_directories(include)

# Eigen library v. 3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Boost
find_package(Boost COMPONENTS iostreams filesystem regex system REQUIRED)
#find_package(Boost COMPONENTS lexical_cast)
include_directories(${Boost_INCLUDE_DIRS})  


add_library(grd
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
  src/thirdparty/mcqd.cpp
)

add_executable(extractSkip app/extractSkip.cpp)
target_link_libraries(extractSkip grd boost_system boost_filesystem boost_iostreams)

add_executable(compareSignatures app/compareSignatures.cpp)
target_link_libraries(compareSignatures grd boost_system boost_filesystem boost_iostreams)

add_executable(computeLoopClosureStats app/computeLoopClosureStats.cpp)
target_link_libraries(computeLoopClosureStats grd boost_system boost_filesystem boost_iostreams)

add_executable(testGlarot app/testGlarot.cpp)
target_link_libraries(testGlarot grd boost_system boost_filesystem boost_iostreams)
