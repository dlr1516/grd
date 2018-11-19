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
  src/GlareSignature.cpp
  src/GnuplotVisualizer.cpp
  src/KernelSignature.cpp
  src/ParamMap.cpp
  src/Profiler.cpp
  src/SkipExtractor.cpp
  src/TfUtils.cpp
  src/fileutils.cpp
  src/thirdparty/mcqd.cpp
)

#add_executable(extractSkip app/extractSkip.cpp)
#target_link_libraries(extractSkip cglarot boost_system boost_filesystem boost_iostreams)

#add_executable(drawSkip app/drawSkip.cpp)
#target_link_libraries(drawSkip cglarot boost_system boost_filesystem boost_iostreams)

#add_executable(compareSignatures app/compareSignatures.cpp)
#target_link_libraries(compareSignatures cglarot boost_system boost_filesystem boost_iostreams)

#add_executable(testApproximations app/testApproximations.cpp)
#target_link_libraries(testApproximations cglarot boost_system boost_filesystem boost_iostreams)

#add_executable(testDistributionSignature app/testDistributionSignature.cpp)
#target_link_libraries(testDistributionSignature cglarot boost_system boost_filesystem boost_iostreams)

#add_executable(drawDistributionSignature app/drawDistributionSignature.cpp)
#target_link_libraries(drawDistributionSignature cglarot boost_system boost_filesystem boost_iostreams)

#add_executable(computeLoopClosureStats app/computeLoopClosureStats.cpp)
#target_link_libraries(computeLoopClosureStats cglarot boost_system boost_filesystem boost_iostreams)



#add_executable(testCorrespondenceGraph app/testCorrespondenceGraph.cpp)
#target_link_libraries(testCorrespondenceGraph cglarot boost_system boost_filesystem boost_iostreams)

#add_executable(testGlarot app/testGlarot.cpp)
#target_link_libraries(testGlarot cglarot boost_system boost_filesystem boost_iostreams)
