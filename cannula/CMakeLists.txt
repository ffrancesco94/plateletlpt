cmake_minimum_required(VERSION 2.8)

PROJECT(platelets)

set(CMAKE_CXX_STANDARD 14)

set(VTK_DIR "/scratch/fiusco/VTK-8.2.0/VTK-build/")
set(BOOST_ROOT "/scratch/fiusco/boost/")
set(EIGEN_DIR "/scratch/fiusco/Eigen/")

find_package(VTK REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem system)
include(${VTK_USE_FILE})
include_directories(${EIGEN_DIR} ${Boost_INCLUDE_DIR} ../lptmodel/)

MESSAGE(STATUS ${Boost_INCLUDE_DIR})

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

if (NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()
message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
message("CMAKE_C_FLAGS is ${CMAKE_C_FLAGS}")


include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "-std=c++14 -msse3")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()
#set(CMAKE_BUILD_TYPE Debug)

message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

add_executable(platelets MACOSX_BUNDLE ../lptmodel/BBox ../lptmodel/BVH ../lptmodel/RayTracer ../lptmodel/vtkhelpers ../lptmodel/Model ../lptmodel/CoordinateSystem ../lptmodel/Injector ../lptmodel/InputFileList ../lptmodel/Absorber ../lptmodel/ActivationModel ../lptmodel/Particle ../lptmodel/ParticleForces platelets_cannula)

if(VTK_LIBRARIES)
  target_link_libraries(platelets ${VTK_LIBRARIES} ${Boost_LIBRARIES})
else()
  target_link_libraries(platelets vtkHybrid vtkWidgets ${Boost_LIBRARIES})
endif()
