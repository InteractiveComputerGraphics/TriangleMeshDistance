cmake_minimum_required(VERSION 3.0)

set( EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../extern/Eigen )
if( NOT EIGEN3_INCLUDE_DIR )
    message( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif()

set(SOURCE_FILES
    test_main.cpp
    tests.cpp
	../TriangleMeshDistance/include/tmd/TriangleMeshDistance.h
	)

add_executable(tests ${SOURCE_FILES})
target_include_directories(tests PUBLIC ../TriangleMeshDistance/include )
target_include_directories(tests PUBLIC ${EIGEN3_INCLUDE_DIR} )
