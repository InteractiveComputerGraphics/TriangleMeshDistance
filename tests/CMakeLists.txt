cmake_minimum_required(VERSION 3.15)

set(SOURCE_FILES
    test_main.cpp
    tests.cpp
	../TriangleMeshDistance/include/tmd/TriangleMeshDistance.h
	)

add_executable(tests ${SOURCE_FILES})
target_include_directories(tests PUBLIC ../TriangleMeshDistance/include )
