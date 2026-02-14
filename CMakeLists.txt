cmake_minimum_required(VERSION 3.15)
project(TriangleMeshDistance)
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(tests)

if(NOT DEFINED CMAKE_INSTALL_INCLUDEDIR)
    SET(CMAKE_INSTALL_INCLUDEDIR include CACHE
        PATH "Output directory for header files")
endif()

add_library(trianglemeshdistance INTERFACE)

target_include_directories(trianglemeshdistance
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/TriangleMeshDistance/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(
    DIRECTORY TriangleMeshDistance/include/
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

install(
    TARGETS trianglemeshdistance
    EXPORT TriangleMeshDistanceTargets
)

install(EXPORT TriangleMeshDistanceTargets
    FILE trianglemeshdistance-config.cmake
    NAMESPACE trianglemeshdistance::
    DESTINATION share/trianglemeshdistance
)
