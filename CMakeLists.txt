cmake_minimum_required(VERSION 3.14)

project(DVMpp
    VERSION 1.0.0
    DESCRIPTION "2D Discrete Vortex Method Code"
    LANGUAGES CXX)

# Options
option(BUILD_TESTING "Build the testing tree" ON)
option(USE_FMM "Enable Fast Multipole Method for O(N) velocity computation" OFF)
option(BUILD_DOC "Build documentation with Doxygen" ON)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# Add source subdirectory
add_subdirectory(src)

# Add documentation if requested
if(BUILD_DOC)
    add_subdirectory(doc)
endif()

# Add testing infrastructure
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
