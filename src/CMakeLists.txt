# DVMpp source CMakeLists.txt

# Define source files
set(DVM_SOURCE_FILES
    DVMBase.cpp
    main.cpp
    VortexBlobs.cpp
    VortexSheet.cpp
    Probe.cpp
    XmlHandler.cpp
    Random.cpp
    pugixml.cpp
)

# Define header files for IDE support
set(DVM_HEADER_FILES
    BaseTypes.hpp
    DVMBase.hpp
    Exceptions.hpp
    Probe.hpp
    Random.hpp
    VelocityKernel.hpp
    VortexBlobs.hpp
    VortexSheet.hpp
    XmlHandler.hpp
    pugiconfig.hpp
    pugixml.hpp
)

# Add VelocityKernel implementation
list(APPEND DVM_SOURCE_FILES VelocityKernel.cpp)

# Create executable
add_executable(DVMpp ${DVM_SOURCE_FILES} ${DVM_HEADER_FILES})

# Require C++17
target_compile_features(DVMpp PRIVATE cxx_std_17)

# Set compiler options
target_compile_options(DVMpp PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-O3 -march=native -Wall -Wextra -pipe>
    $<$<CXX_COMPILER_ID:MSVC>:/O2 /W4>
)

# Find and link Armadillo
find_package(Armadillo REQUIRED)
if(ARMADILLO_VERSION_STRING VERSION_LESS "7.6")
    message(FATAL_ERROR "Armadillo version must be at least 7.6\nFound: ${ARMADILLO_VERSION_STRING}")
endif()
target_include_directories(DVMpp PRIVATE ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(DVMpp PRIVATE ${ARMADILLO_LIBRARIES})

# Find and link OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(DVMpp PRIVATE OpenMP::OpenMP_CXX)

# Find and link BLAS
find_package(BLAS REQUIRED)
target_link_libraries(DVMpp PRIVATE ${BLAS_LIBRARIES})

# Find and link Boost
find_package(Boost 1.56.0 REQUIRED COMPONENTS program_options)
target_include_directories(DVMpp PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(DVMpp PRIVATE Boost::program_options)
target_compile_definitions(DVMpp PRIVATE HAS_BOOST)

# Release mode optimizations
target_compile_definitions(DVMpp PRIVATE
    $<$<CONFIG:Release>:ARMA_NO_DEBUG>
)

# FMM support
if(USE_FMM)
    include(FetchContent)
    FetchContent_Declare(exafmm
        GIT_REPOSITORY https://github.com/exafmm/exafmm-t.git
        GIT_TAG master
    )
    FetchContent_MakeAvailable(exafmm)
    target_link_libraries(DVMpp PRIVATE exafmm)
    target_compile_definitions(DVMpp PRIVATE USE_FMM)
    message(STATUS "FMM support enabled via exafmm-t")
endif()

# Create a library target for testing (without main)
set(DVM_LIB_SOURCES
    DVMBase.cpp
    VortexBlobs.cpp
    VortexSheet.cpp
    Probe.cpp
    XmlHandler.cpp
    Random.cpp
    VelocityKernel.cpp
    pugixml.cpp
)

add_library(DVMpp_lib STATIC ${DVM_LIB_SOURCES})
target_compile_features(DVMpp_lib PUBLIC cxx_std_17)
target_include_directories(DVMpp_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(DVMpp_lib PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(DVMpp_lib PUBLIC
    ${ARMADILLO_LIBRARIES}
    OpenMP::OpenMP_CXX
    ${BLAS_LIBRARIES}
    Boost::program_options
)
target_compile_definitions(DVMpp_lib PUBLIC HAS_BOOST)

if(USE_FMM)
    target_link_libraries(DVMpp_lib PUBLIC exafmm)
    target_compile_definitions(DVMpp_lib PUBLIC USE_FMM)
endif()
