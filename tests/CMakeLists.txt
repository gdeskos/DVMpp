# DVMpp Test Suite CMakeLists.txt

include(FetchContent)

# Fetch GoogleTest
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

# Enable testing
enable_testing()
include(GoogleTest)

# Add test subdirectories
add_subdirectory(unit)
add_subdirectory(integration)
add_subdirectory(regression)

# Copy test fixtures to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/fixtures
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
