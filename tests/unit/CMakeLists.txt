# Unit Tests CMakeLists.txt

# Test executable for VortexBlobs
add_executable(test_VortexBlobs test_VortexBlobs.cpp)
target_link_libraries(test_VortexBlobs PRIVATE
    DVMpp_lib
    GTest::gtest_main
    GTest::gmock
)
gtest_discover_tests(test_VortexBlobs)

# Test executable for VortexSheet
add_executable(test_VortexSheet test_VortexSheet.cpp)
target_link_libraries(test_VortexSheet PRIVATE
    DVMpp_lib
    GTest::gtest_main
    GTest::gmock
)
gtest_discover_tests(test_VortexSheet)

# Test executable for Random
add_executable(test_Random test_Random.cpp)
target_link_libraries(test_Random PRIVATE
    DVMpp_lib
    GTest::gtest_main
)
gtest_discover_tests(test_Random)

# Test executable for XmlHandler
add_executable(test_XmlHandler test_XmlHandler.cpp)
target_link_libraries(test_XmlHandler PRIVATE
    DVMpp_lib
    GTest::gtest_main
)
target_compile_definitions(test_XmlHandler PRIVATE
    TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../fixtures"
)
gtest_discover_tests(test_XmlHandler)

# Test executable for Probe
add_executable(test_Probe test_Probe.cpp)
target_link_libraries(test_Probe PRIVATE
    DVMpp_lib
    GTest::gtest_main
)
gtest_discover_tests(test_Probe)

# Test executable for VelocityKernel
add_executable(test_VelocityKernel test_VelocityKernel.cpp)
target_link_libraries(test_VelocityKernel PRIVATE
    DVMpp_lib
    GTest::gtest_main
)
gtest_discover_tests(test_VelocityKernel)

# Test executable for BackgroundMesh
add_executable(test_BackgroundMesh test_BackgroundMesh.cpp)
target_link_libraries(test_BackgroundMesh PRIVATE
    DVMpp_lib
    GTest::gtest_main
)
gtest_discover_tests(test_BackgroundMesh)
