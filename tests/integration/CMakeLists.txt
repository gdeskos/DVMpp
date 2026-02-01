# Integration Tests CMakeLists.txt

# Test executable for Biot-Savart integration
add_executable(test_BiotSavart test_BiotSavart.cpp)
target_link_libraries(test_BiotSavart PRIVATE
    DVMpp_lib
    GTest::gtest_main
)
target_compile_definitions(test_BiotSavart PRIVATE
    TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../fixtures"
)
gtest_discover_tests(test_BiotSavart)
