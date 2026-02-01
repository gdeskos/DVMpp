# Regression Tests CMakeLists.txt

# Test executable for cylinder flow regression
add_executable(test_CylinderFlow test_CylinderFlow.cpp)
target_link_libraries(test_CylinderFlow PRIVATE
    DVMpp_lib
    GTest::gtest_main
)
target_compile_definitions(test_CylinderFlow PRIVATE
    TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../fixtures"
    GOLDEN_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/golden_data"
    EXAMPLES_DIR="${PROJECT_SOURCE_DIR}/examples"
)
gtest_discover_tests(test_CylinderFlow)
