set(DVM_SOURCE_FILES BaseTypes.cpp DVMBase.cpp main.cpp pugixml.cpp)

add_executable(DVM++ ${DVM_SOURCE_FILES})

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(DVM++ ${ARMADILLO_LIBRARIES})

include(FindOpenMP)
if (NOT OpenMP_CXX_FLAGS)
	message(FATAL_ERROR "The compiler does not have OpenMP support")
endif()

include(FindBLAS)
if (NOT BLAS_FOUND)
	message(FATAL_ERROR "Cannot find BLAS support")
endif()
target_link_libraries(DVM++ ${BLAS_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -Wall -Wextra -pipe ${OpenMP_CXX_FLAGS}")
