cmake_minimum_required(VERSION 3.10)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_BUILD_TYPE Debug)


# set the project name
project( radialPDE   VERSION 3.0)



add_subdirectory(${PROJECT_BINARY_DIR}/../external/googletest)
add_subdirectory(${PROJECT_BINARY_DIR}/../external/eigen)


set(SOURCES  src/geometry.cpp src/kineticEnergy.cpp src/potential.cpp src/functional.cpp)


add_library(main_lib STATIC ${SOURCES})
set_target_properties(main_lib PROPERTIES RULE_MESSAGES ON)

target_include_directories(main_lib PUBLIC
                           "${PROJECT_BINARY_DIR}/../src"
			   )


target_include_directories(main_lib SYSTEM PUBLIC ${PROJECT_BINARY_DIR}/../external/eigen ${PROJECT_BINARY_DIR}/../external/json/single_include		           )
			 
target_link_libraries (main_lib eigen)

add_executable(main src/main.cpp)
target_link_libraries(main main_lib)
target_link_libraries(main )


set(TEST_SOURCES tst/base.cpp )
add_executable(tst tst/main.cpp ${TEST_SOURCES})
target_include_directories(tst PUBLIC
                           "${PROJECT_BINARY_DIR}/../src"
                           )
target_link_libraries(tst PUBLIC ${PROJECT_BINARY_DIR}/../external gtest eigen main_lib ${MPI_LIBRARIES})
