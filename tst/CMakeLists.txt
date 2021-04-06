
set (TEST_SOURCES test-wavefunction.cpp ) 

add_executable(tst main.cpp ${TEST_SOURCES})
target_link_libraries(tst PRIVATE  gtest  PRIVATE AMReX::amrex AMReX::Flags_CXX )
target_link_libraries(tst PUBLIC  gp_lib)
target_include_directories(tst PUBLIC
                           "${PROJECT_SOURCE_DIR}/src"
                           )
target_link_libraries(tst PUBLIC MPI::MPI_CXX)
target_link_libraries(tst PUBLIC stdc++fs)
target_include_directories(tst PUBLIC
${PROJECT_SOURCE_DIR}/external/json/include
			   )
if(OpenMP_CXX_FOUND)
    target_link_libraries(tst PUBLIC OpenMP::OpenMP_CXX)
endif()
           