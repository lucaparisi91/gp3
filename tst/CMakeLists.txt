
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


add_executable(tstOMP openMPTests/main.cpp ../src/timers.cpp )
target_include_directories(tstOMP PUBLIC
                           "${PROJECT_SOURCE_DIR}/src"
                           )
target_link_libraries(tstOMP PUBLIC MPI::MPI_CXX)
target_link_libraries(tstOMP PUBLIC stdc++fs)
if(OpenMP_CXX_FOUND)
    target_link_libraries(tstOMP PUBLIC OpenMP::OpenMP_CXX)
endif()
set_target_properties(tstOMP PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS_ADDITIONAL} )
