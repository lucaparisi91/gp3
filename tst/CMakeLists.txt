
set (TEST_SOURCES test-wavefunction.cpp ) 


add_executable(tst main.cpp ${TEST_SOURCES})
target_link_libraries(tst PRIVATE  gtest  ${MPI_LIBRARIES} amrex )
target_link_libraries(tst PUBLIC  gp_lib)
target_include_directories(tst PUBLIC
                           "${PROJECT_SOURCE_DIR}/src"
                           )
target_link_libraries(tst PUBLIC MPI::MPI_CXX)
target_include_directories(gp_lib PUBLIC
${PROJECT_SOURCE_DIR}/external/json/include
			   )