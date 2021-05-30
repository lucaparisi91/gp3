
set (MPI_TOOLS_SR openMPTests/mpiTools.cpp )


add_library(mpi_tools ${MPI_TOOLS_SR}  )
target_link_libraries(mpi_tools PUBLIC MPI::MPI_CXX)
target_include_directories(mpi_tools PUBLIC
${PROJECT_SOURCE_DIR}/external/json/include
			   )
set_target_properties(mpi_tools PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS_ADDITIONAL} )



set (TEST_SOURCES test-wavefunction.cpp openMPTests/parallelTests.cpp ) 

add_executable(tst main.cpp ${TEST_SOURCES})
target_link_libraries(tst PRIVATE  gtest  PRIVATE AMReX::amrex AMReX::Flags_CXX )
target_link_libraries(tst PUBLIC  gp_lib mpi_tools)
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
target_link_libraries(tstOMP PUBLIC MPI::MPI_CXX mpi_tools)
target_link_libraries(tstOMP PUBLIC stdc++fs)
if(OpenMP_CXX_FOUND)
    target_link_libraries(tstOMP PUBLIC OpenMP::OpenMP_CXX)
endif()
set_target_properties(tstOMP PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS_ADDITIONAL} )

set (TEST_PARALLEL_SOURCES openMPTests/parallelTests.cpp )

add_executable(tstParallel openMPTests/mainTestsParallel.cpp ${TEST_PARALLEL_SOURCES})
target_link_libraries(tstParallel PRIVATE gtest )
target_link_libraries(tstParallel PUBLIC  mpi_tools)
target_include_directories(tstParallel PUBLIC
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