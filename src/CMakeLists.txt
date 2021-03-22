

set(GP_FILES wavefunction.cpp geometry.cpp tools.cpp functional.cpp stepper.cpp normalization.cpp)


add_library(gp_lib ${GP_FILES} )
target_link_libraries(gp_lib PRIVATE amrex gfortran  )
target_link_libraries(gp_lib PUBLIC MPI::MPI_CXX)
target_include_directories(gp_lib PUBLIC
${PROJECT_SOURCE_DIR}/external/json/include
			   )