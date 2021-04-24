
#library of gp routines
set(GP_FILES wavefunction.cpp geometry.cpp tools.cpp functional.cpp stepper.cpp normalization.cpp gpDriver.cpp timers.cpp functionalFactory.cpp)

add_library(gp_lib ${GP_FILES} )
target_link_libraries(gp_lib PRIVATE AMReX::amrex AMReX::Flags_CXX)
target_link_libraries(gp_lib PUBLIC MPI::MPI_CXX)
target_include_directories(gp_lib PUBLIC
${PROJECT_SOURCE_DIR}/external/json/include
			   )
set_target_properties(gp_lib PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS_ADDITIONAL} )


# exectuable gp
add_executable(gp main.cpp )
set_target_properties(gp PROPERTIES LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS_ADDITIONAL} )
target_link_libraries(gp PUBLIC gp_lib MPI::MPI_CXX)
target_link_libraries(gp PRIVATE AMReX::amrex AMReX::Flags_CXX  )
target_link_libraries(gp PUBLIC stdc++fs)
target_include_directories(gp PUBLIC
${PROJECT_SOURCE_DIR}/external/json/include
			   )		   
set_target_properties(gp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS_ADDITIONAL} )

if(OpenMP_CXX_FOUND)
    target_link_libraries(gp PUBLIC OpenMP::OpenMP_CXX)
endif()
