# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/luca/software/cmake/cmake-3.17.0-rc3-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /home/luca/software/cmake/cmake-3.17.0-rc3-Linux-x86_64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/luca/source/GP3/1D

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/luca/source/GP3/1D/build

# Include any dependencies generated for this target.
include external/eigen/bench/spbench/CMakeFiles/spsolver.dir/depend.make

# Include the progress variables for this target.
include external/eigen/bench/spbench/CMakeFiles/spsolver.dir/progress.make

# Include the compile flags for this target's objects.
include external/eigen/bench/spbench/CMakeFiles/spsolver.dir/flags.make

external/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o: external/eigen/bench/spbench/CMakeFiles/spsolver.dir/flags.make
external/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o: ../external/eigen/bench/spbench/sp_solver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luca/source/GP3/1D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o"
	cd /home/luca/source/GP3/1D/build/external/eigen/bench/spbench && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/spsolver.dir/sp_solver.cpp.o -c /home/luca/source/GP3/1D/external/eigen/bench/spbench/sp_solver.cpp

external/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/spsolver.dir/sp_solver.cpp.i"
	cd /home/luca/source/GP3/1D/build/external/eigen/bench/spbench && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/luca/source/GP3/1D/external/eigen/bench/spbench/sp_solver.cpp > CMakeFiles/spsolver.dir/sp_solver.cpp.i

external/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/spsolver.dir/sp_solver.cpp.s"
	cd /home/luca/source/GP3/1D/build/external/eigen/bench/spbench && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/luca/source/GP3/1D/external/eigen/bench/spbench/sp_solver.cpp -o CMakeFiles/spsolver.dir/sp_solver.cpp.s

# Object files for target spsolver
spsolver_OBJECTS = \
"CMakeFiles/spsolver.dir/sp_solver.cpp.o"

# External object files for target spsolver
spsolver_EXTERNAL_OBJECTS =

external/eigen/bench/spbench/spsolver: external/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o
external/eigen/bench/spbench/spsolver: external/eigen/bench/spbench/CMakeFiles/spsolver.dir/build.make
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libcholmod.so
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libamd.so
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libcolamd.so
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libcamd.so
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libccolamd.so
external/eigen/bench/spbench/spsolver: external/eigen/blas/libeigen_blas_static.a
external/eigen/bench/spbench/spsolver: external/eigen/lapack/libeigen_lapack_static.a
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libumfpack.so
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libcolamd.so
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libamd.so
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libcholmod.so
external/eigen/bench/spbench/spsolver: external/eigen/blas/libeigen_blas_static.a
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/librt.so
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libcamd.so
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libccolamd.so
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libumfpack.so
external/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/librt.so
external/eigen/bench/spbench/spsolver: external/eigen/bench/spbench/CMakeFiles/spsolver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/luca/source/GP3/1D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable spsolver"
	cd /home/luca/source/GP3/1D/build/external/eigen/bench/spbench && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/spsolver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/eigen/bench/spbench/CMakeFiles/spsolver.dir/build: external/eigen/bench/spbench/spsolver

.PHONY : external/eigen/bench/spbench/CMakeFiles/spsolver.dir/build

external/eigen/bench/spbench/CMakeFiles/spsolver.dir/clean:
	cd /home/luca/source/GP3/1D/build/external/eigen/bench/spbench && $(CMAKE_COMMAND) -P CMakeFiles/spsolver.dir/cmake_clean.cmake
.PHONY : external/eigen/bench/spbench/CMakeFiles/spsolver.dir/clean

external/eigen/bench/spbench/CMakeFiles/spsolver.dir/depend:
	cd /home/luca/source/GP3/1D/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/luca/source/GP3/1D /home/luca/source/GP3/1D/external/eigen/bench/spbench /home/luca/source/GP3/1D/build /home/luca/source/GP3/1D/build/external/eigen/bench/spbench /home/luca/source/GP3/1D/build/external/eigen/bench/spbench/CMakeFiles/spsolver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/eigen/bench/spbench/CMakeFiles/spsolver.dir/depend

