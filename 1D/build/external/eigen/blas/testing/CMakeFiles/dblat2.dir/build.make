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
include external/eigen/blas/testing/CMakeFiles/dblat2.dir/depend.make

# Include the progress variables for this target.
include external/eigen/blas/testing/CMakeFiles/dblat2.dir/progress.make

# Include the compile flags for this target's objects.
include external/eigen/blas/testing/CMakeFiles/dblat2.dir/flags.make

external/eigen/blas/testing/CMakeFiles/dblat2.dir/dblat2.f.o: external/eigen/blas/testing/CMakeFiles/dblat2.dir/flags.make
external/eigen/blas/testing/CMakeFiles/dblat2.dir/dblat2.f.o: ../external/eigen/blas/testing/dblat2.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luca/source/GP3/1D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object external/eigen/blas/testing/CMakeFiles/dblat2.dir/dblat2.f.o"
	cd /home/luca/source/GP3/1D/build/external/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/luca/source/GP3/1D/external/eigen/blas/testing/dblat2.f -o CMakeFiles/dblat2.dir/dblat2.f.o

external/eigen/blas/testing/CMakeFiles/dblat2.dir/dblat2.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/dblat2.dir/dblat2.f.i"
	cd /home/luca/source/GP3/1D/build/external/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/luca/source/GP3/1D/external/eigen/blas/testing/dblat2.f > CMakeFiles/dblat2.dir/dblat2.f.i

external/eigen/blas/testing/CMakeFiles/dblat2.dir/dblat2.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/dblat2.dir/dblat2.f.s"
	cd /home/luca/source/GP3/1D/build/external/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/luca/source/GP3/1D/external/eigen/blas/testing/dblat2.f -o CMakeFiles/dblat2.dir/dblat2.f.s

# Object files for target dblat2
dblat2_OBJECTS = \
"CMakeFiles/dblat2.dir/dblat2.f.o"

# External object files for target dblat2
dblat2_EXTERNAL_OBJECTS =

external/eigen/blas/testing/dblat2: external/eigen/blas/testing/CMakeFiles/dblat2.dir/dblat2.f.o
external/eigen/blas/testing/dblat2: external/eigen/blas/testing/CMakeFiles/dblat2.dir/build.make
external/eigen/blas/testing/dblat2: external/eigen/blas/libeigen_blas.so
external/eigen/blas/testing/dblat2: external/eigen/blas/testing/CMakeFiles/dblat2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/luca/source/GP3/1D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable dblat2"
	cd /home/luca/source/GP3/1D/build/external/eigen/blas/testing && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dblat2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/eigen/blas/testing/CMakeFiles/dblat2.dir/build: external/eigen/blas/testing/dblat2

.PHONY : external/eigen/blas/testing/CMakeFiles/dblat2.dir/build

external/eigen/blas/testing/CMakeFiles/dblat2.dir/clean:
	cd /home/luca/source/GP3/1D/build/external/eigen/blas/testing && $(CMAKE_COMMAND) -P CMakeFiles/dblat2.dir/cmake_clean.cmake
.PHONY : external/eigen/blas/testing/CMakeFiles/dblat2.dir/clean

external/eigen/blas/testing/CMakeFiles/dblat2.dir/depend:
	cd /home/luca/source/GP3/1D/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/luca/source/GP3/1D /home/luca/source/GP3/1D/external/eigen/blas/testing /home/luca/source/GP3/1D/build /home/luca/source/GP3/1D/build/external/eigen/blas/testing /home/luca/source/GP3/1D/build/external/eigen/blas/testing/CMakeFiles/dblat2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/eigen/blas/testing/CMakeFiles/dblat2.dir/depend

