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
include external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/depend.make

# Include the progress variables for this target.
include external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/progress.make

# Include the compile flags for this target's objects.
include external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/flags.make

external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.o: external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/flags.make
external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.o: external/eigen/doc/snippets/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp
external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.o: ../external/eigen/doc/snippets/SelfAdjointEigenSolver_compute_MatrixType.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luca/source/GP3/1D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.o"
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.o -c /home/luca/source/GP3/1D/build/external/eigen/doc/snippets/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp

external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.i"
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/luca/source/GP3/1D/build/external/eigen/doc/snippets/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp > CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.i

external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.s"
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/luca/source/GP3/1D/build/external/eigen/doc/snippets/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp -o CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.s

# Object files for target compile_SelfAdjointEigenSolver_compute_MatrixType
compile_SelfAdjointEigenSolver_compute_MatrixType_OBJECTS = \
"CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.o"

# External object files for target compile_SelfAdjointEigenSolver_compute_MatrixType
compile_SelfAdjointEigenSolver_compute_MatrixType_EXTERNAL_OBJECTS =

external/eigen/doc/snippets/compile_SelfAdjointEigenSolver_compute_MatrixType: external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/compile_SelfAdjointEigenSolver_compute_MatrixType.cpp.o
external/eigen/doc/snippets/compile_SelfAdjointEigenSolver_compute_MatrixType: external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/build.make
external/eigen/doc/snippets/compile_SelfAdjointEigenSolver_compute_MatrixType: external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/luca/source/GP3/1D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable compile_SelfAdjointEigenSolver_compute_MatrixType"
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/link.txt --verbose=$(VERBOSE)
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && ./compile_SelfAdjointEigenSolver_compute_MatrixType >/home/luca/source/GP3/1D/build/external/eigen/doc/snippets/SelfAdjointEigenSolver_compute_MatrixType.out

# Rule to build all files generated by this target.
external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/build: external/eigen/doc/snippets/compile_SelfAdjointEigenSolver_compute_MatrixType

.PHONY : external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/build

external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/clean:
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && $(CMAKE_COMMAND) -P CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/cmake_clean.cmake
.PHONY : external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/clean

external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/depend:
	cd /home/luca/source/GP3/1D/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/luca/source/GP3/1D /home/luca/source/GP3/1D/external/eigen/doc/snippets /home/luca/source/GP3/1D/build /home/luca/source/GP3/1D/build/external/eigen/doc/snippets /home/luca/source/GP3/1D/build/external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/eigen/doc/snippets/CMakeFiles/compile_SelfAdjointEigenSolver_compute_MatrixType.dir/depend

