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
include external/eigen/test/CMakeFiles/integer_types_4.dir/depend.make

# Include the progress variables for this target.
include external/eigen/test/CMakeFiles/integer_types_4.dir/progress.make

# Include the compile flags for this target's objects.
include external/eigen/test/CMakeFiles/integer_types_4.dir/flags.make

external/eigen/test/CMakeFiles/integer_types_4.dir/integer_types.cpp.o: external/eigen/test/CMakeFiles/integer_types_4.dir/flags.make
external/eigen/test/CMakeFiles/integer_types_4.dir/integer_types.cpp.o: ../external/eigen/test/integer_types.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luca/source/GP3/1D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/eigen/test/CMakeFiles/integer_types_4.dir/integer_types.cpp.o"
	cd /home/luca/source/GP3/1D/build/external/eigen/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/integer_types_4.dir/integer_types.cpp.o -c /home/luca/source/GP3/1D/external/eigen/test/integer_types.cpp

external/eigen/test/CMakeFiles/integer_types_4.dir/integer_types.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/integer_types_4.dir/integer_types.cpp.i"
	cd /home/luca/source/GP3/1D/build/external/eigen/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/luca/source/GP3/1D/external/eigen/test/integer_types.cpp > CMakeFiles/integer_types_4.dir/integer_types.cpp.i

external/eigen/test/CMakeFiles/integer_types_4.dir/integer_types.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/integer_types_4.dir/integer_types.cpp.s"
	cd /home/luca/source/GP3/1D/build/external/eigen/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/luca/source/GP3/1D/external/eigen/test/integer_types.cpp -o CMakeFiles/integer_types_4.dir/integer_types.cpp.s

# Object files for target integer_types_4
integer_types_4_OBJECTS = \
"CMakeFiles/integer_types_4.dir/integer_types.cpp.o"

# External object files for target integer_types_4
integer_types_4_EXTERNAL_OBJECTS =

external/eigen/test/integer_types_4: external/eigen/test/CMakeFiles/integer_types_4.dir/integer_types.cpp.o
external/eigen/test/integer_types_4: external/eigen/test/CMakeFiles/integer_types_4.dir/build.make
external/eigen/test/integer_types_4: external/eigen/test/CMakeFiles/integer_types_4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/luca/source/GP3/1D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable integer_types_4"
	cd /home/luca/source/GP3/1D/build/external/eigen/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/integer_types_4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/eigen/test/CMakeFiles/integer_types_4.dir/build: external/eigen/test/integer_types_4

.PHONY : external/eigen/test/CMakeFiles/integer_types_4.dir/build

external/eigen/test/CMakeFiles/integer_types_4.dir/clean:
	cd /home/luca/source/GP3/1D/build/external/eigen/test && $(CMAKE_COMMAND) -P CMakeFiles/integer_types_4.dir/cmake_clean.cmake
.PHONY : external/eigen/test/CMakeFiles/integer_types_4.dir/clean

external/eigen/test/CMakeFiles/integer_types_4.dir/depend:
	cd /home/luca/source/GP3/1D/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/luca/source/GP3/1D /home/luca/source/GP3/1D/external/eigen/test /home/luca/source/GP3/1D/build /home/luca/source/GP3/1D/build/external/eigen/test /home/luca/source/GP3/1D/build/external/eigen/test/CMakeFiles/integer_types_4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/eigen/test/CMakeFiles/integer_types_4.dir/depend

