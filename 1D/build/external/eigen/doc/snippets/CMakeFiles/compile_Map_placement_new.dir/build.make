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
include external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/depend.make

# Include the progress variables for this target.
include external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/progress.make

# Include the compile flags for this target's objects.
include external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/flags.make

external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.o: external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/flags.make
external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.o: external/eigen/doc/snippets/compile_Map_placement_new.cpp
external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.o: ../external/eigen/doc/snippets/Map_placement_new.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luca/source/GP3/1D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.o"
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.o -c /home/luca/source/GP3/1D/build/external/eigen/doc/snippets/compile_Map_placement_new.cpp

external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.i"
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/luca/source/GP3/1D/build/external/eigen/doc/snippets/compile_Map_placement_new.cpp > CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.i

external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.s"
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/luca/source/GP3/1D/build/external/eigen/doc/snippets/compile_Map_placement_new.cpp -o CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.s

# Object files for target compile_Map_placement_new
compile_Map_placement_new_OBJECTS = \
"CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.o"

# External object files for target compile_Map_placement_new
compile_Map_placement_new_EXTERNAL_OBJECTS =

external/eigen/doc/snippets/compile_Map_placement_new: external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/compile_Map_placement_new.cpp.o
external/eigen/doc/snippets/compile_Map_placement_new: external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/build.make
external/eigen/doc/snippets/compile_Map_placement_new: external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/luca/source/GP3/1D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable compile_Map_placement_new"
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compile_Map_placement_new.dir/link.txt --verbose=$(VERBOSE)
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && ./compile_Map_placement_new >/home/luca/source/GP3/1D/build/external/eigen/doc/snippets/Map_placement_new.out

# Rule to build all files generated by this target.
external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/build: external/eigen/doc/snippets/compile_Map_placement_new

.PHONY : external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/build

external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/clean:
	cd /home/luca/source/GP3/1D/build/external/eigen/doc/snippets && $(CMAKE_COMMAND) -P CMakeFiles/compile_Map_placement_new.dir/cmake_clean.cmake
.PHONY : external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/clean

external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/depend:
	cd /home/luca/source/GP3/1D/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/luca/source/GP3/1D /home/luca/source/GP3/1D/external/eigen/doc/snippets /home/luca/source/GP3/1D/build /home/luca/source/GP3/1D/build/external/eigen/doc/snippets /home/luca/source/GP3/1D/build/external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/eigen/doc/snippets/CMakeFiles/compile_Map_placement_new.dir/depend

