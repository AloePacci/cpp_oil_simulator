# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/aloepacci/trabajo/basic_simulator

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/aloepacci/trabajo/basic_simulator/basic_sim/build

# Include any dependencies generated for this target.
include basic_sim/CMakeFiles/BasicSimulator.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include basic_sim/CMakeFiles/BasicSimulator.dir/compiler_depend.make

# Include the progress variables for this target.
include basic_sim/CMakeFiles/BasicSimulator.dir/progress.make

# Include the compile flags for this target's objects.
include basic_sim/CMakeFiles/BasicSimulator.dir/flags.make

basic_sim/CMakeFiles/BasicSimulator.dir/main.cpp.o: basic_sim/CMakeFiles/BasicSimulator.dir/flags.make
basic_sim/CMakeFiles/BasicSimulator.dir/main.cpp.o: ../main.cpp
basic_sim/CMakeFiles/BasicSimulator.dir/main.cpp.o: basic_sim/CMakeFiles/BasicSimulator.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aloepacci/trabajo/basic_simulator/basic_sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object basic_sim/CMakeFiles/BasicSimulator.dir/main.cpp.o"
	cd /home/aloepacci/trabajo/basic_simulator/basic_sim/build/basic_sim && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT basic_sim/CMakeFiles/BasicSimulator.dir/main.cpp.o -MF CMakeFiles/BasicSimulator.dir/main.cpp.o.d -o CMakeFiles/BasicSimulator.dir/main.cpp.o -c /home/aloepacci/trabajo/basic_simulator/basic_sim/main.cpp

basic_sim/CMakeFiles/BasicSimulator.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BasicSimulator.dir/main.cpp.i"
	cd /home/aloepacci/trabajo/basic_simulator/basic_sim/build/basic_sim && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aloepacci/trabajo/basic_simulator/basic_sim/main.cpp > CMakeFiles/BasicSimulator.dir/main.cpp.i

basic_sim/CMakeFiles/BasicSimulator.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BasicSimulator.dir/main.cpp.s"
	cd /home/aloepacci/trabajo/basic_simulator/basic_sim/build/basic_sim && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aloepacci/trabajo/basic_simulator/basic_sim/main.cpp -o CMakeFiles/BasicSimulator.dir/main.cpp.s

# Object files for target BasicSimulator
BasicSimulator_OBJECTS = \
"CMakeFiles/BasicSimulator.dir/main.cpp.o"

# External object files for target BasicSimulator
BasicSimulator_EXTERNAL_OBJECTS =

basic_sim/BasicSimulator: basic_sim/CMakeFiles/BasicSimulator.dir/main.cpp.o
basic_sim/BasicSimulator: basic_sim/CMakeFiles/BasicSimulator.dir/build.make
basic_sim/BasicSimulator: utils/map/libmap.a
basic_sim/BasicSimulator: utils/basic_simulator/libbasic_simulator.a
basic_sim/BasicSimulator: utils/map/libmap.a
basic_sim/BasicSimulator: basic_sim/CMakeFiles/BasicSimulator.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/aloepacci/trabajo/basic_simulator/basic_sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable BasicSimulator"
	cd /home/aloepacci/trabajo/basic_simulator/basic_sim/build/basic_sim && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BasicSimulator.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
basic_sim/CMakeFiles/BasicSimulator.dir/build: basic_sim/BasicSimulator
.PHONY : basic_sim/CMakeFiles/BasicSimulator.dir/build

basic_sim/CMakeFiles/BasicSimulator.dir/clean:
	cd /home/aloepacci/trabajo/basic_simulator/basic_sim/build/basic_sim && $(CMAKE_COMMAND) -P CMakeFiles/BasicSimulator.dir/cmake_clean.cmake
.PHONY : basic_sim/CMakeFiles/BasicSimulator.dir/clean

basic_sim/CMakeFiles/BasicSimulator.dir/depend:
	cd /home/aloepacci/trabajo/basic_simulator/basic_sim/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/aloepacci/trabajo/basic_simulator /home/aloepacci/trabajo/basic_simulator/basic_sim /home/aloepacci/trabajo/basic_simulator/basic_sim/build /home/aloepacci/trabajo/basic_simulator/basic_sim/build/basic_sim /home/aloepacci/trabajo/basic_simulator/basic_sim/build/basic_sim/CMakeFiles/BasicSimulator.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : basic_sim/CMakeFiles/BasicSimulator.dir/depend

