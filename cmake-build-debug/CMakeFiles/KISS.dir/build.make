# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

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
CMAKE_COMMAND = /snap/clion/149/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/149/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/vsm/CLionProjects/KISS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/vsm/CLionProjects/KISS/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/KISS.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/KISS.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/KISS.dir/flags.make

CMakeFiles/KISS.dir/main.cpp.o: CMakeFiles/KISS.dir/flags.make
CMakeFiles/KISS.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vsm/CLionProjects/KISS/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/KISS.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/KISS.dir/main.cpp.o -c /home/vsm/CLionProjects/KISS/main.cpp

CMakeFiles/KISS.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KISS.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vsm/CLionProjects/KISS/main.cpp > CMakeFiles/KISS.dir/main.cpp.i

CMakeFiles/KISS.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KISS.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vsm/CLionProjects/KISS/main.cpp -o CMakeFiles/KISS.dir/main.cpp.s

# Object files for target KISS
KISS_OBJECTS = \
"CMakeFiles/KISS.dir/main.cpp.o"

# External object files for target KISS
KISS_EXTERNAL_OBJECTS =

KISS: CMakeFiles/KISS.dir/main.cpp.o
KISS: CMakeFiles/KISS.dir/build.make
KISS: CMakeFiles/KISS.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/vsm/CLionProjects/KISS/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable KISS"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/KISS.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/KISS.dir/build: KISS

.PHONY : CMakeFiles/KISS.dir/build

CMakeFiles/KISS.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/KISS.dir/cmake_clean.cmake
.PHONY : CMakeFiles/KISS.dir/clean

CMakeFiles/KISS.dir/depend:
	cd /home/vsm/CLionProjects/KISS/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vsm/CLionProjects/KISS /home/vsm/CLionProjects/KISS /home/vsm/CLionProjects/KISS/cmake-build-debug /home/vsm/CLionProjects/KISS/cmake-build-debug /home/vsm/CLionProjects/KISS/cmake-build-debug/CMakeFiles/KISS.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/KISS.dir/depend
