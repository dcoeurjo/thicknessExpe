# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/williamaufort/Documents/Depot/thickness/localErrorVol

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/williamaufort/Documents/Depot/thickness/localErrorVol/build

# Include any dependencies generated for this target.
include CMakeFiles/localErrorVol.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/localErrorVol.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/localErrorVol.dir/flags.make

CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o: CMakeFiles/localErrorVol.dir/flags.make
CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o: ../localErrorVol.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/williamaufort/Documents/Depot/thickness/localErrorVol/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o -c /home/williamaufort/Documents/Depot/thickness/localErrorVol/localErrorVol.cpp

CMakeFiles/localErrorVol.dir/localErrorVol.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/localErrorVol.dir/localErrorVol.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/williamaufort/Documents/Depot/thickness/localErrorVol/localErrorVol.cpp > CMakeFiles/localErrorVol.dir/localErrorVol.cpp.i

CMakeFiles/localErrorVol.dir/localErrorVol.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/localErrorVol.dir/localErrorVol.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/williamaufort/Documents/Depot/thickness/localErrorVol/localErrorVol.cpp -o CMakeFiles/localErrorVol.dir/localErrorVol.cpp.s

CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o.requires:
.PHONY : CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o.requires

CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o.provides: CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o.requires
	$(MAKE) -f CMakeFiles/localErrorVol.dir/build.make CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o.provides.build
.PHONY : CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o.provides

CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o.provides.build: CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o

# Object files for target localErrorVol
localErrorVol_OBJECTS = \
"CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o"

# External object files for target localErrorVol
localErrorVol_EXTERNAL_OBJECTS =

localErrorVol: CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o
localErrorVol: CMakeFiles/localErrorVol.dir/build.make
localErrorVol: /usr/local/lib/libmpfr.so
localErrorVol: /usr/lib/x86_64-linux-gnu/libgmp.so
localErrorVol: /usr/local/lib/libCGAL_Core.so
localErrorVol: /usr/local/lib/libCGAL.so
localErrorVol: /usr/local/lib/libboost_thread.so
localErrorVol: /usr/local/lib/libboost_system.so
localErrorVol: /usr/lib/x86_64-linux-gnu/libpthread.so
localErrorVol: /usr/local/lib/libCGAL_Core.so
localErrorVol: /usr/local/lib/libCGAL.so
localErrorVol: /usr/local/lib/libboost_thread.so
localErrorVol: /usr/local/lib/libboost_system.so
localErrorVol: /usr/lib/x86_64-linux-gnu/libpthread.so
localErrorVol: CMakeFiles/localErrorVol.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable localErrorVol"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/localErrorVol.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/localErrorVol.dir/build: localErrorVol
.PHONY : CMakeFiles/localErrorVol.dir/build

CMakeFiles/localErrorVol.dir/requires: CMakeFiles/localErrorVol.dir/localErrorVol.cpp.o.requires
.PHONY : CMakeFiles/localErrorVol.dir/requires

CMakeFiles/localErrorVol.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/localErrorVol.dir/cmake_clean.cmake
.PHONY : CMakeFiles/localErrorVol.dir/clean

CMakeFiles/localErrorVol.dir/depend:
	cd /home/williamaufort/Documents/Depot/thickness/localErrorVol/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/williamaufort/Documents/Depot/thickness/localErrorVol /home/williamaufort/Documents/Depot/thickness/localErrorVol /home/williamaufort/Documents/Depot/thickness/localErrorVol/build /home/williamaufort/Documents/Depot/thickness/localErrorVol/build /home/williamaufort/Documents/Depot/thickness/localErrorVol/build/CMakeFiles/localErrorVol.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/localErrorVol.dir/depend

