# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/share/cmake-3.5.0-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /usr/share/cmake-3.5.0-Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mariano/BrainToolsCARM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mariano/BrainToolsCARM/QtCreator

# Include any dependencies generated for this target.
include CMakeFiles/CoregTool.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/CoregTool.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/CoregTool.dir/flags.make

CMakeFiles/CoregTool.dir/coregmain.o: CMakeFiles/CoregTool.dir/flags.make
CMakeFiles/CoregTool.dir/coregmain.o: ../coregmain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mariano/BrainToolsCARM/QtCreator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/CoregTool.dir/coregmain.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CoregTool.dir/coregmain.o -c /home/mariano/BrainToolsCARM/coregmain.cpp

CMakeFiles/CoregTool.dir/coregmain.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CoregTool.dir/coregmain.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mariano/BrainToolsCARM/coregmain.cpp > CMakeFiles/CoregTool.dir/coregmain.i

CMakeFiles/CoregTool.dir/coregmain.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CoregTool.dir/coregmain.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mariano/BrainToolsCARM/coregmain.cpp -o CMakeFiles/CoregTool.dir/coregmain.s

CMakeFiles/CoregTool.dir/coregmain.o.requires:

.PHONY : CMakeFiles/CoregTool.dir/coregmain.o.requires

CMakeFiles/CoregTool.dir/coregmain.o.provides: CMakeFiles/CoregTool.dir/coregmain.o.requires
	$(MAKE) -f CMakeFiles/CoregTool.dir/build.make CMakeFiles/CoregTool.dir/coregmain.o.provides.build
.PHONY : CMakeFiles/CoregTool.dir/coregmain.o.provides

CMakeFiles/CoregTool.dir/coregmain.o.provides.build: CMakeFiles/CoregTool.dir/coregmain.o


CMakeFiles/CoregTool.dir/brainregistration.o: CMakeFiles/CoregTool.dir/flags.make
CMakeFiles/CoregTool.dir/brainregistration.o: ../brainregistration.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mariano/BrainToolsCARM/QtCreator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/CoregTool.dir/brainregistration.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CoregTool.dir/brainregistration.o -c /home/mariano/BrainToolsCARM/brainregistration.cpp

CMakeFiles/CoregTool.dir/brainregistration.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CoregTool.dir/brainregistration.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mariano/BrainToolsCARM/brainregistration.cpp > CMakeFiles/CoregTool.dir/brainregistration.i

CMakeFiles/CoregTool.dir/brainregistration.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CoregTool.dir/brainregistration.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mariano/BrainToolsCARM/brainregistration.cpp -o CMakeFiles/CoregTool.dir/brainregistration.s

CMakeFiles/CoregTool.dir/brainregistration.o.requires:

.PHONY : CMakeFiles/CoregTool.dir/brainregistration.o.requires

CMakeFiles/CoregTool.dir/brainregistration.o.provides: CMakeFiles/CoregTool.dir/brainregistration.o.requires
	$(MAKE) -f CMakeFiles/CoregTool.dir/build.make CMakeFiles/CoregTool.dir/brainregistration.o.provides.build
.PHONY : CMakeFiles/CoregTool.dir/brainregistration.o.provides

CMakeFiles/CoregTool.dir/brainregistration.o.provides.build: CMakeFiles/CoregTool.dir/brainregistration.o


CMakeFiles/CoregTool.dir/brainio.o: CMakeFiles/CoregTool.dir/flags.make
CMakeFiles/CoregTool.dir/brainio.o: ../brainio.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mariano/BrainToolsCARM/QtCreator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/CoregTool.dir/brainio.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CoregTool.dir/brainio.o -c /home/mariano/BrainToolsCARM/brainio.cpp

CMakeFiles/CoregTool.dir/brainio.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CoregTool.dir/brainio.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mariano/BrainToolsCARM/brainio.cpp > CMakeFiles/CoregTool.dir/brainio.i

CMakeFiles/CoregTool.dir/brainio.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CoregTool.dir/brainio.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mariano/BrainToolsCARM/brainio.cpp -o CMakeFiles/CoregTool.dir/brainio.s

CMakeFiles/CoregTool.dir/brainio.o.requires:

.PHONY : CMakeFiles/CoregTool.dir/brainio.o.requires

CMakeFiles/CoregTool.dir/brainio.o.provides: CMakeFiles/CoregTool.dir/brainio.o.requires
	$(MAKE) -f CMakeFiles/CoregTool.dir/build.make CMakeFiles/CoregTool.dir/brainio.o.provides.build
.PHONY : CMakeFiles/CoregTool.dir/brainio.o.provides

CMakeFiles/CoregTool.dir/brainio.o.provides.build: CMakeFiles/CoregTool.dir/brainio.o


# Object files for target CoregTool
CoregTool_OBJECTS = \
"CMakeFiles/CoregTool.dir/coregmain.o" \
"CMakeFiles/CoregTool.dir/brainregistration.o" \
"CMakeFiles/CoregTool.dir/brainio.o"

# External object files for target CoregTool
CoregTool_EXTERNAL_OBJECTS =

CoregTool: CMakeFiles/CoregTool.dir/coregmain.o
CoregTool: CMakeFiles/CoregTool.dir/brainregistration.o
CoregTool: CMakeFiles/CoregTool.dir/brainio.o
CoregTool: CMakeFiles/CoregTool.dir/build.make
CoregTool: CMakeFiles/CoregTool.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mariano/BrainToolsCARM/QtCreator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable CoregTool"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CoregTool.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CoregTool.dir/build: CoregTool

.PHONY : CMakeFiles/CoregTool.dir/build

CMakeFiles/CoregTool.dir/requires: CMakeFiles/CoregTool.dir/coregmain.o.requires
CMakeFiles/CoregTool.dir/requires: CMakeFiles/CoregTool.dir/brainregistration.o.requires
CMakeFiles/CoregTool.dir/requires: CMakeFiles/CoregTool.dir/brainio.o.requires

.PHONY : CMakeFiles/CoregTool.dir/requires

CMakeFiles/CoregTool.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/CoregTool.dir/cmake_clean.cmake
.PHONY : CMakeFiles/CoregTool.dir/clean

CMakeFiles/CoregTool.dir/depend:
	cd /home/mariano/BrainToolsCARM/QtCreator && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mariano/BrainToolsCARM /home/mariano/BrainToolsCARM /home/mariano/BrainToolsCARM/QtCreator /home/mariano/BrainToolsCARM/QtCreator /home/mariano/BrainToolsCARM/QtCreator/CMakeFiles/CoregTool.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/CoregTool.dir/depend
