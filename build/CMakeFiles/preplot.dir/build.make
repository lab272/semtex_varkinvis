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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/hmb/develop-git/semtex-xxt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/hmb/develop-git/semtex-xxt/build

# Include any dependencies generated for this target.
include CMakeFiles/preplot.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/preplot.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/preplot.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/preplot.dir/flags.make

CMakeFiles/preplot.dir/utility/preplot.c.o: CMakeFiles/preplot.dir/flags.make
CMakeFiles/preplot.dir/utility/preplot.c.o: ../utility/preplot.c
CMakeFiles/preplot.dir/utility/preplot.c.o: CMakeFiles/preplot.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmb/develop-git/semtex-xxt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/preplot.dir/utility/preplot.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/preplot.dir/utility/preplot.c.o -MF CMakeFiles/preplot.dir/utility/preplot.c.o.d -o CMakeFiles/preplot.dir/utility/preplot.c.o -c /Users/hmb/develop-git/semtex-xxt/utility/preplot.c

CMakeFiles/preplot.dir/utility/preplot.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/preplot.dir/utility/preplot.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmb/develop-git/semtex-xxt/utility/preplot.c > CMakeFiles/preplot.dir/utility/preplot.c.i

CMakeFiles/preplot.dir/utility/preplot.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/preplot.dir/utility/preplot.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmb/develop-git/semtex-xxt/utility/preplot.c -o CMakeFiles/preplot.dir/utility/preplot.c.s

# Object files for target preplot
preplot_OBJECTS = \
"CMakeFiles/preplot.dir/utility/preplot.c.o"

# External object files for target preplot
preplot_EXTERNAL_OBJECTS =

preplot: CMakeFiles/preplot.dir/utility/preplot.c.o
preplot: CMakeFiles/preplot.dir/build.make
preplot: src/libsrc.a
preplot: femlib/libfem.a
preplot: veclib/libvec.a
preplot: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/libblas.tbd
preplot: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/liblapack.tbd
preplot: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/libblas.tbd
preplot: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/liblapack.tbd
preplot: CMakeFiles/preplot.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/hmb/develop-git/semtex-xxt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable preplot"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/preplot.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/preplot.dir/build: preplot
.PHONY : CMakeFiles/preplot.dir/build

CMakeFiles/preplot.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/preplot.dir/cmake_clean.cmake
.PHONY : CMakeFiles/preplot.dir/clean

CMakeFiles/preplot.dir/depend:
	cd /Users/hmb/develop-git/semtex-xxt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/hmb/develop-git/semtex-xxt /Users/hmb/develop-git/semtex-xxt /Users/hmb/develop-git/semtex-xxt/build /Users/hmb/develop-git/semtex-xxt/build /Users/hmb/develop-git/semtex-xxt/build/CMakeFiles/preplot.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/preplot.dir/depend

