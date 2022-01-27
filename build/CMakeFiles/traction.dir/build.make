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
include CMakeFiles/traction.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/traction.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/traction.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/traction.dir/flags.make

CMakeFiles/traction.dir/utility/traction.cpp.o: CMakeFiles/traction.dir/flags.make
CMakeFiles/traction.dir/utility/traction.cpp.o: ../utility/traction.cpp
CMakeFiles/traction.dir/utility/traction.cpp.o: CMakeFiles/traction.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmb/develop-git/semtex-xxt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/traction.dir/utility/traction.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/traction.dir/utility/traction.cpp.o -MF CMakeFiles/traction.dir/utility/traction.cpp.o.d -o CMakeFiles/traction.dir/utility/traction.cpp.o -c /Users/hmb/develop-git/semtex-xxt/utility/traction.cpp

CMakeFiles/traction.dir/utility/traction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/traction.dir/utility/traction.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hmb/develop-git/semtex-xxt/utility/traction.cpp > CMakeFiles/traction.dir/utility/traction.cpp.i

CMakeFiles/traction.dir/utility/traction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/traction.dir/utility/traction.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hmb/develop-git/semtex-xxt/utility/traction.cpp -o CMakeFiles/traction.dir/utility/traction.cpp.s

# Object files for target traction
traction_OBJECTS = \
"CMakeFiles/traction.dir/utility/traction.cpp.o"

# External object files for target traction
traction_EXTERNAL_OBJECTS =

traction: CMakeFiles/traction.dir/utility/traction.cpp.o
traction: CMakeFiles/traction.dir/build.make
traction: src/libsrc.a
traction: femlib/libfem.a
traction: veclib/libvec.a
traction: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/libblas.tbd
traction: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/liblapack.tbd
traction: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/libblas.tbd
traction: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/liblapack.tbd
traction: CMakeFiles/traction.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/hmb/develop-git/semtex-xxt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable traction"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/traction.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/traction.dir/build: traction
.PHONY : CMakeFiles/traction.dir/build

CMakeFiles/traction.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/traction.dir/cmake_clean.cmake
.PHONY : CMakeFiles/traction.dir/clean

CMakeFiles/traction.dir/depend:
	cd /Users/hmb/develop-git/semtex-xxt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/hmb/develop-git/semtex-xxt /Users/hmb/develop-git/semtex-xxt /Users/hmb/develop-git/semtex-xxt/build /Users/hmb/develop-git/semtex-xxt/build /Users/hmb/develop-git/semtex-xxt/build/CMakeFiles/traction.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/traction.dir/depend

