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
include CMakeFiles/eneq.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/eneq.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/eneq.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/eneq.dir/flags.make

CMakeFiles/eneq.dir/utility/eneq.cpp.o: CMakeFiles/eneq.dir/flags.make
CMakeFiles/eneq.dir/utility/eneq.cpp.o: ../utility/eneq.cpp
CMakeFiles/eneq.dir/utility/eneq.cpp.o: CMakeFiles/eneq.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmb/develop-git/semtex-xxt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/eneq.dir/utility/eneq.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/eneq.dir/utility/eneq.cpp.o -MF CMakeFiles/eneq.dir/utility/eneq.cpp.o.d -o CMakeFiles/eneq.dir/utility/eneq.cpp.o -c /Users/hmb/develop-git/semtex-xxt/utility/eneq.cpp

CMakeFiles/eneq.dir/utility/eneq.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/eneq.dir/utility/eneq.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hmb/develop-git/semtex-xxt/utility/eneq.cpp > CMakeFiles/eneq.dir/utility/eneq.cpp.i

CMakeFiles/eneq.dir/utility/eneq.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/eneq.dir/utility/eneq.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hmb/develop-git/semtex-xxt/utility/eneq.cpp -o CMakeFiles/eneq.dir/utility/eneq.cpp.s

# Object files for target eneq
eneq_OBJECTS = \
"CMakeFiles/eneq.dir/utility/eneq.cpp.o"

# External object files for target eneq
eneq_EXTERNAL_OBJECTS =

eneq: CMakeFiles/eneq.dir/utility/eneq.cpp.o
eneq: CMakeFiles/eneq.dir/build.make
eneq: src/libsrc.a
eneq: femlib/libfem.a
eneq: veclib/libvec.a
eneq: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/libblas.tbd
eneq: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/liblapack.tbd
eneq: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/libblas.tbd
eneq: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/liblapack.tbd
eneq: CMakeFiles/eneq.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/hmb/develop-git/semtex-xxt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable eneq"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/eneq.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/eneq.dir/build: eneq
.PHONY : CMakeFiles/eneq.dir/build

CMakeFiles/eneq.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/eneq.dir/cmake_clean.cmake
.PHONY : CMakeFiles/eneq.dir/clean

CMakeFiles/eneq.dir/depend:
	cd /Users/hmb/develop-git/semtex-xxt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/hmb/develop-git/semtex-xxt /Users/hmb/develop-git/semtex-xxt /Users/hmb/develop-git/semtex-xxt/build /Users/hmb/develop-git/semtex-xxt/build /Users/hmb/develop-git/semtex-xxt/build/CMakeFiles/eneq.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/eneq.dir/depend

