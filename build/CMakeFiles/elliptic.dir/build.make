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
include CMakeFiles/elliptic.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/elliptic.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/elliptic.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/elliptic.dir/flags.make

CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.o: CMakeFiles/elliptic.dir/flags.make
CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.o: ../elliptic/helmholtz.cpp
CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.o: CMakeFiles/elliptic.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmb/develop-git/semtex-xxt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.o -MF CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.o.d -o CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.o -c /Users/hmb/develop-git/semtex-xxt/elliptic/helmholtz.cpp

CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hmb/develop-git/semtex-xxt/elliptic/helmholtz.cpp > CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.i

CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hmb/develop-git/semtex-xxt/elliptic/helmholtz.cpp -o CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.s

CMakeFiles/elliptic.dir/elliptic/drive.cpp.o: CMakeFiles/elliptic.dir/flags.make
CMakeFiles/elliptic.dir/elliptic/drive.cpp.o: ../elliptic/drive.cpp
CMakeFiles/elliptic.dir/elliptic/drive.cpp.o: CMakeFiles/elliptic.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmb/develop-git/semtex-xxt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/elliptic.dir/elliptic/drive.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/elliptic.dir/elliptic/drive.cpp.o -MF CMakeFiles/elliptic.dir/elliptic/drive.cpp.o.d -o CMakeFiles/elliptic.dir/elliptic/drive.cpp.o -c /Users/hmb/develop-git/semtex-xxt/elliptic/drive.cpp

CMakeFiles/elliptic.dir/elliptic/drive.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/elliptic.dir/elliptic/drive.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hmb/develop-git/semtex-xxt/elliptic/drive.cpp > CMakeFiles/elliptic.dir/elliptic/drive.cpp.i

CMakeFiles/elliptic.dir/elliptic/drive.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/elliptic.dir/elliptic/drive.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hmb/develop-git/semtex-xxt/elliptic/drive.cpp -o CMakeFiles/elliptic.dir/elliptic/drive.cpp.s

# Object files for target elliptic
elliptic_OBJECTS = \
"CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.o" \
"CMakeFiles/elliptic.dir/elliptic/drive.cpp.o"

# External object files for target elliptic
elliptic_EXTERNAL_OBJECTS =

elliptic: CMakeFiles/elliptic.dir/elliptic/helmholtz.cpp.o
elliptic: CMakeFiles/elliptic.dir/elliptic/drive.cpp.o
elliptic: CMakeFiles/elliptic.dir/build.make
elliptic: src/libsrc.a
elliptic: femlib/libfem.a
elliptic: veclib/libvec.a
elliptic: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/libblas.tbd
elliptic: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/liblapack.tbd
elliptic: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/libblas.tbd
elliptic: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/usr/lib/liblapack.tbd
elliptic: CMakeFiles/elliptic.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/hmb/develop-git/semtex-xxt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable elliptic"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/elliptic.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/elliptic.dir/build: elliptic
.PHONY : CMakeFiles/elliptic.dir/build

CMakeFiles/elliptic.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/elliptic.dir/cmake_clean.cmake
.PHONY : CMakeFiles/elliptic.dir/clean

CMakeFiles/elliptic.dir/depend:
	cd /Users/hmb/develop-git/semtex-xxt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/hmb/develop-git/semtex-xxt /Users/hmb/develop-git/semtex-xxt /Users/hmb/develop-git/semtex-xxt/build /Users/hmb/develop-git/semtex-xxt/build /Users/hmb/develop-git/semtex-xxt/build/CMakeFiles/elliptic.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/elliptic.dir/depend

