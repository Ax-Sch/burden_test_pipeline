# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build

# Include any dependencies generated for this target.
include tests/CMakeFiles/metaltest.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/metaltest.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/metaltest.dir/flags.make

tests/CMakeFiles/metaltest.dir/metaltest.cpp.o: tests/CMakeFiles/metaltest.dir/flags.make
tests/CMakeFiles/metaltest.dir/metaltest.cpp.o: ../tests/metaltest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/metaltest.dir/metaltest.cpp.o"
	cd /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/metaltest.dir/metaltest.cpp.o -c /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/tests/metaltest.cpp

tests/CMakeFiles/metaltest.dir/metaltest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/metaltest.dir/metaltest.cpp.i"
	cd /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/tests/metaltest.cpp > CMakeFiles/metaltest.dir/metaltest.cpp.i

tests/CMakeFiles/metaltest.dir/metaltest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/metaltest.dir/metaltest.cpp.s"
	cd /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/tests/metaltest.cpp -o CMakeFiles/metaltest.dir/metaltest.cpp.s

# Object files for target metaltest
metaltest_OBJECTS = \
"CMakeFiles/metaltest.dir/metaltest.cpp.o"

# External object files for target metaltest
metaltest_EXTERNAL_OBJECTS =

tests/metaltest: tests/CMakeFiles/metaltest.dir/metaltest.cpp.o
tests/metaltest: tests/CMakeFiles/metaltest.dir/build.make
tests/metaltest: libsrc/libgoncalo.a
tests/metaltest: /usr/lib/x86_64-linux-gnu/libz.so
tests/metaltest: tests/CMakeFiles/metaltest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable metaltest"
	cd /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/metaltest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/metaltest.dir/build: tests/metaltest

.PHONY : tests/CMakeFiles/metaltest.dir/build

tests/CMakeFiles/metaltest.dir/clean:
	cd /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/metaltest.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/metaltest.dir/clean

tests/CMakeFiles/metaltest.dir/depend:
	cd /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05 /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/tests /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/tests /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/tests/CMakeFiles/metaltest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/metaltest.dir/depend

