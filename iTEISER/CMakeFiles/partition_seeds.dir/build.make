# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/hani/Dropbox/cmake/iTEISER

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/hani/Dropbox/cmake/iTEISER

# Include any dependencies generated for this target.
include CMakeFiles/partition_seeds.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/partition_seeds.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/partition_seeds.dir/flags.make

CMakeFiles/partition_seeds.dir/partition_seeds.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/partition_seeds.c.o: partition_seeds.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/partition_seeds.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/partition_seeds.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/partition_seeds.c

CMakeFiles/partition_seeds.dir/partition_seeds.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/partition_seeds.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/partition_seeds.c > CMakeFiles/partition_seeds.dir/partition_seeds.c.i

CMakeFiles/partition_seeds.dir/partition_seeds.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/partition_seeds.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/partition_seeds.c -o CMakeFiles/partition_seeds.dir/partition_seeds.c.s

CMakeFiles/partition_seeds.dir/partition_seeds.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/partition_seeds.c.o.requires

CMakeFiles/partition_seeds.dir/partition_seeds.c.o.provides: CMakeFiles/partition_seeds.dir/partition_seeds.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/partition_seeds.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/partition_seeds.c.o.provides

CMakeFiles/partition_seeds.dir/partition_seeds.c.o.provides.build: CMakeFiles/partition_seeds.dir/partition_seeds.c.o

CMakeFiles/partition_seeds.dir/dataio.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/dataio.c.o: dataio.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/dataio.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/dataio.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/dataio.c

CMakeFiles/partition_seeds.dir/dataio.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/dataio.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/dataio.c > CMakeFiles/partition_seeds.dir/dataio.c.i

CMakeFiles/partition_seeds.dir/dataio.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/dataio.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/dataio.c -o CMakeFiles/partition_seeds.dir/dataio.c.s

CMakeFiles/partition_seeds.dir/dataio.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/dataio.c.o.requires

CMakeFiles/partition_seeds.dir/dataio.c.o.provides: CMakeFiles/partition_seeds.dir/dataio.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/dataio.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/dataio.c.o.provides

CMakeFiles/partition_seeds.dir/dataio.c.o.provides.build: CMakeFiles/partition_seeds.dir/dataio.c.o

CMakeFiles/partition_seeds.dir/hashtable.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/hashtable.c.o: hashtable.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/hashtable.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/hashtable.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/hashtable.c

CMakeFiles/partition_seeds.dir/hashtable.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/hashtable.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/hashtable.c > CMakeFiles/partition_seeds.dir/hashtable.c.i

CMakeFiles/partition_seeds.dir/hashtable.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/hashtable.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/hashtable.c -o CMakeFiles/partition_seeds.dir/hashtable.c.s

CMakeFiles/partition_seeds.dir/hashtable.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/hashtable.c.o.requires

CMakeFiles/partition_seeds.dir/hashtable.c.o.provides: CMakeFiles/partition_seeds.dir/hashtable.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/hashtable.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/hashtable.c.o.provides

CMakeFiles/partition_seeds.dir/hashtable.c.o.provides.build: CMakeFiles/partition_seeds.dir/hashtable.c.o

CMakeFiles/partition_seeds.dir/readFASTA.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/readFASTA.c.o: readFASTA.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/readFASTA.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/readFASTA.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/readFASTA.c

CMakeFiles/partition_seeds.dir/readFASTA.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/readFASTA.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/readFASTA.c > CMakeFiles/partition_seeds.dir/readFASTA.c.i

CMakeFiles/partition_seeds.dir/readFASTA.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/readFASTA.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/readFASTA.c -o CMakeFiles/partition_seeds.dir/readFASTA.c.s

CMakeFiles/partition_seeds.dir/readFASTA.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/readFASTA.c.o.requires

CMakeFiles/partition_seeds.dir/readFASTA.c.o.provides: CMakeFiles/partition_seeds.dir/readFASTA.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/readFASTA.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/readFASTA.c.o.provides

CMakeFiles/partition_seeds.dir/readFASTA.c.o.provides.build: CMakeFiles/partition_seeds.dir/readFASTA.c.o

CMakeFiles/partition_seeds.dir/read_write_motif.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/read_write_motif.c.o: read_write_motif.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/read_write_motif.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/read_write_motif.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/read_write_motif.c

CMakeFiles/partition_seeds.dir/read_write_motif.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/read_write_motif.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/read_write_motif.c > CMakeFiles/partition_seeds.dir/read_write_motif.c.i

CMakeFiles/partition_seeds.dir/read_write_motif.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/read_write_motif.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/read_write_motif.c -o CMakeFiles/partition_seeds.dir/read_write_motif.c.s

CMakeFiles/partition_seeds.dir/read_write_motif.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/read_write_motif.c.o.requires

CMakeFiles/partition_seeds.dir/read_write_motif.c.o.provides: CMakeFiles/partition_seeds.dir/read_write_motif.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/read_write_motif.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/read_write_motif.c.o.provides

CMakeFiles/partition_seeds.dir/read_write_motif.c.o.provides.build: CMakeFiles/partition_seeds.dir/read_write_motif.c.o

CMakeFiles/partition_seeds.dir/information.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/information.c.o: information.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/information.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/information.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/information.c

CMakeFiles/partition_seeds.dir/information.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/information.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/information.c > CMakeFiles/partition_seeds.dir/information.c.i

CMakeFiles/partition_seeds.dir/information.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/information.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/information.c -o CMakeFiles/partition_seeds.dir/information.c.s

CMakeFiles/partition_seeds.dir/information.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/information.c.o.requires

CMakeFiles/partition_seeds.dir/information.c.o.provides: CMakeFiles/partition_seeds.dir/information.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/information.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/information.c.o.provides

CMakeFiles/partition_seeds.dir/information.c.o.provides.build: CMakeFiles/partition_seeds.dir/information.c.o

CMakeFiles/partition_seeds.dir/mi_library.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/mi_library.c.o: mi_library.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/mi_library.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/mi_library.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/mi_library.c

CMakeFiles/partition_seeds.dir/mi_library.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/mi_library.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/mi_library.c > CMakeFiles/partition_seeds.dir/mi_library.c.i

CMakeFiles/partition_seeds.dir/mi_library.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/mi_library.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/mi_library.c -o CMakeFiles/partition_seeds.dir/mi_library.c.s

CMakeFiles/partition_seeds.dir/mi_library.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/mi_library.c.o.requires

CMakeFiles/partition_seeds.dir/mi_library.c.o.provides: CMakeFiles/partition_seeds.dir/mi_library.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/mi_library.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/mi_library.c.o.provides

CMakeFiles/partition_seeds.dir/mi_library.c.o.provides.build: CMakeFiles/partition_seeds.dir/mi_library.c.o

CMakeFiles/partition_seeds.dir/teiser_functions.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/teiser_functions.c.o: teiser_functions.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/teiser_functions.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/teiser_functions.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/teiser_functions.c

CMakeFiles/partition_seeds.dir/teiser_functions.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/teiser_functions.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/teiser_functions.c > CMakeFiles/partition_seeds.dir/teiser_functions.c.i

CMakeFiles/partition_seeds.dir/teiser_functions.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/teiser_functions.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/teiser_functions.c -o CMakeFiles/partition_seeds.dir/teiser_functions.c.s

CMakeFiles/partition_seeds.dir/teiser_functions.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/teiser_functions.c.o.requires

CMakeFiles/partition_seeds.dir/teiser_functions.c.o.provides: CMakeFiles/partition_seeds.dir/teiser_functions.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/teiser_functions.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/teiser_functions.c.o.provides

CMakeFiles/partition_seeds.dir/teiser_functions.c.o.provides.build: CMakeFiles/partition_seeds.dir/teiser_functions.c.o

CMakeFiles/partition_seeds.dir/statistics.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/statistics.c.o: statistics.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/statistics.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/statistics.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/statistics.c

CMakeFiles/partition_seeds.dir/statistics.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/statistics.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/statistics.c > CMakeFiles/partition_seeds.dir/statistics.c.i

CMakeFiles/partition_seeds.dir/statistics.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/statistics.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/statistics.c -o CMakeFiles/partition_seeds.dir/statistics.c.s

CMakeFiles/partition_seeds.dir/statistics.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/statistics.c.o.requires

CMakeFiles/partition_seeds.dir/statistics.c.o.provides: CMakeFiles/partition_seeds.dir/statistics.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/statistics.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/statistics.c.o.provides

CMakeFiles/partition_seeds.dir/statistics.c.o.provides.build: CMakeFiles/partition_seeds.dir/statistics.c.o

CMakeFiles/partition_seeds.dir/sequences.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/sequences.c.o: sequences.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_10)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/sequences.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/sequences.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/sequences.c

CMakeFiles/partition_seeds.dir/sequences.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/sequences.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/sequences.c > CMakeFiles/partition_seeds.dir/sequences.c.i

CMakeFiles/partition_seeds.dir/sequences.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/sequences.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/sequences.c -o CMakeFiles/partition_seeds.dir/sequences.c.s

CMakeFiles/partition_seeds.dir/sequences.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/sequences.c.o.requires

CMakeFiles/partition_seeds.dir/sequences.c.o.provides: CMakeFiles/partition_seeds.dir/sequences.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/sequences.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/sequences.c.o.provides

CMakeFiles/partition_seeds.dir/sequences.c.o.provides.build: CMakeFiles/partition_seeds.dir/sequences.c.o

CMakeFiles/partition_seeds.dir/matchmaker.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/matchmaker.c.o: matchmaker.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_11)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/matchmaker.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/matchmaker.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/matchmaker.c

CMakeFiles/partition_seeds.dir/matchmaker.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/matchmaker.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/matchmaker.c > CMakeFiles/partition_seeds.dir/matchmaker.c.i

CMakeFiles/partition_seeds.dir/matchmaker.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/matchmaker.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/matchmaker.c -o CMakeFiles/partition_seeds.dir/matchmaker.c.s

CMakeFiles/partition_seeds.dir/matchmaker.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/matchmaker.c.o.requires

CMakeFiles/partition_seeds.dir/matchmaker.c.o.provides: CMakeFiles/partition_seeds.dir/matchmaker.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/matchmaker.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/matchmaker.c.o.provides

CMakeFiles/partition_seeds.dir/matchmaker.c.o.provides.build: CMakeFiles/partition_seeds.dir/matchmaker.c.o

CMakeFiles/partition_seeds.dir/folding_energy.c.o: CMakeFiles/partition_seeds.dir/flags.make
CMakeFiles/partition_seeds.dir/folding_energy.c.o: folding_energy.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_12)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/partition_seeds.dir/folding_energy.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/partition_seeds.dir/folding_energy.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/folding_energy.c

CMakeFiles/partition_seeds.dir/folding_energy.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/partition_seeds.dir/folding_energy.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/folding_energy.c > CMakeFiles/partition_seeds.dir/folding_energy.c.i

CMakeFiles/partition_seeds.dir/folding_energy.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/partition_seeds.dir/folding_energy.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/folding_energy.c -o CMakeFiles/partition_seeds.dir/folding_energy.c.s

CMakeFiles/partition_seeds.dir/folding_energy.c.o.requires:
.PHONY : CMakeFiles/partition_seeds.dir/folding_energy.c.o.requires

CMakeFiles/partition_seeds.dir/folding_energy.c.o.provides: CMakeFiles/partition_seeds.dir/folding_energy.c.o.requires
	$(MAKE) -f CMakeFiles/partition_seeds.dir/build.make CMakeFiles/partition_seeds.dir/folding_energy.c.o.provides.build
.PHONY : CMakeFiles/partition_seeds.dir/folding_energy.c.o.provides

CMakeFiles/partition_seeds.dir/folding_energy.c.o.provides.build: CMakeFiles/partition_seeds.dir/folding_energy.c.o

# Object files for target partition_seeds
partition_seeds_OBJECTS = \
"CMakeFiles/partition_seeds.dir/partition_seeds.c.o" \
"CMakeFiles/partition_seeds.dir/dataio.c.o" \
"CMakeFiles/partition_seeds.dir/hashtable.c.o" \
"CMakeFiles/partition_seeds.dir/readFASTA.c.o" \
"CMakeFiles/partition_seeds.dir/read_write_motif.c.o" \
"CMakeFiles/partition_seeds.dir/information.c.o" \
"CMakeFiles/partition_seeds.dir/mi_library.c.o" \
"CMakeFiles/partition_seeds.dir/teiser_functions.c.o" \
"CMakeFiles/partition_seeds.dir/statistics.c.o" \
"CMakeFiles/partition_seeds.dir/sequences.c.o" \
"CMakeFiles/partition_seeds.dir/matchmaker.c.o" \
"CMakeFiles/partition_seeds.dir/folding_energy.c.o"

# External object files for target partition_seeds
partition_seeds_EXTERNAL_OBJECTS =

bin/partition_seeds: CMakeFiles/partition_seeds.dir/partition_seeds.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/dataio.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/hashtable.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/readFASTA.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/read_write_motif.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/information.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/mi_library.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/teiser_functions.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/statistics.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/sequences.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/matchmaker.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/folding_energy.c.o
bin/partition_seeds: CMakeFiles/partition_seeds.dir/build.make
bin/partition_seeds: CMakeFiles/partition_seeds.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable bin/partition_seeds"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/partition_seeds.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/partition_seeds.dir/build: bin/partition_seeds
.PHONY : CMakeFiles/partition_seeds.dir/build

CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/partition_seeds.c.o.requires
CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/dataio.c.o.requires
CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/hashtable.c.o.requires
CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/readFASTA.c.o.requires
CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/read_write_motif.c.o.requires
CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/information.c.o.requires
CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/mi_library.c.o.requires
CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/teiser_functions.c.o.requires
CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/statistics.c.o.requires
CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/sequences.c.o.requires
CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/matchmaker.c.o.requires
CMakeFiles/partition_seeds.dir/requires: CMakeFiles/partition_seeds.dir/folding_energy.c.o.requires
.PHONY : CMakeFiles/partition_seeds.dir/requires

CMakeFiles/partition_seeds.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/partition_seeds.dir/cmake_clean.cmake
.PHONY : CMakeFiles/partition_seeds.dir/clean

CMakeFiles/partition_seeds.dir/depend:
	cd /Users/hani/Dropbox/cmake/iTEISER && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles/partition_seeds.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/partition_seeds.dir/depend

