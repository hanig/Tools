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
include CMakeFiles/print_motifs.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/print_motifs.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/print_motifs.dir/flags.make

CMakeFiles/print_motifs.dir/print_motifs.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/print_motifs.c.o: print_motifs.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/print_motifs.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/print_motifs.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/print_motifs.c

CMakeFiles/print_motifs.dir/print_motifs.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/print_motifs.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/print_motifs.c > CMakeFiles/print_motifs.dir/print_motifs.c.i

CMakeFiles/print_motifs.dir/print_motifs.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/print_motifs.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/print_motifs.c -o CMakeFiles/print_motifs.dir/print_motifs.c.s

CMakeFiles/print_motifs.dir/print_motifs.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/print_motifs.c.o.requires

CMakeFiles/print_motifs.dir/print_motifs.c.o.provides: CMakeFiles/print_motifs.dir/print_motifs.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/print_motifs.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/print_motifs.c.o.provides

CMakeFiles/print_motifs.dir/print_motifs.c.o.provides.build: CMakeFiles/print_motifs.dir/print_motifs.c.o

CMakeFiles/print_motifs.dir/dataio.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/dataio.c.o: dataio.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/dataio.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/dataio.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/dataio.c

CMakeFiles/print_motifs.dir/dataio.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/dataio.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/dataio.c > CMakeFiles/print_motifs.dir/dataio.c.i

CMakeFiles/print_motifs.dir/dataio.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/dataio.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/dataio.c -o CMakeFiles/print_motifs.dir/dataio.c.s

CMakeFiles/print_motifs.dir/dataio.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/dataio.c.o.requires

CMakeFiles/print_motifs.dir/dataio.c.o.provides: CMakeFiles/print_motifs.dir/dataio.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/dataio.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/dataio.c.o.provides

CMakeFiles/print_motifs.dir/dataio.c.o.provides.build: CMakeFiles/print_motifs.dir/dataio.c.o

CMakeFiles/print_motifs.dir/hashtable.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/hashtable.c.o: hashtable.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/hashtable.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/hashtable.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/hashtable.c

CMakeFiles/print_motifs.dir/hashtable.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/hashtable.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/hashtable.c > CMakeFiles/print_motifs.dir/hashtable.c.i

CMakeFiles/print_motifs.dir/hashtable.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/hashtable.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/hashtable.c -o CMakeFiles/print_motifs.dir/hashtable.c.s

CMakeFiles/print_motifs.dir/hashtable.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/hashtable.c.o.requires

CMakeFiles/print_motifs.dir/hashtable.c.o.provides: CMakeFiles/print_motifs.dir/hashtable.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/hashtable.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/hashtable.c.o.provides

CMakeFiles/print_motifs.dir/hashtable.c.o.provides.build: CMakeFiles/print_motifs.dir/hashtable.c.o

CMakeFiles/print_motifs.dir/readFASTA.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/readFASTA.c.o: readFASTA.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/readFASTA.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/readFASTA.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/readFASTA.c

CMakeFiles/print_motifs.dir/readFASTA.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/readFASTA.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/readFASTA.c > CMakeFiles/print_motifs.dir/readFASTA.c.i

CMakeFiles/print_motifs.dir/readFASTA.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/readFASTA.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/readFASTA.c -o CMakeFiles/print_motifs.dir/readFASTA.c.s

CMakeFiles/print_motifs.dir/readFASTA.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/readFASTA.c.o.requires

CMakeFiles/print_motifs.dir/readFASTA.c.o.provides: CMakeFiles/print_motifs.dir/readFASTA.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/readFASTA.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/readFASTA.c.o.provides

CMakeFiles/print_motifs.dir/readFASTA.c.o.provides.build: CMakeFiles/print_motifs.dir/readFASTA.c.o

CMakeFiles/print_motifs.dir/read_write_motif.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/read_write_motif.c.o: read_write_motif.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/read_write_motif.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/read_write_motif.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/read_write_motif.c

CMakeFiles/print_motifs.dir/read_write_motif.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/read_write_motif.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/read_write_motif.c > CMakeFiles/print_motifs.dir/read_write_motif.c.i

CMakeFiles/print_motifs.dir/read_write_motif.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/read_write_motif.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/read_write_motif.c -o CMakeFiles/print_motifs.dir/read_write_motif.c.s

CMakeFiles/print_motifs.dir/read_write_motif.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/read_write_motif.c.o.requires

CMakeFiles/print_motifs.dir/read_write_motif.c.o.provides: CMakeFiles/print_motifs.dir/read_write_motif.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/read_write_motif.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/read_write_motif.c.o.provides

CMakeFiles/print_motifs.dir/read_write_motif.c.o.provides.build: CMakeFiles/print_motifs.dir/read_write_motif.c.o

CMakeFiles/print_motifs.dir/information.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/information.c.o: information.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/information.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/information.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/information.c

CMakeFiles/print_motifs.dir/information.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/information.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/information.c > CMakeFiles/print_motifs.dir/information.c.i

CMakeFiles/print_motifs.dir/information.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/information.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/information.c -o CMakeFiles/print_motifs.dir/information.c.s

CMakeFiles/print_motifs.dir/information.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/information.c.o.requires

CMakeFiles/print_motifs.dir/information.c.o.provides: CMakeFiles/print_motifs.dir/information.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/information.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/information.c.o.provides

CMakeFiles/print_motifs.dir/information.c.o.provides.build: CMakeFiles/print_motifs.dir/information.c.o

CMakeFiles/print_motifs.dir/mi_library.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/mi_library.c.o: mi_library.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/mi_library.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/mi_library.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/mi_library.c

CMakeFiles/print_motifs.dir/mi_library.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/mi_library.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/mi_library.c > CMakeFiles/print_motifs.dir/mi_library.c.i

CMakeFiles/print_motifs.dir/mi_library.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/mi_library.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/mi_library.c -o CMakeFiles/print_motifs.dir/mi_library.c.s

CMakeFiles/print_motifs.dir/mi_library.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/mi_library.c.o.requires

CMakeFiles/print_motifs.dir/mi_library.c.o.provides: CMakeFiles/print_motifs.dir/mi_library.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/mi_library.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/mi_library.c.o.provides

CMakeFiles/print_motifs.dir/mi_library.c.o.provides.build: CMakeFiles/print_motifs.dir/mi_library.c.o

CMakeFiles/print_motifs.dir/teiser_functions.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/teiser_functions.c.o: teiser_functions.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/teiser_functions.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/teiser_functions.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/teiser_functions.c

CMakeFiles/print_motifs.dir/teiser_functions.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/teiser_functions.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/teiser_functions.c > CMakeFiles/print_motifs.dir/teiser_functions.c.i

CMakeFiles/print_motifs.dir/teiser_functions.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/teiser_functions.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/teiser_functions.c -o CMakeFiles/print_motifs.dir/teiser_functions.c.s

CMakeFiles/print_motifs.dir/teiser_functions.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/teiser_functions.c.o.requires

CMakeFiles/print_motifs.dir/teiser_functions.c.o.provides: CMakeFiles/print_motifs.dir/teiser_functions.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/teiser_functions.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/teiser_functions.c.o.provides

CMakeFiles/print_motifs.dir/teiser_functions.c.o.provides.build: CMakeFiles/print_motifs.dir/teiser_functions.c.o

CMakeFiles/print_motifs.dir/statistics.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/statistics.c.o: statistics.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/statistics.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/statistics.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/statistics.c

CMakeFiles/print_motifs.dir/statistics.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/statistics.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/statistics.c > CMakeFiles/print_motifs.dir/statistics.c.i

CMakeFiles/print_motifs.dir/statistics.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/statistics.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/statistics.c -o CMakeFiles/print_motifs.dir/statistics.c.s

CMakeFiles/print_motifs.dir/statistics.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/statistics.c.o.requires

CMakeFiles/print_motifs.dir/statistics.c.o.provides: CMakeFiles/print_motifs.dir/statistics.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/statistics.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/statistics.c.o.provides

CMakeFiles/print_motifs.dir/statistics.c.o.provides.build: CMakeFiles/print_motifs.dir/statistics.c.o

CMakeFiles/print_motifs.dir/sequences.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/sequences.c.o: sequences.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_10)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/sequences.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/sequences.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/sequences.c

CMakeFiles/print_motifs.dir/sequences.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/sequences.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/sequences.c > CMakeFiles/print_motifs.dir/sequences.c.i

CMakeFiles/print_motifs.dir/sequences.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/sequences.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/sequences.c -o CMakeFiles/print_motifs.dir/sequences.c.s

CMakeFiles/print_motifs.dir/sequences.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/sequences.c.o.requires

CMakeFiles/print_motifs.dir/sequences.c.o.provides: CMakeFiles/print_motifs.dir/sequences.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/sequences.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/sequences.c.o.provides

CMakeFiles/print_motifs.dir/sequences.c.o.provides.build: CMakeFiles/print_motifs.dir/sequences.c.o

CMakeFiles/print_motifs.dir/matchmaker.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/matchmaker.c.o: matchmaker.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_11)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/matchmaker.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/matchmaker.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/matchmaker.c

CMakeFiles/print_motifs.dir/matchmaker.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/matchmaker.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/matchmaker.c > CMakeFiles/print_motifs.dir/matchmaker.c.i

CMakeFiles/print_motifs.dir/matchmaker.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/matchmaker.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/matchmaker.c -o CMakeFiles/print_motifs.dir/matchmaker.c.s

CMakeFiles/print_motifs.dir/matchmaker.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/matchmaker.c.o.requires

CMakeFiles/print_motifs.dir/matchmaker.c.o.provides: CMakeFiles/print_motifs.dir/matchmaker.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/matchmaker.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/matchmaker.c.o.provides

CMakeFiles/print_motifs.dir/matchmaker.c.o.provides.build: CMakeFiles/print_motifs.dir/matchmaker.c.o

CMakeFiles/print_motifs.dir/folding_energy.c.o: CMakeFiles/print_motifs.dir/flags.make
CMakeFiles/print_motifs.dir/folding_energy.c.o: folding_energy.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_12)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/print_motifs.dir/folding_energy.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/print_motifs.dir/folding_energy.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/folding_energy.c

CMakeFiles/print_motifs.dir/folding_energy.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/print_motifs.dir/folding_energy.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/folding_energy.c > CMakeFiles/print_motifs.dir/folding_energy.c.i

CMakeFiles/print_motifs.dir/folding_energy.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/print_motifs.dir/folding_energy.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/folding_energy.c -o CMakeFiles/print_motifs.dir/folding_energy.c.s

CMakeFiles/print_motifs.dir/folding_energy.c.o.requires:
.PHONY : CMakeFiles/print_motifs.dir/folding_energy.c.o.requires

CMakeFiles/print_motifs.dir/folding_energy.c.o.provides: CMakeFiles/print_motifs.dir/folding_energy.c.o.requires
	$(MAKE) -f CMakeFiles/print_motifs.dir/build.make CMakeFiles/print_motifs.dir/folding_energy.c.o.provides.build
.PHONY : CMakeFiles/print_motifs.dir/folding_energy.c.o.provides

CMakeFiles/print_motifs.dir/folding_energy.c.o.provides.build: CMakeFiles/print_motifs.dir/folding_energy.c.o

# Object files for target print_motifs
print_motifs_OBJECTS = \
"CMakeFiles/print_motifs.dir/print_motifs.c.o" \
"CMakeFiles/print_motifs.dir/dataio.c.o" \
"CMakeFiles/print_motifs.dir/hashtable.c.o" \
"CMakeFiles/print_motifs.dir/readFASTA.c.o" \
"CMakeFiles/print_motifs.dir/read_write_motif.c.o" \
"CMakeFiles/print_motifs.dir/information.c.o" \
"CMakeFiles/print_motifs.dir/mi_library.c.o" \
"CMakeFiles/print_motifs.dir/teiser_functions.c.o" \
"CMakeFiles/print_motifs.dir/statistics.c.o" \
"CMakeFiles/print_motifs.dir/sequences.c.o" \
"CMakeFiles/print_motifs.dir/matchmaker.c.o" \
"CMakeFiles/print_motifs.dir/folding_energy.c.o"

# External object files for target print_motifs
print_motifs_EXTERNAL_OBJECTS =

bin/print_motifs: CMakeFiles/print_motifs.dir/print_motifs.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/dataio.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/hashtable.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/readFASTA.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/read_write_motif.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/information.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/mi_library.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/teiser_functions.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/statistics.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/sequences.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/matchmaker.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/folding_energy.c.o
bin/print_motifs: CMakeFiles/print_motifs.dir/build.make
bin/print_motifs: CMakeFiles/print_motifs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable bin/print_motifs"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/print_motifs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/print_motifs.dir/build: bin/print_motifs
.PHONY : CMakeFiles/print_motifs.dir/build

CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/print_motifs.c.o.requires
CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/dataio.c.o.requires
CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/hashtable.c.o.requires
CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/readFASTA.c.o.requires
CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/read_write_motif.c.o.requires
CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/information.c.o.requires
CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/mi_library.c.o.requires
CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/teiser_functions.c.o.requires
CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/statistics.c.o.requires
CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/sequences.c.o.requires
CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/matchmaker.c.o.requires
CMakeFiles/print_motifs.dir/requires: CMakeFiles/print_motifs.dir/folding_energy.c.o.requires
.PHONY : CMakeFiles/print_motifs.dir/requires

CMakeFiles/print_motifs.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/print_motifs.dir/cmake_clean.cmake
.PHONY : CMakeFiles/print_motifs.dir/clean

CMakeFiles/print_motifs.dir/depend:
	cd /Users/hani/Dropbox/cmake/iTEISER && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles/print_motifs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/print_motifs.dir/depend

