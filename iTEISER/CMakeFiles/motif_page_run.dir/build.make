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
include CMakeFiles/motif_page_run.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/motif_page_run.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/motif_page_run.dir/flags.make

CMakeFiles/motif_page_run.dir/motif_page_run.c.o: CMakeFiles/motif_page_run.dir/flags.make
CMakeFiles/motif_page_run.dir/motif_page_run.c.o: motif_page_run.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/motif_page_run.dir/motif_page_run.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/motif_page_run.dir/motif_page_run.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/motif_page_run.c

CMakeFiles/motif_page_run.dir/motif_page_run.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/motif_page_run.dir/motif_page_run.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/motif_page_run.c > CMakeFiles/motif_page_run.dir/motif_page_run.c.i

CMakeFiles/motif_page_run.dir/motif_page_run.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/motif_page_run.dir/motif_page_run.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/motif_page_run.c -o CMakeFiles/motif_page_run.dir/motif_page_run.c.s

CMakeFiles/motif_page_run.dir/motif_page_run.c.o.requires:
.PHONY : CMakeFiles/motif_page_run.dir/motif_page_run.c.o.requires

CMakeFiles/motif_page_run.dir/motif_page_run.c.o.provides: CMakeFiles/motif_page_run.dir/motif_page_run.c.o.requires
	$(MAKE) -f CMakeFiles/motif_page_run.dir/build.make CMakeFiles/motif_page_run.dir/motif_page_run.c.o.provides.build
.PHONY : CMakeFiles/motif_page_run.dir/motif_page_run.c.o.provides

CMakeFiles/motif_page_run.dir/motif_page_run.c.o.provides.build: CMakeFiles/motif_page_run.dir/motif_page_run.c.o

CMakeFiles/motif_page_run.dir/dataio.c.o: CMakeFiles/motif_page_run.dir/flags.make
CMakeFiles/motif_page_run.dir/dataio.c.o: dataio.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/motif_page_run.dir/dataio.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/motif_page_run.dir/dataio.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/dataio.c

CMakeFiles/motif_page_run.dir/dataio.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/motif_page_run.dir/dataio.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/dataio.c > CMakeFiles/motif_page_run.dir/dataio.c.i

CMakeFiles/motif_page_run.dir/dataio.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/motif_page_run.dir/dataio.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/dataio.c -o CMakeFiles/motif_page_run.dir/dataio.c.s

CMakeFiles/motif_page_run.dir/dataio.c.o.requires:
.PHONY : CMakeFiles/motif_page_run.dir/dataio.c.o.requires

CMakeFiles/motif_page_run.dir/dataio.c.o.provides: CMakeFiles/motif_page_run.dir/dataio.c.o.requires
	$(MAKE) -f CMakeFiles/motif_page_run.dir/build.make CMakeFiles/motif_page_run.dir/dataio.c.o.provides.build
.PHONY : CMakeFiles/motif_page_run.dir/dataio.c.o.provides

CMakeFiles/motif_page_run.dir/dataio.c.o.provides.build: CMakeFiles/motif_page_run.dir/dataio.c.o

CMakeFiles/motif_page_run.dir/hashtable.c.o: CMakeFiles/motif_page_run.dir/flags.make
CMakeFiles/motif_page_run.dir/hashtable.c.o: hashtable.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/motif_page_run.dir/hashtable.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/motif_page_run.dir/hashtable.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/hashtable.c

CMakeFiles/motif_page_run.dir/hashtable.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/motif_page_run.dir/hashtable.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/hashtable.c > CMakeFiles/motif_page_run.dir/hashtable.c.i

CMakeFiles/motif_page_run.dir/hashtable.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/motif_page_run.dir/hashtable.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/hashtable.c -o CMakeFiles/motif_page_run.dir/hashtable.c.s

CMakeFiles/motif_page_run.dir/hashtable.c.o.requires:
.PHONY : CMakeFiles/motif_page_run.dir/hashtable.c.o.requires

CMakeFiles/motif_page_run.dir/hashtable.c.o.provides: CMakeFiles/motif_page_run.dir/hashtable.c.o.requires
	$(MAKE) -f CMakeFiles/motif_page_run.dir/build.make CMakeFiles/motif_page_run.dir/hashtable.c.o.provides.build
.PHONY : CMakeFiles/motif_page_run.dir/hashtable.c.o.provides

CMakeFiles/motif_page_run.dir/hashtable.c.o.provides.build: CMakeFiles/motif_page_run.dir/hashtable.c.o

CMakeFiles/motif_page_run.dir/readFASTA.c.o: CMakeFiles/motif_page_run.dir/flags.make
CMakeFiles/motif_page_run.dir/readFASTA.c.o: readFASTA.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/motif_page_run.dir/readFASTA.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/motif_page_run.dir/readFASTA.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/readFASTA.c

CMakeFiles/motif_page_run.dir/readFASTA.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/motif_page_run.dir/readFASTA.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/readFASTA.c > CMakeFiles/motif_page_run.dir/readFASTA.c.i

CMakeFiles/motif_page_run.dir/readFASTA.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/motif_page_run.dir/readFASTA.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/readFASTA.c -o CMakeFiles/motif_page_run.dir/readFASTA.c.s

CMakeFiles/motif_page_run.dir/readFASTA.c.o.requires:
.PHONY : CMakeFiles/motif_page_run.dir/readFASTA.c.o.requires

CMakeFiles/motif_page_run.dir/readFASTA.c.o.provides: CMakeFiles/motif_page_run.dir/readFASTA.c.o.requires
	$(MAKE) -f CMakeFiles/motif_page_run.dir/build.make CMakeFiles/motif_page_run.dir/readFASTA.c.o.provides.build
.PHONY : CMakeFiles/motif_page_run.dir/readFASTA.c.o.provides

CMakeFiles/motif_page_run.dir/readFASTA.c.o.provides.build: CMakeFiles/motif_page_run.dir/readFASTA.c.o

CMakeFiles/motif_page_run.dir/read_write_motif.c.o: CMakeFiles/motif_page_run.dir/flags.make
CMakeFiles/motif_page_run.dir/read_write_motif.c.o: read_write_motif.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/motif_page_run.dir/read_write_motif.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/motif_page_run.dir/read_write_motif.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/read_write_motif.c

CMakeFiles/motif_page_run.dir/read_write_motif.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/motif_page_run.dir/read_write_motif.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/read_write_motif.c > CMakeFiles/motif_page_run.dir/read_write_motif.c.i

CMakeFiles/motif_page_run.dir/read_write_motif.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/motif_page_run.dir/read_write_motif.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/read_write_motif.c -o CMakeFiles/motif_page_run.dir/read_write_motif.c.s

CMakeFiles/motif_page_run.dir/read_write_motif.c.o.requires:
.PHONY : CMakeFiles/motif_page_run.dir/read_write_motif.c.o.requires

CMakeFiles/motif_page_run.dir/read_write_motif.c.o.provides: CMakeFiles/motif_page_run.dir/read_write_motif.c.o.requires
	$(MAKE) -f CMakeFiles/motif_page_run.dir/build.make CMakeFiles/motif_page_run.dir/read_write_motif.c.o.provides.build
.PHONY : CMakeFiles/motif_page_run.dir/read_write_motif.c.o.provides

CMakeFiles/motif_page_run.dir/read_write_motif.c.o.provides.build: CMakeFiles/motif_page_run.dir/read_write_motif.c.o

CMakeFiles/motif_page_run.dir/information.c.o: CMakeFiles/motif_page_run.dir/flags.make
CMakeFiles/motif_page_run.dir/information.c.o: information.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/motif_page_run.dir/information.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/motif_page_run.dir/information.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/information.c

CMakeFiles/motif_page_run.dir/information.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/motif_page_run.dir/information.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/information.c > CMakeFiles/motif_page_run.dir/information.c.i

CMakeFiles/motif_page_run.dir/information.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/motif_page_run.dir/information.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/information.c -o CMakeFiles/motif_page_run.dir/information.c.s

CMakeFiles/motif_page_run.dir/information.c.o.requires:
.PHONY : CMakeFiles/motif_page_run.dir/information.c.o.requires

CMakeFiles/motif_page_run.dir/information.c.o.provides: CMakeFiles/motif_page_run.dir/information.c.o.requires
	$(MAKE) -f CMakeFiles/motif_page_run.dir/build.make CMakeFiles/motif_page_run.dir/information.c.o.provides.build
.PHONY : CMakeFiles/motif_page_run.dir/information.c.o.provides

CMakeFiles/motif_page_run.dir/information.c.o.provides.build: CMakeFiles/motif_page_run.dir/information.c.o

CMakeFiles/motif_page_run.dir/mi_library.c.o: CMakeFiles/motif_page_run.dir/flags.make
CMakeFiles/motif_page_run.dir/mi_library.c.o: mi_library.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/motif_page_run.dir/mi_library.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/motif_page_run.dir/mi_library.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/mi_library.c

CMakeFiles/motif_page_run.dir/mi_library.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/motif_page_run.dir/mi_library.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/mi_library.c > CMakeFiles/motif_page_run.dir/mi_library.c.i

CMakeFiles/motif_page_run.dir/mi_library.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/motif_page_run.dir/mi_library.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/mi_library.c -o CMakeFiles/motif_page_run.dir/mi_library.c.s

CMakeFiles/motif_page_run.dir/mi_library.c.o.requires:
.PHONY : CMakeFiles/motif_page_run.dir/mi_library.c.o.requires

CMakeFiles/motif_page_run.dir/mi_library.c.o.provides: CMakeFiles/motif_page_run.dir/mi_library.c.o.requires
	$(MAKE) -f CMakeFiles/motif_page_run.dir/build.make CMakeFiles/motif_page_run.dir/mi_library.c.o.provides.build
.PHONY : CMakeFiles/motif_page_run.dir/mi_library.c.o.provides

CMakeFiles/motif_page_run.dir/mi_library.c.o.provides.build: CMakeFiles/motif_page_run.dir/mi_library.c.o

CMakeFiles/motif_page_run.dir/statistics.c.o: CMakeFiles/motif_page_run.dir/flags.make
CMakeFiles/motif_page_run.dir/statistics.c.o: statistics.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/motif_page_run.dir/statistics.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/motif_page_run.dir/statistics.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/statistics.c

CMakeFiles/motif_page_run.dir/statistics.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/motif_page_run.dir/statistics.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/statistics.c > CMakeFiles/motif_page_run.dir/statistics.c.i

CMakeFiles/motif_page_run.dir/statistics.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/motif_page_run.dir/statistics.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/statistics.c -o CMakeFiles/motif_page_run.dir/statistics.c.s

CMakeFiles/motif_page_run.dir/statistics.c.o.requires:
.PHONY : CMakeFiles/motif_page_run.dir/statistics.c.o.requires

CMakeFiles/motif_page_run.dir/statistics.c.o.provides: CMakeFiles/motif_page_run.dir/statistics.c.o.requires
	$(MAKE) -f CMakeFiles/motif_page_run.dir/build.make CMakeFiles/motif_page_run.dir/statistics.c.o.provides.build
.PHONY : CMakeFiles/motif_page_run.dir/statistics.c.o.provides

CMakeFiles/motif_page_run.dir/statistics.c.o.provides.build: CMakeFiles/motif_page_run.dir/statistics.c.o

CMakeFiles/motif_page_run.dir/sequences.c.o: CMakeFiles/motif_page_run.dir/flags.make
CMakeFiles/motif_page_run.dir/sequences.c.o: sequences.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/motif_page_run.dir/sequences.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/motif_page_run.dir/sequences.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/sequences.c

CMakeFiles/motif_page_run.dir/sequences.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/motif_page_run.dir/sequences.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/sequences.c > CMakeFiles/motif_page_run.dir/sequences.c.i

CMakeFiles/motif_page_run.dir/sequences.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/motif_page_run.dir/sequences.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/sequences.c -o CMakeFiles/motif_page_run.dir/sequences.c.s

CMakeFiles/motif_page_run.dir/sequences.c.o.requires:
.PHONY : CMakeFiles/motif_page_run.dir/sequences.c.o.requires

CMakeFiles/motif_page_run.dir/sequences.c.o.provides: CMakeFiles/motif_page_run.dir/sequences.c.o.requires
	$(MAKE) -f CMakeFiles/motif_page_run.dir/build.make CMakeFiles/motif_page_run.dir/sequences.c.o.provides.build
.PHONY : CMakeFiles/motif_page_run.dir/sequences.c.o.provides

CMakeFiles/motif_page_run.dir/sequences.c.o.provides.build: CMakeFiles/motif_page_run.dir/sequences.c.o

CMakeFiles/motif_page_run.dir/matchmaker.c.o: CMakeFiles/motif_page_run.dir/flags.make
CMakeFiles/motif_page_run.dir/matchmaker.c.o: matchmaker.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_10)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/motif_page_run.dir/matchmaker.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/motif_page_run.dir/matchmaker.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/matchmaker.c

CMakeFiles/motif_page_run.dir/matchmaker.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/motif_page_run.dir/matchmaker.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/matchmaker.c > CMakeFiles/motif_page_run.dir/matchmaker.c.i

CMakeFiles/motif_page_run.dir/matchmaker.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/motif_page_run.dir/matchmaker.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/matchmaker.c -o CMakeFiles/motif_page_run.dir/matchmaker.c.s

CMakeFiles/motif_page_run.dir/matchmaker.c.o.requires:
.PHONY : CMakeFiles/motif_page_run.dir/matchmaker.c.o.requires

CMakeFiles/motif_page_run.dir/matchmaker.c.o.provides: CMakeFiles/motif_page_run.dir/matchmaker.c.o.requires
	$(MAKE) -f CMakeFiles/motif_page_run.dir/build.make CMakeFiles/motif_page_run.dir/matchmaker.c.o.provides.build
.PHONY : CMakeFiles/motif_page_run.dir/matchmaker.c.o.provides

CMakeFiles/motif_page_run.dir/matchmaker.c.o.provides.build: CMakeFiles/motif_page_run.dir/matchmaker.c.o

CMakeFiles/motif_page_run.dir/folding_energy.c.o: CMakeFiles/motif_page_run.dir/flags.make
CMakeFiles/motif_page_run.dir/folding_energy.c.o: folding_energy.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles $(CMAKE_PROGRESS_11)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/motif_page_run.dir/folding_energy.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/motif_page_run.dir/folding_energy.c.o   -c /Users/hani/Dropbox/cmake/iTEISER/folding_energy.c

CMakeFiles/motif_page_run.dir/folding_energy.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/motif_page_run.dir/folding_energy.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/hani/Dropbox/cmake/iTEISER/folding_energy.c > CMakeFiles/motif_page_run.dir/folding_energy.c.i

CMakeFiles/motif_page_run.dir/folding_energy.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/motif_page_run.dir/folding_energy.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/hani/Dropbox/cmake/iTEISER/folding_energy.c -o CMakeFiles/motif_page_run.dir/folding_energy.c.s

CMakeFiles/motif_page_run.dir/folding_energy.c.o.requires:
.PHONY : CMakeFiles/motif_page_run.dir/folding_energy.c.o.requires

CMakeFiles/motif_page_run.dir/folding_energy.c.o.provides: CMakeFiles/motif_page_run.dir/folding_energy.c.o.requires
	$(MAKE) -f CMakeFiles/motif_page_run.dir/build.make CMakeFiles/motif_page_run.dir/folding_energy.c.o.provides.build
.PHONY : CMakeFiles/motif_page_run.dir/folding_energy.c.o.provides

CMakeFiles/motif_page_run.dir/folding_energy.c.o.provides.build: CMakeFiles/motif_page_run.dir/folding_energy.c.o

# Object files for target motif_page_run
motif_page_run_OBJECTS = \
"CMakeFiles/motif_page_run.dir/motif_page_run.c.o" \
"CMakeFiles/motif_page_run.dir/dataio.c.o" \
"CMakeFiles/motif_page_run.dir/hashtable.c.o" \
"CMakeFiles/motif_page_run.dir/readFASTA.c.o" \
"CMakeFiles/motif_page_run.dir/read_write_motif.c.o" \
"CMakeFiles/motif_page_run.dir/information.c.o" \
"CMakeFiles/motif_page_run.dir/mi_library.c.o" \
"CMakeFiles/motif_page_run.dir/statistics.c.o" \
"CMakeFiles/motif_page_run.dir/sequences.c.o" \
"CMakeFiles/motif_page_run.dir/matchmaker.c.o" \
"CMakeFiles/motif_page_run.dir/folding_energy.c.o"

# External object files for target motif_page_run
motif_page_run_EXTERNAL_OBJECTS =

bin/motif_page_run: CMakeFiles/motif_page_run.dir/motif_page_run.c.o
bin/motif_page_run: CMakeFiles/motif_page_run.dir/dataio.c.o
bin/motif_page_run: CMakeFiles/motif_page_run.dir/hashtable.c.o
bin/motif_page_run: CMakeFiles/motif_page_run.dir/readFASTA.c.o
bin/motif_page_run: CMakeFiles/motif_page_run.dir/read_write_motif.c.o
bin/motif_page_run: CMakeFiles/motif_page_run.dir/information.c.o
bin/motif_page_run: CMakeFiles/motif_page_run.dir/mi_library.c.o
bin/motif_page_run: CMakeFiles/motif_page_run.dir/statistics.c.o
bin/motif_page_run: CMakeFiles/motif_page_run.dir/sequences.c.o
bin/motif_page_run: CMakeFiles/motif_page_run.dir/matchmaker.c.o
bin/motif_page_run: CMakeFiles/motif_page_run.dir/folding_energy.c.o
bin/motif_page_run: CMakeFiles/motif_page_run.dir/build.make
bin/motif_page_run: CMakeFiles/motif_page_run.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable bin/motif_page_run"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/motif_page_run.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/motif_page_run.dir/build: bin/motif_page_run
.PHONY : CMakeFiles/motif_page_run.dir/build

CMakeFiles/motif_page_run.dir/requires: CMakeFiles/motif_page_run.dir/motif_page_run.c.o.requires
CMakeFiles/motif_page_run.dir/requires: CMakeFiles/motif_page_run.dir/dataio.c.o.requires
CMakeFiles/motif_page_run.dir/requires: CMakeFiles/motif_page_run.dir/hashtable.c.o.requires
CMakeFiles/motif_page_run.dir/requires: CMakeFiles/motif_page_run.dir/readFASTA.c.o.requires
CMakeFiles/motif_page_run.dir/requires: CMakeFiles/motif_page_run.dir/read_write_motif.c.o.requires
CMakeFiles/motif_page_run.dir/requires: CMakeFiles/motif_page_run.dir/information.c.o.requires
CMakeFiles/motif_page_run.dir/requires: CMakeFiles/motif_page_run.dir/mi_library.c.o.requires
CMakeFiles/motif_page_run.dir/requires: CMakeFiles/motif_page_run.dir/statistics.c.o.requires
CMakeFiles/motif_page_run.dir/requires: CMakeFiles/motif_page_run.dir/sequences.c.o.requires
CMakeFiles/motif_page_run.dir/requires: CMakeFiles/motif_page_run.dir/matchmaker.c.o.requires
CMakeFiles/motif_page_run.dir/requires: CMakeFiles/motif_page_run.dir/folding_energy.c.o.requires
.PHONY : CMakeFiles/motif_page_run.dir/requires

CMakeFiles/motif_page_run.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/motif_page_run.dir/cmake_clean.cmake
.PHONY : CMakeFiles/motif_page_run.dir/clean

CMakeFiles/motif_page_run.dir/depend:
	cd /Users/hani/Dropbox/cmake/iTEISER && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER /Users/hani/Dropbox/cmake/iTEISER/CMakeFiles/motif_page_run.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/motif_page_run.dir/depend

