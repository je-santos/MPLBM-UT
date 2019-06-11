##########################################################################
## Makefile for Permeability and Shan-Chen 2 phase
##
## The present Makefile is a pure configuration file, in which
## you can select compilation options. Compilation dependencies
## are managed automatically through the Python library SConstruct.
##
## If you don't have Python, or if compilation doesn't work for other
## reasons, consult the Palabos user's guide for instructions on manual
## compilation.
##########################################################################

# USE: multiple arguments are separated by spaces.
#   For example: projectFiles = file1.cpp file2.cpp
#                optimFlags   = -O -finline-functions

# Leading directory of the Palabos source code
palabosRoot   = ../Palabos
# Name of source files in current directory to compile and link with Palabos
projectFiles = 2-phase_LBM/ShanChen.cpp 1-phase_LBM/permeability.cpp

# Set optimization flags on/off
optimize     = true
# Set debug mode and debug flags on/off
debug        = false
# Set profiling flags on/off
profile      = false
# Set MPI-parallel mode on/off (parallelism in cluster-like environment)
MPIparallel  = true
# Set SMP-parallel mode on/off (shared-memory parallelism)
SMPparallel  = false
# Decide whether to include calls to the POSIX API. On non-POSIX systems,
#   including Windows, this flag must be false, unless a POSIX environment is
#   emulated (such as with Cygwin).
usePOSIX     = true
# This flag must be defined true if you are using the external library
# cvmlcpp. But first, you need to download cvmlcpp and put it into the
# directory "externalLibraries" of the Palabos root.
useCVMLCPP   = false

# Path to external libraries (other than Palabos)
libraryPaths =
# Path to inlude directories (other than Palabos)
includePaths =
# Dynamic and static libraries (other than Palabos)
libraries    =

# Compiler to use without MPI parallelism
serialCXX    = g++
# Compiler to use with MPI parallelism
parallelCXX  = mpicxx
# General compiler flags (e.g. -Wall to turn on all warnings on g++)
compileFlags =
# General linker flags (don't put library includes into this flag)
linkFlags    =
# Compiler flags to use when optimization mode is on
optimFlags   = -O3
# Compiler flags to use when debug mode is on
debugFlags   = -g
# Compiler flags to use when profile mode is on
profileFlags = -pg


##########################################################################
# All code below this line is just about forwarding the options
# to SConstruct. It is recommended not to modify anything there.
##########################################################################

SCons     = $(palabosRoot)/scons/scons.py -f $(palabosRoot)/SConstruct

SConsArgs = palabosRoot=$(palabosRoot) \
            projectFiles="$(projectFiles)" \
            optimize=$(optimize) \
            debug=$(debug) \
            profile=$(profile) \
            MPIparallel=$(MPIparallel) \
            SMPparallel=$(SMPparallel) \
            usePOSIX=$(usePOSIX) \
	    useCVMLCPP=$(useCVMLCPP) \
            serialCXX=$(serialCXX) \
            parallelCXX=$(parallelCXX) \
            compileFlags="$(compileFlags)" \
            linkFlags="$(linkFlags)" \
            optimFlags="$(optimFlags)" \
            debugFlags="$(debugFlags)" \
	    profileFlags="$(profileFlags)" \
	    libraryPaths="$(libraryPaths)" \
	    includePaths="$(includePaths)" \
	    libraries="$(libraries)"

compile:
	python $(SCons) $(SConsArgs)

clean:
	python $(SCons) -c $(SConsArgs)
	/bin/rm -vf `find $(palabosRoot) -name '*~'`
