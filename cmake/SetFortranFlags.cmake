######################################################
# Determine and set the Fortran compiler flags we want
######################################################

####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
INCLUDE(SetCompileFlag)

# Make sure the build type is uppercase
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

IF(BT STREQUAL "RELEASE")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or RELWITHDEBINFO."
      FORCE)
ELSEIF(BT STREQUAL "DEBUG")
    SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or RELWITHDEBINFO."
      FORCE)
ELSEIF(BT STREQUAL "RELWITHDEBINFO")
    SET (CMAKE_BUILD_TYPE RELWITHDEBINFO CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or RELWITHDEBINFO."
      FORCE)
ELSEIF(NOT BT)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Specifies the build type on single-configuration generators, options are Debug, Release, RelWithDebInfo and MinSizeRel."
      FORCE)
    MESSAGE(STATUS "CMAKE_BUILD_TYPE not given, defaulting to Release")
ELSE()
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE not valid, choices are Debug, Release, RelWithDebInfo and MinSizeRel")
ENDIF(BT STREQUAL "RELEASE")

#########################################################
# If the compiler flags have already been set, return now
#########################################################

IF(CMAKE_Fortran_FLAGS_DEBUG AND CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_MINSIZEREL AND CMAKE_Fortran_FLAGS_RELWITHDEBINFO AND SET_FORTRAN_FLAGS)
    RETURN ()
ENDIF(CMAKE_Fortran_FLAGS_DEBUG AND CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_MINSIZEREL AND CMAKE_Fortran_FLAGS_RELWITHDEBINFO AND SET_FORTRAN_FLAGS)

SET(SET_FORTRAN_FLAGS ON CACHE BOOL "If set, Fortran flags are already configured." FORCE)
message(STATUS "Setting Fortran flags")
########################################################################
# Determine the appropriate flags for this compiler for each build type.
# For each option type, a list of possible flags is given that work
# for various compilers.  The first flag that works is chosen.
# If none of the flags work, nothing is added (unless the REQUIRED
# flag is given in the call).  This way unknown compiles are supported.
#######################################################################

#####################
### GENERAL FLAGS ###
#####################

# Enable pre-processing
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran
                 REQUIRED
                 "-fpp" # Intel
                 "-cpp" # GNU
                )

#TODO: This is strongly discouraged and cfftnd.f90 should be fixed
#      to avoid the use of this flag
# Allow argument mismatch
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran
                 "-fallow-argument-mismatch" # GNU
                )

# Allow any line length in free format
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran
                 "-ffree-line-length-none" # GNU
                )

# Puts automatic arrays and arrays created for temporary computations on the heap instead of the stack
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran
                 "-heap-arrays" # Intel
                )

# Set record length units to bytes
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran
                 "-assume byterecl" # Intel
                )

###################
### DEBUG FLAGS ###
###################

# NOTE: debugging symbols (-g) are already on by default

# Traceback
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" Fortran
                 "-traceback"   # Intel
                 "-fbacktrace"  # GNU
                )

# Disable optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" Fortran
                 "-O0" # All compilers
                )

# Turn on all warnings
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" Fortran
                 "-warn all,nointerfaces" # Intel
                 "-Wall"                  # GNU
                )

# Check array bounds
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" Fortran
                 "-check all"     # Intel
                 "-fcheck=all"    # GNU
                )

# Generate complete debugging information
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" Fortran
                 "-debug extended" # Intel
                )

# Disable floating-point exception
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" Fortran
                 "-fpe0"              # Intel
                 "-mno-fp-exceptions" # GNU
                )

# Improves floating-point precision and consistency
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" Fortran
                 "-mp1" # Intel
                )

#####################
### RELEASE FLAGS ###
#####################

# NOTE: agressive optimizations (-O3) are already turned on by default

# There is some bug where -march=native doesn't work on Mac
IF(APPLE)
    SET(GNUNATIVE "-mtune=native")
ELSE()
    SET(GNUNATIVE "-march=native")
ENDIF()
# Optimize for the host's architecture
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" Fortran
                 "-xHost"     # Intel
                 ${GNUNATIVE} # GNU
                )

# Unroll loops
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" Fortran
                 "-unroll"        # Intel
                 "-funroll-loops" # GNU
                )

# Inline functions
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" Fortran
                 "-inline"            # Intel
                 "-finline-functions" # GNU
                )

# Interprocedural (link-time) optimizations
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" Fortran
#                 "-ipo"  # Intel
#                 "-flto" # GNU
#                )

# Single-file optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" Fortran
                 "-ip"  # Intel
                )

# Vectorized code report
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" Fortran
#                 "-vec-report0"  # Intel
#                )

############################
### RELWITHDEBINFO FLAGS ###
############################

# NOTE: optimizations (-O2) and debugging symbols (-g) are already on by default

# Traceback
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" Fortran
                 "-traceback"   # Intel
                 "-fbacktrace"  # GNU
                )

########################
### MINSIZEREL FLAGS ###
########################

# NOTE: optimizations that do not increase code size (-Os) are already turned on by default
