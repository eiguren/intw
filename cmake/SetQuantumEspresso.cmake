#############################################################################
#
#############################################################################

message("-- Checking QE.")

# List of supported QE versions
set(SUPPORTED_QE
  "6.7+"
  )


# Find QE build dir
file(GLOB files_and_dirs RELATIVE ${QE_HOME} ${QE_HOME}/*)
set(dirs "")
foreach(file_or_dir ${files_and_dirs})
  if(IS_DIRECTORY ${QE_HOME}/${file_or_dir})
    list(APPEND dirs ${file_or_dir})
  endif()
endforeach()

find_file(QE_CMAKE CMakeCache.txt PATHS ${QE_HOME} PATH_SUFFIXES ${dirs} NO_DEFAULT_PATH)
if("${QE_CMAKE}" STREQUAL "QE_CMAKE-NOTFOUND")
  message("-- QE built with autotools.")
  set(QE_BUILD_DIR ${QE_HOME})
else()
  message("-- QE built with CMake.")
  get_filename_component(QE_BUILD_DIR ${QE_CMAKE} DIRECTORY)
endif()
message("-- Found QE_BUILD_DIR: ${QE_BUILD_DIR}")


# Some checks
message("-- Checking QE build.")
find_program(QE_PW_EXE pw.x PATHS ${QE_BUILD_DIR}/bin NO_DEFAULT_PATH)
if("${QE_PW_EXE}" STREQUAL "QE_PW_EXE-NOTFOUND")
  message(FATAL_ERROR "QE_PW_EXE NOT FOUND! QE is not compiled or QE_HOME is incorrect.")
endif()
find_program(QE_PH_EXE ph.x PATHS ${QE_BUILD_DIR}/bin NO_DEFAULT_PATH)
if("${QE_PH_EXE}" STREQUAL "QE_PH_EXE-NOTFOUND")
  message(FATAL_ERROR "QE_PH_EXE NOT FOUND! QE is not compiled or QE_HOME is incorrect.")
endif()


# Find QE version number
find_file(QE_VERSION_FILE qe_version.h PATHS ${QE_HOME}/include NO_DEFAULT_PATH)
if("${QE_VERSION_FILE}" STREQUAL "QE_VERSION_FILE-NOTFOUND")
  # 6.6 version have the version number in a different file
  find_file(QE_VERSION_FILE version.h PATHS ${QE_HOME}/include NO_DEFAULT_PATH)
endif()
if("${QE_VERSION_FILE}" STREQUAL "QE_VERSION_FILE-NOTFOUND")
  # Older versions have the version number in a different file
  find_file(QE_VERSION_FILE version.f90 PATHS ${QE_HOME}/Modules NO_DEFAULT_PATH)
endif()

if("${QE_VERSION_FILE}" STREQUAL "QE_VERSION_FILE-NOTFOUND")
  message(FATAL_ERROR "QE_VERSION_FILE NOT FOUND!")
endif()

file(STRINGS ${QE_VERSION_FILE} QE_VERSION_LINE REGEX "version_number")
string(REGEX MATCH "'([0-9]+)(\\.[0-9]+)?(\\.[0-9]+)?(MaX)?(\\+)?'" QE_VERSION ${QE_VERSION_LINE})
string(REPLACE "'" "" QE_VERSION "${QE_VERSION}")
message("-- Found QE_VERSION: ${QE_VERSION}")


# Check if the given QE version is supported by intw
list(FIND SUPPORTED_QE ${QE_VERSION} QE_SUPPORTED)
if("${QE_SUPPORTED}" EQUAL "-1")
  message("QE especified by QE_HOME is not suported!")
  message("Supported QE versions: ${SUPPORTED_QE}")
  message(FATAL_ERROR)
endif()

