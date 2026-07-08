cmake_minimum_required(VERSION 3.16)
project(wtgivlocal_opcode C)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Find csdl.h
find_path(CSOUND_INCLUDE_DIR
  NAMES csdl.h
  PATHS
    /Library/Frameworks/CsoundLib64.framework/Headers
    /Library/Frameworks/CsoundLib.framework/Headers
    /usr/local/include
    /opt/homebrew/include
  PATH_SUFFIXES csound
)

if(NOT CSOUND_INCLUDE_DIR)
  message(FATAL_ERROR
    "Could not find csdl.h.\n"
    "On macOS it is usually located at:\n"
    "  /Library/Frameworks/CsoundLib64.framework/Headers")
endif()

message(STATUS "Found Csound headers in: ${CSOUND_INCLUDE_DIR}")

# Opcode plugin: MODULE on macOS produces a .so (usually fine for the Csound plugin dir)
add_library(wtgivlocal MODULE wtgivlocal.c)

target_include_directories(wtgivlocal PRIVATE ${CSOUND_INCLUDE_DIR})

# No "lib" prefix
set_target_properties(wtgivlocal PROPERTIES PREFIX "")

# If you need .dylib instead of .so, uncomment:
# set_target_properties(wtgivlocal PROPERTIES SUFFIX ".dylib")