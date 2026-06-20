cmake_minimum_required(VERSION 3.16)
project(wtgivlocal_opcode C)

# Build type di default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Trova csdl.h
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
    "Impossibile trovare csdl.h.\n"
    "Su macOS di solito è in:\n"
    "  /Library/Frameworks/CsoundLib64.framework/Headers")
endif()

message(STATUS "Found Csound headers in: ${CSOUND_INCLUDE_DIR}")

# Plugin opcode: MODULE su macOS produce .so (di solito ok per Csound plugin dir)
add_library(wtgivlocal MODULE wtgivlocal.c)

target_include_directories(wtgivlocal PRIVATE ${CSOUND_INCLUDE_DIR})

# Niente prefisso "lib"
set_target_properties(wtgivlocal PROPERTIES PREFIX "")

# Se ti serve .dylib invece di .so, decommenta:
# set_target_properties(wtgivlocal PROPERTIES SUFFIX ".dylib")