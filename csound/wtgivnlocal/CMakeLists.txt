cmake_minimum_required(VERSION 3.16)
project(wtgivnlocal_opcode C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

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
  message(FATAL_ERROR "Impossibile trovare csdl.h (Csound headers).")
endif()

message(STATUS "Found Csound headers in: ${CSOUND_INCLUDE_DIR}")

add_library(wtgivnlocal MODULE wtgivnlocal.c)
target_include_directories(wtgivnlocal PRIVATE ${CSOUND_INCLUDE_DIR})
set_target_properties(wtgivnlocal PROPERTIES PREFIX "")