cmake_minimum_required(VERSION 3.16)
project(wtgivkura C)

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
  message(FATAL_ERROR "Could not find csdl.h (Csound headers).")
endif()

add_library(wtgivkura MODULE wtgivkura.c)
target_include_directories(wtgivkura PRIVATE ${CSOUND_INCLUDE_DIR})
set_target_properties(wtgivkura PROPERTIES PREFIX "")