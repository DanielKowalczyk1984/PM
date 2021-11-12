find_path(
  GUROBI_INCLUDE_DIR
  NAMES gurobi_c.h
  HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
  PATH_SUFFIXES include
)

find_library(
  GUROBI_LIBRARY
  NAMES gurobi gurobi91
  HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
  PATH_SUFFIXES lib
)

if(MSVC)
  # determine Visual Studio year
  if(MSVC_TOOLSET_VERSION EQUAL 142)
    set(MSVC_YEAR "2019")
  elseif(MSVC_TOOLSET_VERSION EQUAL 141)
    set(MSVC_YEAR "2017")
  elseif(MSVC_TOOLSET_VERSION EQUAL 140)
    set(MSVC_YEAR "2015")
  endif()
  string(FIND "${CMAKE_MSVC_RUNTIME_LIBRARY}" "DLL" FOUND_DLL)

  if(${FOUND_DLL} EQUAL -1)
    set(M_FLAG "mt")
  else()
    set(M_FLAG "md")
  endif()

  find_library(
    GUROBI_CXX_LIBRARY
    NAMES gurobi_c++${M_FLAG}${MSVC_YEAR}
    HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES lib
  )

  find_library(
    GUROBI_CXX_DEBUG_LIBRARY
    NAMES gurobi_c++${M_FLAG}d${MSVC_YEAR}
    HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES lib
  )

else()
  find_library(
    GUROBI_CXX_LIBRARY
    NAMES gurobi_c++
    HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES lib
  )
  set(GUROBI_CXX_DEBUG_LIBRARY ${GUROBI_CXX_LIBRARY})
endif()

set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}")

# Version detection
file(READ "${GUROBI_INCLUDE_DIR}/gurobi_c.h" GUROBI_C_H_CONTENTS)
string(REGEX MATCH "#define GRB_VERSION_MAJOR *([0-9]+)" _dummy "${GUROBI_C_H_CONTENTS}")
set(GRB_VERSION_MAJOR "${CMAKE_MATCH_1}")
string(REGEX MATCH "#define GRB_VERSION_MINOR *([0-9]+)" _dummy"${GUROBI_C_H_CONTENTS}")
set(GRB_VERSION_MINOR "${CMAKE_MATCH_1}")
string(REGEX MATCH "#define GRB_VERSION_TECHNICAL *([0-9]+)" _dummy "${GUROBI_C_H_CONTENTS}")
set(GRB_VERSION_MICRO "${CMAKE_MATCH_1}")
set(GRB_VERSION "${GRB_VERSION_MAJOR}.${GRB_VERSION_MINOR}.${GRB_VERSION_MICRO}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Gurobi
  REQUIRED_VARS GUROBI_INCLUDE_DIR GUROBI_LIBRARIES
  VERSION_VAR GRB_VERSION
)
