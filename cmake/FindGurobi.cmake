if(GUROBI_INCLUDE_DIR)
  set(GUROBI_FOUND TRUE)
  set(GUROBI_LIBRARIES "${GUROBI_LIBRARY};${GUROBI_CXX_LIBRARY}")
else(GUROBI_INCLUDE_DIR)
  find_path(
    GUROBI_INCLUDE_DIR
    NAMES gurobi_c++.h gurobi_c.h
    PATHS "$ENV{GUROBI_HOME}/include" "/opt/gurobi903/linux64/include"
  )

  find_library(
    GUROBI_LIBRARY
    NAMES gurobi81
          gurobi90
          gurobi91
          gurobi80
          gurobi75
          gurobi65
          gurobi60
          gurobi
          gurobi45
          gurobi46
          gurobi50
          gurobi51
          gurobi52
          gurobi55
          gurobi56
          gurobi60
          gurobi563
    PATHS "$ENV{GUROBI_HOME}/lib" "/opt/gurobi903/linux64/lib"
  )

  find_library(
    GUROBI_CXX_LIBRARY
    NAMES gurobi_c++ libgurobi75 libgurobi56 libgurobi
    PATHS "$ENV{GUROBI_HOME}/lib" "/opt/gurobi903/linux64/lib"
  )

  set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}")

  # Version detection
  file(READ "${GUROBI_INCLUDE_DIR}/gurobi_c.h" GUROBI_C_H_CONTENTS)
  string(REGEX MATCH "#define GRB_VERSION_MAJOR *([0-9]+)" _dummy "${GUROBI_C_H_CONTENTS}")
  set(GRB_VERSION_MAJOR "${CMAKE_MATCH_1}")
  string(REGEX MATCH "#define GRB_VERSION_MINOR *([0-9]+)" _dummy "${GUROBI_C_H_CONTENTS}")
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

  mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARIES)

endif(GUROBI_INCLUDE_DIR)
