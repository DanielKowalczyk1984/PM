if(NOT EXISTS "${CMAKE_SOURCE_DIR}/ThirdParty/coinbrew")
  message(STATUS "Downloading coinbrew from https://raw.githubusercontent.com/coin-or/coinbrew/...")
  file(DOWNLOAD "https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew"
       "${CMAKE_SOURCE_DIR}/ThirdParty/coinbrew"
  )
  file(CHMOD "${CMAKE_SOURCE_DIR}/ThirdParty/coinbrew" PERMISSIONS OWNER_READ OWNER_WRITE
       OWNER_EXECUTE
  )
  message(STATUS "Coinbrew downloaded successfully.")
endif()

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/ThirdParty/dist")
  set(TMP_GRB_LIB "-L$ENV{GUROBI_HOME}/lib -lgurobi91 -lpthread -lm")
  execute_process(
    COMMAND ./coinbrew build Osi --with-gurobi-lib=${TMP_GRB_LIB}
            --with-gurobi-incdir=$ENV{GUROBI_HOME}/include --tests none
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ThirdParty
  )
  message("Build Osi completed.")
endif()

find_path(
  OSI_INCLUDE_DIR
  NAMES OsiGrbSolverInterface.hpp
  PATHS "${CMAKE_SOURCE_DIR}/ThirdParty/dist/include/coin"
)

find_library(
  OSI_LIBRARY
  NAMES Osi
  PATHS "${CMAKE_SOURCE_DIR}/ThirdParty/dist/lib"
)

find_library(
  COINUtils_LIBRARY
  NAMES CoinUtils
  PATHS "${CMAKE_SOURCE_DIR}/ThirdParty/dist/lib"
)

find_library(
  OSI_GRB_LIBRARY
  NAMES OsiGrb
  PATHS "${CMAKE_SOURCE_DIR}/ThirdParty/dist/lib"
)

# Version detection
file(READ "${OSI_INCLUDE_DIR}/OsiConfig.h" OSI_CONFIG_H_CONTENTS)
string(REGEX MATCH "#define OSI_VERSION_MAJOR *([0-9]+)" _dummy "${OSI_CONFIG_H_CONTENTS}")
set(OSI_VERSION_MAJOR "${CMAKE_MATCH_1}")
string(REGEX MATCH "#define OSI_VERSION_MINOR *([0-9]+)" _dummy "${OSI_CONFIG_H_CONTENTS}")
set(OSI_VERSION_MINOR "${CMAKE_MATCH_1}")
string(REGEX MATCH "#define OSI_VERSION_RELEASE *([0-9]+)" _dummy "${OSI_CONFIG_H_CONTENTS}")
set(OSI_VERSION_RELEASE "${CMAKE_MATCH_1}")
set(OSI_VERSION "${OSI_VERSION_MAJOR}.${OSI_VERSION_MINOR}.${OSI_VERSION_RELEASE}")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
  OSI
  REQUIRED_VARS OSI_INCLUDE_DIR OSI_LIBRARY COINUtils_LIBRARY OSI_GRB_LIBRARY
  VERSION_VAR OSI_VERSION
)

mark_as_advanced(OSI_INCLUDE_DIR OSI_LIBRARY COINUtils_LIBRARY OSI_GRB_LIBRARY)
