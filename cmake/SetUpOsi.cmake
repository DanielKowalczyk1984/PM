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

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/ThirdParty/coin-or-x64-linux-release")
  set(TMP_GRB_LIB "-L$ENV{GUROBI_HOME}/lib -lgurobi91 -lpthread -lm")
  execute_process(
    COMMAND ./coinbrew fetch Cgl@master WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ThirdParty
  )
  execute_process(
    COMMAND ./coinbrew build Cgl --with-gurobi-lib=${TMP_GRB_LIB}
            --with-gurobi-incdir=$ENV{GUROBI_HOME}/include --tests none
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ThirdParty
  )
  message(STATUS "Build Osi completed.")
endif()

set(coin_modules "Osi;CoinUtils;OsiGrb;Cgl")

if(${CMAKE_HOST_UNIX})
  find_path(
    OSI_INCLUDE_DIR
    NAMES OsiGrbSolverInterface.hpp
    PATHS "${CMAKE_SOURCE_DIR}/ThirdParty/coin-or-x64-linux-release/include/coin-or" REQUIRED
  )

  foreach(coin_lib ${coin_modules})
    string(TOUPPER "${coin_lib}" coin_lib_toupper)
    find_library(
      ${coin_lib_toupper}_LIBRARY
      NAMES ${coin_lib}
      PATHS "${CMAKE_SOURCE_DIR}/ThirdParty/coin-or-x64-linux-release/lib"
    )

    find_library(
      ${coin_lib_toupper}_LIBRARY_DEBUG
      NAMES ${coin_lib}
      PATHS "${CMAKE_SOURCE_DIR}/ThirdParty/coin-or-x64-linux-debug/lib"
    )

  endforeach()
elseif(${CMAKE_HOST_WIN32})
  find_path(
    OSI_INCLUDE_DIR
    NAMES OsiGrbSolverInterface.hpp
    PATHS "${CMAKE_SOURCE_DIR}/ThirdParty/coin-or-x64-MD/include/coin-or" REQUIRED
  )

  foreach(coin_lib ${coin_modules})
    string(TOUPPER "${coin_lib}" coin_lib_toupper)
    find_library(
      ${coin_lib_toupper}_LIBRARY
      NAMES ${coin_lib}
      PATHS "${CMAKE_SOURCE_DIR}/ThirdParty/coin-or-x64-MD/lib"
    )

    find_library(
      ${coin_lib_toupper}_LIBRARY_DEBUG
      NAMES ${coin_lib}
      PATHS "${CMAKE_SOURCE_DIR}/ThirdParty/coin-or-x64-MDd/lib"
    )

  endforeach()

endif()

include(FindPackageHandleStandardArgs)
foreach(coin_lib ${coin_modules})
  string(TOUPPER "${coin_lib}" coin_lib_toupper)

  find_package_handle_standard_args(${coin_lib_toupper} REQUIRED_VARS ${coin_lib_toupper}_LIBRARY)

  find_package_handle_standard_args(
    ${coin_lib_toupper}_DEBUG REQUIRED_VARS ${coin_lib_toupper}_LIBRARY_DEBUG
  )

endforeach()

mark_as_advanced(OSI_INCLUDE_DIR OSI_LIBRARY COINUTILS_LIBRARY OSIGRB_LIBRARY CGL_LIBRARY)
