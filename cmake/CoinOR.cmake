if(${PROJECT_NAME}_ENABLE_COIN_OR)

  if(NOT EXISTS "${CMAKE_BINARY_DIR}/coin_or/coinbrew")
    message(
      STATUS
        "Downloading coinbrew from https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew"
    )

    file(DOWNLOAD "https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew"
         "${CMAKE_BINARY_DIR}/coin_or/coinbrew"
    )

    file(CHMOD ${CMAKE_BINARY_DIR}/coin_or/coinbrew FILE_PERMISSIONS OWNER_READ OWNER_WRITE
         OWNER_EXECUTE
    )
    message(STATUS "coinbrew downloaded successfully.")

    execute_process(
      COMMAND ./coinbrew fetch osi@master WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/coin_or
    )

    set(COIN_GUROBI_CFLAGS "-I$ENV{GUROBI_HOME}/include")
    set(COIN_GUROBI_LFLAGS "-L$ENV{GUROBI_HOME}/lib -lgurobi91 -lgurobi_c++ -lm -lpthread")

    execute_process(
      COMMAND ./coinbrew build osi -t none --with-gurobi --with-gurobi-lflags=${COIN_GUROBI_LFLAGS}
              --with-gurobi-cflags=${COIN_GUROBI_CFLAGS}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/coin_or
    )

  endif()

  find_library(
    COIN_OSI_LIBRARY
    NAMES Osi OsiGrb CoinUtils
    HINTS ${CMAKE_BINARY_DIR}/coin_or/dist/lib
  )

  find_library(
    COIN_OSI_GRB_LIBRARY
    NAMES OsiGrb
    HINTS ${CMAKE_BINARY_DIR}/coin_or/dist/lib
  )
  set(COIN_INCLUDE_DIR "${CMAKE_BINARY_DIR}/coin_or/dist/include/coin-or")

endif(${PROJECT_NAME}_ENABLE_COIN_OR)
