if(${PROJECT_NAME}_ENABLE_CONAN)
  #
  # Setup Conan requires and options here:
  #

  set(CONAN_REQUIRES "boost/1.73.0
                      fmt/7.0.3
                      glib/2.67.0")

  set(PM_CONAN_OPTIONS "pcre:with_utf=True
                        pcre:with_unicode_properties=True")

  #
  # If `conan.cmake` (from https://github.com/conan-io/cmake-conan) does not
  # exist, download it.
  #
  if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(
      STATUS
        "Downloading conan.cmake from https://github.com/conan-io/cmake-conan..."
    )
    file(DOWNLOAD
         "https://github.com/conan-io/cmake-conan/raw/v0.15/conan.cmake"
         "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Cmake-Conan downloaded succesfully.")
  endif()

  include(${CMAKE_BINARY_DIR}/conan.cmake)

  conan_cmake_run(
    REQUIRES ${CONAN_REQUIRES}
    OPTIONS ${PM_CONAN_OPTIONS}
    BASIC_SETUP CMAKE_TARGETS # Individual targets to link to
    BUILD missing)

  verbose_message("Conan is setup and all requires have been installed.")
endif()
