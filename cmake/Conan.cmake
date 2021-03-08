if(${PROJECT_NAME}_ENABLE_CONAN)
  #
  # Setup Conan requires and options here:
  #

  set(CONAN_REQUIRES "boost/1.73.0" "fmt/7.0.3" "docopt.cpp/0.6.3" "glib/2.67.0")

  set(CONAN_OPTIONS "pcre:with_utf=True" "pcre:with_unicode_properties=True")

  #
  # If `conan.cmake` (from https://github.com/conan-io/cmake-conan) does not exist, download it.
  #
  if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan...")
    file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/v0.16.1/conan.cmake"
         "${CMAKE_BINARY_DIR}/conan.cmake"
    )
    message(STATUS "Cmake-Conan downloaded successfully.")
  endif()

  include(${CMAKE_BINARY_DIR}/conan.cmake)

  conan_cmake_configure(REQUIRES ${CONAN_REQUIRES} OPTIONS ${CONAN_OPTIONS} GENERATORS cmake)
  conan_cmake_autodetect(settings)
  conan_cmake_install(PATH_OR_REFERENCE . BUILD missing SETTINGS ${settings})
  include(${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup(TARGETS)

  verbose_message("Conan is setup and all requires have been installed. ")
endif()
