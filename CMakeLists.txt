cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(
  PM
  VERSION 0.1.0
  LANGUAGES CXX C
  DESCRIPTION "Branch-and-Price algorithm for parallel machine scheduling problems"
)

# Define a local module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_SOURCE_DIR})

#
# Prevent building in the source directory
#

include(ConfigSafeGuards)

#
# Set project options
#

include(StandardSettings)
# find_program(IWYU_PATH NAMES include-what-you-use iwyu) if(NOT IWYU_PATH) message(FATAL_ERROR
# "Could not find the program include-what-you-use") endif() set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE
# ${IWYU_PATH}) set(CMAKE_C_INCLUDE_WHAT_YOU_USE ${IWYU_PATH})
include(Utils)
include(StaticAnalyzers)
message(STATUS "Started CMake for ${PROJECT_NAME} v${PROJECT_VERSION}...\n")

#
# Enable package managers
#

include(Conan)
include(Vcpkg)
include(CPM)

#
# code coverage
#

find_package(codecov)

#
# Third Party packages
#

find_package(Gurobi)
include(CoinOR)
add_subdirectory(or-utils)
add_subdirectory(branch-and-bound)

#
# Create library, setup header and source files
#

# Find all headers and implementation files
include(SourcesAndHeaders)

if(${PROJECT_NAME}_BUILD_EXECUTABLE)
  add_executable(${PROJECT_NAME} ${exe_sources})
  add_library(${PROJECT_NAME}_LIB ${headers} ${sources})
elseif(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
  add_library(${PROJECT_NAME} INTERFACE)
else()
  add_library(${PROJECT_NAME} ${headers} ${sources})
endif()

verbose_message("Found the following header files:")
verbose_message(${headers})
verbose_message("Found the following source files:")

if(${PROJECT_NAME}_BUILD_EXECUTABLE)
  verbose_message(${exe_sources})
else()
  verbose_message(${sources})
endif()
message(STATUS "Added all header and implementation files.\n")

#
# Set the project standard and warnings
#

if(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
  target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)
else()
  target_compile_features(${PROJECT_NAME}_LIB PUBLIC cxx_std_20)
endif()

include(CompilerWarnings)
set_project_warnings(${PROJECT_NAME})
set_project_warnings(${PROJECT_NAME}_LIB)

verbose_message("Applied compiler warnings. Using standard ${CXX_STANDARD}.\n")

# RelWIthDebInfo Flags
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -mtune=native -fPIC")
set(CMAKE_C_FLAGS_RELWITHDEBINFO
    "${CMAKE_C_FLAGS_RELWITHDEBINGFO} -std=c99 -mtune=native -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes"
)

# Release Flags
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mtune=native  -fPIC")
set(CMAKE_C_FLAGS_RELEASE
    "${CMAKE_C_FLAGS_RELEASE} -std=c99 -D_XOPEN_SOURCE=500 -mtune=native -Wnested-externs -Wstrict-prototypes"
)

# Debug Flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -mtune=native -O0 -fPIC")
set(CMAKE_C_FLAGS_DEBUG
    "${CMAKE_C_FLAGS_DEBUG} -O0 -std=c99 -mtune=native -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes"
)

#
# Set the build/user include directories
#

# Allow usage of header files in the `src` directory, but only for utilities
if(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
  target_include_directories(
    ${PROJECT_NAME} INTERFACE $<INSTALL_INTERFACE:include>
                              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )
else()
  target_include_directories(
    ${PROJECT_NAME}
    PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
  )

  if(${PROJECT_NAME}_BUILD_EXECUTABLE)
    target_include_directories(
      ${PROJECT_NAME}_LIB
      PUBLIC $<INSTALL_INTERFACE:include>
             $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
             ${GUROBI_INCLUDE_DIR}
             ${CONAN_INCLUDE_DIRS_FMT}
             ${CONAN_INCLUDE_DIRS_BOOST}
             ${CONAN_INCLUDE_DIRS_DOCOPT.CPP}
             ${CONAN_INCLUDE_DIRS_RANGE-V3}
      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
    )

    target_link_libraries(
      ${PROJECT_NAME}_LIB
      ${GUROBI_LIBRARY}
      ${GUROBI_CXX_LIBRARY}
      CONAN_PKG::fmt
      CONAN_PKG::boost
      CONAN_PKG::docopt.cpp
      CONAN_PKG::range-v3
      ORUtils::ORUtils
      branch-and-bound::branch-and-bound
      -lm
      -lpthread
    )

    target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_LIB)
  endif()
endif()

message(STATUS "Finished setting up include directories.")

#
# Provide alias to library for
#

if(${PROJECT_NAME}_BUILD_EXECUTABLE)
  add_executable(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
else()
  add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
endif()

verbose_message("Project is now aliased as ${PROJECT_NAME}::${PROJECT_NAME}.\n")

#
# Format the project using the `clang-format` target (i.e: cmake --build build --target
# clang-format)
#

add_clang_format_target()

#
# Quick `ConfigVersion.cmake` creation
#

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

# include Doxygen
include(Doxygen)
