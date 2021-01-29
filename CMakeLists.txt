cmake_minimum_required(VERSION 3.15)

project(
  PM
  VERSION 0.1.0
  LANGUAGES CXX C
  DESCRIPTION "Branch-and-Price algorithm for parallel machine scheduling problems")

#
# Set project options
#

# Define a local module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake
                      ${CMAKE_SOURCE_DIR})

include(StandardSettings)
include(Utils)
include(StaticAnalyzers)
message(STATUS "Started CMake for ${PROJECT_NAME} v${PROJECT_VERSION}...\n")

#
# Prevent building in the source directory
#

include(ConfigSafeGuards)

#
# Enable package managers
#

include(Conan)
include(Vcpkg)

#
# code coverage
#

find_package(codecov)

#
# Create library, setup header and source files
#

# Find all headers and implementation files
include(cmake/SourcesAndHeaders.cmake)

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

# Find Gurobi
find_package(Gurobi)

add_subdirectory(or-utils)
add_subdirectory(branch-and-bound)

#
# Set the project standard and warnings
#

if(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
  target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)
else()
  target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)
endif()

include(CompilerWarnings)
set_project_warnings(${PROJECT_NAME})

verbose_message("Applied compiler warnings. Using standard ${CXX_STANDARD}.\n")

# Set the flags for compiling add_definitions(-Wall -Wno-deprecated -pedantic
# -Wmissing-declarations -Wpointer-arith -Wcast-qual -Wcast-align
# -Wwrite-strings)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O5 -std=c++20 -m64")
set(CMAKE_C_FLAGS_RELEASE
    "${CMAKE_C_FLAGS_RELEASE} -std=c99 -m64 -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes"
)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++20 -O0 -m64")
set(CMAKE_C_FLAGS_DEBUG
    "${CMAKE_C_FLAGS_DEBUG} -O0 -std=c99 -m64 -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes"
)

#
# Set the build/user include directories
#

# Allow usage of header files in the `src` directory, but only for utilities
if(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
  target_include_directories(
    ${PROJECT_NAME}
    INTERFACE $<INSTALL_INTERFACE:include>
              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
else()
  target_include_directories(
    ${PROJECT_NAME}
    PUBLIC $<INSTALL_INTERFACE:include>
           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

  if(${PROJECT_NAME}_BUILD_EXECUTABLE)
    target_include_directories(
      ${PROJECT_NAME}_LIB
      PUBLIC $<INSTALL_INTERFACE:include>
             $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
             ${GUROBI_INCLUDE_DIR}
             CONAN_PKG::fmt
             CONAN_PKG::glib
             CONAN_PKG::boost
             ORUtils::ORUtils
             branch-and-bound
      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

    target_link_libraries(
      ${PROJECT_NAME}_LIB
      ${GUROBI_LIBRARY}
      ${GUROBI_CXX_LIBRARY}
      CONAN_PKG::fmt
      CONAN_PKG::glib
      CONAN_PKG::boost
      ORUtils::ORUtils
      branch-and-bound
      -lm
      -lpthread)

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
# Format the project using the `clang-format` target (i.e: cmake --build build
# --target clang-format)
#

add_clang_format_target()

#
# Quick `ConfigVersion.cmake` creation
#

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

# Include CppCcheck
include(CppcheckTargets)
include(customCppCheck)

# Include dev-tools
include(dev-tools)

# include Doxygen
include(cmake/Doxygen.cmake)
# set_property(TARGET ${PROJECT_NAME}
# PROPERTY OUTPUT_NAME "${PROJECT_NAME}$<CONFIG>${_output}" )
