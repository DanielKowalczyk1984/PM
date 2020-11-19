cmake_minimum_required(VERSION 3.13...3.16 FATAL_ERROR)
project(
  or_utils
  VERSION 1.0.0
  LANGUAGES C)

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

# ##############################################################################
# Create target and set properties

add_library(
  or_utils
  src/alloc.c
  src/binomial-heap.c
  src/copy.c
  src/graph.c
  src/heap.c
  src/ksubset.c
  src/sortus.c
  src/time.c
  src/time_bis.c
  src/util.c
  src/lp.c)

# Add an alias so that library can be used inside the build tree, e.g. when
# testing
add_library(ORUtils::ORUtils ALIAS or_utils)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake
                      ${CMAKE_SOURCE_DIR})
find_package(Gurobi)

# Set target properties
target_include_directories(
  or_utils
  PUBLIC $<INSTALL_INTERFACE:include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${GUROBI_INCLUDE_DIR})

target_link_libraries(or_utils ${GUROBI_LIBRARY} -lm -lpthread)

# target_compile_features(or_utils PRIVATE c_variadic_macro)
target_compile_options(
  or_utils
  PRIVATE
    $<$<OR:$<C_COMPILER_ID:Clang>,$<C_COMPILER_ID:AppleClang>,$<C_COMPILER_ID:GNU>>:-Wall
    -Wextra
    -Wpedantic>)

# ##############################################################################
# Installation instructions

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/ORUtils)

install(
  TARGETS or_utils
  EXPORT or_utils-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# This is required so that the exported target has the name JSONUtils and not
# jsonutils
set_target_properties(or_utils PROPERTIES EXPORT_NAME ORUtils)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export the targets to a script
install(
  EXPORT or_utils-targets
  FILE ORUtilsTargets.cmake
  NAMESPACE ORUtils::
  DESTINATION ${INSTALL_CONFIGDIR})

# Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/ORUtilsConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake/ORUtilsConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/ORUtilsConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CONFIGDIR})

# Install the config, configversion and custom find modules
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ORUtilsConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/ORUtilsConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR})

# ##############################################################################
# Exporting from the build tree

export(
  EXPORT or_utils-targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/ORUtilsTargets.cmake
  NAMESPACE ORUtils::)

# Register package in user's package registry
export(PACKAGE ORUtils)

# ##############################################################################
# Add test enable_testing() add_subdirectory(test)
