CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(GSBPA)
include(ExternalProject)

# add json export compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Define a local module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Find Glib 
find_package(Glib)

# Find Gurobi
find_package(Gurobi)

# Find Boost
find_package(Boost)


# Include EGlib
ExternalProject_Add(project_EGlib
	URL http://www.dii.uchile.cl/~daespino/SOurce/EGlib.tar.bz2
	PREFIX ${CMAKE_SOURCE_DIR}/third_party/EGlib
	CONFIGURE_COMMAND autoheader -v  && autoconf -v  && ./configure && sed -i "s/^MAINOPT.*//" Makefile.common && sed -i "s/\(MAINOPT\)//g" Makefile.common
	BUILD_IN_SOURCE 1
	BUILD_COMMAND make library
	INSTALL_COMMAND ""
)

ExternalProject_Get_Property(project_EGlib source_dir)
set_target_properties(project_EGlib PROPERTIES EXCLUDE_FROM_ALL TRUE)
set(EGlib_INCLUDE_DIR ${source_dir}/include)
set(EGlib_LIBRARY_DIR ${source_dir}/lib/libEGlib.a)

ExternalProject_Add(project_TdZdd
    GIT_REPOSITORY https://github.com/DanielKowalczyk1984/TdZdd.git
    GIT_TAG "master"
    PREFIX ${CMAKE_SOURCE_DIR}/third_party/TdZDD
    BUILD_COMMAND ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND ""
    EXCLUDE_FROM_ALL
)

ExternalProject_Get_Property(project_TdZdd source_dir)
set_target_properties(project_TdZdd PROPERTIES EXCLUDE_FROM_ALL TRUE)
set(TdZDD_INCLUDE_DIR ${source_dir}/include)

# Set the flags for compiling
add_definitions(-Wall -Wno-deprecated -pedantic -Wmissing-declarations -Wpointer-arith  -Wcast-qual -Wcast-align -Wwrite-strings)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 m64")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -std=c99 -m64 -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -O0 -m64")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -std=c99 -m64 -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes")

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
file(GLOB SOURCE_FILES  "src/*.c" "src/*.cc")

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE ${TdZDD_INCLUDE_DIR} ${EGlib_INCLUDE_DIR} ${GUROBI_INCLUDE_DIR} ${GLIB_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${COMMON_INCLUDES})
target_link_libraries(${PROJECT_NAME} ${EGlib_LIBRARY_DIR} ${GUROBI_LIBRARY} ${GLIB_LIBRARIES} ${Boost_LIBRARIES} -lm -lpthread)

# Include CppCcheck
include(CppcheckTargets)
include(customCppCheck)

# Include dev-tools
include(dev-tools)