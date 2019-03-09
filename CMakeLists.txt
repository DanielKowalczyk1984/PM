CMAKE_MINIMUM_REQUIRED(VERSION 3.13)
project(PM)
include(FetchContent)
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

# Download the tdzdd library
FetchContent_Declare(
    project_TdZdd
    GIT_REPOSITORY https://github.com/DanielKowalczyk1984/TdZdd.git
    GIT_TAG "master" 
)

FetchContent_GetProperties(project_TdZdd)
if(NOT project_TdZdd_POPULATED)
    FetchContent_Populate(project_TdZdd)
    set(TdZDD_INCLUDE_DIR ${project_tdzdd_SOURCE_DIR}/include)
endif()

# Set the flags for compiling
add_definitions(-Wall -Wno-deprecated -pedantic -Wmissing-declarations -Wpointer-arith  -Wcast-qual -Wcast-align -Wwrite-strings)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++14 -m64")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -std=c99 -m64 -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 -O0 -m64")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -std=c99 -m64 -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes")

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
file(GLOB SOURCE_FILES "src/*.cc" "src/*.cpp" "src/*.c" )

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${TdZDD_INCLUDE_DIR} ${GUROBI_INCLUDE_DIR} ${GLIB_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${COMMON_INCLUDES})
target_link_libraries(${PROJECT_NAME} PUBLIC ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} ${GLIB_LIBRARIES} ${Boost_LIBRARIES} -lm -lpthread)

# Include CppCcheck
include(CppcheckTargets)
include(customCppCheck)

# Include dev-tools
include(dev-tools)
