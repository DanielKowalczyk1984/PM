CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(GSBPA)
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

# Set the flags for compiling
add_definitions(-Wall -Wno-deprecated -pedantic -Wmissing-declarations -Wpointer-arith  -Wundef -Wcast-qual -Wcast-align -Wwrite-strings)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 m64")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -std=c99 -m64 -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -O0 -m64")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -std=c99 -m64 -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes")

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
FILE(GLOB SOURCE_FILES  "src/*.c" "src/*.cc")

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE ${GUROBI_INCLUDE_DIR} ${GLIB_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${COMMON_INCLUDES})
target_link_libraries(${PROJECT_NAME} ${GUROBI_LIBRARY} ${GLIB_LIBRARIES} ${Boost_LIBRARIES} -lm -lpthread)

# Include CppCcheck
include(CppcheckTargets)
include(customCppCheck)

# Include dev-tools
include(dev-tools)