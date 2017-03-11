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

# Find CppCheck
find_package(CppCheck)

# static analysis. Should be before adding subprojects
set (ALL_ANALYSIS_TARGETS)

# Set the flags for compiling
add_definitions(-Wall -Wno-deprecated -pedantic -Wmissing-declarations -Wpointer-arith  -Wundef -Wcast-qual -Wcast-align -Wwrite-strings)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 m64")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -std=c99 -m64 -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -O0 -m64")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -std=c99 -m64 -D_XOPEN_SOURCE=500 -Wnested-externs -Wstrict-prototypes")

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
FILE(GLOB SOURCE_FILES  "src/*.c" "src/*.cc")
message(${GUROBI_LIBRARY})
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE ${GUROBI_INCLUDE_DIR} ${GLIB_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${COMMON_INCLUDES})
target_link_libraries(${PROJECT_NAME} ${GUROBI_LIBRARY} ${GLIB_LIBRARIES} ${Boost_LIBRARIES} -lm -lpthread)

include(${CMAKE_SOURCE_DIR}/cmake/analysis.cmake)
add_analysis(${PROJECT_NAME} SOURCE_FILES)

# Add the "make analysis" target
if( CPPCHECK_FOUND )
    add_custom_target(analysis)
    ADD_DEPENDENCIES(analysis ${ALL_ANALYSIS_TARGETS})
    set_target_properties(analysis PROPERTIES EXCLUDE_FROM_ALL TRUE)
    message("analysis analysis targets are ${ALL_ANALYSIS_TARGETS}")
endif()

