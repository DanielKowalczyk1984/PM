# Additional targets to perform clang-format/clang-tidy
# Get all project files
get_target_property(ALL_CXX_SOURCE_FILES ${PROJECT_NAME} SOURCES)
file(GLOB_RECURSE INCLUDES "${CMAKE_SOURCE_DIR}/include/*.h*")

# Adding clang-format target if executable is found
find_program(CLANG_FORMAT "clang-format")
if(CLANG_FORMAT)
  add_custom_target(
    clang-format
    COMMAND /usr/bin/clang-format
    -i
    -style=file
    ${ALL_CXX_SOURCE_FILES}
    )
    set_target_properties(clang-format PROPERTIES EXCLUDE_FROM_ALL TRUE)
    add_custom_target(
      clang-format-inc
      COMMAND /usr/bin/clang-format
      -i
      -style=file
      ${INCLUDES}
      )
    set_target_properties(clang-format-inc PROPERTIES EXCLUDE_FROM_ALL TRUE)
    foreach(_file ${ALL_CXX_SOURCE_FILES})
        get_filename_component(_name ${_file} NAME_WE)
        add_custom_target(clang-format-${_name})
        add_custom_command(
            TARGET
            clang-format-${_name}
            COMMAND /usr/bin/clang-format
            -i
            -style=file
            ${_file}
            )
        set_target_properties(clang-format-${_name} PROPERTIES EXCLUDE_FROM_ALL TRUE)
    endforeach()
    foreach(_file ${INCLUDES})
        get_filename_component(_name ${_file} NAME_WE)
        add_custom_target(clang-format-inc-${_name})
        add_custom_command(
            TARGET
            clang-format-inc-${_name}
            COMMAND /usr/bin/clang-format
            -i
            -style=file
            ${_file}
            )
        set_target_properties(clang-format-inc-${_name} PROPERTIES EXCLUDE_FROM_ALL TRUE)
    endforeach()
endif()
