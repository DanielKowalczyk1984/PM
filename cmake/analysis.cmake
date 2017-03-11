# Make sure cppcheck binary is available
if( NOT CPPCHECK_FOUND )
    find_package(CppCheck)
endif()

# add a target for CppCheck
# _target - The name of the project that this is for. Will generate ${_target}_analysis 
# _sources - The name of the variable holding the sources list. 
#            This is the name of the variable not the actual list
#
# Macro instead of function to make the PARENT_SCOPE stuff easier
macro(add_analysis _target _sources)
    if( CPPCHECK_FOUND )

        # Get the compile_commands
        set(cppcheck_project --project=${CMAKE_BINARY_DIR}/compile_commands.json)

        # Add to the all target to have a high level "make analysis"
        LIST(APPEND ALL_ANALYSIS_TARGETS "${_target}_analysis")
        set(ALL_ANALYSIS_TARGETS "${ALL_ANALYSIS_TARGETS}" )

        # This is used to make the command run correctly on the command line.
        # The COMMAND argument expects a list and this does the change
        # I need to check which version works with 2.7
        set(CPPCHECK_ARG "--enable=information")
        if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VESION} GREATER 2.7)
            separate_arguments(tmp_args UNIX_COMMAND ${CPPCHECK_ARG})
        else ()
            # cmake 2.6 has different arguments 
            string(REPLACE " " ";" tmp_args ${CPPCHECK_ARG})         
        endif ()

        # add a custom _target_analysis target
        add_custom_target(${_target}_analysis)
        set_target_properties(${_target}_analysis PROPERTIES EXCLUDE_FROM_ALL TRUE)
        
        # add the cppcheck command to the target
        add_custom_command(TARGET ${_target}_analysis PRE_BUILD
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
            COMMAND ${CPPCHECK_BIN} ${tmp_args} ${cppcheck_project} ${${_sources}}
            DEPENDS ${${_sources}}
            COMMENT "Running cppcheck: ${_target}"
            VERBATIM)
        message("adding cppcheck analysys target for ${_target}")
    endif()

endmacro()