list(APPEND CPPCHECK_TEMPLATE_ARG "--enable=warning,performance,style")
add_cppcheck(${PROJECT_NAME} FORCE)