list(APPEND CPPCHECK_TEMPLATE_ARG "--enable=warning,performance,style")
add_cppcheck(${PROJECT_NAME} FORCE)
add_cppcheck_indvidual(${PROJECT_NAME} FORCE)
