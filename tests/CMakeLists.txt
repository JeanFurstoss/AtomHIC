### CMakeLists.txt --- 
##
######################################################################
## 
### Commentary: Test oxford library
## 
######################################################################

#First test
add_test(test ${EXECUTABLE_OUTPUT_PATH}/oxfordNewStudent)

#Second test
add_test(student1 ${EXECUTABLE_OUTPUT_PATH}/oxfordNewStudent "Jim" "Clara")

#Third test
add_test(student2 ${EXECUTABLE_OUTPUT_PATH}/oxfordNewStudent "William" "Lincoln")

# Specific target that launches tests in a verbose mode (windows)
if(WIN32)
  add_custom_target(RUN_TESTS_VERBOSE
    ${CMAKE_COMMAND} -E env CTEST_OUTPUT_ON_FAILURE=1
    ${CMAKE_CTEST_COMMAND} -C Release --verbose
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
