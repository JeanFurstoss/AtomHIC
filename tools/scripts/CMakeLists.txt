### CMakeLists.txt --- 
######################################################################
## 
### Commentary: 
## 
######################################################################

## #################################################################
## Copy AtomHic_test script in binary directory
## #################################################################


set(EXTENSION sh)

set(TEST_SCRIPT AtomHic_test.${EXTENSION})

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_SCRIPT}
  ${EXECUTABLE_OUTPUT_PATH}/${TEST_SCRIPT}
  COPYONLY IMMEDIATE)
