### cimlibTest.cmake --- 
## 
######################################################################
## 
### Commentary: 
## 
######################################################################

######################################################################
## Test
######################################################################

set(EXECUTABLE ${EXECUTABLE_OUTPUT_PATH}/${PROJECT_NAME})

######################################################################
## Executable used to compare result files (and arguments)
######################################################################

if("${NUMERICAL_DIFFERENCE_EXECUTABLE}" STREQUAL "")
  # Default Compare
  set(NUMERICAL_DIFFERENCE_EXECUTABLE "${EXECUTABLE_OUTPUT_PATH}/compare")
endif()

if("${NUMERICAL_DIFFERENCE_TOLERANCE}" STREQUAL "")
  set(NUMERICAL_DIFFERENCE_TOLERANCE 0.01)
endif()

######################################################################
## Test name
######################################################################

file(RELATIVE_PATH file_path "${CMAKE_SOURCE_DIR}/tests" "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")

get_filename_component(dir_path ${file_path} PATH)
set(TEST_NAME)
while(NOT "${dir_path}" STREQUAL "")
  get_filename_component(NAME ${dir_path} NAME)
  if("${TEST_NAME}" STREQUAL "")
    set(TEST_NAME ${NAME})
  else()
    set(TEST_NAME ${NAME}_${TEST_NAME})
  endif()
  get_filename_component(dir_path ${dir_path} PATH)
endwhile()

######################################################################
## Copy input files in working directory of test
######################################################################

file(GLOB_RECURSE INPUT_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.dat *.txt *.xsf *.lmp *.cfg *.ath)

list(REMOVE_ITEM INPUT_FILES "CMakeLists.txt")

list(LENGTH INPUT_FILES NUMBER_OF_INPUT_FILES)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ForTest)
	file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ForTest DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

######################################################################
## Configure test
######################################################################

file(GLOB_RECURSE REFERENCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} Ref_*)

if(WIN32)
  set(EXTENSION bat)
else()
  set(EXTENSION sh)
endif()

unset(RESULT_FILES)
foreach(FILE ${REFERENCE_FILES})
  if (${NUMBER_OF_INPUT_FILES} GREATER 0)
    # Remove reference file from input files
    list(REMOVE_ITEM INPUT_FILES ${FILE})
    list(LENGTH INPUT_FILES NUMBER_OF_INPUT_FILES)
  endif()
  unset(DIR)
  get_filename_component(DIR ${FILE} DIRECTORY)
  if(DIR)
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/CMakeLists.txt")
      continue()
    endif()
    set(DIR "${DIR}/")
  endif()
  get_filename_component(FILE ${FILE} NAME)
  string(REPLACE Ref_ "" FILE ${FILE})
  if(RESULT_FILES)
    set(RESULT_FILES "${RESULT_FILES} ${DIR}${FILE}")
  else()
    set(RESULT_FILES "${DIR}${FILE}")
  endif()
endforeach()

# Copy input files
if (${NUMBER_OF_INPUT_FILES} GREATER 0)
  foreach(FILE ${INPUT_FILES})
    file(COPY ${FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  endforeach()
endif()

######################################################################
## Create test
######################################################################

add_test(NAME ${PROJECT_NAME}
	COMMAND ${EXECUTABLE_OUTPUT_PATH}/AtomHic_test.${EXTENSION} -e ${EXECUTABLE} -i ${CMAKE_CURRENT_SOURCE_DIR} -c ${NUMERICAL_DIFFERENCE_EXECUTABLE} -r "${RESULT_FILES}" -t ${NUMERICAL_DIFFERENCE_TOLERANCE}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

######################################################################
## Create tests list
######################################################################

if(WIN32)
  set(LIST_OF_TESTS "${LIST_OF_TESTS};${TEST_NAME}" CACHE STRING "List of possible tests" FORCE)
endif()

