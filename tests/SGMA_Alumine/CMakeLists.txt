project(test_SGMA_Alumine)

## #################################################################
## Sources
## #################################################################

set(${PROJECT_NAME}_SOURCES main.cc)

## #################################################################
## Build rules
## #################################################################

#The file(GLOB...) allows for wildcard additions:
file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cpp")
#add_subdirectory(ForTest)
#include_directories(ForTest)

#Generate executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_link_libraries(${PROJECT_NAME} AtomHic)

include(AtomHicTest)
