### CMakeLists.txt --- 
##
######################################################################
## 
### Commentary: Build of AtomHic library
## 
######################################################################

#Project name
project(AtomHic)
 
#The file(GLOB...) allows for wildcard additions:
file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cpp")

## #################################################################
## Build rules
## #################################################################

#Choose library type statique ou dynamique
if(BUILD_SHARED_LIBRARY)
  set(LIBRARY_TYPE SHARED)
else()
  set(LIBRARY_TYPE STATIC)
endif()
  
#Generate the shared or static library from the sources
add_library(${PROJECT_NAME} ${LIBRARY_TYPE} ${SOURCES} ${HEADERS})

## #################################################################
## Export header file
## #################################################################

generate_export_header(
  ${PROJECT_NAME} EXPORT_FILE_NAME 
  "${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}Export.h")

## #################################################################
## Configure file
## #################################################################

configure_file (
  "${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.h.in"
  "${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}Config.h")
