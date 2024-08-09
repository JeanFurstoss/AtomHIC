### CMakeLists.txt --- 
##
######################################################################
## 
### Commentary: Build of AtomHic library
## 
######################################################################
### GMM
#include(FetchContent)
#FetchContent_Declare(
#  gmm
#  GIT_REPOSITORY https://github.com/andreacasalino/Gaussian-Mixture-Model.git
#  #  GIT_TAG        58b3f207cf6df73a0b7cdc8bde933b06d20ab502
#)
#FetchContent_MakeAvailable(gmm)


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

### GMM
#target_link_libraries(${PROJECT_NAME} PUBLIC
#GaussianMixtureModels
#)

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
