### universityBuild.cmake --- 
## 
####################################################################
## 
### Commentary: 
## 
####################################################################

## #################################################################
## File containing functions for generation of 
## export macros for libraries
## #################################################################

include(GenerateExportHeader)

## #################################################################
## Hide variables
## #################################################################

mark_as_advanced(CMAKE_BUILD_TYPE)
mark_as_advanced(CMAKE_CONFIGURATION_TYPES)
mark_as_advanced(CMAKE_INSTALL_PREFIX)

## #################################################################
## Include directories
## #################################################################

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/include)

## #################################################################
## Configure libraries and executables path
## #################################################################

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Only on windows
if(WIN32)
  ## for multi config types (MSVC based)
  foreach(CONFIG_TYPE ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG_TYPE} CONFIG_TYPE_UPPERCASE)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG_TYPE_UPPERCASE} 
        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG_TYPE_UPPERCASE}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  endforeach()
endif()
