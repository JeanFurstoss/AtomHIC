### CMakeLists.txt --- 
##
######################################################################
## 
### Commentary: Main CMake file
##              Creation of library, tool, tests
## 
######################################################################

#Required version
cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#Project name
project(AtomHic)

#Choose generation of static or shared library
option(BUILD_SHARED_LIBRARY true)

#Enable or disable generation of tools
option(BUILD_TOOLS true)

#Enable or disable testing
option(BUILD_TESTS true)

# build type : Release or Debug
set(CMAKE_BUILD_TYPE "Debug")
#set(CMAKE_BUILD_TYPE "Release")

#Secondary cmake path (in order to not increase size of this main file)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
include(AtomHicBuild)
include(AtomHicDependencies)

#find_package(MPI)
#include_directories(${MPI_INCLUDE_PATH})
#set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

####################################################################
# Add directories containing sources, tools and tests and database #
####################################################################

#Database
add_subdirectory(data)

#Sources
add_subdirectory(src)

#Tools
if(BUILD_TOOLS)
  add_subdirectory(tools)
endif()

#Tests  
if(BUILD_TESTS)
  # Activate testing
  enable_testing()
  add_subdirectory(tests)
endif()

####################################################################
