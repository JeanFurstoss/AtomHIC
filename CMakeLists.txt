#**********************************************************************************
#*   CMakeLists.txt                                                               *
#**********************************************************************************
#* This file contains the building instructions of the AtomHic library            *
#**********************************************************************************
#* (C) Jan 2025 - Jean Furstoss                                                   *
#*     Université de Poitiers, Institut PPRIME                                    *
#*     UPR CNRS 3346, 86360 Chasseuneuil-du-Poitou, France                        *
#*     jean.furstoss@univ-poitiers.fr                                             *
#* Last modification: J. Furstoss - 28 Janv 2025                                  *
#**********************************************************************************
#* This program is free software: you can redistribute it and/or modify           *
#* it under the terms of the GNU General Public License as published by           *
#* the Free Software Foundation, either version 3 of the License, or              *
#* (at your option) any later version.                                            *
#*                                                                                *
#* This program is distributed in the hope that it will be useful,                *
#* but WITHOUT ANY WARRANTY; without even the implied warranty of                 *
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  *
#* GNU General Public License for more details.                                   *
#*                                                                                *
#* You should have received a copy of the GNU General Public License              *
#* along with this program.  If not, see <http://www.gnu.org/licenses/>.          *
#**********************************************************************************

#Required version
cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/MakeLibrary.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FetchEigen.cmake)

#Project name
project(AtomHic)

#Choose generation of static or shared library
option(BUILD_SHARED_LIBRARY "Build as shared library" ON)

#Enable or disable generation of tools
option(BUILD_TOOLS "Enable or disable generation of tools" ON)

#Enable or disable testing
option(BUILD_TESTS "Enable or disable testing" ON)

# build type : Release or Debug
#set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_BUILD_TYPE "Release")

#Secondary cmake path (in order to not increase size of this main file)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
include(AtomHicBuild)
include(AtomHicDependencies)

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
    message(STATUS "Tools generation is enabled.")
    add_subdirectory(tools)
else()
    message(WARNING "Tools generation is disabled. Use -DBUILD_TOOLS=ON to enable it.")
endif()

#Tests  
if(BUILD_TESTS)
    # Activate testing
    enable_testing()
    add_subdirectory(tests)
endif()
####################################################################
