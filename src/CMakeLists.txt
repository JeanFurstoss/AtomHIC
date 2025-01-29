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
