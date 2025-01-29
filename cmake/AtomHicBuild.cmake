#**********************************************************************************
#*   CMakeLists.txt                                                               *
#**********************************************************************************
#* File containing functions for generation of export macros for libraries	  *
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
#* What is still needed to do here:                                               *
#*	- manage windows compîlation			                          *
#**********************************************************************************

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
