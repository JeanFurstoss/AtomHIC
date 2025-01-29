#**********************************************************************************
#*   CMakeLists.txt                                                               *
#**********************************************************************************
#* This file contains the building instructions of an AtomHic executable          *
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

file(GLOB SUBDIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)

# Ajouter chaque sous-dossier trouvé
foreach(subdir ${SUBDIRS})
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${subdir})
        add_subdirectory(${subdir})
    endif()
endforeach()

