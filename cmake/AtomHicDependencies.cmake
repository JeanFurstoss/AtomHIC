#**********************************************************************************
#*   CMakeLists.txt                                                               *
#**********************************************************************************
#* Dependencies to external libraries	  					  *
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

## Python is not necessary in this case, it's only to show how 
## to use "find_package" function
find_package(PythonLibs QUIET)

if(PYTHONLIBS_FOUND)
  # Include python directories (headers)
  include_directories(${PYTHON_INCLUDE_DIRS})
  # To link an executable or a library to this library (lib, dll, so, a)
  # target_link_libraries(oxfordNewStudent ${PYTHON_LIBRARIES})
endif(PYTHONLIBS_FOUND)

## Example with a library specified by the user
## Assume that user sets environment variable for his local version of MPI
# MPI_INCLUDE_DIR=/home/toto/tools/openmpi/include
# MPI_LIBRARY_DIR=/home/toto/tools/openmpi/lib
## Specify include directory
# include_directories($ENV{MPI_INCLUDE_DIR})
## Specify library path
# link_directories($ENV{MPI_LIBRARY_DIR})
## To link an executable or a library to this library
# target_link_libraries(AtomHic mpi)

####################################################################
