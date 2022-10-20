### universityDependencies.cmake --- 
## 
######################################################################
## 
### Commentary: Dependencies to external libraries
## 
######################################################################

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
# target_link_libraries(oxfordNewStudent mpi)

####################################################################
