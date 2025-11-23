### CMakeLists.txt --- 
##
######################################################################
## 
### Commentary: Main program
## 
######################################################################

#Project name
project(SampleGBShiftDoF)
 
#The file(GLOB...) allows for wildcard additions:
file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cpp")

#Generate executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

#Link executable to oxford library
target_link_libraries(${PROJECT_NAME} AtomHic)
