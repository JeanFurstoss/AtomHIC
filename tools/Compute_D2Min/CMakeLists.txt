### CMakeLists.txt --- 
##
######################################################################
## 
### Commentary: Main program
## 
######################################################################

#Project name
project(Compute_D2Min)
 
#The file(GLOB...) allows for wildcard additions:
file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cpp")

#Generate executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_link_libraries(${PROJECT_NAME} AtomHic)
