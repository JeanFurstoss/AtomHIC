### CMakeLists.txt --- 
######################################################################
## 
### Commentary: 
## 
######################################################################

project(compare)

######################################################################
## Sources
######################################################################

# Sources
set(${PROJECT_NAME}_SOURCES main.cc compareExact.cc compareSum.cc)

# Include
set(${PROJECT_NAME}_HEADERS compare.h)

## #################################################################
## Build rules
## #################################################################

add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS})

target_link_libraries(${PROJECT_NAME} AtomHic)

