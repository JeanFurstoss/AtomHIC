### CMakeLists.txt --- 
## 
######################################################################
## 
### Commentary: Create "doc" target to build documentation with doxygen
## 
######################################################################

find_package(Doxygen QUIET)

if(DOXYGEN_FOUND)

  # Copy doxygen configuration file
  configure_file(doxygen.cfg.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen.cfg @ONLY)

  add_custom_target(doc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${DOXYGEN_EXECUTABLE} doxygen.cfg
    )

else(DOXYGEN_FOUND)

  message(WARNING "Doxygen must be installed to generate documentation")

endif(DOXYGEN_FOUND)
