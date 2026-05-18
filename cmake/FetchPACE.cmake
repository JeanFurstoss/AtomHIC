function(LINK_PACE TARGET_NAME MODE)


if(PACE_INSTALL_FOLDER)
	message(STATUS "using PACE from local folder at ${PACE_INSTALL_FOLDER} for target ${TARGET_NAME}")

	target_include_directories(${TARGET_NAME} ${MODE} ${PACE_INSTALL_FOLDER})
else()
	if (TARGET pace)
	else()
		message(STATUS "fetching PACE from git")
		set(CMAKE_DISABLE_FIND_PACKAGE_yaml-cpp TRUE)
		include(FetchContent)
		FetchContent_Declare(
		  pace
		  GIT_REPOSITORY https://github.com/ICAMS/lammps-user-pace.git 
		  GIT_TAG        v.2025.12.3 
		)
		FetchContent_MakeAvailable(pace)
	endif()
	
	set_target_properties(pace PROPERTIES POSITION_INDEPENDENT_CODE ON)
	
	target_link_libraries(${TARGET_NAME} ${MODE} pace)
endif()

endfunction()

