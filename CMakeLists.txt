cmake_minimum_required(VERSION 3.12)


project(templateMatchingAlignment VERSION 0.1.0)

#set(FORCE_STATIC 1)
if(DEFINED FORCE_STATIC )
	message(STATUS "NOTICE: **** Static build for templateMatchingAlignment *****")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
	add_definitions(-DIO_DATA_STATIC -DLXGEO_FACTORY_STATIC)
else()
	message(STATUS "NOTICE: **** Dynamic build templateMatchingAlignment *****")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
endif()

set(CMAKE_BINARY_DIR "../out/${CMAKE_BUILD_TYPE}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)

add_subdirectory("$ENV{LX_GEO_REPOS_ROOT}/LxGeoCppLibs/libs/shared" LxGeoFactoryShared)
add_subdirectory("$ENV{LX_GEO_REPOS_ROOT}/LxGeoCppLibs/libs/IO_DATA" IO_DATA)
add_subdirectory(app)