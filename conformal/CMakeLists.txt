#Requied CMake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

#We need to add to how to find TBB which sits in the cmake dir
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

#Name of our cmake project
project(ConformalMap LANGUAGES CXX)

#Require c++111
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# enable -fPic and include the path to the cur dir
set(CMAKE_POSITION_INDEPENDENT_CODE True)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

#Compiler flags
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Ofast")

#Required packages
find_package( Boost 1.40 COMPONENTS program_options unit_test_framework REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

#include the giee directories when compiling and linking
find_package(SEPlib REQUIRED)
if(SEPlib_FOUND)
	message("-- Found SEPlib")
	include_directories(${SEPlib_INCLUDE_DIRS})
	LINK_DIRECTORIES(${SEPlib_LIBRARY_DIRS})
endif(SEPlib_FOUND)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

include_directories(src)

#add the rules for building the C++ code
add_subdirectory(src)
add_subdirectory(src/python)

link_directories(${PROJECT_SOURCE_DIR}/lib)
# add_subdirectory(test)

#add the rules for python binding
# add_subdirectory(python)

















# project(sepVector)
# cmake_minimum_required(VERSION 2.8)
# set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

# # configure a header file to pass some of the CMake settings
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
# find_package(Boost REQUIRED)


# find_package(NumPy  REQUIRED)
# set(sepVector_MAJOR_VERSION 0)
# set(sepVector_MINOR_VERSION 2)
# set(sepVector_VERSION_STRING
#   ${sepVector_MAJOR_VERSION}.${sepVector_MINOR_VERSION})

#  SET(sepVector_DIRNAME_lib     "lib")
#  SET(sepVector_DIRNAME_include "include")

# include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# FIND_PACKAGE(hypercube)
# if(hypercube_FOUND)
#   message("Found hypercube")
#   INCLUDE_DIRECTORIES(${hypercube_DIRNAME_include})
#   link_directories(hypercube_DIRNAME_lib)
# else(hypercube_FOUND)
#    message("Unable to find hypercube")

#     # we have a submodule setup for foo, assume it is under external/foo
#     # now we need to clone this submodule
#     execute_process(COMMAND git submodule update --init --recursive -- external/hypercube
#                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
#     message("before add subdirectory")
#     add_subdirectory(external/hypercube)
#     INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/external/hypercube/lib)

# endif(hypercube_FOUND)





# # add the binary tree to the search path for include files
# # so that we will find TutorialConfig.h
# include_directories("${PROJECT_BINARY_DIR}")


# include(CheckCXXCompilerFlag)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
# CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
# CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
# if(COMPILER_SUPPORTS_CXX11)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# else()
#         message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
# endif()

# enable_language(C CXX)






#  SET(sepVector_DIRNAME_lib     "lib")
#  SET(sepVector_DIRNAME_include "include")
#  SET(sepVector_DIRNAME_cmake "cmake")
#   INCLUDE(CMakePackageConfigHelpers)

#   CONFIGURE_PACKAGE_CONFIG_FILE(
#   "${CMAKE_CURRENT_SOURCE_DIR}/sepVectorConfig.cmake.in"
#   "${CMAKE_CURRENT_BINARY_DIR}/sepVectorConfig.cmake"
#   INSTALL_DESTINATION
#     ${sepVector_DIRNAME_cmake}
#   PATH_VARS
#     sepVector_DIRNAME_include
#     sepVector_DIRNAME_lib
#   )

#   INSTALL(
#   FILES
#     ${CMAKE_CURRENT_BINARY_DIR}/sepVectorConfig.cmake
#   DESTINATION
#     cmake
#   )


#   INSTALL(
#   EXPORT
#     sepVectorExports
#   DESTINATION
#    cmake
#   )

# add_subdirectory(lib)
# add_subdirectory(python)
# add_subdirectory(env-modules)

# include(InstallRequiredSystemLibraries)
# set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
# set(CPACK_PACKAGE_VERSION_MAJOR "${sepVector_MAJOR_VERSION}")
# set(CPACK_PACKAGE_VERSION_MINOR "${sepVector_MINOR_VERSION}")
# set(CPACK_PACKAGE_VERSION_PATCH "0")
# set(CPACK_PACKAGE_NAME "sepVector")
# set(CPACK_PACKAGE_CONTACT "Robert Clapp")
# include(CPack)
