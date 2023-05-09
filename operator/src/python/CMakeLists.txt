# project(WEM-python)
cmake_minimum_required(VERSION 3.0)

set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

#require python3
find_package(PythonInterp 3 REQUIRED)

#need the lib directory for include
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include")

#need pybind11
find_package(pybind11 REQUIRED)

#need python include directory
include_directories(${PYTHON_INCLUDE_DIR})

#create pybind11 application
include_directories(${pybind11_INCLUDE_DIR})
include_directories(${TBB_INCLUDE_DIRS})
pybind11_add_module(pyOp MODULE pyOp.cpp)
target_link_libraries(pyOp PRIVATE Operator)


#install the libraries
install( TARGETS pyOp DESTINATION ${PROJECT_SOURCE_DIR}/lib/python)

#We will wrap the python/C++ interface with a strict python
install(FILES Operator.py DESTINATION ${PROJECT_SOURCE_DIR}/lib/python)
