
set(SOURCES
	ConformalMap.cpp ComplexLog.cpp
)

set(INCLUDES
	MapChain.h Rotation.h
	ConformalMap.h ComplexLog.h Identity.h ShiftAndScale.h
)

include_directories(${TBB_INCLUDE_DIRS})
add_library(ConformalMap ${SOURCES} ${INCLUDES})

target_link_libraries(ConformalMap
						genericCpp sepVector hypercube 
						)

install(TARGETS ConformalMap DESTINATION ${PROJECT_SOURCE_DIR}/lib)
