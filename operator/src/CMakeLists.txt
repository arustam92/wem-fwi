
set(SOURCES FFT1_c.cpp FFT1_r.cpp Taper.cpp)

set(INCLUDES
			StackAndSpread.h
			Spline2D.h Pad.h Interpolation3D.h Spline3D.h Hilbert.h
			LeakyIntegration.h LanczosInterpolation3D.h LanczosInterpolation2D.h
			FFT1.h Taper.h Weight.h Operator.h)

set(ISPC_SOURCES lanczos_kernel.ispc)
set(ISPC_OBJS)
set(ISPC_INCLUDES)

foreach(src ${ISPC_SOURCES})
	string(REPLACE "ispc" "o" obj ${src})
	string(REPLACE "ispc" "h" inc ${src})
	add_custom_command(OUTPUT ${PROJECT_SOURCE_DIR}/build/src/${obj}
	COMMAND ispc --target=avx2-i32x16 ${CMAKE_CURRENT_SOURCE_DIR}/${src} -O3 --opt=fast-math
	-o ${PROJECT_SOURCE_DIR}/build/src/${obj}
	-h ${CMAKE_CURRENT_SOURCE_DIR}/${inc}
	DEPENDS ${src})
	list(APPEND ISPC_OBJS ${PROJECT_SOURCE_DIR}/build/src/${obj})
	list(APPEND ISPC_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/${inc})
endforeach()

include_directories(${TBB_INCLUDE_DIRS})
add_library(Operator ${ISPC_OBJS} ${ISPC_INCLUDES} ${SOURCES} ${INCLUDES})

target_link_libraries(Operator
						genericCpp sepVector hypercube buffers
						jsonCpp sep3d sep ${TBB_LIBRARIES} fftw3f m
						)
install(TARGETS Operator DESTINATION ${PROJECT_SOURCE_DIR}/lib)
