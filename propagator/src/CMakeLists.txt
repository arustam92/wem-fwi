
set(SOURCES
	Injection.cpp
	ext_born_refl.cpp ext_born_tomo.cpp BornRefl.cpp BornTomo.cpp
	ext_born_prop.cpp ext_wem_prop.cpp extWEM.cpp extBornFull.cpp
	 		born_prop.cpp  wem_prop.cpp
	 		WEMnl.cpp SplitStep.cpp SSF.cpp
			LinOneWay.cpp Born_up.cpp BornFull.cpp Born_refl.cpp Born_down.cpp PSPI.cpp
			OneWay.cpp Phshift_cmplx.cpp RefSampler.cpp Scatter.cpp ic.cpp
			Selector.cpp )

set(INCLUDES
	Injection.h
	PropParam.h
	extBornFull.h extWEM.h
			 BornTomo.h  BornRefl.h OneStep.h SplitStep.h
			Reflect.h dReflect.h illumination.h Acquisition.h
			LinOneWay.h Born_up.h BornFull.h Born_refl.h Born_down.h
			ic.h LinOneWay.h
			OneWay.h Phshift.h RefSampler.h Scatter.h WEM.h Device.h
			Selector.h)

set(ISPC_SOURCES phase_shift.ispc)
set(ISPC_OBJS)
set(ISPC_INCLUDES)

foreach(src ${ISPC_SOURCES})
	string(REPLACE "ispc" "o" obj ${src})
	string(REPLACE "ispc" "h" inc ${src})
	add_custom_command(OUTPUT ${PROJECT_SOURCE_DIR}/build/src/${obj}
	COMMAND ispc --target=avx2-i32x16 --pic
	${CMAKE_CURRENT_SOURCE_DIR}/${src} -O3 --opt=fast-math
	-o ${PROJECT_SOURCE_DIR}/build/src/${obj}
	-h ${CMAKE_CURRENT_SOURCE_DIR}/${inc}
	DEPENDS ${src})
	list(APPEND ISPC_OBJS ${PROJECT_SOURCE_DIR}/build/src/${obj})
	list(APPEND ISPC_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/${inc})
endforeach()

include_directories(${TBB_INCLUDE_DIRS})

add_library(Propagator ${ISPC_OBJS} ${SOURCES} ${INCLUDES})

target_link_libraries(Propagator Operator ConformalMap sepVector genericCpp ${TBB_LIBRARIES} fftw3f)

install(TARGETS Propagator DESTINATION ${PROJECT_SOURCE_DIR}/lib)
