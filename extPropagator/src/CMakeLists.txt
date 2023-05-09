# add SCATE


set(SOURCES WEM_nl.cpp SplitStep.cpp SSF.cpp Born_full_cmplx.cpp
			LinOneWay.cpp Born_up.cpp Born_full.cpp Born_refl.cpp Born_down.cpp PSPI.cpp
			OneWay.cpp Phshift.cpp RefSampler.cpp Scatter.cpp WEM.cpp ic.cpp
			Selector.cpp )

set(INCLUDES OneStep.h SplitStep.h Reflect.h dReflect.h illumination.h Acquisition.h
			LinOneWay.h Born_up.h Born_full.h Born_refl.h Born_down.h
			ic.h LinOneWay.h
			OneWay.h Phshift.h RefSampler.h Scatter.h WEM.h Device.h
			Selector.h)

include_directories(${TBB_INCLUDE_DIRS})

if(SNAPSHOT)
	message("-- Building with SNAPSHOT")
	add_definitions(-DSNAPSHOT)
endif(SNAPSHOT)

add_library(Propagator ${SOURCES} ${INCLUDES})

target_link_libraries(Propagator Operator
						genericCpp sepVector hypercube buffers
						jsonCpp sep3d sep ${TBB_LIBRARIES} fftw3f m ${LIBVTUNE} dl
						)

install(TARGETS Propagator DESTINATION ${PROJECT_SOURCE_DIR}/lib)
