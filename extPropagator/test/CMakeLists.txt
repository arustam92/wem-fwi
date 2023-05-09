
set(SOURCES dotBorn.cpp
			check_adjBorn.cpp checkBorn.cpp check_adjWEM.cpp
			checkWEM.cpp dotPh.cpp dotWEM.cpp)

include_directories(${TBB_INCLUDE_DIRS})

foreach(item ${SOURCES})
	string( REPLACE ".cpp" ".x" exec ${item} )
    add_executable( ${exec} ${item} )
    # Make sure YourLib is linked to each app
    target_link_libraries( ${exec} Propagator Operator genericCpp sepVector hypercube buffers jsonCpp sep3d sep ${TBB_LIBRARIES} fftw3f m)
    install(TARGETS ${exec} DESTINATION ${PROJECT_SOURCE_DIR}/bin)
endforeach(item)
