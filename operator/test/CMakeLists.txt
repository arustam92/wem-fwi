
# set(SOURCES checkInj.cpp checkSort.cpp)
include_directories(${TBB_INCLUDE_DIRS})


foreach(item ${SOURCES})
	string( REPLACE ".cpp" ".x" exec ${item} )
    add_executable( ${exec} ${item} )
    # Make sure YourLib is linked to each app
    target_link_libraries( ${exec} Operator genericCpp sepVector hypercube ${TBB_LIBRARIES} fftw3f m)
    install(TARGETS ${exec} DESTINATION ${PROJECT_SOURCE_DIR}/bin)
endforeach(item)
