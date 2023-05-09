include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(SOURCES complexHyper.cc complex1DReg.cc complex2DReg.cc complex3DReg.cc)

set(INCLUDES complexHyper.h complex1DReg.h complex2DReg.h complex3DReg.h)


add_library(complexHyper ${SOURCES})
install(TARGETS complexHyper DESTINATION lib)