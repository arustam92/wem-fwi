# - Find FFTW
# Find the native FFTW includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (FFTW_INCLUDES)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h HINTS ${FFTW_ROOT}/include)
find_path (FFTWF_INCLUDES fftw3.f03 HINTS ${FFTW_ROOT}/include )

find_library (FFTW_LIB NAMES libfftw3.a HINTS ${FFTW_ROOT}/lib)

find_library (FFTWF_LIB NAMES libfftw3f.a HINTS ${FFTW_ROOT}/lib)

if(NOT FFTW_LIB) 
  MESSAGE("COULD NOT FIND FFTW, TRY SETTING FFTW_ROOT")
endif(NOT FFTW_LIB)
#
set(FFTW_LIBRARIES ${FFTWF_LIB} ${FFTW_LIB})
set(FFTWF_LIBRARIES ${FFTWF_LIB} ${FFTW_LIB})

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTWF_LIBRARIES FFTW_LIBRARIES FFTWF_INCLUDES FFTW_INCLUDES)

mark_as_advanced (FFTWF_LIBRARIES FFTW_LIBRARIES FFTWF_INCLUDES FFTW_INCLUDES )
