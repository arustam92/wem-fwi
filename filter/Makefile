S = ./src
B = ./bin
T = ./test


UCXINCLUDES=-I/sep/gbarnier/genericIO/include -I/opt/SEP/SEP8.1/include
UCXLIBS=-L/sep/gbarnier/genericIO/lib -L/opt/SEP/SEP8.1/lib -lgenericCpp -lsepVector -lhypercube -lsep3d -lsep -L/homes/sep/arustam/repository/complex/lib -lcomplexHyper

MYINC = -I/homes/sep/arustam/repository/base/include -I/homes/sep/arustam/repository/complex/include -I/opt/FFTW/include
MYLIB = -L/homes/sep/arustam/repository/base/lib -lBase -L/opt/FFTW/lib -lfftw3f -lm



$(B)/%.x: $(S)/%.cpp
	g++ -std=c++11 -fopenmp ${UCXINCLUDES} ${MYINC} $< /homes/sep/arustam/repository/base/build/Signal.o ${UCXLIBS} ${MYLIB} -o$@	


include $(T)/shaping.m
include $(T)/kolmo.m


	
