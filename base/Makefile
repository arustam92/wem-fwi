TARGET = libBase.a
.PHONY: default
	default: $(TARGET);

SRCDIR:=./src
BIN:=./bin
BUILDDIR:=./build

CC:=g++ -fopenmp -std=c++11
SRCEXT:=cpp

REPO := /homes/sep/arustam/repository.dev

SEPVECINC:=-I/sep/gbarnier/genericIO/include -I/opt/SEP/SEP8.1/include -I${REPO}/complex/include
SEPVECLIBS:=-L/sep/gbarnier/genericIO/lib -L/opt/SEP/SEP8.1/lib -L${REPO}/complex/lib  -lgenericCpp -lsepVector -lhypercube -lsep3d -lsep  -lcomplexHyper

MYINC:=-I include
MYLIB:=lib

FFTWINC:=-I/opt/FFTW/include
FFTWLIB:=-L/opt/FFTW/lib -lfftw3f -lm

SOURCES:=$(shell find $(SRCDIR) -type f -name "*.$(SRCEXT)")
OBJECTS:=$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))


$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(SEPVECINC) $(MYINC) -c -o $@ $<"; $(CC) $(SEPVECINC) $(MYINC) $(FFTWINC) $< $(SEPVECLIBS) $(FFTWLIB) -c -o$@

${TARGET}: ${OBJECTS}
	@mkdir -p ${MYLIB}
	@echo "Creating $@"...; 
	ar crf ${MYLIB}/$@ $?

main.x: main.cpp
	$(CC) $(SEPVECINC) $(MYINC) $(FFTWINC) $< $(SEPVECLIBS) $(FFTWLIB) -L./lib -lBase -o$@
	
