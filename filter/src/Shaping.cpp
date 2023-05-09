#include <iostream>
#include "float1DReg.h"
#include "float2DReg.h"
#include "ioModes.h"
#include <string>

#include "Signal.h"


using namespace SEP;

int main(int argc, char **argv) {

	ioModes modes(argc,argv); 
	std::shared_ptr <genericIO> ptrIo = modes.getDefaultIO();
	std::shared_ptr <paramObj> ptrPar = ptrIo->getParamObj(); 

	int poslag = ptrPar->getInt("poslag");
	int neglag = ptrPar->getInt("neglag",0);
	float eps = ptrPar->getFloat("eps",0);
	int gap = ptrPar->getInt("gap",1);		

	// prepare input

	std::shared_ptr<genericRegFile> ptrData = ptrIo->getRegFile("in",usageIn);	
	std::shared_ptr <hypercube> hypData = ptrData->getHyper();
	std::shared_ptr <Signal> data(new Signal(hypData));
	ptrData->readFloatStream(data);

	std::shared_ptr<genericRegFile> ptrWav = ptrIo->getRegFile("wave",usageIn);	
	std::shared_ptr <hypercube> hypWav = ptrWav->getHyper();
	std::shared_ptr <Signal> wave(new Signal(hypWav));
	ptrWav->readFloatStream(wave);

	//prepare output
	std::vector<axis> ax = hypData->getAxes();
	axis f1(poslag+neglag+1, -neglag*ax[0].d, ax[0].d);
	std::shared_ptr<hypercube> hypFilter (new hypercube(f1, ax[1]));	
	std::shared_ptr <Signal> filter (new Signal(hypFilter));

	// actual computations
	double start = omp_get_wtime();
	std::shared_ptr <Signal> acorr = data->xcorr(data);
	std::shared_ptr <Signal> corr = wave->xcorr(data);
	filter->levinson(acorr,corr,eps);
	double end = omp_get_wtime();

	std::cerr << "Time "<< end-start << std::endl;

	// write output
	std::shared_ptr<genericRegFile> ptrOut = ptrIo->getRegFile("out",usageOut);
	ptrOut->setHyper(filter->getHyper()); 	// Copy the pointer to the hypercube object
	ptrOut->writeDescription(); // write description on SEP header file 
	ptrOut->writeFloatStream(filter);	

	return 0;
}