#include "float1DReg.h"
#include "float2DReg.h"
#include "complex2DReg.h"
#include "ioModes.h"
#include <string>

#include "Signal.h"
#include <FFT1.h>


using namespace SEP;

int main(int argc, char **argv) {

	ioModes modes(argc,argv); 
	std::shared_ptr <genericIO> ptrIo = modes.getDefaultIO();
	std::shared_ptr <paramObj> ptrPar = ptrIo->getParamObj(); 

	float eps = ptrPar->getFloat("eps",0);
	bool wantwave = ptrPar->getBool("wantwave",true);	

	// prepare input
	std::shared_ptr<genericRegFile> ptrData = ptrIo->getRegFile("in",usageIn);	
	std::shared_ptr <hypercube> hypData = ptrData->getHyper();
	std::shared_ptr <Signal> data(new Signal(hypData));
	ptrData->readFloatStream(data);

	//prepare output

	std::shared_ptr<float2DReg> out = data->kolmogorov(eps, wantwave);


	// write output
	std::shared_ptr<genericRegFile> ptrOut = ptrIo->getRegFile("out",usageOut);
	ptrOut->setHyper(out->getHyper()); 	// Copy the pointer to the hypercube object
	ptrOut->writeDescription(); // write description on SEP header file 
	ptrOut->writeFloatStream(out);	

	return 0;
}