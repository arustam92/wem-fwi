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

	// prepare input

	std::shared_ptr<genericRegFile> ptrData = ptrIo->getRegFile("in",usageIn);	
	std::shared_ptr <hypercube> hypData = ptrData->getHyper();
	std::shared_ptr <Signal> data(new Signal(hypData));
	ptrData->readFloatStream(data);

	std::shared_ptr<genericRegFile> ptrFil = ptrIo->getRegFile("filter",usageIn);	
	std::shared_ptr <hypercube> hypFil = ptrFil->getHyper();
	std::shared_ptr <float2DReg> filter(new float2DReg(hypFil));
	ptrFil->readFloatStream(filter);

	//prepare output
	std::shared_ptr <Signal> out = data->conv(filter);

	// write output
	std::shared_ptr<genericRegFile> ptrOut = ptrIo->getRegFile("out",usageOut);
	ptrOut->setHyper(out->getHyper()); 	// Copy the pointer to the hypercube object
	ptrOut->writeDescription(); // write description on SEP header file 
	ptrOut->writeFloatStream(out);	

	return 0;
}