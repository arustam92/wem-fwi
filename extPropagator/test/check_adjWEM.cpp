#include <WEM.h>

#include <iostream>
#include "float1DReg.h"
#include "float3DReg.h"
#include "complex1DReg.h"
#include "float2DReg.h"
#include <complex3DReg.h>
#include "ioModes.h"


using namespace SEP;

int main(int argc, char **argv)
{

	ioModes modes(argc,argv);
	std::shared_ptr <genericIO> ptrIo = modes.getDefaultIO();
	std::shared_ptr <paramObj> ptrPar = ptrIo->getParamObj();

	// prepare input

	std::shared_ptr<genericRegFile> ptrVel = ptrIo->getRegFile("in",usageIn);
	std::shared_ptr <hypercube> hypVel = ptrVel->getHyper();
	std::shared_ptr <float2DReg> slow(new float2DReg(hypVel));
	ptrVel->readFloatStream(slow);

	std::shared_ptr<genericRegFile> ptrDat = ptrIo->getRegFile("data",usageIn);
	std::shared_ptr <hypercube> hypDat = ptrDat->getHyper();
	std::shared_ptr <float3DReg> data(new float3DReg(hypDat));
	ptrDat->readFloatStream(data);

	// int nrec = ptrPar->getInt("n_src");
	// axis ax(nrec,0,1);
	axis ax_t = hypDat->getAxis(1);
	std::shared_ptr <float1DReg> wave(new float1DReg(ax_t));

	WEM wem (slow,ptrPar);

	wem.adjoint(wave,data,0);

	// write output
	std::shared_ptr<genericRegFile> ptrOut = ptrIo->getRegFile("out",usageOut);
	ptrOut->setHyper(wave->getHyper()); 	// Copy the pointer to the hypercube object
	ptrOut->writeDescription(); // write description on SEP header file
	ptrOut->writeFloatStream(wave);



	return 0;
}
