#include <Born_full.h>

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

	std::shared_ptr<genericRegFile> ptrDvel = ptrIo->getRegFile("dvel",usageIn);
	std::shared_ptr <hypercube> hypDvel = ptrDvel->getHyper();
	std::shared_ptr <float2DReg> d_slow(new float2DReg(hypDvel));
	ptrDvel->readFloatStream(d_slow);

	std::shared_ptr<genericRegFile> ptrWav = ptrIo->getRegFile("wave",usageIn);
	std::shared_ptr <hypercube> hypWav = ptrWav->getHyper();
	std::shared_ptr <float1DReg> wave(new float1DReg(hypWav));
	ptrWav->readFloatStream(wave);

	int nrec = ptrPar->getInt("nr");
	int nshots = ptrPar->getInt("ns");
	axis ax_s(nshots,0,ptrPar->getFloat("dsx"));
	axis ax_r(nrec,0,ptrPar->getFloat("drx"));
	axis ax_t(hypWav->getAxis(1).n,hypWav->getAxis(1).o,hypWav->getAxis(1).d);
	std::shared_ptr <float3DReg> data(new float3DReg(ax_t,ax_r,ax_s));

	Born_full born (wave,slow,ptrPar);

	if (!ptrPar->getString("wfld").empty()) {
		std::shared_ptr<complex3DReg> wfld_w(new complex3DReg(hypWav->getAxis(1),hypVel->getAxis(1),hypVel->getAxis(2)));
		born.setWfld(wfld_w);
	}

	born.forward(d_slow,data,0);

	// write output
	std::shared_ptr<genericRegFile> ptrOut = ptrIo->getRegFile("out",usageOut);
	ptrOut->setHyper(data->getHyper()); 	// Copy the pointer to the hypercube object
	ptrOut->writeDescription(); // write description on SEP header file
	ptrOut->writeFloatStream(data);

	if (!ptrPar->getString("wfld").empty()) {
		std::shared_ptr<float3DReg> wfl = born.getWfld();
		std::shared_ptr<genericRegFile> ptrWfld = ptrIo->getRegFile("wfld",usageOut);
		ptrWfld->setHyper(wfl->getHyper()); 	// Copy the pointer to the hypercube object
		ptrWfld->writeDescription(); // write description on SEP header file
		ptrWfld->writeFloatStream(wfl);
	}


	return 0;
}
