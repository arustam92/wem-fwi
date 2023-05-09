#include <WEM.h>

#include <iostream>
#include "float1DReg.h"
#include "float3DReg.h"
#include "complex1DReg.h"
#include "float2DReg.h"
#include <complex3DReg.h>
#include "ioModes.h"
#include <PSPI.h>
#include <SSF.h>

using namespace SEP;

int main(int argc, char **argv)
{

	ioModes modes(argc,argv);
	std::shared_ptr <genericIO> ptrIo = modes.getDefaultIO();
	std::shared_ptr <paramObj> ptrPar = ptrIo->getParamObj();

	// prepare input

	std::shared_ptr<genericRegFile> ptrVel = ptrIo->getRegFile("in",usageIn);
	std::shared_ptr <hypercube> hypVel = ptrVel->getHyper();
	std::shared_ptr <complex2DReg> slow(new complex2DReg(hypVel));
	ptrVel->readComplexStream(slow);

	std::shared_ptr<genericRegFile> ptrWav = ptrIo->getRegFile("wave",usageIn);
	std::shared_ptr <hypercube> hypWav = ptrWav->getHyper();
	std::shared_ptr <float1DReg> wave(new float1DReg(hypWav));
	ptrWav->readFloatStream(wave);

	int nshots = ptrPar->getInt("ns");
	int nrec = ptrPar->getInt("nr");
	axis ax_s(nshots,0,ptrPar->getFloat("dsx"));
	axis ax_r(nrec,0,ptrPar->getFloat("drx"));
	axis ax_t(hypWav->getAxis(1).n,hypWav->getAxis(1).o,hypWav->getAxis(1).d);
	std::shared_ptr <float3DReg> data(new float3DReg(ax_t,ax_r,ax_s));

	int nref = ptrPar->getInt("nref",1);
	std::shared_ptr<RefSampler> ref (new RefSampler(slow,nref));

	{
		std::shared_ptr<complex1DReg> m (new complex1DReg(hypVel->getAxis(1)));
		std::shared_ptr<complex1DReg> d (new complex1DReg(hypVel->getAxis(1)));
		PSPI pspi(slow,ptrPar,ref);
		pspi.setDomainRange(m,d);
		pspi.setFreq(10.);
		pspi.setDepth(10);
		std::cerr << "PSPI" << std::endl;
		pspi.dotTest();

	}
	{
		std::shared_ptr<complex1DReg> m (new complex1DReg(hypVel->getAxis(1)));
		std::shared_ptr<complex1DReg> d (new complex1DReg(hypVel->getAxis(1)));
		SSF ssf(slow,ptrPar,ref);
		ssf.setDomainRange(m,d);
		ssf.setFreq(10.);
		ssf.setDepth(10);
		std::cerr << "SSF" << std::endl;
		ssf.dotTest();

	}

	{
		std::shared_ptr<complex2DReg> m (new complex2DReg(hypVel->getAxis(1),hypVel->getAxis(2)));
		std::shared_ptr<complex2DReg> d (new complex2DReg(hypVel->getAxis(1),hypVel->getAxis(2)));
		Down down (slow,ptrPar,ref);
		down.setFreq(10.);
		down.setDomainRange(m,d);
		std::cerr << "Down" << std::endl;
		down.dotTest();
	}

	{
		std::shared_ptr<complex2DReg> m (new complex2DReg(hypVel->getAxis(1),hypVel->getAxis(2)));
		std::shared_ptr<complex2DReg> d (new complex2DReg(hypVel->getAxis(1),hypVel->getAxis(2)));
		Up up (slow,ptrPar,ref);
		up.setFreq(10.);
		up.setDomainRange(m,d);
		std::cerr << "Up" << std::endl;
		up.dotTest();
	}

	{
		WEM wem (slow,ptrPar);
		wem.setDomainRange(wave,data);
		std::cerr << "WEM" << std::endl;
		wem.dotTest();
	}


	return 0;
}
