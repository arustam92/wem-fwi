#include <BornFull.h>

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
	std::shared_ptr <complex2DReg> slow(new complex2DReg(hypVel));
	ptrVel->readComplexStream(slow);

	std::shared_ptr<genericRegFile> ptrDvel = ptrIo->getRegFile("dvel",usageIn);
	std::shared_ptr <hypercube> hypDvel = ptrDvel->getHyper();
	std::shared_ptr <complex2DReg> d_slow(new complex2DReg(hypDvel));
	ptrDvel->readComplexStream(d_slow);

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

	std::shared_ptr<complex2DReg> m (new complex2DReg(hypVel));
	std::shared_ptr<float3DReg> d (new float3DReg(data->getHyper()));

	{
		dReflect dr(slow);
		std::shared_ptr<complex2DReg> m = slow->clone();
		std::shared_ptr<complex2DReg> d (new complex2DReg(hypVel->getAxis(1),hypVel->getAxis(2)));
		dr.setDomainRange(m,d);
		std::cerr << "dReflect" << '\n';
		dr.dotTest();
	}

{
	int nref = ptrPar->getInt("nref",1);
	std::shared_ptr<RefSampler> ref(new RefSampler(slow,nref));

	Scatter sc (slow,1);
	std::shared_ptr<complex1DReg> scm(new complex1DReg(slow->getHyper()->getAxis(1)));
	std::shared_ptr<complex1DReg> scd(new complex1DReg(slow->getHyper()->getAxis(1)));
	sc.setDomainRange(scm,scd);
	sc.setDepth(10);
	sc.setFreq(50.);
	std::cerr << "Scatter" << std::endl;
	sc.dotTest();

		std::shared_ptr<complex2DReg> bg_wfld(new complex2DReg(slow->getHyper()->getAxis(1),slow->getHyper()->getAxis(2)));
		bg_wfld->random();

		std::shared_ptr<complex2DReg> m(new complex2DReg(slow->getHyper()->getAxis(1),slow->getHyper()->getAxis(2)));
		std::shared_ptr<complex2DReg> d(new complex2DReg(slow->getHyper()->getAxis(1),slow->getHyper()->getAxis(2)));
		IC ic(bg_wfld);
		ic.setDomainRange(m,d);
		std::cerr << "IC:" << std::endl;
		ic.dotTest();

}

{	int nref = ptrPar->getInt("nref",1);
	std::shared_ptr<RefSampler> ref(new RefSampler(slow,nref));
	std::shared_ptr<Down> down(new Down(slow,ptrPar,ref));
	down->setFreq(50.);
	std::shared_ptr<complex2DReg> bg_wfld(new complex2DReg(slow->getHyper()->getAxis(1),slow->getHyper()->getAxis(2)));
	bg_wfld->random();
	LinDown linDown (slow,ptrPar,bg_wfld,down);
	linDown.setFreq(50.);
	std::shared_ptr<complex2DReg> m(new complex2DReg(slow->getHyper()));
	std::shared_ptr<complex2DReg> d(new complex2DReg(hypVel->getAxis(1),hypVel->getAxis(2)));
	linDown.setDomainRange(m,d);
	std::cerr << "LinDown" << std::endl;
	linDown.dotTest();
}

{	int nref = ptrPar->getInt("nref",1);
	std::shared_ptr<RefSampler> ref(new RefSampler(slow,nref));
	std::shared_ptr<Up> up(new Up(slow,ptrPar,ref));
	up->setFreq(50.);
	std::shared_ptr<complex2DReg> bg_wfld(new complex2DReg(slow->getHyper()->getAxis(1),slow->getHyper()->getAxis(2)));
	bg_wfld->random();
	LinUp linUp (slow,ptrPar,bg_wfld,up);
	linUp.setFreq(50.);
	std::shared_ptr<complex2DReg> m(new complex2DReg(slow->getHyper()));
	std::shared_ptr<complex2DReg> d(new complex2DReg(hypVel->getAxis(1),hypVel->getAxis(2)));
	linUp.setDomainRange(m,d);
	std::cerr << "LinUp" << std::endl;
	linUp.dotTest();
}

	BornFull born (wave,slow,ptrPar);
	born.setDomainRange(m,d);
	std::cerr << "Born" << std::endl;
	born.dotTest();

	return 0;
}
