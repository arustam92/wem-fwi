#include <extWEM.h>

#include <iostream>
#include "float1DReg.h"
#include "float3DReg.h"
#include "complex1DReg.h"
#include "float2DReg.h"
#include <complex3DReg.h>
#include "ioModes.h"

#include <cstring>

using namespace SEP;
int main(int argc, char *argv[])
{
		int n3,n2,n1,nt;
		n2=100;
		n1=100;
		n3 = 10;
		nt = 100;

		axis a1(n1,0,0.01);
		axis a2(n2,0,0.05);
		axis a3(n3,0,1);

		ioModes modes(argc,argv);
		std::shared_ptr <genericIO> ptrIo = modes.getDefaultIO();
		std::shared_ptr <paramObj> ptrPar = ptrIo->getParamObj();

		auto slow = std::make_shared<complex3DReg>(a1, a2, a3);
		slow->set(1.f);

		std::vector<std::shared_ptr<complex3DReg>> model;
		model.push_back(slow);
		model.push_back(slow);

		auto wave = std::make_shared<float1DReg>(nt);

		int nshots = ptrPar->getInt("ns");
		int nrec = ptrPar->getInt("nr");
		axis ax_s(nshots,0,ptrPar->getFloat("dsx"));
		axis ax_r(nrec,0,ptrPar->getFloat("drx"));
		axis ax_t(wave->getHyper()->getAxis(1).n,wave->getHyper()->getAxis(1).o,wave->getHyper()->getAxis(1).d);
		std::shared_ptr <float3DReg> data(new float3DReg(ax_t,ax_r,ax_s));

		extWEM extwem (wave, slow->getHyper(),ptrPar);
		std::cerr << "WEM" << std::endl;
		extwem.forward(model, data, false);

	return 0;
}
