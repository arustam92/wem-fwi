#include <axis.h>
#include <float2DReg.h>
#include <complex1DReg.h>
#include <complex2DReg.h>
#include <hypercube.h>
#include <cmath>
#include <Phshift.h>
#include <Scatter.h>
#include <Selector.h>

#include <cstring>

// typedef boost::multi_array<std::complex<float>, 3> float2D;
// typedef boost::multi_array<std::complex<float>, 2> complex2D;
using namespace SEP;
int main(int argc, char const *argv[])
{		
		int n2,n1;
		n2=500;
		n1=500;

		axis a1(n1,0,0.001);
		axis a2(n2,0,5);

		std::shared_ptr<hypercube> hyp (new hypercube(a1,a2));

		std::shared_ptr<complex1DReg> inp1 (new complex1DReg(a1));
		std::shared_ptr<complex1DReg> out1 (new complex1DReg(a1));
		std::shared_ptr<complex2DReg> input (new complex2DReg(hyp));
		std::shared_ptr<complex2DReg> output (new complex2DReg(hyp));
		std::shared_ptr<float2DReg> weight_m (new float2DReg(hyp));
		std::shared_ptr<float2DReg> weight_d = weight_m->clone();
		std::shared_ptr<float2DReg> w = weight_m->clone();


		{
			std::vector<int> loc(100);
			
			for (int i=0; i<100; ++i) loc[i] = (int)std::rand()%n1;
			std::cerr << "Selector" << std::endl;
			Selector select;
			select.setLocation(loc);
			select.setDomainRange(inp1,out1);
			select.dotTest();
		}

		{ // Injection
			Phshift ph(5,inp1->getHyper());
			ph.setFreq(5);
			ph.setSlow(0.002);
			ph.setDomainRange(inp1,out1);
			std::cerr << "Phshift" << std::endl;
			ph.dotTest();
		}

	return 0;
}