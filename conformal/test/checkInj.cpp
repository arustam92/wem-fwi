#include <axis.h>
#include <float1DReg.h>
#include <float2DReg.h>
#include <complex1DReg.h>
#include <complex2DReg.h>
#include <hypercube.h>
#include <cmath>
#include <Injection.h>
#include <Weight.h>
#include <FFT1.h>
#include <Taper.h>
#include <numeric>
#include <cstring>

// typedef boost::multi_array<std::complex<float>, 3> float2D;
// typedef boost::multi_array<std::complex<float>, 2> complex2D;
using namespace SEP;
int main(int argc, char const *argv[])
{
		int n2,n1,n3;
		n2=500;
		n1=500;
		n3=1;

		axis a1(n1,0,0.001);
		axis a2(n2,0,5);
		axis a3(n3,0,5);

		std::shared_ptr<hypercube> hyp (new hypercube(a1,a2));
		std::shared_ptr<hypercube> hyp3 (new hypercube(a1,a2,a3));

		std::shared_ptr<complex1DReg> inp1 (new complex1DReg(a1));
		std::shared_ptr<complex1DReg> out1 (new complex1DReg(a1));
		std::shared_ptr<complex2DReg> input (new complex2DReg(hyp));
		std::shared_ptr<complex2DReg> output (new complex2DReg(hyp));
		std::shared_ptr<complex3DReg> input3 (new complex3DReg(hyp3));
		std::shared_ptr<float2DReg> weight_m (new float2DReg(hyp));
		std::shared_ptr<float2DReg> weight_d = weight_m->clone();
		// std::shared_ptr<float2DReg> w = weight_m->clone();

		{ // Injection
			std::vector<int> iz (n1,0);
			std::vector<int> ix(n2);
			std::iota(ix.begin(),ix.end(),0);
			Injection inj (iz,ix);
			inj.setStep(1);
			inj.setShot(0);
			inj.Operator<complex3DReg,complex2DReg>::setDomainRange(input3,output);
			std::cerr << "Injection<complex3DReg,complex2DReg>" << std::endl;
			inj.Operator<complex3DReg,complex2DReg>::dotTest();
		}
		{ // Injection
			std::vector<int> iz (n1,0);
			std::vector<int> ix(n2);
			std::iota(ix.begin(),ix.end(),0);
			Injection inj (iz,ix);
			inj.setStep(1);
			inj.setShot(0);
			inj.Operator<complex1DReg,complex2DReg>::setDomainRange(inp1,output);
			std::cerr << "Injection<complex1DReg,complex2DReg>" << std::endl;
			inj.Operator<complex1DReg,complex2DReg>::dotTest();
		}


		// {
		// 	// Weight
		// 	w->random();
		// 	Weight weight_op (w);
		// 	weight_op.Operator<float2DReg,float2DReg>::setDomainRange(weight_m,weight_d);
		// 	std::cerr << "Weight<float2DReg,float2DReg>" << std::endl;
		// 	weight_op.Operator<float2DReg,float2DReg>::dotTest();

		// 	weight_op.Operator<complex2DReg,complex2DReg>::setDomainRange(input,output);
		// 	std::cerr << "Weight<complex2DReg,complex2DReg>" << std::endl;
		// 	weight_op.Operator<complex2DReg,complex2DReg>::dotTest();
		// }

		{ // FFT
			std::shared_ptr<complex1DReg> m1(new complex1DReg(a1));
			std::shared_ptr<complex1DReg> d1(new complex1DReg(a1));

			inp1->random();
			out1->random();

			FFT1 fft_1(inp1,"out",1,1);
			std::cerr << "FFT<complex1DReg,complex1DReg>" << std::endl;

			fft_1.forward(inp1,d1,0);
			fft_1.adjoint(m1,out1,0);

			std::cerr << inp1->dot(m1) << std::endl;
			std::cerr << std::conj(out1->dot(d1)) << std::endl;
			std::cerr << std::conj(inp1->dot(m1))/out1->dot(d1)-1. << std::endl;
		}

		{
			std::cerr << "FFT<float2DReg,complex2DReg>" << std::endl;
			std::shared_ptr<float3DReg> w (new float3DReg(hyp3));
			std::shared_ptr<float3DReg> w1 (new float3DReg(hyp3));
			std::shared_ptr<complex3DReg> outC (new complex3DReg(hyp3));
			std::shared_ptr<complex3DReg> output (new complex3DReg(hyp3));

			w->random();
			outC->random();

			FFT1 fft_r1(w,output,1,2);
			fft_r1.forward(w,output,0);

			FFT1 fft_r2(w1,output,1,2);
			fft_r2.adjoint(w1,output,0);

			w1->scaleAdd(w,1,-1);
			std::cerr << w1->dot(w1)<< std::endl;
		}

		{ // Taper
			Taper tap(a1.n,10);
			tap.setDomainRange(inp1,out1);
			std::cerr << "Taper" << std::endl;
			tap.dotTest();
		}


		// { // dReflect
		// 	std::shared_ptr<float2DReg> w (new float2DReg(hyp));
		// 	for (int i=0; i<n1/2; i++) {
		// 		for (int j=0; j<n2; j++) {
		// 			(*w->_mat)[i][j] = 0.001;
		// 		}
		// 	}
		// 	for (int i=n1/2; i<n1; i++) {
		// 		for (int j=0; j<n2; j++) {
		// 			(*w->_mat)[i][j] = 0.002;
		// 		}
		// 	}
		// 	dReflect drefl(w);
		// 	drefl.setDomainRange(weight_m,weight_d);
		// 	std::cerr << "dReflect" << std::endl;
		// 	drefl.dotTest();

		// 	std::shared_ptr<complex2DReg> input (new complex2DReg(hyp));
		// 	std::shared_ptr<complex2DReg> output (new complex2DReg(hyp));
		// 	Reflect refl (w);
		// 	refl.setDomainRange(input,output);
		// 	std::cerr << "Reflect" << std::endl;
		// 	refl.dotTest();
		// }



	return 0;
}
