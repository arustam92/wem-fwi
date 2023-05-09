#include <axis.h>
#include <float2DReg.h>
#include <complex2DReg.h>
#include <hypercube.h>
#include <cmath>
#include <FFT1.h>

// typedef boost::multi_array<std::complex<float>, 3> float2D;
// typedef boost::multi_array<std::complex<float>, 2> complex2D;
using namespace SEP;
int main(int argc, char const *argv[])
{		
		int n2,n1;
		n2=5;
		n1=5;

		axis a1(n1,0,0.001);
		axis a2(n2,0,5);

		std::shared_ptr<hypercube> hyp (new hypercube(a1,a2));

		std::shared_ptr<complex2DReg> input (new complex2DReg(hyp));
		std::shared_ptr<complex2DReg> output (new complex2DReg(hyp));

	
			for (int j=0; j<n2; j++) {
				(*input->_mat)[j][0] = 1;
				for (int k=1; k<n1; k++) {
					
					(*input->_mat)[j][k] = 0;
		}
		}

		std::shared_ptr<Operator> fft1 (new FFT1(input,1,1));
		fft1->forward(input,output,0);

		for (int j=0; j<n2; j++) {
				for (int k=0; k<n1; k++) {
					
					std::cout << (*input->_mat)[j][k] <<(*output->_mat)[j][k] << std::endl; 
		}
		std::cout<<"\t"<<std::endl;
		}

		std::cout << output->getHyper()->getAxis(1).o << std::endl;


	return 0;
}