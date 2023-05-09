#include <axis.h>
#include <float1DReg.h>
#include <float2DReg.h>
#include <complex1DReg.h>
#include <complex2DReg.h>
#include <hypercube.h>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <cstdlib>

// typedef boost::multi_array<std::complex<float>, 3> float2D;
// typedef boost::multi_array<std::complex<float>, 2> complex2D;
using namespace SEP;
int main(int argc, char const *argv[])
{		
		int n2,n1;
		n2=5;
		n1=1000;

		axis a1(n1,0,0.001);
		axis a2(n2,0,5);

		std::shared_ptr<hypercube> hyp (new hypercube(a1,a2));

		std::shared_ptr<float2DReg> input (new float2DReg(a1,a2));
		std::shared_ptr<float1DReg> output (new float1DReg(a1));

		std::vector<std::vector<int>> v;


		for (int i=0; i<n2; i++) {			
			for (int k=0; k<n1; k++) {
						
						(*input->_mat)[i][k] =std::rand()%10; 
						std::cout  << (*input->_mat)[i][k]  << std::endl;
				}
				std::cout  << "ZZZZZZZZZZZZZZZZZZZZZ\n" << std::endl;

		}

		std::sort(input->_mat->data(),input->_mat->data()+n1);

		for (int i=0; i<n2; i++) {			
			for (int k=0; k<n1; k++) {
						
						std::cout  << (*input->_mat)[i][k]  << std::endl;
				}
				std::cout  << "\n" << std::endl;

		}

		float *i = std::lower_bound(input->_mat->data(),input->_mat->data()+n1,5);
		int off = i-input->_mat->data();
		std::cerr << (*input->_mat)[0][off] << std::endl;


	return 0;
}