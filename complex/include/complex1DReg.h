#pragma once
#include <complexHyper.h>
#include <cstdint>
#include <iostream>
#include "boost/multi_array.hpp"

typedef boost::multi_array<std::complex<float>, 1> complex1D;

namespace SEP 
{
	class complex1DReg : public complexHyper {

 		public:

			/******************************* Constructors *******************************/
			complex1DReg() { ; }

			complex1DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }

			complex1DReg(const int n)
			{
				std::vector<SEP::axis> a(1, SEP::axis(n));
				std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
				initNoData(hyp);
			}

  			complex1DReg(const SEP::axis &a)
  			{
    			std::vector<SEP::axis> as(1, a);
    			std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(as));
    			initNoData(hyp);
  			}

  			complex1DReg(std::shared_ptr<SEP::hypercube> hyper, const complex1D &vals)
  			{
				initData(hyper, vals);
			}

			complex1DReg(const int n, complex1D &vals)
			{
				std::vector<SEP::axis> a(1, SEP::axis(n));
    			std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    			initData(hyp, vals);
  			}

  			complex1DReg(const SEP::axis &a, const complex1D &vals)
  			{
    			std::vector<SEP::axis> as(1, a);
    			std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(as));
    			initData(hyp, vals);
  			}
  
  			/******************************* Member functions ***************************/
  			void allocate()
  			{
    			std::vector<int> ns = getHyper()->getNs();
    			_mat.reset(new complex1D(boost::extents[ns[0]]));
    			setData(_mat->data());
    		}

  			std::shared_ptr<complex1DReg> clone() const;
  			std::shared_ptr<complex1DReg> cloneSpace() const;

  			virtual void cleanMemory() { setSpace(); }

  			/******************************* Member variables ***************************/			
  			std::shared_ptr<complex1D> _mat;

 			protected:
  				void initNoData(std::shared_ptr<SEP::hypercube> hyp);
  				//  void initData(std::shared_ptr<SEP::hypercube> hyp, const float *vals);
  				void initData(std::shared_ptr<SEP::hypercube> hyp, const complex1D &vals);
};

}  // namespace SEP
