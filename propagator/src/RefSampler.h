#pragma once
#include <float1DReg.h>
#include <float2DReg.h>
#include "complex2DReg.h"
#include "complex3DReg.h"
#include "complex1DReg.h"
#include "boost/multi_array.hpp"

namespace SEP {


class RefSampler

	{
	public:

		RefSampler(const std::shared_ptr<complex2DReg>& slow, int nref);

		inline std::complex<float>& getRefSlow(int iz, int iref) {return _slow_ref[iz][iref];}
		inline std::vector<int>& getRefLoc(int iz, int iref) { return index[iz][iref];}

	private:

		complex2D _slow_ref;

		void uniform_sample();
		void getIndexes();

		bool compareSlow (int a, int b);
		bool compareSlowRef (int a, int b);

		const std::shared_ptr<complex2DReg>& _slow;
		float1D _slow_slice;
		float1D _slow_ref_slice;
		std::vector<int> is_sort, isref_sort;
		boost::multi_array<std::vector<int>,2> index;

		int _nx, _nref, _nz;

	};
}
