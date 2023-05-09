#pragma once
#include <float1DReg.h>
#include <float2DReg.h>
#include "boost/multi_array.hpp"

namespace SEP {


class RefSampler

	{
	public:

		RefSampler(std::shared_ptr<float2DReg> slow, int nref);

		inline float getRefSlow(int iz, int iref) {return _slow_ref[iz][iref];}
		inline std::vector<int> getRefLoc(int iz, int iref) { return index[iz][iref];}

	private:

		float3D _slow_ref;
		boost::multi_array<std::vector<int>,3> index;
		
		void uniform_sample(std::shared_ptr<float2DReg> slow, int nref);
		void getIndexes();

		bool compareSlow (int a, int b);
		bool compareSlowRef (int a, int b);

		std::shared_ptr<float3DReg> _slow;
		float1D _slow_slice;
		float1D _slow_ref_slice;
		std::vector<int> is_sort, isref_sort;


		int _nx, _nref, _nz, _nw;

	};
}
