#include <RefSampler.h>
#include <numeric>
#include <algorithm>
#include <functional>
#include <vector>

using namespace SEP;
using namespace std::placeholders;

RefSampler::RefSampler(const std::shared_ptr<complex2DReg>& slow, int nref) : _slow(slow) {
			_nref=nref;
			_nx = _slow->getHyper()->getAxis(1).n;
			_nz = _slow->getHyper()->getAxis(2).n;

			_slow_ref.resize(boost::extents[_nz][nref]);
			index.resize(boost::extents[_nz][nref]);
			_slow_slice.resize(boost::extents[_nx]);
			_slow_ref_slice.resize(boost::extents[_nx]);

			is_sort.assign(_nx,0);
			isref_sort.assign(nref,0);

			uniform_sample();
			getIndexes();
		};


void RefSampler::uniform_sample() {

		int o = _nx/(2*_nref);
		for (int iz=0; iz < _nz; iz++) {
			for (int iref=0; iref<_nref; iref++) {
				int off = o+iref*_nx/_nref;
				_slow_ref[iz][iref] = (*_slow->_mat)[iz][off];
				if (_slow_ref[iz][iref].real() <= 0)
				throw(SEPException(std::string("Ref velocity ") + std::to_string(_slow_ref[iz][iref].real()) +  std::string(" at the layer ") + std::to_string(iz)));
			}
		}
}

void RefSampler::getIndexes() {

	for (int iz=0; iz<_nz; iz++) {

		for (int ix=0; ix<_nx; ++ix) _slow_slice[ix] = (*_slow->_mat)[iz][ix].real();

		// std::copy_n(_slow->_mat->data()+iz*_nx,_nx,_slow_slice.data());

		for (int iref=0; iref<_nref; iref++) {
			_slow_ref_slice[iref] = _slow_ref[iz][iref].real();
		}

		std::iota(is_sort.begin(),is_sort.end(),0);
		std::iota(isref_sort.begin(),isref_sort.end(),0);

		// sort indexes of slowness slice according to the slowness values
		std::sort(is_sort.begin(),is_sort.end(),std::bind(&RefSampler::compareSlow,this,_1,_2));
		// sort indexes of REF slowness slice according to the REF slowness values
		std::sort(isref_sort.begin(),isref_sort.end(),std::bind(&RefSampler::compareSlowRef,this,_1,_2));
		// sort slowness slice
		std::sort(_slow_slice.data(),_slow_slice.data()+_nx);
		// sort reference slowness slice
		std::sort(_slow_ref_slice.data(),_slow_ref_slice.data()+_nx);

		int off = 0;
		while (_slow_slice[off]-_slow_ref_slice[0] < _slow_ref_slice[1]-_slow_slice[off]) off++;

		index[iz][isref_sort[0]] = std::vector<int> (is_sort.begin(),is_sort.begin()+off);

		int off_beg = off;
		// reference slowness between first+1 and last-1
		for (int iref=1; iref<_nref-1; iref++) {
			while (_slow_slice[off]-_slow_ref_slice[iref] <	_slow_ref_slice[iref+1]-_slow_slice[off]) off++;

			index[iz][isref_sort[iref]] = std::vector<int> (is_sort.begin()+off_beg,is_sort.begin()+off);

			off_beg=off;
		}

		// and the last reference slowness
		index[iz][isref_sort[_nref-1]] = std::vector<int> (is_sort.begin()+off,is_sort.end());
	}

}

bool RefSampler::compareSlow (int a, int b) {
	return _slow_slice[a] < _slow_slice[b];
}
bool RefSampler::compareSlowRef (int a, int b) {
	return _slow_ref_slice[a] < _slow_ref_slice[b];
}
