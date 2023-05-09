
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <BornFull.h>
#include <chrono>

using namespace std::chrono;

using namespace SEP;

BornFull::BornFull (std::shared_ptr<complex1DReg> wave, std::shared_ptr<complex2DReg> slow, std::shared_ptr<paramObj> par, bool verb,
std::shared_ptr<float2DReg> scoord, std::shared_ptr<float2DReg> rcoord) :
Acquisition(slow->getHyper(),par,scoord,rcoord) {

	param = par;
	_slow = slow;
	_verb_ = verb;

	nref = par->getInt("nref",1);
	setBgSlow(slow);

	wave_f = wave;

	float df = wave_f->getHyper()->getAxis(1).d;
	fmin = 0;
	fmax = wave_f->getHyper()->getAxis(1).n;

	freq.resize(boost::extents[fmax]);
	for (int i=0; i<freq.size(); i++) {
		freq[i] = (i*df) + wave_f->getHyper()->getAxis(1).o;
	}

	if (param->getBool("illum",false)) {
		c_ptr_adjoint = &BornFull::c_adj_illum;
	}
	else {
		c_ptr_adjoint = &BornFull::c_adj_reg;
	}

	nshots = _src->getZ().size();

}

void BornFull::forward (std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) data->scale(0.);

	if(_verb_) std::cerr << "Forward complex born ... " << std::endl;
		auto start = high_resolution_clock::now();

	(*this.*fwd_propagate)(model,data);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	if(_verb_) std::cerr << duration.count()/1e6 << " s" << std::endl;

}

void BornFull::c_adj_illum (std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) model->scale(0.);

	if(_verb_) std::cerr << "Adjoint complex born illumination... " << std::endl;
		auto start = high_resolution_clock::now();

	illum = std::make_shared<Illumination> (_slow->getHyper());

	(*this.*adj_propagate)(model,data);

	illum->compensate(model);

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	if(_verb_)	std::cerr << duration.count()/1e6 << " s" << std::endl;
}

void BornFull::c_adj_reg (std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) model->scale(0.);

	if(_verb_) std::cerr << "Adjoint complex born ... " << std::endl;
		auto start = high_resolution_clock::now();

		(*this.*adj_propagate)(model,data);

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	if(_verb_)	std::cerr << duration.count()/1e6 << " s" << std::endl;
}
