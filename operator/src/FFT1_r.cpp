	#include <float2DReg.h>
#include <FFT1.h>


using namespace SEP;

FFT1::FFT1 (std::shared_ptr<floatHyper> model, std::shared_ptr<complexHyper> data,  int rank, int a) {

	_rank = rank;
	_axis = a;

	int ndim = model->getHyper()->getNdimG1();
	size = model->getHyper()->getAxis(_axis).n;

	howmany = 1;
	// for (int i=ndim;i>_rank;i--) howmany *= model->getHyper()->getAxis(i).n;
	for (int i=1;i<=ndim;i++) howmany *= model->getHyper()->getAxis(i).n;
	howmany /= size;
	/*
	j-th element of k-th transform is
	out = in + j*stride + k*idist
	obnl
	*/
 	idist = a > 1 ? 1 : size ;
	// if (a == 1) idist = size;
	// if (a == 3) idist = 1;
 	odist = idist;

	istride = 1;
	for (int i = a; i > 1; i--) istride *= model->getHyper()->getAxis(i-1).n;
	ostride = istride;

	dataFFT = std::make_shared<fft_data>(size*howmany);

	for (int i=0; i<size*howmany; i++) {
		dataFFT->getModel()[i][0] = model->getVals()[i];
		dataFFT->getModel()[i][1] = 0.;
	}

	// ERROR HEREEEE
	dataFFT->setData(data->getVals());


 	{
		tbb::mutex::scoped_lock lock(mutex);
 	_fwd_plan = fftwf_plan_many_dft(_rank, &size, howmany,
	                             dataFFT->getModel(), NULL, istride, idist,
	                             dataFFT->getData(), NULL, ostride, odist,
	                             FFTW_FORWARD,FFTW_ESTIMATE);

 	_adj_plan = fftwf_plan_many_dft(_rank, &size, howmany,
                             dataFFT->getData(), NULL, istride, idist,
                             dataFFT->getModel(), NULL, ostride, odist,
                             FFTW_BACKWARD,FFTW_ESTIMATE);

	float pi = 4.f * std::atan(1.f);
	float dk = 2*pi/(model->getHyper()->getAxis(1).d*size);
	k.resize(boost::extents[size]);
	shift.resize(boost::extents[size]);
	for (int ik=0; ik<k.size()/2; ik++) 
		k[ik] = ik*dk;
	for (int ik=1; ik<k.size()/2; ik++) 
		k[k.size()-ik] = -k[ik];

	float o1 = model->getHyper()->getAxis(1).o;
	std::complex<float> I = {0.f,1.f};
	for (int ik=0; ik < k.size(); ++ik) {
		shift[ik] = std::exp(I*o1*k[ik]);
	}
	// shift[k.size()/2] = 1.f;
 	}
}

void FFT1::forward (const std::shared_ptr<floatHyper> model, std::shared_ptr<complexHyper> data, bool add) {

	if(!add) data->scale(0.);
	fftwf_execute(_fwd_plan);

	// for (int ik=0; ik < k.size(); ++ik)
	// 	data->getVals()[ik] *= shift[ik];
	data->scale(1./std::sqrt(size));

}

void FFT1::adjoint (std::shared_ptr<floatHyper> model, const std::shared_ptr<complexHyper> data, bool add) {

	if(!add) model->scale(0.);

	// for (int ik=0; ik < k.size(); ++ik)
	// 	data->getVals()[ik] /= shift[ik];
	fftwf_execute(_adj_plan);

	float scale = 1.0/std::sqrt(size);
	for (int i=0; i < size*howmany; i++) model->getVals()[i] += dataFFT->getModel()[i][0] * scale;

}
