#include <FFT1.h>

using namespace SEP;

FFT1::FFT1 (std::shared_ptr<complexHyper> model, std::string mode, int rank, int a) {

	_rank = rank;
	_axis = a;

	int ndim = model->getHyper()->getNdimG1();
	size = model->getHyper()->getAxis(_axis).n;

	howmany = 1;
	for (int i=ndim;i>_rank;i--) howmany *= model->getHyper()->getAxis(i).n;

	/*
	j-th element of k-th transform is
	out = in + j*stride + k*idist
	*/
 	idist = size;
 	odist = size;
	istride = 1;
	ostride = 1;

 	// fftwf_complex *_in = reinterpret_cast<fftwf_complex*>(model->getVals());
 	// fftwf_complex *_out = reinterpret_cast<fftwf_complex*>(data->getVals());
	fftwf_complex *in = (fftwf_complex*) fftwf_malloc(size*howmany*sizeof(fftwf_complex));
	fftwf_complex *out = (fftwf_complex*) fftwf_malloc(size*howmany*sizeof(fftwf_complex));


 	{
		tbb::mutex::scoped_lock lock(mutex);
	 	if (mode=="in") {
	 		_fwd_plan = fftwf_plan_dft_1d(size, in, in, FFTW_FORWARD,FFTW_ESTIMATE);
 	 		_adj_plan = fftwf_plan_dft_1d(size, in, in, FFTW_BACKWARD,FFTW_ESTIMATE);
	 	}
	 	if (mode=="out") {
	 		_fwd_plan = fftwf_plan_dft_1d(size, in, out, FFTW_FORWARD,FFTW_ESTIMATE);
 	 		_adj_plan = fftwf_plan_dft_1d(size, out, in, FFTW_BACKWARD,FFTW_ESTIMATE);
	 	}

 	}
 	// {

 	// _fwd_plan = fftwf_plan_many_dft(_rank, &size, howmany,
  //                            in, NULL, istride, idist,
  //                            out, NULL, ostride, odist,
  //                            FFTW_FORWARD,FFTW_ESTIMATE);
 	// _adj_plan = fftwf_plan_many_dft(_rank, &size, howmany,
  //                            out, NULL, istride, idist,
  //                            in, NULL, ostride, odist,
  //                            FFTW_BACKWARD,FFTW_ESTIMATE);
 	// }
 	fftwf_free(in);
 	fftwf_free(out);

	float pi = 4.f * std::atan(1.f);
	float dk = 2*pi/(model->getHyper()->getAxis(1).d*size);
	k.resize(boost::extents[size]);
	shift.resize(boost::extents[size]);
	for (int ik=0; ik<k.size()/2; ik++) {
		k[ik] = ik*dk;
		k[size-ik-1] = -k[ik];
		// k[ik+nk/2] = (-nk/2 + ik)*dk;
	}
	float o1 = model->getHyper()->getAxis(1).o;
	std::complex<float> I = {0.f,1.f};
	for (int ik=0; ik < k.size(); ++ik) {
		shift[ik] = std::exp(I*o1*k[ik]);
	}
}

void FFT1::forward (std::shared_ptr<complexHyper> model, std::shared_ptr<complexHyper> data, bool add) {

	if (!add) data->scale(0.);

	_in = (reinterpret_cast<fftwf_complex*>(model->getVals()));
	_out = (reinterpret_cast<fftwf_complex*>(data->getVals()));

	fftwf_execute_dft(_fwd_plan,_in,_out);

	for (int ik=0; ik < k.size(); ++ik)
		data->getVals()[ik] *= shift[ik];

	data->scale(1.0/std::sqrt(size));

}

void FFT1::adjoint (std::shared_ptr<complexHyper> model, std::shared_ptr<complexHyper> data, bool add) {

	if (!add) model->scale(0.);

	_in = (reinterpret_cast<fftwf_complex*>(model->getVals()));
	_out = (reinterpret_cast<fftwf_complex*>(data->getVals()));

	for (int ik=0; ik < k.size(); ++ik)
		data->getVals()[ik] /= shift[ik];

	fftwf_execute_dft(_adj_plan, _out, _in);

	model->scale(1.0/std::sqrt(size));

}
