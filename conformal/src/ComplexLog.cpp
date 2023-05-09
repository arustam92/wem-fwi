#include <ComplexLog.h>

using namespace SEP;

std::shared_ptr<complex2DReg> ComplexLog::forward(std::shared_ptr<complex2DReg> inSpace) {
  if (!inSpace->getHyper()->checkSame(getInHyper())) throw(SEPException(std::string("Input space is wrong")));

  std::shared_ptr<complex2DReg> outSpace = std::make_shared<complex2DReg>(getOutHyper());
  axis ax1 = inSpace->getHyper()->getAxis(1);
  axis ax2 = inSpace->getHyper()->getAxis(2);

  int iz, ix;
  std::complex<float> zmap, z0, vz1, vz2;
  std::complex<float> v01=0,v10=0,v11=0;
  float az, ax;
  complex1D _buffer_;
  _buffer_.resize(boost::extents[getOutHyper()->getAxis(1).n]);
  for (int i2=0; i2 < ax2.n; ++i2) {
    for (int i1=0; i1 < ax1.n; ++i1) {
      zmap = getInMap()[i1 + i2*ax1.n];
      iz = (zmap.real() - ax2.o) / ax2.d;
      ix = (zmap.imag() - ax1.o) / ax1.d;
      // nearest neighbor
      z0 = getInCoord()[ix + iz*ax1.n];
      az = (zmap.real()-z0.real()) / ax2.d;
      ax = (zmap.imag()-z0.imag()) / ax1.d;
      if (iz > 0 && iz < ax2.n && ix > 0 && ix < ax1.n) {
        if (iz + 1 < ax2.n) v10 = (*inSpace->_mat)[iz+1][ix];
        else v10 = 0;
        if (ix + 1 < ax1.n) v01 = (*inSpace->_mat)[iz][ix+1];
        else v01 = 0;
        if (iz+1 < ax2.n && ix+1 < ax1.n) v11 = (*inSpace->_mat)[iz+1][ix+1];
        else v11 = 0;
        vz1 = (1-az)*(*inSpace->_mat)[iz][ix] + az*v10;
        vz2 = (1-az)*v01 + az*v11;
        _buffer_[i1] = (1-ax)*vz1 + ax*vz2;
      }
      (*outSpace->_mat)[i2][i1] = _buffer_[i1];
    }
  }
  return outSpace;
};

std::shared_ptr<complex2DReg> ComplexLog::inverse(std::shared_ptr<complex2DReg> outSpace) {
  if (!outSpace->getHyper()->checkSame(getOutHyper())) throw(SEPException(std::string("Output space is wrong")));

  std::shared_ptr<complex2DReg> inSpace = std::make_shared<complex2DReg>(getInHyper());
  axis ax1 = outSpace->getHyper()->getAxis(1);
  axis ax2 = outSpace->getHyper()->getAxis(2);

  int iz, ix;
  std::complex<float> umap;
  for (int i2=0; i2 < ax2.n; ++i2) {
    for (int i1=0; i1 < ax1.n; ++i1) {
      umap = getOutMap()[i1 + i2*ax1.n];
      float amp = std::abs(getInCoord()[i1 + i2*ax1.n]);
      if (amp > std::exp(_eps_) && amp < .99f) {
        iz = (umap.real() - ax2.o) / ax2.d;
        ix = (umap.imag() - ax1.o) / ax1.d;
        (*inSpace->_mat)[i2][i1] = (*outSpace->_mat)[iz][ix];
      }

    }
  }
  return inSpace;
};

void ComplexLog::forward(const std::shared_ptr<complex2DReg>& inSpace, std::shared_ptr<complex2DReg>& outSpace) {
  if (!inSpace->getHyper()->checkSame(getInHyper())) throw(SEPException(std::string("Input space is wrong")));

  outSpace->setHyper(getOutHyper());

  axis ax1 = inSpace->getHyper()->getAxis(1);
  axis ax2 = inSpace->getHyper()->getAxis(2);

  int iz, ix;
  std::complex<float> zmap, z0, vz1, vz2;
  std::complex<float> v01=0,v10=0,v11=0;
  float az, ax;
  complex1D _buffer_;
  _buffer_.resize(boost::extents[getOutHyper()->getAxis(1).n]);
  for (int i2=0; i2 < ax2.n; ++i2) {
    for (int i1=0; i1 < ax1.n; ++i1) {
      zmap = getInMap()[i1 + i2*ax1.n];
      iz = (zmap.real() - ax2.o) / ax2.d;
      ix = (zmap.imag() - ax1.o) / ax1.d;
      // nearest neighbor
      z0 = getInCoord()[ix + iz*ax1.n];
      az = (zmap.real()-z0.real()) / ax2.d;
      ax = (zmap.imag()-z0.imag()) / ax1.d;
      if (iz >= 0 && iz < ax2.n && ix >= 0 && ix < ax1.n) {
        // if (iz + 1 < ax2.n) v10 = (*inSpace->_mat)[iz+1][ix];
        // else v10 = 0;
        // if (ix + 1 < ax1.n) v01 = (*inSpace->_mat)[iz][ix+1];
        // else v01 = 0;
        // if (iz+1 < ax2.n && ix+1 < ax1.n) v11 = (*inSpace->_mat)[iz+1][ix+1];
        // else v11 = 0;
        // vz1 = (1-az)*(*inSpace->_mat)[iz][ix] + az*v10;
        // vz2 = (1-az)*v01 + az*v11;
        // _buffer_[i1] = (1-ax)*vz1 + ax*vz2;
        _buffer_[i1] = (*inSpace->_mat)[iz][ix];
      }
      (*outSpace->_mat)[i2][i1] = _buffer_[i1];
    }
  }
};

void ComplexLog::inverse(std::shared_ptr<complex2DReg>& inSpace, const std::shared_ptr<complex2DReg>& outSpace) {
  if (!outSpace->getHyper()->checkSame(getOutHyper())) throw(SEPException(std::string("Output space is wrong")));

  inSpace->setHyper(getInHyper());

  axis ax1 = outSpace->getHyper()->getAxis(1);
  axis ax2 = outSpace->getHyper()->getAxis(2);

  int iz, ix;
  std::complex<float> umap;
  for (int i2=0; i2 < ax2.n; ++i2) {
    for (int i1=0; i1 < ax1.n; ++i1) {
      umap = getOutMap()[i1 + i2*ax1.n];
      float amp = std::abs(getInCoord()[i1 + i2*ax1.n]);
      if (amp > std::exp(_eps_) && amp < .99f) {
        iz = (umap.real() - ax2.o) / ax2.d;
        ix = (umap.imag() - ax1.o) / ax1.d;
        (*inSpace->_mat)[i2][i1] = (*outSpace->_mat)[iz][ix];
      }

    }
  }
};

std::shared_ptr<ConformalMap::Point> ComplexLog::forward(std::vector<int>& re, std::vector<int>& im) {

  std::shared_ptr<Point> point = std::make_shared<Point>();
  int n1 = getInHyper()->getAxis(1).n;
  axis outax1 = getOutHyper()->getAxis(1);
  axis outax2 = getOutHyper()->getAxis(2);
  std::complex<float> c, z;
  if (re.size() > 1) {
    for (int i=0; i < re.size(); ++i) {
      int ind = im[i] + re[i]*n1;
      c = getOutMap()[ind];
      z = getInCoord()[ind];
      // if (std::abs(z) < std::exp(_eps_)) {
      //   point->real.resize(outax2.n);
      //   point->imag.resize(outax1.n);
      //   std::iota(point->imag.begin(),point->imag.end(),0);
      //   std::fill(point->real.begin(),point->real.end(),0);
      //   iz = 0;
      //   ix =
      // }
      // else {
        int iz = (c.real() - outax2.o) / outax2.d;
        iz = std::max(0,iz);
        int ix = (c.imag() - outax1.o) / outax1.d;
        point->real.push_back(iz);
        point->imag.push_back(ix);
      // }
      }
  }
  else
    point = forward(re[0],im[0]);

  return point;
};

std::shared_ptr<ConformalMap::Point> ComplexLog::forward(int& re, int& im) {

  std::shared_ptr<Point> point = std::make_shared<Point>();
  int n1 = getInHyper()->getAxis(1).n;
  axis outax1 = getOutHyper()->getAxis(1);
  axis outax2 = getOutHyper()->getAxis(2);
  std::complex<float> c, z;
    int ind = im + re*n1;
    c = getOutMap()[ind];
    z = getInCoord()[ind];
    // if (std::abs(z) <= std::exp(_eps_)) {
      point->real.resize(outax1.n);
      point->imag.resize(outax1.n);
      std::iota(point->imag.begin(),point->imag.end(),0);
      std::fill(point->real.begin(),point->real.end(),0);
    // }
    // else {
    //   int iz = (c.real() - outax2.o) / outax2.d;
    //   int ix = (c.imag() - outax1.o) / outax1.d;
    //   point->real.push_back(iz);
    //   point->imag.push_back(ix);
    // }
  return point;
};
