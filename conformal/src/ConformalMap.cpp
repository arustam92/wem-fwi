#include <ConformalMap.h>

using namespace SEP;

std::shared_ptr<complex2DReg> ConformalMap::forward(std::shared_ptr<complex2DReg> inSpace) {
  if (!inSpace->getHyper()->checkSame(getInHyper())) throw(std::string("Input space is wrong"));

  std::shared_ptr<complex2DReg> outSpace = std::make_shared<complex2DReg>(getOutHyper());
  axis ax1 = inSpace->getHyper()->getAxis(1);
  axis ax2 = inSpace->getHyper()->getAxis(2);

  int iz, ix;
  std::complex<float> zmap;
  for (int i2=0; i2 < ax2.n; ++i2) {
    for (int i1=0; i1 < ax1.n; ++i1) {
      zmap = getInMap()[i1 + i2*ax1.n];
      iz = (zmap.real() - ax2.o) / ax2.d;
      ix = (zmap.imag() - ax1.o) / ax1.d;
      if (iz >= 0 && iz < ax2.n && ix >= 0 && ix < ax1.n)
      (*outSpace->_mat)[i2][i1] = (*inSpace->_mat)[iz][ix];
    }
  }
  return outSpace;
};

std::shared_ptr<complex2DReg> ConformalMap::inverse(std::shared_ptr<complex2DReg> outSpace) {
  if (!outSpace->getHyper()->checkSame(getOutHyper())) throw(std::string("Output space is wrong"));

  std::shared_ptr<complex2DReg> inSpace = std::make_shared<complex2DReg>(getInHyper());
  axis ax1 = outSpace->getHyper()->getAxis(1);
  axis ax2 = outSpace->getHyper()->getAxis(2);

  int iz, ix;
  std::complex<float> umap;
  for (int i2=0; i2 < ax2.n; ++i2) {
    for (int i1=0; i1 < ax1.n; ++i1) {
      umap = getOutMap()[i1 + i2*ax1.n];
      iz = (umap.real() - ax2.o) / ax2.d;
      ix = (umap.imag() - ax1.o) / ax1.d;
      if (iz >= 0 && iz < ax2.n && ix >= 0 && ix < ax1.n)
      (*inSpace->_mat)[i2][i1] = (*outSpace->_mat)[iz][ix];
    }
  }
  return inSpace;
};

void ConformalMap::forward(const std::shared_ptr<complex2DReg>& inSpace, std::shared_ptr<complex2DReg>& outSpace) {
  if (!inSpace->getHyper()->checkSame(getInHyper())) throw(std::string("Input space is wrong"));

  outSpace->setHyper(getOutHyper());

  axis ax1 = inSpace->getHyper()->getAxis(1);
  axis ax2 = inSpace->getHyper()->getAxis(2);

  int iz, ix;
  std::complex<float> zmap;
  for (int i2=0; i2 < ax2.n; ++i2) {
    for (int i1=0; i1 < ax1.n; ++i1) {
      zmap = getInMap()[i1 + i2*ax1.n];
      iz = (zmap.real() - ax2.o) / ax2.d;
      ix = (zmap.imag() - ax1.o) / ax1.d;
      if (iz >= 0 && iz < ax2.n && ix >= 0 && ix < ax1.n)
      (*outSpace->_mat)[i2][i1] = (*inSpace->_mat)[iz][ix];
    }
  }
};

void ConformalMap::inverse(std::shared_ptr<complex2DReg>& inSpace, const std::shared_ptr<complex2DReg>& outSpace) {
  if (!outSpace->getHyper()->checkSame(getOutHyper())) throw(std::string("Output space is wrong"));

  inSpace->setHyper(getInHyper());

  axis ax1 = outSpace->getHyper()->getAxis(1);
  axis ax2 = outSpace->getHyper()->getAxis(2);

  int iz, ix;
  std::complex<float> umap;
  for (int i2=0; i2 < ax2.n; ++i2) {
    for (int i1=0; i1 < ax1.n; ++i1) {
      umap = getOutMap()[i1 + i2*ax1.n];
      iz = (umap.real() - ax2.o) / ax2.d;
      ix = (umap.imag() - ax1.o) / ax1.d;
      if (iz >= 0 && iz < ax2.n && ix >= 0 && ix < ax1.n)
      (*inSpace->_mat)[i2][i1] = (*outSpace->_mat)[iz][ix];
    }
  }
};
