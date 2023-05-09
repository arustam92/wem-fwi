#pragma once
#include <complexHyper.h>
#include <float2DReg.h>
#include "boost/multi_array.hpp"
typedef boost::multi_array<std::complex<float>, 2> complex2D;
namespace SEP {
class complex2DReg : public complexHyper {
 public:
  complex2DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  complex2DReg(const int n1, const int n2) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  complex2DReg(const SEP::axis &a1, const SEP::axis &a2) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  complex2DReg(std::shared_ptr<SEP::hypercube> hyper, const complex2D &vals) {
    initData(hyper, vals);
  }
  complex2DReg(int n1, int n2, complex2D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  complex2DReg(SEP::axis &a1, SEP::axis &a2, complex2D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new complex2D(boost::extents[ns[1]][ns[0]]));

    setData(_mat->data());
    ;
  }
  std::shared_ptr<complex2DReg> clone() const;
  std::shared_ptr<complex2DReg> cloneSpace() const;
  virtual void cleanMemory() { setSpace(); }
  std::shared_ptr<complex2D> _mat;

  std::shared_ptr<float2DReg> real();
  std::shared_ptr<float2DReg> imag();

 private:
  void initNoData(std::shared_ptr<SEP::hypercube> hyp);
  void initData(std::shared_ptr<SEP::hypercube> hyp, const complex2D &vals);
};
}  // namespace SEP
