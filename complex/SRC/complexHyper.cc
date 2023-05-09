#include <complexHyper.h>
#include <hypercube.h>
#include <iostream>
#include <random>
using namespace SEP;

void complexHyper::add(const std::shared_ptr<complexHyper> vec2) 
{
	assert(checkSame(vec2));
  	std::shared_ptr<complexHyper> vec2H = std::dynamic_pointer_cast<complexHyper>(vec2);

  	for (long long i = 0; i < _hyper->getN123(); i++) _vals[i] += vec2H->_vals[i];
  	calcCheckSum();
}

void complexHyper::scaleAdd(const std::complex<float> sc1, std::shared_ptr<complexHyper> vec2,
                          const std::complex<float> sc2) {
  assert(checkSame(vec2));
  std::shared_ptr<complexHyper> vec2H =
      std::dynamic_pointer_cast<complexHyper>(vec2);

  for (long long i = 0; i < _hyper->getN123(); i++)
    _vals[i] = _vals[i] * sc1 + sc2 * vec2H->_vals[i];
  calcCheckSum();
}

// void floatHyper::signum() {
//   assert(!spaceOnly());
//   for (long long i = 0; i < _hyper->getN123(); i++) {
//     if(_vals[i]>1e-20) _vals[i]=1;
//     else if(_vals[i]<-1e-20) _vals[i]=-1;
//     else _vals[i]=0;
//   }
//   calcCheckSum();
// }

void complexHyper::scale(std::complex<float> sc) {
  assert(!spaceOnly());
  for (long long i = 0; i < _hyper->getN123(); i++) _vals[i] = _vals[i] * sc;
  calcCheckSum();
}

void complexHyper::random() {
  assert(!spaceOnly());
  double re,im;
  for (long long i = 0; i < _hyper->getN123(); i++) {
    re = ((double)rand() / (RAND_MAX)) - .5;
    im = ((double)rand() / (RAND_MAX)) - .5;
    _vals[i] = {re,im};
  }
    // _vals[i] = ((double)rand() / (RAND_MAX)) - .5;
  calcCheckSum();
}

double complexHyper::L2Obj() const {
  double x = 0;
  std::cerr<<"running l2 norm"<<std::endl;
  for (long long i = 0; i < _hyper->getN123(); i++) {
    x += (_vals[i] * std::conj(_vals[i])).real();
  }
  return x;
}

double complexHyper::L1Obj() const {
  std::cerr<<"running l1 norm"<<std::endl;
  double x = 0;
  for (long long i = 0; i < _hyper->getN123(); i++) {
    x += std::abs(_vals[i]);
  }
  return x;
}

double complexHyper::dot(const std::shared_ptr<complexHyper> vec2) const {
  assert(checkSame(vec2));
  std::shared_ptr<complexHyper> vec2H =
      std::dynamic_pointer_cast<complexHyper>(vec2);

  double dt = 0.;
  for (long long i = 0; i < _hyper->getN123(); i++) {
    dt += double((_vals[i] * std::conj(vec2H->_vals[i])).real());
  }

  return dt;
}

void complexHyper::createMask(const float zero, const float err) {
  for (long long i = 0; i < _hyper->getN123(); i++) {
    if (std::abs(_vals[i]) - std::abs(zero) > err)

      _vals[i] = {0.,0};
    else
      _vals[i] = {1,1};
  }
  calcCheckSum();
}

void complexHyper::infoStream(const int lev, std::stringstream &x) {
  _hyper->infoStream(x);
  if (spaceOnly())
    x << "Only space\n";
  else {
    x << "Allocated\n";
    long long npts = std::min((const long long)lev, getHyper()->getN123());
    for (long long i = 0; i < npts; i++)
      x << std::to_string(i) << std::string(" (") << std::to_string(_vals[i].real()) << "," << std::to_string(_vals[i].imag())
        << std::endl;
  }
}
// void floatHyper::softClip(const float scale) {
//   float sc2 = scale * scale;
//   for (int i = 0; i < _hyper->getN123(); i++)
//     _vals[i] = scale * _vals[i] / sqrtf(1. + sc2 * _vals[i] * _vals[i]);
//   calcCheckSum();
// }

float complexHyper::absMax() const {
  float val = std::abs(_vals[0]);
  for (int i = 1; i < _hyper->getN123(); i++)
    val = std::max(val, std::abs(_vals[i]));
  return val;
}
// float floatHyper::max() const {
//   float val = fabsf(_vals[0]);
//   for (int i = 1; i < _hyper->getN123(); i++)
//     val = std::max(val, _vals[i]);
//   return val;
// }
// float floatHyper::min() const {
//   float val = fabsf(_vals[0]);
//   for (int i = 1; i < _hyper->getN123(); i++)
//     val = std::min(val, _vals[i]);
//   return val;
// }
void complexHyper::calcCheckSum() {
  uint32_t sum1 = 0, sum2 = 0;
  uint32_t *data = (uint32_t *)_vals;
  uint32_t mx = 4294967295;
  for (long long i = 0; i < _hyper->getN123(); i++) {
    sum1 = (sum1 + data[i]) % mx;
    sum2 = (sum2 + sum1) % mx;
  }
  setCheckSum(sum2 * 2 ^ 32 + sum1);
}

bool complexHyper::checkSame(const std::shared_ptr<complexHyper> vec2) const {
  if (!vec2) {
    std::cerr << "Not allocated vec2" << std::endl;
    return false;
  }
//  if (_hyper == vec2->getHyper()) return true;
  return true;
  std::cerr<<_hyper->getAxis(1).n<<" "<<vec2->getHyper()->getAxis(1).n<<std::endl;
  std::cerr << "Not from the same Hypercube" << std::endl;

  return false;
}
