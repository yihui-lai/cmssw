#ifndef CMS_PLUTO_MIXED_DOUBLE_DALITZ_H
#define CMS_PLUTO_MIXED_DOUBLE_DALITZ_H
#include "TLorentzVector.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace cmspluto {
// Petri, arXiv:1010.2378, Eq. (3.36), divided by its angular bound 2.
// b1,b2 are beta squared; c1,c2 are cos(theta). No identical-fermion
// exchange diagram exists for the mixed-flavour double-Dalitz channel.
inline double mixedAngular(double b1, double b2, double c1, double c2, double phi) {
  const double a = b1 * (1-c1*c1), b = b2 * (1-c2*c2);
  return 1 - (a+b)/2 + a*b*std::pow(std::sin(phi),2)/2;
}

// LO pointlike pseudoscalar -> gamma* gamma* -> mu mu ee, F(s1,s2)=1.
// dGamma/(ds1 ds2 dcos1 dcos2 dphi) is proportional to
// lambda(1,s1/M^2,s2/M^2)^(3/2) beta1 beta2 angular / (s1 s2).
// Logarithmic s proposals cancel both photon-pole factors. Every remaining
// factor is bounded by one, giving a fixed, non-adaptive rejection envelope.
template<class Flat>
std::array<TLorentzVector,4> mixedPointlike(const TLorentzVector& parent,
                                         double m1, double m2, Flat flat) {
  const double M = parent.M(), M2 = M*M;
  if (!(m1 > 0 && m2 > 0 && M > 2*(m1+m2)))
    throw std::runtime_error("Mixed double-Dalitz masses are below threshold");
  const double low1=4*m1*m1, low2=4*m2*m2;
  const double log1=std::log(std::pow(M-2*m2,2)/low1);
  const double log2=std::log(std::pow(M-2*m1,2)/low2);
  const double pi=std::acos(-1.);
  for (unsigned attempt=0; attempt<1000000; ++attempt) {
    const double s1=low1*std::exp(log1*flat()), s2=low2*std::exp(log2*flat());
    if (std::sqrt(s1)+std::sqrt(s2) >= M) continue;
    const double b1=1-low1/s1, b2=1-low2/s2;
    const double x=s1/M2, y=s2/M2;
    const double lambda=std::max(0.,std::pow(1-x-y,2)-4*x*y);
    const double c1=2*flat()-1, c2=2*flat()-1, phi=2*pi*flat();
    const double weight=std::pow(lambda,1.5)*std::sqrt(b1*b2)*mixedAngular(b1,b2,c1,c2,phi);
    if (!std::isfinite(weight) || weight < 0 || weight > 1)
      throw std::runtime_error("Invalid mixed double-Dalitz rejection weight");
    if (flat() >= weight) continue;
    const double k=M*std::sqrt(lambda)/2;
    const double e1=(M2+s1-s2)/(2*M), e2=(M2+s2-s1)/(2*M);
    const double p1=std::sqrt(s1*b1)/2, p2=std::sqrt(s2*b2)/2;
    const double t1=std::sqrt(1-c1*c1), t2=std::sqrt(1-c2*c2);
    // An independent azimuth about the pair axis supplies the third overall
    // orientation angle in addition to the isotropic pair direction below.
    const double az=2*pi*flat();
    std::array<TLorentzVector,4> out = {
      TLorentzVector(p1*t1*std::cos(az),p1*t1*std::sin(az),p1*c1,std::sqrt(s1)/2),
      TLorentzVector(-p1*t1*std::cos(az),-p1*t1*std::sin(az),-p1*c1,std::sqrt(s1)/2),
      TLorentzVector(p2*t2*std::cos(az+phi),p2*t2*std::sin(az+phi),p2*c2,std::sqrt(s2)/2),
      TLorentzVector(-p2*t2*std::cos(az+phi),-p2*t2*std::sin(az+phi),-p2*c2,std::sqrt(s2)/2)};
    const double theta=std::acos(2*flat()-1), orientPhi=2*pi*flat();
    for (unsigned i=0;i<4;++i) {
      out[i].Boost(0,0,i<2 ? k/e1 : -k/e2);
      out[i].RotateY(theta); out[i].RotateZ(orientPhi);
      out[i].Boost(parent.BoostVector());
    }
    return out;
  }
  throw std::runtime_error("Mixed double-Dalitz rejection limit exhausted");
}
}
#endif
