#ifndef CMS_PLUTO_PIPI_DILEPTON_H
#define CMS_PLUTO_PIPI_DILEPTON_H
#include "TLorentzVector.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace cmspluto {
// Magnetic term of Petri arXiv:1010.2378, section 3.5.2, Eq. (3.85).
// Electric/CP-violating amplitudes are absent; M(s_pipi,s_ll) is constant.
inline double pipiMagneticAngular(double betaL2, double cosPi, double cosL, double phi) {
  return (1-cosPi*cosPi)*(1-betaL2*(1-cosL*cosL)*std::pow(std::sin(phi),2));
}

template<class Flat>
std::array<TLorentzVector,4> pipiMagneticPointlike(const TLorentzVector& parent,
                                                double ml, double mpi, Flat flat) {
  const double M=parent.M(), M2=M*M, smin=4*mpi*mpi, qmin=4*ml*ml;
  if (!(ml>0 && mpi>0 && M>2*(ml+mpi)))
    throw std::runtime_error("PiPi dilepton masses are below threshold");
  const double smax=std::pow(M-2*ml,2), qmax=std::pow(M-2*mpi,2);
  const double logq=std::log(qmax/qmin), pi=std::acos(-1.);
  const auto lambda=[](double x,double y){return std::max(0.,std::pow(1-x-y,2)-4*x*y);};
  // Each factor has a monotonic bound over the physical domain. Dividing by
  // their product avoids inefficient near-threshold eta -> mumu pi pi trials
  // while retaining a rigorous fixed envelope, independent of sampled events.
  const double bound=std::pow(lambda(smin/M2,qmin/M2),1.5)
    *(smax/M2)*std::pow(1-smin/smax,1.5)*std::sqrt(1-qmin/qmax);
  if (!(bound>0) || !std::isfinite(bound)) throw std::runtime_error("Invalid PiPi envelope");
  for (unsigned attempt=0;attempt<1000000;++attempt) {
    const double s=smin+(smax-smin)*flat(), q=qmin*std::exp(logq*flat());
    if (std::sqrt(s)+std::sqrt(q)>=M) continue;
    const double bp=1-smin/s, bl=1-qmin/q, lam=lambda(s/M2,q/M2);
    const double cp=2*flat()-1, cl=2*flat()-1, phi=2*pi*flat();
    // Phase space times magnetic |A|^2: lambda^(3/2) s beta_pi^3
    // beta_l angular / q. The logarithmic q proposal cancels 1/q.
    const double weight=std::pow(lam,1.5)*(s/M2)*std::pow(bp,1.5)*std::sqrt(bl)
      *pipiMagneticAngular(bl,cp,cl,phi)/bound;
    if (!std::isfinite(weight) || weight<0 || weight>1)
      throw std::runtime_error("PiPi magnetic rejection envelope violated");
    if (flat()>=weight) continue;
    const double k=M*std::sqrt(lam)/2;
    const double el=(M2+q-s)/(2*M), ep=(M2+s-q)/(2*M);
    const double pl=std::sqrt(q*bl)/2, pp=std::sqrt(s*bp)/2;
    const double tl=std::sqrt(1-cl*cl), tp=std::sqrt(1-cp*cp), az=2*pi*flat();
    // Output ordering is l-, l+, pi+, pi- as used by the external adapter.
    std::array<TLorentzVector,4> out={
      TLorentzVector(pl*tl*std::cos(az+phi),pl*tl*std::sin(az+phi),pl*cl,std::sqrt(q)/2),
      TLorentzVector(-pl*tl*std::cos(az+phi),-pl*tl*std::sin(az+phi),-pl*cl,std::sqrt(q)/2),
      TLorentzVector(pp*tp*std::cos(az),pp*tp*std::sin(az),pp*cp,std::sqrt(s)/2),
      TLorentzVector(-pp*tp*std::cos(az),-pp*tp*std::sin(az),-pp*cp,std::sqrt(s)/2)};
    const double theta=std::acos(2*flat()-1), orientPhi=2*pi*flat();
    for (unsigned i=0;i<4;++i) {
      out[i].Boost(0,0,i<2 ? k/el : -k/ep);
      out[i].RotateY(theta); out[i].RotateZ(orientPhi);
      out[i].Boost(parent.BoostVector());
    }
    return out;
  }
  throw std::runtime_error("PiPi dilepton rejection limit exhausted");
}
}
#endif
