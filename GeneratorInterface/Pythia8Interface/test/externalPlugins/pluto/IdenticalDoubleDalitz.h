#ifndef CMS_PLUTO_IDENTICAL_DOUBLE_DALITZ_H
#define CMS_PLUTO_IDENTICAL_DOUBLE_DALITZ_H
#include "MixedDoubleDalitz.h"
#include <complex>

namespace cmspluto {
namespace fourlepton {
using Complex=std::complex<double>;
using Spinor=std::array<Complex,4>;
using Current=std::array<Complex,3>;

// Dirac basis, ubar*u=2m and vbar*v=-2m. Both independent spin states are
// summed explicitly; the choice of antiparticle spin basis does not affect it.
inline Spinor spinor(const TLorentzVector& p,double m,unsigned spin,bool antiparticle) {
  const double root=std::sqrt(p.E()+m);
  Spinor out{};
  const Complex a=spin==0 ? Complex(p.Pz(),0) : Complex(p.Px(),-p.Py());
  const Complex b=spin==0 ? Complex(p.Px(),p.Py()) : Complex(-p.Pz(),0);
  if (antiparticle) { out[0]=a/root; out[1]=b/root; out[2+spin]=root; }
  else { out[spin]=root; out[2]=a/root; out[3]=b/root; }
  return out;
}

inline Current current(const Spinor& u,const Spinor& v) {
  const Complex a=std::conj(u[0]), b=std::conj(u[1]);
  const Complex c=std::conj(u[2]), d=std::conj(u[3]), I(0,1);
  return {a*v[3]+b*v[2]+c*v[1]+d*v[0],
          I*(-a*v[3]+b*v[2]-c*v[1]+d*v[0]),
          a*v[2]-b*v[3]+c*v[0]-d*v[1]};
}

inline Complex contraction(const Current& a,const Current& b,const TVector3& q) {
  return q.X()*(a[1]*b[2]-a[2]*b[1])
        +q.Y()*(a[2]*b[0]-a[0]*b[2])
        +q.Z()*(a[0]*b[1]-a[1]*b[0]);
}

// Rest-frame epsilon contraction of the two conserved vector currents.
// Ordering is l-(0), l+(1), l-(2), l+(3). The physical amplitude is D-X.
// Global e^2 and the constant pseudoscalar form factor cancel in the sampler.
class Amplitudes {
public:
  Amplitudes(std::array<TLorentzVector,4> p,double m) {
    TLorentzVector parent;
    for (const auto& daughter:p) parent+=daughter;
    const double M=parent.M();
    if (!(M>4*m && m>0)) throw std::runtime_error("Invalid four-lepton masses");
    for (auto& daughter:p) daughter.Boost(-parent.BoostVector());
    const auto qD=p[0]+p[1], rD=p[2]+p[3];
    const auto qX=p[0]+p[3], rX=p[2]+p[1];
    dScale_=M/(qD.M2()*rD.M2()); xScale_=M/(qX.M2()*rX.M2());
    dAxis_=qD.Vect(); xAxis_=qX.Vect();
    std::array<std::array<Spinor,2>,4> spinors;
    for (unsigned i=0;i<4;++i) for (unsigned s=0;s<2;++s)
      spinors[i][s]=spinor(p[i],m,s,i%2==1);
    for (unsigned s=0;s<2;++s) for (unsigned t=0;t<2;++t) {
      d1_[s][t]=current(spinors[0][s],spinors[1][t]);
      d2_[s][t]=current(spinors[2][s],spinors[3][t]);
      x1_[s][t]=current(spinors[0][s],spinors[3][t]);
      x2_[s][t]=current(spinors[2][s],spinors[1][t]);
    }
  }
  std::array<Complex,2> operator()(unsigned a,unsigned b,unsigned c,unsigned d) const {
    return {dScale_*contraction(d1_[a][b],d2_[c][d],dAxis_),
            xScale_*contraction(x1_[a][d],x2_[c][b],xAxis_)};
  }
private:
  using Currents=std::array<std::array<Current,2>,2>;
  Currents d1_,d2_,x1_,x2_;
  TVector3 dAxis_,xAxis_;
  double dScale_,xScale_;
};

struct Scores { double direct=0, exchanged=0, coherent=0; };
inline Scores scores(const std::array<TLorentzVector,4>& p,double m) {
  const Amplitudes amplitudes(p,m);
  Scores result;
  for (unsigned a=0;a<2;++a) for (unsigned b=0;b<2;++b)
    for (unsigned c=0;c<2;++c) for (unsigned d=0;d<2;++d) {
      const auto values=amplitudes(a,b,c,d);
      result.direct+=std::norm(values[0]); result.exchanged+=std::norm(values[1]);
      result.coherent+=std::norm(values[0]-values[1]);
    }
  return result;
}
}

// Equal mixture of direct and exchanged pairings followed by bounded
// coherent rejection. |D-X|^2 <= 2(|D|^2+|X|^2) for the spin sum, so no
// empirically tuned or growing envelope is needed. Identical-particle
// factorials affect rates, not this normalized exclusive shape.
template<class Flat>
std::array<TLorentzVector,4> identicalPointlike(const TLorentzVector& parent,double m,Flat flat) {
  const TLorentzVector rest(0,0,0,parent.M());
  for (unsigned attempt=0;attempt<100000;++attempt) {
    auto out=mixedPointlike(rest,m,m,flat);
    if (flat()<0.5) std::swap(out[1],out[3]);
    const auto score=fourlepton::scores(out,m);
    const double denominator=2*(score.direct+score.exchanged);
    if (!(denominator>0) || !std::isfinite(denominator))
      throw std::runtime_error("Invalid four-lepton proposal density");
    const double weight=score.coherent/denominator;
    if (!std::isfinite(weight) || weight<0 || weight>1+1e-12)
      throw std::runtime_error("Four-lepton coherent envelope violated");
    if (flat()>=std::min(1.,weight)) continue;
    for (auto& p:out) p.Boost(parent.BoostVector());
    return out;
  }
  throw std::runtime_error("Identical double-Dalitz rejection limit exhausted");
}
}
#endif
