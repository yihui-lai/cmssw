#include "IdenticalDoubleDalitz.h"
#include <cassert>
#include <iostream>
#include <random>

double directAnalytic(std::array<TLorentzVector,4> p,double m) {
  TLorentzVector parent;
  for (const auto& v:p) parent+=v;
  for (auto& v:p) v.Boost(-parent.BoostVector());
  const auto q1=p[0]+p[1],q2=p[2]+p[3];
  const double s1=q1.M2(),s2=q2.M2(),M2=parent.M2();
  const auto axis=q1.Vect().Unit();
  p[0].Boost(-q1.BoostVector()); p[2].Boost(-q2.BoostVector());
  const auto n1=p[0].Vect().Unit(),n2=p[2].Vect().Unit();
  const double c1=n1.Dot(axis),c2=n2.Dot(axis),triple=n1.Cross(n2).Dot(axis);
  const double b1=1-4*m*m/s1,b2=1-4*m*m/s2;
  const double lambda=std::pow(M2-s1-s2,2)-4*s1*s2;
  return lambda/(s1*s2)*(2-b1*(1-c1*c1)-b2*(1-c2*c2)+b1*b2*triple*triple);
}
bool close(double a,double b) { return std::abs(a-b)<1e-7*std::max({1.,std::abs(a),std::abs(b)}); }
int main() {
  std::mt19937_64 rng(9876);
  auto flat=[&](){return std::generate_canonical<double,53>(rng);};
  double worst=0;
  bool nonzeroInterference=false;
  for (double M:{0.547862,0.95778}) for (double m:{0.000511,0.105658}) {
    const TLorentzVector parent(0,0,0,M);
    for (int n=0;n<12;++n) {
      auto p=cmspluto::mixedPointlike(parent,m,m,flat);
      const auto score=cmspluto::fourlepton::scores(p,m);
      const double analytic=directAnalytic(p,m);
      worst=std::max(worst,std::abs(score.direct-analytic)/std::max(1.,std::abs(analytic)));
      assert(close(score.direct,analytic));
      assert(score.coherent>=0 && score.coherent<=2*(score.direct+score.exchanged)*(1+1e-12));
      if (std::abs(score.coherent-score.direct-score.exchanged)>1e-5*(score.direct+score.exchanged))
        nonzeroInterference=true;
      std::swap(p[1],p[3]);
      const auto exchanged=cmspluto::fourlepton::scores(p,m);
      assert(close(score.direct,exchanged.exchanged));
      assert(close(score.exchanged,exchanged.direct));
      assert(close(score.coherent,exchanged.coherent));
      std::swap(p[0],p[2]);
      assert(close(score.coherent,cmspluto::fourlepton::scores(p,m).coherent));
      for (auto& v:p) v.Boost(0.2,0.1,0.3);
      assert(close(score.coherent,cmspluto::fourlepton::scores(p,m).coherent));
      const TLorentzVector boostedParent(0,0,5,std::hypot(5.,M));
      const auto generated=cmspluto::identicalPointlike(boostedParent,m,flat);
      TLorentzVector sum;
      for (const auto& v:generated) { sum+=v; assert(std::abs(v.M()-m)<1e-8); }
      assert((sum-boostedParent).Vect().Mag()+std::abs(sum.E()-boostedParent.E())<1e-9);
    }
  }
  assert(nonzeroInterference);
  // Same-momentum, same-spin negative fermions must cancel at amplitude level.
  const double m=0.105658,p=0.1,t=0.03;
  const TLorentzVector negative(p,0,0,std::hypot(p,m));
  const std::array<TLorentzVector,4> pauli={negative,
    TLorentzVector(-p,t,0,std::sqrt(p*p+t*t+m*m)),negative,
    TLorentzVector(-p,-t,0,std::sqrt(p*p+t*t+m*m))};
  const cmspluto::fourlepton::Amplitudes amplitudes(pauli,m);
  for (unsigned s=0;s<2;++s) for (unsigned a=0;a<2;++a) for (unsigned b=0;b<2;++b) {
    const auto value=amplitudes(s,a,s,b);
    assert(std::norm(value[0]-value[1])<1e-20*std::max(1.,std::norm(value[0])+std::norm(value[1])));
  }
  std::cout << "Identical double Dalitz: direct analytic agreement, exchange/boost symmetry, Pauli cancellation, envelope and closure PASS; worst direct relative discrepancy=" << worst << '\n';
}
