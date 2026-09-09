#include "MixedDoubleDalitz.h"
#include <cassert>
#include <iostream>
#include <random>

int main() {
  // Independent integral of the angular polynomial: uniform cosines and phi.
  // Three-point Gauss-Legendre exactly integrates its cosine dependence.
  const double nodes[]={-std::sqrt(3./5),0,std::sqrt(3./5)};
  const double ws[]={5./18,4./9,5./18};
  for (double b1 : {0.,0.4,1.}) for (double b2 : {0.,0.7,1.}) {
    double avg=0;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) for (int k=0;k<8;++k) {
      const double w=cmspluto::mixedAngular(b1,b2,nodes[i],nodes[j],k*std::acos(-1.)/4);
      assert(w>=0 && w<=1);
      avg+=ws[i]*ws[j]*w/8;
    }
    assert(std::abs(avg-(1-b1/3)*(1-b2/3)) < 1e-14);
  }
  std::mt19937_64 rng(12345);
  auto flat=[&](){return std::generate_canonical<double,53>(rng);};
  for (double mass : {0.547862,0.95778}) for (double pz : {0.,5.}) {
    TLorentzVector parent(0,0,pz,std::hypot(pz,mass));
    for (int n=0;n<20;++n) {
      auto out=cmspluto::mixedPointlike(parent,0.105658,0.000511,flat);
      TLorentzVector sum;
      for (unsigned i=0;i<4;++i) {
        sum+=out[i];
        assert(std::abs(out[i].M()-(i<2 ? 0.105658 : 0.000511)) < 1e-8);
      }
      assert((sum-parent).Vect().Mag()+std::abs(sum.E()-parent.E()) < 1e-9);
      assert((out[0]+out[1]).M()>=2*0.105658);
      assert((out[2]+out[3]).M()>=2*0.000511);
    }
  }
  std::cout << "Mixed double Dalitz: angular integral, envelope, thresholds, masses and closure PASS\n";
}
