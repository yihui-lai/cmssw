#include "PiPiDilepton.h"
#include <cassert>
#include <iostream>
#include <random>

int main() {
  const double pi=std::acos(-1.);
  const double nodes[]={-std::sqrt(3./5),0,std::sqrt(3./5)}, ws[]={5./18,4./9,5./18};
  for (double beta2 : {0.,0.5,1.}) {
    double avg=0;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) for (int k=0;k<8;++k) {
      const double w=cmspluto::pipiMagneticAngular(beta2,nodes[i],nodes[j],k*pi/4);
      assert(w>=0 && w<=1);
      avg+=ws[i]*ws[j]*w/8;
    }
    assert(std::abs(avg-(2./3)*(1-beta2/3))<1e-14);
  }
  // A massless lepton along the magnetic-current direction has zero weight.
  assert(std::abs(cmspluto::pipiMagneticAngular(1,0,0,pi/2))<1e-14);
  assert(cmspluto::pipiMagneticAngular(1,0,0,0)==1);
  std::mt19937_64 rng(12345);
  auto flat=[&](){return std::generate_canonical<double,53>(rng);};
  for (double mass : {0.547862,0.95778}) for (double ml : {0.000511,0.105658})
    for (double pz : {0.,5.}) {
      TLorentzVector parent(0,0,pz,std::hypot(pz,mass));
      for (int n=0;n<10;++n) {
        auto out=cmspluto::pipiMagneticPointlike(parent,ml,0.139570,flat);
        TLorentzVector sum;
        for (unsigned i=0;i<4;++i) {
          sum+=out[i];
          assert(std::abs(out[i].M()-(i<2 ? ml : 0.139570))<1e-8);
        }
        assert((sum-parent).Vect().Mag()+std::abs(sum.E()-parent.E())<1e-9);
        assert((out[0]+out[1]).M()>=2*ml);
        assert((out[2]+out[3]).M()>=2*0.139570);
      }
    }
  std::cout << "PiPi magnetic: angular integral and sign, fixed envelope, thresholds, masses and closure PASS\n";
}
