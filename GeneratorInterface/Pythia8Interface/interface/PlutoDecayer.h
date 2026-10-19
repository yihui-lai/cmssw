#ifndef GeneratorInterface_Pythia8Interface_PlutoDecayer_h
#define GeneratorInterface_Pythia8Interface_PlutoDecayer_h

#include <string>
#include <vector>
#include "Pythia8/Pythia.h"

// Forced-exclusive eta/eta' -> lepton(s)(+2pi/gamma) LO decayer, registered
// with Pythia8 through the classic DecayHandler/setDecayPtr mechanism (no
// Pythia version dependency, unlike the newer Init:plugins loader). The
// decay densities are constant-form-factor models from Petri,
// arXiv:1010.2378; see GeneratorInterface/Pythia8Interface/doc/PlutoDecayer.md
// for the physics assumptions. This does not depend on the external Pluto
// library: only ROOT's TLorentzVector is used.
class PlutoDecayer : public Pythia8::DecayHandler {
public:
  PlutoDecayer(Pythia8::Pythia* pythiaPtr, Pythia8::Settings* settingsPtr);

  bool decay(std::vector<int>& idProd,
             std::vector<double>& mProd,
             std::vector<Pythia8::Vec4>& pProd,
             int iDec,
             const Pythia8::Event& event) override;

private:
  Pythia8::Pythia* pythia_;
  bool filter_;
  int parentId_;
  std::string mode_;
  std::string model_;
  std::vector<int> daughterIds_;
  std::vector<double> daughterMasses_;
  double rhoMass_ = 0;
  double rhoGamma_ = 0;
  double rhopMass_ = 0;
  double rhopGamma_ = 0;
  double llpipiBound_ = 0;  // only set for mode=2e2pi/2mu2pi
  double tffBound_ = 0;     // rejection bound of the eta (Pade) or eta' (Pade+VMD) transition form
                            // factor; only set for mode=2mugamma/2egamma/2mu2e
  int aId_ = 0;             // PDG id of the BSM boson A (dark photon or ALP); only set for mode=2muA.
                            // A must already be a registered Pythia particle (mass/width/lifetime/decay
                            // table set by the fragment, e.g. via "<id>:all = ..."); this code only
                            // produces the eta/eta' -> mu+ mu- A kinematics, A's own decay is handled
                            // by Pythia's native decay machinery, not by PlutoDecayer.
  std::string aType_;       // "darkphoton" or "alp"; only set for mode=2muA -- see PlutoBSMDalitz.h
  double aBound_ = 0;       // rejection bound for the 2muA production weight; only set for mode=2muA
  double fourBound_ = 1.;   // sup|F(s1)F(s2)|^2 over the two dilepton invariant masses; only set for
                            // mode=4e/4mu -- see PlutoIdenticalDoubleDalitz.h's identicalPointlike
};

#endif
