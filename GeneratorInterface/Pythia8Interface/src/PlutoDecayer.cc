#include "GeneratorInterface/Pythia8Interface/interface/PlutoDecayer.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoMixedDoubleDalitz.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoPiPiDilepton.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoIdenticalDoubleDalitz.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoSingleDalitz.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaPrimeLLPiPi.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoLLPiPiChPT.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaTFF.h"

#include "TGenPhaseSpace.h"
#include "TRandom.h"

#include <cmath>
#include <mutex>
#include <stdexcept>

using namespace Pythia8;

namespace {
  // Only the phaseSpace control model touches ROOT's global gRandom (via
  // TGenPhaseSpace). Route it through Pythia's own, CMSSW-managed random
  // stream for reproducibility, and guard the global with a mutex since
  // multiple streams/threads may run PlutoDecayer concurrently.
  std::mutex gRandomMutex;
  class PythiaRandom : public TRandom {
  public:
    explicit PythiaRandom(Pythia8::Rndm& rng) : rng_(rng) {}
    Double_t Rndm() override { return rng_.flat(); }
    void RndmArray(Int_t n, Double_t* out) override {
      for (int i = 0; i < n; ++i)
        out[i] = Rndm();
    }
    void RndmArray(Int_t n, Float_t* out) override {
      for (int i = 0; i < n; ++i)
        out[i] = Rndm();
    }
    void SetSeed(ULong_t = 0) override {}  // Seed/state belong to CMSSW's engine.
  private:
    Pythia8::Rndm& rng_;
  };
  class RandomScope {
  public:
    explicit RandomScope(TRandom* rng) : old_(gRandom) { gRandom = rng; }
    ~RandomScope() { gRandom = old_; }

  private:
    TRandom* old_;
  };

  template <std::size_t N>
  std::vector<TLorentzVector> toVector(const std::array<TLorentzVector, N>& in) {
    return std::vector<TLorentzVector>(in.begin(), in.end());
  }
}  // namespace

PlutoDecayer::PlutoDecayer(Pythia8::Pythia* pythiaPtr, Pythia8::Settings* settingsPtr) : pythia_(pythiaPtr) {
  filter_ = settingsPtr->flag("Pluto:filter");
  parentId_ = settingsPtr->mode("Pluto:parent");
  mode_ = settingsPtr->word("Pluto:mode");
  model_ = settingsPtr->word("Pluto:model");

  if (!filter_)
    return;

  if (!settingsPtr->flag("Pluto:allowForcedDecay"))
    throw std::runtime_error(
        "Set Pluto:allowForcedDecay = on to acknowledge this is a forced "
        "exclusive sample, not a branching-ratio-normalized one; apply the "
        "branching normalization separately.");
  if (parentId_ != 221 && parentId_ != 331)
    throw std::runtime_error("Pluto:parent must be 221 (eta) or 331 (eta')");
  if (mode_ == "2mu")
    daughterIds_ = {13, -13};
  else if (mode_ == "2mugamma")
    daughterIds_ = {13, -13, 22};
  else if (mode_ == "4e")
    daughterIds_ = {11, -11, 11, -11};
  else if (mode_ == "2mu2e")
    daughterIds_ = {13, -13, 11, -11};
  else if (mode_ == "4mu")
    daughterIds_ = {13, -13, 13, -13};
  else if (mode_ == "2e2pi")
    daughterIds_ = {11, -11, 211, -211};
  else if (mode_ == "2mu2pi")
    daughterIds_ = {13, -13, 211, -211};
  else
    throw std::runtime_error("Unknown Pluto:mode \"" + mode_ + "\"");
  if (model_ != "pointlike" && model_ != "phaseSpace")
    throw std::runtime_error("Unknown Pluto:model \"" + model_ + "\" (use pointlike or phaseSpace)");

  double threshold = 0;
  for (int id : daughterIds_) {
    const double m = pythia_->particleData.m0(id);
    daughterMasses_.push_back(m);
    threshold += m;
  }
  if (pythia_->particleData.m0(parentId_) <= threshold)
    throw std::runtime_error("Pluto channel is below threshold");

  rhoMass_ = pythia_->particleData.m0(113);
  rhoGamma_ = pythia_->particleData.mWidth(113);
  rhopMass_ = pythia_->particleData.m0(100113);
  rhopGamma_ = pythia_->particleData.mWidth(100113);
  if (model_ == "pointlike" && (mode_ == "2e2pi" || mode_ == "2mu2pi"))
    llpipiBound_ = gen::pluto::llpipiWeightBound(pythia_->particleData.m0(parentId_),
                                                  daughterMasses_[2],
                                                  daughterMasses_[0],
                                                  parentId_ == 331,
                                                  rhoMass_,
                                                  rhoGamma_,
                                                  rhopMass_,
                                                  rhopGamma_);
  if (model_ == "pointlike" && parentId_ == 221 && mode_ == "2mugamma")
    etaPadeBound_ = gen::pluto::etaPadeSingleBound(pythia_->particleData.m0(221), daughterMasses_[0]);
  else if (model_ == "pointlike" && parentId_ == 221 && mode_ == "2mu2e")
    etaPadeBound_ = gen::pluto::etaPadeDoubleBound(
        pythia_->particleData.m0(221), daughterMasses_[0], daughterMasses_[2]);
}

bool PlutoDecayer::decay(std::vector<int>& idProd,
                          std::vector<double>& mProd,
                          std::vector<Pythia8::Vec4>& pProd,
                          int /*iDec*/,
                          const Pythia8::Event& /*event*/) {
  if (!filter_)
    return false;
  if (idProd.size() != 1 || mProd.size() != 1 || pProd.size() != 1 || idProd.front() != parentId_)
    return false;

  const auto& p = pProd.front();
  TLorentzVector parent(p.px(), p.py(), p.pz(), p.e());
  double threshold = 0;
  for (double m : daughterMasses_)
    threshold += m;
  if (parent.M() <= threshold)
    throw std::runtime_error("Actual Pluto parent mass is below decay threshold");

  auto flat = [this]() { return pythia_->rndm.flat(); };
  std::vector<TLorentzVector> result;
  if (model_ == "phaseSpace") {
    std::lock_guard<std::mutex> lock(gRandomMutex);
    PythiaRandom rng(pythia_->rndm);
    RandomScope scope(&rng);
    TGenPhaseSpace generator;
    if (!generator.SetDecay(parent, static_cast<int>(daughterMasses_.size()), daughterMasses_.data()))
      throw std::runtime_error("Pluto phase-space initialization failed");
    bool accepted = false;
    for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
      const double weight = generator.Generate();
      if (!std::isfinite(weight) || weight < 0 || weight > 1.)
        throw std::runtime_error("Invalid Pluto phase-space rejection weight");
      if (flat() < weight) {
        accepted = true;
        break;
      }
    }
    if (!accepted)
      throw std::runtime_error("Pluto phase-space generation exhausted its attempt limit");
    for (unsigned i = 0; i < daughterMasses_.size(); ++i)
      result.push_back(*generator.GetDecay(i));
  } else if (mode_ == "2mu") {
    // Trivial isotropic two-body decay; no rejection sampling needed.
    const double M = parent.M(), m = daughterMasses_[0];
    const double pMag = std::sqrt(std::max(0., M * M / 4 - m * m));
    const double c = 2 * flat() - 1, s = std::sqrt(1 - c * c), pi = std::acos(-1.), phi = 2 * pi * flat();
    TLorentzVector lMinus(pMag * s * std::cos(phi), pMag * s * std::sin(phi), pMag * c, M / 2);
    TLorentzVector lPlus(-lMinus.Px(), -lMinus.Py(), -lMinus.Pz(), M / 2);
    lMinus.Boost(parent.BoostVector());
    lPlus.Boost(parent.BoostVector());
    result = {lMinus, lPlus};
  } else if (mode_ == "2mugamma") {
    // eta: actual data-fitted eta TFF (arXiv:1504.07742 Appendix A), not a
    // generic resonance guess. eta': rho0-pole VMD (matches EvtGen), still
    // missing the omega contribution documented as a known gap.
    result = parentId_ == 221
                 ? toVector(gen::pluto::singleDalitzEtaPade(parent, daughterMasses_[0], etaPadeBound_, flat))
                 : toVector(gen::pluto::singleDalitz(parent, daughterMasses_[0], rhoMass_, rhoGamma_, flat));
  } else if (mode_ == "4e" || mode_ == "4mu") {
    // Amplitude-level form factor dressing of each diagram's two virtual
    // photons (direct and exchange evaluated at their own, different,
    // invariant masses) before they interfere -- see PlutoIdenticalDoubleDalitz.h.
    // eta uses the data-fitted eta TFF; eta' uses the rho0-pole model.
    result = toVector(gen::pluto::identicalPointlike(
        parent, daughterMasses_[0], flat, rhoMass_, rhoGamma_, parentId_ == 221));
  } else if (mode_ == "2mu2e") {
    // Factorized double-virtual form factor (arXiv:1511.04916 Eq. 8); no
    // exchange interference exists here, so this reweighting is unambiguous.
    // eta uses the data-fitted eta TFF (PlutoEtaTFF.h); eta' the rho0-pole.
    result = parentId_ == 221
                 ? toVector(gen::pluto::mixedPointlikeEtaPade(
                       parent, daughterMasses_[0], daughterMasses_[2], etaPadeBound_, flat))
                 : toVector(gen::pluto::mixedPointlikeResonant(
                       parent, daughterMasses_[0], daughterMasses_[2], rhoMass_, rhoGamma_, flat));
  } else {
    // eta/eta' 2e2pi/2mu2pi: mass-dependence from Zillinger/Kubis/Sanchez-Puertas,
    // arXiv:2210.14925 (P(s)*Omega(s)*Fbar(s_l), Omega approximated -- see
    // PlutoLLPiPiChPT.h), combined with pipiMagneticPointlike's angular shape.
    result = toVector(gen::pluto::llpipiChPT(parent,
                                              daughterMasses_[0],
                                              daughterMasses_[2],
                                              parentId_ == 331,
                                              rhoMass_,
                                              rhoGamma_,
                                              rhopMass_,
                                              rhopGamma_,
                                              llpipiBound_,
                                              flat));
  }

  TLorentzVector sum;
  for (const auto& daughter : result)
    sum += daughter;
  const auto residual = sum - parent;
  if (std::abs(residual.E()) + residual.Vect().Mag() > 1e-7 * std::max(1., parent.E()))
    throw std::runtime_error("Pluto decay failed four-momentum closure");

  for (unsigned i = 0; i < daughterIds_.size(); ++i) {
    idProd.push_back(daughterIds_[i]);
    mProd.push_back(result[i].M());
    pProd.emplace_back(result[i].Px(), result[i].Py(), result[i].Pz(), result[i].E());
  }
  return true;
}
