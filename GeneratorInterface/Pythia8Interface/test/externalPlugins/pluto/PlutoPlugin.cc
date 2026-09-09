#include "Pythia8/Pythia.h"
#include "Pythia8/Plugins.h"
#include "PChannel.h"
#include "PDistributionManager.h"
#include "PParticle.h"
#include "PUtils.h"
#include "TGenPhaseSpace.h"
#include "MixedDoubleDalitz.h"
#include "PiPiDilepton.h"
#include "IdenticalDoubleDalitz.h"
#include <array>
#include <cmath>
#include <memory>
#include <mutex>
#include <stdexcept>

namespace {
  std::mutex plutoMutex;
  class PythiaRandom : public TRandom {
  public:
    explicit PythiaRandom(Pythia8::Rndm& rng) : rng_(rng) {}
    Double_t Rndm() override { return rng_.flat(); }
    void RndmArray(Int_t n, Double_t* out) override { for (int i = 0; i < n; ++i) out[i] = Rndm(); }
    void RndmArray(Int_t n, Float_t* out) override { for (int i = 0; i < n; ++i) out[i] = Rndm(); }
    void SetSeed(ULong_t = 0) override {} // Seed/state belong to CMSSW's engine.
  private:
    Pythia8::Rndm& rng_;
  };
  class RandomScope {
  public:
    explicit RandomScope(TRandom* rng) : oldRoot_(gRandom) {
      gRandom = rng;
      try { oldPluto_ = makePUtilsREngine()->SetExternalRandom(rng); }
      catch (...) { gRandom = oldRoot_; throw; }
    }
    ~RandomScope() { makePUtilsREngine()->SetExternalRandom(oldPluto_); gRandom = oldRoot_; }
  private:
    TRandom* oldRoot_;
    TRandom* oldPluto_;
  };
}

// Serial Pythia external-decay plugin. Forced exclusive samples only. The
// native model names describe the upstream eta implementations. Other models
// are explicitly named: phaseSpace controls or pointlike LO decay baselines.
class CMSPluto : public Pythia8::DecayHandler {
public:
  CMSPluto(Pythia8::Pythia* pythia, Pythia8::Settings* settings, Pythia8::Logger*)
      : pythia_(pythia), random_(pythia->rndm), parentId_(settings->mode("Pluto:parent")),
        mode_(settings->word("Pluto:mode")), model_(settings->word("Pluto:model")) {
    if (!settings->flag("Pluto:allowForcedDecay"))
      throw std::runtime_error("Set Pluto:allowForcedDecay = on for an exclusive forced sample; apply its branching normalization separately");
    if (parentId_ != 221 && parentId_ != 331) throw std::runtime_error("Pluto:parent must be 221 or 331");
    if (mode_ == "4e") ids_ = {11,-11,11,-11};
    else if (mode_ == "2mu2e") ids_ = {13,-13,11,-11};
    else if (mode_ == "4mu") ids_ = {13,-13,13,-13};
    else if (mode_ == "2e2pi") ids_ = {11,-11,211,-211};
    else if (mode_ == "2mu2pi") ids_ = {13,-13,211,-211};
    else throw std::runtime_error("Unknown Pluto:mode");
    for (unsigned i=0; i<4; ++i) masses_[i] = pythia_->particleData.m0(ids_[i]);
    const double threshold = masses_[0]+masses_[1]+masses_[2]+masses_[3];
    if (pythia_->particleData.m0(parentId_) <= threshold) throw std::runtime_error("Pluto channel is below threshold");
    if (model_ != "native" && model_ != "phaseSpace" && model_ != "mixedPointlike" && model_ != "pipiMagneticPointlike" && model_ != "identicalPointlike")
      throw std::runtime_error("Unknown Pluto:model");
    if (model_ == "identicalPointlike" && mode_ != "4e" && mode_ != "4mu")
      throw std::runtime_error("identicalPointlike requires 4e or 4mu");
    if (model_ == "pipiMagneticPointlike" && mode_ != "2e2pi" && mode_ != "2mu2pi")
      throw std::runtime_error("pipiMagneticPointlike requires 2e2pi or 2mu2pi");
    if (model_ == "mixedPointlike" && mode_ != "2mu2e")
      throw std::runtime_error("mixedPointlike is only for 2mu2e; identical-lepton exchange is not implemented");
    if (model_ == "native" && (parentId_ != 221 || (mode_ != "4e" && mode_ != "2e2pi")))
      throw std::runtime_error("No validated native Pluto model for this channel. phaseSpace is a kinematic control, not rare-decay dynamics");
    if (model_ == "native" && mode_ == "2e2pi")
      throw std::runtime_error("Native eta->ee pi pi is quarantined due to upstream sampling and angular defects. Select pipiMagneticPointlike for the explicit magnetic baseline or phaseSpace for a kinematic control");
    std::lock_guard<std::mutex> lock(plutoMutex);
    RandomScope scope(&random_);
    if (model_ == "native") initializeNative();
  }

  std::vector<int> handledParticles() override { return {parentId_}; }

  bool decay(std::vector<int>& ids, std::vector<double>& masses, std::vector<Pythia8::Vec4>& momenta,
             int, const Pythia8::Event&) override {
    if (ids.size()!=1 || masses.size()!=1 || momenta.size()!=1 || ids.front()!=parentId_)
      throw std::runtime_error("Invalid Pythia input to CMSPluto");
    std::lock_guard<std::mutex> lock(plutoMutex);
    RandomScope scope(&random_);
    const auto& p = momenta.front();
    TLorentzVector parent(p.px(),p.py(),p.pz(),p.e());
    const double threshold = masses_[0]+masses_[1]+masses_[2]+masses_[3];
    if (parent.M() <= threshold) throw std::runtime_error("Actual parent mass is below decay threshold");
    std::array<TLorentzVector,4> result;
    if (model_ == "identicalPointlike") {
      result = cmspluto::identicalPointlike(parent, masses_[0], [&]() { return random_.Rndm(); });
    } else if (model_ == "pipiMagneticPointlike") {
      result = cmspluto::pipiMagneticPointlike(parent, masses_[0], masses_[2], [&]() { return random_.Rndm(); });
    } else if (model_ == "mixedPointlike") {
      result = cmspluto::mixedPointlike(parent, masses_[0], masses_[2], [&]() { return random_.Rndm(); });
    } else if (model_ == "phaseSpace") {
      TGenPhaseSpace generator;
      if (!generator.SetDecay(parent,4,masses_.data())) throw std::runtime_error("Phase-space initialization failed");
      bool accepted = false;
      for (unsigned attempt=0; attempt<1000000; ++attempt) {
        const double weight = generator.Generate();
        if (!std::isfinite(weight) || weight < 0 || weight > 1.) throw std::runtime_error("Invalid phase-space rejection weight");
        if (random_.Rndm() < weight) { accepted = true; break; }
      }
      if (!accepted) throw std::runtime_error("Phase-space generation exhausted its attempt limit");
      for (unsigned i=0; i<4; ++i) result[i] = *generator.GetDecay(i);
    } else {
      bool accepted = false;
      for (unsigned attempt=0; attempt<1000000; ++attempt) {
        parent_->SetPxPyPzE(0,0,0,parent.M());
        for (auto& particle : particles_) { particle->SetActive(); particle->SetW(1.); }
        bool good = true;
        for (auto& channel : channels_) {
          if (channel->Decay() != 0) { good = false; break; }
        }
        if (!good) continue;
        for (auto& channel : channels_)
          for (int i=0; i<channel->GetNumNotFinalized(); ++i)
            if (!channel->GetDistributionNotFinalized(i)->EndOfChain()) good = false;
        if (good) { accepted = true; break; }
      }
      if (!accepted) throw std::runtime_error("Native Pluto decay exhausted its attempt limit");
      for (unsigned i=0; i<4; ++i) { result[i] = *final_[i]; result[i].Boost(parent.BoostVector()); }
    }
    TLorentzVector sum;
    for (const auto& daughter : result) sum += daughter;
    const auto residual = sum-parent;
    if (std::abs(residual.E())+residual.Vect().Mag() > 1e-7*std::max(1.,parent.E()))
      throw std::runtime_error("Pluto decay failed four-momentum closure");
    for (unsigned i=0; i<4; ++i) {
      ids.push_back(ids_[i]); masses.push_back(result[i].M());
      momenta.emplace_back(result[i].Px(),result[i].Py(),result[i].Pz(),result[i].E());
    }
    return true;
  }

private:
  PParticle* particle(const char* name) {
    particles_.push_back(std::make_unique<PParticle>(name));
    return particles_.back().get();
  }
  void channel(std::initializer_list<PParticle*> particles) {
    auto array = std::make_unique<PParticle*[]>(particles.size());
    std::copy(particles.begin(),particles.end(),array.get());
    channels_.push_back(std::make_unique<PChannel>(array.get(),particles.size()-1));
    arrays_.push_back(std::move(array));
  }
  void initializeNative() {
    makeDistributionManager()->Exec("eta_decays");
    parent_ = particle("eta");
    final_[0] = particle("e-"); final_[1] = particle("e+");
    final_[2] = particle(mode_ == "4e" ? "e-" : "pi+");
    final_[3] = particle(mode_ == "4e" ? "e+" : "pi-");
    if (mode_ == "4e") {
      auto* pair1 = particle("dilepton"); auto* pair2 = particle("dilepton");
      channel({parent_,pair1,pair2});
      channel({pair1,final_[1],final_[0]});
      channel({pair2,final_[3],final_[2]});
    } else channel({parent_,final_[0],final_[1],final_[2],final_[3]});
    // Match PReaction's chain setup: attach models, expose all descendants,
    // then resolve tentative (granddaughter-dependent) distributions.
    for (auto& c : channels_) makeDistributionManager()->Attach(c.get());
    for (auto& c : channels_) c->SetDaughters();
    for (auto& c : channels_) c->Init();
  }
  Pythia8::Pythia* pythia_;
  PythiaRandom random_;
  int parentId_;
  std::string mode_, model_;
  std::array<int,4> ids_;
  std::array<double,4> masses_;
  PParticle* parent_ = nullptr;
  std::array<PParticle*,4> final_{};
  std::vector<std::unique_ptr<PParticle>> particles_;
  std::vector<std::unique_ptr<PParticle*[]>> arrays_;
  std::vector<std::unique_ptr<PChannel>> channels_;
};

void registerPluto(Pythia8::Settings* s) {
  s->addMode("Pluto:parent",221,true,true,221,331);
  s->addWord("Pluto:mode","4e");
  s->addWord("Pluto:model","native");
  s->addFlag("Pluto:allowForcedDecay",false);
}
using namespace Pythia8; // Required by Pythia's plugin registration macros.
PYTHIA8_PLUGIN_CLASS(DecayHandler, CMSPluto, true, true, false)
PYTHIA8_PLUGIN_SETTINGS(registerPluto)
PYTHIA8_PLUGIN_VERSIONS(PYTHIA_VERSION_INTEGER)
PYTHIA8_PLUGIN_PARALLEL(false)
