// Cheap smoke test for the BSM samplers: four-momentum closure (parent =
// sum of daughters) and A's mass being exactly the configured value (not
// sampled -- A is always on-shell by construction). Run after any physics
// change to the BSM headers as a fast sanity check before the full
// bsm_cms_gen.cc shape validation.
#include "GeneratorInterface/Pythia8Interface/interface/EtaBSMDalitz.h"
#include "GeneratorInterface/Pythia8Interface/interface/EtaBSMPiPiA.h"
#include <cstdio>
#include <random>

using namespace gen::pluto;

namespace {
const double M = 0.95778, MA = 0.400115, MMU = 0.105658, MPI = 0.13957;
const double RHO_M = 0.7754, RHO_G = 0.1462, RHOP_M = 1.465, RHOP_G = 0.400;
}  // namespace

int main() {
  std::mt19937_64 rng(1);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  TLorentzVector parent(0, 0, 0, M);

  {
    double bound = darkPhotonSingleBound(M, MMU, MA, true);
    double maxClose = 0, maxMassDiff = 0;
    for (int i = 0; i < 20000; i++) {
      auto out = singleDalitzDarkPhoton(parent, MMU, MA, true, bound, flat);
      TLorentzVector sum = out[0] + out[1] + out[2];
      maxClose = std::max(maxClose, std::abs((sum - parent).E()));
      maxMassDiff = std::max(maxMassDiff, std::abs(out[2].M() - MA));
    }
    printf("llA darkphoton: max|E_sum-M|=%.3e  max|m_A sampled - set|=%.3e\n", maxClose, maxMassDiff);
  }
  {
    double bound = darkPhotonPiPiABound(M, MPI, MA, MMU, true, RHO_M, RHO_G, RHOP_M, RHOP_G);
    double maxClose = 0, maxMassDiff = 0;
    for (int i = 0; i < 20000; i++) {
      auto out = piPiADarkPhoton(parent, MPI, MA, MMU, true, RHO_M, RHO_G, RHOP_M, RHOP_G, bound, flat);
      TLorentzVector sum = out[0] + out[1] + out[2];
      maxClose = std::max(maxClose, std::abs((sum - parent).E()));
      maxMassDiff = std::max(maxMassDiff, std::abs(out[2].M() - MA));
    }
    printf("pipiA darkphoton: max|E_sum-M|=%.3e  max|m_A sampled - set|=%.3e\n", maxClose, maxMassDiff);
  }
  return 0;
}
