// Generates the MC side of the PlutoDecayer-vs-paper validation.
//
// Draws events straight from the samplers in
// GeneratorInterface/Pythia8Interface/interface/ -- i.e. whatever physics is
// CURRENTLY wired into PlutoDecayer -- and histograms the same observables
// as Figs. 4/8/9/10/11 of arXiv:1511.04916 and Fig. 7 of arXiv:0705.0954.
// No CMSSW build or cmsRun needed (header-only + ROOT's TLorentzVector).
//
// Models under test:
//   eta  2mugamma/2mu2e/4e/4mu : data-fitted Pade TFF (arXiv:1504.07742)
//   eta' 2mugamma/2mu2e        : Pade + rho/omega/phi VMD TFF (arXiv:1307.2061, 1511.04916)
//   eta' 4e/4mu                : rho-pole model
//   eta/eta' 2e2pi/2mu2pi      : llpipiChPT (arXiv:2210.14925)
//
// Usage: ./compare_papers [output.csv]   (default: pluto_mc.csv)
#include "GeneratorInterface/Pythia8Interface/interface/PlutoMixedDoubleDalitz.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoIdenticalDoubleDalitz.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaTFF.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoLLPiPiChPT.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoSingleDalitz.h"

#include <array>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace {
const double ME = 0.000511, MMU = 0.105658, MPI = 0.13957;
const double METAH = 0.547862, METAP = 0.95778;
const double RHO_M = 0.7754, RHO_G = 0.1462;
const double RHOP_M = 1.465, RHOP_G = 0.400;
const int NBINS = 60;

struct Hist {
  double lo, hi;
  std::vector<double> counts;
  Hist(double l, double h) : lo(l), hi(h), counts(NBINS, 0.0) {}
  void fill(double x, double w = 1.0) {
    if (x < lo || x >= hi) return;
    int b = static_cast<int>((x - lo) / (hi - lo) * NBINS);
    if (b >= 0 && b < NBINS) counts[b] += w;
  }
  void dump(FILE* f, const char* name) const {
    for (int i = 0; i < NBINS; ++i) {
      double c = lo + (hi - lo) * (i + 0.5) / NBINS;
      fprintf(f, "%s,%.6f,%.6f\n", name, c, counts[i]);
    }
  }
};
}  // namespace

int main(int argc, char** argv) {
  std::mt19937_64 rng(13571113);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  const int N = 40000;
  const char* outPath = argc > 1 ? argv[1] : "pluto_mc.csv";
  FILE* f = fopen(outPath, "w");
  if (!f) {
    fprintf(stderr, "cannot open %s\n", outPath);
    return 1;
  }
  fprintf(f, "series,x,count\n");

  // --- eta -> 2mu2e (eta Pade TFF) ---
  {
    TLorentzVector parent(0, 0, 0, METAH);
    const double bound = gen::pluto::etaPadeDoubleBound(METAH, MMU, ME);
    Hist hMuMu(4 * MMU * MMU, std::pow(METAH - 2 * ME, 2)), hEE(4 * ME * ME, std::pow(METAH - 2 * MMU, 2));
    for (int i = 0; i < N; ++i) {
      auto out = gen::pluto::mixedPointlikeEtaPade(parent, MMU, ME, bound, flat);
      hMuMu.fill((out[0] + out[1]).M2());
      hEE.fill((out[2] + out[3]).M2());
    }
    hMuMu.dump(f, "eta2mu2e_mumu");
    hEE.dump(f, "eta2mu2e_ee");
    fprintf(stderr, "eta 2mu2e done\n"); fflush(stderr);
  }

  // --- eta -> 4mu (eta Pade TFF, amplitude-level) ---
  {
    TLorentzVector parent(0, 0, 0, METAH);
    Hist hPair(4 * MMU * MMU, std::pow(METAH - 2 * MMU, 2));
    for (int i = 0; i < N; ++i) {
      auto out = gen::pluto::identicalPointlike(parent, MMU, flat, RHO_M, RHO_G, true);
      hPair.fill((out[0] + out[1]).M2());
      hPair.fill((out[2] + out[3]).M2());  // both pairs are physically equivalent observables
    }
    hPair.dump(f, "eta4mu_pair");
    fprintf(stderr, "eta 4mu done\n"); fflush(stderr);
  }

  // --- eta -> 4e (eta Pade TFF, amplitude-level) ---
  {
    TLorentzVector parent(0, 0, 0, METAH);
    Hist hPair(4 * ME * ME, std::pow(METAH - 2 * ME, 2));
    for (int i = 0; i < N; ++i) {
      auto out = gen::pluto::identicalPointlike(parent, ME, flat, RHO_M, RHO_G, true);
      hPair.fill((out[0] + out[1]).M2());
      hPair.fill((out[2] + out[3]).M2());
    }
    hPair.dump(f, "eta4e_pair");
    fprintf(stderr, "eta 4e done\n"); fflush(stderr);
  }

  // --- eta -> l+l-gamma (2mugamma), eta Pade TFF, dilepton mass ---
  {
    TLorentzVector parent(0, 0, 0, METAH);
    for (const auto& lep : {std::pair<const char*, double>{"eta_2mugamma_mumu", MMU},
                             std::pair<const char*, double>{"eta_2mugamma_ee", ME}}) {
      const double bound = gen::pluto::etaPadeSingleBound(METAH, lep.second);
      Hist h(2 * lep.second, METAH);
      for (int i = 0; i < N; ++i) {
        auto out = gen::pluto::singleDalitzEtaPade(parent, lep.second, bound, flat);
        h.fill((out[0] + out[1]).M());
      }
      h.dump(f, lep.first);
      fprintf(stderr, "%s done\n", lep.first); fflush(stderr);
    }
  }

  // --- eta' -> 2mu2e (two-regime Pade + rho/omega/phi VMD TFF) ---
  {
    TLorentzVector parent(0, 0, 0, METAP);
    const double bound = gen::pluto::etaPrimeDoubleBound(METAP, MMU, ME);
    Hist hMuMu(4 * MMU * MMU, std::pow(METAP - 2 * ME, 2)), hEE(4 * ME * ME, std::pow(METAP - 2 * MMU, 2));
    for (int i = 0; i < N; ++i) {
      auto out = gen::pluto::mixedPointlikeEtaPrime(parent, MMU, ME, bound, flat);
      hMuMu.fill((out[0] + out[1]).M2());
      hEE.fill((out[2] + out[3]).M2());
    }
    hMuMu.dump(f, "etaprime2mu2e_mumu");
    hEE.dump(f, "etaprime2mu2e_ee");
    fprintf(stderr, "etaprime 2mu2e done\n"); fflush(stderr);
  }

  // --- eta' -> 4mu (rho-pole resonant, amplitude-level) ---
  {
    TLorentzVector parent(0, 0, 0, METAP);
    Hist hPair(4 * MMU * MMU, std::pow(METAP - 2 * MMU, 2));
    for (int i = 0; i < N; ++i) {
      auto out = gen::pluto::identicalPointlike(parent, MMU, flat, RHO_M, RHO_G, false);
      hPair.fill((out[0] + out[1]).M2());
      hPair.fill((out[2] + out[3]).M2());
    }
    hPair.dump(f, "etaprime4mu_pair");
    fprintf(stderr, "etaprime 4mu done\n"); fflush(stderr);
  }

  // --- eta' -> 4e (rho-pole resonant, amplitude-level) ---
  {
    TLorentzVector parent(0, 0, 0, METAP);
    Hist hPair(4 * ME * ME, std::pow(METAP - 2 * ME, 2));
    for (int i = 0; i < N; ++i) {
      auto out = gen::pluto::identicalPointlike(parent, ME, flat, RHO_M, RHO_G, false);
      hPair.fill((out[0] + out[1]).M2());
      hPair.fill((out[2] + out[3]).M2());
    }
    hPair.dump(f, "etaprime4e_pair");
    fprintf(stderr, "etaprime 4e done\n"); fflush(stderr);
  }

  // --- eta' -> l+l-gamma (2mugamma/2egamma), Pade + rho/omega/phi VMD TFF, dilepton mass ---
  {
    TLorentzVector parent(0, 0, 0, METAP);
    for (const auto& lep : {std::pair<const char*, double>{"etaprime_2mugamma_mumu", MMU},
                             std::pair<const char*, double>{"etaprime_2mugamma_ee", ME}}) {
      const double bound = gen::pluto::etaPrimeSingleBound(METAP, lep.second);
      Hist h(2 * lep.second, METAP);
      for (int i = 0; i < N; ++i) {
        auto out = gen::pluto::singleDalitzEtaPrime(parent, lep.second, bound, flat);
        h.fill((out[0] + out[1]).M());
      }
      h.dump(f, lep.first);
      fprintf(stderr, "%s done\n", lep.first); fflush(stderr);
    }
  }

  // --- eta/eta' -> l+l- pi+pi- (llpipiChPT), dilepton mass, k^2-weighted to
  // match arXiv:0705.0954 Fig. 7's k^2 * dGamma/d(sqrt(k^2)) convention ---
  {
    struct Ch {
      const char* name;
      double mParent;
      double mLep;
      bool isEtaPrime;
    };
    const std::array<Ch, 4> chans = {{{"eta_pipiee", METAH, ME, false},
                                       {"eta_pipimumu", METAH, MMU, false},
                                       {"etaprime_pipiee", METAP, ME, true},
                                       {"etaprime_pipimumu", METAP, MMU, true}}};
    for (const auto& ch : chans) {
      TLorentzVector parent(0, 0, 0, ch.mParent);
      const double bound =
          gen::pluto::llpipiWeightBound(ch.mParent, MPI, ch.mLep, ch.isEtaPrime, RHO_M, RHO_G, RHOP_M, RHOP_G);
      Hist hLL(2 * ch.mLep, ch.mParent - 2 * MPI);
      for (int i = 0; i < N; ++i) {
        auto out = gen::pluto::llpipiChPT(
            parent, ch.mLep, MPI, ch.isEtaPrime, RHO_M, RHO_G, RHOP_M, RHOP_G, bound, flat);
        const double mll = (out[0] + out[1]).M();
        hLL.fill(mll, mll * mll);
      }
      hLL.dump(f, ch.name);
      fprintf(stderr, "%s done\n", ch.name); fflush(stderr);
    }
  }

  fclose(f);
  printf("done\n");
  return 0;
}
