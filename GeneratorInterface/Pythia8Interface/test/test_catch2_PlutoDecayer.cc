#include "catch.hpp"

#include "GeneratorInterface/Pythia8Interface/interface/PlutoIdenticalDoubleDalitz.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoMixedDoubleDalitz.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoPiPiDilepton.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoSingleDalitz.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaPrimeLLPiPi.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoLLPiPiChPT.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaTFF.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaPrimeTFF.h"

#include <random>

namespace {
  bool close(double a, double b) { return std::abs(a - b) < 1e-7 * std::max({1., std::abs(a), std::abs(b)}); }

  // Independent cross-check of the identical-lepton double-Dalitz "direct"
  // score, built straight from the pair invariant masses/angles rather than
  // from the spinor machinery under test.
  double directAnalytic(std::array<TLorentzVector, 4> p, double m) {
    TLorentzVector parent;
    for (const auto& v : p)
      parent += v;
    for (auto& v : p)
      v.Boost(-parent.BoostVector());
    const auto q1 = p[0] + p[1], q2 = p[2] + p[3];
    const double s1 = q1.M2(), s2 = q2.M2(), M2 = parent.M2();
    const auto axis = q1.Vect().Unit();
    p[0].Boost(-q1.BoostVector());
    p[2].Boost(-q2.BoostVector());
    const auto n1 = p[0].Vect().Unit(), n2 = p[2].Vect().Unit();
    const double c1 = n1.Dot(axis), c2 = n2.Dot(axis), triple = n1.Cross(n2).Dot(axis);
    const double b1 = 1 - 4 * m * m / s1, b2 = 1 - 4 * m * m / s2;
    const double lambda = std::pow(M2 - s1 - s2, 2) - 4 * s1 * s2;
    return lambda / (s1 * s2) * (2 - b1 * (1 - c1 * c1) - b2 * (1 - c2 * c2) + b1 * b2 * triple * triple);
  }
}  // namespace

TEST_CASE("PlutoIdenticalDoubleDalitz: 4e/4mu direct score, exchange symmetry, Pauli cancellation, closure",
          "[PlutoDecayer]") {
  std::mt19937_64 rng(9876);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  bool nonzeroInterference = false;
  for (double M : {0.547862, 0.95778})
    for (double m : {0.000511, 0.105658}) {
      const TLorentzVector parent(0, 0, 0, M);
      for (int n = 0; n < 12; ++n) {
        auto p = gen::pluto::mixedPointlike(parent, m, m, flat);
        const auto score = gen::pluto::fourlepton::scores(p, m);
        REQUIRE(close(score.direct, directAnalytic(p, m)));
        REQUIRE(score.coherent >= 0);
        REQUIRE(score.coherent <= 2 * (score.direct + score.exchanged) * (1 + 1e-12));

        if (std::abs(score.coherent - score.direct - score.exchanged) > 1e-5 * (score.direct + score.exchanged))
          nonzeroInterference = true;

        // Swapping the two like-charge leptons swaps direct <-> exchanged,
        // and leaves the physical (coherent) density unchanged.
        std::swap(p[1], p[3]);
        const auto exchanged = gen::pluto::fourlepton::scores(p, m);
        REQUIRE(close(score.direct, exchanged.exchanged));
        REQUIRE(close(score.exchanged, exchanged.direct));
        REQUIRE(close(score.coherent, exchanged.coherent));
        std::swap(p[0], p[2]);
        REQUIRE(close(score.coherent, gen::pluto::fourlepton::scores(p, m).coherent));

        // Boost invariance of the physical density.
        for (auto& v : p)
          v.Boost(0.2, 0.1, 0.3);
        REQUIRE(close(score.coherent, gen::pluto::fourlepton::scores(p, m).coherent));

        // Generated events: correct masses and exact four-momentum closure.
        const TLorentzVector boostedParent(0, 0, 5, std::hypot(5., M));
        const auto generated = gen::pluto::identicalPointlike(boostedParent, m, flat);
        TLorentzVector sum;
        for (const auto& v : generated) {
          sum += v;
          REQUIRE(std::abs(v.M() - m) < 1e-8);
        }
        REQUIRE((sum - boostedParent).Vect().Mag() + std::abs(sum.E() - boostedParent.E()) < 1e-9);
      }
    }
  REQUIRE(nonzeroInterference);

  // Same-momentum, same-spin negative fermions must cancel at amplitude
  // level (Pauli exclusion), independent of the sampler.
  const double m = 0.105658, p = 0.1, t = 0.03;
  const TLorentzVector negative(p, 0, 0, std::hypot(p, m));
  const std::array<TLorentzVector, 4> pauli = {
      {negative, TLorentzVector(-p, t, 0, std::sqrt(p * p + t * t + m * m)), negative,
       TLorentzVector(-p, -t, 0, std::sqrt(p * p + t * t + m * m))}};
  const gen::pluto::fourlepton::Amplitudes amplitudes(pauli, m);
  for (unsigned s = 0; s < 2; ++s)
    for (unsigned a = 0; a < 2; ++a)
      for (unsigned b = 0; b < 2; ++b) {
        const auto value = amplitudes(s, a, s, b);
        REQUIRE(std::norm(value[0] - value[1]) < 1e-20 * std::max(1., std::norm(value[0]) + std::norm(value[1])));
      }
}

TEST_CASE("PlutoIdenticalDoubleDalitz: rho0-pole-dressed 4e/4mu, bound survives, closure, actually changes shape",
          "[PlutoDecayer]") {
  const double rhoMass = 0.7754, rhoGamma = 0.1462;
  const double m = 0.105658;

  // The algebraic bound |D-X|^2 <= 2(|D|^2+|X|^2) must still hold with the
  // form factor applied -- checked directly at a handful of kinematic
  // points, not just relied on as an algebraic argument.
  std::mt19937_64 rng(112233);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  bool everDiffered = false;
  for (double M : {0.547862, 0.95778}) {
    const TLorentzVector parent(0, 0, 0, M);
    for (int n = 0; n < 20; ++n) {
      auto p = gen::pluto::mixedPointlike(parent, m, m, flat);
      const auto plain = gen::pluto::fourlepton::scores(p, m);
      const auto dressed = gen::pluto::fourlepton::scores(p, m, rhoMass, rhoGamma);
      REQUIRE(dressed.coherent >= 0);
      REQUIRE(dressed.coherent <= 2 * (dressed.direct + dressed.exchanged) * (1 + 1e-9));
      if (std::abs(dressed.coherent - plain.coherent) > 1e-6 * std::max(1., plain.coherent))
        everDiffered = true;
    }
  }
  REQUIRE(everDiffered);  // confirms the dressing actually does something, not a no-op

  // rhoMass<=0 must reproduce the undressed behavior exactly (regression
  // guard on the default-argument path used by every existing call site).
  {
    const TLorentzVector parent(0, 0, 0, 0.547862);
    auto p = gen::pluto::mixedPointlike(parent, m, m, flat);
    const auto a = gen::pluto::fourlepton::scores(p, m);
    const auto b = gen::pluto::fourlepton::scores(p, m, 0, 0);
    REQUIRE(a.direct == b.direct);
    REQUIRE(a.exchanged == b.exchanged);
    REQUIRE(a.coherent == b.coherent);
  }

  for (double M : {0.547862, 0.95778}) {
    const TLorentzVector parent(0, 0, 0, M);
    for (int n = 0; n < 20; ++n) {
      auto out = gen::pluto::identicalPointlike(parent, m, flat, rhoMass, rhoGamma);
      TLorentzVector sum;
      for (const auto& v : out) {
        sum += v;
        REQUIRE(std::abs(v.M() - m) < 1e-8);
      }
      REQUIRE((sum - parent).Vect().Mag() + std::abs(sum.E() - parent.E()) < 1e-9);
    }
  }
}

TEST_CASE("PlutoMixedDoubleDalitz: 2mu2e angular integral, envelope, thresholds, closure", "[PlutoDecayer]") {
  // Three-point Gauss-Legendre exactly integrates the angular polynomial's
  // cosine dependence; cross-checks the closed-form average independently.
  const double nodes[] = {-std::sqrt(3. / 5), 0, std::sqrt(3. / 5)};
  const double ws[] = {5. / 18, 4. / 9, 5. / 18};
  for (double b1 : {0., 0.4, 1.})
    for (double b2 : {0., 0.7, 1.}) {
      double avg = 0;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          for (int k = 0; k < 8; ++k) {
            const double w = gen::pluto::mixedAngular(b1, b2, nodes[i], nodes[j], k * std::acos(-1.) / 4);
            REQUIRE(w >= 0);
            REQUIRE(w <= 1);
            avg += ws[i] * ws[j] * w / 8;
          }
      REQUIRE(std::abs(avg - (1 - b1 / 3) * (1 - b2 / 3)) < 1e-14);
    }

  std::mt19937_64 rng(12345);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  for (double mass : {0.547862, 0.95778})
    for (double pz : {0., 5.}) {
      TLorentzVector parent(0, 0, pz, std::hypot(pz, mass));
      for (int n = 0; n < 20; ++n) {
        auto out = gen::pluto::mixedPointlike(parent, 0.105658, 0.000511, flat);
        TLorentzVector sum;
        for (unsigned i = 0; i < 4; ++i) {
          sum += out[i];
          REQUIRE(std::abs(out[i].M() - (i < 2 ? 0.105658 : 0.000511)) < 1e-8);
        }
        REQUIRE((sum - parent).Vect().Mag() + std::abs(sum.E() - parent.E()) < 1e-9);
        REQUIRE((out[0] + out[1]).M() >= 2 * 0.105658);
        REQUIRE((out[2] + out[3]).M() >= 2 * 0.000511);
      }
    }
}

TEST_CASE("PlutoPiPiDilepton: 2e2pi/2mu2pi angular integral and sign, envelope, closure", "[PlutoDecayer]") {
  const double pi = std::acos(-1.);
  const double nodes[] = {-std::sqrt(3. / 5), 0, std::sqrt(3. / 5)}, ws[] = {5. / 18, 4. / 9, 5. / 18};
  for (double beta2 : {0., 0.5, 1.}) {
    double avg = 0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 8; ++k) {
          const double w = gen::pluto::pipiMagneticAngular(beta2, nodes[i], nodes[j], k * pi / 4);
          REQUIRE(w >= 0);
          REQUIRE(w <= 1);
          avg += ws[i] * ws[j] * w / 8;
        }
    REQUIRE(std::abs(avg - (2. / 3) * (1 - beta2 / 3)) < 1e-14);
  }
  // A massless lepton along the magnetic-current direction has zero weight.
  REQUIRE(std::abs(gen::pluto::pipiMagneticAngular(1, 0, 0, pi / 2)) < 1e-14);
  REQUIRE(gen::pluto::pipiMagneticAngular(1, 0, 0, 0) == 1);

  std::mt19937_64 rng(12345);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  for (double mass : {0.547862, 0.95778})
    for (double ml : {0.000511, 0.105658})
      for (double pz : {0., 5.}) {
        TLorentzVector parent(0, 0, pz, std::hypot(pz, mass));
        for (int n = 0; n < 10; ++n) {
          auto out = gen::pluto::pipiMagneticPointlike(parent, ml, 0.139570, flat);
          TLorentzVector sum;
          for (unsigned i = 0; i < 4; ++i) {
            sum += out[i];
            REQUIRE(std::abs(out[i].M() - (i < 2 ? ml : 0.139570)) < 1e-8);
          }
          REQUIRE((sum - parent).Vect().Mag() + std::abs(sum.E() - parent.E()) < 1e-9);
          REQUIRE((out[0] + out[1]).M() >= 2 * ml);
          REQUIRE((out[2] + out[3]).M() >= 2 * 0.139570);
        }
      }
}

TEST_CASE("PlutoMixedDoubleDalitz: 2mu2e resonant form factor, envelope, closure", "[PlutoDecayer]") {
  const double rhoMass = 0.7754, rhoGamma = 0.1462;
  std::mt19937_64 rng(3141);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  for (double mass : {0.547862, 0.95778})
    for (double pz : {0., 5.}) {
      TLorentzVector parent(0, 0, pz, std::hypot(pz, mass));
      for (int n = 0; n < 20; ++n) {
        auto out = gen::pluto::mixedPointlikeResonant(parent, 0.105658, 0.000511, rhoMass, rhoGamma, flat);
        TLorentzVector sum;
        for (unsigned i = 0; i < 4; ++i) {
          sum += out[i];
          REQUIRE(std::abs(out[i].M() - (i < 2 ? 0.105658 : 0.000511)) < 1e-8);
        }
        REQUIRE((sum - parent).Vect().Mag() + std::abs(sum.E() - parent.E()) < 1e-9);
      }
    }
}

TEST_CASE("PlutoSingleDalitz: 2mugamma VMD pole bound, massless photon, closure", "[PlutoDecayer]") {
  const double rhoMass = 0.7754, rhoGamma = 0.1462;  // PDG rho0
  // The form factor used inside the sampler, FF(q^2)/FF(mRho^2), must be
  // bounded by 1 and peak (=1) exactly at the pole -- checked algebraically,
  // independent of the rejection sampler.
  for (double q2 : {0., 0.1, rhoMass * rhoMass, 0.5, 0.9}) {
    const double poleBound = rhoMass * rhoMass * rhoGamma * rhoGamma;
    const double diff = rhoMass * rhoMass - q2;
    const double ff = poleBound / (diff * diff + poleBound);
    REQUIRE(ff > 0);
    REQUIRE(ff <= 1 + 1e-12);
  }
  {
    const double poleBound = rhoMass * rhoMass * rhoGamma * rhoGamma;
    REQUIRE(std::abs(poleBound / (0 + poleBound) - 1) < 1e-12);
  }

  std::mt19937_64 rng(2468);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  for (double mass : {0.547862, 0.95778})
    for (double ml : {0.000511, 0.105658})
      for (double pz : {0., 5.}) {
        TLorentzVector parent(0, 0, pz, std::hypot(pz, mass));
        for (int n = 0; n < 20; ++n) {
          auto out = gen::pluto::singleDalitz(parent, ml, rhoMass, rhoGamma, flat);
          TLorentzVector sum;
          for (unsigned i = 0; i < 3; ++i)
            sum += out[i];
          REQUIRE(std::abs(out[0].M() - ml) < 1e-8);
          REQUIRE(std::abs(out[1].M() - ml) < 1e-8);
          REQUIRE(std::abs(out[2].M()) < 1e-6);  // photon; E^2-p^2 cancellation limits precision
          REQUIRE((sum - parent).Vect().Mag() + std::abs(sum.E() - parent.E()) < 1e-9);
        }
      }
}

TEST_CASE("PlutoEtaPrimeLLPiPi: F0 finiteness/positivity, envelope, closure (eta' only)", "[PlutoDecayer]") {
  const double rhoMass = 0.7754, rhoGamma = 0.1462;
  const double mEtaPrime = 0.95778, ml = 0.105658, mpi = 0.139570;

  // F0 must stay finite and non-negative across the physical domain,
  // including right at the rho pole in the dipion channel.
  const double sPiPiPole = rhoMass * rhoMass;
  for (double sLL : {4 * ml * ml + 1e-4, 0.1, 0.3})
    for (double sPiPi : {4 * mpi * mpi + 1e-4, sPiPiPole, 0.5}) {
      const double f0 = gen::pluto::etaPrimeF0(sLL, sPiPi, ml, mpi, rhoMass, rhoGamma);
      REQUIRE(std::isfinite(f0));
      REQUIRE(f0 >= 0);
    }

  const double bound = gen::pluto::etaPrimeF0Bound(mEtaPrime, ml, mpi, rhoMass, rhoGamma, 100);
  REQUIRE(bound > 0);

  std::mt19937_64 rng(13579);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  const TLorentzVector parent(0, 0, 0, mEtaPrime);
  for (int n = 0; n < 30; ++n) {
    auto out = gen::pluto::etaPrimeLLPiPi(parent, ml, mpi, rhoMass, rhoGamma, bound, flat);
    TLorentzVector sum;
    for (unsigned i = 0; i < 4; ++i) {
      sum += out[i];
      REQUIRE(std::abs(out[i].M() - (i < 2 ? ml : mpi)) < 1e-8);
    }
    REQUIRE((sum - parent).Vect().Mag() + std::abs(sum.E() - parent.E()) < 1e-9);
  }
}

TEST_CASE("PlutoLLPiPiChPT: arXiv:2210.14925-based weight, both parents, envelope, closure",
          "[PlutoDecayer]") {
  const double rhoMass = 0.7754, rhoGamma = 0.1462;
  const double rhopMass = 1.465, rhopGamma = 0.400;
  const double mEta = 0.547862, mEtaPrime = 0.95778, mpi = 0.139570;

  // The polynomial/Omega/Fbar pieces must stay finite and non-negative
  // (|complex|^2) across the physical domain for both parents.
  for (bool isEtaPrime : {false, true}) {
    const double M = isEtaPrime ? mEtaPrime : mEta;
    for (double ml : {0.000511, 0.105658}) {
      const double lowS = 4 * mpi * mpi, highS = std::pow(M - 2 * ml, 2);
      const double lowSl = 4 * ml * ml, highSl = std::pow(M - 2 * mpi, 2);
      for (double s : {lowS + 1e-4, (lowS + highS) / 2, highS - 1e-4}) {
        for (double sl : {lowSl + 1e-4, (lowSl + highSl) / 2, highSl - 1e-4}) {
          const double w = gen::pluto::llpipiWeight(s, sl, mpi, ml, isEtaPrime, rhoMass, rhoGamma, rhopMass, rhopGamma);
          REQUIRE(std::isfinite(w));
          REQUIRE(w >= 0);
        }
      }
    }
  }

  std::mt19937_64 rng(24680);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  for (bool isEtaPrime : {false, true}) {
    const double M = isEtaPrime ? mEtaPrime : mEta;
    const TLorentzVector parent(0, 0, 0, M);
    for (double ml : {0.000511, 0.105658}) {
      const double bound = gen::pluto::llpipiWeightBound(M, mpi, ml, isEtaPrime, rhoMass, rhoGamma, rhopMass, rhopGamma, 80);
      REQUIRE(bound > 0);
      for (int n = 0; n < 20; ++n) {
        auto out = gen::pluto::llpipiChPT(parent, ml, mpi, isEtaPrime, rhoMass, rhoGamma, rhopMass, rhopGamma, bound, flat);
        TLorentzVector sum;
        for (unsigned i = 0; i < 4; ++i) {
          sum += out[i];
          REQUIRE(std::abs(out[i].M() - (i < 2 ? ml : mpi)) < 1e-8);
        }
        REQUIRE((sum - parent).Vect().Mag() + std::abs(sum.E() - parent.E()) < 1e-9);
      }
    }
  }
}

TEST_CASE("PlutoEtaTFF: reproduces the paper's own b_eta self-check, finite and monotonic", "[PlutoDecayer]") {
  // Eq. A.2 of arXiv:1504.07742: b_eta = (t1*r1 - t2)*m_eta^2/t1, quoted
  // there as 0.5749 from these exact Table 5 coefficients -- note the
  // normalization scale in the paper's own Eq. 1 is m_eta^2 (the parent
  // meson's own mass), not m_pi^2; an earlier version of this test used
  // m_pi^2 by mistake and failed by exactly that ratio, which is what
  // led to finding and fixing this test (the underlying
  // etaPadeFormFactor code was already correct once /Q2 was fixed). h is
  // chosen away from catastrophic-cancellation territory (h too small
  // makes (Norm(h)-1) lose precision against the h/m_eta^2 denominator).
  const double mEta = 0.547862;
  const double h = 1e-4;
  const double slope = (gen::pluto::etaPadeFormFactorNorm(h) - 1.0) / (h / (mEta * mEta));
  REQUIRE(std::abs(slope - 0.5749) < 0.01);

  REQUIRE(gen::pluto::etaPadeFormFactorNorm(0.) == 1.0);
  double prev = 0.;
  for (double q2 = 1e-4; q2 < mEta * mEta; q2 += 0.001) {
    const double v = gen::pluto::etaPadeFormFactorNorm(q2);
    REQUIRE(std::isfinite(v));
    REQUIRE(v >= prev);  // confirmed-by-eye monotonic rise across eta's whole range
    prev = v;
  }
}

TEST_CASE("PlutoSingleDalitz/PlutoMixedDoubleDalitz eta Pade branch: envelope, closure, masses",
          "[PlutoDecayer]") {
  const double mEta = 0.547862, mMu = 0.105658, mE = 0.000511;
  std::mt19937_64 rng(987123);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };

  const double bound1 = gen::pluto::etaPadeSingleBound(mEta, mMu);
  REQUIRE(bound1 > 0);
  const TLorentzVector parent(0, 0, 0, mEta);
  for (int n = 0; n < 20; ++n) {
    auto out = gen::pluto::singleDalitzEtaPade(parent, mMu, bound1, flat);
    TLorentzVector sum;
    for (unsigned i = 0; i < 3; ++i)
      sum += out[i];
    REQUIRE(std::abs(out[0].M() - mMu) < 1e-8);
    REQUIRE(std::abs(out[1].M() - mMu) < 1e-8);
    REQUIRE(std::abs(out[2].M()) < 1e-6);
    REQUIRE((sum - parent).Vect().Mag() + std::abs(sum.E() - parent.E()) < 1e-9);
  }

  const double bound2 = gen::pluto::etaPadeDoubleBound(mEta, mMu, mE);
  REQUIRE(bound2 > 0);
  for (int n = 0; n < 20; ++n) {
    auto out = gen::pluto::mixedPointlikeEtaPade(parent, mMu, mE, bound2, flat);
    TLorentzVector sum;
    for (unsigned i = 0; i < 4; ++i) {
      sum += out[i];
      REQUIRE(std::abs(out[i].M() - (i < 2 ? mMu : mE)) < 1e-8);
    }
    REQUIRE((sum - parent).Vect().Mag() + std::abs(sum.E() - parent.E()) < 1e-9);
  }
}

TEST_CASE("PlutoIdenticalDoubleDalitz useEtaPade branch: bound survives, closure, differs from constant-F",
          "[PlutoDecayer]") {
  const double mEta = 0.547862, mMu = 0.105658;
  std::mt19937_64 rng(456789);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  const TLorentzVector parent(0, 0, 0, mEta);

  bool everDiffered = false;
  for (int n = 0; n < 20; ++n) {
    auto p = gen::pluto::mixedPointlike(parent, mMu, mMu, flat);
    const auto plain = gen::pluto::fourlepton::scores(p, mMu);
    const auto padeScore = gen::pluto::fourlepton::scores(p, mMu, 0, 0, true);
    REQUIRE(padeScore.coherent >= 0);
    REQUIRE(padeScore.coherent <= 2 * (padeScore.direct + padeScore.exchanged) * (1 + 1e-9));
    if (std::abs(padeScore.coherent - plain.coherent) > 1e-6 * std::max(1., plain.coherent))
      everDiffered = true;
  }
  REQUIRE(everDiffered);

  for (int n = 0; n < 20; ++n) {
    auto out = gen::pluto::identicalPointlike(parent, mMu, flat, 0, 0, true);
    TLorentzVector sum;
    for (const auto& v : out) {
      sum += v;
      REQUIRE(std::abs(v.M() - mMu) < 1e-8);
    }
    REQUIRE((sum - parent).Vect().Mag() + std::abs(sum.E() - parent.E()) < 1e-9);
  }
}

TEST_CASE("PlutoEtaPrimeTFF: normalization, slope, seam continuity, omega peak, bounds", "[PlutoDecayer]") {
  const double mEtaP = 0.95778;
  auto F = [](double q2) { return std::sqrt(gen::pluto::etaPrimeFormFactorNormSq(q2)); };

  REQUIRE(F(0.) == 1.0);
  // Slope at the origin (GeV^-2); arXiv:1307.2061 Table IV Pade gives ~1.42.
  const double h = 1e-5;
  REQUIRE(std::abs((F(h) - 1.0) / h - 1.42) < 0.05);

  // Padé -> VMD seam is magnitude-matched.
  const double seam = gen::pluto::etaprime_tff::matchQ2;
  REQUIRE(std::abs(gen::pluto::etaPrimeFormFactorNormSq(seam + 1e-9) /
                       gen::pluto::etaPrimeFormFactorNormSq(seam) - 1.) < 1e-6);

  // Omega peak: the VMD region must have its maximum at the omega mass and
  // dominate the continuum below the seam.
  double best = 0, bestMass = 0;
  for (double m = 0.3; m < 0.95; m += 0.0005) {
    const double v = gen::pluto::etaPrimeFormFactorNormSq(m * m);
    REQUIRE(std::isfinite(v));
    if (v > best) {
      best = v;
      bestMass = m;
    }
  }
  REQUIRE(std::abs(bestMass - 0.7825) < 0.005);
  REQUIRE(best > 5 * gen::pluto::etaPrimeFormFactorNormSq(0.6 * 0.6));

  // Samplers: closure, masses, and the rejection bound is never violated
  // (an exceeded bound throws) over many attempts, including the omega peak.
  const double mMu = 0.105658, mE = 0.000511;
  std::mt19937_64 rng(2468101);
  auto flat = [&]() { return std::generate_canonical<double, 53>(rng); };
  const TLorentzVector parent(0, 0, 0, mEtaP);

  const double bound1 = gen::pluto::etaPrimeSingleBound(mEtaP, mMu);
  const double bound1e = gen::pluto::etaPrimeSingleBound(mEtaP, mE);
  const double bound2 = gen::pluto::etaPrimeDoubleBound(mEtaP, mMu, mE);
  REQUIRE(bound1 >= best);
  REQUIRE(bound2 >= best);
  int inOmegaWindow = 0;
  const int n = 4000;
  for (int i = 0; i < n; ++i) {
    auto out = gen::pluto::singleDalitzEtaPrime(parent, mMu, bound1, flat);
    TLorentzVector sum = out[0] + out[1] + out[2];
    REQUIRE((sum - parent).Vect().Mag() + std::abs(sum.E() - parent.E()) < 1e-9);
    REQUIRE(std::abs(out[0].M() - mMu) < 1e-8);
    const double mll = (out[0] + out[1]).M();
    if (mll > 0.77 && mll < 0.79)
      ++inOmegaWindow;
    REQUIRE_NOTHROW(gen::pluto::singleDalitzEtaPrime(parent, mE, bound1e, flat));
    auto four = gen::pluto::mixedPointlikeEtaPrime(parent, mMu, mE, bound2, flat);
    TLorentzVector sum4;
    for (unsigned k = 0; k < 4; ++k) {
      sum4 += four[k];
      REQUIRE(std::abs(four[k].M() - (k < 2 ? mMu : mE)) < 1e-8);
    }
    REQUIRE((sum4 - parent).Vect().Mag() + std::abs(sum4.E() - parent.E()) < 1e-9);
  }
  // The omega/rho region carries a few percent of the eta' -> mu mu gamma rate.
  REQUIRE(inOmegaWindow > 0.02 * n);
}
