#ifndef GeneratorInterface_Pythia8Interface_PlutoBSMDalitz_h
#define GeneratorInterface_Pythia8Interface_PlutoBSMDalitz_h
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaTFF.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaPrimeTFF.h"
#include "TLorentzVector.h"
#include <array>
#include <cmath>
#include <stdexcept>

// BSM generalization of PlutoSingleDalitz.h's P -> gamma l+ l- Dalitz decay
// to P -> A l+ l-, where A is a new, on-shell, non-zero-mass neutral boson
// (a dark photon or an ALP) instead of a real photon. A is inserted into
// the Pythia event record as an explicit daughter (PDG id and mass/width/
// lifetime set by the caller's fragment via Pythia's own particle-data
// settings, e.g. "<id>:all = ..."); its own decay to a *different* lepton
// flavor is handled afterward by Pythia's native decay machinery, not by
// this code -- PlutoDecayer only needs to produce correct P -> A l+ l-
// kinematics and hand A off with the right PDG id.
//
// Physics content and its limits (read before trusting this past leading
// order in m_A/M):
//   - Kinematics (momenta, energies, the P -> A + (l+l- system) two-body
//     split via the Kallen/triangle function, generalizing the massless-
//     photon case's k = M(1-x)/2) are EXACT for any m_A < M - 2*m_l.
//   - The l+l- *shape* (dGamma/dx, the (1+cos^2(theta)+...) angular
//     correlation) is carried over unchanged from the massless-photon
//     Dalitz formula already used elsewhere in this package. For a dark
//     photon this is physically motivated at leading order (A couples to
//     leptons exactly like a photon via kinetic mixing, so the l+l-
//     producing vertex and its momentum dependence do not care whether
//     the recoiling boson is real, virtual, or a massive on-shell dark
//     photon) but is NOT a re-derivation of the P-gamma-gamma* vertex
//     with a massive external leg: it is a good approximation only while
//     m_A is not a large fraction of M (the m_A windows requested for
//     eta/eta' here reach up to ~0.6-0.8 of M, where this approximation
//     is not guaranteed to be accurate -- treat normalization/shape
//     uncertainties accordingly until a genuine massive-vector-exchange
//     amplitude is derived or sourced from the literature).
//   - For the ALP hypothesis (pseudoscalar, J^PC = 0^-+), the above
//     vector-exchange shape does not apply at all (the photon-like
//     coupling structure that justifies it is specific to a particle
//     with the photon's quantum numbers). Lacking a sourced ALP
//     production amplitude, alpWeight() below uses genuine, unweighted
//     three-body phase space (the "flat" I floor for the production
//     step, with no angular correlation) -- a deliberate, documented
//     approximation, not a derived ALP amplitude.
namespace gen {
  namespace pluto {
    namespace bsmdalitz {
      // Two-body momentum (in the parent rest frame) of a P -> A + X split,
      // where X is the l+l- system with invariant mass q = sqrt(q2). Uses
      // the factorized (q+mA)(q-mA) form rather than the quartic Kallen
      // polynomial for numerical stability near threshold.
      inline double splitMomentum(double M, double q, double mA) {
        const double a = M * M - (q + mA) * (q + mA);
        const double b = M * M - (q - mA) * (q - mA);
        if (!(a > 0 && b > 0))
          return -1.;  // outside the physical (q, mA) region for this M
        return std::sqrt(a * b) / (2 * M);
      }
    }  // namespace bsmdalitz

    // Core kinematics, shared by the dark-photon and ALP weight functions
    // below: samples x = q^2/M^2 (q = l+l- invariant mass) log-uniformly
    // between threshold and the P -> A + (l+l-) edge, builds the l+l- pair
    // in its own rest frame (angle c = cos(theta) relative to the boost
    // axis), boosts by the exact two-body velocity -p/E_ll (generalizing
    // PlutoSingleDalitz.h's massless-photon -k/(M-k)), and places A with
    // the correct (p, E_A). weight(x, c) must return the physics weight
    // already divided by a safe upper bound, i.e. in [0,1].
    template <class WeightXC, class Flat>
    std::array<TLorentzVector, 3> bsmSingleDalitzGeneric(
        const TLorentzVector& parent, double mLep, double mA, WeightXC weight, Flat flat) {
      using namespace bsmdalitz;
      const double M = parent.M(), M2 = M * M;
      if (!(mLep > 0 && mA > 0 && M > 2 * mLep + mA))
        throw std::runtime_error("BSM single Dalitz masses are below threshold");
      const double xlo = 4 * mLep * mLep / M2, xhi = (M - mA) * (M - mA) / M2;
      if (!(xhi > xlo))
        throw std::runtime_error("BSM single Dalitz: empty kinematic range (m_A too close to threshold)");
      const double logRange = std::log(xhi / xlo), pi = std::acos(-1.);
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double x = xlo * std::exp(logRange * flat());
        const double beta2 = 1 - 4 * mLep * mLep / (x * M2);
        if (!(beta2 > 0))
          continue;
        const double beta = std::sqrt(beta2);
        const double c = 2 * flat() - 1;
        const double w = weight(x, c);
        if (!std::isfinite(w) || w < 0 || w > 1)
          throw std::runtime_error("Invalid BSM single Dalitz rejection weight -- raise the bound's safety margin");
        if (flat() >= w)
          continue;
        const double q = M * std::sqrt(x);
        const double p = splitMomentum(M, q, mA);
        if (!(p > 0))
          continue;  // numerically at the edge; resample
        const double eA = std::sqrt(p * p + mA * mA);
        const double eLL = M - eA;
        const double betaZ = -p / eLL;
        const double pl = q * beta / 2;
        const double s = std::sqrt(1 - c * c), phi = 2 * pi * flat();
        TLorentzVector a(0, 0, p, eA);
        TLorentzVector lMinus(pl * s * std::cos(phi), pl * s * std::sin(phi), pl * c, q / 2);
        TLorentzVector lPlus(-pl * s * std::cos(phi), -pl * s * std::sin(phi), -pl * c, q / 2);
        lMinus.Boost(0, 0, betaZ);
        lPlus.Boost(0, 0, betaZ);
        const double theta = std::acos(2 * flat() - 1), orientPhi = 2 * pi * flat();
        std::array<TLorentzVector, 3> out = {{lMinus, lPlus, a}};
        for (auto& part : out) {
          part.RotateY(theta);
          part.RotateZ(orientPhi);
          part.Boost(parent.BoostVector());
        }
        return out;
      }
      throw std::runtime_error("BSM single Dalitz rejection limit exhausted");
    }

    // Dark-photon-hypothesis weight: the same Kroll-Wada-type beta(x) *
    // (1+c^2+(1-beta^2)(1-c^2))/2 angular shape and transition-form-factor
    // dressing already used for the real-photon Dalitz decay (PlutoSingleDalitz.h),
    // evaluated at q2 = x*M^2 -- see the file header for the leading-order
    // caveat. isEtaPrime selects the eta (Pade) vs eta' (Pade+VMD) TFF.
    //
    // Grid-search bound + ready-to-use weight functor for the dark-photon
    // hypothesis. formFactorBound must come from darkPhotonSingleBound()
    // below (computed once per job, same pattern as etaPadeSingleBound/
    // etaPrimeSingleBound).
    inline double darkPhotonFormFactorNormSq(double q2, bool isEtaPrime) {
      return isEtaPrime ? etaPrimeFormFactorNormSq(q2) : std::pow(etaPadeFormFactorNorm(q2), 2);
    }

    inline double darkPhotonSingleBound(double mParent, double mLep, double mA, bool isEtaPrime, unsigned grid = 2000) {
      const double M2 = mParent * mParent;
      const double xlo = 4 * mLep * mLep / M2, xhi = (mParent - mA) * (mParent - mA) / M2;
      if (!(xhi > xlo))
        throw std::runtime_error("darkPhotonSingleBound: empty kinematic range");
      double maxW = 0.;
      for (unsigned i = 0; i <= grid; ++i) {
        const double x = xlo * std::exp(std::log(xhi / xlo) * i / grid);
        maxW = std::max(maxW, darkPhotonFormFactorNormSq(x * M2, isEtaPrime));
      }
      if (!(maxW > 0.))
        throw std::runtime_error("darkPhotonSingleBound: grid search found no valid point");
      return 1.2 * maxW;  // the angular factor's own max is exactly 1 (see header derivation), folded in at call time
    }

    template <class Flat>
    std::array<TLorentzVector, 3> singleDalitzDarkPhoton(
        const TLorentzVector& parent, double mLep, double mA, bool isEtaPrime, double formFactorBound, Flat flat) {
      const double M2 = parent.M() * parent.M();
      auto weight = [&](double x, double c) {
        const double beta2 = 1 - 4 * mLep * mLep / (x * M2);
        const double beta = std::sqrt(std::max(beta2, 0.));
        const double ang = (1 + c * c + (1 - beta2) * (1 - c * c)) / 2;
        // (1-x)^3 carried over unchanged from the massless-photon formula (PlutoSingleDalitz.h) --
        // this is the inherited-shape approximation documented at the top of this file.
        return beta * std::pow(1 - x, 3) * ang * darkPhotonFormFactorNormSq(x * M2, isEtaPrime) / formFactorBound;
      };
      return bsmSingleDalitzGeneric(parent, mLep, mA, weight, flat);
    }

    // ALP-hypothesis weight: genuine, unweighted three-body phase space for
    // the production step (no angular correlation, no transition form
    // factor) -- see the file header for why this is the documented
    // approximation rather than a derived ALP amplitude. beta(x) is the
    // l+l- pair's own two-body phase-space factor (unaffected by what it
    // recoils against); p(x)/pMax is the P -> A + (l+l-) two-body
    // phase-space factor.
    inline double alpSingleBound(double mParent, double mLep, double mA, unsigned grid = 2000) {
      const double M2 = mParent * mParent;
      const double xlo = 4 * mLep * mLep / M2, xhi = (mParent - mA) * (mParent - mA) / M2;
      if (!(xhi > xlo))
        throw std::runtime_error("alpSingleBound: empty kinematic range");
      double maxW = 0.;
      for (unsigned i = 0; i <= grid; ++i) {
        const double x = xlo + (xhi - xlo) * i / grid;  // linear grid: no 1/x pole to cancel here (no TFF dressing)
        const double beta2 = 1 - 4 * mLep * mLep / (x * M2);
        if (!(beta2 > 0))
          continue;
        const double q = mParent * std::sqrt(x);
        const double p = bsmdalitz::splitMomentum(mParent, q, mA);
        if (!(p > 0))
          continue;
        maxW = std::max(maxW, std::sqrt(beta2) * p);
      }
      if (!(maxW > 0.))
        throw std::runtime_error("alpSingleBound: grid search found no valid point");
      return 1.2 * maxW;
    }

    template <class Flat>
    std::array<TLorentzVector, 3> singleDalitzALP(
        const TLorentzVector& parent, double mLep, double mA, double bound, Flat flat) {
      const double M = parent.M(), M2 = M * M;
      auto weight = [&](double x, double /*c*/) {
        const double beta2 = 1 - 4 * mLep * mLep / (x * M2);
        if (!(beta2 > 0))
          return 0.;
        const double q = M * std::sqrt(x);
        const double p = bsmdalitz::splitMomentum(M, q, mA);
        if (!(p > 0))
          return 0.;
        return std::sqrt(beta2) * p / bound;
      };
      return bsmSingleDalitzGeneric(parent, mLep, mA, weight, flat);
    }
  }  // namespace pluto
}  // namespace gen
#endif
