#ifndef GeneratorInterface_Pythia8Interface_EtaBSMPiPiA_h
#define GeneratorInterface_Pythia8Interface_EtaBSMPiPiA_h
#include "GeneratorInterface/Pythia8Interface/interface/EtaLLPiPiChPT.h"
#include "GeneratorInterface/Pythia8Interface/interface/EtaBSMDalitz.h"
#include "TLorentzVector.h"
#include <array>
#include <cmath>
#include <stdexcept>

// BSM generalization of EtaLLPiPiChPT.h's eta/eta' -> l+ l- pi+ pi- to
// eta/eta' -> A pi+ pi-, A a new on-shell boson (dark photon or ALP),
// paralleling EtaBSMDalitz.h's generalization of the single-Dalitz decay.
// A is inserted as an explicit event-record daughter; its own decay (to a
// *different* lepton flavor than this header's name might suggest -- here
// A -> mu+ mu-) is handled by Pythia's native decay machinery, not here.
//
// Where the dark-photon weight below comes from: in llpipiChPT's amplitude,
// the dilepton invariant mass q = s_l only enters through two factors --
// bl = 1 - 4 m_l^2/q (the dilepton pair's own two-body phase space) and
// Fbar(q) (the coherent rho/rho(1450) photon-to-dilepton form factor). For
// A on-shell, q = m_A^2 is FIXED rather than sampled, so BOTH of those
// factors become q-INDEPENDENT CONSTANTS for a given m_A: they rescale the
// overall rate (by ~epsilon^2 x |Fbar(m_A^2)|^2, physically correct for a
// dark photon whose coupling to leptons is a photon-like constant rescale)
// but do not affect the SHAPE of anything sampled here. The pion-pair mass
// dependence P(s)*Omega(s) is untouched -- it describes the eta/eta' -> pi
// pi (virtual photon) vertex, which does not care what the recoiling boson
// is, under the same photon-like-coupling assumption already used for the
// 2muA channel (same leading-order caveat: not re-derived for a genuinely
// massive exchange, see EtaBSMDalitz.h). The lepton-pair-orientation
// variables (cl, phi) of the original 4-body amplitude have no analogue
// here (A has not split into 2 bodies at this level) -- rather than
// analytically marginalize them out of the weight, they are kept as
// "phantom" sampled variables that enter the weight exactly as in the
// original formula but are discarded afterward, which reproduces their
// exact marginal effect without needing to re-derive it by hand.
namespace gen {
  namespace pluto {
    namespace bsmpipi {
      inline double lambda1(double x, double y) { return std::max(0., std::pow(1 - x - y, 2) - 4 * x * y); }
    }  // namespace bsmpipi

    // Dark-photon hypothesis: reuses llpipiWeight(s, mA2, ...) (P(s)*Omega(s)*Fbar(mA2),
    // Fbar's own q-dependence now frozen at mA2) for the pion-pair mass shape, and the
    // llpipiChPT angular/phase-space structure (with the dilepton-plane variables cl, phi
    // kept as phantom, unused-for-output sampling variables -- see file header).
    inline double bsmPiPiAWeightDarkPhoton(double s, double cp, double cl, double phi, double M2, double mA2,
                                            double mPi, double mLepForA, bool isEtaPrime, double rhoMass,
                                            double rhoGamma, double rhopMass, double rhopGamma) {
      const double bp = 1 - 4 * mPi * mPi / s;
      if (!(bp > 0))
        return 0.;
      const double lam = bsmpipi::lambda1(s / M2, mA2 / M2);
      const double bl = 1 - 4 * mLepForA * mLepForA / mA2;  // constant in s; kept for an exact formula match
      if (!(bl > 0))
        throw std::runtime_error("bsmPiPiAWeightDarkPhoton: m_A below 2*mLepForA threshold");
      const double massWeight = llpipiWeight(s, mA2, mPi, mLepForA, isEtaPrime, rhoMass, rhoGamma, rhopMass, rhopGamma);
      return std::pow(lam, 1.5) * (s / M2) * std::pow(bp, 1.5) * std::sqrt(bl) * (1 - cp * cp) *
             (1 - bl * (1 - cl * cl) * std::pow(std::sin(phi), 2)) * massWeight;
    }

    inline double darkPhotonPiPiABound(double mParent, double mPi, double mA, double mLepForA, bool isEtaPrime,
                                        double rhoMass, double rhoGamma, double rhopMass, double rhopGamma,
                                        unsigned grid = 4000) {
      const double M2 = mParent * mParent, mA2 = mA * mA;
      const double smin = 4 * mPi * mPi, smax = std::pow(mParent - mA, 2);
      if (!(smax > smin))
        throw std::runtime_error("darkPhotonPiPiABound: empty kinematic range");
      double maxW = 0.;
      for (unsigned i = 0; i <= grid; ++i) {
        const double s = smin + (smax - smin) * i / grid;
        // cp, cl, phi factors are each individually bounded by <= 1 (same folding-in
        // convention used throughout this package); grid search only the s-dependent part.
        maxW = std::max(maxW, bsmPiPiAWeightDarkPhoton(s, 0., 0., 0., M2, mA2, mPi, mLepForA, isEtaPrime, rhoMass,
                                                         rhoGamma, rhopMass, rhopGamma));
      }
      if (!(maxW > 0.))
        throw std::runtime_error("darkPhotonPiPiABound: grid search found no valid point");
      return 1.2 * maxW;
    }

    template <class Flat>
    std::array<TLorentzVector, 3> piPiADarkPhoton(const TLorentzVector& parent, double mPi, double mA,
                                                    double mLepForA, bool isEtaPrime, double rhoMass, double rhoGamma,
                                                    double rhopMass, double rhopGamma, double bound, Flat flat) {
      const double M = parent.M(), M2 = M * M, mA2 = mA * mA;
      if (!(mPi > 0 && mA > 0 && M > 2 * mPi + mA))
        throw std::runtime_error("piPiADarkPhoton masses are below threshold");
      const double smin = 4 * mPi * mPi, smax = std::pow(M - mA, 2);
      const double pi = std::acos(-1.);
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double s = smin + (smax - smin) * flat();
        const double cp = 2 * flat() - 1, cl = 2 * flat() - 1, phi = 2 * pi * flat();  // cl, phi: phantom, see header
        const double w =
            bsmPiPiAWeightDarkPhoton(s, cp, cl, phi, M2, mA2, mPi, mLepForA, isEtaPrime, rhoMass, rhoGamma, rhopMass, rhopGamma) /
            bound;
        if (!std::isfinite(w) || w < 0 || w > 1)
          throw std::runtime_error("piPiADarkPhoton rejection weight out of range -- raise the bound's safety margin");
        if (flat() >= w)
          continue;
        const double p = bsmdalitz::splitMomentum(M, std::sqrt(s), mA);
        if (!(p > 0))
          continue;
        const double eA = std::sqrt(p * p + mA2), ePP = M - eA;
        const double bp = 1 - 4 * mPi * mPi / s;
        const double pp = std::sqrt(s * bp) / 2;
        const double tp = std::sqrt(1 - cp * cp), az = 2 * pi * flat();
        TLorentzVector a(0, 0, p, eA);
        TLorentzVector piPlus(pp * tp * std::cos(az), pp * tp * std::sin(az), pp * cp, std::sqrt(s) / 2);
        TLorentzVector piMinus(-pp * tp * std::cos(az), -pp * tp * std::sin(az), -pp * cp, std::sqrt(s) / 2);
        piPlus.Boost(0, 0, -p / ePP);
        piMinus.Boost(0, 0, -p / ePP);
        const double theta = std::acos(2 * flat() - 1), orientPhi = 2 * pi * flat();
        std::array<TLorentzVector, 3> out = {{piPlus, piMinus, a}};
        for (auto& part : out) {
          part.RotateY(theta);
          part.RotateZ(orientPhi);
          part.Boost(parent.BoostVector());
        }
        return out;
      }
      throw std::runtime_error("piPiADarkPhoton rejection limit exhausted");
    }

    // ALP hypothesis: genuine, unweighted phase space for the production step
    // (no P(s)*Omega(s) dressing, no (1-cp^2) angular correlation) -- same
    // documented approximation as alpSingleBound/singleDalitzALP.
    inline double alpPiPiABound(double mParent, double mPi, double mA, unsigned grid = 4000) {
      const double smin = 4 * mPi * mPi, smax = std::pow(mParent - mA, 2);
      if (!(smax > smin))
        throw std::runtime_error("alpPiPiABound: empty kinematic range");
      double maxW = 0.;
      for (unsigned i = 0; i <= grid; ++i) {
        const double s = smin + (smax - smin) * i / grid;
        const double bp = 1 - 4 * mPi * mPi / s;
        if (!(bp > 0))
          continue;
        const double p = bsmdalitz::splitMomentum(mParent, std::sqrt(s), mA);
        if (!(p > 0))
          continue;
        maxW = std::max(maxW, std::sqrt(bp) * p);
      }
      if (!(maxW > 0.))
        throw std::runtime_error("alpPiPiABound: grid search found no valid point");
      return 1.2 * maxW;
    }

    template <class Flat>
    std::array<TLorentzVector, 3> piPiAALP(
        const TLorentzVector& parent, double mPi, double mA, double bound, Flat flat) {
      const double M = parent.M();
      const double smin = 4 * mPi * mPi, smax = std::pow(M - mA, 2);
      if (!(mPi > 0 && mA > 0 && M > 2 * mPi + mA))
        throw std::runtime_error("piPiAALP masses are below threshold");
      const double pi = std::acos(-1.);
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double s = smin + (smax - smin) * flat();
        const double bp = 1 - 4 * mPi * mPi / s;
        if (!(bp > 0))
          continue;
        const double p = bsmdalitz::splitMomentum(M, std::sqrt(s), mA);
        if (!(p > 0))
          continue;
        const double w = std::sqrt(bp) * p / bound;
        if (!std::isfinite(w) || w < 0 || w > 1)
          throw std::runtime_error("piPiAALP rejection weight out of range -- raise the bound's safety margin");
        if (flat() >= w)
          continue;
        const double eA = std::sqrt(p * p + mA * mA), ePP = M - eA;
        const double pp = std::sqrt(s * bp) / 2;
        const double cp = 2 * flat() - 1, tp = std::sqrt(1 - cp * cp), az = 2 * pi * flat();
        TLorentzVector a(0, 0, p, eA);
        TLorentzVector piPlus(pp * tp * std::cos(az), pp * tp * std::sin(az), pp * cp, std::sqrt(s) / 2);
        TLorentzVector piMinus(-pp * tp * std::cos(az), -pp * tp * std::sin(az), -pp * cp, std::sqrt(s) / 2);
        piPlus.Boost(0, 0, -p / ePP);
        piMinus.Boost(0, 0, -p / ePP);
        const double theta = std::acos(2 * flat() - 1), orientPhi = 2 * pi * flat();
        std::array<TLorentzVector, 3> out = {{piPlus, piMinus, a}};
        for (auto& part : out) {
          part.RotateY(theta);
          part.RotateZ(orientPhi);
          part.Boost(parent.BoostVector());
        }
        return out;
      }
      throw std::runtime_error("piPiAALP rejection limit exhausted");
    }
  }  // namespace pluto
}  // namespace gen
#endif
