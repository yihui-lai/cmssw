#ifndef GeneratorInterface_Pythia8Interface_PlutoMixedDoubleDalitz_h
#define GeneratorInterface_Pythia8Interface_PlutoMixedDoubleDalitz_h
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaTFF.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaPrimeTFF.h"
#include "TLorentzVector.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace gen {
  namespace pluto {
    // Petri, arXiv:1010.2378, Eq. (3.36), divided by its angular bound 2.
    // b1,b2 are beta squared; c1,c2 are cos(theta). No identical-fermion
    // exchange diagram exists for the mixed-flavour double-Dalitz channel.
    inline double mixedAngular(double b1, double b2, double c1, double c2, double phi) {
      const double a = b1 * (1 - c1 * c1), b = b2 * (1 - c2 * c2);
      return 1 - (a + b) / 2 + a * b * std::pow(std::sin(phi), 2) / 2;
    }

    // LO pointlike pseudoscalar -> gamma* gamma* -> mu mu ee, F(s1,s2)=1.
    // dGamma/(ds1 ds2 dcos1 dcos2 dphi) is proportional to
    // lambda(1,s1/M^2,s2/M^2)^(3/2) beta1 beta2 angular / (s1 s2).
    // Logarithmic s proposals cancel both photon-pole factors. Every remaining
    // factor is bounded by one, giving a fixed, non-adaptive rejection envelope.
    template <class Flat>
    std::array<TLorentzVector, 4> mixedPointlike(const TLorentzVector& parent, double m1, double m2, Flat flat) {
      const double M = parent.M(), M2 = M * M;
      if (!(m1 > 0 && m2 > 0 && M > 2 * (m1 + m2)))
        throw std::runtime_error("Mixed double-Dalitz masses are below threshold");
      const double low1 = 4 * m1 * m1, low2 = 4 * m2 * m2;
      const double log1 = std::log(std::pow(M - 2 * m2, 2) / low1);
      const double log2 = std::log(std::pow(M - 2 * m1, 2) / low2);
      const double pi = std::acos(-1.);
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double s1 = low1 * std::exp(log1 * flat()), s2 = low2 * std::exp(log2 * flat());
        if (std::sqrt(s1) + std::sqrt(s2) >= M)
          continue;
        const double b1 = 1 - low1 / s1, b2 = 1 - low2 / s2;
        const double x = s1 / M2, y = s2 / M2;
        const double lambda = std::max(0., std::pow(1 - x - y, 2) - 4 * x * y);
        const double c1 = 2 * flat() - 1, c2 = 2 * flat() - 1, phi = 2 * pi * flat();
        const double weight = std::pow(lambda, 1.5) * std::sqrt(b1 * b2) * mixedAngular(b1, b2, c1, c2, phi);
        if (!std::isfinite(weight) || weight < 0 || weight > 1)
          throw std::runtime_error("Invalid mixed double-Dalitz rejection weight");
        if (flat() >= weight)
          continue;
        const double k = M * std::sqrt(lambda) / 2;
        const double e1 = (M2 + s1 - s2) / (2 * M), e2 = (M2 + s2 - s1) / (2 * M);
        const double p1 = std::sqrt(s1 * b1) / 2, p2 = std::sqrt(s2 * b2) / 2;
        const double t1 = std::sqrt(1 - c1 * c1), t2 = std::sqrt(1 - c2 * c2);
        // An independent azimuth about the pair axis supplies the third overall
        // orientation angle in addition to the isotropic pair direction below.
        const double az = 2 * pi * flat();
        std::array<TLorentzVector, 4> out = {{
            TLorentzVector(p1 * t1 * std::cos(az), p1 * t1 * std::sin(az), p1 * c1, std::sqrt(s1) / 2),
            TLorentzVector(-p1 * t1 * std::cos(az), -p1 * t1 * std::sin(az), -p1 * c1, std::sqrt(s1) / 2),
            TLorentzVector(p2 * t2 * std::cos(az + phi), p2 * t2 * std::sin(az + phi), p2 * c2, std::sqrt(s2) / 2),
            TLorentzVector(
                -p2 * t2 * std::cos(az + phi), -p2 * t2 * std::sin(az + phi), -p2 * c2, std::sqrt(s2) / 2)}};
        const double theta = std::acos(2 * flat() - 1), orientPhi = 2 * pi * flat();
        for (unsigned i = 0; i < 4; ++i) {
          out[i].Boost(0, 0, i < 2 ? k / e1 : -k / e2);
          out[i].RotateY(theta);
          out[i].RotateZ(orientPhi);
          out[i].Boost(parent.BoostVector());
        }
        return out;
      }
      throw std::runtime_error("Mixed double-Dalitz rejection limit exhausted");
    }

    // Same density as mixedPointlike, but with the constant F=1 replaced by
    // a factorized double-virtual form factor F(s1,s2) = Frho(s1)*Frho(s2),
    // Frho the single rho0-pole VMD factor also used in PlutoSingleDalitz.h.
    // Escribano & Gonzalez-Solis, arXiv:1511.04916, Eq. (8) ("standard
    // factorisation ansatz"), motivate this product form for the
    // double-virtual TFF from single-virtual data; that paper's Table 4/6
    // shows this correction is large for lepton pairs reaching sizeable
    // virtuality -- e.g. eta -> e+e-mu+mu- increases by ~50% and
    // eta' -> e+e-mu+mu- by roughly a factor 2 relative to the constant-F
    // (their "QED") baseline that mixedPointlike implements. This is only
    // valid for the non-identical-lepton (2mu2e) channel: applying the same
    // factorized reweighting to the identical-lepton case would need to be
    // done separately for the direct and exchange kinematics before they
    // interfere, which is not implemented here -- see PlutoIdenticalDoubleDalitz.h.
    template <class Flat>
    std::array<TLorentzVector, 4> mixedPointlikeResonant(
        const TLorentzVector& parent, double m1, double m2, double rhoMass, double rhoGamma, Flat flat) {
      const double M = parent.M(), M2 = M * M;
      if (!(m1 > 0 && m2 > 0 && M > 2 * (m1 + m2)))
        throw std::runtime_error("Mixed double-Dalitz masses are below threshold");
      const double low1 = 4 * m1 * m1, low2 = 4 * m2 * m2;
      const double log1 = std::log(std::pow(M - 2 * m2, 2) / low1);
      const double log2 = std::log(std::pow(M - 2 * m1, 2) / low2);
      const double pi = std::acos(-1.);
      const double poleBound = rhoMass * rhoMass * rhoGamma * rhoGamma;  // each Frho <= 1, so the
                                                                          // overall bound is unchanged
      const auto formFactor = [&](double s) {
        const double diff = rhoMass * rhoMass - s;
        return poleBound / (diff * diff + poleBound);
      };
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double s1 = low1 * std::exp(log1 * flat()), s2 = low2 * std::exp(log2 * flat());
        if (std::sqrt(s1) + std::sqrt(s2) >= M)
          continue;
        const double b1 = 1 - low1 / s1, b2 = 1 - low2 / s2;
        const double x = s1 / M2, y = s2 / M2;
        const double lambda = std::max(0., std::pow(1 - x - y, 2) - 4 * x * y);
        const double c1 = 2 * flat() - 1, c2 = 2 * flat() - 1, phi = 2 * pi * flat();
        const double weight = std::pow(lambda, 1.5) * std::sqrt(b1 * b2) * mixedAngular(b1, b2, c1, c2, phi) *
                               formFactor(s1) * formFactor(s2);
        if (!std::isfinite(weight) || weight < 0 || weight > 1)
          throw std::runtime_error("Invalid mixed double-Dalitz (resonant) rejection weight");
        if (flat() >= weight)
          continue;
        const double k = M * std::sqrt(lambda) / 2;
        const double e1 = (M2 + s1 - s2) / (2 * M), e2 = (M2 + s2 - s1) / (2 * M);
        const double p1 = std::sqrt(s1 * b1) / 2, p2 = std::sqrt(s2 * b2) / 2;
        const double t1 = std::sqrt(1 - c1 * c1), t2 = std::sqrt(1 - c2 * c2);
        const double az = 2 * pi * flat();
        std::array<TLorentzVector, 4> out = {{
            TLorentzVector(p1 * t1 * std::cos(az), p1 * t1 * std::sin(az), p1 * c1, std::sqrt(s1) / 2),
            TLorentzVector(-p1 * t1 * std::cos(az), -p1 * t1 * std::sin(az), -p1 * c1, std::sqrt(s1) / 2),
            TLorentzVector(p2 * t2 * std::cos(az + phi), p2 * t2 * std::sin(az + phi), p2 * c2, std::sqrt(s2) / 2),
            TLorentzVector(
                -p2 * t2 * std::cos(az + phi), -p2 * t2 * std::sin(az + phi), -p2 * c2, std::sqrt(s2) / 2)}};
        const double theta = std::acos(2 * flat() - 1), orientPhi = 2 * pi * flat();
        for (unsigned i = 0; i < 4; ++i) {
          out[i].Boost(0, 0, i < 2 ? k / e1 : -k / e2);
          out[i].RotateY(theta);
          out[i].RotateZ(orientPhi);
          out[i].Boost(parent.BoostVector());
        }
        return out;
      }
      throw std::runtime_error("Mixed double-Dalitz (resonant) rejection limit exhausted");
    }

    // eta (221) only: same as mixedPointlikeResonant, but each pair's
    // rho-pole VMD factor is replaced by the actual data-fitted eta TFF
    // (PlutoEtaTFF.h, arXiv:1504.07742), applied as a factorized product
    // over the two (generally different) dilepton invariant masses --
    // arXiv:1511.04916 Eq. 8's "standard factorisation ansatz", the same
    // pattern already used for the rho-pole version, just with a better
    // (data-validated, eta-specific) single-virtual input.
    template <class Flat>
    std::array<TLorentzVector, 4> mixedPointlikeEtaPade(
        const TLorentzVector& parent, double m1, double m2, double weightBound, Flat flat) {
      const double M = parent.M(), M2 = M * M;
      if (!(m1 > 0 && m2 > 0 && M > 2 * (m1 + m2)))
        throw std::runtime_error("Mixed double-Dalitz masses are below threshold");
      const double low1 = 4 * m1 * m1, low2 = 4 * m2 * m2;
      const double log1 = std::log(std::pow(M - 2 * m2, 2) / low1);
      const double log2 = std::log(std::pow(M - 2 * m1, 2) / low2);
      const double pi = std::acos(-1.);
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double s1 = low1 * std::exp(log1 * flat()), s2 = low2 * std::exp(log2 * flat());
        if (std::sqrt(s1) + std::sqrt(s2) >= M)
          continue;
        const double b1 = 1 - low1 / s1, b2 = 1 - low2 / s2;
        const double x = s1 / M2, y = s2 / M2;
        const double lambda = std::max(0., std::pow(1 - x - y, 2) - 4 * x * y);
        const double c1 = 2 * flat() - 1, c2 = 2 * flat() - 1, phi = 2 * pi * flat();
        const double formFactor =
            std::pow(etaPadeFormFactorNorm(s1), 2) * std::pow(etaPadeFormFactorNorm(s2), 2) / weightBound;
        const double weight =
            std::pow(lambda, 1.5) * std::sqrt(b1 * b2) * mixedAngular(b1, b2, c1, c2, phi) * formFactor;
        if (!std::isfinite(weight) || weight < 0 || weight > 1)
          throw std::runtime_error(
              "Mixed double-Dalitz (eta Pade) rejection weight violated -- raise the bound's safety margin");
        if (flat() >= weight)
          continue;
        const double k = M * std::sqrt(lambda) / 2;
        const double e1 = (M2 + s1 - s2) / (2 * M), e2 = (M2 + s2 - s1) / (2 * M);
        const double p1 = std::sqrt(s1 * b1) / 2, p2 = std::sqrt(s2 * b2) / 2;
        const double t1 = std::sqrt(1 - c1 * c1), t2 = std::sqrt(1 - c2 * c2);
        const double az = 2 * pi * flat();
        std::array<TLorentzVector, 4> out = {{
            TLorentzVector(p1 * t1 * std::cos(az), p1 * t1 * std::sin(az), p1 * c1, std::sqrt(s1) / 2),
            TLorentzVector(-p1 * t1 * std::cos(az), -p1 * t1 * std::sin(az), -p1 * c1, std::sqrt(s1) / 2),
            TLorentzVector(p2 * t2 * std::cos(az + phi), p2 * t2 * std::sin(az + phi), p2 * c2, std::sqrt(s2) / 2),
            TLorentzVector(
                -p2 * t2 * std::cos(az + phi), -p2 * t2 * std::sin(az + phi), -p2 * c2, std::sqrt(s2) / 2)}};
        const double theta = std::acos(2 * flat() - 1), orientPhi = 2 * pi * flat();
        for (unsigned i = 0; i < 4; ++i) {
          out[i].Boost(0, 0, i < 2 ? k / e1 : -k / e2);
          out[i].RotateY(theta);
          out[i].RotateZ(orientPhi);
          out[i].Boost(parent.BoostVector());
        }
        return out;
      }
      throw std::runtime_error("Mixed double-Dalitz (eta Pade) rejection limit exhausted");
    }

    // Same density as mixedPointlike for an arbitrary factorized form factor:
    // formFactor(s1,s2) must return F(s1,s2)-squared already divided by a
    // rejection bound, i.e. in [0,1]. (The variants above predate this and
    // carry their own copies of the same body; left untouched.)
    template <class FormFactor, class Flat>
    std::array<TLorentzVector, 4> mixedPointlikeGeneric(
        const TLorentzVector& parent, double m1, double m2, FormFactor formFactor, Flat flat) {
      const double M = parent.M(), M2 = M * M;
      if (!(m1 > 0 && m2 > 0 && M > 2 * (m1 + m2)))
        throw std::runtime_error("Mixed double-Dalitz masses are below threshold");
      const double low1 = 4 * m1 * m1, low2 = 4 * m2 * m2;
      const double log1 = std::log(std::pow(M - 2 * m2, 2) / low1);
      const double log2 = std::log(std::pow(M - 2 * m1, 2) / low2);
      const double pi = std::acos(-1.);
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double s1 = low1 * std::exp(log1 * flat()), s2 = low2 * std::exp(log2 * flat());
        if (std::sqrt(s1) + std::sqrt(s2) >= M)
          continue;
        const double b1 = 1 - low1 / s1, b2 = 1 - low2 / s2;
        const double x = s1 / M2, y = s2 / M2;
        const double lambda = std::max(0., std::pow(1 - x - y, 2) - 4 * x * y);
        const double c1 = 2 * flat() - 1, c2 = 2 * flat() - 1, phi = 2 * pi * flat();
        const double weight =
            std::pow(lambda, 1.5) * std::sqrt(b1 * b2) * mixedAngular(b1, b2, c1, c2, phi) * formFactor(s1, s2);
        if (!std::isfinite(weight) || weight < 0 || weight > 1)
          throw std::runtime_error(
              "Mixed double-Dalitz rejection weight violated -- raise the form-factor bound's safety margin");
        if (flat() >= weight)
          continue;
        const double k = M * std::sqrt(lambda) / 2;
        const double e1 = (M2 + s1 - s2) / (2 * M), e2 = (M2 + s2 - s1) / (2 * M);
        const double p1 = std::sqrt(s1 * b1) / 2, p2 = std::sqrt(s2 * b2) / 2;
        const double t1 = std::sqrt(1 - c1 * c1), t2 = std::sqrt(1 - c2 * c2);
        const double az = 2 * pi * flat();
        std::array<TLorentzVector, 4> out = {{
            TLorentzVector(p1 * t1 * std::cos(az), p1 * t1 * std::sin(az), p1 * c1, std::sqrt(s1) / 2),
            TLorentzVector(-p1 * t1 * std::cos(az), -p1 * t1 * std::sin(az), -p1 * c1, std::sqrt(s1) / 2),
            TLorentzVector(p2 * t2 * std::cos(az + phi), p2 * t2 * std::sin(az + phi), p2 * c2, std::sqrt(s2) / 2),
            TLorentzVector(
                -p2 * t2 * std::cos(az + phi), -p2 * t2 * std::sin(az + phi), -p2 * c2, std::sqrt(s2) / 2)}};
        const double theta = std::acos(2 * flat() - 1), orientPhi = 2 * pi * flat();
        for (unsigned i = 0; i < 4; ++i) {
          out[i].Boost(0, 0, i < 2 ? k / e1 : -k / e2);
          out[i].RotateY(theta);
          out[i].RotateZ(orientPhi);
          out[i].Boost(parent.BoostVector());
        }
        return out;
      }
      throw std::runtime_error("Mixed double-Dalitz rejection limit exhausted");
    }

    // eta' (331) only: factorized double-virtual form factor
    // |F(s1)|^2 |F(s2)|^2 (arXiv:1511.04916 Eq. 8's factorisation ansatz)
    // built from the two-regime Pade + rho/omega/phi VMD single-virtual
    // form factor of PlutoEtaPrimeTFF.h. weightBound comes from
    // etaPrimeDoubleBound, computed once by the caller.
    template <class Flat>
    std::array<TLorentzVector, 4> mixedPointlikeEtaPrime(
        const TLorentzVector& parent, double m1, double m2, double weightBound, Flat flat) {
      return mixedPointlikeGeneric(
          parent,
          m1,
          m2,
          [weightBound](double s1, double s2) {
            return etaPrimeFormFactorNormSq(s1) * etaPrimeFormFactorNormSq(s2) / weightBound;
          },
          flat);
    }
  }  // namespace pluto
}  // namespace gen
#endif
