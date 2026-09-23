#ifndef GeneratorInterface_Pythia8Interface_PlutoSingleDalitz_h
#define GeneratorInterface_Pythia8Interface_PlutoSingleDalitz_h
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaTFF.h"
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaPrimeTFF.h"
#include "TLorentzVector.h"
#include <array>
#include <cmath>
#include <stdexcept>

namespace gen {
  namespace pluto {
    // Kroll-Wada single Dalitz decay P -> gamma l+ l-, with a rho0-pole
    // vector-meson-dominance form factor matching EvtGen's
    // EvtEta2MuMuGamma model: |F(q^2)|^2 ~ 1/((mRho^2-q^2)^2 + mRho^2 GammaRho^2),
    // rather than a constant. With x=q^2/M^2 the dilepton mass fraction and
    // theta the l- angle in the dilepton rest frame (relative to the photon
    // direction; the density only depends on cos^2(theta), so the reference
    // sign is immaterial), up to constants independent of the decay
    // variables:
    //   dGamma/(dx dcos(theta)) ~ (1/x) beta(x) (1-x)^3 |F(x M^2)|^2
    //     * (1 + cos(theta)^2 + (1-beta(x)^2) sin(theta)^2),
    // beta(x) = sqrt(1 - 4 m^2/(x M^2)). A logarithmic x proposal cancels
    // the photon pole; the angular/phase-space factors are bounded by 1 (up
    // to the shared /2), and |F(q^2)|^2 is bounded algebraically by its
    // value at the pole, 1/(mRho^2 GammaRho^2), giving a fixed rejection
    // envelope with no grid search needed (unlike PlutoEtaPrimeLLPiPi.h's
    // two-resonance interference, this is a single resonance).
    template <class Flat>
    std::array<TLorentzVector, 3> singleDalitz(
        const TLorentzVector& parent, double m, double rhoMass, double rhoGamma, Flat flat) {
      const double M = parent.M(), M2 = M * M;
      if (!(m > 0 && M > 2 * m))
        throw std::runtime_error("Single Dalitz masses are below threshold");
      const double rhoMass2 = rhoMass * rhoMass;
      const double poleBound = rhoMass2 * rhoGamma * rhoGamma;  // 1/formFactorBound
      const double low = 4 * m * m / M2, high = 1.;
      const double logRange = std::log(high / low), pi = std::acos(-1.);
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double x = low * std::exp(logRange * flat());
        const double beta2 = 1 - 4 * m * m / (x * M2);
        if (!(beta2 > 0))
          continue;
        const double beta = std::sqrt(beta2);
        const double c = 2 * flat() - 1;
        const double q2 = x * M2;
        const double poleDiff = rhoMass2 - q2;
        const double formFactor = poleBound / (poleDiff * poleDiff + poleBound);  // in (0,1]
        const double weight =
            beta * std::pow(1 - x, 3) * (1 + c * c + (1 - beta2) * (1 - c * c)) / 2 * formFactor;
        if (!std::isfinite(weight) || weight < 0 || weight > 1)
          throw std::runtime_error("Invalid single Dalitz rejection weight");
        if (flat() >= weight)
          continue;
        const double q = M * std::sqrt(x);
        const double k = M * (1 - x) / 2;  // photon energy/momentum in the parent frame
        const double pl = q * beta / 2;
        const double s = std::sqrt(1 - c * c), phi = 2 * pi * flat();
        TLorentzVector gamma(0, 0, k, k);
        TLorentzVector lMinus(pl * s * std::cos(phi), pl * s * std::sin(phi), pl * c, q / 2);
        TLorentzVector lPlus(-pl * s * std::cos(phi), -pl * s * std::sin(phi), -pl * c, q / 2);
        lMinus.Boost(0, 0, -k / (M - k));
        lPlus.Boost(0, 0, -k / (M - k));
        const double theta = std::acos(2 * flat() - 1), orientPhi = 2 * pi * flat();
        std::array<TLorentzVector, 3> out = {{lMinus, lPlus, gamma}};
        for (auto& p : out) {
          p.RotateY(theta);
          p.RotateZ(orientPhi);
          p.Boost(parent.BoostVector());
        }
        return out;
      }
      throw std::runtime_error("Single Dalitz rejection limit exhausted");
    }

    // eta (221) only: same kinematic/angular structure as singleDalitz, but
    // the constant/rho-pole form factor is replaced by the actual
    // data-fitted eta TFF (PlutoEtaTFF.h, arXiv:1504.07742 Appendix A) --
    // verified monotonic and finite across eta's whole kinematic range, so
    // formFactorBound (from etaPadeSingleBound, computed once by the
    // caller) is a safe rejection envelope.
    template <class Flat>
    std::array<TLorentzVector, 3> singleDalitzEtaPade(
        const TLorentzVector& parent, double m, double formFactorBound, Flat flat) {
      const double M = parent.M(), M2 = M * M;
      if (!(m > 0 && M > 2 * m))
        throw std::runtime_error("Single Dalitz masses are below threshold");
      const double low = 4 * m * m / M2, high = 1.;
      const double logRange = std::log(high / low), pi = std::acos(-1.);
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double x = low * std::exp(logRange * flat());
        const double beta2 = 1 - 4 * m * m / (x * M2);
        if (!(beta2 > 0))
          continue;
        const double beta = std::sqrt(beta2);
        const double c = 2 * flat() - 1;
        const double q2 = x * M2;
        const double formFactor = std::pow(etaPadeFormFactorNorm(q2), 2) / formFactorBound;
        const double weight =
            beta * std::pow(1 - x, 3) * (1 + c * c + (1 - beta2) * (1 - c * c)) / 2 * formFactor;
        if (!std::isfinite(weight) || weight < 0 || weight > 1)
          throw std::runtime_error("Invalid single Dalitz (eta Pade) rejection weight -- raise the bound's safety margin");
        if (flat() >= weight)
          continue;
        const double q = M * std::sqrt(x);
        const double k = M * (1 - x) / 2;
        const double pl = q * beta / 2;
        const double s = std::sqrt(1 - c * c), phi = 2 * pi * flat();
        TLorentzVector gamma(0, 0, k, k);
        TLorentzVector lMinus(pl * s * std::cos(phi), pl * s * std::sin(phi), pl * c, q / 2);
        TLorentzVector lPlus(-pl * s * std::cos(phi), -pl * s * std::sin(phi), -pl * c, q / 2);
        lMinus.Boost(0, 0, -k / (M - k));
        lPlus.Boost(0, 0, -k / (M - k));
        const double theta = std::acos(2 * flat() - 1), orientPhi = 2 * pi * flat();
        std::array<TLorentzVector, 3> out = {{lMinus, lPlus, gamma}};
        for (auto& p : out) {
          p.RotateY(theta);
          p.RotateZ(orientPhi);
          p.Boost(parent.BoostVector());
        }
        return out;
      }
      throw std::runtime_error("Single Dalitz (eta Pade) rejection limit exhausted");
    }

    // Same kinematic/angular structure as singleDalitz/singleDalitzEtaPade
    // for an arbitrary form factor: formFactorSq(q2) must return
    // |F(q^2)/F(0)|^2 already divided by a rejection bound, i.e. in [0,1].
    // (The two variants above predate this and carry their own copies of
    // the same body; they're left untouched rather than migrated.)
    template <class FormFactorSq, class Flat>
    std::array<TLorentzVector, 3> singleDalitzGeneric(
        const TLorentzVector& parent, double m, FormFactorSq formFactorSq, Flat flat) {
      const double M = parent.M(), M2 = M * M;
      if (!(m > 0 && M > 2 * m))
        throw std::runtime_error("Single Dalitz masses are below threshold");
      const double low = 4 * m * m / M2, high = 1.;
      const double logRange = std::log(high / low), pi = std::acos(-1.);
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double x = low * std::exp(logRange * flat());
        const double beta2 = 1 - 4 * m * m / (x * M2);
        if (!(beta2 > 0))
          continue;
        const double beta = std::sqrt(beta2);
        const double c = 2 * flat() - 1;
        const double weight =
            beta * std::pow(1 - x, 3) * (1 + c * c + (1 - beta2) * (1 - c * c)) / 2 * formFactorSq(x * M2);
        if (!std::isfinite(weight) || weight < 0 || weight > 1)
          throw std::runtime_error("Invalid single Dalitz rejection weight -- raise the form-factor bound's safety margin");
        if (flat() >= weight)
          continue;
        const double q = M * std::sqrt(x);
        const double k = M * (1 - x) / 2;
        const double pl = q * beta / 2;
        const double s = std::sqrt(1 - c * c), phi = 2 * pi * flat();
        TLorentzVector gamma(0, 0, k, k);
        TLorentzVector lMinus(pl * s * std::cos(phi), pl * s * std::sin(phi), pl * c, q / 2);
        TLorentzVector lPlus(-pl * s * std::cos(phi), -pl * s * std::sin(phi), -pl * c, q / 2);
        lMinus.Boost(0, 0, -k / (M - k));
        lPlus.Boost(0, 0, -k / (M - k));
        const double theta = std::acos(2 * flat() - 1), orientPhi = 2 * pi * flat();
        std::array<TLorentzVector, 3> out = {{lMinus, lPlus, gamma}};
        for (auto& p : out) {
          p.RotateY(theta);
          p.RotateZ(orientPhi);
          p.Boost(parent.BoostVector());
        }
        return out;
      }
      throw std::runtime_error("Single Dalitz rejection limit exhausted");
    }

    // eta' (331) only: the two-regime Pade + rho/omega/phi VMD transition
    // form factor of PlutoEtaPrimeTFF.h. formFactorBound comes from
    // etaPrimeSingleBound, computed once by the caller.
    template <class Flat>
    std::array<TLorentzVector, 3> singleDalitzEtaPrime(
        const TLorentzVector& parent, double m, double formFactorBound, Flat flat) {
      return singleDalitzGeneric(
          parent, m, [formFactorBound](double q2) { return etaPrimeFormFactorNormSq(q2) / formFactorBound; }, flat);
    }
  }  // namespace pluto
}  // namespace gen
#endif
