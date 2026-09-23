#ifndef GeneratorInterface_Pythia8Interface_PlutoLLPiPiChPT_h
#define GeneratorInterface_Pythia8Interface_PlutoLLPiPiChPT_h
#include "GeneratorInterface/Pythia8Interface/interface/PlutoEtaPrimeLLPiPi.h"
#include "TLorentzVector.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace gen {
  namespace pluto {
    // eta/eta' -> l+ l- pi+ pi-, mass-dependence from Zillinger, Kubis,
    // Sanchez-Puertas, arXiv:2210.14925, Section 3 (their Standard-Model
    // amplitude, not the CP-violating BSM part of that paper). Unlike
    // PlutoPiPiDilepton.h (constant F=1, eta only) and
    // PlutoEtaPrimeLLPiPi.h (EvtGen-matched, eta' only), this is the first
    // model here validated against *measured* branching ratios for both
    // parents (their Table 1: eta/eta' -> pi+pi-e+e- agree with KLOE/BESIII
    // data at the few-percent level) and it covers eta, which had no
    // resonance-aware reference before.
    //
    // Their amplitude factorizes as f_1(s,s_l) = P(s) * Omega_1^1(s) * Fbar(s_l),
    // s = pion-pair mass^2, s_l = dilepton mass^2:
    //   - P(s): a polynomial fit to eta/eta' -> pi+pi-gamma data (Eq. 3.4),
    //     ported here with their exact coefficients.
    //   - Fbar(s_l): a coherent two-resonance (rho, rho(1450)) dilepton form
    //     factor (Eq. 3.5), ported exactly; rho(1450) mass/width are PDG
    //     values, not given explicitly in the paper's text.
    //   - Omega_1^1(s): the P-wave pi-pi Omnes function, a dispersive
    //     integral over the measured pi-pi phase shift. The paper does not
    //     give a closed form for it -- its own authors obtained it from
    //     unpublished work (acknowledgments: "Akdag and Isken"). APPROXIMATED
    //     here by a single energy-dependent-width rho Breit-Wigner (the same
    //     ansatz already used for Fbar's own resonances), which the paper
    //     itself motivates by noting the pi-pi system is "fully dominated by
    //     the rho resonance at the energies of interest" -- NOT the rigorous
    //     dispersive result.
    //
    // The angular (helicity-angle) sampling below is NOT re-derived from
    // this paper's epsilon-tensor amplitude -- that would need the full
    // squared matrix element, a separate and nontrivial derivation this
    // pass does not attempt. It reuses pipiMagneticPointlike's
    // already-cross-checked-against-EvtGen angular structure unchanged, and
    // only replaces the mass-dependence (which was constant, F=1, before).
    // This is a hybrid, not a first-principles port of the full amplitude.
    inline std::complex<double> llpipiPolynomial(double s, bool isEtaPrime) {
      if (!isEtaPrime) {
        constexpr double A = 17.9, alpha = 1.52;  // GeV^-3, GeV^-2
        return A * (1 + alpha * s);
      }
      constexpr double A = 16.7, alpha = 1.00, beta = -0.55, kappa2 = 6.72e-3;
      constexpr double mOmega = 0.78266, gammaOmega = 0.00868;  // PDG
      const std::complex<double> omegaPole =
          kappa2 / std::complex<double>(mOmega * mOmega - s, -mOmega * gammaOmega);
      return A * (1.0 + alpha * s + beta * s * s + omegaPole);
    }

    // Stand-in for the true Omnes function Omega_1^1(s); see header comment.
    inline std::complex<double> llpipiOmegaApprox(double s, double mPi, double rhoMass, double rhoGamma) {
      const double width = etaPrimeRhoWidth(s, mPi, rhoMass, rhoGamma);
      return 1.0 / std::complex<double>(rhoMass * rhoMass - s, -std::sqrt(s) * width);
    }

    // Eq. 3.5: coherent rho + rho(1450) dilepton form factor.
    inline std::complex<double> llpipiFbar(
        double sl, double mLep, double rhoMass, double rhoGamma, double rhopMass, double rhopGamma) {
      const double widthRho = etaPrimeRhoWidth(sl, mLep, rhoMass, rhoGamma);
      const double widthRhop = etaPrimeRhoWidth(sl, mLep, rhopMass, rhopGamma);
      const std::complex<double> dRho(rhoMass * rhoMass - sl, -std::sqrt(sl) * widthRho);
      const std::complex<double> dRhop(rhopMass * rhopMass - sl, -std::sqrt(sl) * widthRhop);
      return (rhoMass * rhoMass * rhopMass * rhopMass) / (dRho * dRhop);
    }

    inline double llpipiWeight(double s, double sl, double mPi, double mLep, bool isEtaPrime,
                                double rhoMass, double rhoGamma, double rhopMass, double rhopGamma) {
      const auto p = llpipiPolynomial(s, isEtaPrime);
      const auto om = llpipiOmegaApprox(s, mPi, rhoMass, rhoGamma);
      const auto fb = llpipiFbar(sl, mLep, rhoMass, rhoGamma, rhopMass, rhopGamma);
      return std::norm(p) * std::norm(om) * std::norm(fb);
    }

    // Grid-search the weight's maximum over the physical (s,sl) domain, same
    // technique and same caveat as PlutoEtaPrimeLLPiPi.h's etaPrimeF0Bound:
    // this rational/polynomial product has no simple closed-form bound.
    inline double llpipiWeightBound(double mParent, double mPi, double mLep, bool isEtaPrime, double rhoMass,
                                     double rhoGamma, double rhopMass, double rhopGamma, unsigned grid = 300) {
      const double lowS = 4 * mPi * mPi, highS = std::pow(mParent - 2 * mLep, 2);
      const double lowSl = 4 * mLep * mLep, highSl = std::pow(mParent - 2 * mPi, 2);
      double maxW = 0.;
      for (unsigned i = 0; i <= grid; ++i) {
        const double s = lowS + (highS - lowS) * i / grid;
        for (unsigned j = 0; j <= grid; ++j) {
          const double sl = lowSl + (highSl - lowSl) * j / grid;
          if (std::sqrt(s) + std::sqrt(sl) >= mParent)
            continue;
          maxW = std::max(
              maxW, llpipiWeight(s, sl, mPi, mLep, isEtaPrime, rhoMass, rhoGamma, rhopMass, rhopGamma));
        }
      }
      if (!(maxW > 0.))
        throw std::runtime_error("LLPiPi ChPT weight grid search found no valid phase space point");
      return 1.2 * maxW;
    }

    template <class Flat>
    std::array<TLorentzVector, 4> llpipiChPT(const TLorentzVector& parent, double ml, double mpi, bool isEtaPrime,
                                              double rhoMass, double rhoGamma, double rhopMass, double rhopGamma,
                                              double weightBound, Flat flat) {
      const double M = parent.M(), M2 = M * M, smin = 4 * mpi * mpi, qmin = 4 * ml * ml;
      if (!(ml > 0 && mpi > 0 && M > 2 * (ml + mpi)))
        throw std::runtime_error("LLPiPi ChPT masses are below threshold");
      const double smax = std::pow(M - 2 * ml, 2), qmax = std::pow(M - 2 * mpi, 2);
      const double logq = std::log(qmax / qmin), pi = std::acos(-1.);
      const auto lambda = [](double x, double y) { return std::max(0., std::pow(1 - x - y, 2) - 4 * x * y); };
      const double phaseBound = std::pow(lambda(smin / M2, qmin / M2), 1.5) * (smax / M2) *
                                 std::pow(1 - smin / smax, 1.5) * std::sqrt(1 - qmin / qmax);
      const double bound = phaseBound * weightBound;
      if (!(bound > 0) || !std::isfinite(bound))
        throw std::runtime_error("Invalid LLPiPi ChPT envelope");
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double s = smin + (smax - smin) * flat(), q = qmin * std::exp(logq * flat());
        if (std::sqrt(s) + std::sqrt(q) >= M)
          continue;
        const double bp = 1 - smin / s, bl = 1 - qmin / q, lam = lambda(s / M2, q / M2);
        const double cp = 2 * flat() - 1, cl = 2 * flat() - 1, phi = 2 * pi * flat();
        const double massWeight = llpipiWeight(s, q, mpi, ml, isEtaPrime, rhoMass, rhoGamma, rhopMass, rhopGamma);
        const double weight = std::pow(lam, 1.5) * (s / M2) * std::pow(bp, 1.5) * std::sqrt(bl) * (1 - cp * cp) *
                               (1 - bl * (1 - cl * cl) * std::pow(std::sin(phi), 2)) * massWeight / bound;
        if (!std::isfinite(weight) || weight < 0 || weight > 1)
          throw std::runtime_error("LLPiPi ChPT rejection envelope violated (raise the grid-search safety margin)");
        if (flat() >= weight)
          continue;
        const double k = M * std::sqrt(lam) / 2;
        const double el = (M2 + q - s) / (2 * M), ep = (M2 + s - q) / (2 * M);
        const double pl = std::sqrt(q * bl) / 2, pp = std::sqrt(s * bp) / 2;
        const double tl = std::sqrt(1 - cl * cl), tp = std::sqrt(1 - cp * cp), az = 2 * pi * flat();
        std::array<TLorentzVector, 4> out = {{
            TLorentzVector(pl * tl * std::cos(az + phi), pl * tl * std::sin(az + phi), pl * cl, std::sqrt(q) / 2),
            TLorentzVector(-pl * tl * std::cos(az + phi), -pl * tl * std::sin(az + phi), -pl * cl, std::sqrt(q) / 2),
            TLorentzVector(pp * tp * std::cos(az), pp * tp * std::sin(az), pp * cp, std::sqrt(s) / 2),
            TLorentzVector(-pp * tp * std::cos(az), -pp * tp * std::sin(az), -pp * cp, std::sqrt(s) / 2)}};
        const double theta = std::acos(2 * flat() - 1), orientPhi = 2 * pi * flat();
        for (unsigned i = 0; i < 4; ++i) {
          out[i].Boost(0, 0, i < 2 ? k / el : -k / ep);
          out[i].RotateY(theta);
          out[i].RotateZ(orientPhi);
          out[i].Boost(parent.BoostVector());
        }
        return out;
      }
      throw std::runtime_error("LLPiPi ChPT rejection limit exhausted");
    }
  }  // namespace pluto
}  // namespace gen
#endif
