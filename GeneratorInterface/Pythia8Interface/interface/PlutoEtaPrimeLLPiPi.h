#ifndef GeneratorInterface_Pythia8Interface_PlutoEtaPrimeLLPiPi_h
#define GeneratorInterface_Pythia8Interface_PlutoEtaPrimeLLPiPi_h
#include "TLorentzVector.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace gen {
  namespace pluto {
    // eta' -> l+ l- pi+ pi-, matching EvtGen's EvtEtaLLPiPi model (Zhang
    // Zhen-Yu et al., Chinese Phys. C 36, 926 (2012), Eqs. 3, 6-9). Unlike
    // pipiMagneticPointlike's constant form factor, this includes the
    // coherent two-resonance (rho0 in the dilepton channel and rho0 in the
    // dipion channel) amplitude F0. EvtGen provides no eta (221) version of
    // this model, so eta keeps the pipiMagneticPointlike baseline; this is
    // eta'-only. rhoMass/rhoGamma are taken from the caller's particle data
    // table, matching EvtGen's EvtPDL::getMeanMass/getWidth("rho0").
    //
    // The angular part reproduces pipiMagneticPointlike's shape exactly (up
    // to the s_pipi*beta_pipi^2 normalization folded into the phase-space
    // weight below) -- that agreement is a cross-check that the earlier,
    // independently-derived constant-form-factor baseline has the right
    // "magnetic" angular structure. Only F0 is new here.
    inline double etaPrimeRhoWidth(double s, double m, double rhoMass, double rhoGamma) {
      const double rhoMass2 = rhoMass * rhoMass;
      if (s < 4 * m * m)
        return 0.;
      const double num = 1.0 - 4 * m * m / s;
      const double den = 1.0 - 4 * m * m / rhoMass2;
      const double ratio = den > 0. ? num / den : 0.;
      return rhoGamma * (s / rhoMass2) * std::pow(std::max(0., ratio), 1.5);
    }

    // Eq. 7, with the fixed c1=1, c2=0 (hence c3=1) used by EvtGen's
    // default construction: par1=-0.5, parLL=0, parPiPi=1.5.
    inline double etaPrimeF0(double sLL, double sPiPi, double mLep, double mPi, double rhoMass, double rhoGamma) {
      constexpr double par1 = -0.5, parPiPi = 1.5;
      const double rhoMass2 = rhoMass * rhoMass;
      const double widthL = etaPrimeRhoWidth(sLL, mLep, rhoMass, rhoGamma);
      const double widthPi = etaPrimeRhoWidth(sPiPi, mPi, rhoMass, rhoGamma);
      const double dL = rhoMass2 - sLL, dPi = rhoMass2 - sPiPi;
      const double wL = rhoMass * widthL, wPi = rhoMass * widthPi;
      const double denomL = dL * dL + wL * wL, denomPi = dPi * dPi + wPi * wPi;
      if (!(denomL > 0.) || !(denomPi > 0.))
        return 0.;
      const double rho4 = rhoMass2 * rhoMass2, denomProd = denomL * denomPi;
      const double re = par1 + parPiPi * rho4 * (dPi * dL - wL * wPi) / denomProd;
      const double im = parPiPi * rho4 * (wPi * dL + wL * dPi) / denomProd;
      return re * re + im * im;
    }

    // Grid-search the maximum of etaPrimeF0 over the physical (sLL,sPiPi)
    // domain, for use as a rejection-sampling envelope. Unlike the fixed
    // algebraic bounds elsewhere in this file, F0's rational-function
    // maximum near the rho pole has no simple closed form, so this scans a
    // grid and applies a safety margin; call once (e.g. at decayer
    // construction), not per event.
    inline double etaPrimeF0Bound(
        double mEta, double mLep, double mPi, double rhoMass, double rhoGamma, unsigned grid = 400) {
      const double M2 = mEta * mEta;
      const double lowLL = 4 * mLep * mLep, highLL = std::pow(mEta - 2 * mPi, 2);
      const double lowPiPi = 4 * mPi * mPi, highPiPi = std::pow(mEta - 2 * mLep, 2);
      double maxF0 = 0.;
      for (unsigned i = 0; i <= grid; ++i) {
        const double sLL = lowLL + (highLL - lowLL) * i / grid;
        for (unsigned j = 0; j <= grid; ++j) {
          const double sPiPi = lowPiPi + (highPiPi - lowPiPi) * j / grid;
          if (std::sqrt(sLL) + std::sqrt(sPiPi) >= mEta)
            continue;
          const double lambda = std::pow(1 - sLL / M2 - sPiPi / M2, 2) - 4 * sLL * sPiPi / (M2 * M2);
          if (lambda <= 0.)
            continue;
          maxF0 = std::max(maxF0, etaPrimeF0(sLL, sPiPi, mLep, mPi, rhoMass, rhoGamma));
        }
      }
      if (!(maxF0 > 0.))
        throw std::runtime_error("eta' LLPiPi form-factor grid search found no valid phase space point");
      return 1.2 * maxF0;  // safety margin around the grid resolution
    }

    template <class Flat>
    std::array<TLorentzVector, 4> etaPrimeLLPiPi(const TLorentzVector& parent,
                                                  double ml,
                                                  double mpi,
                                                  double rhoMass,
                                                  double rhoGamma,
                                                  double f0Bound,
                                                  Flat flat) {
      const double M = parent.M(), M2 = M * M, smin = 4 * mpi * mpi, qmin = 4 * ml * ml;
      if (!(ml > 0 && mpi > 0 && M > 2 * (ml + mpi)))
        throw std::runtime_error("eta' LLPiPi masses are below threshold");
      const double smax = std::pow(M - 2 * ml, 2), qmax = std::pow(M - 2 * mpi, 2);
      const double logq = std::log(qmax / qmin), pi = std::acos(-1.);
      const auto lambda = [](double x, double y) { return std::max(0., std::pow(1 - x - y, 2) - 4 * x * y); };
      const double bound = std::pow(lambda(smin / M2, qmin / M2), 1.5) * (smax / M2) *
                            std::pow(1 - smin / smax, 1.5) * std::sqrt(1 - qmin / qmax) * f0Bound;
      if (!(bound > 0) || !std::isfinite(bound))
        throw std::runtime_error("Invalid eta' LLPiPi envelope");
      for (unsigned attempt = 0; attempt < 1000000; ++attempt) {
        const double s = smin + (smax - smin) * flat(), q = qmin * std::exp(logq * flat());
        if (std::sqrt(s) + std::sqrt(q) >= M)
          continue;
        const double bp = 1 - smin / s, bl = 1 - qmin / q, lam = lambda(s / M2, q / M2);
        const double cp = 2 * flat() - 1, cl = 2 * flat() - 1, phi = 2 * pi * flat();
        const double f0 = etaPrimeF0(q, s, ml, mpi, rhoMass, rhoGamma);
        const double weight = std::pow(lam, 1.5) * (s / M2) * std::pow(bp, 1.5) * std::sqrt(bl) *
                               (1 - cp * cp) * (1 - bl * (1 - cl * cl) * std::pow(std::sin(phi), 2)) * f0 / bound;
        if (!std::isfinite(weight) || weight < 0 || weight > 1)
          throw std::runtime_error("eta' LLPiPi rejection envelope violated (raise the grid-search safety margin)");
        if (flat() >= weight)
          continue;
        const double k = M * std::sqrt(lam) / 2;
        const double el = (M2 + q - s) / (2 * M), ep = (M2 + s - q) / (2 * M);
        const double pl = std::sqrt(q * bl) / 2, pp = std::sqrt(s * bp) / 2;
        const double tl = std::sqrt(1 - cl * cl), tp = std::sqrt(1 - cp * cp), az = 2 * pi * flat();
        // Output ordering is l-, l+, pi+, pi-, matching pipiMagneticPointlike.
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
      throw std::runtime_error("eta' LLPiPi rejection limit exhausted");
    }
  }  // namespace pluto
}  // namespace gen
#endif
