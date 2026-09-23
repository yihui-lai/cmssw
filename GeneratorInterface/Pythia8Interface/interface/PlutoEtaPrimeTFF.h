#ifndef GeneratorInterface_Pythia8Interface_PlutoEtaPrimeTFF_h
#define GeneratorInterface_Pythia8Interface_PlutoEtaPrimeTFF_h
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace gen {
  namespace pluto {
    // eta' (331) single-virtual transition form factor over its whole
    // time-like range 0 <= q^2 <= M_eta'^2, following the prescription of
    // Escribano & Gonzalez-Solis, arXiv:1511.04916, Section 2 ("eta' ->
    // gamma gamma*"): a two-regime model, because eta' -- unlike eta --
    // reaches past the pole of its own Pade approximant.
    //
    //   sqrt(q^2) <= 0.70 GeV : eta'-specific Pade approximant P_1^6, fitted
    //     to space-like data (Escribano/Masjuan/Sanchez-Puertas,
    //     arXiv:1307.2061, Table IV, eta' column). Its pole sits at
    //     q^2 = 1/r1 = 0.6935 GeV^2 (sqrt = 0.833 GeV, inside the paper's
    //     quoted (0.83,0.86) GeV) -- inside eta''s accessible range up to
    //     M_eta' = 0.958 GeV, which is exactly why a single Pade formula (as
    //     used for eta in PlutoEtaTFF.h, whose pole is above M_eta) cannot
    //     cover eta'.
    //   sqrt(q^2) >  0.70 GeV : coherent rho + omega + phi vector-meson
    //     dominance (Landsberg's finite-width VMD, 1511.04916 Eq. (VMD)),
    //       F(q^2) = [sum_V w_V]^-1 sum_V w_V M_V^2/(M_V^2 - q^2 - i M_V G_V(q^2)),
    //     with the paper's energy-dependent rho width and constant omega/phi
    //     widths, rescaled by one real constant so |F| is continuous at the
    //     0.70 GeV seam (the paper's own matching prescription). This is
    //     where the sharp omega peak near sqrt(q^2) = 0.78 GeV comes from --
    //     the piece the previous single-rho-pole eta' model was missing.
    //
    // WHAT IS FROM THE PAPERS vs. WHAT IS DERIVED HERE (be honest about the
    // difference when trusting this):
    //  - Pade coefficients: exact, from arXiv:1307.2061 Table IV. That table
    //    quotes them to only 1-3 significant figures (t2 = 0.007 has ONE), and
    //    the authors' footnote says full precision is available on request --
    //    a ~7% uncertainty on t2 alone moves the slope by well under 1%, and
    //    the Pade only sets the shape below 0.70 GeV, so this is accepted.
    //  - Masses/widths: PDG 2016 values (the paper's own source, its ref [1]),
    //    the energy-dependent rho width formula from the paper's footnote 5.
    //  - VMD weights w_V = g_{V eta' gamma} / (2 g_{V gamma}): the paper takes
    //    these "from the PDG", but neither rho->eta' gamma nor omega->eta'
    //    gamma is kinematically allowed (M_rho, M_omega < M_eta'), so they
    //    cannot be read off a measured branching ratio. They are DERIVED here
    //    from the quark-flavor-basis mixing model instead: eta' = sin(phi)
    //    eta_q + cos(phi) eta_s with phi = 39.6 deg (the radiative-decay-
    //    based value quoted in arXiv:1307.2061, Sec. III; FKS gives 39.3),
    //    ideal rho/omega/phi flavor content (rho=(uu-dd)/sqrt2,
    //    omega=(uu+dd)/sqrt2, phi=ss), and photon coupling by quark charge:
    //      w_V ~ (sum_q e_q V_q P_q) * (sum_q e_q V_q)
    //      => w_rho : w_omega : w_phi = sin(phi)/(2 sqrt2) : sin(phi)/(18 sqrt2) : cos(phi)/9
    //    (overall scale is irrelevant: the prefactor normalizes F(0) = 1).
    //    Cross-checked, not merely asserted: the same construction for eta
    //    (eta = cos(phi) eta_q - sin(phi) eta_s) predicts a slope
    //    b = 1.87 GeV^-2 against the data-fitted Pade slope 1.92 GeV^-2
    //    (PlutoEtaTFF.h), and for eta' b = 1.48 vs. the Table IV Pade's 1.42
    //    -- ~2% and ~4% agreement. That slope test constrains the
    //    (rho+omega) : phi split and the mixing angle, but NOT the rho : omega
    //    split (9:1 here, the (e_u-e_d)^2 : (e_u+e_d)^2 isospin ratio), since
    //    rho and omega are nearly degenerate; the omega weight is what sets
    //    the height of the narrow peak, and is checked separately against the
    //    digitized eta' -> l+l-gamma spectra (test/pluto_validation/csv/).
    //
    // NOT covered: only |F|^2 is provided (everything the single-Dalitz and
    // the factorized eta' -> l+l-l'+l'- weights need). The 4e/4mu identical-
    // fermion channels need the complex amplitude at the amplitude level,
    // and the Pade (real) / VMD (complex, ~39 deg at the seam from the rho
    // alone) pieces are only magnitude-matched by the paper's prescription,
    // so a phase-continuity choice would have to be made -- deliberately not
    // done here.
    namespace etaprime_tff {
      constexpr double kPi = 3.14159265358979323846;
      // arXiv:1307.2061 Table IV, eta' column (P_1^6): F(Q^2) = P(Q^2)/(1 + r1 Q^2),
      // P(Q^2) = t1 + t2 Q^2 + ... + t6 Q^10, Q^2 = -q^2 [GeV^2, GeV^-1 for F].
      constexpr double t1 = 0.343, t2 = 0.007, t3 = 0.986e-3, t4 = 0.744e-4, t5 = 0.252e-5, t6 = 0.290e-7, r1 = 1.442;
      constexpr double matchSqrtS = 0.70;  // GeV, the paper's Pade -> VMD matching point
      constexpr double matchQ2 = matchSqrtS * matchSqrtS;
      // PDG 2016 (the paper's ref [1]).
      constexpr double mRho = 0.77526, gRho = 0.1491;
      constexpr double mOmega = 0.78265, gOmega = 0.00849;
      constexpr double mPhi = 1.019461, gPhi = 0.004266;
      constexpr double mPi = 0.13957;
      constexpr double mixingAngleDeg = 39.6;

      struct Weights {
        double rho, omega, phi;
      };
      inline const Weights& weights() {
        static const Weights w = [] {
          const double phi = mixingAngleDeg * kPi / 180., s = std::sin(phi), c = std::cos(phi), r2 = std::sqrt(2.);
          return Weights{s / (2 * r2), s / (18 * r2), c / 9};
        }();
        return w;
      }

      // Paper's footnote: Gamma_rho(q^2) = Gamma_rho (q^2/M_rho^2) sigma^3(q^2)/sigma^3(M_rho^2).
      // Zero below the two-pion threshold, where the rho cannot decay.
      inline double rhoWidth(double q2) {
        const double threshold = 4 * mPi * mPi;
        if (q2 <= threshold)
          return 0.;
        const double sigma = std::sqrt(1 - threshold / q2), sigma0 = std::sqrt(1 - threshold / (mRho * mRho));
        return gRho * (q2 / (mRho * mRho)) * std::pow(sigma / sigma0, 3);
      }

      inline double padeNorm(double q2) {
        const double Q2 = -q2;
        const double num = t1 + Q2 * (t2 + Q2 * (t3 + Q2 * (t4 + Q2 * (t5 + Q2 * t6))));
        const double den = 1 + r1 * Q2;
        if (!(den > 0))
          throw std::runtime_error("eta' Pade form factor evaluated at or past its fitted pole");
        return num / den / t1;
      }

      inline std::complex<double> vmdNorm(double q2) {
        const auto& w = weights();
        const auto term = [q2](double weight, double mass, double width) {
          return weight * mass * mass / std::complex<double>(mass * mass - q2, -mass * width);
        };
        return (term(w.rho, mRho, rhoWidth(q2)) + term(w.omega, mOmega, gOmega) + term(w.phi, mPhi, gPhi)) /
               (w.rho + w.omega + w.phi);
      }

      // Real rescale of the VMD piece so |F| is continuous at the seam.
      inline double matchScale() {
        static const double r = padeNorm(matchQ2) / std::abs(vmdNorm(matchQ2));
        return r;
      }
    }  // namespace etaprime_tff

    // |F(q^2)/F(0)|^2 for eta', q^2 in GeV^2, time-like, 0 <= q^2 <= M_eta'^2.
    inline double etaPrimeFormFactorNormSq(double q2) {
      using namespace etaprime_tff;
      if (q2 <= matchQ2) {
        const double f = padeNorm(q2);
        return f * f;
      }
      const double r = matchScale();
      return r * r * std::norm(vmdNorm(q2));
    }

    // Maximum of |F|^2 on [lo, hi]: dense grid (fine enough to resolve the
    // omega's ~0.0066 GeV^2 full width many times over) followed by a
    // golden-section polish around the best grid cell, so the narrow peak
    // is neither missed nor undersampled.
    inline double etaPrimeMaxNormSq(double lo, double hi, unsigned grid = 20000) {
      double best = 0., bestQ2 = lo;
      const double step = (hi - lo) / grid;
      for (unsigned i = 0; i <= grid; ++i) {
        const double q2 = lo + step * i, f = etaPrimeFormFactorNormSq(q2);
        if (f > best) {
          best = f;
          bestQ2 = q2;
        }
      }
      double a = std::max(lo, bestQ2 - step), b = std::min(hi, bestQ2 + step);
      const double gr = 0.6180339887498949;
      for (unsigned it = 0; it < 80; ++it) {
        const double c = b - gr * (b - a), d = a + gr * (b - a);
        if (etaPrimeFormFactorNormSq(c) > etaPrimeFormFactorNormSq(d))
          b = d;
        else
          a = c;
      }
      return std::max(best, etaPrimeFormFactorNormSq(0.5 * (a + b)));
    }

    // Rejection-envelope bound for singleDalitzEtaPrime: |F|^2 is the only
    // factor without an algebraic bound of 1 there.
    inline double etaPrimeSingleBound(double mParent, double mLep, unsigned grid = 20000) {
      const double maxW = etaPrimeMaxNormSq(4 * mLep * mLep, mParent * mParent, grid);
      if (!(maxW > 0.))
        throw std::runtime_error("eta' single bound search found no valid point");
      return 1.05 * maxW;
    }

    // Rejection-envelope bound for the factorized double-virtual weight
    // |F(s1)|^2 |F(s2)|^2 over the joint physical domain sqrt(s1)+sqrt(s2)<M.
    // Exact on the grid: for each s1, the best partner s2 is the running
    // maximum of |F|^2 over every s2 the kinematics still allows.
    inline double etaPrimeDoubleBound(double mParent, double mLep1, double mLep2, unsigned grid = 20000) {
      const double low1 = 4 * mLep1 * mLep1, high1 = std::pow(mParent - 2 * mLep2, 2);
      const double low2 = 4 * mLep2 * mLep2, high2 = std::pow(mParent - 2 * mLep1, 2);
      std::vector<double> runMax2(grid + 1);
      double running = 0.;
      for (unsigned j = 0; j <= grid; ++j) {
        running = std::max(running, etaPrimeFormFactorNormSq(low2 + (high2 - low2) * j / grid));
        runMax2[j] = running;
      }
      double best = 0.;
      for (unsigned i = 0; i <= grid; ++i) {
        const double s1 = low1 + (high1 - low1) * i / grid;
        const double s2max = std::pow(mParent - std::sqrt(s1), 2);
        if (s2max <= low2)
          continue;
        const double pos = std::ceil((std::min(s2max, high2) - low2) / (high2 - low2) * grid);
        best = std::max(best, etaPrimeFormFactorNormSq(s1) * runMax2[std::min<unsigned>(grid, static_cast<unsigned>(pos))]);
      }
      if (!(best > 0.))
        throw std::runtime_error("eta' double bound search found no valid point");
      return 1.05 * best;
    }
  }  // namespace pluto
}  // namespace gen
#endif
