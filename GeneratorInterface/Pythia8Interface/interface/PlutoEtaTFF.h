#ifndef GeneratorInterface_Pythia8Interface_PlutoEtaTFF_h
#define GeneratorInterface_Pythia8Interface_PlutoEtaTFF_h
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace gen {
  namespace pluto {
    // eta (221)-only single-virtual transition form factor, data-fitted
    // rational (Pade) approximant from Escribano, Masjuan, Sanchez-Puertas,
    // arXiv:1504.07742, Appendix A, Table 5 (their best P_1^7(Q^2) fit to
    // Q^2*F_etagammagamma*(Q^2), Q^2 = -q^2 in their space-like-positive
    // convention). This replaces the generic PDG-rho-pole guess used
    // elsewhere in this package, for eta specifically: it is the actual
    // data-validated shape of eta's own TFF, not a resonance-model stand-in.
    // eta' keeps the rho-pole model (this paper does not cover eta').
    //
    // Verified against the paper's own self-check (Eq. A.2): these exact
    // coefficients reproduce b_eta = 0.5749, matching the value quoted
    // there, before being used for anything here.
    //
    // The fitted pole (Q^2 = -1/r1 = -0.5106 GeV^2, i.e. time-like
    // q^2 = +0.5106 GeV^2) sits beyond eta's kinematic reach
    // (q^2_max = M_eta^2 = 0.300 GeV^2), so this stays real, finite and
    // -- checked numerically -- monotonically increasing across eta's
    // entire accessible range: no divergence, unlike a naive Taylor-series
    // extrapolation of b_eta/c_eta alone would give this far from q^2=0.
    inline double etaPadeFormFactor(double q2) {
      constexpr double t1 = 0.27349, t2 = 1.1771e-2, t3 = -1.1048e-3, t4 = 2.8861e-5, t5 = 2.2974e-6,
                       t6 = -1.5096e-7, t7 = 2.3655e-9, r1 = 1.9584;
      const double Q2 = -q2;  // their space-like-positive convention
      const double num = Q2 * (t1 + Q2 * (t2 + Q2 * (t3 + Q2 * (t4 + Q2 * (t5 + Q2 * (t6 + Q2 * t7))))));
      const double den = 1 + r1 * Q2;
      if (!(den > 0))
        throw std::runtime_error("eta Pade form factor evaluated past its fitted pole");
      if (std::abs(q2) < 1e-12)
        return t1;
      return num / den / Q2;  // NOT /q2: num/den = Q^2*F(Q^2), so dividing by
                               // Q2 (not q2=-Q2) recovers F(Q2). An earlier
                               // version divided by q2, introducing a sign
                               // bug caught by this file's own unit test
                               // (which reproduces arXiv:1504.07742 Eq. A.2's
                               // b_eta=0.5749 self-check) -- verify against
                               // that check again if this line is touched.
    }

    inline double etaPadeFormFactorNorm(double q2) { return etaPadeFormFactor(q2) / etaPadeFormFactor(0.); }

    // Grid-search bound for the single-virtual-photon weight |F(q^2)/F(0)|^2
    // (used by singleDalitz's eta branch). The function is confirmed
    // monotonically increasing over eta's domain, but a grid search (same
    // technique as the other resonant headers) is used rather than trusting
    // that analytically.
    inline double etaPadeSingleBound(double mParent, double mLep, unsigned grid = 400) {
      const double low = 4 * mLep * mLep, high = mParent * mParent;
      double maxW = 0.;
      for (unsigned i = 0; i <= grid; ++i) {
        const double q2 = low + (high - low) * i / grid;
        maxW = std::max(maxW, std::pow(etaPadeFormFactorNorm(q2), 2));
      }
      if (!(maxW > 0.))
        throw std::runtime_error("eta Pade single bound grid search found no valid point");
      return 1.2 * maxW;
    }

    // Grid-search bound for the double-virtual (factorized) weight
    // |F(s1)/F(0)|^2 * |F(s2)/F(0)|^2 over the joint physical domain
    // sqrt(s1)+sqrt(s2)<mParent (used by mixedPointlikeResonant's eta branch).
    inline double etaPadeDoubleBound(double mParent, double mLep1, double mLep2, unsigned grid = 300) {
      const double low1 = 4 * mLep1 * mLep1, high1 = std::pow(mParent - 2 * mLep2, 2);
      const double low2 = 4 * mLep2 * mLep2, high2 = std::pow(mParent - 2 * mLep1, 2);
      double maxW = 0.;
      for (unsigned i = 0; i <= grid; ++i) {
        const double s1 = low1 + (high1 - low1) * i / grid;
        for (unsigned j = 0; j <= grid; ++j) {
          const double s2 = low2 + (high2 - low2) * j / grid;
          if (std::sqrt(s1) + std::sqrt(s2) >= mParent)
            continue;
          maxW = std::max(maxW, std::pow(etaPadeFormFactorNorm(s1), 2) * std::pow(etaPadeFormFactorNorm(s2), 2));
        }
      }
      if (!(maxW > 0.))
        throw std::runtime_error("eta Pade double bound grid search found no valid point");
      return 1.2 * maxW;
    }
  }  // namespace pluto
}  // namespace gen
#endif
