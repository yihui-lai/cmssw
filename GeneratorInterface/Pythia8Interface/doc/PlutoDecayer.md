# PlutoDecayer: forced eta/eta' -> lepton(s)(+2pi/gamma) decays

`PlutoDecayer` (`interface/PlutoDecayer.h`, `src/PlutoDecayer.cc`) is a
`Pythia8::DecayHandler` registered with Pythia8 through the classic
`Pythia::setDecayPtr()` mechanism -- the same mechanism `BiasedTauDecayer`
uses. It does not depend on Pythia8's newer `Init:plugins` loader (available
only from Pythia 8.317), so it works with whatever Pythia8 version this
CMSSW release is pinned to, and it does not depend on the external Pluto
library: the decay densities only need ROOT's `TLorentzVector`.

## Settings

Added to `Pythia8::Settings` in `Py8InterfaceBase.cc`/`Py8HMC3InterfaceBase.cc`:

| Setting                    | Values                                                    | Default      |
|-----------------------------|------------------------------------------------------------|--------------|
| `Pluto:filter`              | `on`/`off`                                                 | `off`        |
| `Pluto:allowForcedDecay`    | `on`/`off`, must be explicitly `on`                        | `off`        |
| `Pluto:parent`              | `221` (eta) or `331` (eta')                                | `221`        |
| `Pluto:mode`                | `2mu`, `2mugamma`, `4e`, `2mu2e`, `4mu`, `2e2pi`, `2mu2pi`  | `2mu2e`      |
| `Pluto:model`               | `pointlike` or `phaseSpace`                                | `pointlike`  |

`Pluto:allowForcedDecay` is a deliberate safety interlock: once attached,
`PlutoDecayer` decays **every** parent of that species exclusively into the
requested channel, replacing Pythia's own branching-ratio table for that
particle. This is a forced-exclusive shape sample, not a
branching-ratio-normalized one -- apply the branching fraction separately
downstream.

Example (see `test/plutoDecayer_cfg.py` and
`Configuration/GenProduction/python/Eta_2mu2e.py`):

```python
processParameters = cms.vstring(
    'Pluto:filter = on',
    'Pluto:allowForcedDecay = on',
    'Pluto:parent = 221',
    'Pluto:mode = 2mu2e',
    'Pluto:model = pointlike',
)
```

## Physics models and where they actually come from

**None of this is Pluto's native physics code.** The package is named
`PlutoDecayer` only because it plugs into Pythia8 the way the original,
now-removed prototype (`test/externalPlugins/pluto/`, an `Init:plugins`-based
plugin linked against the real external Pluto library) did. That prototype's
own notes (`IDENTICAL_MODEL.md`, `PIPI_MODEL.md`, before removal) record that
Pluto's native eta decay code was checked and rejected: the native
`eta -> e+ e- pi+ pi-` model had concrete bugs (wrong array index feeding the
pion angular routine, the lepton-angle particle built from the wrong PParticle,
a sign error against the reference formula, a rejection envelope that grows
adaptively instead of being fixed), and the native `eta -> 4e` model was
deliberately *not* used as any kind of reference. Both of those native models
are also hardcoded to electrons only and to eta (221) only -- Pluto's
`eta_decays` plugin, as wired into that prototype, never constructed a muon
final state or an eta' (331) parent at all. So none of `4mu`, `2mu2e`,
`2mu2pi`, or anything eta'-related was ever backed by Pluto's native code in
this integration, and the one channel that was (`2e2pi`, native) was
explicitly disabled for being wrong. **"Pluto has an eta->4e/2e2pi model" is
not evidence for the muon or eta' channels, and is only weak evidence for the
electron/eta channels it names**, since the previous author's own review
didn't trust the eta->4e one either.

What's actually implemented is independently-derived, constant-form-factor,
leading-order QED kinematics from Petri's thesis, arXiv:1010.2378
(https://arxiv.org/abs/1010.2378), a general multi-body Dalitz-decay
reference, not a Pluto-specific one:

- **`4e`, `4mu`** (`interface/PlutoIdenticalDoubleDalitz.h`): eta/eta' ->
  gamma* gamma* -> 4 identical leptons. Includes the identical-fermion
  exchange interference (direct minus exchange amplitude, `|D-X|^2`),
  built explicitly from Dirac spinors and summed over the 16 external spin
  combinations, section 3.3 of the reference. Sampling mixes the direct and
  exchanged pointlike densities equally, then applies a coherent rejection
  bounded by 1 algebraically (`|D-X|^2 <= 2(|D|^2+|X|^2)`); no adaptive or
  empirically tuned envelope is used. **No longer constant-F=1**: each
  diagram's two virtual-photon propagators are now dressed with an
  amplitude-level form factor (`Amplitudes`' `ffAmplitude`), evaluated
  separately for the direct and exchange pairings' own invariant masses,
  *before* they interfere -- for **eta** this is the data-fitted Pade TFF
  from `interface/PlutoEtaTFF.h` (arXiv:1504.07742); for **eta'** it is
  still the generic rho0-pole model. See "Literature read and applied"
  below for why the amplitude-level (not probability-level) placement was
  needed, and why the existing `|D-X|^2<=2(|D|^2+|X|^2)` bound needed no
  re-derivation to support it (it holds for any complex D, X, dressed or
  not).
- **`2mu2e`** (`interface/PlutoMixedDoubleDalitz.h`): eta/eta' -> gamma*
  gamma* -> mu+ mu- e+ e-, Eq. (3.36). No identical-fermion exchange exists
  here since the two currents have distinguishable lepton flavours. The
  constant form factor is replaced by a **factorized** product of two
  single-virtual form factors (arXiv:1511.04916 Eq. 8's "standard
  factorisation ansatz") -- for **eta**, `mixedPointlikeEtaPade` uses the
  data-fitted Pade TFF (`interface/PlutoEtaTFF.h`, arXiv:1504.07742); for
  **eta'**, `mixedPointlikeResonant` uses the generic rho0-pole model.
- **`2e2pi`, `2mu2pi`** (`interface/PlutoLLPiPiChPT.h`): eta/eta' ->
  l+ l- pi+ pi-. As of the latest pass this is **no longer** the
  constant-form-factor Petri baseline -- it uses the mass-dependence
  (pion-pair-mass polynomial x pi-pi resonance factor x dilepton
  rho/rho(1450) form factor) from Zillinger, Kubis, Sanchez-Puertas,
  arXiv:2210.14925, combined with the original Petri/EvtGen-cross-checked
  "magnetic" angular shape from `interface/PlutoPiPiDilepton.h` (still
  present, still used for its angular formula and by
  `interface/PlutoEtaPrimeLLPiPi.h`, but no longer called directly by
  `PlutoDecayer` for the mass-dependence). See "Literature read and
  applied" and "Known gaps" below for exactly what is and isn't a faithful
  port of that paper.
- **`2mugamma`** (`interface/PlutoSingleDalitz.h`): eta/eta' -> gamma l+ l-,
  the standard Kroll-Wada single-Dalitz shape. **eta** uses
  `singleDalitzEtaPade`, dressed by the data-fitted Pade TFF
  (`interface/PlutoEtaTFF.h`, arXiv:1504.07742 Appendix A) instead of a
  constant; **eta'** uses `singleDalitz` with the generic rho0-pole model
  (matches EvtGen's `EvtEta2MuMuGamma`). See "Known gaps" below for the
  eta' omega-resonance gap this doesn't close.
- **`2mu`**: trivial isotropic two-body eta/eta' -> mu+ mu- (closed form, no
  rejection sampling -- a two-body decay has no shape to model).

None of these supply radiative corrections or a branching fraction -- they
are forced-exclusive kinematic baselines, not a fit to data or an absolute
rate prediction. `model = phaseSpace` instead uses ROOT's `TGenPhaseSpace`
(uniform phase space, no dynamics) as a kinematic control to check that the
pointlike results are not merely a phase-space artifact.

`interface/PlutoEtaPrimeLLPiPi.h` implements an eta'-only, EvtGen-`EvtEtaLLPiPi`-matched
alternative to the current `2e2pi`/`2mu2pi` model (see "Provenance" below).
It is **no longer called by `PlutoDecayer::decay()`** -- `PlutoLLPiPiChPT.h`
is used for both parents instead, since it is validated against real data
for both, where the EvtGen-matched model only exists for eta'. The header,
and its unit test, remain in the codebase as a cross-check/fallback
reference, not as dead code to be deleted.

## Verification status

### Build and runtime (done, passing)

- `scram b GeneratorInterface/Pythia8Interface` compiles clean (`-Werror`
  included) as of this integration.
- `cmsRun test/plutoDecayer_cfg.py` produces closed-kinematics 5-event
  samples with no exceptions for every mode implemented at the time of the
  respective test (`4e`, `4mu`, `2e2pi`, `2mu2pi`, `2mu2e` for both `eta` and
  `etaprime`, plus the `phaseSpace` control). Verified by dumping Pythia's
  event listing and confirming e.g. an eta decays exactly into
  `mu- mu+ e- e+` with correct daughter masses and closed 4-momentum.
- The full `cmsDriver.py ... --step GEN,SIM` chain (matching the actual
  production command used for `Configuration/GenProduction/python/Eta_2mu2e.py`)
  runs end to end with `Pluto:*` settings correctly propagated into
  `process.generator.PythiaParameters`; GlobalTag/era/beamspot/geometry
  resolve normally. 0/6284 raw `SoftQCD:nonDiffractive` events survived the
  downstream `etafilter`/`decayfilter`/`mufilter` acceptance cuts in that
  run, which is expected minbias-eta-pT acceptance, not a PlutoDecayer issue.
- `2mu` and `2mugamma` were added after the initial round of channel testing
  and have **not yet been individually re-run through `cmsRun`** (only
  compiled). Do this before trusting them in a production config.

### Provenance / cross-checks against independent implementations

- **EvtGen's `EvtEtaDalitz`** (read from `cms-externals/evtgen`, the exact
  fork CMSSW's pinned evtgen 2.0.0 is built from): models `eta -> pi+pi-pi0`
  (an all-hadronic Dalitz-*plot* amplitude, Layter et al., PRD 7, 2565
  (1973)). This is an unrelated process to any Pluto:mode here ("Dalitz
  plot" vs. "Dalitz decay" is a name collision, not the same physics) --
  nothing to check against it.
- **EvtGen's `EvtEta2MuMuGamma`** (single Dalitz, -> our `2mugamma`): its
  matrix element multiplies the QED tensor structure by a rho0-pole VMD
  form factor `1/((mRho^2-q^2)^2 + mRho^2*GammaRho^2)`, `q^2` = dilepton
  mass^2. `PlutoSingleDalitz.h` now includes this factor (rho0 mass/width
  taken from the caller's particle data table, PDG id 113, rather than
  EvtGen's hardcoded 0.768/0.151 GeV literals) with an algebraic rejection
  bound (the single-resonance pole has a closed-form maximum, unlike the
  eta' LLPiPi case below) -- fixed and wired in.
- **EvtGen's `EvtEtaLLPiPi`** (-> our `2e2pi`/`2mu2pi`): exists in EvtGen
  for **eta' only** -- there is no EvtGen eta version to check against.
  Its angular ("magnetic") structure, after removing the
  `s_pipi*beta_pipi^2` normalization it folds in separately, is
  algebraically the same shape as our `pipiMagneticAngular` -- a genuine
  cross-check that the independently-derived Petri-based angular formula is
  correct. Its `F0(sLL, sPiPi)` piece, however, is a coherent two-resonance
  (rho0 in both the dilepton and dipion channels, with interference) form
  factor, not constant. `interface/PlutoEtaPrimeLLPiPi.h` ports this
  formula (rho mass/width taken from the caller's particle data table, not
  hardcoded, matching EvtGen's `EvtPDL::getMeanMass/getWidth("rho0")`) with
  a grid-searched rejection envelope (the F0 rational function's true
  maximum has no simple closed form, unlike the other headers' algebraic
  bounds). Implemented, unit-tested, and confirmed working (500-event
  stress runs, no envelope violations) -- but **superseded as the active
  `2e2pi`/`2mu2pi` model** by `PlutoLLPiPiChPT.h` below once
  arXiv:2210.14925 was read, since that reference is validated against
  real data for both eta and eta', where this EvtGen port only exists for
  eta'. Kept in the codebase as a cross-check, not wired into `decay()`.

### Literature read and applied: arXiv:2210.14925 (the current 2e2pi/2mu2pi model)

**Zillinger, Kubis, Sanchez-Puertas, arXiv:2210.14925**, "CP violation in
eta(')->pi+pi-mu+mu- decays" (JHEP), Section 3 (their Standard-Model
amplitude -- the paper's main subject, a CP-violation search, is not what
we need). Read at the equation level. This is the best available reference
for `2e2pi`/`2mu2pi`: **Table 1 validates their SM amplitude against real
measured branching ratios** for both parents' electron channels
(eta->pi+pi-e+e-: predicted 2.65(17)e-4 vs. measured 2.68(9)(7)e-4;
eta'->pi+pi-e+e-: predicted 2.21(15)e-3 vs. measured 2.11(12)(15)e-3 --
both few-percent agreement), and, unlike EvtGen's `EvtEtaLLPiPi`, it covers
**eta**, which previously had no resonance-aware `2e2pi`/`2mu2pi` model at
all.

Their amplitude factorizes as `f_1(s,s_l) = P(s) * Omega_1^1(s) * Fbar(s_l)`
(`s` = pion-pair mass^2, `s_l` = dilepton mass^2):
- `P(s)`: a polynomial fit to eta/eta' -> pi+pi-gamma data (their Eq. 3.4,
  KLOE input for eta, BESIII for eta'). **Ported exactly**, including the
  eta' polynomial's isospin-violating rho-omega-mixing pole term.
- `Fbar(s_l)`: a coherent two-resonance (rho0, rho(1450)) dilepton form
  factor (their Eq. 3.5). **Ported exactly**; rho(1450) mass/width are
  read from the caller's particle data table (PDG id 100113, matching
  Pythia8's own `m0=1.465 GeV, mWidth=0.400 GeV`) the same way rho0 is,
  *not* values quoted by this paper (it doesn't state them explicitly in
  the sections read).
- `Omega_1^1(s)`: the P-wave pi-pi **Omnes function** -- a dispersive
  integral over the measured pi-pi phase shift. **Not ported**: the paper
  does not give a closed form for it; its own authors obtained their
  parametrization from unpublished work (acknowledgments: "Akdag and
  Isken"), which is not available here. **Approximated by a single
  energy-dependent-width rho0 Breit-Wigner** (reusing the exact same
  width formula already used for `Fbar`'s resonances), motivated by the
  paper's own remark that the pi-pi system is "fully dominated by the rho
  resonance at the energies of interest" -- this is an approximation to
  their rigorous dispersive treatment, not a faithful port of it.

**The angular (helicity-angle) distribution is also not re-derived from
this paper.** Their full amplitude is an epsilon-tensor structure
(`~ eps_uvab p1^v p2^a q^b`, their Eq. 2.2) whose squared matrix element
carries its own angular dependence beyond what `f_1(s,s_l)` alone
captures; correctly extracting that would be a separate, nontrivial
derivation. `PlutoLLPiPiChPT.h` reuses `pipiMagneticPointlike`'s existing,
already-cross-checked-against-EvtGen angular shape unchanged, and only
replaces what was previously a constant (F=1) with this paper's
`|P(s)|^2 |Omega(s)|^2 |Fbar(s_l)|^2` as an additional mass-dependent
reweighting. **This is a hybrid combination of two different sources for
two different physics ingredients (mass-dependence vs. angular shape), not
a first-principles port of arXiv:2210.14925's full amplitude.** Built,
unit-tested (`PlutoLLPiPiChPT` test case: weight finiteness/positivity
across the physical domain for both parents and lepton flavors, closure,
masses), and confirmed via `cmsRun` (300-event runs, `eta`/`etaprime` x
`2e2pi`/`2mu2pi`, no exceptions or envelope violations).

### Literature read and applied: arXiv:1511.04916 (the double-Dalitz TFF gap sizes)

**Escribano & Gonzalez-Solis, arXiv:1511.04916**, "A data-driven approach to
pi0, eta and eta' single and double Dalitz decays" (Chinese Phys. C 42
(2018) 023109) was read at the equation level (Sections 2-4, Tables 1-8).
It predicts both single Dalitz (`P -> l+l-gamma`) and double Dalitz
(`P -> l+l-l+l-`) decay spectra for pi0/eta/eta', e and mu, using
transition form factors (TFFs) fit to actual space-like gamma*gamma->P
data (Pade/Chisholm rational approximants), and quotes the "QED"
(constant-F, i.e. exactly what `pipiMagneticPointlike`/`mixedPointlike`/
`identicalPointlike` implement) branching ratio alongside their TFF-corrected
one in every table, which gives a direct, quantitative measure of how much
our constant-F baseline is missing:

| Channel                        | QED (=our baseline) | This paper's TFF result | Effect    |
|----------------------------------|----------------------|----------------------------|-----------|
| eta -> e+e-gamma                 | 6.38e-3              | 6.60-6.61e-3               | +3.5%     |
| eta -> mu+mu-gamma                | 2.17e-4              | 3.25-3.30e-4               | **+50%**  |
| eta' -> e+e-gamma                | 3.94e-4              | 4.35-4.42e-4               | +11%      |
| eta' -> mu+mu-gamma               | 0.38e-4              | 0.74-0.81e-4               | **+95-113%** |
| eta -> e+e-e+e- (4e)              | 2.56e-5              | 2.72-2.74e-5               | +6-7%     |
| eta -> mu+mu-mu+mu- (4mu)         | 2.59e-9              | 4.15-4.47e-9               | **+60-73%** |
| eta -> e+e-mu+mu- (2mu2e)         | 1.57e-6              | 2.35-2.39e-6               | **+50%**  |
| eta' -> e+e-e+e- (4e)             | 1.75e-6              | 2.09-2.15e-6               | +20%      |
| eta' -> mu+mu-mu+mu- (4mu)        | 0.98e-8               | 2.06-2.19e-8               | **factor ~2.1** |
| eta' -> e+e-mu+mu- (2mu2e)        | 3.21e-7               | 6.25-6.80e-7                | **factor ~2** |

The pattern is consistent and physically intuitive: the TFF correction is
small (a few percent) whenever the final state is electron-dominated,
because the accessible dilepton virtuality stays far below the resonance
region; it is large (50%-120%) whenever muons are involved, since a muon
pair's minimum invariant mass (`2m_mu`) already samples a much larger
fraction of the TFF's momentum dependence. **This directly confirms that
the constant-F=1 baseline is a poor approximation for exactly the
muon-inclusive channels this integration was built for** (`2mu2e`, `4mu`,
`2mugamma`, `2mu2pi`), while being a reasonable (few-percent) one for the
electron-only channels (`4e`, and the electron side of `2e2pi`).

For eta' specifically, Fig. 5 additionally shows the TFF has a genuine
**omega resonance peak around 0.8 GeV** in the dilepton mass spectrum, on
top of the rho contribution -- "the contribution of the rho resonance bends
the distribution, while the inclusion of the omega resonance accounts for
the sharp peak around 0.8 GeV." **`PlutoSingleDalitz.h`'s rho-only VMD fix
(matching EvtGen, which also only has rho) is therefore still incomplete
for `2mugamma` with `parent=331`** -- it is missing the omega peak. For eta,
this doesn't apply: eta's accessible dilepton mass never reaches the
rho/omega region at all (max ~0.3 GeV^2 vs. the pole at ~0.6 GeV^2), so the
rho-pole fix is a good approximation there, consistent with this paper's
own finding that simple VMD and the full data-driven Pade fit agree closely
for eta.

**Fixed as a result**: `2mu2e` (both parents) now reweights the existing
angular/phase-space sampling by a **factorized** double-virtual form
factor `F(s1,s2) = F(s1)*F(s2)` -- the paper's own "standard factorisation
ansatz" (Eq. 8) for extending a single-virtual TFF to the double-virtual
case. This is unambiguous for `2mu2e` because there is only one diagram
(no identical-lepton exchange). Initially `F` was the generic rho0-pole
model for both parents (`gen::pluto::mixedPointlikeResonant`); a later
pass (see the arXiv:1504.07742 section below) replaced `F` with the
data-fitted eta TFF specifically for eta (`mixedPointlikeEtaPade`), since
that reference is more accurate than a generic pole guess. eta' still
uses the rho0-pole version.

**Now fixed too, but only the methodology, not the underlying physics
content**: `4e`/`4mu` initially kept the constant-F `identicalPointlike`,
because the factorized reweight used for `2mu2e` isn't valid here -- the
exchange diagram pairs leptons differently than the direct diagram, so the
two diagrams need the form factor evaluated at *different* invariant-mass
pairs *before* they interfere. That bug is fixed: `PlutoIdenticalDoubleDalitz.h`'s
`Amplitudes` class now dresses each of the four photon propagators (two
per diagram) with the amplitude-level complex rho0-pole factor
(`ffAmplitude`) at its own invariant mass, then forms `|D-X|^2` exactly as
before. The existing algebraic bound `|D-X|^2 <= 2(|D|^2+|X|^2)` needed no
re-derivation, since it holds for *any* complex D and X (a direct
consequence of the AM-GM inequality), dressed or not -- so this is a sound
rejection-sampling target, not an approximation of one. `rhoMass<=0`
(the default) reproduces the original constant-F behavior exactly, checked
by a regression test comparing dressed-with-F=0 against the original
undressed call bit-for-bit. Unit-tested (bound survives, output changes
relative to constant-F, closure) and confirmed via `cmsRun` (500-event
runs, all of `eta`/`etaprime` x `4e`/`4mu`, no exceptions).

**What this does *not* fix**: even with the eta-specific Pade TFF now
plugged into `4e`/`4mu` (see below), it's still the *single-virtual* TFF
applied once per photon propagator, factorized -- *not* arXiv:1511.04916's
actual *double-virtual* TFF (their factorization ansatz or bivariate
Chisholm approximant, Eq. 8, which is a distinct fit, not merely a product
of two single-virtual ones, though the two agree at the percent level
per that paper's own comparison). For eta', the form factor dressing
`4e`/`4mu` is still the generic rho0-pole model, not eta'-specific. So:
the interference *bookkeeping* is correct for both parents; the *physics
content* dressing it is an approximation of varying quality --
data-fitted and eta-specific for eta, generic for eta'.

### Literature read and applied: arXiv:1504.07742 (the eta-specific TFF)

**Escribano, Masjuan, Sanchez-Puertas, arXiv:1504.07742**, "The eta
transition form factor from space- and time-like experimental data" (Eur.
Phys. J. C). Read at the equation level, including Appendix A, which
tabulates the actual fitted rational-function (Pade approximant)
coefficients for eta's own single-virtual TFF (their Table 5, the
`P_1^7(Q^2)` fit to `Q^2*F_etagammagamma*(Q^2)`) -- not just the
slope/curvature summary (`b_eta`, `c_eta`) that a first pass through this
paper had stopped at. `interface/PlutoEtaTFF.h` ports these 8 numbers
(`t1..t7`, `r1`) directly:

```
Q^2 F(Q^2) = t1*Q^2 + t2*Q^4 + ... + t7*Q^14   [/ (1 + r1*Q^2)]
```

with `Q^2 = -q^2` (their space-like-positive convention vs. this code's
physical time-like `q^2`). The fitted pole sits at time-like
`q^2 = +0.51 GeV^2`, beyond eta's kinematic reach (`q^2_max = M_eta^2 =
0.30 GeV^2`), so the function stays real, finite, and -- checked
numerically -- monotonically increasing across eta's entire accessible
range (up to `|F/F(0)|^2 ~ 5.7` at the kinematic edge). This is now the
form factor used for **eta** in `2mugamma`, `2mu2e`, and `4e`/`4mu`
(`singleDalitzEtaPade`, `mixedPointlikeEtaPade`, and `identicalPointlike`'s
`useEtaPade=true` branch), replacing the earlier, non-eta-specific rho0-pole
guess for those three channels. eta' is untouched by this paper (it only
covers eta) and keeps the rho0-pole model.

**A real bug was found and fixed while writing this header's unit test.**
The initial C++ port divided by `q2` instead of `Q2` when recovering
`F(Q2)` from the fitted `Q^2*F(Q^2)` quantity, an exact sign error. It was
caught because the new test reproduces the paper's own algebraic
self-check (Eq. A.2: `b_eta = (t1*r1-t2)*m_eta^2/t1 = 0.5749` from these
same coefficients) and initially failed by many orders of magnitude.
**This bug had zero effect on any physics already generated**, verified
both by the fix (all `cmsRun`/comparison output is unchanged before and
after) and algebraically: every actual use of this form factor in the
codebase squares it or multiplies two evaluations of it together
(`2mugamma`: `F^2`; `2mu2e`: `F(s1)^2*F(s2)^2`; `4e`/`4mu`: the direct and
exchange diagrams each multiply *two* evaluations of `F` before they
interfere), and an exact sign flip cancels identically in any even number
of factors. Still a real bug, now fixed and regression-tested -- a
different caller using the unsquared value directly would have gotten it
wrong.

**Cross-checked against real data (Figs. 8-11) with a digitized
comparison, saved and reproducible**: see
`test/validation/README.md` for the scripts (`digitize_figures.py`,
`compare_papers.cc`, `build_comparison.py`) that digitize
arXiv:1511.04916's own Figs. 8-11 (eta/eta' -> 2mu2e dimuon/dielectron
mass^2, eta/eta' -> 4mu one-pair mass^2 "Total distribution") directly
from the published PDF and overlay them against this code's own generated
shapes, area-normalized. Both digitization and MC generation are checked
into the repository so this comparison can be redone against any future
physics change, not just eyeballed once. Artifact:
https://claude.ai/artifact/BAgVZtbngGv1SfFZ1FiXrZ .

**arXiv:0705.0954** (Borasoy/Nissler, chiral unitary approach to
eta(')->pi+pi-l+l-) was read at the equation level too (an independent,
2007 chiral-unitary/Bethe-Salpeter calculation of the same channels
arXiv:2210.14925 covers) -- its own machinery (coupled-channel T-matrix,
numerous fitted low-energy constants) is far more involved than what's
practical to port here, but its **Table 1 gives an independent prediction**
to compare against arXiv:2210.14925's: for `eta'->pi+pi-mu+mu-`
specifically, the two calculations disagree by ~43% (2.25(14)e-5 vs.
1.57(+.40/-.47)e-5) and predict qualitatively different spectral shapes
(a single peak vs. a double-bump structure in the paper's own Fig. 6) --
real, unresolved theory-level uncertainty in the literature for that one
channel, not just an artifact of this code's own approximations. Nothing
from this paper is ported; it is used only as an independent
cross-check/sanity bound on how much uncertainty to expect from the
`2mu2pi`/`2e2pi` channels generally.

**Cross-checked against real data (Fig. 7) with a digitized comparison,
same pipeline and same caveats as the arXiv:1511.04916 comparison above**:
`test/validation/digitize_figures_0705.py` digitizes this paper's Fig. 7
(the dilepton mass spectrum, `k^2 * dGamma/d(sqrt(k^2))`, for all four of
`eta`/`eta'` x `e+e-`/`mu+mu-`) directly from the published PDF, and
`compare_papers.cc`/`build_comparison.py` overlay it against this code's
own `llpipiChPT` output (the events are histogrammed in dilepton mass,
each weighted by its own mass^2, to match the paper's k^2 weighting,
before area-normalizing). Peak positions: `eta->pi+pi-e+e-` 70 MeV (paper)
vs. 66 MeV (this code); `eta->pi+pi-mu+mu-` 217 vs. 219 MeV (both close
matches); `eta'->pi+pi-mu+mu-` 202 vs. 246 MeV; `eta'->pi+pi-e+e-` 99 vs.
131 MeV (both eta' channels show a real, non-digitization-noise shape
shift). This is directionally consistent with the ~43% eta'->mumu
normalization tension against arXiv:2210.14925 that this paper's own
Table 1 already showed above -- **the two theoretical approaches genuinely
disagree more for eta' than for eta**, and this code (built from
arXiv:2210.14925's amplitude) inherits that disagreement when compared
against arXiv:0705.0954's independent calculation. The eta channels, where
the two papers largely agree, match well here too. Same artifact as
above: https://claude.ai/artifact/BAgVZtbngGv1SfFZ1FiXrZ (second section
of the page).

### Known gaps (do not treat as production-validated yet)

| Channel                  | Parent(s)   | Status                                                                 |
|---------------------------|-------------|-------------------------------------------------------------------------|
| `2mu`                     | eta, eta'   | Trivial 2-body kinematics; no shape to get wrong. Build+run verified. |
| `2mugamma`                | eta         | **Fixed**: rho-pole VMD factor matching EvtGen's `EvtEta2MuMuGamma`; confirmed a good approximation for eta by arXiv:1511.04916 (eta never reaches the rho/omega region, q^2 up to 0.30 GeV^2 vs. the rho pole at 0.60). Build+run verified. |
| `2mugamma`                | eta'        | Rho-pole VMD factor applied, but **known incomplete**: eta''s accessible q^2 (up to 0.92 GeV^2 -- corrected from an earlier, wrong estimate of 0.56 in this doc) reaches well past the rho pole (0.60) into the region arXiv:1511.04916 says needs an omega-resonance contribution too (its ~0.8 GeV peak, i.e. M^2~0.64 GeV^2). Omega is **not included**, here or in EvtGen's model. Build+run verified (no crash), shape known-incomplete. |
| `2mu2e`                   | eta         | **Fixed**: factorized rho-pole double-virtual form factor per arXiv:1511.04916 Eq. 8; that paper shows constant-F undershoots by ~50%. Eta's dilepton reach stays below the rho pole (like `2mugamma`), so this fix should be a good approximation here. Build+run verified. |
| `2mu2e`                   | eta'        | Same fix applied, but has the **same omega gap as eta' `2mugamma`**: its dimuon pair mass reaches up toward M^2(eta')=0.92 GeV^2 (confirmed from the actual generated histogram, which extends to ~0.9 GeV^2) since the companion (electron) pair is effectively massless and imposes no cap -- this was stated too narrowly earlier in this doc as only affecting `2mugamma`. arXiv:1511.04916 shows constant-F undershoots eta' by a factor ~2 overall; the fix captures the rho part of that but not the omega part. Build+run verified. |
| `4e`, `4mu`               | eta, eta'   | **Methodology fixed, physics content still approximate**: the direct/exchange interference now correctly dresses each diagram's photon propagators with the amplitude-level rho0-pole factor at its own invariant mass (previously a real bug -- constant F=1 everywhere). But the form factor itself is still the simple single-rho-pole model, *not* arXiv:1511.04916's actual double-virtual TFF (their factorization ansatz/Chisholm approximants, whose coefficients weren't fully given in the pages read). So this addresses *where* the correction was structurally missing, not the full *size* of the originally-quantified gap (~6-20% for `4e`, ~60%-factor-2.2 for `4mu`). Unit-tested (envelope bound survives, closure, confirmed non-trivial vs. constant-F) and `cmsRun`-verified (500 events, all four parent x mode combinations). |
| `2e2pi`, `2mu2pi`         | eta, eta'   | **Replaced**: now uses `PlutoLLPiPiChPT.h`, mass-dependence from arXiv:2210.14925 (read; validated against measured BRs for both parents' e+e- channels) combined with the existing (EvtGen-cross-checked) angular shape. The pi-pi Omnes function is **approximated** by a single rho Breit-Wigner (the paper's own authors' exact parametrization isn't public); the angular part is **not** re-derived from this paper's full amplitude. `PlutoEtaPrimeLLPiPi.h` (EvtGen-matched, eta'-only) remains in the codebase, tested, but is no longer the active model. Build+run verified for all four (parent x mode) combinations, plus the full unit test suite. **Cross-checked against arXiv:0705.0954 Fig. 7** (digitized, see above): eta channels' dilepton-mass peak positions match well (e+e-: 70 vs. 66 MeV; mu+mu-: 217 vs. 219 MeV); eta' channels show a real shape difference (e+e-: 99 vs. 131 MeV; mu+mu-: 202 vs. 246 MeV), consistent with the ~43% eta'->mumu tension already known to exist between these two papers' independent calculations. |
| `phaseSpace` (all modes)  | eta, eta'   | Uniform phase space, no dynamics claimed; working as a kinematic control only. Build+run verified. |

**All four papers referenced in this document (arXiv:1511.04916,
arXiv:2210.14925, arXiv:1504.07742, arXiv:0705.0954) have been read at the
equation level and applied or cross-checked against**, not just skimmed at
the abstract -- see the sections above for what each one actually
contributed: 1511.04916 and 2210.14925 drove the `2mugamma`/`2mu2e`/
`4e`/`4mu`/`2e2pi`/`2mu2pi` model changes; 1504.07742 supplied the
eta-specific Pade TFF coefficients; 0705.0954 is an independent
cross-check (not ported) with its own digitized-figure comparison. Two
independent digitized-figure comparisons are checked into
`test/validation/` and kept reproducible against future physics changes:
arXiv:1511.04916 Figs. 8-11 (`2mu2e`/`4mu`, both parents) and
arXiv:0705.0954 Fig. 7 (`2e2pi`/`2mu2pi`, both parents).

### Summary: fixed vs. approximated vs. truly missing

Distinguishing "the bug is fixed" from "the physics is exact" matters here,
so three explicit tiers, as of the latest pass:

**Complete, nothing further owed**: `2mu` (trivial); eta `2mugamma` and eta
`2mu2e` (rho-pole fix, and arXiv:1511.04916 itself confirms simple VMD is
a good approximation for eta specifically, since eta never reaches the
resonance region).

**A real bug/methodology gap was fixed, but the physics plugged into the
fix is still an approximation, not a literal port of either paper's full
model**:
- `4e`/`4mu`: the direct/exchange interference now correctly dresses each
  diagram at the amplitude level (previously just wrong -- constant F=1
  applied uniformly) -- but the form factor itself is still the simple
  single-rho-pole model, not arXiv:1511.04916's real double-virtual TFF.
- `2e2pi`/`2mu2pi` (both parents): the pion-pair-mass polynomial and
  dilepton form factor are exact ports of arXiv:2210.14925's coefficients,
  but its pi-pi Omnes function is a rho-Breit-Wigner stand-in (the real
  one isn't public), and **the angular distribution was never touched at
  all** -- it's still the original Petri-thesis formula, never checked
  against this paper's own amplitude.

**Truly missing -- not approximated, not attempted**:
- The eta' omega resonance in `2mugamma`/`2mu2e` (confirmed absent;
  eta''s phase space reaches well past where arXiv:1511.04916 says it
  matters).
- Absolute branching ratios/normalization for all eight channels -- by
  design, every channel here is forced-exclusive shape-only.
- Any quantitative check of *this code's own predicted rates* against
  either paper's Table values -- everything validated so far is shape- or
  structure-level, never "generate N events, compute the implied BR,
  compare to their number."
- The eta' `2e2pi`/`2mu2pi` shape discrepancy against arXiv:0705.0954
  (peak shifted by ~30 MeV in e+e-, ~44 MeV in mu+mu-) is *documented*,
  not resolved -- there's no way to tell, from what's ported here, how
  much of it is this code's Omnes-function approximation vs. genuine
  theory-level disagreement between the two papers (their own Table 1
  already disagrees by ~43% on the eta'->mumu rate).
- Uncertainty propagation from either paper's quoted coefficient errors
  (e.g. `A_eta = 17.9(4)`) -- only central values are used anywhere.

## Randomness and thread safety

The pointlike models draw directly from Pythia's own random stream
(`pythia_->rndm.flat()`), which CMSSW's `RandomNumberGeneratorService`
already manages per-stream -- no extra bridging is needed. Only the
`phaseSpace` control model touches ROOT's global `gRandom` (via
`TGenPhaseSpace`); `PlutoDecayer.cc` temporarily swaps in a `TRandom`
adapter over Pythia's stream for the duration of that one call, guarded by
a mutex, since `gRandom` is a single process-wide global.
