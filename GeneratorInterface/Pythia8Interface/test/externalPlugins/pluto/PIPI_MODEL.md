# Pion-pair dilepton model

`pipiMagneticPointlike` supports eta and eta-prime, with either electron or
muon pairs and charged pions. It is adapter-side code, not a claim that public
Pluto supplies a correct ready-made muon model.

The reference is the magnetic term of Petri, arXiv:1010.2378, Eq. (3.85):
https://arxiv.org/abs/1010.2378 . Electric/CP-violating terms are omitted, and
the magnetic form factor is constant. This is a named leading-order baseline,
not a fitted rho/omega or pion-rescattering model. It supplies decay shapes,
not branching fractions or a physical channel mixture.

With `s = m_pipi^2`, `q = m_ll^2`, and beta the pair-rest-frame daughter speed,
the sampling density, up to constants independent of the decay variables, is

```
dGamma/(ds dq dcos(theta_pi) dcos(theta_l) dphi)
  ~ lambda(1,s/M^2,q/M^2)^(3/2) * s * beta_pi^3 * beta_l / q
    * sin(theta_pi)^2 * [1-beta_l^2*sin(theta_l)^2*sin(phi)^2].
```

Uniform s and logarithmic q proposals remove the photon-pole factor. Each
remaining factor has a fixed bound over the allowed mass domain. The code
uses their product as a rejection envelope, including the narrow eta ->
mumu pi pi phase space; it does not increase its maximum after seeing events.
The overall decay orientation is isotropic. Both lepton and pion masses are
taken from Pythia by the adapter.

## Why the native electron model remains disabled

The fetched native Pluto source has multiple independent issues:

- `PChannelModel`'s five-argument wrapper creates `m[5]`; the pion angular
  routine reads `mass[5]`, beyond that array, and assigns the other angle
  from the wrong slot.
- It constructs the lepton-angle particle from `pip`, not `ep`.
- Its rejection maximum grows after encountering larger weights.
- Its magnetic angular bracket has a plus sign where Eq. (3.85) has minus.
- `PKinematics::lambda` expects masses and internally squares them. The
  angular routine instead passes squared masses, changing the mass dependence.

The new named model avoids these paths. It does not silently modify or enable
the old native model. Independent form-factor and normalization choices are
still needed for precision physics, particularly for eta-prime.
