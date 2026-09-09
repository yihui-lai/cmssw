# Identical-lepton double Dalitz

`identicalPointlike` implements eta/eta-prime -> 4e or 4mu at leading order,
including finite masses and the identical-fermion exchange interference.
The common pseudoscalar transition form factor is constant, F=1. This is a
new adapter-side model; native Pluto eta->4e is not used as an interference
reference. No radiative corrections or branching normalization are supplied.

The amplitudes follow Petri, arXiv:1010.2378, section 3.3:
https://arxiv.org/abs/1010.2378 . With output ordered l-(0), l+(1), l-(2), l+(3),
the two photon pairings are (01)(23) and (03)(21). The physical amplitude is
`D-X`, not an incoherent sum. The code evaluates the vector currents with
Dirac spinors and explicitly sums the 16 external spin combinations. In the
parent rest frame the epsilon contraction reduces to the pair three-momentum
dotted into the cross product of the spatial currents, multiplied by the
parent mass and divided by both photon virtualities.

Sampling uses an equal mixture of the direct and exchanged pointlike
densities, which have equal integrals by relabelling. It then accepts with

```
sum_spins |D-X|^2 / [2 * sum_spins (|D|^2+|X|^2)].
```

This is bounded by one algebraically. It does not estimate a maximum from
generated events. The identical-particle factorials cancel from the normalized
exclusive shape; this does not calculate a branching fraction.

Small checks (`tests/test_identical_dalitz.cc`) cover:

- Direct spin sum versus the independent angular expression, Eq. (3.36).
- Swapping either pair of identical particles and Lorentz-boost invariance.
- Amplitude-level Pauli cancellation for equal-momentum, equal-spin fermions.
- Positive coherent density, fixed rejection bound, masses and closure.

The worst direct-term relative discrepancy in the tested points was
6.61e-11. Both parents and both lepton flavours pass five-event plugin tests
with repeated and changed seeds. These are analytic/technical checks, not a
large statistical physics comparison or fitted transition-form-factor result.

For the mixed 2mu2e channel the exchange pairing is absent because the two
currents have distinguishable lepton flavours. Thus lepton universality does
not justify using the identical-channel interference correction for 2mu2e.
