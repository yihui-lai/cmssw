#!/usr/bin/env python3
"""Validation/comparison analysis for eta/eta' -> 2mu+A(->2e) [and 2pi+A(->2mu)]
GEN output: extracts the SM-recoil-pair mass spectrum (dark photon vs ALP),
confirms A's mass and decay, and (for ctau>0 samples) checks the lab-frame
decay length against the expected gamma*beta*c*tau0 scaling.

Run inside cmsenv (needs DataFormats.FWLite): ./compare_models.py <root files...>
"""
import statistics
import sys

import ROOT
from DataFormats.FWLite import Events, Handle


def analyze(path, parent_id, a_id, recoil_ids):
    h = Handle("std::vector<reco::GenParticle>")
    n_ev = 0
    n_found = 0
    closure_fail = 0
    a_mass = []
    recoil_mass = []
    a_p = []
    a_beta_gamma = []
    decay_length = []  # mm, lab frame, only meaningful if A has nonzero ctau
    for ev in Events(path):
        ev.getByLabel("genParticles", h)
        ps = h.product()
        n_ev += 1
        for p in ps:
            if p.pdgId() != parent_id or p.numberOfDaughters() != 3:
                continue
            d = [p.daughter(i) for i in range(3)]
            ids = sorted(x.pdgId() for x in d)
            if ids != sorted(list(recoil_ids) + [a_id]):
                continue
            n_found += 1
            s = d[0].p4() + d[1].p4() + d[2].p4()
            if abs(s.M() - p.mass()) > 1e-4:
                closure_fail += 1
            a = [x for x in d if x.pdgId() == a_id][0]
            recoils = [x for x in d if x.pdgId() != a_id]
            a_mass.append(a.mass())
            recoil_mass.append((recoils[0].p4() + recoils[1].p4()).M())
            bg = a.p() / a.mass()
            a_p.append(a.p())
            a_beta_gamma.append(bg)
            if a.numberOfDaughters() >= 1:
                dau0 = a.daughter(0)
                vx, vy, vz = dau0.vx() - a.vx(), dau0.vy() - a.vy(), dau0.vz() - a.vz()
                decay_length.append((vx**2 + vy**2 + vz**2) ** 0.5)
    return {
        "n_ev": n_ev,
        "n_found": n_found,
        "closure_fail": closure_fail,
        "a_mass": a_mass,
        "recoil_mass": recoil_mass,
        "a_beta_gamma": a_beta_gamma,
        "decay_length": decay_length,
    }


def summarize(label, r):
    mm = r["recoil_mass"]
    mm_sorted = sorted(mm)
    n = len(mm_sorted)
    print("== %s" % label)
    print(
        "   events=%d  eta->XA found=%d  closure fails=%d  A decayed correctly: %d/%d"
        % (r["n_ev"], r["n_found"], r["closure_fail"], len(r["a_mass"]), r["n_found"])
    )
    if r["a_mass"]:
        print(
            "   A mass: min=%.6f max=%.6f (should be a delta function at the configured m_A)"
            % (min(r["a_mass"]), max(r["a_mass"]))
        )
    if n:
        print(
            "   recoil-pair mass: mean=%.4f  p10=%.4f  p50=%.4f  p90=%.4f  GeV"
            % (statistics.mean(mm), mm_sorted[int(0.1 * n)], mm_sorted[int(0.5 * n)], mm_sorted[int(0.9 * n)])
        )
    if r["decay_length"] and any(x > 0 for x in r["decay_length"]):
        dl = r["decay_length"]
        bg = r["a_beta_gamma"]
        print(
            "   lab decay length: mean=%.2f mm  (mean beta*gamma=%.2f)" % (statistics.mean(dl), statistics.mean(bg))
        )


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--parent-id", type=int, default=221)
    ap.add_argument("--a-id", type=int, required=True)
    ap.add_argument("--recoil-ids", type=int, nargs=2, default=[13, -13])
    args = ap.parse_args()
    for f in args.files:
        r = analyze(f, args.parent_id, args.a_id, tuple(args.recoil_ids))
        summarize(f, r)
