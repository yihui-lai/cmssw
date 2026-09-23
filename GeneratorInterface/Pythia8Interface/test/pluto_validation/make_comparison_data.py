#!/usr/bin/env python3
"""Merge the reference curves in csv/ with PlutoDecayer's MC histograms.

Reads  : csv/*.csv        (curves digitized from the papers' figures)
         <mc.csv>         (written by ./compare_papers)
Writes : <out.json>       {"paper": {series: [[x,y],..]}, "mc": {series: [[x,y],..]}}

Each MC histogram is multiplied by ONE constant so its trapezoidal area
equals the reference curve's area: PlutoDecayer predicts shapes only (forced
exclusive decays, no branching fractions), so area-matching is the only
meaningful way to overlay them. The scale factors are printed.

Usage: make_comparison_data.py [mc.csv] [out.json]
"""
import collections
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_DIR = os.path.join(HERE, "csv")

# (series name, csv file, x column, y column, x -> GeV or GeV^2 conversion)
CHANNELS = [
    # arXiv:1511.04916 (Escribano & Gonzalez-Solis)
    ("eta_2mugamma_ee", "Fig4_eta_to_gamma_ll_thiswork_ee.csv", "sqrt_s_GeV", "1e6_x_dGamma_dsqrt_s", 1.0),
    ("eta_2mugamma_mumu", "Fig4_eta_to_gamma_ll_thiswork_mumu.csv", "sqrt_s_GeV", "1e6_x_dGamma_dsqrt_s", 1.0),
    ("etaprime_2mugamma_ee", "Fig_etaprime_to_gamma_ll_thiswork_ee.csv", "sqrt_s_GeV", "1e5_x_dGamma_dsqrt_s", 1.0),
    ("etaprime_2mugamma_mumu", "Fig_etaprime_to_gamma_ll_thiswork_mumu.csv", "sqrt_s_GeV", "1e5_x_dGamma_dsqrt_s", 1.0),
    ("eta2mu2e_ee", "Fig8_eta_to_ee_mumu_vs_Mee2.csv", "M2_ll_GeV2", "1e8_x_dGamma_dM2_ll", 1.0),
    ("eta2mu2e_mumu", "Fig8_eta_to_ee_mumu_vs_Mmumu2.csv", "M2_ll_GeV2", "1e8_x_dGamma_dM2_ll", 1.0),
    ("etaprime2mu2e_ee", "Fig10_etaprime_to_ee_mumu_vs_Mee2.csv", "M2_ll_GeV2", "1e7_x_dGamma_dM2_ll", 1.0),
    ("etaprime2mu2e_mumu", "Fig10_etaprime_to_ee_mumu_vs_Mmumu2.csv", "M2_ll_GeV2", "1e7_x_dGamma_dM2_ll", 1.0),
    ("eta4e_pair", "Fig9a_eta_to_4e_total.csv", "M2_ee_GeV2", "1e7_x_dGamma_dM2", 1.0),
    ("etaprime4e_pair", "Fig11a_etaprime_to_4e_total.csv", "M2_ee_GeV2", "1e6_x_dGamma_dM2", 1.0),
    ("eta4mu_pair", "Fig9b_eta_to_4mu_total.csv", "M2_mumu_GeV2", "1e11_x_dGamma_dM2", 1.0),
    ("etaprime4mu_pair", "Fig11b_etaprime_to_4mu_total.csv", "M2_mumu_GeV2", "1e6_x_dGamma_dM2", 1.0),
    # arXiv:0705.0954 (Borasoy & Nissler) Fig. 7: x in MeV; only the "_central"
    # (minimal-chi^2 fit) column of each error band is used.
    ("eta_pipiee", "Fig7_eta_to_pipi_ee.csv", "sqrt_k2_MeV", "k2_dGamma_dsqrtk2_keV2_central", 1e-3),
    ("eta_pipimumu", "Fig7_eta_to_pipi_mumu.csv", "sqrt_k2_MeV", "k2_dGamma_dsqrtk2_keV2_central", 1e-3),
    ("etaprime_pipiee", "Fig7_etaprime_to_pipi_ee.csv", "sqrt_k2_MeV", "k2_dGamma_dsqrtk2_keV2_central", 1e-3),
    ("etaprime_pipimumu", "Fig7_etaprime_to_pipi_mumu.csv", "sqrt_k2_MeV", "k2_dGamma_dsqrtk2_keV2_central", 1e-3),
]


def load_reference(fname, xcol, ycol, xscale):
    pts = []
    with open(os.path.join(CSV_DIR, fname)) as f:
        for row in csv.DictReader(f):
            try:
                pts.append([float(row[xcol]) * xscale, float(row[ycol])])
            except (KeyError, ValueError):
                continue
    return sorted(pts)


def area(pts):
    return sum((y0 + y1) / 2 * (x1 - x0) for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]))


def main():
    mc_path = sys.argv[1] if len(sys.argv) > 1 else "pluto_mc.csv"
    out_path = sys.argv[2] if len(sys.argv) > 2 else "comparison_data.json"
    mc = collections.defaultdict(list)
    with open(mc_path) as f:
        for row in csv.DictReader(f):
            mc[row["series"]].append([float(row["x"]), float(row["count"])])

    out = {"paper": {}, "mc": {}}
    for name, fname, xcol, ycol, xscale in CHANNELS:
        ref = load_reference(fname, xcol, ycol, xscale)
        hist = sorted(mc[name])
        scale = area(ref) / area(hist) if area(hist) > 0 else 1.0
        out["paper"][name] = ref
        out["mc"][name] = [[x, c * scale] for x, c in hist]
        print(f"{name:26s} ref points={len(ref):5d}  MC scale factor={scale:.4g}")
    with open(out_path, "w") as f:
        json.dump(out, f)
    print("wrote", out_path)


if __name__ == "__main__":
    main()
