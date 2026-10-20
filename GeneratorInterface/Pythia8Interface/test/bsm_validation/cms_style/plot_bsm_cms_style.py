#!/usr/bin/env python3
"""CMS-style (mplhep) validation plots for the BSM eta/eta' -> l+l-A / pi+pi-A
channels, comparing EtaDecayer's MC sampler output against a numerical
marginalization of its OWN weight function (see bsm_cms_gen.cc) -- there is
no external published curve for these private BSM channels, so this is a
sampler self-consistency check, not a comparison to data or theory papers.

Reads bsm_mc.csv (written by bsm_cms_gen), writes <series>.png/.pdf here.
Run with the HOST system python3 (has mplhep); NOT inside the cms apptainer.
"""
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as mh

mh.style.use("CMS")

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "bsm_mc.csv")

XLABEL = {
    "mll": r"$M(\ell^{+}\ell^{-})$ [GeV]",
    "mpipi": r"$M(\pi^{+}\pi^{-})$ [GeV]",
    "costheta": r"$\cos\theta^{*}$",
}

CHANNEL_TITLE = {
    "eta_llA_darkphoton": r"$\eta \to \ell^{+}\ell^{-}A$, $A$=dark photon",
    "eta_llA_alp": r"$\eta \to \ell^{+}\ell^{-}A$, $A$=ALP",
    "eta_pipiA_darkphoton": r"$\eta \to \pi^{+}\pi^{-}A$, $A$=dark photon",
    "eta_pipiA_alp": r"$\eta \to \pi^{+}\pi^{-}A$, $A$=ALP",
    "etaprime_llA_darkphoton": r"$\eta' \to \ell^{+}\ell^{-}A$, $A$=dark photon",
    "etaprime_llA_alp": r"$\eta' \to \ell^{+}\ell^{-}A$, $A$=ALP",
    "etaprime_pipiA_darkphoton": r"$\eta' \to \pi^{+}\pi^{-}A$, $A$=dark photon",
    "etaprime_pipiA_alp": r"$\eta' \to \pi^{+}\pi^{-}A$, $A$=ALP",
}


def load():
    data = defaultdict(lambda: defaultdict(list))
    with open(CSV) as fh:
        for row in csv.DictReader(fh):
            data[row["series"]][row["kind"]].append((float(row["x"]), float(row["value"])))
    for series in data:
        for kind in data[series]:
            data[series][kind].sort()
    return data


def observable_of(series):
    for obs in ("costheta", "mll", "mpipi"):
        if series.endswith("_" + obs):
            return obs
    raise ValueError(series)


def channel_of(series):
    obs = observable_of(series)
    return series[: -(len(obs) + 1)]


def plot_one(series, mc, theory, outdir):
    obs = observable_of(series)
    channel = channel_of(series)
    xs_mc = [p[0] for p in mc]
    ys_mc = [p[1] for p in mc]
    n_mc = sum(ys_mc)
    width = xs_mc[1] - xs_mc[0] if len(xs_mc) > 1 else 1.0

    xs_th = [p[0] for p in theory]
    ys_th = [p[1] for p in theory]
    # theory curve already normalized to unit area in dumpCurve's own binning;
    # rescale both to the same area-under-histogram convention (bin width x counts).
    n_th = sum(ys_th)
    scale_mc_to_density = 1.0 / (n_mc * width) if n_mc > 0 else 1.0
    scale_th_to_density = 1.0 / (n_th * (xs_th[1] - xs_th[0])) if n_th > 0 else 1.0

    y_mc_density = [y * scale_mc_to_density for y in ys_mc]
    yerr_mc_density = [((y ** 0.5) if y > 0 else 1.0) * scale_mc_to_density for y in ys_mc]
    y_th_density = [y * scale_th_to_density for y in ys_th]

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.errorbar(xs_mc, y_mc_density, yerr=yerr_mc_density, fmt="o", color="black",
                markersize=4, label="MC (EtaDecayer sampler)", zorder=3)
    ax.plot(xs_th, y_th_density, "-", color="#CC3311", linewidth=2.5,
            label="Theory (own weight, numerically marginalized)", zorder=2)
    ax.set_xlabel(XLABEL[obs])
    ax.set_ylabel("Normalized density")
    ymax = max(max(y_mc_density, default=0), max(y_th_density, default=0))
    ax.set_ylim(bottom=0, top=ymax * 1.3)
    ax.legend(fontsize="x-small", loc="lower right")
    ax.text(0.04, 0.96, CHANNEL_TITLE[channel], transform=ax.transAxes, ha="left", va="top",
            fontsize="small")
    mh.label.exp_label(ax=ax, exp="CMS", data=True, llabel="Simulation Private Work", rlabel="")

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{series}.{ext}"), bbox_inches="tight", dpi=200, transparent=True)
    plt.close(fig)


def main():
    data = load()
    outdir = HERE
    made = []
    for series in sorted(data):
        mc = data[series].get("mc", [])
        theory = data[series].get("theory", [])
        if not mc or not theory:
            print(f"SKIP {series}: missing mc or theory rows")
            continue
        plot_one(series, mc, theory, outdir)
        made.append(series)
    print(f"wrote {len(made)} plots (x2 for png+pdf) to {outdir}")


if __name__ == "__main__":
    main()
