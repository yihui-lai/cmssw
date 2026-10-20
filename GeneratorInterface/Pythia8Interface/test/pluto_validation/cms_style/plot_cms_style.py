#!/usr/bin/env python3
"""CMS-style (mplhep) generator-level validation plots: MC (real EtaDecayer
accept-reject samplers) vs. theory (each sampler's own documented physical
density, evaluated independently via flat Monte Carlo integration -- see
gen_sm_validation.cc's module docstring for the exact method per channel).

Run with the HOST system python3 (mplhep lives there, not in the cmsenv
apptainer used to build/run gen_sm_validation):
    python3 plot_cms_style.py sm_validation.csv
"""
import sys
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as mh

mh.style.use("CMS")

NBINS = 40

GREEK = {"eta": r"$\eta$", "etaprime": r"$\eta^{\prime}$"}


def nice_label(tag):
    parent, rest = tag.split("_", 1)
    p = GREEK[parent]
    chan = {
        "2mu": r"$\mu^+\mu^-$",
        "2mugamma": r"$\mu^+\mu^-\gamma$",
        "2egamma": r"$e^+e^-\gamma$",
        "2mu2e": r"$\mu^+\mu^-e^+e^-$",
        "4mu": r"$\mu^+\mu^-\mu^+\mu^-$",
        "4e": r"$e^+e^-e^+e^-$",
        "2epi": r"$e^+e^-\pi^+\pi^-$",
        "2mupi": r"$\mu^+\mu^-\pi^+\pi^-$",
    }.get(rest, rest)
    return f"{p} → {chan}"


def area_normalize(mc_vals, mc_counts, th_x, th_w):
    mc_area = np.sum(mc_counts)
    # theory histogram over the same binning range for area matching
    th_counts, _ = np.histogram(th_x, bins=len(mc_counts), range=(mc_vals[0], mc_vals[-1]), weights=th_w)
    th_area = np.sum(th_counts)
    if th_area <= 0:
        return th_counts
    return th_counts * (mc_area / th_area)


def make_plot(df, tag, obs, xlabel, out_prefix, xrange=None):
    mc = df[(df.series == f"{tag}_mc_{obs}")]["x"].to_numpy()
    th = df[(df.series == f"{tag}_th_{obs}")]
    if len(mc) == 0 or len(th) == 0:
        print(f"  [skip] {tag}_{obs}: no data (mc={len(mc)}, th={len(th)})")
        return False
    th_x = th["x"].to_numpy()
    th_w = th["w"].to_numpy()

    lo, hi = (xrange if xrange is not None else (mc.min(), mc.max()))
    edges = np.linspace(lo, hi, NBINS + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mc_counts, _ = np.histogram(mc, bins=edges)
    mc_err = np.sqrt(np.maximum(mc_counts, 1))

    th_counts, _ = np.histogram(th_x, bins=edges, weights=th_w)
    th_area = np.sum(th_counts)
    mc_area = np.sum(mc_counts)
    if th_area > 0:
        th_counts_scaled = th_counts * (mc_area / th_area)
    else:
        print(f"  [skip] {tag}_{obs}: theory curve has zero area in range")
        return False

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.errorbar(centers, mc_counts, yerr=mc_err, fmt="o", color="black", markersize=4, label="MC (EtaDecayer)")
    ax.plot(centers, th_counts_scaled, "-", color="#CC0000", linewidth=2, label="Theory (own weight, MC-integrated)")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Events / bin (area-normalized)")
    ax.legend(fontsize="x-small", loc="upper right", bbox_to_anchor=(1.0, 0.88))
    mh.label.exp_label(ax=ax, exp="CMS", data=True, llabel="Simulation Private Work", rlabel="")
    fig.text(0.95, 0.98, nice_label(tag), ha="right", va="top", fontsize="small", transform=ax.transAxes)
    fig.savefig(f"{out_prefix}.png", bbox_inches="tight", dpi=200, transparent=True)
    fig.savefig(f"{out_prefix}.pdf", bbox_inches="tight", transparent=True)
    plt.close(fig)

    # quick numeric verdict: chi2-like residual in populated bins
    mask = mc_counts > 0
    resid = np.abs(mc_counts[mask] - th_counts_scaled[mask]) / mc_err[mask]
    verdict = "OK" if np.median(resid) < 3 else "CHECK"
    print(f"  [{verdict}] {out_prefix}  median|residual|/err={np.median(resid):.2f}  max={np.max(resid):.2f}")
    return True


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "sm_validation.csv"
    print(f"loading {path} ...")
    df = pd.read_csv(path, dtype={"series": "category", "x": "float32", "w": "float32"})
    print(f"loaded {len(df)} rows, {df.series.nunique()} series")

    plans = []
    # 2mu trivial
    for tag in ("eta_2mu", "etaprime_2mu"):
        plans.append((tag, "costheta", r"$\cos\theta^{*}(\mu^-)$", None))
    # single Dalitz
    for tag in ("eta_2mugamma", "eta_2egamma", "etaprime_2mugamma", "etaprime_2egamma"):
        plans.append((tag, "mass", r"$M(\ell^+\ell^-)$ [GeV]", None))
        plans.append((tag, "costheta", r"$\cos\theta^{*}(\ell^-)$", None))
    # mixed double Dalitz
    for tag in ("eta_2mu2e", "etaprime_2mu2e"):
        for sub in ("mumu", "ee"):
            plans.append((f"{tag}_{sub}", "mass", rf"$M$ [GeV] ({sub} pair, {tag})", None))
            plans.append((f"{tag}_{sub}", "costheta", rf"$\cos\theta^{{*}}$ ({sub} pair, {tag})", None))
    # identical double Dalitz
    for tag in ("eta_4mu", "eta_4e", "etaprime_4mu", "etaprime_4e"):
        plans.append((tag, "mass", r"$M(\ell^+\ell^-)$ pair [GeV]", None))
        plans.append((tag, "costheta", r"$\cos\theta^{*}(\ell^-)$ pair", None))
    # ChPT 2pi
    for tag in ("eta_2epi", "eta_2mupi", "etaprime_2epi", "etaprime_2mupi"):
        plans.append((tag, "mll", r"$M(\ell^+\ell^-)$ [GeV]", None))
        plans.append((tag, "mpipi", r"$M(\pi^+\pi^-)$ [GeV]", None))
        plans.append((tag, "costhetal", r"$\cos\theta^{*}(\ell^-)$", None))
        plans.append((tag, "costhetapi", r"$\cos\theta^{*}(\pi^-)$", None))

    n_ok = 0
    for tag, obs, xlabel, xr in plans:
        label_tag = tag
        if "_mumu" in tag or "_ee" in tag:
            label_tag = tag.rsplit("_", 1)[0]
        print(f"plotting {tag}_{obs} ...")
        ok = make_plot(df, tag, obs, xlabel, f"{tag}_{obs}", xrange=xr)
        n_ok += int(ok)
    print(f"\n{n_ok}/{len(plans)} plots produced")


if __name__ == "__main__":
    main()
