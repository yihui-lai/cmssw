#!/usr/bin/env python3
"""Render comparison_data.json (from make_comparison_data.py) as PNG overlays.

One PNG per paper: paper_compare_1511.png (12 panels) and
paper_compare_0705.png (4 panels). Reference curve = line, PlutoDecayer MC
(area-scaled to the reference) = dots. Only needs PIL (no matplotlib).

Usage: plot_comparison_png.py [comparison_data.json] [output_dir]
"""
import json
import math
import os
import sys
from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))


def _find_font(candidates):
    for c in candidates:
        if os.path.exists(c):
            return c
    return None


_FONTS = {
    "reg": ["/usr/share/fonts/google-noto/NotoSans-Regular.ttf", "/usr/share/fonts/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf"],
    "bold": ["/usr/share/fonts/google-noto/NotoSans-Bold.ttf", "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
             "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf"],
    "mono": ["/usr/share/fonts/google-noto/NotoSansMono-Regular.ttf", "/usr/share/fonts/dejavu/DejaVuSansMono.ttf",
             "/usr/share/fonts/liberation-mono/LiberationMono-Regular.ttf"],
}
REG, BOLD, MONO = "reg", "bold", "mono"


def F(kind, size):
    path = _find_font(_FONTS[kind])
    return ImageFont.truetype(path, size) if path else ImageFont.load_default()


P1511 = [
    ("η → e⁺e⁻γ  dilepton mass", "eta_2mugamma_ee", "M(e⁺e⁻) [GeV]", "log", 0),
    ("η → μ⁺μ⁻γ  dilepton mass", "eta_2mugamma_mumu", "M(μ⁺μ⁻) [GeV]", "log", 0),
    ("η′ → e⁺e⁻γ  dilepton mass", "etaprime_2mugamma_ee", "M(e⁺e⁻) [GeV]", "log", 1),
    ("η′ → μ⁺μ⁻γ  dilepton mass", "etaprime_2mugamma_mumu", "M(μ⁺μ⁻) [GeV]", "log", 1),
    ("η → μ⁺μ⁻e⁺e⁻  dimuon M²", "eta2mu2e_mumu", "M²(μ⁺μ⁻) [GeV²]", "log", 0),
    ("η → μ⁺μ⁻e⁺e⁻  dielectron M²", "eta2mu2e_ee", "M²(e⁺e⁻) [GeV²]", "log", 0),
    ("η′ → μ⁺μ⁻e⁺e⁻  dimuon M²", "etaprime2mu2e_mumu", "M²(μ⁺μ⁻) [GeV²]", "log", 1),
    ("η′ → μ⁺μ⁻e⁺e⁻  dielectron M²", "etaprime2mu2e_ee", "M²(e⁺e⁻) [GeV²]", "log", 1),
    ("η → 4e  pair M²", "eta4e_pair", "M²(e⁺e⁻) [GeV²]", "log", 0),
    ("η′ → 4e  pair M²", "etaprime4e_pair", "M²(e⁺e⁻) [GeV²]", "log", 1),
    ("η → 4μ  pair M²", "eta4mu_pair", "M²(μ⁺μ⁻) [GeV²]", "log", 0),
    ("η′ → 4μ  pair M²", "etaprime4mu_pair", "M²(μ⁺μ⁻) [GeV²]", "log", 1),
]
P0705 = [
    ("η → π⁺π⁻e⁺e⁻  dilepton mass", "eta_pipiee", "M(e⁺e⁻) [GeV]", "linear", 0),
    ("η → π⁺π⁻μ⁺μ⁻  dilepton mass", "eta_pipimumu", "M(μ⁺μ⁻) [GeV]", "linear", 0),
    ("η′ → π⁺π⁻e⁺e⁻  dilepton mass", "etaprime_pipiee", "M(e⁺e⁻) [GeV]", "linear", 1),
    ("η′ → π⁺π⁻μ⁺μ⁻  dilepton mass", "etaprime_pipimumu", "M(μ⁺μ⁻) [GeV]", "linear", 1),
]
REF_COL = "#1f1f1f"
MC_COL = ["#d1495b", "#00798c"]  # eta, eta'
NOTES = {
    "etaprime_2mugamma_mumu": "Padé + ρ/ω/φ VMD TFF",
    "etaprime_2mugamma_ee": "Padé + ρ/ω/φ VMD TFF",
    "etaprime2mu2e_mumu": "Padé + ρ/ω/φ VMD TFF",
    "etaprime2mu2e_ee": "Padé + ρ/ω/φ VMD TFF",
    "etaprime4e_pair": "ρ-pole model",
    "etaprime4mu_pair": "ρ-pole model",
}


def nice_ticks(lo, hi, n=5):
    span = hi - lo
    step0 = span / n
    mag = 10 ** math.floor(math.log10(step0))
    step = min((s * mag for s in (1, 2, 2.5, 5, 10) if s * mag >= step0), default=step0)
    t = math.ceil(lo / step) * step
    out = []
    while t <= hi + 1e-12:
        out.append(round(t, 12))
        t += step
    return out


def panel(img, box, spec, data, sub):
    title, key, xlabel, scale, ci = spec
    x0, y0, x1, y1 = box
    d = ImageDraw.Draw(img)
    d.text((x0, y0), title.replace("→", "->"), font=F(BOLD, 15), fill="#111")
    note = NOTES.get(key, "")
    if note:
        d.text((x1, y0 + 2), note, font=F(REG, 12), fill="#b0303f" if note.startswith("Padé") else "#777", anchor="ra")
    ml, mr, mt, mb = 62, 10, 26, 40
    px0, py0, px1, py1 = x0 + ml, y0 + mt, x1 - mr, y1 - mb
    ref = data["paper"].get(key, [])
    mc = data["mc"].get(key, [])
    pts = [p for p in ref if p[1] > 0] + [p for p in mc if p[1] > 0]
    xmax = max(p[0] for p in ref)
    if scale == "linear":
        xmin = 0.0
    else:
        xmin = min(min(p[0] for p in mc), min(p[0] for p in ref if p[1] > 0))
    ys = [p[1] for p in pts if xmin <= p[0] <= xmax]
    if scale == "log":
        ymax = max(ys)
        ymin = max(min(ys), ymax * 1e-6)
        l0, l1 = math.floor(math.log10(ymin)), math.ceil(math.log10(ymax))
        fy = lambda v: (math.log10(max(v, 10.0**l0)) - l0) / (l1 - l0)
    else:
        l0, l1 = 0.0, max(ys) * 1.08
        fy = lambda v: v / l1
    fx = lambda v: (v - xmin) / (xmax - xmin)
    X = lambda v: px0 + fx(v) * (px1 - px0)
    Y = lambda v: py1 - fy(v) * (py1 - py0)
    d.rectangle([px0, py0, px1, py1], outline="#888")
    for t in nice_ticks(xmin, xmax):
        if xmin <= t <= xmax:
            d.line([X(t), py1, X(t), py1 + 4], fill="#888")
            d.text((X(t), py1 + 6), f"{t:g}", font=F(MONO, 11), fill="#444", anchor="ma")
    if scale == "log":
        for e in range(int(l0), int(l1) + 1):
            yy = Y(10.0**e)
            d.line([px0 - 4, yy, px0, yy], fill="#888")
            d.line([px0, yy, px1, yy], fill="#e4e4e4")
            d.text((px0 - 6, yy), f"1e{e}", font=F(MONO, 11), fill="#444", anchor="rm")
    else:
        for t in nice_ticks(0, l1):
            yy = Y(t)
            d.line([px0 - 4, yy, px0, yy], fill="#888")
            d.line([px0, yy, px1, yy], fill="#e4e4e4")
            d.text((px0 - 6, yy), f"{t:g}", font=F(MONO, 11), fill="#444", anchor="rm")
    d.text(((px0 + px1) / 2, y1 - 4), xlabel, font=F(REG, 12), fill="#333", anchor="ms")
    # reference line
    line = [(X(x), Y(y)) for x, y in sorted(ref) if xmin <= x <= xmax and (scale != "log" or y > 0)]
    if len(line) > 1:
        d.line(line, fill=REF_COL, width=2)
    # MC dots
    for x, y in mc:
        if xmin <= x <= xmax and (scale != "log" or y > 0):
            cx, cy = X(x), Y(y)
            d.ellipse([cx - 3.2, cy - 3.2, cx + 3.2, cy + 3.2], fill=MC_COL[ci], outline="white")


def render(panels, data, out, ncol, title):
    pw, ph = 470, 300
    nrow = math.ceil(len(panels) / ncol)
    W, H = ncol * pw + 30, nrow * ph + 90
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    d.text((15, 12), title, font=F(BOLD, 20), fill="#111")
    d.line([15, 52, 45, 52], fill=REF_COL, width=3)
    d.text((52, 44), "reference (digitized, csv/)", font=F(REG, 13), fill="#222")
    d.ellipse([280, 46, 288, 54], fill=MC_COL[0])
    d.ellipse([294, 46, 302, 54], fill=MC_COL[1])
    d.text((310, 44), "PlutoDecayer MC (η / η′), area-scaled to reference", font=F(REG, 13), fill="#222")
    for i, spec in enumerate(panels):
        r, c = divmod(i, ncol)
        box = (15 + c * pw, 75 + r * ph, 15 + (c + 1) * pw - 12, 75 + (r + 1) * ph - 6)
        panel(img, box, spec, data, i)
    img.save(out)
    print("wrote", out)


if __name__ == "__main__":
    data_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "output", "comparison_data.json")
    outdir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(data_path))
    data = json.load(open(data_path))
    render(P1511, data, os.path.join(outdir, "paper_compare_1511.png"), 3, "PlutoDecayer vs arXiv:1511.04916 (Escribano & González-Solís)")
    render(P0705, data, os.path.join(outdir, "paper_compare_0705.png"), 2, "PlutoDecayer vs arXiv:0705.0954 (Borasoy & Nissler)")
