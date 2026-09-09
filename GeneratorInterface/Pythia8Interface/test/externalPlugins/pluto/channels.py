"""Explicit forced-channel catalogue; branching ratios are not supplied here.

Native model availability is a source-level fact, not physics validation.
The phaseSpace controls use ROOT phase space and do not model rare decays.
"""

PARENTS = {"eta": 221, "etaprime": 331}
DAUGHTERS = {
    "4e": (11, -11, 11, -11),
    "2mu2e": (13, -13, 11, -11),
    "4mu": (13, -13, 13, -13),
    "2e2pi": (11, -11, 211, -211),
    "2mu2pi": (13, -13, 211, -211),
}

POINTLIKE_MODELS = {
    "4e": "identicalPointlike",
    "2mu2e": "mixedPointlike",
    "4mu": "identicalPointlike",
    "2e2pi": "pipiMagneticPointlike",
    "2mu2pi": "pipiMagneticPointlike",
}

def pointlike_commands(parent, mode):
    """Explicitly choose the finite-mass LO constant-form-factor baseline.

    Includes identical-lepton exchange for 4e/4mu and magnetic-only pion modes.
    This name makes the physics assumption explicit; it does not silently
    replace native VMD models or assign branching fractions.
    """
    return commands(parent, mode, model=POINTLIKE_MODELS[mode])

def commands(parent, mode, model="native"):
    parent_id = PARENTS[parent]
    if mode not in DAUGHTERS:
        raise ValueError(f"Unknown decay mode: {mode}")
    if model not in ("native", "phaseSpace", "mixedPointlike", "pipiMagneticPointlike", "identicalPointlike"):
        raise ValueError(f"Unknown model: {model}")
    if model == "identicalPointlike" and mode not in ("4e", "4mu"):
        raise ValueError("identicalPointlike requires two identical lepton pairs")
    if model == "mixedPointlike" and mode != "2mu2e":
        raise ValueError("mixedPointlike only implements the distinguishable 2mu2e channel")
    if model == "pipiMagneticPointlike" and mode not in ("2e2pi", "2mu2pi"):
        raise ValueError("pipiMagneticPointlike requires two leptons and two charged pions")
    if model == "native" and (parent != "eta" or mode != "4e"):
        raise ValueError("A native rare-decay model still needs implementation/validation for this channel")
    return [
        "Init:plugins = {libCMSPluto.so::CMSPluto}",
        f"Pluto:parent = {parent_id}",
        f"Pluto:mode = {mode}",
        f"Pluto:model = {model}",
        "Pluto:allowForcedDecay = on",
    ]
