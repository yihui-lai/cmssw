"""Forced finite-mass, constant-form-factor eta/eta-prime decay smoke test.

Exercises PlutoDecayer (GeneratorInterface/Pythia8Interface/interface/PlutoDecayer.h),
registered with Pythia8 through the classic DecayHandler/setDecayPtr mechanism
(works against whatever Pythia8 is pinned for this CMSSW release -- no
Init:plugins, no external Pluto library).
"""
from pathlib import Path
import runpy
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.maxEvents = 5
options.outputFile = "pluto-decayer.root"
options.register("parent", "eta", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "eta or etaprime")
options.register("mode", "2mu2e", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "4e, 2mu2e, 4mu, 2e2pi or 2mu2pi")
options.register("model", "pointlike", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "pointlike or phaseSpace")
options.parseArguments()

parent_id = {"eta": 221, "etaprime": 331}[options.parent]
commands = [
    "Pluto:filter = on",
    "Pluto:allowForcedDecay = on",
    f"Pluto:parent = {parent_id}",
    f"Pluto:mode = {options.mode}",
    f"Pluto:model = {options.model}",
    "HardQCD:all = on", "PhaseSpace:pTHatMin = 50", "PartonLevel:MPI = off",
    "TimeShower:QEDshowerByL = off",
]

test_dir = Path(__file__).resolve().parent
make_process = runpy.run_path(str(test_dir / "pluginSmokeCommon.py"))["make_process"]
process = make_process(commands, options.maxEvents, options.outputFile)
