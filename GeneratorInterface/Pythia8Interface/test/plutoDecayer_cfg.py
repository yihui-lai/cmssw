"""Forced finite-mass, constant-form-factor eta/eta-prime decay smoke test."""
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
options.parseArguments()
test_dir = Path(__file__).resolve().parent
catalogue = runpy.run_path(str(test_dir / "externalPlugins/pluto/channels.py"))
commands = catalogue["pointlike_commands"](options.parent, options.mode)
commands += ["HardQCD:all = on", "PhaseSpace:pTHatMin = 50", "PartonLevel:MPI = off",
             "TimeShower:QEDshowerByL = off"]
make_process = runpy.run_path(str(test_dir / "pluginSmokeCommon.py"))["make_process"]
process = make_process(commands, options.maxEvents, options.outputFile)
