#!/bin/bash
# Runs etaDecayer_cfg.py, forcing every eta and then one eta per event, and checks in both cases that
# eta -> mu- mu+ e- e+ decays were made by EtaDecayer (with Pythia's own eta decay table this final state,
# BR ~ 1e-6, would essentially never appear in 10 minimum-bias events).
set -o pipefail
function die { echo "$1: status $2"; exit $2; }
for onePerEvent in False True; do
  log=etaDecayer_onePerEvent${onePerEvent}.log
  cmsRun ${SCRAM_TEST_PATH}/etaDecayer_cfg.py onePerEvent=${onePerEvent} 2>&1 | tee ${log} || die "cmsRun etaDecayer_cfg.py onePerEvent=${onePerEvent}" $?
  grep -Eq "TrigReport Events total = 10 passed = [1-9][0-9]* " ${log} || die "no eta -> mu mu e e decay found (onePerEvent=${onePerEvent})" 1
done
