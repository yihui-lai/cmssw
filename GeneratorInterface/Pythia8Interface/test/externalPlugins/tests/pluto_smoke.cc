#include "Pythia8/Pythia.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <dlfcn.h>

int main(int argc, char** argv) {
  std::cout << std::unitbuf;
  if (argc < 4) return 2;
  const int parentId = std::stoi(argv[1]);
  const std::string mode = argv[2];
  // cmsRun already loads ROOT. Reproduce that environment in this standalone
  // executable, and clear ROOT's optional-symbol lookup error before Pythia's
  // loader checks dlerror(). A non-null dlopen handle is required here.
  if (!dlopen("libCore.so", RTLD_NOW | RTLD_GLOBAL)) throw std::runtime_error(dlerror());
  dlerror();
  Pythia8::Pythia pythia;
  auto command = [&](std::string value) { if (!pythia.readString(value)) throw std::runtime_error(value); };
  command("Init:plugins = {libCMSPluto.so::CMSPluto}");
  command("Pluto:parent = " + std::string(argv[1]));
  command("Pluto:mode = " + mode);
  command("Pluto:model = " + std::string(argv[3]));
  command("Pluto:allowForcedDecay = on");
  command("ProcessLevel:all = off");
  command("ParticleDecays:limitTau0 = on");
  command("ParticleDecays:tau0Max = 10");
  command("TimeShower:QEDshowerByL = off");
  command("Random:setSeed = on");
  command("Random:seed = " + std::string(argc > 4 ? argv[4] : "12345"));
  if (!pythia.init()) return 3;
  std::vector<int> expected;
  if (mode == "4e") expected = {-11,-11,11,11};
  if (mode == "2mu2e") expected = {-13,-11,11,13};
  if (mode == "4mu") expected = {-13,-13,13,13};
  if (mode == "2e2pi") expected = {-211,-11,11,211};
  if (mode == "2mu2pi") expected = {-211,-13,13,211};
  for (int event=0; event<5; ++event) {
    pythia.event.reset();
    const double mass = pythia.particleData.m0(parentId);
    const double pz = event%2 ? 5. : 0.;
    const double energy = std::hypot(pz,mass);
    const int mother = pythia.event.append(parentId,1,0,0,0,0,0,0,0,0,pz,energy,mass);
    if (!pythia.moreDecays()) return 4;
    if (pythia.event[mother].isFinal()) return 5;
    std::vector<int> found;
    Pythia8::Vec4 sum;
    std::cout << "CHECK " << parentId << " " << mode << " " << event;
    for (const auto& daughter : pythia.event) {
      if (!daughter.isFinal()) continue;
      if (daughter.mother1() != mother) return 6;
      found.push_back(daughter.id());
      sum += daughter.p();
      std::cout << " " << daughter.id() << ":" << std::setprecision(17) << daughter.px();
    }
    std::cout << "\n";
    std::sort(found.begin(),found.end());
    if (found != expected) return 7;
    if (std::abs(sum.px())+std::abs(sum.py())+std::abs(sum.pz()-pz)+std::abs(sum.e()-energy)>1e-7*energy) return 8;
  }
}
