//
//
// Package:    Phase2TrackerDigitizer
// Class:      Phase2TrackerDigitizer
//
// *\class SiPhase2TrackerDigitizer Phase2TrackerDigitizer.cc SimTracker/SiPhase2Digitizer/src/Phase2TrackerDigitizer.cc
//
// Author: Suchandra Dutta, Suvankar Roy Chowdhury, Subir Sarkar
// Date: January 29, 2016
// Description: <one line class summary>
//
// Implementation:
//     <Notes on implementation>
//

#include <memory>
#include <set>
#include <iostream>
#include "SimTracker/SiPhase2Digitizer/plugins/Phase2TrackerDigitizer.h"
#include "SimTracker/SiPhase2Digitizer/plugins/Phase2TrackerDigitizerAlgorithm.h"
#include "SimTracker/SiPhase2Digitizer/plugins/SSDigitizerAlgorithm.h"
#include "SimTracker/SiPhase2Digitizer/plugins/PSSDigitizerAlgorithm.h"
#include "SimTracker/SiPhase2Digitizer/plugins/PSPDigitizerAlgorithm.h"
#include "SimTracker/SiPhase2Digitizer/plugins/PixelDigitizerAlgorithm.h"
#include "SimTracker/SiPhase2Digitizer/plugins/Pixel3DDigitizerAlgorithm.h"
#include "SimTracker/SiPhase2Digitizer/plugins/DigitizerUtility.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "SimGeneral/MixingModule/interface/PileUpEventPrincipal.h"

// Random Number
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CLHEP/Random/RandomEngine.h"

namespace cms {

  Phase2TrackerDigitizer::Phase2TrackerDigitizer(const edm::ParameterSet& iConfig,
                                                 edm::ProducerBase& mixMod,
                                                 edm::ConsumesCollector& iC)
      : first_(true),
        hitsProducer_(iConfig.getParameter<std::string>("hitsProducer")),
        trackerContainers_(iConfig.getParameter<std::vector<std::string> >("ROUList")),
        geometryType_(iConfig.getParameter<std::string>("GeometryType")),
        isOuterTrackerReadoutAnalog(iConfig.getParameter<bool>("isOTreadoutAnalog")),
        premixStage1_(iConfig.getParameter<bool>("premixStage1")),
        makeDigiSimLinks_(
            iConfig.getParameter<edm::ParameterSet>("AlgorithmCommon").getUntrackedParameter<bool>("makeDigiSimLinks")) {
    //edm::LogInfo("Phase2TrackerDigitizer") << "Initialize Digitizer Algorithms";
    const std::string alias1("simSiPixelDigis");
    mixMod.produces<edm::DetSetVector<PixelDigi> >("Pixel").setBranchAlias(alias1);
    if (makeDigiSimLinks_) {
      mixMod.produces<edm::DetSetVector<PixelDigiSimLink> >("Pixel").setBranchAlias(alias1);
    }

    if (!iConfig.getParameter<bool>("isOTreadoutAnalog")) {
      const std::string alias2("simSiTrackerDigis");
      if (premixStage1_) {
        // Premixing exploits the ADC field of PixelDigi to store the collected charge
        // But we still want everything else to be treated like for Phase2TrackerDigi
        mixMod.produces<edm::DetSetVector<PixelDigi> >("Tracker").setBranchAlias(alias2);
      } else {
        mixMod.produces<edm::DetSetVector<Phase2TrackerDigi> >("Tracker").setBranchAlias(alias2);
      }
      if (makeDigiSimLinks_) {
        mixMod.produces<edm::DetSetVector<PixelDigiSimLink> >("Tracker").setBranchAlias(alias2);
      }
    }
    // creating algorithm objects and pushing them into the map
    algomap_[AlgorithmType::InnerPixel] = std::make_unique<PixelDigitizerAlgorithm>(iConfig);
    algomap_[AlgorithmType::InnerPixel3D] = std::make_unique<Pixel3DDigitizerAlgorithm>(iConfig);
    algomap_[AlgorithmType::PixelinPS] = std::make_unique<PSPDigitizerAlgorithm>(iConfig);
    algomap_[AlgorithmType::StripinPS] = std::make_unique<PSSDigitizerAlgorithm>(iConfig);
    algomap_[AlgorithmType::TwoStrip] = std::make_unique<SSDigitizerAlgorithm>(iConfig);
  }

  void Phase2TrackerDigitizer::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& iSetup) {
    iSetup.get<IdealMagneticFieldRecord>().get(pSetup_);
    iSetup.get<TrackerTopologyRcd>().get(tTopoHand);

    if (theTkDigiGeomWatcher.check(iSetup)) {
      iSetup.get<TrackerDigiGeometryRecord>().get(geometryType_, pDD_);
      //reset cache
      ModuleTypeCache().swap(moduleTypeCache_);
      detectorUnits_.clear();
      for (auto const& det_u : pDD_->detUnits()) {
        unsigned int detId_raw = det_u->geographicalId().rawId();
        DetId detId = DetId(detId_raw);
        if (DetId(detId).det() == DetId::Detector::Tracker) {
          const Phase2TrackerGeomDetUnit* pixdet = dynamic_cast<const Phase2TrackerGeomDetUnit*>(det_u);
          assert(pixdet);
          detectorUnits_.insert(std::make_pair(detId_raw, pixdet));
        }
      }
    }
  }

  Phase2TrackerDigitizer::~Phase2TrackerDigitizer() {}
  void Phase2TrackerDigitizer::accumulatePixelHits(edm::Handle<std::vector<PSimHit> > hSimHits,
                                                   size_t globalSimHitIndex,
                                                   const unsigned int tofBin) {
    if (hSimHits.isValid()) {
      std::set<unsigned int> detIds;
      std::vector<PSimHit> const& simHits = *hSimHits.product();
      int indx = 0;
      for (auto it = simHits.begin(), itEnd = simHits.end(); it != itEnd; ++it, ++globalSimHitIndex) {
        unsigned int detId_raw = (*it).detUnitId();
        if (detectorUnits_.find(detId_raw) == detectorUnits_.end())
          continue;
        if (detIds.insert(detId_raw).second) {
          // The insert succeeded, so this detector element has not yet been processed.
          AlgorithmType algotype = getAlgoType(detId_raw);
          const Phase2TrackerGeomDetUnit* phase2det = detectorUnits_[detId_raw];
          // access to magnetic field in global coordinates
          GlobalVector bfield = pSetup_->inTesla(phase2det->surface().position());
          LogDebug("PixelDigitizer") << "B-field(T) at " << phase2det->surface().position()
                                     << "(cm): " << pSetup_->inTesla(phase2det->surface().position());
          if (algomap_.find(algotype) != algomap_.end()) {
            algomap_[algotype]->accumulateSimHits(it, itEnd, globalSimHitIndex, tofBin, phase2det, bfield);
          } else
            edm::LogInfo("Phase2TrackerDigitizer") << "Unsupported algorithm: ";
        }
        indx++;
      }
    }
  }

  void Phase2TrackerDigitizer::initializeEvent(edm::Event const& e, edm::EventSetup const& iSetup) {
    edm::Service<edm::RandomNumberGenerator> rng;
    if (!rng.isAvailable()) {
      throw cms::Exception("Configuration")
          << "Phase2TrackerDigitizer requires the RandomNumberGeneratorService\n"
             "which is not present in the configuration file.  You must add the service\n"
             "in the configuration file or remove the modules that require it.";
    }

    // Must initialize all the algorithms
    for (auto const& el : algomap_) {
      if (first_)
        el.second->init(iSetup);

      el.second->initializeEvent(rng->getEngine(e.streamID()));
    }
    first_ = false;
    // Make sure that the first crossing processed starts indexing the sim hits from zero.
    // This variable is used so that the sim hits from all crossing frames have sequential
    // indices used to create the digi-sim link (if configured to do so) rather than starting
    // from zero for each crossing.
    crossingSimHitIndexOffset_.clear();
  }
  void Phase2TrackerDigitizer::accumulate(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
    accumulate_local<edm::Event>(iEvent, iSetup);
  }

  void Phase2TrackerDigitizer::accumulate(PileUpEventPrincipal const& iEvent,
                                          edm::EventSetup const& iSetup,
                                          edm::StreamID const&) {
    accumulate_local<PileUpEventPrincipal>(iEvent, iSetup);
  }

  template <class T>
  void Phase2TrackerDigitizer::accumulate_local(T const& iEvent, edm::EventSetup const& iSetup) {
    for (auto const& v : trackerContainers_) {
      edm::Handle<std::vector<PSimHit> > simHits;
      edm::InputTag tag(hitsProducer_, v);
      iEvent.getByLabel(tag, simHits);

      //edm::EDGetTokenT< std::vector<PSimHit> > simHitToken_(consumes< std::vector<PSimHit>(tag));
      //iEvent.getByToken(simHitToken_, simHits);

      unsigned int tofBin = PixelDigiSimLink::LowTof;
      if (v.find(std::string("HighTof")) != std::string::npos)
        tofBin = PixelDigiSimLink::HighTof;
      accumulatePixelHits(simHits, crossingSimHitIndexOffset_[tag.encode()], tofBin);
      // Now that the hits have been processed, I'll add the amount of hits in this crossing on to
      // the global counter. Next time accumulateStripHits() is called it will count the sim hits
      // as though they were on the end of this collection.
      // Note that this is only used for creating digi-sim links (if configured to do so).
      if (simHits.isValid())
        crossingSimHitIndexOffset_[tag.encode()] += simHits->size();
    }
  }

  // For premixing
  void Phase2TrackerDigitizer::loadAccumulator(const std::map<unsigned int, std::map<int, float> >& accumulator) {
    for (const auto& detMap : accumulator) {
      AlgorithmType algoType = getAlgoType(detMap.first);
      auto& algo = *(algomap_.at(algoType));
      algo.loadAccumulator(detMap.first, detMap.second);
    }
  }

  void Phase2TrackerDigitizer::finalizeEvent(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    //Decide if we want analog readout for Outer Tracker.
    addPixelCollection(iEvent, iSetup, isOuterTrackerReadoutAnalog);
    if (!isOuterTrackerReadoutAnalog) {
      if (premixStage1_) {
        addOuterTrackerCollection<PixelDigi>(iEvent, iSetup);
      } else {
        addOuterTrackerCollection<Phase2TrackerDigi>(iEvent, iSetup);
      }
    }
  }
  Phase2TrackerDigitizer::AlgorithmType Phase2TrackerDigitizer::getAlgoType(unsigned int detId_raw) {
    DetId detId(detId_raw);

    AlgorithmType algotype = AlgorithmType::Unknown;

    //get mType either from the geometry or from our cache (faster)
    TrackerGeometry::ModuleType mType = TrackerGeometry::ModuleType::UNKNOWN;
    auto itr = moduleTypeCache_.find(detId_raw);
    if (itr != moduleTypeCache_.end()) {
      mType = itr->second;
    } else {
      mType = pDD_->getDetectorType(detId);
      moduleTypeCache_.emplace(detId_raw, mType);
    }
    
    // Helper var to switch between InnerPixel or InnerPixel3D digitizers
    bool is_inner = false;

    switch (mType) {
      case TrackerGeometry::ModuleType::Ph1PXB:
        is_inner = true;
        //algotype = AlgorithmType::InnerPixel;
        break;
      case TrackerGeometry::ModuleType::Ph1PXF:
        is_inner = true;
        //algotype = AlgorithmType::InnerPixel;
        break;
      case TrackerGeometry::ModuleType::Ph2PXB:
        is_inner = true;
        //algotype = AlgorithmType::InnerPixel;
        break;
      case TrackerGeometry::ModuleType::Ph2PXF:
        is_inner = true;
        //algotype = AlgorithmType::InnerPixel;
        break;
      case TrackerGeometry::ModuleType::Ph2PSP:
        algotype = AlgorithmType::PixelinPS;
        break;
      case TrackerGeometry::ModuleType::Ph2PSS:
        algotype = AlgorithmType::StripinPS;
        break;
      case TrackerGeometry::ModuleType::Ph2SS:
        algotype = AlgorithmType::TwoStrip;
        break;
      default:
        edm::LogError("Phase2TrackerDigitizer") << "ERROR - Wrong Detector Type, No Algorithm available ";
    }

    if(is_inner)
    {
        const TrackerTopology* tTopo = tTopoHand.product();
        if(tTopo->getITPixelLayerNumber(detId) <= 2)
        {
            algotype = AlgorithmType::InnerPixel3D;
        }
        else
        {
            algotype = AlgorithmType::InnerPixel;
        }
    }

    return algotype;
  }
  void Phase2TrackerDigitizer::addPixelCollection(edm::Event& iEvent,
                                                  const edm::EventSetup& iSetup,
                                                  const bool ot_analog) {
    const TrackerTopology* tTopo = tTopoHand.product();
    std::vector<edm::DetSet<PixelDigi> > digiVector;
    std::vector<edm::DetSet<PixelDigiSimLink> > digiLinkVector;
    for (auto const& det_u : pDD_->detUnits()) {
      DetId detId_raw = DetId(det_u->geographicalId().rawId());
      AlgorithmType algotype = getAlgoType(detId_raw);
      if (algomap_.find(algotype) == algomap_.end())
        continue;

      //Decide if we want analog readout for Outer Tracker.
      if (!ot_analog 
              && (algotype != AlgorithmType::InnerPixel && algotype != AlgorithmType::InnerPixel3D) )
      {
          continue;
      }
      std::map<int, DigitizerUtility::DigiSimInfo> digi_map;
      algomap_[algotype]->digitize(dynamic_cast<const Phase2TrackerGeomDetUnit*>(det_u), digi_map, tTopo);
      edm::DetSet<PixelDigi> collector(det_u->geographicalId().rawId());
      edm::DetSet<PixelDigiSimLink> linkcollector(det_u->geographicalId().rawId());
      for (auto const& digi_p : digi_map) {
        DigitizerUtility::DigiSimInfo info = digi_p.second;
        std::pair<int, int> ip = PixelDigi::channelToPixel(digi_p.first);
        collector.data.emplace_back(ip.first, ip.second, info.sig_tot);
        for (auto const& sim_p : info.simInfoList) {
          linkcollector.data.emplace_back(digi_p.first,
                                          sim_p.second->trackId(),
                                          sim_p.second->hitIndex(),
                                          sim_p.second->tofBin(),
                                          sim_p.second->eventId(),
                                          sim_p.first);
        }
      }
      if (!collector.data.empty())
        digiVector.push_back(std::move(collector));
      if (!linkcollector.data.empty())
        digiLinkVector.push_back(std::move(linkcollector));
    }

    // Step C: create collection with the cache vector of DetSet
    std::unique_ptr<edm::DetSetVector<PixelDigi> > output(new edm::DetSetVector<PixelDigi>(digiVector));
    std::unique_ptr<edm::DetSetVector<PixelDigiSimLink> > outputlink(
        new edm::DetSetVector<PixelDigiSimLink>(digiLinkVector));

    // Step D: write output to file
    iEvent.put(std::move(output), "Pixel");
    if (makeDigiSimLinks_) {
      iEvent.put(std::move(outputlink), "Pixel");
    }
  }
}  // namespace cms
namespace {
  void addToCollector(edm::DetSet<PixelDigi>& collector, const int channel, const DigitizerUtility::DigiSimInfo& info) {
    // For premixing stage1 the channel must be decoded with PixelDigi
    // so that when the row and column are inserted to PixelDigi the
    // coded channel stays the same (so that it can then be decoded
    // with Phase2TrackerDigi in stage2).
    std::pair<int, int> ip = PixelDigi::channelToPixel(channel);
    collector.data.emplace_back(ip.first, ip.second, info.sig_tot);
  }
  void addToCollector(edm::DetSet<Phase2TrackerDigi>& collector,
                      const int channel,
                      const DigitizerUtility::DigiSimInfo& info) {
    std::pair<int, int> ip = Phase2TrackerDigi::channelToPixel(channel);
    collector.data.emplace_back(ip.first, ip.second, info.ot_bit);
  }
}  // namespace
namespace cms {
  template <typename DigiType>
  void Phase2TrackerDigitizer::addOuterTrackerCollection(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    const TrackerTopology* tTopo = tTopoHand.product();
    std::vector<edm::DetSet<DigiType> > digiVector;
    std::vector<edm::DetSet<PixelDigiSimLink> > digiLinkVector;
    for (auto const& det_u : pDD_->detUnits()) {
      DetId detId_raw = DetId(det_u->geographicalId().rawId());
      AlgorithmType algotype = getAlgoType(detId_raw);

      if (algomap_.find(algotype) == algomap_.end() 
              || algotype == AlgorithmType::InnerPixel 
              || algotype == AlgorithmType::InnerPixel3D )
      {
          continue;
      }

      std::map<int, DigitizerUtility::DigiSimInfo> digi_map;
      algomap_[algotype]->digitize(dynamic_cast<const Phase2TrackerGeomDetUnit*>(det_u), digi_map, tTopo);
      edm::DetSet<DigiType> collector(det_u->geographicalId().rawId());
      edm::DetSet<PixelDigiSimLink> linkcollector(det_u->geographicalId().rawId());

      for (auto const& digi_p : digi_map) {
        DigitizerUtility::DigiSimInfo info = digi_p.second;
        addToCollector(collector, digi_p.first, info);
        for (auto const& sim_p : info.simInfoList) {
          linkcollector.data.emplace_back(digi_p.first,
                                          sim_p.second->trackId(),
                                          sim_p.second->hitIndex(),
                                          sim_p.second->tofBin(),
                                          sim_p.second->eventId(),
                                          sim_p.first);
        }
      }

      if (!collector.data.empty())
        digiVector.push_back(std::move(collector));
      if (!linkcollector.data.empty())
        digiLinkVector.push_back(std::move(linkcollector));
    }

    // Step C: create collection with the cache vector of DetSet
    std::unique_ptr<edm::DetSetVector<DigiType> > output(new edm::DetSetVector<DigiType>(digiVector));
    std::unique_ptr<edm::DetSetVector<PixelDigiSimLink> > outputlink(
        new edm::DetSetVector<PixelDigiSimLink>(digiLinkVector));

    // Step D: write output to file
    iEvent.put(std::move(output), "Tracker");
    if (makeDigiSimLinks_) {
      iEvent.put(std::move(outputlink), "Tracker");
    }
  }
}  // namespace cms

#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixModFactory.h"

using cms::Phase2TrackerDigitizer;
DEFINE_DIGI_ACCUMULATOR(Phase2TrackerDigitizer);
