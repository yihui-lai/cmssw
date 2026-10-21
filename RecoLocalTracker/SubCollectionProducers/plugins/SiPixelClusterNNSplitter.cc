// SiPixelClusterNNSplitter: splits merged pixel clusters in the cores of high-pT jets (Phase-2 Inner Tracker) with a
// neural network (ONNX).
//
// Candidates: clusters of at least two pixels, within jetDeltaR of a calorimeter jet above jetPtMin (directions from
// the vertex), with charge above minChargeRatio times the charge expected for one particle crossing the sensor
// along the vertex-to-cluster line.
//
// Network inputs per candidate (all charges in ke):
//   img [9, 16, 8]: the cluster's charge in a 16 x 8 pixel window centred on its charge centroid; then, for each of
//                   the barrel layers 1-4, the charge of the clusters near the jets on that layer, moved along their
//                   line from the vertex to the candidate's radius and binned in the candidate's window; then the
//                   same four layers at their own scale, each centred on that layer's brightest charge blob within
//                   3 rows x 6 columns of the candidate's centroid.
//   ctx [8]:        barrel flag, layer / 12, radius / 10 cm, (z - z_vertex) / 10 cm, eta, cluster size in y / 8, and
//                   the two incidence angles of the vertex-to-cluster line in the module frame (atan of the
//                   tangents, clipped at 20).
// Network outputs: frac [4, 16, 8], the share of each window pixel's charge for up to 4 particles; exist [4], 1 for
// the predicted number of particles; sigma_hit [4, 2], the uncertainty of each piece's position (window pixels).
// The pieces with at least minSubClusterCharge become clusters (fewer than two: the cluster is kept as it is);
// pixels outside the window go to the piece with the most charge. The predicted uncertainty is attached to each
// piece (setSplitClusterErrorX/Y, um) and used by PixelCPEGeneric with useSplitClusterErrors = true.
//
// Vertex: z from the first L1 track-trigger vertex (l1Vertices); if that is not configured or empty, z of the beam
// spot. The network is trained with the vertex z smeared by 0.6 mm, the L1 vertex resolution.
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/L1Trigger/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

class SiPixelClusterNNSplitter : public edm::stream::EDProducer<edm::GlobalCache<cms::Ort::ONNXRuntime>> {
public:
  explicit SiPixelClusterNNSplitter(const edm::ParameterSet& cfg, const cms::Ort::ONNXRuntime*);
  static std::unique_ptr<cms::Ort::ONNXRuntime> initializeGlobalCache(const edm::ParameterSet& cfg);
  static void globalEndJob(const cms::Ort::ONNXRuntime*) {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::Event& ev, const edm::EventSetup& es) override;

private:
  static constexpr int H = 16, W = 8, K = 4, NL = 4;
  static constexpr int NIMG = 1 + 2 * NL, NCTX = 8;
  static constexpr double kOneParticleCharge = 12000.;  // electrons, normal incidence

  const edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> clusterToken_;
  const edm::EDGetTokenT<edm::View<reco::Candidate>> jetToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<l1t::Vertex>> l1VertexToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const double jetPtMin_, jetDeltaR_, minChargeRatio_, existThreshold_, minSubClusterCharge_, minPixelCharge_;
  const int batchSize_;
};

SiPixelClusterNNSplitter::SiPixelClusterNNSplitter(const edm::ParameterSet& cfg, const cms::Ort::ONNXRuntime*)
    : clusterToken_(consumes(cfg.getParameter<edm::InputTag>("pixelClusters"))),
      jetToken_(consumes(cfg.getParameter<edm::InputTag>("jets"))),
      beamSpotToken_(consumes(cfg.getParameter<edm::InputTag>("beamSpot"))),
      geomToken_(esConsumes()),
      topoToken_(esConsumes()),
      jetPtMin_(cfg.getParameter<double>("jetPtMin")),
      jetDeltaR_(cfg.getParameter<double>("jetDeltaR")),
      minChargeRatio_(cfg.getParameter<double>("minChargeRatio")),
      existThreshold_(cfg.getParameter<double>("existThreshold")),
      minSubClusterCharge_(cfg.getParameter<double>("minSubClusterCharge")),
      minPixelCharge_(cfg.getParameter<double>("minPixelCharge")),
      batchSize_(cfg.getParameter<int>("batchSize")) {
  const auto l1 = cfg.getParameter<edm::InputTag>("l1Vertices");
  if (!l1.label().empty())
    l1VertexToken_ = consumes(l1);
  produces<edmNew::DetSetVector<SiPixelCluster>>();
}

std::unique_ptr<cms::Ort::ONNXRuntime> SiPixelClusterNNSplitter::initializeGlobalCache(const edm::ParameterSet& cfg) {
  return std::make_unique<cms::Ort::ONNXRuntime>(cfg.getParameter<edm::FileInPath>("model").fullPath());
}

void SiPixelClusterNNSplitter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pixelClusters", edm::InputTag("siPixelClustersUnsplit"));
  desc.add<edm::FileInPath>(
      "model",
      edm::FileInPath("RecoLocalTracker/SubCollectionProducers/data/SiPixelClusterNNSplitter_Phase2_v3b.onnx"));
  desc.add<edm::InputTag>("jets", edm::InputTag("ak4CaloJetsForNNSplit"))
      ->setComment("calorimeter jets without vertex correction (no dependence on tracking)");
  desc.add<double>("jetPtMin", 100.);
  desc.add<double>("jetDeltaR", 0.1)->setComment("split only clusters within this dR of a jet");
  desc.add<edm::InputTag>("l1Vertices", edm::InputTag("l1tVertexFinder", "L1Vertices"))
      ->setComment("vertex z from the first L1 vertex; empty tag or no vertex: z of the beam spot");
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  desc.add<double>("minChargeRatio", 1.5)->setComment("candidate: charge above this times one particle's charge");
  desc.add<double>("existThreshold", 0.5);
  desc.add<double>("minSubClusterCharge", 3000.)->setComment("electrons");
  desc.add<double>("minPixelCharge", 500.)->setComment("electrons");
  desc.add<int>("batchSize", 2048)->setComment("candidates per inference call (bounds the memory)");
  descriptions.add("siPixelClusterNNSplitter", desc);
}

void SiPixelClusterNNSplitter::produce(edm::Event& ev, const edm::EventSetup& es) {
  const auto& clusters = ev.get(clusterToken_);
  const auto& geom = es.getData(geomToken_);
  const auto& topo = es.getData(topoToken_);

  GlobalPoint pv(0, 0, ev.get(beamSpotToken_).z0());
  if (!l1VertexToken_.isUninitialized()) {
    const auto& l1 = ev.get(l1VertexToken_);
    if (!l1.empty())
      pv = GlobalPoint(0, 0, l1[0].z0());
  }

  std::vector<std::pair<float, float>> jets;  // (eta, phi)
  for (const auto& j : ev.get(jetToken_))
    if (j.pt() > jetPtMin_)
      jets.emplace_back(j.eta(), j.phi());
  auto nearJet = [&jets](const GlobalVector& dir, double dr) {
    for (const auto& [eta, phi] : jets)
      if (reco::deltaR(dir.eta(), dir.phi(), eta, phi) < dr)
        return true;
    return false;
  };
  auto out = std::make_unique<edmNew::DetSetVector<SiPixelCluster>>();
  if (jets.empty()) {
    for (const auto& dets : clusters) {
      edmNew::DetSetVector<SiPixelCluster>::FastFiller filler(*out, dets.detId());
      for (const auto& clu : dets)
        filler.push_back(clu);
      if (filler.empty())
        filler.abort();
    }
    ev.put(std::move(out));
    return;
  }

  // pixels of the barrel layers 1-4 near the jets, seen from the vertex (inputs of the neighbour-layer maps)
  struct BPix {
    GlobalVector rel;
    float rad, phi, cot, q;
    int layer;
  };
  std::vector<BPix> bpix;
  for (const auto& dets : clusters) {
    const DetId id(dets.detId());
    const auto* det = dynamic_cast<const PixelGeomDetUnit*>(geom.idToDetUnit(id));
    if (det == nullptr || id.subdetId() != PixelSubdetector::PixelBarrel)
      continue;
    const int layer = topo.layer(id);
    if (layer < 1 || layer > NL)
      continue;
    for (const auto& clu : dets) {
      const GlobalVector dir =
          det->surface().toGlobal(det->specificTopology().localPosition(MeasurementPoint(clu.x(), clu.y()))) - pv;
      if (!nearJet(dir, jetDeltaR_ + 0.1))
        continue;
      for (int i = 0; i < clu.size(); ++i) {
        const auto p = clu.pixel(i);
        const GlobalVector rel =
            det->surface().toGlobal(det->specificTopology().localPosition(MeasurementPoint(p.x + 0.5, p.y + 0.5))) - pv;
        const float rad = rel.perp();
        bpix.push_back({rel, rad, float(rel.phi()), float(rel.z() / rad), p.adc / 1000.f, layer - 1});
      }
    }
  }

  // pass 1: candidates and network inputs
  struct Cand {
    int r0, c0;  // window origin (row, column)
  };
  std::vector<std::vector<int>> candIdx;  // per module: candidate index of each cluster, -1 if none
  std::vector<Cand> cands;
  cms::Ort::FloatArrays inputs(2);
  for (const auto& dets : clusters) {
    const DetId id(dets.detId());
    const auto* det = dynamic_cast<const PixelGeomDetUnit*>(geom.idToDetUnit(id));
    const bool barrel = id.subdetId() == PixelSubdetector::PixelBarrel;
    std::vector<int> idx(dets.size(), -1);
    int ic = 0;
    for (const auto& clu : dets) {
      const int me = ic++;
      if (det == nullptr || clu.size() < 2)
        continue;
      const GlobalPoint pos =
          det->surface().toGlobal(det->specificTopology().localPosition(MeasurementPoint(clu.x(), clu.y())));
      const GlobalVector dir = pos - pv;
      if (!nearJet(dir, jetDeltaR_))
        continue;
      const double rad = pos.perp(), dz = pos.z() - pv.z();
      const double eta = std::asinh(dz / rad);
      const double path = barrel ? std::cosh(eta) : std::cosh(eta) / std::max(std::abs(std::sinh(eta)), 0.3);
      if (clu.charge() < minChargeRatio_ * kOneParticleCharge * std::min(path, 6.))
        continue;

      double sq = 0, sr = 0, sc = 0;
      for (int i = 0; i < clu.size(); ++i) {
        const auto p = clu.pixel(i);
        sq += p.adc;
        sr += p.x * double(p.adc);
        sc += p.y * double(p.adc);
      }
      const Cand c{int(std::lround(sr / sq)) - H / 2, int(std::lround(sc / sq)) - W / 2};
      std::vector<float> img(NIMG * H * W, 0.f);
      for (int i = 0; i < clu.size(); ++i) {
        const auto p = clu.pixel(i);
        const int r = p.x - c.r0, cc = p.y - c.c0;
        if (r >= 0 && r < H && cc >= 0 && cc < W)
          img[r * W + cc] = p.adc / 1000.f;
      }

      // neighbour-layer maps: projected to the candidate's radius (channels 1-4) and at each layer's scale (5-8)
      const double crad = dir.perp(), cphi = dir.phi(), ccot = dir.z() / crad;
      const auto pitch = det->specificTopology().pitch();
      struct NPix {
        double fr, fc, sc;
        float q;
        int layer;
      };
      std::vector<NPix> sel;
      for (const auto& b : bpix) {
        const double dphi = reco::deltaPhi(double(b.phi), cphi);
        const double srad = std::min(double(b.rad), crad) / 2;  // the native maps reach further on inner layers
        if (std::abs(dphi) * srad >= (H / 2 + 2) * pitch.first ||
            std::abs(b.cot - ccot) * srad >= (W / 2 + 2) * pitch.second)
          continue;
        const GlobalPoint q(pv + b.rel * float(crad / b.rad));  // moved to the candidate's radius
        const MeasurementPoint mp = det->specificTopology().measurementPosition(det->surface().toLocal(q));
        const double fr = mp.x() - 0.5 - c.r0, fc = mp.y() - 0.5 - c.c0;
        const int r = int(std::lround(fr)), cc = int(std::lround(fc));
        if (r >= 0 && r < H && cc >= 0 && cc < W)
          img[((1 + b.layer) * H + r) * W + cc] += b.q;
        sel.push_back({fr, fc, b.rad / crad, b.q, b.layer});
      }
      const double ci = sr / sq - c.r0, cj = sc / sq - c.c0;
      for (int l = 0; l < NL; ++l) {
        std::vector<const NPix*> box;
        for (const auto& n : sel)
          if (n.layer == l && std::abs(n.fr - ci) <= 3 && std::abs(n.fc - cj) <= 6)
            box.push_back(&n);
        double li = ci, lj = cj, best = -1;
        for (const auto* a : box) {  // most charge within 2 rows x 1 column at the layer's own scale
          double qa = 0, wr = 0, wc = 0;
          for (const auto* o : box)
            if (std::abs(a->fr - o->fr) * a->sc <= 2 && std::abs(a->fc - o->fc) * a->sc <= 1) {
              qa += o->q;
              wr += o->fr * o->q;
              wc += o->fc * o->q;
            }
          if (qa > best) {
            best = qa;
            li = wr / qa;
            lj = wc / qa;
          }
        }
        for (const auto& n : sel) {
          if (n.layer != l)
            continue;
          const int r = int(std::lround(ci + (n.fr - li) * n.sc)), cc = int(std::lround(cj + (n.fc - lj) * n.sc));
          if (r >= 0 && r < H && cc >= 0 && cc < W)
            img[((1 + NL + l) * H + r) * W + cc] += n.q;
        }
      }
      inputs[0].insert(inputs[0].end(), img.begin(), img.end());

      const LocalVector ld = det->surface().toLocal(dir);
      const float ta = ld.z() != 0 ? std::clamp(ld.x() / ld.z(), -20.f, 20.f) : 20.f;
      const float tb = ld.z() != 0 ? std::clamp(ld.y() / ld.z(), -20.f, 20.f) : 20.f;
      const float ctx[NCTX] = {barrel ? 1.f : 0.f,
                               float(topo.layer(id) / 12.),
                               float(rad / 10.),
                               float(dz / 10.),
                               float(eta),
                               float(clu.sizeY() / 8.),
                               std::atan(ta),
                               std::atan(tb)};
      inputs[1].insert(inputs[1].end(), ctx, ctx + NCTX);
      idx[me] = cands.size();
      cands.push_back(c);
    }
    candIdx.push_back(std::move(idx));
  }

  // inference in batches (one call over all candidates would need several GB at PU200)
  const std::vector<std::string> outNames{"frac", "exist", "sigma_hit"};
  cms::Ort::FloatArrays outputs(outNames.size());
  const int64_t n = cands.size();
  for (int64_t b0 = 0; b0 < n; b0 += batchSize_) {
    const int64_t nb = std::min<int64_t>(batchSize_, n - b0);
    cms::Ort::FloatArrays in{{inputs[0].begin() + b0 * NIMG * H * W, inputs[0].begin() + (b0 + nb) * NIMG * H * W},
                             {inputs[1].begin() + b0 * NCTX, inputs[1].begin() + (b0 + nb) * NCTX}};
    auto o = globalCache()->run({"img", "ctx"}, in, {{nb, NIMG, H, W}, {nb, NCTX}}, outNames, nb);
    for (size_t k = 0; k < outNames.size(); ++k)
      outputs[k].insert(outputs[k].end(), o[k].begin(), o[k].end());
  }

  // pass 2: output clusters
  int im = 0;
  for (const auto& dets : clusters) {
    const auto& idx = candIdx[im++];
    const auto* det = dynamic_cast<const PixelGeomDetUnit*>(geom.idToDetUnit(DetId(dets.detId())));
    const auto pitch = det ? det->specificTopology().pitch() : std::pair<float, float>(0.f, 0.f);  // cm
    edmNew::DetSetVector<SiPixelCluster>::FastFiller filler(*out, dets.detId());
    int ic = 0;
    for (const auto& clu : dets) {
      const int ci = idx[ic++];
      if (ci < 0) {
        filler.push_back(clu);
        continue;
      }
      const float* frac = &outputs[0][size_t(ci) * K * H * W];
      const float* exist = &outputs[1][size_t(ci) * K];
      const auto& c = cands[ci];
      std::vector<int> slots;  // predicted particles
      for (int k = 0; k < K; ++k)
        if (exist[k] > existThreshold_)
          slots.push_back(k);
      auto shares = [&](int i, std::vector<double>& s) {  // charge shares of pixel i over `slots`
        const auto p = clu.pixel(i);
        const int r = p.x - c.r0, cc = p.y - c.c0;
        s.assign(slots.size(), 0.);
        if (r < 0 || r >= H || cc < 0 || cc >= W)
          return false;
        double tot = 0;
        for (size_t j = 0; j < slots.size(); ++j)
          tot += (s[j] = frac[(slots[j] * H + r) * W + cc]);
        if (tot <= 0)
          return false;
        for (auto& x : s)
          x /= tot;
        return true;
      };
      // drop pieces with too little charge (twice, the shares change when a piece is dropped)
      std::vector<double> s;
      for (int pass = 0; pass < 2 && slots.size() >= 2; ++pass) {
        std::vector<double> q(slots.size(), 0.);
        for (int i = 0; i < clu.size(); ++i)
          if (shares(i, s))
            for (size_t j = 0; j < slots.size(); ++j)
              q[j] += clu.pixel(i).adc * s[j];
        std::vector<int> keep;
        for (size_t j = 0; j < slots.size(); ++j)
          if (q[j] >= minSubClusterCharge_)
            keep.push_back(slots[j]);
        if (keep.size() == slots.size())
          break;
        slots = keep;
      }
      if (slots.size() < 2) {
        filler.push_back(clu);
        continue;
      }
      // the dominant piece takes the pixels outside the window
      std::vector<double> q(slots.size(), 0.);
      for (int i = 0; i < clu.size(); ++i)
        if (shares(i, s))
          for (size_t j = 0; j < slots.size(); ++j)
            q[j] += clu.pixel(i).adc * s[j];
      const size_t dom = std::max_element(q.begin(), q.end()) - q.begin();
      std::vector<SiPixelCluster> subs(slots.size());
      std::vector<bool> empty(slots.size(), true);
      for (int i = 0; i < clu.size(); ++i) {
        const auto p = clu.pixel(i);
        if (!shares(i, s)) {
          s.assign(slots.size(), 0.);
          s[dom] = 1.;
        }
        for (size_t j = 0; j < slots.size(); ++j) {
          const int adc = std::min(65535, int(p.adc * s[j] + 0.5));
          if (adc < minPixelCharge_)
            continue;
          const SiPixelCluster::PixelPos ppos(p.x, p.y);
          if (empty[j]) {
            subs[j] = SiPixelCluster(ppos, adc);
            empty[j] = false;
          } else {
            subs[j].add(ppos, adc);
          }
        }
      }
      for (size_t j = 0; j < slots.size(); ++j) {
        if (empty[j])
          continue;
        const float* sig = &outputs[2][(size_t(ci) * K + slots[j]) * 2];  // window pixels -> um
        subs[j].setSplitClusterErrorX(sig[0] * pitch.first * 1e4f);
        subs[j].setSplitClusterErrorY(sig[1] * pitch.second * 1e4f);
        filler.push_back(subs[j]);
      }
    }
    if (filler.empty())
      filler.abort();
  }
  ev.put(std::move(out));
}

DEFINE_FWK_MODULE(SiPixelClusterNNSplitter);
