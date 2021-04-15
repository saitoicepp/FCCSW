
#include "DetInterface/IGeoSvc.h"
#include "DD4hep/Detector.h"
#include "DD4hep/BitFieldCoder.h"
#include "DDSegmentation/CartesianGridXZ.h"

#include "datamodel/PositionedTrackHitCollection.h"
#include "datamodel/TrackHitCollection.h"
#include "datamodel/TrackHitCollection.h"


#include <cmath>
#include <random>
#include <map>

#include "FastSimpleDigi.h"
#include "RecTracker/TrackingUtils.h"

// todo: review unit handling
#ifdef HAVE_GEANT4_UNITS
#define MM_2_CM 1.0
#define CM_2_MM 1.0
#else
#define MM_2_CM 0.1
#define CM_2_MM 10.0
#endif


DECLARE_COMPONENT(FastSimpleDigi)



FastSimpleDigi::FastSimpleDigi(const std::string& name, ISvcLocator* svcLoc) :
    GaudiAlgorithm(name, svcLoc), m_geoSvc("GeoSvc", name) {

  declareProperty("smearedHits", m_smearedTrackHits, "Smeared Tracker hits (Output)");
  declareProperty("trackHits", m_trackHits, "Tracker hits (Input)");
}

StatusCode FastSimpleDigi::initialize() {
  info() << "initialize" << endmsg;

  

  StatusCode sc = GaudiAlgorithm::initialize();
  if (sc.isFailure()) return sc;

  return sc;
}

StatusCode FastSimpleDigi::execute() {

  auto lcdd = m_geoSvc->lcdd();
  auto readout = lcdd->readout(m_readoutName);
  auto m_decoder = readout.idSpec().decoder();
  auto segmentationXZ = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXZ*>(readout.segmentation().segmentation());
  if (nullptr == segmentationXZ) {
    error() << "Could not retrieve segmentation!" << endmsg;
    return StatusCode::FAILURE;
  }
  m_segGridSizeX = segmentationXZ->gridSizeX();
  m_segGridSizeZ = segmentationXZ->gridSizeZ();
  auto m_volman = lcdd->volumeManager();
  // get hits from event store
  const fcc::TrackHitCollection* hits = m_trackHits.get();
  std::vector<fcc::TrackHit> sortedHits;
  rec::sortTrackHits(hits, sortedHits, m_decoder);

  auto edmPositions = m_smearedTrackHits.createAndPut();
  std::map<std::pair<int, int>, std::vector<size_t>> hit_map; // (x,z) => index of sortedHits
  //for (const auto& hit : sortedHits) {
  for (size_t i_hit = 0; i_hit < sortedHits.size(); i_hit++) {
    const auto& hit = sortedHits.at(i_hit);
    const fcc::BareHit& hitCore = hit.core();
    dd4hep::DDSegmentation::CellID cID = hit.core().cellId;
    /// the local coordinates on the module
    // add segmentation info and smearing here
    int x = m_decoder->get(cID, "x");
    int z = m_decoder->get(cID, "z");
    auto key = std::make_pair(x, z);
    if(hit_map.count(key) == 0) hit_map[key] = std::vector<size_t>();
    hit_map[key].push_back(i_hit);
  }

  for(auto itr : hit_map){
    auto x = itr.first.first;
    auto z = itr.first.second;
    fcc::BareHit hitCore;
    hitCore.energy = 0.;
    hitCore.time = 0.;
    for(auto i_hit : itr.second){
        const auto& hit = sortedHits.at(i_hit);
        hitCore.cellId = hit.core().cellId; // override
        hitCore.energy += hit.core().energy; // sum
        hitCore.time += hit.core().time; // average
        hitCore.bits = hit.core().bits; // override
    }
    hitCore.time /= itr.second.size();

    dd4hep::DDSegmentation::CellID cID = hitCore.cellId;
    std::array<double, 3> localPos = {x * m_segGridSizeX, 0, z * m_segGridSizeZ};
    // global coordinates, will be filled by the transform
    std::array<double, 3> globalPos = {0, 0, 0};
    // direct lookup of transformation in the volume manager is broken in dd4hep
    auto detelement = m_volman.lookupDetElement(cID);
    const auto& localToGlobal = detelement.nominal().worldTransformation();
    localToGlobal.LocalToMaster(localPos.data(), globalPos.data());
    auto position = fcc::Point();
    // default dd4hep units differ from fcc ones
    position.x = globalPos[0] * CM_2_MM;
    position.y = globalPos[1] * CM_2_MM;
    position.z = globalPos[2] * CM_2_MM;
    edmPositions->create(position, hitCore);
  }
  return StatusCode::SUCCESS;
}

StatusCode FastSimpleDigi::finalize() {
  StatusCode sc = GaudiAlgorithm::finalize();
  return sc;
}
