#ifndef RECTRACKER_FASTSIMPLEDIGI_H
#define RECTRACKER_FASTSIMPLEDIGI_H

//#include "DD4hep/Detector.h"
//#include "DD4hep/BitField64.h"

// GAUDI
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

// FCCSW
#include "FWCore/DataHandle.h"

class IGeoSvc;

namespace fcc {
class TrackHitCollection;
class PositionedTrackHitCollection;
}




class FastSimpleDigi : public GaudiAlgorithm {
public:
  FastSimpleDigi(const std::string& name, ISvcLocator* svcLoc);

  ~FastSimpleDigi() = default;

  StatusCode initialize() override final;

  StatusCode execute() override final;

  StatusCode finalize() override final;

private:
  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;

  float m_segGridSizeZ;
  float m_segGridSizeX;
  //dd4hep::DDSegmentation::BitField64* m_decoder;
  //dd4hep::VolumeManager m_volman;

  DataHandle<fcc::TrackHitCollection> m_trackHits{"trackHits", Gaudi::DataHandle::Reader, this};
  DataHandle<fcc::PositionedTrackHitCollection> m_smearedTrackHits{"smearedHits", Gaudi::DataHandle::Writer, this};
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "", "Name of the readout (hits collection) to save"};

};

#endif /* RECTRACKER_FASTSIMPLEDIGI_H */
