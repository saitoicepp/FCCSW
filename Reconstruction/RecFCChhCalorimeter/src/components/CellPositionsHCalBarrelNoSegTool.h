#ifndef RECCALORIMETER_CELLPOSITIONSHCALBARRELNOSEGTOOL_H
#define RECCALORIMETER_CELLPOSITIONSHCALBARRELNOSEGTOOL_H

// GAUDI
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/ServiceHandle.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "DetInterface/IGeoSvc.h"
#include "DetSegmentation/GridPhiEta.h"
#include "FWCore/DataHandle.h"
#include "RecInterface/ICellPositionsTool.h"

// DD4hep
#include "DD4hep/LCDD.h"
#include "DD4hep/Readout.h"
#include "DD4hep/Volumes.h"
#include "DDSegmentation/Segmentation.h"
#include "TGeoManager.h"

class IGeoSvc;
namespace DD4hep {
namespace DDSegmentation {
class Segmentation;
}
}

/** @class CellPositionsHCalBarrelNoSegTool
 *
 */

class CellPositionsHCalBarrelNoSegTool : public GaudiTool, virtual public ICellPositionsTool {
public:
  CellPositionsHCalBarrelNoSegTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsHCalBarrelNoSegTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  virtual void getPositions(const fcc::CaloHitCollection& aCells, fcc::PositionedCaloHitCollection& outputColl) final;

  virtual DD4hep::Geometry::Position xyzPosition(const uint64_t& aCellId) const final;

  virtual int layerId(const uint64_t& aCellId) final;

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "name of the readout"};

  DD4hep::DDSegmentation::BitField64* m_decoder;

  DD4hep::Geometry::VolumeManager m_volman;
};
#endif /* RECCALORIMETER_CELLPOSITIONSHCALBARRELNOSEGTOOL_H */
