#ifndef RECTRACKER_FCCHHSEEDINGGRAPHTOOL_H
#define RECTRACKER_FCCHHSEEDINGGRAPHTOOL_H

// from Gaudi
#include "GaudiAlg/GaudiTool.h"

// FCCSW
#include "FWCore/DataHandle.h"
#include "RecInterface/ILayerGraphTool.h"

namespace tricktrack {
class CMGraph;
}

class FCChhSeedingGraphTool : public GaudiTool, virtual public ILayerGraphTool {
public:
  FCChhSeedingGraphTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~FCChhSeedingGraphTool() = default;
  virtual StatusCode initialize() override final;
  virtual StatusCode finalize() override final;
  virtual tricktrack::CMGraph graph() override final;
private:
  Gaudi::Property<std::vector<std::vector<std::string>>> m_layerPaths {this, "layerPaths", {{
  
      {"b01", "b02", "b03", "b04", "b05", "b06", "b11", "b12", "b13", "b14", "b15", "b16"},
  
  }} };
};

#endif /* RECTRACKER_FCCHHSEEDINGGRAPHTOOL_H */
