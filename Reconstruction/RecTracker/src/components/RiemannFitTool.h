#ifndef RECTRACKER_RIEMANNFITTOOL_H
#define RECTRACKER_RIEMANNFITTOOL_H

// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/PhysicalConstants.h"

// FCCSW
#include "FWCore/DataHandle.h"
#include "RecInterface/ITrackFittingTool.h"
#include "datamodel/PositionedTrackHitCollection.h"
#include "datamodel/TrackCollection.h"
#include "datamodel/TrackHitCollection.h"
#include "datamodel/TrackStateCollection.h"

#include <map>

/** @class RiemannFitTool
 * Track fitting tool implementation using tricktracks' Riemannfit
 */
class RiemannFitTool : public GaudiTool, virtual public ITrackFittingTool {
public:
  RiemannFitTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~RiemannFitTool() = default;
  virtual StatusCode initialize() override final;
  virtual StatusCode finalize() override final;
  virtual std::pair<fcc::TrackCollection*, fcc::TrackStateCollection*>
  fitTracks(const fcc::PositionedTrackHitCollection* theHits,
            std::multimap<unsigned int, unsigned int> seedmap) override final;

  bool trackSelection(const fcc::Track& track, const fcc::TrackState& trackstate);

private:
  Gaudi::Property<double> m_Bz{this, "Bz", 4. * Gaudi::Units::tesla, "Field strength along Z"};
  Gaudi::Property<double> m_hitRes{this, "hitRes", 1e-8 * Gaudi::Units::mm, "Resolution of local hit coordinates"};
  Gaudi::Property<bool> m_doFit{this, "doFit", true, "flag to actually perform the fit"};
  Gaudi::Property<bool> m_calcErrors{this, "calcErrors", true, "flag to actually calculate errors"};
  Gaudi::Property<bool> m_fitOnlyPrimary{this, "fitOnlyPrimary", false, "flag to only fit the particle with trackID 1"};
  Gaudi::Property<bool> m_calcMultipleScattering{this, "MultipleScatteringErrors", false,
                                                 "flag to toggle estimation of multiple scattering errors"};
  Gaudi::Property<bool> m_saveOnlyValidFit{this, "saveOnlyValidFit", false, "flag to save only tracks that succeeds fitting"};

  Gaudi::Property<bool> m_applyTrackSelection{this, "applyTrackSelection", false, "flag to apply track selection to reduce saved tracks"};
  Gaudi::Property<double> m_cut_track_pt{this, "cut_track_pt", 0. * Gaudi::Units::GeV, "Track pT"};
  Gaudi::Property<double> m_cut_track_eta{this, "cut_track_eta", 10., "Track eta"};
  Gaudi::Property<double> m_cut_track_d0{this, "cut_track_d0", 1e+5 * Gaudi::Units::mm, "Track d0"};
  Gaudi::Property<double> m_cut_track_z0{this, "cut_track_z0", 1e+5 * Gaudi::Units::mm, "Track z0"};
  Gaudi::Property<double> m_cut_track_chi2{this, "cut_track_chi2", 0. , "Track chi2"};
  Gaudi::Property<double> m_cut_track_chi2ndf{this, "cut_track_chi2ndf", 0. , "Track chi2 / ndf"};
};

#endif /* RECTRACKER_RIEMANNFITTOOL_H */
