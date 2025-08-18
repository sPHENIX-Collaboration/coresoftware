/*!
 *  \file PrelimDistortionCorrectionAuAu.h
 *  \brief	   Refits AuAu TPC seeds with distortion corrections and current calibrations
 *  \author Tony Frawley
 */

#ifndef TRACKRECO_PRELIMDISTORTIONCORRECTIONAUAU_H
#define TRACKRECO_PRELIMDISTORTIONCORRECTIONAUAU_H

#include "ALICEKF.h"
#include "nanoflann.hpp"

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <tpc/TpcDistortionCorrection.h>
#include <trackbase/TrkrDefs.h>

#include <Eigen/Core>

// STL includes
#include <memory>
#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class PHField;
class TpcDistortionCorrectionContainer;
class TrkrClusterContainer;
class SvtxTrackMap;
class TrackSeedContainer;
class TrackSeed_v2;

class PrelimDistortionCorrectionAuAu : public SubsysReco
{
 public:

  //! constructor
  PrelimDistortionCorrectionAuAu(const std::string &name = "PrelimDistortionCorrectionAuAu");

  //! default destructor
  ~PrelimDistortionCorrectionAuAu() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  // noop
  void set_field_dir(double)
  {}

  // noop
  void useConstBField(bool)
  {}

  // noop
  void setConstBField(float)
  {}

  void useFixedClusterError(bool opt){_use_fixed_clus_err = opt;}
  void setFixedClusterError(int i, double val){_fixed_clus_err.at(i) = val;}
  void use_truth_clusters(bool truth) { _use_truth_clusters = truth; }
  void set_pp_mode(bool mode) {_pp_mode = mode;}

  void setNeonFraction(double frac) { Ne_frac = frac; };
  void setArgonFraction(double frac) { Ar_frac = frac; };
  void setCF4Fraction(double frac) { CF4_frac = frac; };
  void setNitrogenFraction(double frac) { N2_frac = frac; };
  void setIsobutaneFraction(double frac) { isobutane_frac = frac; };

 private:

  //! put refitted seeds on map
  using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;
  void publishSeeds(std::vector<TrackSeed_v2>& seeds, const PositionMap &positions) const;

  /// tpc distortion correction utility class
  TpcDistortionCorrection m_distortionCorrection;

  bool _use_truth_clusters = false;

  /// fetch node pointers
  int get_nodes(PHCompositeNode *topNode);

  size_t _min_clusters_per_track = 3;
  double _max_sin_phi = 1.;
  bool _pp_mode = false;

  TrkrClusterContainer *_cluster_map = nullptr;

  TrackSeedContainer *_track_map = nullptr;

  //! magnetic field map
  PHField* _field_map = nullptr;

  /// acts geometry
  ActsGeometry *_tgeometry = nullptr;

  //!@name distortion correction containers
  //@{
  /** used in input to correct CM clusters before calculating residuals */
  TpcDistortionCorrectionContainer *m_dcc_module_edge{nullptr};
  TpcDistortionCorrectionContainer *m_dcc_static{nullptr};
  TpcDistortionCorrectionContainer *m_dcc_average{nullptr};
  //@}

  std::unique_ptr<ALICEKF> fitter;

  bool _use_fixed_clus_err = false;
  std::array<double,3> _fixed_clus_err = {.1,.1,.1};

  double Ne_frac = 0.00;
  double Ar_frac = 0.75;
  double CF4_frac = 0.20;
  double N2_frac = 0.00;
  double isobutane_frac = 0.05;

};

#endif
