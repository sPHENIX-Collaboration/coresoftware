/*!
 *  \file PrelimDistortionCorrection.h
 *  \brief		Makes preliminary distortion corrections when crossing is unknown
 *  \author Tony Frawley
 */

#ifndef TRACKRECO_PRELIMDISTORTIONCORRECTION_H
#define TRACKRECO_PRELIMDISTORTIONCORRECTION_H

#include "ALICEKF.h"
#include "nanoflann.hpp"

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <tpc/TpcDistortionCorrection.h>
#include <trackbase/TrkrDefs.h>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

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
class TrkrClusterIterationMapv1;
class SvtxTrackMap;
class TrackSeedContainer;
class TrackSeed_v2;

using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;

class PrelimDistortionCorrection : public SubsysReco
{
 public:
  PrelimDistortionCorrection(const std::string &name = "PrelimDistortionCorrection");
  ~PrelimDistortionCorrection() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_field_dir(const double rescale)
  {
    _fieldDir = 1;
    if(rescale > 0)
      { _fieldDir = -1; }
  }
  void set_max_window(double s){_max_dist = s;}
  void useConstBField(bool opt){_use_const_field = opt;}
  void setConstBField(float b) { _const_field = b; }
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
  void publishSeeds(std::vector<TrackSeed_v2>& seeds, PositionMap &positions);

  /// tpc distortion correction utility class
  TpcDistortionCorrection m_distortionCorrection;

  bool _use_truth_clusters = false;

  /// fetch node pointers
  int get_nodes(PHCompositeNode *topNode);
  std::vector<double> radii;
  std::vector<double> _vertex_x;
  std::vector<double> _vertex_y;
  std::vector<double> _vertex_z;
  std::vector<double> _vertex_xerr;
  std::vector<double> _vertex_yerr;
  std::vector<double> _vertex_zerr;
  std::vector<double> _vertex_ids;
  //double _Bz = 1.4*_Bzconst;
  double _max_dist = .05;
  size_t _min_clusters_per_track = 3;
  double _fieldDir = -1;
  double _max_sin_phi = 1.;
  bool _pp_mode = false;

  TrkrClusterContainer *_cluster_map = nullptr;

  TrackSeedContainer *_track_map = nullptr;

  std::unique_ptr<PHField> _field_map = nullptr;

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

  bool _use_const_field = false;
  float _const_field = 1.4;
  bool _use_fixed_clus_err = false;
  std::array<double,3> _fixed_clus_err = {.1,.1,.1};

  double Ne_frac = 0.00;
  double Ar_frac = 0.75;
  double CF4_frac = 0.20;
  double N2_frac = 0.00;
  double isobutane_frac = 0.05;

};

#endif
