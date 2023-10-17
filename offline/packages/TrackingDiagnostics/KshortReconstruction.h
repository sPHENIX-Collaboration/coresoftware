#ifndef KSHORTRECONSTRUCTION_H
#define KSHORTRECONSTRUCTION_H

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>

#include <Acts/Definitions/Algebra.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunused-value"
#include <Acts/Propagator/Propagator.hpp>
#pragma GCC diagnostic pop

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Utilities/Result.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <ActsExamples/EventData/Trajectories.hpp>
#pragma GCC diagnostic pop

#include <Eigen/Dense>

class TFile;
class TH1D;
class TNtuple;

using BoundTrackParam = const Acts::BoundTrackParameters;
using BoundTrackParamResult = Acts::Result<BoundTrackParam>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
using Trajectory = ActsExamples::Trajectories;

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;

class KshortReconstruction : public SubsysReco
{
 public:
  KshortReconstruction(const std::string& name = "KshortReconstruction");
  virtual ~KshortReconstruction() {}

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* /*topNode*/) override;

  void setPtCut(double ptcut) { invariant_pt_cut = ptcut; }
  void setTrackQualityCut(double cut) { _qual_cut = cut; }
  void setPairDCACut(double cut) { pair_dca_cut = cut; }
  void setTrackDCACut(double cut) { track_dca_cut = cut; }
  void setRequireMVTX(bool set) { _require_mvtx = set; }
  void setDecayMass(Float_t decayMassSet) { decaymass = decayMassSet; }  //(muons decaymass = 0.1057) (pions = 0.13957) (electron = 0.000511)
  void set_output_file(const std::string &outputfile) { filepath = outputfile; }

 private:
  void fillNtp(SvtxTrack* track1, SvtxTrack* track2, Acts::Vector3 dcavals1, Acts::Vector3 dcavals2, Acts::Vector3 pca_rel1, Acts::Vector3 pca_rel2, double pair_dca, double invariantMass, double invariantPt, float rapidity, float pseudorapidity, Eigen::Vector3d projected_pos1, Eigen::Vector3d projected_pos2, Eigen::Vector3d projected_mom1, Eigen::Vector3d projected_mom2, Acts::Vector3 pca_rel1_proj, Acts::Vector3 pca_rel2_proj, double pair_dca_proj);

  void fillHistogram(Eigen::Vector3d mom1, Eigen::Vector3d mom2, TH1D* massreco, double& invariantMass, double& invariantPt, float& rapidity, float& pseudorapidity);

  // void findPcaTwoTracks(SvtxTrack *track1, SvtxTrack *track2, Acts::Vector3& pca1, Acts::Vector3& pca2, double& dca);
  void findPcaTwoTracks(Acts::Vector3 pos1, Acts::Vector3 pos2, Acts::Vector3 mom1, Acts::Vector3 mom2, Acts::Vector3& pca1, Acts::Vector3& pca2, double& dca);

  int getNodes(PHCompositeNode* topNode);

  Acts::Vector3 calculateDca(SvtxTrack* track, const Acts::Vector3& momentum, Acts::Vector3 position);

  bool projectTrackToCylinder(SvtxTrack* track, double Radius, Eigen::Vector3d& pos, Eigen::Vector3d& mom);
  bool projectTrackToPoint(SvtxTrack* track, Eigen::Vector3d PCA, Eigen::Vector3d& pos, Eigen::Vector3d& mom);

  Acts::Vector3 getVertex(SvtxTrack* track);

  TNtuple* ntp_reco_info = nullptr;
  ActsGeometry* _tGeometry = nullptr;
  SvtxTrackMap* m_svtxTrackMap = nullptr;
  SvtxVertexMap* m_vertexMap = nullptr;

  std::string filepath = "";
  Float_t decaymass = 0.13957;  // pion decay mass
  bool _require_mvtx = true;
  double _qual_cut = 10.0;
  double pair_dca_cut = 0.05;  // kshort relative cut 500 microns
  double track_dca_cut = 0.01;
  double invariant_pt_cut = 0.1;
  TFile* fout = nullptr;
  TH1D* recomass = nullptr;
};

#endif  // KSHORTRECONSTRUCTION_H
