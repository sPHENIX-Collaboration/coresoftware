// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  SecondaryVertexFinder
 *  \brief		Class for determining primary vertices usin g fitted tracks
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef SECONDARYVERTEXFINDER_H
#define SECONDARYVERTEXFINDER_H

#include <ActsPropagator.h>

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>

#include <Acts/Definitions/Algebra.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Acts/Propagator/Propagator.hpp>
#pragma GCC diagnostic pop

#include <Acts/Utilities/Result.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <ActsExamples/EventData/Trajectories.hpp>

#include <string>
#include <vector>
#include <map>
#include <set>

#include <Eigen/Dense>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class TrkrCluster;
class TrackSeed;
class ActsGeometry;
class TrkrClusterContainer;
class TNtuple;
class TH2D;
class TH1D;
class TFile;
class TLorentzVector;

using BoundTrackParam = const Acts::BoundTrackParameters;
using BoundTrackParamResult = Acts::Result<BoundTrackParam>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
using Trajectory = ActsExamples::Trajectories;

class SecondaryVertexFinder : public SubsysReco
{
 public:

  SecondaryVertexFinder(const std::string &name = "SecondaryVertexFinder");

  ~SecondaryVertexFinder() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void setTrackDcaCut(const double cutxy, const double cutz) {_track_dcaxy_cut = cutxy; _track_dcaz_cut = cutz;}
 void setTwoTrackDcaCut(const double cut) {_two_track_dcacut = cut;}
 void setMinPathCut(const double cut) {_min_path_cut = cut;}
 void setTrackQualityCut(double cut) {_qual_cut = cut;}
 void setRequireMVTX(bool set) {_require_mvtx = set;}
 void setOutfileName(const std::string& filename) {outfile = filename;}
 void setDecayParticleMass(double mass) {_decaymass = mass;}
 void set_write_electrons_node(bool flag) {_write_electrons_node = flag;}
 void set_write_ntuple(bool flag) {_write_ntuple = flag;}

 private:

  int GetNodes(PHCompositeNode* topNode);
  int CreateOutputNode(PHCompositeNode* topNode);

bool passConversionElectronCuts(TLorentzVector tsum, 
				SvtxTrack* tr1, SvtxTrack* tr2, float pair_dca, 
				Eigen::Vector3d PCA, Eigen::Vector3d VTX);

  bool  hasSiliconSeed(SvtxTrack* tr); 
  void outputTrackDetails(SvtxTrack *tr);
  void get_dca(SvtxTrack* track, float& dca3dxy, float& dca3dz, float& dca3dxysigma, float& dca3dzsigma);
  bool circle_circle_intersection(double r0, double x0, double y0, double r1, double x1, double y1, std::vector<double>& intersectionXY);
 
  void findPcaTwoLines(Eigen::Vector3d pos1, Eigen::Vector3d mom1, Eigen::Vector3d pos2, Eigen::Vector3d mom2,
		       double &dca, Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2);
  void getCircleXYTrack(SvtxTrack *track, double& R, Eigen::Vector2d& center);
  double getZFromIntersectionXY(SvtxTrack *track, double& R, Eigen::Vector2d& center, Eigen::Vector2d intersection);
  bool projectTrackToPoint(SvtxTrack* track, Eigen::Vector3d& PCA, Eigen::Vector3d& pos, Eigen::Vector3d& mom);

  void fillNtp(SvtxTrack *track1, SvtxTrack *track2, double dca3dxy1, double dca3dz1, double dca3dxy2, double dca3dz2,  Eigen::Vector3d vpos1,  Eigen::Vector3d vmom1, Eigen::Vector3d vpos2, Eigen::Vector3d vmom2, Acts::Vector3 pca_rel1, Acts::Vector3 pca_rel2, double pair_dca, double invariantMass, double invariantPt, double path, int has_silicon_1, int has_siilicon_2);

  SvtxTrackMap *_track_map{nullptr};
  SvtxTrackMap *_track_map_electrons{nullptr};
  SvtxTrack *_track{nullptr};  
  SvtxVertexMap *_svtx_vertex_map{nullptr};
  ActsGeometry *_tGeometry{nullptr};

 bool _require_mvtx = false;
 bool _write_electrons_node = true;
 bool _write_ntuple = false;

 double _decaymass = 0.000511;  // conversion electrons, default 

 // these are minimal cuts used to make the ntuple
 // They can be tightened later when analyzing the ntuple

 // single track cuts
 double _track_dcaxy_cut = 0.020;  
 double _track_dcaz_cut = 0.020;  
 double _qual_cut = 4.0;

 //track_pair cuts
 double _two_track_dcacut = 0.5;  // 5000 microns 
 double _max_intersection_radius = 40.0;  // discard intersections at greater than 40 cm radius
 double _projected_track_z_cut = 1.0;

 // decay vertex cuts
 double _min_path_cut = 0.2;
 double _costheta_cut = 0.9985;

 // specific conversion electron cuts
 double _conversion_pair_dcacut = 0.2;  // 2000 microns 
 unsigned int _min_tpc_clusters = 40;
 double _deta_cut = 0.05;
 double _invariant_pt_cut = 0.1;
 double _max_mass_cut = 0.03;

  TH2D *recomass{nullptr};
  TH2D *hdecaypos{nullptr};
  TH1D *hdecay_radius{nullptr};
  TNtuple *ntp{nullptr};
  std::string outfile;
  
};

#endif //SECONDARYVERTEXFINDER_H
