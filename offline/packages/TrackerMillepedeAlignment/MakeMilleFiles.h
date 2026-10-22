// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  MakeMilleFiles.h
 *  \brief		Class for moving corrected TPC clusters to the nearest TPC readout layer radius
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef MAKEMILLEFILES_H
#define MAKEMILLEFILES_H

#include "AlignmentDefs.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxAlignmentState.h>
#include <trackbase_historic/SvtxAlignmentStateMap.h>

#include <fun4all/SubsysReco.h>

#include <ActsExamples/EventData/Trajectories.hpp>

#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

class PHCompositeNode;
class PHG4TpcGeomContainer;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class SvtxVertex;
class TrkrCluster;
class TrkrClusterContainer;
class Mille;
class ActsPropagator;

class TFile;
class TH1D;
class TH2D;
class TProfile;
class TNtuple;

using Trajectory = ActsExamples::Trajectories;

class MakeMilleFiles : public SubsysReco
{
 public:
  MakeMilleFiles(const std::string& name = "MakeMilleFiles");

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_binary(bool bin) { _binary = bin; }

  void set_track_map_name(const std::string& name) {m_track_map_name = name;}
  void set_state_map_name(const std::string& name) {m_state_map_name = name;}
  void set_constraintfile_name(const std::string& file) { m_constraintFileName = file; }
  void set_datafile_name(const std::string& file) { data_outfilename = file; }
  void set_steeringfile_name(const std::string& file) { steering_outfilename = file; }
  void set_tfile_name(const std::string& file) { m_tfile_name = file; }

  void set_mvtx_grouping(int group) { mvtx_group = (AlignmentDefs::mvtxGrp) group; }
  void set_intt_grouping(int group) { intt_group = (AlignmentDefs::inttGrp) group; }
  void set_tpc_grouping(int group) { tpc_group = (AlignmentDefs::tpcGrp) group; }
  void set_mms_grouping(int group) { mms_group = (AlignmentDefs::mmsGrp) group; }
  void use_event_vertex(bool useit) { m_useEventVertex = useit; }
  void set_layer_fixed(unsigned int layer);
  void set_mvtx_layer_fixed(unsigned int layer, unsigned int clamshell);
  void set_layer_gparam_fixed(unsigned int layer, unsigned int param);
  void set_layer_lparam_fixed(unsigned int layer, unsigned int param);

  void set_pt_cut(float pt) { m_minPt = pt; }
  void set_min_nmvtx(unsigned int nmvtx) { m_minNmvtx = nmvtx; }
  void set_min_nintt(unsigned int nintt) { m_minNintt = nintt; }
  void set_track_ntp_diagnostics(bool enable = true)
  {
    m_trackNtpDiagnostics = enable;
  }
  void set_min_vertex_tracks(unsigned int ntracks) { m_minVertexTracks = ntracks; }
  void set_vertex_match_cut_cm(float xcut, float ycut, float zcut)
  {
    m_vertexMatchCutCm = {xcut, ycut, zcut};
  }
  void set_max_tracks_above_min_pt_per_vertex(unsigned int ntracks)
  {
    m_maxTracksAboveMinPtPerVertex = ntracks;
  }

  void set_ignore_tpc() {m_ignore_tpc = true;}

  void set_layers_fixed(unsigned int minlayer, unsigned int maxlayer);
  void set_error_inflation_factor(unsigned int layer, float factor)
  {
    m_layerMisalignment.insert(std::make_pair(layer, factor));
  }
  void set_tpc_sector_fixed(unsigned int region, unsigned int sector,
                            unsigned int side)
  {
    // make a combined subsector index
    unsigned int subsector = region * 24 + side * 12 + sector;
    fixed_sectors.insert(subsector);
  }
  void set_vtx_sigma(float xysig, float zsig)
  {
    m_vtxSigma(0) = xysig;
    m_vtxSigma(1) = zsig;
  }
  void set_vtx_sigma_inflation(float factor) { m_vtxSigmaInflation = factor; }

 private:
  Mille* _mille;

  int GetNodes(PHCompositeNode* topNode);
  Acts::Vector3 getEventVertex(const SvtxTrack* track, const SvtxVertex*& svtxVertexOut) const;

  bool is_layer_fixed(unsigned int layer);

  bool is_layer_param_fixed(unsigned int layer, unsigned int param, std::set<std::pair<unsigned int, unsigned int>>& param_fixed);

  bool is_tpc_sector_fixed(unsigned int layer, unsigned int sector, unsigned int side);
  bool is_mvtx_layer_fixed(unsigned int layer, unsigned int clamshell);
  void addTrackToMilleFile(SvtxAlignmentStateMap::StateVec& statevec);
  void getGlobalVtxDerivativesXY(SvtxTrack* track,
                                 const Acts::Vector3& vertex,
                                 float glblvtx_derivative[SvtxAlignmentState::NRES][3]);
  bool getLocalVtxDerivativesXY(SvtxTrack* track,
                                ActsPropagator& propagator,
                                const Acts::Vector3& vertex,
                                float lclvtx_derivative[SvtxAlignmentState::NRES][SvtxAlignmentState::NLOC],
                                Acts::BoundSquareMatrix* perigeeCov = nullptr);
  bool getLocalVtxDerivativesXYpsuedo(SvtxTrack* track,
                                      const SvtxVertex* vertex,
                                      const Acts::BoundSquareMatrix& perigeeCov,
                                      float localderpsuedo[SvtxAlignmentState::NLOC][SvtxAlignmentState::NLOC],
                                      float localmeaserrpsuedo[SvtxAlignmentState::NLOC]);
  Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC> regulariseCovariance(const Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>& inputCov,
                                                 double conditionCutOff,
                                                 double removeHugeLeading,
                                                 double stabilisationDiag);
  Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC> getInverseComplement(Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>& target,
                                           const Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>& existing_sol);
  Acts::Vector3 localToGlobalVertex(SvtxTrack* track,
                                    const Acts::Vector3& vertex,
                                    const Acts::Vector3& localx) const;
  void getProjectionVtxXY(SvtxTrack* track, const Acts::Vector3& vertex,
                          Acts::Vector3& projx, Acts::Vector3& projy);

  std::string data_outfilename = ("mille_output_data_file.bin");
  std::string steering_outfilename = ("steer.txt");

  bool m_useEventVertex = true;
  bool _binary = true;
  float m_minPt = 0.0;
  unsigned int m_minNmvtx = 0;
  unsigned int m_minNintt = 0;
  bool m_trackNtpDiagnostics = false;
  unsigned int m_minVertexTracks = 0;
  unsigned int m_maxTracksAboveMinPtPerVertex =
      std::numeric_limits<unsigned int>::max();
  Acts::Vector3 m_vertexMatchCutCm = {0.2, 0.2, 0.2};

  Acts::Vector2 m_vtxSigma = {0.1, 0.1};
  float m_vtxSigmaInflation = 1.0;

  std::map<unsigned int, float> m_layerMisalignment;
  std::set<unsigned int> fixed_sectors;
  // set default groups to lowest level
  AlignmentDefs::mvtxGrp mvtx_group = AlignmentDefs::mvtxGrp::snsr;
  AlignmentDefs::inttGrp intt_group = AlignmentDefs::inttGrp::chp;
  AlignmentDefs::tpcGrp tpc_group = AlignmentDefs::tpcGrp::htst;
  AlignmentDefs::mmsGrp mms_group = AlignmentDefs::mmsGrp::tl;

  std::set<unsigned int> fixed_layers;
  std::set<std::pair<unsigned int, unsigned int>> fixed_mvtx_layers;
  std::set<std::pair<unsigned int, unsigned int>> fixed_layer_gparams, fixed_layer_lparams;

  bool m_ignore_tpc = false;
  
  std::string m_constraintFileName = "mp2con.txt";
  std::ofstream m_constraintFile;

  std::string m_track_map_name{"SvtxTrackMap"};
  SvtxTrackMap* _track_map{nullptr};

  SvtxVertexMap* m_vertex_map{nullptr};

  std::string m_state_map_name{"SvtxAlignmentStateMap"};
  SvtxAlignmentStateMap* _state_map{nullptr};

  ActsGeometry* _tGeometry{nullptr};
  TrkrClusterContainer* _cluster_map{nullptr};
  ClusterErrorPara _ClusErrPara;

  std::string m_tfile_name;
  TFile* m_file{nullptr};
  TH1D* m_trackCutFlow{nullptr};
  TH1D* h_track_pt_after_pruning{nullptr};
  TH1D* h_track_pt_after_makemille_pt_cut{nullptr};
  TH1D* h_vertex_ntracks{nullptr};
  TH1D* h_vertex_ntracks_raw{nullptr};
  TH1D* h_absent_vertex_id{nullptr};
  TH1D* h_track_crossing_with_vertex{nullptr};
  TH1D* h_track_crossing_without_vertex{nullptr};
  TH1D* h_track_vertex_crossing_diff{nullptr};
  TH2D* h_frac_no_vertex_vs_ntracks{nullptr};
  TProfile* p_frac_no_vertex_vs_ntracks{nullptr};
  TH1D* h_ntracks_nocrossing{nullptr};
  TH2D* h_frac_no_vertex_vs_ntracks_nocrossing{nullptr};
  TProfile* p_frac_no_vertex_vs_ntracks_nocrossing{nullptr};
  TH2D* h_vertex_ntracks_vs_ntracks_crossing{nullptr};
  TProfile* p_vertex_ntracks_vs_ntracks_crossing{nullptr};
  TH2D* h_vertex_ntracks_vs_ntracks_crossing_used{nullptr};
  TProfile* p_vertex_ntracks_vs_ntracks_crossing_used{nullptr};
  TH2D* h_vertex_ntracks_vs_ntracks_crossing_attached{nullptr};
  TProfile* p_vertex_ntracks_vs_ntracks_crossing_attached{nullptr};
  TNtuple* m_ntuple{nullptr};
  TNtuple* track_ntp{nullptr};
};

#endif  // MAKEMILLEFILES_H
