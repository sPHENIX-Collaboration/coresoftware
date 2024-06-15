#ifndef TRACKRECO_PHGENFITTRKFITTER_H
#define TRACKRECO_PHGENFITTRKFITTER_H

/*!
 *  \file		PHGenFitTrkFitter.h
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include <fun4all/SubsysReco.h>
#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>
#include <trackbase_historic/ActsTransformations.h>

#if defined(__CLING__)
// needed, it crashes on Ubuntu using singularity with local cvmfs install
// shared pointer later on uses this, forward declaration does not cut it
#include <phgenfit/Track.h>
#else
namespace PHGenFit
{
  class Track;
} /* namespace PHGenFit */
#endif

#include <TMatrixFfwd.h>  // for TMatrixF
#include <TVector3.h>     // for TVector3

#include <map>
#include <memory>  // for shared_ptr
#include <set>
#include <string>
#include <vector>

class TClonesArray;

namespace genfit
{
  class GFRaveVertex;
  class GFRaveVertexFactory;
  class Track;
} /* namespace genfit */

class SvtxTrack;
namespace PHGenFit
{
  class Fitter;
} /* namespace PHGenFit */

class ActsGeometry;
class PHCompositeNode;
class SvtxTrackMap;
class TpcDistortionCorrectionContainer;
class TrkrClusterContainer;
class TrackSeedContainer;

//! \brief		Refit SvtxTracks with PHGenFit.
class PHGenFitTrkFitter : public SubsysReco
{
 public:

  //! Default constructor
  PHGenFitTrkFitter(const std::string& name = "PHGenFitTrkFitter");

  //! Initialization, called for initialization
  int Init(PHCompositeNode*) override;

  //! Initialization Run, called for initialization of a run
  int InitRun(PHCompositeNode*) override;

  //! Process Event, called for each event
  int process_event(PHCompositeNode*) override;

  //! End, write and close files
  int End(PHCompositeNode*) override;

  const std::string& get_track_fitting_alg_name() const
  {
    return _track_fitting_alg_name;
  }

  void set_track_fitting_alg_name(const std::string& trackFittingAlgName)
  {
    _track_fitting_alg_name = trackFittingAlgName;
  }

  int get_primary_pid_guess() const
  {
    return _primary_pid_guess;
  }

  void set_primary_pid_guess(int primaryPidGuess)
  {
    _primary_pid_guess = primaryPidGuess;
  }

  double get_fit_min_pT() const
  {
    return _fit_min_pT;
  }

  void set_fit_min_pT(double cutMinPT)
  {
    _fit_min_pT = cutMinPT;
  }

  double get_vertex_min_ndf() const
  {
    return _vertex_min_ndf;
  }

  void set_vertex_min_ndf(double vertexMinPT)
  {
    _vertex_min_ndf = vertexMinPT;
  }

  //!@name disabled layers interface
  //@{

  //! mark layer as disbled
  void disable_layer(int layer, bool disabled = true);

  //! set disabled layers
  void set_disabled_layers(const std::set<int>&);

  //! clear disabled layers
  void clear_disabled_layers();

  //! get disabled layers
  const std::set<int>& get_disabled_layers() const;

  //@}

  /// fit only tracks with silicon+MM hits
  void set_fit_silicon_mms(bool);

  /// require micromegas in SiliconMM fits
  void set_use_micromegas(bool value)
  {
    m_use_micromegas = value;
  }

 private:
  //! Event counter
  int _event = 0;

  //! Get all the nodes
  int GetNodes(PHCompositeNode*);

  //! Create New nodes
  int CreateNodes(PHCompositeNode*);

  /// get global position for a given cluster
  /**
   * uses ActsTransformation to convert cluster local position into global coordinates
   */
  Acts::Vector3 getGlobalPosition(TrkrDefs::cluskey, TrkrCluster*, short int crossing);

  /*
   * fit track with SvtxTrack as input seed.
   * \param intrack Input SvtxTrack
   * \param invertex Input Vertex, if fit track as a primary vertex
   */
  std::shared_ptr<PHGenFit::Track> ReFitTrack(PHCompositeNode*, const SvtxTrack* intrack);

  //! Make SvtxTrack from PHGenFit::Track and SvtxTrack
  std::shared_ptr<SvtxTrack> MakeSvtxTrack(const SvtxTrack* svtxtrack, const std::shared_ptr<PHGenFit::Track>& genfit_track );

  bool pos_cov_XYZ_to_RZ(
      const TVector3& n,
      const TMatrixF& pos_in,
      const TMatrixF& cov_in,
      TMatrixF& pos_out,
      TMatrixF& cov_out) const;

  //! disabled layers
  /** clusters belonging to disabled layers are not included in track fit */
  std::set<int> _disabled_layers;

  /// Boolean to use normal tracking geometry navigator or the
  /// Acts::DirectedNavigator with a list of sorted silicon+MM surfaces
  bool m_fit_silicon_mms = false;

  /// requires micromegas present when fitting silicon-MM surfaces
  bool m_use_micromegas = true;

  //! KalmanFitterRefTrack, KalmanFitter, DafSimple, DafRef
  std::string _track_fitting_alg_name = "DafRef";

  int _primary_pid_guess = 211;
  double _fit_min_pT = 0.1;
  double _vertex_min_ndf = 20;

  /*
  need to use shared_ptr and not unique_ptr because root5 cint
  requires the existence of a copy constructor, which the unique_ptr forbids
  */
  std::shared_ptr<PHGenFit::Fitter> _fitter;

  /// acts geometry
  ActsGeometry* m_tgeometry = nullptr;

  //! Input Node pointers
  TrkrClusterContainer* m_clustermap = nullptr;

  // track seeds
  TrackSeedContainer* m_seedMap = nullptr;
  TrackSeedContainer* m_tpcSeeds = nullptr;
  TrackSeedContainer* m_siliconSeeds = nullptr;

  //! Output Node pointers
  SvtxTrackMap* m_trackMap = nullptr;

  // crossing z correction
  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;

  // distortion corrections
  TpcDistortionCorrectionContainer* m_dcc_static = nullptr;
  TpcDistortionCorrectionContainer* m_dcc_average = nullptr;
  TpcDistortionCorrectionContainer* m_dcc_fluctuation = nullptr;

  /// tpc distortion correction utility class
  TpcDistortionCorrection m_distortionCorrection;

};

#endif
