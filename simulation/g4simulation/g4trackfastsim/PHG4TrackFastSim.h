// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 *  \file		PHG4TrackFastSim.h
 *  \brief		Kalman Filter based on smeared truth PHG4Hit
 *  \details	Kalman Filter based on smeared truth PHG4Hit
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef G4TRACKFASTSIM_PHG4TRACKFASTSIM_H
#define G4TRACKFASTSIM_PHG4TRACKFASTSIM_H

#include <fun4all/SubsysReco.h>

#include <TMatrixDSymfwd.h>  // for TMatrixDSym
#include <TVector3.h>

// rootcint barfs with this header so we need to hide it
#if !defined(__CINT__) || defined(__CLING__)
// needed, it crashes on Ubuntu using singularity with local cvmfs install
// shared pointer later on uses this, forward declaration does not cut it
#include <gsl/gsl_rng.h>
#include <phgenfit/Track.h>
#else
namespace PHGenFit
{
class Track;
} /* namespace PHGenFit */
#endif

#include <climits>  // for UINT_MAX
#include <iostream>
#include <map>
#include <string>
#include <vector>

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class PHCompositeNode;
class PHG4TruthInfoContainer;

namespace PHGenFit
{
class Fitter;
class Measurement;
class PlanarMeasurement;
class Track;
} /* namespace PHGenFit */
namespace genfit
{
class GFRaveVertex;
class GFRaveVertexFactory;
class Track;
} /* namespace genfit */

class PHG4TrackFastSim : public SubsysReco
{
 public:
  enum DETECTOR_TYPE
  {
    Vertical_Plane,
    Cylinder
  };

  //! Default constructor
  explicit PHG4TrackFastSim(const std::string& name = "PHG4TrackFastSim");

  //! dtor
  virtual ~PHG4TrackFastSim();

  //!Initialization Run, called for initialization of a run
  int InitRun(PHCompositeNode*);

  //!Process Event, called for each event
  int process_event(PHCompositeNode*);

  //!End, write and close files
  int End(PHCompositeNode*);

  bool is_do_evt_display() const
  {
    return m_DoEvtDisplayFlag;
  }

  void set_do_evt_display(bool doEvtDisplay)
  {
    m_DoEvtDisplayFlag = doEvtDisplay;
  }

  const std::string& get_fit_alg_name() const
  {
    return _fit_alg_name;
  }

  void set_fit_alg_name(const std::string& fitAlgName)
  {
    _fit_alg_name = fitAlgName;
  }

  const std::vector<std::string>& get_phg4hits_names() const
  {
    return _phg4hits_names;
  }

  void add_phg4hits(
      const std::string& phg4hitsNames,
      const DETECTOR_TYPE phg4dettype,
      const float radres,
      const float phires,
      const float lonres,
      const float eff,
      const float noise)
  {
    _phg4hits_names.push_back(phg4hitsNames);
    _phg4_detector_type.push_back(phg4dettype);
    _phg4_detector_radres.push_back(radres);
    _phg4_detector_phires.push_back(phires);
    _phg4_detector_lonres.push_back(lonres);
    _phg4_detector_hitfindeff.push_back(eff);
    _phg4_detector_noise.push_back(noise);
  }

  void add_state_name(const std::string& stateName)
  {
    _state_names.push_back(stateName);
  }

  const std::string& get_trackmap_out_name() const
  {
    return _trackmap_out_name;
  }

  void set_trackmap_out_name(const std::string& trackmapOutName)
  {
    _trackmap_out_name = trackmapOutName;
  }

  const std::string& get_sub_top_node_name() const
  {
    return _sub_top_node_name;
  }

  void set_sub_top_node_name(const std::string& subTopNodeName)
  {
    _sub_top_node_name = subTopNodeName;
  }

  bool is_use_vertex_in_fitting() const
  {
    return _use_vertex_in_fitting;
  }

  void set_use_vertex_in_fitting(bool useVertexInFitting)
  {
    _use_vertex_in_fitting = useVertexInFitting;
  }

  double get_vertex_xy_resolution() const
  {
    return _vertex_xy_resolution;
  }

  void set_vertex_xy_resolution(double vertexXyResolution)
  {
    _vertex_xy_resolution = vertexXyResolution;
  }

  double get_vertex_z_resolution() const
  {
    return _vertex_z_resolution;
  }

  void set_vertex_z_resolution(double vertexZResolution)
  {
    _vertex_z_resolution = vertexZResolution;
  }

  int get_primary_assumption_pid() const
  {
    return _primary_assumption_pid;
  }

  void set_primary_assumption_pid(int primaryAssumptionPid)
  {
    _primary_assumption_pid = primaryAssumptionPid;
  }

  void set_primary_tracking(int pTrk)
  {
    _primary_tracking = pTrk;
  }

  //! https://rave.hepforge.org/trac/wiki/RaveMethods
  const std::string& get_vertexing_method() const
  {
    return _vertexing_method;
  }

  //! https://rave.hepforge.org/trac/wiki/RaveMethods
  void set_vertexing_method(const std::string& vertexingMethod)
  {
    _vertexing_method = vertexingMethod;
  }

  double get_vertex_min_ndf() const
  {
    return _vertex_min_ndf;
  }

  void set_vertex_min_ndf(double vertexMinPT)
  {
    _vertex_min_ndf = vertexMinPT;
  }

  void enable_vertexing(const bool& b = true)
  {
    _do_vertexing = b;
  }

  void DisplayEvent() const;

 private:
  /*!
	 * Create needed nodes.
	 */
  int CreateNodes(PHCompositeNode*);

  /*!
	 * Get all the all the required nodes off the node tree.
	 */
  int GetNodes(PHCompositeNode*);

  /*!
	 *
	 */
  int PseudoPatternRecognition(const PHG4Particle* particle,
                               std::vector<PHGenFit::Measurement*>& meas_out, TVector3& seed_pos,
                               TVector3& seed_mom, TMatrixDSym& seed_cov, const bool do_smearing = true);

  PHGenFit::PlanarMeasurement* PHG4HitToMeasurementVerticalPlane(const PHG4Hit* g4hit, const double phi_resolution, const double r_resolution);

  PHGenFit::PlanarMeasurement* PHG4HitToMeasurementCylinder(const PHG4Hit* g4hit, const double phi_resolution, const double z_resolution);

  PHGenFit::Measurement* VertexMeasurement(const TVector3& vtx, double dxy, double dz);

  /*!
	 * Make SvtxTrack from PHGenFit::Track
	 */
  SvtxTrack* MakeSvtxTrack(const PHGenFit::Track* phgf_track_in,
                           const unsigned int truth_track_id = UINT_MAX,
                           const unsigned int nmeas = 0, const TVector3& vtx = TVector3(0.0, 0.0, 0.0));

  typedef std::map<const genfit::Track*, unsigned int> GenFitTrackMap;

  /*
  * Fill SvtxVertexMap from GFRaveVertexes and Tracks
  */
  bool FillSvtxVertexMap(
      const std::vector<genfit::GFRaveVertex*>& rave_vertices,
      const GenFitTrackMap& gf_tracks);

  //! Event counter
  int _event;

  //  DETECTOR_TYPE _detector_type;  // deprecated

  //! Input Node pointers
  PHG4TruthInfoContainer* _truth_container;

  std::vector<PHG4HitContainer*> _phg4hits;
  std::vector<std::string> _phg4hits_names;
  std::vector<DETECTOR_TYPE> _phg4_detector_type;
  std::vector<float> _phg4_detector_radres;
  std::vector<float> _phg4_detector_phires;
  std::vector<float> _phg4_detector_lonres;
  std::vector<float> _phg4_detector_hitfindeff;
  std::vector<float> _phg4_detector_noise;

  //! Output Node pointers

  std::string _sub_top_node_name;
  //	std::string _clustermap_out_name;
  std::string _trackmap_out_name;

  //	SvtxClusterMap* _clustermap_out;

  SvtxTrackMap* _trackmap_out;
  SvtxVertexMap* _vertexmap;

  /*!
	 *	GenFit fitter interface
	 */
  PHGenFit::Fitter* _fitter;
  genfit::GFRaveVertexFactory* _vertex_finder;
  //! https://rave.hepforge.org/trac/wiki/RaveMethods
  std::string _vertexing_method;
  double _vertex_min_ndf;
  bool _do_vertexing;

  /*!
	 * Available choices:
	 * KalmanFitter
	 * KalmanFitterRefTrack
	 * DafSimple
	 * DafRef
	 */
  std::string _fit_alg_name;

  //!
  int _primary_assumption_pid;

  //!
  bool m_DoEvtDisplayFlag;

  /*!
	 * For PseudoPatternRecognition function.
	 */

  bool _use_vertex_in_fitting;

  double _vertex_xy_resolution;
  double _vertex_z_resolution;
  //!
  int _primary_tracking;

  //!
  std::vector<std::string> _state_names;
  std::vector<double> _state_location;

#if !defined(__CINT__) || defined(__CLING__)
  //! random generator that conform with sPHENIX standard
  gsl_rng* m_RandomGenerator;
#endif
};

#endif /*__PHG4TrackFastSim_H__*/
