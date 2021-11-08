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

// #include <phgenfit/Track.h> is needed, it crashes on Ubuntu using
// singularity with local cvmfs install
// shared pointer later on uses this, forward declaration does not cut it
#include <phgenfit/Track.h>

#include <gsl/gsl_rng.h>

#include <climits>  // for UINT_MAX
#include <map>
#include <string>
#include <vector>
#include <memory>

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHParameters;

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
} /* namespace genfit */

class PHG4TrackFastSim : public SubsysReco
{
 public:
  enum DETECTOR_TYPE
  {
    Vertical_Plane = 0,
    Cylinder = 1
  };

  //! Default constructor
  explicit PHG4TrackFastSim(const std::string& name = "PHG4TrackFastSim");

  //! dtor
  ~PHG4TrackFastSim() override;

  //!Initialization Run, called for initialization of a run
  int InitRun(PHCompositeNode*) override;

  //!Process Event, called for each event
  int process_event(PHCompositeNode*) override;

  //!End, write and close files
  int End(PHCompositeNode*) override;

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
    return m_FitAlgoName;
  }

  void set_fit_alg_name(const std::string& fitAlgName)
  {
    m_FitAlgoName = fitAlgName;
  }

  const std::vector<std::string>& get_phg4hits_names() const
  {
    return m_PHG4HitsNames;
  }

  //! adding hits from a PHG4Hit node, which usually belong to one detector or a sub group of detectors
  //! Orders of adding detectors do not matter as the hits are internally assembled in the time order
  //! \param[in] phg4hitsNames node name such as "G4HIT_SVTX"
  //! \param[in] PHG4TrackFastSim::Vertical_Plane or PHG4TrackFastSim::Cylinder
  //! \param[in] radres radial resolution [cm], not used for PHG4TrackFastSim::Cylinder
  //! \param[in] phires azimuthal resolution [cm]
  //! \param[in] lonres z-resolution [cm], not used for PHG4TrackFastSim::Vertical_Plane
  //! \param[in] eff    Efficiency [0-1] for a existing hit to be included in the tracking
  //! \param[in] noise  Noise hit propability [0-1]
  void add_phg4hits(
      const std::string& phg4hitsNames,
      const DETECTOR_TYPE phg4dettype,
      const float radres,
      const float phires,
      const float lonres,
      const float eff,
      const float noise)
  {
    m_PHG4HitsNames.push_back(phg4hitsNames);
    m_phg4_detector_type.push_back(phg4dettype);
    m_phg4_detector_radres.push_back(radres);
    m_phg4_detector_phires.push_back(phires);
    m_phg4_detector_lonres.push_back(lonres);
    m_phg4_detector_hitfindeff.push_back(eff);
    m_phg4_detector_noise.push_back(noise);
  }

  // legacy interface for Babar calorimeter projections
  void add_state_name(const std::string& stateName);

  // add saving of state at plane in z
  void add_zplane_state(const std::string& stateName, const double zplane);

  void add_cylinder_state(const std::string& stateName, const double radius);

  const std::string& get_trackmap_out_name() const
  {
    return m_TrackmapOutNodeName;
  }

  void set_trackmap_out_name(const std::string& trackmapOutName)
  {
    m_TrackmapOutNodeName = trackmapOutName;
  }

  const std::string& get_sub_top_node_name() const
  {
    return m_SubTopnodeName;
  }

  void set_sub_top_node_name(const std::string& subTopNodeName)
  {
    m_SubTopnodeName = subTopNodeName;
  }

  bool is_use_vertex_in_fitting() const
  {
    return m_UseVertexInFittingFlag;
  }

  void set_use_vertex_in_fitting(bool useVertexInFitting)
  {
    m_UseVertexInFittingFlag = useVertexInFitting;
  }

  double get_vertex_xy_resolution() const
  {
    return m_VertexXYResolution;
  }

  void set_vertex_xy_resolution(double vertexXyResolution)
  {
    m_VertexXYResolution = vertexXyResolution;
  }

  double get_vertex_z_resolution() const
  {
    return m_VertexZResolution;
  }

  void set_vertex_z_resolution(double vertexZResolution)
  {
    m_VertexZResolution = vertexZResolution;
  }

  int get_primary_assumption_pid() const
  {
    return m_PrimaryAssumptionPid;
  }

  void set_primary_assumption_pid(int primaryAssumptionPid)
  {
    m_PrimaryAssumptionPid = primaryAssumptionPid;
  }

  void set_primary_tracking(int pTrk)
  {
    m_PrimaryTrackingFlag = pTrk;
  }

  //! https://rave.hepforge.org/trac/wiki/RaveMethods
  const std::string& get_vertexing_method() const
  {
    return m_VertexingMethod;
  }

  //! https://rave.hepforge.org/trac/wiki/RaveMethods
  void set_vertexing_method(const std::string& vertexingMethod)
  {
    m_VertexingMethod = vertexingMethod;
  }

  double get_vertex_min_ndf() const
  {
    return m_VertexMinNdf;
  }

  void set_vertex_min_ndf(double vertexMinNdf)
  {
    m_VertexMinNdf = vertexMinNdf;
  }

  void enable_vertexing(const bool& b = true)
  {
    m_DoVertexingFlag = b;
  }

  void DisplayEvent() const;

  void Smearing(const bool b) { m_SmearingFlag = b; }

 private:
  typedef std::map<const genfit::Track*, unsigned int> GenFitTrackMap;

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
                               std::vector<PHGenFit::Measurement*>& meas_out,
                               SvtxTrack* track_out,
                               TVector3& seed_pos,
                               TVector3& seed_mom,
                               TMatrixDSym& seed_cov,
                               const bool do_smearing = true);

  PHGenFit::PlanarMeasurement* PHG4HitToMeasurementVerticalPlane(const PHG4Hit* g4hit, const double phi_resolution, const double r_resolution);

  PHGenFit::PlanarMeasurement* PHG4HitToMeasurementCylinder(const PHG4Hit* g4hit, const double phi_resolution, const double z_resolution);

  PHGenFit::Measurement* VertexMeasurement(const TVector3& vtx, double dxy, double dz);

  /*!
	 * Make SvtxTrack from PHGenFit::Track
	 */
  bool MakeSvtxTrack(SvtxTrack* track_out, const PHGenFit::Track* phgf_track_in,
                     const unsigned int truth_track_id = UINT_MAX,
                     const unsigned int nmeas = 0, const TVector3& vtx = TVector3(0.0, 0.0, 0.0));

  /*
  * Fill SvtxVertexMap from GFRaveVertexes and Tracks
  */
  bool FillSvtxVertexMap(const std::vector<genfit::GFRaveVertex*>& rave_vertices,
                         const GenFitTrackMap& gf_tracks);

 protected:
  // Pointers first
  //! random generator that conform with sPHENIX standard
  gsl_rng* m_RandomGenerator;

  /*!
	 *	GenFit fitter interface
	 */
  PHGenFit::Fitter* m_Fitter;
  genfit::GFRaveVertexFactory* m_RaveVertexFactory;

  //! Input Node pointers
  PHG4TruthInfoContainer* m_TruthContainer;

  SvtxTrackMap* m_SvtxTrackMapOut;
  SvtxVertexMap* m_SvtxVertexMap;

  std::vector<PHG4HitContainer*> m_PHG4HitContainer;
  std::vector<std::string> m_PHG4HitsNames;
  std::vector<DETECTOR_TYPE> m_phg4_detector_type;
  std::vector<float> m_phg4_detector_radres;
  std::vector<float> m_phg4_detector_phires;
  std::vector<float> m_phg4_detector_lonres;
  std::vector<float> m_phg4_detector_hitfindeff;
  std::vector<float> m_phg4_detector_noise;

  //!
  std::map<std::string, std::pair<int, double>> m_ProjectionsMap;

  std::string m_SubTopnodeName;
  std::string m_TrackmapOutNodeName;
  //! https://rave.hepforge.org/trac/wiki/RaveMethods
  std::string m_VertexingMethod;

  /*!
	 * Available choices:
	 * KalmanFitter
	 * KalmanFitterRefTrack
	 * DafSimple
	 * DafRef
	 */
  std::string m_FitAlgoName;

  double m_VertexMinNdf;
  double m_VertexXYResolution;
  double m_VertexZResolution;

  //! Event counter
  int m_EventCnt;

  int m_PrimaryAssumptionPid;

  bool m_SmearingFlag;

  //!
  bool m_DoEvtDisplayFlag;

  /*!
	 * For PseudoPatternRecognition function.
	 */

  bool m_UseVertexInFittingFlag;

  //!
  int m_PrimaryTrackingFlag;

  bool m_DoVertexingFlag;

  PHParameters* m_Parameter = nullptr;
};

#endif /*__PHG4TrackFastSim_H__*/
