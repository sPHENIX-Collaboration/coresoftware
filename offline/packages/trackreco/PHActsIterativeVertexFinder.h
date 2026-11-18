// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKRECO_PHACTSITERATIVEVERTEXFINDER_H
#define TRACKRECO_PHACTSITERATIVEVERTEXFINDER_H

/*!
 *  \file   PHActsIterativeVertexFinder.h
 *  \brief  Primary vertex finding with the Acts IterativeVertexFinder,
 *          run on the output of PHActsTrkFitter. Drop-in replacement for
 *          PHSimpleVertexFinder: fills the SvtxVertexMap and
 *          TrackVertexCrossingAssocMap and sets SvtxTrack::vertex_id, so that
 *          PHActsVertexPropagator can be run afterwards unchanged.
 */

#include <fun4all/SubsysReco.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/TrackParameters.hpp>

#include <limits>
#include <memory>
#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class TrackSeed;
class TrackSeedContainer;
class TrackVertexCrossingAssoc;

namespace Acts
{
  class BasePropagator;
  class FullBilloirVertexFitter;
  class HelicalTrackLinearizer;
  class IterativeVertexFinder;
  class MagneticFieldProvider;
  class Vertex;
}  // namespace Acts

class PHActsIterativeVertexFinder : public SubsysReco
{
 public:
  /// Algorithm used by the IVF to produce the seed vertex at each iteration
  enum class SeedFinder
  {
    ZScan,         ///< mode of the track z0 distribution (classic IVF seeder), x,y at the beam spot
    TrackDensity,  ///< maximum of the gaussian track density along z, x,y at the beam spot
    ZScanXYFit     ///< ZScan z, then x,y,z fitted from the tracks near it (ActsXYVertexSeeder)
  };

  PHActsIterativeVertexFinder(const std::string &name = "PHActsIterativeVertexFinder");
  ~PHActsIterativeVertexFinder() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  /// node names
  void setTrackMapName(const std::string &name) { m_trackMapName = name; }
  void setVertexMapName(const std::string &name) { m_vertexMapName = name; }

  /// track selection, same meaning as in PHSimpleVertexFinder
  void setTrackQualityCut(double cut) { m_qualCut = cut; }
  void setTrackPtCut(double cut) { m_trackPtCut = cut; }
  void setRequireMVTX(bool set = true) { m_requireMvtx = set; }
  void setNmvtxRequired(unsigned int n) { m_nMvtxRequired = n; }
  void setRequireINTT(bool set = true) { m_requireIntt = set; }
  void setNinttRequired(unsigned int n) { m_nInttRequired = n; }
  /// in pp mode, crossing 0 tracks without a silicon seed are skipped
  void set_pp_mode(bool mode = true) { m_ppMode = mode; }

  /// Beam spot (cm). Its x,y are always used as the transverse position of the
  /// seed vertex. Its width (cm) is only used if the constraint is enabled,
  /// in which case it is also added to the vertex fit
  void setBeamSpot(double x, double y, double z = 0.)
  {
    m_beamSpot = Acts::Vector3(x, y, z);
  }
  void setBeamSpotWidth(double sx, double sy, double sz)
  {
    m_beamSpotWidth = Acts::Vector3(sx, sy, sz);
  }
  void useBeamSpotConstraint(bool use = true) { m_useBeamSpotConstraint = use; }

  /// Acts IterativeVertexFinder configuration
  void setSeedFinder(SeedFinder seeder) { m_seedFinder = seeder; }
  void setMaxVertices(int n) { m_maxVertices = n; }
  void setSignificanceCutSeeding(double cut) { m_significanceCutSeeding = cut; }
  void setMaxChi2CutForSeeding(double cut) { m_maxChi2CutForSeeding = cut; }
  void setReassignTracksAfterFirstFit(bool set = true) { m_reassignTracksAfterFirstFit = set; }
  void setCutOffTrackWeight(double cut) { m_cutOffTrackWeight = cut; }
  /// ZScan seeder: minimum pT (GeV) used in the log(pT/minPt) track weight.
  /// Tracks below it get zero weight, so it acts as a pT cut on the z seed.
  /// The Acts default (0.4 GeV) is tuned for the LHC; in low multiplicity pp
  /// events with only soft tracks it leaves the z scan without input, the
  /// seed falls back to z = 0 and the IVF finds nothing (5+ track crossings)
  void setZScanMinPt(double pt) { m_zScanMinPt = pt; }
  /// ZScanXYFit seeder: half width (cm) of the z window used to pick the
  /// tracks for the transverse seed fit, and the per-track chi2 above which a
  /// track is dropped from that fit
  void setSeedZWindow(double cm) { m_seedZWindow = cm; }
  void setSeedMaxTrackChi2(double cut) { m_seedMaxTrackChi2 = cut; }
  /// minimum number of tracks attached to a vertex to keep it
  void setMinTracksPerVertex(unsigned int n) { m_minTracksPerVertex = n; }
  /// reject vertices with a fit chi2/ndf above this value (default: no cut).
  /// Note the Acts Billoir fit is 4D, so ndf also counts the (trivial) time dof
  void setMaxVertexChi2Ndf(double cut) { m_maxVertexChi2Ndf = cut; }

  /// Outlier pass on every vertex returned by the IVF (default: off). The
  /// track with the largest Billoir chi2 is removed and the vertex refitted
  /// while that chi2 is above the cut. If the vertex would drop below
  /// setMinTracksPerVertex it is rejected and its tracks are left unused.
  /// Needed because the Acts IVF puts all remaining tracks into the fit
  /// without a compatibility check once 4 or fewer are left
  void setOutlierChi2Cut(double cut) { m_outlierChi2Cut = cut; }

  /// Diagnostic: inflate the track covariances given to the IVF (the stored
  /// SvtxTrack covariance is not changed). The full bound covariance is
  /// multiplied by scale^2, i.e. every sigma is multiplied by scale
  void setTrackCovarianceScale(double scale) { m_trackCovScale = scale; }
  /// Diagnostic: add a multiple scattering term in quadrature to the d0 and
  /// z0 sigmas, sigma^2 += (b / pT)^2, with b in cm*GeV (e.g. 23 um*GeV =
  /// 0.0023). Applied after setTrackCovarianceScale
  void setTrackMSTerm(double bD0, double bZ0)
  {
    m_trackMSTermD0 = bD0;
    m_trackMSTermZ0 = bZ0;
  }

  /// Assign tracks not used by any vertex to the closest vertex in z in the
  /// same crossing (as PHSimpleVertexFinder does). The track is not added to
  /// the vertex track list, only its vertex_id is set
  void setAssignUnusedTracks(bool set = true) { m_assignUnusedTracks = set; }

  /// Use a constant field of this value (Tesla) instead of the field map,
  /// e.g. "0" for zero field running. Same convention as PHActsVertexPropagator
  void fieldMap(const std::string &fieldmap) { m_fieldMap = fieldmap; }

 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  void setupFinder();

  bool passTrackSelection(SvtxTrack *track) const;
  unsigned int countSiliconClusters(SvtxTrack *track, uint8_t trkrId) const;
  TrackSeed *getSiliconSeed(SvtxTrack *track) const;

  /// Run the finder on the tracks of one crossing. Returns the number of
  /// vertices written, with ids starting at vertexId
  unsigned int processCrossing(short int crossing, unsigned int vertexId);
  void assignUnusedTracks(short int crossing,
                          const std::vector<unsigned int> &vertexIds);

  /// nodes
  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;
  TrackVertexCrossingAssoc *m_trackVertexCrossingMap = nullptr;
  TrackSeedContainer *m_siliconSeedMap = nullptr;
  ActsGeometry *m_tGeometry = nullptr;

  std::string m_trackMapName = "SvtxTrackMap";
  std::string m_vertexMapName = "SvtxVertexMap";

  /// track selection
  double m_qualCut = 10.;
  double m_trackPtCut = 0.;
  bool m_requireMvtx = true;
  unsigned int m_nMvtxRequired = 2;
  bool m_requireIntt = false;
  unsigned int m_nInttRequired = 1;
  bool m_ppMode = true;

  /// beam spot, cm
  Acts::Vector3 m_beamSpot = Acts::Vector3::Zero();
  Acts::Vector3 m_beamSpotWidth = Acts::Vector3(0.01, 0.01, 10.);
  bool m_useBeamSpotConstraint = false;

  /// finder config
  SeedFinder m_seedFinder = SeedFinder::ZScan;
  int m_maxVertices = 50;
  double m_significanceCutSeeding = 10.;
  double m_maxChi2CutForSeeding = 36.;
  bool m_reassignTracksAfterFirstFit = false;
  double m_cutOffTrackWeight = 0.01;
  double m_zScanMinPt = 0.1;
  double m_seedZWindow = 0.2;
  double m_seedMaxTrackChi2 = 36.;
  unsigned int m_minTracksPerVertex = 2;
  double m_maxVertexChi2Ndf = std::numeric_limits<double>::max();
  double m_outlierChi2Cut = std::numeric_limits<double>::max();
  double m_trackCovScale = 1.;
  double m_trackMSTermD0 = 0.;
  double m_trackMSTermZ0 = 0.;
  bool m_assignUnusedTracks = true;
  std::string m_fieldMap;

  /// Acts tools, built once in InitRun. The linearizer must outlive the
  /// finder because the finder holds a delegate to it
  std::shared_ptr<const Acts::MagneticFieldProvider> m_field;
  std::shared_ptr<const Acts::BasePropagator> m_propagator;
  std::unique_ptr<Acts::HelicalTrackLinearizer> m_linearizer;
  std::unique_ptr<Acts::IterativeVertexFinder> m_finder;
  /// same configuration as the fitter inside the finder, used for refits
  std::unique_ptr<Acts::FullBilloirVertexFitter> m_fitter;

  /// bookkeeping
  unsigned long m_nEvents = 0;
  unsigned long m_nFinderCalls = 0;
  unsigned long m_nFinderFailures = 0;
  unsigned long m_nVertices = 0;
  unsigned long m_nOutliersRemoved = 0;
  unsigned long m_nRefits = 0;
  unsigned long m_nRefitFailures = 0;
  unsigned long m_nVerticesRejectedOutliers = 0;
};

#endif  // TRACKRECO_PHACTSITERATIVEVERTEXFINDER_H
