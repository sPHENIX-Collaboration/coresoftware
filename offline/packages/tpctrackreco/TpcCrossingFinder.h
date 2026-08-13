#ifndef TPCTRACKRECO_TPCCROSSINGFINDER_H
#define TPCTRACKRECO_TPCCROSSINGFINDER_H


#include <fun4all/SubsysReco.h>

#include "TpcCrossingDecision.h"
#include <trackbase/TrkrDefs.h>

#include <array>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

class IdealPadMap;
class PHCompositeNode;
class PHGarfield;
class PHG4TpcGeomContainer;
class SvtxVertexMap;
class Tpc_AssembledTrack;
class Tpc_AssembledTrackContainer;
class TpcCrossingDecisionContainerv1;
class TrkrClusterContainer;
class TrkrHitSetContainer;

class TpcCrossingFinder : public SubsysReco
{
 public:
  explicit TpcCrossingFinder(const std::string& name = "TpcCrossingFinder");
  ~TpcCrossingFinder() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  static constexpr unsigned int NPhiSamples = 3;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setVertexMapNodeName(const std::string& n) { m_vertexMapNodeName = n; }
  void setT0(double v) { m_t0 = v; }
  void setTpcAdcClock(double v) { m_tpcAdcClock = v; }
  void setCrossingPeriodNs(double v) { m_crossingPeriodNs = v; }
  void setReverseDriftStepNs(double v) { m_reverseDriftStepNs = v; }
  void setKEffSide0(double v) { m_kEffSide0 = v; }
  void setKEffSide1(double v) { m_kEffSide1 = v; }
  void setCMVoltageDefault(double v) { m_cmVoltageDefault = v; }
  void setRequireSiliconVertex(bool v) { m_requireSiliconVertex = v; }
  void setResolveAmbiguousWithoutVertex(bool v) { m_resolveAmbiguousWithoutVertex = v; }
  void setPreferTriggeredCrossing(bool v) { m_preferTriggeredCrossing = v; }
  void setTriggeredCrossing(short v) { m_triggeredCrossing = v; }
  void setCollisionZ(double v) { m_collisionZ = v; }
  void setMaxVertexDz(double v) { m_maxVertexDz = v; }
  void setMaxTier2BeamlineZ(double v) { m_maxTier2BeamlineZ = v; }
  void setMaxCandidateVertexZ(double v) { m_maxCandidateVertexZ = v; }
  void setMinBestSecondSeparation(double v) { m_minBestSecondSeparation = v; }
  void setTpcGeometryTolerance(double radial_cm, double z_cm, double central_membrane_cm)
  {
    m_radialTolerance = radial_cm;
    m_zTolerance = z_cm;
    m_centralMembraneTolerance = central_membrane_cm;
  }
  void setTpcHalfLength(double v) { m_tpcHalfLength = v; }
  void setUseSurveyGeometry(bool v) { use_survey_geometry = v; }
  void setMoveTpc(double x, double y, double z) { m_tpcMove = {{x, y, z}}; }
  void setRotateTpc(unsigned int index, double x, double y, double z)
  {
    if (index < m_tpcRotations.size()) m_tpcRotations[index] = {{x, y, z}};
  }
  void setStartZ(double south_z, double north_z)
  {
    m_startZSouth = south_z;
    m_startZNorth = north_z;
  }

 private:
  struct DriftPoint
  {
    float delta_r {0.0F};
    float delta_phi {0.0F};
    float z {0.0F};
  };

  struct DriftPolyline
  {
    double phi {0.0};
    std::vector<DriftPoint> points;
  };

  struct Point
  {
    TrkrDefs::hitsetkey hitsetkey {0};
    TrkrDefs::hitkey hitkey {0};
    unsigned int layer {0};
    unsigned int side {0};
    unsigned int pad {0};
    unsigned int tbin {0};
    double x {0.0};
    double y {0.0};
    double z {0.0};
  };

  struct SiliconVertexHypothesis
  {
    short crossing {0};
    unsigned int vertex_id {0};
    double x {0.0};
    double y {0.0};
    double z {0.0};
    double sigma_z {0.0};
    unsigned int ntracks {0};
  };

  struct ZFitResult
  {
    bool valid {false};
    double slope {0.0};
    double intercept {0.0};
    double chi2 {0.0};
    int ndf {-1};
    double s_at_pca {0.0};
    double minimum_radius {0.0};
    double z_at_pca {0.0};
    double z_at_r0 {0.0};
    std::vector<Point> points;
    std::vector<float> path_length;
  };

  struct Candidate
  {
    short crossing {0};
    bool tpc_valid {false};
    bool has_silicon_vertex {false};
    bool vertex_compatible {false};
    unsigned int silicon_vertex_id {0};
    double tpc_z0 {0.0};
    double silicon_vertex_z {0.0};
    double delta_z {0.0};
    unsigned char rejection_status {0};
    unsigned char confidence_tier {std::numeric_limits<unsigned char>::max()};
    double confidence_score {std::numeric_limits<double>::quiet_NaN()};
    TpcCrossingCandidate qa;
  };

  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  void configure_garfield(PHGarfield* garfield) const;
  bool build_drift_lookup();
  bool sample_drift_lookup(unsigned int layer,
                           unsigned int side,
                           unsigned int pad,
                           unsigned int tbin,
                           short crossing,
                           double& x,
                           double& y,
                           double& z) const;
  bool make_xyz_point(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, short crossing, Point& p) const;
  bool find_time_extrema(const Tpc_AssembledTrack* track,
                         TrkrDefs::hitsetkey& min_hsk,
                         TrkrDefs::hitkey& min_hk,
                         TrkrDefs::hitsetkey& max_hsk,
                         TrkrDefs::hitkey& max_hk) const;
  std::set<short> get_available_crossings() const;
  std::set<short> get_intt_crossings() const;
  std::map<short, std::vector<SiliconVertexHypothesis>> get_vertices_by_crossing() const;
  std::vector<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>> select_representatives(const Tpc_AssembledTrack* track,
                                                                   TrkrDefs::hitsetkey min_hsk,
                                                                   TrkrDefs::hitkey min_hk,
                                                                   TrkrDefs::hitsetkey max_hsk,
                                                                   TrkrDefs::hitkey max_hk) const;
  bool estimate_tpc_z0(std::vector<Point>& points, double& z0) const;
  ZFitResult estimate_tpc_z0_diagnostics(std::vector<Point> points) const;
  bool point_in_tpc(const Point& p) const;
  bool point_in_correct_side(const Point& p) const;
  Candidate test_candidate(const Tpc_AssembledTrack* track,
                           short crossing,
                           TrkrDefs::hitsetkey min_hsk,
                           TrkrDefs::hitkey min_hk,
                           TrkrDefs::hitsetkey max_hsk,
                           TrkrDefs::hitkey max_hk,
                           const std::map<short, std::vector<SiliconVertexHypothesis>>& vertices_by_crossing) const;
  static unsigned int drift_lookup_index(unsigned int layer_index, unsigned int side, unsigned int sector, unsigned int sample);

  std::string m_inputNodeName {"TPC_ASSEMBLEDTRACKS"};
  std::string m_outputNodeName {"TPC_CROSSING_DECISIONS"};
  std::string m_vertexMapNodeName {"SvtxVertexMap"};

  Tpc_AssembledTrackContainer* m_assembledTracks {nullptr};
  TpcCrossingDecisionContainerv1* m_decisions {nullptr};
  TrkrHitSetContainer* m_hits {nullptr};
  TrkrClusterContainer* m_clusterMap {nullptr};
  PHG4TpcGeomContainer* m_geomContainerTpc {nullptr};
  SvtxVertexMap* m_vertexMap {nullptr};
  IdealPadMap* m_idealPadMap {nullptr};
  PHGarfield* m_garfield {nullptr};

  std::array<DriftPolyline, 48 * 2 * 12 * NPhiSamples> m_driftLookup;
  unsigned int m_event {0};
  double m_maxLookupTimeNs {0.0};

  double m_t0 {8};
  double m_tpcAdcClock {56.881262};
  double m_crossingPeriodNs {106.56};
  double m_reverseDriftStepNs {56.881262};
  double m_startZSouth {-102.325};
  double m_startZNorth {102.325};
  double m_tpcHalfLength {105.5};
  double m_radialTolerance {2.0};
  double m_zTolerance {1.0};
  double m_centralMembraneTolerance {1.0};
  double m_cmVoltageDefault {380.0};
  double m_kEffSide0 {0.0};
  double m_kEffSide1 {-1.5};
  double m_maxVertexDz {2.0};
  double m_minBestSecondSeparation {0.4};
  double m_collisionZ {0.0};
  double m_maxTier2BeamlineZ {40.0};
  double m_maxCandidateVertexZ {20.0};
  bool m_requireSiliconVertex {false};
  bool m_resolveAmbiguousWithoutVertex {true};
  bool m_preferTriggeredCrossing {false};
  short m_triggeredCrossing {0};
  bool use_survey_geometry {false};
  std::array<double, 3> m_tpcMove {{0.0, 0.0, 0.0}};
  std::array<std::array<double, 3>, 2> m_tpcRotations {{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
};

#endif  // TPCTRACKRECO_TPCCROSSINGFINDER_H
