#pragma once

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <array>
#include <string>
#include <vector>

class Tpc_AssembledTrackContainer;
class IdealPadMap;
class PHCompositeNode;
class PHGarfield;
class Tpc_PolyClusterContainer;
class TrkrHitSetContainer;

class Tpc_PolyClusterizer : public SubsysReco
{
 public:
  explicit Tpc_PolyClusterizer(const std::string& name = "Tpc_PolyClusterizer");
  ~Tpc_PolyClusterizer() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  static constexpr unsigned int NPhiSamples = 3;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setT0(double v) { m_t0 = v; }
  void setTpcAdcClock(double v) { m_tpcAdcClock = v; }
  void setReverseDriftStepNs(double v) { m_reverseDriftStepNs = v; }
  void setKEffSide0(double v) { m_kEffSide0 = v; }
  void setKEffSide1(double v) { m_kEffSide1 = v; }
  void setCMVoltageDefault(double v) { m_cmVoltageDefault = v; }
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
  struct Point
  {
    TrkrDefs::hitsetkey hitsetkey{0};
    TrkrDefs::hitkey hitkey{0};
    unsigned int layer{0};
    unsigned int side{0};
    unsigned int pad{0};
    unsigned int tbin{0};
    double adc{0.0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
  };

  struct Centroid
  {
    bool ok{false};
    unsigned int layer{0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double rms_x{0.0};
    double rms_y{0.0};
    double rms_z{0.0};
  };

  struct ClusterParameters
  {
    double adc{0.0};
    unsigned int phi_width{0};
    unsigned int time_width{0};
    double phase{0.0};
  };

  struct DriftPoint
  {
    float delta_r{0.0F};
    float delta_phi{0.0F};
    float z{0.0F};
  };

  struct DriftPolyline
  {
    double phi{0.0};
    std::vector<DriftPoint> points;
  };

  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  bool make_xyz_point(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, Point& p) const;
  bool build_drift_lookup();
  bool sample_drift_lookup(unsigned int layer,
                           unsigned int side,
                           unsigned int pad,
                           unsigned int tbin,
                           double& x,
                           double& y,
                           double& z) const;
  ClusterParameters make_cluster_parameters(const std::vector<Point>& points, const Centroid& centroid, int side) const;
  static Centroid make_centroid(const std::vector<Point>& points);
  void configure_garfield(PHGarfield* garfield) const;
  static unsigned int drift_lookup_index(unsigned int layer_index, unsigned int side, unsigned int sector, unsigned int sample);
  std::string m_inputNodeName;
  std::string m_outputNodeName;
  Tpc_AssembledTrackContainer* m_assembledTracks{nullptr};
  Tpc_PolyClusterContainer* m_clusters{nullptr};
  TrkrHitSetContainer* m_hits{nullptr};
  IdealPadMap* m_idealPadMap{nullptr};
  PHGarfield* m_garfield{nullptr};
  std::array<DriftPolyline, 48 * 2 * 12 * NPhiSamples> m_driftLookup;
  unsigned int m_event{0};
  double m_t0{8};
  double m_tpcAdcClock{56.881262};
  double m_reverseDriftStepNs{56.881262};
  double m_startZSouth{-102.325};
  double m_startZNorth{102.325};
  double m_kEffSide0{0.0};
  double m_kEffSide1{-1.5};
  double m_cmVoltageDefault{380.0};
  std::array<double, 3> m_tpcMove{{0.0, 0.0, 0.0}};                                             //{{-0.16775, -0.0337, -0.71365}};
  std::array<std::array<double, 3>, 2> m_tpcRotations{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};  //{{{{0.0, 0.01485 / 10.0, 0.0}}, {{0.0298 / 8.0, 0.0, 0.0}}}};
};
