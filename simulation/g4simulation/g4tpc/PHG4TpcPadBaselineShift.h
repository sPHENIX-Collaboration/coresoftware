// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4TPC_PHG4TpcPadBaselineShift_H
#define G4TPC_PHG4TpcPadBaselineShift_H

#include <fun4all/SubsysReco.h>

#include <climits>
#include <cmath>
#include <limits>
#include <string>

class PHCompositeNode;

class TTree;
class TFile;

class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class PHG4TpcCylinderGeom;

struct ActsSurfaceMaps;
struct ActsTrackingGeometry;

class PHG4TpcPadBaselineShift : public SubsysReco
{
 public:
  PHG4TpcPadBaselineShift(const std::string &name = "PHG4TpcPadBaselineShift");

  virtual ~PHG4TpcPadBaselineShift();
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void setScale(float CScale);
  void setFileName(const std::string &filename);
  void writeTree(int f_writeTree);
  void set_drift_velocity(float vd) { _drift_velocity = vd; }

 private:
  bool is_in_sector_boundary(int phibin, int sector, PHG4TpcCylinderGeom *layergeom);

  TrkrHitSetContainer *m_hits{nullptr};
  TrkrClusterContainer *m_clusterlist{nullptr};
  TrkrClusterHitAssoc *m_clusterhitassoc{nullptr};
  ActsSurfaceMaps *m_surfMaps{nullptr};
  ActsTrackingGeometry *m_tGeometry{nullptr};
  TFile *outfile{nullptr};
  TTree *_rawHits{nullptr};

  //   double pedestal {74.4};
  double SectorFiducialCut{0.5};
  double AdcClockPeriod{53.0};  // ns

  unsigned int MaxTBins{498};

  int _hit_adc{std::numeric_limits<int>::min()};
  int _hit_adc_bls{std::numeric_limits<int>::min()};
  int _hit_layer{std::numeric_limits<int>::min()};
  int _hit_sector{std::numeric_limits<int>::min()};
  int _writeTree{0};
  int NZBinsMax{0};
  //   int NSearch {2};

  float _hit_z{std::numeric_limits<float>::signaling_NaN()};
  float _hit_r{std::numeric_limits<float>::signaling_NaN()};
  float _hit_phi{std::numeric_limits<float>::signaling_NaN()};
  float _hit_e{std::numeric_limits<float>::signaling_NaN()};
  float _CScale{1.};
  float _drift_velocity{8.0e-03};

  //   bool do_hit_assoc {true};

  std::string _filename{"./hitsBLS.root"};
};

#endif  // PHG4TpcPadBaselineShift_H
