// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4TPC_PHG4TpcPadBaselineShift_H
#define G4TPC_PHG4TpcPadBaselineShift_H

#include <fun4all/SubsysReco.h>

#include <climits>
#include <cmath>
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
  //int ResetEvent(PHCompositeNode *topNode) override;

  //int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  //int Reset(PHCompositeNode * /*topNode*/) override;

  //void Print(const std::string &what = "ALL") const override;

  void setScale(float CScale);
  void setFileName(const std::string &filename);
  void writeTree(int f_writeTree);
  void set_drift_velocity(float vd) {_drift_velocity = vd;}

 private:
  bool is_in_sector_boundary(int phibin, int sector, PHG4TpcCylinderGeom *layergeom);
  float _hit_z = NAN;
  float _hit_r = NAN;
  float _hit_phi = NAN;
  float _hit_e = NAN;
  int _hit_adc = INT_MIN;
  int _hit_adc_bls = INT_MIN;
  int _hit_layer = INT_MIN;
  int _hit_sector = INT_MIN;

  TrkrHitSetContainer *m_hits = nullptr;
  TrkrClusterContainer *m_clusterlist = nullptr;
  TrkrClusterHitAssoc *m_clusterhitassoc = nullptr;
  ActsSurfaceMaps *m_surfMaps = nullptr;
  ActsTrackingGeometry *m_tGeometry = nullptr;

  //   bool do_hit_assoc = true;
  //   double pedestal = 74.4;
  int _writeTree = 0;
  double SectorFiducialCut = 0.5;

  //   int NSearch = 2;
  int NZBinsMax = 0;
  float _CScale = 1.;

  double AdcClockPeriod = 53.0;  // ns
  unsigned int MaxTBins = 498;
  float _drift_velocity = 8.0e-03;
 
  TFile *outfile = nullptr;
  std::string _filename = "./hitsBLS.root";

  TTree *_rawHits = nullptr;
};

#endif  // PHG4TpcPadBaselineShift_H
