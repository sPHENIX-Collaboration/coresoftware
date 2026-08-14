#ifndef TPC_LASEREVENTIDENTIFIER_H
#define TPC_LASEREVENTIDENTIFIER_H

#include <fun4all/SubsysReco.h>

class ActsGeometry;
class PHCompositeNode;
class PHG4TpcGeom;
class PHG4TpcGeomContainer;
class LaserEventInfo;
class TF1;
class TFile;
class TH1;
class TrkrHitSet;
class TrkrHitSetContainer;
class TTree;

class LaserEventIdentifier : public SubsysReco
{
 public:
  LaserEventIdentifier(const std::string &name = "LaserEventIdentifier");
  ~LaserEventIdentifier() override;  // = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_max_time_samples(int val) { m_time_samples_max = val; }

  void set_runnumber(int runnum) { m_runnumber = runnum; }

  void set_debug(bool debug) { m_debug = debug; }
  void set_debug_name(const std::string &name) { m_debugFileName = name; }

 private:
  TrkrHitSetContainer *m_hits{nullptr};
  ActsGeometry *m_tGeometry{nullptr};
  PHG4TpcGeomContainer *m_geom_container{nullptr};
  LaserEventInfo *m_laserEventInfo{nullptr};

  TF1 *m_f0{nullptr};
  TF1 *m_f1{nullptr};
  TFile *m_debugFile{nullptr};
  TH1 *m_itHist_0{nullptr};
  TH1 *m_itHist_1{nullptr};
  TTree *m_hitTree{nullptr};

  uint64_t prev_BCO{0};

  double peakWidth0{-999};
  double peakWidth1{-999};

  int m_runnumber{0};
  int m_time_samples_max{425};
  int peakSample0{-999};
  int peakSample1{-999};

  bool isLaserEvent{false};
  bool isGl1LaserEvent{false};
  bool isGl1LaserPileupEvent{false};
  bool m_debug{false};

  std::string m_debugFileName = {"LaserEventIdentifier_debug.root"};
};

#endif
