#ifndef TPC_LASEREVENTIDENTIFIER_H
#define TPC_LASEREVENTIDENTIFIER_H

#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>


#include <TFile.h>
#include <TH1I.h>
#include <TTree.h>

class PHCompositeNode;
class TrkrHitSet;
class TrkrHitSetContainer;
class PHG4TpcCylinderGeom;
class PHG4TpcCylinderGeomContainer;
class LaserEventInfo;
class LaserEventIdentifier : public SubsysReco
{
 public:
  LaserEventIdentifier(const std::string &name = "LaserEventIdentifier");
  ~LaserEventIdentifier() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_max_time_samples(int val) { m_time_samples_max = val; }

  void set_runnumber(int runnum) { m_runnumber = runnum; }

  void set_debug(bool debug) { m_debug = debug; }
  void set_debug_name(const std::string &name) { m_debugFileName = name; }
 private:
  int m_time_samples_max=425;

  TrkrHitSetContainer *m_hits = nullptr;
  ActsGeometry *m_tGeometry = nullptr;
  PHG4TpcCylinderGeomContainer *m_geom_container = nullptr;

  LaserEventInfo *m_laserEventInfo = nullptr ;
  bool m_debug = false;
  std::string m_debugFileName = "LaserEventIdentifier_debug.root";
  TFile *m_debugFile = nullptr;
  TTree *m_hitTree = nullptr;
  TH1I *m_itHist_0 = nullptr;
  TH1I *m_itHist_1 = nullptr;
  bool isLaserEvent = false;
  bool isGl1LaserEvent = false;
  bool isGl1LaserPileupEvent = false;
  int peakSample0 = -999;
  int peakSample1 = -999;
  float peakWidth0 = -999;
  float peakWidth1 = -999;
  int m_runnumber = 0;

  uint64_t prev_BCO = 0;
};

#endif
