#ifndef __TpcCLUSTERIZER_H__
#define __TpcCLUSTERIZER_H__

#include <RVersion.h>
#include <fun4all/SubsysReco.h>
#include <limits.h>
#include <vector>

class PHG4CylinderCellGeom;
class TH1F;
class TProfile2D;
class TStopwatch;
class TrkrHitSetContainer;
class TrkrHit;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;

class TpcClusterizer : public SubsysReco
{
 public:
  TpcClusterizer(const char *name = "TpcClusterizer");
  ~TpcClusterizer();

  int Init(PHCompositeNode *topNode) { return 0; }
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode) { return 0; }

 private:
  bool is_local_maximum(int phi, int z, std::vector<std::vector<double>> &adcval);
  void get_cluster(int phibin, int zbin, int &phiup, int &phidown, int &zup, int &zdown, std::vector<std::vector<double>> &adcval);

  TrkrHitSetContainer *m_hits;
  TrkrClusterContainer *m_clusterlist;
  TrkrClusterHitAssoc *m_clusterhitassoc;

  double zz_shaping_correction;
  double pedestal;

  int NPhiBinsMax;
  int NPhiBinsMin;
  int NZBinsMax;
  int NZBinsMin;

};

#endif
