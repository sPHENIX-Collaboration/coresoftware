#ifndef TPC_TPCCLUSTERIZER_H
#define TPC_TPCCLUSTERIZER_H

#include <fun4all/SubsysReco.h>

#include <vector>
#include <string>

class PHCompositeNode;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;

class TpcClusterizer : public SubsysReco
{
 public:
  TpcClusterizer(const std::string &name = "TpcClusterizer");
  virtual ~TpcClusterizer(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

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
