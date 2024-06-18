// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_TRACKING_TPCCLUSTERQA_H
#define QA_TRACKING_TPCCLUSTERQA_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <vector>

class PHCompositeNode;
class TpcDistortionCorrectionContainer;

class TpcClusterQA : public SubsysReco
{
 public:
  TpcClusterQA(const std::string &name = "TpcClusterQA");

  ~TpcClusterQA() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;

  void makeResidQAHistos (bool value)
  {
    m_residQA = value;
  }
 
 private:
  void createHistos();

  bool m_residQA = false;

  std::string m_trackMapName = "SvtxTrackMap";
  TpcDistortionCorrectionContainer *m_dccStatic{nullptr}, *m_dccAverage{nullptr}, *m_dccFluctuation{nullptr};
  float m_px = std::numeric_limits<float>::quiet_NaN();
  float m_py = std::numeric_limits<float>::quiet_NaN();
  float m_pt = std::numeric_limits<float>::quiet_NaN();    
  int m_ntpc = std::numeric_limits<int>::quiet_NaN();
  std::vector<float> m_clusgz;
  std::vector<int> m_cluslayer;
  std::vector<int> m_clusphisize;
  std::vector<int> m_cluszsize;
  std::vector<int> m_region;

  std::string getHistoPrefix() const;
  std::set<int> m_layers;
  std::multimap<int, int> m_layerRegionMap;
   
  int m_event = 0;
  int m_totalClusters = 0;
  int m_clustersPerSector[24] = {0};
};

#endif  // QA_TRACKING_TPCCLUSTERQA_H
