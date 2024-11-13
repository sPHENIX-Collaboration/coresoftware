// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTREEGEN_H
#define CALOTREEGEN_H

#include <fun4all/SubsysReco.h>

#include <limits>
#include <string>
#include <vector>

class PHCompositeNode;
class RawCluster;
class TTree;
class TFile;
class TH1;

class caloTreeGen : public SubsysReco
{
 public:
  explicit caloTreeGen(const std::string &name = "caloTreeGen");

  ~caloTreeGen() override = default;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
  */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
  */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void doClusters(int clusters, const std::string & clusterNode)
  {
    storeClusters = clusters;
    m_clusterNode = clusterNode;
  }

  void setClusterThresh(float thresh)
  {
    clusterThresh = thresh;
  }

  void doClusterDetails(int fineCluster) { storeClusterDetails = fineCluster; }

  void doEMCal(int emcalOn, const std::string & emcNode)
  {
    storeEMCal = emcalOn;
    m_emcTowerNode = emcNode;
  }

  void setEMCalThresh(float thresh)
  {
    emcalThresh = thresh;
  }

  void doHCals(int hcalsOn, const std::string & ohcNode, const std::string & ihcNode)
  {
    storeHCals = hcalsOn;
    m_ohcTowerNode = ohcNode;
    m_ihcTowerNode = ihcNode;
  }

  void setOHCalThresh(float thresh)
  {
    ohcalThresh = thresh;
  }

  void setIHCalThresh(float thresh)
  {
    ihcalThresh = thresh;
  }

  void doZDC(int zdcOn, const std::string & zdcNode)
  {
    storeZDC = zdcOn;
    m_zdcTowerNode = zdcNode;
  }

  void doTrig(int trigOn, const std::string & trigNode)
  {
    storeTrig = trigOn;
    m_trigNode = trigNode;
  }

 private:
  // private methods
  float getMaxTowerE(RawCluster *cluster);
  std::vector<float> returnClusterTowE(RawCluster *cluster);
  std::vector<int> returnClusterTowPhi(RawCluster *cluster);
  std::vector<int> returnClusterTowEta(RawCluster *cluster);

  // pointers
  TTree *T{nullptr};

  TFile *out{nullptr};

  TH1 *zVertex{nullptr};

  // simple variables
  float m_vertex{std::numeric_limits<float>::quiet_NaN()};
  float totalCaloEEMCal{std::numeric_limits<float>::quiet_NaN()};
  float totalCaloEOHCal{std::numeric_limits<float>::quiet_NaN()};
  float totalCaloEIHCal{std::numeric_limits<float>::quiet_NaN()};
  float totalCaloEZDC{std::numeric_limits<float>::quiet_NaN()};
  float clusterThresh{std::numeric_limits<float>::quiet_NaN()};
  float ohcalThresh{std::numeric_limits<float>::quiet_NaN()};
  float ihcalThresh{std::numeric_limits<float>::quiet_NaN()};
  float emcalThresh{std::numeric_limits<float>::quiet_NaN()};

  int storeClusters{1};
  int storeClusterDetails{1};
  int storeEMCal{1};
  int storeHCals{1};
  int storeZDC{1};
  int storeTrig{1};
  int eventNum{0};
  
  std::string m_clusterNode;
  std::string m_emcTowerNode;
  std::string m_ihcTowerNode;
  std::string m_ohcTowerNode;
  std::string m_trigNode;
  std::string m_zdcTowerNode;

  std::string Outfile{"commissioning.root"};

  // EMCal
  std::vector<float> m_emcTowE;
  std::vector<float> m_emciEta;
  std::vector<float> m_emciPhi;
  std::vector<int> m_emcTime;
  std::vector<float> m_emcChi2;
  std::vector<float> m_emcPed;

  // OHCal
  std::vector<float> m_ohciTowPhi;
  std::vector<float> m_ohciTowEta;
  std::vector<float> m_ohcTowE;
  std::vector<int> m_ohcTime;
  std::vector<float> m_ohcChi2;
  std::vector<float> m_ohcPed;

  // IHCal
  std::vector<float> m_ihciTowPhi;
  std::vector<float> m_ihciTowEta;
  std::vector<float> m_ihcTowE;
  std::vector<int> m_ihcTime;
  std::vector<float> m_ihcChi2;
  std::vector<float> m_ihcPed;

  // ZDC
  std::vector<float> m_zdcTowE;
  std::vector<int> m_zdcSide;

  // SMD
  std::vector<float> m_smdE;
  std::vector<int> m_smdSide;

  // Clusters
  std::vector<float> m_clusterE;
  std::vector<float> m_clusterPhi;
  std::vector<float> m_clusterEta;
  std::vector<float> m_clusterPt;
  std::vector<float> m_clusterChi;
  std::vector<float> m_clusterNtow;
  std::vector<float> m_clusterTowMaxE;
  std::vector<float> m_clusterECore;

  std::vector<std::vector<int> > m_clusTowEta;
  std::vector<std::vector<int> > m_clusTowPhi;
  std::vector<std::vector<float> > m_clusTowE;

  // GL1 information
  std::vector<bool> m_triggerVector;
};

#endif
