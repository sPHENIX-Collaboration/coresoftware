// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef INTTCLUSTERQA_H
#define INTTCLUSTERQA_H

#include <fun4all/SubsysReco.h>

#include <cmath>
#include <map>
#include <set>
#include <string>

class PHCompositeNode;
class TH1;
class TH2;

class InttClusterQA : public SubsysReco
{
 public:
  InttClusterQA(const std::string &name = "InttClusterQA");

  ~InttClusterQA() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;

  void writeSensorInfo(bool value)
  {
    m_sensorInfo = value;
  }

 private:
  void createHistos();

  std::string getHistoPrefix() const;
  std::map<int, int> m_layerLadderMap;
  int m_event = 0;
  int m_totalClusters = 0;
  int m_nclustersPerSensor[4][16][4] = {{{0}}};
  bool m_sensorInfo = false;
  double cluszbin[56] = {-25., -22.57245 - 1.1,             //
                         -22.57245 - 0.9, -22.57245 + 0.9,  //
                         -20.57245 - 0.9, -20.57245 + 0.9,  //
                         -18.57245 - 0.9, -18.57245 + 0.9,  //
                         -16.57245 - 0.9, -16.57245 + 0.9,  //
                         -14.57245 - 0.9, -14.57245 + 0.9,  //
                         -12.57245 - 0.7, -12.57245 + 0.7,  //
                         -10.97245 - 0.7, -10.97245 + 0.7,  //
                         -9.372450 - 0.7, -9.372450 + 0.7,  //
                         -7.772450 - 0.7, -7.772450 + 0.7,  //
                         -6.172450 - 0.7, -6.172450 + 0.7,  //
                         -4.572450 - 0.7, -4.572450 + 0.7,  //
                         -2.972450 - 0.7, -2.972450 + 0.7,  //
                         -1.372450 - 0.7, -1.372450 + 0.7,  //
                         0.4275496 - 0.7, 0.4275496 + 0.7,  //
                         2.0275495 - 0.7, 2.0275495 + 0.7,  //
                         3.6275494 - 0.7, 3.6275494 + 0.7,  //
                         5.2275495 - 0.7, 5.2275495 + 0.7,  //
                         6.8275494 - 0.7, 6.8275494 + 0.7,  //
                         8.4275493 - 0.7, 8.4275493 + 0.7,  //
                         10.027549 - 0.7, 10.027549 + 0.7,  //
                         11.627549 - 0.7, 11.627549 + 0.7,  //
                         13.627549 - 0.9, 13.627549 + 0.9,  //
                         15.627549 - 0.9, 15.627549 + 0.9,  //
                         17.627550 - 0.9, 17.627550 + 0.9,  //
                         19.627550 - 0.9, 19.627550 + 0.9,  //
                         21.627550 - 0.9, 21.627550 + 0.9,  //
                         21.627550 + 1.1, 25.};

  TH1 *h_occupancy{nullptr};
  TH1 *h_clusSize{nullptr};
  TH1 *h_clusPhi_incl{nullptr};
  TH1 *h_clusPhi_l34{nullptr};
  TH1 *h_clusPhi_l56{nullptr};
  TH2 *h_clusZ_clusPhi_l34{nullptr};
  TH2 *h_clusZ_clusPhi_l56{nullptr};
  TH2 *h_cluspersensor[4][16][4] = {{{nullptr}}};
};

#endif  // InttClusterQA_H
