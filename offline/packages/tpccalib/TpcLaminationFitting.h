#ifndef TPCCALIB_TPCLAMINATIONFITTING_H
#define TPCCALIB_TPCLAMINATIONFITTING_H

#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

class PHCompositeNode;

class LaserClusterContainer;
class EventHeader;

class TF1;
class TFile;
class TH1;
class TH2;
class TGraph;
class TNtuple;
class TTree;
class TVector3;

class TpcLaminationFitting : public SubsysReco
{
 public:
  TpcLaminationFitting(const std::string &name = "TpcLaminationFitting");
  ~TpcLaminationFitting() override = default;

  /// output file name for storing the space charge reconstruction matrices
  void setOutputfile(const std::string &outputfile);

  void set_event_sequence(int seq)
  {
    m_event_sequence = seq;
    m_event_index = 100 * seq;
  }

  void set_QAFileName(const std::string &QAFileName)
  {
    m_QAFileName = QAFileName;
  }

  void set_stripePatternFile(const std::string &stripePatternFile)
  {
    m_stripePatternFile = stripePatternFile;
  }

  void set_ppMode(bool mode){ ppMode = mode; }

  void set_fieldOff(bool fieldOff){ m_fieldOff = fieldOff; }

  void set_grid_dimensions(int phibins, int rbins);

  void set_nLayerCut(unsigned int cut) { m_nLayerCut = cut; }

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

 private:
  EventHeader *eventHeader{nullptr};

  int GetNodes(PHCompositeNode *topNode);

  int fitLaminations();
  int InterpolatePhiDistortions();
  int doGlobalRMatching(int side);
  void fill_guarding_bins(TpcDistortionCorrectionContainer *dcc);

  TpcDistortionCorrection m_distortionCorrection;

  LaserClusterContainer *m_correctedCMcluster_map{nullptr};

  TpcDistortionCorrectionContainer *m_dcc_in_module_edge{nullptr};
  TpcDistortionCorrectionContainer *m_dcc_in_static{nullptr};

  TpcDistortionCorrectionContainer *m_dcc_out{nullptr};

  std::string m_outputfile{"CMDistortionCorrections.root"};
  std::string m_QAFileName{""};

  TH2 *m_hLamination[18][2]{{nullptr}};
  TF1 *m_fLamination[18][2]{{nullptr}};
  double m_laminationIdeal[18][2]{{0.0}};
  //double m_laminationCenter[18][2]{{0.0}};
  double m_laminationOffset[18][2]{{0.0}};
  //double m_laminationOffset{0.00337078};
  //double m_laminationOffset{0.002775};
  bool m_laminationGoodFit[18][2]{{false}};
  double m_distanceToFit[18][2]{{0.0}};
  int m_nBinsFit[18][2]{{0}};
  double m_fitRMSE[18][2]{{0.0}};


  TH2 *m_hPetal[2]{nullptr};
  TGraph *m_bestRMatch[2]{nullptr};
  TH2 *m_parameterScan[2]{nullptr};


  TH2 *phiDistortionLamination[2]{nullptr};
  //TH2 *scaleFactorMap[2]{nullptr};

  unsigned int m_nLayerCut{1};
  
  bool m_useHeader{true};

  int m_event_index{0};
  int m_event_sequence{0};

  double m_nClusters{0};
  int m_nEvents{0};
  int m_runnumber{};

  bool ppMode{false};
  double m_ZDC_coincidence{0};
  //std::map<int, float>  m_run_ZDC_map_pp;
  //std::map<int, float>  m_run_ZDC_map_auau;

  std::string m_stripePatternFile = "/sphenix/u/bkimelman/CMStripePattern.root";

  bool m_fieldOff{false};

  TTree *m_laminationTree{nullptr};
  bool m_side{false};
  int m_lamIndex{0};
  double m_lamPhi{0};
  double m_lamShift{0};
  bool m_goodFit{false};
  double m_A{0};
  double m_B{0};
  double m_C{0};
  double m_A_err{0};
  double m_B_err{0};
  double m_C_err{0};
  double m_dist{0};
  double m_rmse{};
  int m_nBins{0};

  int m_phibins{80};
  static constexpr float m_phiMin{0};
  static constexpr float m_phiMax{2. * M_PI};

  int m_rbins{52};
  static constexpr float m_rMin{20};  // cm
  static constexpr float m_rMax{80};  // cm

  /*
  const int nRadii{8};
  const int nStripes[4]{6,6,8,12};
  const int nPads[4]{96,96,128,192};
  const double RValues[4][8] = {{22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.63283566},{31.7648082, 32.89678074, 34.02875328, 35.16072582, 36.29269836, 37.4246709, 38.55664344, 39.68861597},{42.1705532, 44.2119258, 46.2532984, 48.29467608, 50.336069, 52.3774416, 54.4188015, 56.4601868},{59.46048725, 61.6545823, 63.84867738, 66.04277246, 68.23686754, 70.43096262, 72.6250577, 74.81915277}};

  const int keepThisAndAfter[8]{1,0,1,0,1,0,1,0};
  const int keepUntil[4][8]{{4,4,5,4,5,5,5,5},{5,5,6,5,6,5,6,5},{7,7,8,7,8,8,8,8},{11,10,11,11,11,11,12,11}};

  const double phi_petal = M_PI/6.0;
  const int pr_mult = 3;
  const int dw_mult = 8;
  const double diffwidth = 0.06;
  const double adjust = 0.015;
  */

  std::vector<double> m_truthR[2];
  std::vector<double> m_truthPhi[2];

  double m_phiModMin[2]{-M_PI/18, 0.0};
  double m_phiModMax[2]{M_PI/18, M_PI/9};
};

#endif
