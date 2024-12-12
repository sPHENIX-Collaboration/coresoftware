#ifndef TPCCALIB_TPCLAMINATIONFITTING_H
#define TPCCALIB_TPCLAMINATIONFITTING_H

#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

class PHCompositeNode;

class CMFlashDifferenceContainer;
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

  void set_fitFileName(const std::string &fitFileName)
  {
    m_fitFileName = fitFileName;
  }

  void set_grid_dimensions(int phibins, int rbins);

  void set_nLayerCut(unsigned int cut) { m_nLayerCut = cut; }

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

 private:
  EventHeader *eventHeader{nullptr};

  int GetNodes(PHCompositeNode *topNode);

  int fitLaminations();
  int InterpolatePhiDistortions(TH2 *simPhiDistortion[2]);
  void fill_guarding_bins(TpcDistortionCorrectionContainer *dcc);

  TpcDistortionCorrection m_distortionCorrection;

  LaserClusterContainer *m_correctedCMcluster_map{nullptr};
  CMFlashDifferenceContainer *m_cm_flash_diffs{nullptr};

  TpcDistortionCorrectionContainer *m_dcc_in_module_edge{nullptr};
  TpcDistortionCorrectionContainer *m_dcc_in_static{nullptr};

  TpcDistortionCorrectionContainer *m_dcc_out{nullptr};

  std::string m_outputfile{"CMDistortionCorrections.root"};
  std::string m_fitFileName{""};

  TH2 *m_hLamination[18][2]{{nullptr}};
  TF1 *m_fLamination[18][2]{{nullptr}};
  double m_laminationCenter[18][2]{{0.0}};
  bool m_laminationGoodFit[18][2]{{false}};
  double m_distanceToFit[18][2]{{0.0}};
  int m_nBinsFit[18][2]{{0}};

  TH2 *phiDistortionLamination[2]{nullptr};
  TH2 *scaleFactorMap[2]{nullptr};

  unsigned int m_nLayerCut{1};
  
  bool m_useHeader{true};

  int m_event_index{0};
  int m_event_sequence{0};

  double m_nClusters{0};
  int m_nEvents{0};

  TTree *m_laminationTree{nullptr};
  int m_runnumber{0};
  int m_segment{0};
  bool m_side{false};
  int m_lamIndex{0};
  double m_lamPhi{0};
  double m_A{0};
  double m_B{0};
  double m_C{0};
  double m_dist{0};
  int m_nBins{0};

  int m_phibins{24};
  static constexpr float m_phiMin{0};
  static constexpr float m_phiMax{2. * M_PI};

  int m_rbins{12};
  static constexpr float m_rMin{20};  // cm
  static constexpr float m_rMax{80};  // cm
};

#endif
