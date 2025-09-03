// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCCALIB_TPCCENTRALMEMBRANEMATCHING_H
#define TPCCALIB_TPCCENTRALMEMBRANEMATCHING_H

/**
 * \file TpcCentralMembraneMatching.h
 * \brief match reconstructed CM clusters to CM pads, calculate differences, store on the node tree and compute distortion reconstruction maps
 * \author Tony Frawley <frawley@fsunuc.physics.fsu.edu>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

//#include <ffaobjects/SyncObject.h>

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>

class PHCompositeNode;
// class CMFlashClusterContainer;
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

class TpcCentralMembraneMatching : public SubsysReco
{
 public:
  TpcCentralMembraneMatching(const std::string &name = "TpcCentralMembraneMatching");

  ~TpcCentralMembraneMatching() override = default;

  /// set to true to store evaluation histograms and ttrees
  void setSavehistograms(bool value)
  {
    m_savehistograms = value;
  }

  /// output file name for evaluation histograms
  void setHistogramOutputfile(const std::string &outputfile)
  {
    m_histogramfilename = outputfile;
  }

  /// output file name for storing the space charge reconstruction matrices
  void setOutputfile(const std::string &outputfile)
  {
    m_outputfile = outputfile;
  }

  void setDebugOutputFile(const std::string &debugfile)
  {
    m_debugfilename = debugfile;
  }

  void set_nHitsInCluster_minimum(const unsigned int minHits)
  {
    m_nHitsInCuster_minimum = minHits;
  }

  void set_fixShifts(bool fixShifts)
  {
    m_fixShifts = fixShifts;
  }

  void set_fieldOn(bool fieldOn)
  {
    m_fieldOn = fieldOn;
  }

  void set_doFancy(bool fancy)
  {
    m_doFancy = fancy;
  }

  void set_doHadd(bool hadd)
  {
    m_doHadd = hadd;
  }

  void set_averageMode(bool averageMode)
  {
    m_averageMode = averageMode;
  }

  void set_event_sequence(int seq)
  {
    m_event_sequence = seq;
    m_event_index = 100*seq;
  }

  void set_grid_dimensions(int phibins, int rbins);

  //! run initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! end of process
  int End(PHCompositeNode *topNode) override;

 private:
  EventHeader *eventHeader{nullptr};

  int GetNodes(PHCompositeNode *topNode);

  double getPhiRotation_smoothed(TH1 *hitHist, TH1 *clustHist, bool side);

  std::vector<int> doGlobalRMatching(TH2 *r_phi, bool pos);

  void getRegionPhiRotation(bool side);

  int getClusterRMatch(double clusterR, int side);

  //! tpc distortion correction utility class
  TpcDistortionCorrection m_distortionCorrection;

  //! CMFlashClusterContainer *m_corrected_CMcluster_map{nullptr};
  LaserClusterContainer *m_corrected_CMcluster_map{nullptr};
  CMFlashDifferenceContainer *m_cm_flash_diffs{nullptr};

  //!@name distortion correction containers
  //@{
  /** used in input to correct CM clusters before calculating residuals */
  TpcDistortionCorrectionContainer *m_dcc_in_module_edge{nullptr};
  TpcDistortionCorrectionContainer *m_dcc_in_static{nullptr};
  TpcDistortionCorrectionContainer *m_dcc_in_average{nullptr};
  //@}

  //! fluctuation distortion container
  /** used in output to write fluctuation distortions */
  TpcDistortionCorrectionContainer *m_dcc_out{nullptr};

  /// local fluctuation distortion container
  /*
   * this one is used to aggregate multple CM events in a single map
   * it is not stored on the node tree but saved to a dedicated file at the end of the run
   */
  std::unique_ptr<TpcDistortionCorrectionContainer> m_dcc_out_aggregated;

  /// output file, to which aggregated central membrane distortion corrections are stored
  std::string m_outputfile{"CMDistortionCorrections.root"};

  ///@name evaluation histograms
  //@{
  bool m_savehistograms{false};
  std::string m_histogramfilename = "TpcCentralMembraneMatching.root";

  TH2 *hxy_reco{nullptr};
  TH2 *hxy_truth{nullptr};
  TH2 *hdrdphi{nullptr};
  TH2 *hrdr{nullptr};
  TH2 *hrdphi{nullptr};
  TH1 *hdrphi{nullptr};
  TH1 *hdphi{nullptr};
  TH1 *hdr1_single{nullptr};
  TH1 *hdr2_single{nullptr};
  TH1 *hdr3_single{nullptr};
  TH1 *hdr1_double{nullptr};
  TH1 *hdr2_double{nullptr};
  TH1 *hdr3_double{nullptr};
  TH1 *hnclus{nullptr};

  TFile *fout;

  TFile *m_debugfile;
  std::string m_debugfilename{"CMMatcher.root"};

  TH2 *truth_r_phi[2]{nullptr};

  TH2 *reco_r_phi[2]{nullptr};

  TH2 *m_matchResiduals[2]{nullptr};

  //  TNtuple *match_ntup {nullptr};
  TTree *match_tree{nullptr};

  bool m_useHeader{true};
  bool m_averageMode{false};

  int m_event_index{0};
  int m_event_sequence{0};
  bool m_matched{false};
  int m_truthIndex{0};
  float m_truthR{0.0};
  float m_truthPhi{0.0};
  float m_recoR{0.0};
  float m_recoPhi{0.0};
  float m_recoZ{0.0};
  float m_rawR{0.0};
  float m_rawPhi{0.0};
  float m_staticR{0.0};
  float m_staticPhi{0.0};
  bool m_side{false};
  bool m_fitMode{false};
  unsigned int m_adc{0};
  unsigned int m_nhits{0};
  unsigned int m_nLayers{0};
  unsigned int m_nIPhi{0};
  unsigned int m_nIT{0};
  float m_layersSD{0.0};
  float m_IPhiSD{0.0};
  float m_ITSD{0.0};
  float m_layersWeightedSD{0.0};
  float m_IPhiWeightedSD{0.0};
  float m_ITWeightedSD{0.0};
  int m_lowShift{0};
  int m_highShift{0};
  float m_phiRotation{0.0};
  float m_distanceToTruth{0.0};
  float m_NNDistance{0.0};
  float m_NNR{0.0};
  float m_NNPhi{0.0};
  int m_NNIndex{0};

  //@}

  unsigned int m_nHitsInCuster_minimum = 5;

  /// radius cut for matching clusters to pad, for size 2 clusters
  //  double m_rad_cut{0}.5;

  /// phi cut for matching clusters to pad
  /** TODO: this will need to be adjusted to match beam-induced time averaged distortions */
  double m_phi_cut{0.025};

  ///@name distortion correction histograms
  //@{

  /// distortion correction grid size along phi
  int m_phibins{24};

  static constexpr float m_phiMin{0};
  static constexpr float m_phiMax{2. * M_PI};

  /// distortion correction grid size along r
  int m_rbins{12};

  static constexpr float m_rMin{20};  // cm
  static constexpr float m_rMax{80};  // cm

  //@}

  ///@name central membrane pads definitions
  //@{
  static constexpr double mm{1.0};
  static constexpr double cm{10.0};

  static constexpr int nRadii{8};
  static constexpr int nStripes_R1{6};
  static constexpr int nStripes_R2{8};
  static constexpr int nStripes_R3{12};

  static constexpr int nPads_R1{6 * 16};
  static constexpr int nPads_R2{8 * 16};
  static constexpr int nPads_R3{12 * 16};

  /// stripe radii
  static constexpr std::array<double, nRadii> R1_e{{227.0902789 * mm, 238.4100043 * mm, 249.7297296 * mm, 261.049455 * mm, 272.3691804 * mm, 283.6889058 * mm, 295.0086312 * mm, 306.3283566 * mm}};
  static constexpr std::array<double, nRadii> R1{{317.648082 * mm, 328.9678074 * mm, 340.2875328 * mm, 351.6072582 * mm, 362.9269836 * mm, 374.246709 * mm, 385.5664344 * mm, 396.8861597 * mm}};
  static constexpr std::array<double, nRadii> R2{{421.705532 * mm, 442.119258 * mm, 462.532984 * mm, 482.9467608 * mm, 503.36069 * mm, 523.774416 * mm, 544.188015 * mm, 564.601868 * mm}};
  static constexpr std::array<double, nRadii> R3{{594.6048725 * mm, 616.545823 * mm, 638.4867738 * mm, 660.4277246 * mm, 682.3686754 * mm, 704.3096262 * mm, 726.250577 * mm, 748.1915277 * mm}};

  double cx1_e[nStripes_R1][nRadii]{};
  double cx1[nStripes_R1][nRadii]{};
  double cx2[nStripes_R2][nRadii]{};
  double cx3[nStripes_R3][nRadii]{};

  double cy1_e[nStripes_R1][nRadii]{};
  double cy1[nStripes_R1][nRadii]{};
  double cy2[nStripes_R2][nRadii]{};
  double cy3[nStripes_R3][nRadii]{};

  // Check which stripes get removed
  std::array<int, nRadii> nGoodStripes_R1_e = {};
  std::array<int, nRadii> nGoodStripes_R1 = {};
  std::array<int, nRadii> nGoodStripes_R2 = {};
  std::array<int, nRadii> nGoodStripes_R3 = {};

  /// min stripe index
  static constexpr std::array<int, nRadii> keepThisAndAfter{{1, 0, 1, 0, 1, 0, 1, 0}};

  /// max stripe index
  static constexpr std::array<int, nRadii> keepUntil_R1_e{{4, 4, 5, 4, 5, 5, 5, 5}};
  static constexpr std::array<int, nRadii> keepUntil_R1{{5, 5, 6, 5, 6, 5, 6, 5}};
  static constexpr std::array<int, nRadii> keepUntil_R2{{7, 7, 8, 7, 8, 8, 8, 8}};
  static constexpr std::array<int, nRadii> keepUntil_R3{{11, 10, 11, 11, 11, 11, 12, 11}};

  std::array<int, nRadii> nStripesIn_R1_e{};
  std::array<int, nRadii> nStripesIn_R1{};
  std::array<int, nRadii> nStripesIn_R2{};
  std::array<int, nRadii> nStripesIn_R3{};
  std::array<int, nRadii> nStripesBefore_R1_e{};
  std::array<int, nRadii> nStripesBefore_R1{};
  std::array<int, nRadii> nStripesBefore_R2{};
  std::array<int, nRadii> nStripesBefore_R3{};

  static constexpr int nStripesPerPetal{213};
  static constexpr int nPetals{18};
  static constexpr int nTotStripes{nStripesPerPetal * nPetals};

  void CalculateCenters(
      int nPads,
      const std::array<double, nRadii> &R,
      std::array<int, nRadii> &nGoodStripes,
      const std::array<int, nRadii> &keepUntil,
      std::array<int, nRadii> &nStripesIn,
      std::array<int, nRadii> &nStripesBefore,
      double cx[][nRadii], double cy[][nRadii]);

  /// store centers of all central membrane pads
  std::vector<TVector3> m_truth_pos;
  std::vector<int> m_truth_index;

  std::vector<double> m_truth_RPeaks{22.709, 23.841, 24.973, 26.1049, 27.2369, 28.3689, 29.5009, 30.6328, 31.7648, 32.8968, 34.0288, 35.1607, 36.2927, 37.4247, 38.5566, 39.6886, 42.1706, 44.2119, 46.2533, 48.2947, 50.3361, 52.3774, 54.4188, 56.4602, 59.4605, 61.6546, 63.8487, 66.0428, 68.2369, 70.431, 72.6251, 74.8192};

  //@}

  std::vector<double> phiSpacing = {0.0749989, 0.0729917, 0.0711665, 0.0694996, 0.0679712, 0.0665648, 0.0652663, 0.0640638, 0.062947, 0.0619071, 0.0609364, 0.0600281, 0.0591765, 0.0583765, 0.0576234, 0.0569132, 0.0473084, 0.0462573, 0.045299, 0.0444217, 0.0436155, 0.0428722, 0.0421847, 0.0415468, 0.0325076, 0.0319331, 0.031398, 0.0308985, 0.0304311, 0.0299928, 0.029581, 0.0291934};

  bool m_fixShifts{false};
  bool m_fieldOn{true};
  bool m_doFancy{false};
  bool m_doHadd{false};

  std::vector<double> m_reco_RPeaks[2];
  double m_m[2]{};
  double m_b[2]{};
  int m_matchLow[2]{};
  int m_matchHigh[2]{};
  std::vector<int> m_reco_RMatches[2];

  double m_recoRotation[2][3]{{-999, -999, -999}, {-999, -999, -999}};
};

#endif  // PHTPCCENTRALMEMBRANEMATCHER_H
