// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHTPCCENTRALMEMBRANEMATCHER_H
#define PHTPCCENTRALMEMBRANEMATCHER_H

#include <string>

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainer.h>

class PHCompositeNode;
class CMFlashClusterContainer;
class CMFlashDifferenceContainer;

class TF1;
class TNtuple;
class TFile;
class TH1F;
class TH2F;
class TVector3;

class PHTpcCentralMembraneMatcher : public SubsysReco
{
 public:

 PHTpcCentralMembraneMatcher(const std::string &name = "PHTpcCentralMembraneMatcher");

  virtual ~PHTpcCentralMembraneMatcher();

  void set_process(const int proc)  { _process = proc;  }
  void set_histos_on(const bool val) {_histos = val;}

 //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode * topNode);

 protected:
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  CMFlashClusterContainer *_corrected_CMcluster_map{nullptr};
  CMFlashDifferenceContainer *_cm_flash_diffs{nullptr};

  TH2F *hxy_reco;
  TH2F *hxy_truth;
  TH2F *hdrdphi;
  TH2F *hrdr;
  TH2F *hrdphi;

  std::vector<TVector3> reco_pos;
  std::vector<TVector3> truth_pos;
  
  int _process = 0;
  bool _histos = true;

  TFile *fout;

/// check if coords are in a stripe
  //int getSearchResult(double xcheck, double ycheck) const;

  //  int getStripeID(double xcheck, double ycheck) const;

  /// adjust central membrane hits delay with respect to trigger time
  //void setCentralMembraneDelay(int ns) { m_centralMembraneDelay = ns; };

 private:
  /// detector name
  std::string detector = "TPC";


  static constexpr double mm = 1.0;
  static constexpr double cm = 10.0;

  /// inner radius of CM
  static constexpr double begin_CM = 221.4019814 * mm;

  /// outer radius of CM
  static constexpr double end_CM = 759.2138 * mm;

  static constexpr int nRadii = 8;
  static constexpr int nStripes_R1 = 6;
  static constexpr int nStripes_R2 = 8;
  static constexpr int nStripes_R3 = 12;

  static constexpr int nPads_R1 = 6 * 16;
  static constexpr int nPads_R2 = 8 * 16;
  static constexpr int nPads_R3 = 12 * 16;

  /// radius of arc on end of a stripe
  static constexpr double arc_r = 0.5 * mm;

  /// stripe radii
  static constexpr std::array<double, nRadii> R1_e = {{227.0902789 * mm, 238.4100043 * mm, 249.7297296 * mm, 261.049455 * mm, 272.3691804 * mm, 283.6889058 * mm, 295.0086312 * mm, 306.3283566 * mm}};
  static constexpr std::array<double, nRadii> R1 = {{317.648082 * mm, 328.9678074 * mm, 340.2875328 * mm, 351.6072582 * mm, 362.9269836 * mm, 374.246709 * mm, 385.5664344 * mm, 396.8861597 * mm}};
  static constexpr std::array<double, nRadii> R2 = {{421.705532 * mm, 442.119258 * mm, 462.532984 * mm, 482.9467608 * mm, 503.36069 * mm, 523.774416 * mm, 544.188015 * mm, 564.601868 * mm}};
  static constexpr std::array<double, nRadii> R3 = {{594.6048725 * mm, 616.545823 * mm, 638.4867738 * mm, 660.4277246 * mm, 682.3686754 * mm, 704.3096262 * mm, 726.250577 * mm, 748.1915277 * mm}};

  double cx1_e[nStripes_R1][nRadii];
  double cx1[nStripes_R1][nRadii];
  double cx2[nStripes_R2][nRadii];
  double cx3[nStripes_R3][nRadii];

  double cy1_e[nStripes_R1][nRadii];
  double cy1[nStripes_R1][nRadii];
  double cy2[nStripes_R2][nRadii];
  double cy3[nStripes_R3][nRadii];

  //Check which stripes get removed
  std::array<int, nRadii> nGoodStripes_R1_e = {};
  std::array<int, nRadii> nGoodStripes_R1 = {};
  std::array<int, nRadii> nGoodStripes_R2 = {};
  std::array<int, nRadii> nGoodStripes_R3 = {};

  ///min stripe index
  static constexpr std::array<int, nRadii> keepThisAndAfter = {{1, 0, 1, 0, 1, 0, 1, 0}};

  ///max stripe index
  static constexpr std::array<int, nRadii> keepUntil_R1_e = {{4, 4, 5, 4, 5, 5, 5, 5}};
  static constexpr std::array<int, nRadii> keepUntil_R1 = {{5, 5, 6, 5, 6, 5, 6, 5}};
  static constexpr std::array<int, nRadii> keepUntil_R2 = {{7, 7, 8, 7, 8, 8, 8, 8}};
  static constexpr std::array<int, nRadii> keepUntil_R3 = {{11, 10, 11, 11, 11, 11, 12, 11}};

  std::array<int, nRadii> nStripesIn_R1_e = {};
  std::array<int, nRadii> nStripesIn_R1 = {};
  std::array<int, nRadii> nStripesIn_R2 = {};
  std::array<int, nRadii> nStripesIn_R3 = {};
  std::array<int, nRadii> nStripesBefore_R1_e = {};
  std::array<int, nRadii> nStripesBefore_R1 = {};
  std::array<int, nRadii> nStripesBefore_R2 = {};
  std::array<int, nRadii> nStripesBefore_R3 = {};

  static constexpr int nStripesPerPetal = 213;
  static constexpr int nPetals = 18;
  static constexpr int nTotStripes = nStripesPerPetal * nPetals;

  double _rad_cut= 0.2;
  double _phi_cut= 0.01;

void CalculateCenters(
    int nPads,
    const std::array<double, nRadii>& R,
    std::array<int, nRadii>& nGoodStripes,
    const std::array<int, nRadii>& keepUntil,
    std::array<int, nRadii>& nStripesIn,
    std::array<int, nRadii>& nStripesBefore,
    double cx[][nRadii], double cy[][nRadii] );


};

#endif // PHTPCCENTRALMEMBRANEMATCHER_H
