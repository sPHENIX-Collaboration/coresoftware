// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOEMCNOISYTOWER_EMCNOISYTOWERFINDER_H
#define CALOEMCNOISYTOWER_EMCNOISYTOWERFINDER_H

#include <fun4all/SubsysReco.h>
//#include <cdbobjects/CDBTTree.h>

#include <array>
#include <string>
#include <vector>

class Fun4AllHistoManager;
class PHCompositeNode;
class RawCluster;
class TFile;
class TH2;
class TowerInfoContainer;
class TTree;

const int nTowers = 24576;
const int nIB = 384;
const int nSectors = 64;
const int nTowersIB = 64;
const int nTowersSec = 384;

class emcNoisyTowerFinder : public SubsysReco
{
 public:
  explicit emcNoisyTowerFinder(const std::string& name = "emcNoisyTowerFinder", const std::string& outputName = "emcNoisyTowerFinder.root", const std::string& cdbtreename_in = "test.root", float adccut_sg_in = 250, float adccut_k_in = 500, float sigmas_lo_in = 1, float sigmas_hi_in = 4.5, float SG_f_in = 0.55, float Kur_f_in = 0.55, float region_f_in = 0.55);

  ~emcNoisyTowerFinder() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
  */
  int Init(PHCompositeNode* topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
  */
  int InitRun(PHCompositeNode* topNode) override;

  /** Called for each event.
      This is where you do the real work.
  */
  int process_event(PHCompositeNode* topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode* topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode* topNode) override;

  /// Reset
  int Reset(PHCompositeNode* /*topNode*/) override;

  void Print(const std::string& what = "ALL") const override;

  // histogram filler so we don't have to call the End function
  void FillHistograms(const int runnumber, const int segment);

  void CalculateCutOffs(const int runnumber);

  void WriteCDBTree(const int runnumber);

  void set_Normalization(int norm) { do_VariableNormalization = norm; }

 private:
  TTree* T{nullptr};
  TFile* out{nullptr};

  //  CDBTTree *cdbttree{nullptr};

  TFile* fchannels{nullptr};
  TTree* channels{nullptr};

  int fiber_type{0};

  // Fun4AllHistoManager *hm {nullptr};
  std::string Outfile{"commissioning.root"};

  //  TH1F* hEventCounter {nullptr};

  TH2* Fspec{nullptr};
  TH2* Fspec_SG{nullptr};
  TH2* Fspec_K{nullptr};
  TH2* Fspec_sector{nullptr};
  TH2* Fspec_IB{nullptr};

  TH2* Fspeci{nullptr};
  TH2* Fspeci_SG{nullptr};
  TH2* Fspeci_K{nullptr};
  TH2* Fspeci_sector{nullptr};
  TH2* Fspeci_IB{nullptr};

  TH2* Espec{nullptr};
  TH2* Espec_SG{nullptr};
  TH2* Espec_K{nullptr};
  TH2* Espec_sector{nullptr};
  TH2* Espec_IB{nullptr};

  const std::string cdbtreename;

  float adccut_sg{0};
  float adccut_k{0};
  float sigmas_lo{0};
  float sigmas_hi{0};
  float SG_f{0};
  float Kur_f{0};
  float region_f{0};

  int m_hot_channels{0};
  int eventCounter{1};
  int do_VariableNormalization{0};
  std::array<float, nTowers> towerF{0};
  std::array<float, nSectors> sectorF{0};
  std::array<float, nIB> ibF{0};

  std::array<float, nTowers> towerE{0};
  std::array<float, nSectors> sectorE{0};
  std::array<float, nIB> ibE{0};

  std::array<int, nTowers> hottowers{0};
  std::array<int, nIB> hotIB{0};
  //  std::array<int, nSectors> hotsectors {0};
  std::array<int, nTowers> deadtowers{0};

  std::array<int, nTowers> coldtowers{0};
  //  std::array<int, nIB> coldIB{0};
  //  std::array<int, nSectors> coldsectors{0};
  int hot_regions{0};
  int cold_regions{0};

  std::array<int, nTowers> goodevents{0};
  std::array<int, nIB> goodeventsIB{0};
  std::array<int, nSectors> goodeventsSec{0};
};

#endif
