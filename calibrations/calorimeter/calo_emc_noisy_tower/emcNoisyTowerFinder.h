// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TOWERID_H
#define TOWERID_H

#include <fun4all/SubsysReco.h>
//#include <cdbobjects/CDBTTree.h>

#include <string>
#include <vector>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

class TTree;
class PHCompositeNode;
class Fun4AllHistoManager;
class TFile;
class RawCluster;
class TowerInfoContainer;

const int nTowers = 24576;
const int nIB = 384;
const int nSectors = 64;
const int nTowersIB = 64;
const int nTowersSec = 384;

class emcNoisyTowerFinder : public SubsysReco
{
 public:

  explicit emcNoisyTowerFinder(const std::string &name = "emcNoisyTowerFinder", const std::string &outputName = "emcNoisyTowerFinder.root", const std::string &cdbtreename = "test.root", float adccut_sg = 250,float adccut_k = 500, float sigmas_lo = 1, float sigmas_hi = 4.5, float SG_f = 0.55, float Kur_f = 0.55, float region_f = 0.55);

  ~emcNoisyTowerFinder() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
  */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
  */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
  */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;	

  void Print(const std::string &what = "ALL") const override;
  
  //histogram filler so we don't have to call the End function
  void FillHistograms(const int runnumber, const int segment);

  void CalculateCutOffs(const int runnumber);

  void WriteCDBTree(const int runnumber);
  
 private:

  TTree *T = NULL;
  TFile *out = NULL;

//  CDBTTree *cdbttree;
 
  TFile *fchannels;
  TTree *channels;

  int fiber_type = 0;
 
  //Fun4AllHistoManager *hm = nullptr;
  std::string Outfile = "commissioning.root";
  
  TH1F* hEventCounter = NULL;
  
  TH2F* Fspec = NULL;
  TH2F* Fspec_SG = NULL;
  TH2F* Fspec_K = NULL;
  TH2F* Fspec_sector = NULL;
  TH2F* Fspec_IB = NULL;
 
  TH2F* Fspeci = NULL;
  TH2F* Fspeci_SG = NULL;
  TH2F* Fspeci_K = NULL;
  TH2F* Fspeci_sector = NULL;
  TH2F* Fspeci_IB = NULL;

  TH2F* Espec = NULL;
  TH2F* Espec_SG = NULL; 
  TH2F* Espec_K = NULL; 
  TH2F* Espec_sector = NULL;
  TH2F* Espec_IB = NULL;

  const std::string cdbtreename; 
  
  float adccut_sg;
  float adccut_k;
  float sigmas_lo;
  float sigmas_hi;
  float SG_f;
  float Kur_f;
  float region_f;
	
  int m_hot_channels = 0;
  float towerF[nTowers] = {0};
  float sectorF[nSectors] = {0};
  float ibF[nIB] = {0};

  float towerE[nTowers] = {0};
  float sectorE[nSectors] = {0};
  float ibE[nIB] = {0};

  int hottowers[nTowers] = {0};
  int hotIB[nIB] = {0};
  int hotsectors[nSectors] = {0};
  int deadtowers[nTowers] = {0};

	
  int coldtowers[nTowers] = {0};
  int coldIB[nIB] = {0};
  int coldsectors[nSectors] = {0};
  int hot_regions = 0;
  int cold_regions = 0;  


  int goodevents[nTowers] = {0};
  int goodeventsIB[nIB] = {0};
  int goodeventsSec[nSectors] = {0};
};

#endif 
