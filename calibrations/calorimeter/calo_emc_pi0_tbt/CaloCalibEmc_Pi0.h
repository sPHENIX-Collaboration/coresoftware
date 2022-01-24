// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOCALIBEMC_PI0_H
#define CALOCALIBEMC_PI0_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class TH1F;
class TH2F;
class TH3F;
class TH1;
class TTree;
class TString;

class CaloCalibEmc_Pi0 : public SubsysReco
{
 public:

  CaloCalibEmc_Pi0(const std::string &name = "CaloCalibEmc_Pi0",const std::string &fnm = "outJF");

  virtual ~CaloCalibEmc_Pi0() {}

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

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void Loop(TString _filename, int nevts);
  
  void FittingHistos();


 private:
  int _ievent;
  TFile *cal_output;
  std::string _caloname;
  std::string _filename;


 // histos lists
  TH1 *cemc_hist_eta_phi[96][258];
  TH1 *eta_hist[96]= {0};
  TH2F *mass_eta;// = Null;
  TH3F *mass_eta_phi;// = Null;

  TH1F *pairInvMassTotal;

  TTree * _eventTree;
  // TTree variables
  int _eventNumber;
  int _nClusters;
  float _clusterIDs[10000]={0};
  float _clusterEnergies[10000]={0};
  float _clusterPts[10000]={0};
  int _clusterEtas[10000]={0};
  int _clusterPhis[10000]={0};
  

  int maxTowerEta;
  int maxTowerPhi; 
  
  int _maxTowerEtas[10000]={0};
  int _maxTowerPhis[10000]={0};

  float alphaCut;


};

#endif  //   CALOCALIBEMC_PI0_H
