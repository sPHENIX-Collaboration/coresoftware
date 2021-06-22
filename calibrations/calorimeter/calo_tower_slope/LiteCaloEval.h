// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef LITECALOEVAL_H
#define LITECALOEVAL_H

#include <fun4all/SubsysReco.h>
#include <TFile.h>

#include <TNtuple.h>

#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

 
class PHCompositeNode;

class LiteCaloEval : public SubsysReco
{
 public:


  LiteCaloEval(const std::string &name = "LiteCaloEval", const std::string & caloNm = "CEMC", const std::string & fnm = "outJF");

  virtual ~LiteCaloEval();

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

  char user_choice;
  

 private:
  TFile * _tfile = nullptr;
  int _ievent = 0;
  TNtuple * _ntp_tower = nullptr;
  std::string _caloname;
  std::string _filename;
  
  TFile *cal_output = nullptr;

  TH1F *hcal_out_eta_phi[24][64] = {};
  TH1F *hcalout_eta[24] = {};
  TH2F *hcalout_energy_eta = nullptr; // = new TH2F("hcalout_energy_eta", "hcalout energy eta", 10,0,10,24000,-1.1,1.1);
  TH3F *hcalout_e_eta_phi = {};// = new TH3F("hcalout_e_eta_phi", "hcalout e eta phi",50,0,10,24,-1.1,1.1,64,-3.14159,3.14159);
  
  TH1F *hcal_in_eta_phi[24][64] = {};
  TH1F *hcalin_eta[24] = {};
  TH2F *hcalin_energy_eta = nullptr;// = new TH2F("hcalin_energy_eta", "hcalin energy eta", 1000,0,10,240,-1.1,1.1);
  TH3F *hcalin_e_eta_phi = nullptr;// = new TH3F("hcalin_e_eta_phi", "hcalin e eta phi",50,0,10,24,-1.1,1.1,64,-3.14159,3.14159);

  TH1F *cemc_hist_eta_phi[96][258] = {};// = {0};
  TH1F *eta_hist[96] = {};// = {0};
  TH2F *energy_eta_hist = nullptr;// = Null;
  TH3F *e_eta_phi = nullptr;// = Null;
};

#endif // LITECALOEVAL_H
