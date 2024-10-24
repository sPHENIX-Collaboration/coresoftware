#ifndef INTTMIXUPQA_H__
#define INTTMIXUPQA_H__

// std headers
#include <vector>
#include <filesystem>
#include <array>
#include <iostream>
#include <iomanip> // setw, setfill

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <set>

// ROOT headers
#include <TObject.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLine.h>
#include <TLegend.h>
#include <TH2.h>
#include <TROOT.h>
#include <TGraph.h>

// Fun4All headers
#include <fun4all/SubsysReco.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

#include <trackbase/InttEventInfo.h>
#include <trackbase/InttEventInfov1.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>


class PHCompositeNode;

class InttMixupQA : public SubsysReco {


 public:
  InttMixupQA(const std::string &name = "InttMixupQA", const int run_num=0,const int felix_num=0);

  virtual ~InttMixupQA();

  int Init(PHCompositeNode *);
  
  int InitRun(PHCompositeNode *);
  
  /// SubsysReco event processing method
  int process_event(PHCompositeNode *);

  /// SubsysReco end processing method
  int EndRun(PHCompositeNode *);

  int End(PHCompositeNode *);

  int SetOutputDir( std::string dir );

  void SetBcoPeakFileDir (std::string path){bcopeak_dir_=path;} ;

  void SetHotChanFileDir (std::string path){hotchan_dir_=path;};

  //int SetHistBin(std::string type);
 private:

  // general variables
  int run_num_ = 0;
  int felix_num_=0;
  static const int kFelix_num_ = 8; // the number of our FELIX server
  static const int kFee_num_ = 14;  // the number of half-ladders in a single FELIX server
  static const int kChip_num_ = 26; // the number of chip in a half-ladder
  static const int kChan_num_ = 128; // the number of channel in a single chip
  static const int kFirst_pid_ = 3001; // the first pid (packet ID), which means intt0
  static const int divimul=10;


  // variables for the output
  std::string output_dir_ = "./";
  std::string output_basename_ = "InttMixupEventQA_run";
  std::string output_root_= "mixup.root";
  std::string output_pdf_ = "mixup.pdf";
  std::string output_txt_ ="bcopeak.pdf";
  TFile* tf_output_;

  // variables for get bco peak
  std::string bcopeak_dir_ ="./";
  //std::string bcopeak_file;
  std::string bcopeak_file[kFelix_num_];
  bool force_suffix_=false;
  bool is_official_=false;
  //bool is_official_=true;
  std::string suffix_;

  std::string GetFileSuffix(){
    if (force_suffix_==false)
    {
      if(is_official_==true)
       return "_official";
      else
       return "_special";
    }
    return suffix_;
  };

  //variables for hot channel cut
  std::string hotchan_dir_="./";
  std::string hotchan_file;
  
  
  std::set<int> bcopar_[kFelix_num_];
  std::set<int> otbcopar_[kFelix_num_];
  std::map<int,int> hotmap;

  int ievent_ = 0;
  int n=kFelix_num_;

  int prev_bcofull= 0;
  uint64_t long_prev_bcofull=0;
  int pre_allhit[kFelix_num_];
  
  int NmixupEv[kFelix_num_] = {0,0,0,0,0,0,0,0};
  double mixupfraction[kFelix_num_];
  double mixupfraction_sum[kFelix_num_] = {0,0,0,0,0,0,0,0};

  int Nmixup_sum[kFelix_num_] = {0,0,0,0,0,0,0,0};
  int pre_allhit_sum[kFelix_num_] = {0,0,0,0,0,0,0,0};
  double Nmixup_ave[kFelix_num_] = {0,0,0,0,0,0,0,0};
  double pre_allhit_ave[kFelix_num_] = {0,0,0,0,0,0,0,0};
  
  //double thisclonefraction[kFelix_num_];
  double copyfraction[kFelix_num_];
  double copyfraction_sum[kFelix_num_] = {0,0,0,0,0,0,0,0};
  int NmixcopyEv[kFelix_num_] = {0,0,0,0,0,0,0,0};
  
  double mixupfraction_ave[kFelix_num_];
  double Mixevent[kFelix_num_];
  double err[kFelix_num_];  
  double copyfraction_ave[kFelix_num_];
  double Mixcloevent[kFelix_num_];
  double thisclone_ave[kFelix_num_];

  uint64_t first_bcofull;

  int DEFAULT_BCO_VALUE[kFelix_num_];

  TGraph *g_evfraction;
  TGraph *g_cloevfraction ;
  TGraph *g_hitfraction;
  TGraph *g_copyfraction;

  std::map<int, int> premap_hit;
  
  TH1F *h_allmulti_[kFelix_num_];
  TH1F *h_allclone_[kFelix_num_];

  TH1F *h_mixupmulti[kFelix_num_][divimul];
  TH1F *h_divmul[kFelix_num_][divimul];

  TH2F *h_vsprefull_bco[kFelix_num_]; 
  TH2F *h_vsfull_bco[kFelix_num_];
  TH1F *h_prefull_bco[kFelix_num_];
  TH1F *h_full_bco[kFelix_num_];
  TH1F *h_prefull_bco_all[kFelix_num_];
  TH2F *h_vsprefull_bco_all[kFelix_num_];
  TH1F *h_interval;
  TH1F *h_mixinterval;
  TH1F *h_divinter;
  TH1F *h_bcofull_7;
  TH1F *h_bco;

  TH1F *h_mixup[kFelix_num_];
  TH2F *h_prevsNmix[kFelix_num_];
  TH1F *h_nocopyhit[kFelix_num_];
  TH1F *h_copyhit[kFelix_num_];
  TH1F *h_mixcopy[kFelix_num_];
  TH2F *h_mixvscopy[kFelix_num_];

  TFile *tf_bcopeak_[kFelix_num_];
  TH1D *hbco[kFelix_num_];
  TH1D *hbco_sub[kFelix_num_];
  TH2D *h2_bco_felix[kFelix_num_];
  TH2D *h2_bco_felix_sub[kFelix_num_];
  TH1* hbcohist[kFelix_num_];
  TH1* hbcohist2[kFelix_num_];
  //TFile *tf_hotchan_;

  //Mixup fraction hist
  TH2F *h_hitfra[kFelix_num_];
  TH2F *h_bghit[kFelix_num_];
  TH1F *h_NmixEv;
  TH1F *h_AllEv;

  std::ofstream f_felixpeak;
  std::ofstream ffraction;
  std::ofstream fNhit[kFelix_num_];
  std::ofstream f_hotchan;
  TFile *fgraph;

  //for Histgram parameter
  //int bin=1000;
  //int bin=400;
  int bin=8000;
  int bin3=100;
  Long64_t bit=0xFFFFFFFFFF;
  Long64_t bin2 = 10000;
  //std::string p;
  //std::string Au;
  
  void DrawHists();

  //void Mixupfraction();
  
  void GetBcopeak();

  void Readpeak();

  void GetBcomyself();

  void Hotchancut();
};
#endif
