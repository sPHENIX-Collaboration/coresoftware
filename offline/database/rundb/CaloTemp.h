#ifndef RUNDB_CALOTEMP_H
#define RUNDB_CALOTEMP_H

#include <fun4all/SubsysReco.h>

#include <string>

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TTree;
class TProfile2D;
class RunHeader;

class CaloTemp : public SubsysReco
{
 public:
  //! constructor
  CaloTemp(const std::string &name = "CaloTemp", const std::string &filename = "MyNtuple.root");

  //! destructor
  virtual ~CaloTemp();

  //! full initialization
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  void Detector(const std::string &name) { detector = name; }
  void outputTempHists(const bool outputhist) { temphist = outputhist; }
  int getRunTime();
  int getTempHist();

 protected:
  Fun4AllHistoManager *hm{nullptr};
  TFile *outfile{nullptr};
  TProfile2D *h_calo_temp{nullptr};
  TProfile2D *h_calo_temp2{nullptr};
  RunHeader *runheader{nullptr};
  TTree *tree{nullptr};
  int runnumber{-9999};
  int year{-9999};
  int month{-9999};
  int day{-9999};
  int hour{-9999};
  int minute{-9999};
  int second{-9999};
  bool temphist{true};
  std::string detector{"HCALIN"};
  std::string outfilename;
  std::string runtime;
};

#endif
