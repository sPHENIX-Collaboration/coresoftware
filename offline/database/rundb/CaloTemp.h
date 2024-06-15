#ifndef CALOTEMP_H__
#define CALOTEMP_H__

#include <fun4all/SubsysReco.h>
#include <vector>
// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TH2F;
class TH1F;
class TH1;
class TProfile2D;
class RunHeader;

class CaloTemp : public SubsysReco
{
 public:
  //! constructor
  CaloTemp(const std::string &name = "CaloTemp", const std::string &fname = "MyNtuple.root");

  //! destructor
  virtual ~CaloTemp();

  //! full initialization
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  void Detector(const std::string &name) { detector = name; }
  int getRunTime();
  int getTempHist();

 protected:
  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm = nullptr;
  TFile *outfile = nullptr;
  TProfile2D* h_calo_temp = nullptr;
  TProfile2D* h_calo_temp2 = nullptr;
  RunHeader* runheader = nullptr;
  int runnumber {0};
  std::string runtime;

};

#endif
