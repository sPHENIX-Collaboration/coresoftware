#ifndef INTTHITMAPV1_INTTHITMAP_H
#define INTTHITMAPV1_INTTHITMAP_H

#include <intt/InttFeeMap.h>
#include <intt/InttFeeMapv1.h>
#include <intt/InttMap.h>

#include <fun4all/SubsysReco.h>

#include <TObject.h>

#include <cstdint>
#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;
class TH2;

class InttHitMap : public SubsysReco
{
 public:
  InttHitMap(const std::string& name = "InttHitMap", const std::string& fname = "outfile.root", int nevent = 10000);

  virtual ~InttHitMap();

  int Init(PHCompositeNode*);

  int InitRun(PHCompositeNode*);

  /// SubsysReco event processing method
  int process_event(PHCompositeNode*);

  /// SubsysReco end processing method
  int End(PHCompositeNode*);

  bool isBCOcutON_ = false;
  bool isBCOPeak(int felix, int ladder, int bco, uint64_t bcofull);
  void SetBCOcut(const bool flag) { isBCOcutON_ = flag; }

  bool isBeam_ = true;
  void IsBeam(const bool flag) { isBeam_ = flag; }
  int SetBCOFile(const std::string& bcofile);
  int SetFeeMapFile(const std::string& feemapfile);
  InttFeeMapv1 fee_map;
  bool FillHitMap(int felix, int moudle, int barrel, int chip, int chan);
  void SetRunNumber(const int runnum) { runnumber_ = runnum; }
  ///////////////////////////////////

 private:
  TFile* inBCOFile_{nullptr};
  TFile* outFile_{nullptr};
  TTree* outTree_{nullptr};
  int nevents_{0};
  int ievent_{0};
  uint32_t total_event_{0};
  int runnumber_{0};
  TH2* h2_AllMap_[8][14]{};
  TH2* h2_bco_cut_[8]{};
  bool IsCloneHit_[8][14][26][128]{};

  std::string outfname_;
  std::string m_InttRawNodeName{"INTTRAWHIT"};

  struct Half_Chip
  {
    int _felix_id_;
    int _module_id_;
    int _chip_id_;
  };
};

#endif
