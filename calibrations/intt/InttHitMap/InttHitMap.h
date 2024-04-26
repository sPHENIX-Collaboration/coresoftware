#ifndef INTTHITMAP_INTTHITMAP_H
#define INTTHITMAP_INTTHITMAP_H

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
  int SetBCOFile(const std::string &bcofile);
  int SetFeeMapFile(const std::string &feemapfile);
  InttFeeMapv1 fee_map;
  bool FillHitMap(int felix, int moudle, int barrel, int chip, int chan);

  ///////////////////////////////////

 private:
  TFile* inBCOFile_{nullptr};
  TFile* outFile_{nullptr};
  int nevents_{0};
  int ievent_{0};
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
  ////////////////////////////
  // List of Half entry chips//
  ////////////////////////////
  std::vector<Half_Chip> half_chips =  // Chip number 0~25
      {
          // Felix 0
          {0, 7, 14},
          // Felix 2
          {2, 9, 15},
          // Felix 3
          {3, 13, 20},
          {3, 13, 22},
          {3, 13, 24},
          // Felix 7
          {7, 0, 0},
          {7, 0, 1},
          {7, 0, 2},
          {7, 0, 3},
          {7, 0, 4},
          {7, 0, 13},
          {7, 0, 14},
          {7, 0, 15},
          {7, 0, 16},
          {7, 0, 17},
          {7, 0, 6},
          {7, 0, 8},
          {7, 0, 10},
          {7, 0, 12},
          {7, 0, 19},
          {7, 0, 21},
          {7, 0, 23},
          {7, 0, 25},
          {7, 1, 0},
          {7, 1, 1},
      };
};

#endif
