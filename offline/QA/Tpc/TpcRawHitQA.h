// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCRAWHITQA_H
#define TPCRAWHITQA_H

#include <tpc/TpcMap.h>

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class TH1;
class TH2;
class TpcRawHitContainer;

class TpcRawHitQA : public SubsysReco
{
 public:
  TpcRawHitQA(const std::string& name = "TpcRawHitQA");

  ~TpcRawHitQA() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int EndRun(const int runnumber) override;

  void setPresample(bool value)
  {
    longPresample = value;
  }

 private:
  void createHistos();
  std::string getHistoPrefix() const;

  TpcMap M;

  std::vector<TpcRawHitContainer*> rawhitcont_vec;

  TH1* h_nhits_sectors[24]{nullptr};
  TH2* h_nhits_sectors_fees[24]{nullptr};
  TH1* h_nhits_sectors_laser[24]{nullptr};
  TH2* h_nhits_sectors_fees_laser[24]{nullptr};
  TH1* h_nhits_sam[24][3]{{nullptr}};
  TH1* h_adc[24][3]{{nullptr}};
  TH2* h_xy_N{nullptr};
  TH2* h_xy_S{nullptr};

  int FEE_R[26]{2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};
  int FEE_map[26]{4, 5, 0, 2, 1, 11, 9, 10, 8, 7, 6, 0, 1, 3, 7, 6, 5, 4, 3, 2, 0, 2, 1, 3, 5, 4};

  bool longPresample = false;
};

#endif  // TPCRAWHITQA_H
