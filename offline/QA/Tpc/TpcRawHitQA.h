// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCRAWHITQA_H
#define TPCRAWHITQA_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <ffarawobjects/TpcRawHitContainer.h>
#include <ffarawobjects/TpcRawHit.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>

#include <string>
#include <vector>

class TpcRawHitQA : public SubsysReco
{
 public:
  TpcRawHitQA(const std::string& name = "TpcRawHitQA");

  ~TpcRawHitQA() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;
  
 private:
  void createHistos();
  std::string getHistoPrefix() const;

  TpcRawHitContainer* rawhitcont{nullptr};

  TH1* h_nhits_sectors[24]{nullptr};
  TH2* h_nhits_sectors_fees[24]{nullptr};
  TH1* h_bco{nullptr};

};

#endif  // TPCRAWHITQA_H
