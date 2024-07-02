// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MVTXRAWHITQA_H
#define MVTXRAWHITQA_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <ffarawobjects/MvtxRawHitContainer.h>
#include <ffarawobjects/MvtxRawHit.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>

#include <string>
#include <vector>

class MvtxRawHitQA : public SubsysReco
{
 public:
  MvtxRawHitQA(const std::string& name = "MvtxRawHitQA");

  ~MvtxRawHitQA() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;
  
 private:
  void createHistos();
  std::string getHistoPrefix() const;

  MvtxRawHitContainer* rawhitcont{nullptr};

  TH1* h_nhits_per_chip_layper0{nullptr};
  TH1* h_nhits_per_chip_layper1{nullptr};
  TH1* h_nhits_per_chip_layper2{nullptr};
  TH1* h_chipocc_layper0{nullptr};
  TH1* h_chipocc_layper1{nullptr};
  TH1* h_chipocc_layper2{nullptr};
  TH1* h_bco{nullptr};
  TH1* h_strobe_bc{nullptr};
  TH1* h_chip_bc{nullptr};

};

#endif  // MVTXRAWHITQA_H
