// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MVTXRAWHITQA_H
#define MVTXRAWHITQA_H

#include <trackbase/TrkrDefs.h>

#include <ffarawobjects/MvtxRawHitContainer.h>
#include <ffarawobjects/MvtxRawHit.h>

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class TH1;
class TH2;

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
  std::vector<MvtxRawHitContainer*> m_rawhit_containers;

  TH1* h_nhits_layer0{nullptr};
  TH1* h_nhits_layer1{nullptr};
  TH1* h_nhits_layer2{nullptr};
  TH1* h_bco{nullptr};
  TH1* h_strobe_bc{nullptr};
  TH1* h_chip_bc{nullptr};
  TH2* h_nhits_stave_chip_layer0{nullptr};
  TH2* h_nhits_stave_chip_layer1{nullptr};
  TH2* h_nhits_stave_chip_layer2{nullptr};

};

#endif  // MVTXRAWHITQA_H
