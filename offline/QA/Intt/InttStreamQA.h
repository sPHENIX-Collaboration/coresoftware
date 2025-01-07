#ifndef INTTSTREAMQA_H__
#define INTTSTREAMQA_H__

#include <TObject.h>
#include <fun4all/SubsysReco.h>

#include <string>

//#include <intt/InttMapping.h>

/// Class declarations for use in the analysis module
class PHCompositeNode;
class TH1;
class TH2;
class InttRawHitContainer;
/// Definition of this analysis module class
class InttStreamQA : public SubsysReco
{
 public:
  /*!
  @brief Constructor
  @param name Used to initialize Subsysreco
  @param fname It's assigned to fname_
  */
  InttStreamQA(const std::string& name = "InttStreamQA");

  // Destructor
  virtual ~InttStreamQA() = default;

  /// SubsysReco initialize processing method
  int InitRun(PHCompositeNode*);

  /// SubsysReco event processing method
  int process_event(PHCompositeNode*);

  /// SubsysReco end processing method
  int End(PHCompositeNode*);

 private:
  void createHistos();
  std::string getHistoPrefix() const;
  std::vector<InttRawHitContainer*> m_rawhit_containers;

  TH2* h_bco[8]{nullptr};
  TH2* h_hit[8]{nullptr};

  TH1* h_bco_felix[8]{nullptr};         // FPHX bco
                                        //   TH1* h_bco_all{nullptr};
  TH1* h_bcogl1diff_felix[8]{nullptr};  // bcofull - gl1bco for all hits

  TH1* h_bcoreco_diff[8]{nullptr};
  TH1* h_bcoreco_evt_diff[8]{nullptr};
  TH1* h_bcoreco_evt_diff_all{nullptr};

  TH1* h_bcorecointt_diff[8]{nullptr};
  TH1* h_bcointtgl1_diff{nullptr};

  //   TH1* h_bunch[8]{nullptr};
  TH1* h_bunch_all{nullptr};
  TH1* h_bunch_gl1{nullptr};

  //   TH1* h_bunch_strb[8]{nullptr};
  TH2* h_bunch_evt_bcodiff[8]{nullptr};
  TH2* h_bunch_bco[8]{nullptr};

  //   TH1* h_bcoprediff[8]{nullptr};
};

#endif
