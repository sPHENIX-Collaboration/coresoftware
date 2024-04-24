#ifndef TRIGGERVALID_TRIGGERVALID_H
#define TRIGGERVALID_TRIGGERVALID_H

#include <fun4all/SubsysReco.h>
#include <string>
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

class TriggerValid : public SubsysReco
{
 public:
  //! constructor
  TriggerValid(const std::string& name = "TriggerValid", const std::string& fname = "MyNtuple.root");

  //! destructor
  virtual ~TriggerValid();

  //! full initialization
  int Init(PHCompositeNode*);

  //! event processing method
  int process_event(PHCompositeNode*);

  //! end of run method
  int End(PHCompositeNode*);

  int process_towers(PHCompositeNode*);
  int process_primitives(PHCompositeNode*);
  int process_ll1out(PHCompositeNode*);

  void Trigger(const std::string& name) { trigger = name; }
 
  void set_debug(bool debug) { m_debug = debug; }
 
 private:
  int Getpeaktime(TH1* h);

  bool m_debug{0};
  std::string trigger;
  std::string outfilename;
  Fun4AllHistoManager* hm{nullptr};
  TFile* outfile{nullptr};

  TH2F* h_emu_emcal_twr_frequency{nullptr};
  TH2F* h_emu_emcal_2x2_frequency{nullptr};
  TH2F* h_emu_emcal_8x8_frequency{nullptr};
  TH2F* h_emu_ihcal_twr_frequency{nullptr};
  TH2F* h_emu_ihcal_2x2_frequency{nullptr};
  TH2F* h_emu_ohcal_twr_frequency{nullptr};
  TH2F* h_emu_ohcal_2x2_frequency{nullptr};
  TH2F* h_emu_hcal_2x2_frequency{nullptr};

  TH2F* h_emu_jet_frequency{nullptr};
  TH2F* h_emu_photon_frequency{nullptr};

  TH2F* h_emcal_ll1_2x2_frequency{nullptr};
  TH2F* h_emcal_ll1_8x8_frequency{nullptr};
  TH2F* h_hcal_ll1_2x2_frequency{nullptr};
  TH2F* h_jet_ll1_frequency{nullptr};

  TProfile2D* h_emu_emcal_twr_avg_out{nullptr};
  TProfile2D* h_emu_emcal_2x2_avg_out{nullptr};
  TProfile2D* h_emu_emcal_4x4_avg_out{nullptr};
  TProfile2D* h_emu_emcal_8x8_avg_out{nullptr};
  TProfile2D* h_emu_ihcal_twr_avg_out{nullptr};
  TProfile2D* h_emu_ihcal_2x2_avg_out{nullptr};
  TProfile2D* h_emu_ohcal_twr_avg_out{nullptr};
  TProfile2D* h_emu_ohcal_2x2_avg_out{nullptr};
  TProfile2D* h_emu_hcal_2x2_avg_out{nullptr};

  TProfile2D* h_emu_jet_avg_out{nullptr};
  TProfile2D* h_emu_photon_avg_out{nullptr};

  TProfile2D* h_emcal_ll1_2x2_avg_out{nullptr};
  TProfile2D* h_emcal_ll1_4x4_avg_out{nullptr};
  TProfile2D* h_emcal_ll1_8x8_avg_out{nullptr};
  TProfile2D* h_hcal_ll1_2x2_avg_out{nullptr};
  TProfile2D* h_jet_ll1_avg_out{nullptr};

  TH2F* h_emcal_2x2_energy_lutsum{nullptr};
  TH2F* h_emcal_8x8_energy_lutsum{nullptr};

  TH2F* h_hcal_2x2_energy_lutsum{nullptr};

  TH2F* h_jet_energy_lutsum{nullptr};

  TProfile2D* h_match_emcal{nullptr};
  TProfile2D* h_match_emcal_ll1{nullptr};
  TProfile2D* h_match_hcal_ll1{nullptr};
  TProfile2D* h_match_jet_ll1{nullptr};
  int _eventcounter{0};
  int _range{1};

};

#endif
