#ifndef TRIGGER_CALOTRIGGEREMULATOR_H
#define TRIGGER_CALOTRIGGEREMULATOR_H

//#include "LL1Outv1.h"
#include "TriggerDefs.h"

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <map>
#include <vector>

// Forward declarations
class CDBHistos;
class CDBTTree;
class TriggerPrimitive;
class TriggerPrimitiveContainer;
class LL1Out;
class Event;
class TowerInfoContainer;
class CaloPacketContainer;
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TH1;
class TH1I;
class TNtuple;
class TTree;
class TProfile;
class TEfficiency;

class CaloTriggerEmulator : public SubsysReco
{
 public:
  //! constructor
  explicit CaloTriggerEmulator(const std::string &name);

  //! destructor
  ~CaloTriggerEmulator() override = default;

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! full initialization
  int InitRun(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;
  int process_sim();
  int process_offline();

  //! end of run method
  int End(PHCompositeNode *) override;

  //! reset variables
  int ResetEvent(PHCompositeNode *) override;

  //! Get Nodes
  void GetNodes(PHCompositeNode *);

  //! Create Nodes
  void CreateNodes(PHCompositeNode *);

  //! MakePrimitives
  int process_waveforms();

  //! MakeTriggerOutput
  int process_primitives();

  int process_organizer();

  int process_trigger();

  unsigned int getBits(unsigned int sum, TriggerDefs::TriggerId tid);

  int Download_Calibrations();

  //! Set TriggerType
  void setTriggerType(const std::string &name);
  void setTriggerType(TriggerDefs::TriggerId triggerid);

  void setOptMaskFile(const std::string &filename) { m_optmask_file = filename; }

  void setEmcalLUTFile(const std::string &filename) { m_emcal_lutname = filename; }
  void setHcalinLUTFile(const std::string &filename) { m_hcalin_lutname = filename; }
  void setHcaloutLUTFile(const std::string &filename) { m_hcalout_lutname = filename; }

  void useMax(bool max) { m_use_max = max; }

  void useEMCALDefaultLUT(bool def) { m_default_lut_emcal = def; }
  void useHCALINDefaultLUT(bool def) { m_default_lut_hcalin = def; }
  void useHCALOUTDefaultLUT(bool def) { m_default_lut_hcalout = def; }

  void setTriggerSample(int s) { m_trig_sample = s; }
  void setTriggerDelay(int d) { m_trig_sub_delay = d + 1; }

  void useHCALIN(bool use);
  void useHCALOUT(bool use);
  void useEMCAL(bool use);

  void setNSamples(int nsamples) { m_nsamples = nsamples; }
  void setThreshold(int threshold) { m_threshold = threshold; }
  void setJetThreshold(int t1, int t2, int t3, int t4)
  {
    m_threshold_jet[0] = t1;
    m_threshold_jet[1] = t2;
    m_threshold_jet[2] = t3;
    m_threshold_jet[3] = t4;
    return;
  }
  void setPhotonThreshold(int t1, int t2, int t3, int t4)
  {
    m_threshold_photon[0] = t1;
    m_threshold_photon[1] = t2;
    m_threshold_photon[2] = t3;
    m_threshold_photon[3] = t4;
    return;
  }

  void setPairThreshold(int t1, int t2, int t3, int t4)
  {
    m_threshold_pair[0] = t1;
    m_threshold_pair[1] = t2;
    m_threshold_pair[2] = t3;
    m_threshold_pair[3] = t4;
    return;
  }


  bool CheckFiberMasks(TriggerDefs::TriggerPrimKey key);
  void LoadFiberMasks();
  void SetIsData(bool isd) { m_isdata = isd; }
  bool CheckChannelMasks(TriggerDefs::TriggerSumKey key);

  void identify();

 private:
  std::string m_ll1_nodename;
  std::string m_prim_nodename;
  std::string m_waveform_nodename;

  std::string m_optmask_file;

  std::string m_emcal_lutname;
  std::string m_hcalin_lutname;
  std::string m_hcalout_lutname;

  //! Trigger Type
  std::string m_trigger{"NONE"};

  TriggerDefs::TriggerId m_triggerid = TriggerDefs::TriggerId::noneTId;

  bool m_use_max{true};
  bool m_do_hcalin{false};
  bool m_do_hcalout{false};
  bool m_do_emcal{false};


  bool m_default_lut_hcalin{false};
  bool m_default_lut_hcalout{false};
  bool m_default_lut_emcal{false};

  bool m_force_hcalin{false};
  bool m_force_hcalout{false};
  bool m_force_emcal{false};

  //! Waveform conatiner
  Event *m_event{nullptr};
  TowerInfoContainer *m_waveforms_hcalin{nullptr};
  TowerInfoContainer *m_waveforms_hcalout{nullptr};
  TowerInfoContainer *m_waveforms_emcal{nullptr};
  CaloPacketContainer *m_hcal_packets{nullptr};
  CaloPacketContainer *m_emcal_packets{nullptr};

  
  //! LL1 Out
  LL1Out *m_ll1out_photon{nullptr};

  LL1Out *m_ll1out_pair{nullptr};

  LL1Out *m_ll1out_jet{nullptr};

  TriggerPrimitiveContainer *m_primitives_photon{nullptr};

  TriggerPrimitiveContainer *m_primitives_pair{nullptr};

  TriggerPrimitiveContainer *m_primitives_jet{nullptr};

  TriggerPrimitiveContainer *m_primitives_hcalin{nullptr};

  TriggerPrimitiveContainer *m_primitives_hcalout{nullptr};

  TriggerPrimitiveContainer *m_primitives_emcal{nullptr};

  TriggerPrimitiveContainer *m_primitives_hcal_ll1{nullptr};

  TriggerPrimitiveContainer *m_primitives_emcal_ll1{nullptr};

  TriggerPrimitiveContainer *m_primitives_emcal_2x2_ll1{nullptr};

  /* std::map<unsigned int, TH1I> h_mbd_charge_lut; */
  /* std::map<unsigned int, TH1I> h_mbd_time_lut; */
  /* std::map<unsigned int, TH1I> h_mbd_slewing_lut; */

  unsigned int m_l1_hcal_table[4096]{};
  unsigned int m_l1_adc_table[1024]{};
  unsigned int m_l1_8x8_table[1024]{};
  unsigned int m_l1_slewing_table[4096]{};

  std::map<unsigned int, TH1I*> h_emcal_lut{};
  std::map<unsigned int, TH1I*> h_hcalin_lut{};
  std::map<unsigned int, TH1I*> h_hcalout_lut{};

  CDBTTree *cdbttree_adcmask{nullptr};
  CDBHistos *cdbttree_emcal{nullptr};
  CDBHistos *cdbttree_hcalin{nullptr};
  CDBHistos *cdbttree_hcalout{nullptr};

  std::map<unsigned int, std::vector<unsigned int> > m_peak_sub_ped_emcal{};
  std::map<unsigned int, std::vector<unsigned int> > m_peak_sub_ped_hcalin{};
  std::map<unsigned int, std::vector<unsigned int> > m_peak_sub_ped_hcalout{};

  //! Verbosity.
  int m_nevent;
  int m_photon_npassed;
  int m_jet_npassed;
  int m_pair_npassed;

  int m_n_sums;
  int m_n_primitives;
  int m_trig_sub_delay;
  int m_trig_sample{-1};

  unsigned int m_threshold{1};
  unsigned int m_threshold_jet[4] = {0};
  unsigned int m_threshold_pair[4] = {0};
  unsigned int m_threshold_photon[4] = {0};

  int m_isdata{1};
  int m_useoffline{false};
  int m_nsamples = 16;
  int m_idx{8};

  int m_packet_low_hcalout = 8001;
  int m_packet_high_hcalout = 8008;
  int m_packet_low_hcalin = 7001;
  int m_packet_high_hcalin = 7008;
  int m_packet_low_emcal = 6001;
  int m_packet_high_emcal = 6128;

  std::string m_fieldname{""};

  std::vector<unsigned int> m_masks_fiber{};
  std::vector<unsigned int> m_masks_channel{};
  std::map<TriggerDefs::DetectorId, int> m_prim_map{};
  std::map<TriggerDefs::TriggerId, int> m_prim_ll1_map{};
  std::map<TriggerDefs::TriggerId, std::vector<std::string>> m_det_map{};
};

#endif
