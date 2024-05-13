#ifndef TRIGGER_CALOTRIGGEREMULATOR_H
#define TRIGGER_CALOTRIGGEREMULATOR_H

#include "LL1Outv1.h"
#include "TriggerDefs.h"
#include "TriggerPrimitiveContainerv1.h"
#include "TriggerPrimitivev1.h"

#include <fun4all/SubsysReco.h>

#include <TEfficiency.h>
#include <TH2.h>
#include <TProfile.h>
#include <TTree.h>

#include <cstdint>

// Forward declarations
class CDBHistos;
class TriggerPrimitive;
class TriggerPrimitiveContainer;
class LL1Out;

class TowerInfoContainer;
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TProfile;
class TEfficiency;
class TH2D;

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

  unsigned int getBits(unsigned int sum);

  int Download_Calibrations();

  //! Set TriggerType
  void setTriggerType(const std::string &name);
  void setTriggerType(TriggerDefs::TriggerId triggerid);
  void setEmcalLUTFile(const std::string &filename) { m_emcal_lutname = filename; }
  void setHcalinLUTFile(const std::string &filename) { m_hcalin_lutname = filename; }
  void setHcaloutLUTFile(const std::string &filename) { m_hcalout_lutname = filename; }

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
  void setThreshold(int t1, int t2, int t3, int t4)
  {
    m_threshold_calo[0] = t1;
    m_threshold_calo[1] = t2;
    m_threshold_calo[2] = t3;
    m_threshold_calo[3] = t4;
    m_single_threshold = false;
    return;
  }

  bool CheckFiberMasks(TriggerDefs::TriggerPrimKey key);
  bool CheckChannelMasks(TriggerDefs::TriggerSumKey key);

  void identify();

 protected:
  std::string m_ll1_nodename;
  std::string m_prim_nodename;
  std::string m_waveform_nodename;

  std::string m_emcal_lutname;
  std::string m_hcalin_lutname;
  std::string m_hcalout_lutname;

  //! Trigger Type
  std::string m_trigger{"NONE"};

  TriggerDefs::TriggerId m_triggerid = TriggerDefs::TriggerId::noneTId;

  bool m_do_hcalin{false};
  bool m_do_hcalout{false};
  bool m_do_emcal{false};
  bool m_do_mbd{false};

  bool m_default_lut_hcalin{false};
  bool m_default_lut_hcalout{false};
  bool m_default_lut_emcal{false};
  bool m_default_lut_mbd{false};

  bool m_force_hcalin{false};
  bool m_force_hcalout{false};
  bool m_force_emcal{false};
  bool m_force_mbd{false};

  //! Waveform conatiner
  TowerInfoContainer *m_waveforms_hcalin{nullptr};
  TowerInfoContainer *m_waveforms_hcalout{nullptr};
  TowerInfoContainer *m_waveforms_emcal{nullptr};
  TowerInfoContainer *m_waveforms_mbd{nullptr};

  //! LL1 Out
  LL1Out *m_ll1out{nullptr};
  TriggerPrimitiveContainer *m_primitives{nullptr};

  TriggerPrimitiveContainer *m_primitives_hcalin{nullptr};

  TriggerPrimitiveContainer *m_primitives_hcalout{nullptr};

  TriggerPrimitiveContainer *m_primitives_emcal{nullptr};

  TriggerPrimitiveContainer *m_primitives_hcal_ll1{nullptr};

  TriggerPrimitiveContainer *m_primitives_emcal_ll1{nullptr};


  unsigned int m_l1_adc_table[1024]{};
  unsigned int m_l1_adc_table_time[1024]{};
  unsigned int m_l1_slewing_table[4096]{};
  unsigned int m_l1_hcal_table[4096]{};


  std::map<unsigned int, TH1I*> h_emcal_lut;
  std::map<unsigned int, TH1I*> h_hcalin_lut;
  std::map<unsigned int, TH1I*> h_hcalout_lut;

  CDBHistos *cdbttree_emcal{nullptr};
  CDBHistos *cdbttree_hcalin{nullptr};
  CDBHistos *cdbttree_hcalout{nullptr};

  std::string m_fieldname_emcal;
  std::string m_calibName_emcal;

  std::string m_fieldname_hcalin;
  std::string m_calibName_hcalin;

  std::string m_fieldname_hcalout;
  std::string m_calibName_hcalout;

  //! Trigger primitives
  unsigned int m_trig_charge[8]{};
  unsigned int m_trig_nhit;
  unsigned int m_trig_time[4]{};

  //! Trigger primitives
  unsigned int m2_trig_charge[4][8]{};
  unsigned int m2_trig_nhit[4]{};
  unsigned int m2_trig_time[4][4]{};

  //! Trigger ouputs
  std::vector<std::vector<unsigned int> *> m_word_mbd;
  unsigned int m_out_tsum[2]{};
  unsigned int m_out_tavg[2]{};
  unsigned int m_out_trem[2]{};
  unsigned int m_out_nhit[2]{};
  unsigned int m_out_vtx_sub;
  unsigned int m_out_vtx_add;

  unsigned int m_nhit1, m_nhit2, m_timediff1, m_timediff2, m_timediff3;

  std::map<unsigned int, std::vector<unsigned int> > m_peak_sub_ped_emcal;
  std::map<unsigned int, std::vector<unsigned int> > m_peak_sub_ped_mbd;
  std::map<unsigned int, std::vector<unsigned int> > m_peak_sub_ped_hcalin;
  std::map<unsigned int, std::vector<unsigned int> > m_peak_sub_ped_hcalout;

  //! Verbosity.
  int m_nevent;
  int m_npassed;
  int m_n_sums;
  int m_n_primitives;
  int m_trig_sub_delay;
  int m_trig_sample{-1};

  bool m_single_threshold{true};
  unsigned int m_threshold{1};
  unsigned int m_threshold_calo[4] = {0};
  int m_isdata{1};
  int m_nsamples = 31;
  int m_idx{3};

  std::vector<unsigned int> m_masks_fiber;
  std::vector<unsigned int> m_masks_channel;
  std::map<TriggerDefs::DetectorId, int> m_prim_map;
  std::map<TriggerDefs::TriggerId, int> m_prim_ll1_map;
  std::map<TriggerDefs::TriggerId, std::vector<std::string>> m_det_map;
};

#endif
