#ifndef TRIGGER_MBDTRIGGEREMULATOR_H
#define TRIGGER_MBDTRIGGEREMULATOR_H

#include "LL1Outv1.h"
#include "TriggerDefs.h"
#include "TriggerPrimitiveContainerv1.h"
#include "TriggerPrimitivev1.h"

#include <fun4all/SubsysReco.h>

#include <TEfficiency.h>
#include <TH2.h>
#include <TProfile.h>
#include <TTree.h>
#include <TNtuple.h>

#include <cstdint>

// Forward declarations
class CDBHistos;
class TriggerPrimitive;
class TriggerPrimitiveContainer;
class LL1Out;
class Event;
class TowerInfoContainer;
class CaloPacketContainer;
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TProfile;
class TEfficiency;
class TH2D;

class MBDTriggerEmulator : public SubsysReco
{
 public:
  //! constructor
  explicit MBDTriggerEmulator(const std::string &name);

  //! destructor
  ~MBDTriggerEmulator() override = default;

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
  int process_raw();

  //! MakeTriggerOutput
  int process_primitives();

  int process_trigger();

  int Download_Calibrations();

  void useMax(bool max) { m_use_max= max; }

  void setTriggerSample(int s) { m_trig_sample = s; }
  void setTriggerDelay(int d) { m_trig_sub_delay = d + 1; }

  void setNSamples(int nsamples) { m_nsamples = nsamples; }

  void identify();

  void setChargeLUTFile(const std::string &filename) { m_mbd_charge_lutname = filename; }
  void setTimeLUTFile(const std::string &filename) { m_mbd_time_lutname = filename; }
  void setSlewingLUTFile(const std::string &filename) { m_mbd_slewing_lutname = filename; }

 protected:
  std::string m_ll1_nodename;
  std::string m_prim_nodename;
  std::string m_waveform_nodename;

  std::string m_mbd_charge_lutname;
  std::string m_mbd_time_lutname;  
  std::string m_mbd_slewing_lutname;

  TriggerDefs::TriggerId m_triggerid = TriggerDefs::TriggerId::mbdTId;

  bool m_use_max{true};

  bool m_default_lut_mbd{false};

  CaloPacketContainer *m_waveforms_mbd{nullptr};
  
  Event *m_event{nullptr};
  //! LL1 Out
  LL1Out *m_ll1out_mbd{nullptr};

  TriggerPrimitiveContainer *m_primitives_mbd{nullptr};

  /* std::map<unsigned int, TH1I> h_mbd_charge_lut; */
  /* std::map<unsigned int, TH1I> h_mbd_time_lut; */
  /* std::map<unsigned int, TH1I> h_mbd_slewing_lut; */

  std::map<unsigned int, TH1I*> h_mbd_charge_lut{};
  std::map<unsigned int, TH1I*> h_mbd_time_lut{};
  std::map<unsigned int, TH1I*> h_mbd_slewing_lut{};

  CDBHistos *cdbttree_mbd_charge{nullptr};
  CDBHistos *cdbttree_mbd_time{nullptr};
  CDBHistos *cdbttree_mbd_slewing{nullptr};

  //! Trigger primitives
  unsigned int m_trig_charge[8]{};
  unsigned int m_trig_nhit;
  unsigned int m_trig_time[4]{};

  //! Trigger primitives
  unsigned int m2_trig_charge[4][8]{};
  unsigned int m2_trig_nhit[4]{};
  unsigned int m2_trig_time[4][4]{};

  //! Trigger ouputs
  std::vector<std::vector<unsigned int> *> m_word_mbd{};
  unsigned int m_out_tsum[2]{};
  unsigned int m_out_tavg[2]{};
  unsigned int m_out_trem[2]{};
  unsigned int m_out_nhit[2]{};
  unsigned int m_out_vtx_sub;
  unsigned int m_out_vtx_add;

  unsigned int m_nhit1, m_nhit2, m_timediff1, m_timediff2, m_timediff3;

  std::map<unsigned int, std::vector<unsigned int> > m_peak_sub_ped_mbd{};

  //! Verbosity.
  int m_nevent;
  int m_mbd_npassed;

  int m_n_sums;
  int m_n_primitives;
  int m_trig_sub_delay;
  int m_trig_sample{-1};

  unsigned int m_threshold{1};
  int m_isdata{1};
  int m_useoffline{false};
  int m_nsamples = 16;
  int m_idx{8};

};

#endif
