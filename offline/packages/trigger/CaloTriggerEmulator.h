#ifndef _CALOTRIGGEREMULATOR_H__
#define _CALOTRIGGEREMULATOR_H__


#include "TriggerPrimitivev1.h"
#include "TriggerPrimitiveContainerv1.h"
#include "TriggerDefs.h"
#include <cstdint>

#include <fun4all/SubsysReco.h>
#include <TTree.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TH2.h>
#include "LL1Outv1.h"

// Forward declarations
class CDBTTree;
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
  explicit CaloTriggerEmulator(const std::string& name, const std::string& filename);

  //! destructor
  ~CaloTriggerEmulator();

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! full initialization
  int InitRun(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  //! end of run method
  int End(PHCompositeNode *) override;

  //! reset variables
  int  ResetEvent(PHCompositeNode *) override ;

  //! Get Nodes
  void GetNodes(PHCompositeNode *);

  //! Create Nodes
  void CreateNodes(PHCompositeNode *);

  //! MakePrimitives
  int process_waveforms();

  //! MakeTriggerOutput
  int process_primitives();

  int process_trigger();

  int Download_Calibrations();

  //! Set TriggerType
  void setTriggerType(const std::string &name);
  void setEmcalScale(float calib){m_emcal_GeV_per_lut = calib;}
  void setEmcalOffset(float off){m_emcal_lut_offset = off;}
  void setEmcalFloor(float floor){m_emcal_lut_floor = floor;}

  void setNSamples(int nsamples) { m_nsamples = nsamples; }
  void setThreshold(int threshold) { _m_threshold = threshold; }

  bool CheckFiberMasks(TriggerDefs::TriggerPrimKey key);
  bool CheckChannelMasks(TriggerDefs::TriggerSumKey key);

  void identify();
 protected:
  std::string outfilename;
  std::string _ll1_nodename;
  std::string _prim_nodename;
  std::string _waveform_nodename;

  Fun4AllHistoManager *hm = nullptr;
  TFile *outfile = nullptr;

  //!Trigger Type
  std::string _trigger{"NONE"};

  TriggerDefs::TriggerId _triggerid = TriggerDefs::TriggerId::noneTId;

  bool _do_hcalin{false};
  bool _do_hcalout{false};
  bool _do_emcal{false};
  bool _do_mbd{false};

  //! Waveform conatiner
  TowerInfoContainer *_waveforms_hcalin = nullptr;
  TowerInfoContainer *_waveforms_hcalout = nullptr;
  TowerInfoContainer *_waveforms_emcal = nullptr;
  TowerInfoContainer *_waveforms_mbd = nullptr;

  //! LL1 Out
  LL1Out *_ll1out = nullptr;
  TriggerPrimitiveContainer *_primitives = nullptr;

  TriggerPrimitiveContainer *_primitives_hcalin = nullptr;

  TriggerPrimitiveContainer *_primitives_hcalout = nullptr;

  TriggerPrimitiveContainer *_primitives_emcal = nullptr;

  TriggerPrimitive *_primitive = nullptr;

  std::vector<unsigned int> *_sum  = nullptr;
  std::vector<unsigned int> *_bits = nullptr;

  TProfile *avg_primitive  = nullptr;
  TH2D *peak_primitive = nullptr;
  TH2D *primitives = nullptr;
  TH2D *jet_trigger_fire_map = nullptr;
  TH2D *trigger_fire_map = nullptr;
  TH2D *h2_line_up = nullptr;
  TH1D *h_nhit = nullptr;
  TH1D *h_mbd_time = nullptr;
  TH1D *h_mbd_charge = nullptr;

  std::vector<TProfile*> v_avg_primitive_photon;
  std::vector<TH2D*> v_peak_primitive_photon;
  std::vector<TH2D*> v_primitives_photon;
  std::vector<TH2D*> v_trigger_fire_map_photon;

  std::vector<TProfile*> v_avg_primitive_jet;
  std::vector<TH2D*> v_peak_primitive_jet;
  std::vector<TH2D*> v_primitives_jet;
  std::vector<TH2D*> v_trigger_fire_map_jet;

  std::vector<TProfile*> v_avg_primitive_emcal;
  std::vector<TH2D*> v_peak_primitive_emcal;
  std::vector<TH2D*> v_primitives_emcal;
  std::vector<TH2D*> v_trigger_fire_map_emcal;

  std::vector<TProfile*> v_avg_primitive_hcalin;
  std::vector<TH2D*> v_peak_primitive_hcalin;
  std::vector<TH2D*> v_primitives_hcalin;
  std::vector<TH2D*> v_trigger_fire_map_hcalin;

  std::vector<TProfile*> v_avg_primitive_hcalout;
  std::vector<TH2D*> v_peak_primitive_hcalout;
  std::vector<TH2D*> v_primitives_hcalout;
  std::vector<TH2D*> v_trigger_fire_map_hcalout;

  std::vector<TH1D*> v_nhit;
  std::vector<TH2D*> v_line_up;
  std::map<std::string, TH1D*> v_mbd_charge;
  std::map<std::string, TH1D*> v_mbd_time;

  //! Lookup tables
  unsigned int m_l1_adc_table[1024];
  unsigned int m_l1_adc_table_time[1024];
  unsigned int m_l1_slewing_table[4096];
  unsigned int m_l1_hcal_table[4096];
  CDBTTree *cdbttree_emcal = nullptr;
  CDBTTree *cdbttree_hcalin = nullptr;
  CDBTTree *cdbttree_hcalout = nullptr;

  std::string m_fieldname_emcal;
  std::string m_calibName_emcal;

  std::string m_fieldname_hcalin;
  std::string m_calibName_hcalin;

  std::string m_fieldname_hcalout;
  std::string m_calibName_hcalout;

  //! Trigger primitives
  unsigned int m_trig_charge[8];
  unsigned int m_trig_nhit;
  unsigned int m_trig_time[4];


  //! Trigger primitives
  unsigned int m2_trig_charge[4][8];
  unsigned int m2_trig_nhit[4];
  unsigned int m2_trig_time[4][4];

  std::vector<std::vector<unsigned int>*> _sum_mbd;

  //! Trigger ouputs
  std::vector<std::vector<unsigned int>*> _word_mbd;
  unsigned int m_out_tsum[2];
  unsigned int m_out_tavg[2];
  unsigned int m_out_trem[2];
  unsigned int m_out_nhit[2];
  unsigned int m_out_vtx_sub;
  unsigned int m_out_vtx_add;

  unsigned int m_nhit1, m_nhit2, m_timediff1, m_timediff2, m_timediff3;

  float m_emcal_GeV_per_lut{0.0};
  float m_emcal_lut_offset{1.0};
  float m_emcal_lut_floor{4};
  std::map<unsigned int, std::vector<int>*> m_peak_sub_ped_emcal;
  std::map<unsigned int, std::vector<int>*> m_peak_sub_ped_mbd;
  std::map<unsigned int, std::vector<int>*> m_peak_sub_ped_hcalin;
  std::map<unsigned int, std::vector<int>*> m_peak_sub_ped_hcalout;
  //! Verbosity.
  int _nevent;
  int _npassed;
  int _n_sums;
  int _n_primitives;
  int _m_trig_sub_delay;
  unsigned int _m_threshold;
  int m_isdata;
  int m_nsamples = 31;
  int _idx;

  std::vector<unsigned int> _masks_fiber;
  std::vector<unsigned int> _masks_channel;
  std::map<TriggerDefs::DetectorId, int> _m_prim_map;
  std::map<TriggerDefs::TriggerId, int> _m_prim_ll1_map;
  std::map<TriggerDefs::TriggerId, std::vector<std::string>> _m_det_map;
};

#endif
