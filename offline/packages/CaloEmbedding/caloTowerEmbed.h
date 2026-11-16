// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTOWEREMBED_H
#define CALOTOWEREMBED_H

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>  // for TowerInfoContainer, TowerIn...
#include <caloreco/CaloTowerDefs.h>

#include <fun4all/SubsysReco.h>

#include <TFile.h>
#include <TTree.h>

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

class PHCompositeNode;
class RawTowerGeomContainer;

class caloTowerEmbed : public SubsysReco
{
 public:
  caloTowerEmbed(const std::string &name = "caloTowerEmbed");

  ~caloTowerEmbed() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void CreateNodeTree(PHCompositeNode *topNode);

  void set_detector_type(CaloTowerDefs::DetectorSystem dettype)
  {
    m_dettype = dettype;
    return;
  }

  void set_inputNodePrefix(const std::string &name)
  {
    m_inputNodePrefix = name;
    return;
  }

  void set_useRetower(bool a)
  {
    m_useRetower = a;
    return;
  }

  void set_removeBadTowers(bool a)
  {
    m_removeBadTowers = a;
    return;
  }
  void set_nsamples(int _nsamples)
  {
    m_nsamples = _nsamples;
    return;
  }
  void set_embedwaveform(bool embed = true)
  {
    m_embedwaveform = embed;
    return;
  }
  // waveform mode: the software ZS threshold of the CaloTowerBuilder that fits the embedded
  // waveforms. Data HW ZS towers whose sim reaches it are built from the pedestal waveform plus
  // the mean data signal. 
  void set_softwarezerosuppression(float adc)
  {
    m_szs_threshold = adc;
    return;
  }
  // number of |eta| bands the mean data signal of the HW ZS towers is estimated in, per event
  void set_nsignalbands(int n)
  {
    m_nsignalbands = n;
    return;
  }

 private:
  int signal_band(unsigned int channel) const;
  void EstimateDataSignal();

  TowerInfoContainer *_data_towers{nullptr};
  TowerInfoContainer *_sim_towers{nullptr};
  TowerInfoContainer *m_PedestalContainer{nullptr};

  RawTowerGeomContainer *tower_geom{nullptr};

  bool m_useRetower{false};
  bool m_removeBadTowers{false};
  bool m_embedwaveform{false};

  CaloTowerDefs::DetectorSystem m_dettype{CaloTowerDefs::DETECTOR_INVALID};

  std::string m_detector;
  std::string m_inputNodePrefix{"TOWERINFO_CALIB_"};
  std::string m_waveformNodePrefix{"WAVEFORM_"};
  std::string m_templateName;

  // reco pulse template sampled at the data samples, peak (sample m_postsample) = 1
  std::vector<float> m_template;

  int m_eventNumber{-1};
  int m_nsamples{31};
  int m_datasamples{12};
  int m_postsample{6};  // the HW zero suppressed readout keeps samples 0 (pre) and 6 (post)
  float m_pedestal_scale{1.};
  float m_szs_threshold{-1};  // < 0: the detector default, set in InitRun
  int m_nsignalbands{3};
  int m_neta{96};
  std::vector<float> m_bandsignal;  // mean data signal of the HW ZS towers per |eta| band, this event
  long m_npedestalonly{0};
};

#endif  // CALOTOWEREMBED_H
