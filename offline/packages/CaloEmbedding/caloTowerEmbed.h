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

 private:
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

  int m_eventNumber{-1};
  int m_nsamples{31};
  int m_datasamples{12};
  float m_pedestal_scale{1.};
};

#endif  // CALOTOWEREMBED_H
