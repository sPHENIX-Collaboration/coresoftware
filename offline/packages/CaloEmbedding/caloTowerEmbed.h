// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTOWEREMBED_H
#define CALOTOWEREMBED_H

#include <calobase/TowerInfoContainer.h>  // for TowerInfoContainer, TowerIn...

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

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

  void set_useRetower(bool a) { m_useRetower = a; };

  enum DetectorSystem
  {
    CEMC = 0,
    HCALIN = 1,
    HCALOUT = 2,
    EPD = 3
  };

 private:
  TowerInfoContainer *_data_towers[3] {nullptr, nullptr, nullptr};
  TowerInfoContainer *_sim_towers[3] {nullptr, nullptr, nullptr};

  RawTowerGeomContainer *tower_geom {nullptr};
  RawTowerGeomContainer *tower_geomIH {nullptr};
  RawTowerGeomContainer *tower_geomOH {nullptr};

  bool m_useRetower {false};

  //  CDBInterface *cdb = nullptr;
  // CDBTTree *cdbttree = nullptr;
  int m_runNumber {0};
  int m_eventNumber {-1};
};

#endif  // CALOTOWEREMBED_H
