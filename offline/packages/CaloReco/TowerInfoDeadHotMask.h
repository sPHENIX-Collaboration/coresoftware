#ifndef CALORECO_TOWERINFODEADHOTMASK_H
#define CALORECO_TOWERINFODEADHOTMASK_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class RawTowerDeadMap;
class TowerInfoContainer;

class TowerInfoDeadHotMask : public SubsysReco
{
 public:
  explicit TowerInfoDeadHotMask(const std::string &name = "TowerInfoDeadHotMask");

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void detector(const std::string &detector)
  {
    m_detector = detector;
  }

 private:
  void CreateNodeTree(PHCompositeNode *topNode);

  std::string m_detector;

  RawTowerDeadMap *m_deadMap;
  TowerInfoContainer *m_calibTowerInfos;
  RawTowerGeomContainer *m_geometry;
};

#endif
