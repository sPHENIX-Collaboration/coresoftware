#ifndef CALORECO_RAWCLUSTERDEADHOTMASK_H
#define CALORECO_RAWCLUSTERDEADHOTMASK_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class RawTowerDeadMap;
class TowerInfoContainer;

class RawClusterDeadHotMask : public SubsysReco
{
 public:
  explicit RawClusterDeadHotMask(const std::string &name = "RawClusterDeadHotMask");

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  //! half width of the dead area around the center of a dead tower, if cluster center fall into this region it is rejected.
  void deadTowerMaskHalfWidth(double deadTowerMaskHalfWidth)
  {
    m_deadTowerMaskHalfWidth = deadTowerMaskHalfWidth;
  }

  void detector(const std::string &detector)
  {
    m_detector = detector;
  }

  void set_UseTowerInfo(const bool useMode)
  {
    m_UseTowerInfo = useMode;
  }

 private:
  void CreateNodeTree(PHCompositeNode *topNode);

  std::string m_detector;

  //! half width of the dead area around the center of a dead tower, if cluster center fall into this region it is rejected.
  double m_deadTowerMaskHalfWidth;

  RawClusterContainer *m_rawClusters;
  RawClusterContainer *m_towerinfoClusters;
  RawTowerDeadMap *m_deadMap;
  RawTowerContainer *m_calibTowers;
  TowerInfoContainer *m_calibTowerInfos;
  RawTowerGeomContainer *m_geometry;
  bool m_UseTowerInfo = true;
};

#endif
