#ifndef __RawClusterDeadAreaMask_H__
#define __RawClusterDeadAreaMask_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <phparameter/PHParameters.h>
#include <string>

class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class RawTowerDeadMap;

class RawClusterDeadAreaMask : public SubsysReco
{
 public:
  explicit RawClusterDeadAreaMask(const std::string &name = "RawClusterDeadAreaMask");

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  //! half width of the dead area around the center of a dead tower, if cluster center fall into this region it is rejected.
  void deadTowerMaskHalfWidth(double deadTowerMaskHalfWidth)
  {
    m_deadTowerMaskHalfWidth = deadTowerMaskHalfWidth;
  }

  void detector(const std::string &detector)
  {
    m_detector = detector;
  }

 private:
  void CreateNodeTree(PHCompositeNode *topNode);

  std::string m_detector;

  //! half width of the dead area around the center of a dead tower, if cluster center fall into this region it is rejected.
  double m_deadTowerMaskHalfWidth;

  RawClusterContainer *m_rawClusters;
  RawTowerDeadMap *m_deadMap;
  RawTowerContainer *m_calibTowers;
  RawTowerGeomContainer *m_geometry;
};

#endif  // __RawClusterDeadAreaMask_H__
