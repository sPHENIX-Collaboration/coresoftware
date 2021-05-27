#ifndef CALORECO_RAWTOWERDEADTOWERINTERP_H
#define CALORECO_RAWTOWERDEADTOWERINTERP_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;
class RawTowerDeadMap;

//! RawTowerDeadTowerInterp recovers the energy in the known dead towers with interpolation between alive towers near-by
class RawTowerDeadTowerInterp : public SubsysReco
{
 public:
  RawTowerDeadTowerInterp(const std::string &name = "RawTowerDeadTowerInterp");
  ~RawTowerDeadTowerInterp() override
  {
  }

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void
  detector(const std::string &d)
  {
    m_detector = d;
  }

 protected:
  void
  CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer *m_calibTowers;
  RawTowerGeomContainer *m_geometry;
  RawTowerDeadMap *m_deadTowerMap;

  std::string m_detector;

  std::string _calib_tower_node_prefix;
};

#endif
