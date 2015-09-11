#ifndef RawTowerDigitizer_H__
#define RawTowerDigitizer_H__

#include <fun4all/SubsysReco.h>
#include <string>

#include <phool/PHTimeServer.h>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeom;

class RawTowerDigitizer : public SubsysReco
{

public:
  RawTowerDigitizer(const std::string& name = "RawTowerDigitizer");
  virtual
  ~RawTowerDigitizer()
  {
  }

  int
  InitRun(PHCompositeNode *topNode);
  int
  process_event(PHCompositeNode *topNode);
  int
  End(PHCompositeNode *topNode);
  void
  Detector(const std::string &d)
  {
    detector = d;
  }
  void
  EminCut(const double e)
  {
    emin = e;
  }

protected:
  void
  CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer* _sim_towers;
  RawTowerContainer* _raw_towers;
  RawTowerGeom *rawtowergeom;

  std::string detector;
  std::string SimTowerNodeName;
  std::string RawTowerNodeName;
  std::string TowerGeomNodeName;

  int _cell_binning;
  double emin;


  PHTimeServer::timer _timer;

};

#endif /* RawTowerDigitizer_H__ */
