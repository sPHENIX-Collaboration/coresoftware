#ifndef RAWTOWERCOMBINER_H__
#define RAWTOWERCOMBINER_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

#ifndef __CINT__
#include <boost/tuple/tuple.hpp>
#endif

#include <map>
#include <set>
#include <string>
#include <vector>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;

class RawTowerCombiner : public SubsysReco {

 public:
  RawTowerCombiner(const std::string& name="RawTowerCombiner");
  virtual ~RawTowerCombiner(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void AddInputDetector(const std::string &d, const int index, const double sf);
  void checkenergy(const int i = 1) {chkenergyconservation = i;}
  void OutputDetector(const std::string &d) {detector = d;}

 protected:
  bool CompareGeometries(RawTowerGeomContainer *geo1, RawTowerGeomContainer *geo2);

  std::string detector;
  std::string TowerNodeName;
  std::string TowerGeomNodeName;
  int chkenergyconservation;
  int iphibins;
  int ietabins;
  PHTimeServer::timer _timer;
#ifndef __CINT__
  std::map<std::string, boost::tuple<std::string,int, RawTowerContainer *, double> > dettuple;
#endif

};

#endif /* RAWTOWERCOMBINER_H__ */
