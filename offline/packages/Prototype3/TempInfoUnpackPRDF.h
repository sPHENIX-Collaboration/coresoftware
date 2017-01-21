#ifndef __TempInfoUnpackPRDFF__
#define __TempInfoUnpackPRDFF__

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <string>
#include <map>
#include <utility>

class Event;
class Packet;
class RawTowerContainer;
class RawTower;

class TempInfoUnpackPRDF : public SubsysReco
{
public:
  TempInfoUnpackPRDF();
  virtual ~TempInfoUnpackPRDF() {};

  int
  Init(PHCompositeNode *topNode);

  int
  InitRun(PHCompositeNode *topNode);

  int
  process_event(PHCompositeNode *topNode);

  int
  End(PHCompositeNode *topNode);

  void
  CreateNodeTree(PHCompositeNode *topNode);


protected:

  int addPacketInfo(Packet *p, PHCompositeNode *topNode, const time_t  etime, const int evtnr);

  RawTowerContainer* hcalin_temperature;
  RawTowerContainer* hcalout_temperature;
  RawTowerContainer* emcal_temperature;

};


#endif //**TempInfoUnpackPRDFF**//
