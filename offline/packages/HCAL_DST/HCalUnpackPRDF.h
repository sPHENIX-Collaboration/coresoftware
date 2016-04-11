#ifndef __HCALUNPACKPRDFF__
#define __HCALUNPACKPRDFF__

//* Unpacks raw HCAL PRDF files *//
//Abhisek Sen

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>

class Event;
class Packet;
class Packet_hbd_fpgashort;
class RawTowerContainer;
class RawTower;


class HCalUnpackPRDF : public SubsysReco
{
 public:
  HCalUnpackPRDF();

  int Init(PHCompositeNode *topNode);

  int InitRun(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  int End(PHCompositeNode *topNode);
  
  void CreateNodeTree(PHCompositeNode *topNode);

  int GetHBDCh(std::string,int,int);

 private:
  Event* _event;
  Packet_hbd_fpgashort* _packet;
  int _nevents; 
  // HCAL node
  PHCompositeNode * dst_node;
  PHCompositeNode * data_node;
  //Towers
  RawTowerContainer* hcalin_towers;
  RawTowerContainer* hcalout_towers;
  RawTowerContainer* emcal_towers;
};


#endif //**HCALUNPACKPRDFF**//
