#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packetConstants.h>
#include <Event/packet.h>
#include <Event/packet_hbd_fpgashort.h>
#include "RawTower_Prototype2.h"
#include <g4cemc/RawTowerContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>
#include <phool/getClass.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include "FEM.h"
#include <iostream>
#include <string>
#include "HCalUnpackPRDF.h"

using namespace std;

//____________________________________
HCalUnpackPRDF::HCalUnpackPRDF() : SubsysReco( "HCalUnpackPRDF" )
{
  verbosity = 0;
  _nevents = 0;
}


//____________________________________
int HCalUnpackPRDF::Init(PHCompositeNode *topNode)
{
 return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________
int HCalUnpackPRDF::InitRun(PHCompositeNode *topNode)
{
 CreateNodeTree( topNode );
 return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int HCalUnpackPRDF::process_event(PHCompositeNode *topNode)
{
 _nevents++;
 _event = findNode::getClass<Event>( topNode, "PRDF"); 
 if(_event==0)
 {
  cout << "HCalUnpackPRDF::Process_Event - Event not found" << endl;
  return -1;
 }

 if(verbosity)
 {
  cout << PHWHERE << "Process event entered" << std::endl;
 }

 if(verbosity) _event->identify();
 _packet = dynamic_cast<Packet_hbd_fpgashort*>(_event->getPacket(FEM::PACKET_ID));

 if(!_packet)
 {
  //They could be special events at the beginning or end of run
  if(_event->getEvtType()==1)
  {
   cout << "HCalUnpackPRDF::Process_Event - Packet not found" << endl;
   _event->identify();
  }
   return Fun4AllReturnCodes::ABORTEVENT;
 }

 _packet->setNumSamples( FEM::NSAMPLES );
 RawTower_Prototype2 *tower;

 //HCALIN
 for(int ibinz=0; ibinz<FEM::NCH_IHCAL_ROWS; ibinz++)
 {
  for(int ibinphi=0; ibinphi<FEM::NCH_IHCAL_COLUMNS; ibinphi++)
  {
   tower = dynamic_cast<RawTower_Prototype2*>(hcalin_towers->getTower(ibinz,ibinphi));
   if(!tower)
   {
    tower = new RawTower_Prototype2();
    tower->set_energy(0);
    hcalin_towers->AddTower(ibinz,ibinphi,tower);
   }
   int ich = GetHBDCh("HCALIN",ibinz,ibinphi);
   tower->set_HBD_channel_number(ich);
   for(int isamp=0; isamp<FEM::NSAMPLES; isamp++)
   {
    tower->set_signal_samples_hg(isamp,_packet->iValue(ich,isamp));
    tower->set_signal_samples_lg(isamp,_packet->iValue(ich+1,isamp));
   }
  }
 }


 //HCALOUT
 for(int ibinz=0; ibinz<FEM::NCH_OHCAL_ROWS; ibinz++)
 {
  for(int ibinphi=0; ibinphi<FEM::NCH_OHCAL_COLUMNS; ibinphi++)
  {
   tower = dynamic_cast<RawTower_Prototype2*>(hcalout_towers->getTower(ibinz,ibinphi));
   if(!tower)
   {
    tower = new RawTower_Prototype2();
    tower->set_energy(0);
    hcalout_towers->AddTower(ibinz,ibinphi,tower);
   }
   int ich = GetHBDCh("HCALOUT",ibinz,ibinphi);
   tower->set_HBD_channel_number(ich);
   for(int isamp=0; isamp<FEM::NSAMPLES; isamp++)
   {
    tower->set_signal_samples_hg(isamp,_packet->iValue(ich,isamp));
    tower->set_signal_samples_lg(isamp,_packet->iValue(ich+1,isamp));
   }
  }
 }

 //EMCAL
 for(int ibinz=0; ibinz<FEM::NCH_EMCAL_ROWS; ibinz++)
 {
  for(int ibinphi=0; ibinphi<FEM::NCH_EMCAL_COLUMNS; ibinphi++)
  {
   tower = dynamic_cast<RawTower_Prototype2*>(emcal_towers->getTower(ibinz,ibinphi));
   if(!tower)
   {
    tower = new RawTower_Prototype2();
    tower->set_energy(0);
    emcal_towers->AddTower(ibinz,ibinphi,tower);
   }
   int ich = GetHBDCh("EMCAL",ibinz,ibinphi);
   tower->set_HBD_channel_number(ich);
   for(int isamp=0; isamp<FEM::NSAMPLES; isamp++)
   {
    tower->set_signal_samples_hg(isamp,_packet->iValue(ich,isamp));
   }
  }
 }

 if(verbosity)
 {
   cout << "HCALIN Towers: " << endl;
   hcalin_towers->identify();
   RawTowerContainer::ConstRange begin_end = hcalin_towers->getTowers();
   RawTowerContainer::ConstIterator iter;
   for(iter=begin_end.first; iter!=begin_end.second;++iter)
   {
    RawTower_Prototype2 *tower = dynamic_cast<RawTower_Prototype2*>(iter->second);
    tower->identify();
    cout << "Signal Samples: [" << endl;
    for(int isamp=0; isamp<FEM::NSAMPLES;isamp++)
    {
     cout << tower->get_signal_samples_hg(isamp) <<", ";
    }
    cout << " ]" << endl;
   }
  }
 delete _packet;
 return Fun4AllReturnCodes::EVENT_OK;
}

//________________________________________
int HCalUnpackPRDF::GetHBDCh(string caloname,int zbin,int phibin)
{
 if(caloname=="HCALIN")
 {
  return 64+8*zbin+2*phibin;
 } else if(caloname=="HCALOUT")
 {
  return 112+8*zbin+2*phibin;
 } else if(caloname=="EMCAL")
 {
  return 8*zbin+phibin;
 }

 return -9999;
}

//_______________________________________
void HCalUnpackPRDF::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator nodeItr( topNode );
  //DST node
  dst_node = static_cast<PHCompositeNode*>( nodeItr.findFirst("PHCompositeNode", "DST" ));
  if(!dst_node)
  {
    cout << "PHComposite node created: DST" << endl;
    dst_node  = new PHCompositeNode( "DST" );
    topNode->addNode( dst_node );
  }

  //DATA nodes
  data_node =  static_cast<PHCompositeNode*>( nodeItr.findFirst("PHCompositeNode", "DATA_NODE" ));
  if(!data_node)
  {
    cout << "PHComposite node created: DATA_NODE" << endl;
    data_node  = new PHCompositeNode( "DATA_NODE" );
    dst_node->addNode( data_node );
  }


  //HCAL Towers
  hcalin_towers = new RawTowerContainer(RawTowerDefs::HCALIN);
  PHIODataNode<PHObject> *hcalin_towerNode = new PHIODataNode<PHObject>(hcalin_towers, "HCALIN_DATA_TOWERS", "PHObject" );
  data_node->addNode(hcalin_towerNode);

  hcalout_towers = new RawTowerContainer(RawTowerDefs::HCALOUT);
  PHIODataNode<PHObject> *hcalout_towerNode = new PHIODataNode<PHObject>(hcalout_towers, "HCALOUT_DATA_TOWERS", "PHObject" );
  data_node->addNode(hcalout_towerNode);

  //EMCAL towers
  emcal_towers = new RawTowerContainer(RawTowerDefs::CEMC);
  PHIODataNode<PHObject> *emcal_towerNode = new PHIODataNode<PHObject>(emcal_towers, "EMCAL_DATA_TOWERS", "PHObject" );
  data_node->addNode(emcal_towerNode);
}

//___________________________________
int HCalUnpackPRDF::End(PHCompositeNode *topNode)
{
 cout << "HCalUnpackPRDFF::End - Total Events: " << _nevents << endl;
 return Fun4AllReturnCodes::EVENT_OK;
}

