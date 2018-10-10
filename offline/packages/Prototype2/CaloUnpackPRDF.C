#include "RawTower_Prototype2.h"
#include "PROTOTYPE2_FEM.h"
#include "CaloUnpackPRDF.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packetConstants.h>
#include <Event/packet.h>
#include <Event/packet_hbd_fpgashort.h>
#include <calobase/RawTowerContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>
#include <phool/getClass.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <iostream>
#include <string>
#include <cassert>

using namespace std;

//____________________________________
CaloUnpackPRDF::CaloUnpackPRDF() : SubsysReco( "CaloUnpackPRDF" ),
    /*Event**/ _event(NULL),
    /*Packet_hbd_fpgashort**/ _packet(NULL),
    /*int*/ _nevents(0),
    /*PHCompositeNode **/ dst_node(NULL),
    /*PHCompositeNode **/ data_node(NULL),
    /*RawTowerContainer**/ hcalin_towers_lg(NULL),
    /*RawTowerContainer**/ hcalout_towers_lg(NULL),
    /*RawTowerContainer**/ hcalin_towers_hg(NULL),
    /*RawTowerContainer**/ hcalout_towers_hg(NULL),
    /*RawTowerContainer**/ emcal_towers(NULL)
{
}


//____________________________________
int CaloUnpackPRDF::Init(PHCompositeNode *topNode)
{
 return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________
int CaloUnpackPRDF::InitRun(PHCompositeNode *topNode)
{
 CreateNodeTree( topNode );
 return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int CaloUnpackPRDF::process_event(PHCompositeNode *topNode)
{
 _nevents++;
 _event = findNode::getClass<Event>( topNode, "PRDF"); 
 if(_event==0)
 {
  cout << "CaloUnpackPRDF::Process_Event - Event not found" << endl;
  return -1;
 }

 if(Verbosity())
 {
  cout << PHWHERE << "Process event entered" << std::endl;
 }

 if(Verbosity()) _event->identify();
 _packet = dynamic_cast<Packet_hbd_fpgashort*>(_event->getPacket(PROTOTYPE2_FEM::PACKET_ID));

 if(!_packet)
 {
  //They could be special events at the beginning or end of run
  if(_event->getEvtType()==1)
  {
   cout << "CaloUnpackPRDF::Process_Event - Packet not found" << endl;
   _event->identify();
  }
   return Fun4AllReturnCodes::ABORTEVENT;
 }

 _packet->setNumSamples( PROTOTYPE2_FEM::NSAMPLES );
 RawTower_Prototype2 *tower_lg = NULL;
 RawTower_Prototype2 *tower_hg = NULL;

 //HCALIN
 assert(hcalin_towers_lg);
 assert(hcalin_towers_hg);
 for(int ibinz=0; ibinz<PROTOTYPE2_FEM::NCH_IHCAL_ROWS; ibinz++)
 {
  for(int ibinphi=0; ibinphi<PROTOTYPE2_FEM::NCH_IHCAL_COLUMNS; ibinphi++)
  {
   tower_lg = dynamic_cast<RawTower_Prototype2*>(hcalin_towers_lg->getTower(ibinz,ibinphi));
   if(!tower_lg)
   {
       tower_lg = new RawTower_Prototype2();
       tower_lg->set_energy(NAN);
       hcalin_towers_lg->AddTower(ibinz,ibinphi,tower_lg);
   }
   tower_hg = dynamic_cast<RawTower_Prototype2*>(hcalin_towers_hg->getTower(ibinz,ibinphi));
   if(!tower_hg)
   {
       tower_hg = new RawTower_Prototype2();
       tower_hg->set_energy(NAN);
       hcalin_towers_hg->AddTower(ibinz,ibinphi,tower_hg);
   }

   int ich = PROTOTYPE2_FEM::GetHBDCh("HCALIN",ibinz,ibinphi);
   tower_lg->set_HBD_channel_number(ich);
   tower_hg->set_HBD_channel_number(ich);
   for(int isamp=0; isamp<PROTOTYPE2_FEM::NSAMPLES; isamp++)
   {
       tower_hg->set_signal_samples(isamp,_packet->iValue(ich,isamp));
       tower_lg->set_signal_samples(isamp,_packet->iValue(ich+1,isamp));
   }
  }
 }

 //HCALOUT
 assert(hcalout_towers_lg);
 assert(hcalout_towers_hg);
 for(int ibinz=0; ibinz<PROTOTYPE2_FEM::NCH_OHCAL_ROWS; ibinz++)
 {
  for(int ibinphi=0; ibinphi<PROTOTYPE2_FEM::NCH_OHCAL_COLUMNS; ibinphi++)
  {
      tower_lg = dynamic_cast<RawTower_Prototype2*>(hcalout_towers_lg->getTower(ibinz,ibinphi));
      if(!tower_lg)
      {
          tower_lg = new RawTower_Prototype2();
          tower_lg->set_energy(NAN);
          hcalout_towers_lg->AddTower(ibinz,ibinphi,tower_lg);
      }
      tower_hg = dynamic_cast<RawTower_Prototype2*>(hcalout_towers_hg->getTower(ibinz,ibinphi));
      if(!tower_hg)
      {
          tower_hg = new RawTower_Prototype2();
          tower_hg->set_energy(NAN);
          hcalout_towers_hg->AddTower(ibinz,ibinphi,tower_hg);
      }
   int ich = PROTOTYPE2_FEM::GetHBDCh("HCALOUT",ibinz,ibinphi);
   tower_lg->set_HBD_channel_number(ich);
   tower_hg->set_HBD_channel_number(ich);
   for(int isamp=0; isamp<PROTOTYPE2_FEM::NSAMPLES; isamp++)
   {
       tower_hg->set_signal_samples(isamp,_packet->iValue(ich,isamp));
       tower_lg->set_signal_samples(isamp,_packet->iValue(ich+1,isamp));
   }
  }
 }

 //EMCAL
 RawTower_Prototype2 *tower = NULL;
 assert(emcal_towers);
 for(int ibinz=0; ibinz<PROTOTYPE2_FEM::NCH_EMCAL_ROWS; ibinz++)
 {
  for(int ibinphi=0; ibinphi<PROTOTYPE2_FEM::NCH_EMCAL_COLUMNS; ibinphi++)
  {
   tower = dynamic_cast<RawTower_Prototype2*>(emcal_towers->getTower(ibinz,ibinphi));
   if(!tower)
   {
    tower = new RawTower_Prototype2();
    tower->set_energy(NAN);
    emcal_towers->AddTower(ibinz,ibinphi,tower);
   }
   int ich = PROTOTYPE2_FEM::GetHBDCh("EMCAL",ibinz,ibinphi);
   tower->set_HBD_channel_number(ich);
   for(int isamp=0; isamp<PROTOTYPE2_FEM::NSAMPLES; isamp++)
   {
    tower->set_signal_samples(isamp,_packet->iValue(ich,isamp));
   }
  }
 }

 if(Verbosity())
 {
   cout << "HCALIN Towers: " << endl;
   hcalin_towers_hg->identify();
   RawTowerContainer::ConstRange begin_end = hcalin_towers_lg->getTowers();
   RawTowerContainer::ConstIterator iter;
   for(iter=begin_end.first; iter!=begin_end.second;++iter)
   {
    RawTower_Prototype2 *tower = dynamic_cast<RawTower_Prototype2*>(iter->second);
    tower->identify();
    cout << "Signal Samples: [" << endl;
    for(int isamp=0; isamp<PROTOTYPE2_FEM::NSAMPLES;isamp++)
    {
     cout << tower->get_signal_samples(isamp) <<", ";
    }
    cout << " ]" << endl;
   }
  }
 delete _packet;
 return Fun4AllReturnCodes::EVENT_OK;
}


//_______________________________________
void CaloUnpackPRDF::CreateNodeTree(PHCompositeNode *topNode)
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
  data_node =  static_cast<PHCompositeNode*>( nodeItr.findFirst("PHCompositeNode", "RAW_DATA" ));
  if(!data_node)
  {
      if (Verbosity())
    cout << "PHComposite node created: RAW_DATA" << endl;
    data_node  = new PHCompositeNode( "RAW_DATA" );
    dst_node->addNode( data_node );
  }

  PHIODataNode<PHObject> *tower_node = NULL;

  //HCAL Towers
  hcalin_towers_lg = new RawTowerContainer(RawTowerDefs::HCALIN);
  tower_node = new PHIODataNode<PHObject>(hcalin_towers_lg, "TOWER_RAW_LG_HCALIN", "PHObject" );
  data_node->addNode(tower_node);
  hcalin_towers_hg = new RawTowerContainer(RawTowerDefs::HCALIN);
  tower_node = new PHIODataNode<PHObject>(hcalin_towers_hg, "TOWER_RAW_HG_HCALIN", "PHObject" );
  data_node->addNode(tower_node);

  hcalout_towers_lg = new RawTowerContainer(RawTowerDefs::HCALOUT);
  tower_node = new PHIODataNode<PHObject>(hcalout_towers_lg, "TOWER_RAW_LG_HCALOUT", "PHObject" );
  data_node->addNode(tower_node);
  hcalout_towers_hg = new RawTowerContainer(RawTowerDefs::HCALOUT);
  tower_node = new PHIODataNode<PHObject>(hcalout_towers_hg, "TOWER_RAW_HG_HCALOUT", "PHObject" );
  data_node->addNode(tower_node);

  //EMCAL towers
  emcal_towers = new RawTowerContainer(RawTowerDefs::CEMC);
  PHIODataNode<PHObject> *emcal_towerNode = new PHIODataNode<PHObject>(emcal_towers, "TOWER_RAW_CEMC", "PHObject" );
  data_node->addNode(emcal_towerNode);
}

//___________________________________
int CaloUnpackPRDF::End(PHCompositeNode *topNode)
{
 return Fun4AllReturnCodes::EVENT_OK;
}

