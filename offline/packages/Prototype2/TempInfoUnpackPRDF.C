#include "RawTower_Temperature.h"
#include "PROTOTYPE2_FEM.h"
#include "TempInfoUnpackPRDF.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packetConstants.h>
#include <Event/packet.h>
#include <calobase/RawTowerContainer.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>
#include <phool/getClass.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <iostream>
#include <string>
#include <cassert>

using namespace std;

//____________________________________
TempInfoUnpackPRDF::TempInfoUnpackPRDF() :
  SubsysReco("TempInfoUnpackPRDF"), hcalin_temperature(NULL), hcalout_temperature(NULL), emcal_temperature(NULL)
{
}

//____________________________________
int
TempInfoUnpackPRDF::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________
int
TempInfoUnpackPRDF::InitRun(PHCompositeNode *topNode)
{
  CreateNodeTree(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int
TempInfoUnpackPRDF::process_event(PHCompositeNode *topNode)
{
  Event* event = findNode::getClass<Event>(topNode, "PRDF");
  if (!event)
    {
      if (Verbosity() >= VERBOSITY_SOME)
        cout << "TempInfoUnpackPRDF::Process_Event - Event not found" << endl;
      return Fun4AllReturnCodes::DISCARDEVENT;
    }


  Packet *p_hcalin;
  Packet *p_hcalout;
  Packet *p_emcal;

  // if (Verbosity() >= VERBOSITY_SOME)
  //   {
  //     cout << "TempInfoUnpackPRDF::process_event - ";
  //     event->identify();
  //   }


  if (event->getEvtType() == BEGRUNEVENT)
    {
      p_hcalin  = event->getPacket(974);
      p_hcalout = event->getPacket(975);
      p_emcal   = event->getPacket(982);
    }
  else
    {
      p_hcalin  = event->getPacket(1074);
      p_hcalout = event->getPacket(1075);
      p_emcal   = event->getPacket(1082);
    }

  time_t etime= event->getTime();
  int evtnr = event->getEvtSequence();


  if (Verbosity() >= VERBOSITY_SOME && ( p_hcalin || p_hcalout || p_emcal) )
    {
      cout << "TempInfoUnpackPRDF::found temperature packet in Event - ";
      event->identify();
    }



  if ( p_hcalin)
    {
      addPacketInfo (p_hcalin, topNode, etime, evtnr);
      if (Verbosity() > VERBOSITY_SOME) p_hcalin->dump();
      delete p_hcalin;
    }

  if ( p_hcalout)
    {
      addPacketInfo (p_hcalout, topNode, etime, evtnr);
      if (Verbosity() > VERBOSITY_SOME) p_hcalout->dump();
      delete p_hcalout;
    }

  if ( p_emcal)
    {
      addPacketInfo (p_emcal, topNode, etime, evtnr);
      if (Verbosity() > VERBOSITY_SOME) p_emcal->dump();
      delete p_emcal;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________
int  TempInfoUnpackPRDF::addPacketInfo(Packet *p, PHCompositeNode *topNode, const time_t  etime, const int evtnr)
{

  int packetid = p->getIdentifier();

  RawTower_Temperature *tower;


  if ( packetid == 974 || packetid == 1074 ) // Inner Hcal
    {
      for(int ibinz=0; ibinz<PROTOTYPE2_FEM::NCH_IHCAL_ROWS; ibinz++)
	{
	  for(int ibinphi=0; ibinphi<PROTOTYPE2_FEM::NCH_IHCAL_COLUMNS; ibinphi++)
	    {
	      tower = dynamic_cast<RawTower_Temperature*>(hcalin_temperature->getTower(ibinz,ibinphi));
	      if(!tower)
		{
		  tower = new RawTower_Temperature();
		  hcalin_temperature->AddTower(ibinz,ibinphi,tower);
		}
	      tower->add_entry( evtnr, etime, p->iValue(ibinz*PROTOTYPE2_FEM::NCH_IHCAL_COLUMNS + ibinphi)/1000. ); 
	    }
	}
    }

  else  if ( packetid == 975 || packetid == 1075 ) // outer Hcal
    {
      for(int ibinz=0; ibinz<PROTOTYPE2_FEM::NCH_OHCAL_ROWS; ibinz++)
	{
	  for(int ibinphi=0; ibinphi<PROTOTYPE2_FEM::NCH_OHCAL_COLUMNS; ibinphi++)
	    {
	      tower = dynamic_cast<RawTower_Temperature*>(hcalout_temperature->getTower(ibinz,ibinphi));
	      if(!tower)
		{
		  tower = new RawTower_Temperature();
		  hcalout_temperature->AddTower(ibinz,ibinphi,tower);
		}
	      tower->add_entry( evtnr, etime, p->iValue(ibinz*PROTOTYPE2_FEM::NCH_OHCAL_COLUMNS + ibinphi)/1000. ); 
	    }
	}
    }

  else  if ( packetid == 982 || packetid == 1082 ) // emcal
    {
      for(int ibinz=0; ibinz<PROTOTYPE2_FEM::NCH_EMCAL_ROWS; ibinz++)
	{
	  for(int ibinphi=0; ibinphi<PROTOTYPE2_FEM::NCH_EMCAL_COLUMNS; ibinphi++)
	    {
	      tower = dynamic_cast<RawTower_Temperature*>(emcal_temperature->getTower(ibinz,ibinphi));
	      if(!tower)
		{
		  tower = new RawTower_Temperature();
		  emcal_temperature->AddTower(ibinz,ibinphi,tower);
		}
	      // this takes care of the newly found "reverse" mapping. (0,0) is module 7, (0,7) is module 0, and 
	      // the 63 - (...) takes care of the reversed vector.  
	      tower->add_entry( evtnr, etime, p->iValue( 63- (ibinz*PROTOTYPE2_FEM::NCH_EMCAL_COLUMNS + (7-ibinphi) ) /1000. ) ); 
	    }
	}
    }
  return 0;
}






//_______________________________________
void TempInfoUnpackPRDF::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator nodeItr(topNode);
  //DST node
  PHCompositeNode * run_node = static_cast<PHCompositeNode*>(nodeItr.findFirst( "PHCompositeNode", "RUN"));
  if (!run_node)
    {
      run_node = new PHCompositeNode("RUN");
      topNode->addNode(run_node);
      cout << "PHComposite node created: RUN" << endl;
    }

  PHIODataNode<PHObject> *tower_node = NULL;


  //HCAL Towers
  hcalin_temperature = new RawTowerContainer(RawTowerDefs::HCALIN);
  tower_node = new PHIODataNode<PHObject>(hcalin_temperature, "TOWER_TEMPERATURE_HCALIN", "PHObject" );
  run_node->addNode(tower_node);

  hcalout_temperature = new RawTowerContainer(RawTowerDefs::HCALOUT);
  tower_node = new PHIODataNode<PHObject>(hcalout_temperature, "TOWER_TEMPERATURE_HCALOUT", "PHObject" );
  run_node->addNode(tower_node);

  emcal_temperature = new RawTowerContainer(RawTowerDefs::CEMC);
  tower_node = new PHIODataNode<PHObject>(emcal_temperature, "TOWER_TEMPERATURE_EMCAL", "PHObject" );
  run_node->addNode(tower_node);

}

//___________________________________
int TempInfoUnpackPRDF::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

