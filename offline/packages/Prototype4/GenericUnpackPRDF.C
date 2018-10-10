#include "GenericUnpackPRDF.h"
#include "PROTOTYPE4_FEM.h"
#include "RawTower_Prototype4.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>
#include <Event/packetConstants.h>
#include <calobase/RawTowerContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <cassert>
#include <iostream>
#include <string>

using namespace std;

//____________________________________
GenericUnpackPRDF::GenericUnpackPRDF(const string &detector)
  : SubsysReco("GenericUnpackPRDF_" + detector)
  ,  //
  _detector(detector)
  ,
  /*Event**/ _event(NULL)
  ,
  /*PHCompositeNode **/ _towers(NULL)
{
}

//____________________________________
int GenericUnpackPRDF::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________
int GenericUnpackPRDF::InitRun(PHCompositeNode *topNode)
{
  CreateNodeTree(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int GenericUnpackPRDF::process_event(PHCompositeNode *topNode)
{
  _event = findNode::getClass<Event>(topNode, "PRDF");
  if (_event == NULL)
  {
    if (Verbosity() >= VERBOSITY_SOME)
      cout << "GenericUnpackPRDF::Process_Event - Event not found" << endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  if (Verbosity() >= VERBOSITY_SOME)
    _event->identify();

  // search for data event
  if (_event->getEvtType() != DATAEVENT)
    return Fun4AllReturnCodes::EVENT_OK;

  map<int, Packet *> packet_list;

  for (channel_map::const_iterator it = _channel_map.begin();
       it != _channel_map.end(); ++it)
  {
    const int packet_id = it->first.first;
    const int channel = it->first.second;
    const int tower_id = it->second;

    if (packet_list.find(packet_id) == packet_list.end())
    {
      packet_list[packet_id] =_event->getPacket(packet_id);

      if (packet_list[packet_id] and Verbosity() >= VERBOSITY_MORE)
      {
        cout << "GenericUnpackPRDF::process_event - open packet " << packet_id << ":" << endl;
//        packet_list[packet_id]->identify(cout);
        packet_list[packet_id]->dump(cout);
      }
    }

    Packet *packet = packet_list[packet_id];

    if (!packet)
    {
      //          if (Verbosity() >= VERBOSITY_SOME)
      cout << "GenericUnpackPRDF::process_event - failed to locate packet "
           << packet_id << endl;
      continue;
    }
    assert(packet);

    //    packet->setNumSamples(PROTOTYPE4_FEM::NSAMPLES);

    RawTower_Prototype4 *tower =
        dynamic_cast<RawTower_Prototype4 *>(_towers->getTower(tower_id));
    if (!tower)
    {
      tower = new RawTower_Prototype4();
      tower->set_energy(NAN);
      _towers->AddTower(tower_id, tower);
    }
    tower->set_HBD_channel_number(channel);
    for (int isamp = 0; isamp < PROTOTYPE4_FEM::NSAMPLES; isamp++)
    {
      tower->set_signal_samples(isamp, (packet->iValue( isamp, channel)) & PROTOTYPE4_FEM::ADC_DATA_MASK);
    }
  }

  for (map<int, Packet *>::iterator it = packet_list.begin();
       it != packet_list.end(); ++it)
  {
    if (it->second)
      delete it->second;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________
void GenericUnpackPRDF::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator nodeItr(topNode);
  //DST node
  PHCompositeNode *dst_node = static_cast<PHCompositeNode *>(nodeItr.findFirst(
      "PHCompositeNode", "DST"));
  if (!dst_node)
  {
    cout << "PHComposite node created: DST" << endl;
    dst_node = new PHCompositeNode("DST");
    topNode->addNode(dst_node);
  }

  //DATA nodes
  PHCompositeNode *data_node = static_cast<PHCompositeNode *>(nodeItr.findFirst(
      "PHCompositeNode", "RAW_DATA"));
  if (!data_node)
  {
    if (Verbosity())
      cout << "PHComposite node created: RAW_DATA" << endl;
    data_node = new PHCompositeNode("RAW_DATA");
    dst_node->addNode(data_node);
  }

  //output as towers
  _towers = new RawTowerContainer(RawTowerDefs::NONE);
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_towers,
                                                                 "TOWER_RAW_" + _detector, "PHObject");
  data_node->addNode(towerNode);
}

//___________________________________
int GenericUnpackPRDF::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void GenericUnpackPRDF::add_channel(const int packet_id,  //! packet id
                                    const int channel,    //! channel in packet
                                    const int tower_id    //! output tower id
                                    )
{
  channel_typ hbd_channel(packet_id, channel);

  if (_channel_map.find(hbd_channel) != _channel_map.end())
  {
    cout << "GenericUnpackPRDF::add_channel - packet " << packet_id
         << ", channel " << channel << " is already registered as tower "
         << _channel_map.find(hbd_channel)->second << endl;
    exit(12);
  }
  _channel_map[hbd_channel] = tower_id;
}
