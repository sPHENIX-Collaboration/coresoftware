#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packetConstants.h>
#include <Event/packet.h>
#include "RawTower_Prototype2.h"
#include <g4cemc/RawTowerContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>
#include <phool/getClass.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include "PROTOTYPE2_FEM.h"
#include <iostream>
#include <string>
#include <cassert>
#include "RunInfoUnpackPRDF.h"

using namespace std;

//____________________________________
RunInfoUnpackPRDF::RunInfoUnpackPRDF(const string & detector) :
    SubsysReco("RunInfoUnpackPRDF")
{
}

//____________________________________
int
RunInfoUnpackPRDF::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________
int
RunInfoUnpackPRDF::InitRun(PHCompositeNode *topNode)
{
  CreateNodeTree(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int
RunInfoUnpackPRDF::process_event(PHCompositeNode *topNode)
{
  Event* event = findNode::getClass<Event>(topNode, "PRDF");
  if (event == NULL)
    {
      if (Verbosity() >= VERBOSITY_SOME)
        cout << "RunInfoUnpackPRDF::Process_Event - Event not found" << endl;
      return Fun4AllReturnCodes::DISCARDEVENT;
    }

  if (verbosity >= VERBOSITY_SOME)
    {

      cout << "RunInfoUnpackPRDF::process_event - ";
      event->identify();
    }

  if (event->getEvtType() != BEGRUNEVENT)
    return Fun4AllReturnCodes::EVENT_OK;

  map<int, Packet*> packet_list;
//
//  packet_list[packet_id] =
//      dynamic_cast<Packet_hbd_fpgashort*>(_event->getPacket(packet_id));
//  Packet_hbd_fpgashort * packet = packet_list[packet_id];
//
//  tower->set_signal_samples(isamp, packet->iValue(channel, isamp));

  for (typ_channel_map::const_iterator it = channel_map.begin();
      it != channel_map.end(); ++it)
    {
      const string & name = it->first;
      const channel_info & info = it->second;

      if (packet_list.find(info.packet_id) == packet_list.end())
        {
          packet_list[info.packet_id] = event->getPacket(info.packet_id);
        }

      Packet * packet = packet_list[info.packet_id];

      if (!packet)
        {
//          if (Verbosity() >= VERBOSITY_SOME)
          cout << "RunInfoUnpackPRDF::process_event - failed to locate packet "
              << info.packet_id << " from ";
          event->identify();

          continue;
        }

      const int ivalue = packet->iValue(info.offset);

      const double dvalue = ivalue * info.calibration_const;

      if (verbosity >= VERBOSITY_SOME)
        {
          cout << "RunInfoUnpackPRDF::process_event - " << name << " = "
              << dvalue << ", raw = " << ivalue << " @ packet "
              << info.packet_id << ", offset " << info.offset << endl;
        }
    }

  for (map<int, Packet*>::iterator it = packet_list.begin();
      it != packet_list.end(); ++it)
    {
      if (it->second)
        delete it->second;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________
void
RunInfoUnpackPRDF::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator nodeItr(topNode);
  //DST node
  PHCompositeNode * run_node = static_cast<PHCompositeNode*>(nodeItr.findFirst(
      "PHCompositeNode", "RUN"));
  if (!run_node)
    {
      cout << "PHComposite node created: RUN" << endl;
      run_node = new PHCompositeNode("RUN");
      topNode->addNode(run_node);
    }

}

//___________________________________
int
RunInfoUnpackPRDF::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void
RunInfoUnpackPRDF::add_channel(const std::string & name, //! name of the channel
    const int packet_id, //! packet id
    const unsigned int offset, //! offset in packet data
    const double calibration_const //! conversion constant from integer to meaningful value
    )
{
  channel_map.insert(
      make_pair(name, channel_info(packet_id, offset, calibration_const)));
}
