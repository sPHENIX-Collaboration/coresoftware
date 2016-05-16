#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packetConstants.h>
#include <Event/packet.h>
#include "RawTower_Prototype2.h"
#include <pdbcalbase/PdbParameterMap.h>
#include <g4detectors/PHG4Parameters.h>
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
RunInfoUnpackPRDF::RunInfoUnpackPRDF() :
    SubsysReco("RunInfoUnpackPRDF"), runinfo_node_name("RUN_INFO")
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

  if (event->getEvtType() != BEGRUNEVENT)
    return Fun4AllReturnCodes::EVENT_OK;

  if (verbosity >= VERBOSITY_SOME)
    {

      cout << "RunInfoUnpackPRDF::process_event - ";
      event->identify();
    }

  map<int, Packet*> packet_list;

  PHG4Parameters Params("RunInfo");

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

          Params.set_double_param(name, NAN);
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

      Params.set_double_param(name, dvalue);
    }

  for (map<int, Packet*>::iterator it = packet_list.begin();
      it != packet_list.end(); ++it)
    {
      if (it->second)
        delete it->second;
    }

  Params.SaveToNodeTree(topNode, runinfo_node_name);

  if (verbosity >= VERBOSITY_SOME) Params.print();

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

  PdbParameterMap *nodeparams = findNode::getClass<PdbParameterMap>(run_node,
      runinfo_node_name);
  if (not nodeparams)
    {
      run_node->addNode(
          new PHIODataNode<PdbParameterMap>(new PdbParameterMap(),
              runinfo_node_name));
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
