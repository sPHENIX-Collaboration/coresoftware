#include "EventCombiner.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/oncsEvent.h>

#include <TSystem.h>
//____________________________________________________________________________..
EventCombiner::EventCombiner(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int EventCombiner::Init(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfOutputNodeName));
  if (!PrdfNode)
  {
    PHDataNode<Event> *newNode = new PHDataNode<Event>(m_Event, m_PrdfOutputNodeName, "Event");
    std::cout << "Creating new prdfnode 0x" << std::hex << newNode << std::dec << std::endl;
    topNode->addNode(newNode);
  }

  std::cout << "EventCombiner::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventCombiner::process_event(PHCompositeNode *topNode)
{
  std::vector<Event *> subeventeventvec;
  unsigned int total_length = 0;
  for (auto &nam : m_PrdfInputNodeNameSet)
  {
    Event *evt = findNode::getClass<Event>(topNode, nam);
    subeventeventvec.push_back(evt);
    total_length += evt->getEvtLength();
  }
  // safety belts
  int eventno = subeventeventvec[0]->getEvtSequence();
  for (auto &e : subeventeventvec)
  {
    if (e->getEvtSequence() != eventno)
    {
      std::cout << "Event number mismatch, first subevt: " << eventno
                << " current subevt: " << e->getEvtSequence() << std::endl;
    }
  }

  // lets copy them all together
  int nwout;
  int current = 0;
  m_OutArray = new int[total_length];
  // the first Copy is without "DATA", "DATA" skips the first few words which contain event number, length, etc
  subeventeventvec[0]->Copy(m_OutArray, total_length, &nwout);
  current = nwout;
  for (unsigned int icnt = 1; icnt < subeventeventvec.size(); icnt++)
  {
    subeventeventvec[icnt]->Copy(&m_OutArray[current], total_length - current, &nwout, "DATA");
    current += nwout;
    m_OutArray[0] += nwout;
  }

  PHNodeIterator iter(topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfOutputNodeName));

  m_Event = new oncsEvent(m_OutArray);
  PrdfNode->setData(m_Event);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
// the events need to be cleaned up outside of the framework
// in principal this is a prdf input manager which has to clear the data it puts on the tree
int EventCombiner::ResetEvent(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfOutputNodeName));
  PrdfNode->setData(nullptr);  // set pointer in Node to nullptr before deleting it
  delete m_Event;
  m_Event = nullptr;
  delete[] m_OutArray;
  return Fun4AllReturnCodes::EVENT_OK;
}

void EventCombiner::AddPrdfInputNodeName(const std::string &name)
{
  auto result = m_PrdfInputNodeNameSet.insert(name);
  if (!result.second)
  {
    std::cout << "EventCombiner::AddPrdfInputNodeName: Prdf Input Node name "
              << name << " already in list - that will wreak havoc and has to be fixed" << std::endl;
    std::cout << "exiting now" << std::endl;
    gSystem->Exit(1);
  }
  else
  {
    if (Verbosity() > 0)
    {
      std::cout << "Prdf node " << name << " inserted" << std::endl;
    }
  }
  return;
}

void EventCombiner::AddPrdfInputNodeFromManager(const Fun4AllInputManager *in)
{
  AddPrdfInputNodeName(in->GetString("PRDFNODENAME"));
  return;
}
