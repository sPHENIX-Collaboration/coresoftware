#include "EventCombiner.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator

#include <Event/Event.h>
#include <Event/oncsEvent.h>

#include <TSystem.h>
//____________________________________________________________________________..
EventCombiner::EventCombiner(const std::string &name):
 SubsysReco(name)
{
  std::cout << "EventCombiner::EventCombiner(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
EventCombiner::~EventCombiner()
{
  std::cout << "EventCombiner::~EventCombiner() Calling dtor" << std::endl;
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
int EventCombiner::InitRun(PHCompositeNode *topNode)
{
  std::cout << "EventCombiner::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventCombiner::process_event(PHCompositeNode *topNode)
{
  std::vector<Event *> subeventeventvec;
  unsigned int total_length = 0;
  for (auto &nam : m_PrdfInputNodeNameSet)
  {
//    std::cout << "name : " << nam << std::endl;
    Event *evt = findNode::getClass<Event>(topNode,nam);
    subeventeventvec.push_back(evt);
    total_length+= evt->getEvtLength();
//    std::cout << nam << " run number " <<  evt->getRunNumber() << " event no: " <<  evt->getEvtSequence() << std::endl;
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
  std::cout << "total length: " <<  total_length << std::endl;
  m_OutArray = new int[total_length];
  subeventeventvec[0]->Copy( m_OutArray , total_length , &nwout);
  current = nwout;
  for (unsigned int icnt = 1; icnt < subeventeventvec.size(); icnt++)
  {
    subeventeventvec[icnt]->Copy(&m_OutArray[current], total_length-current, &nwout, "DATA");
    std::cout << "current " << current << " total length-current: " <<  total_length-current << ", nwout: " << nwout << std::endl;
    current += nwout;
    m_OutArray[0] +=  nwout;
  }  
  PHNodeIterator iter(topNode);
PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfOutputNodeName));

std::cout << "Looking at prdfnode 0x" << std::hex << PrdfNode << std::dec << std::endl;

    m_Event = new oncsEvent(m_OutArray);
PrdfNode->setData(m_Event);
m_Event->identify();
std::cout << "fetching event from tree" << std::endl;
  Event *evts = findNode::getClass<Event>(topNode,m_PrdfOutputNodeName);
  evts->identify();
  std::cout << "evt at " << std::hex << evts << std::dec << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventCombiner::ResetEvent(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfOutputNodeName));
  PrdfNode->setData(nullptr);  // set pointer in Node to nullptr before deleting it
  delete m_Event;
  m_Event = nullptr;
  delete [] m_OutArray;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventCombiner::EndRun(const int runnumber)
{
  std::cout << "EventCombiner::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventCombiner::End(PHCompositeNode *topNode)
{
  std::cout << "EventCombiner::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventCombiner::Reset(PHCompositeNode *topNode)
{
 std::cout << "EventCombiner::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void EventCombiner::Print(const std::string &what) const
{
  std::cout << "EventCombiner::Print(const std::string &what) const Printing info for " << what << std::endl;
}

void EventCombiner::AddPrdfInputNodeName(const std::string &name)
{
  auto result =  m_PrdfInputNodeNameSet.insert(name);
  if (! result.second)
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
