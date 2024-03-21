//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in LocalPolOffMon.h.
//
// LocalPolOffMon(const std::string &name = "LocalPolOffMon")
// everything is keyed to LocalPolOffMon, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// LocalPolOffMon::~LocalPolOffMon()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int LocalPolOffMon::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int LocalPolOffMon::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int LocalPolOffMon::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int LocalPolOffMon::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int LocalPolOffMon::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int LocalPolOffMon::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int LocalPolOffMon::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void LocalPolOffMon::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

#include "LocalPolOffMon.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

//____________________________________________________________________________..
LocalPolOffMon::LocalPolOffMon(const std::string &name):
 SubsysReco(name)
{
  std::cout << "LocalPolOffMon::LocalPolOffMon(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
LocalPolOffMon::~LocalPolOffMon()
{
  std::cout << "LocalPolOffMon::~LocalPolOffMon() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int LocalPolOffMon::Init(PHCompositeNode * /*topNode*/)
{
  std::cout << "LocalPolOffMon::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LocalPolOffMon::InitRun(PHCompositeNode * /*topNode*/)
{
  std::cout << "LocalPolOffMon::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LocalPolOffMon::process_event(PHCompositeNode * /*topNode*/)
{
  std::cout << "LocalPolOffMon::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LocalPolOffMon::ResetEvent(PHCompositeNode * /*topNode*/)
{
  std::cout << "LocalPolOffMon::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LocalPolOffMon::EndRun(const int runnumber)
{
  std::cout << "LocalPolOffMon::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LocalPolOffMon::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "LocalPolOffMon::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LocalPolOffMon::Reset(PHCompositeNode * /*topNode*/)
{
 std::cout << "LocalPolOffMon::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void LocalPolOffMon::Print(const std::string &what) const
{
  std::cout << "LocalPolOffMon::Print(const std::string &what) const Printing info for " << what << std::endl;
}
