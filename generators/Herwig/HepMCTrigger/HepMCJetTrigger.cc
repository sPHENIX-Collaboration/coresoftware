//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in HepMCJetTrigger.h.
//
// HepMCJetTrigger(const std::string &name = "HepMCJetTrigger")
// everything is keyed to HepMCJetTrigger, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// HepMCJetTrigger::~HepMCJetTrigger()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int HepMCJetTrigger::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int HepMCJetTrigger::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int HepMCJetTrigger::process_event(PHCompositeNode *topNode)
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
// int HepMCJetTrigger::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int HepMCJetTrigger::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int HepMCJetTrigger::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int HepMCJetTrigger::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void HepMCJetTrigger::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

#include "HepMCJetTrigger.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

//____________________________________________________________________________..
HepMCJetTrigger::HepMCJetTrigger(float trigger_thresh, int n_incom, bool up_lim, const std::string &name):
 SubsysReco(name)
{
  std::cout << "HepMCJetTrigger::HepMCJetTrigger(const std::string &name) Calling ctor" << std::endl;
	this->threshold=trigger_thresh;
	this->goal_event_number=n_incom;
	this->set_event_limit=up_lim;
}


//____________________________________________________________________________..
HepMCJetTrigger::~HepMCJetTrigger()
{
  std::cout << "HepMCJetTrigger::~HepMCJetTrigger() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int HepMCJetTrigger::Init(PHCompositeNode *topNode)
{
  std::cout << "HepMCJetTrigger::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HepMCJetTrigger::InitRun(PHCompositeNode *topNode)
{
  std::cout << "HepMCJetTrigger::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HepMCJetTrigger::process_event(PHCompositeNode *topNode)
{
  //std::cout << "HepMCJetTrigger::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
	n_evts++;
  	if(this->set_event_limit == true){ //needed to keep all HepMC output at the same number of events
		if(n_good >= this->goal_event_number) return Fun4AllReturnCodes::ABORTEVENT;
	}
	PHHepMCGenEventMap* phg=findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
	if(!phg){
		return Fun4AllReturnCodes::ABORTEVENT;
	}
	for( PHHepMCGenEventMap::ConstIter eventIter=phg->begin(); eventIter != phg->end(); ++eventIter){
		PHHepMCGenEvent* hepev=eventIter->second;
		if(!hepev) return Fun4AllReturnCodes::ABORTEVENT;
		HepMC::GenEvent* ev=hepev->getEvent();
		if(!ev) return Fun4AllReturnCodes::ABORTEVENT;
		bool good_event=isGoodEvent(ev);
		if (good_event) n_good++;
		if(!good_event) return Fun4AllReturnCodes::ABORTEVENT;
	}		
	return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HepMCJetTrigger::ResetEvent(PHCompositeNode *topNode)
{
  std::cout << "HepMCJetTrigger::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HepMCJetTrigger::EndRun(const int runnumber)
{
  std::cout << "HepMCJetTrigger::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HepMCJetTrigger::End(PHCompositeNode *topNode)
{
  std::cout << "HepMCJetTrigger::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HepMCJetTrigger::Reset(PHCompositeNode *topNode)
{
 std::cout << "HepMCJetTrigger::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void HepMCJetTrigger::Print(const std::string &what) const
{
  std::cout << "HepMCJetTrigger::Print(const std::string &what) const Printing info for " << what << std::endl;
}
bool HepMCJetTrigger::isGoodEvent( HepMC::GenEvent* e1)
{
	//this is really just the call to actually evaluate and return the filter
	if(this->threshold == 0 ) return true;
	std::vector<fastjet::PseudoJet> jets=findAllJets(e1);
	int njetsabove=jetsAboveThreshold(jets);
	if(njetsabove > 0 ) return true;
	else return false;
}

std::vector<fastjet::PseudoJet> HepMCJetTrigger::findAllJets(HepMC::GenEvent* e1)
{
	//do the fast jet clustering, antikt r=-0.4
	fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, 0.4);
	std::vector<fastjet::PseudoJet> input, output;
	for(HepMC::GenEvent::particle_const_iterator iter = e1->particles_begin(); iter != e1->particles_end(); ++iter)
	{
		if(!(*iter)->end_vertex() && (*iter)->status() == 1){
			auto p=(*iter)->momentum();
			fastjet::PseudoJet pj( p.px(), p.py(), p.pz(), p.e());
			pj.set_user_index((*iter)->barcode());
			input.push_back(pj);
		}
	}
	if(input.size() == 0 ) return input;
	fastjet::ClusterSequence js (input, jetdef);
	auto j = js.inclusive_jets();
	for(auto j1:j){
		output.push_back(j1); //just keep in the corect format
	}
	return output;
}
int HepMCJetTrigger::jetsAboveThreshold(std::vector<fastjet::PseudoJet> jets)
{
	//search through for the number of identified jets above the threshold
	int n_good_jets=0;
	for(auto j:jets)
	{
		float pt=j.pt();
		if(pt > this->threshold) n_good_jets++;
	}
	return n_good_jets;
}
			
