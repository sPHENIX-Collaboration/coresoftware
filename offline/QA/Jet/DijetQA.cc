//////////////////////////////////////////////////////////////////
//							    	//
//		Dijet QA module for JetQA		    	//
//		QA module that looks at dijet Ajj, xj		//
//								//
//								//
//								//
//		Author:Skadi 				  	//
//		First commit:  		13 Nov 24		//
//		Mosr recent update:	13 Nov 24		//
//////////////////////////////////////////////////////////////////
#include "DijetQA.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

//____________________________________________________________________________..
DijetQA::DijetQA(const std::string &name):
 SubsysReco(name)
{
  std::cout << "DijetQA::DijetQA(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
DijetQA::~DijetQA()
{
  std::cout << "DijetQA::~DijetQA() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int DijetQA::Init(PHCompositeNode *topNode)
{
  std::cout << "DijetQA::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetQA::InitRun(PHCompositeNode *topNode)
{
  std::cout << "DijetQA::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetQA::process_event(PHCompositeNode *topNode)
{
  std::cout << "DijetQA::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetQA::ResetEvent(PHCompositeNode *topNode)
{
  std::cout << "DijetQA::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetQA::EndRun(const int runnumber)
{
  std::cout << "DijetQA::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetQA::End(PHCompositeNode *topNode)
{
  std::cout << "DijetQA::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetQA::Reset(PHCompositeNode *topNode)
{
 std::cout << "DijetQA::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void DijetQA::Print(const std::string &what) const
{
  std::cout << "DijetQA::Print(const std::string &what) const Printing info for " << what << std::endl;
}
