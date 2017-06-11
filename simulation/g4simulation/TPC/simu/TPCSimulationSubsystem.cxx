#include <iostream>
#include <TH1F.h>
#include <TStopwatch.h>
#include <TPCbase/TPCHitsContainer.h>
#include <TPCbase/TPCDigitsContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>

#include "TPCSimulationSubsystem.h"
#include "TPCSimulation.h"

TPCSimulationSubsystem::TPCSimulationSubsystem() :
  SubsysReco("TPCSimulationSubsystem"),
  fSimulation(new TPCSimulation()),
  fStopwatch(NULL),
  fHist(NULL),
  fTreeFileName("")
{
}

TPCSimulationSubsystem::~TPCSimulationSubsystem()
{
  if(fStopwatch) delete fStopwatch;
  delete fSimulation;
}
//======
int TPCSimulationSubsystem::Init(PHCompositeNode* node)
{
  if(verbosity>0) {
    fStopwatch = new TStopwatch();
    Fun4AllServer *se = Fun4AllServer::instance();
    fHist = new TH1F("TPCSIMTIME","TPCSIMTIME;sec/event",100,0,500);
    se->registerHisto( fHist ); // fHist now owned by me anymore
  }
  if(fTreeFileName.Length()>0) {
    fSimulation->PrepareTree(fTreeFileName);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//======
int TPCSimulationSubsystem::InitRun(PHCompositeNode *node)
{
  fSimulation->SetVerbosity(verbosity);
  PHNodeIterator iter(node);
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if(!dstNode){
    std::cout << PHWHERE << "ERROR!" << std::endl;
    exit(1);
  }
  TPCHitsContainer *hits = findNode::getClass<TPCHitsContainer>( node, "TPCHits" );
  if(!hits) {
    std::cout << PHWHERE << "ERROR!" << std::endl;
    exit(1);
  }
  fSimulation->ConnectHits( hits );
  TPCDigitsContainer *digits = findNode::getClass<TPCDigitsContainer>( node , "TPCDigits" );
  if( !digits ) {
    digits = new TPCDigitsContainer();
    dstNode->addNode( new PHIODataNode<TPCDigitsContainer>( digits, "TPCDigits", "TPCDigitsContainer" ));
  }
  fSimulation->ConnectDigits( digits );
  return Fun4AllReturnCodes::EVENT_OK;
}
//======
int TPCSimulationSubsystem::process_event(PHCompositeNode *topNode)
{
  if(verbosity>1) {
    fStopwatch->Reset();
    fStopwatch->Start();
  }
  //-----
  fSimulation->Hits2Digits();
  //-----
  if(verbosity>0) {
    fHist->Fill( fStopwatch->RealTime() );
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//======
int TPCSimulationSubsystem::End(PHCompositeNode *topNode)
{
  if(fTreeFileName.Length()>0) {
    fSimulation->WriteFile();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}


