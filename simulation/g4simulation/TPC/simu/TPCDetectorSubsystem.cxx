#include <iostream>
#include <string>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4Utils.h>
#include <phool/getClass.h>
#include <Geant4/globals.hh>
#include <TPCbase/TPCHit.h>
#include <TPCbase/TPCHitsContainer.h>
#include "TPCDetector.h"
#include "TPCDetectorSubsystem.h"
#include "TPCSteppingAction.h"
#include "TPCEventAction.h"


#include <sstream>

//=====
TPCDetectorSubsystem::TPCDetectorSubsystem(const std::string&, const int):
  PHG4DetectorSubsystem("TPCDetectorSubsystem"),
  fDetector( NULL ),
  fEventAction( NULL ),
  fSteppingAction( NULL )
{
  //InitializeParameters();
  if(verbosity>0) std::cout << "TPCDetectorSubsystem::ctor" << std::endl;
}
//=====
int TPCDetectorSubsystem::InitRunSubsystem( PHCompositeNode* node )
{
  if(verbosity>0) std::cout << "TPCDetectorSubsystem::InitRunDetectorSubsystem" << std::endl;
  fDetector = new TPCDetector(node);
  fEventAction = new TPCEventAction(node);
  fSteppingAction = new TPCSteppingAction(fDetector);
  fDetector->OverlapCheck(CheckOverlap());
  fDetector->SetVerbosity(verbosity);
  fEventAction->SetVerbosity(verbosity);
  fSteppingAction->SetVerbosity(verbosity);

  PHNodeIterator iter( node );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST" ));
  TPCHitsContainer *hits =  findNode::getClass<TPCHitsContainer>( node , "TPCHits" );
  if( !hits ) {
    hits = new TPCHitsContainer();
    dstNode->addNode( new PHIODataNode<TPCHitsContainer>( hits, "TPCHits", "TPCHitsContainer" ));
  }
  return 0;
}
//=====
int TPCDetectorSubsystem::process_event( PHCompositeNode* node ) {
  if(verbosity>0) std::cout << "TPCDetectorSubsystem::process_event" << std::endl;
  if(fSteppingAction) fSteppingAction->SetInterfacePointers( node );
  return 0;
}
//=====
void TPCDetectorSubsystem::Print() const {
  std::cout << Name() << "TPCDetectorSubsystem::Print()" << std::endl;
  //if (!BeginRunExecuted()) {
  //}
  return;
}
