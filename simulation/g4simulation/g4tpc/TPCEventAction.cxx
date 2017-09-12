#include <iostream>
#include <phool/getClass.h>
#include <Geant4/G4Event.hh>
#include <TPCHitsContainer.h>
#include "TPCEventAction.h"

//=====
TPCEventAction::TPCEventAction( PHCompositeNode *node ):
  fVerbosity(0),
  fNode(node)
{
}
//=====
void TPCEventAction::BeginOfEventAction(const G4Event* evt)
{
  TPCHitsContainer* hits = findNode::getClass<TPCHitsContainer>( fNode, "TPCHits" );
  if(hits) {
    //std::cout << "TPCEventAction::BeginOfEventAction => Reset containers" << std::endl;
    hits->Reset();
  }
}
//=====
void TPCEventAction::EndOfEventAction(const G4Event* evt)
{
}
