#include "PHG4TruthSubsystem.h"
#include "PHG4TruthEventAction.h"
#include "PHG4TruthSteppingAction.h"
#include "PHG4TruthTrackingAction.h"


#include "PHG4TruthInfoContainer.h"

#include "PHG4VtxPointv1.h"
#include "PHG4Particlev2.h"

#include "PHG4InEvent.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4SystemOfUnits.hh>

#include <cassert>
#include <iostream>

using namespace std;

//_______________________________________________________________________
PHG4TruthSubsystem::PHG4TruthSubsystem( const string &name ):
  PHG4Subsystem( name ),
  eventAction_( NULL ),
  steppingAction_( NULL ),
  trackingAction_( NULL ),
  saveOnlyEmbeded_(false)
{}

//_______________________________________________________________________
int PHG4TruthSubsystem::InitRun( PHCompositeNode* topNode )
{

  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create truth information container
  PHG4TruthInfoContainer* truthInfoList =  findNode::getClass<PHG4TruthInfoContainer>( topNode , "G4TruthInfo" );
  if ( !truthInfoList )
    {
      truthInfoList = new PHG4TruthInfoContainer();
      dstNode->addNode( new PHIODataNode<PHObject>( truthInfoList, "G4TruthInfo", "PHObject" ));
    }

  // event action
  eventAction_ = new PHG4TruthEventAction();

  // create stepping action
  //steppingAction_ = new PHG4TruthSteppingAction( eventAction_ );

  // create tracking action
  trackingAction_ = new PHG4TruthTrackingAction( eventAction_ );

  return 0;
}

//_______________________________________________________________________
int
PHG4TruthSubsystem::process_event( PHCompositeNode* topNode )
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if ( eventAction_ )
    {
      eventAction_->SetInterfacePointers( topNode );
    }
  else
    {
      cout << PHWHERE << " No EventAction registered" << endl;
      exit(1);
    }

  if ( trackingAction_ )
    {
      trackingAction_->SetInterfacePointers( topNode );
    }
  else
    {
      cout << PHWHERE << " No TrackingAction registered" << endl;
      exit(1);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TruthSubsystem::process_after_geant(PHCompositeNode * topNode)
{
  if (saveOnlyEmbeded_)
    {
      if (Verbosity()>1)
        {
          cout <<__PRETTY_FUNCTION__<<" - INFO - only save the G4 truth information that is associated with the embedded particle"<<endl;
        }

      PHG4TruthInfoContainer* truthInfoList =  findNode::getClass<PHG4TruthInfoContainer>( topNode , "G4TruthInfo" );
      assert(truthInfoList);

      set<int> savevtxlist;

      // remove particle that is not embedd associated
      PHG4TruthInfoContainer::Range truth_range = truthInfoList->GetParticleRange();
      PHG4TruthInfoContainer::Iterator truthiter = truth_range.first;
      while (truthiter != truth_range.second)
        {
          const int primary_id = (truthiter->second)->get_primary_id();
          if (truthInfoList->isEmbeded(primary_id) <= 0)
            {
              // not a embed associated particle

              truthInfoList->delete_particle(truthiter++);
            }
          else
            {
              // save vertex id for primary particle which leaves no hit
              // in active area
              savevtxlist.insert((truthiter->second)->get_vtx_id());
              ++truthiter;
            }
        }

      // remove vertex that is not embedd associated
      PHG4TruthInfoContainer::VtxRange vtxrange = truthInfoList->GetVtxRange();
      PHG4TruthInfoContainer::VtxIterator vtxiter = vtxrange.first;
      while (vtxiter != vtxrange.second)
        {
          if (savevtxlist.find(vtxiter->first) == savevtxlist.end())
            {
              truthInfoList->delete_vtx(vtxiter++);
            }
          else
            {
              ++vtxiter;
            }
        }
    }

  return 0;
}

int
PHG4TruthSubsystem::ResetEvent(PHCompositeNode *topNode)
{
  trackingAction_->ResetEvent(topNode);
  eventAction_->ResetEvent(topNode);
  return 0;
}

//_______________________________________________________________________
PHG4EventAction* PHG4TruthSubsystem::GetEventAction( void ) const
{ return eventAction_; }


//_______________________________________________________________________
PHG4SteppingAction* PHG4TruthSubsystem::GetSteppingAction( void ) const
{ return steppingAction_; }

PHG4TrackingAction* 
PHG4TruthSubsystem::GetTrackingAction( void ) const
{ 
  return trackingAction_; 
}
