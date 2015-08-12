#include "PHG4ForwardEcalSubsystem.h"
#include "PHG4ForwardEcalDetector.h"
#include "PHG4ForwardEcalSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <fun4all/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;


//_______________________________________________________________________
PHG4ForwardEcalSubsystem::PHG4ForwardEcalSubsystem( const std::string &name, const int lyr ):
  PHG4Subsystem( name ),
  detector_( 0 ),
  steppingAction_( NULL ),
  eventAction_(NULL),
  active(1),
  detector_type(name),
  mappingfile_("")
{

}


//_______________________________________________________________________
int PHG4ForwardEcalSubsystem::Init( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  detector_ = new PHG4ForwardEcalDetector(topNode, Name());
  detector_->SetActive(active);
  detector_->OverlapCheck(overlapcheck);
  detector_->Verbosity(verbosity);
  detector_->SetTowerMappingFile( mappingfile_ );

  if (active)
    {
      // create hit output node
      ostringstream nodename;
      nodename <<  "G4HIT_" << detector_type;

      PHG4HitContainer* scintillator_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
      if (!scintillator_hits)
        {
          scintillator_hits = new PHG4HitContainer();
          PHIODataNode<PHObject> *hitNode = new PHIODataNode<PHObject>(scintillator_hits, nodename.str().c_str(), "PHObject");
          dstNode->addNode(hitNode);
        }

      ostringstream absnodename;
      absnodename << "G4HIT_ABSORBER_" << detector_type;

      PHG4HitContainer* absorber_hits = findNode::getClass<PHG4HitContainer>(topNode, absnodename.str().c_str());
      if (!absorber_hits)
        {
          absorber_hits = new PHG4HitContainer();
          PHIODataNode<PHObject> *abshitNode = new PHIODataNode<PHObject>(absorber_hits, absnodename.str().c_str(), "PHObject");
          dstNode->addNode(abshitNode);
        }

      // create stepping action
      steppingAction_ = new PHG4ForwardEcalSteppingAction(detector_);
    }

  return 0;
}


//_______________________________________________________________________
int
PHG4ForwardEcalSubsystem::process_event( PHCompositeNode * topNode )
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
    {
      steppingAction_->SetInterfacePointers( topNode );
    }
  return 0;
}


//_______________________________________________________________________
PHG4Detector* PHG4ForwardEcalSubsystem::GetDetector( void ) const
{
  return detector_;
}


//_______________________________________________________________________
PHG4SteppingAction* PHG4ForwardEcalSubsystem::GetSteppingAction( void ) const
{
  return steppingAction_;
}
